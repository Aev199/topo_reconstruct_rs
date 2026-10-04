//! Planar surfaces on a shared topological graph. No implicit mechanical ties.
pub mod assembly;
pub mod axes;
pub mod frame;
pub mod graph;
pub mod mesh;
pub mod planes;
pub mod recognize;
pub mod reconcile;

use geo::{Contains, Intersects, Line, LineString, Polygon};
use glam::DVec3;
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

#[derive(Debug, Clone, PartialEq)]
pub enum Error {
    InvalidPrecision,
    InvalidPlane,
    InvalidVertex,
    InvalidRing,
    NonPlanar,
    ShortEdge,
}

#[derive(Debug, Clone, Serialize, serde::Deserialize)]
pub struct PlaneFrame {
    origin: [f64; 3],
    normal: [f64; 3],
    u: [f64; 3],
    v: [f64; 3],
}
impl PlaneFrame {
    pub fn new(origin: [f64; 3], normal: [f64; 3]) -> Result<Self, Error> {
        let o = DVec3::from_array(origin);
        let n = DVec3::from_array(normal);
        if !o.is_finite() || !n.is_finite() || !n.length().is_finite() || n.length() == 0.0 {
            return Err(Error::InvalidPlane);
        }
        let n = n.normalize();
        let helper = if n.z.abs() < 0.9 { DVec3::Z } else { DVec3::X };
        let u = helper.cross(n).normalize();
        Ok(Self {
            origin,
            normal: n.to_array(),
            u: u.to_array(),
            v: n.cross(u).to_array(),
        })
    }
    pub fn origin(&self) -> [f64; 3] {
        self.origin
    }
    /// Unit normal.
    pub fn normal(&self) -> [f64; 3] {
        self.normal
    }
    pub fn distance(&self, point: [f64; 3]) -> f64 {
        (DVec3::from_array(point) - DVec3::from_array(self.origin))
            .dot(DVec3::from_array(self.normal))
    }
    pub fn project(&self, point: [f64; 3]) -> [f64; 2] {
        let d = DVec3::from_array(point) - DVec3::from_array(self.origin);
        [
            d.dot(DVec3::from_array(self.u)),
            d.dot(DVec3::from_array(self.v)),
        ]
    }
    pub fn lift(&self, uv: [f64; 2]) -> [f64; 3] {
        (DVec3::from_array(self.origin)
            + uv[0] * DVec3::from_array(self.u)
            + uv[1] * DVec3::from_array(self.v))
        .to_array()
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct EdgeUse {
    pub edge: usize,
    pub reversed: bool,
}
#[derive(Debug, Clone, Serialize)]
pub struct Surface {
    pub plane: usize,
    /// Exterior first, followed by holes. Geometry lives in plane coordinates.
    pub contours: Vec<Vec<[f64; 2]>>,
    pub boundaries: Vec<Vec<EdgeUse>>,
    /// Model edges lying inside the material (not on a contour): explicit
    /// junction lines shared with other surfaces. They impose topology on
    /// the surface and its mesh without dividing its property region.
    pub embedded_edges: Vec<usize>,
    pub source_elements: Vec<u32>,
}
/// Trial copies (every transactional operation works on one) are cheap:
/// vertices, edges, the edge index and each surface are shared until a
/// copy changes them (`Arc::make_mut`).
#[derive(Debug, Clone, Serialize)]
pub struct Model {
    precision: f64,
    minimum_edge: f64,
    planes: Vec<PlaneFrame>,
    vertices: Arc<Vec<[f64; 3]>>,
    edges: Arc<Vec<[usize; 2]>>,
    surfaces: Vec<Arc<Surface>>,
    /// Edges referenced by no surface (replaced during a split, released by
    /// a trim or merge). Kept, still indexed, to preserve edge ids; refreshed
    /// by `refresh_orphaned_edges`.
    orphaned_edges: Vec<usize>,
    #[serde(skip)]
    edge_index: Arc<BTreeMap<[usize; 2], usize>>,
}
impl Model {
    pub fn new(precision: f64, minimum_edge: f64) -> Result<Self, Error> {
        if !precision.is_finite()
            || !minimum_edge.is_finite()
            || precision <= 0.0
            || minimum_edge <= precision
        {
            return Err(Error::InvalidPrecision);
        }
        Ok(Self {
            precision,
            minimum_edge,
            planes: vec![],
            vertices: vec![].into(),
            edges: vec![].into(),
            surfaces: vec![],
            orphaned_edges: vec![],
            edge_index: BTreeMap::new().into(),
        })
    }
    pub fn add_plane(&mut self, plane: PlaneFrame) -> usize {
        let id = self.planes.len();
        self.planes.push(plane);
        id
    }
    /// Identity is explicit: coordinate proximity alone never creates a connection.
    pub fn add_vertex(&mut self, point: [f64; 3]) -> Result<usize, Error> {
        if !DVec3::from_array(point).is_finite() {
            return Err(Error::InvalidVertex);
        }
        let id = self.vertices.len();
        Arc::make_mut(&mut self.vertices).push(point);
        Ok(id)
    }
    pub fn surfaces(&self) -> &[Arc<Surface>] {
        &self.surfaces
    }
    pub fn edges(&self) -> &[[usize; 2]] {
        &self.edges
    }
    pub fn planes(&self) -> &[PlaneFrame] {
        &self.planes
    }

    /// Transactional insertion. Reject warped input rather than projecting a
    /// shared vertex independently onto each owner's plane.
    pub fn add_surface(
        &mut self,
        plane: usize,
        rings: Vec<Vec<usize>>,
        mut source_elements: Vec<u32>,
    ) -> Result<usize, Error> {
        let contours = self.contours(plane, &rings)?;
        // All validation precedes mutation, including edge interning.
        let boundaries = self.intern(&rings);
        source_elements.sort_unstable();
        source_elements.dedup();
        let id = self.surfaces.len();
        self.surfaces.push(Arc::new(Surface {
            plane,
            contours,
            boundaries,
            embedded_edges: vec![],
            source_elements,
        }));
        Ok(id)
    }

    /// Validated plane contours of vertex rings (exterior first).
    fn contours(&self, plane: usize, rings: &[Vec<usize>]) -> Result<Vec<Vec<[f64; 2]>>, Error> {
        self.contours_changed(plane, rings, &vec![true; rings.len()])
    }

    /// [`Model::contours`] of a valid surface of which only the rings
    /// flagged in `changed` were edited: unchanged rings and the relations
    /// between them were valid and are not tested again.
    fn contours_changed(
        &self,
        plane: usize,
        rings: &[Vec<usize>],
        changed: &[bool],
    ) -> Result<Vec<Vec<[f64; 2]>>, Error> {
        let frame = self.planes.get(plane).ok_or(Error::InvalidPlane)?;
        if rings.is_empty() {
            return Err(Error::InvalidRing);
        }
        let mut contours = vec![];
        for ring in rings {
            if ring.len() < 3 || ring.iter().copied().collect::<BTreeSet<_>>().len() != ring.len() {
                return Err(Error::InvalidRing);
            }
            let mut uv = vec![];
            for (i, &id) in ring.iter().enumerate() {
                let p = *self.vertices.get(id).ok_or(Error::InvalidVertex)?;
                let q = *self
                    .vertices
                    .get(ring[(i + 1) % ring.len()])
                    .ok_or(Error::InvalidVertex)?;
                if frame.distance(p).abs() > self.precision {
                    return Err(Error::NonPlanar);
                }
                let a = frame.project(p);
                let b = frame.project(q);
                let length = (a[0] - b[0]).hypot(a[1] - b[1]);
                if !length.is_finite() || length < self.minimum_edge {
                    return Err(Error::ShortEdge);
                }
                uv.push(a);
            }
            if changed[contours.len()] {
                validate_ring(&uv, self.precision)?;
            }
            contours.push(uv);
        }
        validate_holes(&contours, changed)?;
        Ok(contours)
    }

    fn intern(&mut self, rings: &[Vec<usize>]) -> Vec<Vec<EdgeUse>> {
        rings
            .iter()
            .map(|ring| {
                (0..ring.len())
                    .map(|i| {
                        let a = ring[i];
                        let b = ring[(i + 1) % ring.len()];
                        let key = [a.min(b), a.max(b)];
                        let edge = *Arc::make_mut(&mut self.edge_index)
                            .entry(key)
                            .or_insert_with(|| {
                                let id = self.edges.len();
                                Arc::make_mut(&mut self.edges).push(key);
                                id
                            });
                        EdgeUse {
                            edge,
                            reversed: a > b,
                        }
                    })
                    .collect()
            })
            .collect()
    }

    /// Merge vertex `drop` into `keep`: every edge and contour using `drop`
    /// uses `keep` instead. `keep` must lie on the plane of every surface
    /// using either vertex. Edges collapsing to a point or duplicating
    /// another edge are merged; each affected surface is revalidated. The
    /// model is unchanged on error.
    pub fn merge_vertices(&mut self, drop: usize, keep: usize) -> Result<(), Error> {
        if drop == keep || drop >= self.vertices.len() || keep >= self.vertices.len() {
            return Err(Error::InvalidVertex);
        }
        let target = self.vertices[keep];
        let users: Vec<usize> = (0..self.surfaces.len())
            .filter(|&s| self.surface_edges(s).any(|e| self.edges[e].contains(&drop)))
            .collect();
        let mut trial = self.clone();
        for &s in &users {
            let surface = &trial.surfaces[s];
            if trial.planes[surface.plane].distance(target).abs() > trial.precision {
                return Err(Error::NonPlanar);
            }
            let map = |v: usize| if v == drop { keep } else { v };
            let mut rings = vec![];
            for ring in &surface.boundaries {
                let mut ids: Vec<usize> = ring
                    .iter()
                    .map(|e| {
                        let [a, b] = trial.edges[e.edge];
                        map(if e.reversed { b } else { a })
                    })
                    .collect();
                ids.dedup();
                if ids.len() > 1 && ids.first() == ids.last() {
                    ids.pop();
                }
                rings.push(ids);
            }
            let keys: BTreeSet<[usize; 2]> = surface
                .embedded_edges
                .iter()
                .map(|&e| trial.edges[e].map(map))
                .filter(|[a, b]| a != b)
                .map(|[a, b]| [a.min(b), a.max(b)])
                .collect();
            let ring_keys: BTreeSet<[usize; 2]> = rings
                .iter()
                .flat_map(|r| {
                    (0..r.len()).map(move |i| {
                        let (a, b) = (r[i], r[(i + 1) % r.len()]);
                        [a.min(b), a.max(b)]
                    })
                })
                .collect();
            let embedded: Vec<usize> = keys
                .difference(&ring_keys)
                .map(|&key| {
                    *Arc::make_mut(&mut trial.edge_index)
                        .entry(key)
                        .or_insert_with(|| {
                            Arc::make_mut(&mut trial.edges).push(key);
                            trial.edges.len() - 1
                        })
                })
                .collect();
            trial.rebuild_surface(s, rings, embedded)?;
        }
        *self = trial;
        Ok(())
    }

    /// Remove a vertex joining exactly two collinear edges, replacing them
    /// with one edge in every surface that uses them. The vertex must lie on
    /// the joined segment within `max(precision, deviation)`; contours change
    /// by at most that distance. Fails if the vertex has another edge, the
    /// edges have different users, or a resulting contour is invalid. The
    /// model is unchanged on error.
    pub fn remove_vertex(&mut self, vertex: usize, deviation: f64) -> Result<usize, Error> {
        let all: Vec<usize> = (0..self.surfaces.len()).collect();
        self.remove_vertex_among(vertex, deviation, &all)
    }

    /// [`Model::remove_vertex`] when the caller knows that no surface
    /// outside `candidates` uses the vertex (the users of a contour chain):
    /// only those surfaces are searched.
    pub fn remove_vertex_among(
        &mut self,
        vertex: usize,
        deviation: f64,
        candidates: &[usize],
    ) -> Result<usize, Error> {
        let mut incident: Vec<usize> = candidates
            .iter()
            .flat_map(|&s| {
                self.surface_edges(s)
                    .filter(|&e| self.edges[e].contains(&vertex))
            })
            .collect();
        incident.sort_unstable();
        incident.dedup();
        let [e1, e2] = incident[..] else {
            return Err(Error::InvalidVertex);
        };
        let other = |e: usize| {
            let [a, b] = self.edges[e];
            if a == vertex {
                b
            } else {
                a
            }
        };
        let (p, q) = (other(e1), other(e2));
        let users = |e: usize| -> Vec<usize> {
            let mut list: Vec<usize> = candidates
                .iter()
                .copied()
                .filter(|&s| self.surface_edges(s).any(|x| x == e))
                .collect();
            list.sort_unstable();
            list.dedup();
            list
        };
        let surfaces = users(e1);
        if surfaces != users(e2) || p == q {
            return Err(Error::InvalidVertex);
        }
        let (pp, pq, pv) = (
            DVec3::from_array(self.vertices[p]),
            DVec3::from_array(self.vertices[q]),
            DVec3::from_array(self.vertices[vertex]),
        );
        let d = pq - pp;
        let t = (pv - pp).dot(d) / d.length_squared();
        if !(t > 0. && t < 1.) || pv.distance(pp + d * t) > self.precision.max(deviation) {
            return Err(Error::NonPlanar);
        }
        for &s in &surfaces {
            let embedded = &self.surfaces[s].embedded_edges;
            if embedded.contains(&e1) != embedded.contains(&e2) {
                return Err(Error::InvalidRing);
            }
        }
        let key = [p.min(q), p.max(q)];
        if let Some(&existing) = self.edge_index.get(&key) {
            if surfaces
                .iter()
                .any(|&s| self.surface_edges(s).any(|x| x == existing))
            {
                return Err(Error::InvalidRing);
            }
        }
        // Validate all new contours before mutation.
        let mut updates = vec![];
        for &s in &surfaces {
            let surface = &self.surfaces[s];
            let mut rings = vec![];
            let mut changed = vec![];
            for ring in &surface.boundaries {
                let mut this = false;
                let ids: Vec<usize> = ring
                    .iter()
                    .map(|e| {
                        let [a, b] = self.edges[e.edge];
                        if e.reversed {
                            b
                        } else {
                            a
                        }
                    })
                    .filter(|&v| {
                        let keep = v != vertex;
                        this |= !keep;
                        keep
                    })
                    .collect();
                rings.push(ids);
                changed.push(this);
            }
            let ring_changed = changed.iter().any(|&c| c);
            let contours = self.contours_changed(surface.plane, &rings, &changed)?;
            updates.push((s, rings, contours, ring_changed));
        }
        let merged = *Arc::make_mut(&mut self.edge_index)
            .entry(key)
            .or_insert_with(|| {
                Arc::make_mut(&mut self.edges).push(key);
                self.edges.len() - 1
            });
        for (s, rings, contours, ring_changed) in updates {
            if ring_changed {
                let boundaries = self.intern(&rings);
                let surface = Arc::make_mut(&mut self.surfaces[s]);
                surface.boundaries = boundaries;
                surface.contours = contours;
            }
            let surface = Arc::make_mut(&mut self.surfaces[s]);
            if surface.embedded_edges.contains(&e1) {
                surface.embedded_edges.retain(|&e| e != e1 && e != e2);
                surface.embedded_edges.push(merged);
            }
        }
        Ok(merged)
    }

    /// Replace the material domain of a surface: new rings (exterior first)
    /// and embedded edges, validated like a new surface. Plane and source
    /// provenance are kept. Edges no longer used by any surface are recorded
    /// as unreferenced. The model is unchanged on error.
    pub fn rebuild_surface(
        &mut self,
        surface: usize,
        rings: Vec<Vec<usize>>,
        embedded: Vec<usize>,
    ) -> Result<(), Error> {
        let plane = self.surfaces.get(surface).ok_or(Error::InvalidRing)?.plane;
        let contours = self.contours(plane, &rings)?;
        let frame = &self.planes[plane];
        let ring_keys: BTreeSet<[usize; 2]> = rings
            .iter()
            .flat_map(|r| {
                (0..r.len()).map(move |i| {
                    let (a, b) = (r[i], r[(i + 1) % r.len()]);
                    [a.min(b), a.max(b)]
                })
            })
            .collect();
        for &e in &embedded {
            let key = *self.edges.get(e).ok_or(Error::InvalidRing)?;
            if ring_keys.contains(&key) {
                return Err(Error::InvalidRing);
            }
            let [pa, pb] = key.map(|v| self.vertices[v]);
            if frame.distance(pa).abs() > self.precision
                || frame.distance(pb).abs() > self.precision
            {
                return Err(Error::NonPlanar);
            }
            let (ua, ub) = (frame.project(pa), frame.project(pb));
            let mid = [(ua[0] + ub[0]) / 2., (ua[1] + ub[1]) / 2.];
            if !closed_contains(&contours, mid, self.precision)
                || properly_crosses(&contours, ua, ub, self.precision)
            {
                return Err(Error::InvalidRing);
            }
        }
        let boundaries = self.intern(&rings);
        let target = Arc::make_mut(&mut self.surfaces[surface]);
        target.contours = contours;
        target.boundaries = boundaries;
        target.embedded_edges = embedded;
        Ok(())
    }

    /// Join surface `other` into `keep`: both in one plane (every vertex of
    /// `other` on the plane of `keep`), sharing at least one contour edge.
    /// Shared contour edges leave the contours; one still used by a third
    /// surface (a wall standing on the joint) stays as an embedded edge.
    /// Source elements are united. Returns the new index of every old
    /// surface (`other` is removed). The model is unchanged on error.
    pub fn join_surfaces(
        &mut self,
        keep: usize,
        other: usize,
    ) -> Result<Vec<Option<usize>>, Error> {
        if keep == other || keep >= self.surfaces.len() || other >= self.surfaces.len() {
            return Err(Error::InvalidRing);
        }
        let frame = self.planes[self.surfaces[keep].plane].clone();
        let users_of = |m: &Model, e: usize| -> usize {
            (0..m.surfaces.len())
                .filter(|&s| m.surface_edges(s).any(|x| x == e))
                .count()
        };
        for e in self.surface_edges(other) {
            for v in self.edges[e] {
                if frame.distance(self.vertices[v]).abs() > self.precision {
                    return Err(Error::NonPlanar);
                }
            }
        }
        // Directed contour edges of both, `other` oriented like `keep`.
        let directed = |m: &Model, s: usize| -> Vec<Vec<(usize, usize, usize)>> {
            m.surfaces[s]
                .boundaries
                .iter()
                .map(|ring| {
                    ring.iter()
                        .map(|u| {
                            let [a, b] = m.edges[u.edge];
                            if u.reversed {
                                (b, a, u.edge)
                            } else {
                                (a, b, u.edge)
                            }
                        })
                        .collect()
                })
                .collect()
        };
        let signed = |ring: &[(usize, usize, usize)]| -> f64 {
            ring.iter()
                .map(|&(a, b, _)| {
                    let (p, q) = (
                        frame.project(self.vertices[a]),
                        frame.project(self.vertices[b]),
                    );
                    p[0] * q[1] - p[1] * q[0]
                })
                .sum::<f64>()
        };
        let mine = directed(self, keep);
        let mut theirs = directed(self, other);
        if signed(&mine[0]).signum() != signed(&theirs[0]).signum() {
            for ring in &mut theirs {
                for x in ring.iter_mut() {
                    *x = (x.1, x.0, x.2);
                }
            }
        }
        let mut uses: BTreeMap<usize, Vec<(usize, usize)>> = BTreeMap::new();
        for &(a, b, e) in mine.iter().chain(theirs.iter()).flatten() {
            uses.entry(e).or_default().push((a, b));
        }
        let mut shared = vec![];
        let mut out: BTreeMap<usize, Vec<(usize, usize)>> = BTreeMap::new();
        for (e, list) in &uses {
            match list[..] {
                [(a, b)] => out.entry(a).or_default().push((b, *e)),
                [(a, b), (c, d)] if a == d && b == c => shared.push(*e),
                // The same edge in the same direction: overlapping material.
                _ => return Err(Error::InvalidRing),
            }
        }
        if shared.is_empty() {
            return Err(Error::InvalidRing);
        }
        // Relink the remaining edges into rings; a vertex with two ways out
        // (contours pinching) is refused.
        if out.values().any(|v| v.len() != 1) {
            return Err(Error::InvalidRing);
        }
        let mut rings: Vec<Vec<usize>> = vec![];
        let mut seen = BTreeSet::new();
        for &start in out.keys() {
            if seen.contains(&start) {
                continue;
            }
            let mut ring = vec![];
            let mut v = start;
            loop {
                if !seen.insert(v) {
                    break;
                }
                ring.push(v);
                v = out.get(&v).ok_or(Error::InvalidRing)?[0].0;
            }
            if v != start {
                return Err(Error::InvalidRing);
            }
            rings.push(ring);
        }
        let area = |ring: &[usize]| {
            (0..ring.len())
                .map(|i| {
                    let (p, q) = (
                        frame.project(self.vertices[ring[i]]),
                        frame.project(self.vertices[ring[(i + 1) % ring.len()]]),
                    );
                    p[0] * q[1] - p[1] * q[0]
                })
                .sum::<f64>()
                .abs()
        };
        rings.sort_by(|a, b| area(b).total_cmp(&area(a)));
        let mut embedded: BTreeSet<usize> = self.surfaces[keep]
            .embedded_edges
            .iter()
            .chain(self.surfaces[other].embedded_edges.iter())
            .copied()
            .collect();
        for &e in &shared {
            // Used by `keep` and `other`, and by a third surface.
            if users_of(self, e) > 2 {
                embedded.insert(e);
            }
        }
        let mut trial = self.clone();
        let mut sources = trial.surfaces[keep].source_elements.clone();
        sources.extend(trial.surfaces[other].source_elements.iter().copied());
        sources.sort_unstable();
        sources.dedup();
        trial.rebuild_surface(keep, rings, embedded.into_iter().collect())?;
        Arc::make_mut(&mut trial.surfaces[keep]).source_elements = sources;
        let index = trial.remove_surfaces(&BTreeSet::from([other]));
        *self = trial;
        Ok(index)
    }

    /// Remove surfaces (their edges used by no other surface become
    /// orphaned). Returns the new index of every old surface.
    pub fn remove_surfaces(&mut self, remove: &BTreeSet<usize>) -> Vec<Option<usize>> {
        let mut index = vec![];
        let mut kept = vec![];
        for (s, surface) in std::mem::take(&mut self.surfaces).into_iter().enumerate() {
            if remove.contains(&s) {
                index.push(None);
            } else {
                index.push(Some(kept.len()));
                kept.push(surface);
            }
        }
        self.surfaces = kept;
        self.refresh_orphaned_edges();
        index
    }

    /// Recompute the list of edges no surface uses.
    pub fn refresh_orphaned_edges(&mut self) {
        let used: BTreeSet<usize> = (0..self.surfaces.len())
            .flat_map(|s| self.surface_edges(s).collect::<Vec<_>>())
            .collect();
        self.orphaned_edges = (0..self.edges.len())
            .filter(|e| !used.contains(e))
            .collect();
    }
    pub fn vertices(&self) -> &[[f64; 3]] {
        &self.vertices
    }
    pub fn precision(&self) -> f64 {
        self.precision
    }
    pub fn minimum_edge(&self) -> f64 {
        self.minimum_edge
    }
    /// Boundary and embedded edges of a surface.
    pub fn surface_edges(&self, surface: usize) -> impl Iterator<Item = usize> + '_ {
        let s = &self.surfaces[surface];
        s.boundaries
            .iter()
            .flatten()
            .map(|e| e.edge)
            .chain(s.embedded_edges.iter().copied())
    }
    pub fn edge_between(&self, a: usize, b: usize) -> Option<usize> {
        self.edge_index.get(&[a.min(b), a.max(b)]).copied()
    }

    /// Split one global edge at an explicit vertex. Every surface using the
    /// edge, as a boundary or embedded edge, receives the same two pieces, so
    /// shared identity is preserved. Validation precedes mutation.
    pub fn split_edge(&mut self, edge: usize, vertex: usize) -> Result<usize, Error> {
        self.split_edge_within(edge, vertex, 0.)
    }

    /// `split_edge` for a vertex up to `deviation` off the edge line: the
    /// edge bends through the vertex. Every user surface is revalidated.
    pub fn split_edge_within(
        &mut self,
        edge: usize,
        vertex: usize,
        deviation: f64,
    ) -> Result<usize, Error> {
        let [a, b] = *self.edges.get(edge).ok_or(Error::InvalidVertex)?;
        let p = DVec3::from_array(*self.vertices.get(vertex).ok_or(Error::InvalidVertex)?);
        if vertex == a || vertex == b {
            return Err(Error::InvalidVertex);
        }
        let pa = DVec3::from_array(self.vertices[a]);
        let pb = DVec3::from_array(self.vertices[b]);
        let d = pb - pa;
        let t = (p - pa).dot(d) / d.length_squared();
        if !(t > 0. && t < 1.) || p.distance(pa + d * t) > self.precision.max(deviation) {
            return Err(Error::InvalidVertex);
        }
        if p.distance(pa) < self.minimum_edge || p.distance(pb) < self.minimum_edge {
            return Err(Error::ShortEdge);
        }
        // Off the line, the contours change shape: validate them after the
        // split and restore the model on failure.
        if p.distance(pa + d * t) > self.precision {
            let backup = self.clone();
            let result = self.split_unchecked(edge, vertex, a, b, p);
            let valid = result.is_ok()
                && (0..self.surfaces.len())
                    .filter(|&s| {
                        self.surface_edges(s)
                            .any(|e| self.edges[e].contains(&vertex))
                    })
                    .all(|s| {
                        let c = &self.surfaces[s].contours;
                        c.iter().all(|r| validate_ring(r, self.precision).is_ok())
                            && validate_holes(c, &vec![true; c.len()]).is_ok()
                    });
            if !valid {
                *self = backup;
                return result.and(Err(Error::InvalidRing));
            }
            return result;
        }
        self.split_unchecked(edge, vertex, a, b, p)
    }

    /// Replace `edge` (a, b) by two pieces through `vertex` at `p`.
    fn split_unchecked(
        &mut self,
        edge: usize,
        vertex: usize,
        a: usize,
        b: usize,
        p: DVec3,
    ) -> Result<usize, Error> {
        let users: Vec<usize> = (0..self.surfaces.len())
            .filter(|&s| self.surface_edges(s).any(|e| e == edge))
            .collect();
        for &s in &users {
            // A contour through the vertex already would touch itself there.
            let boundary = || self.surfaces[s].boundaries.iter().flatten();
            if boundary().any(|u| u.edge == edge)
                && boundary().any(|u| self.edges[u.edge].contains(&vertex))
            {
                return Err(Error::InvalidRing);
            }
            if self.planes[self.surfaces[s].plane]
                .distance(p.to_array())
                .abs()
                > self.precision
            {
                return Err(Error::NonPlanar);
            }
        }
        let first = [a.min(vertex), a.max(vertex)];
        let second = [b.min(vertex), b.max(vertex)];
        // A piece may already exist, e.g. when the vertex is joined to an end
        // by another surface. Reuse it: one edge per vertex pair. The original
        // edge keeps its id for the first new piece whenever possible.
        let (existing_first, existing_second) = (
            self.edge_index.get(&first).copied(),
            self.edge_index.get(&second).copied(),
        );
        for (piece, users_of) in [(existing_first, &users), (existing_second, &users)] {
            if let Some(piece) = piece {
                if users_of
                    .iter()
                    .any(|&s| self.surface_edges(s).any(|e| e == piece))
                {
                    // The surface would use both the edge and its piece.
                    return Err(Error::InvalidRing);
                }
            }
        }
        // The original key stays indexed unless its id is reused for a piece,
        // so a later edge between the same vertices reuses this id.
        let rekeyed = existing_first.is_none() || existing_second.is_none();
        if rekeyed {
            Arc::make_mut(&mut self.edge_index).remove(&[a, b]);
        }
        let e_first = match existing_first {
            Some(id) => id,
            None => {
                Arc::make_mut(&mut self.edges)[edge] = first;
                Arc::make_mut(&mut self.edge_index).insert(first, edge);
                edge
            }
        };
        let e_second = match existing_second {
            Some(id) => id,
            None if e_first != edge => {
                Arc::make_mut(&mut self.edges)[edge] = second;
                Arc::make_mut(&mut self.edge_index).insert(second, edge);
                edge
            }
            None => {
                Arc::make_mut(&mut self.edges).push(second);
                Arc::make_mut(&mut self.edge_index).insert(second, self.edges.len() - 1);
                self.edges.len() - 1
            }
        };
        for s in users {
            let plane = self.planes[self.surfaces[s].plane].clone();
            let uv = plane.project(p.to_array());
            let surface = Arc::make_mut(&mut self.surfaces[s]);
            for (ring, contour) in surface.boundaries.iter_mut().zip(&mut surface.contours) {
                if let Some(k) = ring.iter().position(|e| e.edge == edge) {
                    let reversed = ring[k].reversed;
                    let (start, end) = if reversed { (b, a) } else { (a, b) };
                    let piece = |from: usize, to: usize, id: usize| EdgeUse {
                        edge: id,
                        reversed: from > to,
                    };
                    let (e0, e1) = if reversed {
                        (e_second, e_first)
                    } else {
                        (e_first, e_second)
                    };
                    ring[k] = piece(start, vertex, e0);
                    ring.insert(k + 1, piece(vertex, end, e1));
                    contour.insert(k + 1, uv);
                }
            }
            if let Some(k) = surface.embedded_edges.iter().position(|&e| e == edge) {
                surface.embedded_edges.remove(k);
                surface.embedded_edges.extend([e_first, e_second]);
            }
        }
        Ok(e_second)
    }

    /// Move a vertex, revalidating every surface that uses it: planarity,
    /// minimum edge length, ring validity, holes and embedded edges. The
    /// model is unchanged on error.
    pub fn move_vertex(&mut self, vertex: usize, target: [f64; 3]) -> Result<(), Error> {
        if vertex >= self.vertices.len() || !DVec3::from_array(target).is_finite() {
            return Err(Error::InvalidVertex);
        }
        let users: Vec<usize> = (0..self.surfaces.len())
            .filter(|&s| {
                self.surface_edges(s)
                    .any(|e| self.edges[e].contains(&vertex))
            })
            .collect();
        let position = |v: usize| {
            if v == vertex {
                target
            } else {
                self.vertices[v]
            }
        };
        let mut updates = vec![];
        for &s in &users {
            let surface = &self.surfaces[s];
            let plane = &self.planes[surface.plane];
            let mut contours = vec![];
            for ring in &surface.boundaries {
                let mut uv = vec![];
                for e in ring {
                    let [a, b] = self.edges[e.edge];
                    let (from, to) = if e.reversed { (b, a) } else { (a, b) };
                    let (p, q) = (position(from), position(to));
                    if plane.distance(p).abs() > self.precision {
                        return Err(Error::NonPlanar);
                    }
                    if DVec3::from_array(p).distance(DVec3::from_array(q)) < self.minimum_edge {
                        return Err(Error::ShortEdge);
                    }
                    uv.push(plane.project(p));
                }
                validate_ring(&uv, self.precision)?;
                contours.push(uv);
            }
            validate_holes(&contours, &vec![true; contours.len()])?;
            for &e in &surface.embedded_edges {
                let [a, b] = self.edges[e];
                let (p, q) = (position(a), position(b));
                if plane.distance(p).abs() > self.precision
                    || plane.distance(q).abs() > self.precision
                {
                    return Err(Error::NonPlanar);
                }
                let (ua, ub) = (plane.project(p), plane.project(q));
                if (ua[0] - ub[0]).hypot(ua[1] - ub[1]) < self.minimum_edge {
                    return Err(Error::ShortEdge);
                }
                let mid = [(ua[0] + ub[0]) / 2., (ua[1] + ub[1]) / 2.];
                if !closed_contains(&contours, mid, self.precision)
                    || properly_crosses(&contours, ua, ub, self.precision)
                {
                    return Err(Error::InvalidRing);
                }
            }
            updates.push((s, contours));
        }
        Arc::make_mut(&mut self.vertices)[vertex] = target;
        for (s, contours) in updates {
            Arc::make_mut(&mut self.surfaces[s]).contours = contours;
        }
        Ok(())
    }

    /// Record an existing or new edge as an interior junction line of a
    /// surface. The edge must lie on the surface plane and inside its closed
    /// material domain; boundary edges are left unchanged.
    pub fn embed_edge(&mut self, surface: usize, a: usize, b: usize) -> Result<usize, Error> {
        let s = self.surfaces.get(surface).ok_or(Error::InvalidRing)?;
        let plane = &self.planes[s.plane];
        let pa = *self.vertices.get(a).ok_or(Error::InvalidVertex)?;
        let pb = *self.vertices.get(b).ok_or(Error::InvalidVertex)?;
        if a == b {
            return Err(Error::InvalidVertex);
        }
        if plane.distance(pa).abs() > self.precision || plane.distance(pb).abs() > self.precision {
            return Err(Error::NonPlanar);
        }
        let (ua, ub) = (plane.project(pa), plane.project(pb));
        if (ua[0] - ub[0]).hypot(ua[1] - ub[1]) < self.minimum_edge {
            return Err(Error::ShortEdge);
        }
        let key = [a.min(b), a.max(b)];
        if let Some(&edge) = self.edge_index.get(&key) {
            if self.surface_edges(surface).any(|e| e == edge) {
                return Ok(edge);
            }
        }
        let mid = [(ua[0] + ub[0]) / 2., (ua[1] + ub[1]) / 2.];
        if !closed_contains(&s.contours, mid, self.precision)
            || properly_crosses(&s.contours, ua, ub, self.precision)
        {
            return Err(Error::InvalidRing);
        }
        let edge = *Arc::make_mut(&mut self.edge_index)
            .entry(key)
            .or_insert_with(|| {
                Arc::make_mut(&mut self.edges).push(key);
                self.edges.len() - 1
            });
        Arc::make_mut(&mut self.surfaces[surface])
            .embedded_edges
            .push(edge);
        Ok(edge)
    }
}

/// Whether segment `a`-`b` crosses a contour segment at an interior point
/// of both, beyond `eps`. Touching at a vertex or along a line is not a crossing.
pub(crate) fn properly_crosses(
    contours: &[Vec<[f64; 2]>],
    a: [f64; 2],
    b: [f64; 2],
    eps: f64,
) -> bool {
    let side = |p: [f64; 2], q: [f64; 2], r: [f64; 2]| {
        let (dx, dy) = (q[0] - p[0], q[1] - p[1]);
        let length = dx.hypot(dy);
        ((r[0] - p[0]) * dy - (r[1] - p[1]) * dx) / length
    };
    contours.iter().any(|ring| {
        (0..ring.len()).any(|i| {
            let c = ring[i];
            let d = ring[(i + 1) % ring.len()];
            let (s1, s2) = (side(a, b, c), side(a, b, d));
            let (s3, s4) = (side(c, d, a), side(c, d, b));
            ((s1 > eps && s2 < -eps) || (s1 < -eps && s2 > eps))
                && ((s3 > eps && s4 < -eps) || (s3 < -eps && s4 > eps))
        })
    })
}

/// Point membership in a polygon with holes, boundary included within `eps`.
pub(crate) fn closed_contains(contours: &[Vec<[f64; 2]>], p: [f64; 2], eps: f64) -> bool {
    let mut inside = false;
    for ring in contours {
        for i in 0..ring.len() {
            let a = ring[i];
            let b = ring[(i + 1) % ring.len()];
            let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
            let l2 = dx * dx + dy * dy;
            let t = if l2 > 0. {
                (((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / l2).clamp(0., 1.)
            } else {
                0.
            };
            if (p[0] - a[0] - t * dx).hypot(p[1] - a[1] - t * dy) <= eps {
                return true;
            }
            if (a[1] > p[1]) != (b[1] > p[1]) && p[0] < a[0] + (p[1] - a[1]) * dx / dy {
                inside = !inside;
            }
        }
    }
    inside
}

/// Holes inside the exterior, apart from it and from each other. Only
/// pairs with a ring flagged in `changed` are tested (the others were).
fn validate_holes(contours: &[Vec<[f64; 2]>], changed: &[bool]) -> Result<(), Error> {
    // A hole lies inside the exterior without touching it, and holes are
    // apart: with no contact between two rings (exact segment predicate)
    // one ring lies inside the other exactly when one of its vertices
    // does. Boxes rule pairs out first.
    let bbox = |r: &Vec<[f64; 2]>| {
        r.iter()
            .fold(([f64::MAX; 2], [f64::MIN; 2]), |(lo, hi), p| {
                (
                    [lo[0].min(p[0]), lo[1].min(p[1])],
                    [hi[0].max(p[0]), hi[1].max(p[1])],
                )
            })
    };
    let boxes: Vec<_> = contours.iter().map(bbox).collect();
    let ring = |r: &Vec<[f64; 2]>| {
        let points: Vec<_> = r.iter().chain(r.first()).map(|p| (p[0], p[1])).collect();
        Polygon::new(LineString::from(points), vec![])
    };
    let inside = |r: usize, p: [f64; 2]| ring(&contours[r]).contains(&geo::Point::new(p[0], p[1]));
    let disjoint = |a: usize, b: usize| {
        let ((alo, ahi), (blo, bhi)) = (boxes[a], boxes[b]);
        ahi[0] < blo[0] || bhi[0] < alo[0] || ahi[1] < blo[1] || bhi[1] < alo[1]
    };
    let within = |a: usize, b: usize| {
        let ((alo, ahi), (blo, bhi)) = (boxes[a], boxes[b]);
        alo[0] >= blo[0] && alo[1] >= blo[1] && ahi[0] <= bhi[0] && ahi[1] <= bhi[1]
    };
    for i in 1..contours.len() {
        if (changed[0] || changed[i])
            && (!within(i, 0)
                || rings_meet(&contours[0], &contours[i])
                || !inside(0, contours[i][0]))
        {
            return Err(Error::InvalidRing);
        }
        for j in 1..i {
            if (changed[i] || changed[j])
                && !disjoint(i, j)
                && (rings_meet(&contours[i], &contours[j])
                    || inside(i, contours[j][0])
                    || inside(j, contours[i][0]))
            {
                return Err(Error::InvalidRing);
            }
        }
    }
    Ok(())
}

/// Whether any segment of ring `a` meets (crosses or touches) any segment
/// of ring `b`: an x sweep over both, exact segment predicate.
fn rings_meet(a: &[[f64; 2]], b: &[[f64; 2]]) -> bool {
    let segment = |r: &[[f64; 2]], i: usize| {
        let (p, q) = (r[i], r[(i + 1) % r.len()]);
        (
            Line::new((p[0], p[1]), (q[0], q[1])),
            p[0].min(q[0]),
            p[0].max(q[0]),
            p[1].min(q[1]),
            p[1].max(q[1]),
        )
    };
    let mut items: Vec<(f64, bool, usize)> = (0..a.len())
        .map(|i| (segment(a, i).1, false, i))
        .chain((0..b.len()).map(|i| (segment(b, i).1, true, i)))
        .collect();
    items.sort_by(|x, y| x.0.total_cmp(&y.0));
    let mut active: [Vec<usize>; 2] = [vec![], vec![]];
    for (lo, of_b, i) in items {
        let (own, other) = if of_b { (b, a) } else { (a, b) };
        let (line, _, _, ylo, yhi) = segment(own, i);
        let side = usize::from(!of_b);
        active[side].retain(|&j| segment(other, j).2 >= lo);
        for &j in &active[side] {
            let (l, _, _, olo, ohi) = segment(other, j);
            if ohi >= ylo && olo <= yhi && line.intersects(&l) {
                return true;
            }
        }
        active[usize::from(of_b)].push(i);
    }
    false
}

fn validate_ring(points: &[[f64; 2]], epsilon: f64) -> Result<(), Error> {
    let n = points.len();
    let line = |i: usize| {
        Line::new(
            (points[i][0], points[i][1]),
            (points[(i + 1) % n][0], points[(i + 1) % n][1]),
        )
    };
    let o = points[0];
    let area: f64 = (0..n)
        .map(|i| {
            let a = points[i];
            let b = points[(i + 1) % n];
            (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
        })
        .sum();
    if !area.is_finite() || area.abs() <= epsilon * epsilon {
        return Err(Error::InvalidRing);
    }
    // Non-adjacent edges must not meet. Sweep over x: edges whose x ranges
    // do not overlap cannot meet, so only overlapping ones are tested with
    // the exact predicate (the same answers as testing every pair).
    let range = |i: usize, k: usize| {
        let (a, b) = (points[i][k], points[(i + 1) % n][k]);
        (a.min(b), a.max(b))
    };
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| range(a, 0).0.total_cmp(&range(b, 0).0));
    let mut active: Vec<usize> = vec![];
    for &i in &order {
        let (lo, _) = range(i, 0);
        active.retain(|&j| range(j, 0).1 >= lo);
        let (ylo, yhi) = range(i, 1);
        for &j in &active {
            let adjacent = (i + 1) % n == j || (j + 1) % n == i;
            let (jlo, jhi) = range(j, 1);
            if adjacent || jhi < ylo || jlo > yhi {
                continue;
            }
            if line(i).intersects(&line(j)) {
                return Err(Error::InvalidRing);
            }
        }
        active.push(i);
    }
    for i in 0..n {
        let a = points[(i + n - 1) % n];
        let b = points[i];
        let c = points[(i + 1) % n];
        let cross = (a[0] - b[0]) * (c[1] - b[1]) - (a[1] - b[1]) * (c[0] - b[0]);
        let dot = (a[0] - b[0]) * (c[0] - b[0]) + (a[1] - b[1]) * (c[1] - b[1]);
        if cross.abs() <= epsilon * epsilon && dot > 0.0 {
            return Err(Error::InvalidRing);
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    #[test]
    fn hole_validation_agrees_with_the_polygon_relation() {
        use super::{validate_holes, Contains, Intersects, LineString, Polygon};
        // The former definition, through the general polygon relation.
        let reference = |contours: &[Vec<[f64; 2]>]| {
            let polygon = |r: &Vec<[f64; 2]>| {
                let points: Vec<_> = r.iter().chain(r.first()).map(|p| (p[0], p[1])).collect();
                Polygon::new(LineString::from(points), vec![])
            };
            let outer = polygon(&contours[0]);
            for i in 1..contours.len() {
                let hole = polygon(&contours[i]);
                if !outer.contains(&hole) || outer.exterior().intersects(hole.exterior()) {
                    return false;
                }
                for previous in &contours[1..i] {
                    if polygon(previous).intersects(&hole) {
                        return false;
                    }
                }
            }
            true
        };
        // Squares and triangles on a coarse grid: touching, overlapping,
        // nested, outside and apart cases all occur.
        let mut seed = 12345u64;
        let mut next = |n: u64| {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((seed >> 33) % n) as f64
        };
        let outer = vec![[0., 0.], [8., 0.], [8., 6.], [4., 8.], [0., 6.]];
        let mut agree = (0, 0);
        for _ in 0..5000 {
            let mut contours = vec![outer.clone()];
            for _ in 0..1 + next(3) as usize {
                let (x, y, w) = (next(9) - 0.5, next(9) - 0.5, 0.5 + next(3) * 0.5);
                contours.push(if next(2) == 0. {
                    vec![[x, y], [x, y + w], [x + w, y + w], [x + w, y]]
                } else {
                    vec![[x, y], [x + w, y], [x, y + w]]
                });
            }
            let all = vec![true; contours.len()];
            let expected = reference(&contours);
            assert_eq!(
                validate_holes(&contours, &all).is_ok(),
                expected,
                "{contours:?}"
            );
            if expected {
                agree.0 += 1;
            } else {
                agree.1 += 1;
            }
        }
        // Both answers occur often.
        assert!(agree.0 > 500 && agree.1 > 500, "{agree:?}");
    }

    use super::*;
    fn square() -> (Model, usize, Vec<usize>) {
        let mut model = Model::new(1e-8, 0.03).unwrap();
        let plane = model.add_plane(PlaneFrame::new([0.; 3], [0., 0., 1.]).unwrap());
        let ring = [[0., 0., 0.], [2., 0., 0.], [2., 2., 0.], [0., 2., 0.]]
            .map(|p| model.add_vertex(p).unwrap())
            .to_vec();
        (model, plane, ring)
    }
    #[test]
    fn split_through_a_vertex_of_the_same_contour_is_rejected() {
        // A notch reaching within precision of the edge x = 2 (a gap closure
        // moved its tip there): splitting that edge at the tip, an exact
        // split without contour revalidation, would make the contour pass
        // through the tip twice.
        let (mut m, slab, r) = square();
        let tip = m.add_vertex([2. - 5e-9, 1., 0.]).unwrap();
        let (u, w) = (
            m.add_vertex([0., 1.05, 0.]).unwrap(),
            m.add_vertex([0., 0.95, 0.]).unwrap(),
        );
        m.add_surface(slab, vec![vec![r[0], r[1], r[2], r[3], u, tip, w]], vec![1])
            .unwrap();
        let edge = m.surfaces()[0].boundaries[0][1].edge;
        assert_eq!(m.edges[edge], [r[1].min(r[2]), r[1].max(r[2])]);
        let before = serde_json::to_string(&m).unwrap();
        assert_eq!(m.split_edge(edge, tip), Err(Error::InvalidRing));
        assert_eq!(before, serde_json::to_string(&m).unwrap());
    }

    #[test]
    fn perpendicular_surfaces_reuse_one_boundary() {
        let (mut m, slab, r) = square();
        m.add_surface(slab, vec![r.clone()], vec![7]).unwrap();
        let wall = m.add_plane(PlaneFrame::new([0.; 3], [0., 1., 0.]).unwrap());
        let a = m.add_vertex([0., 0., 3.]).unwrap();
        let b = m.add_vertex([2., 0., 3.]).unwrap();
        m.add_surface(wall, vec![vec![r[1], r[0], a, b]], vec![8])
            .unwrap();
        assert_eq!(m.edges().len(), 7);
        let e = &m.surfaces()[0].boundaries[0][0];
        let f = &m.surfaces()[1].boundaries[0][0];
        assert_eq!(e.edge, f.edge);
        assert_ne!(e.reversed, f.reversed);
        for surface in m.surfaces() {
            let frame = &m.planes()[surface.plane];
            for &uv in surface.contours.iter().flatten() {
                assert!(frame.distance(frame.lift(uv)).abs() < 1e-12);
            }
        }
    }
    #[test]
    fn warped_surface_is_rejected_transactionally() {
        let (mut m, p, mut r) = square();
        r[2] = m.add_vertex([2., 2., 0.001]).unwrap();
        let before = serde_json::to_string(&m).unwrap();
        assert_eq!(m.add_surface(p, vec![r], vec![1]), Err(Error::NonPlanar));
        assert_eq!(before, serde_json::to_string(&m).unwrap());
    }
    #[test]
    fn rejects_self_intersections_and_short_edges() {
        let (mut m, p, r) = square();
        assert_eq!(
            m.add_surface(p, vec![vec![r[0], r[2], r[1], r[3]]], vec![]),
            Err(Error::InvalidRing)
        );
        let near = m.add_vertex([0.001, 0., 0.]).unwrap();
        assert_eq!(
            m.add_surface(p, vec![vec![r[0], near, r[2], r[3]]], vec![]),
            Err(Error::ShortEdge)
        );
        assert!(m.edges().is_empty());
    }
    #[test]
    fn holes_must_be_disjoint_and_strictly_inside() {
        let (mut m, p, r) = square();
        let hole = [
            [0.5, 0.5, 0.],
            [1.5, 0.5, 0.],
            [1.5, 1.5, 0.],
            [0.5, 1.5, 0.],
        ]
        .map(|q| m.add_vertex(q).unwrap())
        .to_vec();
        assert!(m
            .add_surface(p, vec![r.clone(), hole.clone()], vec![9, 8, 9])
            .is_ok());
        assert_eq!(m.surfaces()[0].source_elements, vec![8, 9]);
        let before = m.edges().len();
        assert_eq!(
            m.add_surface(p, vec![r.clone(), hole.clone(), hole], vec![]),
            Err(Error::InvalidRing)
        );
        let outside = [[3., 0., 0.], [4., 0., 0.], [4., 1., 0.], [3., 1., 0.]]
            .map(|q| m.add_vertex(q).unwrap())
            .to_vec();
        assert_eq!(
            m.add_surface(p, vec![r, outside], vec![]),
            Err(Error::InvalidRing)
        );
        assert_eq!(m.edges().len(), before);
    }
    #[test]
    fn rigid_transform_preserves_shared_topology() {
        let (base, _, _) = square();
        let rotation = glam::DQuat::from_axis_angle(DVec3::new(1., 2., 3.).normalize(), 0.73);
        let shift = DVec3::new(100., -250., 73.);
        let mut m = Model::new(1e-8, 0.03).unwrap();
        let p = m.add_plane(
            PlaneFrame::new(shift.to_array(), (rotation * DVec3::Z).to_array()).unwrap(),
        );
        let ids = base
            .vertices
            .iter()
            .map(|v| {
                m.add_vertex((rotation * DVec3::from_array(*v) + shift).to_array())
                    .unwrap()
            })
            .collect();
        m.add_surface(p, vec![ids], vec![1]).unwrap();
        assert_eq!(m.edges().len(), 4);
        for &uv in &m.surfaces()[0].contours[0] {
            assert!(m.planes()[p].distance(m.planes()[p].lift(uv)).abs() < 1e-10);
        }
    }
    #[test]
    fn invalid_numeric_inputs_are_rejected() {
        assert!(PlaneFrame::new([0.; 3], [0.; 3]).is_err());
        assert!(PlaneFrame::new([f64::NAN, 0., 0.], [0., 0., 1.]).is_err());
        assert!(Model::new(0., 0.03).is_err());
        assert!(Model::new(1e-8, f64::INFINITY).is_err());
    }
}
