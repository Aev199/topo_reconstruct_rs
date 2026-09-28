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

#[derive(Debug, Clone, PartialEq)]
pub enum Error {
    InvalidPrecision,
    InvalidPlane,
    InvalidVertex,
    InvalidRing,
    NonPlanar,
    ShortEdge,
}

#[derive(Debug, Clone, Serialize)]
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
#[derive(Debug, Clone, Serialize)]
pub struct Model {
    precision: f64,
    minimum_edge: f64,
    planes: Vec<PlaneFrame>,
    vertices: Vec<[f64; 3]>,
    edges: Vec<[usize; 2]>,
    surfaces: Vec<Surface>,
    /// Edges referenced by no surface (replaced by pre-existing pieces during
    /// a split, or released by a trimmed console). Kept to preserve edge ids.
    orphaned_edges: Vec<usize>,
    #[serde(skip)]
    edge_index: BTreeMap<[usize; 2], usize>,
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
            vertices: vec![],
            edges: vec![],
            surfaces: vec![],
            orphaned_edges: vec![],
            edge_index: BTreeMap::new(),
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
        self.vertices.push(point);
        Ok(id)
    }
    pub fn surfaces(&self) -> &[Surface] {
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
        self.surfaces.push(Surface {
            plane,
            contours,
            boundaries,
            embedded_edges: vec![],
            source_elements,
        });
        Ok(id)
    }

    /// Validated plane contours of vertex rings (exterior first).
    fn contours(&self, plane: usize, rings: &[Vec<usize>]) -> Result<Vec<Vec<[f64; 2]>>, Error> {
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
            validate_ring(&uv, self.precision)?;
            contours.push(uv);
        }
        validate_holes(&contours)?;
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
                        let edge = *self.edge_index.entry(key).or_insert_with(|| {
                            let id = self.edges.len();
                            self.edges.push(key);
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

    /// Remove a vertex joining exactly two collinear edges, replacing them
    /// with one edge in every surface that uses them. The vertex must lie on
    /// the joined segment within `max(precision, deviation)`; contours change
    /// by at most that distance. Fails if the vertex has another edge, the
    /// edges have different users, or a resulting contour is invalid. The
    /// model is unchanged on error.
    pub fn remove_vertex(&mut self, vertex: usize, deviation: f64) -> Result<usize, Error> {
        let incident: Vec<usize> = (0..self.edges.len())
            .filter(|&e| {
                self.edges[e].contains(&vertex)
                    && (0..self.surfaces.len()).any(|s| self.surface_edges(s).any(|x| x == e))
            })
            .collect();
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
            (0..self.surfaces.len())
                .filter(|&s| self.surface_edges(s).any(|x| x == e))
                .collect()
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
            let mut ring_changed = false;
            for ring in &surface.boundaries {
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
                        ring_changed |= !keep;
                        keep
                    })
                    .collect();
                rings.push(ids);
            }
            let contours = self.contours(surface.plane, &rings)?;
            updates.push((s, rings, contours, ring_changed));
        }
        let merged = *self.edge_index.entry(key).or_insert_with(|| {
            self.edges.push(key);
            self.edges.len() - 1
        });
        for (s, rings, contours, ring_changed) in updates {
            if ring_changed {
                let boundaries = self.intern(&rings);
                let surface = &mut self.surfaces[s];
                surface.boundaries = boundaries;
                surface.contours = contours;
            }
            let surface = &mut self.surfaces[s];
            if surface.embedded_edges.contains(&e1) {
                surface.embedded_edges.retain(|&e| e != e1 && e != e2);
                surface.embedded_edges.push(merged);
            }
        }
        for e in [e1, e2] {
            if !self.orphaned_edges.contains(&e) {
                self.orphaned_edges.push(e);
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
        let before: BTreeSet<usize> = self.surface_edges(surface).collect();
        let boundaries = self.intern(&rings);
        let target = &mut self.surfaces[surface];
        target.contours = contours;
        target.boundaries = boundaries;
        target.embedded_edges = embedded;
        for e in before {
            if !(0..self.surfaces.len()).any(|s| self.surface_edges(s).any(|f| f == e))
                && !self.orphaned_edges.contains(&e)
            {
                self.orphaned_edges.push(e);
            }
        }
        Ok(())
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
        let [a, b] = *self.edges.get(edge).ok_or(Error::InvalidVertex)?;
        let p = DVec3::from_array(*self.vertices.get(vertex).ok_or(Error::InvalidVertex)?);
        if vertex == a || vertex == b {
            return Err(Error::InvalidVertex);
        }
        let pa = DVec3::from_array(self.vertices[a]);
        let pb = DVec3::from_array(self.vertices[b]);
        let d = pb - pa;
        let t = (p - pa).dot(d) / d.length_squared();
        if !(t > 0. && t < 1.) || p.distance(pa + d * t) > self.precision {
            return Err(Error::InvalidVertex);
        }
        if p.distance(pa) < self.minimum_edge || p.distance(pb) < self.minimum_edge {
            return Err(Error::ShortEdge);
        }
        let users: Vec<usize> = (0..self.surfaces.len())
            .filter(|&s| self.surface_edges(s).any(|e| e == edge))
            .collect();
        for &s in &users {
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
        self.edge_index.remove(&[a, b]);
        let e_first = match existing_first {
            Some(id) => id,
            None => {
                self.edges[edge] = first;
                self.edge_index.insert(first, edge);
                edge
            }
        };
        let e_second = match existing_second {
            Some(id) => id,
            None if e_first != edge => {
                self.edges[edge] = second;
                self.edge_index.insert(second, edge);
                edge
            }
            None => {
                self.edges.push(second);
                self.edge_index.insert(second, self.edges.len() - 1);
                self.edges.len() - 1
            }
        };
        if e_first != edge && e_second != edge {
            // Both pieces existed. The original id is left unreferenced with a
            // key no other edge uses; it is no longer indexed.
            self.orphaned_edges.push(edge);
        }
        for s in users {
            let plane = self.planes[self.surfaces[s].plane].clone();
            let uv = plane.project(p.to_array());
            let surface = &mut self.surfaces[s];
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
            validate_holes(&contours)?;
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
        self.vertices[vertex] = target;
        for (s, contours) in updates {
            self.surfaces[s].contours = contours;
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
        let edge = *self.edge_index.entry(key).or_insert_with(|| {
            self.edges.push(key);
            self.edges.len() - 1
        });
        self.surfaces[surface].embedded_edges.push(edge);
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

fn validate_holes(contours: &[Vec<[f64; 2]>]) -> Result<(), Error> {
    let polygon = |r: &Vec<[f64; 2]>| {
        let points: Vec<_> = r.iter().chain(r.first()).map(|p| (p[0], p[1])).collect();
        Polygon::new(LineString::from(points), vec![])
    };
    let outer = polygon(&contours[0]);
    for i in 1..contours.len() {
        let hole = polygon(&contours[i]);
        if !outer.contains(&hole) || outer.exterior().intersects(hole.exterior()) {
            return Err(Error::InvalidRing);
        }
        for previous in &contours[1..i] {
            if polygon(previous).intersects(&hole) {
                return Err(Error::InvalidRing);
            }
        }
    }
    Ok(())
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
    for i in 0..n {
        for j in i + 1..n {
            if j == i + 1 || (i == 0 && j == n - 1) {
                continue;
            }
            if line(i).intersects(&line(j)) {
                return Err(Error::InvalidRing);
            }
        }
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
