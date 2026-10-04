//! Explicit surface-surface junction lines.
//!
//! The intersection of two non-parallel surfaces, and the overlap of collinear
//! boundaries of coplanar surfaces, becomes a chain of shared model edges. An
//! edge on a contour stays a boundary edge; inside the material it is recorded
//! as an embedded edge. Existing edges are split at explicit vertices computed
//! from the participating planes. No vertex is moved and no material removed.
use crate::reconstruction::{closed_contains, Model, PlaneFrame};
use glam::DVec3;
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};

/// Planes whose normals differ less than this sine are treated as parallel.
const PARALLEL_SINE: f64 = 1e-8;
/// Detection tolerance relative to model precision. Matches the sensitivity of
/// the independent audit; created geometry must still meet model precision.
const DETECTION_FACTOR: f64 = 5.;

#[derive(Debug, Clone, Copy, Serialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum Kind {
    /// On the contour of both surfaces.
    Boundary,
    /// On the contour of exactly one surface.
    TJunction,
    /// Inside the material of both surfaces.
    Crossing,
}

#[derive(Debug, Serialize)]
pub struct Junction {
    pub surfaces: [usize; 2],
    pub kind: Kind,
    pub start: [f64; 3],
    pub end: [f64; 3],
    pub length: f64,
    /// Chain of model edges representing the junction after insertion.
    pub edges: Vec<usize>,
    pub generated_vertices: Vec<usize>,
    pub split_edges: usize,
    pub embedded_edges: usize,
}

#[derive(Debug, Serialize)]
pub struct Issue {
    pub surfaces: Vec<usize>,
    pub start: [f64; 3],
    pub end: [f64; 3],
    pub reason: String,
}

/// A vertex moved onto a junction to avoid a parasitic short edge.
#[derive(Debug, Serialize)]
pub struct Snap {
    pub vertex: usize,
    pub surfaces: [usize; 2],
    pub from: [f64; 3],
    pub to: [f64; 3],
    pub distance: f64,
}

#[derive(Debug, Default, Serialize)]
pub struct Report {
    pub detection_tolerance: f64,
    pub candidate_pairs: usize,
    /// Junction segments that were already represented by shared edges.
    pub already_conforming: usize,
    /// Junction segments whose topology was changed.
    pub junctions: Vec<Junction>,
    /// Vertices created at edge/plane intersections; they have no source node.
    pub generated_vertices: Vec<usize>,
    pub split_edges: usize,
    pub embedded_edges: usize,
    /// Existing vertices moved onto junctions, each by less than the minimum
    /// edge length, keeping all planes they lie on.
    pub snapped_vertices: Vec<Snap>,
    /// Splits resolving embedded edges that touched or crossed other edges of
    /// the same surface without a shared vertex.
    pub crossing_splits: usize,
    pub issues: Vec<Issue>,
}

fn p3(v: [f64; 3]) -> DVec3 {
    DVec3::from_array(v)
}

struct Line {
    origin: DVec3,
    direction: DVec3,
}
impl Line {
    fn at(&self, t: f64) -> DVec3 {
        self.origin + self.direction * t
    }
    fn parameter(&self, p: DVec3) -> f64 {
        (p - self.origin).dot(self.direction)
    }
    fn distance(&self, p: DVec3) -> f64 {
        p.distance(self.at(self.parameter(p)))
    }
}

fn normal(plane: &PlaneFrame) -> DVec3 {
    p3(plane.normal)
}

/// Parameter intervals of `line` inside the closed material of a surface.
fn clip(model: &Model, surface: usize, line: &Line, tol: f64) -> Vec<(f64, f64)> {
    let s = &model.surfaces[surface];
    let plane = &model.planes[s.plane];
    let o = plane.project(line.origin.to_array());
    let e = plane.project((line.origin + line.direction).to_array());
    let d = [e[0] - o[0], e[1] - o[1]];
    let norm = d[0].hypot(d[1]);
    if norm < 0.5 {
        return vec![];
    }
    let d = [d[0] / norm, d[1] / norm];
    let side = |p: [f64; 2]| (p[0] - o[0]) * d[1] - (p[1] - o[1]) * d[0];
    let along = |p: [f64; 2]| (p[0] - o[0]) * d[0] + (p[1] - o[1]) * d[1];
    let mut ts = vec![];
    for ring in &s.contours {
        for i in 0..ring.len() {
            let c = ring[i];
            let f = ring[(i + 1) % ring.len()];
            let (sc, sf) = (side(c), side(f));
            if sc.abs() <= tol {
                ts.push(along(c));
            }
            if sf.abs() <= tol {
                ts.push(along(f));
            }
            if (sc > tol && sf < -tol) || (sc < -tol && sf > tol) {
                let k = sc / (sc - sf);
                ts.push(along([c[0] + (f[0] - c[0]) * k, c[1] + (f[1] - c[1]) * k]));
            }
        }
    }
    ts.sort_by(f64::total_cmp);
    ts.dedup_by(|b, a| (*b - *a).abs() <= tol);
    let mut out: Vec<(f64, f64)> = vec![];
    for w in ts.windows(2) {
        let m = (w[0] + w[1]) / 2.;
        let p = [o[0] + d[0] * m, o[1] + d[1] * m];
        if closed_contains(&s.contours, p, tol) {
            match out.last_mut() {
                Some(last) if (last.1 - w[0]).abs() <= tol => last.1 = w[1],
                _ => out.push((w[0], w[1])),
            }
        }
    }
    out
}

fn on_contour(model: &Model, surface: usize, p: DVec3, tol: f64) -> bool {
    let s = &model.surfaces[surface];
    let uv = model.planes[s.plane].project(p.to_array());
    s.contours.iter().any(|ring| {
        (0..ring.len()).any(|i| {
            let a = ring[i];
            let b = ring[(i + 1) % ring.len()];
            let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
            let l2 = dx * dx + dy * dy;
            let t = (((uv[0] - a[0]) * dx + (uv[1] - a[1]) * dy) / l2).clamp(0., 1.);
            (uv[0] - a[0] - t * dx).hypot(uv[1] - a[1] - t * dy) <= tol
        })
    })
}

/// Parameter of `p` strictly inside edge `[a, b]` and its distance to it.
fn on_segment(p: DVec3, a: DVec3, b: DVec3) -> (f64, f64) {
    let d = b - a;
    let t = (p - a).dot(d) / d.length_squared();
    (t, p.distance(a + d * t.clamp(0., 1.)))
}

fn bounds(model: &Model, surface: usize) -> (DVec3, DVec3) {
    let mut low = DVec3::splat(f64::INFINITY);
    let mut high = DVec3::splat(f64::NEG_INFINITY);
    for e in model.surface_edges(surface) {
        for v in model.edges[e] {
            low = low.min(p3(model.vertices[v]));
            high = high.max(p3(model.vertices[v]));
        }
    }
    (low, high)
}

fn surface_edge_set(model: &Model, surfaces: &[usize]) -> Vec<usize> {
    let set: BTreeSet<_> = surfaces
        .iter()
        .flat_map(|&s| model.surface_edges(s).collect::<Vec<_>>())
        .collect();
    set.into_iter().collect()
}

/// Split every edge of `surfaces` that contains `vertex` in its interior.
fn split_through(
    model: &mut Model,
    surfaces: &[usize],
    vertex: usize,
    tol: f64,
    splits: &mut usize,
) -> Result<(), String> {
    let p = p3(model.vertices[vertex]);
    loop {
        let found = surface_edge_set(model, surfaces).into_iter().find(|&e| {
            let [a, b] = model.edges[e];
            if a == vertex || b == vertex {
                return false;
            }
            let (t, d) = on_segment(p, p3(model.vertices[a]), p3(model.vertices[b]));
            t > 0. && t < 1. && d <= tol
        });
        let Some(edge) = found else {
            return Ok(());
        };
        model.split_edge_within(edge, vertex, tol).map_err(|e| {
            let [a, b] = model.edges[edge].map(|v| p3(model.vertices[v]));
            format!(
                "edge_split_{e:?}: distances to edge ends {:e}, {:e}",
                p.distance(a),
                p.distance(b)
            )
        })?;
        *splits += 1;
    }
}

fn revert(model: &mut Model, moves: &[(usize, DVec3, DVec3)]) {
    for &(v, from, _) in moves.iter().rev() {
        model
            .move_vertex(v, from.to_array())
            .expect("restoring a validated vertex position");
    }
}

struct Segment<'a> {
    surfaces: [usize; 2],
    line: &'a Line,
    start: f64,
    end: f64,
}

/// Fixed inputs of the junction pass.
pub struct Context<'a> {
    /// Vertices that must be nodes of a surface without lying on its edges
    /// (bar contacts, retained nodes of closed openings), per surface.
    pub interior: &'a [BTreeSet<usize>],
    /// Vertices that must not move (bar anchors and interior nodes).
    pub locked: &'a BTreeSet<usize>,
    /// Largest move closing the free end of a junction line (a wall end)
    /// onto another line of the same surface. Other vertices move less than
    /// the minimum edge length.
    pub wall_end_tolerance: f64,
}

/// Move an existing vertex that lies within the minimum edge length of a
/// required junction vertex onto the junction, instead of creating a
/// parasitic short edge. The vertex keeps every plane it lies on. Returns
/// whether the vertex was moved; its new position is read from the model.
fn snap(
    model: &mut Model,
    v: usize,
    required: [f64; 3],
    junction_planes: [&PlaneFrame; 2],
    line: &Line,
    tol: f64,
) -> bool {
    // Preferably onto the required point itself (on the crossed edge and
    // both junction planes); move_vertex checks every plane of the vertex.
    if model.move_vertex(v, required).is_ok() {
        return true;
    }
    let from = p3(model.vertices[v]);
    let mut planes: Vec<PlaneFrame> = (0..model.surfaces.len())
        .filter(|&s| model.surface_edges(s).any(|e| model.edges[e].contains(&v)))
        .map(|s| model.planes[model.surfaces[s].plane].clone())
        .collect();
    planes.extend(junction_planes.iter().map(|p| (*p).clone()));
    let refs: Vec<&PlaneFrame> = planes.iter().collect();
    let Some(target) = super::intersection(from, &refs, model.precision) else {
        return false;
    };
    from.distance(target) < model.minimum_edge
        && line.distance(target) <= tol
        && model.move_vertex(v, target.to_array()).is_ok()
}

fn process(
    model: &mut Model,
    seg: Segment<'_>,
    context: &Context<'_>,
    tol: f64,
    report: &mut Report,
) {
    let interior = context.interior;
    let [a, b] = seg.surfaces;
    let eps = model.precision;
    let line = seg.line;
    let (start, end) = (line.at(seg.start), line.at(seg.end));
    let issue = |report: &mut Report, reason: String| {
        report.issues.push(Issue {
            surfaces: vec![a, b],
            start: start.to_array(),
            end: end.to_array(),
            reason,
        })
    };
    if start.distance(end) < model.minimum_edge {
        issue(report, "junction_shorter_than_minimum_edge".into());
        return;
    }
    let pa = &model.planes[model.surfaces[a].plane].clone();
    let pb = &model.planes[model.surfaces[b].plane].clone();
    let edges = surface_edge_set(model, &[a, b]);
    let mut moves: Vec<(usize, DVec3, DVec3)> = vec![];
    // Mandatory interior mesh nodes on the line join the chain; otherwise the
    // mesher would place an unrelated subdivision vertex next to them.
    let edge_vertices: BTreeSet<usize> = edges
        .iter()
        .flat_map(|&e| model.edges[e])
        .chain(interior[a].iter().copied())
        .chain(interior[b].iter().copied())
        .collect();
    // Existing vertices on the junction, followed by new endpoint vertices.
    let mut chain: Vec<(f64, Option<usize>, DVec3)> = edge_vertices
        .iter()
        .filter_map(|&v| {
            let p = p3(model.vertices[v]);
            let t = line.parameter(p);
            (line.distance(p) <= tol && t >= seg.start - tol && t <= seg.end + tol).then_some((
                t,
                Some(v),
                p,
            ))
        })
        .collect();
    for t in [seg.start, seg.end] {
        if chain.iter().any(|(u, _, _)| (u - t).abs() <= tol) {
            continue;
        }
        let p = line.at(t);
        let containing: Vec<usize> = edges
            .iter()
            .copied()
            .filter(|&e| {
                let [i, j] = model.edges[e];
                let (u, d) = on_segment(p, p3(model.vertices[i]), p3(model.vertices[j]));
                u > 0. && u < 1. && d <= tol
            })
            .collect();
        // Prefer an edge crossing the junction: its intersection with the
        // other surface plane is exact for every surface sharing the edge.
        let mut exact = None;
        for &e in &containing {
            let [i, j] = model.edges[e].map(|v| p3(model.vertices[v]));
            if (j - i).normalize().cross(line.direction).length() < 1e-6 {
                continue;
            }
            let in_a = model.surface_edges(a).any(|x| x == e);
            let other = if in_a { pb } else { pa };
            let (di, dj) = (other.distance(i.to_array()), other.distance(j.to_array()));
            if (di - dj).abs() > 0. {
                exact = Some(i + (j - i) * (di / (di - dj)));
                break;
            }
        }
        let q = exact.or_else(|| {
            containing.first().map(|&e| {
                let [i, j] = model.edges[e].map(|v| p3(model.vertices[v]));
                let d = j - i;
                i + d * ((p - i).dot(d) / d.length_squared())
            })
        });
        let Some(q) = q else {
            issue(report, "junction_endpoint_off_existing_edges".into());
            revert(model, &moves);
            return;
        };
        // A vertex already on the line next to the required endpoint is the
        // endpoint: it moves onto it instead of leaving a tiny gap.
        if let Some(k) = chain
            .iter()
            .enumerate()
            .filter(|(_, c)| {
                c.1.is_some_and(|v| !context.locked.contains(&v))
                    && c.2.distance(q) < model.minimum_edge
            })
            .min_by(|x, y| x.1 .2.distance(q).total_cmp(&y.1 .2.distance(q)))
            .map(|(k, _)| k)
        {
            let v = chain[k].1.unwrap();
            let from = p3(model.vertices[v]);
            if snap(model, v, q.to_array(), [pa, pb], line, tol) {
                let target = p3(model.vertices[v]);
                moves.push((v, from, target));
                chain[k] = (line.parameter(target), Some(v), target);
                continue;
            }
        }
        let near = edge_vertices
            .iter()
            .copied()
            .filter(|v| !context.locked.contains(v) && chain.iter().all(|c| c.1 != Some(*v)))
            .map(|v| (p3(model.vertices[v]).distance(q), v))
            .filter(|&(d, _)| d < model.minimum_edge)
            .min_by(|x, y| x.0.total_cmp(&y.0));
        if let Some((_, v)) = near {
            let from = p3(model.vertices[v]);
            if snap(model, v, q.to_array(), [pa, pb], line, tol) {
                let target = p3(model.vertices[v]);
                moves.push((v, from, target));
                chain.push((line.parameter(target), Some(v), target));
                continue;
            }
        }
        chain.push((line.parameter(q), None, q));
    }
    chain.sort_by(|x, y| x.0.total_cmp(&y.0));
    chain.dedup_by(|y, x| y.1.is_some() && y.1 == x.1);
    for w in chain.windows(2) {
        let gap = w[0].2.distance(w[1].2);
        if gap <= eps {
            // Distinct vertices at one location: possibly an intentional seam
            // (e.g. a hinge). Identity is never merged by coordinates alone.
            issue(report, format!("coincident_distinct_vertices: {gap:e}"));
            revert(model, &moves);
            return;
        }
        if gap < model.minimum_edge {
            issue(
                report,
                format!("junction_vertices_closer_than_minimum_edge: {gap:e}"),
            );
            revert(model, &moves);
            return;
        }
    }
    for &(_, _, q) in &chain {
        let off = pa
            .distance(q.to_array())
            .abs()
            .max(pb.distance(q.to_array()).abs());
        if off > eps {
            issue(report, format!("junction_vertex_off_plane: {off:e}"));
            revert(model, &moves);
            return;
        }
    }
    // Every split must leave two pieces of at least the minimum edge length;
    // checked before mutation so a failed junction leaves no trace.
    for &(_, _, q) in &chain {
        for &e in &edges {
            let [i, j] = model.edges[e].map(|v| p3(model.vertices[v]));
            let (u, d) = on_segment(q, i, j);
            let (di, dj) = (q.distance(i), q.distance(j));
            if u > 0.
                && u < 1.
                && d <= tol
                && di > eps
                && dj > eps
                && di.min(dj) < model.minimum_edge
            {
                issue(
                    report,
                    format!("junction_vertex_near_edge_end: {:e}", di.min(dj)),
                );
                revert(model, &moves);
                return;
            }
        }
    }
    for &(vertex, from, to) in &moves {
        report.snapped_vertices.push(Snap {
            vertex,
            surfaces: [a, b],
            from: from.to_array(),
            to: to.to_array(),
            distance: from.distance(to),
        });
    }
    // All checks passed; the remaining operations are individually valid
    // refinements of shared topology and never move existing vertices.
    let mut generated = vec![];
    let mut ids = vec![];
    for &(_, v, q) in &chain {
        let id = match v {
            Some(v) => v,
            None => {
                let id = model
                    .add_vertex(q.to_array())
                    .expect("finite junction vertex");
                generated.push(id);
                id
            }
        };
        ids.push(id);
    }
    let mut splits = 0;
    for &v in &ids {
        if let Err(reason) = split_through(model, &[a, b], v, tol, &mut splits) {
            issue(report, reason);
            report.split_edges += splits;
            report.generated_vertices.extend(generated);
            return;
        }
    }
    // A chain end just outside a surface (a wall end a few micrometres
    // beyond a slab edge) is put on that surface's contour: the nearest
    // boundary edge bends through it, less than the minimum edge length.
    for &v in [ids[0], ids[ids.len() - 1]].iter() {
        for s in [a, b] {
            let own = model.surface_edges(s).any(|e| model.edges[e].contains(&v));
            let surface = &model.surfaces[s];
            let uv = model.planes[surface.plane].project(model.vertices[v]);
            if own || closed_contains(&surface.contours, uv, eps) {
                continue;
            }
            let p = p3(model.vertices[v]);
            let nearest = surface
                .boundaries
                .iter()
                .flatten()
                .map(|u| u.edge)
                .filter_map(|e| {
                    let [i, j] = model.edges[e];
                    let (t, d) = on_segment(p, p3(model.vertices[i]), p3(model.vertices[j]));
                    (t > 0. && t < 1. && d < model.minimum_edge).then_some((d, e))
                })
                .min_by(|x, y| x.0.total_cmp(&y.0).then(x.1.cmp(&y.1)));
            if let Some((_, e)) = nearest {
                if model.split_edge_within(e, v, model.minimum_edge).is_ok() {
                    splits += 1;
                }
            }
        }
    }
    let mut chain_edges = vec![];
    let mut embedded = 0;
    for w in ids.windows(2) {
        let existing = model.edge_between(w[0], w[1]);
        for s in [a, b] {
            let known = existing.is_some_and(|e| model.surface_edges(s).any(|x| x == e));
            match model.embed_edge(s, w[0], w[1]) {
                Ok(e) => {
                    if !known {
                        embedded += 1;
                    }
                    if s == b {
                        chain_edges.push(e);
                    }
                }
                Err(e) => {
                    issue(report, format!("embed_{e:?}"));
                    report.split_edges += splits;
                    report.embedded_edges += embedded;
                    report.generated_vertices.extend(generated);
                    return;
                }
            }
        }
    }
    report.split_edges += splits;
    report.embedded_edges += embedded;
    if splits == 0 && embedded == 0 && generated.is_empty() && moves.is_empty() {
        report.already_conforming += 1;
        return;
    }
    let mid = line.at((seg.start + seg.end) / 2.);
    let kind = match (
        on_contour(model, a, mid, tol),
        on_contour(model, b, mid, tol),
    ) {
        (true, true) => Kind::Boundary,
        (false, false) => Kind::Crossing,
        _ => Kind::TJunction,
    };
    report.generated_vertices.extend(&generated);
    report.junctions.push(Junction {
        surfaces: [a, b],
        kind,
        start: start.to_array(),
        end: end.to_array(),
        length: start.distance(end),
        edges: chain_edges,
        generated_vertices: generated,
        split_edges: splits,
        embedded_edges: embedded,
    });
}

/// Junction segments of a pair of non-parallel surfaces.
fn angled(model: &Model, a: usize, b: usize, tol: f64) -> Option<(Line, Vec<(f64, f64)>)> {
    let pa = &model.planes[model.surfaces[a].plane];
    let pb = &model.planes[model.surfaces[b].plane];
    let (na, nb) = (normal(pa), normal(pb));
    let u = na.cross(nb);
    if u.length() < PARALLEL_SINE {
        return None;
    }
    let (ha, hb) = (na.dot(p3(pa.origin)), nb.dot(p3(pb.origin)));
    let direction = u.normalize();
    let mut origin = (nb.cross(u) * ha + u.cross(na) * hb) / u.length_squared();
    let (la, ha_) = bounds(model, a);
    let (lb, hb_) = bounds(model, b);
    let centre = (la + ha_ + lb + hb_) / 4.;
    origin += direction * (centre - origin).dot(direction);
    let line = Line { origin, direction };
    let ia = clip(model, a, &line, tol);
    let ib = clip(model, b, &line, tol);
    let mut out = vec![];
    for &(x0, x1) in &ia {
        for &(y0, y1) in &ib {
            let (s, e) = (x0.max(y0), x1.min(y1));
            if e - s > tol {
                out.push((s, e));
            }
        }
    }
    Some((line, out))
}

/// Overlaps of collinear, distinct boundary edges of two coplanar surfaces.
fn coplanar(model: &Model, a: usize, b: usize, tol: f64) -> Vec<(Line, f64, f64)> {
    let pa = &model.planes[model.surfaces[a].plane];
    let pb = &model.planes[model.surfaces[b].plane];
    if normal(pa).cross(normal(pb)).length() >= PARALLEL_SINE
        || model
            .surface_edges(b)
            .flat_map(|e| model.edges[e])
            .any(|v| pa.distance(model.vertices[v]).abs() > tol)
    {
        return vec![];
    }
    let boundary = |s: usize| -> Vec<usize> {
        model.surfaces[s]
            .boundaries
            .iter()
            .flatten()
            .map(|e| e.edge)
            .collect()
    };
    let eb = boundary(b);
    let mut out = vec![];
    for ea in boundary(a) {
        if eb.contains(&ea) {
            continue;
        }
        let [i, j] = model.edges[ea].map(|v| p3(model.vertices[v]));
        let line = Line {
            origin: i,
            direction: (j - i).normalize(),
        };
        let length = i.distance(j);
        for &f in &eb {
            let [k, l] = model.edges[f].map(|v| p3(model.vertices[v]));
            if line.distance(k) > tol || line.distance(l) > tol {
                continue;
            }
            let (tk, tl) = (line.parameter(k), line.parameter(l));
            let (s, e) = (tk.min(tl).max(0.), tk.max(tl).min(length));
            if e - s > tol {
                out.push((
                    Line {
                        origin: line.origin,
                        direction: line.direction,
                    },
                    s,
                    e,
                ));
            }
        }
    }
    out
}

enum Contact {
    /// A vertex lies on the interior of an edge without being part of it.
    Touch { vertex: usize, edge: usize },
    /// A vertex lies within the minimum edge length of an edge interior.
    NearTouch {
        vertex: usize,
        edge: usize,
        distance: f64,
        limit: f64,
    },
    /// Two edges cross at an interior point of both.
    Cross { edges: [usize; 2], point: DVec3 },
}

/// First contact without a shared vertex between an embedded edge of `s` and
/// any other edge of `s`, skipping near touches already found unresolvable.
fn find_contact(
    model: &Model,
    s: usize,
    tol: f64,
    wall_end: f64,
    skip: &BTreeSet<(usize, usize)>,
) -> Option<Contact> {
    let plane = &model.planes[model.surfaces[s].plane];
    let all: Vec<usize> = model.surface_edges(s).collect();
    // Free ends of junction lines: vertices with a single edge in the
    // surface, that edge embedded.
    let mut degree = BTreeMap::<usize, usize>::new();
    for &e in &all {
        for v in model.edges[e] {
            *degree.entry(v).or_default() += 1;
        }
    }
    let free: BTreeSet<usize> = model.surfaces[s]
        .embedded_edges
        .iter()
        .flat_map(|&e| model.edges[e])
        .filter(|v| degree[v] == 1)
        .collect();
    let mut near: Option<Contact> = None;
    let mut nearest = f64::INFINITY;
    for &e in &model.surfaces[s].embedded_edges {
        let [i, j] = model.edges[e];
        let (pi, pj) = (p3(model.vertices[i]), p3(model.vertices[j]));
        for &f in &all {
            if f == e {
                continue;
            }
            let [k, l] = model.edges[f];
            let (pk, pl) = (p3(model.vertices[k]), p3(model.vertices[l]));
            for (v, p, target, a, b) in [
                (k, pk, e, pi, pj),
                (l, pl, e, pi, pj),
                (i, pi, f, pk, pl),
                (j, pj, f, pk, pl),
            ] {
                let [x, y] = model.edges[target];
                if v == x || v == y {
                    continue;
                }
                let (t, d) = on_segment(p, a, b);
                if !(t > 0. && t < 1.) {
                    continue;
                }
                if d <= tol {
                    return Some(Contact::Touch {
                        vertex: v,
                        edge: target,
                    });
                }
                let limit = if free.contains(&v) {
                    wall_end.max(model.minimum_edge)
                } else {
                    model.minimum_edge
                };
                if d < limit && d < nearest && !skip.contains(&(v, target)) {
                    nearest = d;
                    near = Some(Contact::NearTouch {
                        vertex: v,
                        edge: target,
                        distance: d,
                        limit,
                    });
                }
            }
            if [i, j].iter().any(|v| [k, l].contains(v)) {
                continue;
            }
            let uv = |p: DVec3| plane.project(p.to_array());
            let (a, b, c, d) = (uv(pi), uv(pj), uv(pk), uv(pl));
            let cross = |o: [f64; 2], p: [f64; 2], q: [f64; 2]| {
                (p[0] - o[0]) * (q[1] - o[1]) - (p[1] - o[1]) * (q[0] - o[0])
            };
            let (d1, d2) = (cross(a, b, c), cross(a, b, d));
            let (d3, d4) = (cross(c, d, a), cross(c, d, b));
            let (lab, lcd) = (pi.distance(pj), pk.distance(pl));
            if ((d1 > tol * lab && d2 < -tol * lab) || (d1 < -tol * lab && d2 > tol * lab))
                && ((d3 > tol * lcd && d4 < -tol * lcd) || (d3 < -tol * lcd && d4 > tol * lcd))
            {
                return Some(Contact::Cross {
                    edges: [e, f],
                    point: pk + (pl - pk) * (d1 / (d1 - d2)),
                });
            }
        }
    }
    near
}

/// Embedded edges must not touch or cross other edges of the same surface
/// without a shared vertex. Such contacts arise between junctions inserted
/// for different surface pairs, e.g. two walls meeting on a slab, and where a
/// wall ends a few micrometres short of a slab edge. The latter is closed by
/// moving the vertex onto the edge (less than the minimum edge length, on
/// every plane it lies on) instead of leaving a parasitic gap.
fn resolve_crossings(model: &mut Model, context: &Context<'_>, tol: f64, report: &mut Report) {
    let locked = context.locked;
    for s in 0..model.surfaces.len() {
        let mut skip = BTreeSet::new();
        for _ in 0..256 {
            if model.surfaces[s].embedded_edges.is_empty() {
                break;
            }
            let Some(contact) = find_contact(model, s, tol, context.wall_end_tolerance, &skip)
            else {
                break;
            };
            let unresolved = |model: &Model, report: &mut Report, edge: usize, reason: String| {
                let [x, y] = model.edges[edge];
                report.issues.push(Issue {
                    surfaces: vec![s],
                    start: model.vertices[x],
                    end: model.vertices[y],
                    reason,
                });
            };
            match contact {
                Contact::Touch { vertex, edge } => {
                    if let Err(e) = model.split_edge(edge, vertex) {
                        unresolved(model, report, edge, format!("embedded_contact_split_{e:?}"));
                        return;
                    }
                    report.crossing_splits += 1;
                }
                Contact::Cross { edges, point } => {
                    let vertex = model
                        .add_vertex(point.to_array())
                        .expect("finite crossing vertex");
                    report.generated_vertices.push(vertex);
                    for edge in edges {
                        if let Err(e) = model.split_edge(edge, vertex) {
                            unresolved(
                                model,
                                report,
                                edge,
                                format!("embedded_crossing_split_{e:?}"),
                            );
                            return;
                        }
                        report.crossing_splits += 1;
                    }
                }
                Contact::NearTouch {
                    vertex,
                    edge,
                    distance,
                    limit,
                } => {
                    skip.insert((vertex, edge));
                    let from = p3(model.vertices[vertex]);
                    let mut planes: Vec<PlaneFrame> = (0..model.surfaces.len())
                        .filter(|&u| {
                            model
                                .surface_edges(u)
                                .any(|e| e == edge || model.edges[e].contains(&vertex))
                        })
                        .map(|u| model.planes[model.surfaces[u].plane].clone())
                        .collect();
                    let [c, d] = model.edges[edge].map(|v| p3(model.vertices[v]));
                    let foot = c + (d - c) * ((from - c).dot(d - c) / (d - c).length_squared());
                    let refs: Vec<&PlaneFrame> = planes.iter().collect();
                    let target = super::intersection(foot, &refs, model.precision)
                        .filter(|q| q.distance(from) < limit);
                    planes.clear();
                    let moved = !locked.contains(&vertex)
                        && target.is_some_and(|q| model.move_vertex(vertex, q.to_array()).is_ok());
                    if !moved || model.split_edge(edge, vertex).is_err() {
                        if moved {
                            model
                                .move_vertex(vertex, from.to_array())
                                .expect("restoring a validated vertex position");
                        }
                        unresolved(
                            model,
                            report,
                            edge,
                            format!("embedded_near_touch: {distance:e}"),
                        );
                        continue;
                    }
                    let to = p3(model.vertices[vertex]);
                    report.snapped_vertices.push(Snap {
                        vertex,
                        surfaces: [s, s],
                        from: from.to_array(),
                        to: to.to_array(),
                        distance: from.distance(to),
                    });
                    report.crossing_splits += 1;
                }
            }
        }
    }
}

/// Contour edges lying on the material of another surface within the
/// minimum edge length of its plane, but beyond the detection tolerance
/// (e.g. a wall top a micrometre below a slab it was not connected to in
/// the source), are moved onto that plane so the junction is detected. Each
/// vertex keeps every plane it lies on and moves less than the minimum edge.
fn settle(model: &mut Model, context: &Context<'_>, tol: f64, report: &mut Report) {
    let count = model.surfaces.len();
    let limit = model.minimum_edge;
    let boxes: Vec<_> = (0..count).map(|s| bounds(model, s)).collect();
    for a in 0..count {
        for b in 0..count {
            let ((la, ha), (lb, hb)) = (boxes[a], boxes[b]);
            if a == b || (ha + limit).cmplt(lb).any() || (hb + limit).cmplt(la).any() {
                continue;
            }
            let pb = model.planes[model.surfaces[b].plane].clone();
            let pa = &model.planes[model.surfaces[a].plane];
            if normal(pa).cross(normal(&pb)).length() < 1e-6 {
                continue;
            }
            let edges: Vec<usize> = model.surfaces[a]
                .boundaries
                .iter()
                .flatten()
                .map(|u| u.edge)
                .collect();
            for e in edges {
                let [i, j] = model.edges[e];
                let (pi, pj) = (p3(model.vertices[i]), p3(model.vertices[j]));
                let (hi, hj) = (pb.distance(pi.to_array()), pb.distance(pj.to_array()));
                if hi.abs().max(hj.abs()) <= tol || hi.abs().max(hj.abs()) >= limit {
                    continue;
                }
                let mid = pb.project(((pi + pj) / 2.).to_array());
                if !closed_contains(&model.surfaces[b].contours, mid, limit) {
                    continue;
                }
                for (v, h) in [(i, hi), (j, hj)] {
                    if h.abs() <= tol || context.locked.contains(&v) {
                        continue;
                    }
                    let from = p3(model.vertices[v]);
                    let mut planes: Vec<PlaneFrame> = (0..count)
                        .filter(|&s| model.surface_edges(s).any(|x| model.edges[x].contains(&v)))
                        .map(|s| model.planes[model.surfaces[s].plane].clone())
                        .collect();
                    planes.push(pb.clone());
                    let refs: Vec<&PlaneFrame> = planes.iter().collect();
                    let Some(target) = super::intersection(from, &refs, model.precision) else {
                        continue;
                    };
                    if from.distance(target) < limit
                        && model.move_vertex(v, target.to_array()).is_ok()
                    {
                        report.snapped_vertices.push(Snap {
                            vertex: v,
                            surfaces: [a, b],
                            from: from.to_array(),
                            to: target.to_array(),
                            distance: from.distance(target),
                        });
                    }
                }
            }
        }
    }
}

/// Insert all detectable surface junctions into the shared topology.
pub fn insert(model: &mut Model, context: &Context<'_>) -> Report {
    let tol = model.precision * DETECTION_FACTOR;
    let mut report = Report {
        detection_tolerance: tol,
        ..Default::default()
    };
    settle(model, context, tol, &mut report);
    let count = model.surfaces.len();
    let boxes: Vec<_> = (0..count).map(|s| bounds(model, s)).collect();
    for a in 0..count {
        for b in a + 1..count {
            let ((la, ha), (lb, hb)) = (boxes[a], boxes[b]);
            if (ha + tol).cmplt(lb).any() || (hb + tol).cmplt(la).any() {
                continue;
            }
            report.candidate_pairs += 1;
            if let Some((line, segments)) = angled(model, a, b, tol) {
                for (start, end) in segments {
                    let seg = Segment {
                        surfaces: [a, b],
                        line: &line,
                        start,
                        end,
                    };
                    process(model, seg, context, tol, &mut report);
                }
            } else {
                for (line, start, end) in coplanar(model, a, b, tol) {
                    let seg = Segment {
                        surfaces: [a, b],
                        line: &line,
                        start,
                        end,
                    };
                    process(model, seg, context, tol, &mut report);
                }
            }
        }
    }
    resolve_crossings(model, context, tol, &mut report);
    report
}

/// Insert the junction of two given surfaces (the editor's "connect"):
/// their intersection line, T-junction or common boundary becomes shared
/// edges, then embedded edges touching other edges are resolved.
pub fn insert_pair(model: &mut Model, context: &Context<'_>, a: usize, b: usize) -> Report {
    let tol = model.precision * DETECTION_FACTOR;
    let mut report = Report {
        detection_tolerance: tol,
        candidate_pairs: 1,
        ..Default::default()
    };
    if let Some((line, segments)) = angled(model, a, b, tol) {
        for (start, end) in segments {
            let seg = Segment {
                surfaces: [a, b],
                line: &line,
                start,
                end,
            };
            process(model, seg, context, tol, &mut report);
        }
    } else {
        for (line, start, end) in coplanar(model, a, b, tol) {
            let seg = Segment {
                surfaces: [a, b],
                line: &line,
                start,
                end,
            };
            process(model, seg, context, tol, &mut report);
        }
    }
    resolve_crossings(model, context, tol, &mut report);
    report
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::reconstruction::PlaneFrame;
    use glam::DQuat;
    use std::collections::BTreeMap;

    /// A rigid transform, a scale and an optional reversal of every normal.
    pub(crate) struct Placement {
        pub(crate) rotation: DQuat,
        pub(crate) shift: DVec3,
        pub(crate) scale: f64,
        pub(crate) flip: bool,
    }
    impl Placement {
        pub(crate) fn all() -> Vec<Placement> {
            let mut out = vec![];
            for scale in [0.1, 1., 25.] {
                for (rotation, shift, flip) in [
                    (DQuat::IDENTITY, DVec3::ZERO, false),
                    (
                        DQuat::from_axis_angle(DVec3::new(1., 2., 3.).normalize(), 0.61),
                        DVec3::new(812.5, -304.25, 71.),
                        true,
                    ),
                ] {
                    out.push(Placement {
                        rotation,
                        shift,
                        scale,
                        flip,
                    });
                }
            }
            out
        }
        pub(crate) fn point(&self, p: [f64; 3]) -> [f64; 3] {
            (self.rotation * (DVec3::from_array(p) * self.scale) + self.shift).to_array()
        }
    }

    /// Rings of coordinates and the plane normal of one surface.
    pub(crate) type Panel = (Vec<Vec<[f64; 3]>>, [f64; 3]);

    /// Surfaces given as rings of coordinates. Equal coordinates in the input
    /// description denote one explicit vertex; the placement never merges.
    pub(crate) fn build(place: &Placement, surfaces: &[Panel]) -> Model {
        // The production precision: micro-offset cases mean what they say.
        let mut model =
            Model::new(super::super::PRECISION * place.scale, 0.001 * place.scale).unwrap();
        let mut ids = BTreeMap::new();
        for (rings, normal) in surfaces {
            let n = place.rotation * DVec3::from_array(*normal);
            let n = if place.flip { -n } else { n };
            let frame = PlaneFrame::new(place.point(rings[0][0]), n.to_array()).unwrap();
            let plane = model.add_plane(frame);
            let rings = rings
                .iter()
                .map(|ring| {
                    ring.iter()
                        .map(|p| {
                            let key = p.map(|x| (x * 1e6).round() as i64);
                            *ids.entry(key)
                                .or_insert_with(|| model.add_vertex(place.point(*p)).unwrap())
                        })
                        .collect()
                })
                .collect();
            model.add_surface(plane, rings, vec![]).unwrap();
        }
        model
    }

    pub(crate) fn slab(x0: f64, x1: f64) -> Panel {
        (
            vec![vec![[x0, 0., 0.], [x1, 0., 0.], [x1, 4., 0.], [x0, 4., 0.]]],
            [0., 0., 1.],
        )
    }
    pub(crate) fn wall(x0: f64, x1: f64, z0: f64, z1: f64) -> Panel {
        (
            vec![vec![[x0, 2., z0], [x1, 2., z0], [x1, 2., z1], [x0, 2., z1]]],
            [0., 1., 0.],
        )
    }
    pub(crate) fn no_interior(model: &Model) -> Vec<BTreeSet<usize>> {
        vec![BTreeSet::new(); model.surfaces.len()]
    }
    pub(crate) fn run(model: &mut Model) -> Report {
        run_with(model, 0.)
    }
    pub(crate) fn run_with(model: &mut Model, wall_end_tolerance: f64) -> Report {
        let interior = no_interior(model);
        let locked = BTreeSet::new();
        insert(
            model,
            &Context {
                interior: &interior,
                locked: &locked,
                wall_end_tolerance,
            },
        )
    }

    /// The junction between `a` and `b` along the segment `p`-`q` (in the
    /// description frame) is covered by edges shared by both surfaces.
    fn shared_cover(
        model: &Model,
        place: &Placement,
        a: usize,
        b: usize,
        p: [f64; 3],
        q: [f64; 3],
    ) {
        let (p, q) = (p3(place.point(p)), p3(place.point(q)));
        let length = p.distance(q);
        let d = (q - p) / length;
        let eps = model.precision * 5.;
        let shared: BTreeSet<usize> = model
            .surface_edges(a)
            .filter(|&e| model.surface_edges(b).any(|f| f == e))
            .collect();
        let mut parts: Vec<(f64, f64)> = shared
            .iter()
            .filter_map(|&e| {
                let [i, j] = model.edges[e].map(|v| p3(model.vertices[v]));
                let (ti, tj) = ((i - p).dot(d), (j - p).dot(d));
                (i.distance(p + d * ti) <= eps && j.distance(p + d * tj) <= eps)
                    .then_some((ti.min(tj), ti.max(tj)))
            })
            .collect();
        parts.sort_by(|x, y| x.0.total_cmp(&y.0));
        let mut cursor = 0.;
        for (s, e) in parts {
            assert!(s <= cursor + eps, "gap at {cursor} before {s}");
            cursor = f64::max(cursor, e);
        }
        assert!(cursor >= length - eps, "covered to {cursor} of {length}");
    }

    fn assert_idempotent(model: &mut Model) {
        let edges = model.edges.len();
        let vertices = model.vertices.len();
        let again = run(model);
        assert!(
            again.junctions.is_empty() && again.issues.is_empty(),
            "{again:?}"
        );
        assert_eq!((edges, vertices), (model.edges.len(), model.vertices.len()));
    }

    #[test]
    fn wall_vertex_micrometres_before_the_junction_end_becomes_the_end() {
        for place in Placement::all() {
            // The wall top has a vertex 3 um before the slab edge x = 2.
            let wall = (
                vec![vec![
                    [1., 2., -2.],
                    [3., 2., -2.],
                    [3., 2., 0.],
                    [1.999997, 2., 0.],
                    [1., 2., 0.],
                ]],
                [0., 1., 0.],
            );
            let mut m = build(&place, &[slab(0., 2.), wall]);
            let r = run(&mut m);
            assert!(r.issues.is_empty(), "{:?}", r.issues);
            shared_cover(&m, &place, 0, 1, [1., 2., 0.], [2., 2., 0.]);
            assert_idempotent(&mut m);
        }
    }

    #[test]
    fn wall_top_micrometres_below_a_slab_is_settled_and_joined() {
        for place in Placement::all() {
            // 7 um below the slab: beyond the detection tolerance, far
            // below the minimum edge length.
            let mut m = build(&place, &[slab(0., 4.), wall(1., 3., -2., -7e-6)]);
            let r = run(&mut m);
            assert!(r.issues.is_empty(), "{:?}", r.issues);
            assert_eq!(r.snapped_vertices.len(), 2);
            for snap in &r.snapped_vertices {
                assert!(snap.distance < 1e-5 * place.scale);
            }
            shared_cover(&m, &place, 0, 1, [1., 2., 0.], [3., 2., 0.]);
            assert_idempotent(&mut m);
            // A millimetre gap is geometry, not noise: nothing moves.
            let mut m = build(&place, &[slab(0., 4.), wall(1., 3., -2., -0.002)]);
            let r = run(&mut m);
            assert!(r.snapped_vertices.is_empty() && r.junctions.is_empty());
        }
    }

    #[test]
    fn wall_end_just_beyond_a_slab_edge_bends_the_edge_through_it() {
        for place in Placement::all() {
            // The wall corner is held 30 um beyond the slab edge x = 2 by a
            // perpendicular wall, so it cannot move onto the edge.
            let x = 2.00003;
            let cross = (
                vec![vec![[x, 1., -2.], [x, 2., -2.], [x, 3., -2.], [x, 2., 0.]]],
                [1., 0., 0.],
            );
            // The slab edge carries a vertex 7 mm past the wall line.
            let slab = (
                vec![vec![
                    [0., 0., 0.],
                    [2., 0., 0.],
                    [2., 2.007, 0.],
                    [2., 4., 0.],
                    [0., 4., 0.],
                ]],
                [0., 0., 1.],
            );
            let mut m = build(&place, &[slab, wall(1., x, -2., 0.), cross]);
            let r = run(&mut m);
            assert!(r.issues.is_empty(), "{:?}", r.issues);
            let corner = (0..m.vertices.len())
                .find(|&v| {
                    p3(m.vertices[v]).distance(p3(place.point([x, 2., 0.])))
                        < 1e-9 * place.scale.max(1.)
                })
                .unwrap();
            assert!(m.surfaces[0].boundaries[0]
                .iter()
                .any(|e| m.edges[e.edge].contains(&corner)));
            assert_idempotent(&mut m);
        }
    }

    #[test]
    fn t_junction_in_panel_interior_becomes_shared_embedded_edge() {
        for place in Placement::all() {
            let mut m = build(&place, &[slab(0., 4.), wall(1., 3., 0., 2.)]);
            let r = run(&mut m);
            assert!(r.issues.is_empty(), "{:?}", r.issues);
            assert_eq!(r.junctions.len(), 1);
            assert_eq!(r.junctions[0].kind, Kind::TJunction);
            assert!(r.generated_vertices.is_empty());
            // The wall keeps its contour; the slab embeds the same edge.
            assert!(m.surfaces[1].embedded_edges.is_empty());
            assert_eq!(m.surfaces[0].embedded_edges.len(), 1);
            let bottom = m.surfaces[0].embedded_edges[0];
            assert!(m.surfaces[1].boundaries[0].iter().any(|e| e.edge == bottom));
            shared_cover(&m, &place, 0, 1, [1., 2., 0.], [3., 2., 0.]);
            assert_idempotent(&mut m);
        }
    }

    #[test]
    fn interior_crossing_is_split_on_both_panels() {
        for place in Placement::all() {
            let mut m = build(&place, &[slab(0., 4.), wall(1., 3., -1., 1.)]);
            let r = run(&mut m);
            assert!(r.issues.is_empty(), "{:?}", r.issues);
            assert_eq!(r.junctions.len(), 1);
            assert_eq!(r.junctions[0].kind, Kind::Crossing);
            // New vertices on the two vertical wall edges, exact on both planes.
            assert_eq!(r.generated_vertices.len(), 2);
            for &v in &r.generated_vertices {
                for s in &m.surfaces {
                    assert!(m.planes[s.plane].distance(m.vertices[v]).abs() <= m.precision);
                }
            }
            assert_eq!(m.surfaces[1].contours[0].len(), 6);
            assert_eq!(m.surfaces[0].embedded_edges, m.surfaces[1].embedded_edges);
            shared_cover(&m, &place, 0, 1, [1., 2., 0.], [3., 2., 0.]);
            assert_idempotent(&mut m);
        }
    }

    #[test]
    fn junction_leaving_a_panel_splits_its_boundary() {
        for place in Placement::all() {
            // The wall base continues past the slab edge at x = 4.
            let mut m = build(&place, &[slab(0., 4.), wall(2., 6., 0., 2.)]);
            let r = run(&mut m);
            assert!(r.issues.is_empty(), "{:?}", r.issues);
            assert_eq!(r.generated_vertices.len(), 1);
            // Slab contour and wall contour both receive the crossing vertex.
            assert_eq!(m.surfaces[0].contours[0].len(), 5);
            assert_eq!(m.surfaces[1].contours[0].len(), 5);
            shared_cover(&m, &place, 0, 1, [2., 2., 0.], [4., 2., 0.]);
            // The exterior part of the wall base stays a wall-only edge.
            let v = r.generated_vertices[0];
            let outside = m
                .vertices
                .iter()
                .position(|&p| p3(p).distance(p3(place.point([6., 2., 0.]))) < m.precision);
            let outer = m.edge_between(v, outside.unwrap()).unwrap();
            assert!(!m.surface_edges(0).any(|e| e == outer));
            assert_idempotent(&mut m);
        }
    }

    #[test]
    fn junction_ending_on_existing_boundary_vertex_reuses_it() {
        for place in Placement::all() {
            // Wall base from the slab corner region: ends exactly at (0, 2)
            // where the slab has an explicit contour vertex.
            let slab = (
                vec![vec![
                    [0., 0., 0.],
                    [4., 0., 0.],
                    [4., 4., 0.],
                    [0., 4., 0.],
                    [0., 2., 0.],
                ]],
                [0., 0., 1.],
            );
            let mut m = build(&place, &[slab, wall(0., 3., 0., 2.)]);
            let r = run(&mut m);
            assert!(r.issues.is_empty(), "{:?}", r.issues);
            assert!(r.generated_vertices.is_empty());
            assert_eq!(r.split_edges, 0);
            shared_cover(&m, &place, 0, 1, [0., 2., 0.], [3., 2., 0.]);
            assert_idempotent(&mut m);
        }
    }

    #[test]
    fn junction_across_property_regions_is_split_at_their_boundary() {
        for place in Placement::all() {
            // Two coplanar property regions share the edge x = 2.
            let mut m = build(&place, &[slab(0., 2.), slab(2., 4.), wall(1., 3., 0., 2.)]);
            let r = run(&mut m);
            assert!(r.issues.is_empty(), "{:?}", r.issues);
            // One vertex at the region boundary, shared by all three surfaces.
            assert_eq!(r.generated_vertices.len(), 1);
            let v = r.generated_vertices[0];
            for s in 0..3 {
                assert!(m.surface_edges(s).any(|e| m.edges[e].contains(&v)));
            }
            shared_cover(&m, &place, 0, 2, [1., 2., 0.], [2., 2., 0.]);
            shared_cover(&m, &place, 1, 2, [2., 2., 0.], [3., 2., 0.]);
            shared_cover(&m, &place, 0, 1, [2., 0., 0.], [2., 4., 0.]);
            assert_idempotent(&mut m);
        }
    }

    #[test]
    fn nearby_but_separate_surfaces_are_left_unchanged() {
        for place in Placement::all() {
            // A 10 mm gap is a modelling question, not an intersection.
            let mut m = build(&place, &[slab(0., 4.), wall(1., 3., 0.01, 2.)]);
            let before = serde_json::to_string(&m).unwrap();
            let r = run(&mut m);
            assert!(r.junctions.is_empty() && r.issues.is_empty());
            assert_eq!(before, serde_json::to_string(&m).unwrap());
        }
    }

    #[test]
    fn walls_meeting_on_a_slab_share_the_contact_vertex() {
        for place in Placement::all() {
            // Two walls on the slab, the second perpendicular and ending in the
            // middle of the first wall base: an embedded T inside the slab.
            let cross = (
                vec![vec![
                    [2., 2., 0.],
                    [2., 3.5, 0.],
                    [2., 3.5, 2.],
                    [2., 2., 2.],
                ]],
                [1., 0., 0.],
            );
            let mut m = build(&place, &[slab(0., 4.), wall(1., 3., 0., 2.), cross]);
            let r = run(&mut m);
            assert!(r.issues.is_empty(), "{:?}", r.issues);
            shared_cover(&m, &place, 0, 1, [1., 2., 0.], [3., 2., 0.]);
            shared_cover(&m, &place, 0, 2, [2., 2., 0.], [2., 3.5, 0.]);
            shared_cover(&m, &place, 1, 2, [2., 2., 0.], [2., 2., 2.]);
            assert_idempotent(&mut m);
        }
    }

    #[test]
    fn micro_offset_corner_is_snapped_instead_of_a_short_edge() {
        for place in Placement::all() {
            // The wall base ends 50 µm past the slab edge x = 4: a source
            // misalignment. Splitting there would leave a 50 µm edge.
            let mut m = build(&place, &[slab(0., 4.), wall(2., 4.00005, 0., 2.)]);
            let r = run(&mut m);
            assert!(r.issues.is_empty(), "{:?}", r.issues);
            assert_eq!(r.snapped_vertices.len(), 1);
            let snap = &r.snapped_vertices[0];
            assert!(snap.distance <= 5.1e-5 * place.scale, "{}", snap.distance);
            assert!(p3(snap.to).distance(p3(place.point([4., 2., 0.]))) <= m.precision);
            assert!(r.generated_vertices.is_empty());
            for s in 0..2 {
                for e in m.surface_edges(s) {
                    let [i, j] = m.edges[e].map(|v| p3(m.vertices[v]));
                    assert!(i.distance(j) >= m.minimum_edge);
                }
            }
            shared_cover(&m, &place, 0, 1, [2., 2., 0.], [4., 2., 0.]);
            assert_idempotent(&mut m);
        }
    }

    #[test]
    fn wall_ending_just_short_of_slab_edge_is_closed() {
        for place in Placement::all() {
            // The wall base stops 50 µm before the slab edge x = 4: the
            // embedded line would end in a parasitic 50 µm gap.
            let mut m = build(&place, &[slab(0., 4.), wall(2., 3.99995, 0., 2.)]);
            let r = run(&mut m);
            assert!(r.issues.is_empty(), "{:?}", r.issues);
            assert_eq!(r.snapped_vertices.len(), 1);
            assert!(r.snapped_vertices[0].distance <= 5.1e-5 * place.scale);
            // The corner now lies on the slab contour, which is split there.
            let v = r.snapped_vertices[0].vertex;
            assert!(m.surfaces[0].boundaries[0]
                .iter()
                .any(|e| m.edges[e.edge].contains(&v)));
            shared_cover(&m, &place, 0, 1, [2., 2., 0.], [4., 2., 0.]);
            assert_idempotent(&mut m);
        }
    }

    #[test]
    fn wall_end_short_of_a_perpendicular_wall_axis_is_closed_within_tolerance() {
        for place in Placement::all() {
            // Wall A on the slab along y = 2; wall B along x = 3 ends 30 mm
            // short of A's axis. Both stand on the slab.
            let b = (
                vec![vec![
                    [3., 2.03, 0.],
                    [3., 3.5, 0.],
                    [3., 3.5, 2.],
                    [3., 2.03, 2.],
                ]],
                [1., 0., 0.],
            );
            let panels = [slab(0., 4.), wall(1., 3.5, 0., 2.), b];
            // Without a wall-end tolerance the 30 mm gap stays.
            let mut m = build(&place, &panels);
            let r = run(&mut m);
            assert!(r.snapped_vertices.is_empty());
            // With 50 mm the base corner of B closes onto A's axis; its top
            // corner is not on the slab and stays, so B's end leans slightly.
            let mut m = build(&place, &panels);
            let r = run_with(&mut m, 0.05 * place.scale);
            assert!(r.issues.is_empty(), "{:?}", r.issues);
            assert_eq!(r.snapped_vertices.len(), 1);
            let snap = &r.snapped_vertices[0];
            assert!((snap.distance - 0.03 * place.scale).abs() < 1e-9 * place.scale);
            assert!(p3(snap.to).distance(p3(place.point([3., 2., 0.]))) <= m.precision * 10.);
            shared_cover(&m, &place, 0, 2, [3., 2., 0.], [3., 3.5, 0.]);
            shared_cover(&m, &place, 0, 1, [1., 2., 0.], [3.5, 2., 0.]);
        }
    }

    #[test]
    fn coincident_distinct_vertices_are_reported_not_merged() {
        for place in Placement::all() {
            // Two coplanar panels meeting along x = 2 with separate vertices
            // at identical positions: a possible intentional seam.
            let mut m = build(&place, &[slab(0., 2.)]);
            let n = place.rotation * DVec3::Z;
            let n = if place.flip { -n } else { n };
            let plane =
                m.add_plane(PlaneFrame::new(place.point([2., 0., 0.]), n.to_array()).unwrap());
            let ring = [[2., 0., 0.], [4., 0., 0.], [4., 4., 0.], [2., 4., 0.]]
                .map(|p| m.add_vertex(place.point(p)).unwrap())
                .to_vec();
            m.add_surface(plane, vec![ring], vec![]).unwrap();
            let before = serde_json::to_string(&m).unwrap();
            let r = run(&mut m);
            assert!(r
                .issues
                .iter()
                .any(|i| i.reason.starts_with("coincident_distinct_vertices")));
            assert_eq!(before, serde_json::to_string(&m).unwrap());
        }
    }

    #[test]
    fn interior_node_on_junction_joins_the_chain() {
        for place in Placement::all() {
            let mut m = build(&place, &[slab(0., 4.), wall(1., 3., 0., 2.)]);
            let node = m.add_vertex(place.point([2., 2., 0.])).unwrap();
            let mut interior = no_interior(&m);
            interior[0].insert(node);
            let locked = interior[0].clone();
            let r = insert(
                &mut m,
                &Context {
                    interior: &interior,
                    locked: &locked,
                    wall_end_tolerance: 0.,
                },
            );
            assert!(r.issues.is_empty(), "{:?}", r.issues);
            // The wall base is split at the retained slab node.
            assert!(m.surface_edges(1).any(|e| m.edges[e].contains(&node)));
            shared_cover(&m, &place, 0, 1, [1., 2., 0.], [3., 2., 0.]);
        }
    }
}
