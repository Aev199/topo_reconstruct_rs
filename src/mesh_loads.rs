//! Loads of the LIRA model carried onto a Gmsh mesh: every load case stays a
//! case, a pressure goes to the shell elements by the exact part of each
//! element inside its contour (the contour is that of the source), line
//! loads to the bar elements along them (clipped to the elements they
//! cover; along a plate edge, to its nodes), point loads to the nodes of the
//! element they act on (by barycentric weights, so force and moment are
//! kept). The force and the moment of every case are kept; the report
//! compares them with those of the loads on the geometry.
use crate::loads::Load;
use crate::meshing::Mesh;
use crate::reconstruction::assembly::edit::State;
use geo::{Centroid, Area, BooleanOps, Coord, LineString, MultiPolygon, Polygon};
use glam::{DVec2, DVec3};
use serde::Serialize;
use std::collections::BTreeMap;

/// A linear load on the part `t` (fractions of the element from its first node) of a bar element.
#[derive(Debug, Clone, Serialize)]
pub struct BarLoad {
    pub case: u32,
    pub bar: usize,
    pub t: [f64; 2],
    /// kN/m (global) at the two ends of the part.
    pub q: [[f64; 3]; 2],
}

/// Forces in kN, global axes.
#[derive(Debug, Clone, Default, Serialize)]
pub struct MeshLoads {
    /// (case, shell, kN/m2 vector): the pressure acts on the whole element, scaled by
    /// the part of it inside the contour.
    pub pressures: Vec<(u32, usize, [f64; 3])>,
    pub bar_loads: Vec<BarLoad>,
    /// (case, node) -> force and moment.
    pub nodal: BTreeMap<(u32, usize), [f64; 6]>,
    /// Loads that found no element (point loads far from the mesh, pressure
    /// contours without shells).
    pub lost: BTreeMap<u32, [f64; 3]>,
}

impl MeshLoads {
    /// Resultant force per case (kN).
    pub fn resultants(&self, mesh: &Mesh) -> BTreeMap<u32, DVec3> {
        self.resultants_about(mesh, DVec3::ZERO).into_iter().map(|(c, v)| (c, v.0)).collect()
    }

    /// Force, moment about `origin` and the scale of the moment (the sum of |F| |r|) per case.
    pub fn resultants_about(&self, mesh: &Mesh, origin: DVec3) -> BTreeMap<u32, (DVec3, DVec3, f64)> {
        let mut out: BTreeMap<u32, (DVec3, DVec3, f64)> = BTreeMap::new();
        fn add(out: &mut BTreeMap<u32, (DVec3, DVec3, f64)>, origin: DVec3, case: u32, at: DVec3, f: DVec3, m: DVec3) {
            let r = at - origin;
            let e = out.entry(case).or_default();
            e.0 += f;
            e.1 += r.cross(f) + m;
            e.2 += f.length() * r.length() + m.length();
        }
        for &(case, shell, p) in &self.pressures {
            let (area, centre) = shell_area_centroid(mesh, shell);
            add(&mut out, origin, case, centre, DVec3::from_array(p) * area, DVec3::ZERO);
        }
        for l in &self.bar_loads {
            let [a, b] = mesh.bars[l.bar].nodes;
            let (pa, pb) = (point(mesh, a), point(mesh, b));
            let (s, e) = (pa.lerp(pb, l.t[0]), pa.lerp(pb, l.t[1]));
            let (qa, qb) = (DVec3::from_array(l.q[0]), DVec3::from_array(l.q[1]));
            let length = s.distance(e);
            // A linear load: force and moment from its two end values.
            let (d, qd) = (e - s, qb - qa);
            let (ra, rs) = (s - origin, length);
            let force = (qa + qb) / 2. * rs;
            let moment = (ra.cross(qa) + (ra.cross(qd) + d.cross(qa)) / 2. + d.cross(qd) / 3.) * rs;
            let en = out.entry(l.case).or_default();
            en.0 += force;
            en.1 += moment;
            en.2 += (qa.length() + qb.length()) / 2. * rs * ((s + e) / 2. - origin).length();
        }
        for (&(case, node), f) in &self.nodal {
            add(&mut out, origin, case, point(mesh, node), DVec3::new(f[0], f[1], f[2]), DVec3::new(f[3], f[4], f[5]));
        }
        out
    }
}

fn point(mesh: &Mesh, n: usize) -> DVec3 {
    DVec3::from_array(mesh.nodes[n])
}

/// Area and centroid of a shell element (a quad as a fan of triangles).
pub fn shell_area_centroid(mesh: &Mesh, shell: usize) -> (f64, DVec3) {
    let p: Vec<DVec3> = mesh.shells[shell].nodes.iter().map(|&n| point(mesh, n)).collect();
    let (mut area, mut first) = (0., DVec3::ZERO);
    for i in 1..p.len() - 1 {
        let t = (p[i] - p[0]).cross(p[i + 1] - p[0]).length() / 2.;
        area += t;
        first += (p[0] + p[i] + p[i + 1]) / 3. * t;
    }
    if area > 0. { (area, first / area) } else { (0., p[0]) }
}

pub fn shell_area(mesh: &Mesh, shell: usize) -> f64 {
    shell_area_centroid(mesh, shell).0
}

pub fn bar_length(mesh: &Mesh, bar: usize) -> f64 {
    let [a, b] = mesh.bars[bar].nodes;
    point(mesh, a).distance(point(mesh, b))
}

fn inside(polygon: &[DVec2], q: DVec2) -> bool {
    let mut odd = false;
    for i in 0..polygon.len() {
        let (a, b) = (polygon[i], polygon[(i + 1) % polygon.len()]);
        if (a.y > q.y) != (b.y > q.y) && q.x < a.x + (q.y - a.y) / (b.y - a.y) * (b.x - a.x) {
            odd = !odd;
        }
    }
    odd
}

/// Distance of `p` from the segment `a`-`b` and the parameter of its foot.
/// Distance of `p` from the infinite line through `a`, `b`, and the position of its foot along it (0 at `a`, 1 at `b`).
fn line_distance(p: DVec3, a: DVec3, b: DVec3) -> (f64, f64) {
    let d = b - a;
    let t = (p - a).dot(d) / d.length_squared();
    ((p - (a + d * t)).length(), t)
}

fn segment_distance(p: DVec3, a: DVec3, b: DVec3) -> (f64, f64) {
    let d = b - a;
    let t = if d.length_squared() > 0. { ((p - a).dot(d) / d.length_squared()).clamp(0., 1.) } else { 0. };
    (p.distance(a + d * t), t)
}

fn add_force(out: &mut MeshLoads, case: u32, node: usize, f: DVec3, m: DVec3) {
    let entry = out.nodal.entry((case, node)).or_default();
    for k in 0..3 {
        entry[k] += f[k];
        entry[k + 3] += m[k];
    }
}

/// Shell elements by the cells of a grid (a point load finds its element).
struct ShellIndex {
    cell: f64,
    cells: std::collections::HashMap<[i64; 3], Vec<usize>>,
}

impl ShellIndex {
    fn new(mesh: &Mesh) -> Self {
        let mut total = 0.;
        for s in &mesh.shells {
            let (a, b) = (point(mesh, s.nodes[0]), point(mesh, s.nodes[1]));
            total += a.distance(b);
        }
        let cell = if mesh.shells.is_empty() { 1. } else { (total / mesh.shells.len() as f64 * 2.).max(1e-3) };
        let mut cells: std::collections::HashMap<[i64; 3], Vec<usize>> = Default::default();
        for (i, s) in mesh.shells.iter().enumerate() {
            let p: Vec<DVec3> = s.nodes.iter().map(|&n| point(mesh, n)).collect();
            let lo = p.iter().fold(DVec3::splat(f64::MAX), |l, q| l.min(*q));
            let hi = p.iter().fold(DVec3::splat(f64::MIN), |h, q| h.max(*q));
            let (a, b) = (Self::key(cell, lo), Self::key(cell, hi));
            for x in a[0]..=b[0] {
                for y in a[1]..=b[1] {
                    for z in a[2]..=b[2] {
                        cells.entry([x, y, z]).or_default().push(i);
                    }
                }
            }
        }
        ShellIndex { cell, cells }
    }

    fn key(cell: f64, p: DVec3) -> [i64; 3] {
        [0, 1, 2].map(|k| (p[k] / cell).floor() as i64)
    }

    /// The nodes and barycentric weights of `p` in the shell element it lies on
    /// (within `tolerance` of the element's plane), else `None`.
    fn locate(&self, mesh: &Mesh, p: DVec3, tolerance: f64) -> Option<Vec<(usize, f64)>> {
        let mut best: Option<(f64, Vec<(usize, f64)>)> = None;
        for &s in self.cells.get(&Self::key(self.cell, p)).into_iter().flatten() {
            let nodes = &mesh.shells[s].nodes;
            // Triangles of the element (a quad: two).
            let triangles: Vec<[usize; 3]> = if nodes.len() == 3 { vec![[0, 1, 2]] } else { vec![[0, 1, 2], [0, 2, 3]] };
            for t in triangles {
                let (a, b, c) = (point(mesh, nodes[t[0]]), point(mesh, nodes[t[1]]), point(mesh, nodes[t[2]]));
                let n = (b - a).cross(c - a);
                if n.length() < 1e-12 {
                    continue;
                }
                let n = n.normalize();
                let distance = (p - a).dot(n).abs();
                if distance > tolerance {
                    continue;
                }
                let q = p - n * (p - a).dot(n);
                // Barycentric coordinates in the triangle.
                let (v0, v1, v2) = (b - a, c - a, q - a);
                let (d00, d01, d11, d20, d21) = (v0.dot(v0), v0.dot(v1), v1.dot(v1), v2.dot(v0), v2.dot(v1));
                let det = d00 * d11 - d01 * d01;
                if det.abs() < 1e-18 {
                    continue;
                }
                let l1 = (d11 * d20 - d01 * d21) / det;
                let l2 = (d00 * d21 - d01 * d20) / det;
                let l0 = 1. - l1 - l2;
                if l0 >= -1e-9 && l1 >= -1e-9 && l2 >= -1e-9 && best.as_ref().is_none_or(|x| distance < x.0) {
                    best = Some((distance, vec![(nodes[t[0]], l0.max(0.)), (nodes[t[1]], l1.max(0.)), (nodes[t[2]], l2.max(0.))]));
                }
            }
        }
        best.map(|(_, w)| {
            let sum: f64 = w.iter().map(|x| x.1).sum();
            w.into_iter().map(|(n, l)| (n, l / sum)).collect()
        })
    }
}

/// The distinct edges of the shell elements by the cells of a grid (a line
/// load along an element edge finds the edges on its line).
struct EdgeIndex {
    cell: f64,
    edges: Vec<[usize; 2]>,
    cells: std::collections::HashMap<[i64; 3], Vec<usize>>,
}

impl EdgeIndex {
    fn new(mesh: &Mesh) -> Self {
        let mut seen: std::collections::BTreeSet<(usize, usize)> = Default::default();
        let mut edges = vec![];
        for s in &mesh.shells {
            for i in 0..s.nodes.len() {
                let (a, b) = (s.nodes[i], s.nodes[(i + 1) % s.nodes.len()]);
                if seen.insert((a.min(b), a.max(b))) {
                    edges.push([a, b]);
                }
            }
        }
        let total: f64 = edges.iter().map(|e| point(mesh, e[0]).distance(point(mesh, e[1]))).sum();
        let cell = if edges.is_empty() { 1. } else { (total / edges.len() as f64 * 2.).max(1e-3) };
        let mut cells: std::collections::HashMap<[i64; 3], Vec<usize>> = Default::default();
        for (i, e) in edges.iter().enumerate() {
            let (a, b) = (point(mesh, e[0]), point(mesh, e[1]));
            let (lo, hi) = (ShellIndex::key(cell, a.min(b)), ShellIndex::key(cell, a.max(b)));
            for x in lo[0]..=hi[0] {
                for y in lo[1]..=hi[1] {
                    for z in lo[2]..=hi[2] {
                        cells.entry([x, y, z]).or_default().push(i);
                    }
                }
            }
        }
        EdgeIndex { cell, edges, cells }
    }

    /// The edges whose cells the segment `a`-`b` passes (a superset of those on it).
    fn near(&self, a: DVec3, b: DVec3) -> Vec<[usize; 2]> {
        let steps = ((a.distance(b) / (self.cell / 4.)).ceil() as usize).clamp(1, 1_000_000);
        let mut found: std::collections::BTreeSet<usize> = Default::default();
        for k in 0..=steps {
            let p = a.lerp(b, k as f64 / steps as f64);
            let key = ShellIndex::key(self.cell, p);
            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        if let Some(list) = self.cells.get(&[key[0] + dx, key[1] + dy, key[2] + dz]) {
                            found.extend(list.iter().copied());
                        }
                    }
                }
            }
        }
        found.into_iter().map(|i| self.edges[i]).collect()
    }
}

/// Carry `loads` (kN, m; the loads of `loads::transfer` without a combination)
/// onto the mesh. `tolerance` is the distance within which a load line or
/// point belongs to an element.
pub fn transfer(mesh: &Mesh, state: &State, loads: &[Load], tolerance: f64) -> MeshLoads {
    let model = &state.model;
    let mut out = MeshLoads::default();
    // Shells by surface, with their centres in the plane of the surface.
    let mut shells_of: BTreeMap<usize, Vec<(usize, DVec2)>> = BTreeMap::new();
    for (i, shell) in mesh.shells.iter().enumerate() {
        let plane = &model.planes()[model.surfaces()[shell.surface].plane];
        let c = shell.nodes.iter().map(|&n| point(mesh, n)).sum::<DVec3>() / shell.nodes.len() as f64;
        shells_of.entry(shell.surface).or_default().push((i, DVec2::from_array(plane.project(c.to_array()))));
    }
    let edge_index = EdgeIndex::new(mesh);
    let index = ShellIndex::new(mesh);
    for load in loads {
        match load {
            Load::Surface { case, surface, polygons, sigma } => {
                let plane = &model.planes()[model.surfaces()[*surface].plane];
                let rings: Vec<Vec<DVec2>> = polygons
                    .iter()
                    .map(|p| p.iter().map(|&x| DVec2::from_array(plane.project(x))).collect())
                    .collect();
                let boxes: Vec<(DVec2, DVec2)> = rings
                    .iter()
                    .map(|r| r.iter().fold((DVec2::splat(f64::MAX), DVec2::splat(f64::MIN)), |(l, h), p| (l.min(*p), h.max(*p))))
                    .collect();
                let region = MultiPolygon::new(
                    rings
                        .iter()
                        .map(|r| Polygon::new(LineString::from(r.iter().map(|p| Coord { x: p.x, y: p.y }).collect::<Vec<_>>()), vec![]))
                        .collect(),
                );
                let all_convex = rings.iter().all(|r| {
                    let n = r.len();
                    let cross = |i: usize| {
                        let (a, b, c) = (r[i], r[(i + 1) % n], r[(i + 2) % n]);
                        (b - a).perp_dot(c - b)
                    };
                    let tol = 1e-12;
                    (0..n).all(|i| cross(i) >= -tol) || (0..n).all(|i| cross(i) <= tol)
                });
                let mut found = false;
                for &(shell, c) in shells_of.get(surface).map(Vec::as_slice).unwrap_or(&[]) {
                    let corners: Vec<DVec2> = mesh.shells[shell]
                        .nodes
                        .iter()
                        .map(|&n| DVec2::from_array(plane.project(mesh.nodes[n])))
                        .collect();
                    // Elements whose box touches no contour box are outside.
                    let (lo, hi) = corners.iter().fold((DVec2::splat(f64::MAX), DVec2::splat(f64::MIN)), |(l, h), p| (l.min(*p), h.max(*p)));
                    if boxes.iter().all(|b| hi.x < b.0.x || lo.x > b.1.x || hi.y < b.0.y || lo.y > b.1.y) {
                        continue;
                    }
                    let flags: Vec<bool> = corners.iter().map(|q| rings.iter().any(|r| inside(r, *q))).collect();
                    let centre_in = rings.iter().any(|r| inside(r, c));
                    // Whole inside (the pieces of the contour tile the region), else the exact part.
                    let mut nodal = false;
                    // The element's own centre of area (the average of the corners is not).
                    let element_centre = Polygon::new(LineString::from(corners.iter().map(|p| Coord { x: p.x, y: p.y }).collect::<Vec<_>>()), vec![])
                        .centroid()
                        .map_or(c, |p| DVec2::new(p.x(), p.y()));
                    // All corners and the centre inside prove full coverage only for convex contours.
                    let fraction = if flags.iter().all(|f| *f) && centre_in && all_convex {
                        1.
                    } else {
                        let element = Polygon::new(LineString::from(corners.iter().map(|p| Coord { x: p.x, y: p.y }).collect::<Vec<_>>()), vec![]);
                        let area = element.unsigned_area();
                        let covered = element.intersection(&region);
                        let fraction = if area > 0. { (covered.unsigned_area() / area).clamp(0., 1.) } else { 0. };
                        // A constant pressure acts at the element centre: when the covered
                        // part lies off it (more than 1 % of the diagonal) (a small patch in a large element), the load
                        // goes to the nodes by the position of its centroid, which keeps the
                        // moment too.
                        if fraction > 1e-9 && fraction < 1. {
                            if let Some(centroid) = covered.centroid() {
                                let diagonal = (hi - lo).length();
                                if (DVec2::new(centroid.x(), centroid.y()) - element_centre).length() > 0.01 * diagonal {
                                    let at = DVec3::from_array(plane.lift([centroid.x(), centroid.y()]));
                                    if let Some(weights) = index.locate(mesh, at, tolerance.max(1e-6)) {
                                        let force = DVec3::from_array(*sigma) * covered.unsigned_area();
                                        for (n, w) in weights {
                                            add_force(&mut out, *case, n, force * w, DVec3::ZERO);
                                        }
                                        nodal = true;
                                        found = true;
                                    }
                                }
                            }
                        }
                        fraction
                    };
                    if fraction > 1e-9 && !nodal {
                        out.pressures.push((*case, shell, (DVec3::from_array(*sigma) * fraction).to_array()));
                        found = true;
                    }
                }
                if !found {
                    let area: f64 = polygons
                        .iter()
                        .map(|p| {
                            let p: Vec<DVec3> = p.iter().map(|&x| DVec3::from_array(x)).collect();
                            (1..p.len() - 1).map(|i| (p[i] - p[0]).cross(p[i + 1] - p[0])).sum::<DVec3>().length() / 2.
                        })
                        .sum();
                    *out.lost.entry(*case).or_default() = (DVec3::from_array(out.lost.get(case).copied().unwrap_or_default()) + DVec3::from_array(*sigma) * area).to_array();
                }
            }
            Load::Line { case, start, end, q_start, q_end } => {
                let (a, b) = (DVec3::from_array(*start), DVec3::from_array(*end));
                let (qa, qb) = (DVec3::from_array(*q_start), DVec3::from_array(*q_end));
                let length = a.distance(b);
                if length < 1e-12 {
                    continue;
                }
                let mut bars = 0;
                for (i, piece) in mesh.bars.iter().enumerate() {
                    let (p, q) = (point(mesh, piece.nodes[0]), point(mesh, piece.nodes[1]));
                    // The piece lies on the line of the load (it may be longer or shorter than the load).
                    let (da, ta) = line_distance(p, a, b);
                    let (db, tb) = line_distance(q, a, b);
                    if da > tolerance || db > tolerance {
                        continue;
                    }
                    // The piece on the load line (unclamped, so the part beyond the line is cut).
                    let line = b - a;
                    let t_of = |x: DVec3| (x - a).dot(line) / line.length_squared();
                    let (t0, t1) = (t_of(p), t_of(q));
                    let _ = (ta, tb);
                    let (lo, hi) = (t0.min(t1).max(0.), t0.max(t1).min(1.));
                    if (hi - lo) * length < 1e-9 || (t1 - t0).abs() < 1e-12 {
                        continue;
                    }
                    // Fractions of the piece (from its first node) of the covered part.
                    let u = |t: f64| (t - t0) / (t1 - t0);
                    let value = |t: f64| qa.lerp(qb, t);
                    let ends = [(u(lo), value(lo)), (u(hi), value(hi))];
                    let [first, second] = if ends[0].0 <= ends[1].0 { ends } else { [ends[1], ends[0]] };
                    out.bar_loads.push(BarLoad { case: *case, bar: i, t: [first.0.clamp(0., 1.), second.0.clamp(0., 1.)], q: [first.1.to_array(), second.1.to_array()] });
                    bars += 1;
                }
                if bars > 0 {
                    continue;
                }
                // Along plate edges (a free edge, a wall standing on a slab, the cut level of a
                // tower): every edge on the line takes the part of the load over it, as
                // consistent nodal forces of the linear load and the linear shape of the edge, which
                // keep its force and its moment. What no edge lies under is lost.
                let on_line = |n: usize| line_distance(point(mesh, n), a, b);
                let mut covered = DVec3::ZERO;
                for [m, n] in edge_index.near(a, b) {
                    let ((dm, tm), (dn, tn)) = (on_line(m), on_line(n));
                    if dm > tolerance || dn > tolerance || (tn - tm).abs() < 1e-12 {
                        continue;
                    }
                    let (lo, hi) = (tm.min(tn).max(0.), tm.max(tn).min(1.));
                    if (hi - lo) * length < 1e-9 {
                        continue;
                    }
                    // Over [lo, hi] the load is qa + (qb - qa) t; shape of node m is (tn - t) / (tn - tm).
                    let q = |t: f64| qa.lerp(qb, t);
                    let (ql, qh) = (q(lo), q(hi));
                    let width = (hi - lo) * length;
                    let force = (ql + qh) / 2. * width;
                    // ∫ q N_m dl with N_m = (tn - t)/(tn - tm): at the ends of the covered part.
                    let (nl, nh) = ((tn - lo) / (tn - tm), (tn - hi) / (tn - tm));
                    let on_m = width * (ql * (2. * nl + nh) + qh * (nl + 2. * nh)) / 6.;
                    add_force(&mut out, *case, m, on_m, DVec3::ZERO);
                    add_force(&mut out, *case, n, force - on_m, DVec3::ZERO);
                    covered += force;
                }
                let resultant = (qa + qb) / 2. * length;
                let missing = resultant - covered;
                if missing.length() > 1e-9 * resultant.length().max(1e-12) {
                    *out.lost.entry(*case).or_default() = (DVec3::from_array(out.lost.get(case).copied().unwrap_or_default()) + missing).to_array();
                }
            }
            Load::Point { case, at, force, moment } => {
                let p = DVec3::from_array(*at);
                let (f, m) = (DVec3::from_array(*force), DVec3::from_array(*moment));
                // On a bar element: split between its nodes by position.
                let mut best: Option<(f64, usize, f64)> = None;
                for (i, piece) in mesh.bars.iter().enumerate() {
                    let (d, t) = segment_distance(p, point(mesh, piece.nodes[0]), point(mesh, piece.nodes[1]));
                    if d <= tolerance && best.is_none_or(|x| d < x.0) {
                        best = Some((d, i, t));
                    }
                }
                if let Some((_, i, t)) = best {
                    let [a, b] = mesh.bars[i].nodes;
                    // Weights (1 - t, t) keep the force and its first moment on the bar.
                    add_force(&mut out, *case, a, f * (1. - t), m * (1. - t));
                    add_force(&mut out, *case, b, f * t, m * t);
                    continue;
                }
                // On a shell element: the weights of the triangle around it keep force and moment.
                if let Some(weights) = index.locate(mesh, p, tolerance.max(1e-6)) {
                    for (n, w) in weights {
                        add_force(&mut out, *case, n, f * w, m * w);
                    }
                    continue;
                }
                // Off the mesh: the nearest node, with the moment of the offset.
                let nearest = (0..mesh.nodes.len())
                    .map(|n| (point(mesh, n).distance(p), n))
                    .min_by(|x, y| x.0.total_cmp(&y.0));
                match nearest {
                    Some((d, n)) if d <= tolerance.max(1e-9) * 1000. => add_force(&mut out, *case, n, f, m + (p - point(mesh, n)).cross(f)),
                    _ => *out.lost.entry(*case).or_default() = (DVec3::from_array(out.lost.get(case).copied().unwrap_or_default()) + f).to_array(),
                }
            }
        }
    }
    out
}
