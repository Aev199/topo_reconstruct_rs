//! Loads of the LIRA model carried onto a Gmsh mesh: every load case stays a
//! case, a pressure goes to the shell elements whose centres are inside its
//! contour (the contour is that of the source), line loads to the bar
//! elements along them (or, along a plate edge, to its boundary nodes), point
//! loads to the nodes of the element they act on. The resultant of every case
//! is kept; the report compares it with that of the loads on the geometry.
use crate::loads::Load;
use crate::meshing::Mesh;
use crate::reconstruction::assembly::edit::State;
use geo::{Area, BooleanOps, Coord, LineString, MultiPolygon, Polygon};
use glam::{DVec2, DVec3};
use serde::Serialize;
use std::collections::BTreeMap;

/// Forces in kN, global axes.
#[derive(Debug, Clone, Default, Serialize)]
pub struct MeshLoads {
    /// (case, shell, kN/m2 vector).
    pub pressures: Vec<(u32, usize, [f64; 3])>,
    /// (case, bar piece, uniform kN/m vector).
    pub bar_loads: Vec<(u32, usize, [f64; 3])>,
    /// (case, node) -> force and moment.
    pub nodal: BTreeMap<(u32, usize), [f64; 6]>,
    /// Loads that found no element (point loads far from the mesh, pressure
    /// contours without shells).
    pub lost: BTreeMap<u32, [f64; 3]>,
}

impl MeshLoads {
    /// Resultant force per case (kN).
    pub fn resultants(&self, mesh: &Mesh) -> BTreeMap<u32, DVec3> {
        let mut out: BTreeMap<u32, DVec3> = BTreeMap::new();
        for &(case, shell, p) in &self.pressures {
            *out.entry(case).or_default() += DVec3::from_array(p) * shell_area(mesh, shell);
        }
        for &(case, bar, q) in &self.bar_loads {
            *out.entry(case).or_default() += DVec3::from_array(q) * bar_length(mesh, bar);
        }
        for (&(case, _), f) in &self.nodal {
            *out.entry(case).or_default() += DVec3::new(f[0], f[1], f[2]);
        }
        out
    }
}

fn point(mesh: &Mesh, n: usize) -> DVec3 {
    DVec3::from_array(mesh.nodes[n])
}

pub fn shell_area(mesh: &Mesh, shell: usize) -> f64 {
    let p: Vec<DVec3> = mesh.shells[shell].nodes.iter().map(|&n| point(mesh, n)).collect();
    (1..p.len() - 1)
        .map(|i| (p[i] - p[0]).cross(p[i + 1] - p[0]).length() / 2.)
        .sum()
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
    // Bars and shell edges by node, for finding the elements along a line.
    let mut edges_on_boundary: Vec<[usize; 2]> = vec![];
    {
        let mut count: BTreeMap<(usize, usize), usize> = BTreeMap::new();
        for s in &mesh.shells {
            for i in 0..s.nodes.len() {
                let (a, b) = (s.nodes[i], s.nodes[(i + 1) % s.nodes.len()]);
                *count.entry((a.min(b), a.max(b))).or_default() += 1;
            }
        }
        edges_on_boundary.extend(count.into_iter().filter(|(_, n)| *n == 1).map(|((a, b), _)| [a, b]));
    }
    for load in loads {
        match load {
            Load::Surface { case, surface, polygons, sigma } => {
                let plane = &model.planes()[model.surfaces()[*surface].plane];
                let rings: Vec<Vec<DVec2>> = polygons
                    .iter()
                    .map(|p| p.iter().map(|&x| DVec2::from_array(plane.project(x))).collect())
                    .collect();
                let mut found = false;
                let region = MultiPolygon::new(
                    rings
                        .iter()
                        .map(|r| Polygon::new(LineString::from(r.iter().map(|p| Coord { x: p.x, y: p.y }).collect::<Vec<_>>()), vec![]))
                        .collect(),
                );
                for &(shell, c) in shells_of.get(surface).map(Vec::as_slice).unwrap_or(&[]) {
                    let corners: Vec<DVec2> = mesh.shells[shell]
                        .nodes
                        .iter()
                        .map(|&n| DVec2::from_array(plane.project(mesh.nodes[n])))
                        .collect();
                    let flags: Vec<bool> = corners.iter().map(|q| rings.iter().any(|r| inside(r, *q))).collect();
                    let centre_in = rings.iter().any(|r| inside(r, c));
                    // Whole inside, or the part of the element in the contour.
                    let fraction = if flags.iter().all(|f| *f) && centre_in {
                        1.
                    } else if !centre_in && flags.iter().all(|f| !*f) {
                        // Outside (a sliver of the contour poking into the element is ignored).
                        0.
                    } else {
                        let element = Polygon::new(LineString::from(corners.iter().map(|p| Coord { x: p.x, y: p.y }).collect::<Vec<_>>()), vec![]);
                        let area = element.unsigned_area();
                        if area > 0. { (element.intersection(&region).unsigned_area() / area).clamp(0., 1.) } else { 0. }
                    };
                    if fraction > 1e-6 {
                        let p = DVec3::from_array(*sigma) * fraction;
                        out.pressures.push((*case, shell, p.to_array()));
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
                let mut bars = 0;
                for (i, piece) in mesh.bars.iter().enumerate() {
                    let (p, q) = (point(mesh, piece.nodes[0]), point(mesh, piece.nodes[1]));
                    let (da, ta) = segment_distance(p, a, b);
                    let (db, tb) = segment_distance(q, a, b);
                    if da <= tolerance && db <= tolerance && ((ta - tb).abs() * length) > tolerance {
                        // Uniform on the piece, with the value at its middle.
                        let t = (ta + tb) / 2.;
                        out.bar_loads.push((*case, i, (qa + (qb - qa) * t).to_array()));
                        bars += 1;
                    }
                }
                if bars > 0 {
                    continue;
                }
                // Along a plate edge: the boundary nodes on the line take their share.
                let mut chain: Vec<(f64, usize)> = vec![];
                for e in &edges_on_boundary {
                    for &n in e {
                        let (d, t) = segment_distance(point(mesh, n), a, b);
                        if d <= tolerance && !chain.iter().any(|x| x.1 == n) {
                            chain.push((t, n));
                        }
                    }
                }
                chain.sort_by(|x, y| x.0.total_cmp(&y.0));
                if chain.len() < 2 {
                    let f = (qa + qb) / 2. * length;
                    *out.lost.entry(*case).or_default() = (DVec3::from_array(out.lost.get(case).copied().unwrap_or_default()) + f).to_array();
                    continue;
                }
                // Trapezoidal weights of the nodes along the line, scaled to the resultant.
                let resultant = (qa + qb) / 2. * length;
                let mut shares: Vec<(usize, DVec3)> = vec![];
                for w in chain.windows(2) {
                    let (t0, t1) = (w[0].0, w[1].0);
                    let q = qa + (qb - qa) * ((t0 + t1) / 2.);
                    let f = q * ((t1 - t0) * length);
                    shares.push((w[0].1, f / 2.));
                    shares.push((w[1].1, f / 2.));
                }
                let total: DVec3 = shares.iter().map(|s| s.1).sum();
                // The chain may not span the whole line: keep the resultant.
                let scale = if total.length() > 1e-12 { resultant.length() / total.length() } else { 1. };
                for (n, f) in shares {
                    add_force(&mut out, *case, n, f * scale, DVec3::ZERO);
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
                    add_force(&mut out, *case, a, f * (1. - t), m * (1. - t));
                    add_force(&mut out, *case, b, f * t, m * t);
                    continue;
                }
                // The nearest node (a point load is placed at a node of the geometry or on a plate).
                let nearest = (0..mesh.nodes.len())
                    .map(|n| (point(mesh, n).distance(p), n))
                    .min_by(|x, y| x.0.total_cmp(&y.0));
                match nearest {
                    Some((d, n)) if d <= tolerance.max(1e-9) * 1000. => add_force(&mut out, *case, n, f, m),
                    _ => *out.lost.entry(*case).or_default() = (DVec3::from_array(out.lost.get(case).copied().unwrap_or_default()) + f).to_array(),
                }
            }
        }
    }
    out
}
