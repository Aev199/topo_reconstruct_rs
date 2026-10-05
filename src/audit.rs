//! Read-only audit of an assembled geotechnical geometry, for the editor.
//!
//! A port of the independent Python auditors (`scripts/check_v2_global_geometry.py`
//! without its trial-mesh checks, and `scripts/check_plaxis_profile.py`) that
//! runs inside the application after every edit. It reads only the assembled
//! model (vertices, edges, surface contours and edge uses, bar axes and their
//! contact records) and shares no geometric code with the reconstruction:
//! a defect the reconstruction cannot see should not be hidden from its
//! audit by the same helper. The Python auditors stay the second, independent
//! check of saved results.
//!
//! Findings are failures (the strict global audit), review items (near
//! misses, never failing) and PLAXIS items (features far below the target
//! element size).
use crate::reconstruction::assembly::bars::{Axis, Contact};
use crate::reconstruction::Model;
use glam::{DVec2, DVec3};
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet, HashMap};

/// Coarsest numerical precision the audit accepts (model units, 1 um).
pub const MAXIMUM_PRECISION: f64 = 1e-6;

#[derive(Debug, Clone, Serialize, serde::Deserialize, PartialEq)]
pub struct Options {
    /// Near misses within this distance are review items.
    pub near_distance: f64,
    /// Target element size h; PLAXIS items are measured against h/10.
    pub element_size: f64,
    /// Sharp contour corners below this angle (degrees) are PLAXIS items.
    pub sharp_angle: f64,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            near_distance: 0.05,
            element_size: 0.5,
            sharp_angle: 10.,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Class {
    Failure,
    Plaxis,
    Review,
}

#[derive(Debug, Clone, Serialize, PartialEq)]
pub struct Finding {
    pub kind: String,
    pub class: Class,
    pub surfaces: Vec<usize>,
    pub bars: Vec<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub vertex: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub edge: Option<usize>,
    /// Where it is: one point, or the ends of a segment.
    pub points: Vec<[f64; 3]>,
    /// Length, area, distance or angle (degrees), by kind.
    pub value: f64,
    #[serde(skip_serializing_if = "String::is_empty")]
    pub detail: String,
}

impl Finding {
    fn new(kind: &str, class: Class) -> Self {
        Finding {
            kind: kind.into(),
            class,
            surfaces: vec![],
            bars: vec![],
            vertex: None,
            edge: None,
            points: vec![],
            value: 0.,
            detail: String::new(),
        }
    }
    fn surfaces(mut self, s: impl IntoIterator<Item = usize>) -> Self {
        self.surfaces = s.into_iter().collect();
        self
    }
    fn bars(mut self, b: impl IntoIterator<Item = usize>) -> Self {
        self.bars = b.into_iter().collect();
        self
    }
    fn vertex(mut self, v: usize) -> Self {
        self.vertex = Some(v);
        self
    }
    fn edge(mut self, e: usize) -> Self {
        self.edge = Some(e);
        self
    }
    fn at(mut self, p: DVec3) -> Self {
        self.points.push(p.to_array());
        self
    }
    fn value(mut self, v: f64) -> Self {
        self.value = v;
        self
    }
    /// Identity of the site, independent of the order of findings: kind,
    /// objects and position at millimetre resolution.
    pub fn key(&self) -> String {
        let mm = |p: &[f64; 3]| p.map(|x| (x * 1000.).round() as i64);
        format!(
            "{}|s{:?}|b{:?}|v{:?}|e{:?}|{:?}",
            self.kind,
            self.surfaces,
            self.bars,
            self.vertex,
            self.edge,
            self.points.first().map(mm)
        )
    }
}

#[derive(Debug, Clone, Default, Serialize)]
pub struct Report {
    /// Audit tolerance (5 x the precision, at most 5 um).
    pub precision: f64,
    pub options: Option<Options>,
    pub findings: Vec<Finding>,
    pub counts: BTreeMap<String, usize>,
    /// No failure.
    pub passed: bool,
    /// No failure and no PLAXIS item.
    pub plaxis_passed: bool,
    /// Proximities explained by a corner of the geometry (not gaps).
    pub explained_proximities: usize,
}

// ---------------------------------------------------------------- geometry

fn p3(p: [f64; 3]) -> DVec3 {
    DVec3::from_array(p)
}

fn closest_on(p: DVec2, a: DVec2, b: DVec2) -> DVec2 {
    let d = b - a;
    let l = d.length_squared();
    let t = if l > 0. {
        ((p - a).dot(d) / l).clamp(0., 1.)
    } else {
        0.
    };
    a + d * t
}

fn closest_on3(p: DVec3, a: DVec3, b: DVec3) -> DVec3 {
    let d = b - a;
    let l = d.length_squared();
    let t = if l > 0. {
        ((p - a).dot(d) / l).clamp(0., 1.)
    } else {
        0.
    };
    a + d * t
}

/// Parameters and distance of the closest points of segments `p` and `q`.
fn closest_points(p0: DVec3, p1: DVec3, q0: DVec3, q1: DVec3) -> (f64, DVec3) {
    let (d1, d2, r) = (p1 - p0, q1 - q0, p0 - q0);
    let (a, e, f) = (d1.dot(d1), d2.dot(d2), d2.dot(r));
    let (c, b) = (d1.dot(r), d1.dot(d2));
    let denom = a * e - b * b;
    let mut s = if denom > 1e-18 * a * e && denom > 0. {
        ((b * f - c * e) / denom).clamp(0., 1.)
    } else {
        0.
    };
    let mut t = if e > 0. { (b * s + f) / e } else { 0. };
    if t < 0. {
        t = 0.;
        s = if a > 0. { (-c / a).clamp(0., 1.) } else { 0. };
    } else if t > 1. {
        t = 1.;
        s = if a > 0. {
            ((b - c) / a).clamp(0., 1.)
        } else {
            0.
        };
    }
    let (x, y) = (p0 + d1 * s, q0 + d2 * t);
    (x.distance(y), x)
}

/// Polygon with holes in plane coordinates.
struct Polygon {
    rings: Vec<Vec<DVec2>>,
}

impl Polygon {
    fn edges(&self) -> impl Iterator<Item = (DVec2, DVec2)> + '_ {
        self.rings
            .iter()
            .flat_map(|r| (0..r.len()).map(move |i| (r[i], r[(i + 1) % r.len()])))
    }
    fn boundary_distance(&self, p: DVec2) -> f64 {
        self.edges()
            .map(|(a, b)| p.distance(closest_on(p, a, b)))
            .fold(f64::MAX, f64::min)
    }
    /// Even-odd membership (exterior and holes together).
    fn inside(&self, p: DVec2) -> bool {
        let mut odd = false;
        for (a, b) in self.edges() {
            if (a.y > p.y) != (b.y > p.y) && p.x < a.x + (p.y - a.y) / (b.y - a.y) * (b.x - a.x) {
                odd = !odd;
            }
        }
        odd
    }
    /// Distance from the closed material (0 inside or on the boundary).
    fn distance(&self, p: DVec2) -> f64 {
        if self.inside(p) {
            0.
        } else {
            self.boundary_distance(p)
        }
    }
    fn contains_closed(&self, p: DVec2, eps: f64) -> bool {
        self.inside(p) || self.boundary_distance(p) <= eps
    }
    fn area(&self) -> f64 {
        let ring_area = |r: &Vec<DVec2>| {
            (0..r.len())
                .map(|i| r[i].perp_dot(r[(i + 1) % r.len()]))
                .sum::<f64>()
                .abs()
                / 2.
        };
        self.rings.first().map_or(0., ring_area)
            - self.rings.iter().skip(1).map(ring_area).sum::<f64>()
    }
    fn perimeter(&self) -> f64 {
        self.edges().map(|(a, b)| a.distance(b)).sum()
    }
    /// Parameter intervals of segment `a b` inside the closed material.
    fn clip(&self, a: DVec2, b: DVec2, eps: f64) -> Vec<(f64, f64)> {
        let d = b - a;
        let length = d.length();
        if length == 0. {
            return vec![];
        }
        let mut ts = vec![0., 1.];
        for (c, e) in self.edges() {
            let f = e - c;
            let den = d.perp_dot(f);
            if den.abs() > 1e-300 {
                let t = (c - a).perp_dot(f) / den;
                let u = (c - a).perp_dot(d) / den;
                if (-1e-12..=1. + 1e-12).contains(&u) && (0. ..=1.).contains(&t) {
                    ts.push(t);
                }
            }
            // Ends of collinear or touching edges.
            for q in [c, e] {
                let t = (q - a).dot(d) / (length * length);
                if (0. ..=1.).contains(&t) && (a + d * t).distance(q) <= eps {
                    ts.push(t);
                }
            }
        }
        ts.sort_by(f64::total_cmp);
        ts.dedup_by(|x, y| (*x - *y).abs() * length <= 1e-12);
        let mut out: Vec<(f64, f64)> = vec![];
        for w in ts.windows(2) {
            let mid = a + d * ((w[0] + w[1]) / 2.);
            if self.contains_closed(mid, eps) {
                match out.last_mut() {
                    Some(last) if (last.1 - w[0]).abs() * length <= 1e-12 => last.1 = w[1],
                    _ => out.push((w[0], w[1])),
                }
            }
        }
        // A touch at a single point (an end on the boundary).
        if out.is_empty() {
            for t in [0., 1.] {
                if self.contains_closed(a + d * t, eps) {
                    out.push((t, t));
                }
            }
        }
        out
    }
}

/// Whether segments `a b` and `c d` intersect (touching included).
fn segments_meet(a: DVec2, b: DVec2, c: DVec2, d: DVec2, eps: f64) -> bool {
    let side = |p: DVec2, q: DVec2, r: DVec2| (q - p).perp_dot(r - p);
    let (d1, d2, d3, d4) = (side(c, d, a), side(c, d, b), side(a, b, c), side(a, b, d));
    if ((d1 > 0. && d2 < 0.) || (d1 < 0. && d2 > 0.))
        && ((d3 > 0. && d4 < 0.) || (d3 < 0. && d4 > 0.))
    {
        return true;
    }
    a.distance(closest_on(a, c, d)) <= eps
        || b.distance(closest_on(b, c, d)) <= eps
        || c.distance(closest_on(c, a, b)) <= eps
        || d.distance(closest_on(d, a, b)) <= eps
}

/// Why a polygon is not a valid area: self-touching rings, holes outside
/// the exterior or meeting another ring, too small.
fn invalidity(polygon: &Polygon, eps: f64) -> Option<String> {
    let mut segments = vec![];
    for (r, ring) in polygon.rings.iter().enumerate() {
        if ring.len() < 3 {
            return Some(format!("ring {r} has fewer than three vertices"));
        }
        for i in 0..ring.len() {
            let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
            if a.distance(b) <= eps {
                return Some(format!("ring {r} has a zero-length edge"));
            }
            segments.push((r, i, ring.len(), a, b));
        }
    }
    // Sweep over x: any two segments meeting, except neighbours of one ring
    // at their common vertex.
    let mut order: Vec<usize> = (0..segments.len()).collect();
    order.sort_by(|&i, &j| {
        let lo = |k: usize| segments[k].3.x.min(segments[k].4.x);
        lo(i).total_cmp(&lo(j))
    });
    let mut active: Vec<usize> = vec![];
    for &k in &order {
        let (r, i, n, a, b) = segments[k];
        let lo = a.x.min(b.x);
        active.retain(|&m| segments[m].3.x.max(segments[m].4.x) + eps >= lo);
        for &m in &active {
            let (r2, j, _, c, d) = segments[m];
            let neighbours = r == r2 && ((i + 1) % n == j || (j + 1) % n == i);
            if neighbours {
                // Only their common vertex may touch: a fold back onto the
                // other edge is a self-touch.
                let (p, q) = if (i + 1) % n == j { (a, d) } else { (b, c) };
                if p.distance(closest_on(p, c, d)) <= eps && p != c && p != d
                    || q.distance(closest_on(q, a, b)) <= eps && q != a && q != b
                {
                    return Some(format!("ring {r} folds back at vertex {}", (i + 1) % n));
                }
                continue;
            }
            if segments_meet(a, b, c, d, eps) {
                return Some(if r == r2 {
                    format!("ring {r} touches itself")
                } else {
                    format!("rings {r2} and {r} meet")
                });
            }
        }
        active.push(k);
    }
    let exterior = Polygon {
        rings: vec![polygon.rings[0].clone()],
    };
    for (r, hole) in polygon.rings.iter().enumerate().skip(1) {
        if !exterior.inside(hole[0]) {
            return Some(format!("hole {r} lies outside the exterior"));
        }
        for (q, other) in polygon.rings.iter().enumerate().skip(1) {
            if q != r
                && (Polygon {
                    rings: vec![other.clone()],
                })
                .inside(hole[0])
            {
                return Some(format!("hole {r} lies inside hole {q}"));
            }
        }
    }
    if polygon.area() <= eps * eps {
        return Some("no area".into());
    }
    None
}

/// Whether sorted-or-not intervals cover [start, end] within eps.
fn covers(parts: &mut [(f64, f64)], start: f64, end: f64, eps: f64) -> bool {
    parts.sort_by(|x, y| x.0.total_cmp(&y.0));
    let mut cursor = start;
    for &(a, b) in parts.iter() {
        if b < cursor - eps {
            continue;
        }
        if a > cursor + eps {
            return false;
        }
        cursor = cursor.max(b);
        if cursor >= end - eps {
            return true;
        }
    }
    cursor >= end - eps
}

/// Intervals along `start + t d` of the segments lying on that line.
fn on_line(
    pairs: impl Iterator<Item = [usize; 2]>,
    vertices: &[[f64; 3]],
    start: DVec3,
    d: DVec3,
    eps: f64,
) -> Vec<(f64, f64)> {
    pairs
        .filter_map(|[a, b]| {
            let (pa, pb) = (p3(vertices[a]) - start, p3(vertices[b]) - start);
            let (ta, tb) = (pa.dot(d), pb.dot(d));
            ((pa - d * ta).length() <= eps && (pb - d * tb).length() <= eps)
                .then(|| (ta.min(tb), ta.max(tb)))
        })
        .collect()
}

/// Uniform grid over bounding boxes.
struct Grid {
    cell: f64,
    cells: HashMap<[i64; 3], Vec<usize>>,
}

impl Grid {
    fn new(cell: f64) -> Self {
        Grid {
            cell,
            cells: HashMap::new(),
        }
    }
    fn range(&self, lo: DVec3, hi: DVec3) -> ([i64; 3], [i64; 3]) {
        let k = |p: DVec3| (p / self.cell).floor().as_i64vec3().to_array();
        (k(lo), k(hi))
    }
    fn insert(&mut self, item: usize, lo: DVec3, hi: DVec3) {
        let (a, b) = self.range(lo, hi);
        for x in a[0]..=b[0] {
            for y in a[1]..=b[1] {
                for z in a[2]..=b[2] {
                    self.cells.entry([x, y, z]).or_default().push(item);
                }
            }
        }
    }
    fn query(&self, lo: DVec3, hi: DVec3) -> BTreeSet<usize> {
        let (a, b) = self.range(lo, hi);
        let mut out = BTreeSet::new();
        for x in a[0]..=b[0] {
            for y in a[1]..=b[1] {
                for z in a[2]..=b[2] {
                    if let Some(items) = self.cells.get(&[x, y, z]) {
                        out.extend(items.iter().copied());
                    }
                }
            }
        }
        out
    }
}

// ------------------------------------------------------------------ audit

struct Surface {
    origin: DVec3,
    normal: DVec3,
    frame: crate::reconstruction::PlaneFrame,
    polygon: Polygon,
    rings: Vec<Vec<usize>>,
    edge_ids: BTreeSet<usize>,
    vertex_ids: BTreeSet<usize>,
    lo: DVec3,
    hi: DVec3,
}

impl Surface {
    fn project(&self, p: DVec3) -> DVec2 {
        DVec2::from_array(self.frame.project(p.to_array()))
    }
    fn lift(&self, uv: DVec2) -> DVec3 {
        p3(self.frame.lift(uv.to_array()))
    }
    fn height(&self, p: DVec3) -> f64 {
        (p - self.origin).dot(self.normal)
    }
    /// Distance of a 3D point from the closed surface.
    fn distance(&self, p: DVec3) -> f64 {
        self.height(p).hypot(self.polygon.distance(self.project(p)))
    }
}

/// Audit the model with its bar axes and bar-surface contact records.
pub fn run(model: &Model, axes: &[Axis], contacts: &[Contact], options: &Options) -> Report {
    let eps = 5. * model.precision().min(MAXIMUM_PRECISION);
    let near = options.near_distance.max(eps);
    let vertices = model.vertices();
    let edges = model.edges();
    let mut findings = vec![];
    // Bar representation: straight bars through ordered, distinct nodes.
    for (i, defect) in crate::reconstruction::assembly::bars::axis_defects(model, axes) {
        let mut f = Finding::new("broken_bar", Class::Failure).bars([i]);
        if let Some(axis) = axes.get(i) {
            f = f.at(p3(vertices[axis.endpoints[0]]));
        }
        f.detail = defect;
        findings.push(f);
    }
    let surfaces: Vec<Surface> = model
        .surfaces()
        .iter()
        .map(|s| {
            let frame = model.planes()[s.plane].clone();
            let rings: Vec<Vec<usize>> = s
                .boundaries
                .iter()
                .map(|ring| {
                    ring.iter()
                        .map(|u| {
                            let [a, b] = edges[u.edge];
                            if u.reversed {
                                b
                            } else {
                                a
                            }
                        })
                        .collect()
                })
                .collect();
            let edge_ids: BTreeSet<usize> = s
                .boundaries
                .iter()
                .flatten()
                .map(|u| u.edge)
                .chain(s.embedded_edges.iter().copied())
                .collect();
            let vertex_ids: BTreeSet<usize> = edge_ids.iter().flat_map(|&e| edges[e]).collect();
            let (lo, hi) = vertex_ids.iter().fold(
                (DVec3::splat(f64::MAX), DVec3::splat(f64::MIN)),
                |(lo, hi), &v| (lo.min(p3(vertices[v])), hi.max(p3(vertices[v]))),
            );
            // Contours from the vertices themselves, not the stored plane
            // coordinates (which the reconstruction maintains).
            let polygon = Polygon {
                rings: rings
                    .iter()
                    .map(|r| {
                        r.iter()
                            .map(|&v| DVec2::from_array(frame.project(vertices[v])))
                            .collect()
                    })
                    .collect(),
            };
            Surface {
                origin: p3(frame.origin()),
                normal: p3(frame.normal()),
                frame,
                polygon,
                rings,
                edge_ids,
                vertex_ids,
                lo,
                hi,
            }
        })
        .collect();

    // Individual surfaces.
    let mut invalid = BTreeSet::new();
    for (i, s) in surfaces.iter().enumerate() {
        let planarity = s
            .vertex_ids
            .iter()
            .map(|&v| s.height(p3(vertices[v])).abs())
            .fold(0., f64::max);
        let boundary: BTreeSet<usize> = s
            .rings
            .iter()
            .enumerate()
            .flat_map(|(r, _)| model.surfaces()[i].boundaries[r].iter().map(|u| u.edge))
            .collect();
        let reason = if planarity > eps {
            Some(format!("not planar ({planarity:.3e})"))
        } else if let Some(reason) = invalidity(&s.polygon, eps) {
            Some(reason)
        } else {
            model.surfaces()[i].embedded_edges.iter().find_map(|&e| {
                if boundary.contains(&e) {
                    return Some(format!("embedded edge {e} is a contour edge"));
                }
                let [a, b] = edges[e].map(|v| s.project(p3(vertices[v])));
                let inside = s.polygon.clip(a, b, eps);
                let covered = inside.len() == 1 && inside[0].0 <= 1e-9 && inside[0].1 >= 1. - 1e-9;
                (a.distance(b) <= eps || !covered)
                    .then(|| format!("embedded edge {e} leaves the surface"))
            })
        };
        if let Some(reason) = reason {
            invalid.insert(i);
            let mut f = Finding::new("invalid_surface", Class::Failure)
                .surfaces([i])
                .at((s.lo + s.hi) / 2.)
                .value(s.polygon.area());
            f.detail = reason;
            findings.push(f);
        }
    }

    // Surface pairs.
    let reach = DVec3::splat(near);
    let mut order: Vec<usize> = (0..surfaces.len())
        .filter(|i| !invalid.contains(i))
        .collect();
    order.sort_by(|&a, &b| surfaces[a].lo.x.total_cmp(&surfaces[b].lo.x));
    for (k, &i) in order.iter().enumerate() {
        let a = &surfaces[i];
        for &j in &order[k + 1..] {
            let b = &surfaces[j];
            if b.lo.x > a.hi.x + near {
                break;
            }
            if (a.hi + reach).cmplt(b.lo).any() || (b.hi + reach).cmplt(a.lo).any() {
                continue;
            }
            let pair = [i.min(j), i.max(j)];
            let direction = a.normal.cross(b.normal);
            let sine = direction.length();
            let mut segments: Vec<(DVec3, DVec3)> = vec![];
            if sine < 1e-8 {
                let gap = (b.origin - a.origin).dot(a.normal).abs();
                if gap > near {
                    continue;
                }
                let other = Polygon {
                    rings: b
                        .rings
                        .iter()
                        .map(|r| r.iter().map(|&v| a.project(p3(vertices[v]))).collect())
                        .collect(),
                };
                let common = overlap_area(&a.polygon, &other);
                let threshold = eps * a.polygon.perimeter().min(other.perimeter()).max(eps);
                if common > threshold {
                    let class = if gap <= eps {
                        ("coplanar_overlap", Class::Failure)
                    } else {
                        ("near_parallel_faces", Class::Review)
                    };
                    findings.push(
                        Finding::new(class.0, class.1)
                            .surfaces(pair)
                            .at(centre(a, &other))
                            .value(common),
                    );
                    continue;
                }
                if gap > eps {
                    continue;
                }
                // Coplanar boundaries touching along lines.
                for (p, q) in a.polygon.edges() {
                    for (r, t) in other.edges() {
                        if let Some((x, y)) = collinear_overlap(p, q, r, t, eps) {
                            segments.push((a.lift(x), a.lift(y)));
                        }
                    }
                }
            } else {
                let d = direction / sine;
                let offset = (b.origin - a.origin).dot(b.normal);
                let origin = a.origin + d.cross(a.normal) * (offset / sine);
                let span = [a.lo, a.hi, b.lo, b.hi]
                    .iter()
                    .map(|p| p.distance(origin))
                    .fold(0., f64::max)
                    + 1.;
                let (s0, s1) = (origin - d * span, origin + d * span);
                let along = |s: &Surface| -> Vec<(f64, f64)> {
                    s.polygon
                        .clip(s.project(s0), s.project(s1), eps)
                        .into_iter()
                        .map(|(t0, t1)| (-span + 2. * span * t0, -span + 2. * span * t1))
                        .collect()
                };
                for (x0, x1) in along(a) {
                    for (y0, y1) in along(b) {
                        let (start, end) = (x0.max(y0), x1.min(y1));
                        if end - start > eps {
                            segments.push((origin + d * start, origin + d * end));
                        }
                    }
                }
            }
            let shared: Vec<[usize; 2]> = a
                .edge_ids
                .intersection(&b.edge_ids)
                .map(|&e| edges[e])
                .collect();
            for (start, end) in segments {
                let length = start.distance(end);
                if length <= eps {
                    continue;
                }
                let d = (end - start) / length;
                let mut parts = on_line(shared.iter().copied(), vertices, start, d, eps);
                if !covers(&mut parts, 0., length, eps) {
                    let on_boundary = |s: &Surface| {
                        (0..=4).all(|k| {
                            let p = start + (end - start) * (k as f64 / 4.);
                            s.polygon.boundary_distance(s.project(p)) <= eps
                        })
                    };
                    let kind = match (on_boundary(a), on_boundary(b)) {
                        (true, true) => "boundary_junction",
                        (false, false) => "crossing",
                        _ => "t_junction",
                    };
                    let mut f = Finding::new("unrepresented_intersection", Class::Failure)
                        .surfaces(pair)
                        .at(start)
                        .at(end)
                        .value(length);
                    f.kind = format!("unrepresented_{kind}");
                    findings.push(f);
                }
            }
        }
    }

    // Point contacts and PLAXIS gaps: surface vertices and bar nodes near
    // another surface they are not a vertex of.
    let mut owners = BTreeMap::<usize, Vec<usize>>::new();
    for (i, s) in surfaces.iter().enumerate() {
        for &v in &s.vertex_ids {
            owners.entry(v).or_default().push(i);
        }
    }
    let bar_nodes: BTreeSet<usize> = axes
        .iter()
        .flat_map(|a| {
            a.endpoints
                .iter()
                .copied()
                .chain(a.anchors.iter().map(|x| x.vertex))
        })
        .collect();
    let element_tenth = options.element_size / 10.;
    let gap_reach = near.max(element_tenth);
    let mut vertex_grid = Grid::new(gap_reach.max(1e-9) * 4.);
    let candidates: BTreeSet<usize> = owners
        .keys()
        .copied()
        .chain(bar_nodes.iter().copied())
        .collect();
    for &v in &candidates {
        let p = p3(vertices[v]);
        vertex_grid.insert(v, p, p);
    }
    // Adjacency for explained proximities: model edges and bar pieces.
    let used_edges: BTreeSet<usize> = surfaces
        .iter()
        .flat_map(|s| s.edge_ids.iter().copied())
        .collect();
    let mut adjacent = BTreeMap::<usize, BTreeSet<usize>>::new();
    for &e in &used_edges {
        let [a, b] = edges[e];
        adjacent.entry(a).or_default().insert(b);
        adjacent.entry(b).or_default().insert(a);
    }
    let pieces = bar_pieces(axes);
    for &(_, a, b) in &pieces {
        adjacent.entry(a).or_default().insert(b);
        adjacent.entry(b).or_default().insert(a);
    }
    let sin_angle = options.sharp_angle.to_radians().sin();
    let mut explained = 0;
    for (j, b) in surfaces.iter().enumerate() {
        if invalid.contains(&j) {
            continue;
        }
        let r = DVec3::splat(gap_reach);
        for v in vertex_grid.query(b.lo - r, b.hi + r) {
            if b.vertex_ids.contains(&v) {
                continue;
            }
            let p = p3(vertices[v]);
            let d = b.distance(p);
            let of: Vec<usize> = owners.get(&v).cloned().unwrap_or_default();
            if of.iter().any(|o| invalid.contains(o)) {
                continue;
            }
            if !of.is_empty() && d <= eps {
                findings.push(
                    Finding::new("unshared_point_contact", Class::Failure)
                        .surfaces(of.iter().copied().chain([j]))
                        .vertex(v)
                        .at(p)
                        .value(d),
                );
                continue;
            }
            if !of.is_empty() && d <= near {
                findings.push(
                    Finding::new("surface_near_miss", Class::Review)
                        .surfaces(of.iter().copied().chain([j]))
                        .vertex(v)
                        .at(p)
                        .value(d),
                );
            }
            if eps < d && d < element_tenth {
                let corner = adjacent.get(&v).is_some_and(|ws| {
                    ws.iter().any(|&w| {
                        let length = p.distance(p3(vertices[w]));
                        b.vertex_ids.contains(&w)
                            && length >= element_tenth
                            && d >= length * sin_angle
                    })
                });
                if corner {
                    explained += 1;
                } else {
                    findings.push(
                        Finding::new("gap", Class::Plaxis)
                            .surfaces(of.iter().copied().chain([j]))
                            .vertex(v)
                            .at(p)
                            .value(d),
                    );
                }
            }
        }
    }

    audit_bars(
        &surfaces,
        &invalid,
        axes,
        contacts,
        edges,
        vertices,
        eps,
        near,
        &mut findings,
    );
    plaxis_shapes(
        &surfaces,
        &used_edges,
        &pieces,
        edges,
        vertices,
        options,
        &mut findings,
    );

    findings.sort_by(|x, y| {
        x.class
            .cmp(&y.class)
            .then(x.kind.cmp(&y.kind))
            .then(x.surfaces.cmp(&y.surfaces))
            .then(x.bars.cmp(&y.bars))
            .then(x.vertex.cmp(&y.vertex))
            .then(x.edge.cmp(&y.edge))
    });
    let mut counts = BTreeMap::new();
    for f in &findings {
        *counts.entry(f.kind.clone()).or_default() += 1;
    }
    let passed = !findings.iter().any(|f| f.class == Class::Failure);
    Report {
        precision: eps,
        options: Some(options.clone()),
        plaxis_passed: passed && !findings.iter().any(|f| f.class == Class::Plaxis),
        passed,
        counts,
        findings,
        explained_proximities: explained,
    }
}

fn centre(a: &Surface, other: &Polygon) -> DVec3 {
    // A point of the overlap: the first vertex of `other` inside `a`, else
    // the first vertex of `a`.
    other
        .rings
        .first()
        .and_then(|r| r.iter().find(|&&p| a.polygon.inside(p)).copied())
        .map_or(a.lift(a.polygon.rings[0][0]), |p| a.lift(p))
}

/// Common area of two polygons (0 if the boolean operation fails).
fn overlap_area(a: &Polygon, b: &Polygon) -> f64 {
    use geo::{Area, BooleanOps};
    let polygon = |p: &Polygon| {
        let ring = |r: &Vec<DVec2>| {
            geo::LineString::from(
                r.iter()
                    .chain(r.first())
                    .map(|q| (q.x, q.y))
                    .collect::<Vec<_>>(),
            )
        };
        geo::Polygon::new(ring(&p.rings[0]), p.rings[1..].iter().map(ring).collect())
    };
    // Disjoint boxes first: the boolean operation is the slow part.
    let bbox = |p: &Polygon| {
        p.rings[0]
            .iter()
            .fold((DVec2::MAX, DVec2::MIN), |(lo, hi), q| {
                (lo.min(*q), hi.max(*q))
            })
    };
    let ((alo, ahi), (blo, bhi)) = (bbox(a), bbox(b));
    if ahi.cmplt(blo).any() || bhi.cmplt(alo).any() {
        return 0.;
    }
    let (pa, pb) = (polygon(a), polygon(b));
    std::panic::catch_unwind(|| pa.intersection(&pb).unsigned_area()).unwrap_or(0.)
}

/// Common part of two collinear segments, if longer than eps.
fn collinear_overlap(p: DVec2, q: DVec2, r: DVec2, t: DVec2, eps: f64) -> Option<(DVec2, DVec2)> {
    let d = q - p;
    let length = d.length();
    if length <= eps {
        return None;
    }
    let u = d / length;
    let off = |x: DVec2| (x - p).perp_dot(u).abs();
    if off(r) > eps || off(t) > eps {
        return None;
    }
    let (a, b) = ((r - p).dot(u), (t - p).dot(u));
    let (start, end) = (a.min(b).max(0.), a.max(b).min(length));
    (end - start > eps).then(|| (p + u * start, p + u * end))
}

/// Consecutive anchor pairs of every axis: (axis, vertex, vertex).
fn bar_pieces(axes: &[Axis]) -> Vec<(usize, usize, usize)> {
    let mut out = vec![];
    for (i, axis) in axes.iter().enumerate() {
        let mut nodes: Vec<(f64, usize)> = axis.anchors.iter().map(|a| (a.t, a.vertex)).collect();
        nodes.push((0., axis.endpoints[0]));
        nodes.push((1., axis.endpoints[1]));
        nodes.sort_by(|x, y| x.0.total_cmp(&y.0));
        nodes.dedup_by_key(|n| n.1);
        for w in nodes.windows(2) {
            if w[0].1 != w[1].1 {
                out.push((i, w[0].1, w[1].1));
            }
        }
    }
    out
}

#[allow(clippy::too_many_arguments)]
fn audit_bars(
    surfaces: &[Surface],
    invalid: &BTreeSet<usize>,
    axes: &[Axis],
    contacts: &[Contact],
    edges: &[[usize; 2]],
    vertices: &[[f64; 3]],
    eps: f64,
    near: f64,
    findings: &mut Vec<Finding>,
) {
    struct Bar {
        p0: DVec3,
        p1: DVec3,
        anchors: BTreeSet<usize>,
        lo: DVec3,
        hi: DVec3,
    }
    let bars: Vec<Bar> = axes
        .iter()
        .map(|axis| {
            let (p0, p1) = (
                p3(vertices[axis.endpoints[0]]),
                p3(vertices[axis.endpoints[1]]),
            );
            Bar {
                p0,
                p1,
                anchors: axis
                    .anchors
                    .iter()
                    .map(|a| a.vertex)
                    .chain(axis.endpoints)
                    .collect(),
                lo: p0.min(p1),
                hi: p0.max(p1),
            }
        })
        .collect();
    let reach = DVec3::splat(near);
    let extent = bars
        .iter()
        .map(|b| (b.hi - b.lo).max_element())
        .fold(0., f64::max);
    let mut grid = Grid::new((near * 20.).max(extent.min(10.)).max(1e-6));
    for (i, b) in bars.iter().enumerate() {
        grid.insert(i, b.lo, b.hi);
    }
    let mut by_vertex = BTreeMap::<usize, BTreeSet<usize>>::new();
    for (i, b) in bars.iter().enumerate() {
        for &v in &b.anchors {
            by_vertex.entry(v).or_default().insert(i);
        }
    }
    let neighbours: Vec<BTreeSet<usize>> = bars
        .iter()
        .enumerate()
        .map(|(i, b)| {
            let mut n: BTreeSet<usize> = b
                .anchors
                .iter()
                .flat_map(|v| by_vertex[v].iter().copied())
                .collect();
            n.remove(&i);
            n
        })
        .collect();
    for (i, b) in bars.iter().enumerate() {
        let length = b.p0.distance(b.p1);
        if length < near {
            findings.push(
                Finding::new("short_bar", Class::Review)
                    .bars([i])
                    .at(b.p0)
                    .value(length),
            );
        }
        for j in grid.query(b.lo - reach, b.hi + reach) {
            if j <= i {
                continue;
            }
            let c = &bars[j];
            if (b.hi + reach).cmplt(c.lo).any() || (c.hi + reach).cmplt(b.lo).any() {
                continue;
            }
            let (dist, x) = closest_points(b.p0, b.p1, c.p0, c.p1);
            if dist > near {
                continue;
            }
            let (d1, d2) = (b.p1 - b.p0, c.p1 - c.p0);
            let parallel = d1.cross(d2).length() <= 1e-9 * d1.length() * d2.length();
            if dist <= eps {
                if parallel && d1.length() > 0. {
                    let u = d1.normalize();
                    let (ta, tb) = ((c.p0 - b.p0).dot(u), (c.p1 - b.p0).dot(u));
                    let overlap = d1.length().min(ta.max(tb)) - 0f64.max(ta.min(tb));
                    if overlap > eps {
                        findings.push(
                            Finding::new("overlapping_bars", Class::Failure)
                                .bars([i, j])
                                .at(x)
                                .value(overlap),
                        );
                        continue;
                    }
                }
                let shared = b
                    .anchors
                    .intersection(&c.anchors)
                    .any(|&v| p3(vertices[v]).distance(x) <= eps);
                if !shared {
                    findings.push(
                        Finding::new("unshared_bar_intersection", Class::Failure)
                            .bars([i, j])
                            .at(x)
                            .value(dist),
                    );
                }
            } else if !neighbours[i].contains(&j) && neighbours[i].is_disjoint(&neighbours[j]) {
                findings.push(
                    Finding::new("bar_near_miss", Class::Review)
                        .bars([i, j])
                        .at(x)
                        .value(dist),
                );
            }
        }
    }

    // Bars and surfaces.
    let mut point_contacts = BTreeMap::<(usize, usize), Vec<usize>>::new();
    let mut interval_contacts = BTreeMap::<(usize, usize), Vec<(f64, f64)>>::new();
    for c in contacts {
        match c {
            Contact::Point {
                axis,
                surface,
                vertex,
                ..
            } => point_contacts
                .entry((*axis, *surface))
                .or_default()
                .push(*vertex),
            Contact::Interval {
                axis,
                surface,
                start_t,
                end_t,
                ..
            } => interval_contacts
                .entry((*axis, *surface))
                .or_default()
                .push((*start_t, *end_t)),
        }
    }
    let surface_extent = surfaces
        .iter()
        .map(|s| (s.hi - s.lo).max_element())
        .fold(0., f64::max);
    let mut surface_grid = Grid::new((near * 20.).max(surface_extent.min(5.)).max(1e-6));
    for (j, s) in surfaces.iter().enumerate() {
        if !invalid.contains(&j) {
            surface_grid.insert(j, s.lo, s.hi);
        }
    }
    for (i, b) in bars.iter().enumerate() {
        let length = b.p0.distance(b.p1);
        for j in surface_grid.query(b.lo - reach, b.hi + reach) {
            let s = &surfaces[j];
            if (b.hi + reach).cmplt(s.lo).any() || (s.hi + reach).cmplt(b.lo).any() {
                continue;
            }
            let (h0, h1) = (s.height(b.p0), s.height(b.p1));
            let records = point_contacts.get(&(i, j)).cloned().unwrap_or_default();
            let check_point = |x: DVec3, radius: f64, findings: &mut Vec<Finding>| {
                let radius = radius.max(eps);
                let shared = records.iter().chain(b.anchors.iter()).any(|&v| {
                    p3(vertices[v]).distance(x) <= radius
                        && (s.vertex_ids.contains(&v) || records.contains(&v))
                });
                if !shared {
                    findings.push(
                        Finding::new("unshared_bar_surface_intersection", Class::Failure)
                            .surfaces([j])
                            .bars([i])
                            .at(x),
                    );
                }
            };
            if h0.abs() <= eps && h1.abs() <= eps && length > eps {
                let u = (b.p1 - b.p0) / length;
                let parts: Vec<(f64, f64)> = s
                    .polygon
                    .clip(s.project(b.p0), s.project(b.p1), eps)
                    .into_iter()
                    .map(|(t0, t1)| (t0 * length, t1 * length))
                    .collect();
                let lines: Vec<(f64, f64)> = parts
                    .iter()
                    .copied()
                    .filter(|(a, c)| c - a > 10. * eps)
                    .collect();
                let slack = 4. * eps;
                for &(a, c) in &lines {
                    let mut topo =
                        on_line(s.edge_ids.iter().map(|&e| edges[e]), vertices, b.p0, u, eps);
                    topo.extend(
                        interval_contacts
                            .get(&(i, j))
                            .into_iter()
                            .flatten()
                            .map(|&(t0, t1)| (t0 * length, t1 * length)),
                    );
                    if !covers(&mut topo, a + slack, c - slack, slack) {
                        findings.push(
                            Finding::new("bar_in_surface_without_contact", Class::Failure)
                                .surfaces([j])
                                .bars([i])
                                .at(b.p0 + u * a)
                                .at(b.p0 + u * c)
                                .value(c - a),
                        );
                    }
                }
                if !lines.is_empty() {
                    continue;
                }
                for &(a, c) in &parts {
                    check_point(b.p0 + u * ((a + c) / 2.), (c - a) / 2. + eps, findings);
                }
                if parts.is_empty() && records.is_empty() {
                    let d = s.distance(b.p0).min(s.distance(b.p1));
                    if d <= near {
                        findings.push(
                            Finding::new("bar_surface_near_miss", Class::Review)
                                .surfaces([j])
                                .bars([i])
                                .at(b.p0)
                                .value(d),
                        );
                    }
                }
                continue;
            }
            if h0 * h1 < 0. || h0.abs() <= eps || h1.abs() <= eps {
                let t = if (h0 - h1).abs() > 0. {
                    h0 / (h0 - h1)
                } else {
                    0.
                };
                let x = b.p0 + (b.p1 - b.p0) * t.clamp(0., 1.);
                if s.distance(x) <= eps {
                    check_point(x, eps, findings);
                    continue;
                }
            }
            if !records.is_empty() || interval_contacts.contains_key(&(i, j)) {
                continue;
            }
            let d = s.distance(b.p0).min(s.distance(b.p1));
            if d <= near {
                findings.push(
                    Finding::new("bar_surface_near_miss", Class::Review)
                        .surfaces([j])
                        .bars([i])
                        .at(b.p0)
                        .value(d),
                );
            }
        }
    }
}

fn plaxis_shapes(
    surfaces: &[Surface],
    used_edges: &BTreeSet<usize>,
    pieces: &[(usize, usize, usize)],
    edges: &[[usize; 2]],
    vertices: &[[f64; 3]],
    options: &Options,
    findings: &mut Vec<Finding>,
) {
    let tenth = options.element_size / 10.;
    for &e in used_edges {
        let [a, b] = edges[e];
        let length = p3(vertices[a]).distance(p3(vertices[b]));
        if length < tenth {
            let users = surfaces
                .iter()
                .enumerate()
                .filter(|(_, s)| s.edge_ids.contains(&e))
                .map(|(i, _)| i);
            findings.push(
                Finding::new("short_edge", Class::Plaxis)
                    .surfaces(users)
                    .edge(e)
                    .at(p3(vertices[a]))
                    .at(p3(vertices[b]))
                    .value(length),
            );
        }
    }
    for &(i, a, b) in pieces {
        let length = p3(vertices[a]).distance(p3(vertices[b]));
        if length > 0. && length < tenth {
            findings.push(
                Finding::new("short_bar_piece", Class::Plaxis)
                    .bars([i])
                    .at(p3(vertices[a]))
                    .at(p3(vertices[b]))
                    .value(length),
            );
        }
    }
    for (index, s) in surfaces.iter().enumerate() {
        for ring in &s.rings {
            let n = ring.len();
            for k in 0..n {
                let (p, q, r) = (
                    p3(vertices[ring[(k + n - 1) % n]]),
                    p3(vertices[ring[k]]),
                    p3(vertices[ring[(k + 1) % n]]),
                );
                let angle = (p - q).angle_between(r - q).to_degrees();
                if angle < options.sharp_angle {
                    findings.push(
                        Finding::new("sharp_corner", Class::Plaxis)
                            .surfaces([index])
                            .vertex(ring[k])
                            .at(q)
                            .value(angle),
                    );
                }
            }
        }
        // Narrow parts: a contour vertex near a contour edge of the same
        // surface away from its neighbourhood.
        let segments: Vec<[usize; 2]> = s
            .rings
            .iter()
            .flat_map(|r| (0..r.len()).map(move |k| [r[k], r[(k + 1) % r.len()]]))
            .collect();
        let mut neighbours = BTreeMap::<usize, BTreeSet<usize>>::new();
        for &[a, b] in &segments {
            neighbours.entry(a).or_default().extend([a, b]);
            neighbours.entry(b).or_default().extend([a, b]);
        }
        let mut grid = Grid::new(tenth.max(1e-9) * 4.);
        for (k, &[a, b]) in segments.iter().enumerate() {
            let (pa, pb) = (p3(vertices[a]), p3(vertices[b]));
            grid.insert(k, pa.min(pb), pa.max(pb));
        }
        for ring in &s.rings {
            for &v in ring {
                let near: BTreeSet<usize> = neighbours[&v]
                    .iter()
                    .flat_map(|w| neighbours[w].iter().copied())
                    .collect();
                let p = p3(vertices[v]);
                for k in grid.query(p - DVec3::splat(tenth), p + DVec3::splat(tenth)) {
                    let [a, b] = segments[k];
                    if near.contains(&a) || near.contains(&b) {
                        continue;
                    }
                    let d = p.distance(closest_on3(p, p3(vertices[a]), p3(vertices[b])));
                    if d < tenth {
                        findings.push(
                            Finding::new("narrow_face", Class::Plaxis)
                                .surfaces([index])
                                .vertex(v)
                                .at(p)
                                .value(d),
                        );
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reconstruction::{Model, PlaneFrame};

    /// A model of planar panels given as vertex rings; equal coordinates in
    /// the description are one vertex.
    fn model(panels: &[(&[[f64; 3]], [f64; 3])], scale: f64, shift: [f64; 3]) -> Model {
        let mut m = Model::new(1e-6 * scale, 0.001 * scale).unwrap();
        let mut ids = BTreeMap::<[i64; 3], usize>::new();
        let place = |p: [f64; 3]| [0, 1, 2].map(|k| p[k] * scale + shift[k]);
        for (ring, normal) in panels {
            let plane = m.add_plane(PlaneFrame::new(place(ring[0]), *normal).unwrap());
            let r: Vec<usize> = ring
                .iter()
                .map(|&p| {
                    let key = p.map(|x| (x * 1e6).round() as i64);
                    *ids.entry(key)
                        .or_insert_with(|| m.add_vertex(place(p)).unwrap())
                })
                .collect();
            m.add_surface(plane, vec![r], vec![]).unwrap();
        }
        m
    }

    const XY: [[f64; 3]; 4] = [[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.]];

    fn kinds(r: &Report, class: Class) -> Vec<String> {
        r.findings
            .iter()
            .filter(|f| f.class == class)
            .map(|f| f.kind.clone())
            .collect()
    }

    #[test]
    fn shared_boundary_passes_and_t_junction_fails() {
        for (scale, shift) in [(1., [0.; 3]), (10., [100., -50., 3.]), (0.1, [-2., 7., 1.])] {
            let wall = [[0., 0., 0.], [1., 0., 0.], [1., 0., 1.], [0., 0., 1.]];
            let r = run(
                &model(&[(&XY, [0., 0., 1.]), (&wall, [0., 1., 0.])], scale, shift),
                &[],
                &[],
                &Options::default(),
            );
            assert!(r.passed, "{:?}", r.findings);
            let t = [[0., 0.5, 0.], [1., 0.5, 0.], [1., 0.5, 1.], [0., 0.5, 1.]];
            let r = run(
                &model(&[(&XY, [0., 0., 1.]), (&t, [0., 1., 0.])], scale, shift),
                &[],
                &[],
                &Options::default(),
            );
            // The wall's bottom corners lie on the slab without being its
            // vertices: point contacts too.
            assert_eq!(
                kinds(&r, Class::Failure),
                vec![
                    "unrepresented_t_junction",
                    "unshared_point_contact",
                    "unshared_point_contact"
                ]
            );
            let f = &r
                .findings
                .iter()
                .find(|f| f.class == Class::Failure)
                .unwrap();
            assert!((f.value - scale).abs() < 1e-9 * scale);
            let cross = [[0., 0.5, -1.], [1., 0.5, -1.], [1., 0.5, 1.], [0., 0.5, 1.]];
            let r = run(
                &model(&[(&XY, [0., 0., 1.]), (&cross, [0., 1., 0.])], scale, shift),
                &[],
                &[],
                &Options::default(),
            );
            assert_eq!(kinds(&r, Class::Failure)[0], "unrepresented_crossing");
        }
    }

    #[test]
    fn coplanar_overlap_fails_and_touching_neighbours_need_shared_edges() {
        let shifted = [[0.5, 0., 0.], [1.5, 0., 0.], [1.5, 1., 0.], [0.5, 1., 0.]];
        let r = run(
            &model(
                &[(&XY, [0., 0., 1.]), (&shifted, [0., 0., 1.])],
                1.,
                [0.; 3],
            ),
            &[],
            &[],
            &Options::default(),
        );
        assert_eq!(kinds(&r, Class::Failure)[0], "coplanar_overlap");
        assert!((r.findings[0].value - 0.5).abs() < 1e-9);
        // Neighbours sharing the edge x = 1: fine; one starting at y = 0.5
        // with its own vertices on that edge: an unrepresented junction and
        // two point contacts.
        let next = [[1., 0., 0.], [2., 0., 0.], [2., 1., 0.], [1., 1., 0.]];
        assert!(
            run(
                &model(&[(&XY, [0., 0., 1.]), (&next, [0., 0., 1.])], 1., [0.; 3]),
                &[],
                &[],
                &Options::default()
            )
            .passed
        );
        let half = [[1., 0.5, 0.], [2., 0.5, 0.], [2., 1.5, 0.], [1., 1.5, 0.]];
        let r = run(
            &model(&[(&XY, [0., 0., 1.]), (&half, [0., 0., 1.])], 1., [0.; 3]),
            &[],
            &[],
            &Options::default(),
        );
        let failures = kinds(&r, Class::Failure);
        assert!(
            failures.contains(&"unrepresented_boundary_junction".into()),
            "{failures:?}"
        );
        assert!(
            failures.contains(&"unshared_point_contact".into()),
            "{failures:?}"
        );
    }

    #[test]
    fn plaxis_items_gap_short_edge_sharp_corner() {
        // A slab and a second one 20 mm beside it: a gap of each corner.
        let other = [[1.02, 0., 0.], [2., 0., 0.], [2., 1., 0.], [1.02, 1., 0.]];
        let r = run(
            &model(&[(&XY, [0., 0., 1.]), (&other, [0., 0., 1.])], 1., [0.; 3]),
            &[],
            &[],
            &Options::default(),
        );
        assert!(r.passed);
        assert_eq!(r.counts.get("gap"), Some(&4), "{:?}", r.findings);
        assert!(!r.plaxis_passed);
        // A 30 mm edge and a 5 degree corner.
        let thin = [
            [0., 0., 0.],
            [1., 0., 0.],
            [1., 0.03, 0.],
            [0.5, 0.0437, 0.],
        ];
        let r = run(
            &model(&[(&thin, [0., 0., 1.])], 1., [0.; 3]),
            &[],
            &[],
            &Options::default(),
        );
        assert!(r.counts.contains_key("short_edge"), "{:?}", r.counts);
        assert!(r.counts.contains_key("sharp_corner"), "{:?}", r.counts);
    }

    #[test]
    fn invalid_contour_is_a_failure() {
        // A bow tie cannot be built by the model; check the polygon test.
        let p = Polygon {
            rings: vec![vec![
                DVec2::new(0., 0.),
                DVec2::new(1., 1.),
                DVec2::new(1., 0.),
                DVec2::new(0., 1.),
            ]],
        };
        assert!(invalidity(&p, 1e-6).is_some());
        let square = Polygon {
            rings: vec![
                vec![
                    DVec2::new(0., 0.),
                    DVec2::new(4., 0.),
                    DVec2::new(4., 4.),
                    DVec2::new(0., 4.),
                ],
                vec![
                    DVec2::new(1., 1.),
                    DVec2::new(1., 2.),
                    DVec2::new(2., 2.),
                    DVec2::new(2., 1.),
                ],
            ],
        };
        assert!(invalidity(&square, 1e-6).is_none());
        assert!((square.area() - 15.).abs() < 1e-12);
        let mut outside = square;
        outside.rings[1] = vec![DVec2::new(5., 1.), DVec2::new(5., 2.), DVec2::new(6., 2.)];
        assert!(invalidity(&outside, 1e-6).is_some());
    }
}
