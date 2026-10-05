//! Exchange file for PLAXIS 3D: plates on planar hole-free polygons and
//! beams on straight lines, with elastic materials from the LIRA stiffness
//! block. `scripts/plaxis_export.py` builds the model in PLAXIS Input
//! through its Python API (plxscripting).
//!
//! PLAXIS surfaces carry no holes: a surface with openings is cut by lines
//! through the extreme points of each opening into hole-free pieces (one
//! material, common cut lines; PLAXIS intersects them into one mesh).
//! Junction lines and bar-surface crossings need not be exported: PLAXIS
//! intersects the geometry itself.
use crate::parsers::lira::Material;
use crate::reconstruction::assembly::edit::State;
use glam::DVec3;
use hashbrown::HashMap;
use serde::Serialize;
use std::collections::BTreeMap;

/// LIRA force unit (tonne-force) in kN.
pub const TONNE_TO_KN: f64 = 9.80665;

/// The loader script written beside an exchange file.
pub const LOADER: &str = include_str!("../scripts/plaxis_export.py");

/// A Python with plxscripting: the distribution installed with PLAXIS on
/// Windows (newest first), else `python` on the path.
pub fn find_python() -> Option<std::path::PathBuf> {
    let roots = [
        r"C:\ProgramData\Bentley\Geotechnical\PLAXIS Python Distribution V2\python\python.exe",
        r"C:\ProgramData\Bentley\Geotechnical\PLAXIS Python Distribution V1\python\python.exe",
        r"C:\ProgramData\Bentley\Geotechnical\PLAXIS Python Distribution V3\python\python.exe",
    ];
    roots
        .iter()
        .map(std::path::PathBuf::from)
        .find(|p| p.exists())
        .or_else(|| {
            Some(if cfg!(windows) {
                "python".into()
            } else {
                "python3".into()
            })
        })
}

/// Units of the exchange file: lengths stay in model units (m); forces are
/// converted by `force_factor` (LIRA t -> PLAXIS kN by default).
#[derive(Debug, Clone, Serialize)]
pub struct Units {
    pub length: &'static str,
    pub force_factor: f64,
    pub stress: &'static str,
    pub unit_weight: &'static str,
}

/// Which stiffness of a LIRA type goes to PLAXIS.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum StiffnessMode {
    /// The stiffness the LIRA model computes with: the numeric EF, EIy, EIz
    /// of a bar type, the WLKE/PLKE factors of a shell.
    Effective,
    /// From E and the section dimensions only.
    Nominal,
}

/// Export settings.
#[derive(Debug, Clone, Copy)]
pub struct Settings {
    /// Force unit of the model in kN (LIRA t: 9.80665).
    pub force_factor: f64,
    /// Cuts of holed surfaces avoid edges shorter than this.
    pub min_edge: f64,
    pub stiffness: StiffnessMode,
}

#[derive(Debug, Clone, Serialize)]
pub struct PlateMaterial {
    pub name: String,
    pub stiffness: u32,
    /// Young's modulus (kN/m2), Poisson's ratio, thickness (m), unit
    /// weight of the material (kN/m3, 0 when unknown).
    pub e: f64,
    pub nu: f64,
    pub d: f64,
    pub gamma: f64,
    /// How the values were derived from the LIRA type (what was changed or
    /// could not be transferred).
    pub notes: Vec<String>,
}

#[derive(Debug, Clone, Serialize)]
pub struct BeamMaterial {
    pub name: String,
    pub stiffness: u32,
    pub e: f64,
    pub nu: f64,
    /// LIRA section b (along local Y1) x h (along Z1), m.
    pub width: f64,
    pub height: f64,
    /// Area and second moments in PLAXIS terms: the beam's local axis 2 is
    /// the LIRA Z1 axis (the section height, PLAXIS "Height in local
    /// direction 2"), so I3 (bending about axis 3) is LIRA Iy and I2 is Iz.
    pub a: f64,
    pub i2: f64,
    pub i3: f64,
    pub gamma: f64,
    pub notes: Vec<String>,
}

#[derive(Debug, Clone, Serialize)]
pub struct Plate {
    pub surface: usize,
    pub stiffness: u32,
    pub material: Option<String>,
    /// Area of the surface (holes excluded), for checks of the pieces.
    pub area: f64,
    /// Planar simple polygons (3D points), no holes.
    pub polygons: Vec<Vec<[f64; 3]>>,
}

#[derive(Debug, Clone, Serialize)]
pub struct Beam {
    pub bar: usize,
    pub stiffness: u32,
    pub material: Option<String>,
    pub start: [f64; 3],
    pub end: [f64; 3],
    /// Local axis 2 of the beam (the LIRA Z1 axis by the LIRA default rule:
    /// the upward normal of the bar in its vertical plane, global X for a
    /// vertical bar). LIRA rotation angles of sections are not read.
    pub axis2: [f64; 3],
}

#[derive(Debug, Clone, Serialize)]
pub struct Exchange {
    pub format: &'static str,
    pub source: String,
    pub units: Units,
    pub stiffness: StiffnessMode,
    /// What the file does not carry or approximates.
    pub warnings: Vec<String>,
    pub plate_materials: Vec<PlateMaterial>,
    pub beam_materials: Vec<BeamMaterial>,
    pub plates: Vec<Plate>,
    pub beams: Vec<Beam>,
    /// Stiffness types used by the geometry without an exportable material
    /// (PLAXIS objects are created without a material).
    pub missing_materials: Vec<u32>,
    /// Holed surfaces cut into pieces.
    pub cut_surfaces: usize,
    /// Holed surfaces whose cutting failed the area check, exported as
    /// triangles instead.
    pub triangulated_surfaces: usize,
    /// Load case numbers and names.
    pub load_cases: Vec<(u32, String)>,
    pub loads: Vec<crate::loads::Load>,
    pub load_report: Option<crate::loads::Report>,
}

fn signed_area(points: &[[f64; 2]]) -> f64 {
    (0..points.len())
        .map(|i| {
            let (p, q) = (points[i], points[(i + 1) % points.len()]);
            p[0] * q[1] - q[0] * p[1]
        })
        .sum::<f64>()
        / 2.
}

/// A polygon as PLAXIS takes it: PLAXIS fits the plane through the first
/// points ("Define plane: First points"), so three collinear points first
/// (a straight contour run with T-junction vertices) make a planar polygon
/// "not coplanar" and invalid. Vertices on the straight line between their
/// neighbours are dropped (PLAXIS intersects the geometry itself and
/// recreates junction points), and the ring starts at the corner whose
/// triangle with the next two points is the largest.
pub fn plaxis_polygon(ring: Vec<[f64; 3]>) -> Vec<[f64; 3]> {
    let p = |v: [f64; 3]| DVec3::from_array(v);
    let mut ring = ring;
    loop {
        let n = ring.len();
        if n <= 3 {
            break;
        }
        // Distance of each vertex from the chord of its neighbours, relative
        // to the chord: below 1e-9 it is a straight-line vertex.
        let straight = (0..n).find(|&i| {
            let (a, b, c) = (p(ring[(i + n - 1) % n]), p(ring[i]), p(ring[(i + 1) % n]));
            let chord = c - a;
            let len = chord.length();
            len > 0.
                && (b - a).cross(chord).length() / len <= 1e-9 * len.max(1.)
                && (b - a).dot(chord) > 0.
                && (c - b).dot(chord) > 0.
        });
        match straight {
            Some(i) => {
                ring.remove(i);
            }
            None => break,
        }
    }
    let n = ring.len();
    let area = |i: usize| {
        let (a, b, c) = (p(ring[i]), p(ring[(i + 1) % n]), p(ring[(i + 2) % n]));
        (b - a).cross(c - a).length()
    };
    if let Some(start) = (0..n).max_by(|&x, &y| area(x).total_cmp(&area(y))) {
        ring.rotate_left(start);
    }
    ring
}

/// Constrained Delaunay triangles of a polygon with holes (fallback).
fn triangles(contours: &[Vec<[f64; 2]>]) -> Vec<Vec<[f64; 2]>> {
    use spade::{ConstrainedDelaunayTriangulation, Point2, Triangulation};
    let mut cdt = ConstrainedDelaunayTriangulation::<Point2<f64>>::new();
    for ring in contours {
        let handles: Vec<_> = ring
            .iter()
            .filter_map(|p| cdt.insert(Point2::new(p[0], p[1])).ok())
            .collect();
        for i in 0..handles.len() {
            let (a, b) = (handles[i], handles[(i + 1) % handles.len()]);
            if a != b && cdt.can_add_constraint(a, b) {
                cdt.add_constraint(a, b);
            }
        }
    }
    cdt.inner_faces()
        .map(|f| {
            f.vertices()
                .map(|v| [v.position().x, v.position().y])
                .to_vec()
        })
        .filter(|t| {
            let c = [
                (t[0][0] + t[1][0] + t[2][0]) / 3.,
                (t[0][1] + t[1][1] + t[2][1]) / 3.,
            ];
            inside(&contours[0], c) && !contours[1..].iter().any(|h| inside(h, c))
        })
        .map(|mut t| {
            if signed_area(&t) < 0. {
                t.reverse();
            }
            t
        })
        .collect()
}

/// Two segments meet (touching included).
fn segments_cross(a: [f64; 2], b: [f64; 2], c: [f64; 2], d: [f64; 2]) -> bool {
    let orient = |p: [f64; 2], q: [f64; 2], r: [f64; 2]| {
        (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0])
    };
    let (d1, d2) = (orient(a, b, c), orient(a, b, d));
    let (d3, d4) = (orient(c, d, a), orient(c, d, b));
    let within = |p: [f64; 2], q: [f64; 2], r: [f64; 2]| {
        r[0] >= p[0].min(q[0])
            && r[0] <= p[0].max(q[0])
            && r[1] >= p[1].min(q[1])
            && r[1] <= p[1].max(q[1])
    };
    if d1 * d2 < 0. && d3 * d4 < 0. {
        return true;
    }
    (d1 == 0. && within(a, b, c))
        || (d2 == 0. && within(a, b, d))
        || (d3 == 0. && within(c, d, a))
        || (d4 == 0. && within(c, d, b))
}

fn inside(ring: &[[f64; 2]], p: [f64; 2]) -> bool {
    let mut odd = false;
    for i in 0..ring.len() {
        let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
        if (a[1] > p[1]) != (b[1] > p[1])
            && p[0] < a[0] + (p[1] - a[1]) / (b[1] - a[1]) * (b[0] - a[0])
        {
            odd = !odd;
        }
    }
    odd
}

/// A face boundary through one point twice (holes touching at a vertex)
/// as simple loops; unchanged when a loop is a hole (negative).
fn unpinch(face: Vec<[f64; 2]>) -> Vec<Vec<[f64; 2]>> {
    for i in 0..face.len() {
        if let Some(j) = (i + 1..face.len()).find(|&j| face[j] == face[i]) {
            let inner: Vec<[f64; 2]> = face[i..j].to_vec();
            let outer: Vec<[f64; 2]> = face[..i].iter().chain(&face[j..]).copied().collect();
            if signed_area(&inner) <= 0. || signed_area(&outer) <= 0. {
                return vec![face];
            }
            let mut out = unpinch(inner);
            out.extend(unpinch(outer));
            return out;
        }
    }
    vec![face]
}

/// Hole-free pieces of a polygon with holes (plane coordinates). From the
/// lowest and highest point of every hole at its smallest and its largest
/// `u` a cut runs along `v` (down and up) through the material to the
/// first contour it meets: each chain of cuts ends on the outer contour,
/// so the faces are `u`-monotone simple polygons. Cuts whose two pieces
/// share nothing but the cut are then removed: about one piece more per
/// hole, each remaining cut shared by its two pieces.
pub fn hole_free(contours: &[Vec<[f64; 2]>], min_edge: f64) -> Vec<Vec<[f64; 2]>> {
    if contours.len() <= 1 {
        return contours.to_vec();
    }
    // Cuts along v, or along u (axes swapped): the variant with fewer short
    // edges and sharp corners, then fewer pieces.
    let swap = |r: &Vec<[f64; 2]>| r.iter().map(|p| [p[1], p[0]]).rev().collect::<Vec<_>>();
    let along_v = cut_holes(contours, min_edge);
    let swapped: Vec<Vec<[f64; 2]>> = contours.iter().map(swap).collect();
    let along_u: Vec<Vec<[f64; 2]>> = cut_holes(&swapped, min_edge).iter().map(swap).collect();
    let score = |pieces: &[Vec<[f64; 2]>]| {
        let defects = pieces
            .iter()
            .flat_map(|f| {
                let n = f.len();
                (0..n).map(move |k| (f[(k + n - 1) % n], f[k], f[(k + 1) % n]))
            })
            .filter(|&(p, q, r)| {
                let (u, v) = ([p[0] - q[0], p[1] - q[1]], [r[0] - q[0], r[1] - q[1]]);
                let (lu, lv) = (u[0].hypot(u[1]), v[0].hypot(v[1]));
                let cos = (u[0] * v[0] + u[1] * v[1]) / (lu * lv);
                lv < min_edge || cos > 10f64.to_radians().cos()
            })
            .count();
        (defects, pieces.len())
    };
    if score(&along_u) < score(&along_v) {
        along_u
    } else {
        along_v
    }
}

fn cut_holes(contours: &[Vec<[f64; 2]>], min_edge: f64) -> Vec<Vec<[f64; 2]>> {
    let (mut lo, mut hi) = ([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]);
    for p in &contours[0] {
        for k in 0..2 {
            lo[k] = lo[k].min(p[k]);
            hi[k] = hi[k].max(p[k]);
        }
    }
    let eps = 1e-9 * (hi[0] - lo[0]).max(hi[1] - lo[1]);
    let snap = 1e3 * eps;
    let edges: Vec<([f64; 2], [f64; 2])> = contours
        .iter()
        .flat_map(|r| (0..r.len()).map(move |i| (r[i], r[(i + 1) % r.len()])))
        .collect();
    let material =
        |p: [f64; 2]| inside(&contours[0], p) && !contours[1..].iter().any(|h| inside(h, p));
    // A cut within `band` of a vertex goes to it (lines a fraction of a
    // millimetre apart are one line), when the slanted cut stays clear.
    let band = (min_edge / 10.).max(snap);
    let clear = |start: [f64; 2], end: [f64; 2]| {
        edges.iter().all(|&(p, q)| {
            p == end || q == end || p == start || q == start || !segments_cross(start, end, p, q)
        }) && material([(start[0] + end[0]) / 2., (start[1] + end[1]) / 2.])
    };
    // The first contour point a vertical ray from `start` meets: `None`
    // when the ray leaves along a contour edge or meets nothing.
    let ray = |start: [f64; 2], dir: f64| -> Option<(usize, [f64; 2])> {
        if !material([start[0], start[1] + dir * snap]) {
            return None;
        }
        let mut best: Option<(f64, usize, [f64; 2])> = None;
        for (k, &(a, b)) in edges.iter().enumerate() {
            // A vertex next to the ray (within `snap` of its line) is met:
            // a cut never runs a hair's breadth past a contour.
            for end in [a, b] {
                let d = (end[1] - start[1]) * dir;
                let off = (end[0] - start[0]).abs();
                if off <= band
                    && d > snap
                    && best.is_none_or(|x| d < x.0)
                    && (off <= snap || clear(start, end))
                {
                    best = Some((d, k, end));
                }
            }
            if (a[0] - start[0]).abs() <= band && (b[0] - start[0]).abs() <= band {
                // Along the ray: met at its near end (above), or the ray
                // runs on the contour.
                let (da, db) = ((a[1] - start[1]) * dir, (b[1] - start[1]) * dir);
                if da.max(db) > snap && da.min(db) <= snap {
                    return None;
                }
                continue;
            }
            if (a[0] - start[0]) * (b[0] - start[0]) > eps * eps {
                continue;
            }
            let t = if (b[0] - a[0]).abs() <= eps {
                0.
            } else {
                (start[0] - a[0]) / (b[0] - a[0])
            };
            let y = a[1] + t.clamp(0., 1.) * (b[1] - a[1]);
            let d = (y - start[1]) * dir;
            if d > snap && best.is_none_or(|x| d < x.0) {
                best = Some((d, k, [start[0], y]));
            }
        }
        let (_, k, mut hit) = best?;
        for end in [edges[k].0, edges[k].1] {
            if (end[0] - hit[0]).abs() <= snap && (end[1] - hit[1]).abs() <= snap {
                hit = end;
            }
        }
        // A hit closer than `min_edge` to a vertex of its edge goes to that
        // vertex (no short edge), when the slanted cut stays clear of every
        // contour.
        let (a, b) = edges[k];
        if hit != a && hit != b {
            let near = [a, b]
                .into_iter()
                .filter(|e| (e[0] - hit[0]).hypot(e[1] - hit[1]) < min_edge)
                .min_by(|e, f| {
                    (e[0] - hit[0])
                        .hypot(e[1] - hit[1])
                        .total_cmp(&(f[0] - hit[0]).hypot(f[1] - hit[1]))
                });
            if let Some(end) = near {
                if clear(start, end) {
                    hit = end;
                }
            }
        }
        Some((k, hit))
    };
    let mut cuts: Vec<([f64; 2], [f64; 2])> = vec![];
    let mut splits: Vec<Vec<[f64; 2]>> = vec![vec![]; edges.len()];
    for hole in &contours[1..] {
        let umin = hole.iter().map(|p| p[0]).fold(f64::INFINITY, f64::min);
        let umax = hole.iter().map(|p| p[0]).fold(f64::NEG_INFINITY, f64::max);
        for u in [umin, umax] {
            let column: Vec<[f64; 2]> = hole
                .iter()
                .copied()
                .filter(|p| (p[0] - u).abs() <= snap)
                .collect();
            let low =
                column
                    .iter()
                    .copied()
                    .fold([u, f64::INFINITY], |a, p| if p[1] < a[1] { p } else { a });
            let high = column.iter().copied().fold([u, f64::NEG_INFINITY], |a, p| {
                if p[1] > a[1] {
                    p
                } else {
                    a
                }
            });
            for (start, dir) in [(low, -1.), (high, 1.)] {
                let Some((k, hit)) = ray(start, dir) else {
                    continue;
                };
                let (a, b) = edges[k];
                if hit != a && hit != b {
                    splits[k].push(hit);
                }
                let cut = if start <= hit {
                    (start, hit)
                } else {
                    (hit, start)
                };
                if !cuts.contains(&cut) {
                    cuts.push(cut);
                }
            }
        }
    }
    // Plane graph: points, undirected segments.
    let mut points: Vec<[f64; 2]> = vec![];
    let id = |p: [f64; 2], points: &mut Vec<[f64; 2]>| {
        points.iter().position(|q| *q == p).unwrap_or_else(|| {
            points.push(p);
            points.len() - 1
        })
    };
    let mut segments: Vec<(usize, usize)> = vec![];
    for (k, &(a, b)) in edges.iter().enumerate() {
        let mut chain = splits[k].clone();
        let d = [b[0] - a[0], b[1] - a[1]];
        let along = |p: &[f64; 2]| (p[0] - a[0]) * d[0] + (p[1] - a[1]) * d[1];
        chain.sort_by(|p, q| along(p).total_cmp(&along(q)));
        chain.dedup();
        let mut previous = id(a, &mut points);
        for p in chain.into_iter().chain([b]) {
            let next = id(p, &mut points);
            segments.push((previous, next));
            previous = next;
        }
    }
    let first_cut = segments.len();
    for &(a, b) in &cuts {
        let (a, b) = (id(a, &mut points), id(b, &mut points));
        segments.push((a, b));
    }
    // Faces: at each vertex the outgoing half-edges in angular order; a
    // face continues with the edge clockwise next to the reversed one.
    let mut outgoing: Vec<Vec<usize>> = vec![vec![]; points.len()];
    let half: Vec<(usize, usize)> = segments
        .iter()
        .flat_map(|&(a, b)| [(a, b), (b, a)])
        .collect();
    for (h, &(a, _)) in half.iter().enumerate() {
        outgoing[a].push(h);
    }
    let angle = |h: usize| {
        let (a, b) = half[h];
        // + 0.0 turns a -0.0 into 0.0: atan2(-0.0, -x) is -pi, not pi.
        (points[b][1] - points[a][1] + 0.0).atan2(points[b][0] - points[a][0] + 0.0)
    };
    for list in &mut outgoing {
        list.sort_by(|&x, &y| angle(x).total_cmp(&angle(y)));
    }
    let coords = |f: &[usize]| f.iter().map(|&v| points[v]).collect::<Vec<_>>();
    let mut used = vec![false; half.len()];
    let mut faces: Vec<Vec<usize>> = vec![];
    for start in 0..half.len() {
        if used[start] {
            continue;
        }
        let mut face = vec![];
        let mut h = start;
        while !used[h] {
            used[h] = true;
            let (a, b) = half[h];
            face.push(a);
            let list = &outgoing[b];
            let i = list.iter().position(|&x| half[x] == (b, a)).unwrap();
            h = list[(i + list.len() - 1) % list.len()];
        }
        let ring = coords(&face);
        if signed_area(&ring) <= eps * eps {
            continue; // the unbounded face
        }
        // Inside the face: a point just left of its first edge.
        let (p, q) = (ring[0], ring[1]);
        let (dx, dy) = (q[0] - p[0], q[1] - p[1]);
        let len = dx.hypot(dy);
        let sample = [
            (p[0] + q[0]) / 2. - dy / len * snap,
            (p[1] + q[1]) / 2. + dx / len * snap,
        ];
        if material(sample) {
            faces.push(face);
        }
    }
    // Remove cuts between two pieces sharing only the cut's ends: the
    // union of two simple polygons meeting along one edge is simple.
    let position = |f: &[usize], a: usize, b: usize| {
        (0..f.len()).find(|&i| f[i] == a && f[(i + 1) % f.len()] == b)
    };
    for &(a, b) in &segments[first_cut..] {
        let left = (0..faces.len()).find(|&i| position(&faces[i], a, b).is_some());
        let right = (0..faces.len()).find(|&i| position(&faces[i], b, a).is_some());
        let (Some(l), Some(r)) = (left, right) else {
            continue;
        };
        if l == r {
            continue;
        }
        let shared = faces[l].iter().filter(|v| faces[r].contains(v)).count();
        if shared != 2
            || faces[l].len()
                != faces[l]
                    .iter()
                    .collect::<std::collections::BTreeSet<_>>()
                    .len()
        {
            continue;
        }
        let i = position(&faces[l], a, b).unwrap();
        let j = position(&faces[r], b, a).unwrap();
        // Left from b around to a, then right from a around to b (ends
        // excluded).
        let (fl, fr) = (&faces[l], &faces[r]);
        let mut merged: Vec<usize> = (0..fl.len()).map(|k| fl[(i + 1 + k) % fl.len()]).collect();
        merged.extend((1..fr.len() - 1).map(|k| fr[(j + 1 + k) % fr.len()]));
        let (keep, drop) = (l.min(r), l.max(r));
        faces[keep] = merged;
        faces.remove(drop);
    }
    // A cut end left on a straight contour (180 degrees, not a contour
    // vertex) goes.
    let original: std::collections::BTreeSet<usize> = contours
        .iter()
        .flatten()
        .filter_map(|p| points.iter().position(|q| q == p))
        .collect();
    let mut out = vec![];
    for face in faces {
        let n = face.len();
        let kept: Vec<usize> = (0..n)
            .filter(|&k| {
                let v = face[k];
                if original.contains(&v) {
                    return true;
                }
                let (p, q, r) = (
                    points[face[(k + n - 1) % n]],
                    points[v],
                    points[face[(k + 1) % n]],
                );
                let cross = (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0]);
                cross.abs() > eps * ((r[0] - p[0]).hypot(r[1] - p[1])).max(eps)
            })
            .map(|k| face[k])
            .collect();
        out.extend(unpinch(coords(&kept)));
    }
    out
}

/// How a surface became polygons.
pub enum Cut {
    /// One polygon, the contour.
    No,
    /// Cut into hole-free pieces.
    Pieces,
    /// Triangles (the cutting failed its area check).
    Triangles,
}

/// The hole-free planar polygons of surface `s`, as PLAXIS takes them.
pub fn surface_polygons(
    model: &crate::reconstruction::Model,
    s: usize,
    min_edge: f64,
) -> (Vec<Vec<[f64; 3]>>, Cut) {
    let surface = &model.surfaces()[s];
    let plane = &model.planes()[surface.plane];
    if surface.boundaries.len() == 1 {
        // Exact vertex positions of the contour.
        let ring = surface.boundaries[0]
            .iter()
            .map(|u| model.vertices()[model.edges()[u.edge][usize::from(u.reversed)]])
            .collect();
        return (vec![plaxis_polygon(ring)], Cut::No);
    }
    let mut pieces = hole_free(&surface.contours, min_edge);
    let area = signed_area(&surface.contours[0]).abs()
        - surface.contours[1..]
            .iter()
            .map(|h| signed_area(h).abs())
            .sum::<f64>();
    let covered: f64 = pieces.iter().map(|p| signed_area(p)).sum();
    let mut how = Cut::Pieces;
    if (covered - area).abs() > 1e-6 * area {
        // Never export wrong material: triangles always fit.
        how = Cut::Triangles;
        pieces = triangles(&surface.contours);
    }
    let polygons = pieces
        .into_iter()
        .map(|piece| plaxis_polygon(piece.into_iter().map(|uv| plane.lift(uv)).collect()))
        .collect();
    (polygons, how)
}

/// The exchange file of the current geometry; cuts of holed surfaces make
/// no edge shorter than `min_edge` where they can avoid it.
pub fn exchange(
    state: &State,
    materials: &HashMap<u32, Material>,
    settings: Settings,
    source: &str,
) -> Exchange {
    let Settings {
        force_factor,
        min_edge,
        stiffness: mode,
    } = settings;
    let model = &state.model;
    let mut warnings = vec![];
    let mut missing = std::collections::BTreeSet::new();
    let mut plate_materials = BTreeMap::new();
    let mut beam_materials = BTreeMap::new();
    let mut plate_material = |id: u32| -> Option<String> {
        let Some(&Material::Plate {
            e,
            nu,
            thickness,
            density,
            membrane,
            bending,
            shear,
        }) = materials.get(&id)
        else {
            missing.insert(id);
            return None;
        };
        let name = format!("GEI_{id}_h{:.0}", thickness * 1000.);
        plate_materials.entry(id).or_insert_with(|| {
            let mut notes = vec![];
            let (km, kb) = (membrane.unwrap_or(1.), bending.unwrap_or(1.));
            let (mut e, mut d, mut gamma) = (e, thickness, density.unwrap_or(0.));
            if (km - 1.).abs() > 1e-12 || (kb - 1.).abs() > 1e-12 {
                if mode == StiffnessMode::Effective && km > 0. && kb > 0. {
                    // E'd' = E d km and E'd'^3 = E d^3 kb keep the membrane
                    // and bending stiffness; gamma' d' = gamma d the weight.
                    let d2 = thickness * (kb / km).sqrt();
                    e = e * km * thickness / d2;
                    gamma = gamma * thickness / d2;
                    d = d2;
                    notes.push(format!(
                        "WLKE {km}, PLKE {kb}: equivalent thickness {d2:.4} m and E keep the membrane and bending stiffness and the weight"
                    ));
                } else {
                    notes.push(format!("WLKE {km}, PLKE {kb} not applied (nominal stiffness)"));
                }
            }
            if shear.iter().flatten().any(|g| (g - 1.).abs() > 1e-12) {
                let show = |g: Option<f64>| g.map_or("-".to_string(), |g| g.to_string());
                notes.push(format!(
                    "WLKG {}, PLKG {} not transferable: PLAXIS derives G from E and nu",
                    show(shear[0]),
                    show(shear[1])
                ));
            }
            PlateMaterial {
                name: name.clone(),
                stiffness: id,
                e: e * force_factor,
                nu,
                d,
                gamma: gamma * force_factor,
                notes,
            }
        });
        Some(name)
    };
    let mut cut_surfaces = 0;
    let mut triangulated_surfaces = 0;
    let plates: Vec<Plate> = (0..model.surfaces().len())
        .map(|s| {
            let surface = &model.surfaces()[s];
            let (polygons, how) = surface_polygons(model, s, min_edge);
            match how {
                Cut::No => {}
                Cut::Pieces => cut_surfaces += 1,
                Cut::Triangles => {
                    cut_surfaces += 1;
                    triangulated_surfaces += 1;
                }
            }
            let stiffness = state.stiffness[s];
            let area = signed_area(&surface.contours[0]).abs()
                - surface.contours[1..]
                    .iter()
                    .map(|h| signed_area(h).abs())
                    .sum::<f64>();
            Plate {
                surface: s,
                area,
                stiffness,
                material: plate_material(stiffness),
                polygons,
            }
        })
        .collect();
    let mut beams = vec![];
    for (k, axis) in state.axes.iter().enumerate() {
        let [p, q] = axis
            .endpoints
            .map(|v| DVec3::from_array(model.vertices()[v]));
        let mut spans = axis.spans.clone();
        spans.sort_by(|a, b| a.start_t.total_cmp(&b.start_t));
        // Consecutive spans of one stiffness make one beam.
        let mut runs: Vec<(u32, f64, f64)> = vec![];
        for span in &spans {
            match runs.last_mut() {
                Some(r) if r.0 == span.stiffness && (r.2 - span.start_t).abs() < 1e-9 => {
                    r.2 = span.end_t
                }
                _ => runs.push((span.stiffness, span.start_t, span.end_t)),
            }
        }
        let along = (q - p).normalize();
        // LIRA default local axes: Z1 in the vertical plane of the bar,
        // upwards; global X for a vertical bar.
        let horizontal = DVec3::new(along.x, along.y, 0.).length();
        let axis2 = if horizontal < 1e-9 {
            DVec3::X
        } else {
            (DVec3::Z - along * along.z).normalize()
        };
        for (stiffness, t0, t1) in runs {
            if (t1 - t0) * p.distance(q) <= model.precision() {
                warnings.push(format!("bar {k}: zero-length piece [{t0}, {t1}] skipped"));
                continue;
            }
            let material = match materials.get(&stiffness) {
                Some(&Material::Bar {
                    e,
                    nu,
                    width,
                    height,
                    density,
                    stiffness: numeric,
                }) => {
                    let name = format!("S0_{stiffness}_{:.0}x{:.0}", width * 100., height * 100.);
                    beam_materials.entry(stiffness).or_insert_with(|| {
                        let mut notes = vec![];
                        let nominal = (width * height, width * height.powi(3) / 12., height * width.powi(3) / 12.);
                        let (a, iy, iz) = match (mode, numeric) {
                            (StiffnessMode::Effective, Some([ef, eiy, eiz, _])) => {
                                notes.push(format!(
                                    "EA, EIy, EIz from the LIRA type ({:.3}, {:.3}, {:.3} of nominal); GIk not transferable (PLAXIS derives J)",
                                    ef / e / nominal.0,
                                    eiy / e / nominal.1,
                                    eiz / e / nominal.2
                                ));
                                (ef / e, eiy / e, eiz / e)
                            }
                            (StiffnessMode::Nominal, Some(_)) => {
                                notes.push("nominal b x h stiffness; the LIRA type gives other EF/EI".into());
                                nominal
                            }
                            _ => nominal,
                        };
                        if (width - height).abs() > 1e-9 {
                            notes.push("rectangular: orientation from the LIRA default local axes".into());
                        }
                        BeamMaterial {
                            name: name.clone(),
                            stiffness,
                            e: e * force_factor,
                            nu: nu.unwrap_or(0.2),
                            width,
                            height,
                            a,
                            i2: iz,
                            i3: iy,
                            // RO of a bar is per length; PLAXIS weighs gamma x A.
                            gamma: density.unwrap_or(0.) / a * force_factor,
                            notes,
                        }
                    });
                    Some(name)
                }
                _ => {
                    missing.insert(stiffness);
                    None
                }
            };
            beams.push(Beam {
                bar: k,
                stiffness,
                material,
                start: p.lerp(q, t0).to_array(),
                end: p.lerp(q, t1).to_array(),
                axis2: axis2.to_array(),
            });
        }
    }
    warnings.push("LIRA rotation angles of bar sections are not read: the LIRA default local axes are assumed".into());
    warnings.push(
        "loads, supports, releases, eccentricities, soils and stages are not exported".into(),
    );
    Exchange {
        format: "topo-plaxis-1",
        source: source.into(),
        stiffness: mode,
        warnings,
        units: Units {
            length: "m",
            force_factor,
            stress: "kN/m2",
            unit_weight: "kN/m3",
        },
        plate_materials: plate_materials.into_values().collect(),
        beam_materials: beam_materials.into_values().collect(),
        plates,
        beams,
        missing_materials: missing.into_iter().collect(),
        cut_surfaces,
        triangulated_surfaces,
        load_cases: vec![],
        loads: vec![],
        load_report: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn area(r: &[[f64; 2]]) -> f64 {
        (0..r.len())
            .map(|i| {
                let (p, q) = (r[i], r[(i + 1) % r.len()]);
                p[0] * q[1] - q[0] * p[1]
            })
            .sum::<f64>()
            .abs()
            / 2.
    }

    #[test]
    fn holes_are_cut_into_hole_free_pieces_of_the_same_area() {
        // A 10 x 6 slab with two openings, one non-convex.
        let outer = vec![[0., 0.], [10., 0.], [10., 6.], [0., 6.]];
        let a = vec![[1., 1.], [1., 2.], [3., 2.], [3., 1.]];
        let b = vec![
            [5., 2.],
            [5., 5.],
            [8., 5.],
            [8., 4.],
            [6., 4.],
            [6., 3.],
            [8., 3.],
            [8., 2.],
        ];
        for (shift, scale) in [([0., 0.], 1.), ([512.3, -77.1], 0.01), ([-3., 9.], 40.)] {
            let place = |r: &Vec<[f64; 2]>| -> Vec<[f64; 2]> {
                r.iter()
                    .map(|p| [p[0] * scale + shift[0], p[1] * scale + shift[1]])
                    .collect()
            };
            let contours = vec![place(&outer), place(&a), place(&b)];
            let pieces = hole_free(&contours, 0.);
            let total: f64 = pieces.iter().map(|p| area(p)).sum();
            let expected = (60. - 2. - 7.) * scale * scale;
            assert!(
                (total - expected).abs() < 1e-9 * expected,
                "{total} {expected}"
            );
            assert!((2..=3).contains(&pieces.len()), "{}", pieces.len());
            for piece in &pieces {
                let mut unique = piece.clone();
                unique.sort_by(|p, q| p.partial_cmp(q).unwrap());
                unique.dedup();
                assert_eq!(unique.len(), piece.len(), "a slit in {piece:?}");
            }
        }
        // A column of windows (holes one above another, edges on one line)
        // and two holes touching at a corner.
        let column = vec![
            vec![[0., 0.], [6., 0.], [6., 10.], [0., 10.]],
            vec![[2., 1.], [2., 3.], [3., 3.], [3., 1.]],
            vec![[2., 4.], [2., 6.], [3., 6.], [3., 4.]],
            vec![[2., 7.], [2., 9.], [3., 9.], [3., 7.]],
            vec![[4., 1.], [4., 2.], [5., 2.], [5., 1.]],
            vec![[5., 2.], [5., 3.], [5.5, 3.], [5.5, 2.]],
        ];
        let pieces = hole_free(&column, 0.);
        let total: f64 = pieces.iter().map(|p| area(p)).sum();
        assert!((total - (60. - 6. - 1. - 0.5)).abs() < 1e-9, "{total}");
        for piece in &pieces {
            let mut unique = piece.clone();
            unique.sort_by(|p, q| p.partial_cmp(q).unwrap());
            unique.dedup();
            assert_eq!(unique.len(), piece.len(), "a slit or pinch in {piece:?}");
            assert!(signed_area(piece) > 0.);
        }
        // Holes one above another, overlapping in u: their own cuts join
        // them to each other only.
        let stacked = vec![
            vec![[0., 0.], [10., 0.], [10., 10.], [0., 10.]],
            vec![[2., 2.], [6., 2.], [6., 4.], [2., 4.]],
            vec![[4., 5.], [8., 5.], [8., 7.], [4., 7.]],
        ];
        let pieces = hole_free(&stacked, 0.);
        let total: f64 = pieces.iter().map(|p| area(p)).sum();
        assert!((total - 84.).abs() < 1e-9, "{total} {pieces:?}");
        // The fallback triangles cover the material exactly.
        let t: f64 = triangles(&stacked).iter().map(|p| signed_area(p)).sum();
        assert!((t - 84.).abs() < 1e-9, "{t}");
        // Without holes the contour is kept as it is.
        assert_eq!(hole_free(&[outer.clone()], 0.), vec![outer]);
    }

    #[test]
    fn stiffness_modes_axes_and_shell_factors() {
        use crate::reconstruction::assembly::bars::{Anchor, Axis};
        use crate::reconstruction::assembly::junctions::tests::{build, slab, Placement};
        use crate::reconstruction::recognize::SourceSpan;
        let place = &Placement::all()[0];
        let mut model = build(place, &[slab(0., 4.)]);
        let axis = |m: &mut crate::reconstruction::Model, a: [f64; 3], b: [f64; 3]| {
            let (va, vb) = (m.add_vertex(a).unwrap(), m.add_vertex(b).unwrap());
            Axis {
                source_axis: 0,
                endpoints: [va, vb],
                anchors: vec![
                    Anchor {
                        source_node: 0,
                        vertex: va,
                        t: 0.,
                    },
                    Anchor {
                        source_node: 1,
                        vertex: vb,
                        t: 1.,
                    },
                ],
                spans: vec![SourceSpan {
                    element: 1,
                    stiffness: 5,
                    start_t: 0.,
                    end_t: 1.,
                }],
            }
        };
        let beam = axis(&mut model, [0., 1., 3.], [4., 1., 3.]);
        let column = axis(&mut model, [1., 1., 0.], [1., 1., 3.]);
        let state = State {
            model,
            axes: vec![beam, column],
            contacts: vec![],
            stiffness: vec![7],
            patches: vec![0],
            removed: vec![],
            removed_bars: vec![],
            joints: Default::default(),
        };
        let mut materials = HashMap::new();
        // 0.5 x 0.8 m with EIy 0.6 and EIz 0.6 of nominal (skala type 1).
        materials.insert(
            5,
            Material::Bar {
                e: 3e6,
                nu: Some(0.2),
                width: 0.5,
                height: 0.8,
                density: Some(1.),
                stiffness: Some([1.2e6, 38400., 15000., 22572.8]),
            },
        );
        materials.insert(
            7,
            Material::Plate {
                e: 3e6,
                nu: 0.2,
                thickness: 0.2,
                density: Some(2.5),
                membrane: Some(1.),
                bending: Some(0.6),
                shear: [None, None],
            },
        );
        let settings = |stiffness| Settings {
            force_factor: 1.,
            min_edge: 0.05,
            stiffness,
        };
        let x = exchange(
            &state,
            &materials,
            settings(StiffnessMode::Effective),
            "test",
        );
        let m = &x.beam_materials[0];
        assert!((m.a - 0.4).abs() < 1e-12);
        // Height along PLAXIS local axis 2 (LIRA Z1): I3 = EIy / E.
        assert!((m.i3 - 38400. / 3e6).abs() < 1e-12 && (m.i2 - 15000. / 3e6).abs() < 1e-12);
        assert!((m.gamma * m.a - 1.).abs() < 1e-12);
        // Z1 of a horizontal bar points up; of a vertical one along X.
        assert_eq!(x.beams[0].axis2, [0., 0., 1.]);
        assert_eq!(x.beams[1].axis2, [1., 0., 0.]);
        let p = &x.plate_materials[0];
        let (e, d) = (p.e, p.d);
        assert!((e * d - 3e6 * 0.2).abs() < 1e-6, "membrane E d kept");
        assert!(
            (e * d.powi(3) - 0.6 * 3e6 * 0.2f64.powi(3)).abs() < 1e-6,
            "bending E d^3 x PLKE"
        );
        assert!(
            (p.gamma * d - 2.5 * 0.2).abs() < 1e-12,
            "weight per area kept"
        );
        let n = exchange(&state, &materials, settings(StiffnessMode::Nominal), "test");
        let m = &n.beam_materials[0];
        assert!((m.i3 - 0.5 * 0.8f64.powi(3) / 12.).abs() < 1e-12);
        assert!((m.i2 - 0.8 * 0.5f64.powi(3) / 12.).abs() < 1e-12);
        assert!(!m.notes.is_empty());
        assert_eq!(n.plate_materials[0].d, 0.2);
    }

    #[test]
    fn polygons_never_start_with_collinear_points() {
        // A wall whose contour runs through junction vertices on its edges
        // (the PLAXIS "points are not coplanar" case).
        let ring = vec![
            [77.1132, -52.1233, -2.95545],
            [77.1132, -52.1233, -1.33],
            [77.1132, -52.1233, 6.45],
            [79.3073, -46.6925, 6.45],
            [79.4007, -46.4615, 6.45],
            [79.4007, -46.4614, -1.33],
            [79.4007, -46.4614, -2.95545],
        ];
        let out = plaxis_polygon(ring);
        let p = |v: [f64; 3]| DVec3::from_array(v);
        assert!(out.len() < 7);
        let n = out.len();
        for i in 0..n {
            let (a, b, c) = (p(out[(i + n - 1) % n]), p(out[i]), p(out[(i + 1) % n]));
            assert!(
                (b - a).cross(c - b).length() > 1e-6,
                "straight vertex {i} left"
            );
        }
        let (a, b, c) = (p(out[0]), p(out[1]), p(out[2]));
        assert!(
            (b - a).cross(c - a).length() > 1.,
            "first triangle degenerate"
        );
    }
}
