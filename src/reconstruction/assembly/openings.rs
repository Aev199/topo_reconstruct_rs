//! Filling of small openings for a geotechnical model.
//!
//! An opening narrower than the minimum width (the smaller side of its
//! minimum bounding rectangle: a 0.9 m door, a 0.5 m duct) carries no
//! geotechnical meaning but forces a fine mesh around it. It is filled
//! (its ring is dropped from the surface) only when it is free: no other
//! surface or bar is attached to its contour and nothing passes through it.
//! An opening held by structures (a shaft surrounded by walls) is kept.
use super::bars::Axis;
use crate::reconstruction::Model;
use glam::{DVec2, DVec3};
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet, HashMap};

#[derive(Debug, Serialize)]
pub struct FilledOpening {
    pub surface: usize,
    /// Smaller side of the minimum bounding rectangle.
    pub width: f64,
    /// Larger side of the minimum bounding rectangle.
    pub length: f64,
    pub area: f64,
    /// Source nodes of the dropped contour.
    pub source_nodes: Vec<u32>,
}

#[derive(Debug, Serialize)]
pub struct KeptOpening {
    pub surface: usize,
    pub width: f64,
    pub length: f64,
    /// `long`, `attached_surface`, `attached_bar`, `occupied` or
    /// `invalid_<error>`.
    pub reason: String,
}

#[derive(Debug, Default, Serialize)]
pub struct Report {
    pub minimum_width: f64,
    pub maximum_length: f64,
    pub filled: Vec<FilledOpening>,
    pub kept: Vec<KeptOpening>,
}

/// Smaller side of the minimum-area bounding rectangle of `points`
/// (rotating calipers over the convex hull).
pub fn minimum_width(points: &[DVec2]) -> f64 {
    rectangle_sides(points)[0]
}

/// Smaller and larger side of the minimum-area bounding rectangle of
/// `points` (rotating calipers over the convex hull).
pub fn rectangle_sides(points: &[DVec2]) -> [f64; 2] {
    let mut pts: Vec<DVec2> = points.to_vec();
    pts.sort_by(|a, b| a.x.total_cmp(&b.x).then(a.y.total_cmp(&b.y)));
    pts.dedup();
    if pts.len() < 3 {
        return [0., 0.];
    }
    let cross = |o: DVec2, a: DVec2, b: DVec2| (a - o).perp_dot(b - o);
    let mut hull: Vec<DVec2> = vec![];
    for pass in 0..2 {
        let start = hull.len();
        let iter: Box<dyn Iterator<Item = &DVec2>> = if pass == 0 {
            Box::new(pts.iter())
        } else {
            Box::new(pts.iter().rev())
        };
        for &p in iter {
            while hull.len() >= start + 2
                && cross(hull[hull.len() - 2], hull[hull.len() - 1], p) <= 0.
            {
                hull.pop();
            }
            hull.push(p);
        }
        hull.pop();
    }
    let mut best = [f64::MAX; 2];
    let mut best_area = f64::MAX;
    for i in 0..hull.len() {
        let (a, b) = (hull[i], hull[(i + 1) % hull.len()]);
        let d = (b - a).normalize_or_zero();
        if d == DVec2::ZERO {
            continue;
        }
        let n = d.perp();
        let (mut lo, mut hi, mut far) = (f64::MAX, f64::MIN, 0f64);
        for &p in &hull {
            let t = (p - a).dot(d);
            lo = lo.min(t);
            hi = hi.max(t);
            far = far.max((p - a).dot(n).abs());
        }
        let area = (hi - lo) * far;
        if area < best_area {
            best_area = area;
            best = [(hi - lo).min(far), (hi - lo).max(far)];
        }
    }
    best
}

fn inside(polygon: &[DVec2], p: DVec2) -> bool {
    let mut odd = false;
    for i in 0..polygon.len() {
        let (a, b) = (polygon[i], polygon[(i + 1) % polygon.len()]);
        if (a.y > p.y) != (b.y > p.y) && p.x < a.x + (p.y - a.y) / (b.y - a.y) * (b.x - a.x) {
            odd = !odd;
        }
    }
    odd
}

fn crosses(polygon: &[DVec2], a: DVec2, b: DVec2) -> bool {
    (0..polygon.len()).any(|i| {
        let (c, d) = (polygon[i], polygon[(i + 1) % polygon.len()]);
        let side = |p: DVec2, q: DVec2, r: DVec2| (q - p).perp_dot(r - p);
        side(a, b, c) * side(a, b, d) < 0. && side(c, d, a) * side(c, d, b) < 0.
    })
}

fn area(polygon: &[DVec2]) -> f64 {
    (0..polygon.len())
        .map(|i| polygon[i].perp_dot(polygon[(i + 1) % polygon.len()]))
        .sum::<f64>()
        .abs()
        / 2.
}

/// Fill openings narrower than `minimum_width` and no longer than
/// `maximum_length` (a long slot is a structural feature, not a door or a
/// duct) that nothing is attached to or passes through. `locked` vertices
/// (bar anchors, contacts) hold an opening.
pub fn fill(
    model: &mut Model,
    minimum_width_limit: f64,
    maximum_length: f64,
    locked: &BTreeSet<usize>,
    axes: &[Axis],
    source_nodes: &[u32],
) -> Report {
    let mut report = Report {
        minimum_width: minimum_width_limit,
        maximum_length,
        ..Default::default()
    };
    if minimum_width_limit <= 0. {
        return report;
    }
    let at = |m: &Model, v: usize| DVec3::from_array(m.vertices[v]);
    // Users of every edge and vertex; every used edge on a grid.
    let mut edge_users = BTreeMap::<usize, BTreeSet<usize>>::new();
    for s in 0..model.surfaces.len() {
        for e in model.surface_edges(s) {
            edge_users.entry(e).or_default().insert(s);
        }
    }
    let mut vertex_users = BTreeMap::<usize, BTreeSet<usize>>::new();
    for (&e, list) in &edge_users {
        for v in model.edges[e] {
            vertex_users.entry(v).or_default().extend(list);
        }
    }
    let mut bar_vertices: BTreeSet<usize> = locked.clone();
    let mut segments: Vec<[usize; 2]> = edge_users.keys().map(|&e| model.edges[e]).collect();
    for axis in axes {
        bar_vertices.extend(axis.endpoints);
        bar_vertices.extend(axis.anchors.iter().map(|a| a.vertex));
        segments.push(axis.endpoints);
    }
    let cell = minimum_width_limit.max(model.minimum_edge);
    let key = |p: DVec3| (p / cell).floor().as_i64vec3().to_array();
    let mut grid = HashMap::<[i64; 3], Vec<usize>>::new();
    for (k, &[a, b]) in segments.iter().enumerate() {
        let (pa, pb) = (at(model, a), at(model, b));
        let (lo, hi) = (key(pa.min(pb)), key(pa.max(pb)));
        for x in lo[0]..=hi[0] {
            for y in lo[1]..=hi[1] {
                for z in lo[2]..=hi[2] {
                    grid.entry([x, y, z]).or_default().push(k);
                }
            }
        }
    }
    for s in 0..model.surfaces.len() {
        // Holes from the last, so ring indices stay valid after a fill.
        for r in (1..model.surfaces[s].boundaries.len()).rev() {
            let surface = &model.surfaces[s];
            let ring: Vec<usize> = surface.boundaries[r]
                .iter()
                .map(|u| {
                    let [a, b] = model.edges[u.edge];
                    if u.reversed {
                        b
                    } else {
                        a
                    }
                })
                .collect();
            let polygon: Vec<DVec2> = surface.contours[r]
                .iter()
                .map(|&p| DVec2::from_array(p))
                .collect();
            let [width, length] = rectangle_sides(&polygon);
            if width >= minimum_width_limit {
                continue;
            }
            let kept = |reason: &str| KeptOpening {
                surface: s,
                width,
                length,
                reason: reason.into(),
            };
            if length > maximum_length {
                report.kept.push(kept("long"));
                continue;
            }
            if ring.iter().any(|v| bar_vertices.contains(v)) {
                report.kept.push(kept("attached_bar"));
                continue;
            }
            if ring.iter().any(|v| {
                vertex_users
                    .get(v)
                    .is_some_and(|u| u.iter().any(|&t| t != s))
            }) {
                report.kept.push(kept("attached_surface"));
                continue;
            }
            // Anything else meeting the opening: an edge (or bar) through
            // it or lying in it.
            let plane = &model.planes[surface.plane];
            let normal = DVec3::from_array(plane.normal);
            let origin = at(model, ring[0]);
            let height = |p: DVec3| (p - origin).dot(normal);
            let uv = |p: DVec3| DVec2::from_array(plane.project(p.to_array()));
            let points: Vec<DVec3> = ring.iter().map(|&v| at(model, v)).collect();
            let lo = points.iter().fold(DVec3::MAX, |m, &p| m.min(p)) - model.minimum_edge;
            let hi = points.iter().fold(DVec3::MIN, |m, &p| m.max(p)) + model.minimum_edge;
            let own: BTreeSet<[usize; 2]> =
                model.surface_edges(s).map(|e| model.edges[e]).collect();
            let tol = model.minimum_edge;
            let (klo, khi) = (key(lo), key(hi));
            let mut seen = BTreeSet::new();
            let mut occupied = false;
            'cells: for x in klo[0]..=khi[0] {
                for y in klo[1]..=khi[1] {
                    for z in klo[2]..=khi[2] {
                        for &k in grid.get(&[x, y, z]).into_iter().flatten() {
                            if !seen.insert(k) || own.contains(&segments[k]) {
                                continue;
                            }
                            let [a, b] = segments[k];
                            let (pa, pb) = (at(model, a), at(model, b));
                            let (ha, hb) = (height(pa), height(pb));
                            let hit = if ha.abs() <= tol && hb.abs() <= tol {
                                inside(&polygon, uv(pa))
                                    || inside(&polygon, uv(pb))
                                    || crosses(&polygon, uv(pa), uv(pb))
                            } else if ha.abs() <= tol {
                                inside(&polygon, uv(pa))
                            } else if hb.abs() <= tol {
                                inside(&polygon, uv(pb))
                            } else if ha * hb < 0. {
                                inside(&polygon, uv(pa + (pb - pa) * (ha / (ha - hb))))
                            } else {
                                false
                            };
                            if hit {
                                occupied = true;
                                break 'cells;
                            }
                        }
                    }
                }
            }
            if occupied {
                report.kept.push(kept("occupied"));
                continue;
            }
            let rings: Vec<Vec<usize>> = (0..surface.boundaries.len())
                .filter(|&i| i != r)
                .map(|i| {
                    surface.boundaries[i]
                        .iter()
                        .map(|u| {
                            let [a, b] = model.edges[u.edge];
                            if u.reversed {
                                b
                            } else {
                                a
                            }
                        })
                        .collect()
                })
                .collect();
            let embedded = surface.embedded_edges.clone();
            let filled_area = area(&polygon);
            match model.rebuild_surface(s, rings, embedded) {
                Ok(()) => report.filled.push(FilledOpening {
                    surface: s,
                    width,
                    length,
                    area: filled_area,
                    source_nodes: ring
                        .iter()
                        .filter_map(|&v| source_nodes.get(v).copied())
                        .collect(),
                }),
                Err(e) => report.kept.push(kept(&format!("invalid_{e:?}"))),
            }
        }
    }
    report
}

#[cfg(test)]
mod tests {
    use super::super::junctions::tests::{build, wall, Placement};
    use super::*;

    fn slab_with_hole(x0: f64, y0: f64, w: f64, h: f64) -> (Vec<Vec<[f64; 3]>>, [f64; 3]) {
        (
            vec![
                vec![[0., 0., 0.], [6., 0., 0.], [6., 4., 0.], [0., 4., 0.]],
                vec![
                    [x0, y0, 0.],
                    [x0, y0 + h, 0.],
                    [x0 + w, y0 + h, 0.],
                    [x0 + w, y0, 0.],
                ],
            ],
            [0., 0., 1.],
        )
    }

    #[test]
    fn minimum_width_is_the_smaller_side_of_the_bounding_rectangle() {
        let rotated: Vec<DVec2> = [[0., 0.], [3., 0.], [3., 0.8], [0., 0.8]]
            .iter()
            .map(|&p| DVec2::from_angle(0.4).rotate(DVec2::from_array(p)))
            .collect();
        assert!((minimum_width(&rotated) - 0.8).abs() < 1e-9);
    }

    #[test]
    fn only_free_narrow_openings_are_filled() {
        for place in Placement::all() {
            // A 0.9 m x 2 m opening is filled, a 1.2 m x 1.2 m one kept.
            for (w, filled) in [(0.9, 1), (1.2, 0)] {
                let mut m = build(&place, &[slab_with_hole(2., 1., w, 2.)]);
                let r = fill(
                    &mut m,
                    1.0 * place.scale,
                    3.0 * place.scale,
                    &BTreeSet::new(),
                    &[],
                    &[],
                );
                assert_eq!(r.filled.len(), filled, "{r:?}");
                assert_eq!(m.surfaces[0].boundaries.len(), 2 - filled);
                // The 1.2 m opening is wide enough: not even reported.
                assert!(r.kept.is_empty());
            }
            // A 0.5 m x 3.5 m slot is longer than the limit: kept.
            let mut m = build(&place, &[slab_with_hole(2., 0.25, 0.5, 3.5)]);
            let r = fill(
                &mut m,
                1.0 * place.scale,
                3.0 * place.scale,
                &BTreeSet::new(),
                &[],
                &[],
            );
            assert!(r.filled.is_empty());
            assert_eq!(r.kept[0].reason, "long");
            assert!((r.kept[0].length - 3.5 * place.scale).abs() < 1e-9 * place.scale.max(1.));
            // A wall passing through the opening holds it.
            let mut m = build(
                &place,
                &[slab_with_hole(2., 1., 0.9, 2.), wall(2.2, 2.6, -1., 1.)],
            );
            let r = fill(
                &mut m,
                1.0 * place.scale,
                3.0 * place.scale,
                &BTreeSet::new(),
                &[],
                &[],
            );
            assert!(r.filled.is_empty());
            assert_eq!(r.kept[0].reason, "occupied");
        }
    }
}
