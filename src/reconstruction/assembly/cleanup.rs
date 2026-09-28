//! Removal of redundant contour vertices that create short edges.
//!
//! A source node can remain on a straight contour or junction line after the
//! structure it belonged to was aligned elsewhere, for example a slab node a
//! few millimetres from the wall end now standing on the slab edge. Such a
//! vertex joins two collinear edges, carries nothing and forces tiny
//! elements. It is removed when it bounds an edge shorter than the tolerance;
//! the geometry of every surface is unchanged.
use crate::reconstruction::Model;
use glam::DVec3;
use serde::Serialize;
use std::collections::BTreeSet;

#[derive(Debug, Serialize)]
pub struct RemovedVertex {
    pub vertex: usize,
    /// Source node of the vertex, if it has one.
    pub source_node: Option<u32>,
    pub surfaces: Vec<usize>,
    /// The short edge that motivated the removal.
    pub short_edge_length: f64,
    /// Distance of the removed vertex from the joined segment.
    pub deviation: f64,
}

#[derive(Debug, Default, Serialize)]
pub struct Report {
    pub tolerance: f64,
    /// Largest accepted distance of a removed vertex from the straight line
    /// of its neighbours: 1/1000 of the tolerance.
    pub maximum_deviation: f64,
    pub removed: Vec<RemovedVertex>,
    /// Short edges left because neither end is a removable vertex.
    pub remaining_short_edges: usize,
}

/// Accepted straightness defect relative to the cleanup tolerance.
const DEVIATION_RATIO: f64 = 1e-3;

fn length(model: &Model, e: usize) -> f64 {
    let [a, b] = model.edges[e].map(|v| DVec3::from_array(model.vertices[v]));
    a.distance(b)
}

/// Remove redundant collinear vertices bounding edges shorter than
/// `tolerance`. `locked` vertices (bar anchors, retained nodes) stay.
pub fn remove_short_edges(
    model: &mut Model,
    tolerance: f64,
    locked: &BTreeSet<usize>,
    source_nodes: &[u32],
) -> Report {
    let deviation = tolerance * DEVIATION_RATIO;
    let mut report = Report {
        tolerance,
        maximum_deviation: deviation,
        ..Default::default()
    };
    let mut failed = BTreeSet::new();
    loop {
        let used: BTreeSet<usize> = (0..model.surfaces.len())
            .flat_map(|s| model.surface_edges(s).collect::<Vec<_>>())
            .collect();
        let mut short: Vec<(f64, usize)> = used
            .iter()
            .map(|&e| (length(model, e), e))
            .filter(|&(l, _)| l < tolerance)
            .collect();
        short.sort_by(|x, y| x.0.total_cmp(&y.0));
        let degree = |v: usize| {
            used.iter()
                .filter(|&&e| model.edges[e].contains(&v))
                .count()
        };
        let other_edge = |v: usize, e: usize| {
            used.iter()
                .copied()
                .find(|&f| f != e && model.edges[f].contains(&v))
        };
        // Candidate removals in order: shortest edge first; per edge, the end
        // whose other edge is shorter, so one removal repairs both pieces.
        let mut plan: Vec<(f64, usize, usize)> = vec![];
        for &(l, e) in &short {
            let mut ends: Vec<(f64, usize)> = model.edges[e]
                .into_iter()
                .filter(|v| !locked.contains(v) && !failed.contains(v) && degree(*v) == 2)
                .filter_map(|v| other_edge(v, e).map(|f| (length(model, f), v)))
                .collect();
            ends.sort_by(|x, y| x.0.total_cmp(&y.0));
            plan.extend(ends.into_iter().map(|(_, v)| (l, e, v)));
        }
        let mut changed = false;
        for (l, e, v) in plan {
            let surfaces: Vec<usize> = (0..model.surfaces.len())
                .filter(|&s| model.surface_edges(s).any(|x| x == e))
                .collect();
            let offset = {
                let ends: Vec<usize> = model
                    .edges
                    .iter()
                    .filter(|x| x.contains(&v))
                    .flat_map(|x| x.iter().copied().filter(|&w| w != v))
                    .collect();
                let p = DVec3::from_array(model.vertices[v]);
                match ends[..] {
                    [a, b] => {
                        let (a, b) = (
                            DVec3::from_array(model.vertices[a]),
                            DVec3::from_array(model.vertices[b]),
                        );
                        (p - a).cross(b - a).length() / a.distance(b)
                    }
                    _ => 0.,
                }
            };
            if model.remove_vertex(v, deviation).is_ok() {
                report.removed.push(RemovedVertex {
                    vertex: v,
                    source_node: source_nodes.get(v).copied(),
                    surfaces,
                    short_edge_length: l,
                    deviation: offset,
                });
                changed = true;
                break;
            }
            failed.insert(v);
        }
        if !changed {
            report.remaining_short_edges = short.len();
            return report;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::super::junctions::tests::{build, Placement};
    use super::*;

    #[test]
    fn collinear_vertex_near_a_corner_is_removed_without_moving_geometry() {
        for place in Placement::all() {
            // A slab edge carries a stray node 12 mm from the corner.
            let slab = (
                vec![vec![
                    [0., 0., 0.],
                    [4., 0., 0.],
                    [4., 3.988, 0.],
                    [4., 4., 0.],
                    [0., 4., 0.],
                ]],
                [0., 0., 1.],
            );
            let mut m = build(&place, &[slab]);
            let before = m.surfaces[0].contours[0].len();
            let r = remove_short_edges(&mut m, 0.05 * place.scale, &BTreeSet::new(), &[]);
            assert_eq!(r.removed.len(), 1);
            assert_eq!(m.surfaces[0].contours[0].len(), before - 1);
            assert_eq!(r.remaining_short_edges, 0);
            // A locked vertex (e.g. a bar anchor) stays.
            let mut m = build(&place, &[r_slab()]);
            let v = (0..m.vertices.len())
                .find(|&v| {
                    DVec3::from_array(m.vertices[v])
                        .distance(DVec3::from_array(place.point([4., 3.988, 0.])))
                        < 1e-9 * place.scale.max(1.)
                })
                .unwrap();
            let r = remove_short_edges(&mut m, 0.05 * place.scale, &BTreeSet::from([v]), &[]);
            assert!(r.removed.is_empty());
            assert_eq!(r.remaining_short_edges, 1);
        }
    }

    fn r_slab() -> (Vec<Vec<[f64; 3]>>, [f64; 3]) {
        (
            vec![vec![
                [0., 0., 0.],
                [4., 0., 0.],
                [4., 3.988, 0.],
                [4., 4., 0.],
                [0., 4., 0.],
            ]],
            [0., 0., 1.],
        )
    }

    #[test]
    fn a_real_corner_is_never_removed() {
        for place in Placement::all() {
            // A 12 mm notch: its vertices are corners, not collinear.
            let slab = (
                vec![vec![
                    [0., 0., 0.],
                    [4., 0., 0.],
                    [4., 3.988, 0.],
                    [3.9, 3.988, 0.],
                    [3.9, 4., 0.],
                    [0., 4., 0.],
                ]],
                [0., 0., 1.],
            );
            let mut m = build(&place, &[slab]);
            let before = serde_json::to_string(&m).unwrap();
            let r = remove_short_edges(&mut m, 0.05 * place.scale, &BTreeSet::new(), &[]);
            assert!(r.removed.is_empty());
            assert_eq!(before, serde_json::to_string(&m).unwrap());
        }
    }
}
