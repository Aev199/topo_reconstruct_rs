//! Removal of redundant contour vertices that create short edges.
//!
//! A source node can remain on a straight contour or junction line after the
//! structure it belonged to was aligned elsewhere, for example a slab node a
//! few millimetres from the wall end now standing on the slab edge. Such a
//! vertex joins two collinear edges, carries nothing and forces tiny
//! elements. It is removed when it bounds an edge shorter than the tolerance;
//! the geometry of every surface is unchanged.
//!
//! Distinct vertices at one location (duplicated source nodes) and a wall end
//! stopping next to a contour corner are merged into one vertex: the
//! geotechnical model is deliberately simplified, and every merge is
//! recorded with its source nodes.
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

#[derive(Debug, Serialize)]
pub struct MergedVertex {
    pub kind: String,
    /// Vertex that disappeared and the vertex that stays.
    pub dropped: usize,
    pub kept: usize,
    pub dropped_source_node: Option<u32>,
    pub kept_source_node: Option<u32>,
    pub distance: f64,
    /// Movement of the kept vertex onto the planes of both.
    pub kept_movement: f64,
    /// Position of the kept vertex before the merge.
    pub kept_from: [f64; 3],
    pub surfaces: Vec<usize>,
}

#[derive(Debug, Default, Serialize)]
pub struct MergeReport {
    pub tolerance: f64,
    pub merged: Vec<MergedVertex>,
    /// Candidate pairs that could not be merged validly, with the reason.
    pub rejected: Vec<RejectedMerge>,
}

#[derive(Debug, Serialize)]
pub struct RejectedMerge {
    pub vertices: [usize; 2],
    pub source_nodes: [Option<u32>; 2],
    pub distance: f64,
    pub reason: String,
}

fn users(model: &Model, v: usize) -> Vec<usize> {
    (0..model.surfaces.len())
        .filter(|&s| model.surface_edges(s).any(|e| model.edges[e].contains(&v)))
        .collect()
}

/// Move `keep` onto the planes of every surface of both vertices (within
/// `limit`, not at all if locked), then merge `drop` into it.
fn merge(
    model: &mut Model,
    drop: usize,
    keep: usize,
    limit: f64,
    locked: &BTreeSet<usize>,
) -> Result<f64, String> {
    let from = DVec3::from_array(model.vertices[keep]);
    let mut surfaces = users(model, drop);
    surfaces.extend(users(model, keep));
    surfaces.sort_unstable();
    surfaces.dedup();
    let planes: Vec<_> = surfaces
        .iter()
        .map(|&s| model.planes[model.surfaces[s].plane].clone())
        .collect();
    let refs: Vec<_> = planes.iter().collect();
    let target = super::intersection(from, &refs, model.precision).ok_or("inconsistent_planes")?;
    let movement = target.distance(from);
    if movement > limit || target.distance(DVec3::from_array(model.vertices[drop])) > limit {
        return Err(format!("movement_beyond_tolerance: {movement:e}"));
    }
    if locked.contains(&keep) && movement > model.precision {
        return Err(format!("locked_vertex_would_move: {movement:e}"));
    }
    if movement > 0. {
        model
            .move_vertex(keep, target.to_array())
            .map_err(|e| format!("move_{e:?}"))?;
    }
    if let Err(e) = model.merge_vertices(drop, keep) {
        if movement > 0. {
            model
                .move_vertex(keep, from.to_array())
                .expect("restoring a validated vertex position");
        }
        return Err(format!("merge_{e:?}"));
    }
    Ok(movement)
}

/// Merge distinct vertices closer than `tolerance` (duplicated source nodes).
pub fn merge_coincident(
    model: &mut Model,
    tolerance: f64,
    locked: &BTreeSet<usize>,
    source_nodes: &[u32],
) -> MergeReport {
    let mut report = MergeReport {
        tolerance,
        ..Default::default()
    };
    let used: BTreeSet<usize> = (0..model.surfaces.len())
        .flat_map(|s| model.surface_edges(s).collect::<Vec<_>>())
        .flat_map(|e| model.edges[e])
        .collect();
    let cell = tolerance.max(model.precision);
    let key = |p: [f64; 3]| p.map(|x| (x / cell).floor() as i64);
    let mut grid = std::collections::BTreeMap::<[i64; 3], Vec<usize>>::new();
    for &v in &used {
        grid.entry(key(model.vertices[v])).or_default().push(v);
    }
    let mut pairs = vec![];
    for &a in &used {
        let k = key(model.vertices[a]);
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    for &b in grid
                        .get(&[k[0] + dx, k[1] + dy, k[2] + dz])
                        .into_iter()
                        .flatten()
                    {
                        let d = DVec3::from_array(model.vertices[a])
                            .distance(DVec3::from_array(model.vertices[b]));
                        if a < b && d <= tolerance {
                            pairs.push((d, a, b));
                        }
                    }
                }
            }
        }
    }
    pairs.sort_by(|x, y| x.0.total_cmp(&y.0).then(x.1.cmp(&y.1)).then(x.2.cmp(&y.2)));
    let mut gone = BTreeSet::new();
    for (d, a, b) in pairs {
        if gone.contains(&a) || gone.contains(&b) {
            continue;
        }
        // Keep a locked vertex, otherwise the one with more edges.
        let degree = |v: usize| users(model, v).len();
        let (drop, keep) = match (locked.contains(&a), locked.contains(&b)) {
            (true, true) => {
                report.rejected.push(RejectedMerge {
                    vertices: [a, b],
                    source_nodes: [source_nodes.get(a).copied(), source_nodes.get(b).copied()],
                    distance: d,
                    reason: "both_locked".into(),
                });
                continue;
            }
            (true, false) => (b, a),
            (false, true) => (a, b),
            _ if degree(b) > degree(a) => (a, b),
            _ => (b, a),
        };
        let surfaces = {
            let mut s = users(model, drop);
            s.extend(users(model, keep));
            s.sort_unstable();
            s.dedup();
            s
        };
        let kept_from = model.vertices[keep];
        match merge(model, drop, keep, tolerance, locked) {
            Ok(movement) => {
                gone.insert(drop);
                report.merged.push(MergedVertex {
                    kind: "coincident".into(),
                    dropped: drop,
                    kept: keep,
                    dropped_source_node: source_nodes.get(drop).copied(),
                    kept_source_node: source_nodes.get(keep).copied(),
                    distance: d,
                    kept_movement: movement,
                    kept_from,
                    surfaces,
                });
            }
            Err(reason) => report.rejected.push(RejectedMerge {
                vertices: [drop, keep],
                source_nodes: [
                    source_nodes.get(drop).copied(),
                    source_nodes.get(keep).copied(),
                ],
                distance: d,
                reason,
            }),
        }
    }
    report
}

/// Merge the free end of a junction line (a wall end inside a surface) into
/// a nearby vertex of the same surface, within `tolerance`.
pub fn merge_wall_ends(
    model: &mut Model,
    tolerance: f64,
    locked: &BTreeSet<usize>,
    source_nodes: &[u32],
) -> MergeReport {
    let mut report = MergeReport {
        tolerance,
        ..Default::default()
    };
    let mut tried = BTreeSet::new();
    'restart: loop {
        for s in 0..model.surfaces.len() {
            let edges: Vec<usize> = model.surface_edges(s).collect();
            let mut degree = std::collections::BTreeMap::<usize, usize>::new();
            for &e in &edges {
                for v in model.edges[e] {
                    *degree.entry(v).or_default() += 1;
                }
            }
            let free: Vec<usize> = model.surfaces[s]
                .embedded_edges
                .iter()
                .flat_map(|&e| model.edges[e])
                .filter(|v| degree[v] == 1 && !locked.contains(v))
                .collect();
            for v in free {
                let p = DVec3::from_array(model.vertices[v]);
                let neighbour: BTreeSet<usize> = edges
                    .iter()
                    .filter(|&&e| model.edges[e].contains(&v))
                    .flat_map(|&e| model.edges[e])
                    .collect();
                let nearest = degree
                    .keys()
                    .copied()
                    .filter(|w| !neighbour.contains(w) && !tried.contains(&(v, *w)))
                    .map(|w| (DVec3::from_array(model.vertices[w]).distance(p), w))
                    .filter(|&(d, _)| d <= tolerance)
                    .min_by(|x, y| x.0.total_cmp(&y.0).then(x.1.cmp(&y.1)));
                let Some((d, w)) = nearest else {
                    continue;
                };
                tried.insert((v, w));
                let surfaces = {
                    let mut u = users(model, v);
                    u.extend(users(model, w));
                    u.sort_unstable();
                    u.dedup();
                    u
                };
                let kept_from = model.vertices[w];
                match merge(model, v, w, tolerance, locked) {
                    Ok(movement) => {
                        report.merged.push(MergedVertex {
                            kind: "wall_end_at_vertex".into(),
                            dropped: v,
                            kept: w,
                            dropped_source_node: source_nodes.get(v).copied(),
                            kept_source_node: source_nodes.get(w).copied(),
                            distance: d,
                            kept_movement: movement,
                            kept_from,
                            surfaces,
                        });
                        continue 'restart;
                    }
                    Err(reason) => report.rejected.push(RejectedMerge {
                        vertices: [v, w],
                        source_nodes: [source_nodes.get(v).copied(), source_nodes.get(w).copied()],
                        distance: d,
                        reason,
                    }),
                }
            }
        }
        return report;
    }
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
    use super::super::junctions::tests::{build, run, slab, wall, Placement};
    use super::*;
    use crate::reconstruction::PlaneFrame;

    /// Second slab panel sharing the line x = 2 with its own vertices.
    fn duplicate_panel(m: &mut Model, place: &Placement) {
        let n = place.rotation * DVec3::Z;
        let n = if place.flip { -n } else { n };
        let plane = m.add_plane(PlaneFrame::new(place.point([2., 0., 0.]), n.to_array()).unwrap());
        let ring = [[2., 0., 0.], [4., 0., 0.], [4., 4., 0.], [2., 4., 0.]]
            .map(|p| m.add_vertex(place.point(p)).unwrap())
            .to_vec();
        m.add_surface(plane, vec![ring], vec![]).unwrap();
    }

    #[test]
    fn duplicated_nodes_are_merged_into_shared_topology() {
        for place in Placement::all() {
            let mut m = build(&place, &[slab(0., 2.)]);
            duplicate_panel(&mut m, &place);
            let r = merge_coincident(&mut m, 0.001 * place.scale, &BTreeSet::new(), &[]);
            assert_eq!(r.merged.len(), 2, "{:?}", r.rejected);
            let shared = m
                .surface_edges(0)
                .filter(|&e| m.surface_edges(1).any(|f| f == e))
                .count();
            assert_eq!(shared, 1);
            // Idempotent.
            let again = merge_coincident(&mut m, 0.001 * place.scale, &BTreeSet::new(), &[]);
            assert!(again.merged.is_empty() && again.rejected.is_empty());
        }
    }

    #[test]
    fn duplicated_bar_anchors_are_not_merged() {
        let place = &Placement::all()[2];
        let mut m = build(place, &[slab(0., 2.)]);
        duplicate_panel(&mut m, place);
        let locked: BTreeSet<usize> = (0..m.vertices.len()).collect();
        let r = merge_coincident(&mut m, 0.001, &locked, &[]);
        assert!(r.merged.is_empty());
        assert!(r.rejected.iter().all(|x| x.reason == "both_locked"));
    }

    #[test]
    fn wall_end_next_to_a_contour_corner_is_merged_into_it() {
        for place in Placement::all() {
            // Slab edge x = 4 has a vertex 10 mm off the wall axis y = 2; the
            // wall ends 30 mm before that edge.
            let panel = (
                vec![vec![
                    [0., 0., 0.],
                    [4., 0., 0.],
                    [4., 2.01, 0.],
                    [4., 4., 0.],
                    [0., 4., 0.],
                ]],
                [0., 0., 1.],
            );
            let mut m = build(&place, &[panel, wall(1., 3.97, 0., 2.)]);
            let j = run(&mut m);
            assert!(j.issues.is_empty(), "{:?}", j.issues);
            let r = merge_wall_ends(&mut m, 0.05 * place.scale, &BTreeSet::new(), &[]);
            assert_eq!(r.merged.len(), 1, "{:?}", r.rejected);
            let merged = &r.merged[0];
            assert!((merged.kept_movement - 0.01 * place.scale).abs() < 1e-9 * place.scale);
            let corner = DVec3::from_array(m.vertices[merged.kept]);
            assert!(
                corner.distance(DVec3::from_array(place.point([4., 2., 0.])))
                    < 1e-9 * place.scale.max(1.)
            );
            // The wall base now reaches the slab contour at the corner.
            assert!(m.surfaces[1].boundaries[0]
                .iter()
                .any(|e| m.edges[e.edge].contains(&merged.kept)));
            // Without a tolerance nothing moves.
            let mut m = build(&place, &[slab(0., 4.), wall(1., 3.97, 0., 2.)]);
            run(&mut m);
            assert!(merge_wall_ends(&mut m, 0., &BTreeSet::new(), &[])
                .merged
                .is_empty());
        }
    }

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
