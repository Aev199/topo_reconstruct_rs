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
use super::bars::{Axis, Contact};
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

/// Bar axes whose anchor vertices may take part in a merge.
pub struct Bars<'a> {
    pub axes: &'a mut Vec<Axis>,
    pub contacts: &'a mut Vec<Contact>,
}

/// Move `v` to `q`. A bar anchor moves only as the end of its axes: the
/// other end stays and every other anchor is re-placed at its parameter on
/// the new straight line, keeping all its planes. An interior anchor, or an
/// anchor shared with another axis, would bend a bar and is rejected.
fn move_with_axes(
    model: &mut Model,
    axes: &[Axis],
    v: usize,
    q: DVec3,
    limit: f64,
) -> Result<(), String> {
    if q.distance(DVec3::from_array(model.vertices[v])) <= model.precision {
        return Ok(());
    }
    let mut moves = vec![(v, q)];
    for (i, axis) in axes.iter().enumerate() {
        let member = axis.endpoints.contains(&v) || axis.anchors.iter().any(|a| a.vertex == v);
        if !member {
            continue;
        }
        let Some(end) = axis.endpoints.iter().position(|&e| e == v) else {
            return Err("interior_bar_anchor".into());
        };
        let mut ends = axis.endpoints.map(|e| DVec3::from_array(model.vertices[e]));
        ends[end] = q;
        for anchor in &axis.anchors {
            if axis.endpoints.contains(&anchor.vertex) {
                continue;
            }
            let shared = axes.iter().enumerate().any(|(j, other)| {
                j != i
                    && (other.endpoints.contains(&anchor.vertex)
                        || other.anchors.iter().any(|a| a.vertex == anchor.vertex))
            });
            if shared {
                return Err("shared_bar_anchor".into());
            }
            let p = ends[0].lerp(ends[1], anchor.t);
            if p.distance(DVec3::from_array(model.vertices[anchor.vertex])) > limit {
                return Err("bar_anchor_movement_beyond_tolerance".into());
            }
            moves.push((anchor.vertex, p));
        }
    }
    for (u, p) in moves {
        model
            .move_vertex(u, p.to_array())
            .map_err(|e| format!("move_{e:?}"))?;
    }
    Ok(())
}

/// Move both vertices onto the planes of every surface of either (within
/// `limit`), then merge `drop` into `keep`, carrying bar axes along.
fn merge(
    model: &mut Model,
    bars: &mut Bars<'_>,
    drop: usize,
    keep: usize,
    limit: f64,
    fixed: &BTreeSet<usize>,
) -> Result<f64, String> {
    if fixed.contains(&drop) || fixed.contains(&keep) {
        return Err("retained_node".into());
    }
    let from = DVec3::from_array(model.vertices[keep]);
    let mut surfaces = users(model, drop);
    surfaces.extend(users(model, keep));
    // A bar node inside a surface (point contact) stays in its plane.
    for c in bars.contacts.iter() {
        if let Contact::Point {
            vertex, surface, ..
        } = c
        {
            if *vertex == drop || *vertex == keep {
                surfaces.push(*surface);
            }
        }
    }
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
    if bars.axes.iter().any(|a| {
        let has = |v: usize| a.endpoints.contains(&v) || a.anchors.iter().any(|x| x.vertex == v);
        has(drop) && has(keep)
    }) {
        return Err("bar_would_collapse".into());
    }
    let mut trial = model.clone();
    move_with_axes(&mut trial, bars.axes, keep, target, limit)?;
    // Ends of one edge: the edge collapses in the merge itself.
    if model.edge_between(drop, keep).is_some() {
        if bars
            .axes
            .iter()
            .any(|a| a.endpoints.contains(&drop) || a.anchors.iter().any(|x| x.vertex == drop))
        {
            return Err("bar_node_on_collapsed_edge".into());
        }
    } else {
        move_with_axes(&mut trial, bars.axes, drop, target, limit)?;
    }
    trial
        .merge_vertices(drop, keep)
        .map_err(|e| format!("merge_{e:?}"))?;
    *model = trial;
    for axis in bars.axes.iter_mut() {
        for e in axis.endpoints.iter_mut() {
            if *e == drop {
                *e = keep;
            }
        }
        for a in axis.anchors.iter_mut() {
            if a.vertex == drop {
                a.vertex = keep;
            }
        }
    }
    for c in bars.contacts.iter_mut() {
        if let Contact::Point { vertex, .. } = c {
            if *vertex == drop {
                *vertex = keep;
            }
        }
    }
    Ok(movement)
}

/// Merge distinct vertices closer than `tolerance` (duplicated source nodes).
pub fn merge_coincident(
    model: &mut Model,
    bars: &mut Bars<'_>,
    tolerance: f64,
    fixed: &BTreeSet<usize>,
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
        // Keep the vertex with more surfaces; bar anchors are carried along.
        let degree = |v: usize| users(model, v).len();
        let (drop, keep) = if degree(b) > degree(a) {
            (a, b)
        } else {
            (b, a)
        };
        let surfaces = {
            let mut s = users(model, drop);
            s.extend(users(model, keep));
            s.sort_unstable();
            s.dedup();
            s
        };
        let kept_from = model.vertices[keep];
        match merge(model, bars, drop, keep, tolerance, fixed) {
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

/// Collapse surface edges shorter than `tolerance` whose two ends are both
/// needed corners (for example the ends of a lower and an upper wall a few
/// millimetres apart on one slab line): the ends merge into one vertex,
/// shortest edge first, keeping every plane and never bending a bar. Such
/// an edge cannot hold elements of any useful size.
pub fn collapse_short_edges(
    model: &mut Model,
    bars: &mut Bars<'_>,
    tolerance: f64,
    fixed: &BTreeSet<usize>,
    source_nodes: &[u32],
) -> MergeReport {
    let mut report = MergeReport {
        tolerance,
        ..Default::default()
    };
    let mut tried = BTreeSet::new();
    loop {
        let used: BTreeSet<usize> = (0..model.surfaces.len())
            .flat_map(|s| model.surface_edges(s).collect::<Vec<_>>())
            .collect();
        let shortest = used
            .iter()
            .map(|&e| {
                let [a, b] = model.edges[e];
                let d = DVec3::from_array(model.vertices[a])
                    .distance(DVec3::from_array(model.vertices[b]));
                (d, a.min(b), a.max(b))
            })
            .filter(|&(d, a, b)| d < tolerance && !tried.contains(&(a, b)))
            .min_by(|x, y| x.0.total_cmp(&y.0).then(x.1.cmp(&y.1)).then(x.2.cmp(&y.2)));
        let Some((d, a, b)) = shortest else {
            return report;
        };
        tried.insert((a, b));
        // Keep the vertex of more surfaces.
        let (drop, keep) = if users(model, a).len() > users(model, b).len() {
            (b, a)
        } else {
            (a, b)
        };
        let surfaces = {
            let mut s = users(model, drop);
            s.extend(users(model, keep));
            s.sort_unstable();
            s.dedup();
            s
        };
        let kept_from = model.vertices[keep];
        match merge(model, bars, drop, keep, tolerance, fixed) {
            Ok(movement) => report.merged.push(MergedVertex {
                kind: "short_edge".into(),
                dropped: drop,
                kept: keep,
                dropped_source_node: source_nodes.get(drop).copied(),
                kept_source_node: source_nodes.get(keep).copied(),
                distance: d,
                kept_movement: movement,
                kept_from,
                surfaces,
            }),
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
}

/// Merge the free end of a junction line (a wall end inside a surface) into
/// a nearby vertex of the same surface, within `tolerance`.
pub fn merge_wall_ends(
    model: &mut Model,
    bars: &mut Bars<'_>,
    tolerance: f64,
    fixed: &BTreeSet<usize>,
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
                .filter(|v| degree[v] == 1 && !fixed.contains(v))
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
                match merge(model, bars, v, w, tolerance, fixed) {
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

/// Identify bar nodes with surface vertices closer than `tolerance`.
///
/// A beam on a slab edge may keep its own node a fraction of a millimetre
/// from the slab contour vertex: two parallel lines then run side by side
/// and the mesh leaves the bar unconnected along them. The surface vertex
/// moves onto the bar node (a bar never bends), and the merge is recorded.
pub fn merge_bar_anchors(
    model: &mut Model,
    bars: &mut Bars<'_>,
    tolerance: f64,
    fixed: &BTreeSet<usize>,
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
    let anchors: BTreeSet<usize> = bars
        .axes
        .iter()
        .flat_map(|a| a.anchors.iter().map(|x| x.vertex).chain(a.endpoints))
        .filter(|v| !used.contains(v))
        .collect();
    let cell = tolerance.max(model.precision);
    let key = |p: [f64; 3]| p.map(|x| (x / cell).floor() as i64);
    let mut grid = std::collections::BTreeMap::<[i64; 3], Vec<usize>>::new();
    for &v in &used {
        grid.entry(key(model.vertices[v])).or_default().push(v);
    }
    let mut pairs = vec![];
    for &a in &anchors {
        let k = key(model.vertices[a]);
        let mut best: Option<(f64, usize)> = None;
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
                        if d <= tolerance && best.is_none_or(|x| (d, b) < x) {
                            best = Some((d, b));
                        }
                    }
                }
            }
        }
        if let Some((d, b)) = best {
            pairs.push((d, a, b));
        }
    }
    pairs.sort_by(|x, y| x.0.total_cmp(&y.0).then(x.1.cmp(&y.1)).then(x.2.cmp(&y.2)));
    let mut gone = BTreeSet::new();
    for (d, keep, drop) in pairs {
        if gone.contains(&drop) {
            continue;
        }
        let surfaces = users(model, drop);
        let kept_from = model.vertices[keep];
        match merge(model, bars, drop, keep, tolerance, fixed) {
            Ok(movement) => {
                gone.insert(drop);
                report.merged.push(MergedVertex {
                    kind: "bar_anchor".into(),
                    dropped: drop,
                    kept: keep,
                    dropped_source_node: source_nodes.get(drop).copied(),
                    kept_source_node: source_nodes.get(keep).copied(),
                    distance: d,
                    kept_movement: movement,
                    kept_from,
                    surfaces,
                })
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

/// Merge the ends of two different bars that stop within `tolerance` of each
/// other (a beam split by a small gap). Both bars stay straight; a bar end
/// near the interior of another bar is left for review.
pub fn merge_bar_ends(
    model: &mut Model,
    bars: &mut Bars<'_>,
    tolerance: f64,
    fixed: &BTreeSet<usize>,
    source_nodes: &[u32],
) -> MergeReport {
    let mut report = MergeReport {
        tolerance,
        ..Default::default()
    };
    let mut tried = BTreeSet::new();
    loop {
        let mut best: Option<(f64, usize, usize)> = None;
        for (i, a) in bars.axes.iter().enumerate() {
            for b in bars.axes.iter().skip(i + 1) {
                for &u in &a.endpoints {
                    for &w in &b.endpoints {
                        if u == w
                            || tried.contains(&(u.min(w), u.max(w)))
                            || a.endpoints.contains(&w)
                            || b.endpoints.contains(&u)
                        {
                            continue;
                        }
                        let d = DVec3::from_array(model.vertices[u])
                            .distance(DVec3::from_array(model.vertices[w]));
                        if d <= tolerance && best.is_none_or(|x| d < x.0) {
                            best = Some((d, u, w));
                        }
                    }
                }
            }
        }
        let Some((d, u, w)) = best else {
            return report;
        };
        tried.insert((u.min(w), u.max(w)));
        // Keep the end with more connections.
        let connections = |v: usize| {
            users(model, v).len()
                + bars
                    .axes
                    .iter()
                    .filter(|a| a.anchors.iter().any(|x| x.vertex == v))
                    .count()
        };
        let (drop, keep) = if connections(u) > connections(w) {
            (w, u)
        } else {
            (u, w)
        };
        let surfaces = {
            let mut s = users(model, drop);
            s.extend(users(model, keep));
            s.sort_unstable();
            s.dedup();
            s
        };
        let kept_from = model.vertices[keep];
        match merge(model, bars, drop, keep, tolerance, fixed) {
            Ok(movement) => report.merged.push(MergedVertex {
                kind: "bar_ends".into(),
                dropped: drop,
                kept: keep,
                dropped_source_node: source_nodes.get(drop).copied(),
                kept_source_node: source_nodes.get(keep).copied(),
                distance: d,
                kept_movement: movement,
                kept_from,
                surfaces,
            }),
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
    use super::super::bars::Anchor;
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
            let r = merge_coincident(
                &mut m,
                &mut no_bars(),
                0.001 * place.scale,
                &BTreeSet::new(),
                &[],
            );
            assert_eq!(r.merged.len(), 2, "{:?}", r.rejected);
            let shared = m
                .surface_edges(0)
                .filter(|&e| m.surface_edges(1).any(|f| f == e))
                .count();
            assert_eq!(shared, 1);
            // Idempotent.
            let again = merge_coincident(
                &mut m,
                &mut no_bars(),
                0.001 * place.scale,
                &BTreeSet::new(),
                &[],
            );
            assert!(again.merged.is_empty() && again.rejected.is_empty());
        }
    }

    fn no_bars() -> Bars<'static> {
        Bars {
            axes: Box::leak(Box::new(vec![])),
            contacts: Box::leak(Box::new(vec![])),
        }
    }

    fn bar(ends: [usize; 2], middle: Option<usize>) -> Axis {
        let mut anchors = vec![Anchor {
            source_node: 0,
            vertex: ends[0],
            t: 0.,
        }];
        if let Some(v) = middle {
            anchors.push(Anchor {
                source_node: 0,
                vertex: v,
                t: 0.5,
            });
        }
        anchors.push(Anchor {
            source_node: 0,
            vertex: ends[1],
            t: 1.,
        });
        Axis {
            source_axis: 0,
            endpoints: ends,
            anchors,
            spans: vec![],
        }
    }

    #[test]
    fn duplicated_bar_ends_merge_and_carry_their_bars() {
        for place in Placement::all() {
            let mut m = build(&place, &[slab(0., 2.)]);
            duplicate_panel(&mut m, &place);
            // A beam on each panel ends at a duplicated node of line x = 2;
            // the second beam has a mid anchor on its panel.
            let at = |m: &Model, p: [f64; 3]| {
                (0..m.vertices.len())
                    .filter(|&v| {
                        DVec3::from_array(m.vertices[v]).distance(DVec3::from_array(place.point(p)))
                            < 1e-9 * place.scale.max(1.)
                    })
                    .collect::<Vec<_>>()
            };
            let (a, b) = {
                let d = at(&m, [2., 4., 0.]);
                (d[0], d[1])
            };
            let far0 = at(&m, [0., 4., 0.])[0];
            let far1 = at(&m, [4., 4., 0.])[0];
            let mid = m.add_vertex(place.point([3., 4., 0.])).unwrap();
            let mut axes = vec![bar([far0, a], None), bar([b, far1], Some(mid))];
            let mut contacts = vec![Contact::Point {
                axis: 1,
                surface: 1,
                vertex: b,
                t: 0.,
                location: super::super::bars::Location::Boundary,
            }];
            let r = merge_coincident(
                &mut m,
                &mut Bars {
                    axes: &mut axes,
                    contacts: &mut contacts,
                },
                0.001 * place.scale,
                &BTreeSet::new(),
                &[],
            );
            assert_eq!(r.merged.len(), 2, "{:?}", r.rejected);
            let merged = r.merged.iter().find(|x| [a, b].contains(&x.kept)).unwrap();
            // Both bars now end at one vertex; the contact follows it.
            assert_eq!(axes[0].endpoints[1], merged.kept);
            assert_eq!(axes[1].endpoints[0], merged.kept);
            assert!(matches!(contacts[0], Contact::Point { vertex, .. } if vertex == merged.kept));
            // Every anchor lies on its straight bar.
            for axis in &axes {
                let [p, q] = axis.endpoints.map(|e| DVec3::from_array(m.vertices[e]));
                for anchor in &axis.anchors {
                    let x = DVec3::from_array(m.vertices[anchor.vertex]);
                    assert!(x.distance(p.lerp(q, anchor.t)) <= m.precision);
                }
            }
        }
    }

    #[test]
    fn beam_split_by_a_small_gap_is_joined_but_a_short_bar_is_kept() {
        for place in Placement::all() {
            let mut m = build(&place, &[slab(0., 4.)]);
            let v = |m: &mut Model, p: [f64; 3]| m.add_vertex(place.point(p)).unwrap();
            // Two beam pieces above the slab with a 20 mm gap.
            let (a0, a1) = (v(&mut m, [0., 1., 1.]), v(&mut m, [1.99, 1., 1.]));
            let (b0, b1) = (v(&mut m, [2.01, 1., 1.]), v(&mut m, [4., 1., 1.]));
            let mut axes = vec![bar([a0, a1], None), bar([b0, b1], None)];
            let mut contacts = vec![];
            let r = merge_bar_ends(
                &mut m,
                &mut Bars {
                    axes: &mut axes,
                    contacts: &mut contacts,
                },
                0.05 * place.scale,
                &BTreeSet::new(),
                &[],
            );
            assert_eq!(r.merged.len(), 1, "{:?}", r.rejected);
            assert_eq!(axes[0].endpoints[1], axes[1].endpoints[0]);
            // The same pieces joined by a 20 mm bar: nothing to merge.
            let mut m = build(&place, &[slab(0., 4.)]);
            let (a0, a1) = (v(&mut m, [0., 1., 1.]), v(&mut m, [1.99, 1., 1.]));
            let (b0, b1) = (v(&mut m, [2.01, 1., 1.]), v(&mut m, [4., 1., 1.]));
            let mut axes = vec![
                bar([a0, a1], None),
                bar([a1, b0], None),
                bar([b0, b1], None),
            ];
            let r = merge_bar_ends(
                &mut m,
                &mut Bars {
                    axes: &mut axes,
                    contacts: &mut contacts,
                },
                0.05 * place.scale,
                &BTreeSet::new(),
                &[],
            );
            assert!(r.merged.is_empty());
            assert_eq!(axes.len(), 3);
        }
    }

    #[test]
    fn beam_node_next_to_a_slab_edge_vertex_is_identified_and_contacts_refreshed() {
        use super::super::bars::{refresh_contacts, Location};
        for place in Placement::all() {
            // The slab edge x = 4 carries a vertex 0.4 mm off the beam axis.
            let panel = (
                vec![vec![
                    [0., 0., 0.],
                    [4., 0., 0.],
                    [4.0004, 2., 0.],
                    [4., 4., 0.],
                    [0., 4., 0.],
                ]],
                [0., 0., 1.],
            );
            let mut m = build(&place, &[panel]);
            let at = |m: &Model, p: [f64; 3]| {
                (0..m.vertices.len())
                    .find(|&v| {
                        DVec3::from_array(m.vertices[v]).distance(DVec3::from_array(place.point(p)))
                            < 1e-9 * place.scale.max(1.)
                    })
                    .unwrap()
            };
            let (a, b) = (at(&m, [4., 0., 0.]), at(&m, [4., 4., 0.]));
            let mid = m.add_vertex(place.point([4., 2., 0.])).unwrap();
            let before = m.vertices[mid];
            let mut axes = vec![bar([a, b], Some(mid))];
            let mut contacts = vec![];
            let r = merge_bar_anchors(
                &mut m,
                &mut Bars {
                    axes: &mut axes,
                    contacts: &mut contacts,
                },
                0.001 * place.scale,
                &BTreeSet::new(),
                &[],
            );
            assert_eq!(r.merged.len(), 1, "{:?}", r.rejected);
            assert_eq!(r.merged[0].kept, mid);
            // The bar stays straight; the slab contour now runs through it.
            assert_eq!(m.vertices[mid], before);
            assert!(m.surfaces[0].boundaries[0]
                .iter()
                .any(|e| m.edges[e.edge].contains(&mid)));
            refresh_contacts(&m, &axes, &mut contacts);
            let points = contacts
                .iter()
                .filter(|c| {
                    matches!(
                        c,
                        Contact::Point {
                            location: Location::Boundary,
                            ..
                        }
                    )
                })
                .count();
            assert_eq!(points, 3, "{contacts:?}");
            assert!(contacts.iter().any(|c| matches!(c,
                Contact::Interval { start_t, end_t, location: Location::Boundary, .. }
                    if *start_t < 1e-9 && *end_t > 1. - 1e-9)));
            // Idempotent.
            let again = merge_bar_anchors(
                &mut m,
                &mut Bars {
                    axes: &mut axes,
                    contacts: &mut contacts,
                },
                0.001 * place.scale,
                &BTreeSet::new(),
                &[],
            );
            assert!(again.merged.is_empty() && again.rejected.is_empty());
        }
    }

    #[test]
    fn wall_ends_millimetres_apart_on_a_slab_collapse_but_a_step_stays() {
        for place in Placement::all() {
            for (upper_end, merged) in [(2.997, 1), (2.98, 0)] {
                // Lower and upper walls in one plane end 3 mm (or 20 mm)
                // apart on the slab line.
                let mut m = build(
                    &place,
                    &[
                        slab(0., 4.),
                        wall(1., 3., -2., 0.),
                        wall(1., upper_end, 0., 2.),
                    ],
                );
                let j = run(&mut m);
                assert!(j.issues.is_empty(), "{:?}", j.issues);
                let tolerance = 0.01 * place.scale;
                let r =
                    collapse_short_edges(&mut m, &mut no_bars(), tolerance, &BTreeSet::new(), &[]);
                assert_eq!(r.merged.len(), merged, "{:?}", r.rejected);
                let shortest = (0..m.surfaces.len())
                    .flat_map(|s| m.surface_edges(s).collect::<Vec<_>>())
                    .map(|e| {
                        let [a, b] = m.edges[e];
                        DVec3::from_array(m.vertices[a]).distance(DVec3::from_array(m.vertices[b]))
                    })
                    .fold(f64::INFINITY, f64::min);
                if merged == 1 {
                    assert!(shortest >= tolerance);
                    // Every surface stays planar and valid; the lower wall
                    // corner is unchanged (it has more surfaces).
                    let corner = place.point([3., 2., 0.]);
                    assert!(m.vertices.iter().any(|&v| DVec3::from_array(v)
                        .distance(DVec3::from_array(corner))
                        < 1e-9 * place.scale.max(1.)));
                } else {
                    assert!((shortest - 0.02 * place.scale).abs() < 1e-6 * place.scale);
                }
            }
        }
    }

    #[test]
    fn column_through_a_slab_without_a_shared_node_gets_a_point_contact() {
        use super::super::bars::{refresh_contacts, Location};
        for place in Placement::all() {
            let m = {
                let mut m = build(&place, &[slab(0., 4.)]);
                let lo = m.add_vertex(place.point([2., 1., -3.])).unwrap();
                let mid = m.add_vertex(place.point([2., 1., 0.])).unwrap();
                let hi = m.add_vertex(place.point([2., 1., 3.])).unwrap();
                (m, [lo, mid, hi])
            };
            let (m, [lo, mid, hi]) = m;
            let axes = vec![bar([lo, hi], Some(mid))];
            let mut contacts = vec![];
            refresh_contacts(&m, &axes, &mut contacts);
            assert_eq!(contacts.len(), 1, "{contacts:?}");
            assert!(matches!(contacts[0],
                Contact::Point { vertex, surface: 0, location: Location::Interior, .. }
                    if vertex == mid));
        }
    }

    #[test]
    fn stale_contacts_gain_the_in_plane_interval_of_a_beam() {
        use super::super::bars::{refresh_contacts, Location};
        for place in Placement::all() {
            let mut m = build(&place, &[slab(0., 4.)]);
            let corner = (0..m.vertices.len())
                .find(|&v| {
                    DVec3::from_array(m.vertices[v])
                        .distance(DVec3::from_array(place.point([0., 0., 0.])))
                        < 1e-9 * place.scale.max(1.)
                })
                .unwrap();
            let inner = m.add_vertex(place.point([2., 2., 0.])).unwrap();
            let axes = vec![bar([corner, inner], None)];
            // Only the interior end was recorded before later vertex moves.
            let mut contacts = vec![Contact::Point {
                axis: 0,
                surface: 0,
                vertex: inner,
                t: 1.,
                location: Location::Interior,
            }];
            refresh_contacts(&m, &axes, &mut contacts);
            assert!(
                contacts.iter().any(|c| matches!(c,
                Contact::Interval { location: Location::Interior, start_t, end_t, .. }
                    if *start_t < 1e-9 && *end_t > 1. - 1e-9)),
                "{contacts:?}"
            );
            assert!(contacts.iter().any(|c| matches!(c,
                Contact::Point { vertex, location: Location::Boundary, .. } if *vertex == corner)));
        }
    }

    #[test]
    fn bar_passing_through_a_duplicate_is_not_bent() {
        let place = &Placement::all()[2];
        let mut m = build(place, &[slab(0., 2.)]);
        duplicate_panel(&mut m, place);
        let d: Vec<usize> = (0..m.vertices.len())
            .filter(|&v| {
                DVec3::from_array(m.vertices[v])
                    .distance(DVec3::from_array(place.point([2., 4., 0.])))
                    < 1e-9
            })
            .collect();
        let far0 = m.add_vertex(place.point([0., 5., 0.])).unwrap();
        let far1 = m.add_vertex(place.point([4., 3., 0.])).unwrap();
        // Both duplicates are interior anchors of bars.
        let mut axes = vec![bar([far0, far1], Some(d[0])), bar([far1, far0], Some(d[1]))];
        let mut contacts = vec![];
        let r = merge_coincident(
            &mut m,
            &mut Bars {
                axes: &mut axes,
                contacts: &mut contacts,
            },
            0.001,
            &BTreeSet::new(),
            &[],
        );
        // Exactly coincident interior anchors may merge (nothing moves);
        // any merge must keep every anchor on its straight bar.
        assert!(!r.merged.is_empty() || !r.rejected.is_empty());
        for axis in &axes {
            let [p, q] = axis.endpoints.map(|e| DVec3::from_array(m.vertices[e]));
            for anchor in &axis.anchors {
                let x = DVec3::from_array(m.vertices[anchor.vertex]);
                assert!(x.distance(p.lerp(q, anchor.t)) <= m.precision);
            }
        }
        // An interior anchor off the other bar's position would bend it.
        let mut m = build(place, &[slab(0., 2.)]);
        duplicate_panel(&mut m, place);
        let d: Vec<usize> = (0..m.vertices.len())
            .filter(|&v| {
                DVec3::from_array(m.vertices[v])
                    .distance(DVec3::from_array(place.point([2., 4., 0.])))
                    < 1e-9
            })
            .collect();
        let far0 = m.add_vertex(place.point([0., 5., 0.])).unwrap();
        let far1 = m.add_vertex(place.point([4., 3., 0.])).unwrap();
        m.move_vertex(d[1], place.point([2., 4.0005, 0.])).unwrap();
        let mut axes = vec![bar([far0, far1], Some(d[0])), bar([far1, far0], Some(d[1]))];
        let r = merge_coincident(
            &mut m,
            &mut Bars {
                axes: &mut axes,
                contacts: &mut contacts,
            },
            0.001,
            &BTreeSet::new(),
            &[],
        );
        assert!(r
            .rejected
            .iter()
            .any(|x| x.reason == "interior_bar_anchor" || x.reason == "bar_would_collapse"));
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
            let r = merge_wall_ends(
                &mut m,
                &mut no_bars(),
                0.05 * place.scale,
                &BTreeSet::new(),
                &[],
            );
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
            assert!(
                merge_wall_ends(&mut m, &mut no_bars(), 0., &BTreeSet::new(), &[])
                    .merged
                    .is_empty()
            );
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
