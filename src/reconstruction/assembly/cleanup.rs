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
use crate::reconstruction::{Model, PlaneFrame};
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

pub(super) fn users(model: &Model, v: usize) -> Vec<usize> {
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
pub(super) fn move_with_axes(
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
pub(super) fn merge(
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

#[derive(Debug, Clone, Serialize)]
pub struct Straightened {
    pub vertex: usize,
    pub source_node: Option<u32>,
    /// Surface whose plane the edge runs along.
    pub along_surface: usize,
    pub distance: f64,
}

/// A contour edge running along another surface's plane from one of its
/// vertices, with the far end off that plane by less than `limit` (a slab
/// edge continuing a wall base line, 4 um off after the frame solve) is
/// straightened: the far end moves onto that plane, keeping its own planes.
/// Otherwise the junction line and the contour edge are two lines
/// micrometres apart. Only near-collinear edges qualify (offset over length
/// below `angle`); fixed vertices never move.
pub(super) fn straighten_edges(
    model: &mut Model,
    axes: &[Axis],
    limit: f64,
    angle: f64,
    fixed: &BTreeSet<usize>,
    source_nodes: &[u32],
) -> Vec<Straightened> {
    let mut owners = std::collections::BTreeMap::<usize, BTreeSet<usize>>::new();
    for s in 0..model.surfaces.len() {
        for e in model.surface_edges(s).collect::<Vec<_>>() {
            for v in model.edges[e] {
                owners.entry(v).or_default().insert(s);
            }
        }
    }
    let mut out = vec![];
    let contour: Vec<(usize, usize, usize)> = (0..model.surfaces.len())
        .flat_map(|s| {
            model.surfaces[s]
                .boundaries
                .iter()
                .flatten()
                .map(move |u| (s, u.edge))
                .collect::<Vec<_>>()
        })
        .flat_map(|(s, e)| {
            let [a, b] = model.edges[e];
            [(s, a, b), (s, b, a)]
        })
        .collect();
    for (s, a, b) in contour {
        if fixed.contains(&b) {
            continue;
        }
        let (pa, pb) = (
            DVec3::from_array(model.vertices[a]),
            DVec3::from_array(model.vertices[b]),
        );
        let length = pa.distance(pb);
        let own_b = owners.get(&b).cloned().unwrap_or_default();
        let normal_s = DVec3::from_array(model.planes[model.surfaces[s].plane].normal);
        for &w in owners.get(&a).into_iter().flatten() {
            if own_b.contains(&w) {
                continue;
            }
            let plane_w = model.planes[model.surfaces[w].plane].clone();
            if DVec3::from_array(plane_w.normal).cross(normal_s).length() < 1e-9 {
                continue;
            }
            let d = plane_w.distance(pb.to_array()).abs();
            if d <= model.precision || d >= limit || d > length * angle {
                continue;
            }
            let mut planes: Vec<PlaneFrame> = own_b
                .iter()
                .map(|&u| model.planes[model.surfaces[u].plane].clone())
                .collect();
            planes.push(plane_w);
            let refs: Vec<&PlaneFrame> = planes.iter().collect();
            let Some(q) = super::intersection(pb, &refs, model.precision) else {
                continue;
            };
            if q.distance(pb) >= limit || move_with_axes(model, axes, b, q, limit).is_err() {
                continue;
            }
            out.push(Straightened {
                vertex: b,
                source_node: source_nodes.get(b).copied(),
                along_surface: w,
                distance: q.distance(pb),
            });
            break;
        }
    }
    out
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
        // Either end may stay: try the other when one is rejected.
        let before = [model.vertices[drop], model.vertices[keep]];
        let first = merge(model, bars, drop, keep, tolerance, fixed);
        let (drop, keep, kept_from, result) = match first {
            Ok(m) => (drop, keep, before[1], Ok(m)),
            Err(reason) => match merge(model, bars, keep, drop, tolerance, fixed) {
                Ok(m) => (keep, drop, before[0], Ok(m)),
                Err(_) => (drop, keep, before[1], Err(reason)),
            },
        };
        match result {
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

#[derive(Debug, Serialize)]
pub struct CollapsedBar {
    /// Axis index before the collapse; `removed` if the whole axis vanished.
    pub axis: usize,
    pub source_axis: usize,
    pub removed: bool,
    pub dropped: usize,
    pub kept: usize,
    pub dropped_source_node: Option<u32>,
    pub kept_source_node: Option<u32>,
    pub length: f64,
    /// Source bar elements of the collapsed piece.
    pub elements: Vec<u32>,
}

/// A bar through exactly the vertices of another bar with the same
/// stiffness sequence (overlapping source bars): represented once.
#[derive(Debug, Serialize)]
pub struct DuplicateBar {
    pub source_axis: usize,
    pub kept_source_axis: usize,
    pub elements: Vec<u32>,
}

#[derive(Debug, Default, Serialize)]
pub struct BarCollapseReport {
    pub tolerance: f64,
    pub collapsed: Vec<CollapsedBar>,
    pub rejected: Vec<RejectedMerge>,
    pub duplicates: Vec<DuplicateBar>,
}

/// Remove every bar whose anchors are, in order along it, the vertices of
/// an earlier kept bar (lowest source axis) with the same stiffness sequence.
/// Contacts of a removed bar go with it; the kept bar carries its own.
pub fn remove_duplicate_bars(bars: &mut Bars<'_>) -> Vec<DuplicateBar> {
    let key = |axis: &Axis| {
        let mut anchors = axis.anchors.clone();
        anchors.sort_by(|x, y| x.t.total_cmp(&y.t));
        let mut vertices: Vec<usize> = anchors.iter().map(|a| a.vertex).collect();
        let mut spans = axis.spans.clone();
        spans.sort_by(|x, y| x.start_t.total_cmp(&y.start_t));
        let mut stiffness: Vec<u32> = spans.iter().map(|s| s.stiffness).collect();
        if vertices.first() > vertices.last() {
            vertices.reverse();
            stiffness.reverse();
        }
        (vertices, stiffness)
    };
    let mut order: Vec<usize> = (0..bars.axes.len()).collect();
    order.sort_by_key(|&k| bars.axes[k].source_axis);
    let mut kept = std::collections::BTreeMap::new();
    let mut removed = vec![false; bars.axes.len()];
    let mut duplicates = vec![];
    for k in order {
        let axis = &bars.axes[k];
        match kept.entry(key(axis)) {
            std::collections::btree_map::Entry::Vacant(e) => {
                e.insert(axis.source_axis);
            }
            std::collections::btree_map::Entry::Occupied(e) => {
                removed[k] = true;
                duplicates.push(DuplicateBar {
                    source_axis: axis.source_axis,
                    kept_source_axis: *e.get(),
                    elements: axis.spans.iter().map(|s| s.element).collect(),
                });
            }
        }
    }
    if duplicates.is_empty() {
        return duplicates;
    }
    let mut index = vec![None; bars.axes.len()];
    let mut next = 0;
    for (k, gone) in removed.iter().enumerate() {
        if !gone {
            index[k] = Some(next);
            next += 1;
        }
    }
    let mut k = 0;
    bars.axes.retain(|_| {
        k += 1;
        !removed[k - 1]
    });
    bars.contacts.retain_mut(|c| {
        let (Contact::Point { axis, .. } | Contact::Interval { axis, .. }) = c;
        match index[*axis] {
            Some(i) => {
                *axis = i;
                true
            }
            None => false,
        }
    });
    duplicates
}

/// Axis `i` without `drop`, whose piece to `keep` collapses; `None` when the
/// whole axis collapses. Spans of the piece are returned as removed.
fn shorten(axis: &Axis, drop: usize, keep: usize) -> (Option<Axis>, Vec<u32>) {
    let t = |v: usize| axis.anchors.iter().find(|a| a.vertex == v).map(|a| a.t);
    let (Some(td), Some(tk)) = (t(drop), t(keep)) else {
        return (Some(axis.clone()), vec![]);
    };
    let (lo, hi) = (td.min(tk), td.max(tk));
    let tol = 1e-9;
    let mut removed = vec![];
    let mut spans = vec![];
    for span in &axis.spans {
        if span.start_t >= lo - tol && span.end_t <= hi + tol {
            removed.push(span.element);
            continue;
        }
        let mut span = span.clone();
        if (span.start_t - td).abs() <= tol {
            span.start_t = tk;
        }
        if (span.end_t - td).abs() <= tol {
            span.end_t = tk;
        }
        spans.push(span);
    }
    let anchors: Vec<_> = axis
        .anchors
        .iter()
        .filter(|a| a.vertex != drop)
        .cloned()
        .collect();
    if anchors.len() < 2 || spans.is_empty() {
        return (None, removed);
    }
    let mut endpoints = axis.endpoints;
    for e in endpoints.iter_mut() {
        if *e == drop {
            *e = keep;
        }
    }
    // Parameters follow the new ends.
    let t0 = t(endpoints[0]).unwrap();
    let t1 = t(endpoints[1]).unwrap();
    let map = |x: f64| ((x - t0) / (t1 - t0)).clamp(0., 1.);
    let mut result = axis.clone();
    result.endpoints = endpoints;
    result.anchors = anchors
        .into_iter()
        .map(|mut a| {
            a.t = map(a.t);
            a
        })
        .collect();
    result.spans = spans
        .into_iter()
        .map(|mut s| {
            s.start_t = map(s.start_t);
            s.end_t = map(s.end_t);
            s
        })
        .collect();
    (Some(result), removed)
}

/// Collapse bar pieces shorter than `tolerance` (short source bars joining
/// beams, and close consecutive nodes of one axis): the node of fewer
/// connections merges into the other, the piece's source elements are
/// reported, and a bar whose whole length collapses disappears. Other bars
/// and surfaces keep their planes and are never bent.
pub fn collapse_short_bars(
    model: &mut Model,
    bars: &mut Bars<'_>,
    tolerance: f64,
    fixed: &BTreeSet<usize>,
    source_nodes: &[u32],
) -> BarCollapseReport {
    let mut report = BarCollapseReport {
        tolerance,
        ..Default::default()
    };
    let mut tried = BTreeSet::new();
    loop {
        let mut best: Option<(f64, usize, usize, usize)> = None;
        for (i, axis) in bars.axes.iter().enumerate() {
            let mut anchors = axis.anchors.clone();
            anchors.sort_by(|x, y| x.t.total_cmp(&y.t));
            for w in anchors.windows(2) {
                let (u, v) = (w[0].vertex, w[1].vertex);
                let d = DVec3::from_array(model.vertices[u])
                    .distance(DVec3::from_array(model.vertices[v]));
                let key = (u.min(v), u.max(v));
                if d < tolerance
                    && !tried.contains(&key)
                    && best.is_none_or(|b| (d, i, key.0) < (b.0, b.1, b.2.min(b.3)))
                {
                    best = Some((d, i, u, v));
                }
            }
        }
        let Some((d, i, u, v)) = best else {
            return report;
        };
        tried.insert((u.min(v), u.max(v)));
        let connections = |x: usize| {
            users(model, x).len()
                + bars
                    .axes
                    .iter()
                    .filter(|a| a.anchors.iter().any(|y| y.vertex == x))
                    .count()
        };
        let (drop, keep) =
            if (connections(u), std::cmp::Reverse(u)) > (connections(v), std::cmp::Reverse(v)) {
                (v, u)
            } else {
                (u, v)
            };
        // Every bar with this piece between consecutive anchors loses it
        // (duplicate pieces of coincident bars collapse together).
        let consecutive = |axis: &Axis| {
            let mut anchors = axis.anchors.clone();
            anchors.sort_by(|x, y| x.t.total_cmp(&y.t));
            anchors.windows(2).any(|w| {
                let pair = (w[0].vertex, w[1].vertex);
                pair == (drop, keep) || pair == (keep, drop)
            })
        };
        let affected: Vec<usize> = (0..bars.axes.len())
            .filter(|&k| k == i || consecutive(&bars.axes[k]))
            .collect();
        let mut collapsed = vec![];
        let mut axes = vec![];
        let mut index = vec![None; bars.axes.len()];
        for (k, axis) in bars.axes.iter().enumerate() {
            if !affected.contains(&k) {
                index[k] = Some(axes.len());
                axes.push(axis.clone());
                continue;
            }
            let (shortened, elements) = shorten(axis, drop, keep);
            collapsed.push((k, axis.source_axis, shortened.is_none(), elements));
            if let Some(axis) = shortened {
                index[k] = Some(axes.len());
                axes.push(axis);
            }
        }
        let mut contacts: Vec<Contact> = bars
            .contacts
            .iter()
            .filter_map(|c| {
                let mut c = c.clone();
                let (Contact::Point { axis, .. } | Contact::Interval { axis, .. }) = &mut c;
                *axis = index[*axis]?;
                Some(c)
            })
            .collect();
        let mut trial = model.clone();
        let result = merge(
            &mut trial,
            &mut Bars {
                axes: &mut axes,
                contacts: &mut contacts,
            },
            drop,
            keep,
            tolerance,
            fixed,
        );
        match result {
            Ok(_) => {
                *model = trial;
                *bars.axes = axes;
                *bars.contacts = contacts;
                for (axis, source_axis, removed, elements) in collapsed {
                    report.collapsed.push(CollapsedBar {
                        axis,
                        source_axis,
                        removed,
                        dropped: drop,
                        kept: keep,
                        dropped_source_node: source_nodes.get(drop).copied(),
                        kept_source_node: source_nodes.get(keep).copied(),
                        length: d,
                        elements,
                    });
                }
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
        // The closest pair of ends of two bars (not ends of one bar), first
        // in axis order among equals; ends are found through a grid.
        let cell = tolerance.max(model.precision);
        let key = |v: usize| model.vertices[v].map(|x| (x / cell).floor() as i64);
        let mut grid = std::collections::BTreeMap::<[i64; 3], Vec<(usize, usize, usize)>>::new();
        for (i, a) in bars.axes.iter().enumerate() {
            for (k, &u) in a.endpoints.iter().enumerate() {
                grid.entry(key(u)).or_default().push((i, k, u));
            }
        }
        let mut best: Option<(f64, [usize; 4], usize, usize)> = None;
        for (i, a) in bars.axes.iter().enumerate() {
            for (k, &u) in a.endpoints.iter().enumerate() {
                let c = key(u);
                for dx in -1..=1 {
                    for dy in -1..=1 {
                        for dz in -1..=1 {
                            for &(j, l, w) in grid
                                .get(&[c[0] + dx, c[1] + dy, c[2] + dz])
                                .into_iter()
                                .flatten()
                            {
                                let b = &bars.axes[j];
                                if j <= i
                                    || u == w
                                    || tried.contains(&(u.min(w), u.max(w)))
                                    || a.endpoints.contains(&w)
                                    || b.endpoints.contains(&u)
                                {
                                    continue;
                                }
                                let d = DVec3::from_array(model.vertices[u])
                                    .distance(DVec3::from_array(model.vertices[w]));
                                let order = [i, j, k, l];
                                if d <= tolerance && best.is_none_or(|x| (d, order) < (x.0, x.1)) {
                                    best = Some((d, order, u, w));
                                }
                            }
                        }
                    }
                }
            }
        }
        let best = best.map(|(d, _, u, w)| (d, u, w));
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

    fn spanned(ends: [usize; 2], nodes: &[(usize, f64)], first: u32) -> Axis {
        Axis {
            source_axis: 0,
            endpoints: ends,
            anchors: nodes
                .iter()
                .map(|&(vertex, t)| Anchor {
                    source_node: 0,
                    vertex,
                    t,
                })
                .collect(),
            spans: nodes
                .windows(2)
                .enumerate()
                .map(|(k, w)| crate::reconstruction::recognize::SourceSpan {
                    element: first + k as u32,
                    stiffness: 1,
                    start_t: w[0].1,
                    end_t: w[1].1,
                })
                .collect(),
        }
    }

    #[test]
    fn short_bars_collapse_and_report_their_elements() {
        for place in Placement::all() {
            let mut m = build(&place, &[slab(0., 4.)]);
            let v = |m: &mut Model, p: [f64; 3]| m.add_vertex(place.point(p)).unwrap();
            // Two beams joined by a 20 mm bar above the slab.
            let (a0, a1) = (v(&mut m, [0., 1., 1.]), v(&mut m, [1.99, 1., 1.]));
            let (b0, b1) = (v(&mut m, [2.01, 1., 1.]), v(&mut m, [4., 1., 1.]));
            // A long beam with two nodes 30 mm apart at x = 2.
            let (c0, c1, c2, c3) = (
                v(&mut m, [0., 3., 1.]),
                v(&mut m, [2., 3., 1.]),
                v(&mut m, [2.03, 3., 1.]),
                v(&mut m, [4., 3., 1.]),
            );
            let mut axes = vec![
                spanned([a0, a1], &[(a0, 0.), (a1, 1.)], 10),
                spanned([a1, b0], &[(a1, 0.), (b0, 1.)], 20),
                spanned([b0, b1], &[(b0, 0.), (b1, 1.)], 30),
                spanned([c0, c3], &[(c0, 0.), (c1, 0.5), (c2, 0.5075), (c3, 1.)], 40),
            ];
            let mut contacts = vec![];
            let r = collapse_short_bars(
                &mut m,
                &mut Bars {
                    axes: &mut axes,
                    contacts: &mut contacts,
                },
                0.05 * place.scale,
                &BTreeSet::new(),
                &[],
            );
            assert_eq!(r.collapsed.len(), 2, "{:?}", r.rejected);
            // The joining bar is gone, its element reported; the beams meet.
            assert_eq!(axes.len(), 3);
            let joining = r.collapsed.iter().find(|c| c.removed).unwrap();
            assert_eq!(joining.elements, vec![20]);
            assert_eq!(axes[0].endpoints[1], axes[1].endpoints[0]);
            // The long beam keeps three nodes, straight, with two elements.
            let long = &axes[2];
            assert_eq!(long.anchors.len(), 3);
            assert_eq!(long.spans.len(), 2);
            let inner = r.collapsed.iter().find(|c| !c.removed).unwrap();
            assert_eq!(inner.elements, vec![41]);
            let [p, q] = long.endpoints.map(|e| DVec3::from_array(m.vertices[e]));
            for anchor in &long.anchors {
                let x = DVec3::from_array(m.vertices[anchor.vertex]);
                assert!(x.distance(p.lerp(q, anchor.t)) <= m.precision);
            }
            let covered: f64 = long.spans.iter().map(|s| s.end_t - s.start_t).sum();
            assert!((covered - 1.).abs() < 1e-12);
            // Nothing shorter than the tolerance is left; a 60 mm piece stays.
            assert!(collapse_short_bars(
                &mut m,
                &mut Bars {
                    axes: &mut axes,
                    contacts: &mut contacts,
                },
                0.05 * place.scale,
                &BTreeSet::new(),
                &[],
            )
            .collapsed
            .is_empty());
        }
    }

    #[test]
    fn duplicate_short_pieces_between_the_same_nodes_collapse_together() {
        // Two beams from different far nodes both end with a 7 mm piece
        // between the same two nodes (overlapping source bars): the piece
        // collapses on both at once instead of being rejected as a collapse.
        for place in Placement::all() {
            let mut m = build(&place, &[slab(0., 4.)]);
            let v = |m: &mut Model, p: [f64; 3]| m.add_vertex(place.point(p)).unwrap();
            let (a, b) = (v(&mut m, [2., 1., 1.]), v(&mut m, [2.007, 1., 1.]));
            let (f0, f1) = (v(&mut m, [0., 1., 1.]), v(&mut m, [4., 1., 1.]));
            let mut axes = vec![
                spanned([f0, b], &[(f0, 0.), (a, 2. / 2.007), (b, 1.)], 10),
                spanned([f1, a], &[(f1, 0.), (b, 0.9965), (a, 1.)], 20),
            ];
            let mut contacts = vec![];
            let r = collapse_short_bars(
                &mut m,
                &mut Bars {
                    axes: &mut axes,
                    contacts: &mut contacts,
                },
                0.05 * place.scale,
                &BTreeSet::new(),
                &[],
            );
            assert!(r.rejected.is_empty(), "{:?}", r.rejected);
            assert_eq!(r.collapsed.len(), 2);
            assert_eq!(axes.len(), 2);
            assert!(axes
                .iter()
                .all(|a| a.anchors.len() == 2 && a.spans.len() == 1));
            assert_eq!(axes[0].endpoints[1], axes[1].endpoints[1]);
        }
    }

    #[test]
    fn bar_through_the_vertices_of_another_bar_is_kept_once() {
        let mut m = build(&Placement::all()[0], &[slab(0., 4.)]);
        let a = m.add_vertex([0., 1., 1.]).unwrap();
        let b = m.add_vertex([2., 1., 1.]).unwrap();
        let c = m.add_vertex([4., 1., 1.]).unwrap();
        let mut first = spanned([a, b], &[(a, 0.), (b, 1.)], 10);
        first.source_axis = 5;
        // Reversed, and with a source axis numbered lower.
        let mut second = spanned([b, a], &[(b, 0.), (a, 1.)], 20);
        second.source_axis = 3;
        let mut other = spanned([b, c], &[(b, 0.), (c, 1.)], 30);
        other.source_axis = 4;
        let mut stiffer = spanned([a, b], &[(a, 0.), (b, 1.)], 40);
        stiffer.source_axis = 6;
        stiffer.spans[0].stiffness = 2;
        let mut axes = vec![first, second, other, stiffer];
        let mut contacts = vec![];
        let duplicates = remove_duplicate_bars(&mut Bars {
            axes: &mut axes,
            contacts: &mut contacts,
        });
        assert_eq!(duplicates.len(), 1);
        assert_eq!(
            (duplicates[0].source_axis, duplicates[0].kept_source_axis),
            (5, 3)
        );
        assert_eq!(duplicates[0].elements, vec![10]);
        let left: Vec<usize> = axes.iter().map(|a| a.source_axis).collect();
        assert_eq!(left, vec![3, 4, 6]);
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
    fn a_stale_contact_of_a_node_moved_off_the_plane_is_dropped() {
        use super::super::bars::{refresh_contacts, Location};
        for place in Placement::all() {
            let mut m = build(&place, &[slab(0., 4.)]);
            let lo = m.add_vertex(place.point([2., 1., -3.])).unwrap();
            let mid = m.add_vertex(place.point([2., 1., 0.])).unwrap();
            let hi = m.add_vertex(place.point([2., 1., 3.])).unwrap();
            let axes = vec![bar([lo, hi], Some(mid))];
            let mut contacts = vec![Contact::Point {
                axis: 0,
                surface: 0,
                vertex: mid,
                t: 0.5,
                location: Location::Interior,
            }];
            // The node leaves the slab plane (1 mm): the record is stale.
            m.move_vertex(mid, place.point([2., 1., 0.001])).unwrap();
            refresh_contacts(&m, &axes, &mut contacts);
            assert!(contacts.is_empty(), "{contacts:?}");
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

    #[test]
    fn slab_edge_continuing_a_wall_base_is_straightened_onto_the_wall_plane() {
        use super::super::junctions::tests::{build, Placement};
        for place in Placement::all() {
            for (off, straightened) in [(4e-6, true), (0.002, false)] {
                // A wall in x = 0 over y 0..1 stands on a slab whose edge
                // continues the wall base to y = 3, its far end `off` away.
                let wall = (
                    vec![vec![[0., 0., 0.], [0., 1., 0.], [0., 1., 1.], [0., 0., 1.]]],
                    [1., 0., 0.],
                );
                let slab = (
                    vec![vec![
                        [0., 0., 0.],
                        [0., 1., 0.],
                        [off, 3., 0.],
                        [4., 3., 0.],
                        [4., 0., 0.],
                    ]],
                    [0., 0., 1.],
                );
                let mut m = build(&place, &[wall, slab]);
                let far = (0..m.vertices.len())
                    .find(|&v| {
                        DVec3::from_array(m.vertices[v])
                            .distance(DVec3::from_array(place.point([off, 3., 0.])))
                            < 1e-9 * place.scale
                    })
                    .unwrap();
                let r = straighten_edges(
                    &mut m,
                    &[],
                    0.001 * place.scale,
                    0.02,
                    &BTreeSet::new(),
                    &[],
                );
                assert_eq!(r.len(), usize::from(straightened), "off {off}: {r:?}");
                let wall_plane = &m.planes[m.surfaces[0].plane];
                let d = wall_plane.distance(m.vertices[far]).abs();
                if straightened {
                    assert!(d <= m.precision, "{d}");
                } else {
                    assert!(d > 0.0019 * place.scale);
                }
            }
        }
    }
}
