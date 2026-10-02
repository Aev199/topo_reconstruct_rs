//! Assemble two-endpoint axes against fixed surface geometry. Contacts are
//! geometric constraints for later meshing, never inferred mechanical ties.
use super::{frame, intersection, movement_budget, MeshData, Model, PlaneFrame, Policy};
use crate::reconstruction::recognize::SourceSpan;
use geo::{
    line_intersection::{line_intersection, LineIntersection},
    Contains, Line, LineString, Point, Polygon,
};
use glam::{DVec2, DVec3};
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};
mod repair;
pub use repair::{Change as BoundaryChange, Repair as BoundaryRepair};

#[derive(Debug, Clone, Serialize)]
pub struct Anchor {
    /// Source node of the anchor; `NO_SOURCE_NODE` for a generated crossing
    /// of the bar with a surface (reported in `Report::imprinted`).
    pub source_node: u32,
    pub vertex: usize,
    pub t: f64,
}
#[derive(Debug, Clone, Serialize)]
pub struct Axis {
    pub source_axis: usize,
    pub endpoints: [usize; 2],
    pub anchors: Vec<Anchor>,
    pub spans: Vec<SourceSpan>,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Location {
    Boundary,
    Interior,
}
#[derive(Debug, Clone, Serialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum Contact {
    Point {
        axis: usize,
        surface: usize,
        vertex: usize,
        t: f64,
        location: Location,
    },
    Interval {
        axis: usize,
        surface: usize,
        start_t: f64,
        end_t: f64,
        location: Location,
    },
}
#[derive(Debug, Serialize)]
pub struct Issue {
    pub source_axis: usize,
    pub source_elements: Vec<u32>,
    pub source_nodes: Vec<u32>,
    pub reason: String,
}
#[derive(Debug, Default, Serialize)]
pub struct Report {
    pub boundary_repairs: Vec<BoundaryRepair>,
    pub all_axes_built: bool,
    pub axes: Vec<Axis>,
    pub contacts: Vec<Contact>,
    pub issues: Vec<Issue>,
    pub rejected_shared_anchors: BTreeMap<u32, String>,
    pub maximum_additional_movement: f64,
    /// No intersection-driven edge subdivision or mesh generation at this stage.
    pub mesh_constraints_complete: bool,
    /// Surface vertices on a bar axis made anchors of that bar.
    pub imprinted: Vec<Imprint>,
}

/// Anchor source node of a generated bar-surface crossing.
pub const NO_SOURCE_NODE: u32 = u32::MAX;

#[derive(Debug, Clone, Serialize)]
pub struct Imprint {
    pub axis: usize,
    pub vertex: usize,
    /// `NO_SOURCE_NODE` for a crossing.
    pub source_node: u32,
    pub t: f64,
    /// `surface_vertex` (a surface vertex on the axis), `crossing` (the
    /// bar passes through the surface), `slid_anchor` (a bar node slid
    /// along its bar onto the crossing next to it), `bar_crossing` (a
    /// vertex shared with a crossing bar) or `shared_node` (a node of a
    /// crossing bar slid onto the crossing).
    pub kind: String,
    pub surface: Option<usize>,
    /// Movement of the vertex (onto the axis, or along it).
    pub movement: f64,
}

/// A bar passing through a surface (its axis crosses the surface plane
/// inside the material, away from its anchors) shares a generated vertex
/// with it: the vertex becomes an anchor of the bar and, when it falls on a
/// contour or junction edge of the surface, splits that edge. Crossings
/// within the minimum edge length of a surface vertex are left alone (near
/// touches, handled elsewhere). A crossing within it of an interior bar node
/// (a beam node a few micrometres off a wall foot line) moves that node
/// along its bar onto the crossing: the bar stays straight and shares the
/// node with the surface.
pub fn imprint_crossings(model: &mut Model, axes: &mut [Axis]) -> Vec<Imprint> {
    let precision = model.precision;
    let minimum = model.minimum_edge;
    let boxes: Vec<(DVec3, DVec3)> = (0..model.surfaces.len())
        .map(|s| {
            model
                .surface_edges(s)
                .flat_map(|e| model.edges[e])
                .map(|v| DVec3::from_array(model.vertices[v]))
                .fold(
                    (DVec3::splat(f64::INFINITY), DVec3::splat(f64::NEG_INFINITY)),
                    |(lo, hi), p| (lo.min(p), hi.max(p)),
                )
        })
        .collect();
    let mut out = vec![];
    for i in 0..axes.len() {
        let [ea, eb] = axes[i].endpoints;
        let a = DVec3::from_array(model.vertices[ea]);
        let b = DVec3::from_array(model.vertices[eb]);
        let d = b - a;
        let length = d.length();
        if length <= minimum {
            continue;
        }
        let (lo, hi) = (a.min(b), a.max(b));
        for s in 0..model.surfaces.len() {
            let (slo, shi) = boxes[s];
            if (hi + precision).cmplt(slo).any() || (lo - precision).cmpgt(shi).any() {
                continue;
            }
            let plane = model.planes[model.surfaces[s].plane].clone();
            let (da, db) = (plane.distance(a.to_array()), plane.distance(b.to_array()));
            if da.abs() <= precision || db.abs() <= precision || da.signum() == db.signum() {
                continue;
            }
            let t = da / (da - db);
            let p = a + d * t;
            let p = p - DVec3::from_array(plane.normal) * plane.distance(p.to_array());
            let uv = plane.project(p.to_array());
            let Some(_) = location(uv, &model.surfaces[s].contours, precision) else {
                continue;
            };
            let near_vertex = model
                .surface_edges(s)
                .flat_map(|e| model.edges[e])
                .any(|v| DVec3::from_array(model.vertices[v]).distance(p) < minimum);
            if near_vertex {
                continue;
            }
            let near: Vec<usize> = (0..axes[i].anchors.len())
                .filter(|&k| (axes[i].anchors[k].t - t).abs() * length < minimum)
                .collect();
            let mut trial = model.clone();
            let (v, slid, movement) = match near[..] {
                [] => {
                    let Ok(v) = trial.add_vertex(p.to_array()) else {
                        continue;
                    };
                    (v, None, 0.)
                }
                [k] => {
                    let v = axes[i].anchors[k].vertex;
                    let movement = DVec3::from_array(model.vertices[v]).distance(p);
                    if !slides_onto(&mut trial, axes, i, v, p) {
                        continue;
                    }
                    (v, Some(k), movement)
                }
                _ => continue,
            };
            // On a contour or junction edge of the surface: split it.
            let on_edge: Vec<usize> = trial
                .surface_edges(s)
                .filter(|&e| {
                    let [x, y] = trial.edges[e];
                    let (px, py) = (
                        DVec3::from_array(trial.vertices[x]),
                        DVec3::from_array(trial.vertices[y]),
                    );
                    let q = py - px;
                    let u = (p - px).dot(q) / q.length_squared();
                    u > 0. && u < 1. && (px + q * u).distance(p) <= precision
                })
                .collect();
            if on_edge
                .iter()
                .any(|&e| trial.split_edge_within(e, v, precision).is_err())
            {
                continue;
            }
            *model = trial;
            let (source_node, kind) = if let Some(k) = slid {
                // Element spans bounded by the node follow it.
                let old = axes[i].anchors[k].t;
                for span in &mut axes[i].spans {
                    for u in [&mut span.start_t, &mut span.end_t] {
                        if *u == old {
                            *u = t;
                        }
                    }
                }
                axes[i].anchors[k].t = t;
                (axes[i].anchors[k].source_node, "slid_anchor")
            } else {
                axes[i].anchors.push(Anchor {
                    source_node: NO_SOURCE_NODE,
                    vertex: v,
                    t,
                });
                (NO_SOURCE_NODE, "crossing")
            };
            axes[i].anchors.sort_by(|x, y| x.t.total_cmp(&y.t));
            out.push(Imprint {
                axis: i,
                vertex: v,
                source_node,
                t,
                kind: kind.into(),
                surface: Some(s),
                movement,
            });
        }
    }
    out
}

/// Two bars crossing at a point away from their nodes (beams lying in one
/// slab, one passing through the span of another) share a generated vertex
/// there: it becomes an anchor of both (and of a third bar through the
/// same point). A node of one bar within the minimum edge of the crossing
/// slides onto it when it may (see `slides_with_ends`) and is shared; other
/// crossings next to source nodes are left alone (a crossing next to
/// another generated crossing is not), as are bars already sharing a node
/// and parallel bars.
pub fn imprint_bar_crossings(model: &mut Model, axes: &mut [Axis]) -> Vec<Imprint> {
    let precision = model.precision;
    let minimum = model.minimum_edge;
    // Endpoints never move here.
    let points: Vec<[DVec3; 2]> = axes
        .iter()
        .map(|a| a.endpoints.map(|v| DVec3::from_array(model.vertices[v])))
        .collect();
    let ends = |_: &[Axis], i: usize| points[i];
    let lengths: Vec<f64> = points.iter().map(|[a, b]| a.distance(*b)).collect();
    let cell = lengths.iter().sum::<f64>() / lengths.len().max(1) as f64;
    if !(cell > precision) {
        return vec![];
    }
    let key = |p: DVec3| ((p / cell).floor()).as_ivec3().to_array();
    let mut grid = BTreeMap::<[i32; 3], Vec<usize>>::new();
    for i in 0..axes.len() {
        let [a, b] = ends(axes, i);
        let (lo, hi) = (key(a.min(b) - precision), key(a.max(b) + precision));
        for x in lo[0]..=hi[0] {
            for y in lo[1]..=hi[1] {
                for z in lo[2]..=hi[2] {
                    grid.entry([x, y, z]).or_default().push(i);
                }
            }
        }
    }
    let mut pairs = BTreeSet::new();
    for members in grid.values() {
        for (k, &i) in members.iter().enumerate() {
            for &j in &members[k + 1..] {
                pairs.insert((i.min(j), i.max(j)));
            }
        }
    }
    let mut out = vec![];
    // Generated vertices by cell (a point near a cell border is found from
    // the neighbouring pair's own cell only if both round alike; the
    // precision is far below the cell size).
    let mut created = BTreeMap::<[i32; 3], Vec<usize>>::new();
    for (i, j) in pairs {
        let shared = axes[i]
            .anchors
            .iter()
            .any(|x| axes[j].anchors.iter().any(|y| y.vertex == x.vertex));
        if shared {
            continue;
        }
        let ([a0, a1], [b0, b1]) = (ends(axes, i), ends(axes, j));
        let (d1, d2, r) = (a1 - a0, b1 - b0, a0 - b0);
        let (a, b, e) = (d1.length_squared(), d1.dot(d2), d2.length_squared());
        let denominator = a * e - b * b;
        if denominator <= 1e-12 * a * e {
            continue;
        }
        let (c, f) = (d1.dot(r), d2.dot(r));
        let (s, t) = ((b * f - c * e) / denominator, (a * f - b * c) / denominator);
        if !(s > 0. && s < 1. && t > 0. && t < 1.) {
            continue;
        }
        let (pa, pb) = (a0 + d1 * s, b0 + d2 * t);
        if pa.distance(pb) > precision {
            continue;
        }
        let p = (pa + pb) * 0.5;
        // A node of a bar next to the crossing: a vertex generated for a
        // crossing with a third bar at the same point is shared; a source
        // node there is a near touch. A crossing with a third bar elsewhere
        // (bars fanning out from a node, crossed next to it) does not
        // block this one: unshared crossings fail the mesh, a short piece
        // is only a review item.
        let near = |k: usize, u: f64| {
            axes[k]
                .anchors
                .iter()
                .filter(|x| (x.t - u).abs() * lengths[k] < minimum)
                .find(|x| {
                    x.source_node != NO_SOURCE_NODE
                        || DVec3::from_array(model.vertices[x.vertex]).distance(p) <= precision
                })
                .map(|x| (x.source_node == NO_SOURCE_NODE).then_some(x.vertex))
        };
        let (v, bars) = match (near(i, s), near(j, t)) {
            (None, None) => {
                // Another pair through the same point already generated it.
                let existing = created
                    .get(&key(p))
                    .into_iter()
                    .flatten()
                    .copied()
                    .find(|&w| DVec3::from_array(model.vertices[w]).distance(p) <= precision);
                let v = match existing {
                    Some(w) => w,
                    None => {
                        let Ok(v) = model.add_vertex(p.to_array()) else {
                            continue;
                        };
                        created.entry(key(p)).or_default().push(v);
                        v
                    }
                };
                (v, vec![(i, s), (j, t)])
            }
            (Some(Some(v)), None) => (v, vec![(j, t)]),
            (None, Some(Some(v))) => (v, vec![(i, s)]),
            // A node of one bar next to the crossing (a slab corner where a
            // beam leaves the slab, micrometres off a crossing beam) slides
            // along its bar onto the crossing, keeping its planes, and is
            // shared: the bars meet at one vertex instead of a near touch.
            (Some(None), None) | (None, Some(None)) => {
                let (k, u, other, w) = if near(i, s).is_some() {
                    (i, s, j, t)
                } else {
                    (j, t, i, s)
                };
                let Some(n) = axes[k]
                    .anchors
                    .iter()
                    .position(|x| (x.t - u).abs() * lengths[k] < minimum)
                else {
                    continue;
                };
                let v = axes[k].anchors[n].vertex;
                let movement = DVec3::from_array(model.vertices[v]).distance(p);
                if movement > minimum || !slides_with_ends(model, axes, k, v, p, minimum) {
                    continue;
                }
                // Element spans bounded by the node follow it.
                let old = axes[k].anchors[n].t;
                for span in &mut axes[k].spans {
                    for x in [&mut span.start_t, &mut span.end_t] {
                        if *x == old {
                            *x = u;
                        }
                    }
                }
                axes[k].anchors[n].t = u;
                let source_node = axes[k].anchors[n].source_node;
                axes[other].anchors.push(Anchor {
                    source_node,
                    vertex: v,
                    t: w,
                });
                axes[other].anchors.sort_by(|x, y| x.t.total_cmp(&y.t));
                for (axis, t, kind, movement) in [
                    (k, u, "slid_anchor", movement),
                    (other, w, "shared_node", 0.),
                ] {
                    out.push(Imprint {
                        axis,
                        vertex: v,
                        source_node,
                        t,
                        kind: kind.into(),
                        surface: None,
                        movement,
                    });
                }
                continue;
            }
            _ => continue,
        };
        for (k, u) in bars {
            axes[k].anchors.push(Anchor {
                source_node: NO_SOURCE_NODE,
                vertex: v,
                t: u,
            });
            axes[k].anchors.sort_by(|x, y| x.t.total_cmp(&y.t));
            out.push(Imprint {
                axis: k,
                vertex: v,
                source_node: NO_SOURCE_NODE,
                t: u,
                kind: "bar_crossing".into(),
                surface: None,
                movement: 0.,
            });
        }
    }
    out
}

/// Move the interior bar node `v` of axis `i` along the axis to `p`: only a
/// node of no other axis and no surface contour, and only when every plane
/// whose material holds it (the slab a beam lies in) holds `p` too.
fn slides_onto(model: &mut Model, axes: &[Axis], i: usize, v: usize, p: DVec3) -> bool {
    if axes[i].endpoints.contains(&v)
        || axes.iter().enumerate().any(|(j, a)| {
            j != i && (a.endpoints.contains(&v) || a.anchors.iter().any(|x| x.vertex == v))
        })
        || (0..model.surfaces.len())
            .any(|s| model.surface_edges(s).any(|e| model.edges[e].contains(&v)))
    {
        return false;
    }
    keeps_planes(model, v, p) && model.move_vertex(v, p.to_array()).is_ok()
}

/// Every plane whose material holds vertex `v` holds `p` too.
fn keeps_planes(model: &Model, v: usize, p: DVec3) -> bool {
    let precision = model.precision;
    let q = DVec3::from_array(model.vertices[v]);
    model.surfaces.iter().all(|surface| {
        let plane = &model.planes[surface.plane];
        plane.distance(q.to_array()).abs() > precision
            || location(plane.project(q.to_array()), &surface.contours, precision).is_none()
            || plane.distance(p.to_array()).abs() <= precision
    })
}

/// Move the interior node `v` of axis `k` along it to `p` (within `limit`),
/// keeping its planes: a contour vertex may move, and bars ending at `v`
/// follow it as ends (see `move_with_axes`); `v` must be an interior node
/// of no other bar. Nothing changes on failure.
fn slides_with_ends(
    model: &mut Model,
    axes: &[Axis],
    k: usize,
    v: usize,
    p: DVec3,
    limit: f64,
) -> bool {
    let interior = |a: &Axis| !a.endpoints.contains(&v) && a.anchors.iter().any(|x| x.vertex == v);
    if !interior(&axes[k])
        || axes.iter().enumerate().any(|(j, a)| j != k && interior(a))
        || !keeps_planes(model, v, p)
    {
        return false;
    }
    let mut others = axes.to_vec();
    others[k].anchors.retain(|x| x.vertex != v);
    let mut trial = model.clone();
    if super::cleanup::move_with_axes(&mut trial, &others, v, p, limit).is_err() {
        return false;
    }
    *model = trial;
    true
}

/// A surface vertex lying on a bar axis inside its span (a bar running
/// along a slab edge from a slab corner, a bar passing a wall corner) is a
/// shared node: it becomes an anchor of the bar, so the bar and the surface
/// meet at one vertex instead of an unshared contact. Only vertices with a
/// source node are imprinted (their identity is the source node).
pub fn imprint_surface_vertices(
    model: &mut Model,
    axes: &mut [Axis],
    source_nodes: &[u32],
) -> Vec<Imprint> {
    let precision = model.precision;
    // A vertex nearly on the axis (within the minimum edge, e.g. a slab
    // edge 1 um off a beam along it) is moved onto the axis when it stays
    // on all its planes; bar nodes themselves never move here.
    let snap = model.minimum_edge;
    let used: BTreeSet<usize> = (0..model.surfaces.len())
        .flat_map(|s| model.surface_edges(s).collect::<Vec<_>>())
        .flat_map(|e| model.edges[e])
        .filter(|&v| v < source_nodes.len())
        .collect();
    let bar_nodes: BTreeSet<usize> = axes
        .iter()
        .flat_map(|a| a.anchors.iter().map(|x| x.vertex))
        .collect();
    let cell = 1.0_f64.max(precision);
    let key = |p: DVec3| {
        (
            (p.x / cell).floor() as i64,
            (p.y / cell).floor() as i64,
            (p.z / cell).floor() as i64,
        )
    };
    let mut grid = BTreeMap::<(i64, i64, i64), Vec<usize>>::new();
    for &v in &used {
        grid.entry(key(DVec3::from_array(model.vertices[v])))
            .or_default()
            .push(v);
    }
    let mut out = vec![];
    for i in 0..axes.len() {
        let a = DVec3::from_array(model.vertices[axes[i].endpoints[0]]);
        let b = DVec3::from_array(model.vertices[axes[i].endpoints[1]]);
        let d = b - a;
        let length = d.length();
        if length <= precision {
            continue;
        }
        let (lo, hi) = (key(a.min(b) - snap), key(a.max(b) + snap));
        let mut found = vec![];
        for x in lo.0..=hi.0 {
            for y in lo.1..=hi.1 {
                for z in lo.2..=hi.2 {
                    for &v in grid.get(&(x, y, z)).into_iter().flatten() {
                        if axes[i].anchors.iter().any(|x| x.vertex == v) {
                            continue;
                        }
                        let p = DVec3::from_array(model.vertices[v]);
                        let t = (p - a).dot(d) / (length * length);
                        let off = (a + d * t).distance(p);
                        if t * length > snap
                            && (1. - t) * length > snap
                            && (off <= precision || (off <= snap && !bar_nodes.contains(&v)))
                        {
                            found.push((t, v));
                        }
                    }
                }
            }
        }
        found.sort_by(|x, y| x.0.total_cmp(&y.0).then(x.1.cmp(&y.1)));
        found.dedup_by(|x, y| x.1 == y.1);
        for (t, v) in found {
            // Not next to an existing anchor (that would be a coincident node).
            if axes[i]
                .anchors
                .iter()
                .any(|x| (x.t - t).abs() * length <= snap)
            {
                continue;
            }
            let q = a + d * t;
            let movement = DVec3::from_array(model.vertices[v]).distance(q);
            if movement > precision && model.move_vertex(v, q.to_array()).is_err() {
                continue;
            }
            axes[i].anchors.push(Anchor {
                source_node: source_nodes[v],
                vertex: v,
                t,
            });
            out.push(Imprint {
                axis: i,
                vertex: v,
                source_node: source_nodes[v],
                t,
                kind: "surface_vertex".into(),
                surface: None,
                movement,
            });
        }
        axes[i].anchors.sort_by(|x, y| x.t.total_cmp(&y.t));
    }
    out
}

fn polygon(contours: &[Vec<[f64; 2]>]) -> Polygon<f64> {
    let ring = |r: &Vec<[f64; 2]>| {
        LineString::from(
            r.iter()
                .chain(r.first())
                .map(|p| (p[0], p[1]))
                .collect::<Vec<_>>(),
        )
    };
    Polygon::new(
        ring(&contours[0]),
        contours.iter().skip(1).map(ring).collect(),
    )
}
fn location(uv: [f64; 2], contours: &[Vec<[f64; 2]>], precision: f64) -> Option<Location> {
    let p = DVec2::from_array(uv);
    for ring in contours {
        for i in 0..ring.len() {
            let a = DVec2::from_array(ring[i]);
            let d = DVec2::from_array(ring[(i + 1) % ring.len()]) - a;
            let q = a + d * ((p - a).dot(d) / d.length_squared()).clamp(0., 1.);
            if p.distance(q) <= precision {
                return Some(Location::Boundary);
            }
        }
    }
    polygon(contours)
        .contains(&Point::new(p.x, p.y))
        .then_some(Location::Interior)
}

/// Largest turn of a bar whose two ends moved within `budgets` (sum of
/// both), never less than the recognition `angle`.
fn allowed_turn(angle: f64, budgets: f64, length: f64) -> f64 {
    angle.max((budgets / length).min(1.).asin())
}

/// Clip a coplanar axis to the material, including holes and nonconvex outlines.
fn intervals(
    a: [f64; 2],
    b: [f64; 2],
    contours: &[Vec<[f64; 2]>],
    precision: f64,
) -> Vec<(f64, f64, Location)> {
    let a = DVec2::from_array(a);
    let b = DVec2::from_array(b);
    let d = b - a;
    let length = d.length();
    if length <= precision {
        return vec![];
    }
    let line = Line::new((a.x, a.y), (b.x, b.y));
    let parameter =
        |x: geo::Coord<f64>| ((DVec2::new(x.x, x.y) - a).dot(d) / d.length_squared()).clamp(0., 1.);
    let mut cuts = vec![0., 1.];
    for ring in contours {
        for i in 0..ring.len() {
            let p = ring[i];
            // A contour vertex on the axis is a cut even when the adjacent
            // edges are only nearly collinear with it (no exact collinear
            // overlap is reported): the axis may leave the boundary there.
            let v = DVec2::from_array(p);
            let t = (v - a).dot(d) / d.length_squared();
            if t > 0. && t < 1. && (a + d * t).distance(v) <= precision {
                cuts.push(t);
            }
            let q = ring[(i + 1) % ring.len()];
            match line_intersection(line, Line::new((p[0], p[1]), (q[0], q[1]))) {
                Some(LineIntersection::SinglePoint { intersection, .. }) => {
                    cuts.push(parameter(intersection))
                }
                Some(LineIntersection::Collinear { intersection }) => {
                    cuts.push(parameter(intersection.start));
                    cuts.push(parameter(intersection.end));
                }
                None => {}
            }
        }
    }
    cuts.sort_by(f64::total_cmp);
    cuts.dedup_by(|a, b| (*a - *b).abs() * length <= precision);
    let mut result: Vec<(f64, f64, Location)> = vec![];
    for w in cuts.windows(2) {
        if let Some(kind) = location(
            (a + d * ((w[0] + w[1]) / 2.)).to_array(),
            contours,
            precision,
        ) {
            if let Some(last) = result.last_mut() {
                if last.2 == kind && (last.1 - w[0]).abs() * length <= precision {
                    last.1 = w[1];
                    continue;
                }
            }
            result.push((w[0], w[1], kind));
        }
    }
    result
}

struct Proposal {
    points: BTreeMap<u32, DVec3>,
    parameters: BTreeMap<u32, f64>,
    spans: Vec<SourceSpan>,
}

fn propose(
    mesh: &MeshData,
    source: &frame::Report,
    axis: &frame::Axis,
    locked: &BTreeMap<u32, DVec3>,
    unavailable: &BTreeMap<u32, &str>,
    owners: &BTreeMap<u32, Vec<usize>>,
    supports: &[PlaneFrame],
    policy: &Policy,
) -> Result<Proposal, (&'static str, Vec<u32>)> {
    let fail = |reason, nodes| Err((reason, nodes));
    let mut anchors = axis.anchors.clone();
    anchors.sort_by(|a, b| a.t.total_cmp(&b.t));
    let ids: Vec<_> = anchors.iter().map(|a| source.node_ids[a.node]).collect();
    let bad: Vec<_> = ids
        .iter()
        .filter(|n| unavailable.contains_key(n))
        .copied()
        .collect();
    if !bad.is_empty() {
        return fail("unavailable_shared_anchor", bad);
    }
    let [i, j] = axis.endpoints;
    let ends = [source.node_ids[i], source.node_ids[j]];
    if !ends.iter().all(|n| ids.contains(n)) || axis.spans.is_empty() {
        return fail("incomplete_axis_definition", ends.to_vec());
    }
    let reference = mesh.nodes[&ends[1]] - mesh.nodes[&ends[0]];
    let reference_length = reference.length();
    if reference_length <= policy.precision {
        return fail("degenerate_reference_axis", ends.to_vec());
    }
    let up = DVec3::from_array(source.policy.up).normalize();
    let fixed: Vec<_> = ids
        .iter()
        .filter_map(|n| locked.get(n).map(|&p| (*n, p)))
        .collect();
    let ca = DVec3::from_array(source.candidate_points[i]);
    let cb = DVec3::from_array(source.candidate_points[j]);
    let short = reference_length < source.policy.minimum_length;
    let raw = if fixed.len() >= 2 {
        fixed.last().unwrap().1 - fixed[0].1
    } else {
        cb - ca
    };
    if raw.length() <= policy.precision {
        return fail("coincident_axis_anchors", ends.to_vec());
    }
    let cosine = reference.normalize().dot(up).abs();
    let direction = if short {
        reference.normalize()
    } else if fixed.len() >= 2 {
        // Once both endpoints are fixed by already accepted geometry, their
        // exact chord is the only line that satisfies both endpoint
        // constraints.  This is especially important for a branch at a
        // structural joint: snapping a nearly vertical/horizontal branch to
        // the preferred global direction can create a false collinearity
        // rejection.  The angular check below still rejects a real conflict.
        raw.normalize()
    } else if axis.constructive_segment {
        raw.normalize()
    } else if cosine >= source.policy.angle.cos() {
        up * reference.dot(up).signum()
    } else if cosine <= source.policy.angle.sin() {
        (raw - up * raw.dot(up)).normalize_or_zero()
    } else {
        raw.normalize()
    };
    // With both ends fixed by accepted geometry the chord turns by what
    // their accepted movements allow: a short bar whose ends moved within
    // their budgets may turn beyond the recognition angle.
    let allowed = if fixed.len() >= 2 {
        let budgets = movement_budget(mesh, source, i) + movement_budget(mesh, source, j);
        allowed_turn(source.policy.angle, budgets, reference_length)
    } else {
        source.policy.angle
    };
    if direction.length_squared() < 0.5 || direction.dot(reference.normalize()) < allowed.cos() {
        return fail("axis_direction_conflict", ends.to_vec());
    }
    let mut origin = if let Some((_, p)) = fixed.first() {
        *p
    } else {
        (ca + cb) * 0.5
    };
    if fixed.is_empty() {
        let mut common: Option<BTreeSet<usize>> = None;
        for n in &ids {
            let owned: BTreeSet<_> = owners.get(n).into_iter().flatten().copied().collect();
            common = Some(match common {
                None => owned,
                Some(s) => s.intersection(&owned).copied().collect(),
            });
        }
        let parallel: Vec<_> = common
            .unwrap_or_default()
            .into_iter()
            .map(|i| &supports[i])
            .filter(|p| DVec3::from_array(p.normal).dot(direction).abs() <= 1e-10)
            .collect();
        origin = intersection(origin, &parallel, policy.precision)
            .ok_or(("axis_surface_incidence_conflict", ids.clone()))?;
    }
    let mut points = BTreeMap::new();
    for a in &anchors {
        let n = source.node_ids[a.node];
        let candidate = DVec3::from_array(source.candidate_points[a.node]);
        let q = if short {
            let (base, t) = if let Some((n, p)) = fixed.first() {
                (
                    *p,
                    anchors
                        .iter()
                        .find(|a| source.node_ids[a.node] == *n)
                        .unwrap()
                        .t,
                )
            } else {
                ((ca + cb) * 0.5, 0.5)
            };
            // A short axis keeps its source vector, or the vector the frame
            // flattened onto its plane.
            let vector = if source.policy.geotechnical {
                cb - ca
            } else {
                reference
            };
            base + vector * (a.t - t)
        } else {
            origin + direction * (candidate - origin).dot(direction)
        };
        let q = if let Some(&fixed) = locked.get(&n) {
            let off = if short {
                q.distance(fixed)
            } else {
                (fixed - origin).cross(direction).length()
            };
            if off > policy.precision {
                return fail("shared_anchors_not_collinear", vec![n]);
            }
            fixed
        } else if !short {
            // Interior source FE nodes may slide along the geometric axis.
            // Intersect that line with all owning planes; never pin their
            // arbitrary tangential coordinates from the approximate solver.
            let mut t = (candidate - origin).dot(direction);
            let owned: Vec<_> = owners
                .get(&n)
                .into_iter()
                .flatten()
                .map(|&i| &supports[i])
                .collect();
            for plane in &owned {
                // Avoid dividing a satisfied coplanar residual by a tiny slope.
                if plane.distance((origin + direction * t).to_array()).abs() <= policy.precision {
                    continue;
                }
                let slope = DVec3::from_array(plane.normal).dot(direction);
                if slope.abs() > 1e-10 {
                    t = -plane.distance(origin.to_array()) / slope;
                    break;
                }
            }
            let p = origin + direction * t;
            if owned
                .iter()
                .any(|plane| plane.distance(p.to_array()).abs() > policy.precision)
            {
                return fail("axis_surface_incidence_conflict", vec![n]);
            }
            p
        } else {
            q
        };
        if owners
            .get(&n)
            .into_iter()
            .flatten()
            .any(|&i| supports[i].distance(q.to_array()).abs() > policy.precision)
        {
            return fail("axis_surface_incidence_conflict", vec![n]);
        }
        if !q.is_finite()
            || q.distance(mesh.nodes[&n]) > movement_budget(mesh, source, a.node) + policy.precision
            || q.distance(candidate) > policy.junction_movement_limit
        {
            return fail("axis_movement_budget", vec![n]);
        }
        points.insert(n, q);
    }
    let (a, b) = (points[&ends[0]], points[&ends[1]]);
    let d = b - a;
    let length = d.length();
    // Each end may move by its budget, so a short axis may shorten by both
    // (a 10 mm bar whose ends lie on surface planes a few micrometres
    // closer); it must keep a length and its orientation.
    let slack = [i, j]
        .map(|e| movement_budget(mesh, source, e))
        .iter()
        .sum::<f64>();
    if length <= policy.precision
        || length + policy.precision + slack < source.policy.minimum_length.min(reference_length)
        || d.dot(reference) <= 0.
    {
        return fail("axis_length_or_orientation", ends.to_vec());
    }
    let parameters: BTreeMap<_, _> = points
        .iter()
        .map(|(&n, &p)| (n, (p - a).dot(d) / d.length_squared()))
        .collect();
    for w in ids.windows(2) {
        if (parameters[&w[1]] - parameters[&w[0]]) * length <= policy.precision {
            return fail("collapsed_or_reordered_property_span", w.to_vec());
        }
    }
    if points
        .values()
        .any(|p| (*p - a).cross(d / length).length() > policy.precision)
    {
        return fail("axis_straightness", ids);
    }
    let mut spans = vec![];
    for span in &axis.spans {
        let start = anchors
            .iter()
            .find(|a| a.t == span.start_t)
            .ok_or(("missing_span_anchor", vec![]))?;
        let end = anchors
            .iter()
            .find(|a| a.t == span.end_t)
            .ok_or(("missing_span_anchor", vec![]))?;
        spans.push(SourceSpan {
            element: span.element,
            stiffness: span.stiffness,
            start_t: parameters[&source.node_ids[start.node]],
            end_t: parameters[&source.node_ids[end.node]],
        });
    }
    Ok(Proposal {
        points,
        parameters,
        spans,
    })
}

pub(super) fn assemble(
    mesh: &MeshData,
    source: &frame::Report,
    model: &mut Model,
    vertex_source_nodes: &mut Vec<u32>,
    supports: &[PlaneFrame],
    representatives: &[usize],
    policy: &Policy,
    removed_elements: &BTreeSet<u32>,
) -> Report {
    let mut report = Report::default();
    let mut vertices: BTreeMap<_, _> = vertex_source_nodes
        .iter()
        .enumerate()
        .map(|(i, &n)| (n, i))
        .collect();
    let mut locked: BTreeMap<_, _> = vertices
        .iter()
        .map(|(&n, &v)| (n, DVec3::from_array(model.vertices[v])))
        .collect();
    let mut unavailable = BTreeMap::new();
    let mut incidence = BTreeMap::<u32, BTreeSet<usize>>::new();
    for (i, axis) in source.axes.iter().enumerate() {
        for a in &axis.anchors {
            incidence
                .entry(source.node_ids[a.node])
                .or_default()
                .insert(i);
        }
    }
    let mut owner_planes = BTreeMap::<u32, Vec<usize>>::new();
    for (i, surface) in source.surfaces.iter().enumerate() {
        for &node in &surface.nodes {
            let n = source.node_ids[node];
            if incidence.contains_key(&n) {
                owner_planes.entry(n).or_default().push(representatives[i]);
            }
        }
    }
    let element_surface: BTreeMap<_, _> = model
        .surfaces
        .iter()
        .enumerate()
        .flat_map(|(i, s)| s.source_elements.iter().map(move |&e| (e, i)))
        .collect();
    let source_shells: BTreeSet<_> = source
        .surfaces
        .iter()
        .flat_map(|s| s.source_elements.iter().copied())
        .collect();
    let mut owner_surfaces = BTreeMap::<u32, BTreeSet<usize>>::new();
    for e in &mesh.elements {
        // Elements removed as degenerate slivers are absent, not unbuilt.
        if !source_shells.contains(&e.id) || removed_elements.contains(&e.id) {
            continue;
        }
        for &n in &e.nodes {
            if !incidence.contains_key(&n) {
                continue;
            }
            if let Some(&s) = element_surface.get(&e.id) {
                owner_surfaces.entry(n).or_default().insert(s);
            } else {
                unavailable.insert(n, "unbuilt_surface_region");
            }
        }
    }
    let lookup: BTreeMap<_, _> = source
        .node_ids
        .iter()
        .enumerate()
        .map(|(i, &n)| (n, i))
        .collect();
    for (&n, axes) in &incidence {
        if locked.contains_key(&n) {
            continue;
        }
        let candidate = DVec3::from_array(source.candidate_points[lookup[&n]]);
        if let Some(owned) = owner_planes.get(&n) {
            let planes: Vec<_> = owned.iter().map(|&p| &supports[p]).collect();
            if let Some(q) = intersection(candidate, &planes, policy.precision) {
                if axes.len() > 1 {
                    locked.insert(n, q);
                }
            } else {
                unavailable.insert(n, "incompatible_surface_supports");
            }
        } else if axes.len() > 1 {
            locked.insert(n, candidate);
        }
    }
    report.boundary_repairs = repair::repair(
        model,
        &mut locked,
        &repair::Context {
            mesh,
            source,
            vertices: &vertices,
            incidence: &incidence,
            unavailable: &unavailable,
            owner_planes: &owner_planes,
            owner_surfaces: &owner_surfaces,
            supports,
            policy,
        },
    );
    // Shared junctions are fixed before any axis is considered. Accepted axes
    // cannot move a neighbor or depend on processing order.
    for (source_axis, axis) in source.axes.iter().enumerate() {
        let result = (|| {
            let proposal = propose(
                mesh,
                source,
                axis,
                &locked,
                &unavailable,
                &owner_planes,
                supports,
                policy,
            )?;
            for (&n, &p) in &proposal.points {
                for &s in owner_surfaces.get(&n).into_iter().flatten() {
                    let surface = &model.surfaces[s];
                    let plane = &model.planes[surface.plane];
                    if plane.distance(p.to_array()).abs() > policy.precision
                        || location(
                            plane.project(p.to_array()),
                            &surface.contours,
                            policy.precision,
                        )
                        .is_none()
                    {
                        return Err(("anchor_outside_source_surface", vec![n]));
                    }
                }
            }
            Ok(proposal)
        })();
        let proposal = match result {
            Ok(p) => p,
            Err((reason, nodes)) => {
                report.issues.push(Issue {
                    source_axis,
                    source_elements: axis.spans.iter().map(|s| s.element).collect(),
                    source_nodes: nodes,
                    reason: reason.into(),
                });
                continue;
            }
        };
        let output_axis = report.axes.len();
        for (&n, &p) in &proposal.points {
            report.maximum_additional_movement = report
                .maximum_additional_movement
                .max(p.distance(DVec3::from_array(source.candidate_points[lookup[&n]])));
            if !vertices.contains_key(&n) {
                let v = model
                    .add_vertex(p.to_array())
                    .expect("validated finite axis point");
                vertices.insert(n, v);
                vertex_source_nodes.push(n);
            }
        }
        let endpoints = axis.endpoints.map(|i| vertices[&source.node_ids[i]]);
        let (a, b) = (model.vertices[endpoints[0]], model.vertices[endpoints[1]]);
        let mut touched = BTreeSet::new();
        for anchor in &axis.anchors {
            let n = source.node_ids[anchor.node];
            for &s in owner_surfaces.get(&n).into_iter().flatten() {
                let surface = &model.surfaces[s];
                let plane = &model.planes[surface.plane];
                report.contacts.push(Contact::Point {
                    axis: output_axis,
                    surface: s,
                    vertex: vertices[&n],
                    t: proposal.parameters[&n],
                    location: location(
                        plane.project(proposal.points[&n].to_array()),
                        &surface.contours,
                        policy.precision,
                    )
                    .unwrap(),
                });
                touched.insert(s);
            }
        }
        for s in touched {
            let surface = &model.surfaces[s];
            let plane = &model.planes[surface.plane];
            if plane.distance(a).abs() <= policy.precision
                && plane.distance(b).abs() <= policy.precision
            {
                for (start_t, end_t, location) in intervals(
                    plane.project(a),
                    plane.project(b),
                    &surface.contours,
                    policy.precision,
                ) {
                    report.contacts.push(Contact::Interval {
                        axis: output_axis,
                        surface: s,
                        start_t,
                        end_t,
                        location,
                    });
                }
            }
        }
        report.axes.push(Axis {
            source_axis,
            endpoints,
            spans: proposal.spans,
            anchors: axis
                .anchors
                .iter()
                .map(|a| {
                    let n = source.node_ids[a.node];
                    Anchor {
                        source_node: n,
                        vertex: vertices[&n],
                        t: proposal.parameters[&n],
                    }
                })
                .collect(),
        });
    }
    report.rejected_shared_anchors = unavailable
        .into_iter()
        .map(|(n, r)| (n, r.into()))
        .collect();
    report.all_axes_built = report.issues.is_empty();
    report
}

/// Recompute bar-surface contacts from the final geometry.
///
/// Contacts are first derived from source ownership, before vertex merges,
/// junction snapping and console trimming move vertices and rebuild
/// contours. A contact computed then can miss an in-plane interval or keep a
/// stale location, and the mesh would leave the bar unconnected along it.
/// Pairs are those already in contact plus every surface having an anchor as
/// a contour vertex: a point contact exists for each anchor on the surface
/// material, and an interval for each part of an axis lying in its plane.
pub fn refresh_contacts(model: &Model, axes: &[Axis], contacts: &mut Vec<Contact>) {
    let precision = model.precision;
    let mut owners = BTreeMap::<usize, BTreeSet<usize>>::new();
    for s in 0..model.surfaces.len() {
        for e in model.surface_edges(s) {
            for v in model.edges[e] {
                owners.entry(v).or_default().insert(s);
            }
        }
    }
    let mut previous = BTreeSet::new();
    let mut pairs = BTreeSet::new();
    for c in contacts.iter() {
        match *c {
            Contact::Point {
                axis,
                surface,
                vertex,
                ..
            } => {
                previous.insert((axis, surface, vertex));
                pairs.insert((axis, surface));
            }
            Contact::Interval { axis, surface, .. } => {
                pairs.insert((axis, surface));
            }
        }
    }
    for (i, axis) in axes.iter().enumerate() {
        for anchor in &axis.anchors {
            for &s in owners.get(&anchor.vertex).into_iter().flatten() {
                pairs.insert((i, s));
            }
        }
    }
    // A bar node on the plane and material of a surface that does not use
    // it (a column through a slab without a shared source node) is a point
    // contact of that surface too.
    let boxes: Vec<(DVec3, DVec3)> = (0..model.surfaces.len())
        .map(|s| {
            model
                .surface_edges(s)
                .flat_map(|e| model.edges[e])
                .map(|v| DVec3::from_array(model.vertices[v]))
                .fold(
                    (DVec3::splat(f64::INFINITY), DVec3::splat(f64::NEG_INFINITY)),
                    |(lo, hi), p| (lo.min(p), hi.max(p)),
                )
        })
        .collect();
    for (i, axis) in axes.iter().enumerate() {
        for anchor in &axis.anchors {
            let p = model.vertices[anchor.vertex];
            let q = DVec3::from_array(p);
            for (s, &(lo, hi)) in boxes.iter().enumerate() {
                if (q + precision).cmplt(lo).any() || (q - precision).cmpgt(hi).any() {
                    continue;
                }
                let surface = &model.surfaces[s];
                let plane = &model.planes[surface.plane];
                if plane.distance(p).abs() <= precision
                    && location(plane.project(p), &surface.contours, precision).is_some()
                {
                    previous.insert((i, s, anchor.vertex));
                    pairs.insert((i, s));
                }
            }
        }
    }
    let mut result = vec![];
    for (i, s) in pairs {
        let (Some(axis), Some(surface)) = (axes.get(i), model.surfaces.get(s)) else {
            continue;
        };
        let plane = &model.planes[surface.plane];
        for anchor in &axis.anchors {
            let p = model.vertices[anchor.vertex];
            // Geometry decides: a node off the plane is no contact, whatever
            // an earlier record or ownership said.
            if plane.distance(p).abs() > precision {
                continue;
            }
            if let Some(location) = location(plane.project(p), &surface.contours, precision) {
                result.push(Contact::Point {
                    axis: i,
                    surface: s,
                    vertex: anchor.vertex,
                    t: anchor.t,
                    location,
                });
            }
        }
        let (a, b) = (
            model.vertices[axis.endpoints[0]],
            model.vertices[axis.endpoints[1]],
        );
        if plane.distance(a).abs() <= precision && plane.distance(b).abs() <= precision {
            for (start_t, end_t, location) in intervals(
                plane.project(a),
                plane.project(b),
                &surface.contours,
                precision,
            ) {
                result.push(Contact::Interval {
                    axis: i,
                    surface: s,
                    start_t,
                    end_t,
                    location,
                });
            }
        }
    }
    *contacts = result;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        input::ElementData,
        reconstruction::{assembly, planes, recognize},
    };

    #[test]
    fn short_bar_may_turn_as_far_as_its_end_budgets_allow() {
        let angle = 0.02;
        // A long bar: the recognition angle governs.
        assert_eq!(allowed_turn(angle, 0.01, 10.), angle);
        // A 164 mm bar whose ends may move 8.2 mm each: about 5.7 degrees.
        let turn = allowed_turn(angle, 2. * 0.05 * 0.164, 0.164);
        assert!((turn - 0.1f64.asin()).abs() < 1e-12 && turn > 2.3f64.to_radians());
        // Budgets beyond the length never exceed a right angle.
        assert_eq!(allowed_turn(angle, 1., 0.1), std::f64::consts::FRAC_PI_2);
    }

    #[test]
    fn bar_passing_through_a_slab_shares_a_crossing_vertex() {
        use super::super::junctions::tests::{build, slab, Placement};
        for place in Placement::all() {
            // Through the material, through the edge x = 0, and past it.
            for (x, crossing, on_edge) in [(1., true, false), (0., true, true), (-1., false, false)]
            {
                let mut m = build(&place, &[slab(0., 4.)]);
                let before = m.surface_edges(0).count();
                let a = m.add_vertex(place.point([x, 2., -1.])).unwrap();
                let b = m.add_vertex(place.point([x, 2., 1.])).unwrap();
                let mut axes = vec![Axis {
                    source_axis: 0,
                    endpoints: [a, b],
                    anchors: vec![
                        Anchor {
                            source_node: 1,
                            vertex: a,
                            t: 0.,
                        },
                        Anchor {
                            source_node: 2,
                            vertex: b,
                            t: 1.,
                        },
                    ],
                    spans: vec![],
                }];
                let r = imprint_crossings(&mut m, &mut axes);
                assert_eq!(r.len(), usize::from(crossing), "x {x}: {r:?}");
                if crossing {
                    let c = &r[0];
                    assert_eq!((c.kind.as_str(), c.surface), ("crossing", Some(0)));
                    assert!((c.t - 0.5).abs() < 1e-9);
                    assert_eq!(axes[0].anchors[1].source_node, NO_SOURCE_NODE);
                    let p = DVec3::from_array(m.vertices[c.vertex]);
                    assert!(
                        p.distance(DVec3::from_array(place.point([x, 2., 0.])))
                            < 1e-9 * place.scale
                    );
                    // On the edge, the contour edge is split at the crossing.
                    assert_eq!(m.surface_edges(0).count(), before + usize::from(on_edge));
                }
            }
        }
    }

    #[test]
    fn beam_node_next_to_a_wall_foot_slides_onto_it() {
        use super::super::junctions::tests::{build, run, slab, wall, Placement};
        for place in Placement::all() {
            // A beam lying in the slab crosses the foot of a wall standing on
            // it; its node 5 um before the foot line slides onto it and
            // splits the foot line. 2 mm before it, a crossing is generated.
            for (offset, slid) in [(5e-6, true), (0.002, false)] {
                let mut m = build(&place, &[slab(0., 4.), wall(1., 3., 0., 2.)]);
                run(&mut m);
                let mut node = |y: f64| m.add_vertex(place.point([2., y, 0.])).unwrap();
                let (a, b, c) = (node(0.5), node(2. - offset), node(3.5));
                let t = (1.5 - offset) / 3.;
                let mut axes = vec![Axis {
                    source_axis: 0,
                    endpoints: [a, c],
                    anchors: [(1, a, 0.), (2, b, t), (3, c, 1.)]
                        .map(|(source_node, vertex, t)| Anchor {
                            source_node,
                            vertex,
                            t,
                        })
                        .to_vec(),
                    spans: [(1, 0., t), (2, t, 1.)]
                        .map(|(element, start_t, end_t)| SourceSpan {
                            element,
                            stiffness: 1,
                            start_t,
                            end_t,
                        })
                        .to_vec(),
                }];
                let r = imprint_crossings(&mut m, &mut axes);
                assert_eq!(r.len(), 1, "{r:?}");
                let foot = DVec3::from_array(place.point([2., 2., 0.]));
                let v = r[0].vertex;
                assert!(DVec3::from_array(m.vertices[v]).distance(foot) < 1e-9 * place.scale);
                assert!((r[0].t - 0.5).abs() < 1e-9);
                // The foot line (wall contour, slab junction) is split at the
                // shared vertex.
                for surface in [0, 1] {
                    assert!(m.surface_edges(surface).any(|e| m.edges[e].contains(&v)));
                }
                if slid {
                    assert_eq!(
                        (r[0].kind.as_str(), r[0].source_node, v),
                        ("slid_anchor", 2, b)
                    );
                    assert!((r[0].movement - offset * place.scale).abs() < 1e-9 * place.scale);
                    assert_eq!(axes[0].anchors.len(), 3);
                    // The element spans bounded by the node follow it.
                    let spans: Vec<_> =
                        axes[0].spans.iter().map(|s| (s.start_t, s.end_t)).collect();
                    assert_eq!(spans, vec![(0., r[0].t), (r[0].t, 1.)]);
                } else {
                    assert_eq!(r[0].kind, "crossing");
                    assert_eq!(axes[0].anchors.len(), 4);
                }
            }
        }
    }

    #[test]
    fn crossing_bars_share_a_generated_vertex() {
        use super::super::junctions::tests::Placement;
        for place in Placement::all() {
            // Bar A along y = 1; bar B along x = 1 crosses it at (1, 1, 0).
            // 2 mm above A, B misses it; crossing 0.5 mm from A's end, the
            // near touch is left alone.
            for (b0, b1, crossing) in [
                ([1., 0., 0.], [1., 2., 0.], true),
                ([1., 0., 0.002], [1., 2., 0.002], false),
                ([1.9995, 0., 0.], [1.9995, 2., 0.], false),
            ] {
                let mut m = Model::new(1e-7 * place.scale, 0.001 * place.scale).unwrap();
                let mut bar = |p: [f64; 3], q: [f64; 3], source_axis: usize| {
                    let [a, b] = [p, q].map(|x| m.add_vertex(place.point(x)).unwrap());
                    Axis {
                        source_axis,
                        endpoints: [a, b],
                        anchors: [(a, 0.), (b, 1.)]
                            .map(|(vertex, t)| Anchor {
                                source_node: vertex as u32,
                                vertex,
                                t,
                            })
                            .to_vec(),
                        spans: vec![],
                    }
                };
                let mut axes = vec![bar([0., 1., 0.], [2., 1., 0.], 0), bar(b0, b1, 1)];
                if crossing {
                    // A third bar through the same point shares the vertex.
                    axes.push(bar([0., 0., 0.], [2., 2., 0.], 2));
                }
                let r = imprint_bar_crossings(&mut m, &mut axes);
                assert_eq!(r.len(), if crossing { 3 } else { 0 }, "{r:?}");
                if crossing {
                    let v = r[0].vertex;
                    assert!(r.iter().all(|c| c.vertex == v));
                    let p = DVec3::from_array(m.vertices[v]);
                    let expected = DVec3::from_array(place.point([1., 1., 0.]));
                    assert!(p.distance(expected) < 1e-9 * place.scale);
                    for axis in &axes {
                        assert_eq!(axis.anchors.len(), 3);
                        assert_eq!(axis.anchors[1].vertex, v);
                        assert!((axis.anchors[1].t - 0.5).abs() < 1e-9);
                    }
                }
            }
        }
    }

    #[test]
    fn crossings_next_to_each_other_are_both_shared() {
        use super::super::junctions::tests::Placement;
        for place in Placement::all() {
            // Bars B and C fan out from (0, 0, 0) at about 12 degrees; bar A
            // along y = 0.004 crosses both 4 mm from their node, 0.8 mm
            // apart: each crossing gets its own vertex.
            let mut m = Model::new(1e-7 * place.scale, 0.001 * place.scale).unwrap();
            let mut bar = |p: [f64; 3], q: [f64; 3], source_axis: usize| {
                let [a, b] = [p, q].map(|x| m.add_vertex(place.point(x)).unwrap());
                Axis {
                    source_axis,
                    endpoints: [a, b],
                    anchors: [(a, 0.), (b, 1.)]
                        .map(|(vertex, t)| Anchor {
                            source_node: vertex as u32,
                            vertex,
                            t,
                        })
                        .to_vec(),
                    spans: vec![],
                }
            };
            let mut axes = vec![
                bar([-1., 0.004, 0.], [1., 0.004, 0.], 0),
                bar([0., 0., 0.], [1., 1., 0.], 1),
                bar([0., 0., 0.], [0.8, 1., 0.], 2),
            ];
            let r = imprint_bar_crossings(&mut m, &mut axes);
            assert_eq!(r.len(), 4, "{r:?}");
            assert_eq!(axes[0].anchors.len(), 4);
            assert!(axes[1..].iter().all(|a| a.anchors.len() == 3));
        }
    }

    #[test]
    fn node_next_to_a_bar_crossing_slides_onto_it_and_is_shared() {
        use super::super::junctions::tests::{build, slab, Placement};
        for place in Placement::all() {
            // Bar A along y = 4 has a node at the slab corner (4, 4, 0); bar
            // B along x = 4.0004 crosses A 0.4 mm from it (beyond the point
            // slack, within the minimum edge). The corner slides along A onto
            // the crossing, staying in the slab plane, and becomes a node of
            // B; A's element spans follow it, and bar C ending at the corner
            // follows it with its end.
            let mut m = build(&place, &[slab(0., 4.)]);
            let at = |m: &Model, p: [f64; 3]| {
                (0..m.vertices.len()).find(|&v| {
                    DVec3::from_array(m.vertices[v]).distance(DVec3::from_array(place.point(p)))
                        < 1e-9 * place.scale
                })
            };
            let corner = at(&m, [4., 4., 0.]).unwrap();
            let mut v = |p: [f64; 3]| m.add_vertex(place.point(p)).unwrap();
            let (a0, a1) = (v([2., 4., 0.]), v([6., 4., 0.]));
            let (b0, b1) = (v([4.0004, 3., 0.]), v([4.0004, 5., 0.]));
            let c1 = v([5., 6., 0.]);
            let anchor = |vertex: usize, t: f64| Anchor {
                source_node: vertex as u32,
                vertex,
                t,
            };
            let span = |start_t: f64, end_t: f64, element: u32| SourceSpan {
                element,
                stiffness: 0,
                start_t,
                end_t,
            };
            let mut axes = vec![
                Axis {
                    source_axis: 0,
                    endpoints: [a0, a1],
                    anchors: vec![anchor(a0, 0.), anchor(corner, 0.5), anchor(a1, 1.)],
                    spans: vec![span(0., 0.5, 1), span(0.5, 1., 2)],
                },
                Axis {
                    source_axis: 1,
                    endpoints: [b0, b1],
                    anchors: vec![anchor(b0, 0.), anchor(b1, 1.)],
                    spans: vec![span(0., 1., 3)],
                },
                Axis {
                    source_axis: 2,
                    endpoints: [corner, c1],
                    anchors: vec![anchor(corner, 0.), anchor(c1, 1.)],
                    spans: vec![span(0., 1., 4)],
                },
            ];
            let r = imprint_bar_crossings(&mut m, &mut axes);
            let kinds: Vec<_> = r.iter().map(|c| (c.axis, c.kind.as_str())).collect();
            assert_eq!(kinds, vec![(0, "slid_anchor"), (1, "shared_node")]);
            assert!((r[0].movement - 0.0004 * place.scale).abs() < 1e-9 * place.scale);
            let p = DVec3::from_array(m.vertices[corner]);
            assert!(
                p.distance(DVec3::from_array(place.point([4.0004, 4., 0.]))) < 1e-9 * place.scale
            );
            assert!(axes[1]
                .anchors
                .iter()
                .any(|x| x.vertex == corner && (x.t - 0.5).abs() < 1e-9));
            let t = axes[0].anchors[1].t;
            assert!((t - 2.0004 / 4.).abs() < 1e-9);
            let spans: Vec<_> = axes[0].spans.iter().map(|s| (s.start_t, s.end_t)).collect();
            assert_eq!(spans, vec![(0., t), (t, 1.)]);
            assert_eq!(axes[2].endpoints[0], corner);
        }
    }

    #[test]
    fn surface_corner_on_a_bar_axis_becomes_a_bar_anchor() {
        use super::super::junctions::tests::{build, slab, Placement};
        for place in Placement::all() {
            for (offset, imprinted) in [(0., true), (1e-5, true), (0.002, false)] {
                // A bar from x = -1 to x = 2 passes the slab corner (0, 0, 0)
                // and runs along its edge. 10 um off, the corner moves onto
                // the axis (within the minimum edge); 2 mm off, it stays.
                let mut m = build(&place, &[slab(0., 4.)]);
                let corner = (0..m.vertices.len())
                    .find(|&v| {
                        DVec3::from_array(m.vertices[v])
                            .distance(DVec3::from_array(place.point([0., 0., 0.])))
                            < 1e-9 * place.scale
                    })
                    .unwrap();
                let a = m.add_vertex(place.point([-1., offset, 0.])).unwrap();
                let b = m.add_vertex(place.point([2., offset, 0.])).unwrap();
                let mut axes = vec![Axis {
                    source_axis: 0,
                    endpoints: [a, b],
                    anchors: vec![
                        Anchor {
                            source_node: 1,
                            vertex: a,
                            t: 0.,
                        },
                        Anchor {
                            source_node: 2,
                            vertex: b,
                            t: 1.,
                        },
                    ],
                    spans: vec![],
                }];
                let source_nodes: Vec<u32> =
                    (0..m.vertices.len() as u32).map(|v| 100 + v).collect();
                let r = imprint_surface_vertices(&mut m, &mut axes, &source_nodes);
                assert_eq!(r.len(), usize::from(imprinted), "{r:?}");
                if imprinted {
                    let anchors: Vec<_> = axes[0].anchors.iter().map(|x| x.vertex).collect();
                    assert_eq!(anchors, vec![a, corner, b]);
                    let (pa, pb) = (
                        DVec3::from_array(m.vertices[a]),
                        DVec3::from_array(m.vertices[b]),
                    );
                    let pc = DVec3::from_array(m.vertices[corner]);
                    assert!(
                        (pc - pa).cross(pb - pa).length() / (pb - pa).length() < 1e-9 * place.scale
                    );
                    assert!((axes[0].anchors[1].t - 1. / 3.).abs() < 1e-9);
                    assert_eq!(axes[0].anchors[1].source_node, 100 + corner as u32);
                }
            }
        }
    }

    #[test]
    fn axis_leaving_a_boundary_at_a_nearly_collinear_vertex_is_cut_there() {
        // The axis runs along the bottom edge and leaves the material at the
        // vertex (1, 1e-9), within precision of it but not exactly on it.
        let square = vec![vec![[0., 0.], [0.5, 0.], [1., 1e-9], [1., 1.], [0., 1.]]];
        let r = intervals([0., 0.], [1.5, 0.], &square, 1e-8);
        assert_eq!(r.len(), 1, "{r:?}");
        assert_eq!(r[0].2, Location::Boundary);
        assert!(
            r[0].0.abs() < 1e-12 && (r[0].1 - 2. / 3.).abs() < 1e-8,
            "{r:?}"
        );
    }

    fn frame(mesh: &MeshData, up: DVec3) -> frame::Report {
        let axes = recognize::recognize(
            mesh,
            &recognize::Policy {
                angle: 0.02,
                line_tolerance: 0.01,
                numerical_precision: 1e-8,
            },
        )
        .unwrap();
        let planes = planes::recognize(
            mesh,
            &planes::Policy {
                angle: 0.02,
                distance: 0.01,
                precision: 1e-8,
            },
        )
        .unwrap();
        frame::solve(
            mesh,
            &axes,
            &planes,
            &frame::Policy {
                up: up.to_array(),
                angle: 0.02,
                maximum_movement: 0.15,
                relative_movement: 0.05,
                minimum_length: 0.03,
                residual_tolerance: 1e-7,
                iterations: 100,
                panel_tolerance: 0.,
                geotechnical: false,
                over_constrained_panels: false,
            },
        )
        .unwrap()
    }
    fn policy(scale: f64) -> Policy {
        Policy {
            closure_tolerance: 0.001 * scale,
            junction_movement_limit: 0.05 * scale,
            precision: 1e-7 * scale,
            minimum_edge: 0.001 * scale,
        }
    }
    fn slab_beam_column() -> MeshData {
        let mut mesh = MeshData::default();
        for y in 0..3 {
            for x in 0..3 {
                mesh.nodes
                    .insert(1 + x + 3 * y, DVec3::new(x as f64, y as f64, 0.));
            }
        }
        for y in 0..2 {
            for x in 0..2 {
                let a = 1 + x + 3 * y;
                mesh.elements.push(ElementData {
                    id: mesh.elements.len() as u32 + 1,
                    elem_type: 44,
                    stiff_id: 1,
                    nodes: vec![a, a + 1, a + 4, a + 3],
                });
            }
        }
        mesh.nodes.insert(10, DVec3::new(1., 1., -2.));
        mesh.elements.extend([
            ElementData {
                id: 5,
                elem_type: 10,
                stiff_id: 10,
                nodes: vec![4, 5],
            },
            ElementData {
                id: 6,
                elem_type: 10,
                stiff_id: 20,
                nodes: vec![5, 6],
            },
            ElementData {
                id: 7,
                elem_type: 10,
                stiff_id: 30,
                nodes: vec![10, 5],
            },
        ]);
        mesh
    }
    #[test]
    fn shared_material_boundary_moves_with_axis_under_transforms() {
        for scale in [0.1, 1., 10.] {
            for transformed in [false, true] {
                let mut mesh = slab_beam_column();
                mesh.elements.retain(|e| e.id != 7);
                for e in &mut mesh.elements {
                    if e.id == 3 || e.id == 4 {
                        e.stiff_id = 2;
                    }
                }
                let q = if transformed {
                    glam::DQuat::from_axis_angle(DVec3::new(1., 2., 3.).normalize(), 0.6)
                } else {
                    glam::DQuat::IDENTITY
                };
                let shift = DVec3::new(20., -30., 50.);
                let number = |n| if transformed { 100 - n } else { n };
                mesh.nodes = mesh
                    .nodes
                    .iter()
                    .map(|(&n, &p)| (number(n), q * (p * scale) + shift))
                    .collect();
                for e in &mut mesh.elements {
                    e.id = number(e.id);
                    for n in &mut e.nodes {
                        *n = number(*n);
                    }
                    e.nodes.reverse();
                }
                if transformed {
                    mesh.elements.reverse();
                }
                let mut f = frame(&mesh, q * DVec3::Z);
                let i = f.node_ids.iter().position(|&n| n == number(5)).unwrap();
                f.candidate_points[i] =
                    (q * (DVec3::new(1., 1.002, 0.) * scale) + shift).to_array();
                let r = assembly::assemble(&mesh, &f, &policy(scale)).unwrap();
                assert!(r.all_surface_patches_built);
                assert!(
                    r.axis_assembly.all_axes_built,
                    "{:?}",
                    r.axis_assembly.issues
                );
                assert_eq!(r.axis_assembly.boundary_repairs.len(), 1);
                let repair = &r.axis_assembly.boundary_repairs[0];
                assert!(repair.accepted, "{}", repair.reason);
                assert_eq!(repair.surfaces.len(), 2);
                assert_eq!(repair.changes.len(), 1);
                assert_eq!(repair.changes[0].source_node, number(5));
                let p = DVec3::from_array(repair.changes[0].after);
                assert!(p.distance(mesh.nodes[&number(5)]) < policy(scale).precision);
                assert_eq!(
                    r.vertex_source_nodes
                        .iter()
                        .filter(|&&n| n == number(5))
                        .count(),
                    1
                );
                for s in &r.preview.surfaces {
                    assert!(
                        (geo::Area::unsigned_area(&polygon(&s.contours)) - 2. * scale * scale)
                            .abs()
                            < 1e-6 * scale * scale
                    );
                }
            }
        }
    }

    #[test]
    fn boundary_repair_rolls_back_when_neighbor_would_self_touch() {
        let mut mesh = slab_beam_column();
        mesh.elements.retain(|e| e.id != 7);
        for e in &mut mesh.elements {
            if e.id == 3 || e.id == 4 {
                e.stiff_id = 2;
            }
        }
        for n in [7, 8, 9] {
            mesh.nodes.get_mut(&n).unwrap().y = 1.004;
        }
        let mut f = frame(&mesh, DVec3::Z);
        for (i, &n) in f.node_ids.iter().enumerate() {
            if n == 5 {
                f.candidate_points[i][1] = 0.998;
            }
            if n == 8 {
                f.candidate_points[i][1] = 1.;
            }
        }
        let r = assembly::assemble(&mesh, &f, &policy(1.)).unwrap();
        assert!(r.all_surface_patches_built, "{:?}", r.issues);
        assert_eq!(r.axis_assembly.boundary_repairs.len(), 1);
        let repair = &r.axis_assembly.boundary_repairs[0];
        assert!(!repair.accepted);
        assert_eq!(repair.reason, "invalid_neighbor_contour");
        assert!(repair.changes.is_empty());
        let v = r.vertex_source_nodes.iter().position(|&n| n == 5).unwrap();
        assert!((r.preview.vertices[v][1] - 0.998).abs() < 1e-7);
        assert!(!r.axis_assembly.all_axes_built);
    }
    #[test]
    fn interior_joint_shares_one_vertex_with_constructive_segments_and_properties() {
        for scale in [0.1, 1., 10.] {
            for transformed in [false, true] {
                let mut mesh = slab_beam_column();
                let q = if transformed {
                    glam::DQuat::from_axis_angle(DVec3::new(1., 2., 3.).normalize(), 0.6)
                } else {
                    glam::DQuat::IDENTITY
                };
                let shift = DVec3::new(20., -30., 50.);
                let id = |n| if transformed { 100 - n } else { n };
                mesh.nodes = mesh
                    .nodes
                    .iter()
                    .map(|(&n, &p)| (id(n), q * (p * scale) + shift))
                    .collect();
                for e in &mut mesh.elements {
                    e.id = id(e.id);
                    for n in &mut e.nodes {
                        *n = id(*n);
                    }
                    e.nodes.reverse();
                }
                if transformed {
                    mesh.elements.reverse();
                }
                let f = frame(&mesh, q * DVec3::Z);
                let result = assembly::assemble(&mesh, &f, &policy(scale)).unwrap();
                let r = &result.axis_assembly;
                assert!(result.all_surface_patches_built);
                assert!(r.all_axes_built, "{:?}", r.issues);
                assert_eq!(r.axes.len(), 3);
                let beams: Vec<_> = r
                    .axes
                    .iter()
                    .filter(|a| a.spans.len() == 1 && matches!(a.spans[0].stiffness, 10 | 20))
                    .collect();
                assert_eq!(beams.len(), 2);
                assert_eq!(
                    beams
                        .iter()
                        .map(|a| a.spans[0].stiffness)
                        .collect::<BTreeSet<_>>(),
                    BTreeSet::from([10, 20])
                );
                let middle = r
                    .axes
                    .iter()
                    .flat_map(|a| a.anchors.iter())
                    .find(|a| a.source_node == id(5))
                    .unwrap();
                assert!(beams
                    .iter()
                    .all(|beam| beam.endpoints.contains(&middle.vertex)));
                let column = r.axes.iter().find(|a| a.spans[0].stiffness == 30).unwrap();
                assert!(column.endpoints.contains(&middle.vertex));
                assert_eq!(
                    result
                        .vertex_source_nodes
                        .iter()
                        .filter(|&&n| n == id(5))
                        .count(),
                    1
                );
                assert!(r.contacts.iter().any(|c|matches!(c,Contact::Point {vertex,location:Location::Interior,..} if *vertex==middle.vertex)));
                assert!(r.contacts.iter().any(
                    |c| matches!(c,Contact::Interval {start_t,end_t,location:Location::Interior,..}
                if start_t.abs()<1e-7 && (*end_t-1.).abs()<1e-7)
                ));
                let sources: BTreeSet<_> = r
                    .axes
                    .iter()
                    .flat_map(|a| a.spans.iter().map(|s| s.element))
                    .collect();
                assert_eq!(sources, BTreeSet::from([id(5), id(6), id(7)]));
                for axis in &r.axes {
                    let a = DVec3::from_array(result.preview.vertices[axis.endpoints[0]]);
                    let b = DVec3::from_array(result.preview.vertices[axis.endpoints[1]]);
                    for anchor in &axis.anchors {
                        assert!(
                            DVec3::from_array(result.preview.vertices[anchor.vertex])
                                .distance(a.lerp(b, anchor.t))
                                < policy(scale).precision
                        );
                    }
                }
            }
        }
    }
    #[test]
    fn interior_fe_nodes_slide_on_axis_and_property_parameters_follow() {
        let mut mesh = MeshData::default();
        for y in 0..3 {
            for x in 0..5 {
                mesh.nodes
                    .insert(1 + x + 5 * y, DVec3::new(x as f64, y as f64, 0.));
            }
        }
        for y in 0..2 {
            for x in 0..4 {
                let a = 1 + x + 5 * y;
                mesh.elements.push(ElementData {
                    id: mesh.elements.len() as u32 + 1,
                    elem_type: 44,
                    stiff_id: 1,
                    nodes: vec![a, a + 1, a + 6, a + 5],
                });
            }
        }
        mesh.elements.extend([
            ElementData {
                id: 9,
                elem_type: 10,
                stiff_id: 10,
                nodes: vec![7, 8],
            },
            ElementData {
                id: 10,
                elem_type: 10,
                stiff_id: 20,
                nodes: vec![8, 9],
            },
        ]);
        let mut f = frame(&mesh, DVec3::Z);
        for (i, &n) in f.node_ids.iter().enumerate() {
            if (7..=9).contains(&n) {
                f.candidate_points[i][2] = 0.00002;
            }
            if n == 8 {
                f.candidate_points[i][0] += 0.02;
                f.candidate_points[i][1] += 0.002;
            }
        }
        let r = assembly::assemble(&mesh, &f, &policy(1.)).unwrap();
        assert!(
            r.axis_assembly.all_axes_built,
            "{:?}",
            r.axis_assembly.issues
        );
        assert_eq!(r.axis_assembly.axes.len(), 1);
        let axis = &r.axis_assembly.axes[0];
        let anchor = axis.anchors.iter().find(|a| a.source_node == 8).unwrap();
        let point = DVec3::from_array(r.preview.vertices[anchor.vertex]);
        assert!(point.distance(DVec3::new(2.02, 1., 0.)) < 1e-7);
        assert!((anchor.t - 0.51).abs() < 1e-7);
        assert!((axis.spans[0].end_t - anchor.t).abs() < 1e-7);
        assert!((axis.spans[1].start_t - anchor.t).abs() < 1e-7);
        assert_eq!(axis.spans[0].stiffness, 10);
        assert_eq!(axis.spans[1].stiffness, 20);
    }

    #[test]
    fn fixed_middle_joint_accepts_bent_chain_as_constructive_segments() {
        let mesh = slab_beam_column();
        let mut f = frame(&mesh, DVec3::Z);
        let n = f.node_ids.iter().position(|&n| n == 5).unwrap();
        f.candidate_points[n][1] += 0.002;
        let before = f.candidate_points.clone();
        let r = assembly::assemble(&mesh, &f, &policy(1.)).unwrap();
        assert_eq!(f.candidate_points, before);
        assert!(r.all_surface_patches_built);
        assert!(
            r.axis_assembly.all_axes_built,
            "{:?}",
            r.axis_assembly.issues
        );
        assert_eq!(r.axis_assembly.axes.len(), 3);
        assert!(r.axis_assembly.issues.is_empty());
        assert_eq!(
            r.axis_assembly
                .axes
                .iter()
                .flat_map(|a| a.spans.iter().map(|s| s.element))
                .collect::<BTreeSet<_>>(),
            BTreeSet::from([5, 6, 7])
        );
        let vertex = r.vertex_source_nodes.iter().position(|&n| n == 5).unwrap();
        assert!(
            DVec3::from_array(r.preview.vertices[vertex]).distance(DVec3::new(1., 1.002, 0.))
                < 1e-7
        );
    }

    #[test]
    fn fixed_branch_uses_exact_endpoint_chord_before_preferred_direction() {
        let mesh = MeshData {
            nodes: BTreeMap::from([(1, DVec3::ZERO), (2, DVec3::new(1e-7, 0., -4.))])
                .into_iter()
                .collect(),
            elements: vec![ElementData {
                id: 1,
                elem_type: 10,
                stiff_id: 5,
                nodes: vec![1, 2],
            }],
        };
        let axis = frame::Axis {
            constructive_segment: false,
            endpoints: [0, 1],
            anchors: vec![
                frame::Anchor { node: 0, t: 0. },
                frame::Anchor { node: 1, t: 1. },
            ],
            spans: vec![SourceSpan {
                element: 1,
                stiffness: 5,
                start_t: 0.,
                end_t: 1.,
            }],
        };
        let source = frame::Report {
            sliding_parameters: None,
            virtual_incidences: Default::default(),
            budget_cache: Default::default(),
            regularized_directions: false,
            nonlinear_steps: 0,
            candidate_parameters_valid: true,
            policy: frame::Policy {
                up: DVec3::Z.to_array(),
                angle: 0.02,
                maximum_movement: 0.15,
                relative_movement: 0.05,
                minimum_length: 0.03,
                residual_tolerance: 1e-7,
                iterations: 10,
                panel_tolerance: 0.,
                geotechnical: false,
                over_constrained_panels: false,
            },
            accepted: true,
            candidate_constraints_satisfied: true,
            violating_equations: 0,
            largest_constraint_failures: vec![],
            movement_failures: vec![],
            axis_failures: vec![],
            reason: "test".into(),
            iterations: 0,
            candidate_max_residual: 0.,
            candidate_maximum_movement: 0.,
            candidate_over_budget_node_ids: vec![],
            maximum_movement: 0.15,
            node_ids: vec![1, 2],
            reference_points: vec![[0., 0., 0.], [1e-7, 0., -4.]],
            candidate_points: vec![[0., 0., 0.], [1e-7, 0., -4.]],
            candidate_planes: vec![],
            points: vec![[0., 0., 0.], [1e-7, 0., -4.]],
            axes: vec![axis.clone()],
            surfaces: vec![],
            plane_families: vec![],
            equation_count: 0,
            short_axis_indices: vec![],
            short_axis_movement: 0.,
        };
        let locked = BTreeMap::from([(1, DVec3::ZERO), (2, DVec3::new(1e-7, 0., -4.))]);
        let proposal = propose(
            &mesh,
            &source,
            &axis,
            &locked,
            &BTreeMap::new(),
            &BTreeMap::new(),
            &[],
            &policy(1.),
        )
        .expect("the exact fixed endpoint chord is a valid branch");
        assert_eq!(proposal.points[&1], DVec3::ZERO);
        assert_eq!(proposal.points[&2], DVec3::new(1e-7, 0., -4.));
    }

    #[test]
    fn short_bar_may_shorten_within_its_end_budgets() {
        // A 10 mm bar whose ends are fixed by surface planes 10 um closer
        // (geotechnical frame): each end may move by the plane distance, so
        // the bar is built; ends that coincide or swap are still refused.
        let ends = [DVec3::new(2., -1., 3.), DVec3::new(2.01, -1., 3.)];
        let mesh = MeshData {
            nodes: BTreeMap::from([(1, ends[0]), (2, ends[1])])
                .into_iter()
                .collect(),
            elements: vec![ElementData {
                id: 1,
                elem_type: 10,
                stiff_id: 5,
                nodes: vec![1, 2],
            }],
        };
        let axis = frame::Axis {
            constructive_segment: false,
            endpoints: [0, 1],
            anchors: vec![
                frame::Anchor { node: 0, t: 0. },
                frame::Anchor { node: 1, t: 1. },
            ],
            spans: vec![SourceSpan {
                element: 1,
                stiffness: 5,
                start_t: 0.,
                end_t: 1.,
            }],
        };
        let mut f = frame(&mesh, DVec3::Z);
        f.axes = vec![axis.clone()];
        f.short_axis_movement = 0.01;
        f.policy.geotechnical = true;
        let run = |b: DVec3| {
            // The geotechnical frame keeps the length class of a short axis,
            // not its exact length.
            let mut f = f.clone();
            f.candidate_points = vec![ends[0].to_array(), b.to_array()];
            let locked = BTreeMap::from([(1, ends[0]), (2, b)]);
            propose(
                &mesh,
                &f,
                &axis,
                &locked,
                &BTreeMap::new(),
                &BTreeMap::new(),
                &[],
                &policy(1.),
            )
            .map(|_| ())
            .map_err(|(reason, _)| reason)
        };
        assert_eq!(run(DVec3::new(2.00999, -1., 3.)), Ok(()));
        assert!(run(ends[0]).is_err());
        assert!(run(DVec3::new(1.99, -1., 3.)).is_err());
    }

    #[test]
    fn coplanar_intervals_respect_holes_and_boundary_roles() {
        let contours = vec![
            vec![[0., 0.], [4., 0.], [4., 4.], [0., 4.]],
            vec![[1., 1.], [3., 1.], [3., 3.], [1., 3.]],
        ];
        assert_eq!(
            intervals([0., 2.], [4., 2.], &contours, 1e-8),
            vec![
                (0., 0.25, Location::Interior),
                (0.75, 1., Location::Interior)
            ]
        );
        assert_eq!(
            intervals([0., 0.], [4., 0.], &contours, 1e-8),
            vec![(0., 1., Location::Boundary)]
        );
        assert_eq!(intervals([1.1, 2.], [2.9, 2.], &contours, 1e-8), vec![]);
        assert_eq!(location([2., 2.], &contours, 1e-8), None);
    }
    #[test]
    fn short_axis_preserves_vector_and_over_budget_attempt_creates_no_vertices() {
        let mut mesh = MeshData::default();
        mesh.nodes.insert(1, DVec3::ZERO);
        mesh.nodes.insert(2, DVec3::new(0.006, 0.008, 0.));
        mesh.elements.push(ElementData {
            id: 1,
            elem_type: 10,
            stiff_id: 9,
            nodes: vec![1, 2],
        });
        let mut f = frame(&mesh, DVec3::Z);
        for p in &mut f.candidate_points {
            p[2] += 0.0001;
        }
        let r = assembly::assemble(&mesh, &f, &policy(1.)).unwrap();
        assert!(
            r.axis_assembly.all_axes_built,
            "{:?}",
            r.axis_assembly.issues
        );
        let ends = r.axis_assembly.axes[0].endpoints;
        let d = DVec3::from_array(r.preview.vertices[ends[1]])
            - DVec3::from_array(r.preview.vertices[ends[0]]);
        assert!(d.distance(mesh.nodes[&2] - mesh.nodes[&1]) < 1e-8);
        for p in &mut f.candidate_points {
            p[2] += 0.01;
        }
        let refused = assembly::assemble(&mesh, &f, &policy(1.)).unwrap();
        assert!(!refused.axis_assembly.all_axes_built);
        assert!(refused.preview.vertices.is_empty());
        assert!(refused.axis_assembly.axes.is_empty());
        assert_eq!(
            refused.axis_assembly.issues[0].reason,
            "axis_movement_budget"
        );
    }
}
