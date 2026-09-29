//! Closure of small gaps between structures, for PLAXIS.
//!
//! PLAXIS 3D intersects imported geometry itself; a gap or offset far below
//! the target element size (a wall top 25 mm below a slab, a slab edge 30 mm
//! short of a wall) either breaks the intersection or forces a very fine
//! local mesh. A vertex of one surface closer than the tolerance to another
//! surface, without being one of its vertices, is closed onto it:
//!
//! 1. it moves onto that surface's plane, keeping every plane it lies on;
//! 2. if it is still outside that surface, it merges into a nearby contour
//!    vertex of it, or moves onto the nearest contour edge (which is split),
//!    or, if it cannot move, that edge bends through it.
//!
//! Every step is transactional and validated (planarity, valid contours,
//! bars never bent), and every closure is reported with its source node.
//! By default only gaps within one plane are closed: the vertex lies in the
//! other surface's plane (within the minimum edge) and only its distance to
//! that surface's contour is closed. An offset across the plane (a wall top
//! 25 mm below a slab) is closed too (`offsets`, the default), except
//! between parallel structures: a vertex of a surface parallel to the other
//! one (two slabs at different levels) never moves across to it.
//! Junction insertion afterwards represents the new contacts.
use super::cleanup::{self, Bars};
use super::intersection;
use crate::reconstruction::{closed_contains, Model, PlaneFrame};
use glam::{DVec2, DVec3};
use serde::Serialize;
use std::collections::BTreeSet;

#[derive(Debug, Serialize)]
pub struct Closure {
    /// "settled", "merged", "onto_edge" or "bent_edge".
    pub kind: String,
    pub vertex: usize,
    pub source_node: Option<u32>,
    pub surface: usize,
    pub gap: f64,
    /// Movement of the vertex (or of the merged vertex).
    pub movement: f64,
}

#[derive(Debug, Serialize)]
pub struct Rejected {
    pub vertex: usize,
    pub source_node: Option<u32>,
    pub surface: usize,
    pub gap: f64,
    pub reason: String,
}

#[derive(Debug, Default, Serialize)]
pub struct Report {
    pub tolerance: f64,
    pub closed: Vec<Closure>,
    pub rejected: Vec<Rejected>,
}

struct Gap {
    distance: f64,
    height: f64,
    outside: f64,
}

fn surface_vertices(model: &Model, s: usize) -> BTreeSet<usize> {
    model
        .surface_edges(s)
        .flat_map(|e| model.edges[e])
        .collect()
}

/// Distance from a plane point to the closed material of a surface.
fn outside(model: &Model, s: usize, uv: [f64; 2]) -> f64 {
    let contours = &model.surfaces[s].contours;
    if closed_contains(contours, uv, model.precision) {
        return 0.;
    }
    let p = DVec2::from_array(uv);
    contours
        .iter()
        .flat_map(|ring| {
            (0..ring.len()).map(move |i| {
                let a = DVec2::from_array(ring[i]);
                let d = DVec2::from_array(ring[(i + 1) % ring.len()]) - a;
                let t = ((p - a).dot(d) / d.length_squared()).clamp(0., 1.);
                p.distance(a + d * t)
            })
        })
        .fold(f64::INFINITY, f64::min)
}

fn gap(model: &Model, v: usize, s: usize, members: &BTreeSet<usize>) -> Option<Gap> {
    if members.contains(&v) {
        return None;
    }
    let plane = &model.planes[model.surfaces[s].plane];
    let p = model.vertices[v];
    let height = plane.distance(p);
    let outside = outside(model, s, plane.project(p));
    let distance = height.hypot(outside);
    Some(Gap {
        distance,
        height,
        outside,
    })
}

fn bounds(model: &Model, s: usize) -> (DVec3, DVec3) {
    model
        .surface_edges(s)
        .flat_map(|e| model.edges[e])
        .map(|v| DVec3::from_array(model.vertices[v]))
        .fold(
            (DVec3::splat(f64::INFINITY), DVec3::splat(f64::NEG_INFINITY)),
            |(lo, hi), p| (lo.min(p), hi.max(p)),
        )
}

/// Planes of every surface using `v`.
fn planes_of(model: &Model, v: usize) -> Vec<PlaneFrame> {
    cleanup::users(model, v)
        .into_iter()
        .map(|s| model.planes[model.surfaces[s].plane].clone())
        .collect()
}

/// Close the gap between vertex `v` and surface `s` on a trial model.
fn close_one(
    model: &Model,
    bars: &mut Bars<'_>,
    v: usize,
    s: usize,
    tolerance: f64,
    fixed: &BTreeSet<usize>,
) -> Result<(Model, String, f64), String> {
    let mut trial = model.clone();
    let from = DVec3::from_array(model.vertices[v]);
    let plane = model.planes[model.surfaces[s].plane].clone();
    let members = surface_vertices(model, s);
    // 1. Onto the plane, keeping every own plane.
    if plane.distance(from.to_array()).abs() > model.precision {
        if fixed.contains(&v) {
            return Err("retained_node".into());
        }
        let mut planes = planes_of(model, v);
        planes.push(plane.clone());
        let refs: Vec<&PlaneFrame> = planes.iter().collect();
        let target = intersection(from, &refs, model.precision).ok_or("inconsistent_planes")?;
        if target.distance(from) > tolerance {
            return Err(format!(
                "movement_beyond_tolerance: {:e}",
                target.distance(from)
            ));
        }
        cleanup::move_with_axes(&mut trial, bars.axes, v, target, tolerance)?;
    }
    let p = DVec3::from_array(trial.vertices[v]);
    // Settled onto a vertex of the surface: one vertex.
    let coincident = members
        .iter()
        .map(|&w| (DVec3::from_array(trial.vertices[w]).distance(p), w))
        .filter(|&(d, _)| d < trial.minimum_edge)
        .min_by(|x, y| x.0.total_cmp(&y.0).then(x.1.cmp(&y.1)));
    if let Some((_, w)) = coincident {
        let mut merged = trial.clone();
        if cleanup::merge(&mut merged, bars, v, w, tolerance, fixed).is_ok() {
            let reach = DVec3::from_array(merged.vertices[w]).distance(from);
            return Ok((merged, "merged".into(), reach));
        }
    }
    if outside(&trial, s, plane.project(p.to_array())) <= trial.precision {
        let movement = p.distance(from);
        return Ok((trial, "settled".into(), movement));
    }
    // 2a. Into a nearby contour vertex of the surface.
    let nearest = members
        .iter()
        .map(|&w| (DVec3::from_array(trial.vertices[w]).distance(p), w))
        .filter(|&(d, _)| d <= tolerance)
        .min_by(|x, y| x.0.total_cmp(&y.0).then(x.1.cmp(&y.1)));
    if let Some((_, w)) = nearest {
        let mut merged = trial.clone();
        if let Ok(movement) = cleanup::merge(&mut merged, bars, v, w, tolerance, fixed) {
            let reach = DVec3::from_array(merged.vertices[w]).distance(from);
            return Ok((merged, "merged".into(), movement.max(reach)));
        }
    }
    // 2b. Onto a nearby contour edge: the vertex slides along its own planes
    //     and the surface plane to where they cross the edge line, and the
    //     edge is split there. Else the edge bends through the vertex.
    let own = planes_of(&trial, v);
    let normal = DVec3::from_array(plane.normal);
    let mut best: Option<(f64, usize, DVec3)> = None;
    let mut nearest_edge: Option<(f64, usize)> = None;
    for e in trial.surfaces[s]
        .boundaries
        .iter()
        .flatten()
        .map(|u| u.edge)
    {
        let [a, b] = trial.edges[e];
        let (pa, pb) = (
            DVec3::from_array(trial.vertices[a]),
            DVec3::from_array(trial.vertices[b]),
        );
        let d = pb - pa;
        let t = ((p - pa).dot(d) / d.length_squared()).clamp(0., 1.);
        let near = p.distance(pa + d * t);
        if near > tolerance {
            continue;
        }
        if t > 0. && t < 1. && nearest_edge.is_none_or(|x| (near, e) < x) {
            nearest_edge = Some((near, e));
        }
        let Ok(side) = PlaneFrame::new(pa.to_array(), normal.cross(d).to_array()) else {
            continue;
        };
        let mut planes = own.clone();
        planes.push(plane.clone());
        planes.push(side);
        let refs: Vec<&PlaneFrame> = planes.iter().collect();
        let Some(c) = intersection(p, &refs, trial.precision) else {
            continue;
        };
        let u = (c - pa).dot(d) / d.length_squared();
        let movement = c.distance(from);
        if u > 0.
            && u < 1.
            && movement <= tolerance
            && best.is_none_or(|x| (movement, e) < (x.0, x.1))
        {
            best = Some((movement, e, c));
        }
    }
    if let (Some((movement, e, c)), false) = (best, fixed.contains(&v)) {
        let mut moved = trial.clone();
        if cleanup::move_with_axes(&mut moved, bars.axes, v, c, tolerance).is_ok()
            && moved.split_edge(e, v).is_ok()
        {
            return Ok((moved, "onto_edge".into(), movement));
        }
    }
    let Some((_, e)) = nearest_edge else {
        return Err("no_contour_edge_within_tolerance".into());
    };
    trial
        .split_edge_within(e, v, tolerance)
        .map_err(|error| format!("bend_{error:?}"))?;
    let movement = p.distance(from);
    Ok((trial, "bent_edge".into(), movement))
}

/// Close gaps narrower than `tolerance` between surface vertices (and bar
/// nodes) and other surfaces.
pub fn close(
    model: &mut Model,
    bars: &mut Bars<'_>,
    tolerance: f64,
    offsets: bool,
    fixed: &BTreeSet<usize>,
    source_nodes: &[u32],
) -> Report {
    // Planes parallel within the plane recognition angle (0.02 rad).
    let parallel = |model: &Model, v: usize, s: usize| {
        let n = DVec3::from_array(model.planes[model.surfaces[s].plane].normal);
        cleanup::users(model, v).into_iter().any(|u| {
            DVec3::from_array(model.planes[model.surfaces[u].plane].normal)
                .cross(n)
                .length()
                < 0.02
        })
    };
    let eligible = |model: &Model, v: usize, s: usize, g: &Gap| {
        g.distance > model.precision * 10.
            && g.distance < tolerance
            && (g.height.abs() <= model.minimum_edge || (offsets && !parallel(model, v, s)))
    };
    let mut report = Report {
        tolerance,
        ..Default::default()
    };
    let mut tried = BTreeSet::new();
    for _ in 0..4 {
        let mut candidates: BTreeSet<usize> = (0..model.surfaces.len())
            .flat_map(|s| surface_vertices(model, s))
            .collect();
        candidates.extend(
            bars.axes
                .iter()
                .flat_map(|a| a.anchors.iter().map(|x| x.vertex)),
        );
        let mut gaps = vec![];
        for s in 0..model.surfaces.len() {
            let (lo, hi) = bounds(model, s);
            let members = surface_vertices(model, s);
            for &v in &candidates {
                let p = DVec3::from_array(model.vertices[v]);
                if (p + tolerance).cmplt(lo).any() || (p - tolerance).cmpgt(hi).any() {
                    continue;
                }
                if let Some(g) = gap(model, v, s, &members) {
                    // Point touches (a corner on another surface's contour)
                    // are closed here too: junction insertion only imprints
                    // intersection lines.
                    if eligible(model, v, s, &g) {
                        gaps.push((g.distance, v, s));
                    }
                }
            }
        }
        gaps.sort_by(|x, y| x.0.total_cmp(&y.0).then(x.1.cmp(&y.1)).then(x.2.cmp(&y.2)));
        let mut changed = false;
        for (_, v, s) in gaps {
            if !tried.insert((v, s)) {
                continue;
            }
            // Earlier closures may have changed or closed this gap.
            let members = surface_vertices(model, s);
            let Some(g) = gap(model, v, s, &members) else {
                continue;
            };
            if !eligible(model, v, s, &g) {
                continue;
            }
            let _ = (g.height, g.outside);
            match close_one(model, bars, v, s, tolerance, fixed) {
                Ok((trial, kind, movement)) => {
                    *model = trial;
                    changed = true;
                    report.closed.push(Closure {
                        kind,
                        vertex: v,
                        source_node: source_nodes.get(v).copied(),
                        surface: s,
                        gap: g.distance,
                        movement,
                    });
                }
                Err(reason) => report.rejected.push(Rejected {
                    vertex: v,
                    source_node: source_nodes.get(v).copied(),
                    surface: s,
                    gap: g.distance,
                    reason,
                }),
            }
        }
        if !changed {
            break;
        }
    }
    report
}

#[cfg(test)]
mod tests {
    use super::super::junctions::tests::{build, run, slab, wall, Placement};
    use super::*;

    fn no_bars() -> Bars<'static> {
        Bars {
            axes: Box::leak(Box::new(vec![])),
            contacts: Box::leak(Box::new(vec![])),
        }
    }

    #[test]
    fn wall_top_below_a_slab_is_settled_and_joined() {
        for place in Placement::all() {
            // The wall top is 25 mm below the slab.
            let mut m = build(&place, &[slab(0., 4.), wall(1., 3., -2., -0.025)]);
            // Across the slab plane: kept unless offsets are requested.
            let mut kept = m.clone();
            let r = close(
                &mut kept,
                &mut no_bars(),
                0.05 * place.scale,
                false,
                &BTreeSet::new(),
                &[],
            );
            assert!(r.closed.is_empty());
            let r = close(
                &mut m,
                &mut no_bars(),
                0.05 * place.scale,
                true,
                &BTreeSet::new(),
                &[],
            );
            assert_eq!(r.closed.len(), 2, "{:?}", r.rejected);
            assert!(r.closed.iter().all(|c| c.kind == "settled"));
            for c in &r.closed {
                assert!((c.movement - 0.025 * place.scale).abs() < 1e-9 * place.scale);
            }
            let j = run(&mut m);
            assert!(j.issues.is_empty(), "{:?}", j.issues);
            assert_eq!(j.junctions.len(), 1);
            // A 0.1 m gap, or a zero tolerance (a real joint), stays.
            let mut m = build(&place, &[slab(0., 4.), wall(1., 3., -2., -0.1)]);
            assert!(close(
                &mut m,
                &mut no_bars(),
                0.05 * place.scale,
                true,
                &BTreeSet::new(),
                &[]
            )
            .closed
            .is_empty());
        }
    }

    #[test]
    fn slab_edge_short_of_a_wall_reaches_it() {
        for place in Placement::all() {
            // A wall in the plane x = 2.03 below a slab ending at x = 2.
            let cross = (
                vec![vec![
                    [2.03, 0., -2.],
                    [2.03, 4., -2.],
                    [2.03, 4., 0.],
                    [2.03, 0., 0.],
                ]],
                [1., 0., 0.],
            );
            let mut m = build(&place, &[slab(0., 2.), cross]);
            // The wall top lies in the slab plane: a gap within one plane.
            let r = close(
                &mut m,
                &mut no_bars(),
                0.05 * place.scale,
                false,
                &BTreeSet::new(),
                &[],
            );
            assert_eq!(r.closed.len(), 2, "{:?}", r.rejected);
            // The slab edge now lies on the wall top: one shared edge.
            let shared = m
                .surface_edges(0)
                .filter(|&e| m.surface_edges(1).any(|f| f == e))
                .count();
            assert_eq!(shared, 1);
            let j = run(&mut m);
            assert!(j.issues.is_empty(), "{:?}", j.issues);
        }
    }

    #[test]
    fn slabs_at_different_levels_are_never_brought_together() {
        for place in Placement::all() {
            // Two parallel slabs 25 mm apart in height, overlapping in plan.
            let upper = (
                vec![vec![
                    [1.9, 0., 0.025],
                    [4., 0., 0.025],
                    [4., 4., 0.025],
                    [1.9, 4., 0.025],
                ]],
                [0., 0., 1.],
            );
            let mut m = build(&place, &[slab(0., 2.), upper]);
            let before = m.vertices.clone();
            let r = close(
                &mut m,
                &mut no_bars(),
                0.05 * place.scale,
                true,
                &BTreeSet::new(),
                &[],
            );
            assert!(r.closed.is_empty(), "{:?}", r.closed);
            assert_eq!(m.vertices, before);
        }
    }
}
