//! Straight-axis and surface constraints; these do not assign mechanical ties.
use crate::config::ReconstructionConfig;
use crate::geometry::utils::get_plane_basis;
use crate::models::{MacroBar, MacroPanel, MeshData};
use crate::reconstructors::bars::classify_bar;
use geo::{Intersects, LineString, Point, Polygon};
use glam::DVec3;
use hashbrown::HashMap;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
const EPS: f64 = 1e-8;
const AXIS_RESIDUAL_TOL: f64 = 1e-7;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RejectedSurfaceCrossing {
    pub panel_id: u32,
    pub source_element_ids: Vec<u32>,
    pub point: [f64; 3],
    pub t: f64,
    pub reason: String,
}

/// Measurements of the final geometry, not a claim about mechanical contact.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MissingSurfaceContact {
    pub node_id: u32,
    pub panel_id: u32,
    pub inferred: bool,
    pub source_bar_element_ids: Vec<u32>,
    pub anchors: Vec<SurfaceAnchorMeasurement>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SurfaceAnchorMeasurement {
    pub point: [f64; 3],
    pub plane_distance: f64,
    pub projection_inside: bool,
    pub boundary_distance: f64,
    pub nearest_constraint_distance: Option<f64>,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct CoupledPlaneSearch {
    pub iterations: usize,
    pub requested_incidence_count: usize,
    pub retained_incidence_count: usize,
    pub relaxed_incidence_pairs: Vec<[u32; 2]>,
    pub frozen_panel_ids: Vec<u32>,
    pub final_axis_residual: f64,
    pub final_plane_residual: f64,
    pub failure_reason: Option<String>,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct RotationTrial {
    pub proposed_panels: usize,
    pub attempted_candidates: usize,
    pub accepted: bool,
    pub reason: String,
    pub baseline_missing_incidences: usize,
    pub candidate_missing_incidences: Option<usize>,
    pub max_rotation_radians: f64,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct BarContactSummary {
    pub fitted_axes: usize,
    /// Geometric hypotheses, separate from source FE incidence.
    #[serde(default)]
    pub inferred_panel_ids_by_node: BTreeMap<u32, Vec<u32>>,
    #[serde(default)]
    pub ambiguous_panel_ids_by_node: BTreeMap<u32, Vec<u32>>,
    pub rotation_trial: RotationTrial,
    pub joint_solve_converged: bool,
    pub coupled_planes_converged: bool,
    pub coupled_search: CoupledPlaneSearch,
    pub max_panel_displacement: f64,
    pub max_axis_fit_displacement: f64,
    pub anchor_groups: usize,
    pub surface_crossings: usize,
    pub rejected_surface_crossings: usize,
    pub rejected_crossing_details: Vec<RejectedSurfaceCrossing>,
    pub attached_groups: usize,
    pub partial_groups: usize,
    pub unresolved_panel_ids_by_node: BTreeMap<u32, Vec<u32>>,
    #[serde(default)]
    pub missing_surface_contacts: Vec<MissingSurfaceContact>,
    pub moved_groups: usize,
    pub rejected_groups: usize,
    pub rejected_node_ids: Vec<u32>,
    pub rejection_reasons: BTreeMap<String, Vec<u32>>,
    pub max_displacement: f64,
    pub unresolved_bar_junctions: Vec<u32>,
    pub numerical_residual_nodes: Vec<u32>,
    pub max_joint_gap: f64,
    pub joint_residual_tolerance: f64,
}

fn reject(stats: &mut BarContactSummary, node: u32, reason: &str) {
    stats.rejected_groups += 1;
    stats.rejected_node_ids.push(node);
    stats
        .rejection_reasons
        .entry(reason.into())
        .or_default()
        .push(node);
}

fn contains(panel: &MacroPanel, p: DVec3) -> bool {
    let normal = DVec3::from_array(panel.plane_normal);
    if (normal.dot(p) + panel.plane_d).abs() > EPS {
        return false;
    }
    let (u, v) = get_plane_basis(normal);
    let ring = |points: &Vec<[f64; 3]>| {
        let mut points: Vec<_> = points
            .iter()
            .map(|q| {
                let q = DVec3::from_array(*q);
                (q.dot(u), q.dot(v))
            })
            .collect();
        points.push(points[0]);
        LineString::from(points)
    };
    let polygon = Polygon::new(
        ring(&panel.polygons[0]),
        panel.polygons[1..].iter().map(ring).collect(),
    );
    if polygon.intersects(&Point::new(p.dot(u), p.dot(v))) {
        return true;
    }
    // Roundoff at an existing contour edge is not a geometric gap.
    nearest_boundary(panel, p).distance(p) <= EPS
}
fn nearest_boundary(panel: &MacroPanel, p: DVec3) -> DVec3 {
    let mut best = p;
    let mut distance = f64::INFINITY;
    for ring in &panel.polygons {
        for i in 0..ring.len() {
            let a = DVec3::from_array(ring[i]);
            let b = DVec3::from_array(ring[(i + 1) % ring.len()]);
            let d = b - a;
            if d.length_squared() <= EPS * EPS {
                continue;
            }
            let t = ((p - a).dot(d) / d.length_squared()).clamp(0.0, 1.0);
            let q = a + t * d;
            if q.distance(p) < distance {
                best = q;
                distance = q.distance(p);
            }
        }
    }
    best
}
#[derive(Clone, Copy)]
enum Anchor {
    Start,
    End,
    Interior(usize),
}
fn anchor_point(bar: &MacroBar, anchor: Anchor) -> DVec3 {
    let a = DVec3::from_array(bar.start_point);
    let b = DVec3::from_array(bar.end_point);
    match anchor {
        Anchor::Start => a,
        Anchor::End => b,
        Anchor::Interior(i) => a + bar.constraints[i].t * (b - a),
    }
}
fn project(
    start: DVec3,
    panels: &[MacroPanel],
    indices: &BTreeSet<usize>,
    bars: &[MacroBar],
    anchors: &[(usize, Anchor)],
) -> Option<DVec3> {
    let mut equations: Vec<_> = indices
        .iter()
        .map(|&i| (DVec3::from_array(panels[i].plane_normal), panels[i].plane_d))
        .collect();
    for &(i, anchor) in anchors {
        if matches!(anchor, Anchor::Interior(_)) {
            let a = DVec3::from_array(bars[i].start_point);
            let b = DVec3::from_array(bars[i].end_point);
            let (u, v) = get_plane_basis((b - a).normalize());
            equations.extend([(u, -u.dot(a)), (v, -v.dot(a))]);
        }
    }
    let mut basis: Vec<(DVec3, f64)> = vec![];
    for (mut n, d) in equations {
        let mut rhs = -d - n.dot(start);
        for &(u, d) in &basis {
            let weight = n.dot(u);
            n -= weight * u;
            rhs -= weight * d;
        }
        let length = n.length();
        if length < EPS {
            if rhs.abs() > EPS {
                return None;
            }
        } else {
            basis.push((n / length, rhs / length));
        }
    }
    let p = start + basis.iter().map(|(n, d)| *n * *d).sum::<DVec3>();
    p.is_finite().then_some(p)
}
fn on_axis(bar: &MacroBar, point: DVec3) -> bool {
    let a = DVec3::from_array(bar.start_point);
    let d = DVec3::from_array(bar.end_point) - a;
    let t = (point - a).dot(d) / d.length_squared();
    t >= -EPS && t <= 1.0 + EPS && point.distance(a + t * d) <= EPS
}
fn update_endpoint(bar: &mut MacroBar, anchor: Anchor, point: DVec3, panel_ids: Vec<u32>) {
    let old_a = DVec3::from_array(bar.start_point);
    let old_d = DVec3::from_array(bar.end_point) - old_a;
    match anchor {
        Anchor::Start => {
            bar.start_point = point.to_array();
            bar.start_panel_ids = panel_ids;
        }
        Anchor::End => {
            bar.end_point = point.to_array();
            bar.end_panel_ids = panel_ids;
        }
        Anchor::Interior(i) => {
            bar.constraints[i].panel_ids = panel_ids;
            return;
        }
    }
    let a = DVec3::from_array(bar.start_point);
    let d = DVec3::from_array(bar.end_point) - a;
    let parameter = |p: DVec3| (p - a).dot(d) / d.length_squared();
    for c in &mut bar.constraints {
        c.t = parameter(old_a + c.t * old_d);
    }
    for span in &mut bar.property_spans {
        span.start_t = if span.start_t == 0.0 {
            0.0
        } else {
            parameter(old_a + span.start_t * old_d)
        };
        span.end_t = if span.end_t == 1.0 {
            1.0
        } else {
            parameter(old_a + span.end_t * old_d)
        };
    }
    bar.length = d.length();
    bar.bar_type = classify_bar(a, a + d);
}

// Alternating projections onto shared-node equations and endpoint displacement
// balls. Axis endpoints and interior incidence fractions can move; no extra
// geometric vertices are introduced.
fn reconcile_axes(
    bars: &mut [MacroBar],
    original: &[[DVec3; 2]],
    config: &ReconstructionConfig,
) -> bool {
    reconcile_axes_and_planes(
        bars,
        original,
        config,
        &mut [],
        &BTreeMap::new(),
        &mut CoupledPlaneSearch::default(),
    )
}
fn reconcile_axes_and_planes(
    bars: &mut [MacroBar],
    original: &[[DVec3; 2]],
    config: &ReconstructionConfig,
    panels: &mut [MacroPanel],
    incident: &BTreeMap<u32, BTreeSet<usize>>,
    search: &mut CoupledPlaneSearch,
) -> bool {
    *search = CoupledPlaneSearch {
        failure_reason: Some("iteration_limit".into()),
        ..Default::default()
    };
    let mut joints: BTreeMap<u32, Vec<(usize, f64)>> = BTreeMap::new();
    for (i, bar) in bars.iter().enumerate() {
        joints.entry(bar.start_node_id).or_default().push((i, 0.0));
        joints.entry(bar.end_node_id).or_default().push((i, 1.0));
        for c in &bar.constraints {
            if let Some(node) = c.source_node_id {
                joints.entry(node).or_default().push((i, c.t));
            }
        }
    }
    // Coupled constraint selection must not depend on LIRA node/element IDs.
    let center = original
        .iter()
        .flat_map(|ends| ends.iter())
        .copied()
        .sum::<DVec3>()
        / (original.len() * 2).max(1) as f64;
    let quantize = |x: f64| (x * 1e6).round() as i64;
    let mut ordered_nodes: Vec<_> = joints.keys().copied().collect();
    if !panels.is_empty() {
        for anchors in joints.values_mut() {
            anchors.sort_by_key(|&(i, _)| {
                let d = original[i][1] - original[i][0];
                let m = (original[i][0] + original[i][1]) * 0.5 - center;
                (
                    -quantize(d.length_squared()),
                    quantize(m.length_squared()),
                    quantize(d.dot(m).abs()),
                )
            });
        }
        ordered_nodes.sort_by_key(|node| {
            let anchors = &joints[node];
            let p = anchors
                .iter()
                .map(|&(i, t)| original[i][0] * (1.0 - t) + original[i][1] * t)
                .sum::<DVec3>()
                / anchors.len() as f64;
            quantize(p.distance_squared(center))
        });
    }
    let original_fractions: BTreeMap<_, _> = joints
        .iter()
        .flat_map(|(&node, anchors)| anchors.iter().map(move |&(i, t)| ((node, i), t)))
        .collect();
    let plane_origins: Vec<_> = panels.iter().map(|p| p.plane_d).collect();
    let mut offsets = plane_origins.clone();
    let mut frozen_planes = BTreeSet::new();
    let leaders: Vec<_> = (0..panels.len())
        .map(|i| {
            (0..i)
                .find(|&j| {
                    DVec3::from_array(panels[i].plane_normal)
                        .distance(DVec3::from_array(panels[j].plane_normal))
                        < EPS
                        && (panels[i].plane_d - panels[j].plane_d).abs() < EPS
                })
                .unwrap_or(i)
        })
        .collect();
    let mut surface_equations = Vec::new();
    for &node in &ordered_nodes {
        let anchors = &joints[&node];
        let (i, t) = anchors[0];
        let p = original[i][0] * (1.0 - t) + original[i][1] * t;
        if let Some(ids) = incident.get(&node) {
            for &j in ids {
                let n = DVec3::from_array(panels[j].plane_normal);
                if (n.dot(p) + panels[j].plane_d).abs() <= config.joint_tol {
                    surface_equations.push((node, j));
                }
            }
        }
    }
    search.requested_incidence_count = surface_equations.len();
    search.retained_incidence_count = surface_equations.len();
    let mut points: Vec<_> = bars
        .iter()
        .flat_map(|b| {
            [
                DVec3::from_array(b.start_point),
                DVec3::from_array(b.end_point),
            ]
        })
        .collect();
    for iteration in 0..if panels.is_empty() { 10000 } else { 30000 } {
        search.iterations = iteration + 1;
        // Slide interior incidences along their whole axes. Keeping FE-derived
        // fractions fixed overconstrains a network after geometric straightening.
        for &node in &ordered_nodes {
            let anchors = joints.get_mut(&node).unwrap();
            if anchors.len() < 2 {
                continue;
            }
            let target = anchors
                .iter()
                .map(|&(i, t)| points[2 * i] * (1.0 - t) + points[2 * i + 1] * t)
                .sum::<DVec3>()
                / anchors.len() as f64;
            for (i, t) in anchors {
                if *t > 0.0 && *t < 1.0 {
                    let a = points[2 * *i];
                    let d = points[2 * *i + 1] - a;
                    if d.length_squared() > EPS * EPS {
                        let old_t = original_fractions[&(node, *i)];
                        let origin = original[*i][0] + old_t * (original[*i][1] - original[*i][0]);
                        let center = (origin - a).dot(d) / d.length_squared();
                        let perpendicular = origin.distance_squared(a + center * d);
                        let radius = ((config.joint_tol.powi(2) - perpendicular).max(0.0)
                            / d.length_squared())
                        .sqrt();
                        let low = (center - radius).max(1e-10);
                        let high = (center + radius).min(1.0 - 1e-10);
                        if low <= high {
                            *t = ((target - a).dot(d) / d.length_squared()).clamp(low, high);
                        }
                    }
                }
            }
        }
        let equations: Vec<Vec<(usize, f64)>> = ordered_nodes
            .iter()
            .map(|node| &joints[node])
            .flat_map(|anchors| {
                anchors
                    .iter()
                    .skip(1)
                    .map(|&(j, u)| {
                        let (i, t) = anchors[0];
                        let mut coefficients = BTreeMap::<usize, f64>::new();
                        for (index, value) in [
                            (2 * i, 1.0 - t),
                            (2 * i + 1, t),
                            (2 * j, -(1.0 - u)),
                            (2 * j + 1, -u),
                        ] {
                            *coefficients.entry(index).or_default() += value;
                        }
                        coefficients
                            .into_iter()
                            .filter(|(_, v)| v.abs() > 1e-14)
                            .collect()
                    })
                    .collect::<Vec<_>>()
            })
            .collect();

        for equation in &equations {
            let norm: f64 = equation.iter().map(|(_, w)| w * w).sum();
            if norm == 0.0 {
                continue;
            }
            let residual: DVec3 = equation.iter().map(|&(i, w)| points[i] * w).sum();
            for &(i, w) in equation {
                points[i] -= residual * (w / norm);
            }
        }
        for &(node, j) in &surface_equations {
            let (i, t) = joints[&node][0];
            let n = DVec3::from_array(panels[j].plane_normal);
            let leader = leaders[j];
            let residual =
                n.dot(points[2 * i] * (1.0 - t) + points[2 * i + 1] * t) + offsets[leader];
            let step = residual
                / ((1.0 - t).powi(2)
                    + t * t
                    + if frozen_planes.contains(&leader) {
                        0.0
                    } else {
                        1.0
                    });
            points[2 * i] -= n * (step * (1.0 - t));
            points[2 * i + 1] -= n * (step * t);
            if !frozen_planes.contains(&leader) {
                offsets[leader] -= step;
            }
        }
        for i in 0..offsets.len() {
            let leader = leaders[i];
            offsets[leader] = offsets[leader].clamp(
                plane_origins[leader] - config.joint_tol,
                plane_origins[leader] + config.joint_tol,
            );
        }
        for i in 0..offsets.len() {
            offsets[i] = offsets[leaders[i]];
        }
        for (i, point) in points.iter_mut().enumerate() {
            let origin = original[i / 2][i % 2];
            let delta = *point - origin;
            if delta.length() > config.joint_tol {
                *point = origin + delta.normalize() * config.joint_tol;
            }
        }
        let error = equations
            .iter()
            .map(|e| {
                e.iter()
                    .map(|&(i, w)| points[i] * w)
                    .sum::<DVec3>()
                    .length()
            })
            .fold(0.0, f64::max);
        let plane_error = surface_equations
            .iter()
            .map(|&(node, j)| {
                let (i, t) = joints[&node][0];
                (DVec3::from_array(panels[j].plane_normal)
                    .dot(points[2 * i] * (1.0 - t) + points[2 * i + 1] * t)
                    + offsets[j])
                    .abs()
            })
            .fold(0.0, f64::max);
        search.final_axis_residual = error;
        search.final_plane_residual = plane_error;
        // Incompatible source incidences must not poison every feasible plane.
        // Relax the largest residuals in bounded rounds; all omitted incidences
        // are checked again by the final contact pass and remain in diagnostics.
        if iteration % 500 == 499
            && iteration < 16000
            && (plane_error > EPS || error > AXIS_RESIDUAL_TOL)
            && surface_equations.len() > 1
        {
            let mut bar_errors = vec![0.0_f64; bars.len()];
            for anchors in joints.values() {
                let (i, t) = anchors[0];
                let a = points[2 * i] * (1.0 - t) + points[2 * i + 1] * t;
                let gap = anchors
                    .iter()
                    .map(|&(b, u)| a.distance(points[2 * b] * (1.0 - u) + points[2 * b + 1] * u))
                    .fold(0.0, f64::max);
                for &(b, _) in anchors {
                    bar_errors[b] = bar_errors[b].max(gap);
                }
            }
            let mut ranked: Vec<_> = surface_equations
                .iter()
                .enumerate()
                .map(|(k, &(node, j))| {
                    let (i, t) = joints[&node][0];
                    let gap = (DVec3::from_array(panels[j].plane_normal)
                        .dot(points[2 * i] * (1.0 - t) + points[2 * i + 1] * t)
                        + offsets[j])
                        .abs();
                    let anchor = points[2 * i] * (1.0 - t) + points[2 * i + 1] * t;
                    let joint_gap = joints[&node]
                        .iter()
                        .map(|&(b, u)| {
                            anchor.distance(points[2 * b] * (1.0 - u) + points[2 * b + 1] * u)
                        })
                        .fold(0.0, f64::max);
                    (k, gap.max(joint_gap).max(bar_errors[i]))
                })
                .collect();
            ranked.sort_by(|a, b| b.1.total_cmp(&a.1));
            let remove: BTreeSet<_> = ranked
                .iter()
                .take((ranked.len() / 10).max(1))
                .filter(|(_, gap)| *gap > EPS)
                .map(|(k, _)| *k)
                .collect();
            for &k in &remove {
                let (node, j) = surface_equations[k];
                search.relaxed_incidence_pairs.push([node, panels[j].id]);
            }
            surface_equations = surface_equations
                .into_iter()
                .enumerate()
                .filter(|(k, _)| !remove.contains(k))
                .map(|(_, e)| e)
                .collect();
            search.retained_incidence_count = surface_equations.len();
            if !remove.is_empty() {
                for (i, point) in points.iter_mut().enumerate() {
                    *point = original[i / 2][i % 2];
                }
                offsets.clone_from(&plane_origins);
                for (&node, anchors) in &mut joints {
                    for (i, t) in anchors {
                        *t = original_fractions[&(node, *i)];
                    }
                }
                continue;
            }
        }
        if error < AXIS_RESIDUAL_TOL && plane_error < EPS * 0.1 {
            if bars.iter().enumerate().any(|(i, _)| {
                let d = points[2 * i + 1] - points[2 * i];
                let old = original[i][1] - original[i][0];
                d.length() + EPS < old.length().min(config.min_edge) || d.dot(old) <= 0.0
            }) {
                search.failure_reason = Some("axis_quality".into());
                return false;
            }
            let mut parameters = Vec::new();
            for (i, bar) in bars.iter().enumerate() {
                let values: Vec<_> = bar
                    .constraints
                    .iter()
                    .map(|c| {
                        c.source_node_id
                            .and_then(|node| {
                                joints[&node].iter().find(|(j, _)| *j == i).map(|(_, t)| *t)
                            })
                            .unwrap_or(c.t)
                    })
                    .collect();
                for (c, &t) in bar.constraints.iter().zip(&values) {
                    let old = original[i][0] + c.t * (original[i][1] - original[i][0]);
                    let new = points[2 * i] + t * (points[2 * i + 1] - points[2 * i]);
                    if t <= 0.0 || t >= 1.0 || new.distance(old) > config.joint_tol + EPS {
                        search.failure_reason = Some("anchor_displacement".into());
                        return false;
                    }
                }
                let mut order: Vec<_> = bar
                    .constraints
                    .iter()
                    .zip(&values)
                    .map(|(c, t)| (c.t, *t))
                    .collect();
                order.push((0.0, 0.0));
                order.push((1.0, 1.0));
                order.sort_by(|a, b| a.0.total_cmp(&b.0));
                for pair in order.windows(2) {
                    let old = (pair[1].0 - pair[0].0) * (original[i][1] - original[i][0]).length();
                    let new =
                        (pair[1].1 - pair[0].1) * (points[2 * i + 1] - points[2 * i]).length();
                    if new + EPS < old.min(config.min_edge) {
                        search.failure_reason = Some("constraint_spacing".into());
                        return false;
                    }
                }
                parameters.push(values);
            }
            let corrected = if panels.is_empty() {
                Vec::new()
            } else {
                match crate::geometry::topology::shifted_planes_checked(panels, &offsets, config) {
                    Ok(result) => result,
                    Err(affected) => {
                        let previous = frozen_planes.len();
                        frozen_planes.extend(affected.iter().map(|&i| leaders[i]));
                        search.frozen_panel_ids = (0..panels.len())
                            .filter(|&i| frozen_planes.contains(&leaders[i]))
                            .map(|i| panels[i].id)
                            .collect();
                        if frozen_planes.len() == previous {
                            search.failure_reason = Some("panel_geometry".into());
                            return false;
                        }
                        for (i, point) in points.iter_mut().enumerate() {
                            *point = original[i / 2][i % 2];
                        }
                        offsets.clone_from(&plane_origins);
                        for (&node, anchors) in &mut joints {
                            for (i, t) in anchors {
                                *t = original_fractions[&(node, *i)];
                            }
                        }
                        continue;
                    }
                }
            };
            for (i, bar) in bars.iter_mut().enumerate() {
                for (c, &t) in bar.constraints.iter_mut().zip(&parameters[i]) {
                    if c.kinds.iter().any(|k| k == "property_boundary") {
                        for span in &mut bar.property_spans {
                            if (span.start_t - c.t).abs() < EPS {
                                span.start_t = t;
                            }
                            if (span.end_t - c.t).abs() < EPS {
                                span.end_t = t;
                            }
                        }
                    }
                    c.t = t;
                }
                bar.start_point = points[2 * i].to_array();
                bar.end_point = points[2 * i + 1].to_array();
                bar.length = points[2 * i].distance(points[2 * i + 1]);
                bar.bar_type = classify_bar(points[2 * i], points[2 * i + 1]);
            }
            panels.clone_from_slice(&corrected);
            search.failure_reason = None;
            return true;
        }
    }
    false
}

pub fn attach_bar_endpoints(
    bars: &mut [MacroBar],
    panels: &mut [MacroPanel],
    mesh: &MeshData,
    canonical: &HashMap<u32, u32>,
    config: &ReconstructionConfig,
) -> BarContactSummary {
    let original_bars = bars.to_vec();
    let original_panels = panels.to_vec();
    let mut baseline = attach_bar_endpoints_impl(bars, panels, mesh, canonical, config);
    let missing = |r: &BarContactSummary| {
        r.unresolved_panel_ids_by_node
            .values()
            .map(|ids| ids.len())
            .sum::<usize>()
    };
    let original_missing = missing(&baseline);
    let proposals = super::plane_fit::proposals(bars, &original_panels, mesh, canonical, config);
    let attempted = proposals.len();
    baseline.rotation_trial = RotationTrial {
        baseline_missing_incidences: original_missing,
        reason: "no_safe_rotation".into(),
        ..Default::default()
    };
    'trials: for (mut candidate_panels, count) in proposals {
        let mut trial = RotationTrial {
            baseline_missing_incidences: original_missing,
            ..Default::default()
        };
        trial.proposed_panels = count;
        let mut candidate_bars = original_bars.clone();
        let mut candidate = attach_bar_endpoints_impl(
            &mut candidate_bars,
            &mut candidate_panels,
            mesh,
            canonical,
            config,
        );
        if candidate.inferred_panel_ids_by_node != baseline.inferred_panel_ids_by_node
            || candidate.ambiguous_panel_ids_by_node != baseline.ambiguous_panel_ids_by_node
        {
            continue;
        }
        trial.candidate_missing_incidences = Some(missing(&candidate));
        let mut max_bar = 0.0_f64;
        let mut max_panel = 0.0_f64;
        for (old, new) in original_bars.iter().zip(&candidate_bars) {
            let a = DVec3::from_array(old.start_point);
            let d = DVec3::from_array(old.end_point) - a;
            let b = DVec3::from_array(new.start_point);
            let e = DVec3::from_array(new.end_point) - b;
            max_bar = max_bar.max(a.distance(b)).max((a + d).distance(b + e));
            for c in &old.constraints {
                if let Some(node) = c.source_node_id {
                    let Some(other) = new
                        .constraints
                        .iter()
                        .find(|c| c.source_node_id == Some(node))
                    else {
                        continue 'trials;
                    };
                    max_bar = max_bar.max((a + c.t * d).distance(b + other.t * e));
                }
            }
        }
        for (old, new) in original_panels.iter().zip(&candidate_panels) {
            for (a, b) in old
                .polygons
                .iter()
                .flatten()
                .zip(new.polygons.iter().flatten())
            {
                max_panel = max_panel.max(DVec3::from_array(*a).distance(DVec3::from_array(*b)));
            }
            trial.max_rotation_radians = trial.max_rotation_radians.max(
                DVec3::from_array(old.plane_normal)
                    .dot(DVec3::from_array(new.plane_normal))
                    .clamp(-1.0, 1.0)
                    .acos(),
            );
        }
        let quality = |r: &BarContactSummary| {
            (
                r.unresolved_bar_junctions.len(),
                missing(r),
                r.rejected_groups,
                r.rejected_surface_crossings,
            )
        };
        if max_bar > config.joint_tol + EPS || max_panel > config.joint_tol + EPS {
            trial.reason = "cumulative_displacement".into();
        } else if quality(&candidate) < quality(&baseline) {
            trial.accepted = true;
            trial.reason = "improved_contacts".into();
            candidate.max_displacement = max_bar;
            candidate.max_panel_displacement = max_panel;
            candidate.rotation_trial = trial;
            bars.clone_from_slice(&candidate_bars);
            panels.clone_from_slice(&candidate_panels);
            baseline = candidate;
            continue;
        } else {
            trial.reason = "no_contact_improvement".into();
        }
        if !baseline.rotation_trial.accepted
            && baseline
                .rotation_trial
                .candidate_missing_incidences
                .is_none_or(|v| v > trial.candidate_missing_incidences.unwrap_or(usize::MAX))
        {
            baseline.rotation_trial = trial;
        }
    }
    baseline.rotation_trial.attempted_candidates = attempted;
    baseline.missing_surface_contacts = measure_missing_contacts(bars, panels, &baseline);
    baseline
}

fn measure_missing_contacts(
    bars: &[MacroBar],
    panels: &[MacroPanel],
    stats: &BarContactSummary,
) -> Vec<MissingSurfaceContact> {
    let mut result = vec![];
    for (&node, ids) in &stats.unresolved_panel_ids_by_node {
        let mut points = vec![];
        let mut sources = BTreeSet::new();
        for bar in bars {
            let mut local = vec![];
            if bar.start_node_id == node {
                local.push(DVec3::from_array(bar.start_point));
            }
            if bar.end_node_id == node {
                local.push(DVec3::from_array(bar.end_point));
            }
            for (i, c) in bar.constraints.iter().enumerate() {
                if c.source_node_id == Some(node) {
                    local.push(anchor_point(bar, Anchor::Interior(i)));
                }
            }
            if !local.is_empty() {
                sources.extend(bar.source_element_ids.iter().copied());
            }
            points.extend(local);
        }
        for id in ids {
            let Some(panel) = panels.iter().find(|p| p.id == *id) else {
                continue;
            };
            let n = DVec3::from_array(panel.plane_normal);
            let anchors = points
                .iter()
                .map(|&point| {
                    let signed = n.dot(point) + panel.plane_d;
                    let projection = point - signed * n;
                    SurfaceAnchorMeasurement {
                        point: point.to_array(),
                        plane_distance: signed.abs(),
                        projection_inside: contains(panel, projection),
                        boundary_distance: nearest_boundary(panel, projection).distance(projection),
                        nearest_constraint_distance: panel
                            .constraint_points
                            .iter()
                            .chain(panel.polygons.iter().flatten())
                            .map(|p| DVec3::from_array(*p).distance(projection))
                            .min_by(f64::total_cmp),
                    }
                })
                .collect();
            result.push(MissingSurfaceContact {
                node_id: node,
                panel_id: *id,
                inferred: stats
                    .inferred_panel_ids_by_node
                    .get(&node)
                    .is_some_and(|ids| ids.contains(id)),
                source_bar_element_ids: sources.iter().copied().collect(),
                anchors,
            });
        }
    }
    result
}

/// Infer missing incidence before the coupled solve. A point must project inside
/// the actual surface (holes excluded). Separated parallel alternatives are not
/// enough evidence to choose a construction connection across a possible joint.
fn infer_incidence(
    bars: &[MacroBar],
    panels: &[MacroPanel],
    incident: &mut BTreeMap<u32, BTreeSet<usize>>,
    config: &ReconstructionConfig,
    stats: &mut BarContactSummary,
) {
    let mut anchors: BTreeMap<u32, Vec<DVec3>> = BTreeMap::new();
    for bar in bars {
        anchors
            .entry(bar.start_node_id)
            .or_default()
            .push(DVec3::from_array(bar.start_point));
        anchors
            .entry(bar.end_node_id)
            .or_default()
            .push(DVec3::from_array(bar.end_point));
        for (i, c) in bar.constraints.iter().enumerate() {
            if let Some(node) = c.source_node_id {
                anchors
                    .entry(node)
                    .or_default()
                    .push(anchor_point(bar, Anchor::Interior(i)));
            }
        }
    }
    for (node, points) in anchors {
        if incident.contains_key(&node) {
            continue;
        }
        let candidates: BTreeSet<_> = panels
            .iter()
            .enumerate()
            .filter_map(|(i, panel)| {
                let n = DVec3::from_array(panel.plane_normal);
                points
                    .iter()
                    .all(|p| {
                        let distance = n.dot(*p) + panel.plane_d;
                        distance.abs() <= config.joint_tol && contains(panel, *p - distance * n)
                    })
                    .then_some(i)
            })
            .collect();
        let ambiguous = candidates.iter().any(|&i| {
            candidates.iter().any(|&j| {
                if i >= j {
                    return false;
                }
                let a = DVec3::from_array(panels[i].plane_normal);
                let b = DVec3::from_array(panels[j].plane_normal);
                a.dot(b).abs() >= config.tol_angle.cos()
                    && points.iter().any(|p| {
                        let pa = *p - a * (a.dot(*p) + panels[i].plane_d);
                        let pb = *p - b * (b.dot(*p) + panels[j].plane_d);
                        pa.distance(pb) > EPS
                    })
            })
        });
        let ids = candidates.iter().map(|&i| panels[i].id).collect();
        if ambiguous {
            stats.ambiguous_panel_ids_by_node.insert(node, ids);
            // An explicit empty entry prevents later nearest-surface fallback.
            incident.insert(node, BTreeSet::new());
        } else if !candidates.is_empty() {
            stats.inferred_panel_ids_by_node.insert(node, ids);
            incident.insert(node, candidates);
        }
    }
}

fn attach_bar_endpoints_impl(
    bars: &mut [MacroBar],
    panels: &mut [MacroPanel],
    mesh: &MeshData,
    canonical: &HashMap<u32, u32>,
    config: &ReconstructionConfig,
) -> BarContactSummary {
    let mut stats = BarContactSummary::default();
    let mut source_panel = HashMap::new();
    for (i, p) in panels.iter().enumerate() {
        for &id in &p.source_element_ids {
            source_panel.insert(id, i);
        }
    }
    let mut incident: BTreeMap<u32, BTreeSet<usize>> = BTreeMap::new();
    for element in &mesh.elements {
        if let Some(&panel) = source_panel.get(&element.id) {
            for id in &element.nodes {
                incident
                    .entry(canonical.get(id).copied().unwrap_or(*id))
                    .or_default()
                    .insert(panel);
            }
        }
    }
    infer_incidence(bars, panels, &mut incident, config, &mut stats);
    let original_ends: Vec<_> = bars
        .iter()
        .map(|b| {
            [
                DVec3::from_array(b.start_point),
                DVec3::from_array(b.end_point),
            ]
        })
        .collect();
    let original_constraints: Vec<Vec<_>> = bars
        .iter()
        .map(|b| {
            b.constraints
                .iter()
                .enumerate()
                .map(|(ci, _)| anchor_point(b, Anchor::Interior(ci)))
                .collect()
        })
        .collect();
    // Move an entire straight axis onto incident near-parallel planes before
    // resolving individual anchors. Interior points travel with that axis.
    for i in 0..bars.len() {
        let a = DVec3::from_array(bars[i].start_point);
        let b = DVec3::from_array(bars[i].end_point);
        let direction = (b - a).normalize();
        let indices: BTreeSet<_> = bars[i]
            .source_node_ids
            .iter()
            .filter_map(|node| incident.get(node))
            .flat_map(|ids| ids.iter().copied())
            .filter(|&j| {
                let panel = &panels[j];
                let n = DVec3::from_array(panel.plane_normal);
                n.dot(direction).abs() < config.tol_angle.sin()
                    && (n.dot(a) + panel.plane_d).abs() <= config.joint_tol
                    && (n.dot(b) + panel.plane_d).abs() <= config.joint_tol
            })
            .collect();
        if indices.is_empty() {
            continue;
        }
        if let (Some(new_a), Some(new_b)) = (
            project(a, panels, &indices, bars, &[]),
            project(b, panels, &indices, bars, &[]),
        ) {
            if new_a.distance(a) <= config.joint_tol
                && new_b.distance(b) <= config.joint_tol
                && new_a.distance(new_b) + EPS >= a.distance(b).min(config.min_edge)
                && (new_b - new_a).dot(b - a) > 0.0
            {
                let movement = new_a.distance(a).max(new_b.distance(b));
                if movement > EPS {
                    stats.fitted_axes += 1;
                    stats.max_axis_fit_displacement = stats.max_axis_fit_displacement.max(movement);
                }
                bars[i].start_point = new_a.to_array();
                bars[i].end_point = new_b.to_array();
                bars[i].length = new_a.distance(new_b);
                bars[i].bar_type = classify_bar(new_a, new_b);
            }
        }
    }
    let panel_origins = panels.to_vec();
    stats.coupled_planes_converged = reconcile_axes_and_planes(
        bars,
        &original_ends,
        config,
        panels,
        &incident,
        &mut stats.coupled_search,
    );
    stats.joint_solve_converged =
        stats.coupled_planes_converged || reconcile_axes(bars, &original_ends, config);
    for (before, after) in panel_origins.iter().zip(panels.iter()) {
        for (a, b) in before
            .polygons
            .iter()
            .flatten()
            .zip(after.polygons.iter().flatten())
        {
            stats.max_panel_displacement = stats
                .max_panel_displacement
                .max(DVec3::from_array(*a).distance(DVec3::from_array(*b)));
        }
    }
    if !stats.joint_solve_converged {
        // A failed global proposal must not leave independently shifted axes.
        for (bar, ends) in bars.iter_mut().zip(&original_ends) {
            bar.start_point = ends[0].to_array();
            bar.end_point = ends[1].to_array();
            bar.length = ends[0].distance(ends[1]);
            bar.bar_type = classify_bar(ends[0], ends[1]);
        }
    }
    let mut groups: BTreeMap<u32, Vec<(usize, Anchor)>> = BTreeMap::new();
    for (i, b) in bars.iter().enumerate() {
        groups
            .entry(b.start_node_id)
            .or_default()
            .push((i, Anchor::Start));
        groups
            .entry(b.end_node_id)
            .or_default()
            .push((i, Anchor::End));
        for (ci, c) in b.constraints.iter().enumerate() {
            if let Some(node) = c.source_node_id {
                groups
                    .entry(node)
                    .or_default()
                    .push((i, Anchor::Interior(ci)));
            }
        }
    }
    stats.anchor_groups = groups.len();
    // Use geometric ordering: changing source IDs must not change which joint
    // claims a nearby mesh constraint first.
    let center = mesh.nodes.values().copied().sum::<DVec3>() / mesh.nodes.len().max(1) as f64;
    let mut groups: Vec<_> = groups.into_iter().collect();
    for (_, anchors) in &mut groups {
        anchors.sort_by(|&(i, a), &(j, b)| {
            matches!(b, Anchor::Interior(_))
                .cmp(&matches!(a, Anchor::Interior(_)))
                .then_with(|| bars[j].length.total_cmp(&bars[i].length))
        });
    }
    let rank = |node: u32, anchors: &Vec<(usize, Anchor)>| {
        let p = mesh
            .nodes
            .get(&node)
            .copied()
            .unwrap_or_else(|| anchor_point(&bars[anchors[0].0], anchors[0].1));
        (p.distance_squared(center) * 1e6).round() as i64
    };
    groups.sort_by_key(|(node, anchors)| rank(*node, anchors));
    for (node, ends) in groups {
        let (index, anchor) = ends[0];
        let origin = anchor_point(&bars[index], anchor);
        let nearby = |i: usize| {
            let p = &panels[i];
            let n = DVec3::from_array(p.plane_normal);
            let q = origin - n * (n.dot(origin) + p.plane_d);
            if contains(p, q) {
                q.distance(origin) <= config.joint_tol
            } else {
                nearest_boundary(p, origin).distance(origin) <= config.joint_tol
            }
        };
        let known = incident.get(&node);
        // Preserve source incidence: an endpoint is not silently transferred to an
        // unrelated nearby surface when its intended surfaces cannot be satisfied.
        let indices: BTreeSet<usize> = if let Some(known) = known {
            known.clone()
        } else {
            (0..panels.len()).filter(|&i| nearby(i)).collect()
        };
        if indices.is_empty() && !ends.iter().any(|(_, a)| matches!(a, Anchor::Interior(_))) {
            continue;
        }
        let safe = |point: DVec3| {
            ends.iter().all(|&(i, anchor)| {
                let bar = &bars[i];
                let old_a = DVec3::from_array(bar.start_point);
                let old_b = DVec3::from_array(bar.end_point);
                let initial = match anchor {
                    Anchor::Start => original_ends[i][0],
                    Anchor::End => original_ends[i][1],
                    Anchor::Interior(ci) => original_constraints[i][ci],
                };
                if point.distance(initial) > config.joint_tol + EPS {
                    return false;
                }
                let (a, b) = match anchor {
                    Anchor::Start => (point, old_b),
                    Anchor::End => (old_a, point),
                    Anchor::Interior(_) => return on_axis(bar, point),
                };
                let d = b - a;
                a.distance(b) + EPS >= old_a.distance(old_b).min(config.min_edge)
                    && d.dot(old_b - old_a) > 0.0
                    && bar.constraints.iter().all(|c| {
                        let p = old_a + c.t * (old_b - old_a);
                        let t = (p - a).dot(d) / d.length_squared();
                        t > 0.0 && t < 1.0 && p.distance(a + t * d) <= EPS
                    })
            })
        };
        let accepts_surface = |i: usize, point: DVec3| {
            let panel = &panels[i];
            let distance = nearest_boundary(panel, point).distance(point);
            contains(panel, point)
                && (distance <= EPS || distance + EPS >= config.min_edge)
                && panel
                    .constraint_points
                    .iter()
                    .chain(panel.polygons.iter().flatten())
                    .all(|p| {
                        let distance = DVec3::from_array(*p).distance(point);
                        distance <= EPS || distance + EPS >= config.min_edge
                    })
        };
        let choose = |indices: &BTreeSet<usize>| -> Result<DVec3, &'static str> {
            if indices.iter().any(|&i| !nearby(i)) {
                return Err("surface_outside_tolerance");
            }
            let projected = project(origin, panels, indices, bars, &ends)
                .ok_or("incompatible_axis_and_planes")?;
            let mut candidates = vec![projected];
            candidates.extend(
                indices
                    .iter()
                    .map(|&i| nearest_boundary(&panels[i], projected)),
            );
            for &i in indices {
                candidates.extend(
                    panels[i]
                        .constraint_points
                        .iter()
                        .chain(panels[i].polygons.iter().flatten())
                        .map(|p| DVec3::from_array(*p)),
                );
            }
            candidates.sort_by(|a, b| a.distance(origin).total_cmp(&b.distance(origin)));
            candidates
                .into_iter()
                .find(|&point| safe(point) && indices.iter().all(|&i| accepts_surface(i, point)))
                .ok_or("would_bend_axis_exceed_budget_or_create_short_constraint")
        };
        let (indices, point) = match choose(&indices) {
            Ok(point) => (indices, point),
            Err(reason) => {
                reject(&mut stats, node, reason);
                // Preserve feasible source incidences even when another source
                // surface cannot be reached. Missing incidences remain explicit.
                let mut partial: Vec<_> = indices
                    .iter()
                    .filter_map(|&i| {
                        choose(&BTreeSet::from([i])).ok().map(|point| {
                            let compatible: BTreeSet<_> = indices
                                .iter()
                                .copied()
                                .filter(|&j| nearby(j) && accepts_surface(j, point))
                                .collect();
                            (compatible, point)
                        })
                    })
                    .collect();
                partial.sort_by(|(a, p), (b, q)| {
                    b.len()
                        .cmp(&a.len())
                        .then_with(|| p.distance(origin).total_cmp(&q.distance(origin)))
                });
                if let Some((compatible, point)) = partial.into_iter().next() {
                    stats.partial_groups += 1;
                    stats.unresolved_panel_ids_by_node.insert(
                        node,
                        indices
                            .difference(&compatible)
                            .map(|&i| panels[i].id)
                            .collect(),
                    );
                    (compatible, point)
                } else {
                    stats
                        .unresolved_panel_ids_by_node
                        .insert(node, indices.iter().map(|&i| panels[i].id).collect());
                    // Failure to match a surface must not prevent a safe bar-to-bar join.
                    if ends.len() < 2 {
                        continue;
                    }
                    match choose(&BTreeSet::new()) {
                        Ok(point) => (BTreeSet::new(), point),
                        Err(_) => continue,
                    }
                }
            }
        };
        let panel_ids: Vec<_> = indices.iter().map(|&i| panels[i].id).collect();
        for (i, anchor) in ends {
            if let Anchor::Interior(ci) = anchor {
                let a = DVec3::from_array(bars[i].start_point);
                let d = DVec3::from_array(bars[i].end_point) - a;
                let old_t = bars[i].constraints[ci].t;
                let new_t = (point - a).dot(d) / d.length_squared();
                if bars[i].constraints[ci]
                    .kinds
                    .iter()
                    .any(|k| k == "property_boundary")
                {
                    for span in &mut bars[i].property_spans {
                        if (span.start_t - old_t).abs() < EPS {
                            span.start_t = new_t;
                        }
                        if (span.end_t - old_t).abs() < EPS {
                            span.end_t = new_t;
                        }
                    }
                }
                bars[i].constraints[ci].t = new_t;
            }
            update_endpoint(&mut bars[i], anchor, point, panel_ids.clone());
        }
        let attached = !indices.is_empty();
        for i in indices {
            if panels[i]
                .constraint_points
                .iter()
                .all(|p| DVec3::from_array(*p).distance(point) > EPS)
            {
                panels[i].constraint_points.push(point.to_array());
            }
        }
        if attached {
            stats.attached_groups += 1;
        }
        let movement = origin.distance(point);
        stats.max_displacement = stats.max_displacement.max(movement);
        if movement > EPS {
            stats.moved_groups += 1;
        }
    }
    for (i, bar) in bars.iter().enumerate() {
        stats.max_displacement = stats
            .max_displacement
            .max(DVec3::from_array(bar.start_point).distance(original_ends[i][0]))
            .max(DVec3::from_array(bar.end_point).distance(original_ends[i][1]));
        for (ci, _) in bar.constraints.iter().enumerate() {
            stats.max_displacement = stats
                .max_displacement
                .max(anchor_point(bar, Anchor::Interior(ci)).distance(original_constraints[i][ci]));
        }
    }
    // A through-column remains one LINE. Surface crossings are parametric mesh
    // constraints and do not add geometric endpoints or depend on source FE nodes.
    for bar in bars.iter_mut() {
        let a = DVec3::from_array(bar.start_point);
        let d = DVec3::from_array(bar.end_point) - a;
        for panel in panels.iter_mut() {
            let n = DVec3::from_array(panel.plane_normal);
            let denominator = n.dot(d);
            if denominator.abs() < EPS {
                continue;
            }
            let t = -(n.dot(a) + panel.plane_d) / denominator;
            if t <= EPS || t >= 1.0 - EPS {
                continue;
            }
            let point = a + t * d;
            if !contains(panel, point) {
                continue;
            }
            let boundary_distance = nearest_boundary(panel, point).distance(point);
            let spacing_ok = (boundary_distance <= EPS
                || boundary_distance + EPS >= config.min_edge)
                && panel
                    .constraint_points
                    .iter()
                    .chain(panel.polygons.iter().flatten())
                    .all(|p| {
                        let distance = point.distance(DVec3::from_array(*p));
                        distance <= EPS || distance + EPS >= config.min_edge
                    })
                && bar.constraints.iter().all(|c| {
                    let distance = (c.t - t).abs() * d.length();
                    distance <= EPS || distance + EPS >= config.min_edge
                })
                && t * d.length() + EPS >= config.min_edge
                && (1.0 - t) * d.length() + EPS >= config.min_edge;
            if !spacing_ok {
                stats.rejected_surface_crossings += 1;
                stats
                    .rejected_crossing_details
                    .push(RejectedSurfaceCrossing {
                        panel_id: panel.id,
                        source_element_ids: bar.source_element_ids.clone(),
                        point: point.to_array(),
                        t,
                        reason: "would_create_short_constraint".into(),
                    });
                continue;
            }
            let existing = bar
                .constraints
                .iter_mut()
                .find(|c| (c.t - t).abs() * d.length() <= EPS);
            if let Some(c) = existing {
                if !c.panel_ids.contains(&panel.id) {
                    c.panel_ids.push(panel.id);
                }
                if !c.kinds.iter().any(|s| s == "surface_intersection") {
                    c.kinds.push("surface_intersection".into());
                }
            } else {
                bar.constraints.push(crate::models::BarConstraint {
                    t,
                    source_node_id: None,
                    kinds: vec!["surface_intersection".into()],
                    panel_ids: vec![panel.id],
                });
            }
            if panel
                .constraint_points
                .iter()
                .all(|p| point.distance(DVec3::from_array(*p)) > EPS)
            {
                panel.constraint_points.push(point.to_array());
            }
            stats.surface_crossings += 1;
        }
        bar.constraints.sort_by(|a, b| a.t.total_cmp(&b.t));
        for c in &mut bar.constraints {
            c.panel_ids.sort_unstable();
            c.panel_ids.dedup();
            c.kinds.sort();
            c.kinds.dedup();
        }
    }
    let mut joint_positions: BTreeMap<u32, Vec<DVec3>> = BTreeMap::new();
    for bar in bars.iter() {
        joint_positions
            .entry(bar.start_node_id)
            .or_default()
            .push(DVec3::from_array(bar.start_point));
        joint_positions
            .entry(bar.end_node_id)
            .or_default()
            .push(DVec3::from_array(bar.end_point));
        for (ci, c) in bar.constraints.iter().enumerate() {
            if let Some(node) = c.source_node_id {
                joint_positions
                    .entry(node)
                    .or_default()
                    .push(anchor_point(bar, Anchor::Interior(ci)));
            }
        }
    }
    stats.joint_residual_tolerance = AXIS_RESIDUAL_TOL;
    for (id, positions) in joint_positions {
        let gap = positions
            .iter()
            .flat_map(|p| positions.iter().map(move |q| p.distance(*q)))
            .fold(0.0, f64::max);
        stats.max_joint_gap = stats.max_joint_gap.max(gap);
        if gap > AXIS_RESIDUAL_TOL {
            stats.unresolved_bar_junctions.push(id);
        } else if gap > EPS {
            stats.numerical_residual_nodes.push(id);
        }
    }
    for panel in panels {
        panel.constraint_points.sort_by(|a, b| {
            a[0].total_cmp(&b[0])
                .then(a[1].total_cmp(&b[1]))
                .then(a[2].total_cmp(&b[2]))
        });
    }
    stats
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::{BarType, PanelType};
    fn panel() -> MacroPanel {
        MacroPanel {
            id: 1,
            panel_type: PanelType::Slab,
            stiffness_id: 1,
            plane_normal: [0., 0., 1.],
            plane_d: 0.,
            polygons: vec![vec![[0., 0., 0.], [2., 0., 0.], [2., 2., 0.], [0., 2., 0.]]],
            filled_holes: 0,
            filled_hole_area: 0.,
            fe_count: 1,
            source_element_ids: vec![100],
            connected_panel_ids: vec![],
            constraint_points: vec![],
        }
    }
    fn bar(a: [f64; 3], b: [f64; 3], start: u32, end: u32) -> MacroBar {
        MacroBar {
            bar_type: BarType::Column,
            stiffness_id: Some(1),
            start_point: a,
            end_point: b,
            length: DVec3::from_array(a).distance(DVec3::from_array(b)),
            start_node_id: start,
            end_node_id: end,
            source_element_ids: vec![1],
            source_node_ids: vec![start, end],
            start_panel_ids: vec![],
            end_panel_ids: vec![],
            constraints: vec![],
            property_spans: vec![],
        }
    }
    #[test]
    fn disconnected_axis_is_inferred_from_surface_geometry() {
        let mut bars = vec![bar([0.5, 1., 0.04], [1.5, 1., 0.04], 1, 2)];
        let mut panels = vec![panel()];
        let stats = attach_bar_endpoints(
            &mut bars,
            &mut panels,
            &MeshData::default(),
            &HashMap::new(),
            &ReconstructionConfig::default(),
        );
        assert_eq!(stats.inferred_panel_ids_by_node.len(), 2);
        assert_eq!(bars[0].start_panel_ids, vec![1]);
        assert_eq!(bars[0].end_panel_ids, vec![1]);
        assert!(contains(&panels[0], DVec3::from_array(bars[0].start_point)));
        assert!(contains(&panels[0], DVec3::from_array(bars[0].end_point)));
        assert_eq!(bars[0].source_element_ids, vec![1]);
        assert!(bars[0].constraints.is_empty());
    }

    #[test]
    fn separated_parallel_alternatives_are_not_arbitrarily_joined() {
        let mut second = panel();
        second.id = 2;
        second.plane_d = -0.02;
        for p in &mut second.polygons[0] {
            p[2] = 0.02;
        }
        let mut panels = vec![panel(), second];
        let mut bars = vec![bar([1., 1., 0.01], [1., 1., 1.], 1, 2)];
        let stats = attach_bar_endpoints(
            &mut bars,
            &mut panels,
            &MeshData::default(),
            &HashMap::new(),
            &ReconstructionConfig::default(),
        );
        assert_eq!(stats.ambiguous_panel_ids_by_node[&1], vec![1, 2]);
        assert!(stats.inferred_panel_ids_by_node.is_empty());
        assert!(bars[0].start_panel_ids.is_empty());
        assert_eq!(bars[0].start_point, [1., 1., 0.01]);
    }

    #[test]
    fn inference_respects_budget_holes_and_existing_source_incidence() {
        let mut surface = panel();
        surface.polygons.push(vec![
            [0.5, 0.5, 0.],
            [0.5, 1.5, 0.],
            [1.5, 1.5, 0.],
            [1.5, 0.5, 0.],
        ]);
        let bars = vec![
            bar([1., 1., 0.01], [1., 1., 0.06], 1, 2),
            bar([0.2, 0.2, 0.01], [0.2, 0.2, 0.06], 3, 4),
        ];
        let mut incident = BTreeMap::from([(3, BTreeSet::from([0]))]);
        let mut stats = BarContactSummary::default();
        infer_incidence(
            &bars,
            &[surface],
            &mut incident,
            &ReconstructionConfig::default(),
            &mut stats,
        );
        assert!(stats.inferred_panel_ids_by_node.is_empty());
        assert_eq!(incident, BTreeMap::from([(3, BTreeSet::from([0]))]));
    }

    #[test]
    fn missing_contact_measurements_use_final_geometry_and_keep_provenance() {
        let bars = vec![bar([1., 1., 0.06], [1., 1., 1.], 1, 2)];
        let mut stats = BarContactSummary::default();
        stats.unresolved_panel_ids_by_node.insert(1, vec![1]);
        let details = measure_missing_contacts(&bars, &[panel()], &stats);
        assert_eq!(details.len(), 1);
        assert_eq!(details[0].source_bar_element_ids, vec![1]);
        assert!(!details[0].inferred);
        assert_eq!(details[0].anchors[0].plane_distance, 0.06);
        assert!(details[0].anchors[0].projection_inside);
        assert_eq!(details[0].anchors[0].boundary_distance, 1.);
    }

    #[test]
    fn coupled_solver_moves_plane_and_whole_axis_together() {
        let mut bars = vec![bar([0.5, 1., 0.04], [1.5, 1., 0.04], 1, 2)];
        let original = vec![[
            DVec3::from_array(bars[0].start_point),
            DVec3::from_array(bars[0].end_point),
        ]];
        let mut panels = vec![panel()];
        let incidence = BTreeMap::from([(1, BTreeSet::from([0])), (2, BTreeSet::from([0]))]);
        assert!(reconcile_axes_and_planes(
            &mut bars,
            &original,
            &ReconstructionConfig::default(),
            &mut panels,
            &incidence,
            &mut CoupledPlaneSearch::default()
        ));
        assert!(panels[0].plane_d.abs() > 1e-4);
        assert!(contains(&panels[0], DVec3::from_array(bars[0].start_point)));
        assert!(contains(&panels[0], DVec3::from_array(bars[0].end_point)));
        assert!((bars[0].length - 1.0).abs() < EPS);
        assert!(panels[0].polygons[0].iter().all(|p| p[2].abs() <= 0.05));
    }
    #[test]
    fn incompatible_planes_are_frozen_and_omitted_incidence_is_reported() {
        let mut second = panel();
        second.id = 2;
        second.plane_d = -0.02;
        for p in second.polygons.iter_mut().flatten() {
            p[2] = 0.02;
        }
        let mut panels = vec![panel(), second];
        let mut bars = vec![bar([1., 1., 0.01], [1., 1., 1.], 1, 2)];
        let original = vec![[
            DVec3::from_array(bars[0].start_point),
            DVec3::from_array(bars[0].end_point),
        ]];
        let incident = BTreeMap::from([(1, BTreeSet::from([0, 1]))]);
        let mut search = CoupledPlaneSearch::default();
        assert!(reconcile_axes_and_planes(
            &mut bars,
            &original,
            &ReconstructionConfig::default(),
            &mut panels,
            &incident,
            &mut search
        ));
        assert_eq!(search.frozen_panel_ids, vec![1, 2]);
        assert_eq!(search.requested_incidence_count, 2);
        assert_eq!(search.retained_incidence_count, 1);
        assert_eq!(search.relaxed_incidence_pairs.len(), 1);
        assert_eq!(panels[0].plane_d, 0.0);
        assert_eq!(panels[1].plane_d, -0.02);
        assert!(search.failure_reason.is_none());
    }
    #[test]
    fn plane_rotation_is_bounded_by_vertex_motion_not_only_angle() {
        let panels = vec![panel()];
        let n = DVec3::new(-0.01, 0., 1.).normalize().to_array();
        let config = ReconstructionConfig::default();
        let tilted =
            crate::geometry::topology::transformed_planes_checked(&panels, &[n], &[0.], &config)
                .unwrap();
        for (a, b) in panels[0].polygons[0].iter().zip(&tilted[0].polygons[0]) {
            assert!(DVec3::from_array(*a).distance(DVec3::from_array(*b)) <= 0.05);
            assert!(DVec3::from_array(n).dot(DVec3::from_array(*b)).abs() < EPS);
        }
        let mut long = panels.clone();
        for p in &mut long[0].polygons[0] {
            p[0] *= 100.;
        }
        assert!(
            crate::geometry::topology::transformed_planes_checked(&long, &[n], &[0.], &config)
                .is_err()
        );
        assert_eq!(panels[0].plane_normal, [0., 0., 1.]);
    }
    #[test]
    fn shifted_planes_keep_shared_vertices_and_reject_excessive_motion() {
        let floor = panel();
        let mut wall = panel();
        wall.id = 2;
        wall.plane_normal = [1., 0., 0.];
        wall.plane_d = 0.;
        wall.polygons = vec![vec![[0., 0., 0.], [0., 2., 0.], [0., 2., 2.], [0., 0., 2.]]];
        let panels = vec![floor, wall];
        let config = ReconstructionConfig::default();
        let moved =
            crate::geometry::topology::shifted_planes(&panels, &[-0.01, -0.02], &config).unwrap();
        assert_eq!(moved[0].polygons[0][0], moved[1].polygons[0][0]);
        assert_eq!(moved[0].polygons[0][3], moved[1].polygons[0][1]);
        assert_eq!(moved[0].polygons[0][0], [0.02, 0., 0.01]);
        assert!(crate::geometry::topology::shifted_planes(&panels, &[-0.1, 0.], &config).is_none());
        assert_eq!(panels[0].polygons[0][0], [0., 0., 0.]);
    }
    #[test]
    fn joint_solver_moves_whole_crossing_axes_within_budget() {
        let mut a = bar([0., 0., 0.], [2., 0., 0.], 1, 3);
        let mut b = bar([1., -1., 0.], [1., 1., 0.], 4, 5);
        for axis in [&mut a, &mut b] {
            axis.constraints.push(crate::models::BarConstraint {
                t: 0.5,
                source_node_id: Some(2),
                kinds: vec!["junction".into()],
                panel_ids: vec![],
            });
        }
        let original = vec![
            [
                DVec3::from_array(a.start_point),
                DVec3::from_array(a.end_point),
            ],
            [
                DVec3::from_array(b.start_point),
                DVec3::from_array(b.end_point),
            ],
        ];
        a.start_point[2] = 0.006;
        a.end_point[2] = 0.006;
        b.start_point[2] = -0.006;
        b.end_point[2] = -0.006;
        let mut bars = vec![a, b];
        assert!(reconcile_axes(
            &mut bars,
            &original,
            &ReconstructionConfig::default()
        ));
        assert!(
            anchor_point(&bars[0], Anchor::Interior(0))
                .distance(anchor_point(&bars[1], Anchor::Interior(0)))
                < EPS
        );
        for (i, axis) in bars.iter().enumerate() {
            assert!(DVec3::from_array(axis.start_point).distance(original[i][0]) <= 0.01 + EPS);
            assert_eq!(axis.constraints[0].t, 0.5);
            assert!(axis.length > 1.99);
        }
    }
    #[test]
    fn snap_uses_fifty_mm_budget_and_slides_property_boundary() {
        let mut main = bar([0., 0., 0.], [2., 0., 0.], 1, 3);
        main.constraints.push(crate::models::BarConstraint {
            t: 0.5,
            source_node_id: Some(2),
            kinds: vec!["junction".into(), "property_boundary".into()],
            panel_ids: vec![],
        });
        main.property_spans = vec![
            crate::models::BarPropertySpan {
                start_t: 0.,
                end_t: 0.5,
                stiffness_id: 1,
                source_element_ids: vec![1],
            },
            crate::models::BarPropertySpan {
                start_t: 0.5,
                end_t: 1.,
                stiffness_id: 2,
                source_element_ids: vec![2],
            },
        ];
        let branch = bar([1.04, 0.04, 0.], [1.04, 1., 0.], 2, 4);
        let mut bars = vec![main, branch];
        let original: Vec<_> = bars
            .iter()
            .map(|b| {
                [
                    DVec3::from_array(b.start_point),
                    DVec3::from_array(b.end_point),
                ]
            })
            .collect();
        assert!(reconcile_axes(
            &mut bars,
            &original,
            &ReconstructionConfig::default()
        ));
        assert!((bars[0].constraints[0].t - 0.5).abs() > 1e-4);
        assert!(
            anchor_point(&bars[0], Anchor::Interior(0))
                .distance(DVec3::from_array(bars[1].start_point))
                < EPS
        );
        assert_eq!(bars[0].property_spans[0].end_t, bars[0].constraints[0].t);
        assert_eq!(bars[0].property_spans[1].start_t, bars[0].constraints[0].t);
        for (i, b) in bars.iter().enumerate() {
            assert!(DVec3::from_array(b.start_point).distance(original[i][0]) <= 0.05 + EPS);
        }
    }
    #[test]
    fn infeasible_joint_solver_does_not_commit_partial_changes() {
        let a = bar([0., 0., 0.], [2., 0., 0.], 1, 1);
        let original = vec![[DVec3::ZERO, DVec3::X * 2.]];
        let mut bars = vec![a];
        assert!(!reconcile_axes(
            &mut bars,
            &original,
            &ReconstructionConfig::default()
        ));
        assert_eq!(bars[0].start_point, [0., 0., 0.]);
        assert_eq!(bars[0].end_point, [2., 0., 0.]);
    }
    #[test]
    fn shared_endpoints_move_together_and_create_surface_constraint() {
        let mut bars = vec![
            bar([1., 1., 0.005], [1., 1., 1.], 1, 2),
            bar([1., 1., 0.005], [0., 1., 1.], 1, 3),
        ];
        let mut panels = vec![panel()];
        let stats = attach_bar_endpoints(
            &mut bars,
            &mut panels,
            &MeshData::default(),
            &HashMap::new(),
            &ReconstructionConfig::default(),
        );
        assert!(contains(&panels[0], DVec3::from_array(bars[0].start_point)));
        assert_eq!(bars[0].start_point, bars[1].start_point);
        assert_eq!(bars[0].start_panel_ids, vec![1]);
        assert_eq!(panels[0].constraint_points, vec![bars[0].start_point]);
        assert!(DVec3::from_array(bars[0].start_point).distance(DVec3::new(1., 1., 0.005)) <= 0.01);
        assert!(stats.max_displacement <= 0.01);
    }
    #[test]
    fn opening_is_not_treated_as_solid_surface() {
        let mut panels = vec![panel()];
        panels[0].polygons.push(vec![
            [0.5, 0.5, 0.],
            [0.5, 1.5, 0.],
            [1.5, 1.5, 0.],
            [1.5, 0.5, 0.],
        ]);
        let mut bars = vec![bar([1., 1., 0.005], [1., 1., 1.], 1, 2)];
        let stats = attach_bar_endpoints(
            &mut bars,
            &mut panels,
            &MeshData::default(),
            &HashMap::new(),
            &ReconstructionConfig::default(),
        );
        assert_eq!(stats.attached_groups, 0);
        assert!(panels[0].constraint_points.is_empty());
        assert_eq!(bars[0].start_point, [1., 1., 0.005]);
    }
    #[test]
    fn conflicting_surfaces_preserve_feasible_incidence_and_report_missing_one() {
        let mut second = panel();
        second.id = 2;
        second.plane_d = -0.006;
        for point in &mut second.polygons[0] {
            point[2] = 0.006;
        }
        second.source_element_ids = vec![101];
        let mut panels = vec![panel(), second];
        let mesh = MeshData {
            nodes: HashMap::new(),
            elements: [100, 101]
                .into_iter()
                .map(|id| crate::models::ElementData {
                    id,
                    elem_type: 42,
                    stiff_id: 1,
                    nodes: vec![1, 10, 11],
                })
                .collect(),
        };
        let mut bars = vec![bar([1., 1., 0.003], [1., 1., 1.], 1, 2)];
        let stats = attach_bar_endpoints(
            &mut bars,
            &mut panels,
            &mesh,
            &HashMap::new(),
            &ReconstructionConfig::default(),
        );
        assert_eq!(stats.rejected_groups, 1);
        assert_eq!(stats.partial_groups, 1);
        assert_eq!(stats.attached_groups, 1);
        assert_eq!(bars[0].start_panel_ids.len(), 1);
        assert_eq!(stats.unresolved_panel_ids_by_node[&1].len(), 1);
        assert_ne!(
            bars[0].start_panel_ids[0],
            stats.unresolved_panel_ids_by_node[&1][0]
        );
        assert!(contains(
            &panels[(bars[0].start_panel_ids[0] - 1) as usize],
            DVec3::from_array(bars[0].start_point)
        ));
    }
    #[test]
    fn near_boundary_point_is_snapped_instead_of_creating_a_sliver() {
        let mut panels = vec![panel()];
        let mut bars = vec![bar([1., 0.005, 0.], [1., 0.005, 1.], 1, 2)];
        let stats = attach_bar_endpoints(
            &mut bars,
            &mut panels,
            &MeshData::default(),
            &HashMap::new(),
            &ReconstructionConfig::default(),
        );
        assert_eq!(stats.attached_groups, 1);
        assert_eq!(bars[0].start_point, [1., 0., 0.]);
        assert_eq!(panels[0].constraint_points, vec![[1., 0., 0.]]);
    }
    #[test]
    fn rejected_interior_crossing_has_location_and_source_provenance() {
        let mut panels = vec![panel()];
        let mut bars = vec![bar([1., 0.005, -1.], [1., 0.005, 1.], 1, 2)];
        let stats = attach_bar_endpoints(
            &mut bars,
            &mut panels,
            &MeshData::default(),
            &HashMap::new(),
            &ReconstructionConfig::default(),
        );
        assert_eq!(stats.rejected_surface_crossings, 1);
        let detail = &stats.rejected_crossing_details[0];
        assert_eq!(detail.panel_id, 1);
        assert_eq!(detail.source_element_ids, vec![1]);
        assert_eq!(detail.point, [1., 0.005, 0.]);
        assert_eq!(detail.t, 0.5);
        assert!(panels[0].constraint_points.is_empty());
    }
    #[test]
    fn through_column_has_two_ends_and_an_interior_surface_constraint() {
        let mut panels = vec![panel()];
        let mut bars = vec![bar([1., 1., -1.], [1., 1., 1.], 1, 2)];
        let stats = attach_bar_endpoints(
            &mut bars,
            &mut panels,
            &MeshData::default(),
            &HashMap::new(),
            &ReconstructionConfig::default(),
        );
        assert_eq!(bars.len(), 1);
        assert_eq!(bars[0].start_point, [1., 1., -1.]);
        assert_eq!(bars[0].end_point, [1., 1., 1.]);
        assert_eq!(bars[0].constraints.len(), 1);
        assert_eq!(bars[0].constraints[0].t, 0.5);
        assert_eq!(bars[0].constraints[0].panel_ids, vec![1]);
        assert_eq!(stats.surface_crossings, 1);
    }
    #[test]
    fn branch_endpoint_is_kept_on_an_unsplit_main_axis() {
        let mut main = bar([0., 0., 0.], [2., 0., 0.], 1, 3);
        main.constraints.push(crate::models::BarConstraint {
            t: 0.5,
            source_node_id: Some(2),
            kinds: vec!["junction".into()],
            panel_ids: vec![],
        });
        let branch = bar([1., 0.005, 0.], [1., 1., 0.], 2, 4);
        let mut bars = vec![main, branch];
        attach_bar_endpoints(
            &mut bars,
            &mut [],
            &MeshData::default(),
            &HashMap::new(),
            &ReconstructionConfig::default(),
        );
        assert_eq!(bars.len(), 2);
        assert!(on_axis(&bars[0], DVec3::from_array(bars[1].start_point)));
        assert!((bars[0].length - 2.0).abs() < EPS);
        assert!(DVec3::from_array(bars[0].start_point).distance(DVec3::ZERO) <= 0.01);
        assert!(DVec3::from_array(bars[0].end_point).distance(DVec3::X * 2.0) <= 0.01);
    }
}
