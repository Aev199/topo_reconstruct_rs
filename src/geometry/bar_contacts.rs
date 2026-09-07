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

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct BarContactSummary {
    pub fitted_axes: usize,
    pub joint_solve_converged: bool,
    pub max_axis_fit_displacement: f64,
    pub anchor_groups: usize,
    pub surface_crossings: usize,
    pub rejected_surface_crossings: usize,
    pub attached_groups: usize,
    pub moved_groups: usize,
    pub rejected_groups: usize,
    pub rejected_node_ids: Vec<u32>,
    pub rejection_reasons: BTreeMap<String, Vec<u32>>,
    pub max_displacement: f64,
    pub unresolved_bar_junctions: Vec<u32>,
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
// balls. Unknowns are only the two ends of each whole straight axis.
fn reconcile_axes(
    bars: &mut [MacroBar],
    original: &[[DVec3; 2]],
    config: &ReconstructionConfig,
) -> bool {
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
    let equations: Vec<Vec<(usize, f64)>> = joints
        .values()
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
    let mut points: Vec<_> = bars
        .iter()
        .flat_map(|b| {
            [
                DVec3::from_array(b.start_point),
                DVec3::from_array(b.end_point),
            ]
        })
        .collect();
    for _ in 0..4000 {
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
        if error < EPS * 0.1 {
            if bars.iter().enumerate().any(|(i, _)| {
                let d = points[2 * i + 1] - points[2 * i];
                let old = original[i][1] - original[i][0];
                d.length() + EPS < old.length().min(config.min_edge) || d.dot(old) <= 0.0
            }) {
                return false;
            }
            for (i, bar) in bars.iter_mut().enumerate() {
                bar.start_point = points[2 * i].to_array();
                bar.end_point = points[2 * i + 1].to_array();
                bar.length = points[2 * i].distance(points[2 * i + 1]);
                bar.bar_type = classify_bar(points[2 * i], points[2 * i + 1]);
            }
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
    stats.joint_solve_converged = reconcile_axes(bars, &original_ends, config);
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
                .find(|&point| {
                    safe(point)
                        && indices.iter().all(|&i| {
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
                        })
                })
                .ok_or("would_bend_axis_exceed_budget_or_create_short_constraint")
        };
        let (indices, point) = match choose(&indices) {
            Ok(point) => (indices, point),
            Err(reason) => {
                reject(&mut stats, node, reason);
                // Failure to match a surface must not prevent a safe bar-to-bar join.
                if ends.len() < 2 {
                    continue;
                }
                match choose(&BTreeSet::new()) {
                    Ok(point) => (BTreeSet::new(), point),
                    Err(_) => continue,
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
    stats.unresolved_bar_junctions = joint_positions
        .into_iter()
        .filter(|(_, points)| points.iter().any(|p| p.distance(points[0]) > EPS))
        .map(|(id, _)| id)
        .collect();
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
        assert_eq!(bars[0].start_point, [1., 1., 0.]);
        assert_eq!(bars[0].start_point, bars[1].start_point);
        assert_eq!(bars[0].start_panel_ids, vec![1]);
        assert_eq!(panels[0].constraint_points, vec![[1., 1., 0.]]);
        assert_eq!(stats.moved_groups, 1);
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
    fn conflicting_surfaces_are_reported_without_moving_endpoint() {
        let mut second = panel();
        second.id = 2;
        second.plane_d = -0.006;
        for point in &mut second.polygons[0] {
            point[2] = 0.006;
        }
        let mut panels = vec![panel(), second];
        let mut bars = vec![bar([1., 1., 0.003], [1., 1., 1.], 1, 2)];
        let stats = attach_bar_endpoints(
            &mut bars,
            &mut panels,
            &MeshData::default(),
            &HashMap::new(),
            &ReconstructionConfig::default(),
        );
        assert_eq!(stats.rejected_groups, 1);
        assert_eq!(bars[0].start_point, [1., 1., 0.003]);
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
