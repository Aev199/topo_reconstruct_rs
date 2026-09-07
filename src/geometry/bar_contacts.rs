//! Endpoint-to-surface geometry constraints; these do not assign mechanical ties.
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
    pub endpoint_groups: usize,
    pub attached_groups: usize,
    pub moved_groups: usize,
    pub rejected_groups: usize,
    pub rejected_node_ids: Vec<u32>,
    pub max_displacement: f64,
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
fn project(start: DVec3, panels: &[MacroPanel], indices: &BTreeSet<usize>) -> Option<DVec3> {
    let mut basis: Vec<(DVec3, f64)> = vec![];
    for &index in indices {
        let mut n = DVec3::from_array(panels[index].plane_normal);
        let mut rhs = -panels[index].plane_d - n.dot(start);
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
    let mut groups: BTreeMap<u32, Vec<(usize, bool)>> = BTreeMap::new();
    for (i, b) in bars.iter().enumerate() {
        groups.entry(b.start_node_id).or_default().push((i, true));
        groups.entry(b.end_node_id).or_default().push((i, false));
    }
    stats.endpoint_groups = groups.len();
    for (node, ends) in groups {
        let (index, start) = ends[0];
        let origin = DVec3::from_array(if start {
            bars[index].start_point
        } else {
            bars[index].end_point
        });
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
        if indices.is_empty() {
            continue;
        }
        if indices.iter().any(|&i| !nearby(i)) {
            stats.rejected_groups += 1;
            stats.rejected_node_ids.push(node);
            continue;
        }
        let Some(projected) = project(origin, panels, &indices) else {
            stats.rejected_groups += 1;
            stats.rejected_node_ids.push(node);
            continue;
        };
        let mut candidates = vec![projected];
        candidates.extend(
            indices
                .iter()
                .map(|&i| nearest_boundary(&panels[i], projected)),
        );
        for &i in &indices {
            candidates.extend(
                panels[i]
                    .constraint_points
                    .iter()
                    .chain(panels[i].polygons.iter().flatten())
                    .map(|p| DVec3::from_array(*p)),
            );
        }
        candidates.sort_by(|a, b| a.distance(origin).total_cmp(&b.distance(origin)));
        let point = candidates.into_iter().find(|&q| {
            q.distance(origin) <= config.joint_tol + EPS
                && indices.iter().all(|&i| {
                    let panel = &panels[i];
                    let boundary_distance = nearest_boundary(panel, q).distance(q);
                    contains(panel, q)
                        && (boundary_distance <= EPS || boundary_distance + EPS >= config.min_edge)
                        && panel
                            .constraint_points
                            .iter()
                            .chain(panel.polygons.iter().flatten())
                            .all(|p| {
                                let distance = DVec3::from_array(*p).distance(q);
                                distance <= EPS || distance + EPS >= config.min_edge
                            })
                })
        });
        let Some(point) = point else {
            stats.rejected_groups += 1;
            stats.rejected_node_ids.push(node);
            continue;
        };
        // All bars sharing the original endpoint move together, or none do.
        let safe = ends.iter().all(|&(i, start)| {
            let bar = &bars[i];
            let old_start = DVec3::from_array(bar.start_point);
            let old_end = DVec3::from_array(bar.end_point);
            let (a, b) = if start {
                (point, old_end)
            } else {
                (old_start, point)
            };
            point.distance(if start { old_start } else { old_end }) <= config.joint_tol + EPS
                && a.distance(b) + EPS >= old_start.distance(old_end).min(config.min_edge)
                && (b - a).dot(old_end - old_start) > 0.0
        });
        if !safe {
            stats.rejected_groups += 1;
            stats.rejected_node_ids.push(node);
            continue;
        }
        let panel_ids: Vec<_> = indices.iter().map(|&i| panels[i].id).collect();
        for (i, start) in ends {
            let bar = &mut bars[i];
            if start {
                bar.start_point = point.to_array();
                bar.start_panel_ids = panel_ids.clone();
            } else {
                bar.end_point = point.to_array();
                bar.end_panel_ids = panel_ids.clone();
            }
            let a = DVec3::from_array(bar.start_point);
            let b = DVec3::from_array(bar.end_point);
            bar.length = a.distance(b);
            bar.bar_type = classify_bar(a, b);
        }
        for i in indices {
            if panels[i]
                .constraint_points
                .iter()
                .all(|p| DVec3::from_array(*p).distance(point) > EPS)
            {
                panels[i].constraint_points.push(point.to_array());
            }
        }
        stats.attached_groups += 1;
        let movement = origin.distance(point);
        stats.max_displacement = stats.max_displacement.max(movement);
        if movement > EPS {
            stats.moved_groups += 1;
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
            stiffness_id: 1,
            start_point: a,
            end_point: b,
            length: DVec3::from_array(a).distance(DVec3::from_array(b)),
            start_node_id: start,
            end_node_id: end,
            source_element_ids: vec![1],
            source_node_ids: vec![start, end],
            start_panel_ids: vec![],
            end_panel_ids: vec![],
        }
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
}
