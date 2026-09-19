//! Bounded removal of numerically collapsed openings, with source evidence.
use super::{ring_area, MeshData, PlaneFrame};
use geo::{Contains, Intersects, LineString, Polygon};
use glam::DVec3;
use serde::Serialize;
use std::collections::BTreeMap;

#[derive(Debug, Clone, Serialize)]
pub struct FeaturePolicy {
    /// Maximum source distance from the longest hole chord (model units).
    pub maximum_source_width: f64,
    /// Cumulative filled source area / exterior area, per property region.
    pub maximum_filled_area_ratio: f64,
}
impl Default for FeaturePolicy {
    fn default() -> Self {
        Self {
            maximum_source_width: 0.05,
            maximum_filled_area_ratio: 0.001,
        }
    }
}

#[derive(Debug, Serialize)]
pub struct SimplifiedHole {
    pub patch: usize,
    pub source_elements: Vec<u32>,
    pub source_nodes: Vec<u32>,
    pub source_area: f64,
    pub source_width: f64,
    pub candidate_area: f64,
    pub candidate_width: f64,
    pub numerical_area_threshold: f64,
    pub exterior_source_area: f64,
}

fn dimensions(points: &[DVec3]) -> (f64, f64) {
    let mut ends = (points[0], points[0]);
    let mut length = 0.0_f64;
    for &a in points {
        for &b in points {
            if a.distance(b) > length {
                length = a.distance(b);
                ends = (a, b);
            }
        }
    }
    if length == 0. {
        return (0., 0.);
    }
    let direction = (ends.1 - ends.0) / length;
    let width = points
        .iter()
        .map(|&p| {
            let t = (p - ends.0).dot(direction).clamp(0., length);
            p.distance(ends.0 + t * direction)
        })
        .fold(0.0_f64, f64::max);
    (length, width)
}

pub(super) fn simplify(
    loops: &mut Vec<Vec<u32>>,
    mesh: &MeshData,
    points: &BTreeMap<u32, DVec3>,
    plane: &PlaneFrame,
    precision: f64,
    policy: &FeaturePolicy,
    patch: usize,
    ids: &[u32],
) -> Vec<SimplifiedHole> {
    if loops.len() < 2 {
        return vec![];
    }
    let exterior_area = ring_area(&loops[0], plane, |n| mesh.nodes[&n].to_array());
    let polygon = |uv: &Vec<[f64; 2]>| {
        Polygon::new(
            LineString::from(
                uv.iter()
                    .chain(uv.first())
                    .map(|p| (p[0], p[1]))
                    .collect::<Vec<_>>(),
            ),
            vec![],
        )
    };
    let outer_uv = loops[0]
        .iter()
        .map(|n| plane.project(mesh.nodes[n].to_array()))
        .collect();
    let outer = polygon(&outer_uv);
    let mut filled = 0.;
    let mut changes = vec![];
    let mut keep = vec![loops[0].clone()];
    for ring in loops.iter().skip(1) {
        let original: Vec<_> = ring.iter().map(|n| mesh.nodes[n]).collect();
        let candidate: Vec<_> = ring.iter().map(|n| points[n]).collect();
        let (_, source_width) = dimensions(&original);
        let (extent, candidate_width) = dimensions(&candidate);
        let source_area = ring_area(ring, plane, |n| mesh.nodes[&n].to_array());
        let candidate_area = ring_area(ring, plane, |n| points[&n].to_array());
        let threshold = precision * extent.max(precision);
        // Width prevents cancellation in a crossed contour from looking like
        // a collapsed opening. Source validation prevents filling bad rings.
        let uv: Vec<_> = original
            .iter()
            .map(|p| plane.project(p.to_array()))
            .collect();
        let hole = polygon(&uv);
        let source_valid = super::super::validate_ring(&outer_uv, precision).is_ok()
            && super::super::validate_ring(&uv, precision).is_ok()
            && outer.contains(&hole)
            && !outer.exterior().intersects(hole.exterior());
        if source_valid
            && source_width <= policy.maximum_source_width
            && candidate_width <= precision
            && candidate_area <= threshold
            && exterior_area > 0.
            && filled + source_area <= policy.maximum_filled_area_ratio * exterior_area
        {
            filled += source_area;
            changes.push(SimplifiedHole {
                patch,
                source_elements: ids.to_vec(),
                source_nodes: ring.clone(),
                source_area,
                source_width,
                candidate_area,
                candidate_width,
                numerical_area_threshold: threshold,
                exterior_source_area: exterior_area,
            });
        } else {
            keep.push(ring.clone());
        }
    }
    *loops = keep;
    changes
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture(
        width: f64,
        collapsed: bool,
        transform: bool,
    ) -> (MeshData, BTreeMap<u32, DVec3>, PlaneFrame) {
        let apply = |p: DVec3| {
            if transform {
                glam::DQuat::from_rotation_x(0.7) * p + DVec3::new(31., -14., 8.)
            } else {
                p
            }
        };
        let mut mesh = MeshData::default();
        for (i, p) in [
            [0., 0., 0.],
            [10., 0., 0.],
            [10., 10., 0.],
            [0., 10., 0.],
            [4., 4., 0.],
            [4.4, 4., 0.],
            [4.2, 4. + width, 0.],
        ]
        .into_iter()
        .enumerate()
        {
            mesh.nodes.insert(i as u32, apply(DVec3::from_array(p)));
        }
        let mut candidate: BTreeMap<_, _> = mesh.nodes.iter().map(|(&n, &p)| (n, p)).collect();
        if collapsed {
            candidate.insert(6, apply(DVec3::new(4.2, 4., 0.)));
        }
        let normal = if transform {
            glam::DQuat::from_rotation_x(0.7) * DVec3::Z
        } else {
            DVec3::Z
        };
        (
            mesh,
            candidate,
            PlaneFrame::new(apply(DVec3::ZERO).to_array(), normal.to_array()).unwrap(),
        )
    }

    #[test]
    fn closure_is_bounded_and_rigid_transform_invariant() {
        for transform in [false, true] {
            for (width, collapsed, expected) in [(0.01, true, 1), (0.01, false, 0), (0.2, true, 0)]
            {
                let (mesh, points, plane) = fixture(width, collapsed, transform);
                let mut loops = vec![vec![0, 1, 2, 3], vec![4, 5, 6]];
                let result = simplify(
                    &mut loops,
                    &mesh,
                    &points,
                    &plane,
                    1e-7,
                    &FeaturePolicy::default(),
                    0,
                    &[12],
                );
                assert_eq!(result.len(), expected);
                assert_eq!(loops.len(), 2 - expected);
            }
        }
    }

    #[test]
    fn cumulative_area_budget_and_crossed_source_are_not_bypassed() {
        let (mut mesh, mut points, plane) = fixture(0.01, true, false);
        let mut loops = vec![vec![0, 1, 2, 3], vec![4, 5, 6]];
        let policy = FeaturePolicy {
            maximum_filled_area_ratio: 0.000001,
            ..FeaturePolicy::default()
        };
        assert!(simplify(&mut loops, &mesh, &points, &plane, 1e-7, &policy, 0, &[12]).is_empty());
        mesh.nodes.insert(7, DVec3::new(4.4, 4.01, 0.));
        mesh.nodes.insert(8, DVec3::new(4., 4.01, 0.));
        points.insert(7, DVec3::new(4.4, 4., 0.));
        points.insert(8, DVec3::new(4., 4., 0.));
        let mut crossed = vec![vec![0, 1, 2, 3], vec![4, 7, 5, 8]];
        assert!(simplify(
            &mut crossed,
            &mesh,
            &points,
            &plane,
            1e-7,
            &FeaturePolicy::default(),
            0,
            &[12]
        )
        .is_empty());
    }
}
