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
    /// Maximum width of a free console beyond a junction line that is
    /// trimmed (model units). The default, 0.25 m, is half the default mesh
    /// spacing: such a strip cannot hold elements of the target size.
    pub maximum_console_width: f64,
    /// Maximum offset between parallel walls on opposite sides of one slab
    /// for the upper wall to adopt the plane of the wall carrying it (model
    /// units). The default, 0.05 m, equals the default junction movement
    /// limit that bounds every closing vertex.
    pub maximum_stack_offset: f64,
    /// Largest move closing a wall end onto the axis of another wall, and
    /// the edge length below which redundant collinear vertices are removed
    /// (model units, default 0.05 m).
    pub maximum_wall_end_snap: f64,
    /// Widest void between unconnected boundary nodes of one planar region
    /// that is treated as a crack of the source mesh and left out of the
    /// region contour (model units, default 0.01 m; 0 disables).
    pub maximum_crack_width: f64,
    /// Surface edges shorter than this whose ends are both needed corners
    /// collapse into one vertex (model units, default 0.05 m: a tenth of a
    /// 0.5 m PLAXIS element; 0 disables).
    pub maximum_collapsed_edge: f64,
    /// Widest gap between a surface vertex (or bar node) and another surface
    /// that is closed onto it (model units, default 0.05 m: a tenth of a
    /// 0.5 m PLAXIS element; 0 keeps every gap, e.g. for real joints).
    pub maximum_gap: f64,
    /// Also close gaps across a plane (a wall top below a slab), not only
    /// gaps within one plane (default true). Structures parallel to each
    /// other (slabs at different levels) are never brought together.
    pub close_offset_gaps: bool,
}
impl Default for FeaturePolicy {
    fn default() -> Self {
        Self {
            maximum_source_width: 0.05,
            maximum_filled_area_ratio: 0.001,
            maximum_console_width: 0.25,
            maximum_stack_offset: 0.05,
            maximum_wall_end_snap: 0.05,
            maximum_crack_width: 0.01,
            maximum_collapsed_edge: 0.05,
            maximum_gap: 0.05,
            close_offset_gaps: true,
        }
    }
}

#[derive(Debug, Serialize)]
pub struct SimplifiedHole {
    /// `collapsed_opening` (closed by the frame to zero width),
    /// `crack_void` (an enclosed void of the source mesh no wider than the
    /// crack width: a sliver between non-conforming elements) or `seam` (a
    /// crack at a hanging node closed off an opening or contour;
    /// `source_nodes` are its tips, kept as interior vertices).
    pub reason: String,
    pub patch: usize,
    pub source_elements: Vec<u32>,
    pub source_nodes: Vec<u32>,
    /// Source nodes kept as interior vertices: all of a collapsed opening;
    /// of a crack void or seam (source mesh discretization) only those
    /// shared with another surface or a bar. Set by the assembly.
    pub retained_nodes: Vec<u32>,
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

/// Close the seams of a non-conforming source mesh in an invalid ring: a
/// ring node lying on another edge of the ring (a hanging node, within the
/// crack width of it in source coordinates) splits that edge, and the seam
/// left behind (a node left and returned to: a, b, a) closes. A ring
/// changes only when it is invalid and becomes valid with at least three
/// nodes (a ring that is only a crack is left to `simplify`); the source
/// area closed counts against the filled area budget. The seam tips are
/// reported and kept as interior vertices.
#[allow(clippy::too_many_arguments)]
fn close_seams(
    loops: &mut [Vec<u32>],
    mesh: &MeshData,
    points: &BTreeMap<u32, DVec3>,
    plane: &PlaneFrame,
    precision: f64,
    policy: &FeaturePolicy,
    patch: usize,
    ids: &[u32],
) -> Vec<SimplifiedHole> {
    let width = policy.maximum_crack_width;
    if width <= 0. || loops.is_empty() {
        return vec![];
    }
    let uv = |ring: &[u32]| -> Vec<[f64; 2]> {
        ring.iter()
            .map(|n| plane.project(points[n].to_array()))
            .collect()
    };
    let source = |n: u32| mesh.nodes[&n];
    let source_area = |ring: &[u32]| ring_area(ring, plane, |n| mesh.nodes[&n].to_array());
    let exterior_area = source_area(&loops[0]);
    let mut filled = 0.;
    let mut changes = vec![];
    for ring in loops.iter_mut() {
        if super::super::validate_ring(&uv(ring), precision).is_ok() {
            continue;
        }
        let mut r = ring.clone();
        let mut tips = vec![];
        let mut seam_width = 0.0_f64;
        // Each split puts a node of the ring on an edge: bounded.
        for _ in 0..ring.len() {
            let n = r.len();
            let split = (0..n).find_map(|i| {
                let (a, b) = (r[i], r[(i + 1) % n]);
                let (pa, d) = (source(a), source(b) - source(a));
                let length = d.length();
                r.iter().find_map(|&h| {
                    let t = (source(h) - pa).dot(d) / (length * length);
                    let distance = (pa + d * t).distance(source(h));
                    (h != a
                        && h != b
                        && t * length > precision
                        && (1. - t) * length > precision
                        && distance <= width)
                        .then_some((i, h, distance))
                })
            });
            let Some((i, h, distance)) = split else {
                break;
            };
            seam_width = seam_width.max(distance);
            r.insert(i + 1, h);
            while r.len() >= 3 {
                let n = r.len();
                let Some(j) = (0..n).find(|&j| r[(j + n - 1) % n] == r[(j + 1) % n]) else {
                    break;
                };
                tips.push(r[j]);
                let next = (j + 1) % n;
                r.remove(j.max(next));
                r.remove(j.min(next));
            }
            if r.len() < 3 {
                break;
            }
        }
        let closed = r.len() >= 3
            && !tips.is_empty()
            && r.iter().collect::<std::collections::BTreeSet<_>>().len() == r.len()
            && super::super::validate_ring(&uv(&r), precision).is_ok();
        let area = (source_area(ring) - source_area(&r)).abs();
        if !closed || filled + area > policy.maximum_filled_area_ratio * exterior_area {
            continue;
        }
        filled += area;
        tips.sort_unstable();
        tips.dedup();
        changes.push(SimplifiedHole {
            reason: "seam".into(),
            patch,
            source_elements: ids.to_vec(),
            retained_nodes: tips.clone(),
            source_nodes: tips,
            source_area: area,
            source_width: seam_width,
            candidate_area: 0.,
            candidate_width: 0.,
            numerical_area_threshold: 0.,
            exterior_source_area: exterior_area,
        });
        *ring = r;
    }
    changes
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
    let mut changes = close_seams(loops, mesh, points, plane, precision, policy, patch, ids);
    if loops.len() < 2 {
        return changes;
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
    let mut filled: f64 = changes.iter().map(|c| c.source_area).sum();
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
        let collapsed = source_width <= policy.maximum_source_width
            && candidate_width <= precision
            && candidate_area <= threshold;
        let crack = source_width <= policy.maximum_crack_width;
        if source_valid
            && (collapsed || crack)
            && exterior_area > 0.
            && filled + source_area <= policy.maximum_filled_area_ratio * exterior_area
        {
            filled += source_area;
            changes.push(SimplifiedHole {
                reason: if collapsed {
                    "collapsed_opening"
                } else {
                    "crack_void"
                }
                .into(),
                patch,
                source_elements: ids.to_vec(),
                source_nodes: ring.clone(),
                retained_nodes: ring.clone(),
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
    fn seam_at_a_hanging_node_closes() {
        for transform in [false, true] {
            let apply = |p: [f64; 2]| {
                let p = DVec3::new(p[0], p[1], 0.);
                if transform {
                    glam::DQuat::from_rotation_x(0.7) * p + DVec3::new(31., -14., 8.)
                } else {
                    p
                }
            };
            let mut mesh = MeshData::default();
            // Node 7 hangs 2 mm beside the element edge 5-6 (1.1 m), aligned
            // onto it by the frame: the hole 4-5-6-7 runs up the seam and
            // back.
            for (i, p) in [
                [0., 0.],
                [10., 0.],
                [10., 10.],
                [0., 10.],
                [6., 4.2],
                [5., 4.],
                [5., 5.1],
                [5.002, 4.2],
            ]
            .into_iter()
            .enumerate()
            {
                mesh.nodes.insert(i as u32, apply(p));
            }
            let mut points: BTreeMap<_, _> = mesh.nodes.iter().map(|(&n, &p)| (n, p)).collect();
            points.insert(7, apply([5., 4.2]));
            let normal = if transform {
                glam::DQuat::from_rotation_x(0.7) * DVec3::Z
            } else {
                DVec3::Z
            };
            let plane = PlaneFrame::new(apply([0., 0.]).to_array(), normal.to_array()).unwrap();
            let run = |hole: Vec<u32>| {
                let mut loops = vec![vec![0, 1, 2, 3], hole];
                let r = simplify(
                    &mut loops,
                    &mesh,
                    &points,
                    &plane,
                    1e-7,
                    &FeaturePolicy::default(),
                    0,
                    &[12],
                );
                (loops, r)
            };
            // The seam closes; the triangle 4-5-7 stays an opening.
            let (loops, r) = run(vec![4, 5, 6, 7]);
            assert_eq!(r.len(), 1);
            assert_eq!(r[0].reason, "seam");
            assert_eq!(r[0].source_nodes, vec![6]);
            assert_eq!(loops[1], vec![4, 5, 7]);
            assert!((r[0].source_width - 0.002).abs() < 1e-9);
            // A hole that is only a crack is left to the crack rule.
            let (_, r) = run(vec![5, 6, 7]);
            assert!(r.iter().all(|c| c.reason != "seam"));
        }
    }

    #[test]
    fn closure_is_bounded_and_rigid_transform_invariant() {
        for transform in [false, true] {
            // Collapsed by the frame: filled up to the source width limit.
            // Not collapsed: filled only as a crack void (within the crack
            // width, 10 mm by default); a 20 mm opening stays.
            for (width, collapsed, expected) in [
                (0.01, true, 1),
                (0.005, false, 1),
                (0.02, false, 0),
                (0.2, true, 0),
            ] {
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
                assert_eq!(
                    result.len(),
                    expected,
                    "width {width} collapsed {collapsed}"
                );
                assert_eq!(loops.len(), 2 - expected);
                if expected == 1 {
                    let reason = if collapsed {
                        "collapsed_opening"
                    } else {
                        "crack_void"
                    };
                    assert_eq!(result[0].reason, reason);
                }
            }
        }
    }

    #[test]
    fn cumulative_area_budget_and_crossed_source_are_not_bypassed() {
        let (mut mesh, mut points, plane) = fixture(0.01, true, false);
        let mut loops = vec![vec![0, 1, 2, 3], vec![4, 5, 6]];
        let policy = FeaturePolicy {
            maximum_filled_area_ratio: 0.000001,
            maximum_console_width: 0.,
            maximum_stack_offset: 0.,
            maximum_wall_end_snap: 0.,
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
