//! Alignment of stacked walls to the axis of the wall carrying them.
//!
//! Walls above and below one slab often have aligned outer faces but
//! different thicknesses, so their mid-planes are a few centimetres apart.
//! Two junction lines millimetres apart on the slab then force tiny
//! elements. A parallel wall standing on the opposite side of the same slab,
//! overlapping the lower wall along the junction and offset by at most the
//! stacking tolerance, adopts the support plane of the lowest wall of its
//! stack. Vertices then close onto the new plane through the ordinary
//! support intersection with its movement limits; a pair whose vertices would
//! exceed them is kept and reported.
use super::{frame, intersection, movement_budget, MeshData, Policy};
use crate::reconstruction::PlaneFrame;
use glam::DVec3;
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Serialize)]
pub struct StackedWall {
    /// Source plane patch of the moved (upper) wall.
    pub upper: usize,
    /// Patch directly below it, across `slab`.
    pub lower: usize,
    /// Lowest patch of the stack, whose plane is adopted.
    pub root: usize,
    pub slab: usize,
    /// Distance of the upper wall to the adopted plane.
    pub offset: f64,
    pub source_elements: Vec<u32>,
    pub reason: String,
}

#[derive(Debug, Default, Serialize)]
pub struct Report {
    pub tolerance: f64,
    pub aligned: Vec<StackedWall>,
    pub kept: Vec<StackedWall>,
    pub identified: Vec<Identified>,
}

struct Link {
    lower: usize,
    slab: usize,
    overlap: f64,
    offset: f64,
    /// Contact nodes of the upper and lower wall on the slab.
    contact: [Vec<usize>; 2],
    direction: DVec3,
}

/// An upper-wall node that takes the vertex of a lower-wall node: after the
/// alignment both lie on one junction line, less than the minimum edge
/// length apart. Their offset was created by the alignment, not the source.
#[derive(Debug, Clone, Serialize)]
pub struct Identified {
    pub upper_node: u32,
    pub lower_node: u32,
    /// Distance along the junction line.
    pub distance: f64,
}

fn normal(p: &PlaneFrame) -> DVec3 {
    DVec3::from_array(p.normal)
}

/// Largest distance of a patch's candidate points to a plane.
fn offset(source: &frame::Report, patch: usize, plane: &PlaneFrame) -> f64 {
    source.surfaces[patch]
        .nodes
        .iter()
        .map(|&n| plane.distance(source.candidate_points[n]).abs())
        .fold(0.0_f64, f64::max)
}

pub(super) fn align(
    mesh: &MeshData,
    source: &frame::Report,
    policy: &Policy,
    owners: &BTreeMap<u32, Vec<usize>>,
    representatives: &mut [usize],
    tolerance: f64,
) -> Report {
    let mut report = Report {
        tolerance,
        ..Default::default()
    };
    let up = DVec3::from_array(source.policy.up).normalize();
    let planes = &source.candidate_planes;
    let count = source.surfaces.len();
    let nodes: Vec<BTreeSet<usize>> = source
        .surfaces
        .iter()
        .map(|s| s.nodes.iter().copied().collect())
        .collect();
    let (sin, cos) = source.policy.angle.sin_cos();
    let wall = |i: usize| normal(&planes[i]).dot(up).abs() <= sin;
    let slab = |i: usize| normal(&planes[i]).dot(up).abs() >= cos;
    let point = |n: usize| DVec3::from_array(source.candidate_points[n]);
    let mut links = BTreeMap::<usize, Link>::new();
    for s in (0..count).filter(|&s| slab(s)) {
        let ns = normal(&planes[s]) * normal(&planes[s]).dot(up).signum();
        let slab_points: Vec<DVec3> = source.surfaces[s].nodes.iter().map(|&n| point(n)).collect();
        // Geometric contact: a wall resting on a slab need not share source
        // nodes with it. A node counts if it lies on the slab plane and
        // within three median node spacings of a slab node.
        let nearest = |p: DVec3, skip: Option<usize>| {
            slab_points
                .iter()
                .enumerate()
                .filter(|&(k, _)| Some(k) != skip)
                .map(|(_, q)| p.distance(*q))
                .fold(f64::INFINITY, f64::min)
        };
        let step = slab_points.len().div_ceil(64).max(1);
        let mut spacings: Vec<f64> = (0..slab_points.len())
            .step_by(step)
            .map(|k| nearest(slab_points[k], Some(k)))
            .filter(|d| d.is_finite())
            .collect();
        if spacings.is_empty() {
            continue;
        }
        spacings.sort_by(f64::total_cmp);
        let reach = 3. * spacings[spacings.len() / 2];
        // Walls in line contact with the slab, with their side and extent.
        let mut touching = vec![];
        for w in (0..count).filter(|&w| wall(w)) {
            let (contact, rest): (Vec<usize>, Vec<usize>) = nodes[w].iter().partition(|&&n| {
                nodes[s].contains(&n)
                    || (planes[s].distance(point(n).to_array()).abs() <= policy.closure_tolerance
                        && nearest(point(n), None) <= reach)
            });
            if contact.len() < 2 || rest.is_empty() {
                continue;
            }
            let side = rest
                .iter()
                .map(|&n| planes[s].distance(point(n).to_array()) * ns.dot(normal(&planes[s])))
                .sum::<f64>()
                / rest.len() as f64;
            let t = normal(&planes[w]).cross(up).normalize();
            let (lo, hi) =
                contact
                    .iter()
                    .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &n| {
                        let x = point(n).dot(t);
                        (lo.min(x), hi.max(x))
                    });
            touching.push((w, side, t, lo, hi, contact));
        }
        for (u, side_u, t, lo_u, hi_u, contact_u) in &touching {
            let (u, side_u, t, lo_u, hi_u) = (*u, *side_u, *t, *lo_u, *hi_u);
            if side_u <= 0. {
                continue;
            }
            for (l, side_l, _, _, _, contact_l) in &touching {
                let (l, side_l) = (*l, *side_l);
                if side_l >= 0. || representatives[u] == representatives[l] {
                    continue;
                }
                if normal(&planes[u]).dot(normal(&planes[l])).abs() < cos {
                    continue;
                }
                let d = offset(source, u, &planes[l]);
                if d <= policy.closure_tolerance || d > tolerance {
                    continue;
                }
                // Overlap along the junction, in the upper wall's direction.
                let (lo_l, hi_l) =
                    contact_l
                        .iter()
                        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &n| {
                            let x = point(n).dot(t);
                            (lo.min(x), hi.max(x))
                        });
                let overlap = hi_u.min(hi_l) - lo_u.max(lo_l);
                if overlap <= policy.minimum_edge {
                    continue;
                }
                let better = links
                    .get(&u)
                    .is_none_or(|k| overlap > k.overlap || (overlap == k.overlap && d < k.offset));
                if better {
                    links.insert(
                        u,
                        Link {
                            lower: l,
                            slab: s,
                            overlap,
                            offset: d,
                            contact: [contact_u.clone(), contact_l.clone()],
                            direction: t,
                        },
                    );
                }
            }
        }
    }
    // Lowest walls first, so every root keeps its own plane.
    let height =
        |i: usize| nodes[i].iter().map(|&n| point(n).dot(up)).sum::<f64>() / nodes[i].len() as f64;
    let anchors: BTreeSet<usize> = source
        .axes
        .iter()
        .flat_map(|a| a.anchors.iter().map(|x| x.node))
        .collect();
    let mut uppers: Vec<usize> = links.keys().copied().collect();
    uppers.sort_by(|&a, &b| height(a).total_cmp(&height(b)).then(a.cmp(&b)));
    for u in uppers {
        let link = &links[&u];
        let mut root = link.lower;
        let mut guard = BTreeSet::from([u]);
        while let Some(next) = links.get(&root) {
            if !guard.insert(root) {
                break;
            }
            root = next.lower;
        }
        let target = representatives[root];
        let moved = representatives[u];
        let mut entry = StackedWall {
            upper: u,
            lower: link.lower,
            root,
            slab: link.slab,
            offset: offset(source, u, &planes[target]),
            source_elements: source.surfaces[u].source_elements.clone(),
            reason: String::new(),
        };
        if entry.offset > tolerance {
            entry.reason = "stack_offset_exceeds_tolerance".into();
            report.kept.push(entry);
            continue;
        }
        // Every vertex of the moved support must close within its limits.
        let remap = |s: usize| {
            if representatives[s] == moved {
                target
            } else {
                representatives[s]
            }
        };
        let fits = owners.iter().all(|(&id, supports)| {
            if !supports.iter().any(|&s| representatives[s] == moved) {
                return true;
            }
            let Some(i) = source.node_ids.iter().position(|&x| x == id) else {
                return false;
            };
            let p = point(i);
            let chosen: BTreeSet<usize> = supports.iter().map(|&s| remap(s)).collect();
            let refs: Vec<&PlaneFrame> = chosen.iter().map(|&s| &planes[s]).collect();
            let Some(q) = intersection(p, &refs, policy.precision) else {
                return false;
            };
            p.distance(q) <= policy.junction_movement_limit
                && q.distance(mesh.nodes[&id])
                    <= movement_budget(mesh, source, i) + policy.precision
        });
        if !fits {
            entry.reason = "vertex_movement_limit".into();
            report.kept.push(entry);
            continue;
        }
        for r in representatives.iter_mut() {
            if *r == moved {
                *r = target;
            }
        }
        // One-to-one identification of nearly coincident junction nodes.
        let boundary = |n: usize| owners.contains_key(&source.node_ids[n]) && !anchors.contains(&n);
        let mut pairs: Vec<(f64, usize, usize)> = vec![];
        for &a in link.contact[0].iter().filter(|&&n| boundary(n)) {
            for &b in link.contact[1].iter().filter(|&&n| boundary(n)) {
                if a == b || (nodes[link.slab].contains(&a) && nodes[link.slab].contains(&b)) {
                    continue;
                }
                let d = (point(a) - point(b)).dot(link.direction).abs();
                if d < policy.minimum_edge {
                    pairs.push((d, a, b));
                }
            }
        }
        pairs.sort_by(|x, y| x.0.total_cmp(&y.0));
        let mut used = BTreeSet::new();
        for (d, a, b) in pairs {
            if used.insert(a) && used.insert(b) {
                report.identified.push(Identified {
                    upper_node: source.node_ids[a],
                    lower_node: source.node_ids[b],
                    distance: d,
                });
            }
        }
        entry.reason = "aligned_to_lowest_wall".into();
        report.aligned.push(entry);
    }
    report
}
