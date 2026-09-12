//! Conservative feature recovery on fixed shared supports. No hole filling,
//! copied junctions, dropped incidences, or model-specific exceptions.
use super::{frame, intersection, movement_budget, MeshData, Model, PlaneFrame, Policy};
use glam::DVec3;
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum HoleOutcome {
    Restored,
    CommonSupportsForceLineOrPoint,
    InvalidSourceContours,
    MissingClosedVertex,
    MovementBudget,
    AxisAnchorRequiresJointRepair,
    IncompatibleSupports,
    ContourConflict,
}
#[derive(Debug, Serialize)]
pub struct HoleConstraint {
    pub source_nodes: Vec<u32>,
    /// Representative support indices in the frame report.
    pub common_supports: Vec<usize>,
    pub normal_rank: usize,
}
#[derive(Debug, Serialize)]
pub struct HoleNodeChange {
    pub source_node: u32,
    pub before: [f64; 3],
    pub after: [f64; 3],
}
#[derive(Debug, Serialize)]
pub struct HoleRecovery {
    pub patch: usize,
    pub source_elements: Vec<u32>,
    pub outcome: HoleOutcome,
    pub constraints: Vec<HoleConstraint>,
    /// Applied changes only. Failed attempts leave shared coordinates untouched.
    pub changes: Vec<HoleNodeChange>,
}

pub(super) struct Context<'a> {
    pub mesh: &'a MeshData,
    pub source: &'a frame::Report,
    pub policy: &'a Policy,
    pub regions: &'a [(usize, u32, Vec<u32>)],
    /// Exterior/hole roles are fixed using the immutable source, not re-sorted
    /// after movement (a deformed hole must not become the exterior).
    pub rings: &'a BTreeMap<usize, Vec<Vec<u32>>>,
    pub owners: &'a BTreeMap<u32, Vec<usize>>,
    pub representatives: &'a [usize],
    pub supports: &'a [PlaneFrame],
}

fn normal_rank(indices: &[usize], supports: &[PlaneFrame]) -> usize {
    let mut basis: Vec<DVec3> = vec![];
    for &i in indices {
        let mut n = DVec3::from_array(supports[i].normal);
        for _ in 0..2 {
            for &q in &basis {
                n -= q * q.dot(n);
            }
        }
        if n.length() > 1e-10 && basis.len() < 3 {
            basis.push(n.normalize());
        }
    }
    basis.len()
}

fn valid(
    plane: &PlaneFrame,
    rings: &[Vec<u32>],
    points: &BTreeMap<u32, DVec3>,
    policy: &Policy,
) -> bool {
    let Ok(mut model) = Model::new(policy.precision, policy.minimum_edge) else {
        return false;
    };
    let p = model.add_plane(plane.clone());
    let mut vertices = BTreeMap::new();
    for &n in rings.iter().flatten() {
        let Some(q) = points.get(&n) else {
            return false;
        };
        if !vertices.contains_key(&n) {
            let Ok(v) = model.add_vertex(q.to_array()) else {
                return false;
            };
            vertices.insert(n, v);
        }
    }
    model
        .add_surface(
            p,
            rings
                .iter()
                .map(|r| r.iter().map(|n| vertices[n]).collect())
                .collect(),
            vec![],
        )
        .is_ok()
}

struct Attempt {
    report: usize,
    desired: BTreeMap<u32, DVec3>,
    affected: BTreeSet<usize>,
}

pub(super) fn recover(points: &mut BTreeMap<u32, DVec3>, c: &Context<'_>) -> Vec<HoleRecovery> {
    let lookup: BTreeMap<_, _> = c
        .source
        .node_ids
        .iter()
        .enumerate()
        .map(|(i, &n)| (n, i))
        .collect();
    let mut reports = vec![];
    let mut attempts = vec![];
    for (&region, loops) in c.rings {
        let (patch, _, ids) = &c.regions[region];
        let plane = &c.supports[c.representatives[*patch]];
        if loops.len() < 2 || valid(plane, loops, points, c.policy) {
            continue;
        }
        let constraints: Vec<_> = loops
            .iter()
            .skip(1)
            .map(|ring| {
                let mut common: Option<BTreeSet<usize>> = None;
                for n in ring {
                    let owned: BTreeSet<_> =
                        c.owners[n].iter().map(|&p| c.representatives[p]).collect();
                    common = Some(match common {
                        None => owned,
                        Some(s) => s.intersection(&owned).copied().collect(),
                    });
                }
                let indices: Vec<_> = common.unwrap_or_default().into_iter().collect();
                HoleConstraint {
                    source_nodes: ring.clone(),
                    normal_rank: normal_rank(&indices, c.supports),
                    common_supports: indices,
                }
            })
            .collect();
        let original: BTreeMap<_, _> = loops
            .iter()
            .flatten()
            .map(|&n| {
                (
                    n,
                    DVec3::from_array(plane.lift(plane.project(c.mesh.nodes[&n].to_array()))),
                )
            })
            .collect();
        let mut report = HoleRecovery {
            patch: *patch,
            source_elements: ids.clone(),
            outcome: HoleOutcome::ContourConflict,
            constraints,
            changes: vec![],
        };
        let mut desired = BTreeMap::new();
        let failure = if !valid(plane, loops, &original, c.policy) {
            Some(HoleOutcome::InvalidSourceContours)
        } else if report.constraints.iter().any(|h| h.normal_rank >= 2) {
            // Two independent planes common to the entire hole restrict it to
            // a line (three to a point). Moving plane offsets cannot fix this.
            Some(HoleOutcome::CommonSupportsForceLineOrPoint)
        } else {
            loops.iter().skip(1).flatten().find_map(|&n| {
                let Some(&current) = points.get(&n) else {
                    return Some(HoleOutcome::MissingClosedVertex);
                };
                let owned: Vec<_> = c.owners[&n]
                    .iter()
                    .map(|&p| &c.supports[c.representatives[p]])
                    .collect();
                let reference = c.mesh.nodes[&n];
                let Some(q) = intersection(reference, &owned, c.policy.precision) else {
                    return Some(HoleOutcome::IncompatibleSupports);
                };
                let i = lookup[&n];
                if q.distance(DVec3::from_array(c.source.candidate_points[i]))
                    > c.policy.junction_movement_limit
                    || q.distance(reference)
                        > movement_budget(c.mesh, c.source, i) + c.policy.precision
                {
                    return Some(HoleOutcome::MovementBudget);
                }
                // Moving an axis anchor in isolation could bend the whole axis.
                // Leave such repairs to the joint axis/surface stage.
                if q.distance(current) > c.policy.precision
                    && c.source
                        .axes
                        .iter()
                        .any(|a| a.anchors.iter().any(|a| a.node == i))
                {
                    return Some(HoleOutcome::AxisAnchorRequiresJointRepair);
                }
                desired.insert(n, q);
                None
            })
        };
        if let Some(reason) = failure {
            report.outcome = reason;
        } else {
            let affected = c
                .rings
                .iter()
                .filter_map(|(&r, ls)| {
                    ls.iter()
                        .flatten()
                        .any(|n| desired.contains_key(n))
                        .then_some(r)
                })
                .collect();
            attempts.push(Attempt {
                report: reports.len(),
                desired,
                affected,
            });
        }
        reports.push(report);
    }
    // Proposals whose changed vertices affect a common contour are validated
    // and committed together. Disjoint failures cannot block other components;
    // traversal order / source numbering never decides which proposal wins.
    let mut remaining: BTreeSet<_> = (0..attempts.len()).collect();
    while let Some(&seed) = remaining.first() {
        remaining.remove(&seed);
        let mut group = vec![seed];
        let mut affected = attempts[seed].affected.clone();
        loop {
            let connected: Vec<_> = remaining
                .iter()
                .copied()
                .filter(|&i| !affected.is_disjoint(&attempts[i].affected))
                .collect();
            if connected.is_empty() {
                break;
            }
            for i in connected {
                remaining.remove(&i);
                affected.extend(&attempts[i].affected);
                group.push(i);
            }
        }
        let mut proposed = points.clone();
        for &i in &group {
            proposed.extend(attempts[i].desired.iter().map(|(&n, &p)| (n, p)));
        }
        if !affected.iter().all(|&region| {
            let plane = &c.supports[c.representatives[c.regions[region].0]];
            valid(plane, &c.rings[&region], &proposed, c.policy)
        }) {
            continue;
        }
        for i in group {
            let report = &mut reports[attempts[i].report];
            report.outcome = HoleOutcome::Restored;
            report.changes = attempts[i]
                .desired
                .iter()
                .filter_map(|(&n, &q)| {
                    (q != points[&n]).then_some(HoleNodeChange {
                        source_node: n,
                        before: points[&n].to_array(),
                        after: q.to_array(),
                    })
                })
                .collect();
        }
        *points = proposed;
    }
    reports
}
