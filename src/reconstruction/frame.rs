//! Joint proposal from immutable source nodes, whole axes and support planes.
//! This is a geometric constraint solve, not a mechanical or meshing model.
use super::{planes, recognize, PlaneFrame};
use crate::input::MeshData;
use glam::DVec3;
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Serialize)]
pub struct Policy {
    pub up: [f64; 3],
    pub angle: f64,
    pub maximum_movement: f64,
    pub relative_movement: f64,
    pub minimum_length: f64,
    pub residual_tolerance: f64,
    pub iterations: usize,
}
#[derive(Debug, Clone, Serialize)]
pub struct Anchor {
    pub node: usize,
    pub t: f64,
}
#[derive(Debug, Clone, Serialize)]
pub struct Axis {
    pub endpoints: [usize; 2],
    pub anchors: Vec<Anchor>,
    pub spans: Vec<recognize::SourceSpan>,
}
#[derive(Debug, Clone, Serialize)]
pub struct Surface {
    pub plane: PlaneFrame,
    pub nodes: Vec<usize>,
    pub source_elements: Vec<u32>,
    pub stiffness_regions: BTreeMap<u32, Vec<u32>>,
}
#[derive(Debug, Clone, Serialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum ConstraintOrigin {
    AxisAnchor {
        axis: usize,
        node_id: u32,
        component: usize,
    },
    ShortAxisVector {
        axis: usize,
        component: usize,
    },
    AxisDirection {
        axis: usize,
    },
    PlaneIncidence {
        plane: usize,
        node_id: u32,
    },
}
#[derive(Debug, Clone, Serialize)]
pub struct ConstraintFailure {
    pub origin: ConstraintOrigin,
    pub residual: f64,
}
#[derive(Debug, Clone, Serialize)]
pub struct MovementFailure {
    pub node_id: u32,
    pub movement: f64,
    pub budget: f64,
}
#[derive(Debug, Clone, Serialize)]
pub struct AxisFailure {
    pub axis: usize,
    pub original_length: f64,
    pub candidate_length: f64,
    pub source_elements: Vec<u32>,
}

#[derive(Debug, Clone, Serialize)]
pub struct Report {
    pub policy: Policy,
    pub accepted: bool,
    pub candidate_constraints_satisfied: bool,
    pub violating_equations: usize,
    /// At most 100 largest residuals; indices reference axes/surfaces in this report.
    pub largest_constraint_failures: Vec<ConstraintFailure>,
    pub movement_failures: Vec<MovementFailure>,
    pub axis_failures: Vec<AxisFailure>,
    pub reason: String,
    pub iterations: usize,
    pub candidate_max_residual: f64,
    pub candidate_maximum_movement: f64,
    pub candidate_over_budget_node_ids: Vec<u32>,
    pub maximum_movement: f64,
    pub node_ids: Vec<u32>,
    pub reference_points: Vec<[f64; 3]>,
    /// Unaccepted proposal for topology assembly; never an export-ready model.
    pub candidate_points: Vec<[f64; 3]>,
    pub candidate_planes: Vec<PlaneFrame>,
    pub points: Vec<[f64; 3]>,
    pub axes: Vec<Axis>,
    pub surfaces: Vec<Surface>,
    pub plane_families: Vec<Vec<usize>>,
    pub equation_count: usize,
    pub short_axis_indices: Vec<usize>,
}
struct Equation {
    origin: ConstraintOrigin,
    terms: Vec<(usize, f64)>,
    target: f64,
}
impl Equation {
    fn value(&self, x: &[f64]) -> f64 {
        self.terms.iter().map(|&(i, a)| a * x[i]).sum()
    }
    fn residual(&self, x: &[f64]) -> f64 {
        self.value(x) - self.target
    }
}
fn equation(terms: Vec<(usize, f64)>, origin: ConstraintOrigin) -> Equation {
    let mut combined = BTreeMap::<usize, f64>::new();
    for (i, a) in terms {
        *combined.entry(i).or_default() += a;
    }
    Equation {
        target: 0.0,
        origin,
        terms: combined
            .into_iter()
            .filter(|(_, a)| a.abs() > 1e-14)
            .collect(),
    }
}

/// Near-coplanar patches may meet at a vertex rather than a full source edge.
/// Validate the entire connected family, never chain unbounded plane drift.
fn support_families(mesh: &MeshData, report: &planes::Report) -> Vec<Vec<usize>> {
    let mut owners = BTreeMap::<u32, Vec<usize>>::new();
    for (i, p) in report.patches.iter().enumerate() {
        for &n in &p.source_nodes {
            owners.entry(n).or_default().push(i);
        }
    }
    let mut neighbors = vec![BTreeSet::new(); report.patches.len()];
    for ids in owners.values() {
        for (i, &a) in ids.iter().enumerate() {
            for &b in &ids[i + 1..] {
                let x = &report.patches[a];
                let y = &report.patches[b];
                if DVec3::from_array(x.plane.normal)
                    .dot(DVec3::from_array(y.plane.normal))
                    .abs()
                    >= report.policy.angle.cos()
                    && x.source_nodes.iter().all(|n| {
                        y.plane.distance(mesh.nodes[n].to_array()).abs() <= report.policy.distance
                    })
                    && y.source_nodes.iter().all(|n| {
                        x.plane.distance(mesh.nodes[n].to_array()).abs() <= report.policy.distance
                    })
                {
                    neighbors[a].insert(b);
                    neighbors[b].insert(a);
                }
            }
        }
    }
    let mut seen = BTreeSet::new();
    let mut families = vec![];
    for seed in 0..report.patches.len() {
        if seen.contains(&seed) {
            continue;
        }
        let mut stack = vec![seed];
        let mut group = BTreeSet::new();
        while let Some(i) = stack.pop() {
            if group.insert(i) {
                stack.extend(&neighbors[i]);
            }
        }
        seen.extend(group.iter().copied());
        let nodes: BTreeSet<_> = group
            .iter()
            .flat_map(|&i| report.patches[i].source_nodes.iter().copied())
            .collect();
        let points: Vec<_> = nodes.iter().map(|id| mesh.nodes[id]).collect();
        let fitted = planes::fit(
            &points,
            DVec3::from_array(report.patches[seed].plane.normal),
        );
        let valid = fitted.is_some_and(|(n, _)| {
            let center = points.iter().copied().sum::<DVec3>() / points.len() as f64;
            points
                .iter()
                .all(|p| n.dot(*p - center).abs() <= report.policy.distance)
                && group.iter().all(|&i| {
                    n.dot(DVec3::from_array(report.patches[i].plane.normal))
                        .abs()
                        >= report.policy.angle.cos()
                })
        });
        if valid {
            families.push(group.into_iter().collect());
        } else {
            families.extend(group.into_iter().map(|i| vec![i]));
        }
    }
    families
}

pub fn solve(
    mesh: &MeshData,
    axes: &recognize::Report,
    planes: &planes::Report,
    policy: &Policy,
) -> Result<Report, &'static str> {
    let up = DVec3::from_array(policy.up);
    if !up.is_finite()
        || !up.length().is_finite()
        || up.length() == 0.0
        || policy.iterations == 0
        || !policy.angle.is_finite()
        || policy.angle <= 0.0
        || policy.angle >= std::f64::consts::FRAC_PI_4
        || [
            policy.maximum_movement,
            policy.relative_movement,
            policy.minimum_length,
            policy.residual_tolerance,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v <= 0.0)
    {
        return Err("invalid frame policy");
    }
    let up = up.normalize();
    let helper = if up.z.abs() < 0.9 { DVec3::Z } else { DVec3::X };
    let u = helper.cross(up).normalize();
    let v = up.cross(u);
    let ids: BTreeSet<_> = axes
        .axes
        .iter()
        .flat_map(|a| a.anchors.iter().map(|c| c.node))
        .chain(
            planes
                .patches
                .iter()
                .flat_map(|p| p.source_nodes.iter().copied()),
        )
        .collect();
    let node_ids: Vec<_> = ids.into_iter().collect();
    let reference: Option<Vec<_>> = node_ids
        .iter()
        .map(|id| mesh.nodes.get(id).copied().filter(|p| p.is_finite()))
        .collect();
    let reference = reference.ok_or("missing source node")?;
    let map: BTreeMap<_, _> = node_ids
        .iter()
        .enumerate()
        .map(|(i, &id)| (id, i))
        .collect();
    let center = reference.first().copied().unwrap_or(DVec3::ZERO);
    let mut x: Vec<_> = reference
        .iter()
        .flat_map(|p| (*p - center).to_array())
        .collect();
    let base = x.len();
    let mut equations = vec![];
    let mut budgets = vec![policy.maximum_movement; reference.len()];
    let mut new_axes = vec![];
    for (axis_index, a) in axes.axes.iter().enumerate() {
        if a.endpoint_nodes.iter().any(|n| !map.contains_key(n))
            || a.anchors
                .iter()
                .any(|c| !c.t.is_finite() || c.t < 0.0 || c.t > 1.0)
        {
            return Err("invalid axis anchors");
        }
        let ends = [map[&a.endpoint_nodes[0]], map[&a.endpoint_nodes[1]]];
        let d = reference[ends[1]] - reference[ends[0]];
        let length = d.length();
        if !length.is_finite() || length <= policy.residual_tolerance {
            return Err("source axis shorter than minimum length");
        }
        let cap = policy
            .maximum_movement
            .min(policy.relative_movement * length);
        let mut anchors = vec![];
        for c in &a.anchors {
            let node = map[&c.node];
            budgets[node] = budgets[node].min(cap);
            for k in 0..3 {
                equations.push(equation(
                    vec![
                        (node * 3 + k, 1.),
                        (ends[0] * 3 + k, -(1. - c.t)),
                        (ends[1] * 3 + k, -c.t),
                    ],
                    ConstraintOrigin::AxisAnchor {
                        axis: axis_index,
                        node_id: c.node,
                        component: k,
                    },
                ));
            }
            anchors.push(Anchor { node, t: c.t });
        }
        let cosine = d.normalize().dot(up).abs();
        // An unresolved short feature may translate, but must not be flattened,
        // contracted or assigned a new direction from noisy source coordinates.
        if length < policy.minimum_length {
            for k in 0..3 {
                let mut e = equation(
                    vec![(ends[1] * 3 + k, 1.), (ends[0] * 3 + k, -1.)],
                    ConstraintOrigin::ShortAxisVector {
                        axis: axis_index,
                        component: k,
                    },
                );
                e.target = d[k];
                equations.push(e);
            }
        }
        let directions = if length < policy.minimum_length {
            vec![]
        } else if cosine >= policy.angle.cos() {
            vec![u, v]
        } else if cosine <= policy.angle.sin() {
            vec![up]
        } else {
            vec![]
        };
        for n in directions {
            equations.push(equation(
                (0..3)
                    .flat_map(|k| [(ends[0] * 3 + k, n[k]), (ends[1] * 3 + k, -n[k])])
                    .collect(),
                ConstraintOrigin::AxisDirection { axis: axis_index },
            ));
        }
        new_axes.push(Axis {
            endpoints: ends,
            anchors,
            spans: a.spans.clone(),
        });
    }
    let plane_families = support_families(mesh, planes);
    let mut plane_to_family = vec![0; planes.patches.len()];
    let mut normals = vec![];
    for (family, members) in plane_families.iter().enumerate() {
        let nodes: BTreeSet<_> = members
            .iter()
            .flat_map(|&i| planes.patches[i].source_nodes.iter().copied())
            .collect();
        if nodes.is_empty() {
            return Err("empty support plane");
        }
        let points: Vec<_> = nodes.iter().map(|n| mesh.nodes[n]).collect();
        let old = if members.len() == 1 {
            DVec3::from_array(planes.patches[members[0]].plane.normal)
        } else {
            planes::fit(
                &points,
                DVec3::from_array(planes.patches[members[0]].plane.normal),
            )
            .ok_or("invalid support family")?
            .0
        };
        let cosine = old.dot(up);
        let normal = if cosine.abs() >= policy.angle.cos() {
            up * cosine.signum()
        } else if cosine.abs() <= policy.angle.sin() {
            (old - up * cosine).normalize()
        } else {
            old
        };
        let origin = points.iter().map(|p| *p - center).sum::<DVec3>() / points.len() as f64;
        x.push(normal.dot(origin));
        normals.push(normal);
        for &pi in members {
            plane_to_family[pi] = family;
        }
    }
    let mut surfaces = vec![];
    for (pi, p) in planes.patches.iter().enumerate() {
        let family = plane_to_family[pi];
        let normal = normals[family];
        let nodes: Vec<_> = p.source_nodes.iter().map(|id| map[id]).collect();
        for &node in &nodes {
            equations.push(equation(
                vec![
                    (node * 3, normal.x),
                    (node * 3 + 1, normal.y),
                    (node * 3 + 2, normal.z),
                    (base + family, -1.),
                ],
                ConstraintOrigin::PlaneIncidence {
                    plane: pi,
                    node_id: node_ids[node],
                },
            ));
        }
        surfaces.push(Surface {
            plane: p.plane.clone(),
            nodes,
            source_elements: p.source_elements.clone(),
            stiffness_regions: p.stiffness_regions.clone(),
        });
    }
    equations.retain(|e| !e.terms.is_empty());
    // LSQR operates on A and A^T directly; avoid normal equations A A^T
    // whose conditioning is squared for nearly dependent plane constraints.
    let norm = |v: &[f64]| v.iter().fold(0.0_f64, |a, b| a.hypot(*b));
    let scales: Vec<_> = equations
        .iter()
        .map(|e| e.terms.iter().map(|(_, a)| a * a).sum::<f64>().sqrt())
        .collect();
    let multiply = |v: &[f64]| -> Vec<f64> {
        equations
            .iter()
            .zip(&scales)
            .map(|(e, s)| e.value(v) / s)
            .collect()
    };
    let transpose = |v: &[f64]| -> Vec<f64> {
        let mut result = vec![0.0; x.len()];
        for ((e, s), v) in equations.iter().zip(&scales).zip(v) {
            for &(i, a) in &e.terms {
                result[i] += a * v / s;
            }
        }
        result
    };
    let mut u: Vec<_> = equations
        .iter()
        .zip(&scales)
        .map(|(e, s)| -e.residual(&x) / s)
        .collect();
    let mut beta = norm(&u);
    if beta > 0. {
        for v in &mut u {
            *v /= beta;
        }
    }
    let mut v = transpose(&u);
    let mut alpha = norm(&v);
    if alpha > 0. {
        for p in &mut v {
            *p /= alpha;
        }
    }
    let mut w = v.clone();
    let mut rho_bar = alpha;
    let mut phi_bar = beta;
    let mut correction = vec![0.0; x.len()];
    let original = x.clone();
    let mut residual = equations
        .iter()
        .map(|e| e.residual(&x).abs())
        .fold(0.0_f64, f64::max);
    let mut iterations = 0;
    for step in 0..policy.iterations {
        if residual <= policy.residual_tolerance {
            break;
        }
        let av = multiply(&v);
        for (u, av) in u.iter_mut().zip(av) {
            *u = av - alpha * (*u);
        }
        beta = norm(&u);
        if beta > 0. {
            for p in &mut u {
                *p /= beta;
            }
        }
        let atu = transpose(&u);
        for (v, atu) in v.iter_mut().zip(atu) {
            *v = atu - beta * (*v);
        }
        alpha = norm(&v);
        if alpha > 0. {
            for p in &mut v {
                *p /= alpha;
            }
        }
        let rho = rho_bar.hypot(beta);
        if !rho.is_finite() || rho <= 1e-30 {
            break;
        }
        let c = rho_bar / rho;
        let s = beta / rho;
        let theta = s * alpha;
        rho_bar = -c * alpha;
        let phi = c * phi_bar;
        phi_bar *= s;
        for ((dx, w), v) in correction.iter_mut().zip(&mut w).zip(&v) {
            *dx += (phi / rho) * (*w);
            *w = *v - (theta / rho) * (*w);
        }
        // No mutation of x until closure borrows end; assess true residual.
        let candidate: Vec<_> = original
            .iter()
            .zip(&correction)
            .map(|(a, b)| a + b)
            .collect();
        residual = equations
            .iter()
            .map(|e| e.residual(&candidate).abs())
            .fold(0.0_f64, f64::max);
        iterations = step + 1;
    }
    for ((x, o), d) in x.iter_mut().zip(original).zip(correction) {
        *x = o + d;
    }
    let candidate: Vec<_> = (0..reference.len())
        .map(|i| center + DVec3::new(x[i * 3], x[i * 3 + 1], x[i * 3 + 2]))
        .collect();
    let valid = candidate.iter().all(|p| p.is_finite())
        && new_axes.iter().all(|a| {
            let [i, j] = a.endpoints;
            let d = candidate[j] - candidate[i];
            d.length() + policy.residual_tolerance
                >= policy
                    .minimum_length
                    .min(reference[j].distance(reference[i]))
                && d.dot(reference[j] - reference[i]) > 0.0
        });
    let within_budget = candidate
        .iter()
        .zip(&reference)
        .zip(&budgets)
        .all(|((p, o), cap)| p.distance(*o) <= *cap + 1e-10);
    let accepted = residual <= policy.residual_tolerance && valid && within_budget;
    let mut maximum_movement = 0.0_f64;
    let mut candidate_surfaces = surfaces.clone();
    {
        for (i, p) in candidate_surfaces.iter_mut().enumerate() {
            let origin = DVec3::from_array(p.plane.origin);
            let family = plane_to_family[i];
            let n = normals[family];
            let projected = origin + n * (x[base + family] - n.dot(origin - center));
            p.plane = PlaneFrame::new(projected.to_array(), n.to_array())
                .map_err(|_| "invalid solved plane")?;
        }
    }
    if accepted {
        surfaces = candidate_surfaces.clone();
        maximum_movement = candidate
            .iter()
            .zip(&reference)
            .map(|(a, b)| a.distance(*b))
            .fold(0.0_f64, f64::max);
    }
    let mut failures: Vec<_> = equations
        .iter()
        .filter_map(|e| {
            let residual = e.residual(&x).abs();
            (residual > policy.residual_tolerance).then(|| ConstraintFailure {
                origin: e.origin.clone(),
                residual,
            })
        })
        .collect();
    let violating_equations = failures.len();
    failures.sort_by(|a, b| b.residual.total_cmp(&a.residual));
    failures.truncate(100);
    let movement_failures = candidate
        .iter()
        .zip(&reference)
        .zip(&budgets)
        .enumerate()
        .filter_map(|(i, ((p, o), cap))| {
            let movement = p.distance(*o);
            (movement > *cap + 1e-10).then_some(MovementFailure {
                node_id: node_ids[i],
                movement,
                budget: *cap,
            })
        })
        .collect();
    let axis_failures = new_axes
        .iter()
        .enumerate()
        .filter_map(|(axis, a)| {
            let [i, j] = a.endpoints;
            let old = reference[j] - reference[i];
            let new = candidate[j] - candidate[i];
            (new.length() + policy.residual_tolerance < policy.minimum_length.min(old.length())
                || new.dot(old) <= 0.0)
                .then_some(AxisFailure {
                    axis,
                    original_length: old.length(),
                    candidate_length: new.length(),
                    source_elements: a.spans.iter().map(|s| s.element).collect(),
                })
        })
        .collect();
    Ok(Report {
        policy: policy.clone(),
        accepted,
        candidate_constraints_satisfied: residual <= policy.residual_tolerance,
        violating_equations,
        largest_constraint_failures: failures,
        movement_failures,
        axis_failures,
        reason: if accepted {
            "converged"
        } else if !within_budget {
            "movement_budget_exceeded"
        } else if !valid {
            "invalid_axis"
        } else {
            "constraints_or_budget_not_satisfied"
        }
        .into(),
        iterations,
        candidate_max_residual: residual,
        candidate_maximum_movement: candidate
            .iter()
            .zip(&reference)
            .map(|(a, b)| a.distance(*b))
            .fold(0.0_f64, f64::max),
        candidate_over_budget_node_ids: candidate
            .iter()
            .zip(&reference)
            .zip(&budgets)
            .enumerate()
            .filter_map(|(i, ((p, o), cap))| (p.distance(*o) > *cap + 1e-10).then_some(node_ids[i]))
            .collect(),
        maximum_movement,
        node_ids: node_ids.clone(),
        reference_points: reference.iter().map(|p| p.to_array()).collect(),
        candidate_points: candidate.iter().map(|p| p.to_array()).collect(),
        candidate_planes: candidate_surfaces.iter().map(|s| s.plane.clone()).collect(),
        points: if accepted {
            candidate.iter().map(|p| p.to_array()).collect()
        } else {
            reference.iter().map(|p| p.to_array()).collect()
        },
        axes: new_axes.clone(),
        surfaces,
        short_axis_indices: new_axes
            .iter()
            .enumerate()
            .filter_map(|(i, a)| {
                (reference[a.endpoints[0]].distance(reference[a.endpoints[1]])
                    < policy.minimum_length)
                    .then_some(i)
            })
            .collect(),
        plane_families,
        equation_count: equations.len(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::input::ElementData;
    fn source() -> MeshData {
        let points = [
            [0., 0., 0.02],
            [1., 0., 0.025],
            [2., 0., 0.03],
            [1.01, 0., 2.],
            [0., 2., 0.02],
        ];
        MeshData {
            nodes: points
                .iter()
                .enumerate()
                .map(|(i, p)| (i as u32 + 1, DVec3::from_array(*p)))
                .collect(),
            elements: vec![
                ElementData {
                    id: 1,
                    elem_type: 10,
                    stiff_id: 1,
                    nodes: vec![1, 2],
                },
                ElementData {
                    id: 2,
                    elem_type: 10,
                    stiff_id: 2,
                    nodes: vec![2, 3],
                },
                ElementData {
                    id: 3,
                    elem_type: 10,
                    stiff_id: 3,
                    nodes: vec![2, 4],
                },
                ElementData {
                    id: 4,
                    elem_type: 42,
                    stiff_id: 4,
                    nodes: vec![1, 3, 5],
                },
            ],
        }
    }
    fn policy() -> Policy {
        Policy {
            up: [0., 0., 1.],
            angle: 0.02,
            maximum_movement: 0.15,
            relative_movement: 0.05,
            minimum_length: 0.03,
            residual_tolerance: 1e-8,
            iterations: 2000,
        }
    }
    fn run(m: &MeshData, p: &Policy) -> Report {
        let a = recognize::recognize(
            m,
            &recognize::Policy {
                angle: 0.02,
                line_tolerance: 0.01,
                numerical_precision: 1e-8,
            },
        )
        .unwrap();
        let s = planes::recognize(
            m,
            &planes::Policy {
                angle: 0.02,
                distance: 0.01,
                precision: 1e-8,
            },
        )
        .unwrap();
        solve(m, &a, &s, p).unwrap()
    }
    fn short_feature(both_ends_on_plane: bool) -> MeshData {
        let mut m = MeshData::default();
        for (id, point) in [
            (1, [0., 0., 0.]),
            (2, [1., 0., 0.0002]),
            (3, [0., 1., 0.]),
            (4, [0.006, 0., 0.00001]),
        ] {
            m.nodes.insert(id, DVec3::from_array(point));
        }
        m.elements.push(ElementData {
            id: 10,
            elem_type: 10,
            stiff_id: 77,
            nodes: vec![1, 4],
        });
        m.elements.push(ElementData {
            id: 20,
            elem_type: 42,
            stiff_id: 88,
            nodes: vec![1, 2, 3],
        });
        if both_ends_on_plane {
            m.elements.push(ElementData {
                id: 21,
                elem_type: 42,
                stiff_id: 88,
                nodes: vec![1, 3, 4],
            });
        }
        m
    }

    #[test]
    fn short_axis_translates_with_support_without_losing_vector_or_properties() {
        for transformed in [false, true] {
            let mut m = short_feature(false);
            let mut p = policy();
            if transformed {
                let rotation =
                    glam::DQuat::from_axis_angle(DVec3::new(1., 2., 3.).normalize(), 0.6);
                m.nodes = m
                    .nodes
                    .into_iter()
                    .map(|(id, point)| (id + 100, rotation * point + DVec3::new(20., 30., 40.)))
                    .collect();
                for e in &mut m.elements {
                    e.id += 200;
                    for n in &mut e.nodes {
                        *n += 100;
                    }
                }
                p.up = (rotation * DVec3::Z).to_array();
            }
            let r = run(&m, &p);
            assert!(r.accepted, "{}: {}", r.reason, r.candidate_max_residual);
            assert_eq!(r.short_axis_indices, vec![0]);
            assert_eq!(r.axes.len(), 1);
            let [i, j] = r.axes[0].endpoints;
            let old =
                DVec3::from_array(r.reference_points[j]) - DVec3::from_array(r.reference_points[i]);
            let new = DVec3::from_array(r.points[j]) - DVec3::from_array(r.points[i]);
            assert!(old.distance(new) < 2. * p.residual_tolerance);
            assert!(r.maximum_movement > 1e-6);
            assert!(r.axis_failures.is_empty());
            assert_eq!(r.axes[0].spans.len(), 1);
            assert_eq!(r.axes[0].spans[0].stiffness, 77);
            assert_eq!(
                r.axes[0].spans[0].element,
                if transformed { 210 } else { 10 }
            );
            assert_eq!(
                r.surfaces[0]
                    .stiffness_regions
                    .keys()
                    .copied()
                    .collect::<Vec<_>>(),
                vec![88]
            );
        }
    }

    #[test]
    fn conflicting_short_axis_is_reported_without_collapsing_or_dropping_it() {
        let m = short_feature(true);
        let r = run(&m, &policy());
        assert!(!r.accepted);
        assert!(!r.candidate_constraints_satisfied);
        assert_eq!(r.points, r.reference_points);
        assert_eq!(r.axes.len(), 1);
        assert_eq!(r.axes[0].spans[0].element, 10);
        assert!(r
            .largest_constraint_failures
            .iter()
            .any(|f| matches!(f.origin, ConstraintOrigin::ShortAxisVector { .. })));
    }

    #[test]
    fn coplanar_patches_meeting_at_one_vertex_share_support() {
        let mut m = MeshData::default();
        for (i, p) in [
            [0., 0., 0.],
            [1., 0., 0.],
            [0., 1., 0.],
            [-1., 0., 0.],
            [0., -1., 0.],
        ]
        .iter()
        .enumerate()
        {
            m.nodes.insert(i as u32 + 1, DVec3::from_array(*p));
        }
        m.elements = vec![
            ElementData {
                id: 1,
                elem_type: 42,
                stiff_id: 1,
                nodes: vec![1, 2, 3],
            },
            ElementData {
                id: 2,
                elem_type: 42,
                stiff_id: 2,
                nodes: vec![1, 4, 5],
            },
        ];
        let p = planes::recognize(
            &m,
            &planes::Policy {
                angle: 0.02,
                distance: 0.01,
                precision: 1e-8,
            },
        )
        .unwrap();
        assert_eq!(p.patches.len(), 2);
        assert_eq!(support_families(&m, &p), vec![vec![0, 1]]);
        let a = recognize::recognize(
            &m,
            &recognize::Policy {
                angle: 0.02,
                line_tolerance: 0.01,
                numerical_precision: 1e-8,
            },
        )
        .unwrap();
        let result = solve(&m, &a, &p, &policy()).unwrap();
        assert!(result.accepted);
        assert_eq!(result.plane_families, vec![vec![0, 1]]);
        assert_eq!(result.surfaces[0].stiffness_regions.len(), 1);
        assert_eq!(result.surfaces[1].stiffness_regions.len(), 1);
    }

    #[test]
    fn plane_and_axes_share_interior_anchor_after_regularization() {
        let r = run(&source(), &policy());
        assert!(r.accepted, "{}", r.candidate_max_residual);
        assert_eq!(r.axes.len(), 2);
        for a in &r.axes {
            let x = DVec3::from_array(r.points[a.endpoints[0]]);
            let y = DVec3::from_array(r.points[a.endpoints[1]]);
            for c in &a.anchors {
                assert!(DVec3::from_array(r.points[c.node]).distance(x + c.t * (y - x)) < 2e-8);
            }
        }
        for s in &r.surfaces {
            for &n in &s.nodes {
                assert!(s.plane.distance(r.points[n]).abs() < 2e-8);
            }
        }
        let p = r.points[r.node_ids.iter().position(|n| *n == 2).unwrap()];
        let q = r.points[r.node_ids.iter().position(|n| *n == 4).unwrap()];
        assert!((p[0] - q[0]).abs() < 2e-8 && (p[1] - q[1]).abs() < 2e-8);
    }
    #[test]
    fn failure_returns_original_geometry_not_last_iterate() {
        let mut p = policy();
        p.maximum_movement = 1e-6;
        p.iterations = 100;
        let r = run(&source(), &p);
        assert!(!r.accepted);
        assert_eq!(r.points, r.reference_points);
        assert_eq!(r.maximum_movement, 0.);
        assert!(!r.movement_failures.is_empty());
        assert!(r.movement_failures.iter().all(|f| f.movement > f.budget));
    }
    #[test]
    fn rigid_transform_with_up_preserves_a_feasible_solution() {
        let mut m = source();
        let rotation = glam::DQuat::from_axis_angle(DVec3::new(1., 2., 3.).normalize(), 0.6);
        for v in m.nodes.values_mut() {
            *v = rotation * (*v) + DVec3::new(20., 30., 40.);
        }
        let mut p = policy();
        p.up = (rotation * DVec3::Z).to_array();
        let r = run(&m, &p);
        assert!(r.accepted);
        assert!(r.maximum_movement < 0.15);
    }
}
