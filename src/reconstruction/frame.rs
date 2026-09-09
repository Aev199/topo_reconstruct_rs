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
pub struct Report {
    pub policy: Policy,
    pub accepted: bool,
    pub reason: String,
    pub iterations: usize,
    pub candidate_max_residual: f64,
    pub candidate_maximum_movement: f64,
    pub candidate_over_budget_node_ids: Vec<u32>,
    pub maximum_movement: f64,
    pub node_ids: Vec<u32>,
    pub reference_points: Vec<[f64; 3]>,
    pub points: Vec<[f64; 3]>,
    pub axes: Vec<Axis>,
    pub surfaces: Vec<Surface>,
    pub equation_count: usize,
    pub short_axis_indices: Vec<usize>,
}
struct Equation {
    terms: Vec<(usize, f64)>,
}
impl Equation {
    fn residual(&self, x: &[f64]) -> f64 {
        self.terms.iter().map(|&(i, a)| a * x[i]).sum()
    }
}
fn equation(terms: Vec<(usize, f64)>) -> Equation {
    let mut combined = BTreeMap::<usize, f64>::new();
    for (i, a) in terms {
        *combined.entry(i).or_default() += a;
    }
    Equation {
        terms: combined
            .into_iter()
            .filter(|(_, a)| a.abs() > 1e-14)
            .collect(),
    }
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
    for a in &axes.axes {
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
                equations.push(equation(vec![
                    (node * 3 + k, 1.),
                    (ends[0] * 3 + k, -(1. - c.t)),
                    (ends[1] * 3 + k, -c.t),
                ]));
            }
            anchors.push(Anchor { node, t: c.t });
        }
        let cosine = d.normalize().dot(up).abs();
        let directions = if cosine >= policy.angle.cos() {
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
            ));
        }
        new_axes.push(Axis {
            endpoints: ends,
            anchors,
            spans: a.spans.clone(),
        });
    }
    let mut normals = vec![];
    let mut surfaces = vec![];
    for (pi, p) in planes.patches.iter().enumerate() {
        let old = DVec3::from_array(p.plane.normal);
        let cosine = old.dot(up);
        let normal = if cosine.abs() >= policy.angle.cos() {
            up * cosine.signum()
        } else if cosine.abs() <= policy.angle.sin() {
            (old - up * cosine).normalize()
        } else {
            old
        };
        let nodes: Vec<_> = p.source_nodes.iter().map(|id| map[id]).collect();
        if nodes.is_empty() {
            return Err("empty support plane");
        }
        let origin =
            nodes.iter().map(|&i| reference[i] - center).sum::<DVec3>() / nodes.len() as f64;
        x.push(normal.dot(origin));
        normals.push(normal);
        for &node in &nodes {
            equations.push(equation(vec![
                (node * 3, normal.x),
                (node * 3 + 1, normal.y),
                (node * 3 + 2, normal.z),
                (base + pi, -1.),
            ]));
        }
        surfaces.push(Surface {
            plane: p.plane.clone(),
            nodes,
            source_elements: p.source_elements.clone(),
            stiffness_regions: p.stiffness_regions.clone(),
        });
    }
    equations.retain(|e| !e.terms.is_empty());
    // Minimum-displacement projection onto the common linear constraints.
    // Preconditioned conjugate gradients on A A^T, without materializing a
    // dense matrix. Redundant source incidences are permitted.
    let mut r: Vec<_> = equations.iter().map(|e| e.residual(&x)).collect();
    let diagonal: Vec<f64> = equations
        .iter()
        .map(|e| e.terms.iter().map(|(_, a)| a * a).sum())
        .collect();
    let mut direction: Vec<_> = r.iter().zip(&diagonal).map(|(r, d)| r / d).collect();
    let mut rz: f64 = r.iter().zip(&direction).map(|(a, b)| a * b).sum();
    let mut residual = r.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
    let mut iterations = 0;
    for step in 0..policy.iterations {
        if residual <= policy.residual_tolerance {
            break;
        }
        let mut transpose = vec![0.0; x.len()];
        for (e, p) in equations.iter().zip(&direction) {
            for &(i, a) in &e.terms {
                transpose[i] += a * p;
            }
        }
        let denominator: f64 = transpose.iter().map(|v| v * v).sum();
        if !denominator.is_finite() || denominator <= 1e-30 {
            break;
        }
        let alpha = rz / denominator;
        for (xi, di) in x.iter_mut().zip(&transpose) {
            *xi -= alpha * di;
        }
        for (ri, e) in r.iter_mut().zip(&equations) {
            *ri -= alpha * e.residual(&transpose);
        }
        // True residual decides acceptance; recursive residual alone can drift.
        residual = equations
            .iter()
            .map(|e| e.residual(&x).abs())
            .fold(0.0_f64, f64::max);
        iterations = step + 1;
        let next_rz: f64 = r.iter().zip(&diagonal).map(|(r, d)| r * r / d).sum();
        if !next_rz.is_finite() || rz <= 1e-30 {
            break;
        }
        let beta = next_rz / rz;
        for ((p, r), d) in direction.iter_mut().zip(&r).zip(&diagonal) {
            *p = r / d + beta * (*p);
        }
        rz = next_rz;
    }
    let candidate: Vec<_> = (0..reference.len())
        .map(|i| center + DVec3::new(x[i * 3], x[i * 3 + 1], x[i * 3 + 2]))
        .collect();
    let valid = candidate.iter().all(|p| p.is_finite())
        && new_axes.iter().all(|a| {
            let [i, j] = a.endpoints;
            let d = candidate[j] - candidate[i];
            d.length()
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
    if accepted {
        for (i, p) in surfaces.iter_mut().enumerate() {
            let origin = DVec3::from_array(p.plane.origin);
            let n = normals[i];
            let projected = origin + n * (x[base + i] - n.dot(origin - center));
            p.plane = PlaneFrame::new(projected.to_array(), n.to_array())
                .map_err(|_| "invalid solved plane")?;
        }
        maximum_movement = candidate
            .iter()
            .zip(&reference)
            .map(|(a, b)| a.distance(*b))
            .fold(0.0_f64, f64::max);
    }
    Ok(Report {
        policy: policy.clone(),
        accepted,
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
