//! Damped Gauss–Newton using the existing LSQR kernel, in one shared system.
use super::*;

struct Parameter {
    axis: usize,
    anchor: usize,
    variable: usize,
    scale: f64,
    rows: Vec<(usize, usize)>,
}

fn linearize(equations: &mut [Equation], x: &[f64], axes: &[Axis], parameters: &[Parameter]) {
    for p in parameters {
        let axis = &axes[p.axis];
        let node = axis.anchors[p.anchor].node;
        let [a, b] = axis.endpoints;
        let t = x[p.variable] / p.scale;
        for &(row, k) in &p.rows {
            let d = x[b * 3 + k] - x[a * 3 + k];
            equations[row].terms = vec![
                (node * 3 + k, 1.),
                (a * 3 + k, -(1. - t)),
                (b * 3 + k, -t),
                (p.variable, -d / p.scale),
            ];
            equations[row].target = -d * t;
        }
    }
}

pub(super) fn solve(
    equations: &mut Vec<Equation>,
    x: &mut Vec<f64>,
    axes: &[Axis],
    reference: &[DVec3],
    node_map: &BTreeMap<u32, usize>,
    policy: &Policy,
    maximum_steps: usize,
) -> (f64, usize, usize, Option<Vec<Vec<f64>>>) {
    let mut parameters = Vec::<Parameter>::new();
    let mut lookup = BTreeMap::new();
    for (row, e) in equations.iter().enumerate() {
        if let ConstraintOrigin::AxisAnchor {
            axis,
            node_id,
            component,
        } = e.origin
        {
            let a = &axes[axis];
            let node = node_map[&node_id];
            let length = reference[a.endpoints[0]].distance(reference[a.endpoints[1]]);
            // Short features preserve their complete source vector and parameters.
            if a.endpoints.contains(&node) || length < policy.minimum_length {
                continue;
            }
            let index = *lookup.entry((axis, node)).or_insert_with(|| {
                let anchor = a.anchors.iter().position(|c| c.node == node).unwrap();
                let index = parameters.len();
                parameters.push(Parameter {
                    axis,
                    anchor,
                    variable: x.len(),
                    scale: length,
                    rows: vec![],
                });
                // Length-scaled parameters keep all unknowns in the same units.
                x.push(a.anchors[anchor].t * length);
                index
            });
            parameters[index].rows.push((row, component));
        }
    }
    let mut iterations = 0;
    let mut steps = 0;
    for _ in 0..maximum_steps {
        linearize(equations, x, axes, &parameters);
        let maximum = equations
            .iter()
            .map(|e| e.residual(x).abs())
            .fold(0.0_f64, f64::max);
        if maximum <= policy.residual_tolerance {
            break;
        }
        let cost: f64 = equations.iter().map(|e| e.residual(x).powi(2)).sum();
        let mut proposal = x.clone();
        let (_, count) = project(
            equations,
            &mut proposal,
            policy.iterations,
            policy.residual_tolerance * 0.1,
        );
        iterations += count;
        let mut accepted = false;
        // Re-evaluate the nonlinear equations, not only the tangent system.
        for backtrack in 0..20 {
            let fraction = 0.5_f64.powi(backtrack);
            let trial: Vec<_> = x
                .iter()
                .zip(&proposal)
                .map(|(a, b)| a + fraction * (b - a))
                .collect();
            if trial.iter().any(|v| !v.is_finite()) {
                continue;
            }
            linearize(equations, &trial, axes, &parameters);
            let trial_cost: f64 = equations.iter().map(|e| e.residual(&trial).powi(2)).sum();
            if trial_cost < cost {
                *x = trial;
                steps += 1;
                accepted = true;
                break;
            }
        }
        if !accepted {
            break;
        }
    }
    linearize(equations, x, axes, &parameters);
    let maximum = equations
        .iter()
        .map(|e| e.residual(x).abs())
        .fold(0.0_f64, f64::max);
    let mut ts: Vec<Vec<f64>> = axes
        .iter()
        .map(|a| a.anchors.iter().map(|c| c.t).collect())
        .collect();
    for p in parameters {
        ts[p.axis][p.anchor] = x[p.variable] / p.scale;
    }
    (maximum, iterations, steps, Some(ts))
}
