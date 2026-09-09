//! Joint regularization of already recognized axes. FE chain recognition is separate.
use glam::DVec3;
use serde::Serialize;

#[derive(Debug, Clone, Serialize)]
pub struct Axis {
    pub ends: [usize; 2],
    pub source_elements: Vec<u32>,
    /// Explicitly protected inclinations or other intentional features.
    pub protected: bool,
}
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub enum Direction {
    Vertical,
    Horizontal,
    Inclined,
}
#[derive(Debug, Clone)]
pub struct Policy {
    pub up: [f64; 3],
    pub angular_tolerance: f64,
    pub relative_movement: f64,
    pub maximum_movement: f64,
    pub minimum_length: f64,
}
#[derive(Debug, Clone, PartialEq)]
pub enum Failure {
    InvalidInput,
    BudgetExceeded,
    ProtectedAxis,
    DegenerateAxis,
}
#[derive(Debug, Clone, Serialize)]
pub struct Result {
    pub points: Vec<[f64; 3]>,
    pub directions: Vec<Direction>,
    pub movements: Vec<f64>,
    pub source_elements: Vec<Vec<u32>>,
}

/// Always solve from immutable reference points, not a previously moved state.
/// Shared endpoint indices are the only evidence of a connection in this API.
/// Return a complete proposal or a failure, never a partially regularized graph.
pub fn regularize(
    reference: &[[f64; 3]],
    axes: &[Axis],
    policy: &Policy,
) -> std::result::Result<Result, Failure> {
    let up = DVec3::from_array(policy.up);
    if !up.is_finite()
        || !up.length().is_finite()
        || up.length() == 0.0
        || !policy.angular_tolerance.is_finite()
        || policy.angular_tolerance <= 0.0
        || policy.angular_tolerance >= std::f64::consts::FRAC_PI_4
        || [
            policy.relative_movement,
            policy.maximum_movement,
            policy.minimum_length,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v <= 0.0)
        || reference.iter().any(|p| !DVec3::from_array(*p).is_finite())
    {
        return Err(Failure::InvalidInput);
    }
    let up = up.normalize();
    let helper = if up.z.abs() < 0.9 { DVec3::Z } else { DVec3::X };
    let u = helper.cross(up).normalize();
    let v = up.cross(u);
    let basis = [u, v, up];
    let origin = reference
        .first()
        .copied()
        .map(DVec3::from_array)
        .unwrap_or(DVec3::ZERO);
    let coordinates: Vec<[f64; 3]> = reference
        .iter()
        .map(|p| {
            let p = DVec3::from_array(*p) - origin;
            basis.map(|axis| axis.dot(p))
        })
        .collect();
    let mut parents: [Vec<usize>; 3] = std::array::from_fn(|_| (0..reference.len()).collect());
    let mut directions = vec![];
    let mut budgets = vec![policy.maximum_movement; reference.len()];
    for axis in axes {
        let [a, b] = axis.ends;
        if a >= reference.len() || b >= reference.len() || a == b {
            return Err(Failure::InvalidInput);
        }
        let d = DVec3::from_array(reference[b]) - DVec3::from_array(reference[a]);
        if !d.length().is_finite() || d.length() < policy.minimum_length {
            return Err(Failure::DegenerateAxis);
        }
        let cosine = d.normalize().dot(up).abs();
        let direction = if cosine >= policy.angular_tolerance.cos() {
            Direction::Vertical
        } else if cosine <= policy.angular_tolerance.sin() {
            Direction::Horizontal
        } else {
            Direction::Inclined
        };
        directions.push(direction);
        let cap = policy
            .maximum_movement
            .min(policy.relative_movement * d.length());
        budgets[a] = budgets[a].min(cap);
        budgets[b] = budgets[b].min(cap);
        if axis.protected {
            continue;
        }
        for (component, parent) in parents.iter_mut().enumerate() {
            if (direction == Direction::Vertical && component < 2)
                || (direction == Direction::Horizontal && component == 2)
            {
                let ra = root(parent, a);
                let rb = root(parent, b);
                parent[rb] = ra;
            }
        }
    }
    let mut adjusted = coordinates.clone();
    for (component, parent) in parents.iter().enumerate() {
        let mut groups = std::collections::BTreeMap::<usize, Vec<f64>>::new();
        for (i, p) in coordinates.iter().enumerate() {
            groups
                .entry(root(parent, i))
                .or_default()
                .push(p[component]);
        }
        let averages: std::collections::BTreeMap<_, _> = groups
            .into_iter()
            .map(|(id, mut values)| {
                // Summation order does not depend on input numbering.
                values.sort_by(f64::total_cmp);
                (id, values.iter().sum::<f64>() / values.len() as f64)
            })
            .collect();
        for (i, p) in adjusted.iter_mut().enumerate() {
            p[component] = averages[&root(parent, i)];
        }
    }
    let points: Vec<_> = adjusted
        .iter()
        .map(|p| (origin + p[0] * u + p[1] * v + p[2] * up).to_array())
        .collect();
    let movements: Vec<_> = points
        .iter()
        .zip(reference)
        .map(|(a, b)| DVec3::from_array(*a).distance(DVec3::from_array(*b)))
        .collect();
    if movements
        .iter()
        .zip(budgets)
        .any(|(d, cap)| !d.is_finite() || *d > cap + 1e-10)
    {
        return Err(Failure::BudgetExceeded);
    }
    for axis in axes {
        let [a, b] = axis.ends;
        if axis.protected && (movements[a] > 1e-10 || movements[b] > 1e-10) {
            return Err(Failure::ProtectedAxis);
        }
        let old = DVec3::from_array(reference[b]) - DVec3::from_array(reference[a]);
        let new = DVec3::from_array(points[b]) - DVec3::from_array(points[a]);
        if new.length() < policy.minimum_length || new.dot(old) <= 0.0 {
            return Err(Failure::DegenerateAxis);
        }
    }
    Ok(Result {
        points,
        directions,
        movements,
        source_elements: axes.iter().map(|a| a.source_elements.clone()).collect(),
    })
}
fn root(parent: &[usize], mut i: usize) -> usize {
    while parent[i] != i {
        i = parent[i];
    }
    i
}

#[cfg(test)]
mod tests {
    use super::*;
    fn policy() -> Policy {
        Policy {
            up: [0., 0., 1.],
            angular_tolerance: 0.02,
            relative_movement: 0.05,
            maximum_movement: 0.15,
            minimum_length: 0.03,
        }
    }
    fn axis(a: usize, b: usize) -> Axis {
        Axis {
            ends: [a, b],
            source_elements: vec![9],
            protected: false,
        }
    }
    #[test]
    fn stacked_columns_share_one_vertical_axis() {
        let source = [[0., 0., 0.], [0.01, 0., 3.], [0.02, 0., 6.]];
        let result = regularize(&source, &[axis(0, 1), axis(1, 2)], &policy()).unwrap();
        assert_eq!(result.points[0][0], result.points[1][0]);
        assert_eq!(result.points[1][0], result.points[2][0]);
        assert_eq!(result.points[2][2], 6.);
        assert_eq!(result.source_elements, vec![vec![9], vec![9]]);
    }
    #[test]
    fn horizontal_beam_keeps_plan_direction() {
        let result = regularize(&[[0., 0., 1.], [3., 4., 1.01]], &[axis(0, 1)], &policy()).unwrap();
        assert_eq!(result.points[0][2], result.points[1][2]);
        assert_eq!(&result.points[1][..2], &[3., 4.]);
    }
    #[test]
    fn small_angle_does_not_override_movement_budget() {
        assert_eq!(
            regularize(&[[0., 0., 0.], [1., 0., 100.]], &[axis(0, 1)], &policy()).unwrap_err(),
            Failure::BudgetExceeded
        );
    }
    #[test]
    fn protected_axis_blocks_incompatible_group_movement() {
        let mut protected = axis(1, 2);
        protected.protected = true;
        assert_eq!(
            regularize(
                &[[0., 0., 0.], [0.01, 0., 3.], [2., 0., 5.]],
                &[axis(0, 1), protected],
                &policy()
            )
            .unwrap_err(),
            Failure::ProtectedAxis
        );
    }
    #[test]
    fn idempotent_and_covariant_with_vertical() {
        let source = [[0., 0., 0.], [0.01, 0., 3.], [3.01, 4., 3.01]];
        let axes = [axis(0, 1), axis(1, 2)];
        let first = regularize(&source, &axes, &policy()).unwrap();
        let second = regularize(&first.points, &axes, &policy()).unwrap();
        assert!(second.movements.iter().all(|d| *d < 1e-10));
        let rotation = glam::DQuat::from_axis_angle(DVec3::new(1., 2., 3.).normalize(), 0.73);
        let offset = DVec3::new(100., -20., 30.);
        let transformed: Vec<_> = source
            .iter()
            .map(|p| (rotation * DVec3::from_array(*p) + offset).to_array())
            .collect();
        let mut config = policy();
        config.up = (rotation * DVec3::Z).to_array();
        let result = regularize(&transformed, &axes, &config).unwrap();
        for (a, b) in result.points.iter().zip(first.points) {
            assert!(
                DVec3::from_array(*a).distance(rotation * DVec3::from_array(b) + offset) < 1e-10
            );
        }
    }
}
