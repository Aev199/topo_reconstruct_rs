//! Recover continuous straight axes from FE adjacency without splitting at properties.
use crate::input::{ElementData, MeshData};
#[cfg(test)]
use glam::DVec3;
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Serialize)]
pub struct Policy {
    pub angle: f64,
    pub line_tolerance: f64,
    pub numerical_precision: f64,
}
#[derive(Debug, Clone, Serialize)]
pub struct SourceSpan {
    pub element: u32,
    pub stiffness: u32,
    pub start_t: f64,
    pub end_t: f64,
}
#[derive(Debug, Clone, Serialize)]
pub struct SourceAnchor {
    pub node: u32,
    pub t: f64,
    pub distance_to_axis: f64,
    pub incident_elements: Vec<u32>,
}
#[derive(Debug, Clone, Serialize)]
pub struct RecognizedAxis {
    pub endpoints: [[f64; 3]; 2],
    pub endpoint_nodes: [u32; 2],
    pub spans: Vec<SourceSpan>,
    /// Interior FE nodes remain parametric, including branches and property boundaries.
    pub anchors: Vec<SourceAnchor>,
}
#[derive(Debug, Clone, Serialize)]
pub struct RejectedElement {
    pub element: u32,
    pub reason: String,
}
#[derive(Debug, Clone, Serialize)]
pub struct Report {
    pub policy: Policy,
    pub axes: Vec<RecognizedAxis>,
    pub rejected: Vec<RejectedElement>,
    pub unmerged_components: Vec<Vec<u32>>,
    pub non_bar_elements: Vec<u32>,
}

/// Connected near-collinear FE paths are accepted only when the whole path fits
/// one chord. A failed whole-path fit is partitioned into validated subpaths.
/// Ambiguous optimal partitions retain all possible cut locations.
/// No proximity welding and no interpretation of mechanical connectivity.
pub fn recognize(mesh: &MeshData, policy: &Policy) -> Result<Report, &'static str> {
    if !policy.angle.is_finite()
        || policy.angle <= 0.0
        || policy.angle >= std::f64::consts::FRAC_PI_4
        || !policy.line_tolerance.is_finite()
        || !policy.numerical_precision.is_finite()
        || policy.numerical_precision <= 0.0
        || policy.line_tolerance < policy.numerical_precision
    {
        return Err("invalid recognition policy");
    }
    let mut report = Report {
        policy: policy.clone(),
        axes: vec![],
        rejected: vec![],
        unmerged_components: vec![],
        non_bar_elements: vec![],
    };
    let mut elements = BTreeMap::<u32, &ElementData>::new();
    let mut adjacency = BTreeMap::<u32, Vec<u32>>::new();
    let mut incidence = BTreeMap::<u32, Vec<u32>>::new();
    let mut ids = BTreeSet::new();
    for e in &mesh.elements {
        if !ids.insert(e.id) {
            return Err("duplicate source element ID");
        }
        for &n in &e.nodes {
            incidence.entry(n).or_default().push(e.id);
        }
        if !e.is_bar() {
            report.non_bar_elements.push(e.id);
            continue;
        }
        let valid = e
            .nodes
            .iter()
            .all(|n| mesh.nodes.get(n).is_some_and(|p| p.is_finite()));
        let reason = if !valid {
            Some("missing_or_nonfinite_node")
        } else {
            let length = mesh.nodes[&e.nodes[0]].distance(mesh.nodes[&e.nodes[1]]);
            (!length.is_finite() || length <= policy.numerical_precision)
                .then_some("degenerate_element")
        };
        if let Some(reason) = reason {
            report.rejected.push(RejectedElement {
                element: e.id,
                reason: reason.into(),
            });
            continue;
        }
        elements.insert(e.id, e);
        for &n in &e.nodes {
            adjacency.entry(n).or_default().push(e.id);
        }
    }
    for ids in incidence.values_mut() {
        ids.sort_unstable();
        ids.dedup();
    }
    let other = |e: &ElementData, n| {
        if e.nodes[0] == n {
            e.nodes[1]
        } else {
            e.nodes[0]
        }
    };
    let mut continuation = BTreeMap::<(u32, u32), u32>::new();
    for (&node, attached) in &adjacency {
        for &id in attached {
            let d = (mesh.nodes[&other(elements[&id], node)] - mesh.nodes[&node]).normalize();
            let candidates: Vec<_> = attached
                .iter()
                .copied()
                .filter(|&j| {
                    j != id
                        && d.dot(
                            (mesh.nodes[&other(elements[&j], node)] - mesh.nodes[&node])
                                .normalize(),
                        ) <= -policy.angle.cos()
                })
                .collect();
            if candidates.len() == 1 {
                continuation.insert((id, node), candidates[0]);
            }
        }
    }
    let mutual: BTreeMap<_, _> = continuation
        .iter()
        .filter(|((id, node), next)| continuation.get(&(**next, *node)) == Some(id))
        .map(|(k, v)| (*k, *v))
        .collect();
    let mut visited = BTreeSet::<u32>::new();
    for &seed in elements.keys() {
        if visited.contains(&seed) {
            continue;
        }
        let mut component = BTreeSet::new();
        let mut stack = vec![seed];
        while let Some(id) = stack.pop() {
            if !component.insert(id) {
                continue;
            }
            for &n in &elements[&id].nodes {
                if let Some(&next) = mutual.get(&(id, n)) {
                    stack.push(next);
                }
            }
        }
        visited.extend(component.iter().copied());
        let mut ends: Vec<_> = component
            .iter()
            .flat_map(|&id| elements[&id].nodes.iter().copied().map(move |n| (id, n)))
            .filter(|&(id, n)| !mutual.contains_key(&(id, n)))
            .collect();
        ends.sort_by(|a, b| {
            let p = mesh.nodes[&a.1];
            let q = mesh.nodes[&b.1];
            p.x.total_cmp(&q.x)
                .then(p.y.total_cmp(&q.y))
                .then(p.z.total_cmp(&q.z))
        });
        let mut path = vec![];
        let mut nodes = vec![];
        if ends.len() == 2 {
            let (mut id, mut n) = ends[0];
            nodes.push(n);
            loop {
                path.push(id);
                n = other(elements[&id], n);
                nodes.push(n);
                if let Some(&next) = mutual.get(&(id, n)) {
                    id = next;
                } else {
                    break;
                }
                if path.len() > component.len() {
                    break;
                }
            }
        }
        let build = |path: &[u32], nodes: &[u32]| -> Option<RecognizedAxis> {
            let a = mesh.nodes[nodes.first()?];
            let b = mesh.nodes[nodes.last()?];
            let delta = b - a;
            let length = delta.length();
            if length <= policy.numerical_precision {
                return None;
            }
            let parameters: Vec<_> = nodes
                .iter()
                .map(|n| (mesh.nodes[n] - a).dot(delta) / delta.length_squared())
                .collect();
            let anchors: Vec<_> = nodes
                .iter()
                .zip(&parameters)
                .map(|(&node, &t)| SourceAnchor {
                    node,
                    t,
                    distance_to_axis: mesh.nodes[&node].distance(a + t * delta),
                    incident_elements: incidence[&node].clone(),
                })
                .collect();
            if anchors
                .iter()
                .any(|c| c.distance_to_axis > policy.line_tolerance)
                || parameters.windows(2).any(|p| p[1] <= p[0])
            {
                return None;
            }
            if nodes.windows(2).any(|ns| {
                (mesh.nodes[&ns[1]] - mesh.nodes[&ns[0]])
                    .normalize()
                    .dot(delta / length)
                    < policy.angle.cos()
            }) {
                return None;
            }
            let spans = path
                .iter()
                .enumerate()
                .map(|(i, &id)| SourceSpan {
                    element: id,
                    stiffness: elements[&id].stiff_id,
                    start_t: parameters[i],
                    end_t: parameters[i + 1],
                })
                .collect();
            Some(RecognizedAxis {
                endpoints: [a.to_array(), b.to_array()],
                endpoint_nodes: [nodes[0], *nodes.last()?],
                spans,
                anchors,
            })
        };
        if path.len() == component.len() {
            if let Some(axis) = build(&path, &nodes) {
                report.axes.push(axis);
                continue;
            }
        }
        report
            .unmerged_components
            .push(component.iter().copied().collect());
        if path.len() == component.len() && nodes.len() == path.len() + 1 {
            // Minimum-count valid interval cover. Keep every cut belonging to an
            // optimal cover rather than choosing by source IDs or traversal order.
            let n = path.len();
            let mut valid = vec![vec![false; n + 1]; n];
            for i in 0..n {
                for j in i + 1..=n {
                    valid[i][j] = build(&path[i..j], &nodes[i..=j]).is_some();
                }
            }
            let mut prefix = vec![n + 1; n + 1];
            let mut suffix = vec![n + 1; n + 1];
            prefix[0] = 0;
            suffix[n] = 0;
            for j in 1..=n {
                for i in 0..j {
                    if valid[i][j] {
                        prefix[j] = prefix[j].min(prefix[i] + 1);
                    }
                }
            }
            for i in (0..n).rev() {
                for j in i + 1..=n {
                    if valid[i][j] {
                        suffix[i] = suffix[i].min(suffix[j] + 1);
                    }
                }
            }
            let cuts: Vec<_> = (0..=n)
                .filter(|&i| prefix[i] + suffix[i] == prefix[n])
                .collect();
            for pair in cuts.windows(2) {
                let (i, j) = (pair[0], pair[1]);
                if let Some(axis) = build(&path[i..j], &nodes[i..=j]) {
                    report.axes.push(axis);
                } else {
                    // A subchord need not satisfy its parent's angular bound.
                    for k in i..j {
                        report.axes.push(
                            build(&path[k..k + 1], &nodes[k..=k + 1])
                                .ok_or("failed to retain valid source bar")?,
                        );
                    }
                }
            }
            continue;
        }
        for id in component {
            report.axes.push(
                build(&[id], &elements[&id].nodes).ok_or("failed to retain valid source bar")?,
            );
        }
    }
    report.non_bar_elements.sort_unstable();
    report.rejected.sort_by_key(|e| e.element);
    Ok(report)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn mesh(points: &[[f64; 3]], ends: &[[u32; 2]]) -> MeshData {
        MeshData {
            nodes: points
                .iter()
                .enumerate()
                .map(|(i, p)| (i as u32 + 1, DVec3::from_array(*p)))
                .collect(),
            elements: ends
                .iter()
                .enumerate()
                .map(|(i, n)| ElementData {
                    id: i as u32 + 1,
                    elem_type: 10,
                    stiff_id: i as u32 + 1,
                    nodes: n.to_vec(),
                })
                .collect(),
        }
    }
    fn policy() -> Policy {
        Policy {
            angle: 0.02,
            line_tolerance: 0.01,
            numerical_precision: 1e-8,
        }
    }
    #[test]
    fn branch_and_properties_do_not_split_straight_axis() {
        let m = mesh(
            &[[0., 0., 0.], [1., 0., 0.], [2., 0., 0.], [1., 1., 0.]],
            &[[1, 2], [2, 3], [2, 4]],
        );
        let r = recognize(&m, &policy()).unwrap();
        assert_eq!(r.axes.len(), 2);
        let a = r.axes.iter().find(|a| a.spans.len() == 2).unwrap();
        assert_eq!(a.endpoint_nodes, [1, 3]);
        assert_eq!(a.spans[0].stiffness, 1);
        assert_eq!(a.spans[1].stiffness, 2);
        assert_eq!(a.anchors[1].incident_elements, vec![1, 2, 3]);
        assert_eq!(a.anchors[1].t, 0.5);
    }
    #[test]
    fn ambiguous_branch_is_not_arbitrarily_chosen() {
        let m = mesh(
            &[[0., 0., 0.], [1., 0., 0.], [2., 0., 0.], [2., 0.001, 0.]],
            &[[1, 2], [2, 3], [2, 4]],
        );
        assert_eq!(recognize(&m, &policy()).unwrap().axes.len(), 3);
    }
    #[test]
    fn cumulative_curvature_is_not_hidden_by_small_local_angles() {
        let points: Vec<_> = (0..30)
            .map(|i| {
                let t = i as f64 * 0.01;
                [t.sin() * 10., (1. - t.cos()) * 10., 0.]
            })
            .collect();
        let ends: Vec<_> = (1..30).map(|i| [i, i + 1]).collect();
        let r = recognize(&mesh(&points, &ends), &policy()).unwrap();
        assert!(r.axes.len() > 1);
        let spans: BTreeSet<_> = r
            .axes
            .iter()
            .flat_map(|a| a.spans.iter().map(|s| s.element))
            .collect();
        assert_eq!(spans.len(), 29);
        assert_eq!(r.axes.iter().map(|a| a.spans.len()).sum::<usize>(), 29);
        assert!(r.axes.iter().all(|a| a
            .anchors
            .iter()
            .all(|n| n.distance_to_axis <= policy().line_tolerance)));
        assert_eq!(r.unmerged_components.len(), 1);
    }
    #[test]
    fn failed_whole_path_retains_straight_subpaths_and_short_property_span() {
        let points = [
            [0., 0., 0.],
            [1., 0., 0.],
            [1.006, 0., 0.],
            [2., 0., 0.],
            [3., 0.015, 0.],
            [4., 0.045, 0.],
        ];
        let ends = [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6]];
        let original = mesh(&points, &ends);
        let baseline = recognize(&original, &policy()).unwrap();
        assert_eq!(baseline.unmerged_components.len(), 1);
        let merged = baseline
            .axes
            .iter()
            .find(|a| a.spans.iter().any(|s| s.element == 2))
            .unwrap();
        assert!(merged.spans.len() >= 3);
        assert_eq!(
            merged
                .spans
                .iter()
                .find(|s| s.element == 2)
                .unwrap()
                .stiffness,
            2
        );
        let groups = |r: &Report, offset: u32| -> BTreeSet<BTreeSet<u32>> {
            r.axes
                .iter()
                .map(|a| a.spans.iter().map(|s| s.element - offset).collect())
                .collect()
        };
        for angle in [0.6, 3.0] {
            let rotation = glam::DQuat::from_axis_angle(DVec3::new(1., 2., 3.).normalize(), angle);
            let mut m = original.clone();
            m.nodes = m
                .nodes
                .into_iter()
                .map(|(id, p)| (id + 100, rotation * p + DVec3::new(20., 30., 40.)))
                .collect();
            m.elements.reverse();
            for e in &mut m.elements {
                e.id += 200;
                for n in &mut e.nodes {
                    *n += 100;
                }
                e.nodes.reverse();
            }
            let r = recognize(&m, &policy()).unwrap();
            assert_eq!(groups(&baseline, 0), groups(&r, 200));
            assert_eq!(r.axes.iter().map(|a| a.spans.len()).sum::<usize>(), 5);
        }
    }

    #[test]
    fn invalid_source_is_accounted_for() {
        let m = mesh(&[[0., 0., 0.], [1., 0., 0.]], &[[1, 2], [1, 1], [2, 9]]);
        let r = recognize(&m, &policy()).unwrap();
        assert_eq!(r.axes.len(), 1);
        assert_eq!(r.rejected.len(), 2);
    }
    #[test]
    fn input_order_does_not_change_recognition() {
        let mut m = mesh(
            &[[0., 0., 0.], [1., 0., 0.], [2., 0., 0.]],
            &[[1, 2], [2, 3]],
        );
        let a = serde_json::to_string(&recognize(&m, &policy()).unwrap()).unwrap();
        m.elements.reverse();
        assert_eq!(
            a,
            serde_json::to_string(&recognize(&m, &policy()).unwrap()).unwrap()
        );
    }
}
