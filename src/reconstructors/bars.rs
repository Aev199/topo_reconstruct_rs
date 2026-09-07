use crate::config::ReconstructionConfig;
use crate::models::{BarConstraint, BarPropertySpan, BarType, MacroBar, MeshData};
use glam::DVec3;
use hashbrown::{HashMap, HashSet};

pub struct BarReconstructor;

pub fn classify_bar(start: DVec3, end: DVec3) -> BarType {
    let d = end - start;
    let z = (d.z / d.length()).abs();
    if z > 0.95 {
        BarType::Column
    } else if z < 0.05 {
        BarType::Beam
    } else {
        BarType::Brace
    }
}

impl BarReconstructor {
    pub fn reconstruct(
        mesh: &MeshData,
        canonical: &HashMap<u32, u32>,
        _levels: &[f64],
        config: &ReconstructionConfig,
    ) -> Vec<MacroBar> {
        let canonical_id = |id: u32| canonical.get(&id).copied().unwrap_or(id);
        let elements: Vec<_> = mesh.elements.iter().filter(|e| e.is_bar()).collect();
        let protected: HashSet<_> = mesh
            .elements
            .iter()
            .filter(|e| e.nodes.len() == 1)
            .map(|e| canonical_id(e.nodes[0]))
            .collect();
        let mut adjacency: HashMap<u32, Vec<(u32, u32, u32)>> = HashMap::new();
        for e in &elements {
            let a = canonical_id(e.nodes[0]);
            let b = canonical_id(e.nodes[1]);
            if a == b {
                continue;
            }
            adjacency.entry(a).or_default().push((b, e.id, e.stiff_id));
            adjacency.entry(b).or_default().push((a, e.id, e.stiff_id));
        }
        let element_by_id: HashMap<_, _> = mesh.elements.iter().map(|e| (e.id, e)).collect();
        let mut visited = HashSet::new();
        let mut result = Vec::new();
        for e in elements {
            if visited.contains(&e.id) {
                continue;
            }
            let a = canonical_id(e.nodes[0]);
            let b = canonical_id(e.nodes[1]);
            let (Some(&pa), Some(&pb)) = (mesh.nodes.get(&a), mesh.nodes.get(&b)) else {
                continue;
            };
            if pa.distance(pb) < 1e-6 {
                continue;
            }
            let direction = (pb - pa).normalize();
            visited.insert(e.id);
            let mut trace = |start: u32, previous: u32, sign: f64| {
                let mut chain = Vec::new();
                let mut current = start;
                let mut previous = previous;
                loop {
                    let Some(neighbors) = adjacency.get(&current) else {
                        break;
                    };
                    // A transverse branch or support is metadata, not a break in the line.
                    // Ambiguous overlapping forward continuations are not selected arbitrarily.
                    let candidates: Vec<_> = neighbors
                        .iter()
                        .filter(|&&(node, id, _)| {
                            node != previous
                                && !visited.contains(&id)
                                && mesh.nodes.get(&node).is_some_and(|p| {
                                    let d = *p - mesh.nodes[&current];
                                    d.length() > 1e-6 && d.normalize().dot(direction * sign) > 0.999
                                })
                        })
                        .copied()
                        .collect();
                    if candidates.len() != 1 {
                        break;
                    }
                    let (node, id, _) = candidates[0];
                    visited.insert(id);
                    chain.push((node, id));
                    previous = current;
                    current = node;
                }
                chain
            };
            let forward = trace(b, a, 1.0);
            let backward = trace(a, b, -1.0);
            let mut nodes: Vec<_> = backward.iter().rev().map(|p| p.0).collect();
            nodes.extend([a, b]);
            nodes.extend(forward.iter().map(|p| p.0));
            let mut ids: Vec<_> = backward.iter().rev().map(|p| p.1).collect();
            ids.push(e.id);
            ids.extend(forward.iter().map(|p| p.1));
            // Split only where the actual polyline deviates from a straight chord.
            // A straight chain spanning storeys, supports or property changes stays one line.
            let mut cuts = vec![0, nodes.len() - 1];
            let mut stack = vec![(0, nodes.len() - 1)];
            while let Some((lo, hi)) = stack.pop() {
                if hi <= lo + 1 {
                    continue;
                }
                let a = mesh.nodes[&nodes[lo]];
                let d = mesh.nodes[&nodes[hi]] - a;
                let mut worst = (0.0, lo);
                for i in lo + 1..hi {
                    let p = mesh.nodes[&nodes[i]];
                    let t = if d.length_squared() > 1e-20 {
                        ((p - a).dot(d) / d.length_squared()).clamp(0.0, 1.0)
                    } else {
                        0.0
                    };
                    let distance = p.distance(a + t * d);
                    if distance > worst.0 {
                        worst = (distance, i);
                    }
                }
                if worst.0 > config.simplify_tol {
                    cuts.push(worst.1);
                    stack.push((lo, worst.1));
                    stack.push((worst.1, hi));
                }
            }
            cuts.sort_unstable();
            cuts.dedup();
            for range in cuts.windows(2) {
                let lo = range[0];
                let hi = range[1];
                let start = mesh.nodes[&nodes[lo]];
                let end = mesh.nodes[&nodes[hi]];
                let d = end - start;
                if d.length() < 1e-6 {
                    continue;
                }
                let t = |index: usize| {
                    ((mesh.nodes[&nodes[index]] - start).dot(d) / d.length_squared())
                        .clamp(0.0, 1.0)
                };
                let mut constraints = Vec::new();
                for i in lo + 1..hi {
                    let mut kinds = Vec::new();
                    if adjacency.get(&nodes[i]).is_some_and(|n| n.len() != 2) {
                        kinds.push("junction".into());
                    }
                    if protected.contains(&nodes[i]) {
                        kinds.push("support".into());
                    }
                    if element_by_id[&ids[i - 1]].stiff_id != element_by_id[&ids[i]].stiff_id {
                        kinds.push("property_boundary".into());
                    }
                    if !kinds.is_empty() {
                        constraints.push(BarConstraint {
                            t: t(i),
                            source_node_id: Some(nodes[i]),
                            kinds,
                            panel_ids: vec![],
                        });
                    }
                }
                let mut property_spans: Vec<BarPropertySpan> = Vec::new();
                for i in lo..hi {
                    let stiffness = element_by_id[&ids[i]].stiff_id;
                    if let Some(span) = property_spans
                        .last_mut()
                        .filter(|span| span.stiffness_id == stiffness)
                    {
                        span.end_t = if i + 1 == hi { 1.0 } else { t(i + 1) };
                        span.source_element_ids.push(ids[i]);
                    } else {
                        property_spans.push(BarPropertySpan {
                            start_t: if i == lo { 0.0 } else { t(i) },
                            end_t: if i + 1 == hi { 1.0 } else { t(i + 1) },
                            stiffness_id: stiffness,
                            source_element_ids: vec![ids[i]],
                        });
                    }
                }
                for span in &mut property_spans {
                    span.source_element_ids.sort_unstable();
                }
                let mut source_element_ids = ids[lo..hi].to_vec();
                source_element_ids.sort_unstable();
                result.push(MacroBar {
                    bar_type: classify_bar(start, end),
                    stiffness_id: if property_spans.len() == 1 {
                        Some(property_spans[0].stiffness_id)
                    } else {
                        None
                    },
                    start_point: start.to_array(),
                    end_point: end.to_array(),
                    length: d.length(),
                    start_node_id: nodes[lo],
                    end_node_id: nodes[hi],
                    source_element_ids,
                    source_node_ids: nodes[lo..=hi].to_vec(),
                    start_panel_ids: vec![],
                    end_panel_ids: vec![],
                    constraints,
                    property_spans,
                });
            }
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::ElementData;
    fn element(id: u32, a: u32, b: u32) -> ElementData {
        ElementData {
            id,
            elem_type: 10,
            stiff_id: 1,
            nodes: vec![a, b],
        }
    }
    #[test]
    fn branch_is_preserved_and_each_source_element_is_accounted_for() {
        let mesh = MeshData {
            nodes: HashMap::from_iter([
                (1, DVec3::ZERO),
                (2, DVec3::X),
                (3, DVec3::X * 2.),
                (4, DVec3::new(1., 1., 0.)),
            ]),
            elements: vec![element(1, 1, 2), element(2, 2, 3), element(3, 2, 4)],
        };
        let bars = BarReconstructor::reconstruct(
            &mesh,
            &HashMap::new(),
            &[],
            &ReconstructionConfig::default(),
        );
        assert_eq!(bars.len(), 2);
        let main = bars
            .iter()
            .find(|b| b.source_element_ids.len() == 2)
            .unwrap();
        assert_eq!(main.constraints.len(), 1);
        assert_eq!(main.constraints[0].source_node_id, Some(2));
        assert!((main.constraints[0].t - 0.5).abs() < 1e-12);
        let mut ids: Vec<_> = bars
            .iter()
            .flat_map(|b| b.source_element_ids.iter().copied())
            .collect();
        ids.sort_unstable();
        assert_eq!(ids, vec![1, 2, 3]);
    }
    #[test]
    fn support_is_a_constraint_and_does_not_split_geometry() {
        let mut mesh = MeshData {
            nodes: HashMap::from_iter([(1, DVec3::ZERO), (2, DVec3::X), (3, DVec3::X * 2.)]),
            elements: vec![element(1, 1, 2), element(2, 2, 3)],
        };
        let config = ReconstructionConfig::default();
        let bars = BarReconstructor::reconstruct(&mesh, &HashMap::new(), &[], &config);
        assert_eq!(bars.len(), 1);
        assert_eq!(bars[0].source_element_ids, vec![1, 2]);
        mesh.elements.push(ElementData {
            id: 3,
            elem_type: 56,
            stiff_id: 1,
            nodes: vec![2],
        });
        assert_eq!(
            BarReconstructor::reconstruct(&mesh, &HashMap::new(), &[], &config).len(),
            1
        );
    }
    #[test]
    fn levels_and_properties_do_not_split_a_straight_column() {
        let mut mesh = MeshData {
            nodes: HashMap::from_iter([(1, DVec3::ZERO), (2, DVec3::Z * 3.), (3, DVec3::Z * 6.)]),
            elements: vec![element(1, 1, 2), element(2, 2, 3)],
        };
        mesh.elements[1].stiff_id = 2;
        let bars = BarReconstructor::reconstruct(
            &mesh,
            &HashMap::new(),
            &[0., 3., 6.],
            &ReconstructionConfig::default(),
        );
        assert_eq!(bars.len(), 1);
        let bar = &bars[0];
        assert_eq!(bar.length, 6.0);
        assert_eq!(bar.stiffness_id, None);
        assert_eq!(bar.property_spans.len(), 2);
        assert_eq!(bar.property_spans[0].end_t, 0.5);
        assert_eq!(bar.property_spans[1].start_t, 0.5);
        assert_eq!(bar.source_element_ids, vec![1, 2]);
    }
}
