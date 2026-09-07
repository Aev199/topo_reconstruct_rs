use crate::config::ReconstructionConfig;
use crate::models::{BarType, MacroBar, MeshData};
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
        levels: &[f64],
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
                    // Preserve branches, supports and spring nodes as explicit endpoints.
                    if neighbors.len() != 2 || protected.contains(&current) {
                        break;
                    }
                    let next = neighbors
                        .iter()
                        .find(|&&(node, id, stiff)| {
                            node != previous
                                && !visited.contains(&id)
                                && stiff == e.stiff_id
                                && mesh.nodes.get(&node).is_some_and(|p| {
                                    let d = *p - mesh.nodes[&current];
                                    d.length() > 1e-6 && d.normalize().dot(direction * sign) > 0.999
                                })
                        })
                        .copied();
                    let Some((node, id, _)) = next else {
                        break;
                    };
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
            let column = classify_bar(mesh.nodes[&nodes[0]], mesh.nodes[nodes.last().unwrap()])
                == BarType::Column;
            let mut cuts = vec![0];
            for i in 1..nodes.len() - 1 {
                if column
                    && levels
                        .iter()
                        .any(|z| (mesh.nodes[&nodes[i]].z - z).abs() < config.tol_dist)
                {
                    cuts.push(i);
                }
            }
            cuts.push(nodes.len() - 1);
            for range in cuts.windows(2) {
                let lo = range[0];
                let hi = range[1];
                let start = mesh.nodes[&nodes[lo]];
                let end = mesh.nodes[&nodes[hi]];
                let d = end - start;
                let straight = d.length() > 1e-6
                    && nodes[lo..=hi].iter().all(|id| {
                        let p = mesh.nodes[id];
                        let t = ((p - start).dot(d) / d.length_squared()).clamp(0.0, 1.0);
                        p.distance(start + t * d) <= config.simplify_tol
                    });
                let ranges: Vec<_> = if straight {
                    vec![(lo, hi)]
                } else {
                    (lo..hi).map(|i| (i, i + 1)).collect()
                };
                for (lo, hi) in ranges {
                    let start = mesh.nodes[&nodes[lo]];
                    let end = mesh.nodes[&nodes[hi]];
                    let mut source_element_ids = ids[lo..hi].to_vec();
                    source_element_ids.sort_unstable();
                    result.push(MacroBar {
                        bar_type: classify_bar(start, end),
                        stiffness_id: e.stiff_id,
                        start_point: start.to_array(),
                        end_point: end.to_array(),
                        length: start.distance(end),
                        start_node_id: nodes[lo],
                        end_node_id: nodes[hi],
                        source_element_ids,
                        source_node_ids: nodes[lo..=hi].to_vec(),
                        start_panel_ids: vec![],
                        end_panel_ids: vec![],
                    });
                }
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
        assert_eq!(bars.len(), 3);
        assert!(bars
            .iter()
            .all(|b| b.start_node_id == 2 || b.end_node_id == 2));
        let mut ids: Vec<_> = bars
            .iter()
            .flat_map(|b| b.source_element_ids.iter().copied())
            .collect();
        ids.sort_unstable();
        assert_eq!(ids, vec![1, 2, 3]);
    }
    #[test]
    fn support_stops_merging_but_plain_subdivision_does_not() {
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
            2
        );
    }
}
