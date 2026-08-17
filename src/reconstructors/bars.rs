use crate::config::ReconstructionConfig;
use crate::models::{BarType, ElementData, MacroBar, MeshData};
use glam::DVec3;
use hashbrown::{HashMap, HashSet};

pub struct BarReconstructor;

impl BarReconstructor {
    /// Двунаправленная трассировка стержней с делением колонн по отметкам перекрытий
    pub fn reconstruct(
        mesh_data: &MeshData,
        canonical_nodes: &HashMap<u32, u32>,
        slab_elevations: &[f64],
        config: &ReconstructionConfig,
    ) -> Vec<MacroBar> {
        let bar_elems: Vec<&ElementData> = mesh_data
            .elements
            .iter()
            .filter(|e| e.nodes.len() == 2)
            .collect();

        // Граф смежности: node -> [(neighbor_node, elem_id, stiff_id)]
        let mut adj: HashMap<u32, Vec<(u32, u32, u32)>> = HashMap::new();
        for el in &bar_elems {
            let n1 = canonical_nodes.get(&el.nodes[0]).copied().unwrap_or(el.nodes[0]);
            let n2 = canonical_nodes.get(&el.nodes[1]).copied().unwrap_or(el.nodes[1]);
            adj.entry(n1).or_default().push((n2, el.id, el.stiff_id));
            adj.entry(n2).or_default().push((n1, el.id, el.stiff_id));
        }

        let mut visited = HashSet::new();
        let mut macro_bars = Vec::new();

        for el in &bar_elems {
            if visited.contains(&el.id) {
                continue;
            }

            let stiff = el.stiff_id;
            let n1 = canonical_nodes.get(&el.nodes[0]).copied().unwrap_or(el.nodes[0]);
            let n2 = canonical_nodes.get(&el.nodes[1]).copied().unwrap_or(el.nodes[1]);

            let (Some(&p1), Some(&p2)) = (mesh_data.nodes.get(&n1), mesh_data.nodes.get(&n2)) else {
                continue;
            };

            let v = p2 - p1;
            let v_norm = v.length();
            if v_norm < 1e-6 {
                continue;
            }
            let dir_vec = v / v_norm;
            visited.insert(el.id);

            let mut trace_side = |start_node: u32, prev_node: u32, sign: f64| -> Vec<u32> {
                let mut curr = start_node;
                let mut prev = prev_node;
                let mut chain = Vec::new();

                loop {
                    let mut next_found = None;
                    if let Some(neighbors) = adj.get(&curr) {
                        for &(next_node, next_elem_id, next_stiff) in neighbors {
                            if next_node == prev || visited.contains(&next_elem_id) || next_stiff != stiff {
                                continue;
                            }

                            if let (Some(&p_curr), Some(&p_next)) =
                                (mesh_data.nodes.get(&curr), mesh_data.nodes.get(&next_node))
                            {
                                let v_next = p_next - p_curr;
                                let norm_next = v_next.length();
                                if norm_next < 1e-6 {
                                    continue;
                                }
                                let dir_next = v_next / norm_next;

                                if (dir_vec * sign).dot(dir_next) > 0.999 {
                                    next_found = Some((next_node, next_elem_id));
                                    break;
                                }
                            }
                        }
                    }

                    if let Some((next_node, next_elem_id)) = next_found {
                        visited.insert(next_elem_id);
                        chain.push(next_node);
                        prev = curr;
                        curr = next_node;
                    } else {
                        break;
                    }
                }
                chain
            };

            let forward_chain = trace_side(n2, n1, 1.0);
            let backward_chain = trace_side(n1, n2, -1.0);

            let mut full_chain = Vec::with_capacity(backward_chain.len() + 2 + forward_chain.len());
            for &n in backward_chain.iter().rev() {
                full_chain.push(n);
            }
            full_chain.push(n1);
            full_chain.push(n2);
            for &n in &forward_chain {
                full_chain.push(n);
            }

            let (Some(&p_start), Some(&p_end)) = (
                mesh_data.nodes.get(&full_chain[0]),
                mesh_data.nodes.get(&full_chain[full_chain.len() - 1]),
            ) else {
                continue;
            };

            let total_vec = p_end - p_start;
            let total_len = total_vec.length();
            if total_len < 1e-4 {
                continue;
            }

            let unit_z = (total_vec.z / total_len).abs();
            let is_column = unit_z > 0.95;
            let is_beam = unit_z < 0.05;
            let bar_type = if is_column {
                BarType::Column
            } else if is_beam {
                BarType::Beam
            } else {
                BarType::Brace
            };

            // Деление колонн по отметкам плит
            let mut sub_chains = Vec::new();
            if is_column && !slab_elevations.is_empty() {
                let mut current_sub = vec![full_chain[0]];
                for i in 1..full_chain.len() - 1 {
                    let node = full_chain[i];
                    let z_node = mesh_data.nodes.get(&node).map(|p| p.z).unwrap_or(0.0);
                    current_sub.push(node);

                    if slab_elevations.iter().any(|&sz| (z_node - sz).abs() < config.tol_dist) {
                        if current_sub.len() >= 2 {
                            sub_chains.push(current_sub);
                            current_sub = vec![node];
                        }
                    }
                }
                current_sub.push(full_chain[full_chain.len() - 1]);
                if current_sub.len() >= 2 {
                    sub_chains.push(current_sub);
                }
            } else {
                sub_chains.push(full_chain);
            }

            for chain in sub_chains {
                let (Some(&start_p), Some(&end_p)) = (
                    mesh_data.nodes.get(&chain[0]),
                    mesh_data.nodes.get(&chain[chain.len() - 1]),
                ) else {
                    continue;
                };

                let seg_len = (end_p - start_p).length();
                if seg_len < 1e-4 {
                    continue;
                }

                macro_bars.push(MacroBar {
                    bar_type,
                    stiffness_id: stiff,
                    start_point: [round_4(start_p.x), round_4(start_p.y), round_4(start_p.z)],
                    end_point: [round_4(end_p.x), round_4(end_p.y), round_4(end_p.z)],
                    length: (seg_len * 1000.0).round() / 1000.0,
                });
            }
        }

        macro_bars
    }
}

fn round_4(val: f64) -> f64 {
    (val * 10000.0).round() / 10000.0
}