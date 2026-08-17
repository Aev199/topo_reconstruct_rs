#![allow(dead_code, unused_imports, unused_variables)]

use crate::models::ElementData;
use glam::DVec3;
use hashbrown::{HashMap, HashSet};

/// Извлечение пар канонических узлов (ребер) стен и балок на отметке перекрытия Z
pub fn extract_cutting_edge_nodes_for_slab(
    elements: &[ElementData],
    nodes: &HashMap<u32, DVec3>,
    canonical_nodes: &HashMap<u32, u32>,
    z_slab: f64,
    tol_dist: f64,
    split_by_walls: bool,
    split_by_beams: bool,
) -> HashSet<(u32, u32)> {
    let mut cutting_edges = HashSet::new();

    // 1. Ребра стен на отметке плиты
    if split_by_walls {
        for el in elements.iter().filter(|e| e.nodes.len() >= 3) {
            let pts: Vec<DVec3> = el
                .nodes
                .iter()
                .filter_map(|nid| canonical_nodes.get(nid).and_then(|cid| nodes.get(cid).copied()))
                .collect();

            if pts.len() < 3 {
                continue;
            }

            let v1 = pts[1] - pts[0];
            let v2 = pts[2] - pts[0];
            let norm = v1.cross(v2);
            let norm_len = norm.length();

            // Проверяем, что стена вертикальная (|nz| < 0.15)
            if norm_len < 1e-7 || (norm.z / norm_len).abs() > 0.15 {
                continue;
            }

            let n_len = el.nodes.len();
            for i in 0..n_len {
                let n_a = canonical_nodes.get(&el.nodes[i]).copied().unwrap_or(el.nodes[i]);
                let n_b = canonical_nodes.get(&el.nodes[(i + 1) % n_len]).copied().unwrap_or(el.nodes[(i + 1) % n_len]);

                if let (Some(pa), Some(pb)) = (nodes.get(&n_a), nodes.get(&n_b)) {
                    if (pa.z - z_slab).abs() < tol_dist && (pb.z - z_slab).abs() < tol_dist {
                        let edge = if n_a < n_b { (n_a, n_b) } else { (n_b, n_a) };
                        cutting_edges.insert(edge);
                    }
                }
            }
        }
    }

    // 2. Ребра балок на отметке плиты
    if split_by_beams {
        for el in elements.iter().filter(|e| e.nodes.len() == 2) {
            let n_a = canonical_nodes.get(&el.nodes[0]).copied().unwrap_or(el.nodes[0]);
            let n_b = canonical_nodes.get(&el.nodes[1]).copied().unwrap_or(el.nodes[1]);

            if let (Some(pa), Some(pb)) = (nodes.get(&n_a), nodes.get(&n_b)) {
                if (pa.z - z_slab).abs() < tol_dist && (pb.z - z_slab).abs() < tol_dist {
                    let edge = if n_a < n_b { (n_a, n_b) } else { (n_b, n_a) };
                    cutting_edges.insert(edge);
                }
            }
        }
    }

    cutting_edges
}