use crate::models::ElementData;
use glam::DVec3;
use hashbrown::HashMap;

/// Извлечение 3D-отрезков ребер стен и балок, лежащих на высотной отметке перекрытия Z
pub fn extract_cutting_segments_for_slab(
    elements: &[ElementData],
    nodes: &HashMap<u32, DVec3>,
    canonical_nodes: &HashMap<u32, u32>,
    z_slab: f64,
    tol_dist: f64,
    split_by_walls: bool,
    split_by_beams: bool,
) -> Vec<(DVec3, DVec3)> {
    let mut cutting_segments = Vec::new();

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

            let n_len = pts.len();
            for i in 0..n_len {
                let p_a = pts[i];
                let p_b = pts[(i + 1) % n_len];

                if (p_a.z - z_slab).abs() < tol_dist && (p_b.z - z_slab).abs() < tol_dist {
                    if (p_b - p_a).length_squared() > 1e-6 {
                        cutting_segments.push((p_a, p_b));
                    }
                }
            }
        }
    }

    // 2. Оси балок на отметке плиты
    if split_by_beams {
        for el in elements.iter().filter(|e| e.nodes.len() == 2) {
            let p_a = canonical_nodes
                .get(&el.nodes[0])
                .and_then(|cid| nodes.get(cid).copied());
            let p_b = canonical_nodes
                .get(&el.nodes[1])
                .and_then(|cid| nodes.get(cid).copied());

            if let (Some(pa), Some(pb)) = (p_a, p_b) {
                if (pa.z - z_slab).abs() < tol_dist && (pb.z - z_slab).abs() < tol_dist {
                    if (pb - pa).length_squared() > 1e-6 {
                        cutting_segments.push((pa, pb));
                    }
                }
            }
        }
    }

    cutting_segments
}