#![allow(dead_code, unused_imports, unused_variables, unused_assignments)]

use crate::config::ReconstructionConfig;
use crate::geometry::partition::extract_cutting_edge_nodes_for_slab;
use crate::geometry::utils::{
    clean_polygon_coords_3d, get_plane_basis, remove_short_edges_3d, snap_loop_to_wall_lines,
};
use crate::models::{ElementData, MacroPanel, MeshData, PanelType};
use glam::DVec3;
use hashbrown::{HashMap, HashSet};
use rayon::prelude::*;
use std::collections::VecDeque;

pub struct PanelReconstructor;

#[derive(Clone)]
struct ParsedShell {
    id: u32,
    stiff_id: u32,
    normal: DVec3,
    d: f64,
    centroid: DVec3,
    unique_nodes: Vec<u32>,
}

struct PlaneCluster {
    is_horiz: bool,
    is_vert: bool,
    z: Option<f64>,
    normal: DVec3,
    d: f64,
    elements: Vec<ParsedShell>,
}

impl PanelReconstructor {
    /// Высокоскоростная параллельная реконструкция плит и стен с дотягиванием до стен
    pub fn reconstruct(
        mesh_data: &MeshData,
        canonical_nodes: &HashMap<u32, u32>,
        config: &ReconstructionConfig,
    ) -> Vec<MacroPanel> {
        let shell_elems: Vec<&ElementData> = mesh_data
            .elements
            .iter()
            .filter(|e| e.nodes.len() >= 3)
            .collect();

        // 1. Извлечение опорных линий вертикальных стен для дотягивания плит
        let wall_segments_by_z = Self::extract_wall_datum_lines(mesh_data, canonical_nodes, config.tol_dist);

        // 2. Первичный парсинг элементов и определение их ориентации
        let parsed_shells: Vec<ParsedShell> = shell_elems
            .into_iter()
            .filter_map(|el| {
                let mut unique_nodes = Vec::with_capacity(el.nodes.len());
                for &nid in &el.nodes {
                    let cid = canonical_nodes.get(&nid).copied().unwrap_or(nid);
                    if !unique_nodes.contains(&cid) {
                        unique_nodes.push(cid);
                    }
                }

                if unique_nodes.len() < 3 {
                    return None;
                }

                let pts: Vec<DVec3> = unique_nodes
                    .iter()
                    .filter_map(|cid| mesh_data.nodes.get(cid).copied())
                    .collect();

                if pts.len() < 3 {
                    return None;
                }

                let centroid = pts.iter().copied().sum::<DVec3>() / (pts.len() as f64);
                let v1 = pts[1] - pts[0];
                let v2 = pts[2] - pts[0];
                let mut norm = v1.cross(v2);
                let mut n_len = norm.length();

                if n_len < 1e-7 && pts.len() >= 4 {
                    let v2_alt = pts[3] - pts[0];
                    norm = v1.cross(v2_alt);
                    n_len = norm.length();
                }

                if n_len < 1e-7 {
                    return None;
                }

                norm /= n_len;

                // Для плит нормаль строго +Z
                if norm.z.abs() > 0.85 {
                    norm = DVec3::new(0.0, 0.0, 1.0);
                } else if norm.z < -0.85 {
                    norm = DVec3::new(0.0, 0.0, 1.0);
                }

                let d = -norm.dot(centroid);

                Some(ParsedShell {
                    id: el.id,
                    stiff_id: el.stiff_id,
                    normal: norm,
                    d,
                    centroid,
                    unique_nodes,
                })
            })
            .collect();

        // 3. Кластеризация элементов по плоскостям
        let mut clusters: Vec<PlaneCluster> = Vec::new();
        for el in parsed_shells {
            let nz = el.normal.z;
            let is_horiz = nz.abs() > 0.85;
            let is_vert = nz.abs() < 0.15;

            let mut matched_idx = None;
            for (idx, cl) in clusters.iter().enumerate() {
                if is_horiz && cl.is_horiz {
                    if let Some(cl_z) = cl.z {
                        if (el.centroid.z - cl_z).abs() < config.tol_dist {
                            matched_idx = Some(idx);
                            break;
                        }
                    }
                } else if is_vert && cl.is_vert {
                    let n2d_el = el.normal.truncate().normalize();
                    let n2d_cl = cl.normal.truncate().normalize();
                    if n2d_el.dot(n2d_cl) > 0.98 && (el.d - cl.d).abs() < config.tol_dist {
                        matched_idx = Some(idx);
                        break;
                    }
                } else if !is_horiz && !is_vert && !cl.is_horiz && !cl.is_vert {
                    if el.normal.dot(cl.normal) > 0.98 && (el.d - cl.d).abs() < config.tol_dist {
                        matched_idx = Some(idx);
                        break;
                    }
                }
            }

            if let Some(idx) = matched_idx {
                clusters[idx].elements.push(el);
            } else {
                let cluster_norm = if is_horiz {
                    DVec3::new(0.0, 0.0, 1.0)
                } else if is_vert {
                    let n2d = el.normal.truncate().normalize();
                    DVec3::new(n2d.x, n2d.y, 0.0)
                } else {
                    el.normal
                };

                clusters.push(PlaneCluster {
                    is_horiz,
                    is_vert,
                    z: if is_horiz { Some(el.centroid.z) } else { None },
                    normal: cluster_norm,
                    d: el.d,
                    elements: vec![el],
                });
            }
        }

        // 4. Параллельная реконструкция кластеров (Rayon)
        let panels_results: Vec<Vec<MacroPanel>> = clusters
            .into_par_iter()
            .enumerate()
            .map(|(cl_idx, cl)| {
                let cluster_norm = cl.normal;
                let (u_axis, v_axis) = get_plane_basis(cluster_norm);

                let cutting_edges = if cl.is_horiz {
                    extract_cutting_edge_nodes_for_slab(
                        &mesh_data.elements,
                        &mesh_data.nodes,
                        canonical_nodes,
                        cl.z.unwrap_or(0.0),
                        config.tol_dist,
                        config.split_slabs_by_walls,
                        config.split_slabs_by_beams,
                    )
                } else {
                    HashSet::new()
                };

                let relevant_wall_lines = if cl.is_horiz {
                    let slab_z = cl.z.unwrap_or(0.0);
                    wall_segments_by_z
                        .iter()
                        .filter(|(z, _)| (z - slab_z).abs() < config.tol_dist)
                        .map(|(_, seg)| *seg)
                        .collect::<Vec<_>>()
                } else {
                    Vec::new()
                };

                // Упорядочивание узлов КЭ в едином базисе кластера
                let ordered_elements: Vec<(u32, u32, Vec<u32>)> = cl
                    .elements
                    .iter()
                    .map(|el| {
                        let pts: Vec<DVec3> = el
                            .unique_nodes
                            .iter()
                            .filter_map(|cid| mesh_data.nodes.get(cid).copied())
                            .collect();

                        let elem_centroid = el.centroid;

                        let mut angle_nodes: Vec<(f64, u32)> = pts
                            .iter()
                            .zip(&el.unique_nodes)
                            .map(|(&p, &nid)| {
                                let rel = p - elem_centroid;
                                let u = rel.dot(u_axis);
                                let v = rel.dot(v_axis);
                                (v.atan2(u), nid)
                            })
                            .collect();

                        angle_nodes.sort_by(|a, b| a.0.total_cmp(&b.0));
                        let ordered_nodes: Vec<u32> = angle_nodes.into_iter().map(|(_, nid)| nid).collect();

                        (el.id, el.stiff_id, ordered_nodes)
                    })
                    .collect();

                let elem_by_id: HashMap<u32, &(u32, u32, Vec<u32>)> =
                    ordered_elements.iter().map(|e| (e.0, e)).collect();

                // Построение графа смежности КЭ
                let mut edge_to_elems: HashMap<(u32, u32), Vec<u32>> = HashMap::new();
                for &(elem_id, _, ref nodes) in &ordered_elements {
                    let n_len = nodes.len();
                    for i in 0..n_len {
                        let a = nodes[i];
                        let b = nodes[(i + 1) % n_len];
                        let edge = if a < b { (a, b) } else { (b, a) };
                        edge_to_elems.entry(edge).or_default().push(elem_id);
                    }
                }

                let mut elem_adj: HashMap<u32, Vec<u32>> = HashMap::new();
                for (&edge, e_ids) in &edge_to_elems {
                    if cutting_edges.contains(&edge) {
                        continue;
                    }

                    if e_ids.len() > 1 {
                        for i in 0..e_ids.len() {
                            for j in (i + 1)..e_ids.len() {
                                elem_adj.entry(e_ids[i]).or_default().push(e_ids[j]);
                                elem_adj.entry(e_ids[j]).or_default().push(e_ids[i]);
                            }
                        }
                    }
                }

                // Поиск компонент связности (Connected Components)
                let mut visited_elems: HashSet<u32> = HashSet::new();
                let mut local_panels = Vec::new();
                let mut sub_id = 1;

                for &(start_id, _, _) in &ordered_elements {
                    if visited_elems.contains(&start_id) {
                        continue;
                    }

                    let mut comp_elems = Vec::new();
                    let mut queue = VecDeque::new();

                    queue.push_back(start_id);
                    visited_elems.insert(start_id);

                    while let Some(curr_id) = queue.pop_front() {
                        if let Some(&el_tuple) = elem_by_id.get(&curr_id) {
                            comp_elems.push(el_tuple);
                            if let Some(neighbors) = elem_adj.get(&curr_id) {
                                for &nbr in neighbors {
                                    if visited_elems.insert(nbr) {
                                        queue.push_back(nbr);
                                    }
                                }
                            }
                        }
                    }

                    // 5. Трассировка контуров с дотягиванием до стен
                    let polygons = Self::extract_cycles_from_elements(
                        &comp_elems,
                        mesh_data,
                        u_axis,
                        v_axis,
                        &relevant_wall_lines,
                        config.tol_dist,
                    );

                    if polygons.is_empty() {
                        continue;
                    }

                    let panel_type = if cl.is_horiz {
                        PanelType::Slab
                    } else if cl.is_vert {
                        PanelType::Wall
                    } else {
                        PanelType::InclinedPanel
                    };

                    local_panels.push(MacroPanel {
                        id: (cl_idx * 1000 + sub_id) as u32,
                        panel_type,
                        stiffness_id: comp_elems[0].1,
                        plane_normal: [cluster_norm.x, cluster_norm.y, cluster_norm.z],
                        plane_d: cl.d,
                        polygons,
                        fe_count: comp_elems.len(),
                        connected_panel_ids: vec![],
                    });
                    sub_id += 1;
                }

                local_panels
            })
            .collect();

        panels_results.into_iter().flatten().collect()
    }

    /// Трассировка контуров по Half-Edges с дотягиванием до линий стен
    fn extract_cycles_from_elements(
        elements: &[&(u32, u32, Vec<u32>)],
        mesh_data: &MeshData,
        u_axis: DVec3,
        v_axis: DVec3,
        wall_segments: &[(DVec3, DVec3)],
        snap_tol: f64,
    ) -> Vec<Vec<[f64; 3]>> {
        let mut dir_edge_count: HashMap<(u32, u32), u32> = HashMap::new();
        for (_, _, nodes) in elements {
            let n_len = nodes.len();
            for i in 0..n_len {
                let u = nodes[i];
                let v = nodes[(i + 1) % n_len];
                *dir_edge_count.entry((u, v)).or_insert(0) += 1;
            }
        }

        let mut out_edges: HashMap<u32, Vec<u32>> = HashMap::new();
        let mut boundary_edges = Vec::new();
        for (&(u, v), &count) in &dir_edge_count {
            if count > 0 && !dir_edge_count.contains_key(&(v, u)) {
                boundary_edges.push((u, v));
                out_edges.entry(u).or_default().push(v);
            }
        }

        if boundary_edges.is_empty() {
            return Vec::new();
        }

        let mut visited_half_edges: HashSet<(u32, u32)> = HashSet::new();
        let mut result_polygons = Vec::new();

        for &(start_u, start_v) in &boundary_edges {
            if visited_half_edges.contains(&(start_u, start_v)) {
                continue;
            }

            let mut cycle = vec![start_u];
            let mut curr = start_v;
            visited_half_edges.insert((start_u, start_v));

            let mut closed = false;
            while let Some(neighbors) = out_edges.get(&curr) {
                cycle.push(curr);
                if curr == start_u {
                    closed = true;
                    break;
                }

                let mut next_node = None;
                for &nxt in neighbors {
                    if !visited_half_edges.contains(&(curr, nxt)) {
                        next_node = Some(nxt);
                        break;
                    }
                }

                if let Some(nxt) = next_node {
                    visited_half_edges.insert((curr, nxt));
                    curr = nxt;
                } else {
                    break;
                }
            }

            if closed && cycle.len() >= 4 {
                let loop_nodes = &cycle[..cycle.len() - 1];
                let raw_3d: Vec<DVec3> = loop_nodes
                    .iter()
                    .filter_map(|nid| mesh_data.nodes.get(nid).copied())
                    .collect();

                // 1. Дотягивание точек контура плиты до опорных линий стен
                let snapped_3d = snap_loop_to_wall_lines(&raw_3d, wall_segments, snap_tol);

                // 2. Схлопывание паразитных микросегментов (< 3 см)
                let fused_3d = remove_short_edges_3d(&snapped_3d, 0.03);

                // 3. Вычисление площади для фильтрации вырожденных микропетель
                let mut area_2d = 0.0;
                let n_pts = fused_3d.len();
                for i in 0..n_pts {
                    let p1 = fused_3d[i];
                    let p2 = fused_3d[(i + 1) % n_pts];
                    let x1 = p1.dot(u_axis);
                    let y1 = p1.dot(v_axis);
                    let x2 = p2.dot(u_axis);
                    let y2 = p2.dot(v_axis);
                    area_2d += (x1 * y2) - (x2 * y1);
                }
                area_2d = (area_2d * 0.5).abs();

                if area_2d > 0.10 {
                    // 4. Окончательное выпрямление прямых линий (удаление всех промежуточных узлов)
                    let clean_3d = clean_polygon_coords_3d(&fused_3d, 0.9995);
                    if clean_3d.len() >= 3 {
                        result_polygons.push(clean_3d.into_iter().map(|p| [p.x, p.y, p.z]).collect());
                    }
                }
            }
        }

        result_polygons
    }

    /// Сбор опорных линий (отрезков) всех вертикальных стен с привязкой к отметкам Z
    fn extract_wall_datum_lines(
        mesh_data: &MeshData,
        canonical_nodes: &HashMap<u32, u32>,
        tol_dist: f64,
    ) -> Vec<(f64, (DVec3, DVec3))> {
        let mut datum_lines = Vec::new();
        for el in mesh_data.elements.iter().filter(|e| e.nodes.len() >= 3) {
            let pts: Vec<DVec3> = el
                .nodes
                .iter()
                .filter_map(|nid| canonical_nodes.get(nid).and_then(|cid| mesh_data.nodes.get(cid).copied()))
                .collect();

            if pts.len() < 3 {
                continue;
            }

            let v1 = pts[1] - pts[0];
            let v2 = pts[2] - pts[0];
            let norm = v1.cross(v2);
            let norm_len = norm.length();

            // Только вертикальные стены
            if norm_len < 1e-7 || (norm.z / norm_len).abs() > 0.15 {
                continue;
            }

            let n_len = pts.len();
            for i in 0..n_len {
                let pa = pts[i];
                let pb = pts[(i + 1) % n_len];
                if (pa.z - pb.z).abs() < tol_dist && (pb - pa).length_squared() > 1e-4 {
                    datum_lines.push((pa.z, (pa, pb)));
                }
            }
        }
        datum_lines
    }

    /// Извлечение уникальных высотных отметок плит перекрытий
    pub fn extract_slab_elevations(panels: &[MacroPanel], tol_dist: f64) -> Vec<f64> {
        let mut slab_z: Vec<f64> = Vec::new();
        for p in panels {
            if p.panel_type == PanelType::Slab {
                for poly in &p.polygons {
                    for pt in poly {
                        slab_z.push(pt[2]);
                    }
                }
            }
        }

        slab_z.sort_by(|a, b| a.total_cmp(b));
        let mut unique_levels: Vec<f64> = Vec::new();
        for z in slab_z {
            if unique_levels.is_empty()
                || (z - unique_levels[unique_levels.len() - 1]).abs() > tol_dist
            {
                unique_levels.push((z * 1000.0).round() / 1000.0);
            }
        }
        unique_levels
    }
}