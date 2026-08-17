use crate::config::ReconstructionConfig;
use crate::geometry::utils::{clean_polygon_coords_3d, get_plane_basis};
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
    ordered_nodes: Vec<u32>,
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
    /// Высокоскоростная параллельная реконструкция плит и стен за O(N)
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

        // 1. Парсинг и локальная CCW-сортировка узлов каждого КЭ вокруг СОБСТВЕННОГО центроида
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

                // Собственный центроид данного элемента
                let elem_centroid = pts.iter().copied().sum::<DVec3>() / (pts.len() as f64);
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

                // Каноническая ориентация нормали
                if norm.x.abs() > 1e-4 {
                    if norm.x < 0.0 { norm = -norm; }
                } else if norm.y.abs() > 1e-4 {
                    if norm.y < 0.0 { norm = -norm; }
                } else if norm.z < 0.0 {
                    norm = -norm;
                }

                let d = -norm.dot(elem_centroid);
                let (u_axis, v_axis) = get_plane_basis(norm);

                // Сортировка узлов строго вокруг СОБСТВЕННОГО центра КЭ
                let mut angle_nodes: Vec<(f64, u32)> = pts
                    .iter()
                    .zip(unique_nodes)
                    .map(|(&p, nid)| {
                        let rel = p - elem_centroid;
                        let u = rel.dot(u_axis);
                        let v = rel.dot(v_axis);
                        (v.atan2(u), nid)
                    })
                    .collect();

                angle_nodes.sort_by(|a, b| a.0.total_cmp(&b.0));
                let ordered_nodes: Vec<u32> = angle_nodes.into_iter().map(|(_, nid)| nid).collect();

                Some(ParsedShell {
                    id: el.id,
                    stiff_id: el.stiff_id,
                    normal: norm,
                    d,
                    centroid: elem_centroid,
                    ordered_nodes,
                })
            })
            .collect();

        // 2. Кластеризация элементов по плоскостям
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
                clusters.push(PlaneCluster {
                    is_horiz,
                    is_vert,
                    z: if is_horiz { Some(el.centroid.z) } else { None },
                    normal: el.normal,
                    d: el.d,
                    elements: vec![el],
                });
            }
        }

        // 3. Параллельное извлечение макропанелей через Half-Edges (Rayon)
        let panels_results: Vec<Vec<MacroPanel>> = clusters
            .into_par_iter()
            .enumerate()
            .map(|(cl_idx, cl)| {
                let norm = cl.normal;
                let elem_by_id: HashMap<u32, &ParsedShell> = cl.elements.iter().map(|e| (e.id, e)).collect();

                // Построение графа смежности КЭ
                let mut edge_to_elems: HashMap<(u32, u32), Vec<u32>> = HashMap::new();
                for el in &cl.elements {
                    let n_len = el.ordered_nodes.len();
                    for i in 0..n_len {
                        let a = el.ordered_nodes[i];
                        let b = el.ordered_nodes[(i + 1) % n_len];
                        let edge = if a < b { (a, b) } else { (b, a) };
                        edge_to_elems.entry(edge).or_default().push(el.id);
                    }
                }

                let mut elem_adj: HashMap<u32, Vec<u32>> = HashMap::new();
                for e_ids in edge_to_elems.values() {
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

                for start_el in &cl.elements {
                    if visited_elems.contains(&start_el.id) {
                        continue;
                    }

                    let mut comp_elems = Vec::new();
                    let mut queue = VecDeque::new();

                    queue.push_back(start_el.id);
                    visited_elems.insert(start_el.id);

                    while let Some(curr_id) = queue.pop_front() {
                        if let Some(&el) = elem_by_id.get(&curr_id) {
                            comp_elems.push(el);
                            if let Some(neighbors) = elem_adj.get(&curr_id) {
                                for &nbr in neighbors {
                                    if visited_elems.insert(nbr) {
                                        queue.push_back(nbr);
                                    }
                                }
                            }
                        }
                    }

                    // 4. Трассировка контуров через Half-Edges для компоненты
                    let polygons = Self::extract_cycles_from_elements(&comp_elems, mesh_data);
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
                        stiffness_id: comp_elems[0].stiff_id,
                        plane_normal: [norm.x, norm.y, norm.z],
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

    /// Трассировка замкнутых циклов (периметр и проемы) по направленным полуребрам за O(N)
    fn extract_cycles_from_elements(
        elements: &[&ParsedShell],
        mesh_data: &MeshData,
    ) -> Vec<Vec<[f64; 3]>> {
        // Подсчет направленных полуребер
        let mut dir_edge_count: HashMap<(u32, u32), u32> = HashMap::new();
        for el in elements {
            let n_len = el.ordered_nodes.len();
            for i in 0..n_len {
                let u = el.ordered_nodes[i];
                let v = el.ordered_nodes[(i + 1) % n_len];
                *dir_edge_count.entry((u, v)).or_insert(0) += 1;
            }
        }

        // Фильтрация граничных полуребер (отсутствует противоположное ребро)
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

        // Сборка замкнутых контуров
        let mut visited_half_edges: HashSet<(u32, u32)> = HashSet::new();
        let mut result_polygons = Vec::new();

        for &(start_u, start_v) in &boundary_edges {
            if visited_half_edges.contains(&(start_u, start_v)) {
                continue;
            }

            let mut cycle = vec![start_u];
            let mut curr_u = start_u;
            let mut curr_v = start_v;
            visited_half_edges.insert((curr_u, curr_v));

            let mut closed = false;
            while let Some(neighbors) = out_edges.get(&curr_v) {
                cycle.push(curr_v);
                if curr_v == start_u {
                    closed = true;
                    break;
                }

                let mut next_node = None;
                for &nxt in neighbors {
                    if !visited_half_edges.contains(&(curr_v, nxt)) {
                        next_node = Some(nxt);
                        break;
                    }
                }

                if let Some(nxt) = next_node {
                    visited_half_edges.insert((curr_v, nxt));
                    curr_u = curr_v;
                    curr_v = nxt;
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

                // Удаление промежуточных коллинеарных вершин
                let clean_3d = clean_polygon_coords_3d(&raw_3d, 0.9995);
                if clean_3d.len() >= 3 {
                    result_polygons.push(clean_3d.into_iter().map(|p| [p.x, p.y, p.z]).collect());
                }
            }
        }

        result_polygons
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