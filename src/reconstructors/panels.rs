#![allow(dead_code, unused_imports, unused_variables, unused_assignments)]

use crate::config::ReconstructionConfig;
use crate::geometry::partition::extract_cutting_edge_nodes_for_slab;
use crate::geometry::utils::{
    clean_polygon_coords_3d, get_plane_basis, remove_short_edges_3d, snap_loop_to_wall_lines,
};
use crate::models::{ElementData, MacroPanel, MeshData, PanelType};
use geo::{Area, Contains, Intersects, LineString, Polygon};
use glam::DVec3;
use hashbrown::{HashMap, HashSet};
use rayon::prelude::*;
use std::collections::VecDeque;

pub struct PanelReconstructor;

#[derive(Clone)]
struct ParsedShell {
    id: u32,
    source_ids: Vec<u32>,
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
        let shell_elems: Vec<&ElementData> =
            mesh_data.elements.iter().filter(|e| e.is_shell()).collect();

        // 1. Извлечение опорных линий вертикальных стен для дотягивания плит
        let wall_segments_by_z =
            Self::extract_wall_datum_lines(mesh_data, canonical_nodes, config.tol_dist);

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

                if pts.len() != unique_nodes.len()
                    || pts.len() < 3
                    || pts.iter().any(|p| !p.is_finite())
                {
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

                // Canonical orientation, independent of FE winding.
                let major = if norm.x.abs() >= norm.y.abs() && norm.x.abs() >= norm.z.abs() {
                    norm.x
                } else if norm.y.abs() >= norm.z.abs() {
                    norm.y
                } else {
                    norm.z
                };
                if major < 0.0 {
                    norm = -norm;
                }
                if norm.z.abs() >= config.tol_angle.cos() {
                    norm = DVec3::Z;
                } else if norm.z.abs() <= config.tol_angle.sin() {
                    norm = DVec3::new(norm.x, norm.y, 0.0).normalize();
                }
                // Do not flatten warped elements beyond the allowed displacement.
                if pts
                    .iter()
                    .any(|p| norm.dot(*p - centroid).abs() > config.tol_dist)
                {
                    return None;
                }

                let d = -norm.dot(centroid);

                Some(ParsedShell {
                    id: el.id,
                    source_ids: vec![el.id],
                    stiff_id: el.stiff_id,
                    normal: norm,
                    d,
                    centroid,
                    unique_nodes,
                })
            })
            .collect();

        // Identical faces are coalesced; conflicting properties require review.
        let mut unique: Vec<ParsedShell> = Vec::new();
        let mut signatures: HashMap<Vec<u32>, usize> = HashMap::new();
        let mut conflicts = HashSet::new();
        for el in parsed_shells {
            let mut signature = el.unique_nodes.clone();
            signature.sort_unstable();
            if let Some(&index) = signatures.get(&signature) {
                if unique[index].stiff_id == el.stiff_id {
                    unique[index].source_ids.extend(el.source_ids);
                } else {
                    conflicts.insert(index);
                }
            } else {
                signatures.insert(signature, unique.len());
                unique.push(el);
            }
        }
        let mut parsed_shells: Vec<_> = unique
            .into_iter()
            .enumerate()
            .filter(|(i, _)| !conflicts.contains(i))
            .map(|(_, el)| el)
            .collect();
        // Order by height and distance to the model centre: independent of FE IDs,
        // translation and horizontal rotations except geometrically symmetric ties.
        let centre = if parsed_shells.is_empty() {
            DVec3::ZERO
        } else {
            parsed_shells.iter().map(|el| el.centroid).sum::<DVec3>() / parsed_shells.len() as f64
        };
        let quantum = config.weld_tol * 1e-6;
        let rank = |el: &ParsedShell| {
            let relative = el.centroid - centre;
            (
                (relative.z / quantum).round(),
                (relative.truncate().length() / quantum).round(),
            )
        };
        parsed_shells.sort_by(|a, b| {
            let ra = rank(a);
            let rb = rank(b);
            ra.0.total_cmp(&rb.0)
                .then(ra.1.total_cmp(&rb.1))
                .then(a.stiff_id.cmp(&b.stiff_id))
                .then(a.centroid.x.total_cmp(&b.centroid.x))
                .then(a.centroid.y.total_cmp(&b.centroid.y))
        });
        let source_ids: HashMap<u32, Vec<u32>> = parsed_shells
            .iter()
            .map(|el| (el.id, el.source_ids.clone()))
            .collect();

        // 3. Кластеризация элементов по плоскостям
        let mut clusters: Vec<PlaneCluster> = Vec::new();
        for el in parsed_shells {
            let nz = el.normal.z;
            let is_horiz = nz.abs() >= config.tol_angle.cos();
            let is_vert = nz.abs() <= config.tol_angle.sin();

            let mut matched_idx = None;
            for (idx, cl) in clusters.iter().enumerate() {
                if is_horiz && cl.is_horiz {
                    if let Some(cl_z) = cl.z {
                        if el
                            .unique_nodes
                            .iter()
                            .all(|id| (mesh_data.nodes[id].z - cl_z).abs() <= config.tol_dist)
                        {
                            matched_idx = Some(idx);
                            break;
                        }
                    }
                } else if is_vert && cl.is_vert {
                    let n2d_el = el.normal.truncate().normalize();
                    let n2d_cl = cl.normal.truncate().normalize();
                    if n2d_el.dot(n2d_cl) >= config.tol_angle.cos()
                        && el.unique_nodes.iter().all(|id| {
                            (cl.normal.dot(mesh_data.nodes[id]) + cl.d).abs() <= config.tol_dist
                        })
                    {
                        matched_idx = Some(idx);
                        break;
                    }
                } else if !is_horiz && !is_vert && !cl.is_horiz && !cl.is_vert {
                    if el.normal.dot(cl.normal) >= config.tol_angle.cos()
                        && el.unique_nodes.iter().all(|id| {
                            (cl.normal.dot(mesh_data.nodes[id]) + cl.d).abs() <= config.tol_dist
                        })
                    {
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
                        let ordered_nodes: Vec<u32> =
                            angle_nodes.into_iter().map(|(_, nid)| nid).collect();

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

                    if e_ids.len() == 2 {
                        for i in 0..e_ids.len() {
                            for j in (i + 1)..e_ids.len() {
                                if elem_by_id[&e_ids[i]].1 != elem_by_id[&e_ids[j]].1 {
                                    continue;
                                }
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
                        config,
                        cluster_norm,
                        cl.d,
                    );

                    if polygons.is_empty() {
                        continue;
                    }

                    let (polygons, filled_holes, filled_hole_area) =
                        fill_sliver_holes(polygons, u_axis, v_axis, config.simplify_tol);
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
                        filled_holes,
                        filled_hole_area,
                        fe_count: comp_elems.iter().map(|e| source_ids[&e.0].len()).sum(),
                        source_element_ids: comp_elems
                            .iter()
                            .flat_map(|e| source_ids[&e.0].iter().copied())
                            .collect(),
                        connected_panel_ids: vec![],
                    });
                    sub_id += 1;
                }

                local_panels
            })
            .collect();

        let mut panels: Vec<MacroPanel> = panels_results.into_iter().flatten().collect();
        for (i, panel) in panels.iter_mut().enumerate() {
            panel.id = (i + 1) as u32;
            panel.source_element_ids.sort_unstable();
        }
        panels
    }

    /// Трассировка контуров по Half-Edges с дотягиванием до линий стен
    fn extract_cycles_from_elements(
        elements: &[&(u32, u32, Vec<u32>)],
        mesh_data: &MeshData,
        u_axis: DVec3,
        v_axis: DVec3,
        wall_segments: &[(DVec3, DVec3)],
        config: &ReconstructionConfig,
        normal: DVec3,
        plane_d: f64,
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

        if dir_edge_count.values().any(|&count| count > 1) {
            return Vec::new();
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

        // Ambiguous/non-manifold boundaries must not produce arbitrary loops.
        let mut incoming: HashMap<u32, usize> = HashMap::new();
        for &(_, v) in &boundary_edges {
            *incoming.entry(v).or_default() += 1;
        }
        if out_edges
            .iter()
            .any(|(u, edges)| edges.len() != 1 || incoming.get(u) != Some(&1))
        {
            return Vec::new();
        }
        boundary_edges.sort_unstable();
        let mut visited_half_edges: HashSet<(u32, u32)> = HashSet::new();
        let mut result_polygons = Vec::new();
        let mut original_polygons = Vec::new();

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

                // Every exported ring lies on its declared plane.
                let planar: Vec<DVec3> = raw_3d
                    .iter()
                    .map(|p| *p - normal * (normal.dot(*p) + plane_d))
                    .collect();
                original_polygons.push(planar.iter().map(|p| p.to_array()).collect());
                let snapped = snap_loop_to_wall_lines(&planar, wall_segments, config.simplify_tol);
                let snapped: Vec<DVec3> = snapped
                    .iter()
                    .map(|p| *p - normal * (normal.dot(*p) + plane_d))
                    .collect();
                let cleaned = simplify_ring(&snapped, config.min_edge, config.simplify_tol);
                // Validate all rings together below. No silent deletion of small holes.
                result_polygons.push(cleaned.into_iter().map(|p| [p.x, p.y, p.z]).collect());
            }
        }

        order_valid_rings(result_polygons, u_axis, v_axis)
            .or_else(|| order_valid_rings(original_polygons, u_axis, v_axis))
            .unwrap_or_default()
    }

    /// Сбор опорных линий (отрезков) всех вертикальных стен с привязкой к отметкам Z
    fn extract_wall_datum_lines(
        mesh_data: &MeshData,
        canonical_nodes: &HashMap<u32, u32>,
        tol_dist: f64,
    ) -> Vec<(f64, (DVec3, DVec3))> {
        let mut datum_lines = Vec::new();
        for el in mesh_data.elements.iter().filter(|e| e.is_shell()) {
            let pts: Vec<DVec3> = el
                .nodes
                .iter()
                .filter_map(|nid| {
                    canonical_nodes
                        .get(nid)
                        .and_then(|cid| mesh_data.nodes.get(cid).copied())
                })
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
/// Remove a vertex only within a distance budget, retaining a simple polygon.
fn simplify_ring(pts: &[DVec3], min_edge: f64, deviation: f64) -> Vec<DVec3> {
    let mut result = pts.to_vec();
    loop {
        if result.len() <= 3 {
            break;
        }
        let mut remove = None;
        for i in 0..result.len() {
            let a = result[(i + result.len() - 1) % result.len()];
            let b = result[i];
            let c = result[(i + 1) % result.len()];
            let ac = c - a;
            if ac.length_squared() <= 1e-20 {
                continue;
            }
            let t = ((b - a).dot(ac) / ac.length_squared()).clamp(0.0, 1.0);
            let error = b.distance(a + t * ac);
            let short = a.distance(b).min(b.distance(c)) < min_edge;
            if error <= deviation && (short || error <= 1e-8) {
                let mut candidate = result.clone();
                candidate.remove(i);
                let bounded = pts.iter().all(|p| {
                    (0..candidate.len()).any(|j| {
                        let a = candidate[j];
                        let b = candidate[(j + 1) % candidate.len()];
                        let ab = b - a;
                        if ab.length_squared() <= 1e-20 {
                            return false;
                        }
                        let t = ((*p - a).dot(ab) / ab.length_squared()).clamp(0.0, 1.0);
                        p.distance(a + t * ab) <= deviation + 1e-12
                    })
                });
                if bounded {
                    remove = Some(i);
                    break;
                }
            }
        }
        match remove {
            Some(i) => {
                result.remove(i);
            }
            None => break,
        }
    }
    result
}

pub(crate) fn order_valid_rings(
    mut rings: Vec<Vec<[f64; 3]>>,
    u: DVec3,
    v: DVec3,
) -> Option<Vec<Vec<[f64; 3]>>> {
    let mut polygons = Vec::new();
    for ring in &rings {
        if ring.len() < 3 {
            return None;
        }
        let mut coords: Vec<(f64, f64)> = ring
            .iter()
            .map(|p| {
                let p = DVec3::from_array(*p);
                (p.dot(u), p.dot(v))
            })
            .collect();
        if coords.iter().any(|p| !p.0.is_finite() || !p.1.is_finite()) {
            return None;
        }
        coords.push(coords[0]);
        let line = LineString::from(coords);
        let segments: Vec<_> = line.lines().collect();
        for i in 0..segments.len() {
            if segments[i].start == segments[i].end {
                return None;
            }
            for j in i + 1..segments.len() {
                if j == i + 1 || (i == 0 && j == segments.len() - 1) {
                    continue;
                }
                if segments[i].intersects(&segments[j]) {
                    return None;
                }
            }
        }
        let polygon = Polygon::new(line, vec![]);
        if polygon.unsigned_area() <= 1e-12 {
            return None;
        }
        polygons.push(polygon);
    }
    let outer = (0..polygons.len()).max_by(|&a, &b| {
        polygons[a]
            .unsigned_area()
            .total_cmp(&polygons[b].unsigned_area())
    })?;
    rings.swap(0, outer);
    polygons.swap(0, outer);
    for i in 1..polygons.len() {
        if !polygons[0].contains(&polygons[i])
            || polygons[0].exterior().intersects(polygons[i].exterior())
        {
            return None;
        }
        for j in 1..i {
            if polygons[i].intersects(&polygons[j]) {
                return None;
            }
        }
    }
    for (i, ring) in rings.iter_mut().enumerate() {
        if (polygons[i].signed_area() > 0.0) != (i == 0) {
            ring.reverse();
        }
    }
    Some(rings)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::utils::canonicalize_nodes;

    fn square(x: f64, y: f64, size: f64) -> Vec<[f64; 3]> {
        vec![
            [x, y, 0.0],
            [x + size, y, 0.0],
            [x + size, y + size, 0.0],
            [x, y + size, 0.0],
        ]
    }
    #[test]
    fn small_hole_is_preserved_and_ordered_after_outer() {
        let rings = order_valid_rings(
            vec![square(1.0, 1.0, 0.3), square(0.0, 0.0, 4.0)],
            DVec3::X,
            DVec3::Y,
        )
        .unwrap();
        assert_eq!(rings.len(), 2);
        assert_eq!(rings[0], square(0.0, 0.0, 4.0));
    }
    #[test]
    fn crossed_touching_and_nested_holes_are_rejected() {
        let crossed = vec![[0., 0., 0.], [3., 2., 0.], [0., 2., 0.], [2., 0., 0.]];
        assert!(order_valid_rings(vec![crossed], DVec3::X, DVec3::Y).is_none());
        assert!(order_valid_rings(
            vec![square(0., 0., 4.), square(0., 1., 0.3)],
            DVec3::X,
            DVec3::Y
        )
        .is_none());
        assert!(order_valid_rings(
            vec![
                square(0., 0., 4.),
                square(1., 1., 2.),
                square(1.5, 1.5, 0.3)
            ],
            DVec3::X,
            DVec3::Y
        )
        .is_none());
    }
    fn mesh(slope: f64, different_stiffness: bool) -> MeshData {
        MeshData {
            nodes: HashMap::from_iter([
                (1, DVec3::new(0., 0., 0.)),
                (2, DVec3::new(2., 0., 2. * slope)),
                (3, DVec3::new(2., 2., 2. * slope)),
                (4, DVec3::new(0., 2., 0.)),
            ]),
            elements: vec![
                ElementData {
                    id: 1,
                    elem_type: 42,
                    stiff_id: 1,
                    nodes: vec![1, 2, 3],
                },
                ElementData {
                    id: 2,
                    elem_type: 42,
                    stiff_id: if different_stiffness { 2 } else { 1 },
                    nodes: vec![4, 3, 1],
                },
            ],
        }
    }
    #[test]
    fn reversed_shell_winding_merges_and_preserves_incline() {
        let mesh = mesh(0.3, false);
        let config = ReconstructionConfig::default();
        let panels = PanelReconstructor::reconstruct(
            &mesh,
            &canonicalize_nodes(&mesh.nodes, config.weld_tol),
            &config,
        );
        assert_eq!(panels.len(), 1);
        assert_eq!(panels[0].panel_type, PanelType::InclinedPanel);
        assert_eq!(panels[0].source_element_ids.len(), 2);
        let normal = DVec3::from_array(panels[0].plane_normal);
        for ring in &panels[0].polygons {
            for p in ring {
                assert!((normal.dot(DVec3::from_array(*p)) + panels[0].plane_d).abs() < 1e-10);
            }
        }
    }
    #[test]
    fn stiffness_boundary_is_retained() {
        let mesh = mesh(0.0, true);
        let config = ReconstructionConfig::default();
        let panels = PanelReconstructor::reconstruct(
            &mesh,
            &canonicalize_nodes(&mesh.nodes, config.weld_tol),
            &config,
        );
        assert_eq!(panels.len(), 2);
        assert_ne!(panels[0].stiffness_id, panels[1].stiffness_id);
    }
    #[test]
    fn duplicate_faces_are_not_exported_as_valid_surface() {
        let mut mesh = mesh(0.0, false);
        let mut duplicate = mesh.elements[0].clone();
        duplicate.id = 3;
        mesh.elements.push(duplicate);
        let config = ReconstructionConfig::default();
        let panels = PanelReconstructor::reconstruct(
            &mesh,
            &canonicalize_nodes(&mesh.nodes, config.weld_tol),
            &config,
        );
        assert_eq!(panels.len(), 1);
        assert_eq!(panels[0].source_element_ids.len(), 3);
        assert_eq!(panels[0].polygons[0].len(), 4);
    }
}

#[cfg(test)]
mod healing_tests {
    use super::*;
    use crate::geometry::utils::canonicalize_nodes;

    #[test]
    fn small_crack_is_welded_into_one_planar_surface() {
        let mesh = MeshData {
            nodes: HashMap::from_iter([
                (1, DVec3::new(0., 0., 0.)),
                (2, DVec3::new(2., 0., 0.)),
                (3, DVec3::new(2., 2., 0.)),
                (4, DVec3::new(0.005, 0., 0.002)),
                (5, DVec3::new(2.005, 2., 0.002)),
                (6, DVec3::new(0., 2., 0.002)),
            ]),
            elements: vec![
                ElementData {
                    id: 1,
                    elem_type: 42,
                    stiff_id: 1,
                    nodes: vec![1, 2, 3],
                },
                ElementData {
                    id: 2,
                    elem_type: 42,
                    stiff_id: 1,
                    nodes: vec![4, 5, 6],
                },
            ],
        };
        let mut config = ReconstructionConfig::default();
        config.weld_tol = 0.01;
        let panels = PanelReconstructor::reconstruct(
            &mesh,
            &canonicalize_nodes(&mesh.nodes, config.weld_tol),
            &config,
        );
        assert_eq!(panels.len(), 1);
        assert_eq!(panels[0].fe_count, 2);
        assert!(panels[0].polygons[0]
            .iter()
            .all(|p| (p[2] + panels[0].plane_d).abs() < 1e-12 && p[2].abs() <= 0.002));
    }

    #[test]
    fn narrow_notch_is_removed_within_budget() {
        let ring = vec![
            DVec3::ZERO,
            DVec3::new(1., 0., 0.),
            DVec3::new(1., 0.005, 0.),
            DVec3::new(1.005, 0., 0.),
            DVec3::new(2., 0., 0.),
            DVec3::new(2., 2., 0.),
            DVec3::new(0., 2., 0.),
        ];
        let cleaned = simplify_ring(&ring, 0.03, 0.01);
        assert_eq!(cleaned.len(), 4);
        assert!(order_valid_rings(
            vec![cleaned.iter().map(|p| p.to_array()).collect()],
            DVec3::X,
            DVec3::Y
        )
        .is_some());
    }

    #[test]
    fn annulus_reconstruction_keeps_small_opening() {
        let points = [
            (0., 0.),
            (4., 0.),
            (4., 4.),
            (0., 4.),
            (1., 1.),
            (1.3, 1.),
            (1.3, 1.3),
            (1., 1.3),
        ];
        let mesh = MeshData {
            nodes: points
                .iter()
                .enumerate()
                .map(|(i, &(x, y))| ((i + 1) as u32, DVec3::new(x, y, 0.)))
                .collect(),
            elements: [[1, 2, 6, 5], [2, 3, 7, 6], [3, 4, 8, 7], [4, 1, 5, 8]]
                .iter()
                .enumerate()
                .map(|(i, ids)| ElementData {
                    id: (i + 1) as u32,
                    elem_type: 44,
                    stiff_id: 1,
                    nodes: ids.to_vec(),
                })
                .collect(),
        };
        let config = ReconstructionConfig::default();
        let panels = PanelReconstructor::reconstruct(
            &mesh,
            &canonicalize_nodes(&mesh.nodes, config.weld_tol),
            &config,
        );
        assert_eq!(panels.len(), 1);
        assert_eq!(panels[0].polygons.len(), 2);
        assert_eq!(panels[0].fe_count, 4);
    }
}

/// Remove narrow closed voids within the healing budget, regardless of model IDs.
/// Width is a support-strip width, so a long thin slit is detected without deleting
/// a small but well-proportioned opening. The outer boundary is never removed here.
fn fill_sliver_holes(
    rings: Vec<Vec<[f64; 3]>>,
    u: DVec3,
    v: DVec3,
    tolerance: f64,
) -> (Vec<Vec<[f64; 3]>>, usize, f64) {
    let mut result = Vec::new();
    let mut count = 0;
    let mut area = 0.0;
    for (i, ring) in rings.into_iter().enumerate() {
        let points: Vec<glam::DVec2> = ring
            .iter()
            .map(|p| {
                let p = DVec3::from_array(*p);
                glam::DVec2::new(p.dot(u), p.dot(v))
            })
            .collect();
        let mut width = f64::INFINITY;
        for j in 0..points.len() {
            let edge = points[(j + 1) % points.len()] - points[j];
            if edge.length_squared() <= 1e-20 {
                continue;
            }
            let normal = glam::DVec2::new(-edge.y, edge.x).normalize();
            let min = points
                .iter()
                .map(|p| p.dot(normal))
                .fold(f64::INFINITY, f64::min);
            let max = points
                .iter()
                .map(|p| p.dot(normal))
                .fold(f64::NEG_INFINITY, f64::max);
            width = width.min(max - min);
        }
        if i > 0 && width <= tolerance {
            count += 1;
            // Subtract a local origin to avoid cancellation on large coordinates.
            let origin = points[0];
            area += (0..points.len())
                .map(|j| (points[j] - origin).perp_dot(points[(j + 1) % points.len()] - origin))
                .sum::<f64>()
                .abs()
                * 0.5;
        } else {
            result.push(ring);
        }
    }
    (result, count, area)
}

#[cfg(test)]
mod sliver_tests {
    use super::*;
    #[test]
    fn heals_long_thin_void_but_retains_small_regular_opening() {
        let outer = vec![[0., 0., 0.], [10., 0., 0.], [10., 10., 0.], [0., 10., 0.]];
        let slit = vec![[1., 1., 0.], [8., 1., 0.], [8., 1.001, 0.], [1., 1.001, 0.]];
        let hole = vec![[2., 2., 0.], [2.3, 2., 0.], [2.3, 2.3, 0.], [2., 2.3, 0.]];
        let (rings, count, area) =
            fill_sliver_holes(vec![outer, slit, hole], DVec3::X, DVec3::Y, 0.01);
        assert_eq!(rings.len(), 2);
        assert_eq!(count, 1);
        assert!((area - 0.007).abs() < 1e-12);
    }
}

#[cfg(test)]
mod ordering_tests {
    use super::*;
    use crate::geometry::utils::canonicalize_nodes;
    #[test]
    fn geometry_does_not_depend_on_element_traversal_order() {
        let mut mesh = MeshData {
            nodes: HashMap::from_iter([
                (1, DVec3::new(0., 0., 0.)),
                (2, DVec3::new(2., 0., 0.08)),
                (3, DVec3::new(2., 2., 0.08)),
                (4, DVec3::new(0., 2., 0.)),
            ]),
            elements: vec![
                ElementData {
                    id: 1,
                    elem_type: 42,
                    stiff_id: 1,
                    nodes: vec![1, 2, 3],
                },
                ElementData {
                    id: 2,
                    elem_type: 42,
                    stiff_id: 1,
                    nodes: vec![1, 3, 4],
                },
            ],
        };
        let config = ReconstructionConfig::default();
        let canonical = canonicalize_nodes(&mesh.nodes, config.weld_tol);
        let a = PanelReconstructor::reconstruct(&mesh, &canonical, &config);
        mesh.elements.reverse();
        let b = PanelReconstructor::reconstruct(&mesh, &canonical, &config);
        assert_eq!(
            serde_json::to_value(a).unwrap(),
            serde_json::to_value(b).unwrap()
        );
    }
}
