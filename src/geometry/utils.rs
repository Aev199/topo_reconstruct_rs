#![allow(dead_code, unused_imports, unused_variables)]

use glam::DVec3;
use hashbrown::HashMap;

/// Канонизация узлов: объединение близко расположенных узлов в единые индексы
pub fn canonicalize_nodes(nodes: &HashMap<u32, DVec3>, precision: u32) -> HashMap<u32, u32> {
    let factor = 10f64.powi(precision as i32);
    let mut coord_map: HashMap<(i64, i64, i64), u32> = HashMap::with_capacity(nodes.len());
    let mut canonical_map = HashMap::with_capacity(nodes.len());

    for (&nid, &pt) in nodes {
        let key = (
            (pt.x * factor).round() as i64,
            (pt.y * factor).round() as i64,
            (pt.z * factor).round() as i64,
        );
        let canon_id = *coord_map.entry(key).or_insert(nid);
        canonical_map.insert(nid, canon_id);
    }
    canonical_map
}

/// Построение ортонормированного 2D-базиса (u, v), перпендикулярного вектору нормали
pub fn get_plane_basis(normal: DVec3) -> (DVec3, DVec3) {
    let norm = normal.normalize();
    let u = if norm.z.abs() < 0.9 {
        norm.cross(DVec3::Z).normalize()
    } else {
        norm.cross(DVec3::X).normalize()
    };
    let v = norm.cross(u).normalize();
    (u, v)
}

/// Проекция 2D-точки P на отрезок AB с возвратом спроецированной точки и расстояния
pub fn project_point_to_segment_2d(p: DVec3, a: DVec3, b: DVec3) -> (DVec3, f64) {
    let ab = DVec3::new(b.x - a.x, b.y - a.y, 0.0);
    let len_sq = ab.length_squared();
    if len_sq < 1e-8 {
        let d = (DVec3::new(p.x, p.y, 0.0) - DVec3::new(a.x, a.y, 0.0)).length();
        return (DVec3::new(a.x, a.y, p.z), d);
    }
    let ap = DVec3::new(p.x - a.x, p.y - a.y, 0.0);
    let t = (ap.dot(ab) / len_sq).clamp(0.0, 1.0);
    let proj_2d = DVec3::new(a.x + t * ab.x, a.y + t * ab.y, p.z);
    let dist = (DVec3::new(p.x, p.y, 0.0) - DVec3::new(proj_2d.x, proj_2d.y, 0.0)).length();
    (proj_2d, dist)
}

/// Дотягивание (snapping) вершин контура плиты до опорных линий стен
pub fn snap_loop_to_wall_lines(
    loop_pts: &[DVec3],
    wall_segments: &[(DVec3, DVec3)],
    snap_tol: f64,
) -> Vec<DVec3> {
    if loop_pts.is_empty() || wall_segments.is_empty() {
        return loop_pts.to_vec();
    }

    let mut snapped = Vec::with_capacity(loop_pts.len());
    for &pt in loop_pts {
        let mut best_pt = pt;
        let mut min_dist = snap_tol;

        for &(w_a, w_b) in wall_segments {
            let (proj, dist) = project_point_to_segment_2d(pt, w_a, w_b);
            if dist < min_dist {
                min_dist = dist;
                best_pt = proj;
            }
        }
        snapped.push(best_pt);
    }
    snapped
}

/// Удаление паразитных микросегментов (длиной < min_len)
pub fn remove_short_edges_3d(pts: &[DVec3], min_len: f64) -> Vec<DVec3> {
    if pts.len() < 3 {
        return pts.to_vec();
    }
    let mut result = Vec::with_capacity(pts.len());
    let n = pts.len();
    for i in 0..n {
        let p_curr = pts[i];
        let p_next = pts[(i + 1) % n];
        if (p_next - p_curr).length() >= min_len {
            result.push(p_curr);
        }
    }
    if result.len() >= 3 {
        result
    } else {
        pts.to_vec()
    }
}

/// Удаление промежуточных коллинеарных узлов КЭ на прямых гранях в 3D
pub fn clean_polygon_coords_3d(pts: &[DVec3], tol_collinear: f64) -> Vec<DVec3> {
    if pts.len() < 3 {
        return pts.to_vec();
    }

    let mut current = pts.to_vec();
    let mut changed = true;

    while changed && current.len() >= 3 {
        changed = false;
        let n = current.len();
        let mut next_pts = Vec::with_capacity(n);

        for i in 0..n {
            let p_prev = current[(i + n - 1) % n];
            let p_curr = current[i];
            let p_next = current[(i + 1) % n];

            let v1 = p_curr - p_prev;
            let v2 = p_next - p_curr;
            let l1 = v1.length();
            let l2 = v2.length();

            if l1 < 1e-4 || l2 < 1e-4 {
                changed = true;
                continue;
            }

            let d1 = v1 / l1;
            let d2 = v2 / l2;

            // Если векторы сонаправлены (скалярное произведение близко к 1.0) — точка на прямой
            if d1.dot(d2) > tol_collinear {
                changed = true;
            } else {
                next_pts.push(p_curr);
            }
        }
        current = next_pts;
    }

    if current.len() >= 3 {
        current
    } else {
        Vec::new()
    }
}