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
        let &canon_id = coord_map.entry(key).or_insert(nid);
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
                // Пропускаем p_curr
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