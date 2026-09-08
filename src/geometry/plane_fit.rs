//! Plane proposals fitted to source-linked bar incidences; acceptance is separate.
use crate::{
    config::ReconstructionConfig,
    models::{MacroBar, MacroPanel, MeshData},
};
use glam::DVec3;
use hashbrown::HashMap;
use std::collections::{BTreeMap, BTreeSet};

// Jacobi diagonalization of a symmetric 3x3 covariance matrix. Reject rank-one
// observations: a single straight axis cannot determine a surface orientation.
fn fit(points: &[DVec3], reference: DVec3) -> Option<(DVec3, f64)> {
    if points.len() < 3 {
        return None;
    }
    let center = points.iter().copied().sum::<DVec3>() / points.len() as f64;
    let mut a = [[0.0; 3]; 3];
    let mut v = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    for p in points {
        let d = (*p - center).to_array();
        for i in 0..3 {
            for j in 0..3 {
                a[i][j] += d[i] * d[j];
            }
        }
    }
    let scale = a[0][0] + a[1][1] + a[2][2];
    if scale <= 1e-20 {
        return None;
    }
    for _ in 0..40 {
        let (p, q) = [(0, 1), (0, 2), (1, 2)]
            .into_iter()
            .max_by(|&(p, q), &(i, j)| a[p][q].abs().total_cmp(&a[i][j].abs()))
            .unwrap();
        if a[p][q].abs() < scale * 1e-14 {
            break;
        }
        let angle = 0.5 * (2.0 * a[p][q]).atan2(a[q][q] - a[p][p]);
        let (c, s) = (angle.cos(), angle.sin());
        let (app, aqq, apq) = (a[p][p], a[q][q], a[p][q]);
        a[p][p] = c * c * app - 2.0 * c * s * apq + s * s * aqq;
        a[q][q] = s * s * app + 2.0 * c * s * apq + c * c * aqq;
        a[p][q] = 0.0;
        a[q][p] = 0.0;
        for k in 0..3 {
            if k != p && k != q {
                let (x, y) = (a[k][p], a[k][q]);
                a[k][p] = c * x - s * y;
                a[p][k] = a[k][p];
                a[k][q] = s * x + c * y;
                a[q][k] = a[k][q];
            }
            let (x, y) = (v[k][p], v[k][q]);
            v[k][p] = c * x - s * y;
            v[k][q] = s * x + c * y;
        }
    }
    let mut order = [0, 1, 2];
    order.sort_by(|&i, &j| a[i][i].total_cmp(&a[j][j]));
    if a[order[1]][order[1]] < scale * 1e-10 {
        return None;
    }
    let k = order[0];
    let mut n = DVec3::new(v[0][k], v[1][k], v[2][k]).normalize();
    if n.dot(reference) < 0.0 {
        n = -n;
    }
    Some((n, -n.dot(center)))
}

pub(crate) fn proposals(
    bars: &[MacroBar],
    panels: &[MacroPanel],
    mesh: &MeshData,
    canonical: &HashMap<u32, u32>,
    config: &ReconstructionConfig,
) -> Vec<(Vec<MacroPanel>, usize)> {
    let mut groups: Vec<Vec<usize>> = vec![];
    for i in 0..panels.len() {
        if let Some(group) = groups.iter_mut().find(|g| {
            let j = g[0];
            DVec3::from_array(panels[i].plane_normal)
                .distance(DVec3::from_array(panels[j].plane_normal))
                < 1e-8
                && (panels[i].plane_d - panels[j].plane_d).abs() < 1e-8
        }) {
            group.push(i);
        } else {
            groups.push(vec![i]);
        }
    }
    let mut owners = BTreeMap::new();
    for (g, indices) in groups.iter().enumerate() {
        for &i in indices {
            for &id in &panels[i].source_element_ids {
                owners.insert(id, g);
            }
        }
    }
    let mut node_groups: BTreeMap<u32, BTreeSet<usize>> = BTreeMap::new();
    for e in &mesh.elements {
        if let Some(&g) = owners.get(&e.id) {
            for n in &e.nodes {
                node_groups
                    .entry(canonical.get(n).copied().unwrap_or(*n))
                    .or_default()
                    .insert(g);
            }
        }
    }
    let mut samples = vec![Vec::<DVec3>::new(); groups.len()];
    for b in bars {
        let a = DVec3::from_array(b.start_point);
        let d = DVec3::from_array(b.end_point) - a;
        let anchors = std::iter::once((b.start_node_id, 0.0))
            .chain(std::iter::once((b.end_node_id, 1.0)))
            .chain(
                b.constraints
                    .iter()
                    .filter_map(|c| c.source_node_id.map(|n| (n, c.t))),
            );
        for (n, t) in anchors {
            if let Some(gs) = node_groups.get(&n) {
                for &g in gs {
                    let p = a + t * d;
                    let panel = &panels[groups[g][0]];
                    if (DVec3::from_array(panel.plane_normal).dot(p) + panel.plane_d).abs()
                        <= config.joint_tol
                        && samples[g].iter().all(|q| p.distance(*q) > 1e-7)
                    {
                        samples[g].push(p);
                    }
                }
            }
        }
    }
    let mut normals: Vec<_> = panels.iter().map(|p| p.plane_normal).collect();
    let mut offsets: Vec<_> = panels.iter().map(|p| p.plane_d).collect();
    let mut active = BTreeSet::new();
    let mut ranked = Vec::new();
    for (g, indices) in groups.iter().enumerate() {
        let p = &panels[indices[0]];
        let old = DVec3::from_array(p.plane_normal);
        let Some((n, d)) = fit(&samples[g], old) else {
            continue;
        };
        let before: f64 = samples[g]
            .iter()
            .map(|q| (old.dot(*q) + p.plane_d).powi(2))
            .sum();
        let after: f64 = samples[g].iter().map(|q| (n.dot(*q) + d).powi(2)).sum();
        if n.dot(old) < config.tol_angle.cos() || before - after < 1e-14 || n.distance(old) < 1e-9 {
            continue;
        }
        if indices
            .iter()
            .flat_map(|&i| panels[i].polygons.iter().flatten())
            .any(|q| (n.dot(DVec3::from_array(*q)) + d).abs() > config.joint_tol)
        {
            continue;
        }
        ranked.push((g, before - after));
        for &i in indices {
            normals[i] = n.to_array();
            offsets[i] = d;
            active.insert(i);
        }
    }
    let mut candidates = Vec::new();
    ranked.sort_by(|a, b| b.1.total_cmp(&a.1));
    for &(g, _) in ranked.iter().take(6) {
        let mut ns: Vec<_> = panels.iter().map(|p| p.plane_normal).collect();
        let mut ds: Vec<_> = panels.iter().map(|p| p.plane_d).collect();
        for &i in &groups[g] {
            ns[i] = normals[i];
            ds[i] = offsets[i];
        }
        if let Ok(result) = super::topology::transformed_planes_checked(panels, &ns, &ds, config) {
            candidates.push((result, groups[g].len()));
        }
    }
    while !active.is_empty() {
        match super::topology::transformed_planes_checked(panels, &normals, &offsets, config) {
            Ok(result) => {
                candidates.insert(0, (result, active.len()));
                return candidates;
            }
            Err(affected) => {
                let mut changed = false;
                for group in &groups {
                    if group.iter().any(|i| affected.contains(i)) {
                        for &i in group {
                            changed |= active.remove(&i);
                            normals[i] = panels[i].plane_normal;
                            offsets[i] = panels[i].plane_d;
                        }
                    }
                }
                if !changed {
                    return candidates;
                }
            }
        }
    }
    candidates
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn fits_tilted_plane_and_preserves_orientation() {
        let points = vec![
            DVec3::new(0., 0., 0.),
            DVec3::new(2., 0., 0.02),
            DVec3::new(0., 2., 0.),
            DVec3::new(2., 2., 0.02),
        ];
        let (n, d) = fit(&points, DVec3::Z).unwrap();
        assert!(n.z > 0.99);
        assert!(points.iter().all(|p| (n.dot(*p) + d).abs() < 1e-10));
        let shifted: Vec<_> = points
            .iter()
            .map(|p| *p + DVec3::new(1000., -500., 10.))
            .collect();
        let (m, e) = fit(&shifted, DVec3::Z).unwrap();
        assert!(n.distance(m) < 1e-10);
        assert!(shifted.iter().all(|p| (m.dot(*p) + e).abs() < 1e-10));
    }
    #[test]
    fn rejects_collinear_and_coincident_samples() {
        assert!(fit(&[DVec3::ZERO, DVec3::X, DVec3::X * 2.], DVec3::Z).is_none());
        assert!(fit(&[DVec3::ZERO; 4], DVec3::Z).is_none());
    }
}
