//! Connected planar hypotheses from source shells, before contour construction.
use super::PlaneFrame;
use crate::input::{ElementData, MeshData};
use glam::DVec3;
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Serialize)]
pub struct Policy {
    pub angle: f64,
    pub distance: f64,
    pub precision: f64,
}
#[derive(Debug, Clone, Serialize)]
pub struct Patch {
    pub plane: PlaneFrame,
    pub source_elements: Vec<u32>,
    pub source_nodes: Vec<u32>,
    pub stiffness_regions: BTreeMap<u32, Vec<u32>>,
    pub maximum_deviation: f64,
}
#[derive(Debug, Clone, Serialize)]
pub struct Rejection {
    pub element: u32,
    pub reason: String,
}
#[derive(Debug, Clone, Serialize)]
pub struct Report {
    pub policy: Policy,
    pub patches: Vec<Patch>,
    pub rejected: Vec<Rejection>,
    pub unmerged_components: Vec<Vec<u32>>,
    pub non_shell_elements: Vec<u32>,
}
struct Facet {
    nodes: Vec<u32>,
    normal: DVec3,
    center: DVec3,
}

pub(super) fn fit(points: &[DVec3], reference: DVec3) -> Option<(DVec3, f64)> {
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

/// Recover the perimeter of a supported convex FE facet independently of
/// tensor-product or cyclic input numbering. Shared by recognition and assembly.
pub(super) fn ordered_facet_nodes(
    mesh: &MeshData,
    e: &ElementData,
    frame: &PlaneFrame,
    precision: f64,
) -> Option<Vec<u32>> {
    if !e.is_shell() || e.nodes.iter().collect::<BTreeSet<_>>().len() != e.nodes.len() {
        return None;
    }
    let points: Option<Vec<_>> = e
        .nodes
        .iter()
        .map(|n| mesh.nodes.get(n).copied().filter(|p| p.is_finite()))
        .collect();
    let points = points?;
    let center = points.iter().copied().sum::<DVec3>() / points.len() as f64;
    let normal = DVec3::from_array(frame.normal);
    let local = PlaneFrame::new(center.to_array(), frame.normal).ok()?;
    let mut nodes = e.nodes.clone();
    nodes.sort_by(|a, b| {
        let p = local.project(mesh.nodes[a].to_array());
        let q = local.project(mesh.nodes[b].to_array());
        p[1].atan2(p[0]).total_cmp(&q[1].atan2(q[0]))
    });
    for i in 0..nodes.len() {
        let a = mesh.nodes[&nodes[i]];
        let b = mesh.nodes[&nodes[(i + 1) % nodes.len()]];
        let c = mesh.nodes[&nodes[(i + 2) % nodes.len()]];
        if a.distance(b) <= precision || (b - a).cross(c - b).dot(normal) <= precision * precision {
            return None;
        }
    }
    Some(nodes)
}

pub fn recognize(mesh: &MeshData, policy: &Policy) -> Result<Report, &'static str> {
    if !policy.angle.is_finite()
        || policy.angle <= 0.0
        || policy.angle >= std::f64::consts::FRAC_PI_4
        || !policy.distance.is_finite()
        || !policy.precision.is_finite()
        || policy.precision <= 0.0
        || policy.distance < policy.precision
    {
        return Err("invalid plane policy");
    }
    let mut result = Report {
        policy: policy.clone(),
        patches: vec![],
        rejected: vec![],
        unmerged_components: vec![],
        non_shell_elements: vec![],
    };
    let mut facets = BTreeMap::<u32, Facet>::new();
    let mut elements = BTreeMap::<u32, &ElementData>::new();
    for e in &mesh.elements {
        if elements.insert(e.id, e).is_some() {
            return Err("duplicate source element ID");
        }
        if !e.is_shell() {
            result.non_shell_elements.push(e.id);
            continue;
        }
        let points: Option<Vec<_>> = e
            .nodes
            .iter()
            .map(|n| mesh.nodes.get(n).copied().filter(|p| p.is_finite()))
            .collect();
        let facet = points.and_then(|points| {
            if e.nodes.iter().copied().collect::<BTreeSet<_>>().len() != e.nodes.len() {
                return None;
            }
            let (n, _) = fit(&points, DVec3::Z)?;
            let center = points.iter().copied().sum::<DVec3>() / points.len() as f64;
            if points
                .iter()
                .any(|p| n.dot(*p - center).abs() > policy.distance)
            {
                return None;
            }
            let frame = PlaneFrame::new(center.to_array(), n.to_array()).ok()?;
            let nodes = ordered_facet_nodes(mesh, e, &frame, policy.precision)?;
            Some(Facet {
                nodes,
                normal: n,
                center,
            })
        });
        if let Some(facet) = facet {
            facets.insert(e.id, facet);
        } else {
            result.rejected.push(Rejection {
                element: e.id,
                reason: "invalid_degenerate_or_excessively_warped_shell".into(),
            });
        }
    }
    let mut edges = BTreeMap::<[u32; 2], Vec<u32>>::new();
    for (&id, f) in &facets {
        for i in 0..f.nodes.len() {
            let a = f.nodes[i];
            let b = f.nodes[(i + 1) % f.nodes.len()];
            edges.entry([a.min(b), a.max(b)]).or_default().push(id);
        }
    }
    let mut neighbors = BTreeMap::<u32, BTreeSet<u32>>::new();
    for owners in edges.values() {
        for (i, &a) in owners.iter().enumerate() {
            for &b in &owners[i + 1..] {
                let x = &facets[&a];
                let y = &facets[&b];
                if x.normal.dot(y.normal).abs() >= policy.angle.cos()
                    && y.nodes
                        .iter()
                        .all(|n| x.normal.dot(mesh.nodes[n] - x.center).abs() <= policy.distance)
                    && x.nodes
                        .iter()
                        .all(|n| y.normal.dot(mesh.nodes[n] - y.center).abs() <= policy.distance)
                {
                    neighbors.entry(a).or_default().insert(b);
                    neighbors.entry(b).or_default().insert(a);
                }
            }
        }
    }
    let build = |ids: &BTreeSet<u32>| -> Option<Patch> {
        let nodes: BTreeSet<_> = ids
            .iter()
            .flat_map(|id| facets[id].nodes.iter().copied())
            .collect();
        let points: Vec<_> = nodes.iter().map(|n| mesh.nodes[n]).collect();
        let (normal, _) = fit(&points, facets[ids.first()?].normal)?;
        let center = points.iter().copied().sum::<DVec3>() / points.len() as f64;
        let maximum_deviation = points
            .iter()
            .map(|p| normal.dot(*p - center).abs())
            .fold(0.0_f64, f64::max);
        if maximum_deviation > policy.distance
            || ids
                .iter()
                .any(|id| normal.dot(facets[id].normal).abs() < policy.angle.cos())
        {
            return None;
        }
        let mut stiffness_regions = BTreeMap::<u32, Vec<u32>>::new();
        for &id in ids {
            stiffness_regions
                .entry(elements[&id].stiff_id)
                .or_default()
                .push(id);
        }
        Some(Patch {
            plane: PlaneFrame::new(center.to_array(), normal.to_array()).ok()?,
            source_elements: ids.iter().copied().collect(),
            source_nodes: nodes.into_iter().collect(),
            stiffness_regions,
            maximum_deviation,
        })
    };
    let mut seen = BTreeSet::new();
    for &seed in facets.keys() {
        if seen.contains(&seed) {
            continue;
        }
        let mut group = BTreeSet::new();
        let mut stack = vec![seed];
        while let Some(id) = stack.pop() {
            if !group.insert(id) {
                continue;
            }
            if let Some(ids) = neighbors.get(&id) {
                stack.extend(ids);
            }
        }
        seen.extend(group.iter().copied());
        if let Some(patch) = build(&group) {
            result.patches.push(patch);
        } else {
            result
                .unmerged_components
                .push(group.iter().copied().collect());
            for id in group {
                result
                    .patches
                    .push(build(&BTreeSet::from([id])).ok_or("failed to retain source facet")?);
            }
        }
    }
    result.rejected.sort_by_key(|r| r.element);
    result.non_shell_elements.sort_unstable();
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn policy() -> Policy {
        Policy {
            angle: 0.02,
            distance: 0.01,
            precision: 1e-8,
        }
    }
    fn mesh() -> MeshData {
        MeshData {
            nodes: [
                [0., 0., 0.],
                [1., 0., 0.],
                [0., 1., 0.],
                [1., 1., 0.],
                [2., 0., 0.],
                [2., 1., 0.],
            ]
            .iter()
            .enumerate()
            .map(|(i, p)| (i as u32 + 1, DVec3::from_array(*p)))
            .collect(),
            elements: vec![
                ElementData {
                    id: 1,
                    elem_type: 44,
                    stiff_id: 1,
                    nodes: vec![1, 2, 3, 4],
                },
                ElementData {
                    id: 2,
                    elem_type: 44,
                    stiff_id: 2,
                    nodes: vec![2, 5, 4, 6],
                },
            ],
        }
    }
    #[test]
    fn shared_plane_preserves_property_regions() {
        let r = recognize(&mesh(), &policy()).unwrap();
        assert_eq!(r.patches.len(), 1);
        assert_eq!(r.patches[0].stiffness_regions.len(), 2);
        assert_eq!(r.patches[0].source_elements, vec![1, 2]);
    }
    #[test]
    fn disconnected_close_planes_are_separate() {
        let mut m = mesh();
        m.elements[1].nodes = vec![7, 8, 9, 10];
        for (id, p) in [
            (7, [0., 0., 0.005]),
            (8, [1., 0., 0.005]),
            (9, [0., 1., 0.005]),
            (10, [1., 1., 0.005]),
        ] {
            m.nodes.insert(id, DVec3::from_array(p));
        }
        assert_eq!(recognize(&m, &policy()).unwrap().patches.len(), 2);
    }
    #[test]
    fn folded_surface_is_not_flattened() {
        let mut m = mesh();
        m.nodes.insert(5, DVec3::new(1., 0., 1.));
        m.nodes.insert(6, DVec3::new(1., 1., 1.));
        assert_eq!(recognize(&m, &policy()).unwrap().patches.len(), 2);
    }
    #[test]
    fn warped_and_missing_nodes_are_reported() {
        let mut m = mesh();
        m.nodes.get_mut(&1).unwrap().z = 0.5;
        m.nodes.remove(&6);
        let r = recognize(&m, &policy()).unwrap();
        assert_eq!(r.rejected.len(), 2);
        assert!(r.patches.is_empty());
    }
    #[test]
    fn gradual_curvature_does_not_chain_into_one_plane() {
        let mut m = MeshData::default();
        for i in 0..30_u32 {
            let t = i as f64 * 0.01;
            for j in 0..2_u32 {
                m.nodes.insert(
                    i * 2 + j + 1,
                    DVec3::new(t.sin() * 10., j as f64, (1. - t.cos()) * 10.),
                );
            }
        }
        for i in 0..29_u32 {
            m.elements.push(ElementData {
                id: i + 1,
                elem_type: 44,
                stiff_id: 1,
                nodes: vec![i * 2 + 1, i * 2 + 3, i * 2 + 2, i * 2 + 4],
            });
        }
        let r = recognize(&m, &policy()).unwrap();
        assert_eq!(r.patches.len(), 29);
        assert_eq!(r.unmerged_components.len(), 1);
        assert!(r.rejected.is_empty());
    }

    #[test]
    fn numbering_order_does_not_change_membership() {
        let mut m = mesh();
        let a = recognize(&m, &policy()).unwrap();
        m.elements.reverse();
        let b = recognize(&m, &policy()).unwrap();
        assert_eq!(a.patches[0].source_elements, b.patches[0].source_elements);
    }
}
