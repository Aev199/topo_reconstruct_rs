//! Surface topology preview. Engineering closure tolerance is distinct from
//! numerical planarity. Mechanical ties and mesh readiness are not inferred.
use super::{frame, planes, Model, PlaneFrame};
use crate::input::MeshData;
use glam::DVec3;
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Serialize)]
pub struct Policy {
    pub closure_tolerance: f64,
    pub precision: f64,
    pub minimum_edge: f64,
}
#[derive(Debug, Serialize)]
pub struct Issue {
    pub patch: usize,
    pub reason: String,
}
#[derive(Debug, Serialize)]
pub struct Report {
    pub policy: Policy,
    /// This stage assembles surfaces only; axes/property-region geometry still
    /// require reconciliation. Never treat this preview as an export gate.
    pub export_ready: bool,
    pub all_surface_patches_built: bool,
    pub preview: Model,
    pub vertex_source_nodes: Vec<u32>,
    pub surface_source_patches: Vec<usize>,
    pub issues: Vec<Issue>,
    pub maximum_closure_movement: f64,
    pub rejected_vertices: BTreeMap<u32, String>,
    pub support_representatives: Vec<usize>,
}

fn boundary(
    mesh: &MeshData,
    ids: &[u32],
    plane: &PlaneFrame,
    precision: f64,
) -> Result<Vec<Vec<u32>>, &'static str> {
    let elements: BTreeMap<_, _> = mesh.elements.iter().map(|e| (e.id, e)).collect();
    let mut counts = BTreeMap::<[u32; 2], usize>::new();
    for id in ids {
        let e = elements.get(id).ok_or("missing_source_element")?;
        let nodes =
            planes::ordered_facet_nodes(mesh, e, plane, precision).ok_or("invalid_source_face")?;
        for i in 0..nodes.len() {
            let (a, b) = (nodes[i], nodes[(i + 1) % nodes.len()]);
            *counts.entry([a.min(b), a.max(b)]).or_default() += 1;
        }
    }
    if counts.values().any(|&n| n > 2) {
        return Err("nonmanifold_source_edges");
    }
    let mut adjacency = BTreeMap::<u32, Vec<u32>>::new();
    for ([a, b], count) in counts {
        if count == 1 {
            adjacency.entry(a).or_default().push(b);
            adjacency.entry(b).or_default().push(a);
        }
    }
    if adjacency.is_empty() || adjacency.values().any(|v| v.len() != 2) {
        return Err("ambiguous_boundary");
    }
    let mut remaining: BTreeSet<_> = adjacency.keys().copied().collect();
    let mut rings = vec![];
    while let Some(&start) = remaining.first() {
        let mut ring = vec![];
        let (mut previous, mut current) = (start, start);
        loop {
            if !remaining.remove(&current) {
                return Err("invalid_boundary_cycle");
            }
            ring.push(current);
            let next = adjacency[&current]
                .iter()
                .copied()
                .find(|n| *n != previous)
                .ok_or("invalid_boundary_cycle")?;
            previous = current;
            current = next;
            if current == start {
                break;
            }
        }
        rings.push(ring);
    }
    Ok(rings)
}

/// Minimum movement onto the intersection of fixed support planes, using an
/// orthonormal constraint basis. Dependent inconsistent planes are rejected.
fn intersection(point: DVec3, planes: &[&PlaneFrame], precision: f64) -> Option<DVec3> {
    let mut basis: Vec<(DVec3, f64)> = vec![];
    for plane in planes {
        let mut n = DVec3::from_array(plane.normal);
        let mut rhs = -plane.distance(point.to_array());
        for &(q, t) in &basis {
            let a = n.dot(q);
            n -= a * q;
            rhs -= a * t;
        }
        let length = n.length();
        if length > 1e-10 && basis.len() < 3 {
            basis.push((n / length, rhs / length));
        }
    }
    let result = point + basis.iter().map(|(n, t)| *n * *t).sum::<DVec3>();
    (result.is_finite()
        && planes
            .iter()
            .all(|p| p.distance(result.to_array()).abs() <= precision))
    .then_some(result)
}

pub fn assemble(
    mesh: &MeshData,
    source: &frame::Report,
    policy: &Policy,
) -> Result<Report, &'static str> {
    if !policy.closure_tolerance.is_finite() || policy.closure_tolerance < policy.precision {
        return Err("invalid closure tolerance");
    }
    let mut model =
        Model::new(policy.precision, policy.minimum_edge).map_err(|_| "invalid assembly policy")?;
    if source.candidate_planes.len() != source.surfaces.len()
        || source.candidate_points.len() != source.node_ids.len()
    {
        return Err("incomplete frame proposal");
    }
    let lookup: BTreeMap<_, _> = source
        .node_ids
        .iter()
        .enumerate()
        .map(|(i, &n)| (n, i))
        .collect();
    let mut issues = vec![];
    let mut rings = BTreeMap::new();
    for (i, s) in source.surfaces.iter().enumerate() {
        match boundary(
            mesh,
            &s.source_elements,
            &source.candidate_planes[i],
            policy.precision,
        ) {
            Ok(r) => {
                rings.insert(i, r);
            }
            Err(reason) => issues.push(Issue {
                patch: i,
                reason: reason.into(),
            }),
        }
    }
    // Include every support owning a boundary node, including a patch whose
    // own boundary failed. Never silently disconnect it from valid neighbors.
    let boundary_nodes: BTreeSet<_> = rings.values().flatten().flatten().copied().collect();
    let mut owners = BTreeMap::<u32, Vec<usize>>::new();
    for (i, s) in source.surfaces.iter().enumerate() {
        for &index in &s.nodes {
            let id = *source.node_ids.get(index).ok_or("invalid surface node")?;
            if boundary_nodes.contains(&id) {
                owners.entry(id).or_default().push(i);
            }
        }
    }
    // Connected nearly coplanar source patches can use one support. Validate
    // every candidate point against the chosen support, preventing chain drift.
    let mut support_representatives: Vec<_> = (0..source.surfaces.len()).collect();
    let mut order: Vec<_> = (0..source.surfaces.len()).collect();
    order.sort_by_key(|&i| std::cmp::Reverse(source.surfaces[i].nodes.len()));
    let mut assigned = BTreeSet::new();
    for master in order {
        if !assigned.insert(master) {
            continue;
        }
        let plane = &source.candidate_planes[master];
        let mut connected: BTreeSet<_> = source.surfaces[master].nodes.iter().copied().collect();
        loop {
            let mut changed = false;
            for j in 0..source.surfaces.len() {
                if assigned.contains(&j) {
                    continue;
                }
                let surface = &source.surfaces[j];
                if surface.nodes.iter().any(|n| connected.contains(n))
                    && DVec3::from_array(plane.normal)
                        .dot(DVec3::from_array(source.candidate_planes[j].normal))
                        .abs()
                        >= source.policy.angle.cos()
                    && surface.nodes.iter().all(|&n| {
                        plane.distance(source.candidate_points[n]).abs() <= policy.closure_tolerance
                    })
                {
                    assigned.insert(j);
                    support_representatives[j] = master;
                    connected.extend(surface.nodes.iter().copied());
                    changed = true;
                }
            }
            if !changed {
                break;
            }
        }
    }
    let mut rejected_vertices = BTreeMap::new();
    let mut vertices = BTreeMap::new();
    let mut vertex_source_nodes = vec![];
    let mut maximum_closure_movement = 0.0_f64;
    for (&id, supports) in &owners {
        let i = lookup[&id];
        let p = DVec3::from_array(source.candidate_points[i]);
        let planes: Vec<_> = supports
            .iter()
            .map(|&s| &source.candidate_planes[support_representatives[s]])
            .collect();
        let Some(q) = intersection(p, &planes, policy.precision) else {
            rejected_vertices.insert(id, "inconsistent_supports".into());
            continue;
        };
        let movement = p.distance(q);
        let reference = *mesh.nodes.get(&id).ok_or("missing reference node")?;
        // Cumulative movement from immutable input, including local axis budgets.
        let mut budget = source.policy.maximum_movement;
        for axis in &source.axes {
            if axis.anchors.iter().any(|a| a.node == i) {
                let a = mesh.nodes[&source.node_ids[axis.endpoints[0]]];
                let b = mesh.nodes[&source.node_ids[axis.endpoints[1]]];
                budget = budget.min(source.policy.relative_movement * a.distance(b));
            }
        }
        if movement > policy.closure_tolerance || q.distance(reference) > budget + policy.precision
        {
            rejected_vertices.insert(
                id,
                format!(
                    "movement: closure={movement}, total={}, budget={budget}",
                    q.distance(reference)
                ),
            );
            continue;
        }
        maximum_closure_movement = maximum_closure_movement.max(movement);
        vertices.insert(
            id,
            model
                .add_vertex(q.to_array())
                .map_err(|_| "invalid closed vertex")?,
        );
        vertex_source_nodes.push(id);
    }
    let mut surface_source_patches = vec![];
    for (patch, mut loops) in rings {
        if loops.iter().flatten().any(|n| !vertices.contains_key(n)) {
            issues.push(Issue {
                patch,
                reason: "support_intersection_or_movement_budget".into(),
            });
            continue;
        }
        let plane = &source.candidate_planes[support_representatives[patch]];
        let area = |ring: &Vec<u32>| {
            let uv: Vec<_> = ring
                .iter()
                .map(|n| plane.project(model.vertices[vertices[n]]))
                .collect();
            (0..uv.len())
                .map(|i| {
                    let a = uv[i];
                    let b = uv[(i + 1) % uv.len()];
                    a[0] * b[1] - a[1] * b[0]
                })
                .sum::<f64>()
                .abs()
        };
        loops.sort_by(|a, b| area(b).total_cmp(&area(a)));
        let plane_id = model.add_plane(plane.clone());
        let mapped = loops
            .iter()
            .map(|r| r.iter().map(|n| vertices[n]).collect())
            .collect();
        match model.add_surface(
            plane_id,
            mapped,
            source.surfaces[patch].source_elements.clone(),
        ) {
            Ok(_) => surface_source_patches.push(patch),
            Err(error) => issues.push(Issue {
                patch,
                reason: format!("contour_{error:?}"),
            }),
        }
    }
    Ok(Report {
        policy: policy.clone(),
        export_ready: false,
        all_surface_patches_built: issues.is_empty(),
        preview: model,
        vertex_source_nodes,
        surface_source_patches,
        issues,
        maximum_closure_movement,
        rejected_vertices,
        support_representatives,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn tensor_numbering_preserves_outer_boundary_and_hole() {
        use crate::input::ElementData;
        let mut mesh = MeshData::default();
        for y in 0..4 {
            for x in 0..4 {
                mesh.nodes.insert(
                    y * 4 + x + 1,
                    DVec3::new(x as f64 + 10., y as f64 - 20., 7.),
                );
            }
        }
        for y in 0..3 {
            for x in 0..3 {
                if x == 1 && y == 1 {
                    continue;
                }
                let a = y * 4 + x + 1;
                mesh.elements.push(ElementData {
                    id: mesh.elements.len() as u32 + 1,
                    elem_type: 44,
                    stiff_id: 9,
                    nodes: vec![a, a + 1, a + 5, a + 4],
                });
            }
        }
        let plane = PlaneFrame::new([0., 0., 7.], [0., 0., 1.]).unwrap();
        let ids: Vec<_> = mesh.elements.iter().map(|e| e.id).collect();
        let original = boundary(&mesh, &ids, &plane, 1e-8).unwrap();
        assert_eq!(original.len(), 2);
        assert_eq!(
            original.iter().map(Vec::len).collect::<BTreeSet<_>>(),
            BTreeSet::from([4, 12])
        );
        for e in &mut mesh.elements {
            e.nodes.swap(2, 3);
        }
        assert_eq!(boundary(&mesh, &ids, &plane, 1e-8).unwrap(), original);
        let mut model = Model::new(1e-8, 0.01).unwrap();
        let p = model.add_plane(plane);
        let vertices: BTreeMap<_, _> = mesh
            .nodes
            .iter()
            .map(|(&id, p)| (id, model.add_vertex(p.to_array()).unwrap()))
            .collect();
        let mut rings = original;
        rings.sort_by_key(|r| std::cmp::Reverse(r.len()));
        model
            .add_surface(
                p,
                rings
                    .iter()
                    .map(|r| r.iter().map(|n| vertices[n]).collect())
                    .collect(),
                ids,
            )
            .unwrap();
        assert_eq!(model.surfaces[0].contours.len(), 2);
        assert_eq!(model.surfaces[0].source_elements.len(), 8);
    }

    #[test]
    fn assembles_common_edge_despite_small_unaccepted_frame_residual() {
        use super::super::{planes, recognize};
        use crate::input::ElementData;
        let mut mesh = MeshData::default();
        for (i, p) in [
            [0., 0., 0.],
            [0., 1., 0.],
            [1., 1., 0.],
            [1., 0., 0.],
            [0., 1., 1.],
            [0., 0., 1.],
        ]
        .iter()
        .enumerate()
        {
            mesh.nodes.insert(i as u32 + 1, DVec3::from_array(*p));
        }
        mesh.elements = vec![
            ElementData {
                id: 1,
                elem_type: 44,
                stiff_id: 10,
                nodes: vec![1, 2, 3, 4],
            },
            ElementData {
                id: 2,
                elem_type: 44,
                stiff_id: 20,
                nodes: vec![1, 2, 5, 6],
            },
        ];
        let axes = recognize::recognize(
            &mesh,
            &recognize::Policy {
                angle: 0.02,
                line_tolerance: 0.01,
                numerical_precision: 1e-8,
            },
        )
        .unwrap();
        let planes = planes::recognize(
            &mesh,
            &planes::Policy {
                angle: 0.02,
                distance: 0.01,
                precision: 1e-8,
            },
        )
        .unwrap();
        let mut f = frame::solve(
            &mesh,
            &axes,
            &planes,
            &frame::Policy {
                up: [0., 0., 1.],
                angle: 0.02,
                maximum_movement: 0.15,
                relative_movement: 0.05,
                minimum_length: 0.03,
                residual_tolerance: 1e-7,
                iterations: 100,
            },
        )
        .unwrap();
        f.accepted = false;
        f.candidate_points[0] = [0.0003, 0., 0.0002];
        let p = Policy {
            closure_tolerance: 0.001,
            precision: 1e-7,
            minimum_edge: 0.001,
        };
        let result = assemble(&mesh, &f, &p).unwrap();
        assert!(result.all_surface_patches_built, "{:?}", result.issues);
        assert!(!result.export_ready);
        assert_eq!(result.preview.surfaces.len(), 2);
        assert_eq!(result.preview.vertices.len(), 6);
        assert_eq!(result.preview.edges.len(), 7);
        assert!(result.maximum_closure_movement > 0.0003);
        for surface in &result.preview.surfaces {
            for ring in &surface.boundaries {
                for edge in ring {
                    for v in result.preview.edges[edge.edge] {
                        assert!(
                            result.preview.planes[surface.plane]
                                .distance(result.preview.vertices[v])
                                .abs()
                                < 1e-7
                        );
                    }
                }
            }
        }
        f.candidate_points[0] = [0.002, 0., 0.002];
        let blocked = assemble(&mesh, &f, &p).unwrap();
        assert!(!blocked.all_surface_patches_built);
        assert!(blocked.preview.surfaces.is_empty());
        assert_eq!(blocked.issues.len(), 2);
    }
    #[test]
    fn closes_shared_point_on_two_planes_and_rejects_parallel_gap() {
        let a = PlaneFrame::new([0., 0., 0.], [0., 0., 1.]).unwrap();
        let b = PlaneFrame::new([0., 0., 0.], [1., 0., 0.]).unwrap();
        let p = intersection(DVec3::new(0.0003, 2., 0.0002), &[&a, &b], 1e-9).unwrap();
        assert!(p.distance(DVec3::new(0., 2., 0.)) < 1e-9);
        let c = PlaneFrame::new([0., 0., 0.0001], [0., 0., 1.]).unwrap();
        assert!(intersection(p, &[&a, &c], 1e-9).is_none());
    }
    #[test]
    fn redundant_planes_do_not_move_point_twice() {
        let a = PlaneFrame::new([0., 0., 0.], [0., 0., 1.]).unwrap();
        let p = intersection(DVec3::new(1., 2., 0.0003), &[&a, &a], 1e-9).unwrap();
        assert_eq!(p, DVec3::new(1., 2., 0.));
    }
}
