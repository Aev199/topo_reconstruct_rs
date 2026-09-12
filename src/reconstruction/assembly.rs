//! Surface topology preview. Engineering closure tolerance is distinct from
//! numerical planarity. Mechanical ties and mesh readiness are not inferred.
mod holes;
use super::{frame, planes, Model, PlaneFrame};
use crate::input::MeshData;
use glam::DVec3;
pub use holes::{HoleConstraint, HoleNodeChange, HoleOutcome, HoleRecovery};
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Serialize)]
pub struct Policy {
    /// Maximum normal deviation when reconciling support planes.
    pub closure_tolerance: f64,
    /// Additional displacement allowed when closing shared junctions.
    pub junction_movement_limit: f64,
    pub precision: f64,
    pub minimum_edge: f64,
}
#[derive(Debug, Serialize)]
pub struct Issue {
    pub patch: usize,
    pub source_elements: Vec<u32>,
    pub reason: String,
    pub boundary_source_nodes: Vec<Vec<u32>>,
}
#[derive(Debug, Serialize)]
pub struct RegionSplit {
    pub patch: usize,
    pub stiffness: u32,
    pub source_elements: Vec<u32>,
    pub parts: Vec<Vec<u32>>,
}
#[derive(Debug, Serialize)]
pub struct Report {
    pub policy: Policy,
    /// This stage assembles surface property regions; axes and inter-surface
    /// intersections still require reconciliation. This is not an export gate.
    pub export_ready: bool,
    pub all_surface_patches_built: bool,
    pub preview: Model,
    pub vertex_source_nodes: Vec<u32>,
    pub surface_source_patches: Vec<usize>,
    pub surface_stiffness: Vec<u32>,
    pub pinched_region_splits: Vec<RegionSplit>,
    pub hole_recovery: Vec<HoleRecovery>,
    pub issues: Vec<Issue>,
    pub maximum_closure_movement: f64,
    pub rejected_vertices: BTreeMap<u32, String>,
    pub support_representatives: Vec<usize>,
    pub support_offset_projection_applied: bool,
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

/// Open a pinched boundary by isolating its incident face fans along existing
/// source edges. This partitions material, never fills a hole or moves a node.
/// Every successful recursion strictly reduces the source group; unresolved
/// nonmanifold cases still reach the normal rejecting boundary validator.
fn split_pinched_regions(ids: &[u32], facets: &BTreeMap<u32, Vec<u32>>) -> Vec<Vec<u32>> {
    let mut edges = BTreeMap::<[u32; 2], Vec<u32>>::new();
    for &id in ids {
        let ns = &facets[&id];
        for i in 0..ns.len() {
            let (a, b) = (ns[i], ns[(i + 1) % ns.len()]);
            edges.entry([a.min(b), a.max(b)]).or_default().push(id);
        }
    }
    if edges.values().any(|owners| owners.len() > 2) {
        return vec![ids.to_vec()];
    }
    let mut degree = BTreeMap::<u32, usize>::new();
    for (edge, owners) in &edges {
        if owners.len() == 1 {
            for &n in edge {
                *degree.entry(n).or_default() += 1;
            }
        }
    }
    let pinch: BTreeSet<_> = degree
        .into_iter()
        .filter_map(|(n, d)| (d > 2).then_some(n))
        .collect();
    if pinch.is_empty() {
        return vec![ids.to_vec()];
    }
    let incident: BTreeSet<_> = ids
        .iter()
        .copied()
        .filter(|id| facets[id].iter().any(|n| pinch.contains(n)))
        .collect();
    let mut neighbors = BTreeMap::<u32, Vec<u32>>::new();
    for owners in edges.values() {
        if let [a, b] = owners.as_slice() {
            if incident.contains(a) == incident.contains(b) {
                neighbors.entry(*a).or_default().push(*b);
                neighbors.entry(*b).or_default().push(*a);
            }
        }
    }
    let mut remaining: BTreeSet<_> = ids.iter().copied().collect();
    let mut groups = vec![];
    while let Some(&seed) = remaining.first() {
        let mut stack = vec![seed];
        let mut group = vec![];
        while let Some(id) = stack.pop() {
            if remaining.remove(&id) {
                group.push(id);
                stack.extend(neighbors.get(&id).into_iter().flatten().copied());
            }
        }
        group.sort_unstable();
        groups.push(group);
    }
    if groups.len() == 1 {
        return groups;
    }
    groups
        .iter()
        .flat_map(|g| split_pinched_regions(g, facets))
        .collect()
}

/// Split property regions by edge connectivity; a point contact does not
/// turn two separate material areas into one polygon with an invalid hole.
fn property_regions(
    mesh: &MeshData,
    source: &frame::Report,
    precision: f64,
    splits: &mut Vec<RegionSplit>,
) -> Result<Vec<(usize, u32, Vec<u32>)>, &'static str> {
    let elements: BTreeMap<_, _> = mesh.elements.iter().map(|e| (e.id, e)).collect();
    let mut result = vec![];
    for (patch, surface) in source.surfaces.iter().enumerate() {
        let represented: Vec<_> = surface
            .stiffness_regions
            .values()
            .flatten()
            .copied()
            .collect();
        if represented.len() != represented.iter().collect::<BTreeSet<_>>().len()
            || represented.iter().copied().collect::<BTreeSet<_>>()
                != surface.source_elements.iter().copied().collect()
        {
            return Err("invalid property coverage");
        }
        for (&stiffness, ids) in &surface.stiffness_regions {
            let mut facets = BTreeMap::new();
            let mut edges = BTreeMap::<[u32; 2], Vec<u32>>::new();
            let mut neighbors = BTreeMap::<u32, BTreeSet<u32>>::new();
            for &id in ids {
                let e = elements.get(&id).ok_or("missing property element")?;
                if e.stiff_id != stiffness {
                    return Err("property mismatch");
                }
                neighbors.entry(id).or_default();
                let ns = planes::ordered_facet_nodes(
                    mesh,
                    e,
                    &source.candidate_planes[patch],
                    precision,
                )
                .ok_or("invalid property facet")?;
                facets.insert(id, ns.clone());
                for i in 0..ns.len() {
                    let (a, b) = (ns[i], ns[(i + 1) % ns.len()]);
                    edges.entry([a.min(b), a.max(b)]).or_default().push(id);
                }
            }
            for owners in edges.values() {
                for &a in owners {
                    for &b in owners {
                        if a != b {
                            neighbors.entry(a).or_default().insert(b);
                        }
                    }
                }
            }
            let mut remaining: BTreeSet<_> = ids.iter().copied().collect();
            while let Some(&seed) = remaining.first() {
                let mut stack = vec![seed];
                let mut group = vec![];
                while let Some(id) = stack.pop() {
                    if remaining.remove(&id) {
                        group.push(id);
                        stack.extend(&neighbors[&id]);
                    }
                }
                group.sort_unstable();
                let parts = split_pinched_regions(&group, &facets);
                if parts.len() > 1 {
                    splits.push(RegionSplit {
                        patch,
                        stiffness,
                        source_elements: group,
                        parts: parts.clone(),
                    });
                }
                result.extend(parts.into_iter().map(|part| (patch, stiffness, part)));
            }
        }
    }
    Ok(result)
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

/// Enforce concurrence by projecting support offsets onto the linear
/// compatibility constraints. Normals remain fixed; no per-surface node copies.
fn concurrent_supports(planes: &[PlaneFrame], junctions: &[Vec<usize>]) -> Vec<PlaneFrame> {
    if planes.is_empty() {
        return vec![];
    }
    let center = DVec3::from_array(planes[0].origin);
    let offsets: Vec<_> = planes
        .iter()
        .map(|p| DVec3::from_array(p.normal).dot(DVec3::from_array(p.origin) - center))
        .collect();
    let mut constraints: Vec<Vec<f64>> = vec![];
    for junction in junctions {
        let mut basis: Vec<(DVec3, Vec<f64>)> = vec![];
        for &i in junction {
            let mut n = DVec3::from_array(planes[i].normal);
            let mut row = vec![0.; planes.len()];
            row[i] = 1.;
            for (q, r) in &basis {
                let a = n.dot(*q);
                n -= a * *q;
                for (v, b) in row.iter_mut().zip(r) {
                    *v -= a * b;
                }
            }
            let length = n.length();
            if length > 1e-10 && basis.len() < 3 {
                for v in &mut row {
                    *v /= length;
                }
                basis.push((n / length, row));
            } else {
                // Reorthogonalize to avoid amplifying repeated junction rows.
                for _ in 0..2 {
                    for q in &constraints {
                        let a: f64 = row.iter().zip(q).map(|(a, b)| a * b).sum();
                        for (v, b) in row.iter_mut().zip(q) {
                            *v -= a * b;
                        }
                    }
                }
                let length = row.iter().map(|v| v * v).sum::<f64>().sqrt();
                if length > 1e-10 {
                    for v in &mut row {
                        *v /= length;
                    }
                    constraints.push(row);
                }
            }
        }
    }
    let mut corrected = offsets.clone();
    for row in constraints {
        let error: f64 = row.iter().zip(&offsets).map(|(a, b)| a * b).sum();
        for (d, a) in corrected.iter_mut().zip(row) {
            *d -= error * a;
        }
    }
    planes
        .iter()
        .enumerate()
        .map(|(i, p)| {
            let n = DVec3::from_array(p.normal);
            let mut result = p.clone();
            result.origin =
                (DVec3::from_array(p.origin) + n * (corrected[i] - offsets[i])).to_array();
            result
        })
        .collect()
}

fn movement_budget(mesh: &MeshData, source: &frame::Report, index: usize) -> f64 {
    let mut budget = source.policy.maximum_movement;
    for axis in &source.axes {
        if axis.anchors.iter().any(|a| a.node == index) {
            let a = mesh.nodes[&source.node_ids[axis.endpoints[0]]];
            let b = mesh.nodes[&source.node_ids[axis.endpoints[1]]];
            budget = budget.min(source.policy.relative_movement * a.distance(b));
        }
    }
    budget
}

fn ring_area(ring: &[u32], plane: &PlaneFrame, point: impl Fn(u32) -> [f64; 3]) -> f64 {
    let uv: Vec<_> = ring.iter().map(|&n| plane.project(point(n))).collect();
    (0..uv.len())
        .map(|i| {
            let (a, b, o) = (uv[i], uv[(i + 1) % uv.len()], uv[0]);
            (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
        })
        .sum::<f64>()
        .abs()
        / 2.
}

pub fn assemble(
    mesh: &MeshData,
    source: &frame::Report,
    policy: &Policy,
) -> Result<Report, &'static str> {
    if !policy.closure_tolerance.is_finite()
        || policy.closure_tolerance < policy.precision
        || !policy.junction_movement_limit.is_finite()
        || policy.junction_movement_limit < policy.precision
    {
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
    let mut pinched_region_splits = vec![];
    let regions = property_regions(mesh, source, policy.precision, &mut pinched_region_splits)?;
    for (i, (patch, _, ids)) in regions.iter().enumerate() {
        match boundary(
            mesh,
            ids,
            &source.candidate_planes[*patch],
            policy.precision,
        ) {
            Ok(mut r) => {
                // Fix exterior/hole roles from the immutable source geometry.
                let area = |ring: &Vec<u32>| {
                    ring_area(ring, &source.candidate_planes[*patch], |n| {
                        mesh.nodes[&n].to_array()
                    })
                };
                r.sort_by(|a, b| area(b).total_cmp(&area(a)));
                rings.insert(i, r);
            }
            Err(reason) => issues.push(Issue {
                patch: *patch,
                source_elements: ids.clone(),
                reason: reason.into(),
                boundary_source_nodes: vec![],
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
    let junctions: Vec<Vec<usize>> = owners
        .values()
        .map(|ids| {
            ids.iter()
                .map(|&i| support_representatives[i])
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect()
        })
        .collect();
    let proposal = concurrent_supports(&source.candidate_planes, &junctions);
    let support_offsets_adjusted = source.surfaces.iter().enumerate().all(|(i, s)| {
        s.nodes.iter().all(|&n| {
            proposal[support_representatives[i]]
                .distance(source.candidate_points[n])
                .abs()
                <= policy.closure_tolerance
        })
    });
    let closed_supports = if support_offsets_adjusted {
        proposal
    } else {
        source.candidate_planes.clone()
    };
    let mut rejected_vertices = BTreeMap::new();
    let mut closed_points = BTreeMap::new();
    for (&id, supports) in &owners {
        let i = lookup[&id];
        let p = DVec3::from_array(source.candidate_points[i]);
        let planes: Vec<_> = supports
            .iter()
            .map(|&s| &closed_supports[support_representatives[s]])
            .collect();
        let Some(q) = intersection(p, &planes, policy.precision) else {
            rejected_vertices.insert(id, "inconsistent_supports".into());
            continue;
        };
        let movement = p.distance(q);
        let reference = *mesh.nodes.get(&id).ok_or("missing reference node")?;
        // Cumulative movement from immutable input, including local axis budgets.
        let budget = movement_budget(mesh, source, i);
        if movement > policy.junction_movement_limit
            || q.distance(reference) > budget + policy.precision
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
        closed_points.insert(id, q);
    }
    let hole_recovery = holes::recover(
        &mut closed_points,
        &holes::Context {
            mesh,
            source,
            policy,
            regions: &regions,
            rings: &rings,
            owners: &owners,
            representatives: &support_representatives,
            supports: &closed_supports,
        },
    );
    let mut vertices = BTreeMap::new();
    let mut vertex_source_nodes = vec![];
    let mut maximum_closure_movement = 0.0_f64;
    for (&id, &q) in &closed_points {
        maximum_closure_movement = maximum_closure_movement
            .max(q.distance(DVec3::from_array(source.candidate_points[lookup[&id]])));
        vertices.insert(
            id,
            model
                .add_vertex(q.to_array())
                .map_err(|_| "invalid closed vertex")?,
        );
        vertex_source_nodes.push(id);
    }
    let mut surface_source_patches = vec![];
    let mut surface_stiffness = vec![];
    for (region, loops) in rings {
        let (patch, stiffness, ids) = &regions[region];
        let patch = *patch;
        if loops.iter().flatten().any(|n| !vertices.contains_key(n)) {
            issues.push(Issue {
                patch,
                source_elements: ids.clone(),
                reason: "support_intersection_or_movement_budget".into(),
                boundary_source_nodes: loops.clone(),
            });
            continue;
        }
        let plane = &closed_supports[support_representatives[patch]];
        let area = |ring: &Vec<u32>| ring_area(ring, plane, |n| model.vertices[vertices[&n]]);
        let degenerate_hole = loops
            .iter()
            .skip(1)
            .any(|ring| area(ring) <= policy.precision * policy.precision);
        let plane_id = model.add_plane(plane.clone());
        let mapped = loops
            .iter()
            .map(|r| r.iter().map(|n| vertices[n]).collect())
            .collect();
        match model.add_surface(plane_id, mapped, ids.clone()) {
            Ok(_) => {
                surface_source_patches.push(patch);
                surface_stiffness.push(*stiffness);
            }
            Err(error) => issues.push(Issue {
                patch,
                source_elements: ids.clone(),
                reason: if degenerate_hole {
                    "degenerate_hole_after_closure".into()
                } else {
                    format!("contour_{error:?}")
                },
                boundary_source_nodes: loops,
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
        surface_stiffness,
        pinched_region_splits,
        hole_recovery,
        issues,
        maximum_closure_movement,
        rejected_vertices,
        support_representatives,
        support_offset_projection_applied: support_offsets_adjusted,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    fn planar_frame(mesh: &MeshData, up: DVec3) -> frame::Report {
        use super::super::recognize;
        let axes = recognize::recognize(
            mesh,
            &recognize::Policy {
                angle: 0.02,
                line_tolerance: 0.01,
                numerical_precision: 1e-8,
            },
        )
        .unwrap();
        let surfaces = planes::recognize(
            mesh,
            &planes::Policy {
                angle: 0.02,
                distance: 0.01,
                precision: 1e-8,
            },
        )
        .unwrap();
        frame::solve(
            mesh,
            &axes,
            &surfaces,
            &frame::Policy {
                up: up.to_array(),
                angle: 0.02,
                maximum_movement: 0.15,
                relative_movement: 0.05,
                minimum_length: 0.03,
                residual_tolerance: 1e-7,
                iterations: 100,
            },
        )
        .unwrap()
    }

    #[test]
    fn pinched_opening_is_partitioned_without_filling_or_duplicate_sources() {
        use crate::input::ElementData;
        let policy = Policy {
            closure_tolerance: 0.001,
            junction_movement_limit: 0.05,
            precision: 1e-7,
            minimum_edge: 0.001,
        };
        let mut mesh = MeshData::default();
        for y in 0..4 {
            for x in 0..4 {
                mesh.nodes
                    .insert(1 + x + 4 * y, DVec3::new(x as f64, y as f64, 0.));
            }
        }
        for y in 0..3 {
            for x in 0..3 {
                // An interior void touches the exterior at just one vertex.
                if (x == 0 && y == 0) || (x == 1 && y == 1) {
                    continue;
                }
                let a = 1 + x + 4 * y;
                mesh.elements.push(ElementData {
                    id: mesh.elements.len() as u32 + 1,
                    elem_type: 44,
                    stiff_id: 9,
                    nodes: vec![a, a + 1, a + 4, a + 5],
                });
            }
        }
        let f = planar_frame(&mesh, DVec3::Z);
        let ids: Vec<_> = mesh.elements.iter().map(|e| e.id).collect();
        assert_eq!(
            boundary(&mesh, &ids, &f.candidate_planes[0], 1e-7),
            Err("ambiguous_boundary")
        );
        let base = assemble(&mesh, &f, &policy).unwrap();
        assert!(base.all_surface_patches_built, "{:?}", base.issues);
        assert_eq!(base.pinched_region_splits.len(), 1);
        assert!(base.preview.surfaces.len() > 1);
        let mut sources: Vec<_> = base
            .preview
            .surfaces
            .iter()
            .flat_map(|s| s.source_elements.iter().copied())
            .collect();
        sources.sort_unstable();
        assert_eq!(sources, ids);
        let area: f64 = base
            .preview
            .surfaces
            .iter()
            .map(|s| {
                s.contours
                    .iter()
                    .enumerate()
                    .map(|(j, r)| {
                        let a = (0..r.len())
                            .map(|i| {
                                let (a, b) = (r[i], r[(i + 1) % r.len()]);
                                a[0] * b[1] - a[1] * b[0]
                            })
                            .sum::<f64>()
                            .abs()
                            / 2.;
                        if j == 0 {
                            a
                        } else {
                            -a
                        }
                    })
                    .sum::<f64>()
            })
            .sum();
        assert!((area - 7.).abs() < 1e-9);
        assert_eq!(
            base.vertex_source_nodes.iter().filter(|&&n| n == 6).count(),
            1
        );
        let mut uses = BTreeMap::<usize, usize>::new();
        for e in base
            .preview
            .surfaces
            .iter()
            .flat_map(|s| s.boundaries.iter().flatten())
        {
            *uses.entry(e.edge).or_default() += 1;
        }
        assert!(uses.values().any(|&n| n == 2));
        assert!(base.surface_stiffness.iter().all(|&s| s == 9));
        // Run the whole recognition and assembly after a rigid transform and
        // reversed node/element numbering, not merely the graph helper.
        let q = glam::DQuat::from_axis_angle(DVec3::new(1., 2., 3.).normalize(), 0.6);
        let shift = DVec3::new(30., -40., 70.);
        let mut moved = mesh.clone();
        moved.nodes = mesh
            .nodes
            .iter()
            .map(|(&n, &p)| (100 - n, q * p + shift))
            .collect();
        moved.elements.reverse();
        for e in &mut moved.elements {
            e.id = 100 - e.id;
            for n in &mut e.nodes {
                *n = 100 - *n;
            }
            e.nodes.reverse();
        }
        let transformed = assemble(&moved, &planar_frame(&moved, q * DVec3::Z), &policy).unwrap();
        assert!(transformed.all_surface_patches_built);
        let groups = |r: &Report, reverse: bool| -> BTreeSet<Vec<u32>> {
            r.preview
                .surfaces
                .iter()
                .map(|s| {
                    let mut ids: Vec<_> = s
                        .source_elements
                        .iter()
                        .map(|&n| if reverse { 100 - n } else { n })
                        .collect();
                    ids.sort_unstable();
                    ids
                })
                .collect()
        };
        assert_eq!(groups(&base, false), groups(&transformed, true));
        for (&n, &p) in base.vertex_source_nodes.iter().zip(&base.preview.vertices) {
            let k = transformed
                .vertex_source_nodes
                .iter()
                .position(|&v| v == 100 - n)
                .unwrap();
            assert!(
                DVec3::from_array(transformed.preview.vertices[k])
                    .distance(q * DVec3::from_array(p) + shift)
                    < 1e-7
            );
        }
    }

    fn annulus(triangle: bool) -> MeshData {
        use crate::input::ElementData;
        let (outer, inner) = if triangle {
            (
                vec![[-2., -2., 0.], [2., -2., 0.], [0., 2., 0.]],
                vec![[-0.005, 0., 0.], [0.005, 0., 0.], [0., 0.5, 0.]],
            )
        } else {
            (
                vec![[-2., -2., 0.], [2., -2., 0.], [2., 2., 0.], [-2., 2., 0.]],
                vec![
                    [-0.005, -0.2, 0.],
                    [0.005, -0.2, 0.],
                    [0.005, 0.2, 0.],
                    [-0.005, 0.2, 0.],
                ],
            )
        };
        let n = outer.len() as u32;
        let mut mesh = MeshData::default();
        for (i, p) in outer.into_iter().chain(inner).enumerate() {
            mesh.nodes.insert(i as u32 + 1, DVec3::from_array(p));
        }
        for i in 0..n {
            mesh.elements.push(ElementData {
                id: i + 1,
                elem_type: 44,
                stiff_id: 9,
                nodes: vec![i + 1, (i + 1) % n + 1, (i + 1) % n + n + 1, i + n + 1],
            });
        }
        mesh
    }

    #[test]
    fn holes_restore_across_shapes_scales_rotations_and_renumbering() {
        for triangle in [true, false] {
            for scale in [0.1, 1., 10.] {
                for transformed in [false, true] {
                    let base = annulus(triangle);
                    let count = if triangle { 3 } else { 4 };
                    let rotation = if transformed {
                        glam::DQuat::from_axis_angle(DVec3::new(1., 2., 3.).normalize(), 0.71)
                    } else {
                        glam::DQuat::IDENTITY
                    };
                    let shift = DVec3::new(17., -29., 51.);
                    let number = |n| if transformed { 100 - n } else { n };
                    let mut mesh = base.clone();
                    mesh.nodes = base
                        .nodes
                        .iter()
                        .map(|(&n, &p)| (number(n), rotation * (p * scale) + shift))
                        .collect();
                    for e in &mut mesh.elements {
                        e.id = number(e.id);
                        for n in &mut e.nodes {
                            *n = number(*n);
                        }
                        e.nodes.reverse();
                    }
                    if transformed {
                        mesh.elements.reverse();
                    }
                    let mut f = planar_frame(&mesh, rotation * DVec3::Z);
                    for n in count + 1..=2 * count {
                        let mut p = base.nodes[&n];
                        p.x = 0.;
                        let i = f.node_ids.iter().position(|&i| i == number(n)).unwrap();
                        f.candidate_points[i] = (rotation * (p * scale) + shift).to_array();
                    }
                    let policy = Policy {
                        closure_tolerance: 0.001 * scale,
                        junction_movement_limit: 0.05 * scale,
                        precision: 1e-7 * scale,
                        minimum_edge: 0.001 * scale,
                    };
                    let result = assemble(&mesh, &f, &policy).unwrap();
                    assert!(result.all_surface_patches_built, "{:?}", result.issues);
                    assert_eq!(result.hole_recovery.len(), 1);
                    assert_eq!(result.hole_recovery[0].outcome, HoleOutcome::Restored);
                    assert_eq!(result.preview.surfaces.len(), 1);
                    assert_eq!(result.preview.surfaces[0].contours.len(), 2);
                    assert_eq!(result.surface_stiffness, vec![9]);
                    assert_eq!(
                        result.preview.surfaces[0].source_elements.len(),
                        count as usize
                    );
                    for n in count + 1..=2 * count {
                        let i = result
                            .vertex_source_nodes
                            .iter()
                            .position(|&i| i == number(n))
                            .unwrap();
                        assert!(
                            DVec3::from_array(result.preview.vertices[i])
                                .distance(mesh.nodes[&number(n)])
                                < policy.precision
                        );
                    }
                    assert!(result.maximum_closure_movement <= policy.junction_movement_limit);
                    // Re-running the assembly is deterministic and leaves input intact.
                    let again = assemble(&mesh, &f, &policy).unwrap();
                    assert_eq!(
                        serde_json::to_string(&result).unwrap(),
                        serde_json::to_string(&again).unwrap()
                    );
                }
            }
        }
    }

    #[test]
    fn multiple_holes_in_one_region_restore_together() {
        use crate::input::ElementData;
        let mut mesh = MeshData::default();
        for y in 0..4 {
            for x in 0..6 {
                mesh.nodes.insert(
                    1 + x + 6 * y,
                    DVec3::new(x as f64 * 0.01, y as f64 * 0.01, 0.),
                );
            }
        }
        for y in 0..3 {
            for x in 0..5 {
                if y == 1 && (x == 1 || x == 3) {
                    continue;
                }
                let a = 1 + x + 6 * y;
                mesh.elements.push(ElementData {
                    id: mesh.elements.len() as u32 + 1,
                    elem_type: 44,
                    stiff_id: 9,
                    nodes: vec![a, a + 1, a + 7, a + 6],
                });
            }
        }
        let mut f = planar_frame(&mesh, DVec3::Z);
        for (i, &n) in f.node_ids.iter().enumerate() {
            let p = mesh.nodes[&n];
            if (n >= 8 && n <= 11) || (n >= 14 && n <= 17) {
                f.candidate_points[i][0] = if p.x < 0.025 { 0.015 } else { 0.035 };
            }
        }
        let policy = Policy {
            closure_tolerance: 0.001,
            junction_movement_limit: 0.05,
            precision: 1e-7,
            minimum_edge: 0.001,
        };
        let result = assemble(&mesh, &f, &policy).unwrap();
        assert!(result.all_surface_patches_built, "{:?}", result.issues);
        assert_eq!(result.hole_recovery.len(), 1);
        assert_eq!(result.hole_recovery[0].outcome, HoleOutcome::Restored);
        assert_eq!(result.hole_recovery[0].constraints.len(), 2);
        assert_eq!(result.preview.surfaces[0].contours.len(), 3);
        assert_eq!(result.preview.surfaces[0].source_elements.len(), 13);
    }

    #[test]
    fn blocked_axis_hole_does_not_prevent_independent_recovery() {
        let mut mesh = annulus(true);
        let second = annulus(false);
        mesh.nodes.extend(
            second
                .nodes
                .iter()
                .map(|(&n, &p)| (n + 20, p + DVec3::X * 10.)),
        );
        mesh.elements
            .extend(second.elements.into_iter().map(|mut e| {
                e.id += 20;
                for n in &mut e.nodes {
                    *n += 20;
                }
                e
            }));
        let mut f = planar_frame(&mesh, DVec3::Z);
        for (i, &n) in f.node_ids.iter().enumerate() {
            if (4..=6).contains(&n) {
                f.candidate_points[i][0] = 0.;
            }
            if (25..=28).contains(&n) {
                f.candidate_points[i][0] = 10.;
            }
        }
        let anchor = f.node_ids.iter().position(|&n| n == 4).unwrap();
        let end = f.node_ids.iter().position(|&n| n == 1).unwrap();
        f.axes.push(frame::Axis {
            endpoints: [anchor, end],
            anchors: vec![frame::Anchor {
                node: anchor,
                t: 0.,
            }],
            spans: vec![],
        });
        let policy = Policy {
            closure_tolerance: 0.001,
            junction_movement_limit: 0.05,
            precision: 1e-7,
            minimum_edge: 0.001,
        };
        let result = assemble(&mesh, &f, &policy).unwrap();
        assert_eq!(result.hole_recovery.len(), 2);
        let blocked = result
            .hole_recovery
            .iter()
            .find(|r| r.source_elements.contains(&1))
            .unwrap();
        assert_eq!(blocked.outcome, HoleOutcome::AxisAnchorRequiresJointRepair);
        assert!(blocked.changes.is_empty());
        assert!(result
            .hole_recovery
            .iter()
            .any(|r| r.outcome == HoleOutcome::Restored));
        assert_eq!(result.preview.surfaces.len(), 1);
        assert_eq!(
            result.preview.surfaces[0].source_elements,
            vec![21, 22, 23, 24]
        );
        let index = result
            .vertex_source_nodes
            .iter()
            .position(|&n| n == 4)
            .unwrap();
        assert!(
            DVec3::from_array(result.preview.vertices[index])
                .distance(DVec3::from_array(f.candidate_points[anchor]))
                < policy.precision
        );
    }

    #[test]
    fn common_support_line_is_reported_without_moving_any_vertex() {
        let mesh = annulus(true);
        let f = planar_frame(&mesh, DVec3::Z);
        let policy = Policy {
            closure_tolerance: 0.001,
            junction_movement_limit: 0.05,
            precision: 1e-7,
            minimum_edge: 0.001,
        };
        let regions = vec![(0, 9, vec![1, 2, 3])];
        let rings = BTreeMap::from([(0, vec![vec![1, 2, 3], vec![4, 5, 6]])]);
        let supports = vec![
            PlaneFrame::new([0.; 3], [0., 0., 1.]).unwrap(),
            PlaneFrame::new([0.; 3], [1., 0., 0.]).unwrap(),
        ];
        let owners = (1..=6)
            .map(|n| (n, if n >= 4 { vec![0, 1] } else { vec![0] }))
            .collect();
        let mut points: BTreeMap<_, _> = mesh
            .nodes
            .iter()
            .map(|(&n, &p)| (n, if n >= 4 { DVec3::new(0., p.y, 0.) } else { p }))
            .collect();
        let before = points.clone();
        let reports = holes::recover(
            &mut points,
            &holes::Context {
                mesh: &mesh,
                source: &f,
                policy: &policy,
                regions: &regions,
                rings: &rings,
                owners: &owners,
                representatives: &[0, 1],
                supports: &supports,
            },
        );
        assert_eq!(points, before);
        assert_eq!(reports.len(), 1);
        assert_eq!(
            reports[0].outcome,
            HoleOutcome::CommonSupportsForceLineOrPoint
        );
        assert_eq!(reports[0].constraints[0].normal_rank, 2);
        assert_eq!(reports[0].constraints[0].common_supports, vec![0, 1]);
    }

    #[test]
    fn hole_recovery_cannot_degenerate_a_neighbor_and_rolls_back_all_vertices() {
        let mut mesh = annulus(true);
        mesh.nodes.insert(9, DVec3::new(0.1, -0.2, 0.));
        let f = planar_frame(&mesh, DVec3::Z);
        let policy = Policy {
            closure_tolerance: 0.001,
            junction_movement_limit: 0.05,
            precision: 1e-7,
            minimum_edge: 0.001,
        };
        let regions = vec![(0, 9, vec![1, 2, 3]), (0, 9, vec![4])];
        let rings = BTreeMap::from([
            (0, vec![vec![1, 2, 3], vec![4, 5, 6]]),
            (1, vec![vec![4, 1, 9]]),
        ]);
        let supports = vec![PlaneFrame::new([0.; 3], [0., 0., 1.]).unwrap()];
        let owners = mesh.nodes.keys().map(|&n| (n, vec![0])).collect();
        let mut points: BTreeMap<_, _> = mesh
            .nodes
            .iter()
            .map(|(&n, &p)| {
                (
                    n,
                    if (4..=6).contains(&n) {
                        DVec3::new(0., p.y, 0.)
                    } else {
                        p
                    },
                )
            })
            .collect();
        // Neighbor is valid now; restoring vertex 4 would coincide with its
        // other vertex 9. A local hole-only validator would miss this.
        points.insert(9, mesh.nodes[&4]);
        let before = points.clone();
        let reports = holes::recover(
            &mut points,
            &holes::Context {
                mesh: &mesh,
                source: &f,
                policy: &policy,
                regions: &regions,
                rings: &rings,
                owners: &owners,
                representatives: &[0],
                supports: &supports,
            },
        );
        assert_eq!(reports.len(), 1);
        assert_eq!(reports[0].outcome, HoleOutcome::ContourConflict);
        assert!(reports[0].changes.is_empty());
        assert_eq!(points, before);
    }

    #[test]
    fn collapsed_hole_is_reported_with_source_nodes_and_never_filled() {
        use crate::input::ElementData;
        let mut mesh = MeshData::default();
        for (i, p) in [
            [-2., -2., 0.],
            [2., -2., 0.],
            [0., 2., 0.],
            [-0.005, 0., 0.],
            [0.005, 0., 0.],
            [0., 0.5, 0.],
        ]
        .into_iter()
        .enumerate()
        {
            mesh.nodes.insert(i as u32 + 1, DVec3::from_array(p));
        }
        for i in 0..3_u32 {
            mesh.elements.push(ElementData {
                id: i + 1,
                elem_type: 44,
                stiff_id: 9,
                nodes: vec![i + 1, (i + 1) % 3 + 1, (i + 1) % 3 + 4, i + 4],
            });
        }
        let policy = Policy {
            closure_tolerance: 0.001,
            junction_movement_limit: 0.05,
            precision: 1e-7,
            minimum_edge: 0.001,
        };
        let mut f = planar_frame(&mesh, DVec3::Z);
        let base = assemble(&mesh, &f, &policy).unwrap();
        assert!(base.all_surface_patches_built, "{:?}", base.issues);
        assert_eq!(base.preview.surfaces[0].contours.len(), 2);
        // A regularized candidate flattens the narrow hole into a line with
        // distinct points, so this is not merely a duplicate-node test.
        for (i, &n) in f.node_ids.iter().enumerate() {
            if n == 4 {
                f.candidate_points[i] = [0., 0., 0.];
            }
            if n == 5 {
                f.candidate_points[i] = [0., 0.1, 0.];
            }
        }
        let result = assemble(&mesh, &f, &policy).unwrap();
        assert!(!result.all_surface_patches_built);
        assert!(result.preview.surfaces.is_empty());
        assert!(result.preview.edges.is_empty());
        assert_eq!(result.issues.len(), 1);
        assert_eq!(result.issues[0].reason, "degenerate_hole_after_closure");
        assert_eq!(result.hole_recovery[0].outcome, HoleOutcome::MovementBudget);
        assert_eq!(result.issues[0].source_elements, vec![1, 2, 3]);
        assert!(result.issues[0].boundary_source_nodes.iter().any(|r| r
            .iter()
            .copied()
            .collect::<BTreeSet<_>>()
            == BTreeSet::from([4, 5, 6])));
    }

    #[test]
    fn nonmanifold_source_is_not_hidden_by_pinch_partition() {
        let facets = BTreeMap::from([(1, vec![1, 2, 3]), (2, vec![2, 1, 4]), (3, vec![1, 2, 5])]);
        assert_eq!(
            split_pinched_regions(&[1, 2, 3], &facets),
            vec![vec![1, 2, 3]]
        );
    }

    #[test]
    fn offset_projection_closes_three_walls_without_changing_normals() {
        let planes = vec![
            PlaneFrame::new([0.; 3], [1., 0., 0.]).unwrap(),
            PlaneFrame::new([0.; 3], [0., 1., 0.]).unwrap(),
            PlaneFrame::new([0.0003, 0., 0.], [1., 1., 0.]).unwrap(),
        ];
        assert!(intersection(DVec3::ZERO, &planes.iter().collect::<Vec<_>>(), 1e-9).is_none());
        let closed = concurrent_supports(&planes, &[vec![0, 1, 2], vec![2, 1, 0]]);
        let point = intersection(DVec3::ZERO, &closed.iter().collect::<Vec<_>>(), 1e-9).unwrap();
        assert!(point.length() < 0.001);
        for (a, b) in planes.iter().zip(&closed) {
            assert_eq!(a.normal, b.normal);
            assert!(a.distance(b.origin).abs() < 0.001);
        }
        let rotation = glam::DQuat::from_axis_angle(DVec3::new(1., 2., 3.).normalize(), 0.6);
        let shift = DVec3::new(10., 20., 30.);
        let moved: Vec<_> = planes
            .iter()
            .map(|p| {
                PlaneFrame::new(
                    (rotation * DVec3::from_array(p.origin) + shift).to_array(),
                    (rotation * DVec3::from_array(p.normal)).to_array(),
                )
                .unwrap()
            })
            .collect();
        let result = concurrent_supports(&moved, &[vec![2, 0, 1]]);
        let transformed = intersection(shift, &result.iter().collect::<Vec<_>>(), 1e-9).unwrap();
        assert!(transformed.distance(rotation * point + shift) < 1e-9);
    }

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
            junction_movement_limit: 0.001,
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
        let mut planar = mesh.clone();
        planar.nodes.insert(5, DVec3::new(-1., 1., 0.));
        planar.nodes.insert(6, DVec3::new(-1., 0., 0.));
        let planar_planes = planes::recognize(
            &planar,
            &planes::Policy {
                angle: 0.02,
                distance: 0.01,
                precision: 1e-8,
            },
        )
        .unwrap();
        assert_eq!(planar_planes.patches.len(), 1);
        let planar_frame = frame::solve(&planar, &axes, &planar_planes, &f.policy).unwrap();
        let materialized = assemble(&planar, &planar_frame, &p).unwrap();
        assert!(materialized.all_surface_patches_built);
        assert_eq!(materialized.surface_source_patches, vec![0, 0]);
        assert_eq!(materialized.surface_stiffness, vec![10, 20]);
        assert_eq!(materialized.preview.surfaces.len(), 2);
        assert_eq!(materialized.preview.edges.len(), 7);
        assert_eq!(materialized.preview.surfaces[0].source_elements, vec![1]);
        assert_eq!(materialized.preview.surfaces[1].source_elements, vec![2]);
        f.candidate_points[0] = [0.002, 0., 0.002];
        let blocked = assemble(&mesh, &f, &p).unwrap();
        assert!(!blocked.all_surface_patches_built);
        assert!(blocked.preview.surfaces.is_empty());
        assert_eq!(blocked.issues.len(), 2);
        let mut larger = p.clone();
        larger.junction_movement_limit = 0.01;
        let permitted = assemble(&mesh, &f, &larger).unwrap();
        assert!(permitted.all_surface_patches_built);
        assert!(permitted.maximum_closure_movement > larger.closure_tolerance);
        assert_eq!(permitted.preview.edges.len(), 7);
        let mut limited = f.clone();
        limited.policy.maximum_movement = 1e-6;
        // Move the shared point tangentially so projection does not remove its
        // accumulated displacement from the immutable source model.
        limited.candidate_points[0] = [0.002, 0.0002, 0.002];
        assert!(
            !assemble(&mesh, &limited, &larger)
                .unwrap()
                .all_surface_patches_built
        );
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
