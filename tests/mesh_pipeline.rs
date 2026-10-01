#[path = "../examples/support/mesh_fixture.rs"]
mod fixture;
use glam::{DQuat, DVec3};
use std::collections::{BTreeMap, BTreeSet};
use topo_reconstruct_rs::reconstruction::{assembly, frame, mesh, planes, recognize};

fn run(scale: f64, rotated: bool) -> (assembly::Report, mesh::Report) {
    run_input(fixture::source(), scale, rotated)
}

fn run_input(
    input: topo_reconstruct_rs::input::MeshData,
    scale: f64,
    rotated: bool,
) -> (assembly::Report, mesh::Report) {
    run_with(input, scale, rotated, None)
}

fn run_with(
    mut input: topo_reconstruct_rs::input::MeshData,
    scale: f64,
    rotated: bool,
    features: Option<assembly::FeaturePolicy>,
) -> (assembly::Report, mesh::Report) {
    let rotation = if rotated {
        DQuat::from_axis_angle(DVec3::new(1., 2., 3.).normalize(), 0.6)
    } else {
        DQuat::IDENTITY
    };
    let shift = if rotated {
        DVec3::new(10., 20., 30.)
    } else {
        DVec3::ZERO
    };
    input.nodes = input
        .nodes
        .into_iter()
        .map(|(id, p)| {
            (
                if rotated { 1000 - id } else { id },
                rotation * (p * scale) + shift,
            )
        })
        .collect();
    if rotated {
        input.elements.reverse();
        for e in &mut input.elements {
            for id in &mut e.nodes {
                *id = 1000 - *id;
            }
        }
    }
    let a = recognize::recognize(
        &input,
        &recognize::Policy {
            angle: 0.02,
            line_tolerance: 0.001 * scale,
            numerical_precision: 1e-8 * scale,
        },
    )
    .unwrap();
    let p = planes::recognize(
        &input,
        &planes::Policy {
            angle: 0.02,
            distance: 0.001 * scale,
            precision: 1e-8 * scale,
        },
    )
    .unwrap();
    let f = frame::solve(
        &input,
        &a,
        &p,
        &frame::Policy {
            up: (rotation * DVec3::Z).to_array(),
            angle: 0.02,
            maximum_movement: 0.05 * scale,
            relative_movement: 0.05,
            minimum_length: 0.01 * scale,
            residual_tolerance: 1e-8 * scale,
            iterations: 1000,
            panel_tolerance: 0.,
            geotechnical: false,
            over_constrained_panels: false,
        },
    )
    .unwrap();
    let policy = assembly::Policy {
        closure_tolerance: 0.001 * scale,
        junction_movement_limit: 0.05 * scale,
        precision: 1e-7 * scale,
        minimum_edge: 0.001 * scale,
    };
    let t = match features {
        None => assembly::assemble(&input, &f, &policy),
        Some(features) => assembly::assemble_geotechnical(&input, &f, &policy, &features),
    }
    .unwrap();
    let m = mesh::build(
        &t,
        &mesh::Policy {
            boundary_spacing: 0.5 * scale,
            maximum_area: 0.5 * scale * scale,
            minimum_angle_degrees: 20.,
            maximum_added_vertices_per_surface: 10000,
        },
    )
    .unwrap();
    (t, m)
}

#[test]
fn fe_to_mesh_preserves_opening_properties_and_shared_joints() {
    for scale in [0.1, 1., 10.] {
        for rotated in [false, true] {
            let (topology, mesh) = run(scale, rotated);
            assert!(
                mesh.topology_valid && mesh.quality_passed,
                "scale={scale} rotated={rotated}: {:?}, angle {}",
                mesh.blockers,
                mesh.minimum_angle_degrees
            );
            assert!(mesh.maximum_edge_ratio.is_finite() && mesh.maximum_edge_ratio >= 1.0);
            assert_eq!(topology.preview.surfaces().len(), 4);
            assert_eq!(topology.axis_assembly.axes.len(), 3);
            let interval_contacts = topology
                .axis_assembly
                .contacts
                .iter()
                .filter(|contact| matches!(contact, assembly::bars::Contact::Interval { .. }))
                .count();
            assert!(mesh.constraint_synchronization.shared_edge_count > 0);
            assert_eq!(
                mesh.constraint_synchronization.interval_contact_count,
                interval_contacts
            );
            assert_eq!(
                mesh.constraint_synchronization
                    .synchronized_interval_endpoints,
                interval_contacts * 2
            );
            assert_eq!(
                topology
                    .preview
                    .surfaces()
                    .iter()
                    .filter(|s| s.contours.len() == 2)
                    .count(),
                1
            );
            assert!(
                mesh.external_mesher_ready,
                "{:?}",
                mesh.external_mesher_blockers
            );
            assert!(mesh.external_mesher_blockers.is_empty());
            assert_eq!(
                mesh.triangles
                    .iter()
                    .map(|t| t.stiffness)
                    .collect::<BTreeSet<_>>(),
                BTreeSet::from([10, 20, 30, 40])
            );
            assert_eq!(
                mesh.bars
                    .iter()
                    .map(|b| b.stiffness)
                    .collect::<BTreeSet<_>>(),
                BTreeSet::from([50, 60, 70])
            );
            let mut areas = BTreeMap::<u32, f64>::new();
            for t in &mesh.triangles {
                let [a, b, c] = t.vertices.map(|n| DVec3::from_array(mesh.vertices[n]));
                *areas.entry(t.stiffness).or_default() += (b - a).cross(c - a).length() / 2.;
            }
            for (stiffness, expected) in [(10, 12.), (20, 11.), (30, 12.), (40, 8.)] {
                assert!((areas[&stiffness] / scale.powi(2) - expected).abs() < 1e-7);
            }
            let mut lengths = BTreeMap::<u32, f64>::new();
            for b in &mesh.bars {
                *lengths.entry(b.stiffness).or_default() +=
                    DVec3::from_array(mesh.vertices[b.vertices[0]])
                        .distance(DVec3::from_array(mesh.vertices[b.vertices[1]]));
            }
            for (stiffness, expected) in [(50, 3.), (60, 3.), (70, 2.)] {
                assert!((lengths[&stiffness] / scale - expected).abs() < 1e-7);
            }
            // The beam/column common node is also a triangle vertex, by identity.
            let a = &topology.axis_assembly.axes;
            let common: BTreeSet<_> = a[0]
                .anchors
                .iter()
                .map(|a| a.vertex)
                .collect::<BTreeSet<_>>()
                .intersection(&a[1].anchors.iter().map(|a| a.vertex).collect())
                .copied()
                .collect();
            assert_eq!(common.len(), 1);
            let node = *common.first().unwrap();
            assert!(mesh.triangles.iter().any(|t| t.vertices.contains(&node)));
            assert!(mesh
                .bars
                .iter()
                .any(|b| b.axis == 0 && b.vertices.contains(&node)));
            assert!(mesh
                .bars
                .iter()
                .any(|b| b.axis == 1 && b.vertices.contains(&node)));
            assert!(!mesh.export_ready);
        }
    }
}

#[test]
fn rejects_incomplete_geometry_and_reports_quality_failure() {
    let (mut topology, _) = run(1., false);
    let mut policy = mesh::Policy {
        boundary_spacing: 1.,
        maximum_area: 0.001,
        minimum_angle_degrees: 20.,
        maximum_added_vertices_per_surface: 1,
    };
    let m = mesh::build(&topology, &policy).unwrap();
    assert!(!m.quality_passed);
    assert!(m.external_mesher_ready, "{:?}", m.external_mesher_blockers);
    assert!(m.external_mesher_blockers.is_empty());
    assert!(!m.quality_diagnostics.is_empty());
    assert!(m
        .quality_diagnostics
        .iter()
        .all(|diagnostic| !diagnostic.source_elements.is_empty()));
    topology.all_surface_patches_built = false;
    policy.maximum_area = 0.5;
    assert!(mesh::build(&topology, &policy).is_err());
}

#[test]
fn partial_mesh_keeps_unresolved_source_coverage_as_a_blocker() {
    let (mut topology, _) = run(1., false);
    topology.all_surface_patches_built = false;
    topology.issues.push(assembly::Issue {
        patch: 99,
        source_elements: vec![9001],
        reason: "synthetic unresolved surface".into(),
        boundary_source_nodes: vec![],
    });
    topology.axis_assembly.all_axes_built = false;
    topology.axis_assembly.issues.push(assembly::bars::Issue {
        source_axis: 99,
        source_elements: vec![9002],
        source_nodes: vec![1],
        reason: "synthetic unresolved axis".into(),
    });
    let policy = mesh::Policy {
        boundary_spacing: 0.5,
        maximum_area: 0.5,
        minimum_angle_degrees: 20.,
        maximum_added_vertices_per_surface: 10000,
    };
    assert!(mesh::build(&topology, &policy).is_err());
    let partial = mesh::build_partial(&topology, &policy).unwrap();
    assert!(!partial.source_coverage_complete);
    assert_eq!(partial.unresolved_surface_source_elements, vec![9001]);
    assert_eq!(partial.unresolved_axis_source_elements, vec![9002]);
    assert!(!partial.triangles.is_empty());
    assert!(!partial.external_mesher_ready);
    assert!(partial
        .external_mesher_blockers
        .contains(&"unresolved_surface_assembly".to_string()));
    assert!(partial
        .external_mesher_blockers
        .contains(&"unresolved_axis_assembly".to_string()));
}

#[test]
fn subresolution_hole_is_reported_without_silent_filling() {
    let mut model = topo_reconstruct_rs::reconstruction::Model::new(1e-7, 0.001).unwrap();
    let plane = model.add_plane(
        topo_reconstruct_rs::reconstruction::PlaneFrame::new([0., 0., 0.], [0., 0., 1.]).unwrap(),
    );
    let points = [
        [0., 0., 0.],
        [10., 0., 0.],
        [10., 10., 0.],
        [0., 10., 0.],
        [4., 4., 0.],
        [4.4, 4., 0.],
        [4.2, 4.00000001, 0.],
    ];
    let vertices: Vec<_> = points
        .into_iter()
        .map(|point| model.add_vertex(point).unwrap())
        .collect();
    model
        .add_surface(
            plane,
            vec![vertices[..4].to_vec(), vertices[4..].to_vec()],
            vec![1],
        )
        .unwrap();
    let topology = assembly::Report {
        policy: assembly::Policy {
            closure_tolerance: 0.001,
            junction_movement_limit: 0.05,
            precision: 1e-7,
            minimum_edge: 0.001,
        },
        export_ready: false,
        all_surface_patches_built: true,
        preview: model,
        vertex_source_nodes: (1..=7).collect(),
        surface_source_patches: vec![0],
        surface_stiffness: vec![10],
        pinched_region_splits: vec![],
        hole_recovery: vec![],
        feature_policy: None,
        simplified_holes: vec![],
        axis_assembly: assembly::bars::Report::default(),
        junctions: assembly::junctions::Report::default(),
        consoles: assembly::consoles::Report::default(),
        stacked_walls: assembly::stacking::Report::default(),
        straightened_edges: vec![],
        short_edges: assembly::cleanup::Report::default(),
        coincident_vertices: assembly::cleanup::MergeReport::default(),
        wall_ends: assembly::cleanup::MergeReport::default(),
        bar_ends: assembly::cleanup::MergeReport::default(),
        bar_anchors: assembly::cleanup::MergeReport::default(),
        bar_tees: Default::default(),
        cracks: vec![],
        removed_slivers: vec![],
        short_edge_merges: assembly::cleanup::MergeReport::default(),
        gaps: assembly::gaps::Report::default(),
        short_bars: assembly::cleanup::BarCollapseReport::default(),
        issues: vec![],
        maximum_closure_movement: 0.,
        rejected_vertices: BTreeMap::new(),
        support_representatives: vec![0],
        support_offset_projection_applied: false,
    };
    let report = mesh::build_partial(
        &topology,
        &mesh::Policy {
            boundary_spacing: 0.5,
            maximum_area: 0.5,
            minimum_angle_degrees: 20.,
            maximum_added_vertices_per_surface: 100,
        },
    )
    .unwrap();
    assert!(report.triangles.is_empty());
    assert_eq!(report.mesh_surface_errors.len(), 1);
    assert_eq!(report.mesh_surface_errors[0].surface, 0);
    assert!(report.mesh_surface_errors[0]
        .reason
        .starts_with("hole_area_below_mesh_resolution"));
    assert_eq!(report.unresolved_surface_source_elements, vec![1]);
    assert!(!report.source_coverage_complete);
    assert!(!report.external_mesher_ready);
}

#[test]
fn nonorthogonal_fragment_passes_quality_under_transforms() {
    for scale in [0.1, 1., 10.] {
        for rotated in [false, true] {
            let mut input = fixture::source();
            for p in input.nodes.values_mut() {
                p.x += 0.2 * p.y;
            }
            let (_topology, mesh) = run_input(input, scale, rotated);
            assert!(
                mesh.topology_valid && mesh.quality_passed,
                "scale={scale}, rotated={rotated}, angle={}, {:?}",
                mesh.minimum_angle_degrees,
                mesh.blockers
            );
            assert!(mesh.blockers.is_empty());
            assert!(mesh.maximum_edge_ratio.is_finite() && mesh.maximum_edge_ratio >= 1.0);
            assert!(!mesh.export_ready);
        }
    }
}

#[test]
fn cantilever_keeps_topology_and_refines_around_open_constraint() {
    let mut input = fixture::source();
    input.elements.retain(|e| {
        e.elem_type != 10
            || !e.nodes.iter().all(|n| input.nodes[n].z == 0.)
            || e.nodes.iter().all(|n| input.nodes[n].x <= 2.)
    });
    let (topology, mesh) = run_input(input, 1., true);
    assert!(mesh.topology_valid, "{:?}", mesh.blockers);
    assert!(mesh.quality_passed, "{:?}", mesh.blockers);
    assert!(mesh.maximum_edge_ratio.is_finite() && mesh.maximum_edge_ratio >= 1.0);
    assert!(!mesh.export_ready);
    assert_eq!(topology.axis_assembly.axes.len(), 2);
    assert_eq!(
        mesh.bars
            .iter()
            .map(|b| b.stiffness)
            .collect::<BTreeSet<_>>(),
        BTreeSet::from([50, 70])
    );
    let mut triangle_edges = BTreeMap::<[usize; 2], usize>::new();
    for triangle in &mesh.triangles {
        for i in 0..3 {
            let edge = [triangle.vertices[i], triangle.vertices[(i + 1) % 3]];
            let edge = [edge[0].min(edge[1]), edge[0].max(edge[1])];
            *triangle_edges.entry(edge).or_default() += 1;
        }
    }
    // The free endpoint is internal to the slab, but the beam remains a
    // conforming two-sided constraint rather than a disconnected overlay.
    for bar in &mesh.bars {
        if bar.axis != 0 {
            continue;
        }
        let edge = [
            bar.vertices[0].min(bar.vertices[1]),
            bar.vertices[0].max(bar.vertices[1]),
        ];
        assert_eq!(triangle_edges.get(&edge), Some(&2));
    }
}

/// Slab with a wall standing on an interior line (T-junction with both wall
/// ends inside the slab) and a second wall passing through the slab (interior
/// crossing). Source FE nodes are shared along both lines.
fn junction_source() -> topo_reconstruct_rs::input::MeshData {
    use topo_reconstruct_rs::input::{ElementData, MeshData};
    let mut mesh = MeshData::default();
    let mut nodes = BTreeMap::new();
    let mut add = |coordinates: [(i32, i32, i32); 4], stiffness| {
        let ids = coordinates
            .into_iter()
            .map(|p| {
                *nodes.entry(p).or_insert_with(|| {
                    let id = mesh.nodes.len() as u32 + 1;
                    mesh.nodes
                        .insert(id, DVec3::new(p.0 as f64, p.1 as f64, p.2 as f64));
                    id
                })
            })
            .collect();
        mesh.elements.push(ElementData {
            id: mesh.elements.len() as u32 + 1,
            elem_type: 44,
            stiff_id: stiffness,
            nodes: ids,
        });
    };
    for x in 0..6 {
        for y in 0..4 {
            add(
                [(x, y, 0), (x + 1, y, 0), (x + 1, y + 1, 0), (x, y + 1, 0)],
                10,
            );
        }
    }
    for x in 1..3 {
        for z in 0..2 {
            add(
                [(x, 2, z), (x + 1, 2, z), (x + 1, 2, z + 1), (x, 2, z + 1)],
                20,
            );
        }
    }
    for y in 1..3 {
        for z in -1..1 {
            add(
                [(5, y, z), (5, y + 1, z), (5, y + 1, z + 1), (5, y, z + 1)],
                30,
            );
        }
    }
    mesh
}

/// Every model edge of a junction is represented by mesh edges used by the
/// triangles of both surfaces: one shared subdivision, no overlaid meshes.
fn assert_mesh_conforming(topology: &assembly::Report, mesh: &mesh::Report) {
    let model = &topology.preview;
    let vertices = &mesh.vertices;
    let mut edges = vec![BTreeSet::new(); model.surfaces().len()];
    for t in &mesh.triangles {
        for i in 0..3 {
            let (a, b) = (t.vertices[i], t.vertices[(i + 1) % 3]);
            edges[t.surface].insert([a.min(b), a.max(b)]);
        }
    }
    let point = |v: usize| DVec3::from_array(vertices[v]);
    let eps = topology.policy.precision;
    for junction in &topology.junctions.junctions {
        let [s, r] = junction.surfaces;
        for &e in &junction.edges {
            let [a, b] = model.edges()[e];
            let (pa, pb) = (point(a), point(b));
            let d = pb - pa;
            let mut on: Vec<(f64, usize)> = edges[s]
                .iter()
                .flatten()
                .copied()
                .collect::<BTreeSet<_>>()
                .into_iter()
                .filter_map(|v| {
                    let t = (point(v) - pa).dot(d) / d.length_squared();
                    ((-1e-12..=1. + 1e-12).contains(&t) && point(v).distance(pa + d * t) <= eps)
                        .then_some((t, v))
                })
                .collect();
            on.sort_by(|x, y| x.0.total_cmp(&y.0));
            assert!(on.len() >= 2 && on[0].1 == a && on[on.len() - 1].1 == b);
            for w in on.windows(2) {
                let key = [w[0].1.min(w[1].1), w[0].1.max(w[1].1)];
                assert!(
                    edges[s].contains(&key) && edges[r].contains(&key),
                    "junction {:?} edge {e} not shared by both meshes",
                    junction.surfaces
                );
            }
        }
    }
}

#[test]
fn surface_junctions_are_shared_topology_and_conforming_mesh() {
    for scale in [0.1, 1., 10.] {
        for rotated in [false, true] {
            let (topology, mesh) = run_input(junction_source(), scale, rotated);
            let junctions = &topology.junctions;
            assert!(junctions.issues.is_empty(), "{:?}", junctions.issues);
            let kinds: BTreeSet<_> = junctions
                .junctions
                .iter()
                .map(|j| format!("{:?}", j.kind))
                .collect();
            assert_eq!(
                kinds,
                BTreeSet::from(["Crossing".to_string(), "TJunction".to_string()]),
                "scale={scale} rotated={rotated}"
            );
            // The crossing wall contour already carries the shared FE nodes at
            // slab level, so no vertex is generated and none is moved.
            assert!(junctions.generated_vertices.is_empty());
            assert!(junctions.snapped_vertices.is_empty());
            assert_eq!(
                topology.preview.vertices().len(),
                topology.vertex_source_nodes.len()
            );
            assert!(
                mesh.topology_valid && mesh.quality_passed,
                "scale={scale} rotated={rotated}: {:?} angle {}",
                mesh.blockers,
                mesh.minimum_angle_degrees
            );
            assert!(
                mesh.external_mesher_ready,
                "{:?}",
                mesh.external_mesher_blockers
            );
            assert_mesh_conforming(&topology, &mesh);
            // All three property regions keep their source provenance.
            assert_eq!(
                mesh.triangles
                    .iter()
                    .map(|t| t.stiffness)
                    .collect::<BTreeSet<_>>(),
                BTreeSet::from([10, 20, 30])
            );
            assert!(!mesh.export_ready);
        }
    }
}

/// Slab whose edge lies 0.15 beyond the axis of the wall carrying it, as in a
/// mid-surface model with the slab edge at the outer wall face. Coordinates
/// are in 5 cm units.
fn console_source() -> topo_reconstruct_rs::input::MeshData {
    use topo_reconstruct_rs::input::{ElementData, MeshData};
    let mut mesh = MeshData::default();
    let mut nodes = BTreeMap::new();
    let mut add = |coordinates: [(i32, i32, i32); 4], stiffness| {
        let ids = coordinates
            .into_iter()
            .map(|p| {
                *nodes.entry(p).or_insert_with(|| {
                    let id = mesh.nodes.len() as u32 + 1;
                    mesh.nodes
                        .insert(id, DVec3::new(p.0 as f64, p.1 as f64, p.2 as f64) * 0.05);
                    id
                })
            })
            .collect();
        mesh.elements.push(ElementData {
            id: mesh.elements.len() as u32 + 1,
            elem_type: 44,
            stiff_id: stiffness,
            nodes: ids,
        });
    };
    let xs = [0, 20, 40, 60, 80, 100, 120];
    let ys = [0, 20, 40, 60, 77, 80];
    for x in xs.windows(2) {
        for y in ys.windows(2) {
            add(
                [
                    (x[0], y[0], 0),
                    (x[1], y[0], 0),
                    (x[1], y[1], 0),
                    (x[0], y[1], 0),
                ],
                10,
            );
        }
        for z in [-40, -20] {
            add(
                [
                    (x[0], 77, z),
                    (x[1], 77, z),
                    (x[1], 77, z + 20),
                    (x[0], 77, z + 20),
                ],
                20,
            );
        }
    }
    mesh
}

#[test]
fn geotechnical_assembly_trims_slab_console_to_wall_axis() {
    for scale in [0.1, 1., 10.] {
        for rotated in [false, true] {
            let features = assembly::FeaturePolicy {
                maximum_console_width: 0.25 * scale,
                maximum_gap: 0.05 * scale,
                maximum_collapsed_edge: 0.05 * scale,
                ..Default::default()
            };
            let (topology, mesh) = run_with(console_source(), scale, rotated, Some(features));
            let trimmed = &topology.consoles.trimmed;
            assert_eq!(trimmed.len(), 1, "{:?}", topology.consoles.kept);
            assert!((trimmed[0].width - 0.15 * scale).abs() < 1e-6 * scale);
            assert!((trimmed[0].area - 0.9 * scale * scale).abs() < 1e-6 * scale * scale);
            assert!(topology.junctions.issues.is_empty());
            assert!(
                mesh.topology_valid && mesh.quality_passed,
                "scale={scale} rotated={rotated}: {:?}",
                mesh.blockers
            );
            assert!(mesh.source_coverage_complete && mesh.external_mesher_ready);
            // Slab and wall share the new slab edge in the mesh.
            let mut edges = [BTreeSet::new(), BTreeSet::new()];
            for t in &mesh.triangles {
                let k = usize::from(t.stiffness == 20);
                for i in 0..3 {
                    let (a, b) = (t.vertices[i], t.vertices[(i + 1) % 3]);
                    edges[k].insert([a.min(b), a.max(b)]);
                }
            }
            let common = edges[0].intersection(&edges[1]).count();
            assert!(common >= 12, "shared mesh edges: {common}");
            // The slab keeps all of its source elements, including the console.
            let slab = mesh
                .surface_source_elements
                .iter()
                .find(|s| s.len() == 30)
                .expect("slab provenance");
            assert_eq!(slab.len(), 30);
        }
    }
}

/// A wall below and a wall above one slab, axes 25 mm apart (aligned outer
/// faces, different thicknesses). The upper wall rests on the slab without
/// shared source nodes. Coordinates are in 25 mm units.
fn stacked_source() -> topo_reconstruct_rs::input::MeshData {
    use topo_reconstruct_rs::input::{ElementData, MeshData};
    let mut mesh = MeshData::default();
    let mut nodes = BTreeMap::new();
    let mut add = |coordinates: [(i32, i32, i32); 4], stiffness| {
        let ids = coordinates
            .into_iter()
            .map(|p| {
                *nodes.entry(p).or_insert_with(|| {
                    let id = mesh.nodes.len() as u32 + 1;
                    mesh.nodes
                        .insert(id, DVec3::new(p.0 as f64, p.1 as f64, p.2 as f64) * 0.025);
                    id
                })
            })
            .collect();
        mesh.elements.push(ElementData {
            id: mesh.elements.len() as u32 + 1,
            elem_type: 44,
            stiff_id: stiffness,
            nodes: ids,
        });
    };
    for x in (0..240).step_by(40) {
        for y in (0..160).step_by(40) {
            add(
                [
                    (x, y, 0),
                    (x + 40, y, 0),
                    (x + 40, y + 40, 0),
                    (x, y + 40, 0),
                ],
                10,
            );
        }
        for z in [-80, -40] {
            add(
                [
                    (x, 80, z),
                    (x + 40, 80, z),
                    (x + 40, 80, z + 40),
                    (x, 80, z + 40),
                ],
                20,
            );
        }
        for z in [0, 40] {
            add(
                [
                    (x, 81, z),
                    (x + 40, 81, z),
                    (x + 40, 81, z + 40),
                    (x, 81, z + 40),
                ],
                30,
            );
        }
    }
    mesh
}

/// Two walls on one slab in one line, the second one 25 mm off the first
/// line (a jog). `lap` (25 mm units) is how far the second wall starts
/// before the end of the first: 0 end to end, 1 a 25 mm strip where the ends
/// overlap, 40 side by side over 1 m.
fn wall_line_source(lap: i32) -> topo_reconstruct_rs::input::MeshData {
    use topo_reconstruct_rs::input::{ElementData, MeshData};
    let mut mesh = MeshData::default();
    let mut nodes = BTreeMap::new();
    let mut add = |coordinates: [(i32, i32, i32); 4], stiffness| {
        let ids = coordinates
            .into_iter()
            .map(|p| {
                *nodes.entry(p).or_insert_with(|| {
                    let id = mesh.nodes.len() as u32 + 1;
                    mesh.nodes
                        .insert(id, DVec3::new(p.0 as f64, p.1 as f64, p.2 as f64) * 0.025);
                    id
                })
            })
            .collect();
        mesh.elements.push(ElementData {
            id: mesh.elements.len() as u32 + 1,
            elem_type: 44,
            stiff_id: stiffness,
            nodes: ids,
        });
    };
    for x in (0..240).step_by(40) {
        for y in (0..160).step_by(40) {
            add(
                [
                    (x, y, 0),
                    (x + 40, y, 0),
                    (x + 40, y + 40, 0),
                    (x, y + 40, 0),
                ],
                10,
            );
        }
    }
    let b: Vec<i32> = std::iter::once(120 - lap)
        .chain((160..=240).step_by(40))
        .collect();
    for z in [0, 40] {
        for x in (0..120).step_by(40) {
            add(
                [
                    (x, 80, z),
                    (x + 40, 80, z),
                    (x + 40, 80, z + 40),
                    (x, 80, z + 40),
                ],
                20,
            );
        }
        for w in b.windows(2) {
            add(
                [
                    (w[0], 81, z),
                    (w[1], 81, z),
                    (w[1], 81, z + 40),
                    (w[0], 81, z + 40),
                ],
                30,
            );
        }
    }
    mesh
}

#[test]
fn walls_in_one_line_adopt_one_plane_but_side_by_side_walls_do_not() {
    for scale in [0.1, 1., 10.] {
        for rotated in [false, true] {
            let features = assembly::FeaturePolicy {
                maximum_console_width: 0.25 * scale,
                maximum_stack_offset: 0.05 * scale,
                maximum_gap: 0.05 * scale,
                maximum_collapsed_edge: 0.05 * scale,
                ..Default::default()
            };
            let (topology, mesh) =
                run_with(wall_line_source(0), scale, rotated, Some(features.clone()));
            let aligned = &topology.stacked_walls.aligned;
            assert_eq!(aligned.len(), 1, "{:?}", topology.stacked_walls.kept);
            assert_eq!(aligned[0].reason, "aligned_wall_line");
            assert!(aligned[0].slab.is_none());
            assert!((aligned[0].offset - 0.025 * scale).abs() < 1e-6 * scale);
            assert!(topology.issues.is_empty() && topology.junctions.issues.is_empty());
            assert!(
                mesh.topology_valid && mesh.quality_passed,
                "scale={scale} rotated={rotated}: {:?}",
                mesh.blockers
            );
            // Both walls lie in one plane.
            let model = &topology.preview;
            let walls: Vec<_> = model
                .surfaces()
                .iter()
                .filter(|s| s.source_elements.iter().all(|&e| e > 24))
                .collect();
            assert_eq!(walls.len(), 2);
            let plane = &model.planes()[walls[0].plane];
            for s in &walls {
                for e in s.boundaries.iter().flatten() {
                    for &v in &model.edges()[e.edge] {
                        let d = plane.distance(model.vertices()[v]);
                        assert!(d.abs() < 1e-6 * scale, "wall vertex off the line: {d}");
                    }
                }
            }

            // Side by side over 1 m, or only a 25 mm strip where the ends
            // overlap: aligned, the walls would overlap in one plane.
            for lap in [40, 1] {
                let (topology, _) = run_with(
                    wall_line_source(lap),
                    scale,
                    rotated,
                    Some(features.clone()),
                );
                assert!(topology.stacked_walls.aligned.is_empty(), "lap {lap}");
                assert!(!topology.stacked_walls.kept.is_empty());
                assert!(topology
                    .stacked_walls
                    .kept
                    .iter()
                    .all(|k| k.reason.starts_with("overlapping_parallel_walls")));
            }
        }
    }
}

#[test]
fn stacked_wall_adopts_axis_of_the_wall_below() {
    for scale in [0.1, 1., 10.] {
        for rotated in [false, true] {
            let features = assembly::FeaturePolicy {
                maximum_console_width: 0.25 * scale,
                maximum_stack_offset: 0.05 * scale,
                ..Default::default()
            };
            let (topology, mesh) = run_with(stacked_source(), scale, rotated, Some(features));
            let stacked = &topology.stacked_walls;
            assert_eq!(stacked.aligned.len(), 1, "{:?}", stacked.kept);
            assert!((stacked.aligned[0].offset - 0.025 * scale).abs() < 1e-6 * scale);
            // Seven upper-wall base nodes take the lower wall's vertices.
            assert_eq!(stacked.identified.len(), 7);
            assert!(topology.issues.is_empty() && topology.junctions.issues.is_empty());
            assert!(
                mesh.topology_valid && mesh.quality_passed,
                "scale={scale} rotated={rotated}: {:?}",
                mesh.blockers
            );
            // One junction line: both walls and the slab share its mesh edges.
            let mut edges = BTreeMap::<u32, BTreeSet<[usize; 2]>>::new();
            for t in &mesh.triangles {
                for i in 0..3 {
                    let (a, b) = (t.vertices[i], t.vertices[(i + 1) % 3]);
                    edges
                        .entry(t.stiffness)
                        .or_default()
                        .insert([a.min(b), a.max(b)]);
                }
            }
            let line: BTreeSet<_> = edges[&20].intersection(&edges[&30]).copied().collect();
            assert!(line.len() >= 12, "wall/wall shared edges {}", line.len());
            assert!(line.is_subset(&edges[&10]));
        }
    }
    // The conservative assembly keeps both walls where the source put them.
    let (topology, _) = run_input(stacked_source(), 1., false);
    assert!(topology.stacked_walls.aligned.is_empty());
}

/// A 2 m x 2 m slab and, on its edge, a needle triangle of the same
/// stiffness 2 mm high whose nodes are not slab nodes (a degenerate sliver
/// of the source mesh).
fn sliver_source() -> topo_reconstruct_rs::input::MeshData {
    use topo_reconstruct_rs::input::{ElementData, MeshData};
    let mut mesh = MeshData::default();
    for j in 0..3 {
        for i in 0..3 {
            mesh.nodes
                .insert(1 + i + 3 * j, DVec3::new(i as f64, j as f64, 0.));
        }
    }
    for j in 0..2 {
        for i in 0..2 {
            let a = 1 + i + 3 * j;
            mesh.elements.push(ElementData {
                id: mesh.elements.len() as u32 + 1,
                elem_type: 44,
                stiff_id: 1,
                nodes: vec![a, a + 1, a + 4, a + 3],
            });
        }
    }
    mesh.nodes.insert(20, DVec3::new(0.2, 2., 0.));
    mesh.nodes.insert(21, DVec3::new(0.7, 2., 0.));
    mesh.nodes.insert(22, DVec3::new(0.45, 2.002, 0.));
    mesh.elements.push(ElementData {
        id: 5,
        elem_type: 42,
        stiff_id: 1,
        nodes: vec![20, 21, 22],
    });
    // A column standing on the sliver apex.
    mesh.nodes.insert(23, DVec3::new(0.45, 2.002, 1.));
    mesh.elements.push(ElementData {
        id: 6,
        elem_type: 10,
        stiff_id: 2,
        nodes: vec![22, 23],
    });
    mesh
}

#[test]
fn region_no_wider_than_a_crack_is_removed_with_provenance() {
    for scale in [0.1, 1., 10.] {
        for rotated in [false, true] {
            let features = assembly::FeaturePolicy {
                maximum_crack_width: 0.01 * scale,
                ..Default::default()
            };
            let (topology, mesh) = run_with(sliver_source(), scale, rotated, Some(features));
            assert_eq!(topology.removed_slivers.len(), 1);
            let sliver = &topology.removed_slivers[0];
            assert_eq!(sliver.source_elements, vec![5]);
            assert!((sliver.width - 0.002 * scale).abs() < 1e-6 * scale);
            assert!(sliver.mean_width < sliver.width);
            assert_eq!(topology.surface_stiffness.len(), 1);
            assert!(topology.issues.is_empty());
            // The column on a removed sliver node is still built.
            assert!(topology.axis_assembly.issues.is_empty());
            assert_eq!(topology.axis_assembly.axes.len(), 1);
            assert!(
                mesh.topology_valid && mesh.source_coverage_complete,
                "scale={scale} rotated={rotated}: {:?}",
                mesh.blockers
            );
        }
    }
}
