#[path = "../examples/support/mesh_fixture.rs"]
mod fixture;
use glam::{DQuat, DVec3};
use std::collections::{BTreeMap, BTreeSet};
use topo_reconstruct_rs::reconstruction::{assembly, frame, mesh, planes, recognize};

fn run(scale: f64, rotated: bool) -> (assembly::Report, mesh::Report) {
    run_input(fixture::source(), scale, rotated)
}

fn run_input(
    mut input: topo_reconstruct_rs::input::MeshData,
    scale: f64,
    rotated: bool,
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
        },
    )
    .unwrap();
    let t = assembly::assemble(
        &input,
        &f,
        &assembly::Policy {
            closure_tolerance: 0.001 * scale,
            junction_movement_limit: 0.05 * scale,
            precision: 1e-7 * scale,
            minimum_edge: 0.001 * scale,
        },
    )
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
        surface_junctions: topo_reconstruct_rs::reconstruction::junctions::Report::default(),
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
