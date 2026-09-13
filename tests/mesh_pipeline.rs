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
            assert_eq!(topology.preview.surfaces().len(), 4);
            assert_eq!(topology.axis_assembly.axes.len(), 2);
            assert_eq!(
                topology
                    .preview
                    .surfaces()
                    .iter()
                    .filter(|s| s.contours.len() == 2)
                    .count(),
                1
            );
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
    topology.all_surface_patches_built = false;
    policy.maximum_area = 0.5;
    assert!(mesh::build(&topology, &policy).is_err());
}

#[test]
fn nonorthogonal_fragment_reports_bad_quality_without_losing_topology() {
    let mut input = fixture::source();
    for p in input.nodes.values_mut() {
        p.x += 0.2 * p.y;
    }
    let (_, mesh) = run_input(input, 1., true);
    assert!(mesh.topology_valid, "{:?}", mesh.blockers);
    // Locked constraints may prevent quality refinement even in valid geometry.
    // This must stay an explicit failure, never a mesh-ready success.
    assert!(!mesh.quality_passed);
    assert!(mesh.blockers.iter().any(|b| b == "minimum_angle_not_met"));
    assert!(!mesh.export_ready);
}
