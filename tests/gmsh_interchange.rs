#[path = "../examples/support/mesh_fixture.rs"]
mod fixture;

use glam::{DQuat, DVec3};
use std::collections::BTreeSet;
use topo_reconstruct_rs::reconstruction::{assembly, frame, gmsh, planes, recognize};

fn build(scale: f64, rotated: bool) -> gmsh::Interchange {
    let mut input = fixture::source();
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
        for element in &mut input.elements {
            for node in &mut element.nodes {
                *node = 1000 - *node;
            }
        }
    }

    let axes = recognize::recognize(
        &input,
        &recognize::Policy {
            angle: 0.02,
            line_tolerance: 0.001 * scale,
            numerical_precision: 1e-8 * scale,
        },
    )
    .unwrap();
    let planes = planes::recognize(
        &input,
        &planes::Policy {
            angle: 0.02,
            distance: 0.001 * scale,
            precision: 1e-8 * scale,
        },
    )
    .unwrap();
    let frame = frame::solve(
        &input,
        &axes,
        &planes,
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
    let topology = assembly::assemble(
        &input,
        &frame,
        &assembly::Policy {
            closure_tolerance: 0.001 * scale,
            junction_movement_limit: 0.05 * scale,
            precision: 1e-7 * scale,
            minimum_edge: 0.001 * scale,
        },
    )
    .unwrap();

    gmsh::from_assembly(&topology, 0.75 * scale).unwrap()
}

#[test]
fn interchange_preserves_surface_property_and_source_ownership() {
    let data = build(1.0, false);
    assert_eq!(data.format, gmsh::FORMAT);
    assert!(data.source_coverage_complete, "{:?}", data.blockers);
    assert_eq!(data.surfaces.len(), 4);

    let stiffness: BTreeSet<_> = data.surfaces.iter().map(|s| s.stiffness).collect();
    assert_eq!(stiffness, BTreeSet::from([10, 20, 30, 40]));
    assert!(data
        .surfaces
        .iter()
        .all(|surface| !surface.source_elements.is_empty() && !surface.rings.is_empty()));
    assert!(data
        .surfaces
        .iter()
        .flat_map(|surface| &surface.rings)
        .all(|ring| ring.len() >= 3 && ring.iter().flatten().all(|x| x.is_finite())));
    assert!(data.surfaces.iter().all(|surface| {
        surface.rings.len() == surface.ring_source_nodes.len()
            && surface
                .rings
                .iter()
                .zip(&surface.ring_source_nodes)
                .all(|(ring, source_nodes)| ring.len() == source_nodes.len())
    }));
}

#[test]
fn interchange_preserves_axis_spans_and_explicit_contacts() {
    let data = build(1.0, false);
    assert_eq!(data.axes.len(), 3);
    assert!(!data.contacts.is_empty());

    let bar_stiffness: BTreeSet<_> = data
        .axes
        .iter()
        .flat_map(|axis| axis.property_spans.iter().map(|span| span.stiffness))
        .collect();
    assert_eq!(bar_stiffness, BTreeSet::from([50, 60, 70]));
    assert!(data.axes.iter().all(|axis| {
        DVec3::from_array(axis.endpoints[0]).distance(DVec3::from_array(axis.endpoints[1])) > 0.0
            && !axis.property_spans.is_empty()
            && axis
                .property_spans
                .iter()
                .all(|span| span.start_t >= 0.0 && span.end_t <= 1.0 && span.end_t > span.start_t)
    }));
}

#[test]
fn interchange_is_invariant_in_structure_under_scale_rotation_and_renumbering() {
    let baseline = build(1.0, false);
    let changed = build(10.0, true);

    assert_eq!(baseline.surfaces.len(), changed.surfaces.len());
    assert_eq!(baseline.axes.len(), changed.axes.len());
    assert_eq!(baseline.contacts.len(), changed.contacts.len());
    let surface_signature = |data: &gmsh::Interchange| {
        let mut values = data
            .surfaces
            .iter()
            .map(|surface| (surface.stiffness, surface.source_elements.len(), surface.rings.len()))
            .collect::<Vec<_>>();
        values.sort_unstable();
        values
    };
    assert_eq!(surface_signature(&baseline), surface_signature(&changed));

    let axis_signature = |data: &gmsh::Interchange| {
        let mut values = data
            .axes
            .iter()
            .map(|axis| {
                let mut stiffness = axis
                    .property_spans
                    .iter()
                    .map(|span| span.stiffness)
                    .collect::<Vec<_>>();
                stiffness.sort_unstable();
                stiffness
            })
            .collect::<Vec<_>>();
        values.sort();
        values
    };
    assert_eq!(axis_signature(&baseline), axis_signature(&changed));

    let axis_lengths = |data: &gmsh::Interchange| {
        let mut values = data
            .axes
            .iter()
            .map(|axis| {
                DVec3::from_array(axis.endpoints[0])
                    .distance(DVec3::from_array(axis.endpoints[1]))
            })
            .collect::<Vec<_>>();
        values.sort_by(f64::total_cmp);
        values
    };
    for (baseline_length, changed_length) in axis_lengths(&baseline)
        .into_iter()
        .zip(axis_lengths(&changed))
    {
        assert!((changed_length / baseline_length - 10.0).abs() < 1e-8);
    }
}
