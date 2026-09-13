//! Run FE input through reconstruction and constrained mesh generation.
#[path = "support/mesh_fixture.rs"]
mod fixture;
use topo_reconstruct_rs::{
    parsers::LiraParser,
    reconstruction::{assembly, frame, mesh, planes, recognize},
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let source = if let Some(path) = std::env::args().nth(1) {
        LiraParser::parse(path)?
    } else {
        fixture::source()
    };
    let axes = recognize::recognize(
        &source,
        &recognize::Policy {
            angle: 0.02,
            line_tolerance: 0.001,
            numerical_precision: 1e-8,
        },
    )?;
    let planes = planes::recognize(
        &source,
        &planes::Policy {
            angle: 0.02,
            distance: 0.001,
            precision: 1e-8,
        },
    )?;
    let frame = frame::solve(
        &source,
        &axes,
        &planes,
        &frame::Policy {
            up: [0., 0., 1.],
            angle: 0.02,
            maximum_movement: 0.05,
            relative_movement: 0.05,
            minimum_length: 0.01,
            residual_tolerance: 1e-8,
            iterations: 1000,
        },
    )?;
    let topology = assembly::assemble(
        &source,
        &frame,
        &assembly::Policy {
            closure_tolerance: 0.001,
            junction_movement_limit: 0.05,
            precision: 1e-7,
            minimum_edge: 0.001,
        },
    )?;
    let mesh = mesh::build(
        &topology,
        &mesh::Policy {
            boundary_spacing: 0.5,
            maximum_area: 0.5,
            minimum_angle_degrees: 20.,
            maximum_added_vertices_per_surface: 10000,
        },
    )?;
    let passed = mesh.topology_valid && mesh.quality_passed;
    serde_json::to_writer_pretty(
        std::io::stdout().lock(),
        &serde_json::json!({"frame":frame,"topology":topology,"mesh":mesh}),
    )?;
    if !passed {
        return Err("mesh validation failed; see report".into());
    }
    Ok(())
}
