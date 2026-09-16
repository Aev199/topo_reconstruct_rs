//! Joint geometric proposal. Not a meshing-ready model.
use topo_reconstruct_rs::{
    parsers::LiraParser,
    reconstruction::{assembly, frame, graph, planes, recognize},
};
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let path = std::env::args()
        .nth(1)
        .ok_or("Usage: recognize_frame model.txt [iterations] [sliding_steps]")?;
    let iterations = std::env::args()
        .nth(2)
        .map(|v| v.parse::<usize>())
        .transpose()?
        .unwrap_or(1000);
    let sliding_steps = std::env::args()
        .nth(3)
        .map(|v| v.parse::<usize>())
        .transpose()?
        .unwrap_or(0);
    let mesh = LiraParser::parse(path)?;
    let axes = recognize::recognize(
        &mesh,
        &recognize::Policy {
            angle: 0.02,
            line_tolerance: 0.01,
            numerical_precision: 1e-8,
        },
    )?;
    let planes = planes::recognize(
        &mesh,
        &planes::Policy {
            angle: 0.02,
            distance: 0.01,
            precision: 1e-8,
        },
    )?;
    let policy = frame::Policy {
        up: [0., 0., 1.],
        angle: 0.02,
        maximum_movement: 0.15,
        relative_movement: 0.05,
        minimum_length: 0.03,
        residual_tolerance: 1e-7,
        iterations,
    };
    if sliding_steps > 0 {
        let result = frame::solve_sliding(&mesh, &axes, &planes, &policy, sliding_steps)?;
        serde_json::to_writer_pretty(
            std::io::stdout().lock(),
            &serde_json::json!({
                "proposal_only":true, "constraint_graph":graph::Graph::from_frame(&result),
                "frame":result, "axis_recognition":axes, "plane_recognition":planes,
            }),
        )?;
        return Ok(());
    }
    let result = frame::solve_with_retry(&mesh, &axes, &planes, &policy, 3)?;
    let topology = assembly::assemble(
        &mesh,
        &result,
        &assembly::Policy {
            closure_tolerance: 0.001,
            junction_movement_limit: 0.05,
            precision: 1e-7,
            minimum_edge: 0.001,
        },
    )?;
    serde_json::to_writer_pretty(
        std::io::stdout().lock(),
        &serde_json::json!({"constraint_graph":graph::Graph::from_frame(&result),"frame":result,"topology":topology,"axis_recognition":axes,"plane_recognition":planes}),
    )?;
    Ok(())
}
