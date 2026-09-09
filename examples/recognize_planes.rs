//! Plane hypotheses only, before contour reconstruction or meshing.
use topo_reconstruct_rs::{
    parsers::LiraParser,
    reconstruction::planes::{recognize, Policy},
};
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let path = std::env::args()
        .nth(1)
        .ok_or("Usage: recognize_planes model.txt")?;
    let mesh = LiraParser::parse(path)?;
    let result = recognize(
        &mesh,
        &Policy {
            angle: 0.02,
            distance: 0.01,
            precision: 1e-8,
        },
    )?;
    serde_json::to_writer_pretty(std::io::stdout().lock(), &result)?;
    Ok(())
}
