//! Axis recognition only. No reconstructed surfaces, loads or remeshing yet.
use topo_reconstruct_rs::{
    parsers::LiraParser,
    reconstruction::recognize::{recognize, Policy},
};
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let path = std::env::args()
        .nth(1)
        .ok_or("Usage: recognize_lira model.txt")?;
    let mesh = LiraParser::parse(path)?;
    let policy = Policy {
        angle: 0.02,
        line_tolerance: 0.01,
        numerical_precision: 1e-8,
    };
    let report = recognize(&mesh, &policy)?;
    serde_json::to_writer_pretty(std::io::stdout().lock(), &report)?;
    Ok(())
}
