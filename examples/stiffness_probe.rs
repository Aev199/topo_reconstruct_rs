//! The MIDAS materials, sections and thicknesses of a LIRA model (development).
//! Usage: stiffness_probe MODEL.txt [--lira] [--frame-cache PATH]
use std::collections::BTreeSet;
use topo_reconstruct_rs::editor::Session;
use topo_reconstruct_rs::midas_stiffness::{build, Mode, Options};
use topo_reconstruct_rs::parsers::lira::LiraParser;
use topo_reconstruct_rs::pipeline::{self, Options as PipelineOptions, Profile};

fn main() -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let mut args = std::env::args().skip(1);
    let input = std::path::PathBuf::from(args.next().ok_or("usage: stiffness_probe MODEL.txt")?);
    let (mut mode, mut cache) = (Mode::Converter, None);
    while let Some(a) = args.next() {
        match a.as_str() {
            "--lira" => mode = Mode::Lira,
            "--frame-cache" => cache = args.next().map(std::path::PathBuf::from),
            _ => {}
        }
    }
    let profile = Profile::plaxis();
    let output = pipeline::run(&input, &profile, &PipelineOptions { mesh: false, frame_cache: cache }, &mut |_| {})?;
    let session = Session::new(&output.topology, profile.audit_options());
    let state = session.state();
    let bytes = std::fs::read(&input)?;
    let materials = LiraParser::materials_from(&bytes);
    let profiles = LiraParser::profiles_from(&bytes);
    let plates: BTreeSet<u32> = state.stiffness.iter().copied().collect();
    let bars: BTreeSet<u32> = state.axes.iter().flat_map(|a| a.spans.iter().map(|s| s.stiffness)).collect();
    let s = build(&materials, &profiles, &plates, &bars, Options { mode, density_multiplier: 1. });
    println!("types: {} plate, {} bar; profiles {}; {} materials, {} sections, {} thicknesses", plates.len(), bars.len(), profiles.len(), s.material_lines.len(), s.sections, s.thickness_lines.len());
    println!("missing types: {:?}", s.missing);
    for n in &s.notes {
        println!("note: {n}");
    }
    for l in s.material_lines.iter().chain(&s.section_lines).chain(&s.thickness_lines).take(60) {
        println!("{l}");
    }
    Ok(())
}
