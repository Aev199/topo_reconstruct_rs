//! Floors of a reconstructed model (development): horizontal surfaces by
//! elevation, and how many surfaces cross each level.
use topo_reconstruct_rs::editor::Session;
use topo_reconstruct_rs::pipeline::{self, Options, Profile};

fn main() -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let mut args = std::env::args().skip(1);
    let input = std::path::PathBuf::from(args.next().ok_or("usage: floors_probe MODEL.txt")?);
    let cache = args.next().map(std::path::PathBuf::from);
    let cut_at: Option<usize> = args.next().and_then(|a| a.parse().ok());
    let profile = Profile::plaxis();
    let output = pipeline::run(&input, &profile, &Options { mesh: false, frame_cache: cache }, &mut |_| {})?;
    let mut session = Session::new(&output.topology, profile.audit_options());
    if let Some(k) = cut_at {
        let f = topo_reconstruct_rs::reconstruction::assembly::cutoff::floors(session.state());
        let z = f[k].z;
        let before = session.audit().counts.clone();
        let t = std::time::Instant::now();
        match session.apply(topo_reconstruct_rs::editor::Edit::CutAbove { z }, "probe") {
            Ok(_) => {
                let cut = session.state().cut.clone().unwrap();
                println!("cut at floor {k} (z {z:.3}) in {:?}: {} surfaces removed, {} clipped, {} bars removed, {} trimmed; supports: {} wall lines, {} columns", t.elapsed(), cut.removed_surfaces, cut.clipped_surfaces, cut.removed_bars, cut.trimmed_bars, cut.walls.len(), cut.columns.len());
                let materials = topo_reconstruct_rs::parsers::lira::LiraParser::parse_materials(&input)?;
                match topo_reconstruct_rs::storeys::with_cap(session.state(), &materials, 1.0) {
                    Some((_, _, r)) => println!("cap: {} surfaces, equivalent thickness {:.3} m (Ix {:.1}, Iy {:.1} m4; storey {:.2} m, removed {:.2} m)", r.surfaces, r.equivalent_thickness, r.ix, r.iy, r.storey, r.removed_height),
                    None => println!("cap: nothing to change"),
                }
                println!("audit before {before:?}\naudit after  {:?}\npassed {}", session.audit().counts, session.audit().passed);
            }
            Err(e) => println!("cut failed: {e}"),
        }
    }
    let state = session.state();
    let model = &state.model;
    // Per surface: z range, horizontal?
    let mut slabs: Vec<(f64, usize, f64)> = vec![];
    let mut ranges = vec![];
    for (s, surface) in model.surfaces().iter().enumerate() {
        let plane = &model.planes()[surface.plane];
        let n = plane.normal();
        let zs: Vec<f64> = model.surface_edges(s).flat_map(|e| model.edges()[e]).map(|v| model.vertices()[v][2]).collect();
        let (lo, hi) = zs.iter().fold((f64::MAX, f64::MIN), |(l, h), z| (l.min(*z), h.max(*z)));
        ranges.push((lo, hi, n[2].abs()));
        if n[2].abs() > 0.98 {
            slabs.push((lo, s, 0.));
        }
    }
    slabs.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut levels: Vec<(f64, usize)> = vec![];
    for (z, _, _) in &slabs {
        match levels.last_mut() {
            Some(l) if (z - l.0).abs() < 0.3 => l.1 += 1,
            _ => levels.push((*z, 1)),
        }
    }
    println!("{} surfaces, {} horizontal, {} levels", ranges.len(), slabs.len(), levels.len());
    for f in topo_reconstruct_rs::reconstruction::assembly::cutoff::floors(state) { println!("  floor z {:.3} slabs {} area {:.1}", f.z, f.slabs, f.area); }
    for (z, n) in &levels {
        let crossing = ranges.iter().filter(|r| r.2 < 0.98 && r.0 < z - 0.05 && r.1 > z + 0.05).count();
        let above = ranges.iter().filter(|r| r.0 >= z - 0.05).count();
        println!("  z {z:8.3}: {n:3} slabs, {crossing:3} crossing surfaces, {above} surfaces from here up");
    }
    println!("bars: {} ({} vertical)", state.axes.len(), state.axes.iter().filter(|a| {
        let p = model.vertices()[a.endpoints[0]]; let q = model.vertices()[a.endpoints[1]];
        ((q[0]-p[0]).powi(2) + (q[1]-p[1]).powi(2)).sqrt() < 0.05 * (q[2]-p[2]).abs()
    }).count());
    Ok(())
}
