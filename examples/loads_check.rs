//! Transfer the loads of a LIRA file onto its reconstructed geometry and print
//! per-case source and exported resultants (development check).
//! Usage: loads_check MODEL.txt [--frame-cache PATH]
use glam::DVec3;
use std::collections::BTreeMap;
use topo_reconstruct_rs::editor::Session;
use topo_reconstruct_rs::loads::{self, Load};
use topo_reconstruct_rs::pipeline::{self, Options, Profile};

fn main() -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let mut args = std::env::args().skip(1);
    let input = std::path::PathBuf::from(args.next().ok_or("usage: loads_check MODEL.txt")?);
    let mut cache = None;
    let mut combine = false;
    while let Some(a) = args.next() {
        if a == "--frame-cache" {
            cache = args.next().map(std::path::PathBuf::from);
        } else if a == "--combine" {
            combine = true;
        }
    }
    let profile = Profile::plaxis();
    let output = pipeline::run(&input, &profile, &Options { mesh: false, frame_cache: cache }, &mut |_| {})?;
    let session = Session::new(&output.topology, profile.audit_options());
    let bytes = std::fs::read(&input)?;
    let mesh = topo_reconstruct_rs::parsers::lira::LiraParser::mesh_from(&bytes)?;
    let set = topo_reconstruct_rs::parsers::loads::parse(&bytes);
    // --combine: every case with factor 1 except the self-weight, simplified.
    let combination = combine.then(|| loads::Combination {
        factors: set.cases.iter().filter(|(_, n)| !loads::is_self_weight(n) && !loads::is_stage(n) && !loads::is_dynamic(n)).map(|(c, _)| (*c, 1.)).collect(),
        simplify: Some(loads::Simplify::default()),
    });
    let (list, report) = loads::transfer(
        session.state(),
        &output.topology.vertex_source_nodes,
        &mesh,
        &set,
        loads::Settings { force_factor: 9.80665, snap: profile.edge_collapse, max_groups: 40, combination, cases: None },
    );
    let mut exported: BTreeMap<u32, DVec3> = BTreeMap::new();
    let mut by_kind: BTreeMap<(u32, &str), DVec3> = BTreeMap::new();
    for l in &list {
        match l {
            Load::Surface { case, polygons, sigma, .. } => {
                let area: f64 = polygons.iter().map(|p| {
                    let p: Vec<DVec3> = p.iter().map(|x| DVec3::from_array(*x)).collect();
                    let mut n = DVec3::ZERO;
                    for i in 0..p.len() { n += p[i].cross(p[(i + 1) % p.len()]); }
                    n.length() / 2.
                }).sum();
                *exported.entry(*case).or_default() += DVec3::from_array(*sigma) * area;
            }
            Load::Line { case, start, end, q_start, q_end } => {
                let len = DVec3::from_array(*start).distance(DVec3::from_array(*end));
                *by_kind.entry((*case, "line")).or_default() += (DVec3::from_array(*q_start) + DVec3::from_array(*q_end)) / 2. * len;
            }
            Load::Point { case, force, .. } => {
                *by_kind.entry((*case, "point")).or_default() += DVec3::from_array(*force);
            }
        }
    }
    println!("by kind: {by_kind:?}");
    for c in &report.cases {
        println!("case {} {:?}: source {:?} exported {:?}", c.case, c.name, c.source, c.exported);
    }
    println!("surface-only exported: {exported:?}");
    println!("skipped {:?}\napprox {:?}", report.skipped, report.approximated);
    Ok(())
}
