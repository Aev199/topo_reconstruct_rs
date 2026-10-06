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
    let mut cut_at: Option<usize> = None;
    while let Some(a) = args.next() {
        if a == "--frame-cache" {
            cache = args.next().map(std::path::PathBuf::from);
        } else if a == "--combine" {
            combine = true;
        } else if a == "--cut" {
            cut_at = args.next().and_then(|s| s.parse().ok());
        }
    }
    let profile = Profile::plaxis();
    let output = pipeline::run(&input, &profile, &Options { mesh: false, frame_cache: cache }, &mut |_| {})?;
    let mut session = Session::new(&output.topology, profile.audit_options());
    let uncut = session.state().clone();
    if let Some(k) = cut_at {
        let z = topo_reconstruct_rs::reconstruction::assembly::cutoff::floors(session.state())[k].z;
        session.apply(topo_reconstruct_rs::editor::Edit::CutAbove { z }, "probe")?;
        println!("cut at z {z:.3}");
    }
    let materials = std::sync::Arc::new(topo_reconstruct_rs::parsers::lira::LiraParser::parse_materials(&input)?);
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
        loads::Settings { force_factor: 9.80665, snap: profile.edge_collapse, max_groups: 40, combination, cases: None, materials: Some(materials.clone()), cut_loads: true },
    );
    let reference = if cut_at.is_some() {
        Some(loads::transfer(
            &uncut,
            &output.topology.vertex_source_nodes,
            &mesh,
            &set,
            loads::Settings { force_factor: 9.80665, snap: profile.edge_collapse, max_groups: 40, combination: None, cases: None, materials: Some(materials.clone()), cut_loads: true },
        ).1)
    } else { None };
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
    // How many loads of each kind, and how many distinct line / point places.
    let mut counts: BTreeMap<&str, usize> = BTreeMap::new();
    let mut places: std::collections::BTreeSet<String> = Default::default();
    let mut points: std::collections::BTreeSet<String> = Default::default();
    for l in &list {
        match l {
            Load::Surface { .. } => *counts.entry("surface").or_default() += 1,
            Load::Line { start, end, .. } => {
                *counts.entry("line").or_default() += 1;
                places.insert(format!("{:.3?}{:.3?}", start, end));
            }
            Load::Point { at, .. } => {
                *counts.entry("point").or_default() += 1;
                points.insert(format!("{at:.3?}"));
            }
        }
    }
    // Point loads that lie on no plate polygon and no beam: PLAXIS calls them mesh-independent.
    {
        let exchange = topo_reconstruct_rs::plaxis::exchange(
            session.state(),
            &materials,
            topo_reconstruct_rs::plaxis::Settings { force_factor: 9.80665, min_edge: profile.edge_collapse, stiffness: topo_reconstruct_rs::plaxis::StiffnessMode::Effective },
            "check",
        );
        let on_polygon = |p: DVec3, tol: f64| exchange.plates.iter().flat_map(|pl| pl.polygons.iter()).any(|poly| {
            let q: Vec<DVec3> = poly.iter().map(|x| DVec3::from_array(*x)).collect();
            let n: DVec3 = (1..q.len() - 1).map(|i| (q[i] - q[0]).cross(q[i + 1] - q[0])).sum();
            if n.length() < 1e-12 { return false; }
            let n = n.normalize();
            if (p - q[0]).dot(n).abs() > tol { return false; }
            // Inside by the sum of angles / crossing number in the plane.
            let (u, v) = { let u = (q[1] - q[0]).normalize(); (u, n.cross(u)) };
            let pt = |x: DVec3| (((x - q[0]).dot(u)), ((x - q[0]).dot(v)));
            let (px, py) = pt(p);
            let mut odd = false;
            for i in 0..q.len() {
                let (a, b) = (pt(q[i]), pt(q[(i + 1) % q.len()]));
                if (a.1 > py) != (b.1 > py) && px < a.0 + (py - a.1) / (b.1 - a.1) * (b.0 - a.0) { odd = !odd; }
            }
            if odd { return true; }
            // On the boundary within the tolerance.
            (0..q.len()).any(|i| {
                let (a, b) = (q[i], q[(i + 1) % q.len()]);
                let d = b - a;
                let t = ((p - a).dot(d) / d.length_squared()).clamp(0., 1.);
                (p - (a + d * t)).length() <= tol
            })
        });
        let on_beam = |p: DVec3, tol: f64| exchange.beams.iter().any(|b| {
            let (a, c) = (DVec3::from_array(b.start), DVec3::from_array(b.end));
            let d = c - a;
            let t = ((p - a).dot(d) / d.length_squared()).clamp(0., 1.);
            (p - (a + d * t)).length() <= tol
        });
        let polygons: Vec<Vec<[f64; 3]>> = exchange.plates.iter().flat_map(|p| p.polygons.iter().cloned()).collect();
        let segments: Vec<([f64; 3], [f64; 3])> = exchange.beams.iter().map(|b| (b.start, b.end)).collect();
        let mut line_off = 0;
        for l in &list {
            if let Load::Line { start, end, .. } = l {
                let (a, b) = (DVec3::from_array(*start), DVec3::from_array(*end));
                if !(on_polygon(a, 1e-3) || on_beam(a, 1e-3)) || !(on_polygon(b, 1e-3) || on_beam(b, 1e-3)) {
                    line_off += 1;
                }
            }
        }
        println!("line loads with an end off structure: {line_off}");
        let mut list = list.clone();
        let (moved, farthest) = loads::attach_points(&mut list, &polygons, &segments);
        println!("attached {moved} point loads, farthest shift {farthest:.4} m");
        let mut off = 0;
        for l in &list {
            if let Load::Point { case, at, force, .. } = l {
                let p = DVec3::from_array(*at);
                if !on_polygon(p, 1e-3) && !on_beam(p, 1e-3) {
                    off += 1;
                    if off <= 12 {
                        println!("point off structure: case {case} at {at:?} force {force:?}");
                    }
                }
            }
        }
        println!("point loads off structure: {off} of {}", counts.get("point").copied().unwrap_or(0));
    }
    println!("counts {counts:?}, distinct lines {}, distinct points {}", places.len(), points.len());
    for c in &report.cases {
        let full = reference.as_ref().and_then(|r| r.cases.iter().find(|x| x.case == c.case)).map(|x| x.exported);
        println!("case {} {:?}: source {:?} exported {:?} uncut {:?}", c.case, c.name, c.source, c.exported, full);
    }
    println!("surface-only exported: {exported:?}");
    println!("skipped {:?}\napprox {:?}", report.skipped, report.approximated);
    Ok(())
}
