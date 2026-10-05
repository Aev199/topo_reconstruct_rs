//! Loads of a LIRA file mapped onto the reconstructed geometry: the resultant
//! must survive, whatever the node order, placement or scale of the model.
use glam::{DQuat, DVec3};
use topo_reconstruct_rs::editor::Session;
use topo_reconstruct_rs::loads::{self, Load};
use topo_reconstruct_rs::pipeline::{self, Options, Profile};

/// A 6 x 4 m slab of 1 m plates (tensor-product numbering for the elements
/// of even index, perimeter order for the others), `n` loaded by case:
/// case 1: 2 tf/m2 downward on the left half, 3 on the right half (global);
/// case 2: a line load along the edge nodes 1-2 of element 1 (positive values act against the axes, as in LIRA).
fn source(transform: impl Fn(DVec3) -> DVec3) -> String {
    let (nx, ny) = (6usize, 4usize);
    let node = |i: usize, j: usize| j * (nx + 1) + i + 1;
    let mut elements = String::new();
    let mut count = 0;
    for j in 0..ny {
        for i in 0..nx {
            let (a, b, c, d) = (node(i, j), node(i + 1, j), node(i + 1, j + 1), node(i, j + 1));
            count += 1;
            if count % 2 == 0 {
                elements += &format!("41 1 {a} {b} {c} {d} /\n");
            } else {
                elements += &format!("41 1 {a} {b} {d} {c} /\n");
            }
        }
    }
    let mut nodes = String::new();
    for j in 0..=ny {
        for i in 0..=nx {
            let p = transform(DVec3::new(i as f64, j as f64, 0.));
            nodes += &format!("{:.6} {:.6} {:.6} /\n", p.x, p.y, p.z);
        }
    }
    let mut loads = String::new();
    let mut n = 0;
    for _ in 0..ny {
        for i in 0..nx {
            n += 1;
            let parameters = if i < nx / 2 { 1 } else { 2 };
            loads += &format!("{n} 16 3 {parameters} 1 /\n");
        }
    }
    // Case 2: a line load along edge 1-2 of element 1 (global z, 4 tf/m).
    loads += "1 19 3 3 2 /\n";
    format!(
        "( 0/ 1; SLAB/ 2; 5/\n39;\n1: LEFT RIGHT ;\n2: EDGE ;\n/\n)\n\
         ( 1/\n{elements})\n( 3/\n1 GEI 0.305915 0.17 0.2 RO 0.254929 /\n)\n( 4/\n{nodes})\n\
         ( 6/\n{loads})\n( 7/\n1 2 0 / 2 3 0 / 3 1 2 4 /\n)\n"
    )
}

fn run(text: &str, tag: &str) -> (Vec<Load>, loads::Report) {
    let path = std::env::temp_dir().join(format!("loads_transfer_{tag}_{}.txt", std::process::id()));
    std::fs::write(&path, text).unwrap();
    let profile = Profile::plaxis();
    let output = pipeline::run(&path, &profile, &Options::default(), &mut |_| {}).unwrap();
    let session = Session::new(&output.topology, profile.audit_options());
    let mesh = topo_reconstruct_rs::parsers::lira::LiraParser::mesh_from(text.as_bytes()).unwrap();
    let set = topo_reconstruct_rs::parsers::loads::parse(text.as_bytes());
    let result = loads::transfer(
        session.state(),
        &output.topology.vertex_source_nodes,
        &mesh,
        &set,
        loads::Settings { force_factor: 10., snap: profile.edge_collapse, max_groups: 40 },
    );
    let _ = std::fs::remove_file(path);
    result
}

fn resultant(list: &[Load], case: u32) -> DVec3 {
    list.iter()
        .filter(|l| l.case() == case)
        .map(|l| match l {
            Load::Surface { polygons, sigma, .. } => {
                let area: f64 = polygons
                    .iter()
                    .map(|p| {
                        let p: Vec<DVec3> = p.iter().map(|&x| DVec3::from_array(x)).collect();
                        (0..p.len()).map(|i| p[i].cross(p[(i + 1) % p.len()])).sum::<DVec3>().length() / 2.
                    })
                    .sum();
                DVec3::from_array(*sigma) * area
            }
            Load::Line { start, end, q_start, q_end, .. } => {
                (DVec3::from_array(*q_start) + DVec3::from_array(*q_end)) / 2.
                    * DVec3::from_array(*start).distance(DVec3::from_array(*end))
            }
            Load::Point { force, .. } => DVec3::from_array(*force),
        })
        .sum()
}

fn check(transform: impl Fn(DVec3) -> DVec3, tag: &str, unit: f64) {
    let (list, report) = run(&source(&transform), tag);
    // Left half 3 x 4 m at -2 tf/m2, right half at -3: -60 tf = -600 kN (force factor 10).
    let expected_case1 = DVec3::new(0., 0., -(2. + 3.) * 12. * 10. * unit * unit);
    let r1 = resultant(&list, 1);
    assert!((r1 - expected_case1).length() < 1e-6 * expected_case1.length(), "{tag}: {r1:?} vs {expected_case1:?}");
    // The two halves are two values on one surface: two surface loads of whole sub-regions.
    let surfaces: Vec<&Load> = list.iter().filter(|l| l.case() == 1).collect();
    assert_eq!(surfaces.len(), 2, "{tag}: {surfaces:?}");
    // Edge load: 1 m at 4 tf/m downward.
    let r2 = resultant(&list, 2);
    let expected_case2 = DVec3::new(0., 0., -4. * 10. * unit);
    assert!((r2 - expected_case2).length() < 1e-6 * expected_case2.length() + 1e-9, "{tag}: {r2:?} vs {expected_case2:?}");
    for c in &report.cases {
        let source = DVec3::from_array(c.source);
        let exported = DVec3::from_array(c.exported);
        assert!((source - exported).length() < 1e-6 * source.length().max(1.), "{tag} case {}: {source:?} vs {exported:?}", c.case);
    }
    assert!(report.skipped.is_empty(), "{tag}: {:?}", report.skipped);
}

#[test]
fn slab_loads_keep_their_resultant_in_any_placement() {
    check(|p| p, "plain", 1.);
    let q = DQuat::from_axis_angle(DVec3::new(1., 2., 3.).normalize(), 0.7);
    check(|p| p * 2., "scaled", 2.);
    check(move |p| q * p + DVec3::new(100., -50., 20.), "rotated", 1.);
}
