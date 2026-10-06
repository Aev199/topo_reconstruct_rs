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
fn source(transform: impl Fn(DVec3) -> DVec3, loaded: impl Fn(usize, usize) -> bool) -> String {
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
    for j in 0..ny {
        for i in 0..nx {
            n += 1;
            let parameters = if i < nx / 2 { 1 } else { 2 };
            if loaded(i, j) {
                loads += &format!("{n} 16 3 {parameters} 1 /\n");
            }
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
    run_with(text, tag, None)
}

fn run_with(text: &str, tag: &str, combination: Option<loads::Combination>) -> (Vec<Load>, loads::Report) {
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
        loads::Settings { force_factor: 10., snap: profile.edge_collapse, max_groups: 40, combination, cases: None, materials: None, cut_loads: true, cut_distribution: Default::default() },
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
    let (list, report) = run(&source(&transform, |_, _| true), tag);
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

fn combination(factors: &[(u32, f64)]) -> Option<loads::Combination> {
    Some(loads::Combination {
        factors: factors.iter().copied().collect(),
        simplify: Some(loads::Simplify::default()),
    })
}

#[test]
fn combination_is_one_uniform_load_per_plate_and_drops_cases_without_a_factor() {
    let q = DQuat::from_axis_angle(DVec3::new(1., 2., 3.).normalize(), 0.7);
    let place = move |p: DVec3| q * p + DVec3::new(100., -50., 20.);
    let text = source(place, |_, _| true);
    // Case 1 weighted 1.5, case 2 left out (as the self-weight case would be).
    let (list, report) = run_with(&text, "combined", combination(&[(1, 1.5)]));
    assert_eq!(list.len(), 1, "{list:?}");
    assert!(matches!(list[0], Load::Surface { case: 0, .. }));
    let expected = DVec3::new(0., 0., -600. * 1.5);
    let got = resultant(&list, 0);
    assert!((got - expected).length() < 1e-6 * expected.length(), "{got:?} vs {expected:?}");
    assert_eq!(report.cases.len(), 1);
    // Both cases (weighted 1 and 2): the edge load between two vertices stays a line load.
    let (list, _) = run_with(&text, "combined2", combination(&[(1, 1.), (2, 2.)]));
    let expected = DVec3::new(0., 0., -600. - 2. * 40.);
    let got = resultant(&list, 0);
    assert!((got - expected).length() < 1e-6 * expected.length(), "{got:?} vs {expected:?}");
    assert_eq!(list.iter().filter(|l| matches!(l, Load::Line { .. })).count(), 1);
    assert_eq!(list.len(), 2);
}

#[test]
fn a_small_patch_becomes_a_point_load_at_its_centre() {
    let text = source(|p| p, |i, j| i == 0 && j == 0);
    let (list, _) = run_with(&text, "patch", combination(&[(1, 1.)]));
    assert_eq!(list.len(), 1, "{list:?}");
    let Load::Point { at, force, .. } = &list[0] else { panic!("{list:?}") };
    assert!((DVec3::from_array(*at) - DVec3::new(0.5, 0.5, 0.)).length() < 1e-9, "{at:?}");
    assert!((DVec3::from_array(*force) - DVec3::new(0., 0., -20.)).length() < 1e-9, "{force:?}");
}

/// The moment of a case about the origin and its force, source and exported (the report).
fn case_of(report: &loads::Report, case: u32) -> &loads::CaseReport {
    report.cases.iter().find(|c| c.case == case).unwrap_or_else(|| panic!("case {case}: {:?}", report.cases))
}

fn near(a: [f64; 3], b: [f64; 3], tolerance: f64) -> bool {
    (DVec3::from_array(a) - DVec3::from_array(b)).length() <= tolerance
}

/// Pressures of opposite signs on the halves of a plate (a couple): the combination must
/// keep the moment, not drop the pair because its force is zero.
#[test]
fn a_couple_of_pressures_survives_the_combination() {
    let text = source(|p| p, |_, _| true).replace("( 7/\n1 2 0 / 2 3 0 / 3 1 2 4 /", "( 7/\n1 2 0 / 2 -2 0 / 3 1 2 4 /");
    let (list, report) = run_with(&text, "couple", combination(&[(1, 1.)]));
    let c = case_of(&report, 0);
    assert!(DVec3::from_array(c.source).length() < 1e-6, "{c:?}");
    assert!(!list.is_empty(), "the couple was dropped");
    // 2 tf/m2 x 12 m2 down on the left, up on the right: 240 kN at x = 1.5 and 4.5, My = 720 kN m.
    assert!(near(c.exported_moment, c.source_moment, 1e-6 * c.moment_scale), "{c:?}");
    assert!(DVec3::from_array(c.source_moment).length() > 700., "{c:?}");
    // The exported loads alone give the same moment about the origin of the report.
    let origin = DVec3::from_array(report.origin);
    let mut moment = DVec3::ZERO;
    for l in &list {
        match l {
            Load::Point { at, force, .. } => moment += (DVec3::from_array(*at) - origin).cross(DVec3::from_array(*force)),
            Load::Surface { .. } | Load::Line { .. } => panic!("a couple is not a distributed load: {l:?}"),
        }
    }
    assert!(near(moment.to_array(), c.source_moment, 1e-6 * c.moment_scale), "{moment:?} {c:?}");
}

/// A load of one sign with a smaller one of the other (net force kept, moment sign kept).
#[test]
fn unequal_pressures_of_both_signs_keep_force_and_moment() {
    let text = source(|p| p, |_, _| true).replace("( 7/\n1 2 0 / 2 3 0 / 3 1 2 4 /", "( 7/\n1 2 0 / 2 -1 0 / 3 1 2 4 /");
    let (_, report) = run_with(&text, "unequal", combination(&[(1, 1.)]));
    let c = case_of(&report, 0);
    assert!(near(c.exported, c.source, 1e-6 * DVec3::from_array(c.source).length()), "{c:?}");
    assert!(near(c.exported_moment, c.source_moment, 1e-6 * c.moment_scale), "{c:?}");
}

/// The same square element numbered in tensor order and in perimeter order, with the edge
/// load given by the local node numbers of each: the same load on the same edge.
#[test]
fn an_edge_load_follows_the_local_node_numbers_of_the_element() {
    let build = |elements: &str, line: &str| {
        format!(
            "( 0/ 1; E/ 2; 5/\n39;\n1: C ;\n/\n)\n( 1/\n{elements})\n( 3/\n1 GEI 0.305915 0.17 0.2 RO 0.254929 /\n)\n\
             ( 4/\n0 0 0 /\n1 0 0 /\n0 1 0 /\n1 1 0 /\n)\n( 6/\n{line})\n( 7/\n1 1 3 4 /\n)\n"
        )
    };
    // Nodes 1 (0,0), 2 (1,0), 3 (0,1), 4 (1,1). Tensor order: 1 2 3 4 (an edge 1-3 is the left side);
    // perimeter order: 1 2 4 3 (the left side is the edge 1-4... local nodes 1 and 4).
    let tensor = build("41 1 1 2 3 4 /\n", "1 19 3 1 1 /\n");
    let perimeter = build("41 1 1 2 4 3 /\n", "1 19 3 1 1 /\n");
    // Edge from local node 1 to local node 3 (tensor) = (0,0)-(0,1); in the perimeter list the
    // local node 3 is node 4 = (1,1): the diagonal. The equivalent perimeter row is "1 to 4".
    let perimeter_equivalent = perimeter.replace("( 7/\n1 1 3 4 /", "( 7/\n1 1 4 4 /");
    let (a, ra) = run_with(&tensor, "edge-tensor", None);
    let (b, rb) = run_with(&perimeter_equivalent, "edge-perimeter", None);
    let (fa, fb) = (resultant(&a, 1), resultant(&b, 1));
    // 4 tf/m on 1 m, positive against the axis: 40 kN down.
    assert!((fa - DVec3::new(0., 0., -40.)).length() < 1e-9, "{fa:?}");
    assert!((fb - DVec3::new(0., 0., -40.)).length() < 1e-9, "{fb:?}");
    let (ca, cb) = (case_of(&ra, 1), case_of(&rb, 1));
    assert!(near(ca.exported_moment, ca.source_moment, 1e-6 * ca.moment_scale.max(1.)), "{ca:?}");
    assert!(near(cb.exported_moment, cb.source_moment, 1e-6 * cb.moment_scale.max(1.)), "{cb:?}");
    // Both put the load on the same line (the left side x = 0).
    for list in [&a, &b] {
        let Load::Line { start, end, .. } = &list[0] else { panic!("{list:?}") };
        assert!(start[0].abs() < 1e-9 && end[0].abs() < 1e-9, "{start:?} {end:?}");
    }
}

#[test]
fn a_report_names_lost_force_and_lost_moment() {
    use topo_reconstruct_rs::loads::{compare_resultants, CaseReport, Report, FORCE_TOLERANCE, MOMENT_TOLERANCE};
    let case = |exported: [f64; 3], exported_moment: [f64; 3]| CaseReport {
        case: 1,
        name: "X".into(),
        loads: 1,
        source: [0., 0., -100.],
        exported,
        source_moment: [0., 50., 0.],
        exported_moment,
        moment_scale: 200.,
        not_in_geometry: [0.; 3],
        not_in_geometry_moment: [0.; 3],
        not_in_geometry_abs: 0.,
        not_in_geometry_scale: 0.,
    };
    let report = |c| Report { origin: [0.; 3], cases: vec![c], ..Default::default() };
    assert!(report(case([0., 0., -100.], [0., 50., 0.])).problems(FORCE_TOLERANCE, MOMENT_TOLERANCE).is_empty());
    // Force kept, moment lost (a point load moved to the element centre).
    assert_eq!(report(case([0., 0., -100.], [0., 20., 0.])).problems(FORCE_TOLERANCE, MOMENT_TOLERANCE).len(), 1);
    // Everything exported as zero while the source is not.
    assert_eq!(report(case([0.; 3], [0.; 3])).problems(FORCE_TOLERANCE, MOMENT_TOLERANCE).len(), 1);
    assert!(compare_resultants(([0.; 3], [0.; 3]), ([0.; 3], [0.; 3]), 0., 0.02, 0.05).is_none());
}

#[test]
fn a_couple_on_plates_outside_the_geometry_is_a_problem_although_its_force_sums_to_zero() {
    use topo_reconstruct_rs::loads::{CaseReport, Report, FORCE_TOLERANCE, MOMENT_TOLERANCE};
    let case = CaseReport {
        case: 1,
        name: "PAIR".into(),
        loads: 0,
        source: [0.; 3],
        exported: [0.; 3],
        source_moment: [0.; 3],
        exported_moment: [0.; 3],
        moment_scale: 0.,
        // +240 kN at x = -3 and -240 kN at x = 3: no force, M = 1440 kN m.
        not_in_geometry: [0.; 3],
        not_in_geometry_moment: [0., 1440., 0.],
        not_in_geometry_abs: 480.,
        not_in_geometry_scale: 1440.,
    };
    let report = Report { origin: [0.; 3], cases: vec![case], ..Default::default() };
    assert_eq!(report.problems(FORCE_TOLERANCE, MOMENT_TOLERANCE).len(), 1);
}

/// Two quads tiling 6 x 2 m: the left one is oblique, (0,0) (4,0) (1,2) (0,2).
fn oblique(second_loaded: bool) -> String {
    let rows = if second_loaded { "1 16 3 1 1 /\n2 16 3 1 1 /\n" } else { "1 16 3 1 1 /\n" };
    format!(
        "( 0/ 1; OBLIQUE/ 2; 5/\n39;\n1: P ;\n/\n)\n( 1/\n41 1 1 2 5 4 /\n41 1 2 3 6 5 /\n)\n\
         ( 3/\n1 GEI 0.305915 0.17 0.2 RO 0.254929 /\n)\n( 4/\n0 0 0 /\n4 0 0 /\n6 0 0 /\n0 2 0 /\n1 2 0 /\n6 2 0 /\n)\n\
         ( 6/\n{rows})\n( 7/\n1 2 0 /\n)\n"
    )
}

#[test]
fn the_centre_of_an_oblique_quad_is_its_centre_of_area() {
    // 2 tf/m2 on 5 m2: 100 kN at the centre of area (1.4, 0.8), not at the corner average (1.25, 1).
    let exact = DVec3::new(1.4, 0.8, 0.);
    for combine in [false, true] {
        let text = oblique(false);
        let combination = combine.then(|| loads::Combination { factors: [(1, 1.)].into_iter().collect(), simplify: Some(loads::Simplify::default()) });
        let (list, report) = run_with(&text, if combine { "obl_c" } else { "obl" }, combination);
        let c = &report.cases[0];
        let force = DVec3::new(0., 0., -100.);
        let expected = (exact - DVec3::from_array(report.origin)).cross(force);
        assert!(near(c.source_moment, expected.to_array(), 1e-6), "{combine}: {:?} vs {expected:?}", c.source_moment);
        assert!(near(c.exported_moment, expected.to_array(), 1e-6), "{combine}: {c:?} {list:?}");
        assert!(report.problems(loads::FORCE_TOLERANCE, loads::MOMENT_TOLERANCE).is_empty());
    }
}

#[test]
fn overlapping_line_loads_are_split_and_added_per_case() {
    use loads::{consolidate, Load};
    let line = |case, x0: f64, x1: f64, q: f64| Load::Line { case, start: [x0, 0., 0.], end: [x1, 0., 0.], q_start: [0., 0., q], q_end: [0., 0., q] };
    // Case 1: 0..4 at -10 and 2..6 at -5 (overlap 2..4); case 2: 0..4 at -1, drawn backwards.
    let mut reversed = line(2, 0., 4., -1.);
    if let Load::Line { start, end, .. } = &mut reversed {
        std::mem::swap(start, end);
    }
    let input = vec![line(1, 0., 4., -10.), line(1, 2., 6., -5.), reversed, Load::Point { case: 1, at: [1., 1., 1.], force: [0., 0., -2.], moment: [0.; 3] }, Load::Point { case: 1, at: [1., 1., 1.], force: [0., 0., -3.], moment: [0.; 3] }];
    let force = |list: &[Load], case: u32| resultant(list, case);
    let (f1, f2) = (force(&input, 1), force(&input, 2));
    let out = consolidate(input, 1e-6);
    assert!((force(&out, 1) - f1).length() < 1e-9 && (force(&out, 2) - f2).length() < 1e-9);
    // No two lines of the output overlap: each pair is equal or disjoint along the axis.
    let spans: Vec<(f64, f64, u32)> = out.iter().filter_map(|l| match l { Load::Line { case, start, end, .. } => Some((start[0].min(end[0]), start[0].max(end[0]), *case)), _ => None }).collect();
    for a in &spans {
        for b in &spans {
            let overlap = a.1.min(b.1) - a.0.max(b.0);
            let same = (a.0 - b.0).abs() < 1e-9 && (a.1 - b.1).abs() < 1e-9;
            assert!(overlap <= 1e-9 || same, "{spans:?}");
        }
    }
    // Case 1 has the segments 0..2, 2..4, 4..6; the overlap carries both.
    let mid = out.iter().find_map(|l| match l { Load::Line { case: 1, start, end, q_start, .. } if (start[0] - 2.).abs() < 1e-9 && (end[0] - 4.).abs() < 1e-9 => Some(q_start[2]), _ => None }).unwrap();
    assert!((mid + 15.).abs() < 1e-9, "{mid}");
    // The two equal points are one.
    let points: Vec<&Load> = out.iter().filter(|l| matches!(l, Load::Point { .. })).collect();
    assert_eq!(points.len(), 1);
}

#[test]
fn a_point_load_off_the_structure_moves_onto_it_keeping_force_and_moment() {
    use loads::{attach_points, Load};
    let plate = vec![[0., 0., 0.], [4., 0., 0.], [4., 4., 0.], [0., 4., 0.]];
    let beam = ([6., 0., 0.], [6., 0., 5.]);
    let point = |at: [f64; 3], force: [f64; 3]| Load::Point { case: 1, at, force, moment: [0.; 3] };
    let mut loads = vec![
        point([1., 1., 0.], [0., 0., -10.]),    // on the plate: stays
        point([2., 2., 0.05], [0., 0., -10.]),  // 5 cm above the plate
        point([5., 2., 0.], [0., 0., -10.]),    // beside the plate, off the beam: nearest is the plate edge
        point([6.02, 0., 2.], [3., 0., -10.]),  // 2 cm from the beam
    ];
    let before: Vec<(glam::DVec3, glam::DVec3)> = loads.iter().map(|l| match l {
        Load::Point { at, force, moment, .. } => (glam::DVec3::from_array(*force), glam::DVec3::from_array(*at).cross(glam::DVec3::from_array(*force)) + glam::DVec3::from_array(*moment)),
        _ => unreachable!(),
    }).collect();
    let (moved, farthest) = attach_points(&mut loads, &[plate], &[beam]);
    assert_eq!(moved, 3);
    assert!((farthest - 1.0).abs() < 1e-9, "{farthest}");
    for (l, (force, moment)) in loads.iter().zip(before) {
        let Load::Point { at, force: f, moment: m, .. } = l else { unreachable!() };
        let after = glam::DVec3::from_array(*at).cross(glam::DVec3::from_array(*f)) + glam::DVec3::from_array(*m);
        assert!((glam::DVec3::from_array(*f) - force).length() < 1e-12 && (after - moment).length() < 1e-9, "{l:?}");
    }
    // Everything lies on the structure now.
    let Load::Point { at, .. } = &loads[1] else { unreachable!() };
    assert!(at[2].abs() < 1e-12);
    let Load::Point { at, .. } = &loads[2] else { unreachable!() };
    assert!((at[0] - 4.).abs() < 1e-12);
}
