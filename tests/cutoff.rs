//! Cutting a tower off at a floor: surfaces crossing the level are clipped,
//! what lies above is removed, the rest stays a valid geometry.
use topo_reconstruct_rs::editor::{Edit, Session};
use topo_reconstruct_rs::pipeline::{self, Options, Profile};
use topo_reconstruct_rs::reconstruction::assembly::cutoff::floors;

/// Three 6 x 4 m slabs at z = 0, 3, 6, one 6 m tall wall on y = 0 through all
/// of them (nodes shared with the slab edges) and a column through the
/// slabs at (3, 2).
fn tower_with(loads: &str) -> String {
    let (nx, ny) = (6usize, 4usize);
    let per_slab = (nx + 1) * (ny + 1);
    let slab_node = |k: usize, i: usize, j: usize| k * per_slab + j * (nx + 1) + i + 1;
    let mut nodes = String::new();
    for k in 0..3 {
        for j in 0..=ny {
            for i in 0..=nx {
                nodes += &format!("{i}.0 {j}.0 {}.0 /\n", k * 3);
            }
        }
    }
    // Wall nodes between the slabs (z = 1, 2, 4, 5), y = 0.
    let base = 3 * per_slab;
    let wall = |i: usize, z: usize| -> usize {
        match z {
            0 => slab_node(0, i, 0),
            3 => slab_node(1, i, 0),
            6 => slab_node(2, i, 0),
            _ => base + (if z < 3 { z - 1 } else { z - 2 }) * (nx + 1) + i + 1,
        }
    };
    for z in [1, 2, 4, 5] {
        for i in 0..=nx {
            nodes += &format!("{i}.0 0.0 {z}.0 /\n");
        }
    }
    let mut elements = String::new();
    for k in 0..3 {
        for j in 0..ny {
            for i in 0..nx {
                elements += &format!(
                    "41 1 {} {} {} {} /\n",
                    slab_node(k, i, j), slab_node(k, i + 1, j), slab_node(k, i + 1, j + 1), slab_node(k, i, j + 1)
                );
            }
        }
    }
    for z in 0..6 {
        for i in 0..nx {
            elements += &format!("41 1 {} {} {} {} /\n", wall(i, z), wall(i + 1, z), wall(i + 1, z + 1), wall(i, z + 1));
        }
    }
    elements += &format!("10 2 {} {} /\n10 2 {} {} /\n", slab_node(0, 3, 2), slab_node(1, 3, 2), slab_node(1, 3, 2), slab_node(2, 3, 2));
    format!(
        "( 0/ 1; TOWER/ 2; 5/\n39;\n1: FLOOR ;\n/\n)\n( 1/\n{elements})\n( 3/\n1 GEI 0.305915 0.17 0.2 RO 0.254929 /\n2 S0 0.305915 20 20/\n 0 RO 0.0101972/\n)\n( 4/\n{nodes})\n{loads}"
    )
}

fn tower() -> String {
    tower_with("")
}

/// Loads: 2 tf/m2 on the top slab (elements 49-72, case 1); 1 tf/m2 along global Y
/// on the upper half of the wall (elements 91-108, case 2).
fn loaded_tower() -> String {
    let mut rows = String::new();
    for e in 49..=72 {
        rows += &format!("{e} 16 3 1 1 /\n");
    }
    for e in 91..=108 {
        rows += &format!("{e} 16 2 2 2 /\n");
    }
    tower_with(&format!("( 6/\n{rows})\n( 7/\n1 2 0 / 2 1 0 /\n)\n"))
}

fn session() -> Session {
    session_of(&tower())
}

fn session_of(text: &str) -> Session {
    static NEXT: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    let path = std::env::temp_dir().join(format!("cutoff_{}_{}.txt", std::process::id(), NEXT.fetch_add(1, std::sync::atomic::Ordering::SeqCst)));
    std::fs::write(&path, text).unwrap();
    let profile = Profile::plaxis();
    let output = pipeline::run(&path, &profile, &Options::default(), &mut |_| {}).unwrap();
    let _ = std::fs::remove_file(&path);
    Session::new(&output.topology, profile.audit_options())
}

#[test]
fn floors_are_recognized() {
    let s = session();
    let f = floors(s.state());
    assert_eq!(f.len(), 3, "{f:?}");
    assert!((f[1].z - 3.).abs() < 1e-9 && f[1].slabs == 1 && (f[1].area - 24.).abs() < 1e-6, "{f:?}");
}

#[test]
fn the_cut_keeps_a_valid_lower_part_and_what_the_upper_part_rested_on() {
    let mut s = session();
    let before = s.state().model.surfaces().len();
    assert!(before >= 4, "{before}");
    s.apply(Edit::CutAbove { z: 3. }, "cut").unwrap_or_else(|e| panic!("{e}"));
    let state = s.state();
    let model = &state.model;
    // Nothing above the level, the slab at the level and the lower wall remain.
    assert!(model.surfaces().iter().enumerate().all(|(i, _)| model.surface_edges(i).flat_map(|e| model.edges()[e]).all(|v| model.vertices()[v][2] <= 3. + 1e-9)));
    let area: f64 = model.surfaces().iter().map(|s| {
        let a = |r: &Vec<[f64; 2]>| (0..r.len()).map(|i| r[i][0] * r[(i + 1) % r.len()][1] - r[(i + 1) % r.len()][0] * r[i][1]).sum::<f64>().abs() / 2.;
        a(&s.contours[0]) - s.contours[1..].iter().map(a).sum::<f64>()
    }).sum();
    // Two slabs of 24 and the wall of 6 x 3.
    assert!((area - (48. + 18.)).abs() < 1e-6, "{area}");
    let (walls, columns) = topo_reconstruct_rs::reconstruction::assembly::cutoff::live_supports(state);
    assert_eq!(walls.len(), 1, "{:?}", state.cut);
    let (a, b, _) = walls[0];
    assert!((a.distance(b) - 6.).abs() < 1e-6 && (a.z - 3.).abs() < 1e-9, "{a:?} {b:?}");
    assert_eq!(columns.len(), 1, "{:?}", state.cut);
    assert!((columns[0].0.x - 3.).abs() < 1e-9 && (columns[0].0.z - 3.).abs() < 1e-9);
    // The column below is whole, the one above is gone.
    assert_eq!(state.axes.len(), 1);
    // The geometry is still valid.
    assert!(s.audit().findings.iter().all(|f| f.class != topo_reconstruct_rs::audit::Class::Failure), "{:?}", s.audit().counts);
    // Undo restores the whole model.
    s.undo().unwrap();
    assert_eq!(s.state().model.surfaces().len(), before);
    assert!(s.state().cut.is_none());
}

#[test]
fn what_was_cut_off_rests_on_the_supports_with_force_and_overturning_moment() {
    use glam::DVec3;
    use topo_reconstruct_rs::loads::{self, Load, CUT_WEIGHT_CASE};
    let text = loaded_tower();
    let mut s = session_of(&text);
    s.apply(Edit::CutAbove { z: 3. }, "cut").unwrap();
    let bytes = text.as_bytes();
    let mesh = topo_reconstruct_rs::parsers::lira::LiraParser::mesh_from(bytes).unwrap();
    let set = topo_reconstruct_rs::parsers::loads::parse(bytes);
    let path = std::env::temp_dir().join(format!("cutoff_materials_{}.txt", std::process::id()));
    std::fs::write(&path, bytes).unwrap();
    let materials = std::sync::Arc::new(topo_reconstruct_rs::parsers::lira::LiraParser::parse_materials(&path).unwrap());
    let _ = std::fs::remove_file(&path);
    let nodes = vec![]; // no vertex of the model carries a node load here
    let (list, report) = loads::transfer(
        s.state(),
        &nodes,
        &mesh,
        &set,
        loads::Settings { force_factor: 10., snap: 0.05, max_groups: 40, combination: None, cases: None, materials: Some(materials.clone()), cut_loads: true },
    );
    assert!(!set.rows.is_empty(), "no load rows parsed");
    assert!(s.state().cut.is_some());
    let total = |case: u32| -> DVec3 {
        list.iter().filter(|l| l.case() == case).map(|l| match l {
            Load::Point { force, .. } => DVec3::from_array(*force),
            Load::Line { start, end, q_start, q_end, .. } => (DVec3::from_array(*q_start) + DVec3::from_array(*q_end)) / 2. * DVec3::from_array(*start).distance(DVec3::from_array(*end)),
            Load::Surface { .. } => DVec3::ZERO,
        }).sum()
    };
    // The floor load of the removed slab: 2 tf/m2 x 24 m2 downwards, x 10 kN/tf.
    let floor = total(1);
    assert!((floor - DVec3::new(0., 0., -480.)).length() < 1e-6, "{floor:?}");
    // The wind: 1 tf/m2 x 18 m2 (positive acts against the axis).
    let wind = total(2);
    assert!((wind - DVec3::new(0., -180., 0.)).length() < 1e-6, "{wind:?}");
    // Its moment about (3, 2, 3) is kept: the true one is (0, -2, 1.5) x (0, -180, 0) = (270, 0, 0) kN m.
    let about = DVec3::new(3., 2., 3.);
    let moment: DVec3 = list.iter().filter(|l| l.case() == 2).map(|l| match l {
        Load::Point { at, force, .. } => (DVec3::from_array(*at) - about).cross(DVec3::from_array(*force)),
        Load::Line { start, end, q_start, q_end, .. } => {
            let (a, b) = (DVec3::from_array(*start), DVec3::from_array(*end));
            (((a + b) / 2.) - about).cross((DVec3::from_array(*q_start) + DVec3::from_array(*q_end)) / 2. * a.distance(b))
        }
        Load::Surface { .. } => DVec3::ZERO,
    }).sum();
    assert!((moment.x - 270.).abs() < 1e-6 && moment.y.abs() < 1e-6, "{moment:?}");
    // The weight of the removed slab, wall panel and column piece.
    let weight = total(CUT_WEIGHT_CASE);
    assert!(weight.z < 0. && weight.x.abs() < 1e-9 && weight.y.abs() < 1e-9, "{weight:?}");
    assert!(report.cases.iter().any(|c| c.case == CUT_WEIGHT_CASE));
    // Without the option nothing of it is transferred.
    let (list, _) = loads::transfer(
        s.state(),
        &nodes,
        &mesh,
        &set,
        loads::Settings { force_factor: 10., snap: 0.05, max_groups: 40, combination: None, cases: None, materials: Some(materials), cut_loads: false },
    );
    assert!(list.iter().all(|l| l.case() != CUT_WEIGHT_CASE && !matches!(l, Load::Line { .. })), "{list:?}");
}

#[test]
fn the_cap_slab_gets_the_bending_stiffness_of_what_stood_on_it() {
    use topo_reconstruct_rs::parsers::lira::Material;
    use topo_reconstruct_rs::storeys::{with_cap, CAP_STIFFNESS_BASE};
    let text = tower();
    let mut s = session_of(&text);
    let path = std::env::temp_dir().join(format!("cutoff_cap_{}.txt", std::process::id()));
    std::fs::write(&path, &text).unwrap();
    let materials = topo_reconstruct_rs::parsers::lira::LiraParser::parse_materials(&path).unwrap();
    let _ = std::fs::remove_file(&path);
    // No cut, no cap.
    assert!(with_cap(s.state(), &materials, 1.).is_none());
    s.apply(Edit::CutAbove { z: 3. }, "cut").unwrap();
    let (state, new_materials, report) = with_cap(s.state(), &materials, 1.).expect("a cap");
    assert_eq!(report.surfaces, 1);
    assert!(report.equivalent_thickness > 0.2, "{report:?}");
    // One storey of 3 m below the level, 3 m cut off: the storey factor is 1.
    assert!((report.storey - 3.).abs() < 1e-9 && (report.removed_height - 3.).abs() < 1e-9, "{report:?}");
    let cap = (0..state.model.surfaces().len()).find(|&i| state.stiffness[i] >= CAP_STIFFNESS_BASE).unwrap();
    let Some(Material::Plate { bending: Some(kb), membrane, thickness, .. }) = new_materials.get(&state.stiffness[cap]) else { panic!() };
    // Only the bending stiffness grew: PLKE = (t_eq / t)^3.
    assert!((kb - (report.equivalent_thickness / thickness).powi(3)).abs() < 1e-9 * kb, "{kb}");
    assert!(membrane.map_or(true, |m| (m - 1.).abs() < 1e-12));
    // A larger factor gives a thicker plate; zero leaves the slab as it was.
    let (_, _, doubled) = with_cap(s.state(), &materials, 8.).unwrap();
    assert!((doubled.equivalent_thickness / report.equivalent_thickness - 2.).abs() < 1e-6);
    assert!(with_cap(s.state(), &materials, 0.).is_none());
}

#[test]
fn a_cut_leaves_no_coincident_vertices() {
    // The cut level is shared by the slab, the wall panels and the column.
    let mut s = session();
    for z in [3., 4.5] {
        s.apply(Edit::CutAbove { z }, "cut").unwrap();
        let v = s.state().model.vertices().to_vec();
        for i in 0..v.len() {
            for j in 0..i {
                let d = (0..3).map(|k| (v[i][k] - v[j][k]).powi(2)).sum::<f64>().sqrt();
                assert!(d > 1e-6, "vertices {j} and {i} coincide at {:?} (cut at {z})", v[i]);
            }
        }
        s.undo().unwrap();
    }
}

/// The slab at z = 0 (surface index) of the tower, with a one-element mesh of its corners.
fn slab_mesh(s: &Session) -> (topo_reconstruct_rs::meshing::Mesh, usize) {
    use topo_reconstruct_rs::meshing::{Mesh, Shell};
    let model = &s.state().model;
    let surface = (0..model.surfaces().len())
        .find(|&i| model.surface_edges(i).flat_map(|e| model.edges()[e]).all(|v| model.vertices()[v][2].abs() < 1e-9))
        .expect("the slab at z = 0");
    let mesh = Mesh {
        nodes: vec![[0., 0., 0.], [6., 0., 0.], [6., 4., 0.], [0., 4., 0.]],
        vertex_nodes: vec![],
        shells: vec![Shell { nodes: vec![0, 1, 2, 3], surface, stiffness: 1 }],
        bars: vec![],
    };
    (mesh, surface)
}

#[test]
fn a_patch_inside_one_big_element_keeps_force_and_moment_on_the_mesh() {
    use glam::DVec3;
    use topo_reconstruct_rs::loads::Load;
    use topo_reconstruct_rs::mesh_loads;
    let s = session();
    let (mesh, surface) = slab_mesh(&s);
    // 1 m2 at (5, 1) under 10 kN/m2, and a point load at (4.5, 1).
    let loads = vec![
        Load::Surface { case: 1, surface, polygons: vec![vec![[4.5, 0.5, 0.], [5.5, 0.5, 0.], [5.5, 1.5, 0.], [4.5, 1.5, 0.]]], sigma: [0., 0., -10.] },
        Load::Point { case: 2, at: [4.5, 1., 0.], force: [0., 0., -10.], moment: [0.; 3] },
    ];
    let on_mesh = mesh_loads::transfer(&mesh, s.state(), &loads, 0.02);
    assert!(on_mesh.lost.is_empty(), "{:?}", on_mesh.lost);
    let about = on_mesh.resultants_about(&mesh, DVec3::ZERO);
    let (force, moment, scale) = about[&1];
    assert!((force - DVec3::new(0., 0., -10.)).length() < 1e-6, "{force:?}");
    assert!((moment - DVec3::new(-10., 50., 0.)).length() < 0.05 * scale, "{moment:?} {scale}");
    // Barycentric distribution of a point load keeps its moment.
    let (force, moment, _) = about[&2];
    assert!((force - DVec3::new(0., 0., -10.)).length() < 1e-6, "{force:?}");
    assert!((moment - DVec3::new(-10., 45., 0.)).length() < 1e-6, "{moment:?}");
}

#[test]
fn a_partial_bar_load_is_not_spread_over_the_whole_bar() {
    use glam::DVec3;
    use topo_reconstruct_rs::loads::Load;
    use topo_reconstruct_rs::meshing::{BarPiece, Mesh};
    use topo_reconstruct_rs::mesh_loads;
    let s = session();
    let state = s.state();
    // The column through (3, 2): one bar piece over the whole axis.
    let axis = (0..state.axes.len())
        .find(|&i| {
            let [a, b] = state.axes[i].endpoints.map(|v| state.model.vertices()[v]);
            (a[0] - 3.).abs() < 1e-9 && (b[0] - 3.).abs() < 1e-9
        })
        .expect("the column");
    let [a, b] = state.axes[axis].endpoints.map(|v| state.model.vertices()[v]);
    let mesh = Mesh {
        nodes: vec![a, b],
        vertex_nodes: vec![],
        shells: vec![],
        bars: vec![BarPiece { nodes: [0, 1], axis, stiffness: 2, t: [0., 1.] }],
    };
    let length = DVec3::from_array(a).distance(DVec3::from_array(b));
    let (z0, z1) = (a[2] + 0.25 * (b[2] - a[2]), a[2] + 0.75 * (b[2] - a[2]));
    let loads = vec![Load::Line { case: 1, start: [3., 2., z0], end: [3., 2., z1], q_start: [0., 0., -10.], q_end: [0., 0., -10.] }];
    let on_mesh = mesh_loads::transfer(&mesh, state, &loads, 0.02);
    let (force, moment, scale) = on_mesh.resultants_about(&mesh, DVec3::new(3., 2., a[2].min(b[2])))[&1];
    // Half the length under 10 kN/m.
    assert!((force - DVec3::new(0., 0., -5. * length)).length() < 1e-6, "{force:?}");
    assert!(moment.length() < 1e-6 * scale.max(1.), "{moment:?}");
}

#[test]
fn supports_follow_the_live_geometry_after_edits() {
    use topo_reconstruct_rs::reconstruction::assembly::cutoff::live_supports;
    let mut s = session();
    s.apply(Edit::CutAbove { z: 3. }, "cut").unwrap();
    assert_eq!(live_supports(s.state()).1.len(), 1);
    // The column is deleted: nothing rests on a point any more.
    let bars: Vec<usize> = (0..s.state().axes.len()).collect();
    s.apply(Edit::DeleteBars { bars }, "delete").unwrap();
    assert!(live_supports(s.state()).1.is_empty(), "{:?}", live_supports(s.state()));
    // The wall is still there.
    assert_eq!(live_supports(s.state()).0.len(), 1);
}

#[test]
fn a_deleted_wall_is_no_support_even_though_its_top_edge_stays_on_the_slab() {
    use topo_reconstruct_rs::reconstruction::assembly::cutoff::live_supports;
    let mut s = session();
    s.apply(Edit::CutAbove { z: 3. }, "cut").unwrap();
    assert_eq!(live_supports(s.state()).0.len(), 1);
    let model = &s.state().model;
    let wall = (0..model.surfaces().len())
        .find(|&i| model.planes()[model.surfaces()[i].plane].normal()[2].abs() < 0.5)
        .expect("the wall");
    s.apply(Edit::DeleteSurface { surface: wall }, "delete").unwrap();
    assert!(live_supports(s.state()).0.is_empty(), "{:?}", live_supports(s.state()).0);
}

fn mesh_resultant(loads: Vec<topo_reconstruct_rs::loads::Load>) -> (glam::DVec3, glam::DVec3, f64, topo_reconstruct_rs::mesh_loads::MeshLoads) {
    use glam::DVec3;
    use topo_reconstruct_rs::mesh_loads;
    let s = session();
    let (mesh, _) = slab_mesh(&s);
    let on_mesh = mesh_loads::transfer(&mesh, s.state(), &loads, 0.02);
    let (f, m, scale) = on_mesh.resultants_about(&mesh, DVec3::ZERO).values().next().copied().unwrap_or_default();
    (f, m, scale, on_mesh)
}

#[test]
fn a_load_over_part_of_a_plate_edge_keeps_force_and_moment() {
    use glam::DVec3;
    use topo_reconstruct_rs::loads::Load;
    // Edge y = 0 of the one 6 x 4 element: 10 kN/m over x = 1.5 .. 4.5 only.
    let (f, m, scale, on_mesh) = mesh_resultant(vec![Load::Line { case: 1, start: [1.5, 0., 0.], end: [4.5, 0., 0.], q_start: [0., 0., -10.], q_end: [0., 0., -10.] }]);
    assert!(on_mesh.lost.is_empty(), "{:?}", on_mesh.lost);
    assert!((f - DVec3::new(0., 0., -30.)).length() < 1e-9, "{f:?}");
    assert!((m - DVec3::new(0., 90., 0.)).length() < 1e-9 * scale.max(1.), "{m:?}");
}

#[test]
fn a_triangular_edge_load_keeps_its_moment() {
    use glam::DVec3;
    use topo_reconstruct_rs::loads::Load;
    // q = -2 x along the whole edge: F = -36, M_y = 144.
    let (f, m, _, _) = mesh_resultant(vec![Load::Line { case: 1, start: [0., 0., 0.], end: [6., 0., 0.], q_start: [0., 0., 0.], q_end: [0., 0., -12.] }]);
    assert!((f - DVec3::new(0., 0., -36.)).length() < 1e-9, "{f:?}");
    assert!((m - DVec3::new(0., 144., 0.)).length() < 1e-9, "{m:?}");
}

#[test]
fn a_concave_contour_does_not_cover_what_its_notch_leaves_free() {
    use glam::DVec3;
    use topo_reconstruct_rs::loads::Load;
    let s = session();
    let (_, surface) = slab_mesh(&s);
    // The slab (6 x 4) under 10 kN/m2 except a notch x 2.5 .. 3.5, y 2.5 .. 4 (1.5 m2) entering from the edge.
    let ring = [[-1., -1., 0.], [7., -1., 0.], [7., 5., 0.], [3.5, 5., 0.], [3.5, 2.5, 0.], [2.5, 2.5, 0.], [2.5, 5., 0.], [-1., 5., 0.]];
    let (f, m, scale, _) = mesh_resultant(vec![Load::Surface { case: 1, surface, polygons: vec![ring.to_vec()], sigma: [0., 0., -10.] }]);
    assert!((f - DVec3::new(0., 0., -225.)).length() < 1e-6, "{f:?}");
    // First moments: (67.5, 43.125) m3.
    assert!((m - DVec3::new(-431.25, 675., 0.)).length() < 1e-6 * scale.max(1.), "{m:?}");
}

#[test]
fn a_patch_a_little_off_the_centre_still_keeps_its_moment() {
    use glam::DVec3;
    use topo_reconstruct_rs::loads::Load;
    let s = session();
    let (_, surface) = slab_mesh(&s);
    // x = 0 .. 5.7 of the 6 x 4 element: the centroid is 0.15 m off its centre.
    let ring = vec![[0., 0., 0.], [5.7, 0., 0.], [5.7, 4., 0.], [0., 4., 0.]];
    let (f, m, _, _) = mesh_resultant(vec![Load::Surface { case: 1, surface, polygons: vec![ring], sigma: [0., 0., -10.] }]);
    let area = 5.7 * 4.;
    assert!((f.z + 10. * area).abs() < 1e-5, "{f:?}");
    // M = r x F: M_x = y F = 2 (-10 A), M_y = -x F = 2.85 (10 A).
    assert!((m - DVec3::new(-20. * area, 28.5 * area, 0.)).length() < 1e-5, "{m:?}");
}
