//! The reconstructed geometry meshed with Gmsh (skipped when the library is
//! not available: set TOPO_GMSH_LIB or put it beside the test executable).
use glam::DVec3;
use topo_reconstruct_rs::editor::Session;
use topo_reconstruct_rs::gmsh::Gmsh;
use topo_reconstruct_rs::meshing::mesh_state;
use topo_reconstruct_rs::pipeline::{self, Options, Profile};

/// A 6 x 4 m slab and a 6 x 3 m wall standing on its edge y = 0; the wall and
/// the slab share the nodes of that edge.
fn source() -> String {
    let (nx, ny, nz) = (6usize, 4usize, 3usize);
    let slab = |i: usize, j: usize| j * (nx + 1) + i + 1;
    let base = (nx + 1) * (ny + 1);
    // Wall nodes: the row k = 0 is the slab's edge j = 0.
    let wall = |i: usize, k: usize| if k == 0 { slab(i, 0) } else { base + (k - 1) * (nx + 1) + i + 1 };
    let mut elements = String::new();
    for j in 0..ny {
        for i in 0..nx {
            elements += &format!("41 1 {} {} {} {} /\n", slab(i, j), slab(i + 1, j), slab(i + 1, j + 1), slab(i, j + 1));
        }
    }
    for k in 0..nz {
        for i in 0..nx {
            elements += &format!("41 1 {} {} {} {} /\n", wall(i, k), wall(i + 1, k), wall(i + 1, k + 1), wall(i, k + 1));
        }
    }
    let mut nodes = String::new();
    for j in 0..=ny {
        for i in 0..=nx {
            nodes += &format!("{i}.0 {j}.0 0.0 /\n");
        }
    }
    for k in 1..=nz {
        for i in 0..=nx {
            nodes += &format!("{i}.0 0.0 {k}.0 /\n");
        }
    }
    format!("( 0/ 1; WALL/ 2; 5/\n39;\n1: C ;\n/\n)\n( 1/\n{elements})\n( 3/\n1 GEI 0.305915 0.17 0.2 RO 0.254929 /\n)\n( 4/\n{nodes})\n")
}

#[test]
fn slab_and_wall_get_one_conforming_mesh() {
    let gmsh = match Gmsh::load() {
        Ok(g) => g,
        Err(e) => {
            eprintln!("Gmsh library not available ({e}): test skipped");
            return;
        }
    };
    let path = std::env::temp_dir().join(format!("gmsh_mesh_{}.txt", std::process::id()));
    std::fs::write(&path, source()).unwrap();
    let profile = Profile::plaxis();
    let output = pipeline::run(&path, &profile, &Options::default(), &mut |_| {}).unwrap();
    let _ = std::fs::remove_file(&path);
    let session = Session::new(&output.topology, profile.audit_options());
    assert_eq!(session.state().model.surfaces().len(), 2);
    let mesh = mesh_state(&gmsh, session.state(), 0.5, false).unwrap();
    // Area of the shells is the area of the two panels.
    let area: f64 = mesh
        .shells
        .iter()
        .map(|s| {
            let p: Vec<DVec3> = s.nodes.iter().map(|&n| DVec3::from_array(mesh.nodes[n])).collect();
            (1..p.len() - 1).map(|i| (p[i] - p[0]).cross(p[i + 1] - p[0]).length() / 2.).sum::<f64>()
        })
        .sum();
    assert!((area - (24. + 18.)).abs() < 1e-6, "{area}");
    // No two nodes at one place, none unused: the two panels share the nodes of their common edge.
    let mut seen = std::collections::HashSet::new();
    for p in &mesh.nodes {
        assert!(seen.insert([0, 1, 2].map(|k| (p[k] * 1e5).round() as i64)), "duplicate node {p:?}");
    }
    let edge_nodes = mesh
        .nodes
        .iter()
        .enumerate()
        .filter(|(_, p)| p[1].abs() < 1e-9 && p[2].abs() < 1e-9)
        .map(|(i, _)| i)
        .collect::<Vec<_>>();
    assert!(edge_nodes.len() >= 7);
    let users = |surface: usize| -> std::collections::HashSet<usize> {
        mesh.shells.iter().filter(|s| s.surface == surface).flat_map(|s| s.nodes.iter().copied()).collect()
    };
    assert!(edge_nodes.iter().all(|n| users(0).contains(n) && users(1).contains(n)));
}

#[test]
fn a_pressure_contour_through_elements_keeps_its_resultant_on_the_mesh() {
    use topo_reconstruct_rs::loads::Load;
    use topo_reconstruct_rs::mesh_loads;
    let gmsh = match Gmsh::load() {
        Ok(g) => g,
        Err(e) => {
            eprintln!("Gmsh library not available ({e}): test skipped");
            return;
        }
    };
    let path = std::env::temp_dir().join(format!("gmsh_mesh_loads_{}.txt", std::process::id()));
    std::fs::write(&path, source()).unwrap();
    let profile = Profile::plaxis();
    let output = pipeline::run(&path, &profile, &Options::default(), &mut |_| {}).unwrap();
    let _ = std::fs::remove_file(&path);
    let session = Session::new(&output.topology, profile.audit_options());
    let state = session.state();
    // Element size 0.7 does not divide the contour (x from 0 to 2.5): elements are cut by it.
    let mesh = mesh_state(&gmsh, state, 0.7, false).unwrap();
    let slab = (0..state.model.surfaces().len())
        .find(|&s| state.model.planes()[state.model.surfaces()[s].plane].normal()[2].abs() > 0.9)
        .unwrap();
    let load = Load::Surface {
        case: 1,
        surface: slab,
        polygons: vec![vec![[0., 0., 0.], [2.5, 0., 0.], [2.5, 4., 0.], [0., 4., 0.]]],
        sigma: [0., 0., -3.],
    };
    let on_mesh = mesh_loads::transfer(&mesh, state, &[load], 0.02);
    let resultant = on_mesh.resultants(&mesh)[&1];
    let expected = -3. * 2.5 * 4.;
    assert!((resultant.z - expected).abs() < 1e-6, "{resultant:?} vs {expected}");
}
