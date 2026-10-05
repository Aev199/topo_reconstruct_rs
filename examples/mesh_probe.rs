//! Mesh a reconstructed LIRA model with Gmsh and print what came out
//! (development). Usage: mesh_probe MODEL.txt [--frame-cache PATH] [--size 0.5] [--quads]
//! The Gmsh library: TOPO_GMSH_LIB or beside the executable.
use topo_reconstruct_rs::editor::Session;
use topo_reconstruct_rs::pipeline::{self, Options, Profile};

fn main() -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let mut args = std::env::args().skip(1);
    let input = std::path::PathBuf::from(args.next().ok_or("usage: mesh_probe MODEL.txt")?);
    let (mut cache, mut size, mut quads) = (None, 0.5, false);
    while let Some(a) = args.next() {
        match a.as_str() {
            "--frame-cache" => cache = args.next().map(std::path::PathBuf::from),
            "--size" => size = args.next().and_then(|s| s.parse().ok()).unwrap_or(size),
            "--quads" => quads = true,
            _ => {}
        }
    }
    let profile = Profile::plaxis();
    let output = pipeline::run(&input, &profile, &Options { mesh: false, frame_cache: cache }, &mut |_| {})?;
    let session = Session::new(&output.topology, profile.audit_options());
    let gmsh = topo_reconstruct_rs::gmsh::Gmsh::load()?;
    let started = std::time::Instant::now();
    let mesh = topo_reconstruct_rs::meshing::mesh_state(&gmsh, session.state(), size, quads)?;
    let tri = mesh.shells.iter().filter(|s| s.nodes.len() == 3).count();
    println!(
        "{} nodes, {} shells ({} triangles, {} quads), {} bar elements, {:?}",
        mesh.nodes.len(), mesh.shells.len(), tri, mesh.shells.len() - tri, mesh.bars.len(), started.elapsed()
    );
    // Conformity: no two nodes in one place; every bar node is a shell node
    // or a bar-only node; free shell edges per surface.
    let mut grid = std::collections::HashMap::new();
    let mut duplicates = 0;
    let vertex_node: std::collections::HashSet<usize> = mesh.vertex_nodes.iter().flatten().copied().collect();
    let (mut vv, mut vm, mut mm) = (0, 0, 0);
    for (i, p) in mesh.nodes.iter().enumerate() {
        let key = [0, 1, 2].map(|k| (p[k] * 1e5).round() as i64);
        if let Some(j) = grid.insert(key, i) {
            duplicates += 1;
            match (vertex_node.contains(&i), vertex_node.contains(&j)) {
                (true, true) => vv += 1,
                (false, false) => mm += 1,
                _ => vm += 1,
            }
        }
    }
    let mut users: std::collections::HashMap<usize, std::collections::BTreeSet<usize>> = Default::default();
    for sh in &mesh.shells {
        for &n in &sh.nodes {
            users.entry(n).or_default().insert(sh.surface);
        }
    }
    let mut shown = 0;
    let mut pairs: std::collections::BTreeMap<(Vec<usize>, Vec<usize>), usize> = Default::default();
    {
        let mut g2 = std::collections::HashMap::new();
        for (i, p) in mesh.nodes.iter().enumerate() {
            let key = [0, 1, 2].map(|k| (p[k] * 1e5).round() as i64);
            if let Some(j) = g2.insert(key, i) {
                let a: Vec<usize> = users.get(&j).map(|s| s.iter().copied().collect()).unwrap_or_default();
                let b: Vec<usize> = users.get(&i).map(|s| s.iter().copied().collect()).unwrap_or_default();
                if shown < 4 { shown += 1; println!("dup at {:?}: surfaces {a:?} vs {b:?}", p); }
                *pairs.entry((a, b)).or_default() += 1;
            }
        }
    }
    let mut top: Vec<_> = pairs.into_iter().collect();
    top.sort_by_key(|x| std::cmp::Reverse(x.1));
    {
        let model = &session.state().model;
        let target = glam::DVec3::new(30.1, 40.47, 28.839);
        for (e, [a, b]) in model.edges().iter().enumerate() {
            let (p, q) = (glam::DVec3::from_array(model.vertices()[*a]), glam::DVec3::from_array(model.vertices()[*b]));
            let d = q - p;
            let t = ((target - p).dot(d) / d.length_squared()).clamp(0., 1.);
            if (p + d * t).distance(target) < 1e-3 {
                let surfaces: Vec<usize> = (0..model.surfaces().len()).filter(|&s| model.surface_edges(s).any(|x| x == e)).collect();
                let bars: Vec<usize> = (0..session.state().axes.len()).filter(|&i| session.state().axes[i].anchors.windows(2).any(|w| (w[0].vertex == *a && w[1].vertex == *b) || (w[0].vertex == *b && w[1].vertex == *a))).collect();
                println!("edge {e}: {a} {:?} - {b} {:?} surfaces {surfaces:?} bars {bars:?}", model.vertices()[*a], model.vertices()[*b]);
            }
        }
    }
    {
        let model = &session.state().model;
        let used: Vec<usize> = (0..model.edges().len()).filter(|&e| (0..model.surfaces().len()).any(|s| model.surface_edges(s).any(|x| x == e))).collect();
        let pts = |e: usize| { let [a, b] = model.edges()[e]; (glam::DVec3::from_array(model.vertices()[a]), glam::DVec3::from_array(model.vertices()[b])) };
        let mut overlaps = 0;
        let mut shown = 0;
        for (i, &e) in used.iter().enumerate() {
            let (p, q) = pts(e);
            let d = q - p;
            let len = d.length();
            for &f in &used[i + 1..] {
                let (r, t) = pts(f);
                let (lo, hi) = (p.min(q) - glam::DVec3::splat(1e-6), p.max(q) + glam::DVec3::splat(1e-6));
                if r.max(t).cmplt(lo).any() || r.min(t).cmpgt(hi).any() { continue; }
                // collinear and overlapping
                let u = d / len;
                let dist = |x: glam::DVec3| (x - p - u * (x - p).dot(u)).length();
                if dist(r) < 1e-6 && dist(t) < 1e-6 {
                    let (s0, s1) = ((r - p).dot(u), (t - p).dot(u));
                    let overlap = s0.max(s1).min(len) - s0.min(s1).max(0.);
                    if overlap > 1e-6 {
                        overlaps += 1;
                        if shown < 5 { shown += 1; println!("overlapping edges {e} {:?}-{:?} and {f} {:?}-{:?}: {overlap:.3} m", model.edges()[e][0], model.edges()[e][1], model.edges()[f][0], model.edges()[f][1]); }
                    }
                }
            }
        }
        println!("overlapping model edge pairs: {overlaps}");
    }
    {
        let mut referenced: std::collections::HashSet<usize> = shell_nodes_set(&mesh);
        referenced.extend(mesh.bars.iter().flat_map(|b| b.nodes));
        let orphan = (0..mesh.nodes.len()).filter(|n| !referenced.contains(n)).count();
        println!("nodes in no element: {orphan}");
    }
    println!("most common duplicate pairs: {:?}", &top[..top.len().min(8)]);
    println!("duplicates: vertex-vertex {vv}, vertex-other {vm}, other-other {mm}");
    let mut edge_use: std::collections::HashMap<(usize, usize), usize> = Default::default();
    for s in &mesh.shells {
        for i in 0..s.nodes.len() {
            let (a, b) = (s.nodes[i], s.nodes[(i + 1) % s.nodes.len()]);
            *edge_use.entry((a.min(b), a.max(b))).or_default() += 1;
        }
    }
    let free = edge_use.values().filter(|&&n| n == 1).count();
    let multi = edge_use.values().filter(|&&n| n > 2).count();
    let shell_nodes: std::collections::HashSet<usize> = mesh.shells.iter().flat_map(|s| s.nodes.iter().copied()).collect();
    let bar_free = mesh.bars.iter().flat_map(|b| b.nodes).filter(|n| !shell_nodes.contains(n)).collect::<std::collections::HashSet<_>>().len();
    println!("duplicate node positions {duplicates}, free shell edges {free}, edges of 3+ shells {multi}, bar-only nodes {bar_free}");
    Ok(())
}

fn shell_nodes_set(mesh: &topo_reconstruct_rs::meshing::Mesh) -> std::collections::HashSet<usize> {
    mesh.shells.iter().flat_map(|s| s.nodes.iter().copied()).collect()
}
