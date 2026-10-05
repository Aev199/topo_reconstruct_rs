//! The reconstructed geometry meshed with Gmsh: shell elements on the
//! surfaces (conforming along shared edges), line elements on the bars.
//! Every model vertex used by a surface or a bar is a node of the mesh, so
//! the topology of the reconstruction is the topology of the mesh.
use crate::gmsh::{Error, Gmsh, Input, SurfaceInput};
use crate::reconstruction::assembly::bars::Contact;
use crate::reconstruction::assembly::edit::State;
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Serialize)]
pub struct Shell {
    /// 3 or 4 nodes.
    pub nodes: Vec<usize>,
    pub surface: usize,
    pub stiffness: u32,
}

#[derive(Debug, Clone, Serialize)]
pub struct BarPiece {
    pub nodes: [usize; 2],
    pub axis: usize,
    pub stiffness: u32,
    /// Position along the axis of its two nodes.
    pub t: [f64; 2],
}

#[derive(Debug, Clone, Default, Serialize)]
pub struct Mesh {
    pub nodes: Vec<[f64; 3]>,
    /// Node of a model vertex (`None`: the vertex is in no surface or bar).
    pub vertex_nodes: Vec<Option<usize>>,
    pub shells: Vec<Shell>,
    pub bars: Vec<BarPiece>,
}

/// Mesh the model with elements of about `size` metres.
pub fn mesh_state(gmsh: &Gmsh, state: &State, size: f64, quads: bool) -> Result<Mesh, Error> {
    let model = &state.model;
    // Lines: the model edges in use, then the bar pieces that are no edge.
    let mut line_of_edge: BTreeMap<usize, usize> = BTreeMap::new();
    let mut lines: Vec<[usize; 2]> = vec![];
    let mut use_edge = |e: usize, lines: &mut Vec<[usize; 2]>| -> usize {
        *line_of_edge.entry(e).or_insert_with(|| {
            lines.push(model.edges()[e]);
            lines.len() - 1
        })
    };
    let mut surfaces = vec![];
    for s in 0..model.surfaces().len() {
        let surface = &model.surfaces()[s];
        let mut loops = vec![];
        for ring in &surface.boundaries {
            loops.push(
                ring.iter()
                    .map(|u| {
                        let l = use_edge(u.edge, &mut lines) as i64 + 1;
                        if u.reversed { -l } else { l }
                    })
                    .collect(),
            );
        }
        let embedded_lines = surface.embedded_edges.iter().map(|&e| use_edge(e, &mut lines)).collect();
        surfaces.push(SurfaceInput { loops, embedded_lines });
    }
    // Bar pieces between consecutive anchors.
    struct Piece {
        line: usize,
        axis: usize,
        t: [f64; 2],
        stiffness: u32,
    }
    let mut pieces: Vec<Piece> = vec![];
    for (i, axis) in state.axes.iter().enumerate() {
        for pair in axis.anchors.windows(2) {
            let (a, b) = (pair[0].vertex, pair[1].vertex);
            if a == b {
                continue;
            }
            let line = match model.edge_between(a, b) {
                Some(e) => use_edge(e, &mut lines),
                None => {
                    lines.push([a, b]);
                    lines.len() - 1
                }
            };
            let mid = (pair[0].t + pair[1].t) / 2.;
            let stiffness = axis
                .spans
                .iter()
                .find(|s| s.start_t.min(s.end_t) - 1e-9 <= mid && mid <= s.start_t.max(s.end_t) + 1e-9)
                .or_else(|| axis.spans.first())
                .map_or(0, |s| s.stiffness);
            pieces.push(Piece { line, axis: i, t: [pair[0].t, pair[1].t], stiffness });
        }
    }
    // Bar ends inside a surface (no edge leads to them): embedded points.
    let surface_vertices: Vec<BTreeSet<usize>> = (0..model.surfaces().len())
        .map(|s| model.surface_edges(s).flat_map(|e| model.edges()[e]).collect())
        .collect();
    let mut embedded_points = vec![];
    for c in &state.contacts {
        if let Contact::Point { surface, vertex, .. } = c {
            if *surface < surface_vertices.len() && !surface_vertices[*surface].contains(vertex) {
                embedded_points.push((*vertex, *surface));
            }
        }
    }
    // Only vertices that a line or embedded point uses.
    let mut used: BTreeSet<usize> = lines.iter().flatten().copied().collect();
    used.extend(embedded_points.iter().map(|p| p.0));
    // A vertex on the inside of a line (a bar or a wall meeting the middle of
    // an edge) splits that line, so that every neighbour shares the nodes.
    let tolerance = model.precision().max(1e-9) * 4.;
    let chains = split_lines(&lines, &used.iter().copied().collect::<Vec<_>>(), model.vertices(), tolerance);
    // Lines of the chains; a surface loop or a bar piece uses all segments of its line.
    let mut split_lines_list: Vec<[usize; 2]> = vec![];
    let mut segment_start: Vec<usize> = vec![];
    for chain in &chains {
        segment_start.push(split_lines_list.len());
        for w in chain.windows(2) {
            split_lines_list.push([w[0], w[1]]);
        }
    }
    let segments_of = |line: usize| -> std::ops::Range<usize> {
        let start = segment_start[line];
        start..start + chains[line].len() - 1
    };
    for s in &mut surfaces {
        s.loops = s
            .loops
            .iter()
            .map(|l| {
                l.iter()
                    .flat_map(|&c| {
                        let line = (c.unsigned_abs() as usize) - 1;
                        let mut ids: Vec<i64> = segments_of(line).map(|k| k as i64 + 1).collect();
                        if c < 0 {
                            ids.reverse();
                            ids.iter_mut().for_each(|x| *x = -*x);
                        }
                        ids
                    })
                    .collect()
            })
            .collect();
        s.embedded_lines = s.embedded_lines.iter().flat_map(|&l| segments_of(l)).collect();
    }
    let compact: BTreeMap<usize, usize> = used.iter().enumerate().map(|(i, &v)| (v, i)).collect();
    let input = Input {
        vertices: used.iter().map(|&v| model.vertices()[v]).collect(),
        lines: split_lines_list.iter().map(|l| [compact[&l[0]], compact[&l[1]]]).collect(),
        surfaces,
        embedded_points: embedded_points.iter().map(|&(v, s)| (compact[&v], s)).collect(),
        size,
        quads,
    };
    let out = gmsh.mesh(&input)?;
    let mut mesh = Mesh {
        nodes: out.nodes.clone(),
        vertex_nodes: vec![None; model.vertices().len()],
        ..Default::default()
    };
    for (&v, &c) in &compact {
        mesh.vertex_nodes[v] = Some(out.vertex_nodes[c]);
    }
    for (s, elements) in out.surface_elements.iter().enumerate() {
        for nodes in elements {
            mesh.shells.push(Shell { nodes: nodes.clone(), surface: s, stiffness: state.stiffness[s] });
        }
    }
    for piece in &pieces {
        let range = segments_of(piece.line);
        let chain = &chains[piece.line];
        let [ta, tb] = piece.t;
        // The piece follows its line when its first anchor is the line's first vertex.
        let first = state.axes[piece.axis].anchors.iter().find(|x| (x.t - ta).abs() < 1e-12).map(|x| x.vertex);
        let forward = first == Some(chain[0]) || first.is_none();
        let mut parts: Vec<([usize; 2], f64)> = vec![];
        for k in range {
            for seg in &out.line_elements[k] {
                let (p, q) = (out.nodes[seg[0]], out.nodes[seg[1]]);
                let length = ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt();
                parts.push((*seg, length));
            }
        }
        let total: f64 = parts.iter().map(|p| p.1).sum::<f64>().max(1e-12);
        let mut done = 0.;
        for (nodes, length) in parts {
            let (f0, f1) = (done / total, (done + length) / total);
            done += length;
            let t = if forward { [ta + (tb - ta) * f0, ta + (tb - ta) * f1] } else { [tb + (ta - tb) * f0, tb + (ta - tb) * f1] };
            mesh.bars.push(BarPiece { nodes, axis: piece.axis, stiffness: piece.stiffness, t });
        }
    }
    // Gmsh leaves nodes that belong to no element (copies of curve nodes);
    // keep the referenced ones.
    let mut keep = vec![false; mesh.nodes.len()];
    for n in mesh.shells.iter().flat_map(|s| s.nodes.iter()).chain(mesh.bars.iter().flat_map(|b| b.nodes.iter())) {
        keep[*n] = true;
    }
    let mut renumber = vec![usize::MAX; mesh.nodes.len()];
    let mut nodes = vec![];
    for (i, p) in mesh.nodes.iter().enumerate() {
        if keep[i] {
            renumber[i] = nodes.len();
            nodes.push(*p);
        }
    }
    mesh.nodes = nodes;
    for s in &mut mesh.shells {
        s.nodes.iter_mut().for_each(|n| *n = renumber[*n]);
    }
    for b in &mut mesh.bars {
        b.nodes.iter_mut().for_each(|n| *n = renumber[*n]);
    }
    for v in mesh.vertex_nodes.iter_mut() {
        *v = v.map(|n| renumber[n]).filter(|&n| n != usize::MAX);
    }
    Ok(mesh)
}

/// Every line as the chain of vertices from its first to its last, with the
/// vertices of `vertices_in_use` that lie inside it (within `tolerance`).
fn split_lines(lines: &[[usize; 2]], in_use: &[usize], vertices: &[[f64; 3]], tolerance: f64) -> Vec<Vec<usize>> {
    let cell = (tolerance * 64.).max(1e-3);
    let key = |p: [f64; 3]| [0, 1, 2].map(|k| (p[k] / cell).floor() as i64);
    let mut grid: std::collections::HashMap<[i64; 3], Vec<usize>> = Default::default();
    for &v in in_use {
        grid.entry(key(vertices[v])).or_default().push(v);
    }
    lines
        .iter()
        .map(|&[a, b]| {
            let (p, q) = (glam::DVec3::from_array(vertices[a]), glam::DVec3::from_array(vertices[b]));
            let d = q - p;
            let length = d.length();
            let mut inside: Vec<(f64, usize)> = vec![];
            if length > 2. * tolerance {
                let (lo, hi) = (key((p.min(q) - glam::DVec3::splat(tolerance)).to_array()), key((p.max(q) + glam::DVec3::splat(tolerance)).to_array()));
                // Cells of the bounding box; a long diagonal line visits few of them in practice.
                let cells = (hi[0] - lo[0] + 1) as i128 * (hi[1] - lo[1] + 1) as i128 * (hi[2] - lo[2] + 1) as i128;
                let candidates: Vec<usize> = if cells <= 4096 {
                    let mut c = vec![];
                    for x in lo[0]..=hi[0] {
                        for y in lo[1]..=hi[1] {
                            for z in lo[2]..=hi[2] {
                                c.extend(grid.get(&[x, y, z]).into_iter().flatten().copied());
                            }
                        }
                    }
                    c
                } else {
                    in_use.to_vec()
                };
                let u = d / length;
                for v in candidates {
                    if v == a || v == b {
                        continue;
                    }
                    let x = glam::DVec3::from_array(vertices[v]);
                    let t = (x - p).dot(u);
                    if t > tolerance && t < length - tolerance && (x - p - u * t).length() <= tolerance {
                        inside.push((t, v));
                    }
                }
            }
            inside.sort_by(|x, y| x.0.total_cmp(&y.0));
            inside.dedup_by_key(|x| x.1);
            let mut chain = vec![a];
            chain.extend(inside.into_iter().map(|x| x.1));
            chain.push(b);
            chain
        })
        .collect()
}
