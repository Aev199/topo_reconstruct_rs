//! Reconstruct an auditable real FE fragment selected from a full preview report.
use clap::Parser;
use geo::{LineString, Polygon, TriangulateEarcut};
use glam::DVec3;
use serde_json::{json, Value};
use std::{collections::BTreeSet, fs::File};
use topo_reconstruct_rs::{
    input::MeshData,
    parsers::LiraParser,
    reconstruction::{assembly, frame, mesh, planes, recognize},
};

#[derive(Parser)]
struct Args {
    input: String,
    #[arg(long)]
    report: String,
    /// Surface indices in the supplied full-model topology preview.
    #[arg(long, required = true)]
    surface: Vec<usize>,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    let full = LiraParser::parse(&args.input)?;
    let previous: Value = serde_json::from_reader(File::open(&args.report)?)?;
    let ids = previous["frame"]["node_ids"]
        .as_array()
        .ok_or("missing reference ids")?;
    let points = previous["frame"]["reference_points"]
        .as_array()
        .ok_or("missing reference points")?;
    if ids.len() != points.len() {
        return Err("reference lengths differ".into());
    }
    for (id, p) in ids.iter().zip(points) {
        let id = id.as_u64().ok_or("invalid source id")? as u32;
        let xyz: [f64; 3] = serde_json::from_value(p.clone())?;
        if full
            .nodes
            .get(&id)
            .is_none_or(|v| v.distance(DVec3::from_array(xyz)) > 1e-10)
        {
            return Err("preview does not match source coordinates".into());
        }
    }
    let surfaces = previous["topology"]["preview"]["surfaces"]
        .as_array()
        .ok_or("missing surfaces")?;
    let mut selected = BTreeSet::new();
    for &s in &args.surface {
        let surface = surfaces.get(s).ok_or("surface index out of range")?;
        for e in surface["source_elements"]
            .as_array()
            .ok_or("missing surface provenance")?
        {
            selected.insert(e.as_u64().ok_or("invalid element id")? as u32);
        }
    }
    if selected.is_empty()
        || selected
            .iter()
            .any(|id| !full.elements.iter().any(|e| e.id == *id && e.is_shell()))
    {
        return Err("invalid shell selection".into());
    }
    let shell_nodes: BTreeSet<_> = full
        .elements
        .iter()
        .filter(|e| selected.contains(&e.id))
        .flat_map(|e| e.nodes.iter().copied())
        .collect();
    let axis_policy = recognize::Policy {
        angle: 0.02,
        line_tolerance: 0.01,
        numerical_precision: 1e-8,
    };
    let all_axes = recognize::recognize(&full, &axis_policy)?;
    // Include whole recognized axes, never cut bars at the fragment box.
    for axis in &all_axes.axes {
        if axis.anchors.iter().any(|a| shell_nodes.contains(&a.node)) {
            selected.extend(axis.spans.iter().map(|s| s.element));
        }
    }
    // Invalid touching source bars must also remain visible to recognition.
    for e in &full.elements {
        if e.is_bar() && e.nodes.iter().any(|n| shell_nodes.contains(n)) {
            selected.insert(e.id);
        }
    }
    let elements: Vec<_> = full
        .elements
        .iter()
        .filter(|e| selected.contains(&e.id))
        .cloned()
        .collect();
    let node_ids: BTreeSet<_> = elements
        .iter()
        .flat_map(|e| e.nodes.iter().copied())
        .collect();
    let source = MeshData {
        nodes: node_ids.iter().map(|id| (*id, full.nodes[id])).collect(),
        elements,
    };
    let cut:Vec<_>=full.elements.iter().filter(|e|!selected.contains(&e.id) && e.nodes.iter().any(|n|node_ids.contains(n))).map(|e|json!({"element":e.id,"type":e.elem_type,"shared_nodes":e.nodes.iter().filter(|n|node_ids.contains(n)).collect::<Vec<_>>()})).collect();
    let global_rejections: Vec<_> = previous["topology"]["axis_assembly"]["issues"]
        .as_array()
        .ok_or("missing global axis issues")?
        .iter()
        .filter(|i| {
            i["source_elements"].as_array().is_some_and(|es| {
                es.iter()
                    .any(|e| selected.contains(&(e.as_u64().unwrap_or(0) as u32)))
            })
        })
        .cloned()
        .collect();
    let axes = recognize::recognize(&source, &axis_policy)?;
    let planes = planes::recognize(
        &source,
        &planes::Policy {
            angle: 0.02,
            distance: 0.01,
            precision: 1e-8,
        },
    )?;
    let frame = frame::solve(
        &source,
        &axes,
        &planes,
        &frame::Policy {
            up: [0., 0., 1.],
            angle: 0.02,
            maximum_movement: 0.15,
            relative_movement: 0.05,
            minimum_length: 0.03,
            residual_tolerance: 1e-7,
            iterations: 1000,
        },
    )?;
    let topology = assembly::assemble(
        &source,
        &frame,
        &assembly::Policy {
            closure_tolerance: 0.001,
            junction_movement_limit: 0.05,
            precision: 1e-7,
            minimum_edge: 0.001,
        },
    )?;
    let (mesh, mesh_error) = match mesh::build(
        &topology,
        &mesh::Policy {
            boundary_spacing: 0.5,
            maximum_area: 0.5,
            minimum_angle_degrees: 20.,
            maximum_added_vertices_per_surface: 10000,
        },
    ) {
        Ok(m) => (Some(m), None),
        Err(e) => (None, Some(e)),
    };
    // Display triangulation only: never presented as a new FE mesh.
    let mut display = Vec::new();
    for (i, s) in topology.preview.surfaces().iter().enumerate() {
        let ring = |r: &Vec<[f64; 2]>| {
            LineString::from(
                r.iter()
                    .chain(r.first())
                    .map(|p| (p[0], p[1]))
                    .collect::<Vec<_>>(),
            )
        };
        let polygon = Polygon::new(
            ring(&s.contours[0]),
            s.contours.iter().skip(1).map(ring).collect(),
        );
        let plane = &topology.preview.planes()[s.plane];
        for t in polygon.earcut_triangles() {
            display.push(json!({"surface":i,"vertices":[plane.lift([t.v1().x,t.v1().y]),plane.lift([t.v2().x,t.v2().y]),plane.lift([t.v3().x,t.v3().y])]}));
        }
    }
    let raw_nodes: Vec<_> = node_ids
        .iter()
        .map(|n| json!({"id":n,"point":source.nodes[n].to_array()}))
        .collect();
    let raw_elements: Vec<_> = source
        .elements
        .iter()
        .map(|e| json!({"id":e.id,"type":e.elem_type,"stiffness":e.stiff_id,"nodes":e.nodes}))
        .collect();
    let output = json!({"selection":{"surface_indices":args.surface,"element_ids":selected,"cut_connections":cut,"original_global_axis_rejections":global_rejections,"local_result_is_not_full_model_acceptance":true},"source":{"nodes":raw_nodes,"elements":raw_elements},"frame":frame,"topology":topology,"mesh":mesh,"mesh_error":mesh_error,"display_triangles":display,"export_ready":false});
    serde_json::to_writer_pretty(std::io::stdout().lock(), &output)?;
    Ok(())
}
