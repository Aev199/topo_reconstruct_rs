//! Apply editor edits to a reconstructed model and print how the audit
//! changes (development: checks that the editor's operations close the
//! defects left by the automatic reconstruction).
//!
//! Usage: edit_probe MODEL.txt [--frame-cache PATH] [--write REPORT.json] [--plaxis EXCHANGE.json] [EDIT_JSON ...]
//! Each edit is the JSON of `editor::Edit` (`{"op": "merge_vertices", ...}`);
//! `{"op": "show", "surfaces": [..], "bars": [..], "vertices": [..]}`
//! prints the named objects instead.
use serde_json::Value;
use std::collections::BTreeMap;
use std::path::PathBuf;
use topo_reconstruct_rs::audit::Class;
use topo_reconstruct_rs::editor::{Edit, Session};
use topo_reconstruct_rs::pipeline::{self, Options, Profile};

fn counts(session: &Session) -> BTreeMap<(String, String), usize> {
    let mut out = BTreeMap::new();
    for f in &session.audit().findings {
        if f.class != Class::Review {
            *out.entry((format!("{:?}", f.class), f.kind.clone()))
                .or_default() += 1;
        }
    }
    out
}

fn show(session: &Session, what: &Value) {
    let state = session.state();
    let model = &state.model;
    let list = |key: &str| -> Vec<usize> {
        what.get(key)
            .and_then(Value::as_array)
            .map(|a| {
                a.iter()
                    .filter_map(Value::as_u64)
                    .map(|x| x as usize)
                    .collect()
            })
            .unwrap_or_default()
    };
    for s in list("surfaces") {
        let surface = &model.surfaces()[s];
        println!(
            "surface {s}: stiffness {} plane {} normal {:?} embedded {:?}",
            state.stiffness[s],
            surface.plane,
            model.planes()[surface.plane].normal(),
            surface
                .embedded_edges
                .iter()
                .map(|&e| (e, model.edges()[e]))
                .collect::<Vec<_>>()
        );
        for r in &surface.boundaries {
            let vs: Vec<usize> = r
                .iter()
                .map(|u| model.edges()[u.edge][usize::from(u.reversed)])
                .collect();
            println!(
                "  ring {:?}",
                vs.iter()
                    .zip(r)
                    .map(|(&v, u)| (v, u.edge, model.vertices()[v]))
                    .collect::<Vec<_>>()
            );
        }
    }
    for e in list("edges") {
        let [a, b] = model.edges()[e];
        println!(
            "edge {e}: {a} {:?} - {b} {:?}",
            model.vertices()[a],
            model.vertices()[b]
        );
    }
    for b in list("bars") {
        let a = &state.axes[b];
        println!(
            "bar {b}: ends {:?} {:?} anchors {:?}",
            a.endpoints,
            a.endpoints.map(|v| model.vertices()[v]),
            a.anchors
                .iter()
                .map(|x| (x.vertex, x.t, x.source_node))
                .collect::<Vec<_>>()
        );
    }
    for v in list("vertices") {
        let users: Vec<usize> = (0..model.surfaces().len())
            .filter(|&s| {
                model
                    .surface_edges(s)
                    .any(|e| model.edges()[e].contains(&v))
            })
            .collect();
        let bars: Vec<usize> = (0..state.axes.len())
            .filter(|&b| state.axes[b].anchors.iter().any(|a| a.vertex == v))
            .collect();
        println!(
            "vertex {v}: {:?} surfaces {users:?} bars {bars:?}",
            model.vertices()[v]
        );
    }
}

fn main() -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let mut args = std::env::args().skip(1);
    let input = PathBuf::from(
        args.next()
            .ok_or("usage: edit_probe MODEL.txt [--frame-cache PATH] EDIT...")?,
    );
    let mut cache = None;
    let mut write = None;
    let mut plaxis = None;
    let mut edits = vec![];
    while let Some(a) = args.next() {
        if a == "--plaxis" {
            plaxis = args.next().map(PathBuf::from);
        } else if a == "--write" {
            write = args.next().map(PathBuf::from);
        } else if a == "--frame-cache" {
            cache = args.next().map(PathBuf::from);
        } else {
            edits.push(a);
        }
    }
    let profile = Profile::plaxis();
    let options = Options {
        mesh: false,
        frame_cache: cache,
    };
    let output = pipeline::run(&input, &profile, &options, &mut |_| {})?;
    let mut session = Session::new(&output.topology, profile.audit_options());
    println!("before: {:?}", counts(&session));
    let mut last = String::new();
    for text in edits {
        // `"$v"` names the vertex the previous edit generated.
        let text = text.replace("\"$v\"", &last);
        let value: Value = serde_json::from_str(&text)?;
        if value["op"] == "edge" {
            let (a, b) = (
                value["a"].as_u64().unwrap() as usize,
                value["b"].as_u64().unwrap() as usize,
            );
            let model = &session.state().model;
            let found: Vec<usize> = (0..model.edges().len())
                .filter(|&e| {
                    let [x, y] = model.edges()[e];
                    (x, y) == (a, b) || (x, y) == (b, a)
                })
                .collect();
            println!("edge {a}-{b}: {found:?}");
            continue;
        }
        if value["op"] == "show" {
            show(&session, &value);
            continue;
        }
        let edit: Edit = serde_json::from_value(value)?;
        let before = counts(&session);
        match session.apply(edit, "probe") {
            Ok(_) => {
                let after = counts(&session);
                let keys: std::collections::BTreeSet<_> =
                    before.keys().chain(after.keys()).collect();
                let diff: Vec<String> = keys
                    .into_iter()
                    .filter(|k| before.get(*k) != after.get(*k))
                    .map(|k| {
                        format!(
                            "{}/{}: {} -> {}",
                            k.0,
                            k.1,
                            before.get(k).unwrap_or(&0),
                            after.get(k).unwrap_or(&0)
                        )
                    })
                    .collect();
                let outcome = &session.journal().last().unwrap().outcome;
                if let Some(v) = outcome
                    .rsplit("vertex ")
                    .next()
                    .filter(|_| outcome.contains("vertex "))
                {
                    last = v.trim().to_string();
                }
                println!("{text}: {outcome} | {}", diff.join(", "));
            }
            Err(e) => println!("{text}: REFUSED {e}"),
        }
    }
    if let Some(path) = plaxis {
        let materials = topo_reconstruct_rs::parsers::lira::LiraParser::parse_materials(&input)?;
        let started = std::time::Instant::now();
        let exchange = topo_reconstruct_rs::plaxis::exchange(
            session.state(),
            &materials,
            topo_reconstruct_rs::plaxis::TONNE_TO_KN,
            profile.edge_collapse,
            &input.display().to_string(),
        );
        println!(
            "plaxis: {} plates, {} polygons, {} cut, {} beams, missing {:?}, {:?}",
            exchange.plates.len(),
            exchange
                .plates
                .iter()
                .map(|p| p.polygons.len())
                .sum::<usize>(),
            exchange.cut_surfaces,
            exchange.beams.len(),
            exchange.missing_materials,
            started.elapsed()
        );
        std::fs::write(&path, serde_json::to_vec(&exchange)?)?;
        // Contours of the holed surfaces, for checks of the cutting.
        let model = &session.state().model;
        let holed: std::collections::BTreeMap<usize, _> = (0..model.surfaces().len())
            .filter(|&s| model.surfaces()[s].contours.len() > 1)
            .map(|s| (s, model.surfaces()[s].contours.clone()))
            .collect();
        std::fs::write(
            path.with_extension("contours.json"),
            serde_json::to_vec(&holed)?,
        )?;
    }
    if let Some(path) = write {
        let mut value = serde_json::to_value(&output)?;
        session.write_into(&mut value)?;
        std::fs::write(path, serde_json::to_vec(&value)?)?;
    }
    Ok(())
}
