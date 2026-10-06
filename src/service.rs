//! Command layer of the desktop editor. The Tauri application and the
//! development HTTP bridge (`examples/editor_server.rs`) both forward their
//! calls to [`Service::dispatch`]: one command name and JSON arguments in,
//! JSON out. All state lives here, so the window holds only a view.
use crate::audit;
use crate::editor::{content_hash, Edit, Project, Session};
use crate::pipeline::{self, Options, Profile};
use glam::DVec2;
use serde::Serialize;
use serde_json::{json, Value};
use std::path::{Path, PathBuf};

pub struct Service {
    /// Where solved frames are cached between runs.
    cache_dir: PathBuf,
    profile: Profile,
    input: Option<PathBuf>,
    /// Content hash of the input as it was read for the open session.
    input_hash: Option<String>,
    output: Option<pipeline::Output>,
    session: Option<Session>,
    /// The journal as last saved in (or opened from) a project: edits
    /// differing from it are unsaved.
    saved_journal: Vec<crate::editor::Entry>,
    /// Materials of the opened input bytes (the snapshot the geometry is).
    materials: Option<hashbrown::HashMap<u32, crate::parsers::lira::Material>>,
}

/// Display data of the current geometry.
#[derive(Debug, Serialize)]
pub struct Scene {
    pub vertices: Vec<[f64; 3]>,
    pub surfaces: Vec<SceneSurface>,
    /// Used model edges: (edge id, vertex, vertex, embedded).
    pub edges: Vec<(usize, usize, usize, bool)>,
    /// Bar pieces: (axis, vertex, vertex).
    pub bars: Vec<(usize, usize, usize)>,
    pub bounds: [[f64; 3]; 2],
}

#[derive(Debug, Serialize)]
pub struct SceneSurface {
    pub stiffness: u32,
    pub patch: usize,
    pub normal: [f64; 3],
    pub area: f64,
    pub source_elements: usize,
    /// Display triangles over model vertex indices (flat).
    pub triangles: Vec<usize>,
}

impl Service {
    pub fn new(cache_dir: PathBuf) -> Self {
        Service {
            cache_dir,
            profile: Profile::plaxis(),
            input: None,
            input_hash: None,
            output: None,
            session: None,
            saved_journal: vec![],
            materials: None,
        }
    }

    fn session(&self) -> Result<&Session, String> {
        self.session
            .as_ref()
            .ok_or_else(|| "no model is open".to_string())
    }

    fn session_mut(&mut self) -> Result<&mut Session, String> {
        self.session
            .as_mut()
            .ok_or_else(|| "no model is open".to_string())
    }

    /// Run one command. `progress` receives pipeline stage names.
    pub fn dispatch(
        &mut self,
        command: &str,
        args: Value,
        progress: &mut dyn FnMut(&str),
    ) -> Result<Value, String> {
        let text = |key: &str| -> Result<String, String> {
            args.get(key)
                .and_then(Value::as_str)
                .map(str::to_string)
                .ok_or_else(|| format!("missing argument `{key}`"))
        };
        let to_value = |v: &dyn erased::Ser| v.value();
        match command {
            "profile" => Ok(serde_json::to_value(&self.profile).map_err(|e| e.to_string())?),
            "open_model" => {
                if let Some(p) = args.get("profile") {
                    self.profile = serde_json::from_value(p.clone()).map_err(|e| e.to_string())?;
                }
                self.open_model(Path::new(&text("path")?), progress)?;
                self.summary()
            }
            "summary" => self.summary(),
            "scene" => to_value(&self.scene()?),
            "audit" => to_value(self.session()?.audit()),
            "apply" => {
                let edit: Edit = serde_json::from_value(
                    args.get("edit").cloned().ok_or("missing argument `edit`")?,
                )
                .map_err(|e| format!("invalid edit: {e}"))?;
                let note = args.get("note").and_then(Value::as_str).unwrap_or("");
                self.session_mut()?.apply(edit, note)?;
                self.summary()
            }
            "undo" => {
                self.session_mut()?.undo()?;
                self.summary()
            }
            "redo" => {
                self.session_mut()?.redo()?;
                self.summary()
            }
            "journal" => to_value(&self.session()?.journal().to_vec()),
            "save_project" => {
                let changed = self.save_project(Path::new(&text("path")?))?;
                Ok(json!({"saved": text("path")?, "input_changed_on_disk": changed}))
            }
            "open_project" => self.open_project(Path::new(&text("path")?), progress),
            "export_plaxis" => {
                let factor = args
                    .get("force_factor")
                    .and_then(Value::as_f64)
                    .unwrap_or(crate::plaxis::TONNE_TO_KN);
                let stiffness = match args.get("stiffness").and_then(Value::as_str) {
                    Some("nominal") => crate::plaxis::StiffnessMode::Nominal,
                    _ => crate::plaxis::StiffnessMode::Effective,
                };
                let with_loads = args.get("loads").and_then(Value::as_bool).unwrap_or(true);
                let combination = combination_from(args.get("combination"))?;
                // The cases to read in the mode "by case" (`include_cases`).
                let cases = args.get("include_cases").and_then(Value::as_array).map(|a| {
                    a.iter().filter_map(Value::as_u64).map(|c| c as u32).collect()
                });
                self.export_plaxis(Path::new(&text("path")?), factor, stiffness, with_loads, combination, cases, CutOptions::from(args.get("cut")))
            }
            "load_cases" => self.load_cases(),
            "floors" => self.floors(),
            "export_midas" => self.export_midas(&args),
            "run_plaxis" => self.run_plaxis(&args),
            "export_report" => {
                self.export_report(Path::new(&text("path")?))?;
                Ok(json!({"saved": text("path")?}))
            }
            _ => Err(format!("unknown command `{command}`")),
        }
    }

    fn open_model(&mut self, path: &Path, progress: &mut dyn FnMut(&str)) -> Result<(), String> {
        let bytes = std::fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
        std::fs::create_dir_all(&self.cache_dir).map_err(|e| e.to_string())?;
        let cache = self
            .cache_dir
            .join(format!("{}.frame.json", content_hash(&bytes)));
        let options = Options {
            mesh: false,
            frame_cache: Some(cache),
        };
        // A library panic must not take the window down with it.
        let profile = self.profile.clone();
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            pipeline::run(path, &profile, &options, progress)
        }))
        .map_err(|_| "the reconstruction failed (internal error)".to_string())?
        .map_err(|e| e.to_string())?;
        // The session must be the reconstruction of exactly these bytes.
        let after = std::fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
        if content_hash(&after) != content_hash(&bytes) {
            return Err(format!(
                "{} changed while it was being read; open it again",
                path.display()
            ));
        }
        self.session = Some(Session::new(&result.topology, profile.audit_options()));
        self.output = Some(result);
        self.input = Some(path.to_path_buf());
        self.input_hash = Some(content_hash(&bytes));
        self.materials = Some(crate::parsers::lira::LiraParser::materials_from(&bytes));
        self.saved_journal.clear();
        Ok(())
    }

    fn summary(&self) -> Result<Value, String> {
        let session = self.session()?;
        let audit = session.audit();
        let state = session.state();
        Ok(json!({
            "input": self.input.as_ref().map(|p| p.display().to_string()),
            "surfaces": state.model.surfaces().len(),
            "bars": state.axes.len(),
            "edits": session.journal().len(),
            // Edits not in the last saved project (opening another model
            // would lose them).
            "dirty": session.journal() != self.saved_journal.as_slice(),
            "can_redo": session.can_redo(),
            "audit": {
                "passed": audit.passed,
                "plaxis_passed": audit.plaxis_passed,
                "counts": audit.counts,
                "failures": audit.findings.iter().filter(|f| f.class == audit::Class::Failure).count(),
                "plaxis": audit.findings.iter().filter(|f| f.class == audit::Class::Plaxis).count(),
                "review": audit.findings.iter().filter(|f| f.class == audit::Class::Review).count(),
            },
            "last": session.journal().last(),
            "cut": state.cut.as_ref().map(|c| json!({"z": c.z, "top": c.top})),
        }))
    }

    /// Geometry for the 3D view.
    pub fn scene(&self) -> Result<Scene, String> {
        let state = self.session()?.state();
        let model = &state.model;
        let vertices = model.vertices().to_vec();
        let mut lo = [f64::MAX; 3];
        let mut hi = [f64::MIN; 3];
        let mut edges = std::collections::BTreeMap::<usize, bool>::new();
        let mut surfaces = vec![];
        for (i, s) in model.surfaces().iter().enumerate() {
            for u in s.boundaries.iter().flatten() {
                edges.insert(u.edge, false);
            }
            for &e in &s.embedded_edges {
                edges.entry(e).or_insert(true);
            }
            let rings: Vec<Vec<usize>> = s
                .boundaries
                .iter()
                .map(|ring| {
                    ring.iter()
                        .map(|u| {
                            let [a, b] = model.edges()[u.edge];
                            if u.reversed {
                                b
                            } else {
                                a
                            }
                        })
                        .collect()
                })
                .collect();
            for &v in rings.iter().flatten() {
                for k in 0..3 {
                    lo[k] = lo[k].min(vertices[v][k]);
                    hi[k] = hi[k].max(vertices[v][k]);
                }
            }
            let plane = &model.planes()[s.plane];
            let uv: Vec<Vec<DVec2>> = rings
                .iter()
                .map(|r| {
                    r.iter()
                        .map(|&v| DVec2::from_array(plane.project(vertices[v])))
                        .collect()
                })
                .collect();
            let area = uv
                .iter()
                .enumerate()
                .map(|(k, r)| {
                    let a = (0..r.len())
                        .map(|j| r[j].perp_dot(r[(j + 1) % r.len()]))
                        .sum::<f64>()
                        .abs()
                        / 2.;
                    if k == 0 {
                        a
                    } else {
                        -a
                    }
                })
                .sum();
            surfaces.push(SceneSurface {
                stiffness: state.stiffness.get(i).copied().unwrap_or(0),
                patch: state.patches.get(i).copied().unwrap_or(0),
                normal: plane.normal(),
                area,
                source_elements: s.source_elements.len(),
                triangles: triangulate(&rings, &uv),
            });
        }
        let bars = state
            .axes
            .iter()
            .enumerate()
            .flat_map(|(i, axis)| {
                let mut nodes: Vec<(f64, usize)> =
                    axis.anchors.iter().map(|a| (a.t, a.vertex)).collect();
                nodes.push((0., axis.endpoints[0]));
                nodes.push((1., axis.endpoints[1]));
                nodes.sort_by(|x, y| x.0.total_cmp(&y.0));
                nodes.dedup_by_key(|n| n.1);
                nodes
                    .windows(2)
                    .map(move |w| (i, w[0].1, w[1].1))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        for &(_, a, b) in &bars {
            for v in [a, b] {
                for k in 0..3 {
                    lo[k] = lo[k].min(vertices[v][k]);
                    hi[k] = hi[k].max(vertices[v][k]);
                }
            }
        }
        Ok(Scene {
            edges: edges
                .into_iter()
                .map(|(e, embedded)| {
                    let [a, b] = model.edges()[e];
                    (e, a, b, embedded)
                })
                .collect(),
            vertices,
            surfaces,
            bars,
            // An empty model (every surface deleted) gets a finite unit box.
            bounds: if lo[0] <= hi[0] {
                [lo, hi]
            } else {
                [[-1.; 3], [1.; 3]]
            },
        })
    }

    /// Save the project; the input hash is that of the file the session
    /// was reconstructed from. Returns whether the file on disk has changed
    /// since.
    fn save_project(&mut self, path: &Path) -> Result<bool, String> {
        let input = self.input.as_ref().ok_or("no model is open")?;
        let opened = self.input_hash.clone().ok_or("no model is open")?;
        let changed = std::fs::read(input)
            .map(|bytes| content_hash(&bytes) != opened)
            .unwrap_or(true);
        let project = Project {
            format: 1,
            input: input.display().to_string(),
            input_hash: opened,
            profile: self.profile.clone(),
            journal: self.session()?.journal().to_vec(),
        };
        let text = serde_json::to_string_pretty(&project).map_err(|e| e.to_string())?;
        std::fs::write(path, text).map_err(|e| e.to_string())?;
        self.saved_journal = project.journal;
        Ok(changed)
    }

    fn open_project(
        &mut self,
        path: &Path,
        progress: &mut dyn FnMut(&str),
    ) -> Result<Value, String> {
        let text = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
        let project: Project =
            serde_json::from_str(&text).map_err(|e| format!("invalid project: {e}"))?;
        // The model path is kept as chosen; a relative one is next to the
        // project file.
        let mut input = PathBuf::from(&project.input);
        if input.is_relative() {
            if let Some(dir) = path.parent() {
                input = dir.join(input);
            }
        }
        let bytes = std::fs::read(&input).map_err(|e| format!("{}: {e}", input.display()))?;
        let changed = content_hash(&bytes) != project.input_hash;
        self.profile = project.profile.clone();
        self.open_model(&input, progress)?;
        let (applied, stopped) = self.session_mut()?.replay(&project.journal);
        // A fully replayed project is saved as is; a partly replayed one
        // differs from the file (saving it would drop the rest).
        self.saved_journal = if stopped.is_none() {
            self.session()?.journal().to_vec()
        } else {
            project.journal.clone()
        };
        let mut summary = self.summary()?;
        summary["replayed"] = json!(applied);
        summary["replay_stopped"] = json!(stopped);
        summary["input_changed"] = json!(changed);
        Ok(summary)
    }

    /// Write the PLAXIS exchange file and, beside it, the loader script.
    fn export_plaxis(
        &self,
        path: &Path,
        force_factor: f64,
        stiffness: crate::plaxis::StiffnessMode,
        with_loads: bool,
        combination: Option<crate::loads::Combination>,
        cases: Option<std::collections::BTreeSet<u32>>,
        cut: CutOptions,
    ) -> Result<Value, String> {
        let input = self.input.as_ref().ok_or("no model is open")?;
        let session = self.session()?;
        // Materials of the bytes the geometry was reconstructed from, never
        // of a file changed on disk since.
        let materials = self.materials.as_ref().ok_or("no model is open")?;
        // The storeys that were cut off: a cap slab with their stiffness.
        let cap = if cut.cap { crate::storeys::with_cap(session.state(), materials, cut.factor) } else { None };
        let state = cap.as_ref().map_or(session.state(), |c| &c.0);
        let exchange_materials = cap.as_ref().map_or(materials, |c| &c.1);
        let exchange = crate::plaxis::exchange(
            state,
            exchange_materials,
            crate::plaxis::Settings {
                force_factor,
                min_edge: self.profile.edge_collapse,
                stiffness,
            },
            &input.display().to_string(),
        );
        let mut exchange = exchange;
        // Bars of S1..S6 sections and profiles of the block 13 (the exchange reads S0 bars only).
        if let Ok(bytes) = self.input_bytes() {
            let profiles = crate::parsers::lira::LiraParser::profiles_from(&bytes);
            crate::plaxis::add_section_beams(&mut exchange, exchange_materials, &profiles, stiffness, force_factor);
        }
        if with_loads {
            self.add_loads(&mut exchange, state, force_factor, combination, cases, cut.loads)?;
        }
        let file = std::fs::File::create(path).map_err(|e| e.to_string())?;
        serde_json::to_writer(std::io::BufWriter::new(file), &exchange)
            .map_err(|e| e.to_string())?;
        let script = path.with_file_name("plaxis_export.py");
        std::fs::write(&script, crate::plaxis::LOADER).map_err(|e| e.to_string())?;
        Ok(json!({
            "saved": path.display().to_string(),
            "script": script.display().to_string(),
            "plates": exchange.plates.len(),
            "polygons": exchange.plates.iter().map(|p| p.polygons.len()).sum::<usize>(),
            "cut_surfaces": exchange.cut_surfaces,
            "triangulated_surfaces": exchange.triangulated_surfaces,
            "beams": exchange.beams.len(),
            "plate_materials": exchange.plate_materials.len(),
            "beam_materials": exchange.beam_materials.len(),
            "missing_materials": exchange.missing_materials,
            "audit_passed": session.audit().passed,
            "warnings": exchange.warnings,
            "loads": exchange.loads.len(),
            "load_cases": exchange.load_report.as_ref().map(|r| r.cases.len()).unwrap_or(0),
            "load_report": exchange.load_report,
            "load_problems": exchange.load_report.as_ref().map(|r| r.problems(crate::loads::FORCE_TOLERANCE, crate::loads::MOMENT_TOLERANCE)).unwrap_or_default(),
            "cap": cap.as_ref().map(|c| &c.2),
            "material_notes": exchange.plate_materials.iter().map(|m| (m.name.clone(), m.notes.clone()))
                .chain(exchange.beam_materials.iter().map(|m| (m.name.clone(), m.notes.clone())))
                .filter(|(_, n)| !n.is_empty())
                .collect::<Vec<_>>(),
        }))
    }

    /// The loads of the opened input, mapped onto the current geometry. The
    /// file must still be the one the geometry was reconstructed from.
    fn add_loads(
        &self,
        exchange: &mut crate::plaxis::Exchange,
        state: &crate::reconstruction::assembly::edit::State,
        force_factor: f64,
        combination: Option<crate::loads::Combination>,
        cases: Option<std::collections::BTreeSet<u32>>,
        cut_loads: bool,
    ) -> Result<(), String> {
        let bytes = self.input_bytes()?;
        let mesh = crate::parsers::lira::LiraParser::mesh_from(&bytes).map_err(|e| e.to_string())?;
        let set = crate::parsers::loads::parse(&bytes);
        let output = self.output.as_ref().ok_or("no model is open")?;
        let (loads, report) = crate::loads::transfer(
            state,
            &output.topology.vertex_source_nodes,
            &mesh,
            &set,
            crate::loads::Settings {
                force_factor,
                snap: self.profile.edge_collapse,
                max_groups: 40,
                combination: combination.clone(),
                cases,
                materials: self.materials.clone().map(std::sync::Arc::new),
                cut_loads,
            },
        );
        exchange.load_cases = if combination.is_some() {
            vec![(crate::loads::COMBINATION, "Сочетание".to_string())]
        } else {
            set.cases.clone()
        };
        // PLAXIS deletes point loads that lie on no plate and no beam: move them onto the structure.
        let polygons: Vec<Vec<[f64; 3]>> = exchange.plates.iter().flat_map(|p| p.polygons.iter().cloned()).collect();
        let segments: Vec<([f64; 3], [f64; 3])> = exchange.beams.iter().map(|b| (b.start, b.end)).collect();
        let mut loads = loads;
        let (moved, farthest) = crate::loads::attach_points(&mut loads, &polygons, &segments);
        let mut report = report;
        if moved > 0 {
            report.approximated.insert(format!("точечные нагрузки перенесены на ближайшую плиту или балку (до {:.0} мм, момент сохранён парой)", farthest * 1000.), moved);
        }
        exchange.loads = loads;
        exchange.load_report = Some(report);
        Ok(())
    }

    /// Mesh the geometry with Gmsh, carry the loads of the chosen cases onto
    /// the mesh and write a MIDAS Civil `.mxt` file.
    fn export_midas(&self, args: &Value) -> Result<Value, String> {
        let path = args.get("path").and_then(Value::as_str).ok_or("missing argument `path`")?;
        let size = args.get("size").and_then(Value::as_f64).filter(|s| *s > 0.).unwrap_or(self.profile.element_size);
        let quads = args.get("quads").and_then(Value::as_bool).unwrap_or(false);
        let factor = args.get("force_factor").and_then(Value::as_f64).unwrap_or(crate::plaxis::TONNE_TO_KN);
        let gmsh = crate::gmsh::Gmsh::load().map_err(|e| format!("gmsh_missing: {e}"))?;
        let cut = CutOptions::from(args.get("cut"));
        let materials = self.materials.as_ref().ok_or("no model is open")?;
        let cap = if cut.cap { crate::storeys::with_cap(self.session()?.state(), materials, cut.factor) } else { None };
        let state = cap.as_ref().map_or(self.session()?.state(), |c| &c.0);
        let exchange_materials = cap.as_ref().map_or(materials, |c| &c.1);
        let cut_loads = cut.loads;
        let mesh = crate::meshing::mesh_state(&gmsh, state, size, quads).map_err(|e| format!("meshing_failed: {e}"))?;
        // Cases: the chosen ones, else all but the self-weight, stages and dynamics.
        let bytes = self.input_bytes()?;
        let source = crate::parsers::lira::LiraParser::mesh_from(&bytes).map_err(|e| e.to_string())?;
        let set = crate::parsers::loads::parse(&bytes);
        let selected: std::collections::BTreeSet<u32> = match args.get("include_cases").and_then(Value::as_array) {
            Some(a) => a.iter().filter_map(Value::as_u64).map(|c| c as u32).collect(),
            None => set
                .cases
                .iter()
                .filter(|(_, n)| !crate::loads::is_self_weight(n) && !crate::loads::is_stage(n) && !crate::loads::is_dynamic(n))
                .map(|(c, _)| *c)
                .chain(state.cut.iter().map(|_| crate::loads::CUT_WEIGHT_CASE))
                .collect(),
        };
        let output = self.output.as_ref().ok_or("no model is open")?;
        let (loads, load_report) = crate::loads::transfer(
            state,
            &output.topology.vertex_source_nodes,
            &source,
            &set,
            crate::loads::Settings {
                force_factor: factor,
                snap: self.profile.edge_collapse,
                max_groups: 40,
                combination: None,
                cases: Some(selected.clone()),
                materials: self.materials.clone().map(std::sync::Arc::new),
                cut_loads,
            },
        );
        let on_mesh = crate::mesh_loads::transfer(&mesh, state, &loads, 0.02);
        // Materials, sections and thicknesses by the rules of the converter (or the stiffness LIRA analysed with).
        let section_mode = match args.get("midas_stiffness").and_then(Value::as_str) {
            Some("lira") => crate::midas_stiffness::Mode::Lira,
            _ => crate::midas_stiffness::Mode::Converter,
        };
        let density_multiplier = args.get("density_multiplier").and_then(Value::as_f64).filter(|m| m.is_finite() && *m >= 0.).unwrap_or(1.);
        let profiles = crate::parsers::lira::LiraParser::profiles_from(&bytes);
        let plate_types: std::collections::BTreeSet<u32> = mesh.shells.iter().map(|s| s.stiffness).collect();
        let bar_types: std::collections::BTreeSet<u32> = mesh.bars.iter().map(|b| b.stiffness).collect();
        let stiffness = crate::midas_stiffness::build(
            exchange_materials,
            &profiles,
            &plate_types,
            &bar_types,
            crate::midas_stiffness::Options { mode: section_mode, density_multiplier },
        );
        let mut cases: Vec<(u32, String)> = set.cases.iter().filter(|(c, _)| selected.contains(c)).cloned().collect();
        if state.cut.is_some() && cut_loads && selected.contains(&crate::loads::CUT_WEIGHT_CASE) {
            cases.push((crate::loads::CUT_WEIGHT_CASE, "Вес отброшенных этажей".into()));
        }
        let (text, report) = crate::midas::write_mxt(&mesh, &on_mesh, &stiffness, &cases);
        std::fs::write(path, text).map_err(|e| format!("{path}: {e}"))?;
        // Force and moment of every case at each stage: source, geometry, mesh.
        let about = on_mesh.resultants_about(&mesh, glam::DVec3::from_array(load_report.origin));
        let mut problems = load_report.problems(crate::loads::FORCE_TOLERANCE, crate::loads::MOMENT_TOLERANCE);
        let per_case: Vec<Value> = load_report
            .cases
            .iter()
            .map(|c| {
                let (f, m, scale) = about.get(&c.case).copied().unwrap_or_default();
                if let Some(what) = crate::loads::compare_resultants(
                    (c.exported, c.exported_moment),
                    (f.to_array(), m.to_array()),
                    c.moment_scale.max(scale),
                    crate::loads::FORCE_TOLERANCE,
                    crate::loads::MOMENT_TOLERANCE,
                ) {
                    problems.push(format!("загружение {} «{}» на сетке: {what}", c.case, c.name));
                }
                json!({"case": c.case, "name": c.name, "source": c.source, "geometry": c.exported, "mesh": f.to_array(),
                       "source_moment": c.source_moment, "geometry_moment": c.exported_moment, "mesh_moment": m.to_array(),
                       "moment_scale": c.moment_scale.max(scale), "lost": on_mesh.lost.get(&c.case)})
            })
            .collect();
        Ok(json!({
            "saved": path,
            "gmsh": gmsh.path().display().to_string(),
            "report": report,
            "triangles": mesh.shells.iter().filter(|s| s.nodes.len() == 3).count(),
            "quads": mesh.shells.iter().filter(|s| s.nodes.len() == 4).count(),
            "cases": per_case,
            "load_problems": problems,
            "cap": cap.as_ref().map(|c| &c.2),
            "skipped": load_report.skipped,
            "warnings": ["LIRA rotation angles of bar sections are not read: the program default local axes (beta 0)"],
            "stiffness_notes": stiffness.notes,
            "missing_materials": stiffness.missing,
        }))
    }

    /// The source file's bytes, when it is still the file the geometry was
    /// reconstructed from.
    fn input_bytes(&self) -> Result<Vec<u8>, String> {
        let input = self.input.as_ref().ok_or("no model is open")?;
        let bytes = std::fs::read(input).map_err(|e| format!("{}: {e}", input.display()))?;
        if Some(content_hash(&bytes)) != self.input_hash {
            return Err(format!(
                "input_changed: {} changed after it was opened; open it again to export loads",
                input.display()
            ));
        }
        Ok(bytes)
    }

    /// The load cases of the opened input for the combination dialog.
    fn load_cases(&self) -> Result<Value, String> {
        let set = crate::parsers::loads::parse(&self.input_bytes()?);
        let mut rows: std::collections::BTreeMap<u32, usize> = Default::default();
        for row in &set.rows {
            *rows.entry(row.case).or_default() += 1;
        }
        let names: std::collections::BTreeMap<u32, &String> = set.cases.iter().map(|(n, s)| (*n, s)).collect();
        let cases: Vec<Value> = rows
            .iter()
            .map(|(&case, &count)| {
                let name = names.get(&case).map(|s| s.to_string()).unwrap_or_default();
                json!({"case": case, "name": name, "rows": count, "self_weight": crate::loads::is_self_weight(&name), "dynamic": crate::loads::is_dynamic(&name), "stage": crate::loads::is_stage(&name)})
            })
            .collect();
        let mut cases = cases;
        if self.session()?.state().cut.is_some() {
            cases.push(json!({"case": crate::loads::CUT_WEIGHT_CASE, "name": "Вес отброшенных этажей", "rows": 0,
                              "self_weight": false, "dynamic": false, "stage": false}));
        }
        Ok(json!({"cases": cases}))
    }

    /// The floors of the geometry and the cut made on it, for the dialog.
    fn floors(&self) -> Result<Value, String> {
        let state = self.session()?.state();
        Ok(json!({
            "floors": crate::reconstruction::assembly::cutoff::floors(state),
            "cut": state.cut,
        }))
    }

    /// Run the loader with a Python that has plxscripting (the PLAXIS
    /// distribution): it builds the model in the open PLAXIS Input.
    fn run_plaxis(&self, args: &Value) -> Result<Value, String> {
        let get = |key: &str| args.get(key).and_then(Value::as_str);
        let exchange = PathBuf::from(get("path").ok_or("missing argument `path`")?);
        let script = exchange.with_file_name("plaxis_export.py");
        if !script.exists() {
            std::fs::write(&script, crate::plaxis::LOADER).map_err(|e| e.to_string())?;
        }
        let python = get("python")
            .filter(|p| !p.is_empty())
            .map(PathBuf::from)
            .or_else(crate::plaxis::find_python)
            .ok_or("no Python with plxscripting found: give its path")?;
        let mut command = std::process::Command::new(&python);
        command
            .arg(&script)
            .arg(&exchange)
            .args(["--host", get("host").unwrap_or("localhost")])
            .args([
                "--port",
                &args
                    .get("port")
                    .and_then(Value::as_u64)
                    .unwrap_or(10000)
                    .to_string(),
            ])
            .args(["--password", get("password").unwrap_or("")]);
        #[cfg(windows)]
        {
            // No console window from the desktop application.
            use std::os::windows::process::CommandExt;
            command.creation_flags(0x0800_0000);
        }
        if args.get("new").and_then(Value::as_bool).unwrap_or(false) {
            command.arg("--new");
        }
        if args.get("phases").and_then(Value::as_bool) == Some(false) {
            command.arg("--no-phases");
        }
        if args
            .get("shift_to_origin")
            .and_then(Value::as_bool)
            .unwrap_or(false)
        {
            command.arg("--shift-to-origin");
        }
        let output = command
            .output()
            .map_err(|e| format!("{}: {e}", python.display()))?;
        let stdout = String::from_utf8_lossy(&output.stdout).to_string();
        let stderr = String::from_utf8_lossy(&output.stderr).to_string();
        if !output.status.success() {
            let tail: String = stderr
                .lines()
                .rev()
                .take(8)
                .collect::<Vec<_>>()
                .into_iter()
                .rev()
                .collect::<Vec<_>>()
                .join("\n");
            return Err(format!("plaxis_loader_failed: {tail}"));
        }
        // The report is the JSON object printed last.
        let report = stdout
            .rfind("\n{")
            .map(|i| &stdout[i + 1..])
            .or_else(|| stdout.starts_with('{').then_some(stdout.as_str()))
            .and_then(|t| serde_json::from_str::<Value>(t).ok())
            .unwrap_or(json!({"output": stdout}));
        Ok(json!({"python": python.display().to_string(), "report": report}))
    }

    fn export_report(&self, path: &Path) -> Result<(), String> {
        let output = self.output.as_ref().ok_or("no model is open")?;
        let mut value = serde_json::to_value(output).map_err(|e| e.to_string())?;
        self.session()?
            .write_into(&mut value)
            .map_err(|e| e.to_string())?;
        let file = std::fs::File::create(path).map_err(|e| e.to_string())?;
        serde_json::to_writer(std::io::BufWriter::new(file), &value).map_err(|e| e.to_string())
    }
}

/// Serialization behind a trait object, for the dispatch table.
mod erased {
    pub trait Ser {
        fn value(&self) -> Result<serde_json::Value, String>;
    }
    impl<T: serde::Serialize> Ser for T {
        fn value(&self) -> Result<serde_json::Value, String> {
            serde_json::to_value(self).map_err(|e| e.to_string())
        }
    }
}

/// Display triangulation of a polygon with holes over its ring vertices
/// (constrained Delaunay, triangles inside the material). Falls back to a
/// fan of the exterior if the constraints cannot be inserted.
fn triangulate(rings: &[Vec<usize>], uv: &[Vec<DVec2>]) -> Vec<usize> {
    use spade::{ConstrainedDelaunayTriangulation, Point2, Triangulation};
    let inside = |p: DVec2| {
        let mut odd = false;
        for r in uv {
            for i in 0..r.len() {
                let (a, b) = (r[i], r[(i + 1) % r.len()]);
                if (a.y > p.y) != (b.y > p.y) && p.x < a.x + (p.y - a.y) / (b.y - a.y) * (b.x - a.x)
                {
                    odd = !odd;
                }
            }
        }
        odd
    };
    let fan = || -> Vec<usize> {
        let r = &rings[0];
        (1..r.len().saturating_sub(1))
            .flat_map(|i| [r[0], r[i], r[i + 1]])
            .collect()
    };
    let mut cdt = ConstrainedDelaunayTriangulation::<Point2<f64>>::new();
    let mut handles = vec![];
    let mut owner = std::collections::HashMap::new();
    for (r, ring) in uv.iter().enumerate() {
        let mut hs = vec![];
        for (k, p) in ring.iter().enumerate() {
            let Ok(h) = cdt.insert(Point2::new(p.x, p.y)) else {
                return fan();
            };
            if owner
                .insert(h, rings[r][k])
                .is_some_and(|v| v != rings[r][k])
            {
                return fan();
            }
            hs.push(h);
        }
        handles.push(hs);
    }
    for hs in &handles {
        for i in 0..hs.len() {
            let (a, b) = (hs[i], hs[(i + 1) % hs.len()]);
            if a == b || !cdt.can_add_constraint(a, b) {
                return fan();
            }
            cdt.add_constraint(a, b);
        }
    }
    let mut out = vec![];
    for face in cdt.inner_faces() {
        let v = face.vertices();
        let ps = v.map(|h| h.position());
        let c = DVec2::new(
            (ps[0].x + ps[1].x + ps[2].x) / 3.,
            (ps[0].y + ps[1].y + ps[2].y) / 3.,
        );
        if inside(c) {
            for h in v {
                out.push(owner[&h.fix()]);
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn square_with_a_hole_is_triangulated_inside_only() {
        let rings = vec![vec![0, 1, 2, 3], vec![4, 5, 6, 7]];
        let uv = vec![
            vec![
                DVec2::new(0., 0.),
                DVec2::new(4., 0.),
                DVec2::new(4., 4.),
                DVec2::new(0., 4.),
            ],
            vec![
                DVec2::new(1., 1.),
                DVec2::new(1., 2.),
                DVec2::new(2., 2.),
                DVec2::new(2., 1.),
            ],
        ];
        let t = triangulate(&rings, &uv);
        assert_eq!(t.len() % 3, 0);
        let point = |v: usize| if v < 4 { uv[0][v] } else { uv[1][v - 4] };
        let area: f64 = t
            .chunks(3)
            .map(|c| ((point(c[1]) - point(c[0])).perp_dot(point(c[2]) - point(c[0]))).abs() / 2.)
            .sum();
        assert!((area - 15.).abs() < 1e-9, "{area}");
    }

    #[test]
    fn commands_need_an_open_model() {
        let mut s = Service::new(std::env::temp_dir());
        assert!(s.dispatch("scene", json!({}), &mut |_| {}).is_err());
        assert!(s.dispatch("nonsense", json!({}), &mut |_| {}).is_err());
        assert_eq!(
            s.dispatch("profile", json!({}), &mut |_| {}).unwrap()["element_size"],
            json!(0.5)
        );
    }
}

/// `{"cases": [{"case": 1, "factor": 1.35}, ...], "simplify": true,
/// "center_tolerance": 0.15, "min_fraction": 0.3}` as a load combination.
fn combination_from(value: Option<&Value>) -> Result<Option<crate::loads::Combination>, String> {
    let Some(value) = value.filter(|v| !v.is_null()) else {
        return Ok(None);
    };
    let mut factors = std::collections::BTreeMap::new();
    for entry in value.get("cases").and_then(Value::as_array).ok_or("combination: missing `cases`")? {
        let case = entry.get("case").and_then(Value::as_u64).ok_or("combination: bad case number")?;
        let factor = entry.get("factor").and_then(Value::as_f64).ok_or("combination: bad factor")?;
        if !factor.is_finite() {
            return Err("combination: factor is not a number".into());
        }
        factors.insert(case as u32, factor);
    }
    let defaults = crate::loads::Simplify::default();
    let simplify = value.get("simplify").and_then(Value::as_bool).unwrap_or(true).then(|| crate::loads::Simplify {
        center_tolerance: value.get("center_tolerance").and_then(Value::as_f64).unwrap_or(defaults.center_tolerance),
        min_fraction: value.get("min_fraction").and_then(Value::as_f64).unwrap_or(defaults.min_fraction),
        max_points: value.get("max_points").and_then(Value::as_u64).map_or(defaults.max_points, |n| (n as usize).max(1)),
    });
    Ok(Some(crate::loads::Combination { factors, simplify }))
}

/// What the export does with the storeys that were cut off:
/// `{"loads": true, "cap": true, "factor": 1.0}`.
#[derive(Debug, Clone, Copy)]
pub struct CutOptions {
    /// Their loads and weight go to the supports at the level.
    pub loads: bool,
    /// A cap slab with their stiffness.
    pub cap: bool,
    /// Factor of the equivalent stiffness.
    pub factor: f64,
}

impl CutOptions {
    fn from(value: Option<&Value>) -> CutOptions {
        let get = |k: &str| value.and_then(|v| v.get(k));
        CutOptions {
            loads: get("loads").and_then(Value::as_bool).unwrap_or(true),
            cap: get("cap").and_then(Value::as_bool).unwrap_or(true),
            factor: get("factor").and_then(Value::as_f64).filter(|f| *f >= 0.).unwrap_or(1.),
        }
    }
}
