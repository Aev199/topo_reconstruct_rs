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
