//! Editing session over a reconstructed geometry: every edit is a
//! transactional model operation (`assembly::edit`), logged in a journal
//! with its provenance and re-checked at once by the same audit. Undo
//! replays the journal from the reconstruction; a project stores only the
//! input, the profile and the journal, and every replayed edit is checked
//! against the geometry it was made on.
use crate::audit::{self, Class};
use crate::reconstruction::assembly::edit::State;
use crate::reconstruction::assembly::Report;
use serde::{Deserialize, Serialize};

/// One manual edit. Indices refer to the geometry at the time of the edit.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(tag = "op", rename_all = "snake_case")]
pub enum Edit {
    MoveVertex {
        vertex: usize,
        to: [f64; 3],
    },
    MergeVertices {
        drop: usize,
        keep: usize,
    },
    DeleteSurface {
        surface: usize,
    },
    JoinSurfaces {
        keep: usize,
        other: usize,
    },
    SplitEdge {
        edge: usize,
        at: [f64; 3],
    },
    /// The gap is a defect: close it, moving the vertex at most `tolerance`.
    CloseGap {
        vertex: usize,
        surface: usize,
        tolerance: f64,
    },
    /// The gap is a joint of the structure: keep it open.
    MarkJoint {
        vertex: usize,
        surface: usize,
    },
    DeleteBar {
        bar: usize,
    },
    /// Two crossing bars share a node (within `tolerance` of the crossing).
    ConnectBars {
        a: usize,
        b: usize,
        tolerance: f64,
    },
    /// A bar shares a node with every surface it passes through.
    ConnectBarToSurfaces {
        bar: usize,
    },
    /// Two surfaces get a shared edge along their junction.
    ConnectSurfaces {
        a: usize,
        b: usize,
    },
}

/// What an edit was made on, checked again when the journal is replayed
/// (the indices of a newer reconstruction may name other objects).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default)]
pub struct Check {
    /// Positions of the vertices the edit names.
    pub vertices: Vec<[f64; 3]>,
    /// Source elements of the surfaces it names.
    pub surfaces: Vec<Vec<u32>>,
    /// End positions of the edge it names.
    pub edge: Vec<[f64; 3]>,
    /// End positions of the bars it names.
    #[serde(default)]
    pub bars: Vec<[f64; 3]>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct Entry {
    pub edit: Edit,
    /// Why (the user's note, or the audit finding it answers).
    pub note: String,
    pub check: Check,
    /// What the operation reported.
    pub outcome: String,
}

pub struct Session {
    base: State,
    state: State,
    journal: Vec<Entry>,
    redo: Vec<Entry>,
    options: audit::Options,
    audit: audit::Report,
}

/// Position tolerance of a replay check (model units).
const CHECK_TOLERANCE: f64 = 1e-6;

impl Session {
    pub fn new(report: &Report, options: audit::Options) -> Self {
        let state = State::from_report(report);
        let mut session = Session {
            base: state.clone(),
            state,
            journal: vec![],
            redo: vec![],
            options,
            audit: audit::Report::default(),
        };
        session.reaudit();
        session
    }

    pub fn state(&self) -> &State {
        &self.state
    }
    pub fn journal(&self) -> &[Entry] {
        &self.journal
    }
    pub fn can_redo(&self) -> bool {
        !self.redo.is_empty()
    }
    /// The audit of the current geometry. Gaps marked as joints are review
    /// items (`accepted_joint`), not PLAXIS items.
    pub fn audit(&self) -> &audit::Report {
        &self.audit
    }

    fn reaudit(&mut self) {
        let mut report = audit::run(
            &self.state.model,
            &self.state.axes,
            &self.state.contacts,
            &self.options,
        );
        for f in &mut report.findings {
            if f.kind == "gap" {
                if let (Some(v), Some(&s)) = (f.vertex, f.surfaces.last()) {
                    if self.state.joints.contains(&(v, s)) {
                        f.kind = "accepted_joint".into();
                        f.class = Class::Review;
                    }
                }
            }
        }
        report.counts.clear();
        for f in &report.findings {
            *report.counts.entry(f.kind.clone()).or_default() += 1;
        }
        report.plaxis_passed =
            report.passed && !report.findings.iter().any(|f| f.class == Class::Plaxis);
        self.audit = report;
    }

    fn check_of(state: &State, edit: &Edit) -> Check {
        let vertex = |v: usize| state.model.vertices().get(v).copied();
        let surface = |s: usize| {
            state
                .model
                .surfaces()
                .get(s)
                .map(|x| x.source_elements.clone())
        };
        let bar = |b: usize| -> Vec<[f64; 3]> {
            state
                .axes
                .get(b)
                .map(|a| a.endpoints.iter().filter_map(|&v| vertex(v)).collect())
                .unwrap_or_default()
        };
        let mut check = Check::default();
        match *edit {
            Edit::DeleteBar { bar: b } | Edit::ConnectBarToSurfaces { bar: b } => {
                check.bars = bar(b)
            }
            Edit::ConnectBars { a, b, .. } => {
                check.bars = bar(a);
                check.bars.extend(bar(b));
            }
            Edit::MoveVertex { vertex: v, .. } => check.vertices.extend(vertex(v)),
            Edit::MergeVertices { drop, keep } => check
                .vertices
                .extend([vertex(drop), vertex(keep)].into_iter().flatten()),
            Edit::DeleteSurface { surface: s } => check.surfaces.extend(surface(s)),
            Edit::ConnectSurfaces { a: keep, b: other } | Edit::JoinSurfaces { keep, other } => {
                check
                    .surfaces
                    .extend([surface(keep), surface(other)].into_iter().flatten())
            }
            Edit::SplitEdge { edge, .. } => {
                if let Some(e) = state.model.edges().get(edge) {
                    check.edge = e.iter().filter_map(|&v| vertex(v)).collect();
                }
            }
            Edit::CloseGap {
                vertex: v,
                surface: s,
                ..
            }
            | Edit::MarkJoint {
                vertex: v,
                surface: s,
            } => {
                check.vertices.extend(vertex(v));
                check.surfaces.extend(surface(s));
            }
        }
        check
    }

    fn matches(a: &Check, b: &Check) -> bool {
        let close = |x: &[[f64; 3]], y: &[[f64; 3]]| {
            x.len() == y.len()
                && x.iter().zip(y).all(|(p, q)| {
                    (0..3).map(|k| (p[k] - q[k]).powi(2)).sum::<f64>().sqrt() <= CHECK_TOLERANCE
                })
        };
        close(&a.vertices, &b.vertices)
            && close(&a.edge, &b.edge)
            && close(&a.bars, &b.bars)
            && a.surfaces == b.surfaces
    }

    fn execute(state: &mut State, edit: &Edit) -> Result<String, String> {
        match *edit {
            Edit::MoveVertex { vertex, to } => state.move_vertex(vertex, to),
            Edit::MergeVertices { drop, keep } => state.merge_vertices(drop, keep),
            Edit::DeleteSurface { surface } => state.delete_surface(surface),
            Edit::JoinSurfaces { keep, other } => state.join_surfaces(keep, other),
            Edit::SplitEdge { edge, at } => state.split_edge(edge, at),
            Edit::CloseGap {
                vertex,
                surface,
                tolerance,
            } => state.close_gap(vertex, surface, tolerance),
            Edit::MarkJoint { vertex, surface } => state.mark_joint(vertex, surface),
            Edit::DeleteBar { bar } => state.delete_bar(bar),
            Edit::ConnectBars { a, b, tolerance } => state.connect_bars(a, b, tolerance),
            Edit::ConnectBarToSurfaces { bar } => state.connect_bar_to_surfaces(bar),
            Edit::ConnectSurfaces { a, b } => state.connect_surfaces(a, b),
        }
    }

    /// Apply an edit. On error nothing changes.
    pub fn apply(&mut self, edit: Edit, note: &str) -> Result<&audit::Report, String> {
        let check = Self::check_of(&self.state, &edit);
        let mut trial = self.state.clone();
        let outcome = Self::execute(&mut trial, &edit)?;
        self.state = trial;
        self.journal.push(Entry {
            edit,
            note: note.into(),
            check,
            outcome,
        });
        self.redo.clear();
        self.reaudit();
        Ok(&self.audit)
    }

    /// Undo the last edit (the journal without it is replayed).
    pub fn undo(&mut self) -> Result<(), String> {
        let last = self.journal.pop().ok_or("nothing to undo")?;
        let mut state = self.base.clone();
        for entry in &self.journal {
            Self::execute(&mut state, &entry.edit).map_err(|e| format!("replay failed: {e}"))?;
        }
        self.state = state;
        self.redo.push(last);
        self.reaudit();
        Ok(())
    }

    /// Apply the last undone edit again.
    pub fn redo(&mut self) -> Result<(), String> {
        let entry = self.redo.pop().ok_or("nothing to redo")?;
        let mut trial = self.state.clone();
        let outcome = Self::execute(&mut trial, &entry.edit)?;
        self.state = trial;
        self.journal.push(Entry { outcome, ..entry });
        self.reaudit();
        Ok(())
    }

    /// Replay a saved journal on a fresh reconstruction. Stops at the first
    /// edit whose objects are not where they were (a different
    /// reconstruction) or which fails; returns how many were applied and
    /// why it stopped.
    pub fn replay(&mut self, journal: &[Entry]) -> (usize, Option<String>) {
        for (k, entry) in journal.iter().enumerate() {
            let check = Self::check_of(&self.state, &entry.edit);
            if !Self::matches(&check, &entry.check) {
                return (
                    k,
                    Some(format!(
                        "edit {k}: the geometry differs from the one it was made on"
                    )),
                );
            }
            if let Err(e) = self.apply(entry.edit.clone(), &entry.note) {
                return (k, Some(format!("edit {k}: {e}")));
            }
        }
        (journal.len(), None)
    }

    /// The reconstruction report with the edited geometry: the model, bars,
    /// contacts and per-surface data replaced, removed surfaces and the
    /// journal recorded, the audit refreshed; a stale trial mesh dropped.
    pub fn write_into(&self, output: &mut serde_json::Value) -> Result<(), serde_json::Error> {
        let topology = &mut output["topology"];
        topology["preview"] = serde_json::to_value(&self.state.model)?;
        topology["axis_assembly"]["axes"] = serde_json::to_value(&self.state.axes)?;
        topology["axis_assembly"]["contacts"] = serde_json::to_value(&self.state.contacts)?;
        topology["surface_stiffness"] = serde_json::to_value(&self.state.stiffness)?;
        topology["surface_source_patches"] = serde_json::to_value(&self.state.patches)?;
        // What the edits changed, for the auditors: vertices moved away
        // from the reconstruction (user decisions, logged in `edits`),
        // vertices created by splits, surfaces removed with provenance.
        let base = self.base.model.vertices();
        let current = self.state.model.vertices();
        let moved: Vec<usize> = (0..base.len().min(current.len()))
            .filter(|&v| base[v] != current[v])
            .collect();
        // Bars whose nodes the edits changed (by source axis): their nodes
        // and spans no longer follow the source bar one to one.
        let nodes = |a: &crate::reconstruction::assembly::bars::Axis| {
            (
                a.endpoints,
                a.anchors
                    .iter()
                    .map(|n| (n.vertex, n.t.to_bits()))
                    .collect::<Vec<_>>(),
            )
        };
        let base_axes: std::collections::BTreeMap<usize, _> = self
            .base
            .axes
            .iter()
            .map(|a| (a.source_axis, nodes(a)))
            .collect();
        let edited_bars: Vec<usize> = self
            .state
            .axes
            .iter()
            .filter(|a| base_axes.get(&a.source_axis) != Some(&nodes(a)))
            .map(|a| a.source_axis)
            .collect();
        topology["user_edits"] = serde_json::json!({
            "edited_bars": edited_bars,
            "moved_vertices": moved,
            "generated_vertices": (base.len()..current.len()).collect::<Vec<_>>(),
            "removed_surfaces": self.state.removed,
            "removed_bars": self.state.removed_bars,
            "accepted_joints": self.state.joints,
        });
        output["edits"] = serde_json::to_value(&self.journal)?;
        output["audit"] = serde_json::to_value(&self.audit)?;
        if let Some(map) = output.as_object_mut() {
            map.remove("mesh");
            map.remove("mesh_error");
        }
        Ok(())
    }
}

/// A saved project: what to reconstruct and the edits made on it.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct Project {
    pub format: u32,
    /// Input model path as chosen by the user.
    pub input: String,
    /// Hash of the input file at the time of the edits.
    pub input_hash: String,
    pub profile: crate::pipeline::Profile,
    pub journal: Vec<Entry>,
}

/// Content hash of an input file (FNV-1a 64, stable across platforms).
pub fn content_hash(bytes: &[u8]) -> String {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in bytes {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    format!("{h:016x}")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reconstruction::assembly::edit::State;

    fn session(state: State) -> Session {
        let mut s = Session {
            base: state.clone(),
            state,
            journal: vec![],
            redo: vec![],
            options: audit::Options::default(),
            audit: audit::Report::default(),
        };
        s.reaudit();
        s
    }

    fn slabs(gap: f64) -> State {
        use crate::reconstruction::{Model, PlaneFrame};
        let mut m = Model::new(1e-6, 0.001).unwrap();
        let plane = m.add_plane(PlaneFrame::new([0.; 3], [0., 0., 1.]).unwrap());
        let mut ring = |pts: [[f64; 3]; 4]| -> Vec<usize> {
            pts.iter().map(|&p| m.add_vertex(p).unwrap()).collect()
        };
        let a = ring([[0., 0., 0.], [2., 0., 0.], [2., 2., 0.], [0., 2., 0.]]);
        let b = ring([
            [2. + gap, 0., 0.],
            [5., 0., 0.],
            [5., 2., 0.],
            [2. + gap, 2., 0.],
        ]);
        m.add_surface(plane, vec![a], vec![1]).unwrap();
        m.add_surface(plane, vec![b], vec![2]).unwrap();
        State {
            model: m,
            axes: vec![],
            contacts: vec![],
            stiffness: vec![1, 1],
            patches: vec![0, 1],
            removed: vec![],
            removed_bars: vec![],
            joints: Default::default(),
        }
    }

    #[test]
    fn edits_are_audited_logged_undone_and_replayed() {
        let mut s = session(slabs(0.02));
        // Four gap items: the corners of each slab at the other one.
        let gaps = |s: &Session| s.audit().counts.get("gap").copied().unwrap_or(0);
        assert_eq!(gaps(&s), 4, "{:?}", s.audit().counts);
        // Mark one as a joint: accepted, not a PLAXIS item.
        s.apply(
            Edit::MarkJoint {
                vertex: 4,
                surface: 0,
            },
            "expansion joint",
        )
        .unwrap();
        assert_eq!(gaps(&s), 3);
        assert_eq!(s.audit().counts.get("accepted_joint"), Some(&1));
        // Close another: the corner moves onto the left slab.
        s.apply(
            Edit::CloseGap {
                vertex: 7,
                surface: 0,
                tolerance: 0.05,
            },
            "not a joint",
        )
        .unwrap();
        let closed = gaps(&s);
        assert!(closed < 3, "{:?}", s.audit().counts);
        assert!(s.audit().passed, "{:?}", s.audit().findings);
        assert_eq!(s.journal().len(), 2);
        // A refused edit changes nothing and is not logged.
        assert!(s
            .apply(
                Edit::MoveVertex {
                    vertex: 0,
                    to: [0., 0., 1.]
                },
                ""
            )
            .is_err());
        assert_eq!(s.journal().len(), 2);
        // Undo, redo.
        s.undo().unwrap();
        assert_eq!(gaps(&s), 3);
        assert!(s.can_redo());
        s.redo().unwrap();
        assert_eq!(gaps(&s), closed);
        // Replay on a fresh session: same geometry; on another geometry the
        // replay stops at the first edit.
        let journal = s.journal().to_vec();
        let mut again = session(slabs(0.02));
        assert_eq!(again.replay(&journal), (2, None));
        assert_eq!(
            format!("{:?}", again.state().model),
            format!("{:?}", s.state().model)
        );
        let mut other = session(slabs(0.03));
        let (applied, reason) = other.replay(&journal);
        assert_eq!(applied, 0);
        assert!(reason.unwrap().contains("differs"));
    }

    #[test]
    fn project_round_trips() {
        let p = Project {
            format: 1,
            input: "model.txt".into(),
            input_hash: content_hash(b"abc"),
            profile: crate::pipeline::Profile::plaxis(),
            journal: vec![Entry {
                edit: Edit::DeleteSurface { surface: 3 },
                note: "n".into(),
                check: Check::default(),
                outcome: "deleted".into(),
            }],
        };
        let text = serde_json::to_string(&p).unwrap();
        assert_eq!(serde_json::from_str::<Project>(&text).unwrap(), p);
        assert_eq!(content_hash(b"abc"), "e71fa2190541574b");
    }
}
