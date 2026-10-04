//! The geotechnical reconstruction pipeline: LIRA text model -> recognized
//! axes and planes -> joint frame solve -> topological assembly with the
//! geotechnical rules -> optional trial mesh. One profile of tolerances,
//! tuned for PLAXIS, is the default; the CLI and the desktop application
//! both call [`run`].
use crate::input::{self, MeshData};
use crate::parsers::{LiraParser, Section};
use crate::reconstruction::{assembly, frame, graph, mesh, planes, recognize, reconcile};
use serde::Serialize;
use std::path::{Path, PathBuf};

/// Tolerances of the geotechnical model (model units, metres).
#[derive(Debug, Clone, Serialize, serde::Deserialize, PartialEq)]
pub struct Profile {
    /// Target element size of the downstream mesh (PLAXIS: 0.5 m); the
    /// trial mesh uses it as boundary spacing and its square as maximum
    /// triangle area.
    pub element_size: f64,
    /// Stacked-wall and wall-line alignment; the vertex closure movement
    /// limit is at least this plus the 1 mm closure tolerance.
    pub stack_offset: f64,
    /// Wall-end closure and redundant-vertex short-edge threshold.
    pub wall_end_snap: f64,
    /// Widest console trimmed beyond a junction line.
    pub console_width: f64,
    /// Widest crack of a converted mesh rebuilt inside one structure.
    pub crack_width: f64,
    /// Edges and bars shorter than this between needed corners collapse.
    pub edge_collapse: f64,
    /// Gaps narrower than this between structures close (h/10).
    pub gap_closure: f64,
    /// Close offsets across a plane (a wall top below a slab).
    pub close_offset_gaps: bool,
    /// Largest simplification tolerance (half the plate thickness within
    /// [gap closure, this]).
    pub simplification_cap: f64,
    /// Free openings narrower than this are filled.
    pub min_opening: f64,
    /// Openings longer than this are kept whatever their width.
    pub max_opening_length: f64,
    /// Base iteration budget of the frame solver.
    pub iterations: usize,
}

impl Profile {
    /// PLAXIS 3D: target elements of 0.5 m, features below h/10 closed.
    pub fn plaxis() -> Self {
        Profile {
            element_size: 0.5,
            stack_offset: 0.05,
            wall_end_snap: 0.05,
            console_width: 0.25,
            crack_width: 0.01,
            edge_collapse: 0.05,
            gap_closure: 0.05,
            close_offset_gaps: true,
            simplification_cap: 0.2,
            min_opening: 1.0,
            max_opening_length: 3.0,
            iterations: 1000,
        }
    }

    /// Every length finite and non-negative, the element size and the
    /// iteration budget positive.
    pub fn validate(&self) -> Result<(), String> {
        for (name, value) in [
            ("stack offset", self.stack_offset),
            ("wall end snap", self.wall_end_snap),
            ("console width", self.console_width),
            ("crack width", self.crack_width),
            ("edge collapse", self.edge_collapse),
            ("gap closure", self.gap_closure),
            ("simplification cap", self.simplification_cap),
            ("minimum opening", self.min_opening),
            ("maximum opening length", self.max_opening_length),
        ] {
            if !value.is_finite() || value < 0. {
                return Err(format!("{name} must be a finite non-negative length"));
            }
        }
        if !self.element_size.is_finite() || self.element_size <= 0. {
            return Err("element size must be a finite positive length".into());
        }
        if self.iterations == 0 {
            return Err("iterations must be positive".into());
        }
        Ok(())
    }
}

impl Default for Profile {
    fn default() -> Self {
        Self::plaxis()
    }
}

/// What to run besides the reconstruction itself.
#[derive(Debug, Clone, Default)]
pub struct Options {
    /// Build the trial mesh (validates meshability; slow on large models).
    pub mesh: bool,
    /// Reuse a solved frame stored for the same input and frame policy
    /// (development and the application's project cache).
    pub frame_cache: Option<PathBuf>,
}

/// Rigid links of the analysis model left out of the geometry.
#[derive(Debug, Serialize)]
pub struct RigidLinks {
    pub removed: usize,
    pub elements: Vec<u32>,
}

/// Every stage's report. Serialized, it is the reconstruction JSON read by
/// the auditors in `scripts/`.
#[derive(Debug, Serialize)]
pub struct Output {
    pub profile: Profile,
    pub constraint_graph: graph::Graph,
    pub frame: frame::Report,
    pub topology: assembly::Report,
    pub reconciliation: reconcile::Report,
    pub axis_recognition: recognize::Report,
    pub plane_recognition: planes::Report,
    pub relaxation: Relaxation,
    pub rigid_links: RigidLinks,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mesh: Option<mesh::Report>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mesh_error: Option<String>,
}

/// Frame relaxation step used (0: default rules, 1: panels merged at
/// over-constrained nodes).
#[derive(Debug, Serialize)]
pub struct Relaxation {
    pub frame: usize,
}

pub type Error = Box<dyn std::error::Error + Send + Sync>;

/// Run the pipeline on a LIRA text model. `stage` is called with the name
/// of every finished stage (progress, timing).
pub fn run(
    input: &Path,
    profile: &Profile,
    options: &Options,
    stage: &mut dyn FnMut(&str),
) -> Result<Output, Error> {
    profile.validate()?;
    let mut mesh_data = LiraParser::parse(input)?;
    // Ties between structures in the analysis model (before rigid links are
    // dropped): structures tied together are joined, not kept apart.
    let connections = input::node_links(&mesh_data, &LiraParser::parse_rigid_bodies(input)?);
    // Rigid links of the analysis model (fans spreading a column into a
    // slab, offsets) are not structures: they are left out.
    let rigid_links = input::rigid_links(&mesh_data, &LiraParser::parse_stiffness(input)?);
    {
        let removed: std::collections::HashSet<u32> = rigid_links.iter().copied().collect();
        mesh_data.elements.retain(|e| !removed.contains(&e.id));
    }
    let surface_thickness = LiraParser::parse_sections(input)?
        .into_iter()
        .filter_map(|(id, section)| match section {
            Section::Plate { thickness } => Some((id, thickness)),
            _ => None,
        })
        .collect();
    stage("parse");
    let axes = recognize::recognize(
        &mesh_data,
        &recognize::Policy {
            angle: 0.02,
            line_tolerance: 0.01,
            numerical_precision: 1e-8,
        },
    )?;
    stage("axis_recognition");
    let plane_report = planes::recognize(
        &mesh_data,
        &planes::Policy {
            angle: 0.02,
            distance: 0.01,
            precision: 1e-8,
        },
    )?;
    stage("plane_recognition");
    let (frame_report, relaxation) = solve_frame(
        input,
        &mesh_data,
        &axes,
        &plane_report,
        profile,
        options,
        &rigid_links,
    )?;
    stage("frame");
    let assembly_policy = assembly::Policy {
        closure_tolerance: 0.001,
        // A stacked wall moves by its offset when closing onto the lower axis.
        // A vertex may move by the largest rule tolerance plus the closure
        // tolerance (a node within it counts as on a plane): an aligned wall
        // is offset by at most the tolerance at its junction, and its slight
        // non-parallelism elsewhere is measured against the closure tolerance.
        junction_movement_limit: profile.stack_offset.max(0.05) + 0.001,
        precision: assembly::PRECISION,
        minimum_edge: 0.001,
    };
    let topology = assembly::assemble_geotechnical(
        &mesh_data,
        &frame_report,
        &assembly_policy,
        &assembly::FeaturePolicy {
            maximum_console_width: profile.console_width,
            maximum_stack_offset: profile.stack_offset,
            maximum_wall_end_snap: profile.wall_end_snap,
            maximum_crack_width: profile.crack_width,
            maximum_collapsed_edge: profile.edge_collapse,
            maximum_gap: profile.gap_closure,
            close_offset_gaps: profile.close_offset_gaps,
            surface_thickness,
            maximum_simplification: profile.simplification_cap,
            minimum_opening_width: profile.min_opening,
            maximum_opening_length: profile.max_opening_length,
            connections,
            ..Default::default()
        },
    )?;
    stage("assembly");
    let reconciliation = reconcile::solve(
        &mesh_data,
        &frame_report,
        &topology,
        &reconcile::Policy::default(),
    )?;
    stage("reconciliation");
    let (trial_mesh, mesh_error) = if options.mesh {
        match mesh::build_partial(
            &topology,
            &mesh::Policy {
                boundary_spacing: profile.element_size,
                maximum_area: profile.element_size * profile.element_size * 2.,
                minimum_angle_degrees: 20.,
                maximum_added_vertices_per_surface: 10000,
            },
        ) {
            Ok(m) => (Some(m), None),
            Err(error) => (None, Some(error.to_string())),
        }
    } else {
        (None, None)
    };
    if options.mesh {
        stage("mesh");
    }
    Ok(Output {
        profile: profile.clone(),
        constraint_graph: graph::Graph::from_frame(&frame_report),
        frame: frame_report,
        topology,
        reconciliation,
        axis_recognition: axes,
        plane_recognition: plane_report,
        relaxation: Relaxation { frame: relaxation },
        rigid_links: RigidLinks {
            removed: rigid_links.len(),
            elements: rigid_links,
        },
        mesh: trial_mesh,
        mesh_error,
    })
}

fn frame_policy(profile: &Profile, over_constrained_panels: bool) -> frame::Policy {
    frame::Policy {
        up: [0., 0., 1.],
        angle: 0.02,
        maximum_movement: 0.15,
        relative_movement: 0.05,
        minimum_length: 0.03,
        residual_tolerance: 1e-7,
        iterations: profile.iterations,
        // Nearly parallel walls through common nodes may be joined into
        // one panel within the tolerance for aligning walls in one line.
        panel_tolerance: profile.stack_offset,
        // A short bar lying in a wall or slab follows its plane.
        geotechnical: true,
        over_constrained_panels,
    }
}

/// Solve the frame with the relaxation ladder (a frame the default rules
/// cannot satisfy is solved again with panels merged at over-constrained
/// nodes), through the cache when one is given.
fn solve_frame(
    input: &Path,
    mesh_data: &MeshData,
    axes: &recognize::Report,
    plane_report: &planes::Report,
    profile: &Profile,
    options: &Options,
    rigid_links: &[u32],
) -> Result<(frame::Report, usize), Error> {
    // Geotechnical mode also closes source gaps across planes (a wall end a
    // few centimetres from the walls it abuts) by moving planes in the solve.
    let gap_tolerance = if profile.close_offset_gaps {
        profile.gap_closure
    } else {
        0.
    };
    let solve = |policy: &frame::Policy| {
        // A finite residual with no movement/axis failure is retried with
        // doubling iteration budgets: 1, 2, 4 and 8 times the base (a
        // curved-wall model needed about 4800 iterations).
        frame::solve_closing_gaps(
            mesh_data,
            axes,
            plane_report,
            policy,
            4,
            gap_tolerance,
            0.001,
        )
    };
    let ladder = || -> Result<(frame::Report, usize), Error> {
        let result = solve(&frame_policy(profile, false))?;
        if !result.accepted {
            let relaxed = solve(&frame_policy(profile, true))?;
            if relaxed.accepted {
                return Ok((relaxed, 1));
            }
        }
        Ok((result, 0))
    };
    let Some(path) = &options.frame_cache else {
        return ladder();
    };
    // The key covers what the frame depends on besides the code: the input
    // text, the policies and the gap tolerance.
    use std::hash::{Hash, Hasher};
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    std::fs::read(input)?.hash(&mut hasher);
    serde_json::to_string(&frame_policy(profile, false))?.hash(&mut hasher);
    gap_tolerance.to_bits().hash(&mut hasher);
    // Kept from the former mode switch so that existing caches stay valid.
    false.hash(&mut hasher);
    if !rigid_links.is_empty() {
        rigid_links.hash(&mut hasher);
    }
    let key = format!("{:016x}", hasher.finish());
    let cached = std::fs::read(path).ok().and_then(|bytes| {
        let value: serde_json::Value = serde_json::from_slice(&bytes).ok()?;
        (value["key"] == key.as_str()).then_some(())?;
        let frame = serde_json::from_value(value["frame"].clone()).ok()?;
        Some((frame, value["relaxation"].as_u64()? as usize))
    });
    if let Some(hit) = cached {
        return Ok(hit);
    }
    let (frame_report, relaxation) = ladder()?;
    let value = serde_json::json!({
        "key": key,
        "relaxation": relaxation,
        "frame": frame_report,
    });
    std::fs::write(path, serde_json::to_vec(&value)?)?;
    Ok((frame_report, relaxation))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn plaxis_profile_is_the_default_and_valid() {
        assert_eq!(Profile::default(), Profile::plaxis());
        assert!(Profile::plaxis().validate().is_ok());
        let mut bad = Profile::plaxis();
        bad.gap_closure = f64::NAN;
        assert!(bad.validate().is_err());
        bad = Profile::plaxis();
        bad.element_size = 0.;
        assert!(bad.validate().is_err());
    }
}
