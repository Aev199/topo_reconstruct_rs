//! Deterministic planning of geometry reconciliation operations.
//!
//! This layer deliberately separates a discrete topological decision from a
//! continuous coordinate solve.  It does not silently edit source geometry:
//! it builds a finite set of auditable hypotheses, selects the best feasible
//! one for each problem within a connected conflict component, and leaves
//! unsupported cases unresolved.  Applying a selected operation will be added
//! only behind the same validation gates.

use super::{assembly, frame};
use crate::input::MeshData;
use glam::DVec3;
use serde::Serialize;
use std::collections::BTreeSet;

#[derive(Debug, Clone, Serialize)]
pub struct Policy {
    /// Numerical comparison tolerance, not a repair budget.
    pub precision: f64,
    /// Area/length² ratio below which a hole is treated as geometrically
    /// degenerate for hypothesis generation.  It is configurable because the
    /// same ratio is not appropriate for every solver profile.
    pub hole_degeneracy_ratio: f64,
    /// Maximum movement available to a future applied repair.
    pub maximum_repair_movement: f64,
}

impl Default for Policy {
    fn default() -> Self {
        Self {
            precision: 1e-8,
            hole_degeneracy_ratio: 1e-4,
            maximum_repair_movement: 0.15,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum ProblemKind {
    DegenerateHole,
    NonCollinearSharedAnchor,
    UnavailableSharedAnchor,
    ContinuousConstraintResidual,
    ContinuousMovementBudget,
    ContinuousAxisFailure,
    ContinuousSolveFailure,
    InvalidSurfaceContour,
    OtherSurfaceTopology,
    OtherAxisAssembly,
}

#[derive(Debug, Clone, Serialize)]
pub struct Evidence {
    pub local_scale: f64,
    pub residual: Option<f64>,
    pub movement: Option<f64>,
    pub movement_budget: Option<f64>,
    pub original_length: Option<f64>,
    pub candidate_length: Option<f64>,
    pub normalized_hole_area: Option<f64>,
    pub boundary_count: Option<usize>,
    pub direct_bar_degree: Option<usize>,
    pub direct_axis_count: Option<usize>,
    pub structural_joint: bool,
    pub depends_on_unbuilt_surface: bool,
}

#[derive(Debug, Clone, Serialize)]
pub struct Problem {
    pub id: usize,
    pub component: usize,
    pub kind: ProblemKind,
    pub patch: Option<usize>,
    pub axis: Option<usize>,
    pub source_elements: Vec<u32>,
    pub source_nodes: Vec<u32>,
    pub evidence: Evidence,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum OperationKind {
    SplitAxisAtStructuralJoint,
    PreserveHole,
    CloseHole,
    DeferUntilSurfaceResolved,
    RetryContinuousSolve,
    LeaveUnresolved,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum DecisionStatus {
    Planned,
    Deferred,
    Unresolved,
}

#[derive(Debug, Clone, Copy, Serialize)]
pub struct Score {
    pub hard_violations: usize,
    pub semantic_loss: usize,
    pub normalized_movement: f64,
    pub topology_edits: usize,
    pub mesh_penalty: usize,
}

#[derive(Debug, Clone, Serialize)]
pub struct Hypothesis {
    pub operation: OperationKind,
    pub feasible: bool,
    pub score: Score,
    pub rationale: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct Decision {
    pub problem_id: usize,
    pub component: usize,
    pub status: DecisionStatus,
    pub selected_operation: OperationKind,
    pub hypotheses: Vec<Hypothesis>,
    pub rationale: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct Component {
    pub id: usize,
    pub problem_ids: Vec<usize>,
    pub source_elements: Vec<u32>,
    pub source_nodes: Vec<u32>,
}

#[derive(Debug, Clone, Serialize)]
pub struct Report {
    pub policy: Policy,
    pub problems: Vec<Problem>,
    pub components: Vec<Component>,
    pub decisions: Vec<Decision>,
    pub unresolved_problem_ids: Vec<usize>,
    /// A safe operation exists, but it depends on a later prerequisite such
    /// as constrained-hole meshing or a completed neighboring surface.
    pub deferred_problem_ids: Vec<usize>,
    pub geometry_changed: bool,
    pub export_ready: bool,
}

/// A small transaction primitive for future geometry application.
///
/// Candidate edits are evaluated on a clone and become visible only after an
/// explicit commit.  The source value is never modified in place.
#[derive(Debug, Clone)]
pub struct Transaction<T: Clone> {
    original: T,
    working: T,
}

impl<T: Clone> Transaction<T> {
    pub fn new(original: T) -> Self {
        Self {
            working: original.clone(),
            original,
        }
    }

    pub fn working(&self) -> &T {
        &self.working
    }

    pub fn working_mut(&mut self) -> &mut T {
        &mut self.working
    }

    pub fn rollback(&mut self) {
        self.working = self.original.clone();
    }

    pub fn commit(self) -> T {
        self.working
    }
}

pub fn solve(
    mesh: &MeshData,
    source: &frame::Report,
    topology: &assembly::Report,
    policy: &Policy,
) -> Result<Report, &'static str> {
    validate_policy(policy)?;
    if source.node_ids.len() != source.candidate_points.len() {
        return Err("reconciliation requires complete frame points");
    }

    let mut problems = collect_problems(mesh, source, topology, policy);
    let components = make_components(&mut problems);
    let mut decisions = Vec::new();
    let mut unresolved_problem_ids = Vec::new();
    let mut deferred_problem_ids = Vec::new();

    for problem in &problems {
        let hypotheses = hypotheses(problem);
        let selected = choose(&hypotheses);
        let status = match selected.operation {
            OperationKind::LeaveUnresolved => DecisionStatus::Unresolved,
            OperationKind::DeferUntilSurfaceResolved
            | OperationKind::PreserveHole
            | OperationKind::RetryContinuousSolve => DecisionStatus::Deferred,
            OperationKind::SplitAxisAtStructuralJoint | OperationKind::CloseHole => {
                DecisionStatus::Planned
            }
        };
        if status == DecisionStatus::Unresolved {
            unresolved_problem_ids.push(problem.id);
        }
        if status == DecisionStatus::Deferred {
            deferred_problem_ids.push(problem.id);
        }
        decisions.push(Decision {
            problem_id: problem.id,
            component: problem.component,
            status,
            selected_operation: selected.operation,
            rationale: selected.rationale.clone(),
            hypotheses,
        });
    }

    Ok(Report {
        policy: policy.clone(),
        problems,
        components,
        decisions,
        unresolved_problem_ids,
        deferred_problem_ids,
        // This first layer is a planner.  No source or candidate coordinate is
        // edited until an operation-specific applier and global validator are
        // present.
        geometry_changed: false,
        export_ready: false,
    })
}

fn validate_policy(policy: &Policy) -> Result<(), &'static str> {
    if !policy.precision.is_finite()
        || policy.precision <= 0.0
        || !policy.hole_degeneracy_ratio.is_finite()
        || policy.hole_degeneracy_ratio <= 0.0
        || !policy.maximum_repair_movement.is_finite()
        || policy.maximum_repair_movement <= 0.0
    {
        return Err("invalid reconciliation policy");
    }
    Ok(())
}

fn collect_problems(
    mesh: &MeshData,
    source: &frame::Report,
    topology: &assembly::Report,
    policy: &Policy,
) -> Vec<Problem> {
    let mut problems = Vec::new();
    let unbuilt_nodes: BTreeSet<u32> = topology
        .issues
        .iter()
        .flat_map(|issue| issue.boundary_source_nodes.iter().flatten().copied())
        .collect();

    for issue in &topology.issues {
        let source_nodes = unique(
            issue
                .boundary_source_nodes
                .iter()
                .flatten()
                .copied()
                .collect(),
        );
        let hole = issue.boundary_source_nodes.get(1);
        let normalized_hole_area = hole.and_then(|ring| {
            let points: Vec<_> = ring
                .iter()
                .filter_map(|id| frame_point(source, *id))
                .collect();
            if points.len() < 3 {
                return None;
            }
            let scale = polygon_scale(&points);
            (scale > policy.precision).then(|| polygon_area(&points) / (scale * scale))
        });
        let hole_is_within_policy = normalized_hole_area
            .map(|ratio| ratio.is_finite() && ratio <= policy.hole_degeneracy_ratio)
            .unwrap_or(true);
        let kind = if issue.reason == "degenerate_hole_after_closure" && hole_is_within_policy {
            ProblemKind::DegenerateHole
        } else if issue.reason.starts_with("contour_") {
            ProblemKind::InvalidSurfaceContour
        } else {
            ProblemKind::OtherSurfaceTopology
        };
        problems.push(Problem {
            id: problems.len(),
            component: 0,
            kind,
            patch: Some(issue.patch),
            axis: None,
            source_elements: unique(issue.source_elements.clone()),
            source_nodes: source_nodes.clone(),
            evidence: Evidence {
                local_scale: local_scale(mesh, &source_nodes, &issue.source_elements),
                residual: None,
                movement: None,
                movement_budget: None,
                original_length: None,
                candidate_length: None,
                normalized_hole_area,
                boundary_count: Some(issue.boundary_source_nodes.len()),
                direct_bar_degree: None,
                direct_axis_count: None,
                structural_joint: false,
                depends_on_unbuilt_surface: source_nodes
                    .iter()
                    .any(|node| unbuilt_nodes.contains(node)),
            },
        });
    }

    for issue in &topology.axis_assembly.issues {
        let source_nodes = unique(issue.source_nodes.clone());
        let direct_bar_degree = source_nodes
            .iter()
            .map(|node| incident_bar_count(mesh, *node))
            .max();
        let direct_axis_count = source_nodes
            .iter()
            .map(|node| incident_axis_count(source, *node))
            .max();
        let structural_joint = direct_bar_degree.is_some_and(|degree| degree != 2)
            || direct_axis_count.is_some_and(|count| count > 1);
        let depends_on_unbuilt_surface =
            source_nodes.iter().any(|node| unbuilt_nodes.contains(node));
        let kind = match issue.reason.as_str() {
            "shared_anchors_not_collinear" => ProblemKind::NonCollinearSharedAnchor,
            "unavailable_shared_anchor" => ProblemKind::UnavailableSharedAnchor,
            _ => ProblemKind::OtherAxisAssembly,
        };
        problems.push(Problem {
            id: problems.len(),
            component: 0,
            kind,
            patch: None,
            axis: Some(issue.source_axis),
            source_elements: unique(issue.source_elements.clone()),
            source_nodes: source_nodes.clone(),
            evidence: Evidence {
                local_scale: local_scale(mesh, &source_nodes, &issue.source_elements),
                residual: None,
                movement: None,
                movement_budget: None,
                original_length: None,
                candidate_length: None,
                normalized_hole_area: None,
                boundary_count: None,
                direct_bar_degree,
                direct_axis_count,
                structural_joint,
                depends_on_unbuilt_surface,
            },
        });
    }
    collect_frame_problems(mesh, source, &unbuilt_nodes, &mut problems);
    problems
}

fn collect_frame_problems(
    mesh: &MeshData,
    source: &frame::Report,
    unbuilt_nodes: &BTreeSet<u32>,
    problems: &mut Vec<Problem>,
) {
    let first_frame_problem = problems.len();
    for failure in &source.largest_constraint_failures {
        let (axis, patch, source_nodes) = match &failure.origin {
            frame::ConstraintOrigin::AxisAnchor { axis, node_id, .. } => {
                (Some(*axis), None, vec![*node_id])
            }
            frame::ConstraintOrigin::ShortAxisVector { axis, .. }
            | frame::ConstraintOrigin::AxisDirection { axis } => {
                (Some(*axis), None, axis_source_nodes(source, *axis))
            }
            frame::ConstraintOrigin::PlaneIncidence { plane, node_id } => {
                (None, Some(*plane), vec![*node_id])
            }
        };
        let source_elements = axis
            .and_then(|i| source.axes.get(i))
            .map(axis_source_elements)
            .or_else(|| {
                patch.and_then(|i| {
                    source
                        .surfaces
                        .get(i)
                        .map(|surface| unique(surface.source_elements.clone()))
                })
            })
            .unwrap_or_default();
        let structural_joint = axis.is_some_and(|_| {
            source_nodes
                .iter()
                .map(|node| incident_bar_count(mesh, *node))
                .any(|degree| degree != 2)
        });
        problems.push(Problem {
            id: problems.len(),
            component: 0,
            kind: ProblemKind::ContinuousConstraintResidual,
            patch,
            axis,
            source_elements: source_elements.clone(),
            source_nodes: source_nodes.clone(),
            evidence: Evidence {
                local_scale: local_scale(mesh, &source_nodes, &source_elements),
                residual: Some(failure.residual.abs()),
                movement: None,
                movement_budget: None,
                original_length: None,
                candidate_length: None,
                normalized_hole_area: None,
                boundary_count: None,
                direct_bar_degree: axis.and_then(|_| {
                    source_nodes
                        .iter()
                        .map(|node| incident_bar_count(mesh, *node))
                        .max()
                }),
                direct_axis_count: axis.and_then(|_| {
                    source_nodes
                        .iter()
                        .map(|node| incident_axis_count(source, *node))
                        .max()
                }),
                structural_joint,
                depends_on_unbuilt_surface: source_nodes
                    .iter()
                    .any(|node| unbuilt_nodes.contains(node)),
            },
        });
    }

    let mut reported_movement_nodes = BTreeSet::new();
    for failure in &source.movement_failures {
        reported_movement_nodes.insert(failure.node_id);
        let source_nodes = vec![failure.node_id];
        let source_elements = incident_element_ids(mesh, failure.node_id);
        problems.push(Problem {
            id: problems.len(),
            component: 0,
            kind: ProblemKind::ContinuousMovementBudget,
            patch: None,
            axis: None,
            source_elements: source_elements.clone(),
            source_nodes: source_nodes.clone(),
            evidence: Evidence {
                local_scale: local_scale(mesh, &source_nodes, &source_elements),
                residual: None,
                movement: Some(failure.movement),
                movement_budget: Some(failure.budget),
                original_length: None,
                candidate_length: None,
                normalized_hole_area: None,
                boundary_count: None,
                direct_bar_degree: Some(incident_bar_count(mesh, failure.node_id)),
                direct_axis_count: Some(incident_axis_count(source, failure.node_id)),
                structural_joint: incident_bar_count(mesh, failure.node_id) != 2,
                depends_on_unbuilt_surface: false,
            },
        });
    }
    for &node in &source.candidate_over_budget_node_ids {
        if reported_movement_nodes.contains(&node) {
            continue;
        }
        let source_nodes = vec![node];
        let source_elements = incident_element_ids(mesh, node);
        problems.push(Problem {
            id: problems.len(),
            component: 0,
            kind: ProblemKind::ContinuousMovementBudget,
            patch: None,
            axis: None,
            source_elements: source_elements.clone(),
            source_nodes: source_nodes.clone(),
            evidence: Evidence {
                local_scale: local_scale(mesh, &source_nodes, &source_elements),
                residual: None,
                movement: None,
                movement_budget: None,
                original_length: None,
                candidate_length: None,
                normalized_hole_area: None,
                boundary_count: None,
                direct_bar_degree: Some(incident_bar_count(mesh, node)),
                direct_axis_count: Some(incident_axis_count(source, node)),
                structural_joint: incident_bar_count(mesh, node) != 2,
                depends_on_unbuilt_surface: false,
            },
        });
    }

    for failure in &source.axis_failures {
        let source_nodes = axis_source_nodes(source, failure.axis);
        let source_elements = unique(failure.source_elements.clone());
        problems.push(Problem {
            id: problems.len(),
            component: 0,
            kind: ProblemKind::ContinuousAxisFailure,
            patch: None,
            axis: Some(failure.axis),
            source_elements: source_elements.clone(),
            source_nodes: source_nodes.clone(),
            evidence: Evidence {
                local_scale: local_scale(mesh, &source_nodes, &source_elements),
                residual: None,
                movement: None,
                movement_budget: None,
                original_length: Some(failure.original_length),
                candidate_length: Some(failure.candidate_length),
                normalized_hole_area: None,
                boundary_count: None,
                direct_bar_degree: source_nodes
                    .iter()
                    .map(|node| incident_bar_count(mesh, *node))
                    .max(),
                direct_axis_count: source_nodes
                    .iter()
                    .map(|node| incident_axis_count(source, *node))
                    .max(),
                structural_joint: source_nodes
                    .iter()
                    .map(|node| incident_bar_count(mesh, *node))
                    .any(|degree| degree != 2),
                depends_on_unbuilt_surface: source_nodes
                    .iter()
                    .any(|node| unbuilt_nodes.contains(node)),
            },
        });
    }

    if source.violating_equations > source.largest_constraint_failures.len() {
        problems.push(Problem {
            id: problems.len(),
            component: 0,
            kind: ProblemKind::ContinuousConstraintResidual,
            patch: None,
            axis: None,
            source_elements: vec![],
            source_nodes: vec![],
            evidence: Evidence {
                local_scale: 1.0,
                residual: Some(source.candidate_max_residual),
                movement: None,
                movement_budget: None,
                original_length: None,
                candidate_length: None,
                normalized_hole_area: None,
                boundary_count: None,
                direct_bar_degree: None,
                direct_axis_count: None,
                structural_joint: false,
                depends_on_unbuilt_surface: false,
            },
        });
    }
    if !source.accepted && problems.len() == first_frame_problem {
        problems.push(Problem {
            id: problems.len(),
            component: 0,
            kind: ProblemKind::ContinuousSolveFailure,
            patch: None,
            axis: None,
            source_elements: vec![],
            source_nodes: vec![],
            evidence: Evidence {
                local_scale: 1.0,
                residual: Some(source.candidate_max_residual),
                movement: None,
                movement_budget: None,
                original_length: None,
                candidate_length: None,
                normalized_hole_area: None,
                boundary_count: None,
                direct_bar_degree: None,
                direct_axis_count: None,
                structural_joint: false,
                depends_on_unbuilt_surface: false,
            },
        });
    }
}

fn axis_source_nodes(source: &frame::Report, axis: usize) -> Vec<u32> {
    source
        .axes
        .get(axis)
        .map(|a| {
            unique(
                a.anchors
                    .iter()
                    .filter_map(|anchor| source.node_ids.get(anchor.node).copied())
                    .collect(),
            )
        })
        .unwrap_or_default()
}

fn axis_source_elements(axis: &frame::Axis) -> Vec<u32> {
    unique(axis.spans.iter().map(|span| span.element).collect())
}

fn incident_element_ids(mesh: &MeshData, node: u32) -> Vec<u32> {
    unique(
        mesh.elements
            .iter()
            .filter(|element| element.nodes.contains(&node))
            .map(|element| element.id)
            .collect(),
    )
}

fn hypotheses(problem: &Problem) -> Vec<Hypothesis> {
    match problem.kind {
        ProblemKind::DegenerateHole => vec![
            Hypothesis {
                operation: OperationKind::PreserveHole,
                feasible: true,
                score: Score {
                    hard_violations: 0,
                    semantic_loss: 0,
                    normalized_movement: 0.0,
                    topology_edits: 1,
                    // A constrained-hole mesher is still required before the
                    // result can become mesh-ready.
                    mesh_penalty: 1,
                },
                rationale:
                    "the opening collapses under closure; preserve it and defer constrained meshing"
                        .into(),
            },
            Hypothesis {
                operation: OperationKind::CloseHole,
                feasible: false,
                score: Score {
                    hard_violations: 1,
                    semantic_loss: 1,
                    normalized_movement: 0.0,
                    topology_edits: 1,
                    mesh_penalty: 0,
                },
                rationale: "rejected because closure would erase or collapse a source opening"
                    .into(),
            },
        ],
        ProblemKind::NonCollinearSharedAnchor => {
            let mut result = vec![Hypothesis {
                operation: OperationKind::SplitAxisAtStructuralJoint,
                feasible: problem.evidence.structural_joint,
                score: Score {
                    hard_violations: usize::from(!problem.evidence.structural_joint),
                    semantic_loss: 0,
                    normalized_movement: 0.0,
                    topology_edits: 1,
                    mesh_penalty: 0,
                },
                rationale: if problem.evidence.structural_joint {
                    "incident constructions form a direct structural joint; keep one shared vertex and separate constructive segments".into()
                } else {
                    "not selected because direct-joint evidence is insufficient".into()
                },
            }];
            result.push(Hypothesis {
                operation: OperationKind::LeaveUnresolved,
                feasible: true,
                score: Score {
                    hard_violations: 0,
                    semantic_loss: 1,
                    normalized_movement: 0.0,
                    topology_edits: 0,
                    mesh_penalty: 2,
                },
                rationale: "retain the conflict explicitly if a structural split cannot be proven"
                    .into(),
            });
            result
        }
        ProblemKind::UnavailableSharedAnchor => {
            let mut result = vec![Hypothesis {
                operation: OperationKind::DeferUntilSurfaceResolved,
                feasible: problem.evidence.depends_on_unbuilt_surface,
                score: Score {
                    hard_violations: usize::from(!problem.evidence.depends_on_unbuilt_surface),
                    semantic_loss: 0,
                    normalized_movement: 0.0,
                    topology_edits: 0,
                    mesh_penalty: 1,
                },
                rationale: if problem.evidence.depends_on_unbuilt_surface {
                    "keep the bar provenance and retry after its owning surface component is resolved".into()
                } else {
                    "no explicit surface dependency was found".into()
                },
            }];
            result.push(Hypothesis {
                operation: OperationKind::LeaveUnresolved,
                feasible: true,
                score: Score {
                    hard_violations: 0,
                    semantic_loss: 1,
                    normalized_movement: 0.0,
                    topology_edits: 0,
                    mesh_penalty: 2,
                },
                rationale: "do not drop or invent the unavailable anchor".into(),
            });
            result
        }
        ProblemKind::ContinuousConstraintResidual => vec![
            Hypothesis {
                operation: OperationKind::RetryContinuousSolve,
                feasible: problem
                    .evidence
                    .residual
                    .is_some_and(|residual| residual.is_finite()),
                score: Score {
                    hard_violations: 0,
                    semantic_loss: 0,
                    normalized_movement: 0.0,
                    topology_edits: 0,
                    mesh_penalty: 1,
                },
                rationale:
                    "the residual is finite; retry the continuous solve with its strict gate".into(),
            },
            Hypothesis {
                operation: OperationKind::LeaveUnresolved,
                feasible: true,
                score: Score {
                    hard_violations: 0,
                    semantic_loss: 1,
                    normalized_movement: 0.0,
                    topology_edits: 0,
                    mesh_penalty: 2,
                },
                rationale: "retain the numerical conflict if a retry cannot satisfy the tolerance"
                    .into(),
            },
        ],
        ProblemKind::ContinuousMovementBudget
        | ProblemKind::ContinuousAxisFailure
        | ProblemKind::ContinuousSolveFailure => {
            vec![Hypothesis {
                operation: OperationKind::LeaveUnresolved,
                feasible: true,
                score: Score {
                    hard_violations: 0,
                    semantic_loss: 0,
                    normalized_movement: 0.0,
                    topology_edits: 0,
                    mesh_penalty: 2,
                },
                rationale: "do not exceed the movement/length gate or invent a repair".into(),
            }]
        }
        ProblemKind::InvalidSurfaceContour
        | ProblemKind::OtherSurfaceTopology
        | ProblemKind::OtherAxisAssembly => vec![Hypothesis {
            operation: OperationKind::LeaveUnresolved,
            feasible: true,
            score: Score {
                hard_violations: 0,
                semantic_loss: 0,
                normalized_movement: 0.0,
                topology_edits: 0,
                mesh_penalty: 2,
            },
            rationale: "no safe generic operation is available for this conflict yet".into(),
        }],
    }
}

fn choose(hypotheses: &[Hypothesis]) -> Hypothesis {
    hypotheses
        .iter()
        .filter(|candidate| candidate.feasible)
        .min_by(|a, b| {
            score_cmp(&a.score, &b.score)
                .then_with(|| operation_rank(a.operation).cmp(&operation_rank(b.operation)))
        })
        .cloned()
        .unwrap_or_else(|| Hypothesis {
            operation: OperationKind::LeaveUnresolved,
            feasible: true,
            score: Score {
                hard_violations: 0,
                semantic_loss: 1,
                normalized_movement: 0.0,
                topology_edits: 0,
                mesh_penalty: 2,
            },
            rationale: "all generated operations were rejected".into(),
        })
}

fn score_cmp(a: &Score, b: &Score) -> std::cmp::Ordering {
    a.hard_violations
        .cmp(&b.hard_violations)
        .then_with(|| a.semantic_loss.cmp(&b.semantic_loss))
        .then_with(|| a.normalized_movement.total_cmp(&b.normalized_movement))
        .then_with(|| a.topology_edits.cmp(&b.topology_edits))
        .then_with(|| a.mesh_penalty.cmp(&b.mesh_penalty))
}

fn operation_rank(operation: OperationKind) -> usize {
    match operation {
        OperationKind::SplitAxisAtStructuralJoint => 0,
        OperationKind::PreserveHole => 1,
        OperationKind::DeferUntilSurfaceResolved => 2,
        OperationKind::RetryContinuousSolve => 3,
        OperationKind::CloseHole => 4,
        OperationKind::LeaveUnresolved => 5,
    }
}

fn make_components(problems: &mut [Problem]) -> Vec<Component> {
    let mut parent: Vec<_> = (0..problems.len()).collect();
    let mut rank = vec![0_usize; problems.len()];
    fn find(parent: &mut [usize], mut value: usize) -> usize {
        while parent[value] != value {
            let next = parent[value];
            parent[value] = parent[next];
            value = parent[value];
        }
        value
    }
    fn union(parent: &mut [usize], rank: &mut [usize], a: usize, b: usize) {
        let mut left = find(parent, a);
        let mut right = find(parent, b);
        if left == right {
            return;
        }
        // Keep the representative deterministic while retaining near-linear
        // union/find behavior for large conflict sets.
        if rank[left] < rank[right] || (rank[left] == rank[right] && left > right) {
            std::mem::swap(&mut left, &mut right);
        }
        parent[right] = left;
        if rank[left] == rank[right] {
            rank[left] += 1;
        }
    }

    let mut node_owner = std::collections::BTreeMap::<u32, usize>::new();
    let mut element_owner = std::collections::BTreeMap::<u32, usize>::new();
    for (i, problem) in problems.iter().enumerate() {
        for &node in &problem.source_nodes {
            if let Some(&owner) = node_owner.get(&node) {
                union(&mut parent, &mut rank, i, owner);
            } else {
                node_owner.insert(node, i);
            }
        }
        for &element in &problem.source_elements {
            if let Some(&owner) = element_owner.get(&element) {
                union(&mut parent, &mut rank, i, owner);
            } else {
                element_owner.insert(element, i);
            }
        }
    }
    let mut groups = std::collections::BTreeMap::<usize, Vec<usize>>::new();
    for i in 0..problems.len() {
        let root = find(&mut parent, i);
        groups.entry(root).or_default().push(i);
    }
    let mut groups: Vec<_> = groups.into_values().collect();
    groups.sort_by_key(|group| group[0]);
    groups
        .into_iter()
        .enumerate()
        .map(|(component, ids)| {
            let mut source_elements = vec![];
            let mut source_nodes = vec![];
            for &id in &ids {
                problems[id].component = component;
                source_elements.extend(&problems[id].source_elements);
                source_nodes.extend(&problems[id].source_nodes);
            }
            Component {
                id: component,
                problem_ids: ids,
                source_elements: unique(source_elements),
                source_nodes: unique(source_nodes),
            }
        })
        .collect()
}

fn unique(mut values: Vec<u32>) -> Vec<u32> {
    values.sort_unstable();
    values.dedup();
    values
}

fn frame_point(source: &frame::Report, node: u32) -> Option<DVec3> {
    source
        .node_ids
        .iter()
        .position(|id| *id == node)
        .and_then(|i| source.candidate_points.get(i).copied())
        .map(DVec3::from_array)
}

fn local_scale(mesh: &MeshData, nodes: &[u32], source_elements: &[u32]) -> f64 {
    let wanted: BTreeSet<_> = source_elements.iter().copied().collect();
    let mut lengths = vec![];
    for element in &mesh.elements {
        if !wanted.is_empty() && !wanted.contains(&element.id) {
            continue;
        }
        if !element.nodes.iter().any(|node| nodes.contains(node)) {
            continue;
        }
        for pair in element.nodes.windows(2) {
            if let (Some(a), Some(b)) = (mesh.nodes.get(&pair[0]), mesh.nodes.get(&pair[1])) {
                let length = a.distance(*b);
                if length.is_finite() && length > 0.0 {
                    lengths.push(length);
                }
            }
        }
    }
    lengths.sort_by(f64::total_cmp);
    lengths
        .get(lengths.len() / 2)
        .copied()
        .or_else(|| {
            nodes
                .iter()
                .filter_map(|node| mesh.nodes.get(node))
                .flat_map(|a| {
                    nodes
                        .iter()
                        .filter_map(move |node| mesh.nodes.get(node).map(|b| a.distance(*b)))
                })
                .filter(|length| length.is_finite() && *length > 0.0)
                .max_by(f64::total_cmp)
        })
        .unwrap_or(1.0)
}

fn polygon_scale(points: &[DVec3]) -> f64 {
    points
        .iter()
        .flat_map(|a| points.iter().map(move |b| a.distance(*b)))
        .max_by(f64::total_cmp)
        .unwrap_or(0.0)
}

fn polygon_area(points: &[DVec3]) -> f64 {
    let sum = points
        .iter()
        .enumerate()
        .map(|(i, point)| point.cross(points[(i + 1) % points.len()]))
        .sum::<DVec3>();
    0.5 * sum.length()
}

fn incident_bar_count(mesh: &MeshData, node: u32) -> usize {
    mesh.elements
        .iter()
        .filter(|element| element.is_bar() && element.nodes.contains(&node))
        .count()
}

fn incident_axis_count(source: &frame::Report, node: u32) -> usize {
    source
        .axes
        .iter()
        .filter(|axis| {
            axis.anchors
                .iter()
                .any(|anchor| source.node_ids.get(anchor.node) == Some(&node))
        })
        .count()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn problem(kind: ProblemKind, structural_joint: bool, depends: bool) -> Problem {
        Problem {
            id: 0,
            component: 0,
            kind,
            patch: None,
            axis: None,
            source_elements: vec![10],
            source_nodes: vec![20],
            evidence: Evidence {
                local_scale: 1.0,
                residual: None,
                movement: None,
                movement_budget: None,
                original_length: None,
                candidate_length: None,
                normalized_hole_area: Some(1e-5),
                boundary_count: Some(2),
                direct_bar_degree: Some(if structural_joint { 4 } else { 2 }),
                direct_axis_count: Some(if structural_joint { 3 } else { 1 }),
                structural_joint,
                depends_on_unbuilt_surface: depends,
            },
        }
    }

    #[test]
    fn structural_joint_prefers_split_over_forced_collinearity() {
        let selected = choose(&hypotheses(&problem(
            ProblemKind::NonCollinearSharedAnchor,
            true,
            false,
        )));
        assert_eq!(
            selected.operation,
            OperationKind::SplitAxisAtStructuralJoint
        );
        assert!(selected.feasible);
    }

    #[test]
    fn degenerate_hole_is_preserved_not_filled() {
        let candidates = hypotheses(&problem(ProblemKind::DegenerateHole, false, false));
        let selected = choose(&candidates);
        assert_eq!(selected.operation, OperationKind::PreserveHole);
        assert!(candidates.iter().any(|candidate| {
            candidate.operation == OperationKind::CloseHole && !candidate.feasible
        }));
    }

    #[test]
    fn unavailable_anchor_waits_for_its_surface_dependency() {
        let selected = choose(&hypotheses(&problem(
            ProblemKind::UnavailableSharedAnchor,
            false,
            true,
        )));
        assert_eq!(selected.operation, OperationKind::DeferUntilSurfaceResolved);
    }

    #[test]
    fn finite_continuous_residual_plans_a_retry_without_relaxing_tolerance() {
        let mut problem = problem(ProblemKind::ContinuousConstraintResidual, false, false);
        problem.evidence.residual = Some(1e-7);
        let selected = choose(&hypotheses(&problem));
        assert_eq!(selected.operation, OperationKind::RetryContinuousSolve);
        assert!(selected.feasible);
    }

    #[test]
    fn transaction_commit_and_rollback_are_explicit() {
        let mut transaction = Transaction::new(vec![1_u32]);
        transaction.working_mut().push(2);
        transaction.rollback();
        assert_eq!(transaction.working(), &vec![1]);
        transaction.working_mut().push(3);
        assert_eq!(transaction.commit(), vec![1, 3]);
    }
}
