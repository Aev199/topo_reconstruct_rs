//! Joint proposal from immutable source nodes, whole axes and support planes.
//! This is a geometric constraint solve, not a mechanical or meshing model.
use super::{planes, recognize, PlaneFrame};
use crate::input::MeshData;
use glam::DVec3;
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};
pub mod gaps;
mod segments;
mod sliding;

#[derive(Debug, Clone, Serialize)]
pub struct Policy {
    pub up: [f64; 3],
    pub angle: f64,
    pub maximum_movement: f64,
    pub relative_movement: f64,
    pub minimum_length: f64,
    pub residual_tolerance: f64,
    pub iterations: usize,
    /// Two support families sharing a node that their intersection line
    /// lies farther from than this (nearly parallel walls, one kinked) become
    /// one panel if all their nodes lie within this of one plane. Zero
    /// disables.
    pub panel_tolerance: f64,
    /// Geotechnical repairs (PLAXIS geometry; deviation from the source is
    /// allowed): a short axis (kept as a vector) whose ends lie on one
    /// support family keeps its vector only along that family's plane, and
    /// every axis node may move at least the plane distance.
    pub geotechnical: bool,
}
#[derive(Debug, Clone, Serialize)]
pub struct Anchor {
    pub node: usize,
    pub t: f64,
}
#[derive(Debug, Clone, Serialize)]
pub struct Axis {
    /// Geometric segment of a chain; shared endpoints do not imply releases.
    pub constructive_segment: bool,
    pub endpoints: [usize; 2],
    pub anchors: Vec<Anchor>,
    pub spans: Vec<recognize::SourceSpan>,
}
#[derive(Debug, Clone, Serialize)]
pub struct Surface {
    pub plane: PlaneFrame,
    pub nodes: Vec<usize>,
    pub source_elements: Vec<u32>,
    pub stiffness_regions: BTreeMap<u32, Vec<u32>>,
}
#[derive(Debug, Clone, Serialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum ConstraintOrigin {
    AxisAnchor {
        axis: usize,
        node_id: u32,
        component: usize,
    },
    ShortAxisVector {
        axis: usize,
        component: usize,
    },
    AxisDirection {
        axis: usize,
    },
    PlaneIncidence {
        plane: usize,
        node_id: u32,
    },
    /// A node near a plane it misses by a gap (`gaps`).
    VirtualIncidence {
        plane: usize,
        node_id: u32,
    },
}
#[derive(Debug, Clone, Serialize)]
pub struct ConstraintFailure {
    pub origin: ConstraintOrigin,
    pub residual: f64,
}
#[derive(Debug, Clone, Serialize)]
pub struct MovementFailure {
    pub node_id: u32,
    pub movement: f64,
    pub budget: f64,
}
#[derive(Debug, Clone, Serialize)]
pub struct AxisFailure {
    pub axis: usize,
    pub original_length: f64,
    pub candidate_length: f64,
    pub source_elements: Vec<u32>,
}

#[derive(Debug, Clone, Serialize)]
pub struct Report {
    /// Experimental nonlinear solve; original Axis::anchors and spans stay immutable.
    pub sliding_parameters: Option<Vec<Vec<f64>>>,
    pub nonlinear_steps: usize,
    pub candidate_parameters_valid: bool,
    pub policy: Policy,
    pub accepted: bool,
    pub candidate_constraints_satisfied: bool,
    pub violating_equations: usize,
    /// At most 100 largest residuals; indices reference axes/surfaces in this report.
    pub largest_constraint_failures: Vec<ConstraintFailure>,
    pub movement_failures: Vec<MovementFailure>,
    pub axis_failures: Vec<AxisFailure>,
    pub reason: String,
    pub iterations: usize,
    pub candidate_max_residual: f64,
    pub candidate_maximum_movement: f64,
    pub candidate_over_budget_node_ids: Vec<u32>,
    pub maximum_movement: f64,
    pub node_ids: Vec<u32>,
    pub reference_points: Vec<[f64; 3]>,
    /// Unaccepted proposal for topology assembly; never an export-ready model.
    pub candidate_points: Vec<[f64; 3]>,
    pub candidate_planes: Vec<PlaneFrame>,
    pub points: Vec<[f64; 3]>,
    pub axes: Vec<Axis>,
    pub surfaces: Vec<Surface>,
    pub plane_families: Vec<Vec<usize>>,
    /// Nearly parallel wall normals were unified (the solve with the
    /// recognised normals did not converge).
    pub regularized_directions: bool,
    pub equation_count: usize,
    pub short_axis_indices: Vec<usize>,
    /// Movement allowed to the nodes of a short axis (translation only):
    /// the plane distance tolerance.
    pub short_axis_movement: f64,
    /// Gaps closed in the solve by virtual plane incidences.
    pub virtual_incidences: gaps::Report,
    #[serde(skip)]
    pub budget_cache: BudgetCache,
}
/// Lazily computed per-node movement budgets; a clone starts empty, so a
/// copied report whose policy or axes are then changed recomputes them.
#[derive(Debug, Default)]
pub struct BudgetCache(std::sync::OnceLock<Vec<f64>>);
impl Clone for BudgetCache {
    fn clone(&self) -> Self {
        Self::default()
    }
}
impl Report {
    /// Movement budget of every node: the maximum movement, capped by the
    /// budget of each axis through it (`short_axis_cap`). Computed once.
    pub fn movement_budgets(&self) -> &[f64] {
        self.budget_cache.0.get_or_init(|| {
            let mut budgets = vec![self.policy.maximum_movement; self.node_ids.len()];
            for axis in &self.axes {
                let [a, b] = axis
                    .endpoints
                    .map(|e| DVec3::from_array(self.reference_points[e]));
                let cap = short_axis_cap(&self.policy, a.distance(b), self.short_axis_movement);
                for anchor in &axis.anchors {
                    budgets[anchor.node] = budgets[anchor.node].min(cap);
                }
            }
            budgets
        })
    }
}
struct Equation {
    origin: ConstraintOrigin,
    terms: Vec<(usize, f64)>,
    target: f64,
}
impl Equation {
    fn value(&self, x: &[f64]) -> f64 {
        self.terms.iter().map(|&(i, a)| a * x[i]).sum()
    }
    fn residual(&self, x: &[f64]) -> f64 {
        self.value(x) - self.target
    }
}
fn equation(terms: Vec<(usize, f64)>, origin: ConstraintOrigin) -> Equation {
    let mut combined = BTreeMap::<usize, f64>::new();
    for (i, a) in terms {
        *combined.entry(i).or_default() += a;
    }
    Equation {
        target: 0.0,
        origin,
        terms: combined
            .into_iter()
            .filter(|(_, a)| a.abs() > 1e-14)
            .collect(),
    }
}

/// Smallest extent of planar points across their principal direction.
fn narrowest_extent(points: &[[f64; 2]]) -> f64 {
    if points.len() < 2 {
        return 0.;
    }
    let n = points.len() as f64;
    let c = points
        .iter()
        .fold([0., 0.], |c, p| [c[0] + p[0] / n, c[1] + p[1] / n]);
    let (mut sxx, mut syy, mut sxy) = (0., 0., 0.);
    for p in points {
        let (x, y) = (p[0] - c[0], p[1] - c[1]);
        sxx += x * x;
        syy += y * y;
        sxy += x * y;
    }
    // Minor principal direction of the 2x2 covariance.
    let angle = 0.5 * (2. * sxy).atan2(sxx - syy) + std::f64::consts::FRAC_PI_2;
    let (s, co) = angle.sin_cos();
    let (lo, hi) = points
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), p| {
            let t = (p[0] - c[0]) * co + (p[1] - c[1]) * s;
            (lo.min(t), hi.max(t))
        });
    hi - lo
}

/// Movement budget of the nodes of an axis: a fraction of its length, but
/// an axis shorter than the minimum length keeps its vector in the solve
/// (it may only translate, undistorted), so its nodes may move as far as a
/// node may lie off its plane, like the surfaces carrying it. In
/// geotechnical mode every axis node may move that far (a 33 mm bar end on
/// a wall corner 2 mm off the wall plane).
fn short_axis_cap(policy: &Policy, length: f64, plane_distance: f64) -> f64 {
    let relative = policy.relative_movement * length;
    policy
        .maximum_movement
        .min(if length < policy.minimum_length || policy.geotechnical {
            relative.max(plane_distance)
        } else {
            relative
        })
}

/// Merge support families whose common nodes cannot lie on one line of
/// their intersection direction: spread across it by more than the plane
/// distance (nearly parallel walls, e.g. an upper wall kinked by a few
/// degrees over a straight lower wall, both through the nodes of the slab
/// edge; the solve would pull those nodes onto the kink line, ill
/// conditioned). Plane offsets are free in the solve, so common nodes along
/// one line (a corner) are consistent. They become one panel when all nodes
/// of both lie within `tolerance` of one plane (a curved wall becomes larger
/// planar panels); pairs are taken by the smallest merged deviation first.
fn merge_kinked_families(
    mesh: &MeshData,
    report: &planes::Report,
    mut families: Vec<Vec<usize>>,
    tolerance: f64,
) -> Vec<Vec<usize>> {
    if tolerance <= 0. {
        return families;
    }
    let points = |members: &[usize]| -> Vec<DVec3> {
        members
            .iter()
            .flat_map(|&i| report.patches[i].source_nodes.iter().copied())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .map(|n| mesh.nodes[&n])
            .collect()
    };
    let plane = |members: &[usize]| -> Option<(DVec3, DVec3, f64)> {
        let p = points(members);
        let (n, _) = planes::fit(
            &p,
            DVec3::from_array(report.patches[members[0]].plane.normal),
        )?;
        let c = p.iter().copied().sum::<DVec3>() / p.len() as f64;
        let deviation = p.iter().map(|x| n.dot(*x - c).abs()).fold(0., f64::max);
        Some((n, c, deviation))
    };
    loop {
        let mut owner = vec![0; report.patches.len()];
        for (f, members) in families.iter().enumerate() {
            for &i in members {
                owner[i] = f;
            }
        }
        let fitted: Vec<_> = families.iter().map(|m| plane(m)).collect();
        let mut node_families = BTreeMap::<u32, BTreeSet<usize>>::new();
        for (i, p) in report.patches.iter().enumerate() {
            for &n in &p.source_nodes {
                node_families.entry(n).or_default().insert(owner[i]);
            }
        }
        let mut common = BTreeMap::<(usize, usize), Vec<DVec3>>::new();
        for (n, fs) in &node_families {
            let fs: Vec<usize> = fs.iter().copied().collect();
            for (k, &f) in fs.iter().enumerate() {
                for &g in &fs[k + 1..] {
                    common.entry((f, g)).or_default().push(mesh.nodes[n]);
                }
            }
        }
        let kinked: Vec<(usize, usize)> = common
            .into_iter()
            .filter_map(|((f, g), nodes)| {
                let (a, b) = (fitted[f]?, fitted[g]?);
                let u = a.0.cross(b.0);
                let spread = if u.length() <= 1e-12 {
                    f64::INFINITY
                } else {
                    let u = u.normalize();
                    let across: Vec<DVec3> = nodes.iter().map(|p| *p - u * u.dot(*p)).collect();
                    let c = across.iter().copied().sum::<DVec3>() / across.len() as f64;
                    across.iter().map(|q| q.distance(c)).fold(0., f64::max)
                };
                (spread > report.policy.distance).then_some((f, g))
            })
            .collect();
        // A node on three or more families whose planes do not meet within
        // the plane distance of it (a slab bump of a few single tilted
        // triangles at a wall): every pair of them is a candidate.
        let mut kinked = kinked;
        for (n, fs) in &node_families {
            if fs.len() < 3 {
                continue;
            }
            let planes: Vec<(DVec3, DVec3)> = fs
                .iter()
                .filter_map(|&f| fitted[f].map(|(n, c, _)| (n, c)))
                .collect();
            let p = mesh.nodes[n];
            let mut normal = glam::DMat3::ZERO;
            let mut rhs = DVec3::ZERO;
            for (n, c) in &planes {
                normal += glam::DMat3::from_cols(*n * n.x, *n * n.y, *n * n.z);
                rhs += *n * n.dot(*c - p);
            }
            let ridge = 1e-9 * (normal.x_axis.x + normal.y_axis.y + normal.z_axis.z);
            let delta = (normal + glam::DMat3::from_diagonal(DVec3::splat(ridge))).inverse() * rhs;
            let off = planes
                .iter()
                .map(|(n, c)| n.dot(p + delta - *c).abs())
                .fold(0., f64::max);
            if delta.length() > report.policy.distance || off > report.policy.distance {
                // Nearly parallel families (within the recognition angle)
                // are consistent through their free offsets in the solve.
                let fs: Vec<usize> = fs.iter().copied().collect();
                for (k, &f) in fs.iter().enumerate() {
                    for &g in &fs[k + 1..] {
                        let (Some(a), Some(b)) = (fitted[f], fitted[g]) else {
                            continue;
                        };
                        if a.0.dot(b.0).abs() < report.policy.angle.cos() {
                            kinked.push((f, g));
                        }
                    }
                }
            }
        }
        kinked.sort_unstable();
        kinked.dedup();
        let best = kinked
            .into_iter()
            .filter_map(|(f, g)| {
                let mut union: Vec<usize> =
                    families[f].iter().chain(&families[g]).copied().collect();
                union.sort_unstable();
                let (_, _, deviation) = plane(&union)?;
                (deviation <= tolerance).then_some((deviation, f, g))
            })
            .min_by(|x, y| x.0.total_cmp(&y.0).then((x.1, x.2).cmp(&(y.1, y.2))));
        let Some((_, f, g)) = best else {
            return families;
        };
        let moved = std::mem::take(&mut families[g]);
        families[f].extend(moved);
        families[f].sort_unstable();
        families.remove(g);
    }
}

/// Near-coplanar patches may meet at a vertex rather than a full source edge.
/// Validate the entire connected family, never chain unbounded plane drift.
fn support_families(mesh: &MeshData, report: &planes::Report) -> Vec<Vec<usize>> {
    let mut owners = BTreeMap::<u32, Vec<usize>>::new();
    for (i, p) in report.patches.iter().enumerate() {
        for &n in &p.source_nodes {
            owners.entry(n).or_default().push(i);
        }
    }
    // Narrowest in-plane extent of each patch: the normal of a narrow patch
    // (a triangle 0.1 m wide) is known only to about distance / width.
    let widths: Vec<f64> = report
        .patches
        .iter()
        .map(|p| {
            let points: Vec<[f64; 2]> = p
                .source_nodes
                .iter()
                .map(|n| p.plane.project(mesh.nodes[n].to_array()))
                .collect();
            narrowest_extent(&points)
        })
        .collect();
    let within = |x: &planes::Patch, y: &planes::Patch| {
        x.source_nodes
            .iter()
            .all(|n| y.plane.distance(mesh.nodes[n].to_array()).abs() <= report.policy.distance)
    };
    // Whether patch `i` may lie in a plane of normal `n` although its own
    // fitted normal differs: always within the recognition angle, beyond it
    // only as far as its width leaves the normal undetermined.
    let tilt_allowed = |i: usize, n: DVec3| {
        let own = DVec3::from_array(report.patches[i].plane.normal);
        let sin = own.cross(n).length();
        n.dot(own).abs() >= report.policy.angle.cos() || sin * widths[i] <= report.policy.distance
    };
    let mut neighbors = vec![BTreeSet::new(); report.patches.len()];
    for ids in owners.values() {
        for (i, &a) in ids.iter().enumerate() {
            for &b in &ids[i + 1..] {
                let x = &report.patches[a];
                let y = &report.patches[b];
                let (nx, ny) = (
                    DVec3::from_array(x.plane.normal),
                    DVec3::from_array(y.plane.normal),
                );
                let parallel = nx.dot(ny).abs() >= report.policy.angle.cos();
                // A narrow patch whose nodes lie in its neighbour's plane is
                // coplanar with it even if its fitted normal is off.
                if (parallel && within(x, y) && within(y, x))
                    || (tilt_allowed(a, ny) && within(x, y))
                    || (tilt_allowed(b, nx) && within(y, x))
                {
                    neighbors[a].insert(b);
                    neighbors[b].insert(a);
                }
            }
        }
    }
    let mut seen = BTreeSet::new();
    let mut families = vec![];
    for seed in 0..report.patches.len() {
        if seen.contains(&seed) {
            continue;
        }
        let mut stack = vec![seed];
        let mut group = BTreeSet::new();
        while let Some(i) = stack.pop() {
            if group.insert(i) {
                stack.extend(&neighbors[i]);
            }
        }
        seen.extend(group.iter().copied());
        let valid = |set: &BTreeSet<usize>| {
            let first = *set.first().unwrap();
            let nodes: BTreeSet<_> = set
                .iter()
                .flat_map(|&i| report.patches[i].source_nodes.iter().copied())
                .collect();
            let points: Vec<_> = nodes.iter().map(|id| mesh.nodes[id]).collect();
            planes::fit(
                &points,
                DVec3::from_array(report.patches[first].plane.normal),
            )
            .is_some_and(|(n, _)| {
                let center = points.iter().copied().sum::<DVec3>() / points.len() as f64;
                points
                    .iter()
                    .all(|p| n.dot(*p - center).abs() <= report.policy.distance)
                    && set.iter().all(|&i| tilt_allowed(i, n))
            })
        };
        if valid(&group) {
            families.push(group.into_iter().collect());
            continue;
        }
        // Not one plane (a narrow patch bridging two facets of a curved
        // wall): larger patches first, each joins the neighbouring family
        // it fits best (least tilt) if that family stays one plane.
        let size = |i: usize| report.patches[i].source_nodes.len();
        let mut order: Vec<usize> = group.iter().copied().collect();
        order.sort_by(|&a, &b| size(b).cmp(&size(a)).then(a.cmp(&b)));
        let mut member = BTreeMap::<usize, usize>::new();
        let mut split: Vec<BTreeSet<usize>> = vec![];
        for i in order {
            let own = DVec3::from_array(report.patches[i].plane.normal);
            let mut options: Vec<(f64, usize)> = neighbors[i]
                .iter()
                .filter_map(|j| member.get(j).copied())
                .collect::<BTreeSet<_>>()
                .into_iter()
                .map(|f| {
                    let n =
                        DVec3::from_array(report.patches[*split[f].first().unwrap()].plane.normal);
                    (own.cross(n).length(), f)
                })
                .collect();
            options.sort_by(|x, y| x.0.total_cmp(&y.0).then(x.1.cmp(&y.1)));
            let joined = options.into_iter().find(|&(_, f)| {
                let mut trial = split[f].clone();
                trial.insert(i);
                valid(&trial)
            });
            match joined {
                Some((_, f)) => {
                    split[f].insert(i);
                    member.insert(i, f);
                }
                None => {
                    member.insert(i, split.len());
                    split.push(BTreeSet::from([i]));
                }
            }
        }
        families.extend(split.into_iter().map(|f| f.into_iter().collect()));
    }
    families
}

fn project(equations: &[Equation], x: &mut Vec<f64>, limit: usize, tolerance: f64) -> (f64, usize) {
    // LSQR operates on A and A^T directly; avoid normal equations A A^T
    // whose conditioning is squared for nearly dependent plane constraints.
    let norm = |v: &[f64]| v.iter().fold(0.0_f64, |a, b| a.hypot(*b));
    let scales: Vec<_> = equations
        .iter()
        .map(|e| e.terms.iter().map(|(_, a)| a * a).sum::<f64>().sqrt())
        .collect();
    let multiply = |v: &[f64]| -> Vec<f64> {
        equations
            .iter()
            .zip(&scales)
            .map(|(e, s)| e.value(v) / s)
            .collect()
    };
    let transpose = |v: &[f64]| -> Vec<f64> {
        let mut result = vec![0.0; x.len()];
        for ((e, s), v) in equations.iter().zip(&scales).zip(v) {
            for &(i, a) in &e.terms {
                result[i] += a * v / s;
            }
        }
        result
    };
    let mut u: Vec<_> = equations
        .iter()
        .zip(&scales)
        .map(|(e, s)| -e.residual(&x) / s)
        .collect();
    let mut beta = norm(&u);
    if beta > 0. {
        for v in &mut u {
            *v /= beta;
        }
    }
    let mut v = transpose(&u);
    let mut alpha = norm(&v);
    if alpha > 0. {
        for p in &mut v {
            *p /= alpha;
        }
    }
    let mut w = v.clone();
    let mut rho_bar = alpha;
    let mut phi_bar = beta;
    let mut correction = vec![0.0; x.len()];
    let original = x.clone();
    let mut residual = equations
        .iter()
        .map(|e| e.residual(&x).abs())
        .fold(0.0_f64, f64::max);
    let mut iterations = 0;
    for step in 0..limit {
        if residual <= tolerance {
            break;
        }
        let av = multiply(&v);
        for (u, av) in u.iter_mut().zip(av) {
            *u = av - alpha * (*u);
        }
        beta = norm(&u);
        if beta > 0. {
            for p in &mut u {
                *p /= beta;
            }
        }
        let atu = transpose(&u);
        for (v, atu) in v.iter_mut().zip(atu) {
            *v = atu - beta * (*v);
        }
        alpha = norm(&v);
        if alpha > 0. {
            for p in &mut v {
                *p /= alpha;
            }
        }
        let rho = rho_bar.hypot(beta);
        if !rho.is_finite() || rho <= 1e-30 {
            break;
        }
        let c = rho_bar / rho;
        let s = beta / rho;
        let theta = s * alpha;
        rho_bar = -c * alpha;
        let phi = c * phi_bar;
        phi_bar *= s;
        for ((dx, w), v) in correction.iter_mut().zip(&mut w).zip(&v) {
            *dx += (phi / rho) * (*w);
            *w = *v - (theta / rho) * (*w);
        }
        // No mutation of x until closure borrows end; assess true residual.
        let candidate: Vec<_> = original
            .iter()
            .zip(&correction)
            .map(|(a, b)| a + b)
            .collect();
        residual = equations
            .iter()
            .map(|e| e.residual(&candidate).abs())
            .fold(0.0_f64, f64::max);
        iterations = step + 1;
    }
    for ((x, o), d) in x.iter_mut().zip(original).zip(correction) {
        *x = o + d;
    }
    (residual, iterations)
}

pub fn solve(
    mesh: &MeshData,
    axes: &recognize::Report,
    planes: &planes::Report,
    policy: &Policy,
) -> Result<Report, &'static str> {
    solve_impl(mesh, axes, planes, policy, 0, &[], false)
}

/// Retry only a numerically incomplete continuous solve with a bounded LSQR
/// budget. Tolerances, source coordinates and movement budgets are unchanged;
/// structural, length and budget failures are returned immediately.
pub fn solve_with_retry(
    mesh: &MeshData,
    axes: &recognize::Report,
    planes: &planes::Report,
    policy: &Policy,
    maximum_attempts: usize,
) -> Result<Report, &'static str> {
    retry(mesh, axes, planes, policy, maximum_attempts, &[])
}

/// `solve_with_retry`, then close gaps below `tolerance` between source nodes
/// and non-parallel support planes by virtual incidences. Incidences the
/// solve cannot satisfy within tolerance and movement budgets are dropped
/// (largest residuals and over-budget nodes first); if none remain, the
/// plain solve is returned. An accepted solve that collapses a source
/// element edge at an incidence node below `minimum` drops those incidences.
/// Every applied or dropped incidence is reported.
pub fn solve_closing_gaps(
    mesh: &MeshData,
    axes: &recognize::Report,
    planes: &planes::Report,
    policy: &Policy,
    maximum_attempts: usize,
    tolerance: f64,
    minimum: f64,
) -> Result<Report, &'static str> {
    let base = retry(mesh, axes, planes, policy, maximum_attempts, &[])?;
    if tolerance <= 0. || !base.accepted {
        return Ok(base);
    }
    let mut active = gaps::candidates(
        mesh,
        planes,
        &base.node_ids,
        tolerance,
        minimum,
        policy.angle,
    );
    let mut report = gaps::Report {
        tolerance,
        candidates: active.len(),
        ..Default::default()
    };
    for round in 1..=4 {
        if active.is_empty() {
            break;
        }
        report.rounds = round;
        let result = retry(mesh, axes, planes, policy, maximum_attempts, &active)?;
        // An accepted solve must not collapse an element edge at an
        // incidence node (a narrow panel squeezed onto a plane).
        let collapsed = if result.accepted {
            gaps::collapsed_edges(mesh, &active, &result.node_ids, &result.points, minimum)
        } else {
            BTreeSet::new()
        };
        if result.accepted && !collapsed.is_empty() {
            let (keep, drop): (Vec<_>, Vec<_>) = active
                .into_iter()
                .partition(|i| !collapsed.contains(&i.node_id));
            report.dropped.extend(drop);
            active = keep;
            continue;
        }
        if result.accepted {
            report.applied = active;
            let mut result = result;
            result.virtual_incidences = report;
            return Ok(result);
        }
        // Drop what the solve could not satisfy.
        let failing: BTreeSet<(u32, usize)> = result
            .largest_constraint_failures
            .iter()
            .filter_map(|f| match f.origin {
                ConstraintOrigin::VirtualIncidence { plane, node_id } => Some((node_id, plane)),
                _ => None,
            })
            .collect();
        let over: BTreeSet<u32> = result
            .candidate_over_budget_node_ids
            .iter()
            .copied()
            .collect();
        let (keep, drop): (Vec<_>, Vec<_>) = active
            .into_iter()
            .partition(|i| !failing.contains(&(i.node_id, i.patch)) && !over.contains(&i.node_id));
        if drop.is_empty() {
            report.dropped.extend(keep);
            active = vec![];
            break;
        }
        report.dropped.extend(drop);
        active = keep;
    }
    report.dropped.extend(active);
    let mut base = base;
    base.virtual_incidences = report;
    Ok(base)
}

fn retry(
    mesh: &MeshData,
    axes: &recognize::Report,
    planes: &planes::Report,
    policy: &Policy,
    maximum_attempts: usize,
    extra: &[gaps::Incidence],
) -> Result<Report, &'static str> {
    // Nearly parallel wall normals are unified only when the solve with the
    // recognised normals cannot converge: the regularisation moves nodes of
    // slightly tilted walls, which other models' budgets may not allow.
    let plain = retry_with(mesh, axes, planes, policy, maximum_attempts, extra, false)?;
    if plain.accepted {
        return Ok(plain);
    }
    let snapped = retry_with(mesh, axes, planes, policy, maximum_attempts, extra, true)?;
    Ok(if snapped.accepted { snapped } else { plain })
}

fn retry_with(
    mesh: &MeshData,
    axes: &recognize::Report,
    planes: &planes::Report,
    policy: &Policy,
    maximum_attempts: usize,
    extra: &[gaps::Incidence],
    snap: bool,
) -> Result<Report, &'static str> {
    if maximum_attempts == 0 {
        return Err("frame retry limit must be positive");
    }
    let mut attempt_policy = policy.clone();
    let mut best = None;
    for attempt in 0..maximum_attempts {
        let result = solve_impl(mesh, axes, planes, &attempt_policy, 0, extra, snap)?;
        let improved = best.as_ref().is_none_or(|previous: &Report| {
            result.candidate_max_residual < previous.candidate_max_residual
        });
        if improved {
            best = Some(result.clone());
        }
        let retryable = !result.accepted
            && result.candidate_parameters_valid
            && result.movement_failures.is_empty()
            && result.axis_failures.is_empty()
            && result.violating_equations > 0
            && result.candidate_max_residual.is_finite()
            && result.candidate_max_residual > attempt_policy.residual_tolerance;
        if !retryable || !improved || attempt + 1 == maximum_attempts {
            break;
        }
        let Some(iterations) = attempt_policy.iterations.checked_mul(2) else {
            break;
        };
        if iterations <= attempt_policy.iterations {
            break;
        }
        attempt_policy.iterations = iterations;
    }
    best.ok_or("frame solve produced no result")
}

/// Joint line/point/plane proposal with free interior line parameters.
/// No topology repair or automatic acceptance into the assembly pipeline.
pub fn solve_sliding(
    mesh: &MeshData,
    axes: &recognize::Report,
    planes: &planes::Report,
    policy: &Policy,
    maximum_steps: usize,
) -> Result<Report, &'static str> {
    if maximum_steps == 0 {
        return Err("nonlinear step limit must be positive");
    }
    solve_impl(mesh, axes, planes, policy, maximum_steps, &[], false)
}

fn solve_impl(
    mesh: &MeshData,
    axes: &recognize::Report,
    planes: &planes::Report,
    policy: &Policy,
    maximum_steps: usize,
    extra: &[gaps::Incidence],
    snap: bool,
) -> Result<Report, &'static str> {
    let up = DVec3::from_array(policy.up);
    if !up.is_finite()
        || !up.length().is_finite()
        || up.length() == 0.0
        || policy.iterations == 0
        || !policy.angle.is_finite()
        || policy.angle <= 0.0
        || policy.angle >= std::f64::consts::FRAC_PI_4
        || [
            policy.maximum_movement,
            policy.relative_movement,
            policy.minimum_length,
            policy.residual_tolerance,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v <= 0.0)
    {
        return Err("invalid frame policy");
    }
    let up = up.normalize();
    let helper = if up.z.abs() < 0.9 { DVec3::Z } else { DVec3::X };
    let u = helper.cross(up).normalize();
    let v = up.cross(u);
    let ids: BTreeSet<_> = axes
        .axes
        .iter()
        .flat_map(|a| a.anchors.iter().map(|c| c.node))
        .chain(
            planes
                .patches
                .iter()
                .flat_map(|p| p.source_nodes.iter().copied()),
        )
        .collect();
    let node_ids: Vec<_> = ids.into_iter().collect();
    let reference: Option<Vec<_>> = node_ids
        .iter()
        .map(|id| mesh.nodes.get(id).copied().filter(|p| p.is_finite()))
        .collect();
    let reference = reference.ok_or("missing source node")?;
    let map: BTreeMap<_, _> = node_ids
        .iter()
        .enumerate()
        .map(|(i, &id)| (id, i))
        .collect();
    let center = reference.first().copied().unwrap_or(DVec3::ZERO);
    let mut x: Vec<_> = reference
        .iter()
        .flat_map(|p| (*p - center).to_array())
        .collect();
    let base = x.len();
    let mut equations = vec![];
    let mut budgets = vec![policy.maximum_movement; reference.len()];
    let mut new_axes = vec![];
    let mut deferred_short = vec![];
    let segments = segments::split(mesh, axes, planes);
    for (axis_index, (a, constructive_segment)) in segments.iter().enumerate() {
        if a.endpoint_nodes.iter().any(|n| !map.contains_key(n))
            || a.anchors
                .iter()
                .any(|c| !c.t.is_finite() || c.t < 0.0 || c.t > 1.0)
        {
            return Err("invalid axis anchors");
        }
        let ends = [map[&a.endpoint_nodes[0]], map[&a.endpoint_nodes[1]]];
        let d = reference[ends[1]] - reference[ends[0]];
        let length = d.length();
        if !length.is_finite() || length <= policy.residual_tolerance {
            return Err("source axis shorter than minimum length");
        }
        let cap = short_axis_cap(policy, length, planes.policy.distance);
        let mut anchors = vec![];
        for c in &a.anchors {
            let node = map[&c.node];
            budgets[node] = budgets[node].min(cap);
            for k in 0..3 {
                equations.push(equation(
                    vec![
                        (node * 3 + k, 1.),
                        (ends[0] * 3 + k, -(1. - c.t)),
                        (ends[1] * 3 + k, -c.t),
                    ],
                    ConstraintOrigin::AxisAnchor {
                        axis: axis_index,
                        node_id: c.node,
                        component: k,
                    },
                ));
            }
            anchors.push(Anchor { node, t: c.t });
        }
        let cosine = d.normalize().dot(up).abs();
        // An unresolved short feature may translate, but must not be flattened,
        // contracted or assigned a new direction from noisy source coordinates.
        if length < policy.minimum_length {
            if policy.geotechnical {
                deferred_short.push((axis_index, ends, d));
            } else {
                for k in 0..3 {
                    let mut e = equation(
                        vec![(ends[1] * 3 + k, 1.), (ends[0] * 3 + k, -1.)],
                        ConstraintOrigin::ShortAxisVector {
                            axis: axis_index,
                            component: k,
                        },
                    );
                    e.target = d[k];
                    equations.push(e);
                }
            }
        }
        let directions = if *constructive_segment || length < policy.minimum_length {
            vec![]
        } else if cosine >= policy.angle.cos() {
            vec![u, v]
        } else if cosine <= policy.angle.sin() {
            vec![up]
        } else {
            vec![]
        };
        for n in directions {
            equations.push(equation(
                (0..3)
                    .flat_map(|k| [(ends[0] * 3 + k, n[k]), (ends[1] * 3 + k, -n[k])])
                    .collect(),
                ConstraintOrigin::AxisDirection { axis: axis_index },
            ));
        }
        new_axes.push(Axis {
            constructive_segment: *constructive_segment,
            endpoints: ends,
            anchors,
            spans: a.spans.clone(),
        });
    }
    let plane_families = merge_kinked_families(
        mesh,
        planes,
        support_families(mesh, planes),
        policy.panel_tolerance,
    );
    let mut plane_to_family = vec![0; planes.patches.len()];
    let mut normals = vec![];
    let mut pending: Vec<(DVec3, Vec<DVec3>)> = vec![];
    for (family, members) in plane_families.iter().enumerate() {
        let nodes: BTreeSet<_> = members
            .iter()
            .flat_map(|&i| planes.patches[i].source_nodes.iter().copied())
            .collect();
        if nodes.is_empty() {
            return Err("empty support plane");
        }
        let points: Vec<_> = nodes.iter().map(|n| mesh.nodes[n]).collect();
        let old = if members.len() == 1 {
            DVec3::from_array(planes.patches[members[0]].plane.normal)
        } else {
            planes::fit(
                &points,
                DVec3::from_array(planes.patches[members[0]].plane.normal),
            )
            .ok_or("invalid support family")?
            .0
        };
        let cosine = old.dot(up);
        let snapped = if cosine.abs() >= policy.angle.cos() {
            up * cosine.signum()
        } else if cosine.abs() <= policy.angle.sin() {
            (old - up * cosine).normalize()
        } else {
            old
        };
        // Horizontal slabs and vertical walls, if every node of the family
        // stays within the plane distance (a large slightly tilted family,
        // 1 degree over 9 m, would move nodes by 16 cm).
        let mean = points.iter().copied().sum::<DVec3>() / points.len() as f64;
        let normal = if snapped == old
            || points
                .iter()
                .all(|p| snapped.dot(*p - mean).abs() <= planes.policy.distance)
        {
            snapped
        } else {
            old
        };
        pending.push((normal, points));
        for &pi in members {
            plane_to_family[pi] = family;
        }
    }
    // Walls in one direction share one normal: nearly parallel wall normals
    // (within the recognition angle) take the direction of the largest wall
    // of their group, if every node of the family stays within the plane
    // distance tolerance. Otherwise bars between them are inconsistent with
    // their fixed directions (a 1e-4 rad difference over the spacing of two
    // vertical bars leaves a residual the solve cannot remove).
    let vertical = |n: DVec3| n.dot(up).abs() <= 1e-12;
    let mut order: Vec<usize> = (0..pending.len())
        .filter(|&f| snap && vertical(pending[f].0))
        .collect();
    order.sort_by(|&a, &b| {
        pending[b].1.len().cmp(&pending[a].1.len()).then_with(|| {
            let angle = |n: DVec3| n.y.atan2(n.x).rem_euclid(std::f64::consts::PI);
            angle(pending[a].0).total_cmp(&angle(pending[b].0))
        })
    });
    let mut directions: Vec<DVec3> = vec![];
    for f in order {
        let normal = pending[f].0;
        let Some(&d) = directions
            .iter()
            .find(|d| d.dot(normal).abs() >= policy.angle.cos())
        else {
            directions.push(normal);
            continue;
        };
        let snapped = d * d.dot(normal).signum();
        let points = &pending[f].1;
        let mean = points.iter().copied().sum::<DVec3>() / points.len() as f64;
        if points
            .iter()
            .all(|p| snapped.dot(*p - mean).abs() <= planes.policy.distance)
        {
            pending[f].0 = snapped;
        }
    }
    for (normal, points) in pending {
        let origin = points.iter().map(|p| *p - center).sum::<DVec3>() / points.len() as f64;
        x.push(normal.dot(origin));
        normals.push(normal);
    }
    // Short axes whose ends lie on common support families keep their
    // vector only across those families' normals.
    let mut node_families = BTreeMap::<usize, BTreeSet<usize>>::new();
    for (pi, p) in planes.patches.iter().enumerate() {
        for id in &p.source_nodes {
            node_families
                .entry(map[id])
                .or_default()
                .insert(plane_to_family[pi]);
        }
    }
    let mut flattened_lengths = BTreeMap::<usize, f64>::new();
    for (axis_index, ends, d) in deferred_short {
        let empty = BTreeSet::new();
        let [a, b] = ends.map(|e| node_families.get(&e).unwrap_or(&empty));
        let mut basis: Vec<DVec3> = vec![];
        for f in a.intersection(b) {
            let mut n = normals[*f];
            for q in &basis {
                n -= *q * q.dot(n);
            }
            if n.length() > 1e-6 {
                basis.push(n.normalize());
            }
        }
        let fixed = basis.len();
        // Like a long axis, a vertical or horizontal short axis is snapped
        // to its direction class; a horizontal one keeps only its plan
        // length and level, its transverse plan offset (source noise, e.g.
        // a chain wall - slab - wall of 25 mm links) is free.
        let cosine = d.normalize().dot(up).abs();
        let plan = d - up * d.dot(up);
        let (target, candidates) = if cosine >= policy.angle.cos() {
            (up * d.dot(up), vec![DVec3::X, DVec3::Y, DVec3::Z])
        } else if cosine <= policy.angle.sin() {
            (plan, vec![plan.normalize(), up])
        } else {
            (d, vec![DVec3::X, DVec3::Y, DVec3::Z])
        };
        for axis in candidates {
            let mut n = axis;
            for q in &basis {
                n -= *q * q.dot(n);
            }
            if n.length() > 1e-6 {
                basis.push(n.normalize());
            }
        }
        let kept: DVec3 = basis.iter().skip(fixed).map(|n| *n * n.dot(target)).sum();
        // A bar mostly across the plane is not noise in it: keep its vector.
        let (fixed, target) = if kept.length() < 0.5 * d.length() {
            basis = vec![DVec3::X, DVec3::Y, DVec3::Z];
            (0, d)
        } else {
            flattened_lengths.insert(axis_index, kept.length());
            (fixed, target)
        };
        for (k, n) in basis.iter().enumerate().skip(fixed) {
            let mut e = equation(
                (0..3)
                    .flat_map(|c| [(ends[1] * 3 + c, n[c]), (ends[0] * 3 + c, -n[c])])
                    .collect(),
                ConstraintOrigin::ShortAxisVector {
                    axis: axis_index,
                    component: k - fixed,
                },
            );
            e.target = n.dot(target);
            equations.push(e);
        }
    }
    let mut surfaces = vec![];
    for (pi, p) in planes.patches.iter().enumerate() {
        let family = plane_to_family[pi];
        let normal = normals[family];
        let nodes: Vec<_> = p.source_nodes.iter().map(|id| map[id]).collect();
        for &node in &nodes {
            equations.push(equation(
                vec![
                    (node * 3, normal.x),
                    (node * 3 + 1, normal.y),
                    (node * 3 + 2, normal.z),
                    (base + family, -1.),
                ],
                ConstraintOrigin::PlaneIncidence {
                    plane: pi,
                    node_id: node_ids[node],
                },
            ));
        }
        surfaces.push(Surface {
            plane: p.plane.clone(),
            nodes,
            source_elements: p.source_elements.clone(),
            stiffness_regions: p.stiffness_regions.clone(),
        });
    }
    for incidence in extra {
        let (Some(&node), Some(&family)) = (
            map.get(&incidence.node_id),
            plane_to_family.get(incidence.patch),
        ) else {
            continue;
        };
        let normal = normals[family];
        equations.push(equation(
            vec![
                (node * 3, normal.x),
                (node * 3 + 1, normal.y),
                (node * 3 + 2, normal.z),
                (base + family, -1.),
            ],
            ConstraintOrigin::VirtualIncidence {
                plane: incidence.patch,
                node_id: incidence.node_id,
            },
        ));
    }
    equations.retain(|e| !e.terms.is_empty());
    let (residual, iterations, nonlinear_steps, sliding_parameters) = if maximum_steps == 0 {
        let (residual, iterations) = project(
            &equations,
            &mut x,
            policy.iterations,
            policy.residual_tolerance,
        );
        (residual, iterations, 0, None)
    } else {
        sliding::solve(
            &mut equations,
            &mut x,
            &new_axes,
            &reference,
            &map,
            policy,
            maximum_steps,
        )
    };
    let candidate: Vec<_> = (0..reference.len())
        .map(|i| center + DVec3::new(x[i * 3], x[i * 3 + 1], x[i * 3 + 2]))
        .collect();
    let valid_parameters = sliding_parameters.as_ref().is_none_or(|parameters| {
        parameters.iter().zip(&new_axes).all(|(ts, axis)| {
            let mut ordered: Vec<_> = axis.anchors.iter().zip(ts).collect();
            ordered.sort_by(|(a, _), (b, _)| a.t.total_cmp(&b.t));
            ts.iter().all(|t| t.is_finite() && *t >= 0.0 && *t <= 1.0)
                && ordered.windows(2).all(|w| w[0].1 < w[1].1)
        })
    });
    let valid = valid_parameters
        && candidate.iter().all(|p| p.is_finite())
        && new_axes.iter().enumerate().all(|(k, a)| {
            let [i, j] = a.endpoints;
            let d = candidate[j] - candidate[i];
            // A short axis flattened onto its plane keeps its in-plane length.
            let expected = flattened_lengths
                .get(&k)
                .copied()
                .unwrap_or(reference[j].distance(reference[i]));
            // Acceptance still requires convergence (residual within the
            // tolerance); before it, a change below the residual is numerical.
            d.length() + policy.residual_tolerance.max(residual)
                >= policy.minimum_length.min(expected)
                && d.dot(reference[j] - reference[i]) > 0.0
                && (sliding_parameters.is_none()
                    || d.normalize().dot((reference[j] - reference[i]).normalize())
                        >= policy.angle.cos())
        });
    let within_budget = candidate
        .iter()
        .zip(&reference)
        .zip(&budgets)
        .all(|((p, o), cap)| p.distance(*o) <= *cap + 1e-10);
    let accepted = residual <= policy.residual_tolerance && valid && within_budget;
    let mut maximum_movement = 0.0_f64;
    let mut candidate_surfaces = surfaces.clone();
    {
        for (i, p) in candidate_surfaces.iter_mut().enumerate() {
            let origin = DVec3::from_array(p.plane.origin);
            let family = plane_to_family[i];
            let n = normals[family];
            let projected = origin + n * (x[base + family] - n.dot(origin - center));
            p.plane = PlaneFrame::new(projected.to_array(), n.to_array())
                .map_err(|_| "invalid solved plane")?;
        }
    }
    if accepted {
        surfaces = candidate_surfaces.clone();
        maximum_movement = candidate
            .iter()
            .zip(&reference)
            .map(|(a, b)| a.distance(*b))
            .fold(0.0_f64, f64::max);
    }
    let mut failures: Vec<_> = equations
        .iter()
        .filter_map(|e| {
            let residual = e.residual(&x).abs();
            (residual > policy.residual_tolerance).then(|| ConstraintFailure {
                origin: e.origin.clone(),
                residual,
            })
        })
        .collect();
    let violating_equations = failures.len();
    failures.sort_by(|a, b| b.residual.total_cmp(&a.residual));
    failures.truncate(100);
    let movement_failures = candidate
        .iter()
        .zip(&reference)
        .zip(&budgets)
        .enumerate()
        .filter_map(|(i, ((p, o), cap))| {
            let movement = p.distance(*o);
            (movement > *cap + 1e-10).then_some(MovementFailure {
                node_id: node_ids[i],
                movement,
                budget: *cap,
            })
        })
        .collect();
    let axis_failures = new_axes
        .iter()
        .enumerate()
        .filter_map(|(axis, a)| {
            let [i, j] = a.endpoints;
            let old = reference[j] - reference[i];
            let new = candidate[j] - candidate[i];
            // A short axis flattened onto its plane keeps its in-plane length.
            let expected = flattened_lengths
                .get(&axis)
                .copied()
                .unwrap_or(old.length());
            // Before convergence a length change below the residual is
            // numerical, not a collapse (it must not block the retry).
            (new.length() + policy.residual_tolerance.max(residual)
                < policy.minimum_length.min(expected)
                || new.dot(old) <= 0.0
                || (sliding_parameters.is_some()
                    && new.normalize().dot(old.normalize()) < policy.angle.cos()))
            .then_some(AxisFailure {
                axis,
                original_length: old.length(),
                candidate_length: new.length(),
                source_elements: a.spans.iter().map(|s| s.element).collect(),
            })
        })
        .collect();
    Ok(Report {
        sliding_parameters,
        nonlinear_steps,
        candidate_parameters_valid: valid_parameters,
        policy: policy.clone(),
        accepted,
        candidate_constraints_satisfied: residual <= policy.residual_tolerance,
        violating_equations,
        largest_constraint_failures: failures,
        movement_failures,
        axis_failures,
        reason: if accepted {
            "converged"
        } else if !valid_parameters {
            "invalid_anchor_order"
        } else if !within_budget {
            "movement_budget_exceeded"
        } else if !valid {
            "invalid_axis"
        } else {
            "constraints_or_budget_not_satisfied"
        }
        .into(),
        iterations,
        candidate_max_residual: residual,
        candidate_maximum_movement: candidate
            .iter()
            .zip(&reference)
            .map(|(a, b)| a.distance(*b))
            .fold(0.0_f64, f64::max),
        candidate_over_budget_node_ids: candidate
            .iter()
            .zip(&reference)
            .zip(&budgets)
            .enumerate()
            .filter_map(|(i, ((p, o), cap))| (p.distance(*o) > *cap + 1e-10).then_some(node_ids[i]))
            .collect(),
        maximum_movement,
        node_ids: node_ids.clone(),
        reference_points: reference.iter().map(|p| p.to_array()).collect(),
        candidate_points: candidate.iter().map(|p| p.to_array()).collect(),
        candidate_planes: candidate_surfaces.iter().map(|s| s.plane.clone()).collect(),
        points: if accepted {
            candidate.iter().map(|p| p.to_array()).collect()
        } else {
            reference.iter().map(|p| p.to_array()).collect()
        },
        axes: new_axes.clone(),
        surfaces,
        short_axis_indices: new_axes
            .iter()
            .enumerate()
            .filter_map(|(i, a)| {
                (reference[a.endpoints[0]].distance(reference[a.endpoints[1]])
                    < policy.minimum_length)
                    .then_some(i)
            })
            .collect(),
        short_axis_movement: planes.policy.distance,
        plane_families,
        regularized_directions: snap,
        equation_count: equations.len(),
        virtual_incidences: gaps::Report::default(),
        budget_cache: Default::default(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::input::ElementData;
    fn source() -> MeshData {
        let points = [
            [0., 0., 0.02],
            [1., 0., 0.025],
            [2., 0., 0.03],
            [1.01, 0., 2.],
            [0., 2., 0.02],
        ];
        MeshData {
            nodes: points
                .iter()
                .enumerate()
                .map(|(i, p)| (i as u32 + 1, DVec3::from_array(*p)))
                .collect(),
            elements: vec![
                ElementData {
                    id: 1,
                    elem_type: 10,
                    stiff_id: 1,
                    nodes: vec![1, 2],
                },
                ElementData {
                    id: 2,
                    elem_type: 10,
                    stiff_id: 2,
                    nodes: vec![2, 3],
                },
                ElementData {
                    id: 3,
                    elem_type: 10,
                    stiff_id: 3,
                    nodes: vec![2, 4],
                },
                ElementData {
                    id: 4,
                    elem_type: 42,
                    stiff_id: 4,
                    nodes: vec![1, 3, 5],
                },
            ],
        }
    }

    fn column_through_floors() -> MeshData {
        let mut mesh = MeshData::default();
        for (id, z) in (0..=4).enumerate() {
            mesh.nodes.insert(
                id as u32 + 1,
                DVec3::new(if z == 2 { 0.003 } else { 0. }, 0., z as f64),
            );
        }
        for i in 0..4 {
            mesh.elements.push(ElementData {
                id: i + 1,
                elem_type: 10,
                stiff_id: 100 + i,
                nodes: vec![i + 1, i + 2],
            });
        }
        for (i, z) in [0., 2., 4.].into_iter().enumerate() {
            let axis_node = i * 2 + 1;
            let first = 10 + i as u32 * 3;
            mesh.nodes.insert(first, DVec3::new(1., 0., z));
            mesh.nodes.insert(first + 1, DVec3::new(0., 1., z));
            mesh.elements.push(ElementData {
                id: 10 + i as u32,
                elem_type: 42,
                stiff_id: 200 + i as u32,
                nodes: vec![axis_node as u32, first, first + 1],
            });
        }
        mesh
    }
    fn run_sliding(m: &MeshData, p: &Policy, scale: f64) -> Report {
        let axes = recognize::recognize(
            m,
            &recognize::Policy {
                angle: 0.02,
                line_tolerance: 0.01 * scale,
                numerical_precision: 1e-8 * scale,
            },
        )
        .unwrap();
        let planes = planes::recognize(
            m,
            &planes::Policy {
                angle: 0.02,
                distance: 0.01 * scale,
                precision: 1e-8 * scale,
            },
        )
        .unwrap();
        solve_sliding(m, &axes, &planes, p, 12).unwrap()
    }

    #[test]
    fn joint_slides_along_whole_beam_under_scale_rotation_and_renumbering() {
        let baseline = run_sliding(&source(), &policy(), 1.);
        for scale in [0.1, 1., 10.] {
            for transformed in [false, true] {
                let mut m = source();
                let rotation = if transformed {
                    glam::DQuat::from_axis_angle(DVec3::new(1., 2., 3.).normalize(), 0.6)
                } else {
                    glam::DQuat::IDENTITY
                };
                let shift = if transformed {
                    DVec3::new(20., 30., 40.)
                } else {
                    DVec3::ZERO
                };
                let id = |n: u32| if transformed { 100 - n } else { n };
                m.nodes = m
                    .nodes
                    .into_iter()
                    .map(|(n, v)| (id(n), rotation * (v * scale) + shift))
                    .collect();
                for e in &mut m.elements {
                    for n in &mut e.nodes {
                        *n = id(*n);
                    }
                }
                if transformed {
                    m.elements.reverse();
                }
                let mut p = policy();
                p.up = (rotation * DVec3::Z).to_array();
                p.maximum_movement *= scale;
                p.minimum_length *= scale;
                p.residual_tolerance *= scale;
                let r = run_sliding(&m, &p, scale);
                assert!(r.accepted, "{}: {}", r.reason, r.candidate_max_residual);
                for (&n, point) in baseline.node_ids.iter().zip(&baseline.points) {
                    let i = r.node_ids.iter().position(|&v| v == id(n)).unwrap();
                    let expected = rotation * (DVec3::from_array(*point) * scale) + shift;
                    assert!(DVec3::from_array(r.points[i]).distance(expected) < 1e-6 * scale);
                }
                let ts = r.sliding_parameters.as_ref().unwrap();
                assert_eq!(r.axes.len(), 3);
                assert!(r.axes.iter().all(|a| a.spans.len() == 1));
                assert_eq!(r.axes.iter().filter(|a| a.constructive_segment).count(), 2);
                assert!(r.axes.iter().zip(ts).all(|(a, t)| {
                    a.anchors
                        .iter()
                        .zip(t)
                        .all(|(c, t)| (c.t - t).abs() <= 1e-5)
                }));
                for (axis, ts) in r.axes.iter().zip(ts) {
                    let [a, b] = axis.endpoints.map(|i| DVec3::from_array(r.points[i]));
                    for (anchor, t) in axis.anchors.iter().zip(ts) {
                        assert!(
                            DVec3::from_array(r.points[anchor.node]).distance(a + t * (b - a))
                                < 3. * p.residual_tolerance
                        );
                    }
                }
                for surface in &r.surfaces {
                    for &n in &surface.nodes {
                        assert!(
                            surface.plane.distance(r.points[n]).abs() < 3. * p.residual_tolerance
                        );
                    }
                }
                let spans: BTreeMap<_, _> = r
                    .axes
                    .iter()
                    .flat_map(|a| a.spans.iter().map(|s| (s.element, s.stiffness)))
                    .collect();
                assert_eq!(spans, BTreeMap::from([(1, 1), (2, 2), (3, 3)]));
                let beam: Vec<_> = r
                    .axes
                    .iter()
                    .filter(|a| matches!(a.spans[0].element, 1 | 2))
                    .collect();
                assert_eq!(beam.len(), 2);
                assert!(beam.iter().all(|a| a.anchors.len() == 2));
            }
        }
    }

    #[test]
    fn column_through_transverse_floors_splits_only_at_structural_supports() {
        for scale in [0.1, 1., 10.] {
            for transformed in [false, true] {
                let mut m = column_through_floors();
                let rotation = if transformed {
                    glam::DQuat::from_axis_angle(DVec3::new(1., 2., 3.).normalize(), 0.6)
                } else {
                    glam::DQuat::IDENTITY
                };
                let shift = if transformed {
                    DVec3::new(20., 30., 40.)
                } else {
                    DVec3::ZERO
                };
                let node_id = |n: u32| if transformed { 1000 - n } else { n };
                let element_id = |e: u32| if transformed { 2000 - e } else { e };
                m.nodes = m
                    .nodes
                    .into_iter()
                    .map(|(n, p)| (node_id(n), rotation * (p * scale) + shift))
                    .collect();
                for e in &mut m.elements {
                    e.id = element_id(e.id);
                    for n in &mut e.nodes {
                        *n = node_id(*n);
                    }
                    e.nodes.reverse();
                }
                if transformed {
                    m.elements.reverse();
                }
                let recognized = recognize::recognize(
                    &m,
                    &recognize::Policy {
                        angle: 0.02,
                        line_tolerance: 0.01 * scale,
                        numerical_precision: 1e-8 * scale,
                    },
                )
                .unwrap();
                assert_eq!(recognized.axes.len(), 1);
                let supports = planes::recognize(
                    &m,
                    &planes::Policy {
                        angle: 0.02,
                        distance: 0.01 * scale,
                        precision: 1e-8 * scale,
                    },
                )
                .unwrap();
                assert_eq!(supports.patches.len(), 3);
                let pieces = segments::split(&m, &recognized, &supports);
                assert_eq!(pieces.len(), 2);
                assert!(pieces.iter().all(|(_, constructive)| *constructive));
                assert!(pieces.iter().all(|(axis, _)| {
                    axis.spans.len() == 2
                        && axis.anchors.len() == 3
                        && (axis.anchors[1].t - 0.5).abs() < 1e-8
                }));
                assert_eq!(
                    pieces
                        .iter()
                        .flat_map(|(axis, _)| axis.spans.iter().map(|s| s.element))
                        .map(|e| if transformed { 2000 - e } else { e })
                        .collect::<BTreeSet<_>>(),
                    BTreeSet::from([1, 2, 3, 4])
                );
                assert!(pieces
                    .iter()
                    .any(|(axis, _)| { axis.endpoint_nodes.contains(&node_id(3)) }));
                let result = solve(
                    &m,
                    &recognized,
                    &supports,
                    &Policy {
                        up: (rotation * DVec3::Z).to_array(),
                        angle: 0.02,
                        maximum_movement: 0.15 * scale,
                        relative_movement: 0.05,
                        minimum_length: 0.03 * scale,
                        residual_tolerance: 1e-8 * scale,
                        iterations: 2000,
                        panel_tolerance: 0.,
                        geotechnical: false,
                    },
                )
                .unwrap();
                assert!(result.accepted, "{}", result.candidate_max_residual);
                assert_eq!(result.axes.len(), 2);
                assert!(result.axes.iter().all(|axis| axis.spans.len() == 2));
            }
        }
    }

    #[test]
    fn sliding_failure_keeps_reference_geometry_and_short_features() {
        let mut p = policy();
        p.maximum_movement = 1e-6;
        let r = run_sliding(&source(), &p, 1.);
        assert!(!r.accepted);
        assert_eq!(r.points, r.reference_points);
        assert!(!r.movement_failures.is_empty());
        let r = run_sliding(&short_feature(true), &policy(), 1.);
        assert!(!r.accepted);
        assert_eq!(r.points, r.reference_points);
        assert_eq!(r.axes[0].spans[0].element, 10);
    }

    #[test]
    fn sliding_proposal_cannot_enter_assembly_with_old_property_parameters() {
        let mesh = source();
        let r = run_sliding(&mesh, &policy(), 1.);
        assert!(r.accepted);
        let result = super::super::assembly::assemble(
            &mesh,
            &r,
            &super::super::assembly::Policy {
                closure_tolerance: 0.001,
                junction_movement_limit: 0.05,
                precision: 1e-7,
                minimum_edge: 0.001,
            },
        );
        assert!(matches!(
            result,
            Err("sliding frame is proposal-only until parameter transfer is implemented")
        ));
    }

    #[test]
    fn retry_doubles_only_a_numerical_iteration_budget() {
        let mesh = source();
        let axes = recognize::recognize(
            &mesh,
            &recognize::Policy {
                angle: 0.02,
                line_tolerance: 0.01,
                numerical_precision: 1e-8,
            },
        )
        .unwrap();
        let planes = planes::recognize(
            &mesh,
            &planes::Policy {
                angle: 0.02,
                distance: 0.01,
                precision: 1e-8,
            },
        )
        .unwrap();
        let mut policy = policy();
        policy.iterations = 1;
        let first = solve(&mesh, &axes, &planes, &policy).unwrap();
        assert!(!first.accepted);
        let retried = solve_with_retry(&mesh, &axes, &planes, &policy, 12).unwrap();
        assert!(retried.accepted, "{}", retried.candidate_max_residual);
        assert!(retried.iterations > policy.iterations);
        assert_eq!(retried.policy.residual_tolerance, policy.residual_tolerance);
        assert_eq!(retried.policy.maximum_movement, policy.maximum_movement);
        assert_eq!(retried.policy.minimum_length, policy.minimum_length);
    }

    fn policy() -> Policy {
        Policy {
            up: [0., 0., 1.],
            angle: 0.02,
            maximum_movement: 0.15,
            relative_movement: 0.05,
            minimum_length: 0.03,
            residual_tolerance: 1e-8,
            iterations: 2000,
            panel_tolerance: 0.,
            geotechnical: false,
        }
    }
    fn run(m: &MeshData, p: &Policy) -> Report {
        let a = recognize::recognize(
            m,
            &recognize::Policy {
                angle: 0.02,
                line_tolerance: 0.01,
                numerical_precision: 1e-8,
            },
        )
        .unwrap();
        let s = planes::recognize(
            m,
            &planes::Policy {
                angle: 0.02,
                distance: 0.01,
                precision: 1e-8,
            },
        )
        .unwrap();
        solve(m, &a, &s, p).unwrap()
    }
    fn short_feature(both_ends_on_plane: bool) -> MeshData {
        let mut m = MeshData::default();
        for (id, point) in [
            (1, [0., 0., 0.]),
            (2, [1., 0., 0.0002]),
            (3, [0., 1., 0.]),
            (4, [0.006, 0., 0.00001]),
        ] {
            m.nodes.insert(id, DVec3::from_array(point));
        }
        m.elements.push(ElementData {
            id: 10,
            elem_type: 10,
            stiff_id: 77,
            nodes: vec![1, 4],
        });
        m.elements.push(ElementData {
            id: 20,
            elem_type: 42,
            stiff_id: 88,
            nodes: vec![1, 2, 3],
        });
        if both_ends_on_plane {
            m.elements.push(ElementData {
                id: 21,
                elem_type: 42,
                stiff_id: 88,
                nodes: vec![1, 3, 4],
            });
        }
        m
    }

    #[test]
    fn short_axis_translates_with_support_without_losing_vector_or_properties() {
        for transformed in [false, true] {
            let mut m = short_feature(false);
            let mut p = policy();
            if transformed {
                let rotation =
                    glam::DQuat::from_axis_angle(DVec3::new(1., 2., 3.).normalize(), 0.6);
                m.nodes = m
                    .nodes
                    .into_iter()
                    .map(|(id, point)| (id + 100, rotation * point + DVec3::new(20., 30., 40.)))
                    .collect();
                for e in &mut m.elements {
                    e.id += 200;
                    for n in &mut e.nodes {
                        *n += 100;
                    }
                }
                p.up = (rotation * DVec3::Z).to_array();
            }
            let r = run(&m, &p);
            assert!(r.accepted, "{}: {}", r.reason, r.candidate_max_residual);
            assert_eq!(r.short_axis_indices, vec![0]);
            assert_eq!(r.axes.len(), 1);
            let [i, j] = r.axes[0].endpoints;
            let old =
                DVec3::from_array(r.reference_points[j]) - DVec3::from_array(r.reference_points[i]);
            let new = DVec3::from_array(r.points[j]) - DVec3::from_array(r.points[i]);
            assert!(old.distance(new) < 2. * p.residual_tolerance);
            assert!(r.maximum_movement > 1e-6);
            assert!(r.axis_failures.is_empty());
            assert_eq!(r.axes[0].spans.len(), 1);
            assert_eq!(r.axes[0].spans[0].stiffness, 77);
            assert_eq!(
                r.axes[0].spans[0].element,
                if transformed { 210 } else { 10 }
            );
            assert_eq!(
                r.surfaces[0]
                    .stiffness_regions
                    .keys()
                    .copied()
                    .collect::<Vec<_>>(),
                vec![88]
            );
        }
    }

    #[test]
    fn conflicting_short_axis_is_reported_without_collapsing_or_dropping_it() {
        let m = short_feature(true);
        let r = run(&m, &policy());
        assert!(!r.accepted);
        assert!(!r.candidate_constraints_satisfied);
        assert_eq!(r.points, r.reference_points);
        assert_eq!(r.axes.len(), 1);
        assert_eq!(r.axes[0].spans[0].element, 10);
        assert!(r
            .largest_constraint_failures
            .iter()
            .any(|f| matches!(f.origin, ConstraintOrigin::ShortAxisVector { .. })));
    }

    #[test]
    fn coplanar_patches_meeting_at_one_vertex_share_support() {
        let mut m = MeshData::default();
        for (i, p) in [
            [0., 0., 0.],
            [1., 0., 0.],
            [0., 1., 0.],
            [-1., 0., 0.],
            [0., -1., 0.],
        ]
        .iter()
        .enumerate()
        {
            m.nodes.insert(i as u32 + 1, DVec3::from_array(*p));
        }
        m.elements = vec![
            ElementData {
                id: 1,
                elem_type: 42,
                stiff_id: 1,
                nodes: vec![1, 2, 3],
            },
            ElementData {
                id: 2,
                elem_type: 42,
                stiff_id: 2,
                nodes: vec![1, 4, 5],
            },
        ];
        let p = planes::recognize(
            &m,
            &planes::Policy {
                angle: 0.02,
                distance: 0.01,
                precision: 1e-8,
            },
        )
        .unwrap();
        assert_eq!(p.patches.len(), 2);
        assert_eq!(support_families(&m, &p), vec![vec![0, 1]]);
        let a = recognize::recognize(
            &m,
            &recognize::Policy {
                angle: 0.02,
                line_tolerance: 0.01,
                numerical_precision: 1e-8,
            },
        )
        .unwrap();
        let result = solve(&m, &a, &p, &policy()).unwrap();
        assert!(result.accepted);
        assert_eq!(result.plane_families, vec![vec![0, 1]]);
        let graph = super::super::graph::Graph::from_frame(&result);
        assert_eq!(graph.components.len(), 1);
        assert_eq!(graph.components[0].planes, vec![0, 1]);
        assert_eq!(result.surfaces[0].stiffness_regions.len(), 1);
        assert_eq!(result.surfaces[1].stiffness_regions.len(), 1);
    }

    #[test]
    fn plane_and_axes_share_interior_anchor_after_regularization() {
        let r = run(&source(), &policy());
        let graph = super::super::graph::Graph::from_frame(&r);
        assert_eq!(graph.components.len(), 1);
        assert_eq!(graph.components[0].axes.len(), 3);
        assert_eq!(graph.components[0].planes.len(), 1);
        assert_eq!(graph.components[0].vertices.len(), r.node_ids.len());
        assert!(r.accepted, "{}", r.candidate_max_residual);
        assert_eq!(r.axes.len(), 3);
        for a in &r.axes {
            let x = DVec3::from_array(r.points[a.endpoints[0]]);
            let y = DVec3::from_array(r.points[a.endpoints[1]]);
            for c in &a.anchors {
                assert!(DVec3::from_array(r.points[c.node]).distance(x + c.t * (y - x)) < 2e-8);
            }
        }
        for s in &r.surfaces {
            for &n in &s.nodes {
                assert!(s.plane.distance(r.points[n]).abs() < 2e-8);
            }
        }
        let p = r.points[r.node_ids.iter().position(|n| *n == 2).unwrap()];
        let q = r.points[r.node_ids.iter().position(|n| *n == 4).unwrap()];
        assert!((p[0] - q[0]).abs() < 2e-8 && (p[1] - q[1]).abs() < 2e-8);
    }
    #[test]
    fn failure_returns_original_geometry_not_last_iterate() {
        let mut p = policy();
        p.maximum_movement = 1e-6;
        p.iterations = 100;
        let r = run(&source(), &p);
        assert!(!r.accepted);
        assert_eq!(r.points, r.reference_points);
        assert_eq!(r.maximum_movement, 0.);
        assert!(!r.movement_failures.is_empty());
        assert!(r.movement_failures.iter().all(|f| f.movement > f.budget));
    }
    #[test]
    fn rigid_transform_with_up_preserves_a_feasible_solution() {
        let mut m = source();
        let original_graph = super::super::graph::Graph::from_frame(&run(&m, &policy()));
        let rotation = glam::DQuat::from_axis_angle(DVec3::new(1., 2., 3.).normalize(), 0.6);
        for v in m.nodes.values_mut() {
            *v = rotation * (*v) + DVec3::new(20., 30., 40.);
        }
        let mut p = policy();
        p.up = (rotation * DVec3::Z).to_array();
        let r = run(&m, &p);
        assert!(r.accepted);
        assert!(r.maximum_movement < 0.15);
        assert_eq!(
            serde_json::to_value(original_graph).unwrap(),
            serde_json::to_value(super::super::graph::Graph::from_frame(&r)).unwrap()
        );
    }

    /// Quad strip between `a` and `b` extruded by `h` along `up`, `n`
    /// elements long; returns the next free node ID.
    fn strip(
        m: &mut MeshData,
        first: u32,
        stiff: u32,
        a: DVec3,
        b: DVec3,
        up: DVec3,
        n: u32,
    ) -> u32 {
        for i in 0..=n {
            let p = a + (b - a) * (i as f64 / n as f64);
            m.nodes.insert(first + 2 * i, p);
            m.nodes.insert(first + 2 * i + 1, p + up);
        }
        for i in 0..n {
            let k = first + 2 * i;
            m.elements.push(ElementData {
                id: first + i,
                elem_type: 44,
                stiff_id: stiff,
                nodes: vec![k, k + 2, k + 3, k + 1],
            });
        }
        first + 2 * n + 2
    }

    /// Wall B ends `gap` short of wall A (a corner left open by the source).
    fn open_corner(gap: f64, transform: impl Fn(DVec3) -> DVec3) -> MeshData {
        let mut m = MeshData::default();
        let up = DVec3::Z * 3.;
        let t = |p: DVec3| transform(p);
        let mut raw = MeshData::default();
        let next = strip(&mut raw, 1, 1, DVec3::ZERO, DVec3::new(0., 2., 0.), up, 4);
        strip(
            &mut raw,
            next,
            2,
            DVec3::new(gap, 0., 0.),
            DVec3::new(2., 0., 0.),
            up,
            4,
        );
        m.elements = raw.elements;
        m.nodes = raw.nodes.into_iter().map(|(k, p)| (k, t(p))).collect();
        m
    }

    fn planes_of(m: &MeshData, scale: f64) -> planes::Report {
        planes::recognize(
            m,
            &planes::Policy {
                angle: 0.02,
                distance: 0.01 * scale,
                precision: 1e-8 * scale,
            },
        )
        .unwrap()
    }

    #[test]
    fn open_corner_gap_becomes_virtual_incidence_under_transforms() {
        for scale in [0.5, 1., 4.] {
            for transformed in [false, true] {
                let rotation = if transformed {
                    glam::DQuat::from_axis_angle(DVec3::new(1., 2., 3.).normalize(), 0.4)
                } else {
                    glam::DQuat::IDENTITY
                };
                let shift = DVec3::new(10., -20., 5.) * transformed as u8 as f64;
                let m = open_corner(0.026, |p| rotation * (p * scale) + shift);
                let planes = planes_of(&m, scale);
                let axes = recognize::recognize(
                    &m,
                    &recognize::Policy {
                        angle: 0.02,
                        line_tolerance: 0.01 * scale,
                        numerical_precision: 1e-8 * scale,
                    },
                )
                .unwrap();
                let mut p = policy();
                p.up = (rotation * DVec3::Z).to_array();
                p.maximum_movement *= scale;
                p.minimum_length *= scale;
                let r = solve_closing_gaps(&m, &axes, &planes, &p, 3, 0.05 * scale, 0.001 * scale)
                    .unwrap();
                assert!(r.accepted);
                let v = &r.virtual_incidences;
                // Both end nodes of wall B (bottom and top of its free end).
                assert_eq!(v.applied.len(), 2, "{v:?}");
                assert!(v.dropped.is_empty());
                let wall_a = &planes.patches[v.applied[0].patch];
                assert!(wall_a.source_nodes.contains(&1));
                let normal = DVec3::from_array(wall_a.plane.normal);
                let on_a =
                    DVec3::from_array(r.points[r.node_ids.iter().position(|&n| n == 1).unwrap()]);
                for c in &v.applied {
                    assert!((c.height - 0.026 * scale).abs() < 1e-6 * scale);
                    let i = r.node_ids.iter().position(|&n| n == c.node_id).unwrap();
                    let q = DVec3::from_array(r.points[i]);
                    assert!((q - on_a).dot(normal).abs() < 1e-6 * scale);
                }
            }
        }
    }

    #[test]
    fn wide_gaps_parallel_structures_and_interior_nodes_are_not_candidates() {
        let m = open_corner(0.08, |p| p);
        let planes = planes_of(&m, 1.);
        let ids: Vec<u32> = m.nodes.keys().copied().collect();
        assert!(gaps::candidates(&m, &planes, &ids, 0.05, 0.001, 0.02).is_empty());

        // Two slabs 30 mm apart in level, overlapping in plan: never joined.
        let mut m = MeshData::default();
        let next = strip(
            &mut m,
            1,
            1,
            DVec3::ZERO,
            DVec3::new(2., 0., 0.),
            DVec3::Y * 2.,
            2,
        );
        strip(
            &mut m,
            next,
            2,
            DVec3::new(1.99, 0., 0.03),
            DVec3::new(4., 0., 0.03),
            DVec3::Y * 2.,
            2,
        );
        let planes = planes_of(&m, 1.);
        let ids: Vec<u32> = m.nodes.keys().copied().collect();
        assert!(gaps::candidates(&m, &planes, &ids, 0.05, 0.001, 0.02).is_empty());

        // A wall crossing just above a slab: interior (non-contour) wall
        // nodes are never pulled onto the slab, only its free edge nodes.
        let mut m = MeshData::default();
        let next = strip(
            &mut m,
            1,
            1,
            DVec3::new(-2., -2., 0.),
            DVec3::new(2., -2., 0.),
            DVec3::Y * 4.,
            4,
        );
        strip(
            &mut m,
            next,
            2,
            DVec3::new(-1., 0., 0.02),
            DVec3::new(1., 0., 0.02),
            DVec3::Z * 3.,
            4,
        );
        let planes = planes_of(&m, 1.);
        let ids: Vec<u32> = m.nodes.keys().copied().collect();
        let c = gaps::candidates(&m, &planes, &ids, 0.05, 0.001, 0.02);
        assert_eq!(c.len(), 5, "{c:?}");
        assert!(c.iter().all(|c| m.nodes[&c.node_id].z.abs() < 0.03));
    }

    #[test]
    fn narrow_panel_above_a_slab_extends_to_it_instead_of_collapsing() {
        // A 30 mm high panel 10 mm above a slab: only its bottom closes.
        let mut m = MeshData::default();
        let next = strip(
            &mut m,
            1,
            1,
            DVec3::new(-2., -2., 0.),
            DVec3::new(2., -2., 0.),
            DVec3::Y * 4.,
            2,
        );
        let first = next;
        strip(
            &mut m,
            next,
            2,
            DVec3::new(-1., 0., 0.01),
            DVec3::new(1., 0., 0.01),
            DVec3::Z * 0.03,
            1,
        );
        let planes = planes_of(&m, 1.);
        let axes = recognize::recognize(
            &m,
            &recognize::Policy {
                angle: 0.02,
                line_tolerance: 0.01,
                numerical_precision: 1e-8,
            },
        )
        .unwrap();
        let r = solve_closing_gaps(&m, &axes, &planes, &policy(), 3, 0.05, 0.001).unwrap();
        assert!(r.accepted);
        let v = &r.virtual_incidences;
        let mut applied: Vec<u32> = v.applied.iter().map(|c| c.node_id).collect();
        applied.sort();
        assert_eq!(applied, vec![first, first + 2], "{v:?}");
        let at =
            |n: u32| DVec3::from_array(r.points[r.node_ids.iter().position(|&x| x == n).unwrap()]);
        for k in [first, first + 2] {
            assert!((at(k).z - at(1).z).abs() < 1e-6);
            assert!((at(k + 1).z - at(k).z) > 0.03);
        }

        // The post-solve guard: an edge collapsed at an incidence node.
        let incidence = |n| gaps::Incidence {
            node_id: n,
            patch: 0,
            height: 0.04,
            distance: 0.04,
        };
        let ids: Vec<u32> = m.nodes.keys().copied().collect();
        let mut points: Vec<[f64; 3]> = ids.iter().map(|n| m.nodes[n].to_array()).collect();
        let top = ids.iter().position(|&n| n == first + 1).unwrap();
        points[top] = m.nodes[&first].to_array();
        let c = gaps::collapsed_edges(
            &m,
            &[incidence(first + 1), incidence(first + 3)],
            &ids,
            &points,
            0.001,
        );
        assert_eq!(c.into_iter().collect::<Vec<_>>(), vec![first + 1]);
    }

    /// Wall A in x = 0 (z 0..1) and wall B above it, rotated `tilt` rad about
    /// the vertical, joined by two vertical bars 2 m apart.
    fn tilted_walls(tilt: f64) -> MeshData {
        let mut m = MeshData::default();
        let next = strip(
            &mut m,
            1,
            1,
            DVec3::ZERO,
            DVec3::new(0., 2., 0.),
            DVec3::Z,
            2,
        );
        let b0 = DVec3::new(0., 0., 1.2);
        let b1 = b0 + DVec3::new(2. * tilt, 2., 0.);
        let bars = strip(&mut m, next, 2, b0, b1, DVec3::Z, 2);
        for (k, lower) in [2u32, 6].into_iter().enumerate() {
            m.elements.push(ElementData {
                id: bars + k as u32,
                elem_type: 10,
                stiff_id: 3,
                nodes: vec![lower, next + 4 * k as u32],
            });
        }
        m
    }

    #[test]
    fn nearly_parallel_walls_share_a_direction_only_when_needed() {
        for (tilt, regularized) in [(0., false), (1e-4, true)] {
            let m = tilted_walls(tilt);
            let axes = recognize::recognize(
                &m,
                &recognize::Policy {
                    angle: 0.02,
                    line_tolerance: 0.01,
                    numerical_precision: 1e-8,
                },
            )
            .unwrap();
            let planes = planes_of(&m, 1.);
            assert_eq!(planes.patches.len(), 2);
            let r = solve_with_retry(&m, &axes, &planes, &policy(), 3).unwrap();
            assert!(
                r.accepted,
                "tilt {tilt}: {} {}",
                r.reason, r.candidate_max_residual
            );
            assert_eq!(r.regularized_directions, regularized, "tilt {tilt}");
            // With the recognised normals the tilted pair cannot converge.
            let plain = retry_with(&m, &axes, &planes, &policy(), 3, &[], false).unwrap();
            assert_eq!(plain.accepted, !regularized, "tilt {tilt}");
        }
    }

    #[test]
    fn narrow_patch_with_an_ill_defined_normal_joins_its_coplanar_neighbour() {
        // A wall in x = 0 and a triangle of another stiffness filling a notch
        // 0.105 m wide: its apex is 3 mm off the wall plane, so its fitted
        // normal is tilted 1.7 degrees (beyond the recognition angle), yet all
        // its nodes lie within the plane distance tolerance of the wall.
        // A 1 m wide panel tilted the same way is a different plane.
        for (width, grouped) in [(0.105, true), (1., false)] {
            let mut m = MeshData::default();
            strip(
                &mut m,
                1,
                1,
                DVec3::new(0., -2., 0.),
                DVec3::ZERO,
                DVec3::Z,
                4,
            );
            // Nodes 9 (y=0, z=0) and 10 (y=0, z=1) are the wall's edge.
            let apex = DVec3::new(width * 0.0298, width, 1.);
            m.nodes.insert(100, apex);
            m.elements.push(ElementData {
                id: 100,
                elem_type: 42,
                stiff_id: 2,
                nodes: vec![9, 100, 10],
            });
            let planes = planes_of(&m, 1.);
            assert_eq!(planes.patches.len(), 2, "width {width}");
            let families = support_families(&m, &planes);
            assert_eq!(families.len() == 1, grouped, "width {width}: {families:?}");
        }
    }

    #[test]
    fn narrow_patch_bridging_two_facets_joins_one_family_only() {
        // Two 1 m wide wall facets meeting at a vertical edge at 3 degrees
        // (separate planes), and a narrow triangle of another stiffness
        // touching both near that edge: all three do not fit one plane, so
        // the triangle joins one facet's family instead of all becoming
        // separate families.
        let turn = 3f64.to_radians();
        let mut m = MeshData::default();
        // Facet A: nodes 1..=6, node 5 = (0, 0, 0), node 6 = (0, 0, 3).
        strip(
            &mut m,
            1,
            1,
            DVec3::new(-1., 0., 0.),
            DVec3::ZERO,
            DVec3::Z * 3.,
            2,
        );
        // Facet B shares the edge 5-6.
        let far = DVec3::new(turn.cos(), turn.sin(), 0.);
        for (id, p) in [(7, far * 0.5), (9, far)] {
            m.nodes.insert(id, p);
            m.nodes.insert(id + 1, p + DVec3::Z * 3.);
        }
        for (id, nodes) in [(3, vec![5, 7, 8, 6]), (4, vec![7, 9, 10, 8])] {
            m.elements.push(ElementData {
                id,
                elem_type: 44,
                stiff_id: 2,
                nodes,
            });
        }
        // A 3 mm wide triangle on the shared edge, its fitted normal ill
        // defined (tilted about 10 degrees from both facets).
        m.nodes.insert(11, DVec3::new(-0.003, 0.0005, 0.4));
        m.elements.push(ElementData {
            id: 5,
            elem_type: 42,
            stiff_id: 3,
            nodes: vec![5, 11, 6],
        });
        let planes = planes_of(&m, 1.);
        assert_eq!(planes.patches.len(), 3);
        let families = support_families(&m, &planes);
        assert_eq!(families.len(), 2, "{families:?}");
        let mut all: Vec<usize> = families.iter().flatten().copied().collect();
        all.sort_unstable();
        assert_eq!(all, vec![0, 1, 2]);
    }

    /// A straight lower wall (plane y = 0) and, above the slab line z = 0,
    /// an upper wall that continues it to x = 0 and then kinks by about two
    /// degrees; the kinked facet shares the slab line node at x = 0.5
    /// (8 mm off the lower wall) with the lower wall, 50 cm from the kink.
    fn kinked_upper_wall(transform: impl Fn(DVec3) -> DVec3) -> MeshData {
        let mut m = MeshData::default();
        let points = [
            (1, [-1., 0., -3.]),
            (2, [0., 0., -3.]),
            (3, [0.5, 0., -3.]),
            (4, [-1., 0., 0.]),
            (5, [0., 0., 0.]),
            (6, [0.5, 0.008, 0.]),
            (7, [-1., 0., 3.]),
            (8, [0., 0., 3.]),
            (9, [0.5, 0.021, 3.]),
        ];
        for (id, p) in points {
            m.nodes.insert(id, transform(DVec3::from_array(p)));
        }
        for (id, stiff, nodes) in [
            (1, 1, vec![1, 2, 5, 4]),
            (2, 1, vec![2, 3, 6, 5]),
            (3, 2, vec![4, 5, 8, 7]),
            (4, 2, vec![5, 6, 9, 8]),
        ] {
            m.elements.push(ElementData {
                id,
                elem_type: 44,
                stiff_id: stiff,
                nodes,
            });
        }
        m
    }

    #[test]
    fn kinked_wall_sharing_nodes_off_the_kink_becomes_one_panel() {
        let rotation = glam::DQuat::from_axis_angle(DVec3::Z, 0.7);
        let shift = DVec3::new(30., -12., 4.);
        for rotated in [false, true] {
            let transform = |p: DVec3| if rotated { rotation * p + shift } else { p };
            let mut m = kinked_upper_wall(transform);
            if rotated {
                m.elements.reverse();
            }
            let planes = planes_of(&m, 1.);
            assert!(planes.patches.len() >= 2, "{}", planes.patches.len());
            let mut p = policy();
            let without = run(&m, &p);
            assert!(!without.accepted);
            p.panel_tolerance = 0.05;
            let with = run(&m, &p);
            assert!(with.accepted, "{:?}", with.movement_failures);
            assert_eq!(with.plane_families.len(), 1, "{:?}", with.plane_families);
            assert!(with.maximum_movement <= 0.05);
            // A tolerance below the kink deviation keeps the families apart.
            p.panel_tolerance = 0.002;
            assert!(run(&m, &p).plane_families.len() > 1);
        }
    }

    #[test]
    fn walls_meeting_at_a_corner_line_stay_separate_panels() {
        // Two 1 m facets 3 degrees apart sharing only their vertical edge:
        // the common nodes lie on the intersection line, no merge.
        let turn = 3f64.to_radians();
        let mut m = MeshData::default();
        strip(
            &mut m,
            1,
            1,
            DVec3::new(-1., 0., 0.),
            DVec3::ZERO,
            DVec3::Z * 3.,
            2,
        );
        let far = DVec3::new(turn.cos(), turn.sin(), 0.);
        for (id, p) in [(7, far * 0.5), (9, far)] {
            m.nodes.insert(id, p);
            m.nodes.insert(id + 1, p + DVec3::Z * 3.);
        }
        for (id, nodes) in [(3, vec![5, 7, 8, 6]), (4, vec![7, 9, 10, 8])] {
            m.elements.push(ElementData {
                id,
                elem_type: 44,
                stiff_id: 2,
                nodes,
            });
        }
        let mut p = policy();
        p.panel_tolerance = 0.05;
        let r = run(&m, &p);
        assert!(r.accepted);
        assert_eq!(r.plane_families.len(), 2);
    }

    #[test]
    fn short_bar_translates_with_a_node_moved_onto_its_wall_plane() {
        // A wall node 8 mm off the wall plane carries a 20 mm bar: the node
        // moves about 5 mm onto the plane, the bar translates with it,
        // beyond 5 % of its length but within the plane distance.
        let mut m = MeshData::default();
        for (id, p) in [
            (1, [0., 0., 0.]),
            (2, [0., 1., 0.]),
            (3, [0., 2., 0.]),
            (4, [0., 0., 1.]),
            (5, [0.008, 1., 1.]),
            (6, [0., 2., 1.]),
            (7, [0.028, 1., 1.]),
        ] {
            m.nodes.insert(id, DVec3::from_array(p));
        }
        for (id, nodes) in [(1, vec![1, 2, 5, 4]), (2, vec![2, 3, 6, 5])] {
            m.elements.push(ElementData {
                id,
                elem_type: 44,
                stiff_id: 1,
                nodes,
            });
        }
        m.elements.push(ElementData {
            id: 3,
            elem_type: 10,
            stiff_id: 2,
            nodes: vec![5, 7],
        });
        let r = run(&m, &policy());
        assert!(r.accepted, "{:?}", r.movement_failures);
        assert_eq!(r.short_axis_indices.len(), 1);
        let k = |id: u32| r.node_ids.iter().position(|&n| n == id).unwrap();
        let moved = DVec3::from_array(r.points[k(5)]).distance(m.nodes[&5]);
        assert!(moved > 0.001 && moved <= 0.01, "{moved}");
        let vector = DVec3::from_array(r.points[k(7)]) - DVec3::from_array(r.points[k(5)]);
        assert!((vector - (m.nodes[&7] - m.nodes[&5])).length() < 1e-6);
        assert!((r.movement_budgets()[k(7)] - 0.01).abs() < 1e-12);
    }

    #[test]
    fn short_axis_lying_on_its_plane_follows_it_in_geotechnical_mode() {
        // The conflicting short feature (both ends on one slightly tilted
        // plane, the bar 0.01 mm off it) is accepted when short axes follow
        // their plane; its in-plane vector is kept.
        let m = short_feature(true);
        let mut p = policy();
        p.geotechnical = true;
        let r = run(&m, &p);
        assert!(r.accepted, "{}: {}", r.reason, r.candidate_max_residual);
        assert!(r.axis_failures.is_empty());
        let [i, j] = r.axes[0].endpoints;
        let old =
            DVec3::from_array(r.reference_points[j]) - DVec3::from_array(r.reference_points[i]);
        let new = DVec3::from_array(r.points[j]) - DVec3::from_array(r.points[i]);
        assert!(old.distance(new) < 1e-4, "{}", old.distance(new));
        assert!(new.length() > 0.9 * old.length());
    }

    #[test]
    fn bar_end_on_a_wall_moves_with_the_wall_in_geotechnical_mode() {
        // A 33 mm bar (above the short-axis length) on a wall node 8 mm off
        // the wall plane: its end moves about 5 mm, beyond 5 % of its length;
        // accepted only in geotechnical mode.
        let mut m = MeshData::default();
        for (id, p) in [
            (1, [0., 0., 0.]),
            (2, [0., 1., 0.]),
            (3, [0., 2., 0.]),
            (4, [0., 0., 1.]),
            (5, [0.008, 1., 1.]),
            (6, [0., 2., 1.]),
            (7, [0.041, 1., 1.]),
        ] {
            m.nodes.insert(id, DVec3::from_array(p));
        }
        for (id, nodes) in [(1, vec![1, 2, 5, 4]), (2, vec![2, 3, 6, 5])] {
            m.elements.push(ElementData {
                id,
                elem_type: 44,
                stiff_id: 1,
                nodes,
            });
        }
        m.elements.push(ElementData {
            id: 3,
            elem_type: 10,
            stiff_id: 2,
            nodes: vec![5, 7],
        });
        let mut p = policy();
        let strict = run(&m, &p);
        assert!(!strict.accepted);
        assert!(!strict.movement_failures.is_empty());
        p.geotechnical = true;
        let r = run(&m, &p);
        assert!(r.accepted, "{:?}", r.movement_failures);
        let k = |id: u32| r.node_ids.iter().position(|&n| n == id).unwrap();
        assert!(DVec3::from_array(r.points[k(5)]).distance(m.nodes[&5]) <= 0.01);
    }

    #[test]
    fn flattened_short_axis_is_valid_with_its_in_plane_length() {
        // A 6.3 mm bar between two nodes of one slab, 2 mm across the slab
        // plane: flattened, it is 0.3 mm shorter, which is not a collapse.
        let mut m = short_feature(true);
        m.nodes.insert(4, DVec3::new(0.006, 0., 0.002));
        let mut p = policy();
        p.geotechnical = true;
        let r = run(&m, &p);
        assert!(r.accepted, "{}: {}", r.reason, r.candidate_max_residual);
        let [i, j] = r.axes[0].endpoints;
        let new = DVec3::from_array(r.points[j]) - DVec3::from_array(r.points[i]);
        assert!(
            new.length() < 0.0062 && new.length() > 0.0058,
            "{}",
            new.length()
        );
    }

    /// Wall nodes A and C (C 0.2 mm off the wall plane) joined through a
    /// slab node B by two 25 mm horizontal bars along the wall.
    fn wall_slab_wall_links() -> MeshData {
        let mut m = MeshData::default();
        for (id, p) in [
            (1, [0., 0., 0.]),
            (2, [1., 0., 0.]),
            (3, [1., 0., 1.]),
            (4, [0.35, 0.0002, 1.]),
            (5, [0.3, 0., 1.]),
            (6, [0., 0., 1.]),
            (7, [0.325, 0.0001, 1.]),
            (8, [0.325, 1., 1.]),
            (9, [1., 1., 1.]),
        ] {
            m.nodes.insert(id, DVec3::from_array(p));
        }
        for (id, ty, stiff, nodes) in [
            (1, 42, 1, vec![1, 2, 4]),
            (2, 42, 1, vec![2, 3, 4]),
            (3, 42, 1, vec![1, 4, 5]),
            (4, 42, 1, vec![1, 5, 6]),
            (5, 42, 2, vec![7, 8, 9]),
            (6, 10, 3, vec![5, 7]),
            (7, 10, 3, vec![7, 4]),
        ] {
            m.elements.push(ElementData {
                id,
                elem_type: ty,
                stiff_id: stiff,
                nodes,
            });
        }
        m
    }

    #[test]
    fn short_links_keep_length_and_level_not_transverse_noise() {
        let m = wall_slab_wall_links();
        let mut p = policy();
        p.geotechnical = true;
        let r = run(&m, &p);
        assert!(
            r.accepted,
            "{}: {} {:?}",
            r.reason, r.candidate_max_residual, r.short_axis_indices
        );
        for axis in &r.axes {
            let [i, j] = axis.endpoints;
            let old =
                DVec3::from_array(r.reference_points[j]) - DVec3::from_array(r.reference_points[i]);
            let new = DVec3::from_array(r.points[j]) - DVec3::from_array(r.points[i]);
            assert!((old.length() - new.length()).abs() < 1e-3);
            assert!(new.z.abs() < 1e-6);
        }
    }

    /// A 3 x 3 m slab (1 m quads) whose middle panel is a pyramid: apex
    /// 18 mm up, base corners 4 mm up (noise within the plane distance).
    fn slab_with_a_bump() -> MeshData {
        let mut m = MeshData::default();
        let mut id = 0;
        let mut node = |m: &mut MeshData, p: [f64; 3]| {
            id += 1;
            m.nodes.insert(id, DVec3::from_array(p));
            id
        };
        let mut grid = [[0u32; 4]; 4];
        for (i, row) in grid.iter_mut().enumerate() {
            for (j, n) in row.iter_mut().enumerate() {
                let z = if (1..=2).contains(&i) && (1..=2).contains(&j) {
                    0.004
                } else {
                    0.
                };
                *n = node(&mut m, [i as f64, j as f64, z]);
            }
        }
        let apex = node(&mut m, [1.5, 1.5, 0.018]);
        let mut element = 0;
        let mut push = |m: &mut MeshData, nodes: Vec<u32>| {
            element += 1;
            m.elements.push(ElementData {
                id: element,
                elem_type: if nodes.len() == 4 { 44 } else { 42 },
                stiff_id: 1,
                nodes,
            });
        };
        for i in 0..3 {
            for j in 0..3 {
                let (a, b, c, d) = (
                    grid[i][j],
                    grid[i + 1][j],
                    grid[i + 1][j + 1],
                    grid[i][j + 1],
                );
                if (i, j) == (1, 1) {
                    for (x, y) in [(a, b), (b, c), (c, d), (d, a)] {
                        push(&mut m, vec![x, y, apex]);
                    }
                } else {
                    push(&mut m, vec![a, b, c, d]);
                }
            }
        }
        m
    }

    #[test]
    fn slab_bump_of_tilted_triangles_joins_the_slab_panel() {
        let m = slab_with_a_bump();
        let mut p = policy();
        p.panel_tolerance = 0.05;
        p.geotechnical = true;
        let r = run(&m, &p);
        assert_eq!(r.plane_families.len(), 1, "{:?}", r.plane_families);
        assert!(r.accepted, "{}: {}", r.reason, r.candidate_max_residual);
        assert!(r.maximum_movement <= 0.02, "{}", r.maximum_movement);
    }
}
