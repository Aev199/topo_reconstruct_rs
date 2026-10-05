//! Transfer of LIRA loads to the reconstructed geometry, for the PLAXIS
//! exchange.
//!
//! Loads are read by element and node (`parsers::loads`); the reconstruction
//! changed the geometry they act on, so each kind is mapped by what it
//! means:
//! * loads on shell elements (pressures, trapezoidal pressures, point
//!   forces) become a pressure vector per element (a point force is spread
//!   over its element); the elements of one surface with one value (values
//!   binned when there are too many) are the region of a surface load, its
//!   outline cut out of the element polygons, clipped to the surface;
//! * line loads along shell edges stay line loads on those edges;
//! * loads on bars follow the bar's axis (`source_axis`/spans of the
//!   reconstruction): distributed ones become line loads on the piece the
//!   element is, point forces and moments point loads at their place;
//! * node forces and moments become point loads at the node (or at the
//!   vertex it became).
//!
//! Everything LIRA-specific that cannot be carried (temperature, dynamic
//! and stage loads, prescribed displacements, plate moments) is counted and
//! reported, and every load case reports the resultant force of the source
//! loads next to that of the exported ones.
use crate::input::{ElementData, MeshData};
use crate::parsers::lira::Material;
use crate::parsers::loads::{LoadRow, LoadSet};
use crate::plaxis::{hole_free, plaxis_polygon, surface_polygons};
use crate::reconstruction::assembly::edit::State;
use geo::{Area, BooleanOps, Coord, LineString, MultiPolygon, Polygon};
use glam::{DVec2, DVec3};
use hashbrown::HashMap;
use serde::Serialize;
use std::collections::BTreeMap;

#[derive(Debug, Clone, Serialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum Load {
    /// Force (kN) and moment (kN m) at a point.
    Point {
        case: u32,
        at: [f64; 3],
        force: [f64; 3],
        moment: [f64; 3],
    },
    /// Line load (kN/m), linear from `q_start` to `q_end`.
    Line {
        case: u32,
        start: [f64; 3],
        end: [f64; 3],
        q_start: [f64; 3],
        q_end: [f64; 3],
    },
    /// Surface load (kN/m2, global components) on polygons without holes.
    Surface {
        case: u32,
        surface: usize,
        polygons: Vec<Vec<[f64; 3]>>,
        sigma: [f64; 3],
    },
}

impl Load {
    pub fn case(&self) -> u32 {
        match self {
            Load::Point { case, .. } | Load::Line { case, .. } | Load::Surface { case, .. } => *case,
        }
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct CaseReport {
    pub case: u32,
    pub name: String,
    pub loads: usize,
    /// Resultant force (kN) of the supported source loads and of the
    /// exported ones.
    pub source: [f64; 3],
    pub exported: [f64; 3],
    /// Resultant (kN) of plate loads on elements that are not in the
    /// reconstructed geometry (absorbed, removed or deleted), not counted in
    /// `source`.
    pub not_in_geometry: [f64; 3],
}

#[derive(Debug, Clone, Default, Serialize)]
pub struct Report {
    pub cases: Vec<CaseReport>,
    /// What was not transferred, with the number of load rows.
    pub skipped: BTreeMap<String, usize>,
    /// What was transferred approximately.
    pub approximated: BTreeMap<String, usize>,
}

/// Whether a load case is the self-weight of the structure (PLAXIS applies
/// it itself): "СВ", "СВ_...", "СОБСТВЕННЫЙ ВЕС ...", "self weight".
pub fn is_self_weight(name: &str) -> bool {
    let n = name.trim().to_lowercase();
    n == "св"
        || n.starts_with("св ")
        || n.starts_with("св_")
        || n.starts_with("св.")
        || n.contains("собствен")
        || n.contains("self weight")
        || n.contains("self-weight")
        || n.contains("selfweight")
}

/// Whether a load case looks dynamic (seismic, dynamic): not part of a static
/// settlement combination by default.
pub fn is_dynamic(name: &str) -> bool {
    let n = name.to_lowercase();
    n.contains("сейсм") || n.contains("динам") || n.contains("seism")
}

/// Whether a load case is a construction stage ("СТАДИЯ 3", "СТАДИЯ №3",
/// "Stage 3"): stages are modelled as load cases in LIRA but are not part of
/// the PLAXIS or MIDAS load sets.
pub fn is_stage(name: &str) -> bool {
    let n = name.trim().to_lowercase();
    n.starts_with("стади") || n.starts_with("stage")
}

/// The pseudo-case of the weight of the storeys that were cut off.
pub const CUT_WEIGHT_CASE: u32 = 1_000_000;

/// The case number of the load combination (`Settings::combination`).
pub const COMBINATION: u32 = 0;

/// How the combination is reduced for PLAXIS: its geometry must stay valid,
/// so a plate load is uniform over the whole plate and a bar load uniform
/// over the whole bar, or a point load where that would move the resultant.
#[derive(Debug, Clone, Copy)]
pub struct Simplify {
    /// Largest distance of the centre of pressure from the centre of the
    /// plate or bar for a uniform load, as a fraction of its size.
    pub center_tolerance: f64,
    /// Smallest loaded part of a plate or bar for a uniform load.
    pub min_fraction: f64,
}

impl Default for Simplify {
    fn default() -> Self {
        Simplify { center_tolerance: 0.15, min_fraction: 0.3 }
    }
}

/// Load cases combined into one case with a factor each.
#[derive(Debug, Clone)]
pub struct Combination {
    /// Factor per load case; cases that are absent or 0 are left out (the
    /// self-weight case, which PLAXIS applies itself).
    pub factors: BTreeMap<u32, f64>,
    /// `None`: the loads keep the contours of the source.
    pub simplify: Option<Simplify>,
}

#[derive(Debug, Clone)]
pub struct Settings {
    /// Force unit of the model in kN (LIRA tf: 9.80665).
    pub force_factor: f64,
    /// Snap distance of region outlines to the surface contour.
    pub snap: f64,
    /// Most surface loads of one surface and case: values are binned.
    pub max_groups: usize,
    /// One combined case instead of the cases of the source.
    pub combination: Option<Combination>,
    /// Only these cases of the source are read (stages, the self-weight and
    /// dynamic cases are left out); `None`: all of them.
    pub cases: Option<std::collections::BTreeSet<u32>>,
    /// Materials of the stiffness types: the weight of the storeys cut off
    /// (`CUT_WEIGHT_CASE`) and the shares of their supports need them.
    pub materials: Option<std::sync::Arc<hashbrown::HashMap<u32, Material>>>,
    /// Loads of the storeys cut off are carried to the level of the cut
    /// (`false`: they are dropped).
    pub cut_loads: bool,
}

struct Context<'a> {
    state: &'a State,
    mesh: &'a MeshData,
    set: &'a LoadSet,
    settings: Settings,
    /// Source node -> vertex of the reconstruction.
    vertex_of_node: HashMap<u32, usize>,
    surface_of_element: HashMap<u32, usize>,
    /// Source bar element -> (axis, t of its first node, t of its second).
    bar_of_element: HashMap<u32, (usize, f64, f64)>,
}

fn element(mesh: &MeshData, id: u32) -> Option<&ElementData> {
    mesh.elements
        .get((id as usize).checked_sub(1)?)
        .filter(|e| e.id == id)
}

fn positions(mesh: &MeshData, e: &ElementData) -> Option<Vec<DVec3>> {
    e.nodes.iter().map(|n| mesh.nodes.get(n).copied()).collect()
}

/// Area of a planar polygon (any simple outline, convex or not).
fn polygon_area(p: &[DVec3]) -> f64 {
    (1..p.len().saturating_sub(1))
        .map(|i| (p[i] - p[0]).cross(p[i + 1] - p[0]))
        .sum::<DVec3>()
        .length()
        / 2.
}

/// The perimeter order of a shell element's nodes, as indices. A quad
/// listed in tensor-product order (1 2 / 3 4) or with crossed sides has a
/// bow-tie outline; the order enclosing the largest area is the perimeter
/// (the listed one wins ties, so nodes on a straight side are harmless).
fn shell_order(p: &[DVec3]) -> Vec<usize> {
    if p.len() != 4 {
        return (0..p.len()).collect();
    }
    let area = |o: [usize; 4]| -> f64 {
        (0..4)
            .map(|i| p[o[i]].cross(p[o[(i + 1) % 4]]))
            .sum::<DVec3>()
            .length()
    };
    let mut best = ([0, 1, 2, 3], area([0, 1, 2, 3]));
    for o in [[0, 1, 3, 2], [0, 2, 1, 3]] {
        let a = area(o);
        if a > best.1 * (1. + 1e-9) + 1e-12 {
            best = (o, a);
        }
    }
    best.0.to_vec()
}

fn shell_ring(p: &[DVec3]) -> Vec<DVec3> {
    shell_order(p).into_iter().map(|i| p[i]).collect()
}

/// `shell_ring` with the node numbers.
fn shell_ring_nodes(p: &[DVec3], nodes: &[u32]) -> Vec<(u32, DVec3)> {
    shell_order(p).into_iter().map(|i| (nodes[i], p[i])).collect()
}

/// Local axes of a shell element: X1 from its first to its second node, Z1
/// the normal by the right-hand rule of the perimeter order.
fn shell_axes(ring: &[DVec3]) -> Option<[DVec3; 3]> {
    let newell: DVec3 = (0..ring.len())
        .map(|i| ring[i].cross(ring[(i + 1) % ring.len()]))
        .sum();
    let z = newell.try_normalize()?;
    let x = (ring[1] - ring[0]).try_normalize()?;
    let y = z.cross(x).try_normalize()?;
    Some([y.cross(z), y, z])
}

/// Local axes of a bar element: X1 first to second node; Y1 from document
/// 17 when given, else Z1 the upward normal in the vertical plane of the
/// bar (global X for a vertical bar), Y1 = Z1 x X1.
pub fn bar_axes(a: DVec3, b: DVec3, y_vector: Option<DVec3>) -> Option<[DVec3; 3]> {
    let x = (b - a).try_normalize()?;
    if let Some(y) = y_vector.and_then(|y| (y - x * y.dot(x)).try_normalize()) {
        return Some([x, y, x.cross(y)]);
    }
    let z = if DVec3::new(x.x, x.y, 0.).length() < 1e-9 {
        DVec3::X
    } else {
        (DVec3::Z - x * x.z).normalize()
    };
    Some([x, z.cross(x), z])
}

fn unit(direction: u8) -> Option<DVec3> {
    match direction {
        1 | 4 => Some(DVec3::X),
        2 | 5 => Some(DVec3::Y),
        3 | 6 => Some(DVec3::Z),
        _ => None,
    }
}

/// Frame of a load code: 0 local, 1 global, 2 projection, 4 levelling
/// (taken as global).
fn frame(code: u16) -> u16 {
    code / 10
}

/// Sign of force loads. A positive value in a LIRA file acts against the
/// positive axis it is given along: gravity loads are entered positive and
/// the wind cases named after a direction (X+, Y-) push towards it only with
/// the sign reversed (checked on the fixtures: the reversed resultants point
/// as the case names say, in global and local axes alike). Moments keep
/// their sign, as the existing LIRA to MIDAS converter does.
const FORCE_SIGN: f64 = -1.0;

fn skip(report: &mut Report, what: impl Into<String>) {
    *report.skipped.entry(what.into()).or_default() += 1;
}

pub fn transfer(
    state: &State,
    vertex_nodes: &[u32],
    mesh: &MeshData,
    set: &LoadSet,
    settings: Settings,
) -> (Vec<Load>, Report) {
    let model = &state.model;
    let mut surface_of_element = HashMap::new();
    for (s, surface) in model.surfaces().iter().enumerate() {
        for &e in &surface.source_elements {
            surface_of_element.insert(e, s);
        }
    }
    let mut bar_of_element = HashMap::new();
    for (i, axis) in state.axes.iter().enumerate() {
        let [a, b] = axis.endpoints.map(|v| DVec3::from_array(model.vertices()[v]));
        let d = b - a;
        let t_of = |p: DVec3| (p - a).dot(d) / d.length_squared();
        for span in &axis.spans {
            let Some(e) = element(mesh, span.element) else {
                continue;
            };
            let Some(p) = positions(mesh, e).filter(|p| p.len() == 2) else {
                continue;
            };
            // The element's first node is at the end of the span the node
            // lies nearer to.
            let (t0, t1) = (span.start_t, span.end_t);
            let first_at_start = (t_of(p[0]) - t0).abs() + (t_of(p[1]) - t1).abs()
                <= (t_of(p[0]) - t1).abs() + (t_of(p[1]) - t0).abs();
            bar_of_element.insert(
                span.element,
                if first_at_start { (i, t0, t1) } else { (i, t1, t0) },
            );
        }
    }
    let vertex_of_node = vertex_nodes
        .iter()
        .enumerate()
        .filter(|(v, _)| *v < model.vertices().len())
        .map(|(v, &n)| (n, v))
        .collect();
    let context = Context {
        state,
        mesh,
        set,
        settings,
        vertex_of_node,
        surface_of_element,
        bar_of_element,
    };
    context.run()
}

/// A point force on a plate element (a "stamp" row of LIRA: the stamp's
/// force given per element it covers).
struct Stamp {
    case: u32,
    element: u32,
    value: f64,
    /// Unit direction of the force, sign included.
    direction: DVec3,
    area: f64,
    /// Height of the element's centre along the force.
    level: f64,
}

/// Pressure of a shell element (tf/m2, global).
#[derive(Default, Clone, Copy)]
struct Pressure {
    vector: DVec3,
}

impl Context<'_> {
    fn node_position(&self, node: u32) -> Option<DVec3> {
        match self.vertex_of_node.get(&node) {
            Some(&v) => Some(DVec3::from_array(self.state.model.vertices()[v])),
            None => self.mesh.nodes.get(&node).copied(),
        }
    }

    fn run(&self) -> (Vec<Load>, Report) {
        let factor = self.settings.force_factor;
        let mut report = Report::default();
        let mut loads: Vec<Load> = vec![];
        let mut source: BTreeMap<u32, DVec3> = BTreeMap::new();
        // Pressures per case and shell element (tf/m2).
        let mut pressure: BTreeMap<u32, HashMap<u32, Pressure>> = BTreeMap::new();
        let mut node_loads: BTreeMap<(u32, u32), (DVec3, DVec3)> = BTreeMap::new();
        let mut bar_uniform: BTreeMap<(u32, u32), DVec3> = BTreeMap::new();
        // Linear segments on bar axes: (case, axis) -> (ta, tb, qa, qb).
        let mut bar_lines: BTreeMap<(u32, usize), Vec<(f64, f64, DVec3, DVec3)>> = BTreeMap::new();
        let mut edge_lines: BTreeMap<u32, Vec<(DVec3, DVec3, DVec3, DVec3)>> = BTreeMap::new();
        // Resultants of loads on plates that are not in the geometry.
        let mut lost: BTreeMap<u32, DVec3> = BTreeMap::new();
        let mut stamps: Vec<Stamp> = vec![];
        // Loads on what was cut off: (position, force, moment) in tf, tf m.
        let mut removed: BTreeMap<u32, Vec<(DVec3, DVec3, DVec3)>> = BTreeMap::new();
        let cut = self.state.cut.as_ref();
        let above = |p: DVec3| cut.is_some_and(|c| p.z > c.z + 1e-6);
        let mut shell_geometry: HashMap<u32, Option<(Vec<DVec3>, [DVec3; 3], f64)>> = HashMap::new();

        for row in &self.set.rows {
            if self.settings.cases.as_ref().is_some_and(|c| !c.contains(&row.case)) {
                continue;
            }
            let params = self.set.parameters(row.parameters);
            let code = row.code;
            // Node loads: the target is a node.
            if code == 0 {
                self.node_row(row, params, &mut node_loads, &mut report);
                continue;
            }
            let Some(e) = element(self.mesh, row.target) else {
                skip(&mut report, format!("код {code}: элемент {} не найден", row.target));
                continue;
            };
            match e.elem_type {
                41 | 42 | 44 => {
                    let geometry = shell_geometry.entry(e.id).or_insert_with(|| {
                        let ring = shell_ring(&positions(self.mesh, e)?);
                        let axes = shell_axes(&ring)?;
                        let area = polygon_area(&ring);
                        Some((ring, axes, area))
                    });
                    let Some((p, axes, area)) = geometry.clone() else {
                        skip(&mut report, "пластина с вырожденной геометрией");
                        continue;
                    };
                    self.shell_row(row, params, &p, &axes, area, &mut pressure, &mut edge_lines, &mut stamps, &mut report);
                }
                10 => self.bar_row(row, params, e, &mut bar_uniform, &mut bar_lines, &mut loads, &mut source, &mut removed, &mut report),
                other => skip(&mut report, format!("нагрузка на элемент типа {other}")),
            }
        }

        self.spread_stamps(&stamps, &mut pressure, &mut report);
        if cut.is_some() {
            // What was cut off: plate pressures, edge lines and node loads above the level.
            for (&case, elements) in pressure.iter_mut() {
                let gone: Vec<u32> = elements
                    .keys()
                    .copied()
                    .filter(|&e| element(self.mesh, e).and_then(|el| positions(self.mesh, el)).is_some_and(|p| above(p.iter().sum::<DVec3>() / p.len() as f64)))
                    .collect();
                for e in gone {
                    let vector = elements.remove(&e).map(|p| p.vector).unwrap_or_default();
                    if let Some(p) = element(self.mesh, e).and_then(|el| positions(self.mesh, el)) {
                        let area = polygon_area(&shell_ring(&p));
                        removed.entry(case).or_default().push((p.iter().sum::<DVec3>() / p.len() as f64, vector * area, DVec3::ZERO));
                    }
                }
            }
            for (&case, list) in edge_lines.iter_mut() {
                list.retain(|&(a, b, qa, qb)| {
                    if above((a + b) / 2.) {
                        removed.entry(case).or_default().push(((a + b) / 2., (qa + qb) / 2. * (b - a).length(), DVec3::ZERO));
                        false
                    } else {
                        true
                    }
                });
            }
            let gone: Vec<(u32, u32)> = node_loads.keys().copied().filter(|&(_, n)| self.mesh.nodes.get(&n).is_some_and(|&p| above(p))).collect();
            for key in gone {
                if let (Some((f, m)), Some(&p)) = (node_loads.remove(&key), self.mesh.nodes.get(&key.1)) {
                    removed.entry(key.0).or_default().push((p, f, m));
                }
            }
            if self.settings.cut_loads && self.settings.cases.as_ref().is_none_or(|c| c.contains(&CUT_WEIGHT_CASE)) {
                if let Some(items) = self.removed_weight() {
                    removed.insert(CUT_WEIGHT_CASE, items);
                }
            }
            self.cut_loads(&removed, &mut loads, &mut source, &mut report);
        }
        if let Some(combination) = &self.settings.combination {
            // One case of the weighted sum; cases without a factor are dropped.
            let weight = |case: u32| combination.factors.get(&case).copied().unwrap_or(0.);
            let mut combined: HashMap<u32, Pressure> = HashMap::new();
            for (&case, elements) in &pressure {
                for (&e, p) in elements {
                    combined.entry(e).or_default().vector += p.vector * weight(case);
                }
            }
            pressure = BTreeMap::from([(COMBINATION, combined)]);
            let mut nodes = BTreeMap::new();
            for (&(case, node), &(force, moment)) in &node_loads {
                let entry: &mut (DVec3, DVec3) = nodes.entry((COMBINATION, node)).or_default();
                entry.0 += force * weight(case);
                entry.1 += moment * weight(case);
            }
            node_loads = nodes;
            let mut uniform = BTreeMap::new();
            for (&(case, e), &q) in &bar_uniform {
                if weight(case) != 0. {
                    *uniform.entry((COMBINATION, e)).or_default() += q * weight(case);
                }
            }
            bar_uniform = uniform;
            let mut segments: BTreeMap<(u32, usize), Vec<(f64, f64, DVec3, DVec3)>> = BTreeMap::new();
            for (&(case, axis), list) in &bar_lines {
                let w = weight(case);
                if w == 0. {
                    continue;
                }
                segments
                    .entry((COMBINATION, axis))
                    .or_default()
                    .extend(list.iter().map(|&(ta, tb, qa, qb)| (ta, tb, qa * w, qb * w)));
            }
            bar_lines = segments;
            let mut edges: Vec<(DVec3, DVec3, DVec3, DVec3)> = vec![];
            for (&case, list) in &edge_lines {
                let w = weight(case);
                if w == 0. {
                    continue;
                }
                edges.extend(list.iter().map(|&(a, b, qa, qb)| (a, b, qa * w, qb * w)));
            }
            edge_lines = BTreeMap::from([(COMBINATION, edges)]);
            // Point loads and trapezoids of bars were made while reading.
            loads = std::mem::take(&mut loads)
                .into_iter()
                .filter(|l| weight(l.case()) != 0.)
                .map(|l| match l {
                    Load::Point { case, at, force, moment } => {
                        let w = weight(case);
                        Load::Point {
                            case: COMBINATION,
                            at,
                            force: (DVec3::from_array(force) * w).to_array(),
                            moment: (DVec3::from_array(moment) * w).to_array(),
                        }
                    }
                    Load::Line { case, start, end, q_start, q_end } => {
                        let w = weight(case);
                        Load::Line {
                            case: COMBINATION,
                            start,
                            end,
                            q_start: (DVec3::from_array(q_start) * w).to_array(),
                            q_end: (DVec3::from_array(q_end) * w).to_array(),
                        }
                    }
                    other => other,
                })
                .collect();
            let mut merged_source: BTreeMap<u32, DVec3> = BTreeMap::new();
            for (&case, &v) in &source {
                *merged_source.entry(COMBINATION).or_default() += v * weight(case);
            }
            source = merged_source;
            let mut merged_lost: BTreeMap<u32, DVec3> = BTreeMap::new();
            for (&case, &v) in &lost {
                *merged_lost.entry(COMBINATION).or_default() += v * weight(case);
            }
            lost = merged_lost;
        }

        // ---- node loads
        for (&(case, node), &(force, moment)) in &node_loads {
            let Some(at) = self.node_position(node) else {
                skip(&mut report, "узел вне модели");
                continue;
            };
            loads.push(Load::Point {
                case,
                at: at.to_array(),
                force: (force * factor).to_array(),
                moment: (moment * factor).to_array(),
            });
            *source.entry(case).or_default() += force * factor;
        }

        // ---- bars: uniform loads of whole elements, merged along their axis
        for (&(case, e), &q) in &bar_uniform {
            let Some(&(axis, t0, t1)) = self.bar_of_element.get(&e) else {
                skip(&mut report, "нагрузка на стержень, которого нет в геометрии");
                continue;
            };
            let (ta, tb) = (t0.min(t1), t0.max(t1));
            bar_lines.entry((case, axis)).or_default().push((ta, tb, q, q));
            if let Some(el) = element(self.mesh, e).and_then(|el| positions(self.mesh, el)) {
                *source.entry(case).or_default() += q * (el[1] - el[0]).length() * factor;
            }
        }
        for (&(case, axis), segments) in &mut bar_lines {
            let [a, b] = self.state.axes[axis].endpoints.map(|v| DVec3::from_array(self.state.model.vertices()[v]));
            segments.sort_by(|x, y| x.0.total_cmp(&y.0).then(x.1.total_cmp(&y.1)));
            let mut merged: Vec<(f64, f64, DVec3, DVec3)> = vec![];
            for &(ta, tb, qa, qb) in segments.iter() {
                if let Some(last) = merged.last_mut() {
                    let uniform = |s: &(f64, f64, DVec3, DVec3)| (s.2 - s.3).length() <= 1e-12;
                    if uniform(last)
                        && (qa - qb).length() <= 1e-12
                        && (last.2 - qa).length() <= 1e-9 * qa.length().max(1e-12)
                        && (ta - last.1).abs() <= 1e-6
                    {
                        last.1 = tb;
                        continue;
                    }
                }
                merged.push((ta, tb, qa, qb));
            }
            if let Some(simplify) = self.simplify() {
                self.simplified_bar(case, a, b, &merged, simplify, &mut loads, &mut report);
                continue;
            }
            for (ta, tb, qa, qb) in merged {
                loads.push(Load::Line {
                    case,
                    start: a.lerp(b, ta).to_array(),
                    end: a.lerp(b, tb).to_array(),
                    q_start: (qa * factor).to_array(),
                    q_end: (qb * factor).to_array(),
                });
            }
        }
        // Source resultants of the non-uniform bar segments were added by
        // `bar_row`; the exported side is computed below for all loads.

        // ---- plate edges
        for (&case, edges) in &edge_lines {
            let mut merged: Vec<(DVec3, DVec3, DVec3, DVec3)> = vec![];
            for &(a, b, qa, qb) in edges {
                *source.entry(case).or_default() += (qa + qb) / 2. * (b - a).length() * factor;
                // Consecutive collinear edges of one value make one line.
                if let Some(last) = merged.iter_mut().find(|l| {
                    (l.1 - a).length() <= 1e-6
                        && (l.2 - l.3).length() <= 1e-12
                        && (qa - qb).length() <= 1e-12
                        && (l.2 - qa).length() <= 1e-9 * qa.length().max(1e-12)
                        && (l.1 - l.0).normalize().dot((b - a).normalize()) > 1. - 1e-9
                }) {
                    last.1 = b;
                    continue;
                }
                merged.push((a, b, qa, qb));
            }
            if self.simplify().is_some() {
                self.simplified_edges(case, &merged, &mut loads, &mut report);
                continue;
            }
            for (a, b, qa, qb) in merged {
                loads.push(Load::Line {
                    case,
                    start: a.to_array(),
                    end: b.to_array(),
                    q_start: (qa * factor).to_array(),
                    q_end: (qb * factor).to_array(),
                });
            }
        }

        // ---- shell pressures: surface loads
        for (&case, elements) in &pressure {
            match self.simplify() {
                Some(simplify) => {
                    self.simplified_surface_loads(case, elements, simplify, &mut loads, &mut source, &mut lost, &mut report)
                }
                None => self.surface_loads(case, elements, &mut loads, &mut source, &mut lost, &mut report),
            }
        }

        // ---- reports per case
        let mut exported: BTreeMap<u32, DVec3> = BTreeMap::new();
        let mut counts: BTreeMap<u32, usize> = BTreeMap::new();
        for load in &loads {
            let case = load.case();
            *counts.entry(case).or_default() += 1;
            let sum = match load {
                Load::Point { force, .. } => DVec3::from_array(*force),
                Load::Line { start, end, q_start, q_end, .. } => {
                    let length = DVec3::from_array(*start).distance(DVec3::from_array(*end));
                    (DVec3::from_array(*q_start) + DVec3::from_array(*q_end)) / 2. * length
                }
                Load::Surface { polygons, sigma, .. } => {
                    let area: f64 = polygons
                        .iter()
                        .map(|p| polygon_area(&p.iter().map(|&x| DVec3::from_array(x)).collect::<Vec<_>>()))
                        .sum();
                    DVec3::from_array(*sigma) * area
                }
            };
            *exported.entry(case).or_default() += sum;
        }
        let cut_name = "Вес отброшенных этажей".to_string();
        let mut names: BTreeMap<u32, &String> = self.set.cases.iter().map(|(n, s)| (*n, s)).collect();
        names.insert(CUT_WEIGHT_CASE, &cut_name);
        let cases: std::collections::BTreeSet<u32> = counts
            .keys()
            .chain(source.keys())
            .chain(lost.keys())
            .copied()
            .collect();
        for case in cases {
            report.cases.push(CaseReport {
                case,
                name: names.get(&case).map(|s| s.to_string()).unwrap_or_default(),
                loads: counts.get(&case).copied().unwrap_or(0),
                source: source.get(&case).copied().unwrap_or_default().to_array(),
                exported: exported.get(&case).copied().unwrap_or_default().to_array(),
                not_in_geometry: lost.get(&case).copied().unwrap_or_default().to_array(),
            });
        }
        (loads, report)
    }

    fn simplify(&self) -> Option<Simplify> {
        self.settings.combination.as_ref().and_then(|c| c.simplify)
    }

    /// The loads of one bar axis (segments `(ta, tb, qa, qb)`, kN/m after
    /// the force factor) as a uniform load over the whole bar when their
    /// resultant acts near its middle and covers enough of it, else as point
    /// loads at the centres of the segments: no new points on the bar and
    /// the resultant is kept.
    #[allow(clippy::too_many_arguments)]
    fn simplified_bar(
        &self,
        case: u32,
        a: DVec3,
        b: DVec3,
        segments: &[(f64, f64, DVec3, DVec3)],
        simplify: Simplify,
        loads: &mut Vec<Load>,
        report: &mut Report,
    ) {
        let factor = self.settings.force_factor;
        let length = a.distance(b);
        let mut force = DVec3::ZERO;
        let mut moment_arm = 0.;
        let mut weight = 0.;
        let mut covered = 0.;
        let mut pieces = vec![];
        for &(ta, tb, qa, qb) in segments {
            let f = (qa + qb) / 2. * (tb - ta) * length;
            // The centre of a trapezoid by the magnitudes of its ends.
            let (ma, mb) = (qa.length(), qb.length());
            let t = if ma + mb > 1e-15 { ta + (tb - ta) * (ma + 2. * mb) / (3. * (ma + mb)) } else { (ta + tb) / 2. };
            force += f;
            moment_arm += f.length() * t;
            weight += f.length();
            covered += tb - ta;
            pieces.push((f, t));
        }
        if weight < 1e-15 || length < 1e-9 {
            return;
        }
        let center = moment_arm / weight;
        if covered.min(1.) >= simplify.min_fraction && (center - 0.5).abs() <= simplify.center_tolerance {
            loads.push(Load::Line {
                case,
                start: a.to_array(),
                end: b.to_array(),
                q_start: (force / length * factor).to_array(),
                q_end: (force / length * factor).to_array(),
            });
            if segments.len() > 1 || covered < 1. - 1e-6 {
                *report.approximated.entry("нагрузка на стержень приведена к равномерной по всей длине".into()).or_default() += 1;
            }
        } else {
            for (f, t) in pieces {
                loads.push(Load::Point {
                    case,
                    at: a.lerp(b, t).to_array(),
                    force: (f * factor).to_array(),
                    moment: [0.; 3],
                });
            }
            *report.approximated.entry("нагрузка на часть стержня заменена точечными силами".into()).or_default() += 1;
        }
    }

    /// Line loads along plate edges: kept as lines between model vertices
    /// when both ends lie on one (no new points in the geometry), else the
    /// resultant is applied at the centre of the load.
    fn simplified_edges(
        &self,
        case: u32,
        lines: &[(DVec3, DVec3, DVec3, DVec3)],
        loads: &mut Vec<Load>,
        report: &mut Report,
    ) {
        let factor = self.settings.force_factor;
        let vertices = self.state.model.vertices();
        let snap = self.settings.snap.max(1e-6);
        let mut grid: HashMap<[i64; 3], Vec<usize>> = HashMap::new();
        let cell = |p: DVec3| [0, 1, 2].map(|k| (p[k] / snap).floor() as i64);
        for (i, v) in vertices.iter().enumerate() {
            grid.entry(cell(DVec3::from_array(*v))).or_default().push(i);
        }
        let nearest = |p: DVec3| -> Option<DVec3> {
            let c = cell(p);
            let mut best: Option<(f64, DVec3)> = None;
            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        for &i in grid.get(&[c[0] + dx, c[1] + dy, c[2] + dz]).into_iter().flatten() {
                            let v = DVec3::from_array(vertices[i]);
                            let d = v.distance(p);
                            if d <= snap && best.is_none_or(|x| d < x.0) {
                                best = Some((d, v));
                            }
                        }
                    }
                }
            }
            best.map(|x| x.1)
        };
        // Equal loads on the same line (after snapping) add up.
        let mut by_line: BTreeMap<[i64; 6], (DVec3, DVec3, DVec3, DVec3)> = BTreeMap::new();
        for &(a, b, qa, qb) in lines {
            if let (Some(sa), Some(sb)) = (nearest(a), nearest(b)) {
                if sa.distance(sb) > 1e-9 {
                    let key = |p: DVec3| [0, 1, 2].map(|k| (p[k] * 1e6).round() as i64);
                    let (ka, kb) = (key(sa), key(sb));
                    let (k, flip) = if ka <= kb { ([ka[0], ka[1], ka[2], kb[0], kb[1], kb[2]], false) } else { ([kb[0], kb[1], kb[2], ka[0], ka[1], ka[2]], true) };
                    let (start, end, q0, q1) = if flip { (sb, sa, qb, qa) } else { (sa, sb, qa, qb) };
                    let entry = by_line.entry(k).or_insert((start, end, DVec3::ZERO, DVec3::ZERO));
                    entry.2 += q0;
                    entry.3 += q1;
                    continue;
                }
            }
            // Not between vertices: the resultant at the centre of the load.
            let length = a.distance(b);
            let (ma, mb) = (qa.length(), qb.length());
            let t = if ma + mb > 1e-15 { (ma + 2. * mb) / (3. * (ma + mb)) } else { 0.5 };
            loads.push(Load::Point {
                case,
                at: a.lerp(b, t).to_array(),
                force: (((qa + qb) / 2. * length) * factor).to_array(),
                moment: [0.; 3],
            });
            *report.approximated.entry("нагрузка по кромке плиты вне вершин модели заменена точечной силой".into()).or_default() += 1;
        }
        for (start, end, q0, q1) in by_line.into_values() {
            loads.push(Load::Line {
                case,
                start: start.to_array(),
                end: end.to_array(),
                q_start: (q0 * factor).to_array(),
                q_end: (q1 * factor).to_array(),
            });
        }
    }

    /// Plate pressures of the combination, one load per plate: uniform over
    /// the whole plate when enough of it is loaded and the centre of
    /// pressure is near its centre, else point loads at the centres of
    /// pressure of the connected loaded parts. No new contours; the
    /// resultant of every plate is kept.
    #[allow(clippy::too_many_arguments)]
    fn simplified_surface_loads(
        &self,
        case: u32,
        elements: &HashMap<u32, Pressure>,
        simplify: Simplify,
        loads: &mut Vec<Load>,
        source: &mut BTreeMap<u32, DVec3>,
        lost: &mut BTreeMap<u32, DVec3>,
        report: &mut Report,
    ) {
        let factor = self.settings.force_factor;
        let model = &self.state.model;
        let mut by_surface: BTreeMap<usize, Vec<(u32, DVec3, f64, DVec3)>> = BTreeMap::new();
        for (&e, pressure) in elements {
            let Some(positions) = element(self.mesh, e).and_then(|el| positions(self.mesh, el)) else { continue };
            let ring = shell_ring(&positions);
            let area = polygon_area(&ring);
            let centre = ring.iter().sum::<DVec3>() / ring.len() as f64;
            let Some(&s) = self.surface_of_element.get(&e) else {
                if pressure.vector.length() > 1e-15 {
                    *report.skipped.entry("нагрузка на пластины, которых нет в геометрии".into()).or_default() += 1;
                    *lost.entry(case).or_default() += pressure.vector * area * factor;
                }
                continue;
            };
            if pressure.vector.length() < 1e-15 {
                continue;
            }
            *source.entry(case).or_default() += pressure.vector * area * factor;
            by_surface.entry(s).or_default().push((e, pressure.vector, area, centre));
        }
        for (s, items) in by_surface {
            let force: DVec3 = items.iter().map(|x| x.1 * x.2).sum();
            let polygons = surface_polygons(model, s, self.settings.snap).0;
            let pieces: Vec<Vec<DVec3>> = polygons
                .iter()
                .map(|p| p.iter().map(|&x| DVec3::from_array(x)).collect())
                .collect();
            let surface_area: f64 = pieces.iter().map(|p| polygon_area(p)).sum();
            if surface_area <= 0. || force.length() < 1e-15 {
                continue;
            }
            // Centre of the plate: area-weighted centres of the fans of its pieces.
            let mut centre_sum = DVec3::ZERO;
            for p in &pieces {
                for i in 1..p.len().saturating_sub(1) {
                    let t = (p[i] - p[0]).cross(p[i + 1] - p[0]).length() / 2.;
                    centre_sum += (p[0] + p[i] + p[i + 1]) / 3. * t;
                }
            }
            let plate_centre = centre_sum / surface_area;
            let direction = force.normalize();
            let weights: Vec<f64> = items.iter().map(|x| (x.1.dot(direction)).max(0.) * x.2).collect();
            let weight_sum: f64 = weights.iter().sum();
            let pressure_centre = if weight_sum > 1e-15 {
                items.iter().zip(&weights).map(|(x, w)| x.3 * *w).sum::<DVec3>() / weight_sum
            } else {
                plate_centre
            };
            let loaded: f64 = items.iter().map(|x| x.2).sum();
            let uniform = loaded / surface_area >= simplify.min_fraction
                && pressure_centre.distance(plate_centre) <= simplify.center_tolerance * surface_area.sqrt();
            if uniform {
                loads.push(Load::Surface {
                    case,
                    surface: s,
                    polygons,
                    sigma: (force / surface_area * factor).to_array(),
                });
                if items.len() > 1 {
                    *report.approximated.entry("давление на пластину приведено к равномерному по всей пластине".into()).or_default() += 1;
                }
                continue;
            }
            // Point loads: one per connected loaded part (elements sharing nodes).
            let mut parent: Vec<usize> = (0..items.len()).collect();
            fn find(parent: &mut [usize], i: usize) -> usize {
                let mut r = i;
                while parent[r] != r {
                    r = parent[r];
                }
                let mut j = i;
                while parent[j] != r {
                    let next = parent[j];
                    parent[j] = r;
                    j = next;
                }
                r
            }
            let mut owner: HashMap<u32, usize> = HashMap::new();
            for (i, x) in items.iter().enumerate() {
                for &n in element(self.mesh, x.0).map(|el| el.nodes.as_slice()).unwrap_or(&[]) {
                    match owner.get(&n) {
                        Some(&j) => {
                            let (a, b) = (find(&mut parent, i), find(&mut parent, j));
                            parent[a] = b;
                        }
                        None => {
                            owner.insert(n, i);
                        }
                    }
                }
            }
            let mut parts: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
            for i in 0..items.len() {
                let r = find(&mut parent, i);
                parts.entry(r).or_default().push(i);
            }
            let plane = &model.planes()[model.surfaces()[s].plane];
            let inside = |p: DVec3| -> bool {
                let q = DVec2::from_array(plane.project(p.to_array()));
                pieces.iter().any(|piece| {
                    let ring: Vec<DVec2> = piece.iter().map(|&v| DVec2::from_array(plane.project(v.to_array()))).collect();
                    let mut odd = false;
                    for i in 0..ring.len() {
                        let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
                        if (a.y > q.y) != (b.y > q.y) && q.x < a.x + (q.y - a.y) / (b.y - a.y) * (b.x - a.x) {
                            odd = !odd;
                        }
                    }
                    odd
                })
            };
            for members in parts.values() {
                let f: DVec3 = members.iter().map(|&i| items[i].1 * items[i].2).sum();
                let w: Vec<f64> = members.iter().map(|&i| items[i].1.dot(direction).max(0.) * items[i].2).collect();
                let ws: f64 = w.iter().sum();
                let mut at = if ws > 1e-15 {
                    members.iter().zip(&w).map(|(&i, w)| items[i].3 * *w).sum::<DVec3>() / ws
                } else {
                    items[members[0]].3
                };
                if !inside(at) {
                    // Outside the plate (an L, a hole): the nearest loaded element centre.
                    at = members
                        .iter()
                        .map(|&i| items[i].3)
                        .min_by(|x, y| x.distance(at).total_cmp(&y.distance(at)))
                        .unwrap_or(at);
                }
                loads.push(Load::Point {
                    case,
                    at: at.to_array(),
                    force: (f * factor).to_array(),
                    moment: [0.; 3],
                });
            }
            *report.approximated.entry("давление на часть пластины заменено точечными силами".into()).or_default() += 1;
        }
    }

    /// The weight of the elements above the cut (tf, downwards), at their
    /// centres: plates by thickness and density, bars by their density.
    fn removed_weight(&self) -> Option<Vec<(DVec3, DVec3, DVec3)>> {
        let cut = self.state.cut.as_ref()?;
        let materials = self.settings.materials.as_ref()?;
        let mut items = vec![];
        for e in &self.mesh.elements {
            let Some(p) = positions(self.mesh, e) else { continue };
            let centre = p.iter().sum::<DVec3>() / p.len() as f64;
            if centre.z <= cut.z + 1e-6 {
                continue;
            }
            let weight = match (e.is_shell(), e.is_bar(), materials.get(&e.stiff_id)) {
                (true, _, Some(Material::Plate { thickness, density: Some(rho), .. })) => rho * thickness * polygon_area(&shell_ring(&p)),
                (_, true, Some(Material::Bar { density: Some(ro), .. })) => ro * (p[1] - p[0]).length(),
                _ => continue,
            };
            items.push((centre, DVec3::new(0., 0., -weight), DVec3::ZERO));
        }
        (!items.is_empty()).then_some(items)
    }

    /// Loads of what was cut off (position, force, moment; tf) as loads on
    /// the supports at the level: every item to the support nearest to it in
    /// plan, a line load along a wall, a point load on a column; the
    /// overturning moment of the items about the supports is kept by a
    /// couple of vertical forces on the supports.
    fn cut_loads(
        &self,
        removed: &BTreeMap<u32, Vec<(DVec3, DVec3, DVec3)>>,
        loads: &mut Vec<Load>,
        source: &mut BTreeMap<u32, DVec3>,
        report: &mut Report,
    ) {
        let Some(cut) = self.state.cut.as_ref().filter(|_| self.settings.cut_loads) else { return };
        let factor = self.settings.force_factor;
        // Supports: (start, end, weight).
        let weight_of = |stiffness: u32, wall: bool| -> f64 {
            match self.settings.materials.as_ref().and_then(|m| m.get(&stiffness)) {
                Some(Material::Plate { thickness, .. }) if wall => *thickness,
                Some(Material::Bar { width, height, .. }) if !wall => width * height,
                _ => if wall { 0.2 } else { 0.25 },
            }
        };
        let mut supports: Vec<(DVec3, DVec3, f64, bool)> = vec![];
        for w in &cut.walls {
            let (a, b) = (DVec3::from_array(w.a), DVec3::from_array(w.b));
            supports.push((a, b, weight_of(w.stiffness, true) * a.distance(b), true));
        }
        for c in &cut.columns {
            let p = DVec3::from_array(c.a);
            supports.push((p, p, weight_of(c.stiffness, false), false));
        }
        supports.retain(|s| s.2 > 0.);
        let plan = |p: DVec3| glam::DVec2::new(p.x, p.y);
        for (&case, items) in removed {
            if items.is_empty() {
                continue;
            }
            if supports.is_empty() {
                *report.skipped.entry("нагрузки отброшенных этажей: нет несущих элементов на уровне отсечения".into()).or_default() += items.len();
                continue;
            }
            let total: DVec3 = items.iter().map(|x| x.1).sum();
            *source.entry(case).or_default() += total * factor;
            let mut shares = vec![DVec3::ZERO; supports.len()];
            let weight_sum: f64 = supports.iter().map(|s| s.2).sum();
            let centre = supports.iter().map(|s| (s.0 + s.1) / 2. * s.2).sum::<DVec3>() / weight_sum;
            let reference = DVec3::new(centre.x, centre.y, cut.z);
            let mut moment = DVec3::ZERO;
            for &(p, f, m) in items {
                let nearest = (0..supports.len())
                    .map(|i| {
                        let (a, b) = (plan(supports[i].0), plan(supports[i].1));
                        let d = b - a;
                        let t = if d.length_squared() > 0. { ((plan(p) - a).dot(d) / d.length_squared()).clamp(0., 1.) } else { 0. };
                        (plan(p).distance(a + d * t), i)
                    })
                    .min_by(|x, y| x.0.total_cmp(&y.0))
                    .map(|x| x.1)
                    .unwrap_or(0);
                shares[nearest] += f;
                moment += (p - reference).cross(f) + m;
            }
            // The moment lost by moving every force to its support.
            let mut assigned = DVec3::ZERO;
            for (i, s) in supports.iter().enumerate() {
                let mid = (s.0 + s.1) / 2.;
                assigned += (DVec3::new(mid.x, mid.y, cut.z) - reference).cross(shares[i]);
            }
            let delta = moment - assigned;
            let (mut sxx, mut syy, mut sxy) = (0., 0., 0.);
            for s in &supports {
                let r = (s.0 + s.1) / 2. - centre;
                sxx += s.2 * r.x * r.x;
                syy += s.2 * r.y * r.y;
                sxy += s.2 * r.x * r.y;
            }
            // f_i = w_i (l1 r_y - l2 r_x) with A l = delta (the moment of the vertical
            // couple about X and Y); least squares, so that supports on one line still
            // take the component of the moment they can.
            let (a11, a12, a21, a22) = (syy, -sxy, -sxy, sxx);
            let (n11, n12, n22) = (a11 * a11 + a21 * a21, a11 * a12 + a21 * a22, a12 * a12 + a22 * a22);
            let rhs = (a11 * delta.x + a21 * delta.y, a12 * delta.x + a22 * delta.y);
            let ridge = 1e-9 * (n11 + n22).max(1e-12);
            let det = (n11 + ridge) * (n22 + ridge) - n12 * n12;
            let (l1, l2) = (((n22 + ridge) * rhs.0 - n12 * rhs.1) / det, ((n11 + ridge) * rhs.1 - n12 * rhs.0) / det);
            let mut left = DVec3::ZERO;
            for (i, s) in supports.iter().enumerate() {
                let r = (s.0 + s.1) / 2. - centre;
                shares[i].z += s.2 * (l1 * r.y - l2 * r.x);
                left += DVec3::new(r.y, -r.x, 0.) * (s.2 * (l1 * r.y - l2 * r.x));
            }
            if (delta.x - left.x).hypot(delta.y - left.y) > 1e-6 * (1. + delta.x.hypot(delta.y)) {
                *report.approximated.entry("опрокидывающий момент отброшенных этажей передан не полностью: опоры на одной прямой".into()).or_default() += 1;
            }
            *report.approximated.entry("нагрузки отброшенных этажей переданы на ближайшие в плане несущие элементы на уровне отсечения".into()).or_default() += 1;
            for (i, s) in supports.iter().enumerate() {
                if shares[i].length() < 1e-12 {
                    continue;
                }
                if s.3 {
                    let q = shares[i] / s.0.distance(s.1);
                    loads.push(Load::Line { case, start: s.0.to_array(), end: s.1.to_array(), q_start: (q * factor).to_array(), q_end: (q * factor).to_array() });
                } else {
                    loads.push(Load::Point { case, at: s.0.to_array(), force: (shares[i] * factor).to_array(), moment: [0.; 3] });
                }
            }
        }
    }

    /// Stamp forces become pressures. The rows of one case, direction and
    /// level with equal forces make one stamp: its total force is spread
    /// evenly over the area of the elements it covers (so the pressure of
    /// the stamp is uniform and the resultant is kept).
    fn spread_stamps(
        &self,
        stamps: &[Stamp],
        pressure: &mut BTreeMap<u32, HashMap<u32, Pressure>>,
        report: &mut Report,
    ) {
        const FORCE_TOLERANCE: f64 = 1e-3;
        let quantize = |x: f64, step: f64| (x / step).round() as i64;
        let mut groups: BTreeMap<(u32, [i64; 3], i64), Vec<&Stamp>> = BTreeMap::new();
        for stamp in stamps {
            let key = (
                stamp.case,
                [0, 1, 2].map(|k| quantize(stamp.direction[k], 1e-6)),
                quantize(stamp.level, 1e-3),
            );
            groups.entry(key).or_default().push(stamp);
        }
        for group in groups.values_mut() {
            group.sort_by(|a, b| a.value.total_cmp(&b.value).then(a.element.cmp(&b.element)));
            let mut start = 0;
            while start < group.len() {
                let first = group[start].value;
                let end = (start..group.len())
                    .find(|&i| (group[i].value - first).abs() > FORCE_TOLERANCE)
                    .unwrap_or(group.len());
                let cluster = &group[start..end];
                let force: f64 = cluster.iter().map(|s| s.value).sum();
                let mut elements: BTreeMap<u32, f64> = BTreeMap::new();
                for s in cluster {
                    elements.insert(s.element, s.area);
                }
                let area: f64 = elements.values().sum();
                if area > 0. {
                    let intensity = force / area;
                    for &element in elements.keys() {
                        pressure
                            .entry(cluster[0].case)
                            .or_default()
                            .entry(element)
                            .or_default()
                            .vector += cluster[0].direction * intensity;
                    }
                    if cluster.len() > 1 {
                        *report
                            .approximated
                            .entry("сосредоточенные силы на пластинах (штампы) заменены равномерным давлением на их элементы".into())
                            .or_default() += 1;
                    }
                } else {
                    skip(report, "сосредоточенная сила на пластине нулевой площади");
                }
                start = end;
            }
        }
    }

    fn node_row(
        &self,
        row: &LoadRow,
        params: &[f64],
        node_loads: &mut BTreeMap<(u32, u32), (DVec3, DVec3)>,
        report: &mut Report,
    ) {
        let (Some(&value), Some(axis)) = (params.first(), unit(row.direction)) else {
            skip(report, format!("нагрузка на узел, направление {}", row.direction));
            return;
        };
        // Loads on a node are in the node's local axes (document 12).
        let vector = match self.set.node_axes.get(&row.target) {
            Some([x, y]) => {
                let z = x.cross(*y);
                let local = [*x, *y, z][(row.direction as usize - 1) % 3];
                local * value
            }
            None => axis * value,
        };
        let vector = if row.direction <= 3 { vector * FORCE_SIGN } else { vector };
        let entry = node_loads.entry((row.case, row.target)).or_default();
        if row.direction <= 3 {
            entry.0 += vector;
        } else {
            entry.1 += vector;
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn shell_row(
        &self,
        row: &LoadRow,
        params: &[f64],
        p: &[DVec3],
        axes: &[DVec3; 3],
        area: f64,
        pressure: &mut BTreeMap<u32, HashMap<u32, Pressure>>,
        edge_lines: &mut BTreeMap<u32, Vec<(DVec3, DVec3, DVec3, DVec3)>>,
        stamps: &mut Vec<Stamp>,
        report: &mut Report,
    ) {
        let code = row.code;
        let base = code % 10;
        let fr = frame(code);
        if row.direction == 0 || row.direction > 3 {
            skip(report, format!("момент на пластину (код {code}, направление {})", row.direction));
            return;
        }
        // The direction as a vector, with the factor of projected loads.
        let direction = match fr {
            0 => axes[row.direction as usize - 1],
            _ => unit(row.direction).unwrap(),
        };
        let direction = direction * FORCE_SIGN;
        if fr == 4 {
            *report.approximated.entry("нагрузки «выравнивания» приняты глобальными".into()).or_default() += 1;
        }
        let normal = axes[2];
        let projected = |dir: DVec3| if fr == 2 { dir * normal.dot(dir).abs() } else { dir };
        let add = |pressure: &mut BTreeMap<u32, HashMap<u32, Pressure>>, v: DVec3| {
            pressure.entry(row.case).or_default().entry(row.target).or_default().vector += v;
        };
        match (base, code) {
            // Uniform pressure.
            (6, _) => {
                let Some(&value) = params.first() else {
                    return skip(report, "нагрузка без параметров");
                };
                add(pressure, projected(direction) * value);
            }
            // Four values at the nodes: the mean (the resultant is kept).
            (7, _) => {
                if params.len() < p.len().min(4) {
                    return skip(report, "нагрузка без параметров");
                }
                let n = p.len().min(4);
                let mean = params[..n].iter().sum::<f64>() / n as f64;
                if params[..n].iter().any(|v| (v - mean).abs() > 1e-9 * mean.abs().max(1e-12)) {
                    *report.approximated.entry("трапециевидное давление на пластину заменено средним по элементу".into()).or_default() += 1;
                }
                add(pressure, projected(direction) * mean);
            }
            // A point force on a plate is a stamp: the forces of its elements
            // are spread over them by `spread_stamps`.
            (5, _) if params.len() <= 4 => {
                let Some(&value) = params.first() else {
                    return skip(report, "нагрузка без параметров");
                };
                stamps.push(Stamp {
                    case: row.case,
                    element: row.target,
                    value,
                    direction: direction.normalize_or_zero(),
                    area,
                    level: (p.iter().sum::<DVec3>() / p.len() as f64).dot(direction.normalize_or_zero()),
                });
            }
            (5, _) => skip(report, "произвольная трапециевидная нагрузка на пластину"),
            // Line loads along an element edge.
            (9 | 0, 9 | 19 | 29 | 49 | 10 | 20 | 30 | 50) => {
                let (Some(&a), Some(&b)) = (params.first(), params.get(1)) else {
                    return skip(report, "нагрузка по линии без параметров");
                };
                let (ia, ib) = (a as usize, b as usize);
                if ia == 0 || ib == 0 || ia > p.len() || ib > p.len() || ia == ib {
                    return skip(report, "нагрузка по линии с неверными узлами");
                }
                let (pa, pb) = (p[ia - 1], p[ib - 1]);
                let (qa, qb) = if base == 9 {
                    (params.get(2).copied().unwrap_or(0.), params.get(2).copied().unwrap_or(0.))
                } else {
                    (params.get(2).copied().unwrap_or(0.), params.get(3).copied().unwrap_or(0.))
                };
                let along = (pb - pa).normalize();
                let dir = if fr == 2 { direction * (1. - along.dot(direction).powi(2)).max(0.).sqrt() } else { direction };
                edge_lines.entry(row.case).or_default().push((pa, pb, dir * qa, dir * qb));
            }
            _ => skip(report, format!("нагрузка на пластину, код {code}")),
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn bar_row(
        &self,
        row: &LoadRow,
        params: &[f64],
        e: &ElementData,
        uniform: &mut BTreeMap<(u32, u32), DVec3>,
        lines: &mut BTreeMap<(u32, usize), Vec<(f64, f64, DVec3, DVec3)>>,
        loads: &mut Vec<Load>,
        source: &mut BTreeMap<u32, DVec3>,
        removed: &mut BTreeMap<u32, Vec<(DVec3, DVec3, DVec3)>>,
        report: &mut Report,
    ) {
        let factor = self.settings.force_factor;
        let code = row.code;
        let base = code % 10;
        let fr = frame(code);
        let Some(p) = positions(self.mesh, e).filter(|p| p.len() == 2) else {
            return skip(report, "стержень без узлов");
        };
        let y = self
            .set
            .element_axes
            .get(&e.id)
            .filter(|v| v.len() >= 3)
            .map(|v| DVec3::new(v[0], v[1], v[2]));
        let Some(axes) = bar_axes(p[0], p[1], y) else {
            return skip(report, "стержень нулевой длины");
        };
        let length = (p[1] - p[0]).length();
        if row.direction == 0 || row.direction > 6 {
            return skip(report, format!("нагрузка на стержень, направление {}", row.direction));
        }
        let moment = row.direction > 3;
        let dir = if fr == 0 {
            axes[(row.direction as usize - 1) % 3]
        } else {
            unit(row.direction).unwrap()
        };
        let dir = if moment { dir } else { dir * FORCE_SIGN };
        // A bar cut off with the upper storeys: its load goes to the supports.
        if self.state.cut.as_ref().is_some_and(|c| (p[0] + p[1]).z / 2. > c.z + 1e-6) {
            let along = axes[0];
            let projected = |d: DVec3| if fr == 2 { d * (1. - along.dot(d).powi(2)).max(0.).sqrt() } else { d };
            if moment {
                return skip(report, "момент на стержень, отброшенный вместе с верхними этажами");
            }
            let entry = match base {
                5 => params.first().zip(params.get(1)).map(|(&v, &a)| (p[0] + along * a.clamp(0., length), dir * v)),
                6 => params.first().map(|&v| ((p[0] + p[1]) / 2., projected(dir) * v * length)),
                7 => match (params.first(), params.get(1), params.get(2), params.get(3)) {
                    (Some(&p1), Some(&a1), Some(&p2), Some(&a2)) => {
                        Some((p[0] + along * ((a1 + a2) / 2.).clamp(0., length), projected(dir) * ((p1 + p2) / 2.) * (a2 - a1).abs()))
                    }
                    _ => None,
                },
                _ => None,
            };
            match entry {
                Some((at, force)) => removed.entry(row.case).or_default().push((at, force, DVec3::ZERO)),
                None => skip(report, format!("нагрузка на стержень, код {code}")),
            }
            return;
        }
        let Some(&(axis, t_first, t_second)) = self.bar_of_element.get(&e.id) else {
            return skip(report, "нагрузка на стержень, которого нет в геометрии");
        };
        let [a, b] = self.state.axes[axis].endpoints.map(|v| DVec3::from_array(self.state.model.vertices()[v]));
        let t_at = |s: f64| t_first + (t_second - t_first) * (s / length).clamp(0., 1.);
        let point_at = |s: f64| a.lerp(b, t_at(s));
        let along = axes[0];
        let projected = |d: DVec3| if fr == 2 { d * (1. - along.dot(d).powi(2)).max(0.).sqrt() } else { d };
        match base {
            5 => {
                let (Some(&value), Some(&s)) = (params.first(), params.get(1)) else {
                    return skip(report, "нагрузка без параметров");
                };
                let v = dir * value * factor;
                if !moment {
                    *source.entry(row.case).or_default() += v;
                }
                loads.push(Load::Point {
                    case: row.case,
                    at: point_at(s).to_array(),
                    force: if moment { [0.; 3] } else { v.to_array() },
                    moment: if moment { v.to_array() } else { [0.; 3] },
                });
            }
            6 => {
                if moment {
                    return skip(report, "распределённый момент на стержень");
                }
                let Some(&value) = params.first() else {
                    return skip(report, "нагрузка без параметров");
                };
                *uniform.entry((row.case, e.id)).or_default() += projected(dir) * value;
            }
            7 => {
                if moment {
                    return skip(report, "распределённый момент на стержень");
                }
                let (Some(&p1), Some(&a1), Some(&p2), Some(&a2)) =
                    (params.first(), params.get(1), params.get(2), params.get(3))
                else {
                    return skip(report, "нагрузка без параметров");
                };
                let (ta, tb) = (t_at(a1), t_at(a2));
                let (qa, qb) = (projected(dir) * p1, projected(dir) * p2);
                if (ta - tb).abs() < 1e-9 {
                    return skip(report, "трапеция нулевой длины");
                }
                // Source resultant, counted here (it is not uniform).
                *source.entry(row.case).or_default() += (qa + qb) / 2. * (a2 - a1).abs() * factor;
                if ta < tb {
                    lines.entry((row.case, axis)).or_default().push((ta, tb, qa, qb));
                } else {
                    lines.entry((row.case, axis)).or_default().push((tb, ta, qb, qa));
                }
            }
            _ => skip(report, format!("нагрузка на стержень, код {code}")),
        }
    }

    /// Surface loads of one case from the pressures of its shell elements.
    fn surface_loads(
        &self,
        case: u32,
        elements: &HashMap<u32, Pressure>,
        loads: &mut Vec<Load>,
        source: &mut BTreeMap<u32, DVec3>,
        lost: &mut BTreeMap<u32, DVec3>,
        report: &mut Report,
    ) {
        let factor = self.settings.force_factor;
        let model = &self.state.model;
        // Elements by surface.
        let mut by_surface: BTreeMap<usize, Vec<(u32, DVec3, f64)>> = BTreeMap::new();
        for (&e, pressure) in elements {
            let Some(&s) = self.surface_of_element.get(&e) else {
                // Not in the geometry (absorbed, removed or deleted).
                let area = element(self.mesh, e)
                    .and_then(|el| positions(self.mesh, el))
                    .map_or(0., |p| polygon_area(&shell_ring(&p)));
                *report.skipped.entry("нагрузка на пластины, которых нет в геометрии".into()).or_default() += 1;
                *lost.entry(case).or_default() += pressure.vector * area * factor;
                continue;
            };
            let area = element(self.mesh, e)
                .and_then(|el| positions(self.mesh, el))
                .map_or(0., |p| polygon_area(&shell_ring(&p)));
            *source.entry(case).or_default() += pressure.vector * area * factor;
            by_surface.entry(s).or_default().push((e, pressure.vector, area));
        }
        for (s, mut items) in by_surface {
            if items.iter().all(|(_, v, _)| v.length() < 1e-15) {
                continue;
            }
            items.sort_by_key(|x| x.0);
            let shells_in_surface = model.surfaces()[s]
                .source_elements
                .iter()
                .filter(|&&e| element(self.mesh, e).is_some_and(|el| el.is_shell()))
                .count();
            // Groups of equal value (binned when there are too many).
            let key = |v: DVec3, step: f64| [0, 1, 2].map(|k| (v[k] / step).round() as i64);
            let mut step = 1e-9;
            let mut groups: BTreeMap<[i64; 3], Vec<usize>> = BTreeMap::new();
            for (i, (_, v, _)) in items.iter().enumerate() {
                groups.entry(key(*v, step)).or_default().push(i);
            }
            if groups.len() > self.settings.max_groups {
                let top = items.iter().map(|x| x.1.abs().max_element()).fold(0., f64::max);
                step = top / self.settings.max_groups as f64;
                groups.clear();
                for (i, (_, v, _)) in items.iter().enumerate() {
                    groups.entry(key(*v, step)).or_default().push(i);
                }
                *report.approximated.entry("давления на поверхности объединены в диапазоны значений (равнодействующая сохранена)".into()).or_default() += 1;
            }
            for members in groups.values() {
                // Area-weighted mean value: the resultant of the group.
                let area: f64 = members.iter().map(|&i| items[i].2).sum();
                if area <= 0. {
                    continue;
                }
                let sigma = members.iter().map(|&i| items[i].1 * items[i].2).sum::<DVec3>() / area;
                if sigma.length() < 1e-15 {
                    continue;
                }
                let whole = members.len() == shells_in_surface;
                let polygons = if whole {
                    surface_polygons(model, s, self.settings.snap).0
                } else {
                    let ids: Vec<u32> = members.iter().map(|&i| items[i].0).collect();
                    self.region(s, &ids)
                };
                if polygons.is_empty() {
                    continue;
                }
                loads.push(Load::Surface {
                    case,
                    surface: s,
                    polygons,
                    sigma: (sigma * factor).to_array(),
                });
            }
        }
    }

    /// The polygons (without holes) of the union of shell elements of
    /// surface `s`, clipped to the surface.
    fn region(&self, s: usize, ids: &[u32]) -> Vec<Vec<[f64; 3]>> {
        let model = &self.state.model;
        let surface = &model.surfaces()[s];
        let plane = &model.planes()[surface.plane];
        let uv = |p: DVec3| DVec2::from_array(plane.project(p.to_array()));
        // Directed edges of the element rings, counter-clockwise in the
        // plane; edges used twice in opposite directions cancel.
        let mut edges: HashMap<(u32, u32), ()> = HashMap::new();
        let mut at: HashMap<u32, DVec2> = HashMap::new();
        for &id in ids {
            let Some(el) = element(self.mesh, id) else { continue };
            let Some(p) = positions(self.mesh, el) else { continue };
            let perimeter = shell_ring_nodes(&p, &el.nodes);
            let ring: Vec<DVec2> = perimeter.iter().map(|&(_, x)| uv(x)).collect();
            let signed: f64 = (0..ring.len())
                .map(|i| ring[i].perp_dot(ring[(i + 1) % ring.len()]))
                .sum();
            let mut nodes: Vec<u32> = perimeter.iter().map(|&(n, _)| n).collect();
            if signed < 0. {
                nodes.reverse();
            }
            for (&(n, _), &point) in perimeter.iter().zip(&ring) {
                at.insert(n, point);
            }
            for i in 0..nodes.len() {
                let (a, b) = (nodes[i], nodes[(i + 1) % nodes.len()]);
                if edges.remove(&(b, a)).is_none() {
                    edges.insert((a, b), ());
                }
            }
        }
        // Boundary loops.
        let mut next: HashMap<u32, Vec<u32>> = HashMap::new();
        for &(a, b) in edges.keys() {
            next.entry(a).or_default().push(b);
        }
        let mut starts: Vec<(u32, u32)> = edges.keys().copied().collect();
        starts.sort_unstable();
        let mut loops: Vec<Vec<DVec2>> = vec![];
        for (a, b) in starts {
            // Already used?
            let Some(list) = next.get_mut(&a) else { continue };
            let Some(position) = list.iter().position(|&x| x == b) else { continue };
            list.swap_remove(position);
            let mut ring = vec![at[&a]];
            let (start, mut current) = (a, b);
            let mut guard = 0;
            while current != start && guard < 1_000_000 {
                guard += 1;
                ring.push(at[&current]);
                let Some(list) = next.get_mut(&current) else { break };
                let Some(n) = list.pop() else { break };
                current = n;
            }
            if ring.len() >= 3 && current == start {
                loops.push(ring);
            }
        }
        // Snap the outline onto the contour of the surface.
        let snap = self.settings.snap;
        let contour_segments: Vec<(DVec2, DVec2)> = surface
            .contours
            .iter()
            .flat_map(|r| (0..r.len()).map(move |i| (DVec2::from_array(r[i]), DVec2::from_array(r[(i + 1) % r.len()]))))
            .collect();
        for ring in &mut loops {
            for point in ring.iter_mut() {
                let mut best: Option<(f64, DVec2)> = None;
                for &(a, b) in &contour_segments {
                    let d = b - a;
                    let t = ((*point - a).dot(d) / d.length_squared()).clamp(0., 1.);
                    let q = a + d * t;
                    let dist = point.distance(q);
                    if dist <= snap && best.is_none_or(|x| dist < x.0) {
                        best = Some((dist, q));
                    }
                }
                if let Some((_, q)) = best {
                    *point = q;
                }
            }
        }
        // Outer loops (counter-clockwise) with the holes they contain.
        let signed = |r: &[DVec2]| -> f64 { (0..r.len()).map(|i| r[i].perp_dot(r[(i + 1) % r.len()])).sum::<f64>() / 2. };
        let to_line = |r: &[DVec2]| LineString::from(r.iter().map(|p| Coord { x: p.x, y: p.y }).collect::<Vec<_>>());
        let outers: Vec<&Vec<DVec2>> = loops.iter().filter(|r| signed(r) > 0.).collect();
        let holes: Vec<&Vec<DVec2>> = loops.iter().filter(|r| signed(r) < 0.).collect();
        let mut region = MultiPolygon::new(vec![]);
        for outer in &outers {
            let inside = |hole: &&Vec<DVec2>| {
                let p = hole[0];
                let mut odd = false;
                for i in 0..outer.len() {
                    let (a, b) = (outer[i], outer[(i + 1) % outer.len()]);
                    if (a.y > p.y) != (b.y > p.y) && p.x < a.x + (p.y - a.y) / (b.y - a.y) * (b.x - a.x) {
                        odd = !odd;
                    }
                }
                odd
            };
            let own: Vec<LineString<f64>> = holes.iter().filter(|h| inside(h)).map(|h| to_line(h)).collect();
            region.0.push(Polygon::new(to_line(outer), own));
        }
        let surface_polygon = Polygon::new(
            LineString::from(surface.contours[0].iter().map(|p| Coord { x: p[0], y: p[1] }).collect::<Vec<_>>()),
            surface.contours[1..]
                .iter()
                .map(|r| LineString::from(r.iter().map(|p| Coord { x: p[0], y: p[1] }).collect::<Vec<_>>()))
                .collect(),
        );
        let clipped = surface_polygon.intersection(&region);
        let threshold = (self.settings.snap * self.settings.snap).max(1e-6);
        let mut out = vec![];
        for polygon in clipped {
            if polygon.unsigned_area() < threshold {
                continue;
            }
            let ring = |l: &LineString<f64>| -> Vec<[f64; 2]> {
                let mut points: Vec<[f64; 2]> = l.coords().map(|c| [c.x, c.y]).collect();
                points.pop();
                points
            };
            let mut contours = vec![ring(polygon.exterior())];
            contours.extend(polygon.interiors().iter().map(ring));
            for piece in hole_free(&contours, self.settings.snap) {
                out.push(plaxis_polygon(piece.into_iter().map(|p| plane.lift(p)).collect()));
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shell_perimeter_order_is_found_for_tensor_numbering_and_collinear_nodes() {
        let q = |a: [f64; 2], b, c, d| [a, b, c, d].map(|p: [f64; 2]| DVec3::new(p[0], p[1], 0.));
        // Perimeter order, tensor-product order and a crossed order of a unit square.
        let square = q([0., 0.], [1., 0.], [1., 1.], [0., 1.]);
        assert_eq!(shell_order(&square), vec![0, 1, 2, 3]);
        let tensor = q([0., 0.], [1., 0.], [0., 1.], [1., 1.]);
        assert_eq!(shell_order(&tensor), vec![0, 1, 3, 2]);
        let crossed = q([0., 0.], [1., 1.], [1., 0.], [0., 1.]);
        assert_eq!(shell_order(&crossed), vec![0, 2, 1, 3]);
        // A node in the middle of a straight side (a transition element) keeps the listed order...
        let hanging = q([0., 0.], [2., 0.], [2., 1.], [1., 1.]);
        let line = [DVec3::ZERO, DVec3::new(2., 0., 0.), DVec3::new(2., 1., 0.), DVec3::new(0., 1., 0.)];
        assert_eq!(shell_order(&hanging), vec![0, 1, 2, 3]);
        assert_eq!(shell_order(&line), vec![0, 1, 2, 3]);
        // ... also in tensor numbering.
        let tensor_hanging = [DVec3::ZERO, DVec3::new(1., 0., 0.), DVec3::new(0., 1., 0.), DVec3::new(2., 0., 0.)];
        assert_eq!(shell_order(&tensor_hanging).len(), 4);
        assert!((polygon_area(&shell_ring(&tensor)) - 1.).abs() < 1e-12);
    }

    #[test]
    fn self_weight_cases_are_recognized_by_name() {
        for name in ["СВ", " св ", "СВ_1.1", "СОБСТВЕННЫЙ ВЕС", "Собственный вес плиты перекрытия 11-го этажа", "Self weight"] {
            assert!(is_self_weight(name), "{name}");
        }
        for name in ["СВЕТ", "СНЕГ_1.4|0.5", "ПОЛЫ НОРМ", "СТАДИЯ 1", "ПОКРЫТИЕ_1.3", ""] {
            assert!(!is_self_weight(name), "{name}");
        }
    }

    #[test]
    fn stage_and_dynamic_cases_are_recognized_by_name() {
        for name in ["СТАДИЯ 1", "СТАДИЯ №10", "СТАДИЯ 10 <ФИНАЛЬНАЯ СТАДИЯ>", "Stage 2"] {
            assert!(is_stage(name), "{name}");
        }
        for name in ["ПОЛЫ", "СНЕГ", "ПОЛЕЗНАЯ НАГРУЗКА_ K=1_2"] {
            assert!(!is_stage(name), "{name}");
        }
        assert!(is_dynamic("СЕЙСМИКА ПО Х") && !is_dynamic("ВЕТЕР"));
    }

    #[test]
    fn polygon_area_of_a_concave_outline_is_not_overcounted() {
        // An L: the fan from its first node covers the notch with a negative triangle.
        let l = [(0., 0.), (2., 0.), (2., 1.), (1., 1.), (1., 2.), (0., 2.)].map(|(x, y)| DVec3::new(x, y, 0.));
        assert!((polygon_area(&l) - 3.).abs() < 1e-12);
    }
}
