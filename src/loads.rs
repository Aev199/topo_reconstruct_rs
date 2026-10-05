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

#[derive(Debug, Clone, Copy)]
pub struct Settings {
    /// Force unit of the model in kN (LIRA tf: 9.80665).
    pub force_factor: f64,
    /// Snap distance of region outlines to the surface contour.
    pub snap: f64,
    /// Most surface loads of one surface and case: values are binned.
    pub max_groups: usize,
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
        let mut shell_geometry: HashMap<u32, Option<(Vec<DVec3>, [DVec3; 3], f64)>> = HashMap::new();

        for row in &self.set.rows {
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
                    self.shell_row(row, params, &p, &axes, area, &mut pressure, &mut edge_lines, &mut report);
                }
                10 => self.bar_row(row, params, e, &mut bar_uniform, &mut bar_lines, &mut loads, &mut source, &mut report),
                other => skip(&mut report, format!("нагрузка на элемент типа {other}")),
            }
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
            self.surface_loads(case, elements, &mut loads, &mut source, &mut lost, &mut report);
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
        let names: BTreeMap<u32, &String> = self.set.cases.iter().map(|(n, s)| (*n, s)).collect();
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
            // A point force: spread over its element (the resultant is kept).
            (5, _) if params.len() <= 4 => {
                let Some(&value) = params.first() else {
                    return skip(report, "нагрузка без параметров");
                };
                *report.approximated.entry("сосредоточенная сила на пластине заменена давлением на её элемент".into()).or_default() += 1;
                add(pressure, direction * value / area);
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
        let Some(&(axis, t_first, t_second)) = self.bar_of_element.get(&e.id) else {
            return skip(report, "нагрузка на стержень, которого нет в геометрии");
        };
        if row.direction == 0 || row.direction > 6 {
            return skip(report, format!("нагрузка на стержень, направление {}", row.direction));
        }
        let moment = row.direction > 3;
        let dir = if fr == 0 {
            axes[(row.direction as usize - 1) % 3]
        } else {
            unit(row.direction).unwrap()
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
    fn polygon_area_of_a_concave_outline_is_not_overcounted() {
        // An L: the fan from its first node covers the notch with a negative triangle.
        let l = [(0., 0.), (2., 0.), (2., 1.), (1., 1.), (1., 2.), (0., 2.)].map(|(x, y)| DVec3::new(x, y, 0.));
        assert!((polygon_area(&l) - 3.).abs() < 1e-12);
    }
}
