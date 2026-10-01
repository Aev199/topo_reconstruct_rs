//! Global synchronization of mesh constraints before local surface meshing.
//!
//! Every interval endpoint is materialized once on its axis.  Boundary edge
//! chains are then built from those same vertices and are shared by every
//! surface owning the edge.  This keeps surface constraints and bar chains
//! conforming without welding vertices by proximity.

use super::{diagnostic, parameter, point, sorted, subdivide, Policy, ENDPOINT_SLACK};
use crate::reconstruction::{
    assembly::bars::{Axis, Contact},
    Model,
};
use glam::DVec3;
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Serialize)]
pub struct Report {
    pub shared_edge_count: usize,
    pub interval_contact_count: usize,
    pub synchronized_interval_endpoints: usize,
    pub endpoint_bindings: Vec<EndpointBinding>,
    pub axis_node_count: usize,
    pub edge_node_count: usize,
}

#[derive(Debug, Serialize)]
pub struct EndpointBinding {
    pub axis: usize,
    pub surface: usize,
    pub contact: usize,
    pub role: EndpointRole,
    pub parameter: f64,
    pub vertex: usize,
    pub generated: bool,
}

#[derive(Debug, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum EndpointRole {
    Start,
    End,
}

pub(super) struct Synchronized {
    pub axis_nodes: Vec<Vec<(f64, usize)>>,
    pub edge_nodes: Vec<Vec<(f64, usize)>>,
    pub report: Report,
}

fn axis_geometry(axis: &Axis, vertices: &[[f64; 3]]) -> Result<(DVec3, DVec3, f64), &'static str> {
    let [a, b] = axis.endpoints;
    let (Some(a), Some(b)) = (vertices.get(a), vertices.get(b)) else {
        return Err("invalid axis endpoint vertex");
    };
    let a = point(a);
    let b = point(b);
    let length = a.distance(b);
    if !length.is_finite() || length == 0. {
        return Err("degenerate axis constraint");
    }
    Ok((a, b, length))
}

fn normalized_parameter(t: f64, tolerance: f64) -> Result<f64, &'static str> {
    if !t.is_finite() || t < -tolerance || t > 1. + tolerance {
        return Err("axis interval parameter outside endpoints");
    }
    Ok(t.clamp(0., 1.))
}

fn insert_axis_node(
    chain: &mut Vec<(f64, usize)>,
    t: f64,
    vertex: usize,
    vertices: &[[f64; 3]],
    tolerance: f64,
    axis_start: DVec3,
    axis_end: DVec3,
    precision: f64,
) -> Result<bool, &'static str> {
    let axis = axis_end - axis_start;
    let p = point(
        vertices
            .get(vertex)
            .ok_or("invalid axis constraint vertex")?,
    );
    if p.distance(axis_start + axis * t) > precision {
        return Err("axis constraint vertex is off axis");
    }
    if let Some((u, _existing)) = chain.iter().find(|(_, n)| *n == vertex) {
        if (*u - t).abs() > tolerance {
            return Err("axis constraint vertex has inconsistent parameter");
        }
        return Ok(false);
    }
    if let Some((_, existing)) = chain.iter().find(|(u, _)| (*u - t).abs() <= tolerance) {
        if *existing != vertex {
            diagnostic(|| {
                format!(
                    "axis vertices {existing} {:?} and {vertex} {:?} at parameter {t}",
                    vertices[*existing], vertices[vertex]
                )
            });
            return Err("distinct axis constraint vertices share a parameter");
        }
        return Ok(false);
    }
    chain.push((t, vertex));
    Ok(true)
}

fn validate_surface(surface: usize, model: &Model) -> Result<(), &'static str> {
    model
        .surfaces
        .get(surface)
        .map(|_| ())
        .ok_or("invalid contact surface")
}

/// Local feature size for shared constraint edges. A constraint vertex
/// chain is fixed during surface refinement, so two nearby constraint lines
/// subdivided at the nominal spacing would enclose unavoidable slivers. The
/// size at a point is its distance to the nearest non-adjacent edge of any
/// surface owning the constraint, capped by the nominal spacing. Every owner
/// receives the same graded chain, so conformity is unchanged.
struct Sizing {
    segments: Vec<(usize, [usize; 2], DVec3, DVec3)>,
    by_surface: Vec<Vec<usize>>,
    spacing: f64,
}

impl Sizing {
    fn new(model: &Model, vertices: &[[f64; 3]], spacing: f64) -> Self {
        let segments = model
            .edges
            .iter()
            .enumerate()
            .map(|(e, &[a, b])| (e, [a, b], point(&vertices[a]), point(&vertices[b])))
            .collect();
        let by_surface = (0..model.surfaces.len())
            .map(|s| model.surface_edges(s).collect())
            .collect();
        Self {
            segments,
            by_surface,
            spacing,
        }
    }

    /// Edges of the owning surfaces within the nominal spacing of `edge`.
    fn features(&self, edge: usize, owners: &BTreeSet<usize>) -> Vec<(DVec3, DVec3)> {
        let (_, ends, a, b) = self.segments[edge];
        let low = a.min(b) - DVec3::splat(self.spacing);
        let high = a.max(b) + DVec3::splat(self.spacing);
        let mut seen = BTreeSet::new();
        let mut out = vec![];
        for &s in owners {
            for &f in &self.by_surface[s] {
                let (_, [c, d], pc, pd) = self.segments[f];
                if f == edge || ends.contains(&c) || ends.contains(&d) || !seen.insert(f) {
                    continue;
                }
                if pc.max(pd).cmplt(low).any() || pc.min(pd).cmpgt(high).any() {
                    continue;
                }
                out.push((pc, pd));
            }
        }
        out
    }

    fn size(&self, p: DVec3, features: &[(DVec3, DVec3)]) -> f64 {
        let mut h = self.spacing;
        for &(a, b) in features {
            let d = b - a;
            let t = ((p - a).dot(d) / d.length_squared()).clamp(0., 1.);
            h = h.min(p.distance(a + d * t));
        }
        h.max(self.spacing * MINIMUM_SIZE_RATIO)
    }
}

/// Lower bound of the graded size relative to the nominal spacing.
const MINIMUM_SIZE_RATIO: f64 = 1e-3;

/// Subdivide a chain with point density `1 / size`. For a constant size this
/// is exactly the uniform subdivision; the placement depends only on the
/// geometry of the edge, not on its orientation or vertex ids.
fn graded(
    chain: &[(f64, usize)],
    vertices: &mut Vec<[f64; 3]>,
    size: impl Fn(DVec3) -> f64,
    spacing: f64,
    limit: usize,
    precision: f64,
) -> Result<Vec<(f64, usize)>, &'static str> {
    let mut out = vec![chain[0]];
    for pair in chain.windows(2) {
        let [(ta, a), (tb, b)] = [pair[0], pair[1]];
        let pa = point(&vertices[a]);
        let pb = point(&vertices[b]);
        let length = pa.distance(pb);
        // Sample the density; the size is 1-Lipschitz, so a quarter of the
        // local size resolves it. Symmetric sampling keeps orientation
        // independence.
        let mut samples = vec![(0., 1. / size(pa))];
        let mut s = 0.;
        let mut uniform = samples[0].1 * spacing <= 1. + 1e-12;
        while s < length {
            let h = size(pa.lerp(pb, s / length));
            s = (s + (h / 4.).max(precision)).min(length);
            let density = 1. / size(pa.lerp(pb, s / length));
            uniform &= density * spacing <= 1. + 1e-12;
            samples.push((s, density));
            if samples.len() > 16 * limit.max(1) {
                return Err("constraint subdivision limit");
            }
        }
        let count = if uniform {
            ((length - precision) / spacing).ceil().max(1.)
        } else {
            let mut total = 0.;
            let mut cumulative = vec![0.];
            for w in samples.windows(2) {
                total += (w[1].0 - w[0].0) * (w[0].1 + w[1].1) / 2.;
                cumulative.push(total);
            }
            let count = (total - precision / spacing).ceil().max(1.);
            if count.is_finite() && count <= limit as f64 {
                let mut k = 0;
                for i in 1..count as usize {
                    let target = total * i as f64 / count;
                    while cumulative[k + 1] < target {
                        k += 1;
                    }
                    let (c0, c1) = (cumulative[k], cumulative[k + 1]);
                    let x = samples[k].0
                        + (samples[k + 1].0 - samples[k].0) * (target - c0) / (c1 - c0);
                    let u = x / length;
                    let id = vertices.len();
                    vertices.push(pa.lerp(pb, u).to_array());
                    out.push((ta + (tb - ta) * u, id));
                }
                out.push((tb, b));
                continue;
            }
            count
        };
        if !count.is_finite() || count > limit as f64 {
            return Err("constraint subdivision limit");
        }
        let count = count as usize;
        for i in 1..count {
            let u = i as f64 / count as f64;
            let id = vertices.len();
            vertices.push(pa.lerp(pb, u).to_array());
            out.push((ta + (tb - ta) * u, id));
        }
        out.push((tb, b));
    }
    Ok(out)
}

pub(super) fn synchronize(
    model: &Model,
    axes: &[Axis],
    contacts: &[Contact],
    vertices: &mut Vec<[f64; 3]>,
    policy: &Policy,
    precision: f64,
) -> std::result::Result<Synchronized, &'static str> {
    let mut axis_nodes = Vec::with_capacity(axes.len());
    let mut geometries = Vec::with_capacity(axes.len());
    for axis in axes {
        let geometry = axis_geometry(axis, vertices)?;
        let mut chain: Vec<_> = axis.anchors.iter().map(|a| (a.t, a.vertex)).collect();
        sorted(&mut chain, vertices, precision)?;
        if chain.len() < 2 {
            return Err("empty axis constraint");
        }
        geometries.push(geometry);
        axis_nodes.push(chain);
    }

    let mut endpoint_bindings = Vec::new();
    let mut generated_ends = BTreeMap::<[i64; 3], Vec<usize>>::new();
    let mut interval_contact_count = 0;
    for (contact, item) in contacts.iter().enumerate() {
        match item {
            Contact::Point {
                axis,
                surface,
                vertex,
                t,
                ..
            } => {
                validate_surface(*surface, model)?;
                let (start, end, length) =
                    *geometries.get(*axis).ok_or("invalid point contact axis")?;
                let tolerance = precision / length;
                let t = normalized_parameter(*t, tolerance)?;
                insert_axis_node(
                    &mut axis_nodes[*axis],
                    t,
                    *vertex,
                    vertices,
                    tolerance,
                    start,
                    end,
                    precision,
                )?;
            }
            Contact::Interval {
                axis,
                surface,
                start_t,
                end_t,
                ..
            } => {
                validate_surface(*surface, model)?;
                let (start, end, length) = *geometries
                    .get(*axis)
                    .ok_or("invalid interval contact axis")?;
                let tolerance = precision / length;
                let start_t = normalized_parameter(*start_t, tolerance)?;
                let end_t = normalized_parameter(*end_t, tolerance)?;
                if end_t - start_t <= tolerance {
                    return Err("degenerate interval contact");
                }
                interval_contact_count += 1;
                // An axis node within the endpoint slack of an interval end
                // (a bar crossing computed separately at a slab edge) is
                // that end.
                let snap = ENDPOINT_SLACK * tolerance;
                for (role, t) in [(EndpointRole::Start, start_t), (EndpointRole::End, end_t)] {
                    let existing = axis_nodes[*axis]
                        .iter()
                        .filter(|(u, _)| (*u - t).abs() <= snap)
                        .min_by(|x, y| (x.0 - t).abs().total_cmp(&(y.0 - t).abs()))
                        .copied();
                    let (t, vertex, generated) = if let Some((u, vertex)) = existing {
                        (u, vertex, false)
                    } else {
                        // The end of another bar's interval at the same point
                        // (two bars leaving a slab at one point) is shared.
                        let p = start.lerp(end, t);
                        let radius = ENDPOINT_SLACK * precision;
                        let cell = |p: DVec3| (p / radius).floor().as_i64vec3().to_array();
                        let c = cell(p);
                        let mut shared = None;
                        'search: for dx in -1..=1 {
                            for dy in -1..=1 {
                                for dz in -1..=1 {
                                    for &w in generated_ends
                                        .get(&[c[0] + dx, c[1] + dy, c[2] + dz])
                                        .into_iter()
                                        .flatten()
                                    {
                                        // On this bar too, within precision.
                                        if point(&vertices[w]).distance(p) <= precision {
                                            shared = Some(w);
                                            break 'search;
                                        }
                                    }
                                }
                            }
                        }
                        match shared {
                            Some(w) => (t, w, false),
                            None => {
                                let vertex = vertices.len();
                                vertices.push(p.to_array());
                                generated_ends.entry(c).or_default().push(vertex);
                                (t, vertex, true)
                            }
                        }
                    };
                    insert_axis_node(
                        &mut axis_nodes[*axis],
                        t,
                        vertex,
                        vertices,
                        tolerance,
                        start,
                        end,
                        precision,
                    )?;
                    endpoint_bindings.push(EndpointBinding {
                        axis: *axis,
                        surface: *surface,
                        contact,
                        role,
                        parameter: t,
                        vertex,
                        generated,
                    });
                }
            }
        }
    }

    // Overlapping bars (a source bar along part of a longer one) run through
    // the same two vertices: that piece is subdivided once and its vertices
    // are shared, keyed by the piece's ends.
    let mut pieces = BTreeMap::<[usize; 2], Vec<(f64, usize)>>::new();
    for chain in &mut axis_nodes {
        sorted(chain, vertices, precision)?;
        // A synchronized interval endpoint is part of the global chain before
        // spacing is applied, so every owner receives the same split vertex.
        let mut expanded = vec![chain[0]];
        for w in chain.windows(2) {
            let [(ta, a), (tb, b)] = [w[0], w[1]];
            let key = [a.min(b), a.max(b)];
            if !pieces.contains_key(&key) {
                let piece = subdivide(
                    &[(0., key[0]), (1., key[1])],
                    vertices,
                    policy.boundary_spacing,
                    policy.maximum_added_vertices_per_surface,
                    precision,
                )?;
                pieces.insert(key, piece[1..piece.len() - 1].to_vec());
            }
            let inner = &pieces[&key];
            let along = |u: f64| if a == key[0] { u } else { 1. - u };
            let mut inner: Vec<_> = inner
                .iter()
                .map(|&(u, v)| (ta + (tb - ta) * along(u), v))
                .collect();
            inner.sort_by(|x, y| x.0.total_cmp(&y.0));
            expanded.extend(inner);
            expanded.push((tb, b));
        }
        *chain = expanded;
    }

    let mut owners = vec![BTreeSet::new(); model.edges.len()];
    for surface in 0..model.surfaces.len() {
        for edge in model.surface_edges(surface) {
            let owner = owners
                .get_mut(edge)
                .ok_or("invalid surface edge reference")?;
            owner.insert(surface);
        }
    }
    let shared_edge_count = owners.iter().filter(|owner| owner.len() > 1).count();
    let sizing = Sizing::new(model, vertices, policy.boundary_spacing);
    // Contacts of each surface, in contact order (an edge only collects
    // those of its owners; scanning all contacts per edge was quadratic).
    let mut by_surface = vec![vec![]; model.surfaces.len()];
    for (k, item) in contacts.iter().enumerate() {
        let (Contact::Point { surface, .. } | Contact::Interval { surface, .. }) = item;
        by_surface[*surface].push(k);
    }
    let mut edge_nodes = Vec::with_capacity(model.edges.len());
    for (edge_id, &[a, b]) in model.edges.iter().enumerate() {
        let pa = point(vertices.get(a).ok_or("invalid model edge vertex")?);
        let pb = point(vertices.get(b).ok_or("invalid model edge vertex")?);
        let mut chain = vec![(0., a), (1., b)];
        let mut owned: Vec<usize> = owners[edge_id]
            .iter()
            .flat_map(|&s| by_surface[s].iter().copied())
            .collect();
        owned.sort_unstable();
        for item in owned.iter().map(|&k| &contacts[k]) {
            let (surface, candidates) = match *item {
                Contact::Point {
                    vertex, surface, ..
                } => (surface, vec![vertex]),
                Contact::Interval {
                    axis,
                    surface,
                    start_t,
                    end_t,
                    ..
                } => (
                    surface,
                    axis_nodes[axis]
                        .iter()
                        .filter(|(t, _)| {
                            let slack = ENDPOINT_SLACK * precision / geometries[axis].2;
                            *t >= start_t - slack && *t <= end_t + slack
                        })
                        .map(|(_, vertex)| *vertex)
                        .collect(),
                ),
            };
            if !owners[edge_id].contains(&surface) {
                continue;
            }
            for vertex in candidates {
                if let Some(t) = parameter(point(&vertices[vertex]), pa, pb, precision) {
                    chain.push((t, vertex));
                }
            }
        }
        if owners[edge_id].is_empty() {
            // An orphaned edge (left by a merge or a removed crack void) is
            // meshed only as a simplified-hole constraint; its ends may have
            // been merged onto one point.
            edge_nodes.push(chain);
            continue;
        }
        sorted(&mut chain, vertices, precision)?;
        let features = sizing.features(edge_id, &owners[edge_id]);
        edge_nodes.push(graded(
            &chain,
            vertices,
            |p| sizing.size(p, &features),
            policy.boundary_spacing,
            policy.maximum_added_vertices_per_surface,
            precision,
        )?);
    }

    // A subdivision introduced on a common boundary is also a subdivision of
    // every bar interval lying on that boundary.  Insert it into the one
    // global axis chain, not into a surface-local copy.
    for item in contacts {
        let Contact::Interval {
            axis,
            surface,
            start_t,
            end_t,
            ..
        } = *item
        else {
            continue;
        };
        let [a, b] = axes[axis].endpoints;
        for edge in model.surface_edges(surface) {
            for &(_, vertex) in &edge_nodes[edge] {
                if let Some(t) = parameter(
                    point(&vertices[vertex]),
                    point(&vertices[a]),
                    point(&vertices[b]),
                    precision,
                ) {
                    let length = point(&vertices[a]).distance(point(&vertices[b]));
                    let tolerance = precision / length;
                    if t >= start_t - ENDPOINT_SLACK * tolerance
                        && t <= end_t + ENDPOINT_SLACK * tolerance
                    {
                        insert_axis_node(
                            &mut axis_nodes[axis],
                            t,
                            vertex,
                            vertices,
                            tolerance,
                            point(&vertices[a]),
                            point(&vertices[b]),
                            precision,
                        )?;
                    }
                }
            }
        }
    }
    for chain in &mut axis_nodes {
        sorted(chain, vertices, precision)?;
    }

    let axis_node_count = axis_nodes.iter().map(Vec::len).sum();
    let edge_node_count = edge_nodes.iter().map(Vec::len).sum();
    Ok(Synchronized {
        axis_nodes,
        edge_nodes,
        report: Report {
            shared_edge_count,
            interval_contact_count,
            synchronized_interval_endpoints: endpoint_bindings.len(),
            endpoint_bindings,
            axis_node_count,
            edge_node_count,
        },
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reconstruction::{
        assembly::bars::{Anchor, Location},
        PlaneFrame,
    };

    fn shared_edge_case() -> (Model, Vec<[f64; 3]>, Vec<Axis>, Vec<Contact>) {
        let mut model = Model::new(1e-8, 0.01).unwrap();
        let plane = model.add_plane(PlaneFrame::new([0., 0., 0.], [0., 0., 1.]).unwrap());
        let points = vec![
            [0., 0., 0.],
            [4., 0., 0.],
            [4., 2., 0.],
            [0., 2., 0.],
            [4., -2., 0.],
            [0., -2., 0.],
        ];
        let vertices: Vec<_> = points
            .iter()
            .map(|&p| model.add_vertex(p).unwrap())
            .collect();
        model
            .add_surface(plane, vec![vec![0, 1, 2, 3]], vec![1])
            .unwrap();
        model
            .add_surface(plane, vec![vec![1, 0, 5, 4]], vec![2])
            .unwrap();
        let axis = Axis {
            source_axis: 0,
            endpoints: [vertices[0], vertices[1]],
            anchors: vec![
                Anchor {
                    source_node: 10,
                    vertex: vertices[0],
                    t: 0.,
                },
                Anchor {
                    source_node: 11,
                    vertex: vertices[1],
                    t: 1.,
                },
            ],
            spans: vec![],
        };
        let contacts = vec![
            Contact::Interval {
                axis: 0,
                surface: 0,
                start_t: 0.25,
                end_t: 0.75,
                location: Location::Boundary,
            },
            Contact::Interval {
                axis: 0,
                surface: 1,
                start_t: 0.25,
                end_t: 0.75,
                location: Location::Boundary,
            },
        ];
        (model, points, vec![axis], contacts)
    }

    fn policy() -> Policy {
        Policy {
            boundary_spacing: 10.,
            maximum_area: 1.,
            minimum_angle_degrees: 20.,
            maximum_added_vertices_per_surface: 32,
        }
    }

    #[test]
    fn interval_endpoints_are_created_once_and_shared_by_surfaces_and_bar() {
        let (model, mut vertices, axes, contacts) = shared_edge_case();
        let synced = synchronize(&model, &axes, &contacts, &mut vertices, &policy(), 1e-8).unwrap();

        assert_eq!(synced.report.shared_edge_count, 1);
        assert_eq!(synced.report.interval_contact_count, 2);
        assert_eq!(synced.report.synchronized_interval_endpoints, 4);
        assert_eq!(synced.report.endpoint_bindings.len(), 4);
        assert_eq!(
            synced
                .report
                .endpoint_bindings
                .iter()
                .filter(|binding| binding.generated)
                .count(),
            2
        );
        let first = &synced.report.endpoint_bindings[0..2];
        let second = &synced.report.endpoint_bindings[2..4];
        assert_eq!(first[0].vertex, second[0].vertex);
        assert_eq!(first[1].vertex, second[1].vertex);
        assert_eq!(
            synced.axis_nodes[0]
                .iter()
                .map(|(_, vertex)| *vertex)
                .collect::<Vec<_>>(),
            vec![0, first[0].vertex, first[1].vertex, 1]
        );
        let common_edge = synced.edge_nodes[0]
            .iter()
            .map(|(_, vertex)| *vertex)
            .collect::<Vec<_>>();
        assert_eq!(common_edge, vec![0, first[0].vertex, first[1].vertex, 1]);
        assert_eq!(synced.report.axis_node_count, 4);
        // The nominal spacing exceeds the panels. The 4 m outer edges are 2 m
        // from the opposite shared edge, so graded sizing splits each once.
        assert_eq!(synced.report.edge_node_count, 18);
    }

    #[test]
    fn overlapping_bars_share_the_subdivision_of_their_common_piece() {
        // A bar from x = 0 to 4 and a longer one from x = 4 back to -2 through
        // x = 0 (overlapping source bars): their common piece gets one set of
        // subdivision vertices, whichever direction each bar runs.
        let (model, mut vertices, mut axes, _) = shared_edge_case();
        vertices.push([-2., 0., 0.]);
        let far = vertices.len() - 1;
        let anchor = |source_node, vertex, t| Anchor {
            source_node,
            vertex,
            t,
        };
        axes.push(Axis {
            source_axis: 1,
            endpoints: [1, far],
            anchors: vec![
                anchor(11, 1, 0.),
                anchor(10, 0, 2. / 3.),
                anchor(12, far, 1.),
            ],
            spans: vec![],
        });
        let mut p = policy();
        p.boundary_spacing = 1.;
        let synced = synchronize(&model, &axes, &[], &mut vertices, &p, 1e-8).unwrap();
        let on_piece = |k: usize| -> BTreeSet<usize> {
            synced.axis_nodes[k]
                .iter()
                .map(|&(_, v)| v)
                .filter(|&v| (0. ..=4.).contains(&vertices[v][0]))
                .collect()
        };
        assert_eq!(on_piece(0).len(), 5);
        assert_eq!(on_piece(0), on_piece(1));
        for chain in &synced.axis_nodes {
            assert!(chain.windows(2).all(|w| w[0].0 < w[1].0));
        }
    }

    #[test]
    fn axis_node_next_to_an_interval_end_is_that_end() {
        // A bar node 50 nm before the interval start (within the endpoint
        // slack of 10 precisions): the start binds to it, no vertex is
        // generated.
        let (model, mut vertices, mut axes, contacts) = shared_edge_case();
        vertices.push([1. - 5e-8, 0., 0.]);
        let node = vertices.len() - 1;
        axes[0].anchors.insert(
            1,
            Anchor {
                source_node: 12,
                vertex: node,
                t: (1. - 5e-8) / 4.,
            },
        );
        let synced = synchronize(&model, &axes, &contacts, &mut vertices, &policy(), 1e-8).unwrap();
        let starts: Vec<_> = synced
            .report
            .endpoint_bindings
            .iter()
            .filter(|b| matches!(b.role, EndpointRole::Start))
            .collect();
        assert!(starts.iter().all(|b| b.vertex == node && !b.generated));
        assert!(starts.iter().all(|b| b.parameter == (1. - 5e-8) / 4.));
    }

    #[test]
    fn bars_leaving_a_surface_at_one_point_share_the_interval_end() {
        // Two bars crossing the common edge y = 0 of both panels at (2, 0):
        // one interval end vertex there, used by both bars and the edge.
        let (model, mut vertices, _, _) = shared_edge_case();
        let mut axes = vec![];
        let mut contacts = vec![];
        for (k, (a, b)) in [([2., -1., 0.], [2., 1., 0.]), ([1., -1., 0.], [3., 1., 0.])]
            .into_iter()
            .enumerate()
        {
            vertices.push(a);
            vertices.push(b);
            let (va, vb) = (vertices.len() - 2, vertices.len() - 1);
            axes.push(Axis {
                source_axis: k,
                endpoints: [va, vb],
                anchors: vec![
                    Anchor {
                        source_node: 20 + k as u32,
                        vertex: va,
                        t: 0.,
                    },
                    Anchor {
                        source_node: 30 + k as u32,
                        vertex: vb,
                        t: 1.,
                    },
                ],
                spans: vec![],
            });
            for (surface, start_t, end_t) in [(1, 0., 0.5), (0, 0.5, 1.)] {
                contacts.push(Contact::Interval {
                    axis: k,
                    surface,
                    start_t,
                    end_t,
                    location: Location::Interior,
                });
            }
        }
        let synced = synchronize(&model, &axes, &contacts, &mut vertices, &policy(), 1e-8).unwrap();
        let at = |k: usize| {
            synced.axis_nodes[k]
                .iter()
                .find(|(t, _)| (t - 0.5).abs() < 1e-12)
                .unwrap()
                .1
        };
        assert_eq!(at(0), at(1));
        assert!(synced.edge_nodes[0].iter().any(|&(_, v)| v == at(0)));
    }

    #[test]
    fn malformed_interval_is_rejected_before_any_mesh_constraint_is_added() {
        let (model, mut vertices, axes, mut contacts) = shared_edge_case();
        contacts[0] = Contact::Interval {
            axis: 0,
            surface: 0,
            start_t: 0.75,
            end_t: 0.75,
            location: Location::Interior,
        };
        assert!(matches!(
            synchronize(&model, &axes, &contacts, &mut vertices, &policy(), 1e-8),
            Err("degenerate interval contact")
        ));
        assert_eq!(vertices.len(), 6);
    }
}
