//! Assemble two-endpoint axes against fixed surface geometry. Contacts are
//! geometric constraints for later meshing, never inferred mechanical ties.
use super::{frame, intersection, movement_budget, MeshData, Model, PlaneFrame, Policy};
use crate::reconstruction::recognize::SourceSpan;
use geo::{
    line_intersection::{line_intersection, LineIntersection},
    Contains, Line, LineString, Point, Polygon,
};
use glam::{DVec2, DVec3};
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Serialize)]
pub struct Anchor {
    pub source_node: u32,
    pub vertex: usize,
    pub t: f64,
}
#[derive(Debug, Serialize)]
pub struct Axis {
    pub source_axis: usize,
    pub endpoints: [usize; 2],
    pub anchors: Vec<Anchor>,
    pub spans: Vec<SourceSpan>,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Location {
    Boundary,
    Interior,
}
#[derive(Debug, Serialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum Contact {
    Point {
        axis: usize,
        surface: usize,
        vertex: usize,
        t: f64,
        location: Location,
    },
    Interval {
        axis: usize,
        surface: usize,
        start_t: f64,
        end_t: f64,
        location: Location,
    },
}
#[derive(Debug, Serialize)]
pub struct Issue {
    pub source_axis: usize,
    pub source_elements: Vec<u32>,
    pub source_nodes: Vec<u32>,
    pub reason: String,
}
#[derive(Debug, Default, Serialize)]
pub struct Report {
    pub all_axes_built: bool,
    pub axes: Vec<Axis>,
    pub contacts: Vec<Contact>,
    pub issues: Vec<Issue>,
    pub rejected_shared_anchors: BTreeMap<u32, String>,
    pub maximum_additional_movement: f64,
    /// No intersection-driven edge subdivision or mesh generation at this stage.
    pub mesh_constraints_complete: bool,
}

fn polygon(contours: &[Vec<[f64; 2]>]) -> Polygon<f64> {
    let ring = |r: &Vec<[f64; 2]>| {
        LineString::from(
            r.iter()
                .chain(r.first())
                .map(|p| (p[0], p[1]))
                .collect::<Vec<_>>(),
        )
    };
    Polygon::new(
        ring(&contours[0]),
        contours.iter().skip(1).map(ring).collect(),
    )
}
fn location(uv: [f64; 2], contours: &[Vec<[f64; 2]>], precision: f64) -> Option<Location> {
    let p = DVec2::from_array(uv);
    for ring in contours {
        for i in 0..ring.len() {
            let a = DVec2::from_array(ring[i]);
            let d = DVec2::from_array(ring[(i + 1) % ring.len()]) - a;
            let q = a + d * ((p - a).dot(d) / d.length_squared()).clamp(0., 1.);
            if p.distance(q) <= precision {
                return Some(Location::Boundary);
            }
        }
    }
    polygon(contours)
        .contains(&Point::new(p.x, p.y))
        .then_some(Location::Interior)
}

/// Clip a coplanar axis to the material, including holes and nonconvex outlines.
fn intervals(
    a: [f64; 2],
    b: [f64; 2],
    contours: &[Vec<[f64; 2]>],
    precision: f64,
) -> Vec<(f64, f64, Location)> {
    let a = DVec2::from_array(a);
    let b = DVec2::from_array(b);
    let d = b - a;
    let length = d.length();
    if length <= precision {
        return vec![];
    }
    let line = Line::new((a.x, a.y), (b.x, b.y));
    let parameter =
        |x: geo::Coord<f64>| ((DVec2::new(x.x, x.y) - a).dot(d) / d.length_squared()).clamp(0., 1.);
    let mut cuts = vec![0., 1.];
    for ring in contours {
        for i in 0..ring.len() {
            let p = ring[i];
            let q = ring[(i + 1) % ring.len()];
            match line_intersection(line, Line::new((p[0], p[1]), (q[0], q[1]))) {
                Some(LineIntersection::SinglePoint { intersection, .. }) => {
                    cuts.push(parameter(intersection))
                }
                Some(LineIntersection::Collinear { intersection }) => {
                    cuts.push(parameter(intersection.start));
                    cuts.push(parameter(intersection.end));
                }
                None => {}
            }
        }
    }
    cuts.sort_by(f64::total_cmp);
    cuts.dedup_by(|a, b| (*a - *b).abs() * length <= precision);
    let mut result: Vec<(f64, f64, Location)> = vec![];
    for w in cuts.windows(2) {
        if let Some(kind) = location(
            (a + d * ((w[0] + w[1]) / 2.)).to_array(),
            contours,
            precision,
        ) {
            if let Some(last) = result.last_mut() {
                if last.2 == kind && (last.1 - w[0]).abs() * length <= precision {
                    last.1 = w[1];
                    continue;
                }
            }
            result.push((w[0], w[1], kind));
        }
    }
    result
}

struct Proposal {
    points: BTreeMap<u32, DVec3>,
    parameters: BTreeMap<u32, f64>,
    spans: Vec<SourceSpan>,
}

fn propose(
    mesh: &MeshData,
    source: &frame::Report,
    axis: &frame::Axis,
    locked: &BTreeMap<u32, DVec3>,
    unavailable: &BTreeMap<u32, &str>,
    owners: &BTreeMap<u32, Vec<usize>>,
    supports: &[PlaneFrame],
    policy: &Policy,
) -> Result<Proposal, (&'static str, Vec<u32>)> {
    let fail = |reason, nodes| Err((reason, nodes));
    let mut anchors = axis.anchors.clone();
    anchors.sort_by(|a, b| a.t.total_cmp(&b.t));
    let ids: Vec<_> = anchors.iter().map(|a| source.node_ids[a.node]).collect();
    let bad: Vec<_> = ids
        .iter()
        .filter(|n| unavailable.contains_key(n))
        .copied()
        .collect();
    if !bad.is_empty() {
        return fail("unavailable_shared_anchor", bad);
    }
    let [i, j] = axis.endpoints;
    let ends = [source.node_ids[i], source.node_ids[j]];
    if !ends.iter().all(|n| ids.contains(n)) || axis.spans.is_empty() {
        return fail("incomplete_axis_definition", ends.to_vec());
    }
    let reference = mesh.nodes[&ends[1]] - mesh.nodes[&ends[0]];
    let reference_length = reference.length();
    if reference_length <= policy.precision {
        return fail("degenerate_reference_axis", ends.to_vec());
    }
    let up = DVec3::from_array(source.policy.up).normalize();
    let fixed: Vec<_> = ids
        .iter()
        .filter_map(|n| locked.get(n).map(|&p| (*n, p)))
        .collect();
    let ca = DVec3::from_array(source.candidate_points[i]);
    let cb = DVec3::from_array(source.candidate_points[j]);
    let short = reference_length < source.policy.minimum_length;
    let raw = if fixed.len() >= 2 {
        fixed.last().unwrap().1 - fixed[0].1
    } else {
        cb - ca
    };
    if raw.length() <= policy.precision {
        return fail("coincident_axis_anchors", ends.to_vec());
    }
    let cosine = reference.normalize().dot(up).abs();
    let direction = if short {
        reference.normalize()
    } else if cosine >= source.policy.angle.cos() {
        up * reference.dot(up).signum()
    } else if cosine <= source.policy.angle.sin() {
        (raw - up * raw.dot(up)).normalize_or_zero()
    } else {
        raw.normalize()
    };
    if direction.length_squared() < 0.5
        || direction.dot(reference.normalize()) < source.policy.angle.cos()
    {
        return fail("axis_direction_conflict", ends.to_vec());
    }
    let mut origin = if let Some((_, p)) = fixed.first() {
        *p
    } else {
        (ca + cb) * 0.5
    };
    if fixed.is_empty() {
        let mut common: Option<BTreeSet<usize>> = None;
        for n in &ids {
            let owned: BTreeSet<_> = owners.get(n).into_iter().flatten().copied().collect();
            common = Some(match common {
                None => owned,
                Some(s) => s.intersection(&owned).copied().collect(),
            });
        }
        let parallel: Vec<_> = common
            .unwrap_or_default()
            .into_iter()
            .map(|i| &supports[i])
            .filter(|p| DVec3::from_array(p.normal).dot(direction).abs() <= 1e-10)
            .collect();
        origin = intersection(origin, &parallel, policy.precision)
            .ok_or(("axis_surface_incidence_conflict", ids.clone()))?;
    }
    let mut points = BTreeMap::new();
    for a in &anchors {
        let n = source.node_ids[a.node];
        let candidate = DVec3::from_array(source.candidate_points[a.node]);
        let q = if short {
            let (base, t) = if let Some((n, p)) = fixed.first() {
                (
                    *p,
                    anchors
                        .iter()
                        .find(|a| source.node_ids[a.node] == *n)
                        .unwrap()
                        .t,
                )
            } else {
                ((ca + cb) * 0.5, 0.5)
            };
            base + reference * (a.t - t)
        } else {
            origin + direction * (candidate - origin).dot(direction)
        };
        let q = if let Some(&fixed) = locked.get(&n) {
            let off = if short {
                q.distance(fixed)
            } else {
                (fixed - origin).cross(direction).length()
            };
            if off > policy.precision {
                return fail("shared_anchors_not_collinear", vec![n]);
            }
            fixed
        } else if !short {
            // Interior source FE nodes may slide along the geometric axis.
            // Intersect that line with all owning planes; never pin their
            // arbitrary tangential coordinates from the approximate solver.
            let mut t = (candidate - origin).dot(direction);
            let owned: Vec<_> = owners
                .get(&n)
                .into_iter()
                .flatten()
                .map(|&i| &supports[i])
                .collect();
            for plane in &owned {
                let slope = DVec3::from_array(plane.normal).dot(direction);
                if slope.abs() > 1e-10 {
                    t = -plane.distance(origin.to_array()) / slope;
                    break;
                }
            }
            let p = origin + direction * t;
            if owned
                .iter()
                .any(|plane| plane.distance(p.to_array()).abs() > policy.precision)
            {
                return fail("axis_surface_incidence_conflict", vec![n]);
            }
            p
        } else {
            q
        };
        if owners
            .get(&n)
            .into_iter()
            .flatten()
            .any(|&i| supports[i].distance(q.to_array()).abs() > policy.precision)
        {
            return fail("axis_surface_incidence_conflict", vec![n]);
        }
        if !q.is_finite()
            || q.distance(mesh.nodes[&n]) > movement_budget(mesh, source, a.node) + policy.precision
            || q.distance(candidate) > policy.junction_movement_limit
        {
            return fail("axis_movement_budget", vec![n]);
        }
        points.insert(n, q);
    }
    let (a, b) = (points[&ends[0]], points[&ends[1]]);
    let d = b - a;
    let length = d.length();
    if length + policy.precision < source.policy.minimum_length.min(reference_length)
        || d.dot(reference) <= 0.
    {
        return fail("axis_length_or_orientation", ends.to_vec());
    }
    let parameters: BTreeMap<_, _> = points
        .iter()
        .map(|(&n, &p)| (n, (p - a).dot(d) / d.length_squared()))
        .collect();
    for w in ids.windows(2) {
        if (parameters[&w[1]] - parameters[&w[0]]) * length <= policy.precision {
            return fail("collapsed_or_reordered_property_span", w.to_vec());
        }
    }
    if points
        .values()
        .any(|p| (*p - a).cross(d / length).length() > policy.precision)
    {
        return fail("axis_straightness", ids);
    }
    let mut spans = vec![];
    for span in &axis.spans {
        let start = anchors
            .iter()
            .find(|a| a.t == span.start_t)
            .ok_or(("missing_span_anchor", vec![]))?;
        let end = anchors
            .iter()
            .find(|a| a.t == span.end_t)
            .ok_or(("missing_span_anchor", vec![]))?;
        spans.push(SourceSpan {
            element: span.element,
            stiffness: span.stiffness,
            start_t: parameters[&source.node_ids[start.node]],
            end_t: parameters[&source.node_ids[end.node]],
        });
    }
    Ok(Proposal {
        points,
        parameters,
        spans,
    })
}

pub(super) fn assemble(
    mesh: &MeshData,
    source: &frame::Report,
    model: &mut Model,
    vertex_source_nodes: &mut Vec<u32>,
    supports: &[PlaneFrame],
    representatives: &[usize],
    policy: &Policy,
) -> Report {
    let mut report = Report::default();
    let mut vertices: BTreeMap<_, _> = vertex_source_nodes
        .iter()
        .enumerate()
        .map(|(i, &n)| (n, i))
        .collect();
    let mut locked: BTreeMap<_, _> = vertices
        .iter()
        .map(|(&n, &v)| (n, DVec3::from_array(model.vertices[v])))
        .collect();
    let mut unavailable = BTreeMap::new();
    let mut incidence = BTreeMap::<u32, BTreeSet<usize>>::new();
    for (i, axis) in source.axes.iter().enumerate() {
        for a in &axis.anchors {
            incidence
                .entry(source.node_ids[a.node])
                .or_default()
                .insert(i);
        }
    }
    let mut owner_planes = BTreeMap::<u32, Vec<usize>>::new();
    for (i, surface) in source.surfaces.iter().enumerate() {
        for &node in &surface.nodes {
            let n = source.node_ids[node];
            if incidence.contains_key(&n) {
                owner_planes.entry(n).or_default().push(representatives[i]);
            }
        }
    }
    let element_surface: BTreeMap<_, _> = model
        .surfaces
        .iter()
        .enumerate()
        .flat_map(|(i, s)| s.source_elements.iter().map(move |&e| (e, i)))
        .collect();
    let source_shells: BTreeSet<_> = source
        .surfaces
        .iter()
        .flat_map(|s| s.source_elements.iter().copied())
        .collect();
    let mut owner_surfaces = BTreeMap::<u32, BTreeSet<usize>>::new();
    for e in &mesh.elements {
        if !source_shells.contains(&e.id) {
            continue;
        }
        for &n in &e.nodes {
            if !incidence.contains_key(&n) {
                continue;
            }
            if let Some(&s) = element_surface.get(&e.id) {
                owner_surfaces.entry(n).or_default().insert(s);
            } else {
                unavailable.insert(n, "unbuilt_surface_region");
            }
        }
    }
    let lookup: BTreeMap<_, _> = source
        .node_ids
        .iter()
        .enumerate()
        .map(|(i, &n)| (n, i))
        .collect();
    for (&n, axes) in &incidence {
        if locked.contains_key(&n) {
            continue;
        }
        let candidate = DVec3::from_array(source.candidate_points[lookup[&n]]);
        if let Some(owned) = owner_planes.get(&n) {
            let planes: Vec<_> = owned.iter().map(|&p| &supports[p]).collect();
            if let Some(q) = intersection(candidate, &planes, policy.precision) {
                if axes.len() > 1 {
                    locked.insert(n, q);
                }
            } else {
                unavailable.insert(n, "incompatible_surface_supports");
            }
        } else if axes.len() > 1 {
            locked.insert(n, candidate);
        }
    }
    // Shared junctions are fixed before any axis is considered. Accepted axes
    // cannot move a neighbor or depend on processing order.
    for (source_axis, axis) in source.axes.iter().enumerate() {
        let result = (|| {
            let proposal = propose(
                mesh,
                source,
                axis,
                &locked,
                &unavailable,
                &owner_planes,
                supports,
                policy,
            )?;
            for (&n, &p) in &proposal.points {
                for &s in owner_surfaces.get(&n).into_iter().flatten() {
                    let surface = &model.surfaces[s];
                    let plane = &model.planes[surface.plane];
                    if plane.distance(p.to_array()).abs() > policy.precision
                        || location(
                            plane.project(p.to_array()),
                            &surface.contours,
                            policy.precision,
                        )
                        .is_none()
                    {
                        return Err(("anchor_outside_source_surface", vec![n]));
                    }
                }
            }
            Ok(proposal)
        })();
        let proposal = match result {
            Ok(p) => p,
            Err((reason, nodes)) => {
                report.issues.push(Issue {
                    source_axis,
                    source_elements: axis.spans.iter().map(|s| s.element).collect(),
                    source_nodes: nodes,
                    reason: reason.into(),
                });
                continue;
            }
        };
        let output_axis = report.axes.len();
        for (&n, &p) in &proposal.points {
            report.maximum_additional_movement = report
                .maximum_additional_movement
                .max(p.distance(DVec3::from_array(source.candidate_points[lookup[&n]])));
            if !vertices.contains_key(&n) {
                let v = model
                    .add_vertex(p.to_array())
                    .expect("validated finite axis point");
                vertices.insert(n, v);
                vertex_source_nodes.push(n);
            }
        }
        let endpoints = axis.endpoints.map(|i| vertices[&source.node_ids[i]]);
        let (a, b) = (model.vertices[endpoints[0]], model.vertices[endpoints[1]]);
        let mut touched = BTreeSet::new();
        for anchor in &axis.anchors {
            let n = source.node_ids[anchor.node];
            for &s in owner_surfaces.get(&n).into_iter().flatten() {
                let surface = &model.surfaces[s];
                let plane = &model.planes[surface.plane];
                report.contacts.push(Contact::Point {
                    axis: output_axis,
                    surface: s,
                    vertex: vertices[&n],
                    t: proposal.parameters[&n],
                    location: location(
                        plane.project(proposal.points[&n].to_array()),
                        &surface.contours,
                        policy.precision,
                    )
                    .unwrap(),
                });
                touched.insert(s);
            }
        }
        for s in touched {
            let surface = &model.surfaces[s];
            let plane = &model.planes[surface.plane];
            if plane.distance(a).abs() <= policy.precision
                && plane.distance(b).abs() <= policy.precision
            {
                for (start_t, end_t, location) in intervals(
                    plane.project(a),
                    plane.project(b),
                    &surface.contours,
                    policy.precision,
                ) {
                    report.contacts.push(Contact::Interval {
                        axis: output_axis,
                        surface: s,
                        start_t,
                        end_t,
                        location,
                    });
                }
            }
        }
        report.axes.push(Axis {
            source_axis,
            endpoints,
            spans: proposal.spans,
            anchors: axis
                .anchors
                .iter()
                .map(|a| {
                    let n = source.node_ids[a.node];
                    Anchor {
                        source_node: n,
                        vertex: vertices[&n],
                        t: proposal.parameters[&n],
                    }
                })
                .collect(),
        });
    }
    report.rejected_shared_anchors = unavailable
        .into_iter()
        .map(|(n, r)| (n, r.into()))
        .collect();
    report.all_axes_built = report.issues.is_empty();
    report
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        input::ElementData,
        reconstruction::{assembly, planes, recognize},
    };

    fn frame(mesh: &MeshData, up: DVec3) -> frame::Report {
        let axes = recognize::recognize(
            mesh,
            &recognize::Policy {
                angle: 0.02,
                line_tolerance: 0.01,
                numerical_precision: 1e-8,
            },
        )
        .unwrap();
        let planes = planes::recognize(
            mesh,
            &planes::Policy {
                angle: 0.02,
                distance: 0.01,
                precision: 1e-8,
            },
        )
        .unwrap();
        frame::solve(
            mesh,
            &axes,
            &planes,
            &frame::Policy {
                up: up.to_array(),
                angle: 0.02,
                maximum_movement: 0.15,
                relative_movement: 0.05,
                minimum_length: 0.03,
                residual_tolerance: 1e-7,
                iterations: 100,
            },
        )
        .unwrap()
    }
    fn policy(scale: f64) -> Policy {
        Policy {
            closure_tolerance: 0.001 * scale,
            junction_movement_limit: 0.05 * scale,
            precision: 1e-7 * scale,
            minimum_edge: 0.001 * scale,
        }
    }
    fn slab_beam_column() -> MeshData {
        let mut mesh = MeshData::default();
        for y in 0..3 {
            for x in 0..3 {
                mesh.nodes
                    .insert(1 + x + 3 * y, DVec3::new(x as f64, y as f64, 0.));
            }
        }
        for y in 0..2 {
            for x in 0..2 {
                let a = 1 + x + 3 * y;
                mesh.elements.push(ElementData {
                    id: mesh.elements.len() as u32 + 1,
                    elem_type: 44,
                    stiff_id: 1,
                    nodes: vec![a, a + 1, a + 4, a + 3],
                });
            }
        }
        mesh.nodes.insert(10, DVec3::new(1., 1., -2.));
        mesh.elements.extend([
            ElementData {
                id: 5,
                elem_type: 10,
                stiff_id: 10,
                nodes: vec![4, 5],
            },
            ElementData {
                id: 6,
                elem_type: 10,
                stiff_id: 20,
                nodes: vec![5, 6],
            },
            ElementData {
                id: 7,
                elem_type: 10,
                stiff_id: 30,
                nodes: vec![10, 5],
            },
        ]);
        mesh
    }
    #[test]
    fn interior_joint_shares_one_vertex_without_splitting_axis_or_properties() {
        for scale in [0.1, 1., 10.] {
            for transformed in [false, true] {
                let mut mesh = slab_beam_column();
                let q = if transformed {
                    glam::DQuat::from_axis_angle(DVec3::new(1., 2., 3.).normalize(), 0.6)
                } else {
                    glam::DQuat::IDENTITY
                };
                let shift = DVec3::new(20., -30., 50.);
                let id = |n| if transformed { 100 - n } else { n };
                mesh.nodes = mesh
                    .nodes
                    .iter()
                    .map(|(&n, &p)| (id(n), q * (p * scale) + shift))
                    .collect();
                for e in &mut mesh.elements {
                    e.id = id(e.id);
                    for n in &mut e.nodes {
                        *n = id(*n);
                    }
                    e.nodes.reverse();
                }
                if transformed {
                    mesh.elements.reverse();
                }
                let f = frame(&mesh, q * DVec3::Z);
                let result = assembly::assemble(&mesh, &f, &policy(scale)).unwrap();
                let r = &result.axis_assembly;
                assert!(result.all_surface_patches_built);
                assert!(r.all_axes_built, "{:?}", r.issues);
                assert_eq!(r.axes.len(), 2);
                let beam = r.axes.iter().find(|a| a.spans.len() == 2).unwrap();
                assert_eq!(beam.endpoints.len(), 2);
                assert_eq!(
                    beam.spans
                        .iter()
                        .map(|s| s.stiffness)
                        .collect::<BTreeSet<_>>(),
                    BTreeSet::from([10, 20])
                );
                let middle = beam
                    .anchors
                    .iter()
                    .find(|a| a.source_node == id(5))
                    .unwrap();
                assert!((middle.t - 0.5).abs() < 1e-7);
                let column = r.axes.iter().find(|a| a.spans.len() == 1).unwrap();
                assert!(column.endpoints.contains(&middle.vertex));
                assert_eq!(
                    result
                        .vertex_source_nodes
                        .iter()
                        .filter(|&&n| n == id(5))
                        .count(),
                    1
                );
                assert!(r.contacts.iter().any(|c|matches!(c,Contact::Point {vertex,location:Location::Interior,..} if *vertex==middle.vertex)));
                assert!(r.contacts.iter().any(
                    |c| matches!(c,Contact::Interval {start_t,end_t,location:Location::Interior,..}
                if start_t.abs()<1e-7 && (*end_t-1.).abs()<1e-7)
                ));
                let sources: BTreeSet<_> = r
                    .axes
                    .iter()
                    .flat_map(|a| a.spans.iter().map(|s| s.element))
                    .collect();
                assert_eq!(sources, BTreeSet::from([id(5), id(6), id(7)]));
                for axis in &r.axes {
                    let a = DVec3::from_array(result.preview.vertices[axis.endpoints[0]]);
                    let b = DVec3::from_array(result.preview.vertices[axis.endpoints[1]]);
                    for anchor in &axis.anchors {
                        assert!(
                            DVec3::from_array(result.preview.vertices[anchor.vertex])
                                .distance(a.lerp(b, anchor.t))
                                < policy(scale).precision
                        );
                    }
                }
            }
        }
    }
    #[test]
    fn interior_fe_nodes_slide_on_axis_and_property_parameters_follow() {
        let mut mesh = MeshData::default();
        for y in 0..3 {
            for x in 0..5 {
                mesh.nodes
                    .insert(1 + x + 5 * y, DVec3::new(x as f64, y as f64, 0.));
            }
        }
        for y in 0..2 {
            for x in 0..4 {
                let a = 1 + x + 5 * y;
                mesh.elements.push(ElementData {
                    id: mesh.elements.len() as u32 + 1,
                    elem_type: 44,
                    stiff_id: 1,
                    nodes: vec![a, a + 1, a + 6, a + 5],
                });
            }
        }
        mesh.elements.extend([
            ElementData {
                id: 9,
                elem_type: 10,
                stiff_id: 10,
                nodes: vec![7, 8],
            },
            ElementData {
                id: 10,
                elem_type: 10,
                stiff_id: 20,
                nodes: vec![8, 9],
            },
        ]);
        let mut f = frame(&mesh, DVec3::Z);
        for (i, &n) in f.node_ids.iter().enumerate() {
            if (7..=9).contains(&n) {
                f.candidate_points[i][2] = 0.00002;
            }
            if n == 8 {
                f.candidate_points[i][0] += 0.02;
                f.candidate_points[i][1] += 0.002;
            }
        }
        let r = assembly::assemble(&mesh, &f, &policy(1.)).unwrap();
        assert!(
            r.axis_assembly.all_axes_built,
            "{:?}",
            r.axis_assembly.issues
        );
        assert_eq!(r.axis_assembly.axes.len(), 1);
        let axis = &r.axis_assembly.axes[0];
        let anchor = axis.anchors.iter().find(|a| a.source_node == 8).unwrap();
        let point = DVec3::from_array(r.preview.vertices[anchor.vertex]);
        assert!(point.distance(DVec3::new(2.02, 1., 0.)) < 1e-7);
        assert!((anchor.t - 0.51).abs() < 1e-7);
        assert!((axis.spans[0].end_t - anchor.t).abs() < 1e-7);
        assert!((axis.spans[1].start_t - anchor.t).abs() < 1e-7);
        assert_eq!(axis.spans[0].stiffness, 10);
        assert_eq!(axis.spans[1].stiffness, 20);
    }

    #[test]
    fn fixed_middle_joint_rejects_bent_beam_without_moving_surface_vertices() {
        let mesh = slab_beam_column();
        let mut f = frame(&mesh, DVec3::Z);
        let n = f.node_ids.iter().position(|&n| n == 5).unwrap();
        f.candidate_points[n][1] += 0.002;
        let before = f.candidate_points.clone();
        let r = assembly::assemble(&mesh, &f, &policy(1.)).unwrap();
        assert_eq!(f.candidate_points, before);
        assert!(r.all_surface_patches_built);
        assert!(!r.axis_assembly.all_axes_built);
        assert_eq!(r.axis_assembly.axes.len(), 1);
        assert_eq!(r.axis_assembly.issues.len(), 1);
        assert_eq!(
            r.axis_assembly.issues[0].reason,
            "shared_anchors_not_collinear"
        );
        assert_eq!(
            r.axis_assembly.issues[0]
                .source_elements
                .iter()
                .copied()
                .collect::<BTreeSet<_>>(),
            BTreeSet::from([5, 6])
        );
        assert_eq!(r.axis_assembly.issues[0].source_nodes, vec![5]);
        let vertex = r.vertex_source_nodes.iter().position(|&n| n == 5).unwrap();
        assert!(
            DVec3::from_array(r.preview.vertices[vertex]).distance(DVec3::new(1., 1.002, 0.))
                < 1e-7
        );
    }
    #[test]
    fn coplanar_intervals_respect_holes_and_boundary_roles() {
        let contours = vec![
            vec![[0., 0.], [4., 0.], [4., 4.], [0., 4.]],
            vec![[1., 1.], [3., 1.], [3., 3.], [1., 3.]],
        ];
        assert_eq!(
            intervals([0., 2.], [4., 2.], &contours, 1e-8),
            vec![
                (0., 0.25, Location::Interior),
                (0.75, 1., Location::Interior)
            ]
        );
        assert_eq!(
            intervals([0., 0.], [4., 0.], &contours, 1e-8),
            vec![(0., 1., Location::Boundary)]
        );
        assert_eq!(intervals([1.1, 2.], [2.9, 2.], &contours, 1e-8), vec![]);
        assert_eq!(location([2., 2.], &contours, 1e-8), None);
    }
    #[test]
    fn short_axis_preserves_vector_and_over_budget_attempt_creates_no_vertices() {
        let mut mesh = MeshData::default();
        mesh.nodes.insert(1, DVec3::ZERO);
        mesh.nodes.insert(2, DVec3::new(0.006, 0.008, 0.));
        mesh.elements.push(ElementData {
            id: 1,
            elem_type: 10,
            stiff_id: 9,
            nodes: vec![1, 2],
        });
        let mut f = frame(&mesh, DVec3::Z);
        for p in &mut f.candidate_points {
            p[2] += 0.0001;
        }
        let r = assembly::assemble(&mesh, &f, &policy(1.)).unwrap();
        assert!(
            r.axis_assembly.all_axes_built,
            "{:?}",
            r.axis_assembly.issues
        );
        let ends = r.axis_assembly.axes[0].endpoints;
        let d = DVec3::from_array(r.preview.vertices[ends[1]])
            - DVec3::from_array(r.preview.vertices[ends[0]]);
        assert!(d.distance(mesh.nodes[&2] - mesh.nodes[&1]) < 1e-8);
        for p in &mut f.candidate_points {
            p[2] += 0.01;
        }
        let refused = assembly::assemble(&mesh, &f, &policy(1.)).unwrap();
        assert!(!refused.axis_assembly.all_axes_built);
        assert!(refused.preview.vertices.is_empty());
        assert!(refused.axis_assembly.axes.is_empty());
        assert_eq!(
            refused.axis_assembly.issues[0].reason,
            "axis_movement_budget"
        );
    }
}
