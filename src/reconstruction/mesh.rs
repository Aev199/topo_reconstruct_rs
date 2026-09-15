//! Conforming mesh probe for a completely assembled geometric fragment.
//! Does not repair geometry, invent connections, or transfer loads/releases.
use super::assembly::{self, bars::Contact};
use geo::{Area, LineString, Polygon};
use glam::DVec3;
use serde::Serialize;
use spade::{
    AngleLimit, ConstrainedDelaunayTriangulation, Point2, RefinementParameters, Triangulation,
};
use std::collections::{BTreeMap, BTreeSet};
mod domain;

#[derive(Clone, Debug, Serialize)]
pub struct Policy {
    pub boundary_spacing: f64,
    pub maximum_area: f64,
    pub minimum_angle_degrees: f64,
    pub maximum_added_vertices_per_surface: usize,
}
#[derive(Debug, Serialize)]
pub struct Triangle {
    pub vertices: [usize; 3],
    pub surface: usize,
    pub stiffness: u32,
}
#[derive(Debug, Serialize)]
pub struct Bar {
    pub vertices: [usize; 2],
    pub axis: usize,
    pub source_element: u32,
    pub stiffness: u32,
}
#[derive(Debug, Serialize)]
pub struct Report {
    pub policy: Policy,
    pub vertices: Vec<[f64; 3]>,
    pub triangles: Vec<Triangle>,
    pub bars: Vec<Bar>,
    /// Region provenance, not a claim that each new triangle is one old FE.
    pub surface_source_elements: Vec<Vec<u32>>,
    pub minimum_angle_degrees: f64,
    pub maximum_triangle_area: f64,
    /// Maximum ratio between the longest and shortest edge of a non-degenerate
    /// emitted triangle. This is diagnostic only; the acceptance threshold is
    /// intentionally profile-specific and is not hard-coded here.
    pub maximum_edge_ratio: f64,
    pub topology_valid: bool,
    pub quality_passed: bool,
    pub blockers: Vec<String>,
    pub export_ready: bool,
}
fn key(a: usize, b: usize) -> [usize; 2] {
    if a < b {
        [a, b]
    } else {
        [b, a]
    }
}
fn point(v: &[f64; 3]) -> DVec3 {
    DVec3::from_array(*v)
}
fn parameter(p: DVec3, a: DVec3, b: DVec3, eps: f64) -> Option<f64> {
    let d = b - a;
    if d.length_squared() == 0. {
        return None;
    }
    let t = (p - a).dot(d) / d.length_squared();
    (t >= 0. && t <= 1. && p.distance(a + d * t) <= eps).then_some(t)
}
fn subdivide(
    chain: &[(f64, usize)],
    vertices: &mut Vec<[f64; 3]>,
    spacing: f64,
    limit: usize,
    precision: f64,
) -> Result<Vec<(f64, usize)>, &'static str> {
    let mut out = vec![chain[0]];
    for pair in chain.windows(2) {
        let [(ta, a), (tb, b)] = [pair[0], pair[1]];
        let pa = point(&vertices[a]);
        let pb = point(&vertices[b]);
        // A rigid transform must not turn an exact two-way split into three
        // pieces merely because the computed length is a few ulps longer.
        let count = ((pa.distance(pb) - precision) / spacing).ceil().max(1.);
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
fn sorted(
    chain: &mut Vec<(f64, usize)>,
    vertices: &[[f64; 3]],
    eps: f64,
) -> Result<(), &'static str> {
    chain.sort_by(|a, b| a.0.total_cmp(&b.0));
    chain.dedup_by(|a, b| a.1 == b.1);
    if chain
        .windows(2)
        .any(|w| point(&vertices[w[0].1]).distance(point(&vertices[w[1].1])) <= eps)
    {
        return Err("distinct constraint vertices coincide");
    }
    Ok(())
}

pub fn build(source: &assembly::Report, policy: &Policy) -> Result<Report, &'static str> {
    if !source.all_surface_patches_built
        || !source.issues.is_empty()
        || !source.axis_assembly.all_axes_built
    {
        return Err("mesh probe requires a completely assembled fragment");
    }
    if !policy.boundary_spacing.is_finite()
        || policy.boundary_spacing <= 0.
        || !policy.maximum_area.is_finite()
        || policy.maximum_area <= 0.
        || !policy.minimum_angle_degrees.is_finite()
        || policy.minimum_angle_degrees <= 0.
        || policy.minimum_angle_degrees > 30.
        || policy.maximum_added_vertices_per_surface == 0
    {
        return Err("invalid mesh policy");
    }
    let model = &source.preview;
    let eps = source.policy.precision;
    let mut vertices = model.vertices.clone();
    let axes = &source.axis_assembly.axes;
    let contacts = &source.axis_assembly.contacts;
    let mut axis_nodes = Vec::new();
    for axis in axes {
        let mut chain: Vec<_> = axis.anchors.iter().map(|a| (a.t, a.vertex)).collect();
        sorted(&mut chain, &vertices, eps)?;
        if chain.len() < 2 {
            return Err("empty axis");
        }
        axis_nodes.push(subdivide(
            &chain,
            &mut vertices,
            policy.boundary_spacing,
            policy.maximum_added_vertices_per_surface,
            eps,
        )?);
    }
    let mut owners = vec![BTreeSet::new(); model.edges.len()];
    for (s, surface) in model.surfaces.iter().enumerate() {
        for edge in surface.boundaries.iter().flatten() {
            owners[edge.edge].insert(s);
        }
    }
    let mut edge_nodes = Vec::new();
    for (e, &[a, b]) in model.edges.iter().enumerate() {
        let mut chain = vec![(0., a), (1., b)];
        let pa = point(&vertices[a]);
        let pb = point(&vertices[b]);
        for c in contacts {
            let (surface, candidates) = match *c {
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
                        .filter(|(t, _)| *t >= start_t && *t <= end_t)
                        .map(|(_, n)| *n)
                        .collect(),
                ),
            };
            if !owners[e].contains(&surface) {
                continue;
            }
            for node in candidates {
                if let Some(t) = parameter(point(&vertices[node]), pa, pb, eps) {
                    chain.push((t, node));
                }
            }
        }
        sorted(&mut chain, &vertices, eps)?;
        edge_nodes.push(subdivide(
            &chain,
            &mut vertices,
            policy.boundary_spacing,
            policy.maximum_added_vertices_per_surface,
            eps,
        )?);
    }
    // Shared boundary subdivisions also become bar nodes where incidence is known.
    for c in contacts {
        if let Contact::Interval {
            axis,
            surface,
            start_t,
            end_t,
            ..
        } = *c
        {
            let [a, b] = axes[axis].endpoints;
            for edge in model.surfaces[surface].boundaries.iter().flatten() {
                for &(_, node) in &edge_nodes[edge.edge] {
                    if let Some(t) = parameter(
                        point(&vertices[node]),
                        point(&vertices[a]),
                        point(&vertices[b]),
                        eps,
                    ) {
                        if t >= start_t && t <= end_t {
                            axis_nodes[axis].push((t, node));
                        }
                    }
                }
            }
        }
    }
    for chain in &mut axis_nodes {
        sorted(chain, &vertices, eps)?;
    }
    let mut triangles = Vec::new();
    let mut blockers = Vec::new();
    let mut topology_valid = true;
    let mut minimum_angle = 180.0_f64;
    let mut max_area = 0.0_f64;
    let mut maximum_edge_ratio = 0.0_f64;
    for (s, surface) in model.surfaces.iter().enumerate() {
        let support = &model.planes[surface.plane];
        // A boundary-aligned, local numerical chart avoids loss of significant
        // digits under rigid transforms. Rounding affects only CDT predicates;
        // shared 3D coordinates stay unchanged and are checked after meshing.
        let edge = surface
            .boundaries
            .iter()
            .flatten()
            .max_by(|a, b| {
                let length = |e: usize| {
                    let [i, j] = model.edges[e];
                    point(&vertices[i]).distance(point(&vertices[j]))
                };
                length(a.edge).total_cmp(&length(b.edge))
            })
            .ok_or("empty surface")?;
        let [a, b] = model.edges[edge.edge];
        let origin = point(&vertices[a]);
        let normal = point(&support.normal);
        let direction = point(&vertices[b]) - origin;
        let u = (direction - normal * normal.dot(direction)).normalize();
        let v = normal.cross(u);
        let plane = super::PlaneFrame {
            origin: origin.to_array(),
            normal: support.normal,
            u: u.to_array(),
            v: v.to_array(),
        };
        let project = |p| {
            let uv = plane.project(p);
            uv.map(|x| (x / eps).round() * eps)
        };
        let ring = |r: &Vec<[f64; 2]>| {
            LineString::from(
                r.iter()
                    .chain(r.first())
                    .map(|p| (p[0], p[1]))
                    .collect::<Vec<_>>(),
            )
        };
        let contours: Vec<Vec<_>> = surface
            .contours
            .iter()
            .map(|r| r.iter().map(|uv| project(support.lift(*uv))).collect())
            .collect();
        let polygon = Polygon::new(
            ring(&contours[0]),
            contours.iter().skip(1).map(ring).collect(),
        );
        let mut constraints = BTreeSet::new();
        let mut boundary = BTreeSet::new();
        let mut barriers = BTreeSet::new();
        let mut has_open_internal_constraint = false;
        let mut internal_constraint_edges = BTreeSet::new();
        let mut nodes = BTreeSet::new();
        for edge in surface.boundaries.iter().flatten() {
            for pair in edge_nodes[edge.edge].windows(2) {
                boundary.insert(key(pair[0].1, pair[1].1));
            }
        }
        barriers.extend(&boundary);
        let boundary_nodes: BTreeSet<_> = boundary.iter().flatten().copied().collect();
        constraints.extend(&boundary);
        for c in contacts {
            match *c {
                Contact::Point {
                    surface, vertex, ..
                } if surface == s => {
                    nodes.insert(vertex);
                }
                Contact::Interval {
                    surface,
                    axis,
                    start_t,
                    end_t,
                    ..
                } if surface == s => {
                    let length = point(&vertices[axes[axis].endpoints[0]])
                        .distance(point(&vertices[axes[axis].endpoints[1]]));
                    let tolerance = eps / length;
                    let chain = &axis_nodes[axis];
                    if ![start_t, end_t]
                        .iter()
                        .all(|t| chain.iter().any(|(u, _)| (t - u).abs() <= tolerance))
                    {
                        return Err("contact endpoint requires explicit intersection vertex");
                    }
                    let mut interval_edges = Vec::new();
                    for pair in chain.windows(2) {
                        if pair[0].0 >= start_t - tolerance && pair[1].0 <= end_t + tolerance {
                            let edge = key(pair[0].1, pair[1].1);
                            constraints.insert(edge);
                            interval_edges.push(edge);
                        }
                    }
                    // Segmenting a boundary-to-boundary construction can
                    // make each individual interval look open. Defer the
                    // barrier decision until all explicitly connected
                    // interval edges on this surface are known.
                    internal_constraint_edges.extend(interval_edges);
                }
                _ => {}
            }
        }
        // An internal constraint edge is a material-domain barrier exactly
        // when it belongs to a path between two distinct source boundary
        // vertices. This preserves a continuous boundary-to-boundary chain
        // across constructive segments while leaving true dangling branches
        // in one refinement domain.
        let mut graph = BTreeMap::<usize, BTreeSet<usize>>::new();
        for &[a, b] in &internal_constraint_edges {
            graph.entry(a).or_default().insert(b);
            graph.entry(b).or_default().insert(a);
        }
        let reachable_without = |start: usize, blocked: [usize; 2]| {
            let mut seen = BTreeSet::from([start]);
            let mut pending = vec![start];
            while let Some(node) = pending.pop() {
                for &next in graph.get(&node).into_iter().flatten() {
                    if key(node, next) == blocked || !seen.insert(next) {
                        continue;
                    }
                    pending.push(next);
                }
            }
            seen
        };
        for &edge in &internal_constraint_edges {
            let left = reachable_without(edge[0], edge);
            let right = reachable_without(edge[1], edge);
            let left_boundary: BTreeSet<_> = left.intersection(&boundary_nodes).copied().collect();
            let right_boundary: BTreeSet<_> =
                right.intersection(&boundary_nodes).copied().collect();
            if left_boundary
                .iter()
                .any(|a| right_boundary.iter().any(|b| a != b))
            {
                barriers.insert(edge);
            } else {
                has_open_internal_constraint = true;
            }
        }
        nodes.extend(constraints.iter().flatten().copied());
        let mut cdt = ConstrainedDelaunayTriangulation::<Point2<f64>>::new();
        let mut handles = BTreeMap::new();
        let mut global = BTreeMap::new();
        for n in nodes {
            if plane.distance(vertices[n]).abs() > eps {
                return Err("nonplanar mesh constraint");
            }
            let uv = project(vertices[n]);
            let h = cdt
                .insert(Point2::new(uv[0], uv[1]))
                .map_err(|_| "invalid CDT vertex")?;
            if global.insert(h.index(), n).is_some() {
                return Err("implicit vertex merge in CDT");
            }
            handles.insert(n, h);
        }
        for &[a, b] in &constraints {
            if !cdt.can_add_constraint(handles[&a], handles[&b]) {
                return Err("constraint crossing requires explicit shared vertex");
            }
            cdt.add_constraint(handles[&a], handles[&b]);
        }
        let (mut faces, complete) = domain::refine(
            &cdt,
            &global,
            &boundary,
            &barriers,
            &constraints,
            has_open_internal_constraint,
            &plane,
            &mut vertices,
            policy,
        )?;
        if !complete {
            blockers.push(format!("refinement_limit: surface={s}"));
        }
        // Refinement may split a constrained bar edge by inserting Steiner
        // vertices. Promote those vertices to parametric axis anchors before
        // emitting bars, and canonicalize equal generated points across
        // adjacent surfaces. Otherwise the surface mesh and the bar mesh
        // would describe different subdivisions of the same construction.
        let face_nodes: BTreeSet<_> = faces.iter().flatten().copied().collect();
        let mut remap = BTreeMap::new();
        for c in contacts {
            let Contact::Interval {
                axis,
                surface,
                start_t,
                end_t,
                ..
            } = *c
            else {
                continue;
            };
            if surface != s {
                continue;
            }
            let [a, b] = axes[axis].endpoints;
            let length = point(&vertices[a]).distance(point(&vertices[b]));
            let tolerance = eps / length;
            for node in face_nodes.iter().copied() {
                let Some(t) = parameter(
                    point(&vertices[node]),
                    point(&vertices[a]),
                    point(&vertices[b]),
                    eps,
                ) else {
                    continue;
                };
                if t < start_t - tolerance || t > end_t + tolerance {
                    continue;
                }
                let canonical = axis_nodes[axis]
                    .iter()
                    .find(|(u, _)| (t - *u).abs() <= tolerance)
                    .map(|(_, n)| *n)
                    .unwrap_or_else(|| {
                        axis_nodes[axis].push((t, node));
                        node
                    });
                if canonical != node {
                    remap.insert(node, canonical);
                }
            }
        }
        for triangle in &mut faces {
            for node in triangle {
                if let Some(&canonical) = remap.get(node) {
                    *node = canonical;
                }
            }
        }
        for chain in &mut axis_nodes {
            sorted(chain, &vertices, eps)?;
        }
        let mut counts = BTreeMap::new();
        let mut area = 0.;
        let mut used = BTreeSet::new();
        for ids in faces {
            let p = ids.map(|n| point(&vertices[n]));
            let normal = (p[1] - p[0]).cross(p[2] - p[0]);
            let ar = normal.length() / 2.;
            if !ar.is_finite() || ar <= eps * eps || normal.dot(point(&plane.normal)) <= 0. {
                minimum_angle = 0.;
                blockers.push(format!("degenerate_or_inverted_triangle: surface={s}"));
            }
            let lengths = [
                p[0].distance(p[1]),
                p[1].distance(p[2]),
                p[2].distance(p[0]),
            ];
            let shortest = lengths.iter().copied().fold(f64::INFINITY, f64::min);
            let longest = lengths.iter().copied().fold(0.0_f64, f64::max);
            if shortest.is_finite() && shortest > eps {
                maximum_edge_ratio = maximum_edge_ratio.max(longest / shortest);
            } else {
                blockers.push(format!("short_triangle_edge: surface={s}"));
            }
            area += ar;
            max_area = max_area.max(ar);
            for i in 0..3 {
                let angle = (p[(i + 1) % 3] - p[i])
                    .normalize()
                    .dot((p[(i + 2) % 3] - p[i]).normalize())
                    .clamp(-1., 1.)
                    .acos()
                    .to_degrees();
                minimum_angle = minimum_angle.min(angle);
                *counts
                    .entry(key(ids[i], ids[(i + 1) % 3]))
                    .or_insert(0usize) += 1;
            }
            used.extend(ids);
            triangles.push(Triangle {
                vertices: ids,
                surface: s,
                stiffness: source.surface_stiffness[s],
            });
        }
        if (area - polygon.unsigned_area()).abs()
            > eps * polygon.unsigned_area().sqrt().max(eps) * 10.
            || constraints.iter().any(|edge| {
                counts.get(edge).copied().unwrap_or(0)
                    != if boundary.contains(edge) { 1 } else { 2 }
            })
            || counts
                .iter()
                .any(|(edge, &count)| count != if boundary.contains(edge) { 1 } else { 2 })
            || handles.keys().any(|n| !used.contains(n))
        {
            topology_valid = false;
            blockers.push(format!("coverage_or_constraint_failure: surface={s}"));
        }
    }
    let mut bars = Vec::new();
    for (axis, chain) in axis_nodes.iter().enumerate() {
        for pair in chain.windows(2) {
            let mid = (pair[0].0 + pair[1].0) / 2.;
            let spans: Vec<_> = axes[axis]
                .spans
                .iter()
                .filter(|s| mid >= s.start_t && mid < s.end_t)
                .collect();
            if spans.len() != 1 {
                return Err("ambiguous bar property interval");
            }
            let span = spans[0];
            if pair[0].0 < span.start_t - 1e-10 || pair[1].0 > span.end_t + 1e-10 {
                return Err("bar crosses property boundary");
            }
            bars.push(Bar {
                vertices: [pair[0].1, pair[1].1],
                axis,
                source_element: span.element,
                stiffness: span.stiffness,
            });
        }
    }
    let quality_passed = !triangles.is_empty()
        && minimum_angle + 1e-7 >= policy.minimum_angle_degrees
        && max_area <= policy.maximum_area * (1. + 1e-7)
        && blockers.is_empty();
    if minimum_angle + 1e-7 < policy.minimum_angle_degrees {
        blockers.push("minimum_angle_not_met".into());
    }
    if max_area > policy.maximum_area * (1. + 1e-7) {
        blockers.push("maximum_area_not_met".into());
    }
    Ok(Report {
        policy: policy.clone(),
        vertices,
        triangles,
        bars,
        surface_source_elements: model
            .surfaces
            .iter()
            .map(|s| s.source_elements.clone())
            .collect(),
        minimum_angle_degrees: minimum_angle,
        maximum_triangle_area: max_area,
        maximum_edge_ratio,
        topology_valid,
        quality_passed,
        blockers,
        export_ready: false,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn subdivision_respects_numerical_precision_and_resource_limit() {
        let eps = 1e-7;
        for (excess, expected) in [(eps * 0.5, 3), (eps * 2., 4)] {
            let mut vertices = vec![[0., 0., 0.], [1. + excess, 0., 0.]];
            let chain = subdivide(&[(0., 0), (1., 1)], &mut vertices, 0.5, 10, eps).unwrap();
            assert_eq!(chain.len(), expected);
            assert_eq!(vertices[1], [1. + excess, 0., 0.]);
        }
        let mut vertices = vec![[0., 0., 0.], [1., 0., 0.]];
        assert!(subdivide(&[(0., 0), (1., 1)], &mut vertices, 0.5, 1, eps).is_err());
    }
}
