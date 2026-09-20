//! Explicit surface-surface junction topology.
//!
//! The source FE mesh may contain panels that geometrically intersect but do
//! not share node/edge identities.  This module derives the finite
//! plane-intersection segments from assembled surface contours, materializes
//! their endpoints once, splits existing boundary edges when necessary and
//! registers the same model edge as a junction constraint on every owner.
//! Numerical coincidence alone is never used as an engineering weld: only a
//! verified finite surface-surface intersection creates topology. Re-running the
//! pass is idempotent: already materialized junction vertices and edges are reused.

use super::{EdgeUse, Model, PlaneFrame, Surface};
use glam::DVec3;
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Kind {
    BoundaryJunction,
    TJunction,
    Crossing,
}

#[derive(Debug, Clone, Serialize)]
pub struct Segment {
    pub surfaces: [usize; 2],
    pub kind: Kind,
    pub vertices: [usize; 2],
    pub length: f64,
}

#[derive(Debug, Clone, Serialize, Default)]
pub struct Report {
    pub candidate_pairs: usize,
    pub detected_segments: usize,
    pub generated_vertices: usize,
    pub split_boundary_edges: usize,
    pub shared_constraint_edges: usize,
    pub segments: Vec<Segment>,
}

#[derive(Debug, Clone)]
struct RawSegment {
    surfaces: [usize; 2],
    kind: Kind,
    start: DVec3,
    end: DVec3,
}

fn cross2(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}

fn sub2(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] - b[0], a[1] - b[1]]
}

fn point_on_segment_2d(p: [f64; 2], a: [f64; 2], b: [f64; 2], eps: f64) -> bool {
    let ab = sub2(b, a);
    let ap = sub2(p, a);
    let length = ab[0].hypot(ab[1]);
    if length <= eps {
        return ap[0].hypot(ap[1]) <= eps;
    }
    if cross2(ab, ap).abs() > eps * length {
        return false;
    }
    let dot = ap[0] * ab[0] + ap[1] * ab[1];
    dot >= -eps * length && dot <= length * length + eps * length
}

fn ring_location(ring: &[[f64; 2]], point: [f64; 2], eps: f64) -> (bool, bool) {
    if ring.len() < 3 {
        return (false, false);
    }
    let mut inside = false;
    for i in 0..ring.len() {
        let a = ring[i];
        let b = ring[(i + 1) % ring.len()];
        if point_on_segment_2d(point, a, b, eps) {
            return (true, true);
        }
        let crosses = (a[1] > point[1]) != (b[1] > point[1]);
        if crosses {
            let x = a[0] + (point[1] - a[1]) * (b[0] - a[0]) / (b[1] - a[1]);
            if x > point[0] {
                inside = !inside;
            }
        }
    }
    (false, inside)
}

fn contains(surface: &Surface, plane: &PlaneFrame, point: DVec3, eps: f64) -> bool {
    let uv = plane.project(point.to_array());
    let (outer_boundary, outer_inside) = ring_location(&surface.contours[0], uv, eps);
    if !outer_boundary && !outer_inside {
        return false;
    }
    for hole in surface.contours.iter().skip(1) {
        let (boundary, inside) = ring_location(hole, uv, eps);
        if boundary {
            return true;
        }
        if inside {
            return false;
        }
    }
    true
}

fn on_boundary(surface: &Surface, plane: &PlaneFrame, point: DVec3, eps: f64) -> bool {
    let uv = plane.project(point.to_array());
    surface
        .contours
        .iter()
        .any(|ring| ring_location(ring, uv, eps).0)
}

fn aabb(surface: &Surface, plane: &PlaneFrame) -> (DVec3, DVec3) {
    let mut low = DVec3::splat(f64::INFINITY);
    let mut high = DVec3::splat(f64::NEG_INFINITY);
    for point in surface.contours.iter().flatten() {
        let p = DVec3::from_array(plane.lift(*point));
        low = low.min(p);
        high = high.max(p);
    }
    (low, high)
}

fn line_of_planes(a: &PlaneFrame, b: &PlaneFrame) -> Option<(DVec3, DVec3)> {
    let na = DVec3::from_array(a.normal);
    let nb = DVec3::from_array(b.normal);
    let cross = na.cross(nb);
    let sine = cross.length();
    if !sine.is_finite() || sine < 1e-8 {
        return None;
    }
    let direction = cross / sine;
    let oa = DVec3::from_array(a.origin);
    let ob = DVec3::from_array(b.origin);
    let offset = (ob - oa).dot(nb);
    let origin = oa + direction.cross(na) * (offset / sine);
    Some((origin, direction))
}

fn dedup_parameters(values: &mut Vec<f64>, eps: f64) {
    values.sort_by(|a, b| a.total_cmp(b));
    values.dedup_by(|a, b| (*a - *b).abs() <= eps);
}

fn merge_intervals(mut values: Vec<(f64, f64)>, eps: f64) -> Vec<(f64, f64)> {
    values.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
    let mut result: Vec<(f64, f64)> = Vec::new();
    for (start, end) in values {
        if end - start <= eps {
            continue;
        }
        if let Some(last) = result.last_mut() {
            if start <= last.1 + eps {
                last.1 = last.1.max(end);
                continue;
            }
        }
        result.push((start, end));
    }
    result
}

fn surface_intervals(
    surface: &Surface,
    plane: &PlaneFrame,
    other: &PlaneFrame,
    origin: DVec3,
    direction: DVec3,
    eps: f64,
) -> Vec<(f64, f64)> {
    let mut parameters = Vec::new();
    for ring in &surface.contours {
        for i in 0..ring.len() {
            let p = DVec3::from_array(plane.lift(ring[i]));
            let q = DVec3::from_array(plane.lift(ring[(i + 1) % ring.len()]));
            let dp = other.distance(p.to_array());
            let dq = other.distance(q.to_array());
            if dp.abs() <= eps && dq.abs() <= eps {
                parameters.push((p - origin).dot(direction));
                parameters.push((q - origin).dot(direction));
                continue;
            }
            if dp > eps && dq > eps || dp < -eps && dq < -eps {
                continue;
            }
            let denominator = dp - dq;
            if denominator.abs() <= f64::EPSILON {
                continue;
            }
            let t = (dp / denominator).clamp(0.0, 1.0);
            let point = p.lerp(q, t);
            parameters.push((point - origin).dot(direction));
        }
    }
    dedup_parameters(&mut parameters, eps);
    let mut intervals = Vec::new();
    for pair in parameters.windows(2) {
        let (start, end) = (pair[0], pair[1]);
        if end - start <= eps {
            continue;
        }
        let midpoint = origin + direction * ((start + end) * 0.5);
        if contains(surface, plane, midpoint, eps) {
            intervals.push((start, end));
        }
    }
    merge_intervals(intervals, eps)
}

fn segment_parameter(point: DVec3, a: DVec3, b: DVec3, eps: f64) -> Option<f64> {
    let ab = b - a;
    let length2 = ab.length_squared();
    if length2 <= eps * eps {
        return None;
    }
    let t = (point - a).dot(ab) / length2;
    let tolerance = eps / length2.sqrt();
    if t < -tolerance || t > 1.0 + tolerance {
        return None;
    }
    let t = t.clamp(0.0, 1.0);
    (a.lerp(b, t).distance(point) <= eps).then_some(t)
}

fn segment_intersection_on_plane(
    a0: DVec3,
    a1: DVec3,
    b0: DVec3,
    b1: DVec3,
    plane: &PlaneFrame,
    eps: f64,
) -> Option<DVec3> {
    let p = plane.project(a0.to_array());
    let p1 = plane.project(a1.to_array());
    let q = plane.project(b0.to_array());
    let q1 = plane.project(b1.to_array());
    let r = sub2(p1, p);
    let s = sub2(q1, q);
    let denominator = cross2(r, s);
    let scale = (r[0].hypot(r[1]) * s[0].hypot(s[1])).max(eps);
    if denominator.abs() <= eps * scale {
        return None;
    }
    let qp = sub2(q, p);
    let t = cross2(qp, s) / denominator;
    let u = cross2(qp, r) / denominator;
    let tolerance = eps / r[0].hypot(r[1]).min(s[0].hypot(s[1])).max(eps);
    if t < -tolerance || t > 1.0 + tolerance || u < -tolerance || u > 1.0 + tolerance {
        return None;
    }
    Some(a0.lerp(a1, t.clamp(0.0, 1.0)))
}

fn push_unique_point(points: &mut Vec<DVec3>, point: DVec3, eps: f64) {
    if !points.iter().any(|p| p.distance(point) <= eps) {
        points.push(point);
    }
}

fn raw_segments(model: &Model, eps: f64) -> (usize, Vec<RawSegment>) {
    let bounds: Vec<_> = model
        .surfaces
        .iter()
        .map(|surface| aabb(surface, &model.planes[surface.plane]))
        .collect();
    let mut candidate_pairs = 0;
    let mut result = Vec::new();
    for i in 0..model.surfaces.len() {
        let a = &model.surfaces[i];
        let pa = &model.planes[a.plane];
        for j in i + 1..model.surfaces.len() {
            let b = &model.surfaces[j];
            let pb = &model.planes[b.plane];
            let (alo, ahi) = bounds[i];
            let (blo, bhi) = bounds[j];
            if ahi.x + eps < blo.x
                || bhi.x + eps < alo.x
                || ahi.y + eps < blo.y
                || bhi.y + eps < alo.y
                || ahi.z + eps < blo.z
                || bhi.z + eps < alo.z
            {
                continue;
            }
            let Some((origin, direction)) = line_of_planes(pa, pb) else {
                continue;
            };
            candidate_pairs += 1;
            let ia = surface_intervals(a, pa, pb, origin, direction, eps);
            let ib = surface_intervals(b, pb, pa, origin, direction, eps);
            for &(a0, a1) in &ia {
                for &(b0, b1) in &ib {
                    let start = a0.max(b0);
                    let end = a1.min(b1);
                    if end - start <= eps {
                        continue;
                    }
                    let p0 = origin + direction * start;
                    let p1 = origin + direction * end;
                    let midpoint = p0.lerp(p1, 0.5);
                    let boundary_a = on_boundary(a, pa, midpoint, eps);
                    let boundary_b = on_boundary(b, pb, midpoint, eps);
                    let kind = match (boundary_a, boundary_b) {
                        (true, true) => Kind::BoundaryJunction,
                        (true, false) | (false, true) => Kind::TJunction,
                        (false, false) => Kind::Crossing,
                    };
                    result.push(RawSegment {
                        surfaces: [i, j],
                        kind,
                        start: p0,
                        end: p1,
                    });
                }
            }
        }
    }
    (candidate_pairs, result)
}

fn ordered_ring_vertices(
    ring: &[EdgeUse],
    edges: &[[usize; 2]],
) -> Result<Vec<usize>, &'static str> {
    if ring.len() < 3 {
        return Err("invalid surface boundary");
    }
    let mut result = Vec::with_capacity(ring.len());
    let mut expected = None;
    for edge_use in ring {
        let [low, high] = *edges
            .get(edge_use.edge)
            .ok_or("invalid surface boundary edge")?;
        let (start, end) = if edge_use.reversed {
            (high, low)
        } else {
            (low, high)
        };
        if expected.is_some_and(|previous| previous != start) {
            return Err("disconnected surface boundary");
        }
        result.push(start);
        expected = Some(end);
    }
    if expected != result.first().copied() {
        return Err("open surface boundary");
    }
    Ok(result)
}

fn intern_edge(
    edges: &mut Vec<[usize; 2]>,
    index: &mut BTreeMap<[usize; 2], usize>,
    a: usize,
    b: usize,
) -> Result<usize, &'static str> {
    if a == b {
        return Err("degenerate junction edge");
    }
    let key = [a.min(b), a.max(b)];
    Ok(*index.entry(key).or_insert_with(|| {
        let id = edges.len();
        edges.push(key);
        id
    }))
}

fn vertex_for_point(model: &mut Model, point: DVec3, tolerance: f64) -> (usize, bool) {
    if let Some((index, _)) = model
        .vertices
        .iter()
        .enumerate()
        .filter_map(|(i, p)| {
            let distance = DVec3::from_array(*p).distance(point);
            (distance <= tolerance).then_some((i, distance))
        })
        .min_by(|a, b| a.1.total_cmp(&b.1))
    {
        return (index, false);
    }
    let index = model
        .add_vertex(point.to_array())
        .expect("finite junction intersection");
    (index, true)
}

fn rebuild_boundaries(
    model: &mut Model,
    split_vertices: &[usize],
    eps: f64,
) -> Result<usize, &'static str> {
    let old_edges = model.edges.clone();
    let original_edge_count: usize = model
        .surfaces
        .iter()
        .map(|surface| surface.boundaries.iter().map(Vec::len).sum::<usize>())
        .sum();
    let old_rings: Vec<Vec<Vec<usize>>> = model
        .surfaces
        .iter()
        .map(|surface| {
            surface
                .boundaries
                .iter()
                .map(|ring| ordered_ring_vertices(ring, &old_edges))
                .collect::<Result<Vec<_>, _>>()
        })
        .collect::<Result<_, _>>()?;

    let mut expanded = Vec::with_capacity(old_rings.len());
    for rings in old_rings {
        let mut surface_rings = Vec::with_capacity(rings.len());
        for ring in rings {
            let mut output = Vec::new();
            for i in 0..ring.len() {
                let a = ring[i];
                let b = ring[(i + 1) % ring.len()];
                let pa = DVec3::from_array(model.vertices[a]);
                let pb = DVec3::from_array(model.vertices[b]);
                let mut interior = Vec::new();
                for &vertex in split_vertices {
                    if vertex == a || vertex == b {
                        continue;
                    }
                    let p = DVec3::from_array(model.vertices[vertex]);
                    if let Some(t) = segment_parameter(p, pa, pb, eps) {
                        let distance_a = pa.distance(p);
                        let distance_b = pb.distance(p);
                        if distance_a >= model.minimum_edge && distance_b >= model.minimum_edge {
                            interior.push((t, vertex));
                        }
                    }
                }
                interior.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
                interior.dedup_by_key(|item| item.1);
                output.push(a);
                output.extend(interior.into_iter().map(|(_, vertex)| vertex));
            }
            surface_rings.push(output);
        }
        expanded.push(surface_rings);
    }

    model.edges.clear();
    model.edge_index.clear();
    for surface in &mut model.surfaces {
        surface.boundaries.clear();
        surface.junctions.clear();
    }
    for (surface_index, rings) in expanded.iter().enumerate() {
        let plane = model.surfaces[surface_index].plane;
        let mut boundaries = Vec::with_capacity(rings.len());
        let mut contours = Vec::with_capacity(rings.len());
        for ring in rings {
            let mut uses = Vec::with_capacity(ring.len());
            let mut contour = Vec::with_capacity(ring.len());
            for i in 0..ring.len() {
                let a = ring[i];
                let b = ring[(i + 1) % ring.len()];
                let edge = intern_edge(&mut model.edges, &mut model.edge_index, a, b)?;
                uses.push(EdgeUse {
                    edge,
                    reversed: a > b,
                });
                contour.push(model.planes[plane].project(model.vertices[a]));
            }
            boundaries.push(uses);
            contours.push(contour);
        }
        model.surfaces[surface_index].boundaries = boundaries;
        model.surfaces[surface_index].contours = contours;
    }
    let new_edge_count: usize = model
        .surfaces
        .iter()
        .map(|surface| surface.boundaries.iter().map(Vec::len).sum::<usize>())
        .sum();
    Ok(new_edge_count.saturating_sub(original_edge_count))
}

/// Materialize all finite non-parallel surface intersections as shared topology.
///
/// Parallel/coplanar overlap remains an audit responsibility because deciding
/// whether two coincident panels should merge is semantic, not purely geometric.
pub fn conform(model: &mut Model, precision: f64) -> Result<Report, &'static str> {
    if !precision.is_finite() || precision <= 0.0 {
        return Err("invalid junction precision");
    }
    let (candidate_pairs, raw) = raw_segments(model, precision);
    if raw.is_empty() {
        return Ok(Report {
            candidate_pairs,
            ..Report::default()
        });
    }

    let mut split_points = Vec::new();
    for segment in &raw {
        push_unique_point(&mut split_points, segment.start, precision);
        push_unique_point(&mut split_points, segment.end, precision);
    }
    // Triple/multi-surface junctions can cross in the interior of both pairwise
    // segments. Materialize that shared point before any constraint is added.
    for i in 0..raw.len() {
        for j in i + 1..raw.len() {
            let shared_surface = raw[i]
                .surfaces
                .iter()
                .copied()
                .find(|surface| raw[j].surfaces.contains(surface));
            let Some(surface) = shared_surface else {
                continue;
            };
            let plane = &model.planes[model.surfaces[surface].plane];
            if let Some(point) = segment_intersection_on_plane(
                raw[i].start,
                raw[i].end,
                raw[j].start,
                raw[j].end,
                plane,
                precision,
            ) {
                push_unique_point(&mut split_points, point, precision);
            }
        }
    }

    let old_vertex_count = model.vertices.len();
    let reuse_tolerance = precision.max(model.minimum_edge);
    let mut split_vertices = Vec::with_capacity(split_points.len());
    for point in split_points {
        let (vertex, _) = vertex_for_point(model, point, reuse_tolerance);
        if !split_vertices.contains(&vertex) {
            split_vertices.push(vertex);
        }
    }

    let split_boundary_edges = rebuild_boundaries(model, &split_vertices, precision)?;

    let mut shared_edges = BTreeSet::new();
    let mut segments = Vec::new();
    for item in raw {
        let start = split_vertices
            .iter()
            .copied()
            .filter_map(|vertex| {
                let p = DVec3::from_array(model.vertices[vertex]);
                let distance = p.distance(item.start);
                (distance <= reuse_tolerance).then_some((vertex, distance))
            })
            .min_by(|a, b| a.1.total_cmp(&b.1))
            .map(|x| x.0)
            .ok_or("missing junction start vertex")?;
        let end = split_vertices
            .iter()
            .copied()
            .filter_map(|vertex| {
                let p = DVec3::from_array(model.vertices[vertex]);
                let distance = p.distance(item.end);
                (distance <= reuse_tolerance).then_some((vertex, distance))
            })
            .min_by(|a, b| a.1.total_cmp(&b.1))
            .map(|x| x.0)
            .ok_or("missing junction end vertex")?;
        let pa = DVec3::from_array(model.vertices[start]);
        let pb = DVec3::from_array(model.vertices[end]);
        if pa.distance(pb) <= model.minimum_edge {
            continue;
        }

        let mut chain = Vec::new();
        for &vertex in &split_vertices {
            let point = DVec3::from_array(model.vertices[vertex]);
            if let Some(t) = segment_parameter(point, pa, pb, precision) {
                chain.push((t, vertex));
            }
        }
        chain.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
        chain.dedup_by_key(|item| item.1);
        for pair in chain.windows(2) {
            let a = pair[0].1;
            let b = pair[1].1;
            if DVec3::from_array(model.vertices[a])
                .distance(DVec3::from_array(model.vertices[b]))
                < model.minimum_edge
            {
                continue;
            }
            let edge = intern_edge(&mut model.edges, &mut model.edge_index, a, b)?;
            for &surface in &item.surfaces {
                model.surfaces[surface].junctions.push(edge);
            }
            shared_edges.insert(edge);
        }
        for &surface in &item.surfaces {
            model.surfaces[surface].junctions.sort_unstable();
            model.surfaces[surface].junctions.dedup();
        }
        segments.push(Segment {
            surfaces: item.surfaces,
            kind: item.kind,
            vertices: [start, end],
            length: pa.distance(pb),
        });
    }

    Ok(Report {
        candidate_pairs,
        detected_segments: segments.len(),
        generated_vertices: model.vertices.len() - old_vertex_count,
        split_boundary_edges,
        shared_constraint_edges: shared_edges.len(),
        segments,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reconstruction::PlaneFrame;

    fn add_ring(model: &mut Model, points: &[[f64; 3]]) -> Vec<usize> {
        points
            .iter()
            .map(|point| model.add_vertex(*point).unwrap())
            .collect()
    }

    fn fixture(shift: DVec3, reverse_wall: bool) -> Model {
        let mut model = Model::new(1e-8, 1e-6).unwrap();
        let slab_plane = model.add_plane(
            PlaneFrame::new(
                shift.to_array(),
                [0.0, 0.0, 1.0],
            )
            .unwrap(),
        );
        let slab = [
            [-2.0, -2.0, 0.0],
            [2.0, -2.0, 0.0],
            [2.0, 2.0, 0.0],
            [-2.0, 2.0, 0.0],
        ]
        .map(|p| (DVec3::from_array(p) + shift).to_array());
        let slab_ring = add_ring(&mut model, &slab);
        model.add_surface(slab_plane, vec![slab_ring], vec![1]).unwrap();

        let wall_normal = if reverse_wall {
            [-1.0, 0.0, 0.0]
        } else {
            [1.0, 0.0, 0.0]
        };
        let wall_plane = model.add_plane(PlaneFrame::new(shift.to_array(), wall_normal).unwrap());
        let wall = [
            [0.0, -1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 1.0, 2.0],
            [0.0, -1.0, 2.0],
        ]
        .map(|p| (DVec3::from_array(p) + shift).to_array());
        let wall_ring = add_ring(&mut model, &wall);
        model.add_surface(wall_plane, vec![wall_ring], vec![2]).unwrap();
        model
    }

    #[test]
    fn t_junction_becomes_one_shared_constraint() {
        let mut model = fixture(DVec3::ZERO, false);
        let report = conform(&mut model, 1e-8).unwrap();
        assert_eq!(report.detected_segments, 1);
        assert_eq!(report.segments[0].kind, Kind::TJunction);
        assert!(report.shared_constraint_edges >= 1);
        let common: BTreeSet<_> = model.surfaces[0]
            .junctions
            .iter()
            .copied()
            .collect::<BTreeSet<_>>()
            .intersection(&model.surfaces[1].junctions.iter().copied().collect())
            .copied()
            .collect();
        assert!(!common.is_empty());
        for edge in common {
            let [a, b] = model.edges[edge];
            assert!((model.vertices[a][0]).abs() < 1e-8);
            assert!((model.vertices[b][0]).abs() < 1e-8);
            assert!((model.vertices[a][2]).abs() < 1e-8);
            assert!((model.vertices[b][2]).abs() < 1e-8);
        }
    }

    #[test]
    fn result_is_transform_invariant_and_idempotent() {
        for (shift, reverse) in [
            (DVec3::ZERO, false),
            (DVec3::new(1.2e6, -8.0e5, 42.0), true),
        ] {
            let mut model = fixture(shift, reverse);
            let first = conform(&mut model, 1e-8).unwrap();
            assert_eq!(first.detected_segments, 1);
            let vertices = model.vertices.clone();
            let edges = model.edges.clone();
            let boundaries: Vec<_> = model
                .surfaces
                .iter()
                .map(|s| (s.boundaries.clone(), s.junctions.clone()))
                .collect();
            let second = conform(&mut model, 1e-8).unwrap();
            assert_eq!(second.generated_vertices, 0);
            assert_eq!(model.vertices, vertices);
            assert_eq!(model.edges, edges);
            assert_eq!(
                model
                    .surfaces
                    .iter()
                    .map(|s| (s.boundaries.clone(), s.junctions.clone()))
                    .collect::<Vec<_>>(),
                boundaries
            );
        }
    }

    #[test]
    fn crossing_splits_boundary_edges_and_keeps_shared_ids() {
        let mut model = Model::new(1e-8, 1e-6).unwrap();
        let slab_plane = model.add_plane(PlaneFrame::new([0., 0., 0.], [0., 0., 1.]).unwrap());
        let slab_ring = add_ring(
            &mut model,
            &[
                [-2., -1., 0.],
                [2., -1., 0.],
                [2., 1., 0.],
                [-2., 1., 0.],
            ],
        );
        model.add_surface(slab_plane, vec![slab_ring], vec![1]).unwrap();

        let wall_plane = model.add_plane(PlaneFrame::new([0., 0., 0.], [1., 0., 0.]).unwrap());
        let wall_ring = add_ring(
            &mut model,
            &[
                [0., -2., -2.],
                [0., 2., -2.],
                [0., 2., 2.],
                [0., -2., 2.],
            ],
        );
        model.add_surface(wall_plane, vec![wall_ring], vec![2]).unwrap();

        let report = conform(&mut model, 1e-8).unwrap();
        assert_eq!(report.detected_segments, 1);
        assert_eq!(report.segments[0].kind, Kind::Crossing);
        assert_eq!(report.generated_vertices, 2);
        assert_eq!(report.split_boundary_edges, 2);
        let common = model.surfaces[0]
            .junctions
            .iter()
            .filter(|edge| model.surfaces[1].junctions.contains(edge))
            .count();
        assert!(common >= 1);
    }

    #[test]
    fn nearby_nonintersecting_surfaces_remain_separate() {
        let mut model = Model::new(1e-8, 1e-6).unwrap();
        let p0 = model.add_plane(PlaneFrame::new([0., 0., 0.], [0., 0., 1.]).unwrap());
        let a = add_ring(
            &mut model,
            &[[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.]],
        );
        model.add_surface(p0, vec![a], vec![1]).unwrap();
        let p1 = model.add_plane(PlaneFrame::new([0., 0., 0.01], [0., 0., 1.]).unwrap());
        let b = add_ring(
            &mut model,
            &[
                [0., 0., 0.01],
                [1., 0., 0.01],
                [1., 1., 0.01],
                [0., 1., 0.01],
            ],
        );
        model.add_surface(p1, vec![b], vec![2]).unwrap();
        let report = conform(&mut model, 1e-8).unwrap();
        assert_eq!(report.detected_segments, 0);
        assert!(model.surfaces.iter().all(|s| s.junctions.is_empty()));
    }
}
