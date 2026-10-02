//! Surface topology preview. Engineering closure tolerance is distinct from
//! numerical planarity. Mechanical ties and mesh readiness are not inferred.
pub mod bars;
pub mod cleanup;
pub mod consoles;
pub mod cracks;
mod features;
pub mod gaps;
mod holes;
pub mod junctions;
pub mod stacking;
use super::{frame, planes, Model, PlaneFrame};
use crate::input::MeshData;
pub use features::{FeaturePolicy, SimplifiedHole};
use glam::DVec3;
pub use holes::{HoleConstraint, HoleNodeChange, HoleOutcome, HoleRecovery};
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Serialize)]
pub struct Policy {
    /// Maximum normal deviation when reconciling support planes.
    pub closure_tolerance: f64,
    /// Additional displacement allowed when closing shared junctions.
    pub junction_movement_limit: f64,
    pub precision: f64,
    pub minimum_edge: f64,
}
#[derive(Debug, Serialize)]
pub struct Issue {
    pub patch: usize,
    pub source_elements: Vec<u32>,
    pub reason: String,
    pub boundary_source_nodes: Vec<Vec<u32>>,
}
#[derive(Debug, Serialize)]
pub struct RegionSplit {
    pub patch: usize,
    pub stiffness: u32,
    pub source_elements: Vec<u32>,
    pub parts: Vec<Vec<u32>>,
}
#[derive(Debug, Serialize)]
pub struct RemovedSliver {
    pub patch: usize,
    pub stiffness: u32,
    pub source_elements: Vec<u32>,
    pub area: f64,
    /// Twice the area over the perimeter: the width of a strip, half the
    /// height of a triangle.
    pub mean_width: f64,
    /// Largest smallest extent of a source facet of the region.
    pub width: f64,
}
#[derive(Debug, Serialize)]
pub struct Report {
    pub policy: Policy,
    /// This stage assembles surface property regions; axes and inter-surface
    /// intersections still require reconciliation. This is not an export gate.
    pub export_ready: bool,
    pub all_surface_patches_built: bool,
    pub preview: Model,
    pub vertex_source_nodes: Vec<u32>,
    pub surface_source_patches: Vec<usize>,
    pub surface_stiffness: Vec<u32>,
    pub pinched_region_splits: Vec<RegionSplit>,
    pub hole_recovery: Vec<HoleRecovery>,
    pub feature_policy: Option<FeaturePolicy>,
    pub simplified_holes: Vec<SimplifiedHole>,
    pub axis_assembly: bars::Report,
    /// Surface-surface junction lines inserted as shared edges. Vertices with
    /// index `>= vertex_source_nodes.len()` are generated junction vertices.
    pub junctions: junctions::Report,
    /// Thin consoles beyond junction lines trimmed in geotechnical assembly.
    pub consoles: consoles::Report,
    /// Walls aligned to the plane of the wall carrying them (geotechnical).
    pub stacked_walls: stacking::Report,
    /// Contour vertices moved onto the plane their edge runs along.
    pub straightened_edges: Vec<cleanup::Straightened>,
    /// Redundant collinear vertices removed at short edges (geotechnical).
    pub short_edges: cleanup::Report,
    /// Duplicated vertices merged before junction insertion (geotechnical).
    pub coincident_vertices: cleanup::MergeReport,
    /// Wall ends merged into nearby vertices (geotechnical).
    pub wall_ends: cleanup::MergeReport,
    /// Bar ends merged across small gaps (geotechnical).
    pub bar_ends: cleanup::MergeReport,
    /// Surface vertices identified with nearly coincident bar nodes.
    pub bar_anchors: cleanup::MergeReport,
    /// Bar ends joined to the span of a nearby bar (geotechnical).
    pub bar_tees: cleanup::TeeReport,
    /// Short edges between needed corners collapsed into one vertex.
    pub short_edge_merges: cleanup::MergeReport,
    /// Bar pieces shorter than the edge collapse tolerance, collapsed.
    pub short_bars: cleanup::BarCollapseReport,
    /// Gaps between structures closed for PLAXIS.
    pub gaps: gaps::Report,
    /// Region contours rebuilt across cracks of the source mesh.
    pub cracks: Vec<cracks::Closure>,
    /// Property regions no wider than a crack (degenerate source slivers),
    /// left out of the geometry (geotechnical).
    pub removed_slivers: Vec<RemovedSliver>,
    pub issues: Vec<Issue>,
    pub maximum_closure_movement: f64,
    pub rejected_vertices: BTreeMap<u32, String>,
    pub support_representatives: Vec<usize>,
    pub support_offset_projection_applied: bool,
}

fn boundary(
    mesh: &MeshData,
    elements: &BTreeMap<u32, &crate::input::ElementData>,
    ids: &[u32],
    plane: &PlaneFrame,
    precision: f64,
) -> Result<Vec<Vec<u32>>, &'static str> {
    let mut counts = BTreeMap::<[u32; 2], usize>::new();
    for id in ids {
        let e = elements.get(id).ok_or("missing_source_element")?;
        let nodes =
            planes::ordered_facet_nodes(mesh, e, plane, precision).ok_or("invalid_source_face")?;
        for i in 0..nodes.len() {
            let (a, b) = (nodes[i], nodes[(i + 1) % nodes.len()]);
            *counts.entry([a.min(b), a.max(b)]).or_default() += 1;
        }
    }
    if counts.values().any(|&n| n > 2) {
        return Err("nonmanifold_source_edges");
    }
    let mut adjacency = BTreeMap::<u32, Vec<u32>>::new();
    for ([a, b], count) in counts {
        if count == 1 {
            adjacency.entry(a).or_default().push(b);
            adjacency.entry(b).or_default().push(a);
        }
    }
    if adjacency.is_empty() || adjacency.values().any(|v| v.len() != 2) {
        return Err("ambiguous_boundary");
    }
    let mut remaining: BTreeSet<_> = adjacency.keys().copied().collect();
    let mut rings = vec![];
    while let Some(&start) = remaining.first() {
        let mut ring = vec![];
        let (mut previous, mut current) = (start, start);
        loop {
            if !remaining.remove(&current) {
                return Err("invalid_boundary_cycle");
            }
            ring.push(current);
            let next = adjacency[&current]
                .iter()
                .copied()
                .find(|n| *n != previous)
                .ok_or("invalid_boundary_cycle")?;
            previous = current;
            current = next;
            if current == start {
                break;
            }
        }
        rings.push(ring);
    }
    Ok(rings)
}

/// Open a pinched boundary by isolating its incident face fans along existing
/// source edges. This partitions material, never fills a hole or moves a node.
/// Every successful recursion strictly reduces the source group; unresolved
/// nonmanifold cases still reach the normal rejecting boundary validator.
fn split_pinched_regions(ids: &[u32], facets: &BTreeMap<u32, Vec<u32>>) -> Vec<Vec<u32>> {
    let mut edges = BTreeMap::<[u32; 2], Vec<u32>>::new();
    for &id in ids {
        let ns = &facets[&id];
        for i in 0..ns.len() {
            let (a, b) = (ns[i], ns[(i + 1) % ns.len()]);
            edges.entry([a.min(b), a.max(b)]).or_default().push(id);
        }
    }
    if edges.values().any(|owners| owners.len() > 2) {
        return vec![ids.to_vec()];
    }
    let mut degree = BTreeMap::<u32, usize>::new();
    for (edge, owners) in &edges {
        if owners.len() == 1 {
            for &n in edge {
                *degree.entry(n).or_default() += 1;
            }
        }
    }
    let pinch: BTreeSet<_> = degree
        .into_iter()
        .filter_map(|(n, d)| (d > 2).then_some(n))
        .collect();
    if pinch.is_empty() {
        return vec![ids.to_vec()];
    }
    let incident: BTreeSet<_> = ids
        .iter()
        .copied()
        .filter(|id| facets[id].iter().any(|n| pinch.contains(n)))
        .collect();
    let mut neighbors = BTreeMap::<u32, Vec<u32>>::new();
    for owners in edges.values() {
        if let [a, b] = owners.as_slice() {
            if incident.contains(a) == incident.contains(b) {
                neighbors.entry(*a).or_default().push(*b);
                neighbors.entry(*b).or_default().push(*a);
            }
        }
    }
    let mut remaining: BTreeSet<_> = ids.iter().copied().collect();
    let mut groups = vec![];
    while let Some(&seed) = remaining.first() {
        let mut stack = vec![seed];
        let mut group = vec![];
        while let Some(id) = stack.pop() {
            if remaining.remove(&id) {
                group.push(id);
                stack.extend(neighbors.get(&id).into_iter().flatten().copied());
            }
        }
        group.sort_unstable();
        groups.push(group);
    }
    if groups.len() == 1 {
        // Every element touches a pinch (fans of triangles meeting at their
        // centres): split the elements around the first pinch into fans,
        // elements joined by an edge through that pinch node.
        let p = *pinch.first().unwrap();
        let mut fans = BTreeMap::<u32, BTreeSet<u32>>::new();
        for (edge, owners) in &edges {
            if let [a, b] = owners.as_slice() {
                if edge.contains(&p) {
                    fans.entry(*a).or_default().insert(*b);
                    fans.entry(*b).or_default().insert(*a);
                }
            }
        }
        let mut remaining: BTreeSet<_> = ids.iter().copied().collect();
        let mut parts = vec![];
        while let Some(&seed) = remaining.first() {
            let around = |id: u32| facets[&id].contains(&p);
            let mut stack = vec![seed];
            let mut part = vec![];
            while let Some(id) = stack.pop() {
                if remaining.remove(&id) {
                    part.push(id);
                    if around(id) {
                        stack.extend(fans.get(&id).into_iter().flatten().copied());
                    } else {
                        stack.extend(
                            neighbors
                                .get(&id)
                                .into_iter()
                                .flatten()
                                .copied()
                                .filter(|&n| !around(n)),
                        );
                    }
                }
            }
            part.sort_unstable();
            parts.push(part);
        }
        if parts.len() == 1 {
            return parts;
        }
        return parts
            .iter()
            .flat_map(|g| split_pinched_regions(g, facets))
            .collect();
    }
    groups
        .iter()
        .flat_map(|g| split_pinched_regions(g, facets))
        .collect()
}

/// Split property regions by edge connectivity; a point contact does not
/// turn two separate material areas into one polygon with an invalid hole.
fn property_regions(
    mesh: &MeshData,
    source: &frame::Report,
    precision: f64,
    splits: &mut Vec<RegionSplit>,
) -> Result<Vec<(usize, u32, Vec<u32>)>, &'static str> {
    let elements: BTreeMap<_, _> = mesh.elements.iter().map(|e| (e.id, e)).collect();
    let mut result = vec![];
    for (patch, surface) in source.surfaces.iter().enumerate() {
        let represented: Vec<_> = surface
            .stiffness_regions
            .values()
            .flatten()
            .copied()
            .collect();
        if represented.len() != represented.iter().collect::<BTreeSet<_>>().len()
            || represented.iter().copied().collect::<BTreeSet<_>>()
                != surface.source_elements.iter().copied().collect()
        {
            return Err("invalid property coverage");
        }
        for (&stiffness, ids) in &surface.stiffness_regions {
            let mut facets = BTreeMap::new();
            let mut edges = BTreeMap::<[u32; 2], Vec<u32>>::new();
            let mut neighbors = BTreeMap::<u32, BTreeSet<u32>>::new();
            for &id in ids {
                let e = elements.get(&id).ok_or("missing property element")?;
                if e.stiff_id != stiffness {
                    return Err("property mismatch");
                }
                neighbors.entry(id).or_default();
                let ns = planes::ordered_facet_nodes(
                    mesh,
                    e,
                    &source.candidate_planes[patch],
                    precision,
                )
                .ok_or("invalid property facet")?;
                facets.insert(id, ns.clone());
                for i in 0..ns.len() {
                    let (a, b) = (ns[i], ns[(i + 1) % ns.len()]);
                    edges.entry([a.min(b), a.max(b)]).or_default().push(id);
                }
            }
            for owners in edges.values() {
                for &a in owners {
                    for &b in owners {
                        if a != b {
                            neighbors.entry(a).or_default().insert(b);
                        }
                    }
                }
            }
            let mut remaining: BTreeSet<_> = ids.iter().copied().collect();
            while let Some(&seed) = remaining.first() {
                let mut stack = vec![seed];
                let mut group = vec![];
                while let Some(id) = stack.pop() {
                    if remaining.remove(&id) {
                        group.push(id);
                        stack.extend(&neighbors[&id]);
                    }
                }
                group.sort_unstable();
                let parts = split_pinched_regions(&group, &facets);
                if parts.len() > 1 {
                    splits.push(RegionSplit {
                        patch,
                        stiffness,
                        source_elements: group,
                        parts: parts.clone(),
                    });
                }
                result.extend(parts.into_iter().map(|part| (patch, stiffness, part)));
            }
        }
    }
    Ok(result)
}

/// Minimum movement onto the intersection of fixed support planes, using an
/// orthonormal constraint basis. Dependent inconsistent planes are rejected.
fn intersection(point: DVec3, planes: &[&PlaneFrame], precision: f64) -> Option<DVec3> {
    let mut basis: Vec<(DVec3, f64)> = vec![];
    for plane in planes {
        let mut n = DVec3::from_array(plane.normal);
        let mut rhs = -plane.distance(point.to_array());
        for &(q, t) in &basis {
            let a = n.dot(q);
            n -= a * q;
            rhs -= a * t;
        }
        let length = n.length();
        if length > 1e-10 && basis.len() < 3 {
            basis.push((n / length, rhs / length));
        }
    }
    let result = point + basis.iter().map(|(n, t)| *n * *t).sum::<DVec3>();
    (result.is_finite()
        && planes
            .iter()
            .all(|p| p.distance(result.to_array()).abs() <= precision))
    .then_some(result)
}

/// A support takes its concurrent offset only if all nodes of its surfaces
/// (support, points) stay within `tolerance` of it; one support that does
/// not fit no longer withdraws the projection from the whole model. Returns
/// the supports and whether every one took its proposal.
fn fitting_supports(
    proposal: Vec<PlaneFrame>,
    candidate: &[PlaneFrame],
    surfaces: impl Iterator<Item = (usize, Vec<[f64; 3]>)>,
    tolerance: f64,
) -> (Vec<PlaneFrame>, bool) {
    let mut fits = vec![true; candidate.len()];
    for (r, points) in surfaces {
        if !points
            .iter()
            .all(|p| proposal[r].distance(*p).abs() <= tolerance)
        {
            fits[r] = false;
        }
    }
    let all = fits.iter().all(|&f| f);
    let supports = proposal
        .into_iter()
        .zip(candidate)
        .zip(&fits)
        .map(|((p, c), &f)| if f { p } else { c.clone() })
        .collect();
    (supports, all)
}

/// Enforce concurrence by projecting support offsets onto the linear
/// compatibility constraints. Normals remain fixed; no per-surface node copies.
fn concurrent_supports(planes: &[PlaneFrame], junctions: &[Vec<usize>]) -> Vec<PlaneFrame> {
    if planes.is_empty() {
        return vec![];
    }
    let center = DVec3::from_array(planes[0].origin);
    let offsets: Vec<_> = planes
        .iter()
        .map(|p| DVec3::from_array(p.normal).dot(DVec3::from_array(p.origin) - center))
        .collect();
    let mut constraints: Vec<Vec<f64>> = vec![];
    for junction in junctions {
        let mut basis: Vec<(DVec3, Vec<f64>)> = vec![];
        for &i in junction {
            let mut n = DVec3::from_array(planes[i].normal);
            let mut row = vec![0.; planes.len()];
            row[i] = 1.;
            for (q, r) in &basis {
                let a = n.dot(*q);
                n -= a * *q;
                for (v, b) in row.iter_mut().zip(r) {
                    *v -= a * b;
                }
            }
            let length = n.length();
            if length > 1e-10 && basis.len() < 3 {
                for v in &mut row {
                    *v /= length;
                }
                basis.push((n / length, row));
            } else {
                // Reorthogonalize to avoid amplifying repeated junction rows.
                for _ in 0..2 {
                    for q in &constraints {
                        let a: f64 = row.iter().zip(q).map(|(a, b)| a * b).sum();
                        for (v, b) in row.iter_mut().zip(q) {
                            *v -= a * b;
                        }
                    }
                }
                let length = row.iter().map(|v| v * v).sum::<f64>().sqrt();
                if length > 1e-10 {
                    for v in &mut row {
                        *v /= length;
                    }
                    constraints.push(row);
                }
            }
        }
    }
    let mut corrected = offsets.clone();
    for row in constraints {
        let error: f64 = row.iter().zip(&offsets).map(|(a, b)| a * b).sum();
        for (d, a) in corrected.iter_mut().zip(row) {
            *d -= error * a;
        }
    }
    planes
        .iter()
        .enumerate()
        .map(|(i, p)| {
            let n = DVec3::from_array(p.normal);
            let mut result = p.clone();
            result.origin =
                (DVec3::from_array(p.origin) + n * (corrected[i] - offsets[i])).to_array();
            result
        })
        .collect()
}

/// Stage timing on stderr when `TOPO_TIMING` is set (diagnostics only).
struct Timer {
    enabled: bool,
    last: std::time::Instant,
}
impl Timer {
    fn new() -> Self {
        Self {
            enabled: std::env::var_os("TOPO_TIMING").is_some(),
            last: std::time::Instant::now(),
        }
    }
    fn lap(&mut self, stage: &str) {
        if self.enabled {
            eprintln!(
                "[timing] assembly {stage}: {:.1}s",
                self.last.elapsed().as_secs_f64()
            );
        }
        self.last = std::time::Instant::now();
    }
}

fn movement_budget(_mesh: &MeshData, source: &frame::Report, index: usize) -> f64 {
    source.movement_budgets()[index]
}

/// Apply the crack mouths identified in any region of a patch to the rings
/// of every region of that patch; consecutive repeated nodes collapse.
/// Remove repeated consecutive nodes and back-and-forth spikes `a, b, a`
/// (cyclically) from a ring whose nodes were identified: a mouth applied from
/// another region of the patch can fold a straight boundary through a node
/// that the crack removed there.
fn collapse_backtracks(ring: &mut Vec<u32>) {
    loop {
        ring.dedup();
        while ring.len() > 1 && ring.first() == ring.last() {
            ring.pop();
        }
        let n = ring.len();
        if n < 3 {
            return;
        }
        let Some(i) = (0..n).find(|&i| ring[(i + n - 1) % n] == ring[(i + 1) % n]) else {
            return;
        };
        // Drop the spike tip and the repeated node after it.
        let next = (i + 1) % n;
        let (first, second) = (i.max(next), i.min(next));
        ring.remove(first);
        ring.remove(second);
    }
}

fn share_crack_mouths(
    rings: &mut BTreeMap<usize, Vec<Vec<u32>>>,
    patch: impl Fn(usize) -> usize,
    cracks: &[cracks::Closure],
) {
    let mut mouths = BTreeMap::<usize, BTreeMap<u32, u32>>::new();
    for closure in cracks {
        let map = mouths.entry(closure.patch).or_default();
        for &[dropped, kept] in &closure.identified {
            map.insert(dropped, kept);
        }
    }
    for (&i, loops) in rings.iter_mut() {
        let Some(map) = mouths.get(&patch(i)) else {
            continue;
        };
        let resolve = |mut n: u32| {
            for _ in 0..map.len() {
                match map.get(&n) {
                    Some(&k) => n = k,
                    None => break,
                }
            }
            n
        };
        for ring in loops.iter_mut() {
            if ring.iter().any(|n| map.contains_key(n)) {
                let mut mapped: Vec<u32> = ring.iter().map(|&n| resolve(n)).collect();
                collapse_backtracks(&mut mapped);
                *ring = mapped;
            }
        }
    }
}

/// Area, mean width (twice the area over the boundary length) and largest
/// facet width (the smallest extent of a convex facet, attained across one
/// of its edges) of a property region; `None` for an invalid facet.
fn region_width(
    mesh: &MeshData,
    elements: &BTreeMap<u32, &crate::input::ElementData>,
    ids: &[u32],
    plane: &PlaneFrame,
    precision: f64,
) -> Option<(f64, f64, f64)> {
    let mut area = 0.;
    let mut widest: f64 = 0.;
    let mut edges = BTreeMap::<[u32; 2], usize>::new();
    for id in ids {
        let ns = planes::ordered_facet_nodes(mesh, elements.get(id)?, plane, precision)?;
        area += ring_area(&ns, plane, |n| mesh.nodes[&n].to_array());
        let uv: Vec<[f64; 2]> = ns
            .iter()
            .map(|n| plane.project(mesh.nodes[n].to_array()))
            .collect();
        let width = (0..uv.len())
            .map(|i| {
                let (a, b) = (uv[i], uv[(i + 1) % uv.len()]);
                let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
                let length = dx.hypot(dy);
                uv.iter()
                    .map(|p| ((p[0] - a[0]) * dy - (p[1] - a[1]) * dx).abs() / length)
                    .fold(0., f64::max)
            })
            .fold(f64::INFINITY, f64::min);
        widest = widest.max(width);
        for i in 0..ns.len() {
            let (a, b) = (ns[i], ns[(i + 1) % ns.len()]);
            *edges.entry([a.min(b), a.max(b)]).or_default() += 1;
        }
    }
    let perimeter: f64 = edges
        .iter()
        .filter(|(_, &count)| count == 1)
        .map(|([a, b], _)| mesh.nodes[a].distance(mesh.nodes[b]))
        .sum();
    (perimeter > 0.).then(|| (area, 2. * area / perimeter, widest))
}

fn ring_area(ring: &[u32], plane: &PlaneFrame, point: impl Fn(u32) -> [f64; 3]) -> f64 {
    let uv: Vec<_> = ring.iter().map(|&n| plane.project(point(n))).collect();
    (0..uv.len())
        .map(|i| {
            let (a, b, o) = (uv[i], uv[(i + 1) % uv.len()], uv[0]);
            (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
        })
        .sum::<f64>()
        .abs()
        / 2.
}

pub fn assemble(
    mesh: &MeshData,
    source: &frame::Report,
    policy: &Policy,
) -> Result<Report, &'static str> {
    assemble_impl(mesh, source, policy, None)
}

pub fn assemble_geotechnical(
    mesh: &MeshData,
    source: &frame::Report,
    policy: &Policy,
    features: &FeaturePolicy,
) -> Result<Report, &'static str> {
    if !features.maximum_source_width.is_finite()
        || features.maximum_source_width <= 0.
        || !features.maximum_filled_area_ratio.is_finite()
        || features.maximum_filled_area_ratio <= 0.
        || features.maximum_filled_area_ratio >= 1.
        || !features.maximum_console_width.is_finite()
        || features.maximum_console_width < 0.
        || !features.maximum_stack_offset.is_finite()
        || features.maximum_stack_offset < 0.
        || !features.maximum_wall_end_snap.is_finite()
        || features.maximum_wall_end_snap < 0.
        || !features.maximum_crack_width.is_finite()
        || features.maximum_crack_width < 0.
        || !features.maximum_collapsed_edge.is_finite()
        || features.maximum_collapsed_edge < 0.
        || !features.maximum_gap.is_finite()
        || features.maximum_gap < 0.
    {
        return Err("invalid feature simplification policy");
    }
    assemble_impl(mesh, source, policy, Some(features))
}

fn assemble_impl(
    mesh: &MeshData,
    source: &frame::Report,
    policy: &Policy,
    features: Option<&FeaturePolicy>,
) -> Result<Report, &'static str> {
    if source.sliding_parameters.is_some() {
        return Err("sliding frame is proposal-only until parameter transfer is implemented");
    }
    if !policy.closure_tolerance.is_finite()
        || policy.closure_tolerance < policy.precision
        || !policy.junction_movement_limit.is_finite()
        || policy.junction_movement_limit < policy.precision
    {
        return Err("invalid closure tolerance");
    }
    let mut model =
        Model::new(policy.precision, policy.minimum_edge).map_err(|_| "invalid assembly policy")?;
    if source.candidate_planes.len() != source.surfaces.len()
        || source.candidate_points.len() != source.node_ids.len()
    {
        return Err("incomplete frame proposal");
    }
    let lookup: BTreeMap<_, _> = source
        .node_ids
        .iter()
        .enumerate()
        .map(|(i, &n)| (n, i))
        .collect();
    let mut issues = vec![];
    let mut rings = BTreeMap::new();
    let mut pinched_region_splits = vec![];
    let mut timer = Timer::new();
    let mut regions = property_regions(mesh, source, policy.precision, &mut pinched_region_splits)?;
    // A region no wider than a crack, facet by facet and on average, is a
    // degenerate sliver of the source mesh (needle triangles along a line),
    // not a structure.
    let elements: BTreeMap<_, _> = mesh.elements.iter().map(|e| (e.id, e)).collect();
    let mut removed_slivers = vec![];
    if let Some(width) = features.map(|f| f.maximum_crack_width).filter(|&w| w > 0.) {
        regions.retain(|(patch, stiffness, ids)| {
            let plane = &source.candidate_planes[*patch];
            match region_width(mesh, &elements, ids, plane, policy.precision) {
                Some((area, mean, widest)) if mean <= width && widest <= width => {
                    removed_slivers.push(RemovedSliver {
                        patch: *patch,
                        stiffness: *stiffness,
                        source_elements: ids.clone(),
                        area,
                        mean_width: mean,
                        width: widest,
                    });
                    false
                }
                _ => true,
            }
        });
    }
    timer.lap("property_regions");
    let mut cracks = vec![];
    for (i, (patch, _, ids)) in regions.iter().enumerate() {
        let plane = &source.candidate_planes[*patch];
        let result = boundary(mesh, &elements, ids, plane, policy.precision).map(|r| {
            // A crack of the source mesh is left out of the region contour.
            let closed = features
                .filter(|f| f.maximum_crack_width > 0.)
                .and_then(|f| {
                    cracks::close(
                        mesh,
                        &elements,
                        ids,
                        &r,
                        plane,
                        policy.precision,
                        f.maximum_crack_width,
                        *patch,
                    )
                });
            match closed {
                Some((rebuilt, closure)) => {
                    cracks.push(closure);
                    rebuilt
                }
                None => r,
            }
        });
        match result {
            Ok(mut r) => {
                // Fix exterior/hole roles from the immutable source geometry.
                let area = |ring: &Vec<u32>| {
                    ring_area(ring, &source.candidate_planes[*patch], |n| {
                        mesh.nodes[&n].to_array()
                    })
                };
                r.sort_by(|a, b| area(b).total_cmp(&area(a)));
                rings.insert(i, r);
            }
            Err(reason) => issues.push(Issue {
                patch: *patch,
                source_elements: ids.clone(),
                reason: reason.into(),
                boundary_source_nodes: vec![],
            }),
        }
    }
    // A crack mouth identified in one region is one point for every region
    // of the same patch (parts split off at the crack share its nodes).
    share_crack_mouths(&mut rings, |i| regions[i].0, &cracks);
    // Include every support owning a boundary node, including a patch whose
    // own boundary failed. Never silently disconnect it from valid neighbors.
    timer.lap("region_boundaries");
    let boundary_nodes: BTreeSet<_> = rings.values().flatten().flatten().copied().collect();
    let mut owners = BTreeMap::<u32, Vec<usize>>::new();
    for (i, s) in source.surfaces.iter().enumerate() {
        for &index in &s.nodes {
            let id = *source.node_ids.get(index).ok_or("invalid surface node")?;
            if boundary_nodes.contains(&id) {
                owners.entry(id).or_default().push(i);
            }
        }
    }
    timer.lap("owners");
    // Connected nearly coplanar source patches can use one support. Validate
    // every candidate point against the chosen support, preventing chain drift.
    let mut support_representatives: Vec<_> = (0..source.surfaces.len()).collect();
    let mut order: Vec<_> = (0..source.surfaces.len()).collect();
    order.sort_by_key(|&i| std::cmp::Reverse(source.surfaces[i].nodes.len()));
    let mut assigned = BTreeSet::new();
    for master in order {
        if !assigned.insert(master) {
            continue;
        }
        let plane = &source.candidate_planes[master];
        let mut connected: BTreeSet<_> = source.surfaces[master].nodes.iter().copied().collect();
        loop {
            let mut changed = false;
            for j in 0..source.surfaces.len() {
                if assigned.contains(&j) {
                    continue;
                }
                let surface = &source.surfaces[j];
                if surface.nodes.iter().any(|n| connected.contains(n))
                    && DVec3::from_array(plane.normal)
                        .dot(DVec3::from_array(source.candidate_planes[j].normal))
                        .abs()
                        >= source.policy.angle.cos()
                    && surface.nodes.iter().all(|&n| {
                        plane.distance(source.candidate_points[n]).abs() <= policy.closure_tolerance
                    })
                {
                    assigned.insert(j);
                    support_representatives[j] = master;
                    connected.extend(surface.nodes.iter().copied());
                    changed = true;
                }
            }
            if !changed {
                break;
            }
        }
    }
    timer.lap("regions_and_representatives");
    let mut stacked_walls = match features {
        Some(features) if features.maximum_stack_offset > 0. => stacking::align(
            mesh,
            source,
            policy,
            &owners,
            &rings,
            &mut support_representatives,
            features.maximum_stack_offset,
        ),
        _ => stacking::Report::default(),
    };
    if let Some(features) = features.filter(|f| f.maximum_stack_offset > 0.) {
        stacking::align_lines(
            mesh,
            source,
            policy,
            &owners,
            &rings,
            &mut support_representatives,
            features.maximum_stack_offset,
            &mut stacked_walls,
        );
    }
    // Candidate points of an aligned wall are intentionally off its support.
    let aligned: BTreeSet<usize> = (0..source.surfaces.len())
        .filter(|&i| {
            stacked_walls
                .aligned
                .iter()
                .any(|w| support_representatives[w.upper] == support_representatives[i])
        })
        .collect();
    let junctions: Vec<Vec<usize>> = owners
        .values()
        .map(|ids| {
            ids.iter()
                .map(|&i| support_representatives[i])
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect()
        })
        .collect();
    timer.lap("stacked_walls");
    let proposal = concurrent_supports(&source.candidate_planes, &junctions);
    let (closed_supports, support_offsets_adjusted) = fitting_supports(
        proposal,
        &source.candidate_planes,
        source
            .surfaces
            .iter()
            .enumerate()
            .filter(|(i, _)| !aligned.contains(i))
            .map(|(i, s)| {
                (
                    support_representatives[i],
                    s.nodes
                        .iter()
                        .map(|&n| source.candidate_points[n])
                        .collect::<Vec<_>>(),
                )
            }),
        policy.closure_tolerance,
    );
    timer.lap("concurrent_supports");
    let mut rejected_vertices = BTreeMap::new();
    let mut closed_points = BTreeMap::new();
    for (&id, supports) in &owners {
        let i = lookup[&id];
        let p = DVec3::from_array(source.candidate_points[i]);
        let planes: Vec<_> = supports
            .iter()
            .map(|&s| &closed_supports[support_representatives[s]])
            .collect();
        let Some(q) = intersection(p, &planes, policy.precision) else {
            rejected_vertices.insert(id, "inconsistent_supports".into());
            continue;
        };
        let movement = p.distance(q);
        let reference = *mesh.nodes.get(&id).ok_or("missing reference node")?;
        // Cumulative movement from immutable input, including local axis budgets.
        let budget = movement_budget(mesh, source, i);
        let (limit, k) = stacking::movement_limit(
            policy.junction_movement_limit,
            supports,
            |s| stacked_walls.shifted.contains(&s),
            |s| support_representatives[s],
        );
        if k >= 2 && movement > policy.junction_movement_limit && movement <= limit {
            stacked_walls.raised_limits.push(stacking::RaisedLimit {
                source_node: id,
                aligned_supports: k,
                limit,
            });
        }
        if movement > limit || q.distance(reference) > budget + policy.precision {
            rejected_vertices.insert(
                id,
                format!(
                    "movement: closure={movement}, total={}, budget={budget}",
                    q.distance(reference)
                ),
            );
            continue;
        }
        closed_points.insert(id, q);
    }
    timer.lap("vertex_closure");
    let hole_recovery = holes::recover(
        &mut closed_points,
        &holes::Context {
            mesh,
            source,
            policy,
            regions: &regions,
            rings: &rings,
            owners: &owners,
            representatives: &support_representatives,
            supports: &closed_supports,
        },
    );
    let mut vertices = BTreeMap::new();
    let mut vertex_source_nodes = vec![];
    let mut maximum_closure_movement = 0.0_f64;
    // Identified stacked-wall nodes share the lower node's vertex, if that
    // position is within the upper node's own movement limits.
    let mut identified = BTreeMap::new();
    for pair in &stacked_walls.identified {
        let (Some(&q), Some(&i)) = (
            closed_points.get(&pair.lower_node),
            lookup.get(&pair.upper_node),
        ) else {
            continue;
        };
        if !closed_points.contains_key(&pair.upper_node) {
            continue;
        }
        let p = DVec3::from_array(source.candidate_points[i]);
        let reference = *mesh
            .nodes
            .get(&pair.upper_node)
            .ok_or("missing reference node")?;
        // The shared position must lie on every support of the upper node.
        let on_supports = owners.get(&pair.upper_node).is_none_or(|patches| {
            patches.iter().all(|&patch| {
                closed_supports[support_representatives[patch]]
                    .distance(q.to_array())
                    .abs()
                    <= policy.precision
            })
        });
        let (limit, _) = stacking::movement_limit(
            policy.junction_movement_limit,
            owners.get(&pair.upper_node).map_or(&[][..], |v| &v[..]),
            |s| stacked_walls.shifted.contains(&s),
            |s| support_representatives[s],
        );
        if on_supports
            && p.distance(q) <= limit
            && q.distance(reference) <= movement_budget(mesh, source, i) + policy.precision
        {
            identified.insert(pair.upper_node, pair.lower_node);
            closed_points.insert(pair.upper_node, q);
        }
    }
    // Source nodes on the same supports (a corner of three planes, a point
    // of one junction line) closed within the minimum edge length are one
    // vertex; typically an aligned wall's corner reaching the corner of the
    // wall below. Bar anchors keep their own nodes.
    let bar_nodes: BTreeSet<u32> = source
        .axes
        .iter()
        .flat_map(|a| a.anchors.iter().map(|x| source.node_ids[x.node]))
        .collect();
    let mut by_supports = BTreeMap::<Vec<usize>, Vec<u32>>::new();
    for &id in closed_points.keys() {
        if identified.contains_key(&id) || bar_nodes.contains(&id) {
            continue;
        }
        let Some(supports) = owners.get(&id) else {
            continue;
        };
        let key: BTreeSet<usize> = supports
            .iter()
            .map(|&s| support_representatives[s])
            .collect();
        if key.len() >= 2 {
            by_supports
                .entry(key.into_iter().collect())
                .or_default()
                .push(id);
        }
    }
    for ids in by_supports.values() {
        let mut kept: Vec<u32> = vec![];
        for &id in ids {
            let q = closed_points[&id];
            match kept
                .iter()
                .find(|k| closed_points[k].distance(q) < policy.minimum_edge)
            {
                Some(&k) => {
                    stacked_walls.coincident.push(stacking::Identified {
                        upper_node: id,
                        lower_node: k,
                        distance: closed_points[&k].distance(q),
                    });
                    identified.insert(id, k);
                }
                None => kept.push(id),
            }
        }
    }
    // Consecutive contour nodes closed within the minimum edge length (a
    // degenerate source element edge, 0.7 mm) are one vertex when the
    // supports of one contain the other's: the kept node lies on every
    // support of the dropped one. Bar anchors are never dropped.
    let supports_of = |id: u32| -> BTreeSet<usize> {
        owners
            .get(&id)
            .map(|s| s.iter().map(|&x| support_representatives[x]).collect())
            .unwrap_or_default()
    };
    for ring in rings.values().flatten() {
        for i in 0..ring.len() {
            let root = |mut id: u32| {
                while let Some(&k) = identified.get(&id) {
                    id = k;
                }
                id
            };
            let (a, b) = (root(ring[i]), root(ring[(i + 1) % ring.len()]));
            if a == b {
                continue;
            }
            let (Some(&qa), Some(&qb)) = (closed_points.get(&a), closed_points.get(&b)) else {
                continue;
            };
            if qa.distance(qb) >= policy.minimum_edge {
                continue;
            }
            let (sa, sb) = (supports_of(a), supports_of(b));
            let (keep, drop) = if sb.is_subset(&sa) && !bar_nodes.contains(&b) {
                (a, b)
            } else if sa.is_subset(&sb) && !bar_nodes.contains(&a) {
                (b, a)
            } else {
                continue;
            };
            stacked_walls.coincident.push(stacking::Identified {
                upper_node: drop,
                lower_node: keep,
                distance: qa.distance(qb),
            });
            identified.insert(drop, keep);
        }
    }
    for (&id, &q) in &closed_points {
        if identified.contains_key(&id) {
            continue;
        }
        maximum_closure_movement = maximum_closure_movement
            .max(q.distance(DVec3::from_array(source.candidate_points[lookup[&id]])));
        vertices.insert(
            id,
            model
                .add_vertex(q.to_array())
                .map_err(|_| "invalid closed vertex")?,
        );
        vertex_source_nodes.push(id);
    }
    for &upper in identified.keys() {
        let mut lower = upper;
        while let Some(&k) = identified.get(&lower) {
            lower = k;
        }
        vertices.insert(upper, vertices[&lower]);
    }
    stacked_walls
        .identified
        .retain(|pair| identified.contains_key(&pair.upper_node));
    let mut surface_source_patches = vec![];
    let mut surface_stiffness = vec![];
    let mut simplified_holes = vec![];
    for (region, mut loops) in rings {
        let (patch, stiffness, ids) = &regions[region];
        let patch = *patch;
        if loops.iter().flatten().any(|n| !vertices.contains_key(n)) {
            issues.push(Issue {
                patch,
                source_elements: ids.clone(),
                reason: "support_intersection_or_movement_budget".into(),
                boundary_source_nodes: loops.clone(),
            });
            continue;
        }
        let plane = &closed_supports[support_representatives[patch]];
        let changes = if let Some(features) = features {
            features::simplify(
                &mut loops,
                mesh,
                &closed_points,
                plane,
                policy.precision,
                features,
                patch,
                ids,
            )
        } else {
            vec![]
        };
        let area = |ring: &Vec<u32>| ring_area(ring, plane, |n| model.vertices[vertices[&n]]);
        let degenerate_hole = loops
            .iter()
            .skip(1)
            .any(|ring| area(ring) <= policy.precision * policy.precision);
        let plane_id = model.add_plane(plane.clone());
        // An edge whose two nodes were identified (a stacked wall aligned
        // onto the wall below) collapses: its contour keeps one vertex.
        let mapped = loops
            .iter()
            .map(|r| {
                let mut ring: Vec<usize> = r.iter().map(|n| vertices[n]).collect();
                ring.dedup();
                while ring.len() > 1 && ring.first() == ring.last() {
                    ring.pop();
                }
                ring
            })
            .collect();
        match model.add_surface(plane_id, mapped, ids.clone()) {
            Ok(_) => {
                simplified_holes.extend(changes);
                surface_source_patches.push(patch);
                surface_stiffness.push(*stiffness);
            }
            Err(error) => issues.push(Issue {
                patch,
                source_elements: ids.clone(),
                reason: if degenerate_hole {
                    "degenerate_hole_after_closure".into()
                } else {
                    format!("contour_{error:?}")
                },
                boundary_source_nodes: loops,
            }),
        }
    }
    timer.lap("holes_and_surfaces");
    let mut axis_assembly = bars::assemble(
        mesh,
        source,
        &mut model,
        &mut vertex_source_nodes,
        &closed_supports,
        &support_representatives,
        policy,
        &removed_slivers
            .iter()
            .flat_map(|s| s.source_elements.iter().copied())
            .collect(),
    );
    maximum_closure_movement =
        maximum_closure_movement.max(axis_assembly.maximum_additional_movement);
    // After axis repairs, which may still move boundary anchors: junctions
    // are derived from final coordinates and never move a vertex.
    // Mandatory interior nodes per surface, the vertices that must stay
    // (all bar vertices and retained nodes), and the retained nodes that no
    // merge may move. Recomputed after merges, which renumber vertices.
    // Nodes of a closed crack void or seam (source mesh discretization) stay
    // only where shared with another surface or a bar: alone they are pairs
    // of vertices millimetres apart inside the material.
    for hole in simplified_holes
        .iter_mut()
        .filter(|hole| hole.reason != "collapsed_opening")
    {
        let own: Vec<usize> = model
            .surfaces
            .iter()
            .enumerate()
            .filter(|(_, surface)| surface.source_elements == hole.source_elements)
            .map(|(s, _)| s)
            .collect();
        let shared = |v: usize| {
            (0..model.surfaces.len())
                .filter(|s| !own.contains(s))
                .any(|s| model.surface_edges(s).any(|e| model.edges[e].contains(&v)))
                || axis_assembly
                    .axes
                    .iter()
                    .any(|a| a.endpoints.contains(&v) || a.anchors.iter().any(|x| x.vertex == v))
        };
        hole.retained_nodes = hole
            .source_nodes
            .iter()
            .copied()
            .filter(|n| {
                vertex_source_nodes
                    .iter()
                    .position(|m| m == n)
                    .is_some_and(shared)
            })
            .collect();
    }
    let hole_nodes: Vec<(usize, Vec<usize>)> = simplified_holes
        .iter()
        .flat_map(|hole| {
            let nodes: Vec<usize> = hole
                .retained_nodes
                .iter()
                .filter_map(|n| vertex_source_nodes.iter().position(|m| m == n))
                .collect();
            model
                .surfaces
                .iter()
                .enumerate()
                .filter(|(_, surface)| surface.source_elements == hole.source_elements)
                .map(|(s, _)| (s, nodes.clone()))
                .collect::<Vec<_>>()
        })
        .collect();
    let protected = |model: &Model, axis_assembly: &bars::Report| {
        let mut interior = vec![BTreeSet::new(); model.surfaces.len()];
        for contact in &axis_assembly.contacts {
            if let bars::Contact::Point {
                surface, vertex, ..
            } = *contact
            {
                interior[surface].insert(vertex);
            }
        }
        let mut fixed = BTreeSet::new();
        for (s, nodes) in &hole_nodes {
            interior[*s].extend(nodes.iter().copied());
            fixed.extend(nodes.iter().copied());
        }
        let mut locked: BTreeSet<usize> = interior.iter().flatten().copied().collect();
        for axis in &axis_assembly.axes {
            locked.extend(axis.endpoints);
            locked.extend(axis.anchors.iter().map(|a| a.vertex));
        }
        (interior, locked, fixed)
    };
    let (_, _, fixed) = protected(&model, &axis_assembly);
    timer.lap("bars");
    let coincident = match features {
        Some(_) => cleanup::merge_coincident(
            &mut model,
            &mut cleanup::Bars {
                axes: &mut axis_assembly.axes,
                contacts: &mut axis_assembly.contacts,
            },
            policy.minimum_edge,
            &fixed,
            &vertex_source_nodes,
        ),
        None => cleanup::MergeReport::default(),
    };
    timer.lap("coincident_merges");
    let bar_ends = match features {
        Some(features) if features.maximum_wall_end_snap > 0. => cleanup::merge_bar_ends(
            &mut model,
            &mut cleanup::Bars {
                axes: &mut axis_assembly.axes,
                contacts: &mut axis_assembly.contacts,
            },
            features.maximum_wall_end_snap,
            &fixed,
            &vertex_source_nodes,
        ),
        _ => cleanup::MergeReport::default(),
    };
    let bar_anchors = match features {
        Some(_) => cleanup::merge_bar_anchors(
            &mut model,
            &mut cleanup::Bars {
                axes: &mut axis_assembly.axes,
                contacts: &mut axis_assembly.contacts,
            },
            policy.minimum_edge,
            &fixed,
            &vertex_source_nodes,
        ),
        None => cleanup::MergeReport::default(),
    };
    let (_, _, fixed) = protected(&model, &axis_assembly);
    // Bar pieces too short for a PLAXIS element collapse.
    let short_bars = match features {
        Some(features) if features.maximum_collapsed_edge > 0. => cleanup::collapse_short_bars(
            &mut model,
            &mut cleanup::Bars {
                axes: &mut axis_assembly.axes,
                contacts: &mut axis_assembly.contacts,
            },
            features.maximum_collapsed_edge,
            &fixed,
            &vertex_source_nodes,
        ),
        _ => cleanup::BarCollapseReport::default(),
    };
    let mut short_bars = short_bars;
    short_bars.duplicates = cleanup::remove_duplicate_bars(&mut cleanup::Bars {
        axes: &mut axis_assembly.axes,
        contacts: &mut axis_assembly.contacts,
    });
    let (_, _, fixed) = protected(&model, &axis_assembly);
    // Bar ends a few millimetres off the span of another bar join it.
    let bar_tees = match features {
        Some(features) if features.maximum_wall_end_snap > 0. => cleanup::join_bar_tees(
            &mut model,
            &mut cleanup::Bars {
                axes: &mut axis_assembly.axes,
                contacts: &mut axis_assembly.contacts,
            },
            features.maximum_wall_end_snap,
            &fixed,
            &vertex_source_nodes,
        ),
        _ => cleanup::TeeReport::default(),
    };
    let mut bar_tees = bar_tees;
    if !bar_tees.joined.is_empty() {
        bar_tees.shared_nodes = cleanup::share_overlapping_bars(&model, &mut axis_assembly.axes);
        let mut bars = cleanup::Bars {
            axes: &mut axis_assembly.axes,
            contacts: &mut axis_assembly.contacts,
        };
        short_bars
            .duplicates
            .extend(cleanup::remove_duplicate_bars(&mut bars));
        short_bars
            .duplicates
            .extend(cleanup::remove_contained_bars(&mut bars));
    }
    let (_, _, fixed) = protected(&model, &axis_assembly);
    timer.lap("bar_end_merges");
    let gaps = match features {
        Some(features) if features.maximum_gap > 0. => gaps::close(
            &mut model,
            &mut cleanup::Bars {
                axes: &mut axis_assembly.axes,
                contacts: &mut axis_assembly.contacts,
            },
            features.maximum_gap,
            features.close_offset_gaps,
            &fixed,
            &vertex_source_nodes,
            &vertex_source_nodes
                .iter()
                .map(|n| source.candidate_points[lookup[n]])
                .collect::<Vec<_>>(),
            policy.junction_movement_limit,
        ),
        _ => gaps::Report::default(),
    };
    // A closed gap can bring a vertex onto another one.
    let mut coincident = coincident;
    if !gaps.closed.is_empty() {
        let (_, _, fixed) = protected(&model, &axis_assembly);
        let again = cleanup::merge_coincident(
            &mut model,
            &mut cleanup::Bars {
                axes: &mut axis_assembly.axes,
                contacts: &mut axis_assembly.contacts,
            },
            policy.minimum_edge,
            &fixed,
            &vertex_source_nodes,
        );
        coincident.merged.extend(again.merged);
        coincident.rejected.extend(again.rejected);
    }
    let (interior, locked, fixed) = protected(&model, &axis_assembly);
    timer.lap("gaps");
    let straightened_edges = if features.is_some() {
        cleanup::straighten_edges(
            &mut model,
            &axis_assembly.axes,
            policy.minimum_edge,
            source.policy.angle.sin(),
            &fixed,
            &vertex_source_nodes,
        )
    } else {
        vec![]
    };
    let junctions = junctions::insert(
        &mut model,
        &junctions::Context {
            interior: &interior,
            locked: &locked,
            wall_end_tolerance: features.map_or(0., |f| f.maximum_wall_end_snap),
        },
    );
    timer.lap("junctions");
    let wall_ends = match features {
        Some(features) if features.maximum_wall_end_snap > 0. => cleanup::merge_wall_ends(
            &mut model,
            &mut cleanup::Bars {
                axes: &mut axis_assembly.axes,
                contacts: &mut axis_assembly.contacts,
            },
            features.maximum_wall_end_snap,
            &fixed,
            &vertex_source_nodes,
        ),
        _ => cleanup::MergeReport::default(),
    };
    let short_edge_merges = match features {
        Some(features) if features.maximum_collapsed_edge > 0. => {
            let (_, _, fixed) = protected(&model, &axis_assembly);
            cleanup::collapse_short_edges(
                &mut model,
                &mut cleanup::Bars {
                    axes: &mut axis_assembly.axes,
                    contacts: &mut axis_assembly.contacts,
                },
                features.maximum_collapsed_edge,
                &fixed,
                &vertex_source_nodes,
            )
        }
        _ => cleanup::MergeReport::default(),
    };
    let (interior, locked, _) = protected(&model, &axis_assembly);
    // A wall end left open by the junction pass may be closed by a merge;
    // its near-touch diagnostic is then obsolete.
    let mut junctions = junctions;
    junctions.issues.retain(|issue| {
        !issue.reason.starts_with("embedded_near_touch")
            || !wall_ends.merged.iter().any(|m| {
                [issue.start, issue.end].iter().any(|p| {
                    DVec3::from_array(*p).distance(DVec3::from_array(m.kept_from))
                        <= policy.precision
                })
            })
    });
    timer.lap("wall_end_merges");
    let consoles = match features {
        Some(features) => consoles::trim(
            &mut model,
            &consoles::Context {
                maximum_width: features.maximum_console_width,
                interior: &interior,
                locked: &locked,
                axes: &axis_assembly.axes,
                contacts: &axis_assembly.contacts,
            },
        ),
        None => consoles::Report::default(),
    };
    timer.lap("consoles");
    let short_edges = match features {
        Some(features) if features.maximum_wall_end_snap > 0. => cleanup::remove_short_edges(
            &mut model,
            features.maximum_wall_end_snap,
            &locked,
            &vertex_source_nodes,
        ),
        _ => cleanup::Report::default(),
    };
    timer.lap("short_edges");
    model.refresh_orphaned_edges();
    axis_assembly.imprinted = bars::imprint_surface_vertices(
        &mut model,
        &mut axis_assembly.axes,
        &vertex_source_nodes,
        features.map_or(0., |f| f.maximum_gap),
    );
    if features.is_some() {
        let mut crossings = bars::imprint_crossings(&mut model, &mut axis_assembly.axes);
        crossings.extend(bars::imprint_bar_crossings(
            &mut model,
            &mut axis_assembly.axes,
        ));
        // Generated vertices (a bar crossing is reported once per bar).
        let generated: BTreeSet<usize> = crossings
            .iter()
            .filter(|c| c.kind == "crossing" || c.kind == "bar_crossing")
            .map(|c| c.vertex)
            .collect();
        junctions.generated_vertices.extend(generated);
        axis_assembly.imprinted.extend(crossings);
    }
    bars::refresh_contacts(&model, &axis_assembly.axes, &mut axis_assembly.contacts);
    Ok(Report {
        policy: policy.clone(),
        export_ready: false,
        all_surface_patches_built: issues.is_empty(),
        preview: model,
        vertex_source_nodes,
        surface_source_patches,
        surface_stiffness,
        pinched_region_splits,
        hole_recovery,
        feature_policy: features.cloned(),
        simplified_holes,
        axis_assembly,
        junctions,
        consoles,
        stacked_walls,
        straightened_edges,
        short_edges,
        coincident_vertices: coincident,
        wall_ends,
        bar_ends,
        bar_anchors,
        bar_tees,
        short_edge_merges,
        short_bars,
        removed_slivers,
        gaps,
        cracks,
        issues,
        maximum_closure_movement,
        rejected_vertices,
        support_representatives,
        support_offset_projection_applied: support_offsets_adjusted,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    fn planar_frame(mesh: &MeshData, up: DVec3) -> frame::Report {
        use super::super::recognize;
        let axes = recognize::recognize(
            mesh,
            &recognize::Policy {
                angle: 0.02,
                line_tolerance: 0.01,
                numerical_precision: 1e-8,
            },
        )
        .unwrap();
        let surfaces = planes::recognize(
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
            &surfaces,
            &frame::Policy {
                up: up.to_array(),
                angle: 0.02,
                maximum_movement: 0.15,
                relative_movement: 0.05,
                minimum_length: 0.03,
                residual_tolerance: 1e-7,
                iterations: 100,
                panel_tolerance: 0.,
                geotechnical: false,
                over_constrained_panels: false,
            },
        )
        .unwrap()
    }

    #[test]
    fn contour_edge_below_the_minimum_edge_becomes_one_vertex() {
        use crate::input::ElementData;
        let policy = Policy {
            closure_tolerance: 0.001,
            junction_movement_limit: 0.05,
            precision: 1e-7,
            minimum_edge: 0.001,
        };
        // Two slab quads joined by a degenerate triangle whose top edge (a
        // contour edge) is 0.5 mm long.
        let mut mesh = MeshData::default();
        for (id, p) in [
            (1, [0., 0.]),
            (2, [1., 0.]),
            (3, [1., 1.]),
            (4, [0., 1.]),
            (5, [1.0005, 1.]),
            (6, [2., 0.]),
            (7, [2., 1.]),
        ] {
            mesh.nodes.insert(id, DVec3::new(p[0], p[1], 0.));
        }
        for (id, nodes) in [
            (1, vec![1, 2, 3, 4]),
            (2, vec![2, 5, 3]),
            (3, vec![2, 6, 7, 5]),
        ] {
            mesh.elements.push(ElementData {
                id,
                elem_type: if nodes.len() == 3 { 42 } else { 44 },
                stiff_id: 1,
                nodes,
            });
        }
        let f = planar_frame(&mesh, DVec3::Z);
        let r = assemble(&mesh, &f, &policy).unwrap();
        assert!(r.all_surface_patches_built, "{:?}", r.issues);
        let pairs: Vec<_> = r
            .stacked_walls
            .coincident
            .iter()
            .map(|c| {
                let mut p = [c.upper_node, c.lower_node];
                p.sort_unstable();
                p
            })
            .collect();
        assert_eq!(pairs, vec![[3, 5]]);
    }

    #[test]
    fn interleaved_fans_touching_at_their_centres_are_split_until_no_pinch() {
        // Topology only (relabelled): two triangle fans around nodes 13 and 14
        // of one property region, every element incident to a pinch node.
        let faces: [&[u32]; 23] = [
            &[16, 17, 13],
            &[14, 18, 15, 19],
            &[20, 14, 12],
            &[12, 14, 15],
            &[21, 18, 14],
            &[11, 13, 1],
            &[13, 6, 5],
            &[13, 22, 6],
            &[11, 1, 14],
            &[1, 5, 14],
            &[14, 5, 6],
            &[14, 6, 22],
            &[13, 17, 10],
            &[13, 10, 3],
            &[13, 3, 9],
            &[13, 9, 4],
            &[13, 4, 8],
            &[13, 8, 21],
            &[13, 21, 2],
            &[13, 2, 7],
            &[13, 7, 22],
            &[21, 14, 2],
            &[14, 22, 7],
        ];
        let facets: BTreeMap<u32, Vec<u32>> = faces
            .iter()
            .enumerate()
            .map(|(i, f)| (i as u32 + 1, f.to_vec()))
            .collect();
        let ids: Vec<u32> = facets.keys().copied().collect();
        let parts = split_pinched_regions(&ids, &facets);
        assert!(parts.len() > 1);
        let mut all: Vec<u32> = parts.iter().flatten().copied().collect();
        all.sort_unstable();
        assert_eq!(all, ids);
        for part in &parts {
            let mut counts = BTreeMap::<[u32; 2], usize>::new();
            for id in part {
                let f = &facets[id];
                for i in 0..f.len() {
                    let (a, b) = (f[i], f[(i + 1) % f.len()]);
                    *counts.entry([a.min(b), a.max(b)]).or_default() += 1;
                }
            }
            let mut degree = BTreeMap::<u32, usize>::new();
            for (e, c) in counts {
                if c == 1 {
                    for n in e {
                        *degree.entry(n).or_default() += 1;
                    }
                }
            }
            assert!(
                degree.values().all(|&d| d == 2),
                "pinched part {part:?}: {degree:?}"
            );
        }
    }

    #[test]
    fn pinched_opening_is_partitioned_without_filling_or_duplicate_sources() {
        use crate::input::ElementData;
        let policy = Policy {
            closure_tolerance: 0.001,
            junction_movement_limit: 0.05,
            precision: 1e-7,
            minimum_edge: 0.001,
        };
        let mut mesh = MeshData::default();
        for y in 0..4 {
            for x in 0..4 {
                mesh.nodes
                    .insert(1 + x + 4 * y, DVec3::new(x as f64, y as f64, 0.));
            }
        }
        for y in 0..3 {
            for x in 0..3 {
                // An interior void touches the exterior at just one vertex.
                if (x == 0 && y == 0) || (x == 1 && y == 1) {
                    continue;
                }
                let a = 1 + x + 4 * y;
                mesh.elements.push(ElementData {
                    id: mesh.elements.len() as u32 + 1,
                    elem_type: 44,
                    stiff_id: 9,
                    nodes: vec![a, a + 1, a + 4, a + 5],
                });
            }
        }
        let f = planar_frame(&mesh, DVec3::Z);
        let ids: Vec<_> = mesh.elements.iter().map(|e| e.id).collect();
        assert_eq!(
            boundary(
                &mesh,
                &mesh.elements.iter().map(|e| (e.id, e)).collect(),
                &ids,
                &f.candidate_planes[0],
                1e-7
            ),
            Err("ambiguous_boundary")
        );
        let base = assemble(&mesh, &f, &policy).unwrap();
        assert!(base.all_surface_patches_built, "{:?}", base.issues);
        assert_eq!(base.pinched_region_splits.len(), 1);
        assert!(base.preview.surfaces.len() > 1);
        let mut sources: Vec<_> = base
            .preview
            .surfaces
            .iter()
            .flat_map(|s| s.source_elements.iter().copied())
            .collect();
        sources.sort_unstable();
        assert_eq!(sources, ids);
        let area: f64 = base
            .preview
            .surfaces
            .iter()
            .map(|s| {
                s.contours
                    .iter()
                    .enumerate()
                    .map(|(j, r)| {
                        let a = (0..r.len())
                            .map(|i| {
                                let (a, b) = (r[i], r[(i + 1) % r.len()]);
                                a[0] * b[1] - a[1] * b[0]
                            })
                            .sum::<f64>()
                            .abs()
                            / 2.;
                        if j == 0 {
                            a
                        } else {
                            -a
                        }
                    })
                    .sum::<f64>()
            })
            .sum();
        assert!((area - 7.).abs() < 1e-9);
        assert_eq!(
            base.vertex_source_nodes.iter().filter(|&&n| n == 6).count(),
            1
        );
        let mut uses = BTreeMap::<usize, usize>::new();
        for e in base
            .preview
            .surfaces
            .iter()
            .flat_map(|s| s.boundaries.iter().flatten())
        {
            *uses.entry(e.edge).or_default() += 1;
        }
        assert!(uses.values().any(|&n| n == 2));
        assert!(base.surface_stiffness.iter().all(|&s| s == 9));
        // Run the whole recognition and assembly after a rigid transform and
        // reversed node/element numbering, not merely the graph helper.
        let q = glam::DQuat::from_axis_angle(DVec3::new(1., 2., 3.).normalize(), 0.6);
        let shift = DVec3::new(30., -40., 70.);
        let mut moved = mesh.clone();
        moved.nodes = mesh
            .nodes
            .iter()
            .map(|(&n, &p)| (100 - n, q * p + shift))
            .collect();
        moved.elements.reverse();
        for e in &mut moved.elements {
            e.id = 100 - e.id;
            for n in &mut e.nodes {
                *n = 100 - *n;
            }
            e.nodes.reverse();
        }
        let transformed = assemble(&moved, &planar_frame(&moved, q * DVec3::Z), &policy).unwrap();
        assert!(transformed.all_surface_patches_built);
        let groups = |r: &Report, reverse: bool| -> BTreeSet<Vec<u32>> {
            r.preview
                .surfaces
                .iter()
                .map(|s| {
                    let mut ids: Vec<_> = s
                        .source_elements
                        .iter()
                        .map(|&n| if reverse { 100 - n } else { n })
                        .collect();
                    ids.sort_unstable();
                    ids
                })
                .collect()
        };
        assert_eq!(groups(&base, false), groups(&transformed, true));
        for (&n, &p) in base.vertex_source_nodes.iter().zip(&base.preview.vertices) {
            let k = transformed
                .vertex_source_nodes
                .iter()
                .position(|&v| v == 100 - n)
                .unwrap();
            assert!(
                DVec3::from_array(transformed.preview.vertices[k])
                    .distance(q * DVec3::from_array(p) + shift)
                    < 1e-7
            );
        }
    }

    fn annulus(triangle: bool) -> MeshData {
        use crate::input::ElementData;
        let (outer, inner) = if triangle {
            (
                vec![[-2., -2., 0.], [2., -2., 0.], [0., 2., 0.]],
                vec![[-0.005, 0., 0.], [0.005, 0., 0.], [0., 0.5, 0.]],
            )
        } else {
            (
                vec![[-2., -2., 0.], [2., -2., 0.], [2., 2., 0.], [-2., 2., 0.]],
                vec![
                    [-0.005, -0.2, 0.],
                    [0.005, -0.2, 0.],
                    [0.005, 0.2, 0.],
                    [-0.005, 0.2, 0.],
                ],
            )
        };
        let n = outer.len() as u32;
        let mut mesh = MeshData::default();
        for (i, p) in outer.into_iter().chain(inner).enumerate() {
            mesh.nodes.insert(i as u32 + 1, DVec3::from_array(p));
        }
        for i in 0..n {
            mesh.elements.push(ElementData {
                id: i + 1,
                elem_type: 44,
                stiff_id: 9,
                nodes: vec![i + 1, (i + 1) % n + 1, (i + 1) % n + n + 1, i + n + 1],
            });
        }
        mesh
    }

    #[test]
    fn holes_restore_across_shapes_scales_rotations_and_renumbering() {
        for triangle in [true, false] {
            for scale in [0.1, 1., 10.] {
                for transformed in [false, true] {
                    let base = annulus(triangle);
                    let count = if triangle { 3 } else { 4 };
                    let rotation = if transformed {
                        glam::DQuat::from_axis_angle(DVec3::new(1., 2., 3.).normalize(), 0.71)
                    } else {
                        glam::DQuat::IDENTITY
                    };
                    let shift = DVec3::new(17., -29., 51.);
                    let number = |n| if transformed { 100 - n } else { n };
                    let mut mesh = base.clone();
                    mesh.nodes = base
                        .nodes
                        .iter()
                        .map(|(&n, &p)| (number(n), rotation * (p * scale) + shift))
                        .collect();
                    for e in &mut mesh.elements {
                        e.id = number(e.id);
                        for n in &mut e.nodes {
                            *n = number(*n);
                        }
                        e.nodes.reverse();
                    }
                    if transformed {
                        mesh.elements.reverse();
                    }
                    let mut f = planar_frame(&mesh, rotation * DVec3::Z);
                    for n in count + 1..=2 * count {
                        let mut p = base.nodes[&n];
                        p.x = 0.;
                        let i = f.node_ids.iter().position(|&i| i == number(n)).unwrap();
                        f.candidate_points[i] = (rotation * (p * scale) + shift).to_array();
                    }
                    let policy = Policy {
                        closure_tolerance: 0.001 * scale,
                        junction_movement_limit: 0.05 * scale,
                        precision: 1e-7 * scale,
                        minimum_edge: 0.001 * scale,
                    };
                    let result = assemble(&mesh, &f, &policy).unwrap();
                    assert!(result.all_surface_patches_built, "{:?}", result.issues);
                    assert_eq!(result.hole_recovery.len(), 1);
                    assert_eq!(result.hole_recovery[0].outcome, HoleOutcome::Restored);
                    assert_eq!(result.preview.surfaces.len(), 1);
                    assert_eq!(result.preview.surfaces[0].contours.len(), 2);
                    assert_eq!(result.surface_stiffness, vec![9]);
                    assert_eq!(
                        result.preview.surfaces[0].source_elements.len(),
                        count as usize
                    );
                    for n in count + 1..=2 * count {
                        let i = result
                            .vertex_source_nodes
                            .iter()
                            .position(|&i| i == number(n))
                            .unwrap();
                        assert!(
                            DVec3::from_array(result.preview.vertices[i])
                                .distance(mesh.nodes[&number(n)])
                                < policy.precision
                        );
                    }
                    assert!(result.maximum_closure_movement <= policy.junction_movement_limit);
                    // Re-running the assembly is deterministic and leaves input intact.
                    let again = assemble(&mesh, &f, &policy).unwrap();
                    assert_eq!(
                        serde_json::to_string(&result).unwrap(),
                        serde_json::to_string(&again).unwrap()
                    );
                }
            }
        }
    }

    #[test]
    fn multiple_holes_in_one_region_restore_together() {
        use crate::input::ElementData;
        let mut mesh = MeshData::default();
        for y in 0..4 {
            for x in 0..6 {
                mesh.nodes.insert(
                    1 + x + 6 * y,
                    DVec3::new(x as f64 * 0.01, y as f64 * 0.01, 0.),
                );
            }
        }
        for y in 0..3 {
            for x in 0..5 {
                if y == 1 && (x == 1 || x == 3) {
                    continue;
                }
                let a = 1 + x + 6 * y;
                mesh.elements.push(ElementData {
                    id: mesh.elements.len() as u32 + 1,
                    elem_type: 44,
                    stiff_id: 9,
                    nodes: vec![a, a + 1, a + 7, a + 6],
                });
            }
        }
        let mut f = planar_frame(&mesh, DVec3::Z);
        for (i, &n) in f.node_ids.iter().enumerate() {
            let p = mesh.nodes[&n];
            if (n >= 8 && n <= 11) || (n >= 14 && n <= 17) {
                f.candidate_points[i][0] = if p.x < 0.025 { 0.015 } else { 0.035 };
            }
        }
        let policy = Policy {
            closure_tolerance: 0.001,
            junction_movement_limit: 0.05,
            precision: 1e-7,
            minimum_edge: 0.001,
        };
        let result = assemble(&mesh, &f, &policy).unwrap();
        assert!(result.all_surface_patches_built, "{:?}", result.issues);
        assert_eq!(result.hole_recovery.len(), 1);
        assert_eq!(result.hole_recovery[0].outcome, HoleOutcome::Restored);
        assert_eq!(result.hole_recovery[0].constraints.len(), 2);
        assert_eq!(result.preview.surfaces[0].contours.len(), 3);
        assert_eq!(result.preview.surfaces[0].source_elements.len(), 13);
    }

    #[test]
    fn blocked_axis_hole_does_not_prevent_independent_recovery() {
        let mut mesh = annulus(true);
        let second = annulus(false);
        mesh.nodes.extend(
            second
                .nodes
                .iter()
                .map(|(&n, &p)| (n + 20, p + DVec3::X * 10.)),
        );
        mesh.elements
            .extend(second.elements.into_iter().map(|mut e| {
                e.id += 20;
                for n in &mut e.nodes {
                    *n += 20;
                }
                e
            }));
        let mut f = planar_frame(&mesh, DVec3::Z);
        for (i, &n) in f.node_ids.iter().enumerate() {
            if (4..=6).contains(&n) {
                f.candidate_points[i][0] = 0.;
            }
            if (25..=28).contains(&n) {
                f.candidate_points[i][0] = 10.;
            }
        }
        let anchor = f.node_ids.iter().position(|&n| n == 4).unwrap();
        let end = f.node_ids.iter().position(|&n| n == 1).unwrap();
        f.axes.push(frame::Axis {
            constructive_segment: false,
            endpoints: [anchor, end],
            anchors: vec![frame::Anchor {
                node: anchor,
                t: 0.,
            }],
            spans: vec![],
        });
        let policy = Policy {
            closure_tolerance: 0.001,
            junction_movement_limit: 0.05,
            precision: 1e-7,
            minimum_edge: 0.001,
        };
        let result = assemble(&mesh, &f, &policy).unwrap();
        assert_eq!(result.hole_recovery.len(), 2);
        let blocked = result
            .hole_recovery
            .iter()
            .find(|r| r.source_elements.contains(&1))
            .unwrap();
        assert_eq!(blocked.outcome, HoleOutcome::AxisAnchorRequiresJointRepair);
        assert!(blocked.changes.is_empty());
        assert!(result
            .hole_recovery
            .iter()
            .any(|r| r.outcome == HoleOutcome::Restored));
        assert_eq!(result.preview.surfaces.len(), 1);
        assert_eq!(
            result.preview.surfaces[0].source_elements,
            vec![21, 22, 23, 24]
        );
        let index = result
            .vertex_source_nodes
            .iter()
            .position(|&n| n == 4)
            .unwrap();
        assert!(
            DVec3::from_array(result.preview.vertices[index])
                .distance(DVec3::from_array(f.candidate_points[anchor]))
                < policy.precision
        );
    }

    #[test]
    fn common_support_line_is_reported_without_moving_any_vertex() {
        let mesh = annulus(true);
        let f = planar_frame(&mesh, DVec3::Z);
        let policy = Policy {
            closure_tolerance: 0.001,
            junction_movement_limit: 0.05,
            precision: 1e-7,
            minimum_edge: 0.001,
        };
        let regions = vec![(0, 9, vec![1, 2, 3])];
        let rings = BTreeMap::from([(0, vec![vec![1, 2, 3], vec![4, 5, 6]])]);
        let supports = vec![
            PlaneFrame::new([0.; 3], [0., 0., 1.]).unwrap(),
            PlaneFrame::new([0.; 3], [1., 0., 0.]).unwrap(),
        ];
        let owners = (1..=6)
            .map(|n| (n, if n >= 4 { vec![0, 1] } else { vec![0] }))
            .collect();
        let mut points: BTreeMap<_, _> = mesh
            .nodes
            .iter()
            .map(|(&n, &p)| (n, if n >= 4 { DVec3::new(0., p.y, 0.) } else { p }))
            .collect();
        let before = points.clone();
        let reports = holes::recover(
            &mut points,
            &holes::Context {
                mesh: &mesh,
                source: &f,
                policy: &policy,
                regions: &regions,
                rings: &rings,
                owners: &owners,
                representatives: &[0, 1],
                supports: &supports,
            },
        );
        assert_eq!(points, before);
        assert_eq!(reports.len(), 1);
        assert_eq!(
            reports[0].outcome,
            HoleOutcome::CommonSupportsForceLineOrPoint
        );
        assert_eq!(reports[0].constraints[0].normal_rank, 2);
        assert_eq!(reports[0].constraints[0].common_supports, vec![0, 1]);
    }

    #[test]
    fn hole_recovery_cannot_degenerate_a_neighbor_and_rolls_back_all_vertices() {
        let mut mesh = annulus(true);
        mesh.nodes.insert(9, DVec3::new(0.1, -0.2, 0.));
        let f = planar_frame(&mesh, DVec3::Z);
        let policy = Policy {
            closure_tolerance: 0.001,
            junction_movement_limit: 0.05,
            precision: 1e-7,
            minimum_edge: 0.001,
        };
        let regions = vec![(0, 9, vec![1, 2, 3]), (0, 9, vec![4])];
        let rings = BTreeMap::from([
            (0, vec![vec![1, 2, 3], vec![4, 5, 6]]),
            (1, vec![vec![4, 1, 9]]),
        ]);
        let supports = vec![PlaneFrame::new([0.; 3], [0., 0., 1.]).unwrap()];
        let owners = mesh.nodes.keys().map(|&n| (n, vec![0])).collect();
        let mut points: BTreeMap<_, _> = mesh
            .nodes
            .iter()
            .map(|(&n, &p)| {
                (
                    n,
                    if (4..=6).contains(&n) {
                        DVec3::new(0., p.y, 0.)
                    } else {
                        p
                    },
                )
            })
            .collect();
        // Neighbor is valid now; restoring vertex 4 would coincide with its
        // other vertex 9. A local hole-only validator would miss this.
        points.insert(9, mesh.nodes[&4]);
        let before = points.clone();
        let reports = holes::recover(
            &mut points,
            &holes::Context {
                mesh: &mesh,
                source: &f,
                policy: &policy,
                regions: &regions,
                rings: &rings,
                owners: &owners,
                representatives: &[0],
                supports: &supports,
            },
        );
        assert_eq!(reports.len(), 1);
        assert_eq!(reports[0].outcome, HoleOutcome::ContourConflict);
        assert!(reports[0].changes.is_empty());
        assert_eq!(points, before);
    }

    #[test]
    fn crack_mouth_is_shared_by_every_region_of_its_patch() {
        let closure = cracks::Closure {
            patch: 7,
            source_elements: 10,
            overlapping_pieces_removed: 0,
            cuts: 1,
            identified: vec![[20, 12]],
            removed_nodes: vec![],
            dropped_contours: 0,
            maximum_width: 0.001,
            contour_nodes_before: 5,
            contour_nodes_after: 4,
        };
        let mut rings = BTreeMap::from([
            (0, vec![vec![10, 11, 12, 13]]),
            (1, vec![vec![30, 20, 31, 32], vec![40, 20, 12, 41]]),
            (2, vec![vec![50, 20, 51]]),
        ]);
        // Regions 0 and 1 belong to patch 7, region 2 to another patch.
        share_crack_mouths(&mut rings, |i| if i == 2 { 3 } else { 7 }, &[closure]);
        assert_eq!(rings[&0], vec![vec![10, 11, 12, 13]]);
        assert_eq!(rings[&1], vec![vec![30, 12, 31, 32], vec![40, 12, 41]]);
        assert_eq!(rings[&2], vec![vec![50, 20, 51]]);
    }

    #[test]
    fn collapsed_hole_is_reported_with_source_nodes_and_never_filled() {
        use crate::input::ElementData;
        let mut mesh = MeshData::default();
        for (i, p) in [
            [-2., -2., 0.],
            [2., -2., 0.],
            [0., 2., 0.],
            [-0.005, 0., 0.],
            [0.005, 0., 0.],
            [0., 0.5, 0.],
        ]
        .into_iter()
        .enumerate()
        {
            mesh.nodes.insert(i as u32 + 1, DVec3::from_array(p));
        }
        for i in 0..3_u32 {
            mesh.elements.push(ElementData {
                id: i + 1,
                elem_type: 44,
                stiff_id: 9,
                nodes: vec![i + 1, (i + 1) % 3 + 1, (i + 1) % 3 + 4, i + 4],
            });
        }
        let policy = Policy {
            closure_tolerance: 0.001,
            junction_movement_limit: 0.05,
            precision: 1e-7,
            minimum_edge: 0.001,
        };
        let mut f = planar_frame(&mesh, DVec3::Z);
        let base = assemble(&mesh, &f, &policy).unwrap();
        assert!(base.all_surface_patches_built, "{:?}", base.issues);
        assert_eq!(base.preview.surfaces[0].contours.len(), 2);
        // A regularized candidate flattens the narrow hole into a line with
        // distinct points, so this is not merely a duplicate-node test.
        for (i, &n) in f.node_ids.iter().enumerate() {
            if n == 4 {
                f.candidate_points[i] = [0., 0., 0.];
            }
            if n == 5 {
                f.candidate_points[i] = [0., 0.1, 0.];
            }
        }
        let result = assemble(&mesh, &f, &policy).unwrap();
        assert!(!result.all_surface_patches_built);
        assert!(result.preview.surfaces.is_empty());
        assert!(result.preview.edges.is_empty());
        assert_eq!(result.issues.len(), 1);
        assert_eq!(result.issues[0].reason, "degenerate_hole_after_closure");
        assert_eq!(result.hole_recovery[0].outcome, HoleOutcome::MovementBudget);
        assert_eq!(result.issues[0].source_elements, vec![1, 2, 3]);
        let repaired =
            assemble_geotechnical(&mesh, &f, &policy, &FeaturePolicy::default()).unwrap();
        assert!(repaired.all_surface_patches_built);
        assert_eq!(repaired.preview.surfaces[0].contours.len(), 1);
        assert_eq!(repaired.preview.surfaces[0].source_elements, vec![1, 2, 3]);
        assert_eq!(repaired.simplified_holes.len(), 1);
        assert_eq!(repaired.preview.vertices, result.preview.vertices);
        let trial = crate::reconstruction::mesh::build(
            &repaired,
            &crate::reconstruction::mesh::Policy {
                boundary_spacing: 0.5,
                maximum_area: 0.5,
                minimum_angle_degrees: 20.,
                maximum_added_vertices_per_surface: 1000,
            },
        )
        .unwrap();
        assert!(trial.topology_valid);
        let used: BTreeSet<_> = trial.triangles.iter().flat_map(|t| t.vertices).collect();
        for n in &repaired.simplified_holes[0].source_nodes {
            assert!(used.contains(
                &repaired
                    .vertex_source_nodes
                    .iter()
                    .position(|id| id == n)
                    .unwrap()
            ));
        }

        assert!(result.issues[0].boundary_source_nodes.iter().any(|r| r
            .iter()
            .copied()
            .collect::<BTreeSet<_>>()
            == BTreeSet::from([4, 5, 6])));
    }

    #[test]
    fn nonmanifold_source_is_not_hidden_by_pinch_partition() {
        let facets = BTreeMap::from([(1, vec![1, 2, 3]), (2, vec![2, 1, 4]), (3, vec![1, 2, 5])]);
        assert_eq!(
            split_pinched_regions(&[1, 2, 3], &facets),
            vec![vec![1, 2, 3]]
        );
    }

    #[test]
    fn offset_projection_closes_three_walls_without_changing_normals() {
        let planes = vec![
            PlaneFrame::new([0.; 3], [1., 0., 0.]).unwrap(),
            PlaneFrame::new([0.; 3], [0., 1., 0.]).unwrap(),
            PlaneFrame::new([0.0003, 0., 0.], [1., 1., 0.]).unwrap(),
        ];
        assert!(intersection(DVec3::ZERO, &planes.iter().collect::<Vec<_>>(), 1e-9).is_none());
        let closed = concurrent_supports(&planes, &[vec![0, 1, 2], vec![2, 1, 0]]);
        let point = intersection(DVec3::ZERO, &closed.iter().collect::<Vec<_>>(), 1e-9).unwrap();
        assert!(point.length() < 0.001);
        for (a, b) in planes.iter().zip(&closed) {
            assert_eq!(a.normal, b.normal);
            assert!(a.distance(b.origin).abs() < 0.001);
        }
        let rotation = glam::DQuat::from_axis_angle(DVec3::new(1., 2., 3.).normalize(), 0.6);
        let shift = DVec3::new(10., 20., 30.);
        let moved: Vec<_> = planes
            .iter()
            .map(|p| {
                PlaneFrame::new(
                    (rotation * DVec3::from_array(p.origin) + shift).to_array(),
                    (rotation * DVec3::from_array(p.normal)).to_array(),
                )
                .unwrap()
            })
            .collect();
        let result = concurrent_supports(&moved, &[vec![2, 0, 1]]);
        let transformed = intersection(shift, &result.iter().collect::<Vec<_>>(), 1e-9).unwrap();
        assert!(transformed.distance(rotation * point + shift) < 1e-9);
    }

    #[test]
    fn tensor_numbering_preserves_outer_boundary_and_hole() {
        use crate::input::ElementData;
        let mut mesh = MeshData::default();
        for y in 0..4 {
            for x in 0..4 {
                mesh.nodes.insert(
                    y * 4 + x + 1,
                    DVec3::new(x as f64 + 10., y as f64 - 20., 7.),
                );
            }
        }
        for y in 0..3 {
            for x in 0..3 {
                if x == 1 && y == 1 {
                    continue;
                }
                let a = y * 4 + x + 1;
                mesh.elements.push(ElementData {
                    id: mesh.elements.len() as u32 + 1,
                    elem_type: 44,
                    stiff_id: 9,
                    nodes: vec![a, a + 1, a + 5, a + 4],
                });
            }
        }
        let plane = PlaneFrame::new([0., 0., 7.], [0., 0., 1.]).unwrap();
        let ids: Vec<_> = mesh.elements.iter().map(|e| e.id).collect();
        let original = boundary(
            &mesh,
            &mesh.elements.iter().map(|e| (e.id, e)).collect(),
            &ids,
            &plane,
            1e-8,
        )
        .unwrap();
        assert_eq!(original.len(), 2);
        assert_eq!(
            original.iter().map(Vec::len).collect::<BTreeSet<_>>(),
            BTreeSet::from([4, 12])
        );
        for e in &mut mesh.elements {
            e.nodes.swap(2, 3);
        }
        assert_eq!(
            boundary(
                &mesh,
                &mesh.elements.iter().map(|e| (e.id, e)).collect(),
                &ids,
                &plane,
                1e-8
            )
            .unwrap(),
            original
        );
        let mut model = Model::new(1e-8, 0.01).unwrap();
        let p = model.add_plane(plane);
        let vertices: BTreeMap<_, _> = mesh
            .nodes
            .iter()
            .map(|(&id, p)| (id, model.add_vertex(p.to_array()).unwrap()))
            .collect();
        let mut rings = original;
        rings.sort_by_key(|r| std::cmp::Reverse(r.len()));
        model
            .add_surface(
                p,
                rings
                    .iter()
                    .map(|r| r.iter().map(|n| vertices[n]).collect())
                    .collect(),
                ids,
            )
            .unwrap();
        assert_eq!(model.surfaces[0].contours.len(), 2);
        assert_eq!(model.surfaces[0].source_elements.len(), 8);
    }

    #[test]
    fn assembles_common_edge_despite_small_unaccepted_frame_residual() {
        use super::super::{planes, recognize};
        use crate::input::ElementData;
        let mut mesh = MeshData::default();
        for (i, p) in [
            [0., 0., 0.],
            [0., 1., 0.],
            [1., 1., 0.],
            [1., 0., 0.],
            [0., 1., 1.],
            [0., 0., 1.],
        ]
        .iter()
        .enumerate()
        {
            mesh.nodes.insert(i as u32 + 1, DVec3::from_array(*p));
        }
        mesh.elements = vec![
            ElementData {
                id: 1,
                elem_type: 44,
                stiff_id: 10,
                nodes: vec![1, 2, 3, 4],
            },
            ElementData {
                id: 2,
                elem_type: 44,
                stiff_id: 20,
                nodes: vec![1, 2, 5, 6],
            },
        ];
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
        let mut f = frame::solve(
            &mesh,
            &axes,
            &planes,
            &frame::Policy {
                up: [0., 0., 1.],
                angle: 0.02,
                maximum_movement: 0.15,
                relative_movement: 0.05,
                minimum_length: 0.03,
                residual_tolerance: 1e-7,
                iterations: 100,
                panel_tolerance: 0.,
                geotechnical: false,
                over_constrained_panels: false,
            },
        )
        .unwrap();
        f.accepted = false;
        f.candidate_points[0] = [0.0003, 0., 0.0002];
        let p = Policy {
            closure_tolerance: 0.001,
            junction_movement_limit: 0.001,
            precision: 1e-7,
            minimum_edge: 0.001,
        };
        let result = assemble(&mesh, &f, &p).unwrap();
        assert!(result.all_surface_patches_built, "{:?}", result.issues);
        assert!(!result.export_ready);
        assert_eq!(result.preview.surfaces.len(), 2);
        assert_eq!(result.preview.vertices.len(), 6);
        assert_eq!(result.preview.edges.len(), 7);
        assert!(result.maximum_closure_movement > 0.0003);
        for surface in &result.preview.surfaces {
            for ring in &surface.boundaries {
                for edge in ring {
                    for v in result.preview.edges[edge.edge] {
                        assert!(
                            result.preview.planes[surface.plane]
                                .distance(result.preview.vertices[v])
                                .abs()
                                < 1e-7
                        );
                    }
                }
            }
        }
        let mut planar = mesh.clone();
        planar.nodes.insert(5, DVec3::new(-1., 1., 0.));
        planar.nodes.insert(6, DVec3::new(-1., 0., 0.));
        let planar_planes = planes::recognize(
            &planar,
            &planes::Policy {
                angle: 0.02,
                distance: 0.01,
                precision: 1e-8,
            },
        )
        .unwrap();
        assert_eq!(planar_planes.patches.len(), 1);
        let planar_frame = frame::solve(&planar, &axes, &planar_planes, &f.policy).unwrap();
        let materialized = assemble(&planar, &planar_frame, &p).unwrap();
        assert!(materialized.all_surface_patches_built);
        assert_eq!(materialized.surface_source_patches, vec![0, 0]);
        assert_eq!(materialized.surface_stiffness, vec![10, 20]);
        assert_eq!(materialized.preview.surfaces.len(), 2);
        assert_eq!(materialized.preview.edges.len(), 7);
        assert_eq!(materialized.preview.surfaces[0].source_elements, vec![1]);
        assert_eq!(materialized.preview.surfaces[1].source_elements, vec![2]);
        f.candidate_points[0] = [0.002, 0., 0.002];
        let blocked = assemble(&mesh, &f, &p).unwrap();
        assert!(!blocked.all_surface_patches_built);
        assert!(blocked.preview.surfaces.is_empty());
        assert_eq!(blocked.issues.len(), 2);
        let mut larger = p.clone();
        larger.junction_movement_limit = 0.01;
        let permitted = assemble(&mesh, &f, &larger).unwrap();
        assert!(permitted.all_surface_patches_built);
        assert!(permitted.maximum_closure_movement > larger.closure_tolerance);
        assert_eq!(permitted.preview.edges.len(), 7);
        let mut limited = f.clone();
        limited.policy.maximum_movement = 1e-6;
        // Move the shared point tangentially so projection does not remove its
        // accumulated displacement from the immutable source model.
        limited.candidate_points[0] = [0.002, 0.0002, 0.002];
        assert!(
            !assemble(&mesh, &limited, &larger)
                .unwrap()
                .all_surface_patches_built
        );
    }
    #[test]
    fn closes_shared_point_on_two_planes_and_rejects_parallel_gap() {
        let a = PlaneFrame::new([0., 0., 0.], [0., 0., 1.]).unwrap();
        let b = PlaneFrame::new([0., 0., 0.], [1., 0., 0.]).unwrap();
        let p = intersection(DVec3::new(0.0003, 2., 0.0002), &[&a, &b], 1e-9).unwrap();
        assert!(p.distance(DVec3::new(0., 2., 0.)) < 1e-9);
        let c = PlaneFrame::new([0., 0., 0.0001], [0., 0., 1.]).unwrap();
        assert!(intersection(p, &[&a, &c], 1e-9).is_none());
    }
    #[test]
    fn redundant_planes_do_not_move_point_twice() {
        let a = PlaneFrame::new([0., 0., 0.], [0., 0., 1.]).unwrap();
        let p = intersection(DVec3::new(1., 2., 0.0003), &[&a, &a], 1e-9).unwrap();
        assert_eq!(p, DVec3::new(1., 2., 0.));
    }

    #[test]
    fn identified_ring_loses_back_and_forth_spikes() {
        let mut ring = vec![234, 940, 2813, 7823, 3267, 2363, 7824, 940];
        collapse_backtracks(&mut ring);
        assert_eq!(ring, vec![2813, 7823, 3267, 2363, 7824, 940]);
        let mut ring = vec![1, 2, 3, 2, 4, 5];
        collapse_backtracks(&mut ring);
        assert_eq!(ring, vec![1, 2, 4, 5]);
        let mut ring = vec![1, 2, 2, 3, 1];
        collapse_backtracks(&mut ring);
        assert_eq!(ring, vec![1, 2, 3]);
    }

    #[test]
    fn support_projection_is_kept_where_it_fits() {
        let plane = |z: f64| PlaneFrame::new([0., 0., z], [0., 0., 1.]).unwrap();
        let candidate = vec![plane(0.), plane(1.)];
        // The first support moves 0.5 mm, the second 5 mm.
        let proposal = vec![plane(0.0005), plane(1.005)];
        let surfaces = vec![
            (0, vec![[0., 0., 0.], [1., 0., 0.]]),
            (1, vec![[0., 0., 1.]]),
        ];
        let (supports, all) = fitting_supports(proposal, &candidate, surfaces.into_iter(), 0.001);
        assert!(!all);
        assert!((supports[0].origin[2] - 0.0005).abs() < 1e-12);
        assert_eq!(supports[1].origin, candidate[1].origin);
    }
}
