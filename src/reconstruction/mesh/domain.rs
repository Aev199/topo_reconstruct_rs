//! Refine only closed material regions; beam edges subdivide regions, not holes.
use super::*;
type Cdt = ConstrainedDelaunayTriangulation<Point2<f64>>;

// Spade's acute-angle refinement can create a long sequence of ever-smaller
// triangles around a fixed constraint vertex. Keep the floor relative to the
// requested area so rigid scaling does not change the refinement behavior.
const MIN_REQUIRED_AREA_RATIO: f64 = 1e-3;

fn interior(
    cdt: &Cdt,
    global: &BTreeMap<usize, usize>,
    boundary: &BTreeSet<[usize; 2]>,
) -> Result<BTreeSet<usize>, &'static str> {
    let mut dual = BTreeMap::<usize, Vec<(usize, bool)>>::new();
    for edge in cdt.undirected_edges() {
        let d = edge.as_directed();
        let a = d.face().fix().index();
        let b = d.rev().face().fix().index();
        let toggle = boundary.contains(&key(
            global
                .get(&d.from().fix().index())
                .copied()
                .unwrap_or(usize::MAX - d.from().fix().index()),
            global
                .get(&d.to().fix().index())
                .copied()
                .unwrap_or(usize::MAX - d.to().fix().index()),
        ));
        dual.entry(a).or_default().push((b, toggle));
        dual.entry(b).or_default().push((a, toggle));
    }
    let outer = cdt.outer_face().fix().index();
    let mut values = BTreeMap::from([(outer, false)]);
    let mut pending = vec![outer];
    while let Some(a) = pending.pop() {
        for &(b, toggle) in dual.get(&a).into_iter().flatten() {
            let value = values[&a] ^ toggle;
            if let Some(previous) = values.get(&b) {
                if *previous != value {
                    return Err("inconsistent mesh boundary graph");
                }
            } else {
                values.insert(b, value);
                pending.push(b);
            }
        }
    }
    Ok(values
        .into_iter()
        .filter_map(|(i, inside)| inside.then_some(i))
        .collect())
}

// Insert size-control points using material boundaries only. Internal beam
// constraints are not holes, and exterior faces must not consume this budget.
fn seed_material(
    region: &mut Cdt,
    mapping: &BTreeMap<usize, usize>,
    boundary: &BTreeSet<[usize; 2]>,
    maximum_area: f64,
    budget: usize,
) -> Result<bool, &'static str> {
    let before = region.num_vertices();
    loop {
        let inside = interior(region, mapping, boundary)?;
        let mut candidates: Vec<_> = region
            .inner_faces()
            .filter(|f| inside.contains(&f.fix().index()))
            .filter_map(|f| {
                let p = f.positions();
                let area = ((p[1].x - p[0].x) * (p[2].y - p[0].y)
                    - (p[1].y - p[0].y) * (p[2].x - p[0].x))
                    .abs()
                    * 0.5;
                (area > maximum_area).then(|| {
                    (
                        area,
                        Point2::new(
                            (p[0].x + p[1].x + p[2].x) / 3.,
                            (p[0].y + p[1].y + p[2].y) / 3.,
                        ),
                    )
                })
            })
            .collect();
        if candidates.is_empty() {
            return Ok(true);
        }
        if region.num_vertices() - before >= budget {
            return Ok(false);
        }
        candidates.sort_by(|a, b| b.0.total_cmp(&a.0));
        let previous = region.num_vertices();
        for (_, p) in candidates.into_iter().take(budget - (previous - before)) {
            region.insert(p).map_err(|_| "invalid material seed")?;
        }
        if previous == region.num_vertices() {
            return Ok(false);
        }
    }
}

// Spade's exclusion treats every constraint as a winding boundary. A dangling
// beam can therefore hide valid material from its angle pass. Revisit only
// material faces, with fixed constraints protected against encroachment.
fn refine_material_angles(
    region: &mut Cdt,
    mapping: &BTreeMap<usize, usize>,
    boundary: &BTreeSet<[usize; 2]>,
    policy: &Policy,
    budget: usize,
) -> Result<bool, &'static str> {
    let before = region.num_vertices();
    loop {
        let inside = interior(region, mapping, boundary)?;
        let fixed: Vec<_> = region
            .undirected_edges()
            .filter(|e| e.is_constraint_edge() || e.is_part_of_convex_hull())
            .map(|e| e.positions())
            .collect();
        let mut candidates = vec![];
        for face in region
            .inner_faces()
            .filter(|f| inside.contains(&f.fix().index()))
        {
            let p = face.positions();
            let area = ((p[1].x - p[0].x) * (p[2].y - p[0].y)
                - (p[1].y - p[0].y) * (p[2].x - p[0].x))
                .abs()
                * 0.5;
            if area < policy.maximum_area * MIN_REQUIRED_AREA_RATIO {
                continue;
            }
            let mut angle = 180.0_f64;
            for i in 0..3 {
                let a = p[(i + 1) % 3];
                let b = p[(i + 2) % 3];
                let o = p[i];
                let (ax, ay, bx, by) = (a.x - o.x, a.y - o.y, b.x - o.x, b.y - o.y);
                angle = angle.min(
                    ((ax * bx + ay * by) / (ax.hypot(ay) * bx.hypot(by)))
                        .clamp(-1., 1.)
                        .acos()
                        .to_degrees(),
                );
            }
            if angle + 1e-7 >= policy.minimum_angle_degrees {
                continue;
            }
            let c = face.circumcenter();
            if !c.x.is_finite() || !c.y.is_finite() {
                continue;
            }
            let in_material = match region.locate(c) {
                spade::PositionInTriangulation::OnFace(f) => inside.contains(&f.index()),
                spade::PositionInTriangulation::OnEdge(e) => {
                    let e = region.directed_edge(e);
                    !e.is_constraint_edge()
                        && [e.face(), e.rev().face()]
                            .iter()
                            .all(|f| inside.contains(&f.fix().index()))
                }
                _ => false,
            };
            if in_material
                && !fixed
                    .iter()
                    .any(|[a, b]| (c.x - a.x) * (c.x - b.x) + (c.y - a.y) * (c.y - b.y) <= 0.)
            {
                candidates.push((angle, c));
            }
        }
        if candidates.is_empty() {
            return Ok(true);
        }
        if region.num_vertices() - before >= budget {
            return Ok(false);
        }
        candidates.sort_by(|a, b| a.0.total_cmp(&b.0));
        let previous = region.num_vertices();
        // Circumcenters become stale after insertion. Recompute before the
        // next point, otherwise neighboring candidates can form tiny edges.
        for (_, p) in candidates.into_iter().take(1) {
            if matches!(
                region.locate(p),
                spade::PositionInTriangulation::OnVertex(_)
            ) {
                continue;
            }
            region
                .insert(p)
                .map_err(|_| "invalid material angle vertex")?;
        }
        if previous == region.num_vertices() {
            return Ok(true);
        }
    }
}

pub(super) fn refine(
    cdt: &Cdt,
    global: &BTreeMap<usize, usize>,
    boundary: &BTreeSet<[usize; 2]>,
    barriers: &BTreeSet<[usize; 2]>,
    constraints: &BTreeSet<[usize; 2]>,
    plane: &super::super::PlaneFrame,
    vertices: &mut Vec<[f64; 3]>,
    policy: &Policy,
) -> Result<(Vec<[usize; 3]>, bool), &'static str> {
    // Reserve size control for the whole material before any local angle
    // refinement can exhaust the per-surface budget in an early region.
    let mut cdt = cdt.clone();
    let mut global = global.clone();
    let initial_count = cdt.num_vertices();
    let size_complete = seed_material(
        &mut cdt,
        &global,
        boundary,
        policy.maximum_area * 0.8,
        policy.maximum_added_vertices_per_surface,
    )?;
    let seeded = cdt.num_vertices() - initial_count;
    for vertex in cdt.vertices() {
        global.entry(vertex.fix().index()).or_insert_with(|| {
            let p = vertex.position();
            let n = vertices.len();
            vertices.push(plane.lift([p.x, p.y]));
            n
        });
    }
    let mut remaining = interior(&cdt, &global, boundary)?;
    let mut neighbors = BTreeMap::<usize, Vec<usize>>::new();
    for edge in cdt.undirected_edges() {
        // Only material boundaries split the refinement domain. Internal
        // construction lines (for example, a beam ending inside a slab) are
        // still constrained edges, but they must remain in one 2D domain so
        // the refinement can build a quality fan around their endpoint.
        let d = edge.as_directed();
        let global_edge = key(
            global[&d.from().fix().index()],
            global[&d.to().fix().index()],
        );
        if barriers.contains(&global_edge) {
            continue;
        }
        let a = d.face().fix().index();
        let b = d.rev().face().fix().index();
        if remaining.contains(&a) && remaining.contains(&b) {
            neighbors.entry(a).or_default().push(b);
            neighbors.entry(b).or_default().push(a);
        }
    }
    let faces: BTreeMap<_, _> = cdt
        .inner_faces()
        .map(|f| {
            (
                f.fix().index(),
                f.vertices().map(|v| global[&v.fix().index()]),
            )
        })
        .collect();
    let uv: BTreeMap<_, _> = cdt
        .vertices()
        .map(|v| (global[&v.fix().index()], v.position()))
        .collect();
    let mut result = Vec::new();
    let mut complete = size_complete;
    let mut added = seeded;
    while let Some(&seed) = remaining.first() {
        remaining.remove(&seed);
        let mut pending = vec![seed];
        let mut counts = BTreeMap::<[usize; 2], usize>::new();
        while let Some(a) = pending.pop() {
            let ids = faces[&a];
            for i in 0..3 {
                *counts.entry(key(ids[i], ids[(i + 1) % 3])).or_default() += 1;
            }
            for &b in neighbors.get(&a).into_iter().flatten() {
                if remaining.remove(&b) {
                    pending.push(b);
                }
            }
        }
        let nodes: BTreeSet<_> = counts.keys().flatten().copied().collect();
        let region_boundary: BTreeSet<_> = counts
            .iter()
            .filter_map(|(&e, &count)| (count == 1).then_some(e))
            .collect();
        let construct_region = || -> Result<(Cdt, BTreeMap<usize, usize>), &'static str> {
            let mut region = Cdt::new();
            let mut handles = BTreeMap::new();
            let mut mapping = BTreeMap::new();
            for &n in &nodes {
                let h = region.insert(uv[&n]).map_err(|_| "invalid region vertex")?;
                handles.insert(n, h);
                mapping.insert(h.index(), n);
            }
            for edge in counts
                .keys()
                .filter(|e| region_boundary.contains(*e) || constraints.contains(*e))
            {
                if !region.can_add_constraint(handles[&edge[0]], handles[&edge[1]]) {
                    return Err("invalid region constraint");
                }
                region.add_constraint(handles[&edge[0]], handles[&edge[1]]);
            }
            Ok((region, mapping))
        };
        let refine_region = || {
            let (mut region, mapping) = construct_region()?;
            let before = region.num_vertices();
            let parameters = RefinementParameters::new()
                .keep_constraint_edges()
                .exclude_outer_faces(true)
                // A small relative floor prevents Spade from endlessly
                // chasing acute fans around source constraints. It is a
                // refinement hint only; the mesh gate still measures every
                // emitted triangle against the configured quality policy.
                .with_min_required_area(policy.maximum_area * MIN_REQUIRED_AREA_RATIO)
                .with_max_allowed_area(policy.maximum_area)
                .with_angle_limit(AngleLimit::from_deg(policy.minimum_angle_degrees))
                .with_max_additional_vertices(
                    policy
                        .maximum_added_vertices_per_surface
                        .saturating_sub(added),
                );
            // Spade currently has an internal panic path in `refine` for some
            // valid-looking constrained configurations (it reports "Failed to
            // locate position"). Mesh generation is a diagnostic gate, so an
            // upstream triangulator panic must become an ordinary rejection.
            let refined = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                region.refine(parameters)
            }))
            .map_err(|_| "CDT refinement panicked")?;
            Ok((region, mapping, before, refined))
        };
        let (mut region, mut mapping, before, refined) = refine_region()?;
        let angle_budget = policy
            .maximum_added_vertices_per_surface
            .saturating_sub(added + region.num_vertices() - before);
        let angle_complete = refine_material_angles(
            &mut region,
            &mapping,
            &region_boundary,
            policy,
            angle_budget,
        )?;
        added += region.num_vertices() - before;
        complete &= refined.refinement_complete && angle_complete;
        for vertex in region.vertices() {
            mapping.entry(vertex.fix().index()).or_insert_with(|| {
                let p = vertex.position();
                let n = vertices.len();
                vertices.push(plane.lift([p.x, p.y]));
                n
            });
        }
        // Verify material membership with our boundary-only graph, independently
        // of the library's refinement exclusion (dangling beam edges are allowed).
        let inside = interior(&region, &mapping, &region_boundary)?;
        for face in region.inner_faces() {
            if inside.contains(&face.fix().index()) {
                result.push(face.vertices().map(|v| mapping[&v.fix().index()]));
            }
        }
    }
    Ok((result, complete))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn material_seeding_fills_both_sides_of_open_beam_but_not_hole() {
        for scale in [0.1, 1., 10.] {
            let points = [
                [0., 0.],
                [10., 0.],
                [10., 10.],
                [0., 10.],
                [0., 2.],
                [4., 4.],
                [6., 4.],
                [6., 6.],
                [4., 6.],
                [3., 2.],
            ];
            let mut region = Cdt::new();
            let mut mapping = BTreeMap::new();
            let mut handles = vec![];
            for (n, p) in points.into_iter().enumerate() {
                let h = region
                    .insert(Point2::new(p[0] * scale, p[1] * scale))
                    .unwrap();
                mapping.insert(h.index(), n);
                handles.push(h);
            }
            let mut boundary = BTreeSet::new();
            for ring in [vec![0, 1, 2, 3, 4], vec![5, 6, 7, 8]] {
                for i in 0..ring.len() {
                    let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
                    region.add_constraint(handles[a], handles[b]);
                    boundary.insert(key(a, b));
                }
            }
            region.add_constraint(handles[4], handles[9]);
            let max_area = 0.5 * scale * scale;
            assert!(!seed_material(&mut region, &mapping, &boundary, max_area, 0).unwrap());
            assert_eq!(region.num_vertices(), points.len());
            assert!(seed_material(&mut region, &mapping, &boundary, max_area, 1000).unwrap());
            assert!(region.num_vertices() > points.len());
            for v in region
                .vertices()
                .filter(|v| !mapping.contains_key(&v.fix().index()))
            {
                let p = v.position();
                assert!(p.x > 0. && p.x < 10. * scale && p.y > 0. && p.y < 10. * scale);
                assert!(
                    !(p.x > 4. * scale && p.x < 6. * scale && p.y > 4. * scale && p.y < 6. * scale)
                );
            }
            let inside = interior(&region, &mapping, &boundary).unwrap();
            let mut area = 0.;
            for f in region
                .inner_faces()
                .filter(|f| inside.contains(&f.fix().index()))
            {
                let p = f.positions();
                let a = ((p[1].x - p[0].x) * (p[2].y - p[0].y)
                    - (p[1].y - p[0].y) * (p[2].x - p[0].x))
                    .abs()
                    / 2.;
                assert!(a <= max_area * (1. + 1e-10));
                area += a;
            }
            assert!((area - 96. * scale * scale).abs() < 1e-8 * scale * scale);
            assert_eq!(region.num_constraints(), 10);
        }
    }
}
