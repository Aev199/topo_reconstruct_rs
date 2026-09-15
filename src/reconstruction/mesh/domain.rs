//! Refine only closed material regions; beam edges subdivide regions, not holes.
use super::*;
type Cdt = ConstrainedDelaunayTriangulation<Point2<f64>>;

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
            global[&d.from().fix().index()],
            global[&d.to().fix().index()],
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

pub(super) fn refine(
    cdt: &Cdt,
    global: &BTreeMap<usize, usize>,
    boundary: &BTreeSet<[usize; 2]>,
    barriers: &BTreeSet<[usize; 2]>,
    constraints: &BTreeSet<[usize; 2]>,
    refine_outer_faces: bool,
    plane: &super::super::PlaneFrame,
    vertices: &mut Vec<[f64; 3]>,
    policy: &Policy,
) -> Result<(Vec<[usize; 3]>, bool), &'static str> {
    let mut remaining = interior(cdt, global, boundary)?;
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
    let mut complete = true;
    let mut added = 0;
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
        let mut region = Cdt::new();
        let mut handles = BTreeMap::new();
        let mut mapping = BTreeMap::new();
        for n in nodes {
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
        let before = region.num_vertices();
        let mut parameters = RefinementParameters::new()
            .keep_constraint_edges()
            .with_max_allowed_area(policy.maximum_area)
            .with_angle_limit(AngleLimit::from_deg(policy.minimum_angle_degrees))
            .with_max_additional_vertices(
                policy
                    .maximum_added_vertices_per_surface
                    .saturating_sub(added),
            );
        if !refine_outer_faces {
            parameters = parameters.exclude_outer_faces(true);
        }
        // Spade currently has an internal panic path in `refine` for some
        // valid-looking constrained configurations (it reports "Failed to
        // locate position"). Mesh generation is a diagnostic gate, so an
        // upstream triangulator panic must become an ordinary rejection of
        // this mesh candidate rather than aborting the whole reconstruction.
        let refined =
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| region.refine(parameters)))
                .map_err(|_| "CDT refinement panicked")?;
        added += region.num_vertices() - before;
        complete &= refined.refinement_complete;
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

