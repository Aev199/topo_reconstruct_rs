//! One transactional pass over movable, non-branching boundary anchors.
use super::*;
use geo::{Area, BooleanOps};

fn no_new_coplanar_overlap(
    before: &Model,
    after: &Model,
    changed: &BTreeSet<usize>,
    precision: f64,
) -> bool {
    for i in 0..before.surfaces.len() {
        for j in i + 1..before.surfaces.len() {
            if !changed.contains(&i) && !changed.contains(&j) {
                continue;
            }
            let a = &before.surfaces[i];
            let b = &before.surfaces[j];
            let pa = &before.planes[a.plane];
            let pb = &before.planes[b.plane];
            if DVec3::from_array(pa.normal)
                .cross(DVec3::from_array(pb.normal))
                .length()
                > 1e-10
                || b.contours
                    .iter()
                    .flatten()
                    .any(|&uv| pa.distance(pb.lift(uv)).abs() > precision)
            {
                continue;
            }
            let projected = |model: &Model, s: usize| {
                let surface = &model.surfaces[s];
                let plane = &model.planes[surface.plane];
                polygon(
                    &surface
                        .contours
                        .iter()
                        .map(|r| r.iter().map(|&uv| pa.project(plane.lift(uv))).collect())
                        .collect::<Vec<_>>(),
                )
            };
            let old_a = polygon(&a.contours);
            let old_b = projected(before, j);
            let new_a = polygon(&after.surfaces[i].contours);
            let new_b = projected(after, j);
            let tolerance =
                precision * (old_a.unsigned_area().sqrt() + old_b.unsigned_area().sqrt());
            let old_overlap = old_a.intersection(&old_b);
            let new_overlap = new_a.intersection(&new_b);
            if new_overlap.difference(&old_overlap).unsigned_area() > tolerance {
                return false;
            }
        }
    }
    true
}

#[derive(Debug, Serialize)]
pub struct Change {
    pub source_node: u32,
    pub before: [f64; 3],
    pub after: [f64; 3],
}
#[derive(Debug, Serialize)]
pub struct Repair {
    pub source_axes: Vec<usize>,
    pub surfaces: Vec<usize>,
    pub accepted: bool,
    pub reason: String,
    pub changes: Vec<Change>,
}
pub(super) struct Context<'a> {
    pub mesh: &'a MeshData,
    pub source: &'a frame::Report,
    pub vertices: &'a BTreeMap<u32, usize>,
    pub incidence: &'a BTreeMap<u32, BTreeSet<usize>>,
    pub unavailable: &'a BTreeMap<u32, &'a str>,
    pub owner_planes: &'a BTreeMap<u32, Vec<usize>>,
    pub owner_surfaces: &'a BTreeMap<u32, BTreeSet<usize>>,
    pub supports: &'a [PlaneFrame],
    pub policy: &'a Policy,
}

fn validate_axes(
    c: &Context<'_>,
    model: &Model,
    locked: &BTreeMap<u32, DVec3>,
    axes: &BTreeSet<usize>,
) -> bool {
    axes.iter().all(|&i| {
        let Ok(p) = propose(
            c.mesh,
            c.source,
            &c.source.axes[i],
            locked,
            c.unavailable,
            c.owner_planes,
            c.supports,
            c.policy,
        ) else {
            return false;
        };
        p.points.iter().all(|(n, p)| {
            c.owner_surfaces.get(n).into_iter().flatten().all(|&s| {
                let surface = &model.surfaces[s];
                let plane = &model.planes[surface.plane];
                plane.distance(p.to_array()).abs() <= c.policy.precision
                    && location(
                        plane.project(p.to_array()),
                        &surface.contours,
                        c.policy.precision,
                    )
                    .is_some()
            })
        })
    })
}

pub(super) fn repair(
    model: &mut Model,
    locked: &mut BTreeMap<u32, DVec3>,
    c: &Context<'_>,
) -> Vec<Repair> {
    // No movement of hole boundaries, axis ends, or multi-axis junctions in
    // this rule. They need a genuinely coupled multi-axis solve.
    let mut protected = BTreeSet::new();
    for a in &c.source.axes {
        protected.extend(a.endpoints.map(|i| c.source.node_ids[i]));
    }
    let by_vertex: BTreeMap<_, _> = c.vertices.iter().map(|(&n, &v)| (v, n)).collect();
    for surface in &model.surfaces {
        for edge in surface.boundaries.iter().skip(1).flatten() {
            protected.extend(model.edges[edge.edge].map(|v| by_vertex[&v]));
        }
    }
    let baseline: BTreeSet<_> = (0..c.source.axes.len())
        .filter(|&i| validate_axes(c, model, locked, &BTreeSet::from([i])))
        .collect();
    let mut proposals = vec![];
    for (i, axis) in c.source.axes.iter().enumerate() {
        if baseline.contains(&i) {
            continue;
        }
        let movable: BTreeSet<_> = axis
            .anchors
            .iter()
            .map(|a| c.source.node_ids[a.node])
            .filter(|n| {
                c.vertices.contains_key(n) && !protected.contains(n) && c.incidence[n].len() == 1
            })
            .collect();
        if movable.is_empty() {
            continue;
        }
        let mut relaxed = locked.clone();
        for n in &movable {
            relaxed.remove(n);
        }
        let Ok(p) = propose(
            c.mesh,
            c.source,
            axis,
            &relaxed,
            c.unavailable,
            c.owner_planes,
            c.supports,
            c.policy,
        ) else {
            continue;
        };
        let changed: BTreeMap<_, _> = movable
            .into_iter()
            .filter_map(|n| {
                let q = p.points[&n];
                (q.distance(locked[&n]) > c.policy.precision).then_some((n, q))
            })
            .collect();
        if changed.is_empty() {
            continue;
        }
        let affected: BTreeSet<_> = model
            .surfaces
            .iter()
            .enumerate()
            .filter_map(|(s, surface)| {
                surface
                    .boundaries
                    .iter()
                    .flatten()
                    .any(|e| {
                        model.edges[e.edge]
                            .iter()
                            .any(|v| changed.contains_key(&by_vertex[v]))
                    })
                    .then_some(s)
            })
            .collect();
        proposals.push((i, changed, affected));
    }
    let mut reports = vec![];
    let mut remaining: BTreeSet<_> = (0..proposals.len()).collect();
    while let Some(&seed) = remaining.first() {
        remaining.remove(&seed);
        let mut group = vec![seed];
        let mut surfaces = proposals[seed].2.clone();
        loop {
            let connected: Vec<_> = remaining
                .iter()
                .copied()
                .filter(|&i| !surfaces.is_disjoint(&proposals[i].2))
                .collect();
            if connected.is_empty() {
                break;
            }
            for i in connected {
                remaining.remove(&i);
                surfaces.extend(&proposals[i].2);
                group.push(i);
            }
        }
        let mut trial = model.clone();
        let mut proposed_locked = locked.clone();
        let mut changes = vec![];
        for &i in &group {
            for (&n, &q) in &proposals[i].1 {
                changes.push(Change {
                    source_node: n,
                    before: locked[&n].to_array(),
                    after: q.to_array(),
                });
                trial.vertices[c.vertices[&n]] = q.to_array();
                proposed_locked.insert(n, q);
            }
        }
        let valid = surfaces.iter().all(|&s| {
            let surface = &model.surfaces[s];
            let rings: Vec<Vec<usize>> = surface
                .boundaries
                .iter()
                .map(|ring| {
                    ring.iter()
                        .map(|e| {
                            let ends = model.edges[e.edge];
                            if e.reversed {
                                ends[1]
                            } else {
                                ends[0]
                            }
                        })
                        .collect()
                })
                .collect();
            let mut check = Model::new(c.policy.precision, c.policy.minimum_edge).unwrap();
            check.vertices = trial.vertices.clone();
            let plane = check.add_plane(model.planes[surface.plane].clone());
            if check
                .add_surface(plane, rings, surface.source_elements.clone())
                .is_err()
            {
                return false;
            }
            trial.surfaces[s].contours = check.surfaces[0].contours.clone();
            true
        });
        let axes: BTreeSet<_> = group.iter().map(|&i| proposals[i].0).collect();
        let mut required = baseline.clone();
        required.extend(&axes);
        // Also retain axes repaired in an earlier, independent component.
        for r in &reports {
            let r: &Repair = r;
            if r.accepted {
                required.extend(&r.source_axes);
            }
        }
        let no_overlap =
            valid && no_new_coplanar_overlap(model, &trial, &surfaces, c.policy.precision);
        let accepted = no_overlap && validate_axes(c, &trial, &proposed_locked, &required);
        reports.push(Repair {
            source_axes: axes.into_iter().collect(),
            surfaces: surfaces.into_iter().collect(),
            accepted,
            reason: if !valid {
                "invalid_neighbor_contour"
            } else if !no_overlap {
                "new_coplanar_overlap"
            } else if !accepted {
                "axis_or_incidence_conflict"
            } else {
                "boundary_anchors_regularized"
            }
            .into(),
            changes: if accepted { changes } else { vec![] },
        });
        if accepted {
            *model = trial;
            *locked = proposed_locked;
        }
    }
    reports
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn rejects_new_overlap_between_individually_valid_surfaces() {
        let mut before = Model::new(1e-7, 0.001).unwrap();
        let p = before.add_plane(PlaneFrame::new([0.; 3], [0., 0., 1.]).unwrap());
        for (left, right) in [(0., 1.), (1.01, 2.)] {
            let ring = [
                [left, 0., 0.],
                [right, 0., 0.],
                [right, 1., 0.],
                [left, 1., 0.],
            ]
            .map(|q| before.add_vertex(q).unwrap())
            .to_vec();
            before.add_surface(p, vec![ring], vec![]).unwrap();
        }
        let mut after = before.clone();
        for v in [1, 2] {
            after.vertices[v][0] += 0.02;
        }
        after.surfaces[0].contours[0] = [0, 1, 2, 3]
            .map(|v| after.planes[p].project(after.vertices[v]))
            .to_vec();
        assert!(!no_new_coplanar_overlap(
            &before,
            &after,
            &BTreeSet::from([0]),
            1e-7
        ));
        assert!(no_new_coplanar_overlap(
            &before,
            &before,
            &BTreeSet::from([0]),
            1e-7
        ));
    }
}
