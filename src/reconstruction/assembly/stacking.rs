//! Alignment of stacked walls to the axis of the wall carrying them.
//!
//! Walls above and below one slab often have aligned outer faces but
//! different thicknesses, so their mid-planes are a few centimetres apart.
//! Two junction lines millimetres apart on the slab then force tiny
//! elements. A parallel wall standing on the opposite side of the same slab,
//! overlapping the lower wall along the junction and offset by at most the
//! stacking tolerance, adopts the support plane of the lowest wall of its
//! stack. Vertices then close onto the new plane through the ordinary
//! support intersection with its movement limits; a pair whose vertices would
//! exceed them is kept and reported.
use super::{frame, intersection, movement_budget, MeshData, Policy};
use crate::reconstruction::PlaneFrame;
use glam::DVec3;
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Serialize)]
pub struct StackedWall {
    /// Source plane patch of the moved (upper) wall.
    pub upper: usize,
    /// Patch directly below it, across `slab`.
    pub lower: usize,
    /// Lowest patch of the stack, whose plane is adopted.
    pub root: usize,
    /// Slab between a stacked pair; none for walls in one line.
    pub slab: Option<usize>,
    /// Distance of the upper wall to the adopted plane.
    pub offset: f64,
    pub source_elements: Vec<u32>,
    pub reason: String,
}

#[derive(Debug, Default, Serialize)]
pub struct Report {
    pub tolerance: f64,
    pub aligned: Vec<StackedWall>,
    pub kept: Vec<StackedWall>,
    pub identified: Vec<Identified>,
    /// Source nodes on the same supports (two or more planes) closed within
    /// the minimum edge length of each other: one vertex, the second node
    /// takes the first node's vertex (`upper_node` takes `lower_node`).
    pub coincident: Vec<Identified>,
    /// Source nodes at the junction of two or more aligned structures, whose
    /// closure movement limit is raised (see `movement_limit`).
    pub raised_limits: Vec<RaisedLimit>,
    /// Patches whose plane was replaced by an alignment.
    #[serde(skip)]
    pub shifted: BTreeSet<usize>,
}

#[derive(Debug, Clone, Serialize)]
pub struct RaisedLimit {
    pub source_node: u32,
    /// Distinct aligned structures the node lies on.
    pub aligned_supports: usize,
    pub limit: f64,
}

/// Closure movement limit of a node: each alignment moves a structure by at
/// most the tolerance, so a node where `k` independently aligned structures
/// meet moves by up to their vector sum, bounded by sqrt(k) times the limit
/// (user decision: such corners are closed).
pub(super) fn movement_limit(
    limit: f64,
    supports: &[usize],
    shifted: impl Fn(usize) -> bool,
    group: impl Fn(usize) -> usize,
) -> (f64, usize) {
    let k = supports
        .iter()
        .filter(|&&s| shifted(s))
        .map(|&s| group(s))
        .collect::<BTreeSet<_>>()
        .len();
    (limit * (k.max(1) as f64).sqrt(), k)
}

struct Link {
    lower: usize,
    slab: usize,
    overlap: f64,
    offset: f64,
    /// Contact nodes of the upper and lower wall on the slab.
    contact: [Vec<usize>; 2],
    direction: DVec3,
}

/// An upper-wall node that takes the vertex of a lower-wall node: after the
/// alignment both lie on one junction line, less than the minimum edge
/// length apart. Their offset was created by the alignment, not the source.
#[derive(Debug, Clone, Serialize)]
pub struct Identified {
    pub upper_node: u32,
    pub lower_node: u32,
    /// Distance along the junction line.
    pub distance: f64,
}

fn normal(p: &PlaneFrame) -> DVec3 {
    DVec3::from_array(p.normal)
}

/// Largest distance of candidate points to a plane. Offsets are measured at
/// the contact (the junction on the slab, the touching wall ends): over a
/// whole wall a slight non-parallelism adds to it, and that is bounded by the
/// vertex movement limits instead.
fn offset(source: &frame::Report, nodes: &[usize], plane: &PlaneFrame) -> f64 {
    nodes
        .iter()
        .map(|&n| plane.distance(source.candidate_points[n]).abs())
        .fold(0.0_f64, f64::max)
}

/// Whether every source vertex of the patches represented by `moved` closes
/// within its limits once they adopt the plane of `target`, and no contour
/// edge through a moved vertex collapses below the minimum edge length (an
/// alignment must not destroy a surface).
fn fits(
    mesh: &MeshData,
    source: &frame::Report,
    policy: &Policy,
    owners: &BTreeMap<u32, Vec<usize>>,
    rings: &BTreeMap<usize, Vec<Vec<u32>>>,
    representatives: &[usize],
    shifted: &BTreeSet<usize>,
    moved: usize,
    target: usize,
) -> Result<(), &'static str> {
    let planes = &source.candidate_planes;
    let lookup: BTreeMap<u32, usize> = source
        .node_ids
        .iter()
        .enumerate()
        .map(|(i, &n)| (n, i))
        .collect();
    let remap = |s: usize| {
        if representatives[s] == moved {
            target
        } else {
            representatives[s]
        }
    };
    let mut closed = BTreeMap::<u32, DVec3>::new();
    for (&id, supports) in owners {
        if !supports.iter().any(|&s| representatives[s] == moved) {
            continue;
        }
        let Some(&i) = lookup.get(&id) else {
            return Err("vertex_movement_limit");
        };
        let p = DVec3::from_array(source.candidate_points[i]);
        let chosen: BTreeSet<usize> = supports.iter().map(|&s| remap(s)).collect();
        let refs: Vec<&PlaneFrame> = chosen.iter().map(|&s| &planes[s]).collect();
        let Some(q) = intersection(p, &refs, policy.precision) else {
            return Err("inconsistent_supports");
        };
        let (limit, _) = movement_limit(
            policy.junction_movement_limit,
            supports,
            |s| shifted.contains(&s) || representatives[s] == moved,
            remap,
        );
        if p.distance(q) > limit
            || q.distance(mesh.nodes[&id]) > movement_budget(mesh, source, i) + policy.precision
        {
            return Err("vertex_movement_limit");
        }
        closed.insert(id, q);
    }
    let at = |id: u32| {
        closed.get(&id).copied().or_else(|| {
            lookup
                .get(&id)
                .map(|&i| DVec3::from_array(source.candidate_points[i]))
        })
    };
    // Nodes on the same supports that meet are identified (one vertex).
    let key = |n: u32| -> Option<BTreeSet<usize>> {
        let set: BTreeSet<usize> = owners.get(&n)?.iter().map(|&s| remap(s)).collect();
        (set.len() >= 2).then_some(set)
    };
    let intact = rings.values().flatten().all(|ring| {
        !ring.iter().any(|n| closed.contains_key(n))
            || (0..ring.len()).all(|k| {
                let (a, b) = (ring[k], ring[(k + 1) % ring.len()]);
                if key(a).is_some_and(|x| Some(x) == key(b)) {
                    return true;
                }
                let (Some(pa), Some(pb)) = (at(a), at(b)) else {
                    return true;
                };
                let before = |n: u32| {
                    lookup
                        .get(&n)
                        .map(|&i| DVec3::from_array(source.candidate_points[i]))
                };
                let old = before(a).zip(before(b)).map_or(0., |(x, y)| x.distance(y));
                a == b || old < policy.minimum_edge || pa.distance(pb) >= policy.minimum_edge
            })
    });
    if intact {
        Ok(())
    } else {
        Err("contour_edge_collapse")
    }
}

pub(super) fn align(
    mesh: &MeshData,
    source: &frame::Report,
    policy: &Policy,
    owners: &BTreeMap<u32, Vec<usize>>,
    rings: &BTreeMap<usize, Vec<Vec<u32>>>,
    representatives: &mut [usize],
    tolerance: f64,
) -> Report {
    let mut report = Report {
        tolerance,
        ..Default::default()
    };
    let up = DVec3::from_array(source.policy.up).normalize();
    let planes = &source.candidate_planes;
    let count = source.surfaces.len();
    let nodes: Vec<BTreeSet<usize>> = source
        .surfaces
        .iter()
        .map(|s| s.nodes.iter().copied().collect())
        .collect();
    let (sin, cos) = source.policy.angle.sin_cos();
    let wall = |i: usize| normal(&planes[i]).dot(up).abs() <= sin;
    let slab = |i: usize| normal(&planes[i]).dot(up).abs() >= cos;
    let point = |n: usize| DVec3::from_array(source.candidate_points[n]);
    let mut links = BTreeMap::<usize, Link>::new();
    for s in (0..count).filter(|&s| slab(s)) {
        let ns = normal(&planes[s]) * normal(&planes[s]).dot(up).signum();
        let slab_points: Vec<DVec3> = source.surfaces[s].nodes.iter().map(|&n| point(n)).collect();
        // Geometric contact: a wall resting on a slab need not share source
        // nodes with it. A node counts if it lies on the slab plane and
        // within three median node spacings of a slab node.
        let nearest = |p: DVec3, skip: Option<usize>| {
            slab_points
                .iter()
                .enumerate()
                .filter(|&(k, _)| Some(k) != skip)
                .map(|(_, q)| p.distance(*q))
                .fold(f64::INFINITY, f64::min)
        };
        let step = slab_points.len().div_ceil(64).max(1);
        let mut spacings: Vec<f64> = (0..slab_points.len())
            .step_by(step)
            .map(|k| nearest(slab_points[k], Some(k)))
            .filter(|d| d.is_finite())
            .collect();
        if spacings.is_empty() {
            continue;
        }
        spacings.sort_by(f64::total_cmp);
        let reach = 3. * spacings[spacings.len() / 2];
        // Walls in line contact with the slab, with their side and extent.
        let mut touching = vec![];
        for w in (0..count).filter(|&w| wall(w)) {
            let (contact, rest): (Vec<usize>, Vec<usize>) = nodes[w].iter().partition(|&&n| {
                nodes[s].contains(&n)
                    || (planes[s].distance(point(n).to_array()).abs() <= policy.closure_tolerance
                        && nearest(point(n), None) <= reach)
            });
            if contact.len() < 2 || rest.is_empty() {
                continue;
            }
            let side = rest
                .iter()
                .map(|&n| planes[s].distance(point(n).to_array()) * ns.dot(normal(&planes[s])))
                .sum::<f64>()
                / rest.len() as f64;
            let t = normal(&planes[w]).cross(up).normalize();
            let (lo, hi) =
                contact
                    .iter()
                    .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &n| {
                        let x = point(n).dot(t);
                        (lo.min(x), hi.max(x))
                    });
            touching.push((w, side, t, lo, hi, contact));
        }
        for (u, side_u, t, lo_u, hi_u, contact_u) in &touching {
            let (u, side_u, t, lo_u, hi_u) = (*u, *side_u, *t, *lo_u, *hi_u);
            if side_u <= 0. {
                continue;
            }
            for (l, side_l, _, _, _, contact_l) in &touching {
                let (l, side_l) = (*l, *side_l);
                if side_l >= 0. || representatives[u] == representatives[l] {
                    continue;
                }
                if normal(&planes[u]).dot(normal(&planes[l])).abs() < cos {
                    continue;
                }
                let d = offset(source, contact_u, &planes[l]);
                if d <= policy.closure_tolerance || d > tolerance + policy.closure_tolerance {
                    continue;
                }
                // Overlap along the junction, in the upper wall's direction.
                let (lo_l, hi_l) =
                    contact_l
                        .iter()
                        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &n| {
                            let x = point(n).dot(t);
                            (lo.min(x), hi.max(x))
                        });
                let overlap = hi_u.min(hi_l) - lo_u.max(lo_l);
                if overlap <= policy.minimum_edge {
                    continue;
                }
                let better = links
                    .get(&u)
                    .is_none_or(|k| overlap > k.overlap || (overlap == k.overlap && d < k.offset));
                if better {
                    links.insert(
                        u,
                        Link {
                            lower: l,
                            slab: s,
                            overlap,
                            offset: d,
                            contact: [contact_u.clone(), contact_l.clone()],
                            direction: t,
                        },
                    );
                }
            }
        }
    }
    // Lowest walls first, so every root keeps its own plane.
    let height =
        |i: usize| nodes[i].iter().map(|&n| point(n).dot(up)).sum::<f64>() / nodes[i].len() as f64;
    let anchors: BTreeSet<usize> = source
        .axes
        .iter()
        .flat_map(|a| a.anchors.iter().map(|x| x.node))
        .collect();
    let mut uppers: Vec<usize> = links.keys().copied().collect();
    uppers.sort_by(|&a, &b| height(a).total_cmp(&height(b)).then(a.cmp(&b)));
    for u in uppers {
        let link = &links[&u];
        let mut root = link.lower;
        let mut guard = BTreeSet::from([u]);
        while let Some(next) = links.get(&root) {
            if !guard.insert(root) {
                break;
            }
            root = next.lower;
        }
        let target = representatives[root];
        let moved = representatives[u];
        let mut entry = StackedWall {
            upper: u,
            lower: link.lower,
            root,
            slab: Some(link.slab),
            offset: offset(source, &link.contact[0], &planes[target]),
            source_elements: source.surfaces[u].source_elements.clone(),
            reason: String::new(),
        };
        if entry.offset > tolerance + policy.closure_tolerance {
            entry.reason = "stack_offset_exceeds_tolerance".into();
            report.kept.push(entry);
            continue;
        }
        // Every vertex of the moved support must close within its limits.
        let fits = fits(
            mesh,
            source,
            policy,
            owners,
            rings,
            representatives,
            &report.shifted,
            moved,
            target,
        );
        if let Err(reason) = fits {
            entry.reason = reason.into();
            report.kept.push(entry);
            continue;
        }
        for (j, r) in representatives.iter_mut().enumerate() {
            if *r == moved {
                *r = target;
                report.shifted.insert(j);
            }
        }
        // One-to-one identification of nearly coincident junction nodes.
        let boundary = |n: usize| owners.contains_key(&source.node_ids[n]) && !anchors.contains(&n);
        let mut pairs: Vec<(f64, usize, usize)> = vec![];
        for &a in link.contact[0].iter().filter(|&&n| boundary(n)) {
            for &b in link.contact[1].iter().filter(|&&n| boundary(n)) {
                if a == b || (nodes[link.slab].contains(&a) && nodes[link.slab].contains(&b)) {
                    continue;
                }
                let d = (point(a) - point(b)).dot(link.direction).abs();
                if d < policy.minimum_edge {
                    pairs.push((d, a, b));
                }
            }
        }
        pairs.sort_by(|x, y| x.0.total_cmp(&y.0));
        let mut used = BTreeSet::new();
        for (d, a, b) in pairs {
            if used.insert(a) && used.insert(b) {
                report.identified.push(Identified {
                    upper_node: source.node_ids[a],
                    lower_node: source.node_ids[b],
                    distance: d,
                });
            }
        }
        entry.reason = "aligned_to_lowest_wall".into();
        report.aligned.push(entry);
    }
    report
}

/// Walls in one line: a parallel wall continuing another one end to end (a
/// wall line jogged by a few centimetres, a low riser skewed off the wall
/// it continues) adopts the plane of the larger wall, like a stacked wall.
/// The two must touch (source nodes within `tolerance`), be offset at the
/// contact by at most `tolerance`, and must not overlap in the common
/// plane: parallel walls side by side (a double wall, a joint) are never
/// merged. Slabs are never aligned: slabs at different levels stay apart.
pub(super) fn align_lines(
    mesh: &MeshData,
    source: &frame::Report,
    policy: &Policy,
    owners: &BTreeMap<u32, Vec<usize>>,
    rings: &BTreeMap<usize, Vec<Vec<u32>>>,
    representatives: &mut [usize],
    tolerance: f64,
    report: &mut Report,
) {
    let up = DVec3::from_array(source.policy.up).normalize();
    let planes = &source.candidate_planes;
    let (sin, cos) = source.policy.angle.sin_cos();
    let wall = |i: usize| normal(&planes[i]).dot(up).abs() <= sin;
    let point = |n: usize| DVec3::from_array(source.candidate_points[n]);
    let count = source.surfaces.len();
    let cell = tolerance.max(policy.precision) * 2.;
    let key = |p: DVec3| {
        (
            (p.x / cell).floor() as i64,
            (p.y / cell).floor() as i64,
            (p.z / cell).floor() as i64,
        )
    };
    let mut grid = BTreeMap::<(i64, i64, i64), Vec<(usize, usize)>>::new();
    for w in (0..count).filter(|&w| wall(w)) {
        for &n in &source.surfaces[w].nodes {
            grid.entry(key(point(n))).or_default().push((n, w));
        }
    }
    // Contact nodes of each parallel pair (a, b), a < b.
    let mut contacts = BTreeMap::<(usize, usize), [BTreeSet<usize>; 2]>::new();
    for (&(x, y, z), items) in &grid {
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let Some(others) = grid.get(&(x + dx, y + dy, z + dz)) else {
                        continue;
                    };
                    for &(n, a) in items {
                        for &(m, b) in others {
                            if a >= b
                                || representatives[a] == representatives[b]
                                || normal(&planes[a]).dot(normal(&planes[b])).abs() < cos
                                || point(n).distance(point(m)) > tolerance
                            {
                                continue;
                            }
                            let c = contacts.entry((a, b)).or_default();
                            c[0].insert(n);
                            c[1].insert(m);
                        }
                    }
                }
            }
        }
    }
    let elements: BTreeMap<u32, &crate::input::ElementData> =
        mesh.elements.iter().map(|e| (e.id, e)).collect();
    // Element polygons of a patch in a plane.
    let polygons = |patch: usize, plane: &PlaneFrame| -> Vec<Vec<[f64; 2]>> {
        source.surfaces[patch]
            .source_elements
            .iter()
            .filter_map(|id| {
                let e = elements.get(id)?;
                let ns = crate::reconstruction::planes::ordered_facet_nodes(mesh, e, plane, 1e-12)?;
                Some(
                    ns.iter()
                        .map(|n| plane.project(mesh.nodes[n].to_array()))
                        .collect(),
                )
            })
            .collect()
    };
    // Material overlap: an element centre of one patch strictly inside an
    // element of the other (both projected into the common plane).
    let overlap = |a: usize, b: usize, plane: &PlaneFrame| {
        let (pa, pb) = (polygons(a, plane), polygons(b, plane));
        let inside = |q: [f64; 2], poly: &[[f64; 2]]| {
            let mut c = false;
            for i in 0..poly.len() {
                let (u, v) = (poly[i], poly[(i + 1) % poly.len()]);
                if (u[1] > q[1]) != (v[1] > q[1])
                    && q[0] < u[0] + (q[1] - u[1]) * (v[0] - u[0]) / (v[1] - u[1])
                {
                    c = !c;
                }
            }
            c
        };
        let centre = |poly: &[[f64; 2]]| {
            let k = poly.len() as f64;
            let s = poly
                .iter()
                .fold([0., 0.], |s, p| [s[0] + p[0], s[1] + p[1]]);
            [s[0] / k, s[1] / k]
        };
        let hit = |from: &[Vec<[f64; 2]>], into: &[Vec<[f64; 2]>]| {
            from.iter().any(|p| {
                let q = centre(p);
                into.iter().any(|poly| inside(q, poly))
            })
        };
        hit(&pa, &pb) || hit(&pb, &pa)
    };
    let mut pairs: Vec<(f64, usize, usize, Vec<usize>)> = vec![];
    for (&(a, b), [ca, cb]) in &contacts {
        // The smaller wall moves onto the larger one.
        let size = |i: usize| source.surfaces[i].source_elements.len();
        let (moved, target, nodes) = if (size(a), b) < (size(b), a) {
            (a, b, ca)
        } else {
            (b, a, cb)
        };
        let nodes: Vec<usize> = nodes.iter().copied().collect();
        let d = offset(source, &nodes, &planes[target]);
        if d <= policy.closure_tolerance || d > tolerance + policy.closure_tolerance {
            continue;
        }
        pairs.push((d, moved, target, nodes));
    }
    pairs.sort_by(|x, y| x.0.total_cmp(&y.0).then((x.1, x.2).cmp(&(y.1, y.2))));
    for (_, u, l, nodes) in pairs {
        if representatives[u] == representatives[l] {
            continue;
        }
        let mut entry = StackedWall {
            upper: u,
            lower: l,
            root: representatives[l],
            slab: None,
            offset: offset(source, &nodes, &planes[l]),
            source_elements: source.surfaces[u].source_elements.clone(),
            reason: String::new(),
        };
        if overlap(u, l, &planes[l]) {
            entry.reason = "overlapping_parallel_walls".into();
            report.kept.push(entry);
            continue;
        }
        // The smaller group (a wall with the stack it carries) moves; if it
        // cannot, the other one moves instead.
        let group = |i: usize| -> usize {
            (0..count)
                .filter(|&j| representatives[j] == representatives[i])
                .map(|j| source.surfaces[j].source_elements.len())
                .sum()
        };
        let order = if (group(u), u) <= (group(l), l) {
            [(u, l), (l, u)]
        } else {
            [(l, u), (u, l)]
        };
        let mut reason = "";
        let chosen = order.into_iter().find(|&(a, b)| {
            let r = fits(
                mesh,
                source,
                policy,
                owners,
                rings,
                representatives,
                &report.shifted,
                representatives[a],
                representatives[b],
            );
            if let Err(e) = r {
                reason = e;
            }
            r.is_ok()
        });
        let Some((a, b)) = chosen else {
            entry.reason = reason.into();
            report.kept.push(entry);
            continue;
        };
        let (moved, target) = (representatives[a], representatives[b]);
        for (j, r) in representatives.iter_mut().enumerate() {
            if *r == moved {
                *r = target;
                report.shifted.insert(j);
            }
        }
        entry.upper = a;
        entry.lower = b;
        entry.root = target;
        entry.source_elements = source.surfaces[a].source_elements.clone();
        entry.reason = "aligned_wall_line".into();
        report.aligned.push(entry);
    }
}

#[cfg(test)]
mod tests {
    use super::movement_limit;

    #[test]
    fn limit_grows_with_distinct_aligned_structures_only() {
        let group = |s: usize| [0, 1, 1, 3][s];
        // No aligned support, or one: the plain limit.
        assert_eq!(movement_limit(0.05, &[0, 3], |_| false, group), (0.05, 0));
        assert_eq!(movement_limit(0.05, &[0, 3], |s| s == 0, group), (0.05, 1));
        // Two aligned supports of one group are one structure.
        assert_eq!(movement_limit(0.05, &[1, 2], |s| s > 0, group), (0.05, 1));
        // Two independent aligned structures: their vector sum.
        let (limit, k) = movement_limit(0.05, &[0, 1, 3], |s| s != 3, group);
        assert_eq!(k, 2);
        assert!((limit - 0.05 * 2f64.sqrt()).abs() < 1e-15);
    }
}
