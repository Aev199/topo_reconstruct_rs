//! Virtual incidences: source nodes that miss a support plane by a gap.
//!
//! A source model may leave a structure a few centimetres short of another
//! one it bears on or abuts (a wall end 26 mm from the corner of two other
//! walls, a wall top 25 mm below a slab). Such a node near the material of a
//! non-parallel support plane, but not on it, becomes a virtual incidence:
//! the frame solve then shifts the planes and nodes with the least total
//! movement so that it lies on that plane, as if the node had been shared.
//! Parallel structures (slabs at different levels) are never joined.
use super::super::planes;
use crate::input::MeshData;
use glam::{DVec2, DVec3};
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Serialize)]
pub struct Incidence {
    pub node_id: u32,
    /// Plane patch the node is made incident to.
    pub patch: usize,
    /// Source distance to the patch plane and to its material.
    pub height: f64,
    pub distance: f64,
}

#[derive(Debug, Clone, Default, Serialize)]
pub struct Report {
    pub tolerance: f64,
    pub candidates: usize,
    /// Incidences in the accepted solve.
    pub applied: Vec<Incidence>,
    /// Incidences dropped because the solve could not satisfy them within
    /// tolerances and movement budgets.
    pub dropped: Vec<Incidence>,
    pub rounds: usize,
}

type Key = (i64, i64, i64);

/// Nodes within `tolerance` of the material of a non-parallel patch they do
/// not belong to, off its plane by at least `minimum`.
pub fn candidates(
    mesh: &MeshData,
    planes: &planes::Report,
    node_ids: &[u32],
    tolerance: f64,
    minimum: f64,
    angle: f64,
) -> Vec<Incidence> {
    let elements: BTreeMap<_, _> = mesh.elements.iter().map(|e| (e.id, e)).collect();
    let mut owners = BTreeMap::<u32, BTreeSet<usize>>::new();
    for (i, p) in planes.patches.iter().enumerate() {
        for &n in &p.source_nodes {
            owners.entry(n).or_default().insert(i);
        }
    }
    let normal = |i: usize| DVec3::from_array(planes.patches[i].plane.normal);
    // A gap is at the edge of a structure: only contour nodes of their own
    // patches (free element edges) are candidates, never interior nodes.
    let mut boundary = BTreeSet::<u32>::new();
    for patch in &planes.patches {
        let mut counts = BTreeMap::<[u32; 2], usize>::new();
        for id in &patch.source_elements {
            let Some(e) = elements.get(id) else { continue };
            let Some(ns) = planes::ordered_facet_nodes(mesh, e, &patch.plane, 1e-9) else {
                continue;
            };
            for i in 0..ns.len() {
                let (a, b) = (ns[i], ns[(i + 1) % ns.len()]);
                *counts.entry([a.min(b), a.max(b)]).or_default() += 1;
            }
        }
        boundary.extend(
            counts
                .into_iter()
                .filter(|&(_, c)| c == 1)
                .flat_map(|(e, _)| e),
        );
    }
    // Coarse grid of candidate nodes.
    let cell = tolerance.max(1e-9) * 20.;
    let key = |p: DVec3| -> Key {
        (
            (p.x / cell).floor() as i64,
            (p.y / cell).floor() as i64,
            (p.z / cell).floor() as i64,
        )
    };
    let mut grid = BTreeMap::<Key, Vec<u32>>::new();
    for &n in node_ids {
        if let Some(p) = mesh.nodes.get(&n) {
            grid.entry(key(*p)).or_default().push(n);
        }
    }
    let mut out = vec![];
    for (b, patch) in planes.patches.iter().enumerate() {
        let plane = &patch.plane;
        let nb = normal(b);
        // Element polygons of the patch in its plane.
        let polygons: Vec<Vec<DVec2>> = patch
            .source_elements
            .iter()
            .filter_map(|id| {
                let e = elements.get(id)?;
                let ns = planes::ordered_facet_nodes(mesh, e, plane, 1e-9)?;
                Some(
                    ns.iter()
                        .map(|n| DVec2::from_array(plane.project(mesh.nodes[n].to_array())))
                        .collect(),
                )
            })
            .collect();
        if polygons.is_empty() {
            continue;
        }
        let (lo, hi) = patch.source_nodes.iter().map(|n| mesh.nodes[n]).fold(
            (DVec3::splat(f64::INFINITY), DVec3::splat(f64::NEG_INFINITY)),
            |(lo, hi), p| (lo.min(p), hi.max(p)),
        );
        let (k0, k1) = (key(lo - tolerance), key(hi + tolerance));
        // Polygons indexed by their boxes in the plane.
        let pcell = cell;
        let pkey = |p: DVec2| ((p.x / pcell).floor() as i64, (p.y / pcell).floor() as i64);
        let mut pgrid = BTreeMap::<(i64, i64), Vec<usize>>::new();
        for (k, poly) in polygons.iter().enumerate() {
            let (a, z) = poly.iter().fold(
                (DVec2::splat(f64::INFINITY), DVec2::splat(f64::NEG_INFINITY)),
                |(a, z), q| (a.min(*q), z.max(*q)),
            );
            let (c0, c1) = (pkey(a - tolerance), pkey(z + tolerance));
            for x in c0.0..=c1.0 {
                for y in c0.1..=c1.1 {
                    pgrid.entry((x, y)).or_default().push(k);
                }
            }
        }
        for (_, nodes) in grid.range(k0..=k1) {
            for &n in nodes {
                let p = mesh.nodes[&n];
                if (p + tolerance).cmplt(lo).any() || (p - tolerance).cmpgt(hi).any() {
                    continue;
                }
                let own = owners.get(&n);
                if own.is_some_and(|o| o.contains(&b)) {
                    continue;
                }
                if own.is_some() && !boundary.contains(&n) {
                    continue;
                }
                // Parallel structures are never joined.
                if own.is_some_and(|o| o.iter().any(|&a| normal(a).cross(nb).length() < angle)) {
                    continue;
                }
                let height = plane.distance(p.to_array()).abs();
                if height < minimum || height >= tolerance {
                    continue;
                }
                // A node of the same structure nearer to the plane (on it or
                // between the node and it) next to the projection: the
                // offset is a step of that structure (a contour edge or
                // crack mouth to that node) or its own extent (a narrow
                // panel), not a gap; pulling the node onto the plane would
                // collapse the structure there.
                let q3 = p - nb * plane.distance(p.to_array());
                let step = (-1..=1).any(|dx| {
                    (-1..=1).any(|dy| {
                        (-1..=1).any(|dz| {
                            let c = key(p);
                            grid.get(&(c.0 + dx, c.1 + dy, c.2 + dz))
                                .into_iter()
                                .flatten()
                                .any(|&m| {
                                    let pm = mesh.nodes[&m];
                                    let hm = plane.distance(pm.to_array());
                                    m != n
                                        && pm.distance(q3) < tolerance
                                        && (hm.abs() < minimum
                                            || (hm * plane.distance(p.to_array()) > 0.
                                                && hm.abs() < height))
                                        && owners.get(&m).is_some_and(|om| {
                                            own.is_some_and(|o| o.intersection(om).next().is_some())
                                        })
                                })
                        })
                    })
                });
                if step {
                    continue;
                }
                let q = DVec2::from_array(plane.project(p.to_array()));
                let mut outside = f64::INFINITY;
                for &k in pgrid.get(&pkey(q)).into_iter().flatten() {
                    let poly = &polygons[k];
                    if contains(poly, q) {
                        outside = 0.;
                        break;
                    }
                    for i in 0..poly.len() {
                        let (a, z) = (poly[i], poly[(i + 1) % poly.len()]);
                        let d = z - a;
                        let t = ((q - a).dot(d) / d.length_squared()).clamp(0., 1.);
                        outside = outside.min(q.distance(a + d * t));
                    }
                }
                let distance = height.hypot(outside);
                if distance < tolerance {
                    out.push(Incidence {
                        node_id: n,
                        patch: b,
                        height,
                        distance,
                    });
                }
            }
        }
    }
    // A node near two parallel planes (between the offset walls of a
    // stacked pair) is ambiguous: that case belongs to wall alignment.
    let mut by_node = BTreeMap::<u32, Vec<usize>>::new();
    for c in &out {
        by_node.entry(c.node_id).or_default().push(c.patch);
    }
    let ambiguous: BTreeSet<u32> = by_node
        .iter()
        .filter(|(_, ps)| {
            ps.iter().enumerate().any(|(i, &a)| {
                ps[i + 1..]
                    .iter()
                    .any(|&b| normal(a).cross(normal(b)).length() < angle)
            })
        })
        .map(|(&n, _)| n)
        .collect();
    out.retain(|c| !ambiguous.contains(&c.node_id));
    out.sort_by(|x, y| {
        x.distance
            .total_cmp(&y.distance)
            .then(x.node_id.cmp(&y.node_id))
            .then(x.patch.cmp(&y.patch))
    });
    out
}

/// Incidence nodes at which a source element edge became shorter than
/// `minimum` in the solved `points` (it was not in the source).
pub fn collapsed_edges(
    mesh: &MeshData,
    incidences: &[Incidence],
    node_ids: &[u32],
    points: &[[f64; 3]],
    minimum: f64,
) -> BTreeSet<u32> {
    let nodes: BTreeSet<u32> = incidences.iter().map(|i| i.node_id).collect();
    let index: BTreeMap<u32, usize> = node_ids.iter().enumerate().map(|(i, &n)| (n, i)).collect();
    let solved = |n: u32| index.get(&n).map(|&i| DVec3::from_array(points[i]));
    let mut out = BTreeSet::new();
    for e in &mesh.elements {
        if !e.is_shell() || !e.nodes.iter().any(|n| nodes.contains(n)) {
            continue;
        }
        for &a in &e.nodes {
            for &b in &e.nodes {
                if a >= b || !(nodes.contains(&a) || nodes.contains(&b)) {
                    continue;
                }
                let (Some(pa), Some(pb)) = (solved(a), solved(b)) else {
                    continue;
                };
                if pa.distance(pb) < minimum && mesh.nodes[&a].distance(mesh.nodes[&b]) >= minimum {
                    out.extend([a, b].into_iter().filter(|n| nodes.contains(n)));
                }
            }
        }
    }
    out
}

fn contains(poly: &[DVec2], q: DVec2) -> bool {
    let mut inside = false;
    for i in 0..poly.len() {
        let (a, b) = (poly[i], poly[(i + 1) % poly.len()]);
        if (a.y > q.y) != (b.y > q.y) && q.x < a.x + (q.y - a.y) * (b.x - a.x) / (b.y - a.y) {
            inside = !inside;
        }
    }
    inside
}
