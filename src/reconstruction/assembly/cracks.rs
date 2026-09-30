//! Region contours rebuilt across cracks of a converted mesh.
//!
//! A converted FE mesh may leave neighbouring elements of one wall or slab
//! with distinct nodes a few millimetres apart, or overlapping collinear
//! edges ending at distinct nodes, instead of shared nodes. The material is
//! one region (its elements are edge connected elsewhere), but its contour
//! runs into the crack and back, so it folds onto itself. The source mesh is
//! not welded: the contour of this region is rebuilt without the crack.
//!
//! A crack is an excursion of a contour between two unconnected nodes at
//! most the crack width apart: the contour path between them is more than
//! twice their distance, encloses void (never material: a thin fin of
//! elements keeps its contour) with a mean width within the crack width, and
//! the chord between them crosses no element. The excursion is cut off and
//! its two mouth nodes become one contour point. Before that, contour pieces
//! traversed from both sides (a zero-width crack along non-matching nodes of
//! two meshes) are removed. Enclosed narrow voids of positive width are not
//! cracks here: the hole policy (`features::simplify`) fills those no wider
//! than the crack width and keeps their nodes. Nodes left off
//! the contour keep no contour vertex; every change is reported with its
//! source nodes.
use super::{planes, MeshData, PlaneFrame};
use glam::DVec2;
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Serialize)]
pub struct Closure {
    pub patch: usize,
    pub source_elements: usize,
    /// Contour pieces traversed from both sides (a zero-width crack along
    /// non-matching nodes) removed before cutting excursions.
    pub overlapping_pieces_removed: usize,
    /// Number of crack excursions cut off a contour.
    pub cuts: usize,
    /// Mouth nodes of cut excursions within the crack width: [dropped, kept].
    pub identified: Vec<[u32; 2]>,
    /// Nodes of cut excursions and dropped crack contours.
    pub removed_nodes: Vec<u32>,
    /// Number of hole contours cut down to a line (cracks as a whole).
    pub dropped_contours: usize,
    /// Widest crack mouth identified into one contour point.
    pub maximum_width: f64,
    pub contour_nodes_before: usize,
    pub contour_nodes_after: usize,
}

fn signed_area(points: &[DVec2]) -> f64 {
    let o = points[0];
    (0..points.len())
        .map(|i| (points[i] - o).perp_dot(points[(i + 1) % points.len()] - o))
        .sum::<f64>()
        / 2.
}

/// Excursion `ring[i..=j]` (forward, wrapping): its length, the signed area
/// of the loop it closes with the chord j -> i, and its number of edges.
fn excursion(points: &[DVec2], i: usize, j: usize) -> (f64, f64, usize) {
    let n = points.len();
    let count = (j + n - i) % n;
    let path: Vec<DVec2> = (0..=count).map(|k| points[(i + k) % n]).collect();
    let length = path.windows(2).map(|w| w[0].distance(w[1])).sum();
    (length, signed_area(&path), count)
}

/// Element pieces (convex polygons) in a grid of cells four crack widths wide.
type MaterialIndex = (BTreeMap<(i64, i64), Vec<usize>>, Vec<Vec<DVec2>>);

/// Contour rings of one region rebuilt without cracks, or `None` when the
/// region has no crack. `rings` are the source boundary rings of the region.
pub fn close(
    mesh: &MeshData,
    elements: &BTreeMap<u32, &crate::input::ElementData>,
    ids: &[u32],
    rings: &[Vec<u32>],
    plane: &PlaneFrame,
    precision: f64,
    width: f64,
    patch: usize,
) -> Option<(Vec<Vec<u32>>, Closure)> {
    if width <= precision || rings.is_empty() {
        return None;
    }
    let uv: BTreeMap<u32, DVec2> = rings
        .iter()
        .flatten()
        .map(|&n| {
            let p = mesh.nodes.get(&n)?;
            Some((n, DVec2::from_array(plane.project(p.to_array()))))
        })
        .collect::<Option<_>>()?;
    let points = |ring: &[u32]| ring.iter().map(|n| uv[n]).collect::<Vec<_>>();
    let before: usize = rings.iter().map(Vec::len).sum();
    // Zero-width cracks first: contour edges traversed from both sides.
    let arranged = arrange(rings, &uv, 10. * precision, width);
    let overlaps = arranged.is_some();
    let overlaps_removed = arranged.as_ref().map_or(0, |a| a.1);
    let rings: &[Vec<u32>] = arranged.as_ref().map_or(rings, |a| &a.0);
    // The outer contour encloses material, the others enclose void.
    let outer = (0..rings.len())
        .max_by(|&a, &b| {
            signed_area(&points(&rings[a]))
                .abs()
                .total_cmp(&signed_area(&points(&rings[b])).abs())
        })
        .unwrap();
    let mut rings: Vec<(Vec<u32>, bool)> = rings
        .iter()
        .enumerate()
        .map(|(k, r)| (r.clone(), k == outer))
        .collect();
    let mut owners = BTreeMap::<u32, BTreeSet<u32>>::new();
    let mut facets = BTreeMap::new();
    for id in ids {
        let nodes = planes::ordered_facet_nodes(mesh, elements.get(id)?, plane, precision)?;
        for &n in &nodes {
            owners.entry(n).or_default().insert(*id);
        }
        facets.insert(*id, nodes);
    }
    let share = |a: u32, b: u32| {
        owners
            .get(&a)
            .zip(owners.get(&b))
            .is_some_and(|(x, y)| !x.is_disjoint(y))
    };
    let key = |p: DVec2| ((p.x / width).floor() as i64, (p.y / width).floor() as i64);
    let mut cuts = 0;
    let mut identified = vec![];
    let mut removed = BTreeSet::new();
    let mut dropped_contours = 0;
    let mut maximum_width = 0.0_f64;
    let mut material: Option<MaterialIndex> = None;
    loop {
        // The longest excursion first: a crack mouth before its inner pairs.
        let mut best: Option<(f64, f64, usize, usize, usize)> = None;
        for (r, (ring, is_outer)) in rings.iter().enumerate() {
            let n = ring.len();
            if n < 4 {
                continue;
            }
            let p = points(ring);
            let orientation = signed_area(&p).signum();
            let mut grid = BTreeMap::<(i64, i64), Vec<usize>>::new();
            for (k, &q) in p.iter().enumerate() {
                grid.entry(key(q)).or_default().push(k);
            }
            // Mouths: unconnected contour nodes within the crack width, and
            // the neighbours of a contour node where the contour turns back
            // (a crack of zero width along overlapping collinear edges).
            let mut mouths = BTreeSet::new();
            for i in 0..n {
                let c = key(p[i]);
                for dx in -1..=1 {
                    for dy in -1..=1 {
                        for &j in grid.get(&(c.0 + dx, c.1 + dy)).into_iter().flatten() {
                            if i != j && p[i].distance(p[j]) <= width && !share(ring[i], ring[j]) {
                                mouths.insert((i, j));
                            }
                        }
                    }
                }
                let (h, k) = ((i + n - 1) % n, (i + 1) % n);
                if (p[h] - p[i]).dot(p[k] - p[i]) > 0. {
                    mouths.insert((h, k));
                }
            }
            {
                for (i, j) in mouths {
                    let d = p[i].distance(p[j]);
                    let (length, area, count) = excursion(&p, i, j);
                    // A contour edge between unconnected nodes within the
                    // crack width is the mouth left by a zero-width crack.
                    let mouth = count == 1 && d <= width && !share(ring[i], ring[j]);
                    if !mouth && (count < 2 || count > n - 2 || length <= 2. * d + precision) {
                        continue;
                    }
                    // Void: opposite to a material contour, along a hole
                    // contour; a zero-area loop is void either way.
                    let zero = area.abs() <= precision * length.max(1.);
                    let void = zero
                        || if *is_outer {
                            area.signum() != orientation
                        } else {
                            area.signum() == orientation
                        };
                    // A mouth wider than the crack width only for a crack of
                    // zero width (overlapping collinear edges).
                    if !void || area.abs() > width * length / 2. || (d > width && !zero) {
                        continue;
                    }
                    if best.is_some_and(|b| {
                        length < b.0 || (length == b.0 && (d, r, i, j) >= (b.1, b.2, b.3, b.4))
                    }) {
                        continue;
                    }
                    // No part of the chord may run inside an element.
                    let index = material.get_or_insert_with(|| index(&facets, mesh, plane, width));
                    if crosses_material(p[i], p[j], index, width, precision) {
                        continue;
                    }
                    best = Some((length, d, r, i, j));
                }
            }
        }
        let Some((length, d, r, i, j)) = best else {
            break;
        };
        let area = excursion(&points(&rings[r].0), i, j).1;
        cuts += 1;
        let ring = &mut rings[r].0;
        let n = ring.len();
        let (a, b) = (ring[i], ring[j]);
        let count = (j + n - i) % n;
        let cut: BTreeSet<usize> = (1..=count).map(|k| (i + k) % n).collect();
        if d <= width {
            // The mouth keeps the node of more region elements.
            let degree = |x: u32| owners.get(&x).map_or(0, BTreeSet::len);
            let keep = if (degree(b), std::cmp::Reverse(b)) > (degree(a), std::cmp::Reverse(a)) {
                b
            } else {
                a
            };
            identified.push([if keep == a { b } else { a }, keep]);
            removed.extend(cut.iter().map(|&k| ring[k]));
            removed.insert(a);
            *ring = (0..n)
                .filter(|k| !cut.contains(k))
                .map(|k| if ring[k] == a { keep } else { ring[k] })
                .collect();
            maximum_width = maximum_width.max(d);
        } else {
            // A wide mouth of a thin crack: the contour follows the chord.
            removed.extend(cut.iter().filter(|&&k| k != j).map(|&k| ring[k]));
            *ring = (0..n)
                .filter(|&k| k == j || !cut.contains(&k))
                .map(|k| ring[k])
                .collect();
            maximum_width = maximum_width.max(2. * area.abs() / length);
        }
    }
    // A hole contour cut down to a line was a crack as a whole. Enclosed
    // narrow voids of positive width are left to the hole policy, which
    // fills those within the crack width and keeps their nodes.
    let mut kept = vec![];
    for (ring, is_outer) in rings {
        if ring.len() < 3 {
            if is_outer {
                return None;
            }
            dropped_contours += 1;
            removed.extend(ring.iter().copied());
            continue;
        }
        kept.push(ring);
    }
    if cuts == 0 && dropped_contours == 0 && !overlaps {
        return None;
    }
    let on_contour: BTreeSet<u32> = kept.iter().flatten().copied().collect();
    let removed: Vec<u32> = removed
        .into_iter()
        .filter(|n| !on_contour.contains(n))
        .collect();
    let after = kept.iter().map(Vec::len).sum();
    Some((
        kept,
        Closure {
            patch,
            source_elements: ids.len(),
            overlapping_pieces_removed: overlaps_removed,
            cuts,
            identified,
            removed_nodes: removed,
            dropped_contours,
            maximum_width,
            contour_nodes_before: before,
            contour_nodes_after: after,
        },
    ))
}

/// Element polygons in a grid (cell four crack widths) over their bounding
/// boxes, or along their edges for a large element.
fn index(
    facets: &BTreeMap<u32, Vec<u32>>,
    mesh: &MeshData,
    plane: &PlaneFrame,
    width: f64,
) -> MaterialIndex {
    let cell = 4. * width;
    let key = |p: DVec2| ((p.x / cell).floor() as i64, (p.y / cell).floor() as i64);
    let mut grid = BTreeMap::<(i64, i64), Vec<usize>>::new();
    let mut pieces = vec![];
    for ns in facets.values() {
        let p: Vec<DVec2> = ns
            .iter()
            .map(|n| DVec2::from_array(plane.project(mesh.nodes[n].to_array())))
            .collect();
        for piece in convex_pieces(p) {
            let (lo, hi) = piece.iter().fold(
                (DVec2::splat(f64::INFINITY), DVec2::splat(f64::NEG_INFINITY)),
                |(lo, hi), &q| (lo.min(q), hi.max(q)),
            );
            let (k0, k1) = (key(lo), key(hi));
            let mut cells = BTreeSet::new();
            if (k1.0 - k0.0 + 1) * (k1.1 - k0.1 + 1) <= 4096 {
                for x in k0.0..=k1.0 {
                    for y in k0.1..=k1.1 {
                        cells.insert((x, y));
                    }
                }
            } else {
                // A chord entering a large element crosses one of its edges:
                // the cells along the edges suffice.
                for i in 0..piece.len() {
                    cells.extend(cells_along(
                        piece[i],
                        piece[(i + 1) % piece.len()],
                        cell,
                        key,
                    ));
                }
            }
            for c in cells {
                grid.entry(c).or_default().push(pieces.len());
            }
            pieces.push(piece);
        }
    }
    (grid, pieces)
}

/// A shell facet as convex pieces: itself, or a nonconvex quadrilateral as
/// the two triangles of its inner diagonal.
fn convex_pieces(p: Vec<DVec2>) -> Vec<Vec<DVec2>> {
    let n = p.len();
    let turn = |i: usize| (p[i] - p[(i + n - 1) % n]).perp_dot(p[(i + 1) % n] - p[i]);
    let sign = signed_area(&p).signum();
    let reflex: Vec<usize> = (0..n).filter(|&i| turn(i) * sign < 0.).collect();
    match (n, reflex.as_slice()) {
        (_, []) => vec![p],
        (4, [r]) => {
            let r = *r;
            let (a, b, c) = (p[(r + 1) % 4], p[(r + 2) % 4], p[(r + 3) % 4]);
            vec![vec![p[r], a, b], vec![p[r], b, c]]
        }
        // Not a valid shell facet; its triangle fan is a conservative cover.
        _ => (1..n - 1).map(|i| vec![p[0], p[i], p[i + 1]]).collect(),
    }
}

/// Whether some part of segment a-b lies strictly inside an element, more
/// than `tolerance` from its boundary. Touching edges or vertices (a crack
/// of zero width runs along element edges) is not crossing.
fn crosses_material(
    a: DVec2,
    b: DVec2,
    (grid, pieces): &MaterialIndex,
    width: f64,
    tolerance: f64,
) -> bool {
    let cell = 4. * width;
    let key = |p: DVec2| ((p.x / cell).floor() as i64, (p.y / cell).floor() as i64);
    let mut candidates = BTreeSet::<usize>::new();
    for c in cells_along(a, b, cell, key) {
        for dx in -1..=1 {
            for dy in -1..=1 {
                candidates.extend(grid.get(&(c.0 + dx, c.1 + dy)).into_iter().flatten());
            }
        }
    }
    let d = b - a;
    candidates.into_iter().any(|k| {
        let piece = &pieces[k];
        let sign = signed_area(piece).signum();
        // Clip the segment by every edge's inner half-plane, inset by the
        // tolerance (Cyrus-Beck).
        let (mut t0, mut t1) = (0.0_f64, 1.0_f64);
        for i in 0..piece.len() {
            let (p, q) = (piece[i], piece[(i + 1) % piece.len()]);
            let e = q - p;
            let length = e.length();
            if length == 0. {
                continue;
            }
            // Signed distance inside the edge: positive towards the interior.
            let inside = |x: DVec2| sign * e.perp_dot(x - p) / length - tolerance;
            let (fa, fd) = (inside(a), sign * e.perp_dot(d) / length);
            if fd.abs() < f64::EPSILON * length {
                if fa <= 0. {
                    return false;
                }
                continue;
            }
            let t = -fa / fd;
            if fd > 0. {
                t0 = t0.max(t);
            } else {
                t1 = t1.min(t);
            }
            if t0 >= t1 {
                return false;
            }
        }
        (t1 - t0) * d.length() > tolerance
    })
}

/// Grid cells visited by a segment, sampled at most one cell apart.
fn cells_along(
    a: DVec2,
    b: DVec2,
    cell: f64,
    key: impl Fn(DVec2) -> (i64, i64),
) -> BTreeSet<(i64, i64)> {
    let steps = ((a.distance(b) / cell).ceil() as usize).max(1);
    (0..=steps)
        .map(|k| key(a.lerp(b, k as f64 / steps as f64)))
        .collect()
}

/// Contour rings with zero-width cracks removed, and the number of removed
/// pieces. Every contour edge is split at contour nodes lying on it (within
/// `tolerance`); a piece traversed twice is between material on both sides
/// and is dropped; the remaining pieces are linked into rings. `None` when
/// nothing overlaps or the pieces do not form simple rings.
fn arrange(
    rings: &[Vec<u32>],
    uv: &BTreeMap<u32, DVec2>,
    tolerance: f64,
    cell: f64,
) -> Option<(Vec<Vec<u32>>, usize)> {
    let key = |p: DVec2| ((p.x / cell).floor() as i64, (p.y / cell).floor() as i64);
    let mut grid = BTreeMap::<(i64, i64), Vec<u32>>::new();
    for &n in rings.iter().flatten() {
        grid.entry(key(uv[&n])).or_default().push(n);
    }
    let mut pieces = BTreeMap::<[u32; 2], usize>::new();
    for ring in rings {
        for i in 0..ring.len() {
            let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
            let (pa, pb) = (uv[&a], uv[&b]);
            let d = pb - pa;
            let mut on = BTreeSet::new();
            let steps = ((pa.distance(pb) / cell).ceil() as usize).max(1);
            for s in 0..=steps {
                let c = key(pa.lerp(pb, s as f64 / steps as f64));
                for dx in -1..=1 {
                    for dy in -1..=1 {
                        for &x in grid.get(&(c.0 + dx, c.1 + dy)).into_iter().flatten() {
                            if x == a || x == b {
                                continue;
                            }
                            let t = (uv[&x] - pa).dot(d) / d.length_squared();
                            let e = uv[&x].distance(pa + d * t);
                            if t > 0. && t < 1. && e <= tolerance {
                                on.insert((t.to_bits(), x));
                            }
                        }
                    }
                }
            }
            let mut chain = vec![a];
            // Nonnegative f64 bits order like the values.
            chain.extend(on.into_iter().map(|(_, x)| x));
            chain.push(b);
            for w in chain.windows(2) {
                if w[0] != w[1] {
                    *pieces.entry([w[0].min(w[1]), w[0].max(w[1])]).or_default() += 1;
                }
            }
        }
    }
    let removed = pieces.values().filter(|&&n| n == 2).count();
    if removed == 0 || pieces.values().any(|&n| n > 2) {
        return None;
    }
    let mut adjacency = BTreeMap::<u32, Vec<u32>>::new();
    for ([a, b], n) in pieces {
        if n == 1 {
            adjacency.entry(a).or_default().push(b);
            adjacency.entry(b).or_default().push(a);
        }
    }
    if adjacency.values().any(|v| v.len() != 2) {
        return None;
    }
    let mut remaining: BTreeSet<_> = adjacency.keys().copied().collect();
    let mut result = vec![];
    while let Some(&start) = remaining.first() {
        let mut ring = vec![];
        let (mut previous, mut current) = (start, start);
        loop {
            if !remaining.remove(&current) {
                return None;
            }
            ring.push(current);
            let next = adjacency[&current]
                .iter()
                .copied()
                .find(|n| *n != previous)?;
            previous = current;
            current = next;
            if current == start {
                break;
            }
        }
        if ring.len() < 3 {
            return None;
        }
        result.push(ring);
    }
    (!result.is_empty()).then_some((result, removed))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::input::ElementData;
    use glam::{DQuat, DVec3};

    struct Place {
        rotation: DQuat,
        shift: DVec3,
        scale: f64,
        id_stride: u32,
    }
    fn places() -> Vec<Place> {
        let mut out = vec![];
        for scale in [0.2, 1., 30.] {
            for (rotation, shift, id_stride) in [
                (DQuat::IDENTITY, DVec3::ZERO, 1),
                (
                    DQuat::from_axis_angle(DVec3::new(1., 2., 3.).normalize(), 0.61),
                    DVec3::new(812.5, -304.25, 71.),
                    7,
                ),
            ] {
                out.push(Place {
                    rotation,
                    shift,
                    scale,
                    id_stride,
                });
            }
        }
        out
    }

    /// Quads given by planar coordinates; `shared` coordinates are one
    /// node, a coordinate repeated with a `tag` is a distinct node.
    fn mesh(place: &Place, quads: &[[(f64, f64, u8); 4]]) -> (MeshData, PlaneFrame, Vec<u32>) {
        let mut mesh = MeshData::default();
        let mut ids = BTreeMap::new();
        let mut next = 1_000_u32;
        for (k, quad) in quads.iter().enumerate() {
            let mut nodes = vec![];
            for &(x, y, tag) in quad {
                let key = (x.to_bits(), y.to_bits(), tag);
                let id = *ids.entry(key).or_insert_with(|| {
                    // Scrambled numbering.
                    next = next.wrapping_mul(place.id_stride).wrapping_add(13) % 100_003;
                    next
                });
                let p = place.rotation * (DVec3::new(x, y, 0.) * place.scale) + place.shift;
                mesh.nodes.insert(id, p);
                nodes.push(id);
            }
            mesh.elements.push(ElementData {
                id: 1 + k as u32 * place.id_stride,
                elem_type: 44,
                stiff_id: 1,
                nodes,
            });
        }
        let plane = PlaneFrame::new(
            place.shift.to_array(),
            (place.rotation * DVec3::Z).to_array(),
        )
        .unwrap();
        let ids = mesh.elements.iter().map(|e| e.id).collect();
        (mesh, plane, ids)
    }

    fn rebuild(
        place: &Place,
        quads: &[[(f64, f64, u8); 4]],
    ) -> (Vec<Vec<u32>>, Option<Closure>, MeshData, PlaneFrame) {
        let (mesh, plane, ids) = mesh(place, quads);
        let precision = 1e-7 * place.scale;
        let elements: BTreeMap<_, _> = mesh.elements.iter().map(|e| (e.id, e)).collect();
        let rings = super::super::boundary(&mesh, &elements, &ids, &plane, precision).unwrap();
        match close(
            &mesh,
            &elements,
            &ids,
            &rings,
            &plane,
            precision,
            0.01 * place.scale,
            0,
        ) {
            Some((r, c)) => (r, Some(c), mesh, plane),
            None => (rings, None, mesh, plane),
        }
    }

    fn area(rings: &[Vec<u32>], mesh: &MeshData, plane: &PlaneFrame, scale: f64) -> Vec<f64> {
        let mut a: Vec<f64> = rings
            .iter()
            .map(|r| {
                let p: Vec<DVec2> = r
                    .iter()
                    .map(|n| DVec2::from_array(plane.project(mesh.nodes[n].to_array())))
                    .collect();
                signed_area(&p).abs() / (scale * scale)
            })
            .collect();
        a.sort_by(|x, y| y.total_cmp(x));
        a
    }

    fn quad(x0: f64, x1: f64, y0: f64, y1: f64) -> [(f64, f64, u8); 4] {
        [(x0, y0, 0), (x1, y0, 0), (x1, y1, 0), (x0, y1, 0)]
    }

    #[test]
    fn wedge_crack_in_a_wall_is_left_out_of_the_contour() {
        for place in places() {
            // A 4 x 2 wall; at x = 2 the two lower elements have distinct
            // bottom nodes 1.5 mm apart and share the node at y = 1.
            let mut quads = vec![];
            for x in 0..4 {
                for y in 0..2 {
                    quads.push(quad(x as f64, x as f64 + 1., y as f64, y as f64 + 1.));
                }
            }
            quads[2] = [(1., 0., 0), (2., 0., 0), (2., 1., 0), (1., 1., 0)];
            quads[4] = [(2.0015, 0., 1), (3., 0., 0), (3., 1., 0), (2., 1., 0)];
            let (rings, closure, mesh, plane) = rebuild(&place, &quads);
            let closure = closure.expect("crack closed");
            assert_eq!(rings.len(), 1);
            assert_eq!(closure.identified.len(), 1);
            assert!((closure.maximum_width - 0.0015 * place.scale).abs() < 1e-9 * place.scale);
            // The crack tip is not on the contour; the wall keeps its outline.
            let a = area(&rings, &mesh, &plane, place.scale);
            assert!((a[0] - 8.).abs() < 1e-3, "{a:?}");
            let unique: BTreeSet<_> = rings[0].iter().collect();
            assert_eq!(unique.len(), rings[0].len());
            // A 20 mm slit is wider than the crack width and stays.
            quads[4][0] = (2.02, 0., 1);
            let (_, closure, _, _) = rebuild(&place, &quads);
            assert!(closure.is_none(), "{closure:?}");
        }
    }

    #[test]
    fn non_matching_interface_is_interior_and_openings_stay() {
        for place in places() {
            // Lower and upper strips meet along y = 1 with non-matching
            // nodes between x = 1 and x = 3.
            let mut quads = vec![];
            for x in 0..4 {
                quads.push(quad(x as f64, x as f64 + 1., 0., 1.));
            }
            let top = [0., 1., 1.5, 2.5, 3., 4.];
            for w in top.windows(2) {
                quads.push(quad(w[0], w[1], 1., 2.));
            }
            let (rings, closure, mesh, plane) = rebuild(&place, &quads);
            let closure = closure.expect("interface removed");
            assert!(closure.overlapping_pieces_removed > 0);
            assert_eq!(rings.len(), 1);
            assert!((area(&rings, &mesh, &plane, place.scale)[0] - 8.).abs() < 1e-6);
            // A notch: the element over [1.5, 2.5] x [1, 2] is missing.
            let quads: Vec<_> = quads
                .into_iter()
                .filter(|q| !(q[0].0 == 1.5 && q[0].1 == 1.))
                .collect();
            let (rings, _, mesh, plane) = rebuild(&place, &quads);
            let a = area(&rings, &mesh, &plane, place.scale);
            // The notch left by the missing element stays in the outline.
            assert_eq!(rings.len(), 1);
            assert!((a[0] - 7.).abs() < 1e-6, "{a:?}");
        }
    }

    fn material(pieces: Vec<Vec<DVec2>>, width: f64) -> MaterialIndex {
        let cell = 4. * width;
        let key = |p: DVec2| ((p.x / cell).floor() as i64, (p.y / cell).floor() as i64);
        let mut grid = BTreeMap::<(i64, i64), Vec<usize>>::new();
        let pieces: Vec<Vec<DVec2>> = pieces.into_iter().flat_map(convex_pieces).collect();
        for (k, piece) in pieces.iter().enumerate() {
            let mut cells = BTreeSet::new();
            for i in 0..piece.len() {
                cells.extend(cells_along(
                    piece[i],
                    piece[(i + 1) % piece.len()],
                    cell,
                    key,
                ));
            }
            for c in cells {
                grid.entry(c).or_default().push(k);
            }
        }
        (grid, pieces)
    }

    #[test]
    fn a_chord_crossing_an_element_away_from_its_midpoint_is_rejected() {
        let v = |x: f64, y: f64| DVec2::new(x, y);
        for scale in [0.01, 1., 40.] {
            let width = 0.01 * scale;
            let tol = 1e-7 * scale;
            // A 2 mm element corner pokes into the gap near one chord end;
            // the chord midpoint is in the void.
            let index = material(
                vec![vec![
                    v(0.007, -0.001),
                    v(0.009, -0.001),
                    v(0.009, 0.001),
                    v(0.007, 0.001),
                ]
                .into_iter()
                .map(|p| p * scale)
                .collect()],
                width,
            );
            let (a, b) = (v(0., 0.) * scale, v(0.01, 0.) * scale);
            assert!(crosses_material(a, b, &index, width, tol));
            assert!(crosses_material(b, a, &index, width, tol));
            // Along an element edge or through a corner only: not crossing.
            let (c, d) = (v(0.005, 0.001) * scale, v(0.012, 0.001) * scale);
            assert!(!crosses_material(c, d, &index, width, tol));
            let (e, f) = (v(0.007, 0.003) * scale, v(0.011, -0.001) * scale);
            assert!(!crosses_material(e, f, &index, width, tol));
        }
        // A nonconvex quadrilateral is covered exactly by its two triangles:
        // its notch is void.
        let dart = vec![v(0., 0.), v(4., 2.), v(0., 4.), v(1., 2.)];
        let index = material(vec![dart], 0.5);
        assert!(!crosses_material(v(0.1, 2.), v(0.9, 2.), &index, 0.5, 1e-7));
        assert!(crosses_material(v(1.5, 1.), v(1.5, 3.), &index, 0.5, 1e-7));
    }

    #[test]
    fn thin_fin_of_material_is_never_cut() {
        for place in places() {
            // A block with a 5 mm fin of two elements on its top edge.
            let quads = vec![
                quad(0., 1., 0., 1.),
                quad(1., 1.0025, 0., 1.),
                quad(1.0025, 1.005, 0., 1.),
                quad(1.005, 2., 0., 1.),
                quad(1., 1.0025, 1., 2.),
                quad(1.0025, 1.005, 1., 2.),
            ];
            let (_, closure, _, _) = rebuild(&place, &quads);
            assert!(closure.is_none(), "{closure:?}");
        }
    }
}
