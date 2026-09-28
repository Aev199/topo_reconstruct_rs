//! Trimming of thin consoles beyond structural junction lines.
//!
//! A mid-surface FE model often lets a slab reach the outer face of a wall
//! whose mid-plane lies inside the slab: a console of half the wall
//! thickness. It carries no structure, yet forces elements much smaller than
//! the target mesh size. A surface is divided by its junction lines into
//! faces; a group of thin faces between a junction line and a free contour is
//! removed, so the contour follows the junction line. Nothing attached is
//! removed: a face is kept if its free contour is shared with another
//! surface, or if it holds bars, openings, retained nodes or other lines.
//! Source provenance is unchanged; every trim is reported.
use super::bars::{Axis, Contact};
use crate::reconstruction::{Model, PlaneFrame};
use glam::DVec3;
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};

pub struct Context<'a> {
    /// Maximum console width (model units).
    pub maximum_width: f64,
    /// Mandatory interior nodes per surface (bar contacts, retained nodes).
    pub interior: &'a [BTreeSet<usize>],
    /// Vertices that must stay in their surfaces (bar anchors).
    pub locked: &'a BTreeSet<usize>,
    pub axes: &'a [Axis],
    pub contacts: &'a [Contact],
}

#[derive(Debug, Serialize)]
pub struct TrimmedConsole {
    pub surface: usize,
    pub source_elements: Vec<u32>,
    /// Largest distance of the removed material from its junction lines.
    pub width: f64,
    pub area: f64,
    /// Removed free contour length.
    pub free_length: f64,
    /// Junction edges that became the contour of the surface.
    pub junction_edges: Vec<usize>,
    /// Edges the surface no longer uses.
    pub released_edges: Vec<usize>,
}

#[derive(Debug, Default, Serialize)]
pub struct Report {
    pub maximum_width: f64,
    pub trimmed: Vec<TrimmedConsole>,
    /// Thin faces on a free contour that were kept, by reason.
    pub kept: BTreeMap<String, usize>,
    pub passes: usize,
}

type Uv = [f64; 2];

fn segment_distance(p: Uv, a: Uv, b: Uv) -> f64 {
    let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
    let l2 = dx * dx + dy * dy;
    let t = if l2 > 0. {
        (((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / l2).clamp(0., 1.)
    } else {
        0.
    };
    (p[0] - a[0] - t * dx).hypot(p[1] - a[1] - t * dy)
}

/// Strict interior of a simple polygon.
fn strictly_inside(polygon: &[Uv], p: Uv, eps: f64) -> bool {
    let n = polygon.len();
    if (0..n).any(|i| segment_distance(p, polygon[i], polygon[(i + 1) % n]) <= eps) {
        return false;
    }
    let mut inside = false;
    for i in 0..n {
        let (a, b) = (polygon[i], polygon[(i + 1) % n]);
        if (a[1] > p[1]) != (b[1] > p[1])
            && p[0] < a[0] + (p[1] - a[1]) * (b[0] - a[0]) / (b[1] - a[1])
        {
            inside = !inside;
        }
    }
    inside
}

/// Planar subdivision of one surface by its contour and embedded edges.
struct Faces {
    /// Half-edge `2k` runs along `edges[k]` from its first to its second
    /// vertex; `2k + 1` is its twin.
    edges: Vec<usize>,
    ends: Vec<[usize; 2]>,
    boundary: Vec<bool>,
    next: Vec<usize>,
    face: Vec<usize>,
    cycles: Vec<Vec<usize>>,
    area: Vec<f64>,
    material: Vec<bool>,
    /// Cycles of separate components lying inside each face.
    inner: Vec<Vec<usize>>,
    uv: BTreeMap<usize, Uv>,
}

impl Faces {
    fn tail(&self, h: usize) -> usize {
        let [a, b] = self.ends[h / 2];
        if h.is_multiple_of(2) {
            a
        } else {
            b
        }
    }
    fn head(&self, h: usize) -> usize {
        self.tail(h ^ 1)
    }

    fn new(model: &Model, s: usize) -> Self {
        let surface = &model.surfaces[s];
        let plane: &PlaneFrame = &model.planes[surface.plane];
        let boundary_set: BTreeSet<usize> = surface
            .boundaries
            .iter()
            .flatten()
            .map(|e| e.edge)
            .collect();
        let edges: Vec<usize> = model.surface_edges(s).collect();
        let ends: Vec<[usize; 2]> = edges.iter().map(|&e| model.edges[e]).collect();
        let uv: BTreeMap<usize, Uv> = ends
            .iter()
            .flatten()
            .map(|&v| (v, plane.project(model.vertices[v])))
            .collect();
        let mut faces = Faces {
            boundary: edges.iter().map(|e| boundary_set.contains(e)).collect(),
            edges,
            ends,
            next: vec![],
            face: vec![],
            cycles: vec![],
            area: vec![],
            material: vec![],
            inner: vec![],
            uv,
        };
        let count = faces.edges.len() * 2;
        let mut outgoing = BTreeMap::<usize, Vec<(f64, usize)>>::new();
        for h in 0..count {
            let (a, b) = (faces.uv[&faces.tail(h)], faces.uv[&faces.head(h)]);
            outgoing
                .entry(faces.tail(h))
                .or_default()
                .push(((b[1] - a[1]).atan2(b[0] - a[0]), h));
        }
        for list in outgoing.values_mut() {
            list.sort_by(|x, y| x.0.total_cmp(&y.0));
        }
        faces.next = (0..count)
            .map(|h| {
                let list = &outgoing[&faces.head(h)];
                let i = list.iter().position(|&(_, g)| g == h ^ 1).unwrap();
                list[(i + list.len() - 1) % list.len()].1
            })
            .collect();
        faces.face = vec![usize::MAX; count];
        for start in 0..count {
            if faces.face[start] != usize::MAX {
                continue;
            }
            let id = faces.cycles.len();
            let mut cycle = vec![];
            let mut h = start;
            while faces.face[h] == usize::MAX {
                faces.face[h] = id;
                cycle.push(h);
                h = faces.next[h];
            }
            let area = cycle
                .iter()
                .map(|&h| {
                    let (a, b) = (faces.uv[&faces.tail(h)], faces.uv[&faces.head(h)]);
                    a[0] * b[1] - b[0] * a[1]
                })
                .sum::<f64>()
                / 2.;
            faces.cycles.push(cycle);
            faces.area.push(area);
        }
        let eps = model.precision;
        faces.material = (0..faces.cycles.len())
            .map(|f| {
                if faces.area[f] <= 0. {
                    return false;
                }
                // A point just left of the longest half-edge of the cycle.
                let h = *faces.cycles[f]
                    .iter()
                    .max_by(|&&x, &&y| faces.length(x).total_cmp(&faces.length(y)))
                    .unwrap();
                let (a, b) = (faces.uv[&faces.tail(h)], faces.uv[&faces.head(h)]);
                let l = faces.length(h);
                let d = (l * 1e-3).min(model.minimum_edge * 0.1);
                let p = [
                    (a[0] + b[0]) / 2. - (b[1] - a[1]) / l * d,
                    (a[1] + b[1]) / 2. + (b[0] - a[0]) / l * d,
                ];
                crate::reconstruction::closed_contains(&surface.contours, p, eps)
            })
            .collect();
        // A cycle without positive area is a separate component inside some
        // face: an opening contour or a floating junction line. It belongs
        // to the smallest material face strictly containing it.
        let polygons: Vec<Vec<Uv>> = (0..faces.cycles.len()).map(|f| faces.polygon(f)).collect();
        faces.inner = vec![vec![]; faces.cycles.len()];
        for c in 0..faces.cycles.len() {
            if faces.area[c] > 0. {
                continue;
            }
            let p = faces.uv[&faces.tail(faces.cycles[c][0])];
            let owner = (0..faces.cycles.len())
                .filter(|&f| faces.material[f] && strictly_inside(&polygons[f], p, eps))
                .min_by(|&x, &y| faces.area[x].total_cmp(&faces.area[y]));
            if let Some(owner) = owner {
                for &h in &faces.cycles[c] {
                    faces.face[h] = owner;
                }
                faces.inner[owner].push(c);
            }
        }
        faces
    }

    fn length(&self, h: usize) -> f64 {
        let (a, b) = (self.uv[&self.tail(h)], self.uv[&self.head(h)]);
        (b[0] - a[0]).hypot(b[1] - a[1])
    }

    fn polygon(&self, f: usize) -> Vec<Uv> {
        self.cycles[f]
            .iter()
            .map(|&h| self.uv[&self.tail(h)])
            .collect()
    }

    /// Half-edges of a face: its outer cycle and the cycles inside it.
    fn members(&self, f: usize) -> Vec<usize> {
        self.cycles[f]
            .iter()
            .chain(self.inner[f].iter().flat_map(|&c| self.cycles[c].iter()))
            .copied()
            .collect()
    }

    /// Area of a face net of the components inside it.
    fn net_area(&self, f: usize) -> f64 {
        self.area[f] + self.inner[f].iter().map(|&c| self.area[c]).sum::<f64>()
    }

    /// Strict interior of a face: inside its outer cycle, outside (and not
    /// on) every component inside it.
    fn contains(&self, f: usize, p: Uv, eps: f64) -> bool {
        strictly_inside(&self.polygon(f), p, eps)
            && self.inner[f].iter().all(|&c| {
                let polygon = self.polygon(c);
                let n = polygon.len();
                !strictly_inside(&polygon, p, eps)
                    && (0..n).all(|i| segment_distance(p, polygon[i], polygon[(i + 1) % n]) > eps)
            })
    }
}

/// Why a thin face on a free contour must be kept, if it must.
fn blocked(
    model: &Model,
    s: usize,
    faces: &Faces,
    f: usize,
    context: &Context<'_>,
) -> Option<&'static str> {
    let eps = model.precision;
    let cycle = faces.members(f);
    if cycle.iter().any(|&h| faces.face[h ^ 1] == f) {
        return Some("dangling_line_inside");
    }
    let polygon = faces.polygon(f);
    let plane = &model.planes[model.surfaces[s].plane];
    if faces.inner[f]
        .iter()
        .any(|&c| faces.cycles[c].iter().any(|&h| faces.boundary[h / 2]))
    {
        return Some("opening_inside");
    }
    for &h in &cycle {
        if !faces.boundary[h / 2] {
            continue;
        }
        let e = faces.edges[h / 2];
        if (0..model.surfaces.len()).any(|o| o != s && model.surface_edges(o).any(|x| x == e)) {
            return Some("contour_shared_with_other_surface");
        }
        let v = faces.tail(h);
        if context.locked.contains(&v) {
            return Some("bar_anchor_on_contour");
        }
    }
    let contour_vertices: BTreeSet<usize> = cycle
        .iter()
        .filter(|&&h| faces.boundary[h / 2])
        .flat_map(|&h| [faces.tail(h), faces.head(h)])
        .collect();
    for &v in &context.interior[s] {
        let p = plane.project(model.vertices[v]);
        if faces.contains(f, p, eps) || contour_vertices.contains(&v) {
            return Some("retained_node_inside");
        }
    }
    for contact in context.contacts {
        let Contact::Interval {
            axis,
            surface,
            start_t,
            end_t,
            ..
        } = *contact
        else {
            continue;
        };
        if surface != s {
            continue;
        }
        let [a, b] = context.axes[axis]
            .endpoints
            .map(|v| DVec3::from_array(model.vertices[v]));
        let (p, q) = (a.lerp(b, start_t), a.lerp(b, end_t));
        let samples = [p, q, p.lerp(q, 0.5)].map(|x| plane.project(x.to_array()));
        let (pu, qu) = (samples[0], samples[1]);
        if samples.iter().any(|&x| faces.contains(f, x, eps))
            || crate::reconstruction::properly_crosses(std::slice::from_ref(&polygon), pu, qu, eps)
        {
            return Some("bar_inside");
        }
    }
    None
}

/// Trim console groups of one surface. Returns whether the surface changed.
fn trim_surface(model: &mut Model, s: usize, context: &Context<'_>, report: &mut Report) -> bool {
    let faces = Faces::new(model, s);
    let eps = model.precision;
    let material: Vec<usize> = (0..faces.cycles.len())
        .filter(|&f| faces.material[f])
        .collect();
    if material.len() < 2 {
        return false;
    }
    let mut candidate = BTreeSet::new();
    let mut widths = BTreeMap::new();
    for &f in &material {
        let cycle = faces.members(f);
        let junctions: Vec<usize> = cycle
            .iter()
            .copied()
            .filter(|&h| !faces.boundary[h / 2])
            .collect();
        if junctions.is_empty() || junctions.len() == cycle.len() {
            continue;
        }
        let mut probes: Vec<Uv> = cycle.iter().map(|&h| faces.uv[&faces.tail(h)]).collect();
        for &h in cycle.iter().filter(|&&h| faces.boundary[h / 2]) {
            let (a, b) = (faces.uv[&faces.tail(h)], faces.uv[&faces.head(h)]);
            probes.push([(a[0] + b[0]) / 2., (a[1] + b[1]) / 2.]);
        }
        let width = probes
            .iter()
            .map(|&p| {
                junctions
                    .iter()
                    .map(|&h| {
                        segment_distance(p, faces.uv[&faces.tail(h)], faces.uv[&faces.head(h)])
                    })
                    .fold(f64::INFINITY, f64::min)
            })
            .fold(0.0_f64, f64::max);
        if width > context.maximum_width {
            continue;
        }
        if let Some(reason) = blocked(model, s, &faces, f, context) {
            *report.kept.entry(reason.into()).or_default() += 1;
            continue;
        }
        candidate.insert(f);
        widths.insert(f, width);
    }
    // Connected groups across shared junction edges. A group is a console
    // only if it borders material that stays, along a junction line.
    let mut removed = BTreeSet::new();
    let mut seen = BTreeSet::new();
    for &seed in &candidate {
        if !seen.insert(seed) {
            continue;
        }
        let mut group = vec![seed];
        let mut i = 0;
        while i < group.len() {
            for h in faces.members(group[i]) {
                let other = faces.face[h ^ 1];
                if !faces.boundary[h / 2] && candidate.contains(&other) && seen.insert(other) {
                    group.push(other);
                }
            }
            i += 1;
        }
        let anchored = group.iter().any(|&f| {
            faces.members(f).into_iter().any(|h| {
                let other = faces.face[h ^ 1];
                !faces.boundary[h / 2] && faces.material[other] && !candidate.contains(&other)
            })
        });
        if anchored {
            removed.extend(group);
        } else {
            *report
                .kept
                .entry("not_beyond_a_junction".into())
                .or_default() += group.len();
        }
    }
    if removed.is_empty() {
        return false;
    }
    let kept = |f: usize| faces.material[f] && !removed.contains(&f);
    // Walk the contour of the remaining material.
    let is_contour = |h: usize| kept(faces.face[h]) && !kept(faces.face[h ^ 1]);
    let mut used = BTreeSet::new();
    let mut rings = vec![];
    for start in 0..faces.face.len() {
        if !is_contour(start) || used.contains(&start) {
            continue;
        }
        let mut ring = vec![];
        let mut h = start;
        loop {
            used.insert(h);
            ring.push(faces.tail(h));
            let mut g = faces.next[h];
            let mut guard = 0;
            while !is_contour(g) {
                g = faces.next[g ^ 1];
                guard += 1;
                if guard > faces.face.len() {
                    return false;
                }
            }
            h = g;
            if h == start {
                break;
            }
            if used.contains(&h) {
                return false;
            }
        }
        rings.push(ring);
    }
    let signed = |ring: &Vec<usize>| {
        (0..ring.len())
            .map(|i| {
                let (a, b) = (faces.uv[&ring[i]], faces.uv[&ring[(i + 1) % ring.len()]]);
                a[0] * b[1] - b[0] * a[1]
            })
            .sum::<f64>()
            / 2.
    };
    if rings.iter().filter(|r| signed(r) > 0.).count() != 1 {
        *report
            .kept
            .entry("would_disconnect_surface".into())
            .or_default() += removed.len();
        return false;
    }
    rings.sort_by(|a, b| signed(b).total_cmp(&signed(a)));
    let embedded: Vec<usize> = (0..faces.edges.len())
        .filter(|&k| !faces.boundary[k] && kept(faces.face[2 * k]) && kept(faces.face[2 * k + 1]))
        .map(|k| faces.edges[k])
        .collect();
    let before: BTreeSet<usize> = faces.edges.iter().copied().collect();
    let area: f64 = removed.iter().map(|&f| faces.net_area(f)).sum();
    let free_length: f64 = removed
        .iter()
        .flat_map(|&f| faces.members(f))
        .filter(|&h| faces.boundary[h / 2])
        .map(|h| faces.length(h))
        .sum();
    let width = removed.iter().map(|f| widths[f]).fold(0.0_f64, f64::max);
    let junction_edges: Vec<usize> = (0..faces.face.len())
        .filter(|&h| !faces.boundary[h / 2] && is_contour(h))
        .map(|h| faces.edges[h / 2])
        .collect();
    if let Err(e) = model.rebuild_surface(s, rings, embedded) {
        *report
            .kept
            .entry(format!("invalid_after_trim_{e:?}"))
            .or_default() += removed.len();
        return false;
    }
    let after: BTreeSet<usize> = model.surface_edges(s).collect();
    report.trimmed.push(TrimmedConsole {
        surface: s,
        source_elements: model.surfaces[s].source_elements.clone(),
        width,
        area: area.max(eps * eps),
        free_length,
        junction_edges,
        released_edges: before.difference(&after).copied().collect(),
    });
    true
}

/// Trim thin consoles of all surfaces until no further change. Trimming a
/// slab console can free the contour of a wall stub beneath it, which is
/// then trimmed in a later pass; a stub standing on a wider plate stays.
pub fn trim(model: &mut Model, context: &Context<'_>) -> Report {
    let mut report = Report {
        maximum_width: context.maximum_width,
        ..Default::default()
    };
    for _ in 0..16 {
        report.passes += 1;
        report.kept.clear();
        let mut changed = false;
        for s in 0..model.surfaces.len() {
            changed |= trim_surface(model, s, context, &mut report);
        }
        if !changed {
            break;
        }
    }
    report
}

#[cfg(test)]
mod tests {
    use super::super::junctions::{
        self,
        tests::{build, no_interior, slab, Panel, Placement},
    };
    use super::*;

    fn wall_y(y: f64, x0: f64, x1: f64, z0: f64, z1: f64) -> Panel {
        (
            vec![vec![[x0, y, z0], [x1, y, z0], [x1, y, z1], [x0, y, z1]]],
            [0., 1., 0.],
        )
    }
    fn wall_x(x: f64, y0: f64, y1: f64, z0: f64, z1: f64) -> Panel {
        (
            vec![vec![[x, y0, z0], [x, y1, z0], [x, y1, z1], [x, y0, z1]]],
            [1., 0., 0.],
        )
    }

    fn run(model: &mut Model, place: &Placement, interior: Option<Vec<BTreeSet<usize>>>) -> Report {
        let interior = interior.unwrap_or_else(|| no_interior(model));
        let locked = BTreeSet::new();
        let j = junctions::insert(
            model,
            &junctions::Context {
                interior: &interior,
                locked: &locked,
            },
        );
        assert!(j.issues.is_empty(), "{:?}", j.issues);
        trim(
            model,
            &Context {
                maximum_width: 0.25 * place.scale,
                interior: &interior,
                locked: &locked,
                axes: &[],
                contacts: &[],
            },
        )
    }

    fn area(model: &Model, s: usize) -> f64 {
        let ring = &model.surfaces[s].contours[0];
        (0..ring.len())
            .map(|i| {
                let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
                a[0] * b[1] - b[0] * a[1]
            })
            .sum::<f64>()
            .abs()
            / 2.
    }

    /// The whole segment `p`-`q` is covered by edges of both surfaces.
    fn shared(model: &Model, place: &Placement, a: usize, b: usize, p: [f64; 3], q: [f64; 3]) {
        let (p, q) = (
            DVec3::from_array(place.point(p)),
            DVec3::from_array(place.point(q)),
        );
        let d = (q - p).normalize();
        let eps = model.precision * 5.;
        let mut parts: Vec<(f64, f64)> = model
            .surface_edges(a)
            .filter(|&e| model.surface_edges(b).any(|f| f == e))
            .filter_map(|e| {
                let [i, j] = model.edges[e].map(|v| DVec3::from_array(model.vertices[v]));
                let (ti, tj) = ((i - p).dot(d), (j - p).dot(d));
                (i.distance(p + d * ti) <= eps && j.distance(p + d * tj) <= eps)
                    .then_some((ti.min(tj), ti.max(tj)))
            })
            .collect();
        parts.sort_by(|x, y| x.0.total_cmp(&y.0));
        let mut cursor = 0.;
        for (s, e) in parts {
            assert!(s <= cursor + eps, "gap at {cursor}");
            cursor = f64::max(cursor, e);
        }
        assert!(cursor >= p.distance(q) - eps);
    }

    #[test]
    fn console_beyond_supporting_wall_is_trimmed_to_its_axis() {
        for place in Placement::all() {
            for (z0, z1) in [(-2., 0.), (-2., 2.)] {
                // Wall below the slab (T-junction) or passing through it.
                let mut m = build(&place, &[slab(0., 4.), wall_y(3.85, 0., 4., z0, z1)]);
                let n = m.surfaces.len();
                let r = run(&mut m, &place, None);
                assert_eq!(r.trimmed.len(), 1, "{:?}", r.kept);
                let s2 = place.scale * place.scale;
                assert!((r.trimmed[0].area - 0.6 * s2).abs() < 1e-9 * s2);
                assert!((r.trimmed[0].width - 0.15 * place.scale).abs() < 1e-9 * place.scale);
                assert!((area(&m, 0) - 4. * 3.85 * s2).abs() < 1e-9 * s2);
                assert!(m.surfaces[0].embedded_edges.is_empty());
                assert_eq!(m.surfaces.len(), n);
                shared(&m, &place, 0, 1, [0., 3.85, 0.], [4., 3.85, 0.]);
                // Idempotent.
                let again = run(&mut m, &place, None);
                assert!(again.trimmed.is_empty());
            }
        }
    }

    #[test]
    fn wide_or_attached_consoles_are_kept() {
        for place in Placement::all() {
            // 0.30 wide: wider than the limit.
            let mut m = build(&place, &[slab(0., 4.), wall_y(3.7, 0., 4., -2., 0.)]);
            assert!(run(&mut m, &place, None).trimmed.is_empty());
            // A parapet on the free edge: the console carries a structure.
            let mut m = build(
                &place,
                &[
                    slab(0., 4.),
                    wall_y(3.85, 0., 4., -2., 0.),
                    wall_y(4., 0., 4., 0., 1.),
                ],
            );
            let r = run(&mut m, &place, None);
            assert!(r.trimmed.is_empty());
            assert!(r.kept.contains_key("contour_shared_with_other_surface"));
            // A retained node (e.g. a column) inside the strip.
            let mut m = build(&place, &[slab(0., 4.), wall_y(3.85, 0., 4., -2., 0.)]);
            let node = m.add_vertex(place.point([2., 3.95, 0.])).unwrap();
            let mut interior = no_interior(&m);
            interior[0].insert(node);
            let r = run(&mut m, &place, Some(interior));
            assert!(r.trimmed.is_empty());
            assert!(r.kept.contains_key("retained_node_inside"));
        }
    }

    #[test]
    fn openings_and_floating_lines_survive_a_trim() {
        for place in Placement::all() {
            let mut panel = slab(0., 4.);
            panel
                .0
                .push(vec![[1., 1., 0.], [1., 2., 0.], [2., 2., 0.], [2., 1., 0.]]);
            // A short wall crossing the slab without touching other lines.
            let mut m = build(
                &place,
                &[
                    panel,
                    wall_y(3.85, 0., 4., -2., 0.),
                    wall_x(3., 0.5, 2.5, -1., 1.),
                ],
            );
            let r = run(&mut m, &place, None);
            assert_eq!(r.trimmed.len(), 1, "{:?}", r.kept);
            assert_eq!(m.surfaces[0].contours.len(), 2);
            assert_eq!(m.surfaces[0].embedded_edges.len(), 1);
            shared(&m, &place, 0, 2, [3., 0.5, 0.], [3., 2.5, 0.]);
            shared(&m, &place, 0, 1, [0., 3.85, 0.], [4., 3.85, 0.]);
        }
    }

    #[test]
    fn perimeter_band_around_closed_wall_loop_is_trimmed() {
        for place in Placement::all() {
            // Four walls under the slab meeting at their axes, 0.15 inside
            // the slab edge on every side: the console is an annulus.
            let (a, b) = (0.15, 3.85);
            let mut m = build(
                &place,
                &[
                    slab(0., 4.),
                    wall_y(a, a, b, -2., 0.),
                    wall_y(b, a, b, -2., 0.),
                    wall_x(a, a, b, -2., 0.),
                    wall_x(b, a, b, -2., 0.),
                ],
            );
            let r = run(&mut m, &place, None);
            assert_eq!(r.trimmed.len(), 1, "{:?}", r.kept);
            let s2 = place.scale * place.scale;
            assert!((area(&m, 0) - 3.7 * 3.7 * s2).abs() < 1e-9 * s2);
            assert_eq!(m.surfaces[0].contours.len(), 1);
            assert!(m.surfaces[0].embedded_edges.is_empty());
            shared(&m, &place, 0, 1, [a, a, 0.], [b, a, 0.]);
            shared(&m, &place, 0, 4, [b, a, 0.], [b, b, 0.]);
        }
    }

    #[test]
    fn corner_consoles_and_wall_stubs_are_trimmed_together() {
        for place in Placement::all() {
            // Two perpendicular walls under the slab, both reaching the slab
            // edge: slab consoles, a corner square and two wall stubs.
            let mut m = build(
                &place,
                &[
                    slab(0., 4.),
                    wall_y(3.85, 0., 4., -2., 0.),
                    wall_x(3.85, 0., 4., -2., 0.),
                ],
            );
            let r = run(&mut m, &place, None);
            assert!(!r.trimmed.is_empty(), "{:?}", r.kept);
            let s2 = place.scale * place.scale;
            assert!((area(&m, 0) - 3.85 * 3.85 * s2).abs() < 1e-9 * s2);
            // The stubs beyond the other wall are trimmed in a later pass.
            assert!((area(&m, 1) - 3.85 * 2. * s2).abs() < 1e-9 * s2);
            assert!((area(&m, 2) - 3.85 * 2. * s2).abs() < 1e-9 * s2);
            shared(&m, &place, 0, 1, [0., 3.85, 0.], [3.85, 3.85, 0.]);
            shared(&m, &place, 0, 2, [3.85, 0., 0.], [3.85, 3.85, 0.]);
            shared(&m, &place, 1, 2, [3.85, 3.85, -2.], [3.85, 3.85, 0.]);
        }
    }
}
