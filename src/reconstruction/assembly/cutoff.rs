//! The upper storeys cut off at a floor level (the model of a tall building
//! reduced to its lower part). The result is an ordinary state: surfaces
//! crossing the level are clipped along it (new faces share their edges with
//! the neighbours; an edge crossing the level gets one new vertex shared by
//! everything that uses it), bars are trimmed, what lies above is removed
//! with provenance. What the removed part rested on is kept in `Cut`: the
//! lines of the walls and the points of the columns at the level, where its
//! loads and weight are applied by the export (`loads`).
use super::bars::{Anchor, Axis, NO_SOURCE_NODE};
use super::edit::{RemovedBar, RemovedSurface, State};
use crate::reconstruction::Model;
use geo::{Contains, InteriorPoint, LineString, Point, Polygon};
use glam::{DVec2, DVec3};
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};

/// A line (walls: `a` to `b`) or a point (columns: `a == b`) at the level of
/// the cut on which a removed part stood.
#[derive(Debug, Clone, Serialize, serde::Deserialize, PartialEq)]
pub struct Support {
    /// Vertices of the geometry (the support follows edits that move or
    /// merge them): the ends of a wall's line, one vertex twice for a column.
    pub a: usize,
    pub b: usize,
    /// LIRA stiffness number of the wall or column.
    pub stiffness: u32,
}

/// What the cut left behind.
#[derive(Debug, Clone, Serialize, serde::Deserialize, PartialEq)]
pub struct Cut {
    pub z: f64,
    /// Highest elevation of what was cut off.
    pub top: f64,
    pub walls: Vec<Support>,
    pub columns: Vec<Support>,
    pub removed_surfaces: usize,
    pub clipped_surfaces: usize,
    pub removed_bars: usize,
    pub trimmed_bars: usize,
}

/// A floor: the horizontal surfaces at one elevation.
#[derive(Debug, Clone, Serialize)]
pub struct Floor {
    pub z: f64,
    pub slabs: usize,
    pub area: f64,
    /// A real floor (a good part of the largest one), not a landing or a pit.
    pub major: bool,
}

fn ring_area(r: &[[f64; 2]]) -> f64 {
    (0..r.len())
        .map(|i| r[i][0] * r[(i + 1) % r.len()][1] - r[(i + 1) % r.len()][0] * r[i][1])
        .sum::<f64>()
        / 2.
}

/// The floors of the model, lowest first: groups of horizontal surfaces
/// within `0.15` m of one elevation.
pub fn floors(state: &State) -> Vec<Floor> {
    let model = &state.model;
    let mut slabs: Vec<(f64, f64)> = vec![];
    for (s, surface) in model.surfaces().iter().enumerate() {
        if model.planes()[surface.plane].normal()[2].abs() < 0.98 {
            continue;
        }
        let zs: Vec<f64> = model
            .surface_edges(s)
            .flat_map(|e| model.edges()[e])
            .map(|v| model.vertices()[v][2])
            .collect();
        let (lo, hi) = zs.iter().fold((f64::MAX, f64::MIN), |(l, h), z| (l.min(*z), h.max(*z)));
        if hi - lo > 0.05 {
            continue;
        }
        let area = ring_area(&surface.contours[0]).abs() - surface.contours[1..].iter().map(|h| ring_area(h).abs()).sum::<f64>();
        slabs.push(((lo + hi) / 2., area));
    }
    slabs.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut out: Vec<Floor> = vec![];
    for (z, area) in slabs {
        match out.last_mut() {
            Some(f) if (z - f.z).abs() < 0.15 => {
                f.z = (f.z * f.slabs as f64 + z) / (f.slabs as f64 + 1.);
                f.slabs += 1;
                f.area += area;
            }
            _ => out.push(Floor { z, slabs: 1, area, major: false }),
        }
    }
    let largest = out.iter().map(|f| f.area).fold(0., f64::max);
    for f in &mut out {
        f.major = f.area >= 0.2 * largest;
    }
    out
}

/// One face of a clipped surface.
struct Piece {
    rings: Vec<Vec<usize>>,
    embedded: Vec<usize>,
}

impl State {
    /// Cut off everything above elevation `z`: the floor slabs at `z` stay
    /// (they close the lower part).
    pub fn cut_above(&mut self, z: f64) -> Result<String, String> {
        if !z.is_finite() {
            return Err("cut_level".into());
        }
        // Vertices within the shortest edge of the level are on it (a floor's
        // vertices differ by millimetres); an edge crossing farther gets a vertex.
        let minimum = self.model.minimum_edge();
        let share = (self.model.precision() * 10.).max(1e-6);
        let tol = minimum.max(self.model.precision() * 10.).max(1e-6);
        let mut model = self.model.clone();
        let height = |model: &Model, v: usize| model.vertices()[v][2] - z;
        if !model.vertices().iter().any(|p| p[2] > z + tol) {
            return Err("cut_nothing_above".into());
        }

        // 1. An edge crossing the level gets a vertex there.
        let used: BTreeSet<usize> = (0..model.surfaces().len()).flat_map(|s| model.surface_edges(s).collect::<Vec<_>>()).collect();
        for e in used {
            let [a, b] = model.edges()[e];
            let (za, zb) = (height(&model, a), height(&model, b));
            if (za < -tol && zb > tol) || (za > tol && zb < -tol) {
                let (pa, pb) = (DVec3::from_array(model.vertices()[a]), DVec3::from_array(model.vertices()[b]));
                let mut p = pa + (pb - pa) * (za / (za - zb));
                p.z = z;
                if p.distance(pa) < minimum || p.distance(pb) < minimum {
                    return Err(format!("cut_near_vertex: the edge {a}-{b} crosses the level {:.3} m from its end; cut at a floor level", p.distance(pa).min(p.distance(pb))));
                }
                let v = model.add_vertex(p.to_array()).map_err(|e| format!("cut_{e:?}"))?;
                model.split_edge(e, v).map_err(|e| format!("cut_split_{e:?}"))?;
            }
        }
        let label = |model: &Model, v: usize| -> i8 {
            let h = height(model, v);
            if h > tol { 1 } else if h < -tol { -1 } else { 0 }
        };

        // 2. Surfaces: kept, removed or clipped.
        let mut remove: BTreeSet<usize> = BTreeSet::new();
        let mut removed_records: Vec<RemovedSurface> = vec![];
        let mut added: Vec<(usize, Piece)> = vec![]; // (original surface, face)
        let mut walls: Vec<Support> = vec![];
        let (mut removed_surfaces, mut clipped_surfaces) = (0, 0);
        let mut removed_top = z;
        for s in 0..model.surfaces().len() {
            let vertices: BTreeSet<usize> = model.surface_edges(s).flat_map(|e| model.edges()[e]).collect();
            let labels: Vec<i8> = vertices.iter().map(|&v| label(&model, v)).collect();
            let (below, above) = (labels.iter().any(|&l| l < 0), labels.iter().any(|&l| l > 0));
            if !above {
                continue;
            }
            let vertical = model.planes()[model.surfaces()[s].plane].normal()[2].abs() < 0.5;
            remove.insert(s);
            for &v in &vertices {
                removed_top = removed_top.max(model.vertices()[v][2]);
            }
            if !below {
                // Entirely above: its edges on the level are where it stood.
                removed_surfaces += 1;
                let surface = &model.surfaces()[s];
                removed_records.push(RemovedSurface {
                    reason: "cut".into(),
                    stiffness: self.stiffness[s],
                    patch: self.patches[s],
                    source_elements: surface.source_elements.clone(),
                    into: None,
                });
                if vertical {
                    for ring in &surface.boundaries {
                        for u in ring {
                            let [a, b] = model.edges()[u.edge];
                            if label(&model, a) == 0 && label(&model, b) == 0 {
                                walls.push(Support { a, b, stiffness: self.stiffness[s] });
                            }
                        }
                    }
                }
                continue;
            }
            clipped_surfaces += 1;
            let faces = clip(&model, s, z, tol, &label).map_err(|e| format!("cut_surface_{s}: {e}"))?;
            if vertical {
                // The chords: where the wall goes on above the level.
                for face in &faces.1 {
                    walls.push(Support { a: face.0, b: face.1, stiffness: self.stiffness[s] });
                }
            }
            for piece in faces.0 {
                added.push((s, piece));
            }
        }
        // Support lines found twice (a chord and the edge of a wall above) once.
        let mut seen: BTreeSet<[usize; 3]> = BTreeSet::new();
        walls.retain(|w| seen.insert([w.a.min(w.b), w.a.max(w.b), w.stiffness as usize]));

        // New faces first, then the old surfaces go.
        let mut new_data: Vec<(u32, usize)> = vec![];
        for (s, piece) in added {
            let plane = model.surfaces()[s].plane;
            let elements = model.surfaces()[s].source_elements.clone();
            let id = model.add_surface(plane, piece.rings.clone(), elements).map_err(|e| format!("cut_face_of_{s}: {e:?}"))?;
            if !piece.embedded.is_empty() {
                model.rebuild_surface(id, piece.rings, piece.embedded).map_err(|e| format!("cut_face_of_{s}: {e:?}"))?;
            }
            new_data.push((self.stiffness[s], self.patches[s]));
        }
        let mut stiffness = self.stiffness.clone();
        let mut patches = self.patches.clone();
        for (st, pa) in new_data {
            stiffness.push(st);
            patches.push(pa);
        }
        let index = model.remove_surfaces(&remove);

        // 3. Bars.
        let mut axes: Vec<Axis> = vec![];
        let mut removed_bars: Vec<RemovedBar> = vec![];
        let mut columns: Vec<Support> = vec![];
        let (mut dropped, mut trimmed) = (0, 0);
        for axis in &self.axes {
            let [e0, e1] = axis.endpoints;
            let (l0, l1) = (label(&model, e0), label(&model, e1));
            removed_top = removed_top.max(model.vertices()[e0][2]).max(model.vertices()[e1][2]);
            if l0 <= 0 && l1 <= 0 {
                axes.push(axis.clone());
                continue;
            }
            let stiffness_at = |t: f64| {
                axis.spans
                    .iter()
                    .find(|s| s.start_t.min(s.end_t) - 1e-9 <= t && t <= s.start_t.max(s.end_t) + 1e-9)
                    .or_else(|| axis.spans.first())
                    .map_or(0, |s| s.stiffness)
            };
            if l0 >= 0 && l1 >= 0 {
                // Entirely above: a column stood on the level where an end is on it.
                dropped += 1;
                let mut elements: Vec<u32> = axis.spans.iter().map(|s| s.element).collect();
                elements.sort_unstable();
                elements.dedup();
                removed_bars.push(RemovedBar { reason: "cut".into(), source_axis: axis.source_axis, removed: true, source_elements: elements });
                for (k, l) in [(0, l0), (1, l1)] {
                    if l == 0 {
                        let v = axis.endpoints[k];
                        columns.push(Support { a: v, b: v, stiffness: stiffness_at(k as f64) });
                    }
                }
                continue;
            }
            // Crossing: keep the part below.
            trimmed += 1;
            let (p0, p1) = (DVec3::from_array(model.vertices()[e0]), DVec3::from_array(model.vertices()[e1]));
            let tc = (z - p0.z) / (p1.z - p0.z);
            let length = p0.distance(p1);
            if tc * length < minimum || (1. - tc) * length < minimum {
                return Err(format!("cut_near_vertex: bar {} crosses the level {:.3} m from its end; cut at a floor level", axis.source_axis, (tc * length).min((1. - tc) * length)));
            }
            let near = axis.anchors.iter().find(|a| a.vertex != e0 && a.vertex != e1 && (a.t - tc).abs() * length < minimum);
            let (tc, vertex) = match near {
                Some(a) => (a.t, a.vertex),
                None => {
                    let mut p = p0 + (p1 - p0) * tc;
                    p.z = z;
                    // A vertex of the model already there (the split of a wall edge the bar runs
                    // along, the bar of another storey) is the node: the topology stays shared.
                    let shared = (0..model.vertices().len()).find(|&v| DVec3::from_array(model.vertices()[v]).distance(p) <= share);
                    match shared {
                        Some(v) => (tc, v),
                        None => (tc, model.add_vertex(p.to_array()).map_err(|e| format!("cut_{e:?}"))?),
                    }
                }
            };
            let (lo, hi, ends) = if l0 < 0 { (0., tc, [e0, vertex]) } else { (tc, 1., [vertex, e1]) };
            let remap = |t: f64| (t - lo) / (hi - lo);
            let mut anchors: Vec<Anchor> = axis
                .anchors
                .iter()
                .filter(|a| a.t >= lo - 1e-12 && a.t <= hi + 1e-12 && a.vertex != vertex)
                .map(|a| Anchor { source_node: a.source_node, vertex: a.vertex, t: remap(a.t).clamp(0., 1.) })
                .collect();
            anchors.push(Anchor { source_node: NO_SOURCE_NODE, vertex, t: if l0 < 0 { 1. } else { 0. } });
            anchors.sort_by(|a, b| a.t.total_cmp(&b.t));
            let spans = axis
                .spans
                .iter()
                .filter_map(|s| {
                    let (start, end) = (s.start_t.max(lo).min(hi), s.end_t.max(lo).min(hi));
                    ((end - start).abs() * length > 1e-9).then(|| {
                        let mut c = s.clone();
                        c.start_t = remap(start);
                        c.end_t = remap(end);
                        c
                    })
                })
                .collect();
            let above_end = if l0 < 0 { 1. } else { 0. };
            columns.push(Support { a: vertex, b: vertex, stiffness: stiffness_at(above_end) });
            axes.push(Axis { source_axis: axis.source_axis, endpoints: ends, anchors, spans });
        }

        // Done: the state follows.
        self.model = model;
        self.stiffness = stiffness;
        self.patches = patches;
        self.renumber_after_cut(&index);
        self.removed.extend(removed_records);
        self.removed_bars.extend(removed_bars);
        self.axes = axes;
        self.refresh_after_cut();
        let top = removed_top;
        self.cut = Some(Cut { z, top, walls, columns, removed_surfaces, clipped_surfaces, removed_bars: dropped, trimmed_bars: trimmed });
        Ok(format!(
            "cut at {z:.3}: {removed_surfaces} surfaces removed, {clipped_surfaces} clipped, {dropped} bars removed, {trimmed} trimmed"
        ))
    }
}

/// The faces of surface `s` below the level, and the chords along the level
/// (where the surface goes on above).
#[allow(clippy::type_complexity)]
fn clip(
    model: &Model,
    s: usize,
    z: f64,
    tol: f64,
    label: &dyn Fn(&Model, usize) -> i8,
) -> Result<(Vec<Piece>, Vec<(usize, usize)>), String> {
    let surface = model.surfaces()[s].clone();
    let plane = &model.planes()[surface.plane];
    let uv = |v: usize| DVec2::from_array(plane.project(model.vertices()[v]));
    let key = |a: usize, b: usize| [a.min(b), a.max(b)];
    let rings: Vec<Vec<usize>> = surface
        .boundaries
        .iter()
        .map(|ring| {
            ring.iter()
                .map(|u| {
                    let [a, b] = model.edges()[u.edge];
                    if u.reversed { b } else { a }
                })
                .collect()
        })
        .collect();
    // Graph of the edges below or on the level; holes entirely below are loose.
    let mut edges: BTreeSet<[usize; 2]> = BTreeSet::new();
    let mut loose: Vec<usize> = vec![];
    for (i, ring) in rings.iter().enumerate() {
        if i > 0 && ring.iter().all(|&v| label(model, v) < 0) {
            loose.push(i);
            continue;
        }
        for k in 0..ring.len() {
            let (a, b) = (ring[k], ring[(k + 1) % ring.len()]);
            if label(model, a) <= 0 && label(model, b) <= 0 {
                edges.insert(key(a, b));
            }
        }
    }
    // The cut line in the plane and the vertices on it.
    let n = glam::DVec3::from_array(plane.normal());
    let dir3 = n.cross(DVec3::Z);
    let o = plane.origin();
    let (po, pd) = (plane.project(o), plane.project((DVec3::from_array(o) + dir3).to_array()));
    let dir = (DVec2::new(pd[0] - po[0], pd[1] - po[1])).normalize();
    let on: BTreeSet<usize> = rings
        .iter()
        .flatten()
        .copied()
        .chain(surface.embedded_edges.iter().flat_map(|&e| model.edges()[e]))
        .filter(|&v| label(model, v) == 0)
        .collect();
    let mut on: Vec<usize> = on.into_iter().collect();
    on.sort_by(|&a, &b| uv(a).dot(dir).total_cmp(&uv(b).dot(dir)));
    let tuple = |r: &[[f64; 2]]| LineString::from(r.iter().map(|p| (p[0], p[1])).collect::<Vec<_>>());
    let material = Polygon::new(tuple(&surface.contours[0]), surface.contours[1..].iter().map(|h| tuple(h)).collect());
    let mut chords: Vec<[usize; 2]> = vec![];
    for w in on.windows(2) {
        let (a, b) = (w[0], w[1]);
        if (uv(a) - uv(b)).length() < model.minimum_edge() * 0.5 {
            continue;
        }
        let m = (uv(a) + uv(b)) / 2.;
        let inside = material.contains(&Point::new(m.x, m.y));
        if inside {
            chords.push(key(a, b));
            edges.insert(key(a, b));
        }
    }
    if edges.is_empty() {
        return Err("nothing of the surface is below the level".into());
    }
    // Faces of the planar graph.
    let mut neighbours: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for [a, b] in &edges {
        neighbours.entry(*a).or_default().push(*b);
        neighbours.entry(*b).or_default().push(*a);
    }
    for (v, list) in neighbours.iter_mut() {
        let p = uv(*v);
        list.sort_by(|&x, &y| {
            let (dx, dy) = (uv(x) - p, uv(y) - p);
            dx.y.atan2(dx.x).total_cmp(&dy.y.atan2(dy.x))
        });
    }
    let mut visited: BTreeSet<(usize, usize)> = BTreeSet::new();
    let mut faces: Vec<Vec<usize>> = vec![];
    for (&u0, list) in &neighbours {
        for &v0 in list {
            if visited.contains(&(u0, v0)) {
                continue;
            }
            let (mut u, mut v) = (u0, v0);
            let mut cycle = vec![];
            loop {
                visited.insert((u, v));
                cycle.push(u);
                let around = &neighbours[&v];
                let pos = around.iter().position(|&x| x == u).ok_or("broken planar graph")?;
                let w = around[(pos + around.len() - 1) % around.len()];
                (u, v) = (v, w);
                if (u, v) == (u0, v0) || cycle.len() > edges.len() * 2 + 4 {
                    break;
                }
            }
            faces.push(cycle);
        }
    }
    // Faces of material (counter-clockwise, with material at an inner point).
    let mut pieces: Vec<(Vec<usize>, Polygon<f64>, f64)> = vec![];
    let (mut positive, mut non_material, face_count) = (0, 0, faces.len());
    let mut notes: Vec<String> = vec![];
    for cycle in faces {
        let pts: Vec<[f64; 2]> = cycle.iter().map(|&v| uv(v).to_array()).collect();
        let area = ring_area(&pts);
        if area <= model.precision() * model.precision() || cycle.len() < 3 {
            continue;
        }
        positive += 1;
        let polygon = Polygon::new(tuple(&pts), vec![]);
        // A point of the face away from the holes that lie entirely below the level.
        let holes: Vec<LineString<f64>> = loose
            .iter()
            .filter(|&&i| {
                let q = uv(rings[i][0]);
                polygon.contains(&Point::new(q.x, q.y))
            })
            .map(|&i| tuple(&rings[i].iter().map(|&v| uv(v).to_array()).collect::<Vec<_>>()))
            .collect();
        let Some(p) = Polygon::new(tuple(&pts), holes).interior_point() else { continue };
        if !material.contains(&p) {
            non_material += 1;
            notes.push(format!("cycle of {} vertices, area {area:.3}, inner point ({:.3}, {:.3}), exterior area {:.3}, holes {}", cycle.len(), p.x(), p.y(), ring_area(&surface.contours[0]), surface.contours.len() - 1));
            continue;
        }
        let distinct: BTreeSet<usize> = cycle.iter().copied().collect();
        if distinct.len() != cycle.len() {
            return Err("a face of the cut touches itself (a hole touching the level): cut at another level or edit the surface".into());
        }
        pieces.push((cycle, polygon, area));
    }
    if pieces.is_empty() {
        return Err(format!("no face below the level ({} edges, {face_count} cycles, {positive} counter-clockwise, {non_material} outside the material, {} chords, {} vertices on the level; {notes:?})", edges.len(), chords.len(), on.len()));
    }
    let mut out: Vec<Piece> = pieces.iter().map(|p| Piece { rings: vec![p.0.clone()], embedded: vec![] }).collect();
    // Loose holes go to the smallest face around them, clockwise.
    for i in loose {
        let first = uv(rings[i][0]);
        let host = (0..pieces.len())
            .filter(|&k| pieces[k].1.contains(&Point::new(first.x, first.y)))
            .min_by(|&a, &b| pieces[a].2.total_cmp(&pieces[b].2));
        if let Some(k) = host {
            let mut ring = rings[i].clone();
            let pts: Vec<[f64; 2]> = ring.iter().map(|&v| uv(v).to_array()).collect();
            if ring_area(&pts) > 0. {
                ring.reverse();
            }
            out[k].rings.push(ring);
        }
    }
    // Embedded edges stay with the face around them (not when they became its boundary).
    for &e in &surface.embedded_edges {
        let [a, b] = model.edges()[e];
        if label(model, a) > 0 || label(model, b) > 0 {
            continue;
        }
        let m = (uv(a) + uv(b)) / 2.;
        for (k, p) in pieces.iter().enumerate() {
            let boundary = p.0.iter().enumerate().any(|(i, &v)| key(v, p.0[(i + 1) % p.0.len()]) == key(a, b));
            if !boundary && p.1.contains(&Point::new(m.x, m.y)) {
                out[k].embedded.push(e);
            }
        }
    }
    let chord_lines = chords.iter().map(|[a, b]| (*a, *b)).collect();
    let _ = (z, tol);
    Ok((out, chord_lines))
}

/// The supports of the cut that still exist: a wall's line needs an edge of a
/// surface along it, a column a bar node at its vertex. Positions follow the
/// vertices (edits move them). Lines: (start, end, stiffness); columns:
/// (point, stiffness).
#[allow(clippy::type_complexity)]
pub fn live_supports(state: &State) -> (Vec<(DVec3, DVec3, u32)>, Vec<(DVec3, u32)>) {
    let Some(cut) = state.cut.as_ref() else { return (vec![], vec![]) };
    let model = &state.model;
    let vertices = model.vertices();
    let tol = (model.precision() * 10.).max(1e-6);
    let mut lines = vec![];
    let used: Vec<usize> = (0..model.surfaces().len()).flat_map(|s| model.surface_edges(s).collect::<Vec<_>>()).collect();
    for w in &cut.walls {
        let (Some(a), Some(b)) = (vertices.get(w.a), vertices.get(w.b)) else { continue };
        let (a, b) = (DVec3::from_array(*a), DVec3::from_array(*b));
        let d = b - a;
        if d.length() < 1e-9 {
            continue;
        }
        let on = |p: DVec3| {
            let t = (p - a).dot(d) / d.length_squared();
            (-1e-9..=1. + 1e-9).contains(&t) && (p - a - d * t).length() <= tol
        };
        if used.iter().any(|&e| model.edges()[e].iter().all(|&v| on(DVec3::from_array(vertices[v])))) {
            lines.push((a, b, w.stiffness));
        }
    }
    let mut points = vec![];
    for c in &cut.columns {
        let Some(p) = vertices.get(c.a) else { continue };
        if state.axes.iter().any(|axis| axis.endpoints.contains(&c.a) || axis.anchors.iter().any(|n| n.vertex == c.a)) {
            points.push((DVec3::from_array(*p), c.stiffness));
        }
    }
    (lines, points)
}

impl State {
    /// A vertex merged into another one: the supports of the cut follow.
    pub(super) fn remap_cut_vertex(&mut self, drop: usize, keep: usize) {
        if let Some(cut) = self.cut.as_mut() {
            for s in cut.walls.iter_mut().chain(cut.columns.iter_mut()) {
                if s.a == drop {
                    s.a = keep;
                }
                if s.b == drop {
                    s.b = keep;
                }
            }
            cut.walls.retain(|w| w.a != w.b);
        }
    }
}
