//! Connectivity of the assembled structure: groups of bars and surfaces
//! hanging in the air, bar ends connected to nothing, and links the source
//! model had (bars or a bar and a shell sharing a node) that the
//! reconstruction lost. Every finding carries repair proposals, edits the
//! editor can apply as they are.
use super::{closest_on3, closest_points, p3, Class, Finding, Fix, Grid, Options, Surface};
use crate::reconstruction::assembly::bars::{Axis, Contact};
use glam::DVec3;
use serde::Serialize;
use serde_json::json;
use std::collections::{BTreeMap, BTreeSet};

/// Two source bars (constructive segments) sharing a node.
#[derive(Debug, Clone, Serialize)]
pub struct BarLink {
    /// Source axes (frame axis indices).
    pub a: usize,
    pub b: usize,
    pub node: u32,
    pub at: [f64; 3],
}

/// A source bar and a shell patch sharing a node.
#[derive(Debug, Clone, Serialize)]
pub struct SurfaceLink {
    pub axis: usize,
    pub patch: usize,
    pub node: u32,
    pub at: [f64; 3],
}

/// What the source model connected, for the comparison with the result.
#[derive(Debug, Clone, Default, Serialize)]
pub struct SourceLinks {
    pub bar_bar: Vec<BarLink>,
    pub bar_surface: Vec<SurfaceLink>,
}

/// Source links of a frame: bars through one node, a bar node on a shell.
pub fn source_links(frame: &crate::reconstruction::frame::Report) -> SourceLinks {
    let mut bars_at: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (k, axis) in frame.axes.iter().enumerate() {
        for anchor in &axis.anchors {
            bars_at.entry(anchor.node).or_default().push(k);
        }
    }
    let mut patches_at: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (p, surface) in frame.surfaces.iter().enumerate() {
        for &n in &surface.nodes {
            if bars_at.contains_key(&n) {
                patches_at.entry(n).or_default().push(p);
            }
        }
    }
    let mut links = SourceLinks::default();
    for (&n, bars) in &bars_at {
        let at = frame.reference_points[n];
        let node = frame.node_ids[n];
        // At most the first few bars of a node pair up: a node of 20 bars
        // is one connection, not 190.
        let limited = &bars[..bars.len().min(6)];
        for (i, &a) in limited.iter().enumerate() {
            for &b in &limited[i + 1..] {
                if a != b {
                    links.bar_bar.push(BarLink { a, b, node, at });
                }
            }
        }
        for &p in patches_at.get(&n).into_iter().flatten() {
            for &axis in limited {
                links.bar_surface.push(SurfaceLink {
                    axis,
                    patch: p,
                    node,
                    at,
                });
            }
        }
    }
    links
}

/// What the audit knows besides the geometry.
pub struct Context<'a> {
    pub links: Option<&'a SourceLinks>,
    /// Source patch of every surface.
    pub patches: &'a [usize],
}

struct Union(Vec<usize>);
impl Union {
    fn find(&mut self, mut x: usize) -> usize {
        while self.0[x] != x {
            self.0[x] = self.0[self.0[x]];
            x = self.0[x];
        }
        x
    }
    fn join(&mut self, a: usize, b: usize) {
        let (a, b) = (self.find(a), self.find(b));
        self.0[a] = b;
    }
}

struct Bar {
    p0: DVec3,
    p1: DVec3,
    nodes: BTreeSet<usize>,
    ends: [usize; 2],
}

/// Largest number of objects listed in one finding.
const LISTED: usize = 200;

pub(super) fn analyse(
    surfaces: &[Surface],
    axes: &[Axis],
    contacts: &[Contact],
    vertices: &[[f64; 3]],
    eps: f64,
    options: &Options,
    context: &Context<'_>,
    findings: &mut Vec<Finding>,
) {
    let bars: Vec<Bar> = axes
        .iter()
        .map(|a| Bar {
            p0: p3(vertices[a.endpoints[0]]),
            p1: p3(vertices[a.endpoints[1]]),
            nodes: a
                .anchors
                .iter()
                .map(|n| n.vertex)
                .chain(a.endpoints)
                .collect(),
            ends: a.endpoints,
        })
        .collect();
    let (b_count, s_count) = (bars.len(), surfaces.len());
    let reach = (options.element_size * 0.6).max(eps);
    let group_reach = options.element_size * 4.;

    // ---- components over bars and surfaces
    let mut union = Union((0..b_count + s_count).collect());
    let mut at_vertex: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (i, bar) in bars.iter().enumerate() {
        for &v in &bar.nodes {
            at_vertex.entry(v).or_default().push(i);
        }
    }
    let mut surface_at: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (s, surface) in surfaces.iter().enumerate() {
        for &v in &surface.vertex_ids {
            surface_at.entry(v).or_default().push(s);
        }
    }
    for (v, list) in &at_vertex {
        for &b in &list[1..] {
            union.join(list[0], b);
        }
        for &s in surface_at.get(v).into_iter().flatten() {
            union.join(list[0], b_count + s);
        }
    }
    for list in surface_at.values() {
        for &s in &list[1..] {
            union.join(b_count + list[0], b_count + s);
        }
    }
    for c in contacts {
        let (Contact::Point { axis, surface, .. } | Contact::Interval { axis, surface, .. }) = *c;
        if axis < b_count && surface < s_count {
            union.join(axis, b_count + surface);
        }
    }
    let component: Vec<usize> = (0..b_count + s_count).map(|x| union.find(x)).collect();
    let mut sizes: BTreeMap<usize, (usize, usize)> = BTreeMap::new();
    for (x, &c) in component.iter().enumerate() {
        let e = sizes.entry(c).or_default();
        if x < b_count {
            e.0 += 1;
        } else {
            e.1 += 1;
        }
    }
    let main = sizes
        .iter()
        .max_by_key(|(c, (b, s))| (b + s, std::cmp::Reverse(**c)))
        .map(|(&c, _)| c);

    // ---- bar index for proposals
    let extent = bars
        .iter()
        .map(|b| (b.p1 - b.p0).abs().max_element())
        .fold(0., f64::max);
    let mut grid = Grid::new((group_reach * 2.).max(extent.min(10.)).max(1e-6));
    for (i, b) in bars.iter().enumerate() {
        grid.insert(i, b.p0.min(b.p1), b.p0.max(b.p1));
    }
    let nearest_bar =
        |p: DVec3, skip: &dyn Fn(usize) -> bool, radius: f64| -> Option<(f64, usize, DVec3)> {
            let mut best: Option<(f64, usize, DVec3)> = None;
            for j in grid.query(p - DVec3::splat(radius), p + DVec3::splat(radius)) {
                if skip(j) {
                    continue;
                }
                let q = closest_on3(p, bars[j].p0, bars[j].p1);
                let d = p.distance(q);
                if d <= radius && best.is_none_or(|b| d < b.0) {
                    best = Some((d, j, q));
                }
            }
            best
        };
    let nearest_surface = |p: DVec3, radius: f64| -> Option<(f64, usize, DVec3)> {
        let mut best: Option<(f64, usize, DVec3)> = None;
        for (s, surface) in surfaces.iter().enumerate() {
            if (p + radius).cmplt(surface.lo).any() || (p - radius).cmpgt(surface.hi).any() {
                continue;
            }
            let d = surface.distance(p);
            if d <= radius && best.is_none_or(|b| d < b.0) {
                let foot = surface.lift(surface.project(p));
                best = Some((d, s, foot));
            }
        }
        best
    };
    let component_of_bar = |i: usize| component[i];

    // ---- 1. groups apart from the main structure
    if let Some(main) = main {
        let mut groups: BTreeMap<usize, (Vec<usize>, Vec<usize>)> = BTreeMap::new();
        for (x, &c) in component.iter().enumerate() {
            if c == main {
                continue;
            }
            let e = groups.entry(c).or_default();
            if x < b_count {
                e.0.push(x);
            } else {
                e.1.push(x - b_count);
            }
        }
        for (_, (gb, gs)) in groups {
            let bars_only = gs.is_empty();
            let mut f = Finding::new(
                "floating_group",
                if bars_only {
                    Class::Plaxis
                } else {
                    Class::Review
                },
            );
            f.value = (gb.len() + gs.len()) as f64;
            f.detail = format!("bars {}, surfaces {}", gb.len(), gs.len());
            f.points.push(match (gb.first(), gs.first()) {
                (Some(&b), _) => bars[b].p0.to_array(),
                (_, Some(&s)) => surfaces[s].lo.to_array(),
                _ => [0.; 3],
            });
            f.bars = gb.iter().copied().take(LISTED).collect();
            f.surfaces = gs.iter().copied().take(LISTED).collect();
            // The closest ways to the main structure: bar ends of the group
            // against bars of the main one.
            let mut options_: Vec<(f64, Fix)> = vec![];
            for &b in &gb {
                for end in 0..2 {
                    let p = if end == 0 { bars[b].p0 } else { bars[b].p1 };
                    let skip = |j: usize| component_of_bar(j) != main;
                    if let Some((d, j, _)) = nearest_bar(p, &skip, group_reach) {
                        options_.push((d, fix("connect_bars", json!({"op": "connect_bars", "a": b, "b": j, "tolerance": d + 1e-3}), d)));
                    }
                    if let Some((d, s, foot)) = nearest_surface(p, group_reach) {
                        if component[b_count + s] == main {
                            options_.push((d, fix("move_end_onto_surface", json!({"op": "move_vertex", "vertex": bars[b].ends[end], "to": foot.to_array()}), d)));
                        }
                    }
                }
            }
            options_.sort_by(|x, y| x.0.total_cmp(&y.0));
            f.fixes = options_.into_iter().map(|x| x.1).take(3).collect();
            if bars_only {
                // Deleting all of them (highest index first, as edits
                // renumber the bars).
                let mut bars_desc = gb.clone();
                bars_desc.sort_unstable_by(|a, b| b.cmp(a));
                f.fixes.push(Fix {
                    title: "delete_group".into(),
                    edit: json!({"op": "delete_bars", "bars": bars_desc}),
                    distance: 0.,
                });
            }
            findings.push(f);
        }
    }

    // ---- 2. bar ends connected to nothing
    let mut used_by_surface: BTreeSet<usize> = BTreeSet::new();
    for s in surfaces {
        used_by_surface.extend(s.vertex_ids.iter().copied());
    }
    let point_contact: BTreeSet<usize> = contacts
        .iter()
        .filter_map(|c| match *c {
            Contact::Point { vertex, .. } => Some(vertex),
            _ => None,
        })
        .collect();
    let mut seen = BTreeSet::new();
    for (i, bar) in bars.iter().enumerate() {
        for end in 0..2 {
            let v = bar.ends[end];
            let alone = at_vertex.get(&v).is_some_and(|l| l.len() == 1)
                && !used_by_surface.contains(&v)
                && !point_contact.contains(&v);
            if !alone || !seen.insert((i, v)) {
                continue;
            }
            let p = if end == 0 { bar.p0 } else { bar.p1 };
            // Seen from this bar's own point of view: not its own span.
            let mut found: Vec<(f64, Fix)> = vec![];
            if let Some((d, j, q)) = nearest_bar(p, &|j| j == i, reach) {
                found.push((
                    d,
                    fix(
                        "connect_bars",
                        json!({"op": "connect_bars", "a": i, "b": j, "tolerance": d + 1e-3}),
                        d,
                    ),
                ));
                // A node of that bar next to the closest point.
                if let Some(&w) = bars[j]
                    .nodes
                    .iter()
                    .filter(|&&w| {
                        p3(vertices[w]).distance(p) <= reach && !bars[j].nodes.contains(&v)
                    })
                    .min_by(|&&a, &&b| {
                        p3(vertices[a])
                            .distance(q)
                            .total_cmp(&p3(vertices[b]).distance(q))
                    })
                {
                    let dw = p3(vertices[w]).distance(p);
                    found.push((
                        dw,
                        fix(
                            "merge_end_into_vertex",
                            json!({"op": "merge_vertices", "drop": v, "keep": w}),
                            dw,
                        ),
                    ));
                }
            }
            if let Some((d, s, foot)) = nearest_surface(p, reach) {
                let _ = s;
                found.push((
                    d,
                    fix(
                        "move_end_onto_surface",
                        json!({"op": "move_vertex", "vertex": v, "to": foot.to_array()}),
                        d,
                    ),
                ));
            }
            found.sort_by(|x, y| x.0.total_cmp(&y.0));
            let mut f = Finding::new("free_bar_end", Class::Review)
                .bars([i])
                .vertex(v)
                .at(p);
            f.value = found.first().map_or(0., |x| x.0);
            f.fixes = found.into_iter().map(|x| x.1).take(3).collect();
            findings.push(f);
        }
    }

    // ---- 3. links of the source model lost in the result
    let Some(links) = context.links else {
        return;
    };
    let by_source: BTreeMap<usize, usize> = axes
        .iter()
        .enumerate()
        .map(|(i, a)| (a.source_axis, i))
        .collect();
    let mut reported = BTreeSet::new();
    for link in &links.bar_bar {
        let (Some(&a), Some(&b)) = (by_source.get(&link.a), by_source.get(&link.b)) else {
            continue;
        };
        if a == b || !reported.insert((a.min(b), a.max(b))) {
            continue;
        }
        if bars[a].nodes.intersection(&bars[b].nodes).next().is_some() {
            continue;
        }
        let (d, x) = closest_points(bars[a].p0, bars[a].p1, bars[b].p0, bars[b].p1);
        if d <= eps {
            // Touching without a shared node: reported as an unshared bar
            // intersection already.
            continue;
        }
        let apart = component[a] != component[b];
        let mut f = Finding::new(
            "lost_bar_link",
            if apart { Class::Failure } else { Class::Plaxis },
        )
        .bars([a, b])
        .at(x);
        f.value = d;
        f.detail = format!("source node {}", link.node);
        let mut found = vec![fix(
            "connect_bars",
            json!({"op": "connect_bars", "a": a, "b": b, "tolerance": d + 1e-3}),
            d,
        )];
        // The nearest pair of ends.
        let pairs = [(0usize, 0usize), (0, 1), (1, 0), (1, 1)];
        if let Some(&(ea, eb)) = pairs.iter().min_by(|x, y| {
            let dist = |&(i, j): &(usize, usize)| {
                let p = |bar: &Bar, e: usize| if e == 0 { bar.p0 } else { bar.p1 };
                p(&bars[a], i).distance(p(&bars[b], j))
            };
            dist(x).total_cmp(&dist(y))
        }) {
            let pa = if ea == 0 { bars[a].p0 } else { bars[a].p1 };
            let pb = if eb == 0 { bars[b].p0 } else { bars[b].p1 };
            let de = pa.distance(pb);
            if de <= d * 3. + reach {
                found.push(fix("merge_end_into_vertex", json!({"op": "merge_vertices", "drop": bars[a].ends[ea], "keep": bars[b].ends[eb]}), de));
            }
        }
        f.fixes = found;
        findings.push(f);
    }
    // Bars attached to a shell in the source and to none of its surfaces.
    let mut surfaces_of_patch: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (s, &p) in context.patches.iter().enumerate() {
        surfaces_of_patch.entry(p).or_default().push(s);
    }
    let mut contact_of: BTreeMap<usize, BTreeSet<usize>> = BTreeMap::new();
    for c in contacts {
        let (Contact::Point { axis, surface, .. } | Contact::Interval { axis, surface, .. }) = *c;
        contact_of.entry(axis).or_default().insert(surface);
    }
    let mut reported = BTreeSet::new();
    for link in &links.bar_surface {
        let Some(&i) = by_source.get(&link.axis) else {
            continue;
        };
        let Some(list) = surfaces_of_patch.get(&link.patch) else {
            continue;
        };
        if !reported.insert((i, link.patch)) {
            continue;
        }
        let touching = |s: usize| {
            contact_of.get(&i).is_some_and(|c| c.contains(&s))
                || bars[i]
                    .nodes
                    .iter()
                    .any(|v| surfaces[s].vertex_ids.contains(v))
        };
        if list.iter().any(|&s| touching(s)) {
            continue;
        }
        // The node of this bar the link was made at.
        let node = axes[i]
            .anchors
            .iter()
            .find(|n| n.source_node == link.node)
            .map(|n| n.vertex);
        let Some(node) = node else {
            continue;
        };
        let p = p3(vertices[node]);
        let nearest = list
            .iter()
            .map(|&s| (surfaces[s].distance(p), s))
            .min_by(|x, y| x.0.total_cmp(&y.0));
        let Some((d, s)) = nearest else {
            continue;
        };
        let apart = list.iter().all(|&s| component[b_count + s] != component[i]);
        let mut f = Finding::new(
            "lost_surface_link",
            if apart { Class::Failure } else { Class::Plaxis },
        )
        .bars([i])
        .surfaces([s])
        .vertex(node)
        .at(p);
        f.value = d;
        f.detail = format!("source node {}, patch {}", link.node, link.patch);
        let foot = surfaces[s].lift(surfaces[s].project(p));
        let is_end = axes[i].endpoints.contains(&node);
        if is_end {
            f.fixes.push(fix(
                "move_end_onto_surface",
                json!({"op": "move_vertex", "vertex": node, "to": foot.to_array()}),
                d,
            ));
        }
        f.fixes.push(fix(
            "connect_bar_to_surfaces",
            json!({"op": "connect_bar_to_surfaces", "bar": i}),
            d,
        ));
        findings.push(f);
    }
}

fn fix(title: &str, edit: serde_json::Value, distance: f64) -> Fix {
    Fix {
        title: title.into(),
        edit,
        distance,
    }
}
