//! Manual edits of an assembled geometry (the editor). Every edit is one
//! transactional operation of the model: it validates every surface it
//! touches and leaves the state unchanged on error. Bar axes follow moved
//! and merged vertices; bar-surface contacts are recomputed from the
//! geometry after every edit, as at the end of the assembly.
use super::bars::{self, Anchor, Axis, Contact};
use super::cleanup::{self, Bars};
use super::gaps;
use super::junctions;
use crate::reconstruction::Model;
use glam::DVec3;
use serde::Serialize;
use std::collections::BTreeSet;

/// A surface taken out of the geometry by an edit, with its provenance.
#[derive(Debug, Clone, Serialize, serde::Deserialize, PartialEq)]
pub struct RemovedSurface {
    /// `deleted` or `joined`.
    pub reason: String,
    pub stiffness: u32,
    pub patch: usize,
    pub source_elements: Vec<u32>,
    /// For `joined`: the surface that took its material (current index).
    pub into: Option<usize>,
}

/// A bar (or a piece of one) taken out by an edit, with its provenance.
#[derive(Debug, Clone, Serialize, serde::Deserialize, PartialEq)]
pub struct RemovedBar {
    /// `deleted` (the whole bar) or `collapsed` (a short piece merged away;
    /// `removed` when nothing of the bar was left).
    pub reason: String,
    pub source_axis: usize,
    pub removed: bool,
    /// Source bar elements no longer represented.
    pub source_elements: Vec<u32>,
}

/// What the editor changes: the model, the bars and the per-surface data.
#[derive(Debug, Clone)]
pub struct State {
    pub model: Model,
    pub axes: Vec<Axis>,
    pub contacts: Vec<Contact>,
    pub stiffness: Vec<u32>,
    pub patches: Vec<usize>,
    pub removed: Vec<RemovedSurface>,
    pub removed_bars: Vec<RemovedBar>,
    /// Gaps marked as joints: (vertex, surface) kept open on purpose.
    pub joints: BTreeSet<(usize, usize)>,
    /// Set when the upper storeys were cut off (`cutoff`).
    pub cut: Option<super::cutoff::Cut>,
}

impl State {
    pub fn from_report(report: &super::Report) -> Self {
        State {
            model: report.preview.clone(),
            axes: report.axis_assembly.axes.clone(),
            contacts: report.axis_assembly.contacts.clone(),
            stiffness: report.surface_stiffness.clone(),
            patches: report.surface_source_patches.clone(),
            removed: vec![],
            removed_bars: vec![],
            joints: BTreeSet::new(),
            cut: None,
        }
    }

    fn check_bar(&self, b: usize) -> Result<&Axis, String> {
        self.axes.get(b).ok_or_else(|| format!("no bar {b}"))
    }

    fn ends(&self, b: usize) -> [DVec3; 2] {
        self.axes[b]
            .endpoints
            .map(|v| DVec3::from_array(self.model.vertices()[v]))
    }

    /// Surfaces renumbered after a removal: per-surface data, contacts,
    /// joints and earlier `into` references follow.
    pub(super) fn renumber_after_cut(&mut self, index: &[Option<usize>]) {
        self.renumber(index);
    }

    pub(super) fn refresh_after_cut(&mut self) {
        self.refresh();
    }

    fn renumber(&mut self, index: &[Option<usize>]) {
        let kept = |s: &usize| index[*s].is_some();
        self.stiffness = (0..index.len())
            .filter(kept)
            .map(|s| self.stiffness[s])
            .collect();
        self.patches = (0..index.len())
            .filter(kept)
            .map(|s| self.patches[s])
            .collect();
        self.joints = self
            .joints
            .iter()
            .filter_map(|&(v, s)| index[s].map(|n| (v, n)))
            .collect();
        for r in &mut self.removed {
            r.into = r.into.and_then(|s| index.get(s).copied().flatten());
        }
    }

    fn refresh(&mut self) {
        self.model.refresh_orphaned_edges();
        // Contacts with removed surfaces go; the rest are recomputed.
        self.contacts.retain(|c| {
            let (Contact::Point { surface, .. } | Contact::Interval { surface, .. }) = c;
            *surface < self.model.surfaces().len()
        });
        bars::refresh_contacts(&self.model, &self.axes, &mut self.contacts);
    }

    fn check_vertex(&self, v: usize) -> Result<(), String> {
        (v < self.model.vertices().len())
            .then_some(())
            .ok_or_else(|| format!("no vertex {v}"))
    }

    fn check_surface(&self, s: usize) -> Result<(), String> {
        (s < self.model.surfaces().len())
            .then_some(())
            .ok_or_else(|| format!("no surface {s}"))
    }

    /// Move a vertex; it must stay on the plane of every surface using it.
    /// The bars it ends follow; an interior bar node only slides along its
    /// bars.
    pub fn move_vertex(&mut self, v: usize, to: [f64; 3]) -> Result<String, String> {
        self.check_vertex(v)?;
        let from = DVec3::from_array(self.model.vertices()[v]);
        let mut trial = self.model.clone();
        let interior = self
            .axes
            .iter()
            .any(|a| !a.endpoints.contains(&v) && a.anchors.iter().any(|n| n.vertex == v));
        if interior {
            // An interior bar node slides along its bars only.
            let mut axes = self.axes.clone();
            slide_node(&mut trial, &mut axes, v, DVec3::from_array(to))?;
            self.axes = axes;
        } else {
            cleanup::move_with_axes(&mut trial, &self.axes, v, DVec3::from_array(to), f64::MAX)?;
        }
        self.model = trial;
        self.refresh();
        Ok(format!("moved {:.4}", from.distance(DVec3::from_array(to))))
    }

    /// Merge vertex `drop` into `keep`: both move onto the planes of every
    /// surface of either, bars follow.
    pub fn merge_vertices(&mut self, drop: usize, keep: usize) -> Result<String, String> {
        self.check_vertex(drop)?;
        self.check_vertex(keep)?;
        let consecutive = self.axes.iter().any(|axis| {
            let mut anchors = axis.anchors.clone();
            anchors.sort_by(|x, y| x.t.total_cmp(&y.t));
            anchors.windows(2).any(|w| {
                let pair = (w[0].vertex, w[1].vertex);
                pair == (drop, keep) || pair == (keep, drop)
            })
        });
        if consecutive {
            return self.collapse_bar_piece(drop, keep);
        }
        let mut trial = self.model.clone();
        let (mut axes, mut contacts) = (self.axes.clone(), self.contacts.clone());
        let movement = cleanup::merge(
            &mut trial,
            &mut Bars {
                axes: &mut axes,
                contacts: &mut contacts,
            },
            drop,
            keep,
            f64::MAX,
            &BTreeSet::new(),
        )?;
        (self.model, self.axes, self.contacts) = (trial, axes, contacts);
        self.joints = self
            .joints
            .iter()
            .map(|&(v, s)| (if v == drop { keep } else { v }, s))
            .collect();
        self.refresh();
        Ok(format!("merged, moved {movement:.4}"))
    }

    /// Delete a surface; its source elements are kept as provenance.
    pub fn delete_surface(&mut self, s: usize) -> Result<String, String> {
        self.check_surface(s)?;
        let surface = &self.model.surfaces()[s];
        self.removed.push(RemovedSurface {
            reason: "deleted".into(),
            stiffness: self.stiffness[s],
            patch: self.patches[s],
            source_elements: surface.source_elements.clone(),
            into: None,
        });
        let index = self.model.remove_surfaces(&BTreeSet::from([s]));
        self.renumber(&index);
        self.refresh();
        Ok("deleted".into())
    }

    /// Join surface `other` into `keep`: one plane, one stiffness, at least
    /// one common contour edge.
    pub fn join_surfaces(&mut self, keep: usize, other: usize) -> Result<String, String> {
        self.check_surface(keep)?;
        self.check_surface(other)?;
        if self.stiffness[keep] != self.stiffness[other] {
            return Err("different stiffness: material regions are never joined".into());
        }
        let record = RemovedSurface {
            reason: "joined".into(),
            stiffness: self.stiffness[other],
            patch: self.patches[other],
            source_elements: self.model.surfaces()[other].source_elements.clone(),
            into: Some(keep),
        };
        let index = self
            .model
            .join_surfaces(keep, other)
            .map_err(|e| format!("join_{e:?}"))?;
        self.removed.push(record);
        self.renumber(&index);
        self.refresh();
        Ok("joined".into())
    }

    /// Split a model edge at the point of it nearest to `at`, for every
    /// surface using it.
    pub fn split_edge(&mut self, edge: usize, at: [f64; 3]) -> Result<String, String> {
        let [a, b] = *self
            .model
            .edges()
            .get(edge)
            .ok_or_else(|| format!("no edge {edge}"))?;
        let (pa, pb) = (
            DVec3::from_array(self.model.vertices()[a]),
            DVec3::from_array(self.model.vertices()[b]),
        );
        let d = pb - pa;
        let t = ((DVec3::from_array(at) - pa).dot(d) / d.length_squared()).clamp(0., 1.);
        let mut trial = self.model.clone();
        let v = trial
            .add_vertex((pa + d * t).to_array())
            .map_err(|e| format!("split_{e:?}"))?;
        trial
            .split_edge(edge, v)
            .map_err(|e| format!("split_{e:?}"))?;
        self.model = trial;
        self.refresh();
        Ok(format!("split at {t:.4}, vertex {v}"))
    }

    /// Close the gap between a vertex and a surface ("not a joint"),
    /// moving the vertex by at most `tolerance`.
    pub fn close_gap(&mut self, v: usize, s: usize, tolerance: f64) -> Result<String, String> {
        self.check_vertex(v)?;
        self.check_surface(s)?;
        let mut trial = self.model.clone();
        let (mut axes, mut contacts) = (self.axes.clone(), self.contacts.clone());
        let (kind, movement) = gaps::close_single(
            &mut trial,
            &mut Bars {
                axes: &mut axes,
                contacts: &mut contacts,
            },
            v,
            s,
            tolerance,
        )?;
        (self.model, self.axes, self.contacts) = (trial, axes, contacts);
        self.joints.remove(&(v, s));
        self.refresh();
        Ok(format!("{kind}, moved {movement:.4}"))
    }

    /// Merge two consecutive nodes of a bar: the piece between them
    /// collapses (a bar left without length disappears), other bars and the
    /// surfaces follow as in a merge.
    fn collapse_bar_piece(&mut self, drop: usize, keep: usize) -> Result<String, String> {
        let length = DVec3::from_array(self.model.vertices()[drop])
            .distance(DVec3::from_array(self.model.vertices()[keep]));
        let mut trial = self.model.clone();
        let (mut axes, mut contacts) = (self.axes.clone(), self.contacts.clone());
        let collapsed = cleanup::collapse_piece(
            &mut trial,
            &mut Bars {
                axes: &mut axes,
                contacts: &mut contacts,
            },
            drop,
            keep,
            length,
            f64::MAX,
            &BTreeSet::new(),
            &[],
        )?;
        (self.model, self.axes, self.contacts) = (trial, axes, contacts);
        self.joints = self
            .joints
            .iter()
            .map(|&(v, s)| (if v == drop { keep } else { v }, s))
            .collect();
        for c in &collapsed {
            self.removed_bars.push(RemovedBar {
                reason: "collapsed".into(),
                source_axis: c.source_axis,
                removed: c.removed,
                source_elements: c.elements.clone(),
            });
        }
        self.refresh();
        Ok(format!("bar piece collapsed, length {length:.4}"))
    }

    /// Delete several bars (a floating group) as one edit.
    pub fn delete_bars(&mut self, bars: &[usize]) -> Result<String, String> {
        let mut list = bars.to_vec();
        list.sort_unstable_by(|a, b| b.cmp(a));
        list.dedup();
        for &b in &list {
            self.check_bar(b)?;
        }
        // Highest index first: each deletion renumbers only the bars after it.
        for &b in &list {
            self.delete_bar(b)?;
        }
        Ok(format!("{} bars deleted", list.len()))
    }

    /// Delete a bar; its source elements are kept as provenance.
    pub fn delete_bar(&mut self, b: usize) -> Result<String, String> {
        let axis = self.check_bar(b)?.clone();
        let mut elements: Vec<u32> = axis.spans.iter().map(|s| s.element).collect();
        elements.sort_unstable();
        elements.dedup();
        self.removed_bars.push(RemovedBar {
            reason: "deleted".into(),
            source_axis: axis.source_axis,
            removed: true,
            source_elements: elements,
        });
        self.axes.remove(b);
        self.contacts.retain(|c| {
            let (Contact::Point { axis, .. } | Contact::Interval { axis, .. }) = c;
            *axis != b
        });
        for c in &mut self.contacts {
            let (Contact::Point { axis, .. } | Contact::Interval { axis, .. }) = c;
            if *axis > b {
                *axis -= 1;
            }
        }
        self.refresh();
        Ok("bar deleted".into())
    }

    /// Two bars crossing within `tolerance` share a node at the crossing: a
    /// node of both near it moves onto the crossing (along the bars), or
    /// the crossing becomes a new node of both (a node of one next to it
    /// slides there). Bars are never bent: bars passing each other more
    /// than the precision apart are refused.
    pub fn connect_bars(&mut self, a: usize, b: usize, tolerance: f64) -> Result<String, String> {
        self.check_bar(a)?;
        self.check_bar(b)?;
        if a == b {
            return Err("same_bar".into());
        }
        let ([p0, p1], [q0, q1]) = (self.ends(a), self.ends(b));
        let Some((s, t)) = closest_parameters(p0, p1, q0, q1) else {
            // Parallel bars never cross: an end of one is seated on the
            // side of the other.
            return self.seat_end(a, b, tolerance);
        };
        let (pa, pb) = (p0.lerp(p1, s), q0.lerp(q1, t));
        let distance = pa.distance(pb);
        if distance > tolerance {
            return Err(format!("bars_apart: {distance:e}"));
        }
        let x = (pa + pb) / 2.;
        let node =
            |axis: &Axis| -> BTreeSet<usize> { axis.anchors.iter().map(|n| n.vertex).collect() };
        let common: Vec<usize> = node(&self.axes[a])
            .intersection(&node(&self.axes[b]))
            .copied()
            .filter(|&v| DVec3::from_array(self.model.vertices()[v]).distance(x) <= tolerance)
            .collect();
        if let Some(&v) = common.first() {
            // Already a node of both: it slides onto the crossing along the
            // bar it ends (or onto the exact crossing when inside both).
            let p = DVec3::from_array(self.model.vertices()[v]);
            let target = if self.axes[b].endpoints.contains(&v) {
                pb
            } else if self.axes[a].endpoints.contains(&v) {
                pa
            } else {
                x
            };
            let mut trial = self.model.clone();
            let mut axes = self.axes.clone();
            slide_node(&mut trial, &mut axes, v, target)?;
            (self.model, self.axes) = (trial, axes);
            self.refresh();
            return Ok(format!("shared node {v}, moved {:.6}", p.distance(target)));
        }
        // The node of either bar nearest to the crossing joins the other
        // bar (a bar end within the tolerance of the other bar's span is
        // drawn onto it); without one, the crossing becomes a new node.
        let mut best: Option<(f64, usize, usize, Anchor)> = None;
        for (k, o) in [(a, b), (b, a)] {
            for n in &self.axes[k].anchors {
                let d = DVec3::from_array(self.model.vertices()[n.vertex]).distance(x);
                if d <= tolerance && best.as_ref().is_none_or(|x| d < x.0) {
                    best = Some((d, k, o, n.clone()));
                }
            }
        }
        let mut trial = self.model.clone();
        let mut axes = self.axes.clone();
        let (v, o, source_node, target) = match best {
            Some((_, k, o, n)) => {
                let [o0, o1] = self.ends(o);
                let p = DVec3::from_array(self.model.vertices()[n.vertex]);
                // An end moves onto the other bar; an interior node only
                // along its own bar (onto the exact crossing).
                if self.axes[k].endpoints.contains(&n.vertex) {
                    let d = o1 - o0;
                    let target = o0 + d * ((p - o0).dot(d) / d.length_squared()).clamp(0., 1.);
                    cleanup::move_with_axes(&mut trial, &axes, n.vertex, target, f64::MAX)?;
                    (n.vertex, o, n.source_node, target)
                } else {
                    // An interior node slides along its bar (and any other
                    // bar it is a node of) onto the crossing.
                    if distance > self.model.precision {
                        return Err(format!("skew_bars: {distance:e}"));
                    }
                    slide_node(&mut trial, &mut axes, n.vertex, x)?;
                    (n.vertex, o, n.source_node, x)
                }
            }
            None => {
                if distance > self.model.precision {
                    return Err(format!("skew_bars: {distance:e}"));
                }
                let v = trial
                    .add_vertex(x.to_array())
                    .map_err(|e| format!("{e:?}"))?;
                let t = DVec3::from_array(trial.vertices()[v]);
                let [p0, p1] = self.ends(a);
                let s = (t - p0).dot(p1 - p0) / (p1 - p0).length_squared();
                axes[a].anchors.push(Anchor {
                    source_node: bars::NO_SOURCE_NODE,
                    vertex: v,
                    t: s,
                });
                axes[a].anchors.sort_by(|x, y| x.t.total_cmp(&y.t));
                (v, b, bars::NO_SOURCE_NODE, x)
            }
        };
        let [o0, o1] = self.ends(o);
        let t = (target - o0).dot(o1 - o0) / (o1 - o0).length_squared();
        if !(t > 0. && t < 1.) {
            return Err("crossing_outside_bar".into());
        }
        let minimum = trial.minimum_edge;
        if axes[o]
            .anchors
            .iter()
            .any(|n| DVec3::from_array(trial.vertices()[n.vertex]).distance(target) < minimum)
        {
            return Err("ShortEdge".into());
        }
        axes[o].anchors.push(Anchor {
            source_node,
            vertex: v,
            t,
        });
        axes[o].anchors.sort_by(|x, y| x.t.total_cmp(&y.t));
        self.model = trial;
        self.axes = axes;
        self.refresh();
        Ok(format!("bars share vertex {v}"))
    }

    /// An end of one of two (parallel) bars within `tolerance` of the other
    /// one's span moves onto it and becomes a node of both.
    fn seat_end(&mut self, a: usize, b: usize, tolerance: f64) -> Result<String, String> {
        let minimum = self.model.minimum_edge();
        let mut best: Option<(f64, usize, usize, usize, DVec3)> = None;
        for (k, o) in [(a, b), (b, a)] {
            let [o0, o1] = self.ends(o);
            let d = o1 - o0;
            for end in 0..2 {
                let v = self.axes[k].endpoints[end];
                let p = DVec3::from_array(self.model.vertices()[v]);
                let t = (p - o0).dot(d) / d.length_squared();
                let foot = o0 + d * t;
                let distance = p.distance(foot);
                let clear = self.axes[o].anchors.iter().all(|n| {
                    DVec3::from_array(self.model.vertices()[n.vertex]).distance(foot) >= minimum
                });
                if t > 0.
                    && t < 1.
                    && clear
                    && distance <= tolerance
                    && best.as_ref().is_none_or(|x| distance < x.0)
                {
                    best = Some((distance, k, o, v, foot));
                }
            }
        }
        let Some((distance, k, o, v, foot)) = best else {
            return Err("parallel_bars".into());
        };
        let mut trial = self.model.clone();
        let mut axes = self.axes.clone();
        cleanup::move_with_axes(&mut trial, &axes, v, foot, f64::MAX)?;
        let [o0, o1] = axes[o]
            .endpoints
            .map(|e| DVec3::from_array(trial.vertices()[e]));
        let t = (foot - o0).dot(o1 - o0) / (o1 - o0).length_squared();
        let source_node = axes[k]
            .anchors
            .iter()
            .find(|n| n.vertex == v)
            .map_or(bars::NO_SOURCE_NODE, |n| n.source_node);
        axes[o].anchors.push(Anchor {
            source_node,
            vertex: v,
            t,
        });
        axes[o].anchors.sort_by(|x, y| x.t.total_cmp(&y.t));
        self.model = trial;
        self.axes = axes;
        self.refresh();
        Ok(format!("bars share vertex {v}, moved {distance:.4}"))
    }

    /// A bar passing through surfaces shares a node with each of them where
    /// it crosses (as the automatic assembly does).
    pub fn connect_bar_to_surfaces(&mut self, b: usize) -> Result<String, String> {
        self.check_bar(b)?;
        let mut trial = self.model.clone();
        // All bars: nodes shared with other bars stay guarded.
        let mut one = self.axes.clone();
        let precision = trial.precision;
        let imprinted: Vec<_> = bars::imprint_crossings_of(&mut trial, &mut one, Some(b))
            .into_iter()
            .filter(|i| i.kind != "slid_anchor" || i.movement > precision)
            .collect();
        // A bar lying in a surface's plane gets nodes where it crosses that
        // surface's contour (its piece inside then is an interval contact).
        let in_plane = imprint_in_plane(&mut trial, &mut one, b)?;
        if imprinted.is_empty() && in_plane == 0 {
            return Err("no_crossing_to_share".into());
        }
        self.axes = one;
        self.model = trial;
        self.refresh();
        Ok(format!(
            "bar shares {} node(s) with surfaces",
            imprinted.len() + in_plane
        ))
    }

    /// Two surfaces crossing, meeting in a T or along a common boundary
    /// without a shared edge get one: the junction becomes shared edges
    /// (embedded inside material). Bar nodes and contacts do not move.
    pub fn connect_surfaces(&mut self, a: usize, b: usize) -> Result<String, String> {
        self.check_surface(a)?;
        self.check_surface(b)?;
        if a == b {
            return Err("same_surface".into());
        }
        let mut interior = vec![BTreeSet::new(); self.model.surfaces().len()];
        for c in &self.contacts {
            if let Contact::Point {
                surface, vertex, ..
            } = *c
            {
                interior[surface].insert(vertex);
            }
        }
        let mut locked: BTreeSet<usize> = interior.iter().flatten().copied().collect();
        for axis in &self.axes {
            locked.extend(axis.endpoints);
            locked.extend(axis.anchors.iter().map(|n| n.vertex));
        }
        let mut trial = self.model.clone();
        let report = junctions::insert_pair(
            &mut trial,
            &junctions::Context {
                interior: &interior,
                locked: &locked,
                wall_end_tolerance: 0.,
            },
            a,
            b,
        );
        if let Some(issue) = report.issues.first() {
            return Err(format!("junction_{}", issue.reason));
        }
        if report.junctions.is_empty() {
            return Err(if report.already_conforming > 0 {
                "already_connected".into()
            } else {
                "no_junction".into()
            });
        }
        self.model = trial;
        self.refresh();
        let length: f64 = report.junctions.iter().map(|j| j.length).sum();
        Ok(format!("connected, junction {length:.4}"))
    }

    /// Mark the gap between a vertex and a surface as a joint: kept open.
    pub fn mark_joint(&mut self, v: usize, s: usize) -> Result<String, String> {
        self.check_vertex(v)?;
        self.check_surface(s)?;
        self.joints.insert((v, s));
        Ok("marked as joint".into())
    }
}

/// Nodes of bar `b` where it crosses the contour of a surface whose plane
/// it lies in: an existing contour vertex on the bar becomes a node of it,
/// otherwise a generated vertex splits the contour edge. Crossings within
/// the minimum edge of a node of the bar are left alone. Returns how many
/// nodes were added.
fn imprint_in_plane(model: &mut Model, axes: &mut [Axis], b: usize) -> Result<usize, String> {
    let precision = model.precision;
    let on_plane = 5. * precision;
    let minimum = model.minimum_edge;
    let mut added = 0;
    for s in 0..model.surfaces.len() {
        let [p, q] = axes[b]
            .endpoints
            .map(|v| DVec3::from_array(model.vertices[v]));
        let plane = model.planes[model.surfaces[s].plane].clone();
        if plane.distance(p.to_array()).abs() > on_plane
            || plane.distance(q.to_array()).abs() > on_plane
        {
            continue;
        }
        let d = q - p;
        let length = d.length();
        // Contour edges (with their vertices) crossed by the bar.
        let mut hits: Vec<(f64, usize, Option<usize>)> = vec![];
        for e in model.surface_edges(s).collect::<Vec<_>>() {
            let [u, w] = model.edges[e];
            let (a, c) = (
                DVec3::from_array(model.vertices[u]),
                DVec3::from_array(model.vertices[w]),
            );
            for (v, x) in [(u, a), (w, c)] {
                let t = (x - p).dot(d) / (length * length);
                if t > 0. && t < 1. && x.distance(p + d * t) <= precision {
                    hits.push((t, e, Some(v)));
                }
            }
            let Some((t, r)) = segment_parameters(p, q, a, c) else {
                continue;
            };
            if t > 0.
                && t < 1.
                && r > 0.
                && r < 1.
                && (p + d * t).distance(a + (c - a) * r) <= precision
            {
                let x = a + (c - a) * r;
                if x.distance(a) >= minimum && x.distance(c) >= minimum {
                    hits.push((t, e, None));
                }
            }
        }
        hits.sort_by(|x, y| x.0.total_cmp(&y.0));
        for (t, e, vertex) in hits {
            let near = axes[b]
                .anchors
                .iter()
                .any(|n| (n.t - t).abs() * length < minimum);
            if near {
                continue;
            }
            let v = match vertex {
                Some(v) => v,
                None => {
                    let v = model
                        .add_vertex((p + d * t).to_array())
                        .map_err(|e| format!("{e:?}"))?;
                    model
                        .split_edge_within(e, v, precision)
                        .map_err(|e| format!("split_{e:?}"))?;
                    v
                }
            };
            axes[b].anchors.push(Anchor {
                source_node: bars::NO_SOURCE_NODE,
                vertex: v,
                t,
            });
            axes[b].anchors.sort_by(|x, y| x.t.total_cmp(&y.t));
            added += 1;
        }
    }
    Ok(added)
}

/// Parameters (t on p-q, r on a-c) of the closest points of two segment
/// lines; `None` when parallel.
fn segment_parameters(p: DVec3, q: DVec3, a: DVec3, c: DVec3) -> Option<(f64, f64)> {
    let (u, v, w) = (q - p, c - a, p - a);
    let (aa, bb, cc, dd, ee) = (u.dot(u), u.dot(v), v.dot(v), u.dot(w), v.dot(w));
    let den = aa * cc - bb * bb;
    if den <= 1e-12 * aa * cc {
        return None;
    }
    Some(((bb * ee - cc * dd) / den, (aa * ee - bb * dd) / den))
}

/// Move node `v` of bars to `q` without bending any bar: `q` must be on the
/// line of every bar having `v` (inside the span for an interior node). A
/// bar end slides along its own line; the parameters of its nodes and
/// source spans follow the new end.
fn slide_node(model: &mut Model, axes: &mut [Axis], v: usize, q: DVec3) -> Result<(), String> {
    let precision = model.precision;
    let mut remaps = vec![];
    let mut interior = vec![];
    for (k, axis) in axes.iter().enumerate() {
        if !axis.anchors.iter().any(|a| a.vertex == v) && !axis.endpoints.contains(&v) {
            continue;
        }
        let [a, b] = axis.endpoints.map(|e| DVec3::from_array(model.vertices[e]));
        let d = b - a;
        let t = (q - a).dot(d) / d.length_squared();
        if q.distance(a + d * t) > precision {
            return Err("node_off_bar".into());
        }
        match axis.endpoints.iter().position(|&e| e == v) {
            None if !(t > 0. && t < 1.) => return Err("node_off_bar".into()),
            // An interior node takes its new parameter; source spans
            // bounded by it follow.
            None => interior.push((k, t)),
            // The new end at `t` of the old line: old parameters map
            // affinely onto the new span.
            Some(0) => remaps.push((k, t, 1.)),
            Some(_) => remaps.push((k, 0., t)),
        }
    }
    model
        .move_vertex(v, q.to_array())
        .map_err(|e| format!("move_{e:?}"))?;
    for (k, t) in interior {
        let axis = &mut axes[k];
        let n = axis.anchors.iter().position(|a| a.vertex == v).unwrap();
        let old = axis.anchors[n].t;
        for s in &mut axis.spans {
            for u in [&mut s.start_t, &mut s.end_t] {
                if (*u - old).abs() <= 1e-12 {
                    *u = t;
                }
            }
        }
        axis.anchors[n].t = t;
        if axis.anchors.windows(2).any(|w| w[1].t <= w[0].t) {
            return Err("node_beyond_neighbour".into());
        }
    }
    for (k, t0, t1) in remaps {
        if !(t1 - t0 > 0.) {
            return Err("bar_would_collapse".into());
        }
        // The moved end stays the end (0 or 1); anything else maps without
        // clamping, a node beyond the new end being refused below.
        let moved = if t0 != 0. { 0. } else { 1. };
        let map = |x: f64| {
            if x == moved {
                moved
            } else {
                (x - t0) / (t1 - t0)
            }
        };
        let axis = &mut axes[k];
        for a in &mut axis.anchors {
            a.t = if axis.endpoints[0] == a.vertex {
                0.
            } else if axis.endpoints[1] == a.vertex {
                1.
            } else {
                map(a.t)
            };
        }
        for s in &mut axis.spans {
            (s.start_t, s.end_t) = (map(s.start_t), map(s.end_t));
        }
        if axis.anchors.iter().any(|a| !(0. ..=1.).contains(&a.t))
            || axis.anchors.windows(2).any(|w| w[1].t <= w[0].t)
            || axis
                .spans
                .iter()
                .any(|s| !(s.start_t >= 0. && s.start_t < s.end_t && s.end_t <= 1.))
        {
            return Err("node_beyond_bar_end".into());
        }
    }
    Ok(())
}

/// Parameters of the closest points of two segments' lines, clamped to the
/// segments; `None` for parallel lines.
fn closest_parameters(p0: DVec3, p1: DVec3, q0: DVec3, q1: DVec3) -> Option<(f64, f64)> {
    let (u, v, w) = (p1 - p0, q1 - q0, p0 - q0);
    let (a, b, c, d, e) = (u.dot(u), u.dot(v), v.dot(v), u.dot(w), v.dot(w));
    let den = a * c - b * b;
    if den <= 1e-12 * a * c {
        return None;
    }
    let t = (a * e - b * d) / den;
    // Clamped to one segment, the point on the other is re-projected.
    let t = t.clamp(0., 1.);
    let s = ((q0 + v * t - p0).dot(u) / a).clamp(0., 1.);
    let t = ((p0 + u * s - q0).dot(v) / c).clamp(0., 1.);
    Some((s, t))
}

#[cfg(test)]
mod tests {
    use super::super::junctions::tests::slab as panel;
    use super::super::junctions::tests::{build, Placement};
    use super::*;

    fn state(model: Model, surfaces: usize) -> State {
        State {
            model,
            axes: vec![],
            contacts: vec![],
            stiffness: vec![1; surfaces],
            patches: (0..surfaces).collect(),
            removed: vec![],
            removed_bars: vec![],
            joints: BTreeSet::new(),
            cut: None,
        }
    }

    fn slab(x0: f64, x1: f64) -> (Vec<Vec<[f64; 3]>>, [f64; 3]) {
        (
            vec![vec![[x0, 0., 0.], [x1, 0., 0.], [x1, 2., 0.], [x0, 2., 0.]]],
            [0., 0., 1.],
        )
    }

    #[test]
    fn neighbouring_slabs_join_into_one_and_keep_provenance() {
        for place in Placement::all() {
            let m = build(&place, &[slab(0., 2.), slab(2., 5.)]);
            let mut s = state(m, 2);
            s.join_surfaces(0, 1).unwrap();
            assert_eq!(s.model.surfaces().len(), 1);
            // One rectangle: the two vertices on the common edge stay as
            // collinear contour vertices.
            assert_eq!(s.model.surfaces()[0].boundaries.len(), 1);
            assert_eq!(s.model.surfaces()[0].boundaries[0].len(), 6);
            assert_eq!(s.removed[0].reason, "joined");
            assert_eq!(s.removed[0].into, Some(0));
            assert_eq!(s.stiffness, vec![1]);
        }
    }

    #[test]
    fn joins_refuse_other_stiffness_and_separate_surfaces() {
        for place in Placement::all() {
            let m = build(&place, &[slab(0., 2.), slab(2., 5.), slab(6., 7.)]);
            let mut s = state(m, 3);
            s.stiffness[1] = 2;
            let before = s.model.clone();
            assert!(s.join_surfaces(0, 1).is_err());
            assert!(s.join_surfaces(0, 2).is_err());
            assert_eq!(format!("{:?}", s.model), format!("{before:?}"));
        }
    }

    #[test]
    fn delete_split_move_and_merge() {
        for place in Placement::all() {
            let m = build(&place, &[slab(0., 2.), slab(2., 5.), slab(6., 7.)]);
            let mut s = state(m, 3);
            s.delete_surface(2).unwrap();
            assert_eq!(s.model.surfaces().len(), 2);
            assert_eq!(s.removed[0].reason, "deleted");
            // Split the common edge in the middle: both slabs get the vertex.
            let shared = s
                .model
                .surface_edges(0)
                .find(|e| s.model.surface_edges(1).any(|f| f == *e))
                .unwrap();
            s.split_edge(shared, place.point([2., 1., 0.])).unwrap();
            let v = s.model.vertices().len() - 1;
            for k in 0..2 {
                assert!(s
                    .model
                    .surface_edges(k)
                    .any(|e| s.model.edges()[e].contains(&v)));
            }
            // Move it along the edge; off the plane it is refused.
            s.move_vertex(v, place.point([2., 1.2, 0.])).unwrap();
            assert!(s.move_vertex(v, place.point([2., 1.2, 0.1])).is_err());
            // Merge it into a corner of the edge: the edge pieces collapse.
            let corner = s.model.edges()[shared][0];
            assert!(s.merge_vertices(v, corner).is_ok());
        }
    }

    #[test]
    fn a_gap_marked_as_joint_or_closed() {
        for place in Placement::all() {
            let m = build(&place, &[slab(0., 2.), slab(2.02, 5.)]);
            let mut s = state(m, 2);
            let v = (0..s.model.vertices().len())
                .find(|&v| {
                    DVec3::from_array(s.model.vertices()[v])
                        .distance(DVec3::from_array(place.point([2.02, 0., 0.])))
                        < 1e-9 * place.scale.max(1.)
                })
                .unwrap();
            s.mark_joint(v, 0).unwrap();
            assert!(s.joints.contains(&(v, 0)));
            s.close_gap(v, 0, 0.05 * place.scale).unwrap();
            assert!(s.joints.is_empty());
            let p = DVec3::from_array(s.model.vertices()[v]);
            assert!(
                p.distance(DVec3::from_array(place.point([2., 0., 0.])))
                    < 1e-6 * place.scale.max(1.)
            );
        }
    }

    fn bar(m: &mut Model, place: &Placement, source_axis: usize, a: [f64; 3], b: [f64; 3]) -> Axis {
        let va = m.add_vertex(place.point(a)).unwrap();
        let vb = m.add_vertex(place.point(b)).unwrap();
        Axis {
            source_axis,
            endpoints: [va, vb],
            anchors: vec![
                Anchor {
                    source_node: 2 * source_axis as u32,
                    vertex: va,
                    t: 0.,
                },
                Anchor {
                    source_node: 2 * source_axis as u32 + 1,
                    vertex: vb,
                    t: 1.,
                },
            ],
            spans: vec![crate::reconstruction::recognize::SourceSpan {
                element: source_axis as u32,
                stiffness: 1,
                start_t: 0.,
                end_t: 1.,
            }],
        }
    }

    fn near(s: &State, v: usize, place: &Placement, p: [f64; 3]) -> bool {
        DVec3::from_array(s.model.vertices()[v]).distance(DVec3::from_array(place.point(p)))
            < 1e-6 * place.scale.max(1.)
    }

    /// Every bar node at its parameter on the straight bar.
    fn straight(s: &State) -> bool {
        s.axes.iter().all(|a| {
            let [p, q] = a
                .endpoints
                .map(|v| DVec3::from_array(s.model.vertices()[v]));
            a.anchors.iter().all(|n| {
                p.lerp(q, n.t)
                    .distance(DVec3::from_array(s.model.vertices()[n.vertex]))
                    <= 10. * s.model.precision
            })
        })
    }

    #[test]
    fn crossing_bars_share_a_node_and_pieces_collapse() {
        for place in Placement::all() {
            let mut m = build(&place, &[slab(0., 4.)]);
            let a = bar(&mut m, &place, 0, [0., 1., 1.], [4., 1., 1.]);
            let b = bar(&mut m, &place, 1, [2., 0., 1.], [2., 3., 1.]);
            // A bar ending 20 mm short of bar 0.
            let c = bar(&mut m, &place, 2, [1., -1., 1.], [1., 0.98, 1.]);
            let mut s = state(m, 1);
            s.axes = vec![a, b, c];
            let tol = 0.05 * place.scale;
            // A crossing inside both spans: a new node of both.
            s.connect_bars(0, 1, tol).unwrap();
            let v = s.model.vertices().len() - 1;
            assert!(near(&s, v, &place, [2., 1., 1.]));
            assert!(s.axes[0].anchors.iter().any(|n| n.vertex == v));
            assert!(s.axes[1].anchors.iter().any(|n| n.vertex == v));
            // Connected again: nothing to move.
            assert!(s.connect_bars(0, 1, tol).is_ok());
            // The end of bar 2 is drawn onto bar 0 and becomes its node.
            s.connect_bars(0, 2, tol).unwrap();
            let end = s.axes[2].endpoints[1];
            assert!(near(&s, end, &place, [1., 1., 1.]));
            let mut nodes: Vec<f64> = s.axes[0].anchors.iter().map(|n| n.t).collect();
            assert_eq!(nodes.len(), 4);
            nodes.dedup();
            assert!(nodes.windows(2).all(|w| w[0] < w[1]));
            assert!(straight(&s));
            // Bars far apart are refused.
            assert!(s.connect_bars(1, 2, tol).is_err());
            // A node of two bars moves along neither alone.
            assert!(s.move_vertex(end, place.point([1.2, 1., 1.])).is_err());
            // An interior node moves along its bar, not off it.
            assert!(s.move_vertex(v, place.point([2., 1., 1.5])).is_err());
            // Two consecutive nodes of bar 0 merge: the piece collapses and
            // bar 2 follows its end.
            s.merge_vertices(end, v).unwrap();
            assert_eq!(s.axes[0].anchors.len(), 3);
            assert_eq!(s.axes[2].endpoints[1], v);
            assert_eq!(s.removed_bars[0].reason, "collapsed");
            assert_eq!(s.removed_bars[0].source_axis, 0);
            assert!(straight(&s));
            // A deleted bar keeps its provenance.
            s.delete_bar(2).unwrap();
            assert_eq!(s.axes.len(), 2);
            assert_eq!(s.removed_bars[1].source_elements, vec![2]);
            // A node left on bar 0 alone slides along it, not off it.
            let d = bar(&mut s.model, &place, 3, [3., 0., 1.], [3., 2., 1.]);
            s.axes.push(d);
            s.connect_bars(0, 2, tol).unwrap();
            let n = s.model.vertices().len() - 1;
            s.delete_bar(2).unwrap();
            assert!(s.move_vertex(n, place.point([3., 1.1, 1.])).is_err());
            s.move_vertex(n, place.point([3.2, 1., 1.])).unwrap();
            assert!(straight(&s));
        }
    }

    #[test]
    fn a_wall_through_a_slab_is_connected_on_request() {
        use super::super::junctions::tests::wall;
        for place in Placement::all() {
            let m = build(&place, &[panel(0., 4.), wall(1., 3., -1., 1.)]);
            let mut s = state(m, 2);
            let shared = |s: &State| {
                s.model
                    .surface_edges(0)
                    .filter(|e| s.model.surface_edges(1).any(|f| f == *e))
                    .count()
            };
            assert_eq!(shared(&s), 0);
            s.connect_surfaces(0, 1).unwrap();
            assert!(shared(&s) > 0);
            assert!(!s.model.surfaces()[0].embedded_edges.is_empty());
            assert!(s.connect_surfaces(0, 1).is_err());
        }
    }

    #[test]
    fn a_bar_through_a_slab_shares_a_node_on_request() {
        for place in Placement::all() {
            let mut m = build(&place, &[slab(0., 2.)]);
            let a = bar(&mut m, &place, 0, [1., 1., -1.], [1., 1., 1.]);
            let mut s = state(m, 1);
            s.axes = vec![a];
            s.connect_bar_to_surfaces(0).unwrap();
            let v = s.axes[0].anchors[1].vertex;
            assert!(near(&s, v, &place, [1., 1., 0.]));
            assert!(s.connect_bar_to_surfaces(0).is_err());
        }
    }

    /// A node of bar `k` at `p`, with source spans split there.
    fn node(s: &mut State, place: &Placement, k: usize, p: [f64; 3]) -> usize {
        let v = s.model.add_vertex(place.point(p)).unwrap();
        let [a, b] = s.ends(k);
        let x = DVec3::from_array(s.model.vertices()[v]);
        let t = (x - a).dot(b - a) / (b - a).length_squared();
        let axis = &mut s.axes[k];
        axis.anchors.push(Anchor {
            source_node: 1000 + v as u32,
            vertex: v,
            t,
        });
        axis.anchors.sort_by(|x, y| x.t.total_cmp(&y.t));
        let mut spans = vec![];
        for span in &axis.spans {
            if span.start_t < t && t < span.end_t {
                spans.push(crate::reconstruction::recognize::SourceSpan {
                    end_t: t,
                    ..span.clone()
                });
                spans.push(crate::reconstruction::recognize::SourceSpan {
                    start_t: t,
                    element: span.element + 100,
                    stiffness: span.stiffness + 1,
                    ..span.clone()
                });
            } else {
                spans.push(span.clone());
            }
        }
        axis.spans = spans;
        v
    }

    #[test]
    fn nodes_never_coincide_or_pass_the_bar_end() {
        for place in Placement::all() {
            let mut m = build(&place, &[slab(0., 2.)]);
            let a = bar(&mut m, &place, 0, [0., 3., 1.], [4., 3., 1.]);
            let mut s = state(m, 1);
            s.axes = vec![a];
            let one = node(&mut s, &place, 0, [1., 3., 1.]);
            let two = node(&mut s, &place, 0, [2., 3., 1.]);
            assert!(bars::axis_defects(&s.model, &s.axes).is_empty());
            // Onto its neighbour: two nodes in one point, an empty span.
            assert!(s.move_vertex(two, place.point([1., 3., 1.])).is_err());
            // Along the bar between its neighbours: fine.
            s.move_vertex(two, place.point([1.5, 3., 1.])).unwrap();
            assert!(
                bars::axis_defects(&s.model, &s.axes).is_empty(),
                "{:?}",
                bars::axis_defects(&s.model, &s.axes)
            );
            // A bar end that is a node of another bar along the same line
            // slides; past an interior node of its own bar it is refused.
            let end = s.axes[0].endpoints[0];
            let mut long = bar(&mut s.model, &place, 1, [-1., 3., 1.], [6., 3., 1.]);
            long.anchors.insert(
                1,
                Anchor {
                    source_node: 9,
                    vertex: end,
                    t: 1. / 7.,
                },
            );
            s.axes.push(long);
            assert!(s.move_vertex(end, place.point([1.2, 3., 1.])).is_err());
            s.move_vertex(end, place.point([0.5, 3., 1.])).unwrap();
            assert!(
                bars::axis_defects(&s.model, &s.axes).is_empty(),
                "{:?}",
                bars::axis_defects(&s.model, &s.axes)
            );
            assert!(s.move_vertex(one, place.point([0., 3., 1.])).is_err());
            assert!(bars::axis_defects(&s.model, &s.axes).is_empty());
        }
    }

    #[test]
    fn an_interior_node_slides_onto_the_crossing() {
        for place in Placement::all() {
            let mut m = build(&place, &[slab(0., 2.)]);
            let a = bar(&mut m, &place, 0, [0., 3., 1.], [4., 3., 1.]);
            let b = bar(&mut m, &place, 1, [2., 0., 1.], [2., 6., 1.]);
            let mut s = state(m, 1);
            s.axes = vec![a, b];
            // A node of bar 0 20 mm from the crossing.
            let v = node(&mut s, &place, 0, [2.02, 3., 1.]);
            s.connect_bars(0, 1, 0.05 * place.scale).unwrap();
            assert!(near(&s, v, &place, [2., 3., 1.]));
            assert!(s.axes[1].anchors.iter().any(|n| n.vertex == v));
            assert!(bars::axis_defects(&s.model, &s.axes).is_empty());
        }
    }

    #[test]
    fn connecting_one_bar_never_bends_another() {
        for place in Placement::all() {
            let mut m = build(&place, &[slab(0., 2.)]);
            let up = bar(&mut m, &place, 0, [1., 1., -1.], [1., 1., 1.]);
            let across = bar(&mut m, &place, 1, [1., -1., 0.0005], [1., 3., 0.0005]);
            let mut s = state(m, 1);
            s.axes = vec![up, across];
            // One node of both bars 0.5 mm above the slab.
            let v = node(&mut s, &place, 0, [1., 1., 0.0005]);
            let t = 0.5;
            s.axes[1].anchors.push(Anchor {
                source_node: 7,
                vertex: v,
                t,
            });
            s.axes[1].anchors.sort_by(|x, y| x.t.total_cmp(&y.t));
            assert!(bars::axis_defects(&s.model, &s.axes).is_empty());
            let _ = s.connect_bar_to_surfaces(0);
            assert!(
                bars::axis_defects(&s.model, &s.axes).is_empty(),
                "{:?}",
                bars::axis_defects(&s.model, &s.axes)
            );
        }
    }

    #[test]
    fn a_bar_in_the_slab_plane_gets_nodes_on_its_contour() {
        for place in Placement::all() {
            let mut m = build(&place, &[panel(0., 4.)]);
            let a = bar(&mut m, &place, 0, [-1., 2., 0.], [5., 2., 0.]);
            let mut s = state(m, 1);
            s.axes = vec![a];
            s.refresh();
            assert!(!s
                .contacts
                .iter()
                .any(|c| matches!(c, Contact::Interval { .. })));
            s.connect_bar_to_surfaces(0).unwrap();
            assert_eq!(s.axes[0].anchors.len(), 4);
            assert!(bars::axis_defects(&s.model, &s.axes).is_empty());
            assert!(
                s.contacts
                    .iter()
                    .any(|c| matches!(c, Contact::Interval { .. })),
                "{:?}",
                s.contacts
            );
        }
    }

    #[test]
    fn the_end_of_a_bar_is_seated_on_a_parallel_bar_and_groups_are_deleted() {
        for place in Placement::all() {
            let mut m = build(&place, &[slab(0., 2.)]);
            let a = bar(&mut m, &place, 0, [0., 3., 0.], [0., 3., 4.]);
            // Bar 1 stands beside bar 0, its lower end 0.2 off bar 0's span.
            let b = bar(&mut m, &place, 1, [0.2, 3., 3.], [0.2, 3., 6.]);
            let c = bar(&mut m, &place, 2, [8., 8., 0.], [8., 8., 1.]);
            let d = bar(&mut m, &place, 3, [8., 9., 0.], [8., 9., 1.]);
            let mut s = state(m, 1);
            s.axes = vec![a, b, c, d];
            let tol = 0.25 * place.scale;
            assert!(s.connect_bars(0, 1, 0.1 * place.scale).is_err());
            s.connect_bars(0, 1, tol).unwrap();
            // One bar's end went onto the other bar: a node of both.
            let shared: Vec<usize> = s.axes[0]
                .anchors
                .iter()
                .map(|n| n.vertex)
                .filter(|v| s.axes[1].anchors.iter().any(|n| n.vertex == *v))
                .collect();
            assert_eq!(shared.len(), 1);
            assert!(bars::axis_defects(&s.model, &s.axes).is_empty());
            // A floating group goes at once, in one journal entry.
            s.delete_bars(&[3, 2]).unwrap();
            assert_eq!(s.axes.len(), 2);
            assert_eq!(s.removed_bars.len(), 2);
        }
    }
}
