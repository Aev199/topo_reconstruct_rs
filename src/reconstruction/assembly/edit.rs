//! Manual edits of an assembled geometry (the editor). Every edit is one
//! transactional operation of the model: it validates every surface it
//! touches and leaves the state unchanged on error. Bar axes follow moved
//! and merged vertices; bar-surface contacts are recomputed from the
//! geometry after every edit, as at the end of the assembly.
use super::bars::{self, Axis, Contact};
use super::cleanup::{self, Bars};
use super::gaps;
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

/// What the editor changes: the model, the bars and the per-surface data.
#[derive(Debug, Clone)]
pub struct State {
    pub model: Model,
    pub axes: Vec<Axis>,
    pub contacts: Vec<Contact>,
    pub stiffness: Vec<u32>,
    pub patches: Vec<usize>,
    pub removed: Vec<RemovedSurface>,
    /// Gaps marked as joints: (vertex, surface) kept open on purpose.
    pub joints: BTreeSet<(usize, usize)>,
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
            joints: BTreeSet::new(),
        }
    }

    /// Surfaces renumbered after a removal: per-surface data, contacts,
    /// joints and earlier `into` references follow.
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
    /// The bars it ends follow; an interior bar node does not move.
    pub fn move_vertex(&mut self, v: usize, to: [f64; 3]) -> Result<String, String> {
        self.check_vertex(v)?;
        let from = DVec3::from_array(self.model.vertices()[v]);
        let mut trial = self.model.clone();
        cleanup::move_with_axes(&mut trial, &self.axes, v, DVec3::from_array(to), f64::MAX)?;
        self.model = trial;
        self.refresh();
        Ok(format!("moved {:.4}", from.distance(DVec3::from_array(to))))
    }

    /// Merge vertex `drop` into `keep`: both move onto the planes of every
    /// surface of either, bars follow.
    pub fn merge_vertices(&mut self, drop: usize, keep: usize) -> Result<String, String> {
        self.check_vertex(drop)?;
        self.check_vertex(keep)?;
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

    /// Mark the gap between a vertex and a surface as a joint: kept open.
    pub fn mark_joint(&mut self, v: usize, s: usize) -> Result<String, String> {
        self.check_vertex(v)?;
        self.check_surface(s)?;
        self.joints.insert((v, s));
        Ok("marked as joint".into())
    }
}

#[cfg(test)]
mod tests {
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
            joints: BTreeSet::new(),
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
}
