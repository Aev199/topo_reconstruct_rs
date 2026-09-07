//! Conform boundary vertices and node common edges. This is geometric connectivity,
//! not a declaration of mechanical ties or interfaces.
use crate::config::ReconstructionConfig;
use crate::geometry::utils::get_plane_basis;
use crate::models::MacroPanel;
use crate::reconstructors::panels::order_valid_rings;
use geo::{Area, BooleanOps, LineString, Polygon};
use glam::DVec3;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

const EPS: f64 = 1e-8;

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct TopologySummary {
    pub merged_vertices: usize,
    pub snapped_junctions: usize,
    pub inserted_vertices: usize,
    pub connected_panel_pairs: usize,
    pub unresolved_near_junctions: usize,
    pub max_displacement: f64,
}

struct Graph<'a> {
    panels: &'a [MacroPanel],
    min_edge: f64,
    points: Vec<DVec3>,
    origins: Vec<Vec<DVec3>>,
    owners: Vec<BTreeSet<usize>>,
    parent: Vec<usize>,
    rings: Vec<Vec<Vec<usize>>>,
    baseline_overlap: BTreeMap<(usize, usize), f64>,
}

impl<'a> Graph<'a> {
    fn new(panels: &'a [MacroPanel], min_edge: f64) -> Self {
        let mut graph = Self {
            panels,
            min_edge,
            points: vec![],
            origins: vec![],
            owners: vec![],
            parent: vec![],
            rings: vec![],
            baseline_overlap: BTreeMap::new(),
        };
        for (owner, panel) in panels.iter().enumerate() {
            let mut rings = Vec::new();
            for ring in &panel.polygons {
                let mut ids = Vec::new();
                for p in ring {
                    let id = graph.points.len();
                    let p = DVec3::from_array(*p);
                    graph.points.push(p);
                    graph.origins.push(vec![p]);
                    graph.owners.push(BTreeSet::from([owner]));
                    graph.parent.push(id);
                    ids.push(id);
                }
                rings.push(ids);
            }
            graph.rings.push(rings);
        }
        for i in 0..panels.len() {
            for j in i + 1..panels.len() {
                if graph.coplanar(i, j) {
                    graph.baseline_overlap.insert((i, j), graph.overlap(i, j));
                }
            }
        }
        graph
    }
    fn root(&self, mut i: usize) -> usize {
        while self.parent[i] != i {
            i = self.parent[i];
        }
        i
    }
    fn coords(&self, i: usize) -> Vec<Vec<[f64; 3]>> {
        self.rings[i]
            .iter()
            .map(|r| {
                r.iter()
                    .map(|&id| self.points[self.root(id)].to_array())
                    .collect()
            })
            .collect()
    }
    fn coplanar(&self, i: usize, j: usize) -> bool {
        let a = &self.panels[i];
        let b = &self.panels[j];
        let na = DVec3::from_array(a.plane_normal);
        let nb = DVec3::from_array(b.plane_normal);
        let sign = if na.dot(nb) >= 0.0 { 1.0 } else { -1.0 };
        (na - nb * sign).length() < EPS && (a.plane_d - b.plane_d * sign).abs() < EPS
    }
    fn polygon(&self, i: usize, normal: DVec3) -> Polygon<f64> {
        let (u, v) = get_plane_basis(normal);
        let rings = self.coords(i);
        let convert = |r: &Vec<[f64; 3]>| -> LineString<f64> {
            let mut points: Vec<_> = r
                .iter()
                .map(|p| {
                    let p = DVec3::from_array(*p);
                    (p.dot(u), p.dot(v))
                })
                .collect();
            points.push(points[0]);
            LineString::from(points)
        };
        Polygon::new(convert(&rings[0]), rings[1..].iter().map(convert).collect())
    }
    fn overlap(&self, i: usize, j: usize) -> f64 {
        let normal = DVec3::from_array(self.panels[i].plane_normal);
        self.polygon(i, normal)
            .intersection(&self.polygon(j, normal))
            .unsigned_area()
    }
    fn valid(&self, affected: &BTreeSet<usize>) -> bool {
        for &i in affected {
            let normal = DVec3::from_array(self.panels[i].plane_normal);
            let coords = self.coords(i);
            if coords
                .iter()
                .flatten()
                .any(|p| (normal.dot(DVec3::from_array(*p)) + self.panels[i].plane_d).abs() > EPS)
            {
                return false;
            }
            for (ri, ring) in coords.iter().enumerate() {
                let original = &self.panels[i].polygons[ri];
                let signed_area = |points: &Vec<[f64; 3]>| {
                    let origin = DVec3::from_array(points[0]);
                    (0..points.len())
                        .map(|k| {
                            let a = DVec3::from_array(points[k]) - origin;
                            let b = DVec3::from_array(points[(k + 1) % points.len()]) - origin;
                            a.cross(b).dot(normal)
                        })
                        .sum::<f64>()
                };
                if signed_area(ring) * signed_area(original) <= 0.0 {
                    return false;
                }
                for k in 0..ring.len() {
                    let length = DVec3::from_array(ring[k])
                        .distance(DVec3::from_array(ring[(k + 1) % ring.len()]));
                    let old = DVec3::from_array(original[k])
                        .distance(DVec3::from_array(original[(k + 1) % original.len()]));
                    if length + EPS < old.min(self.min_edge) {
                        return false;
                    }
                }
            }
            let (u, v) = get_plane_basis(normal);
            if order_valid_rings(coords, u, v).is_none() {
                return false;
            }
        }
        // A local repair must not introduce a new overlap with another panel.
        for (&(i, j), &old) in &self.baseline_overlap {
            if (affected.contains(&i) || affected.contains(&j)) && self.overlap(i, j) > old + 1e-9 {
                return false;
            }
        }
        true
    }
    fn bounded(&self, id: usize, p: DVec3, tolerance: f64) -> bool {
        p.is_finite()
            && self.origins[id]
                .iter()
                .all(|o| o.distance(p) <= tolerance + EPS)
    }
    fn project(&self, start: DVec3, owners: &BTreeSet<usize>) -> Option<DVec3> {
        // Orthogonalize plane equations rather than inverting a near-singular matrix.
        let mut basis: Vec<(DVec3, f64)> = vec![];
        for &i in owners {
            let mut n = DVec3::from_array(self.panels[i].plane_normal);
            let mut rhs = -self.panels[i].plane_d - n.dot(start);
            for &(u, d) in &basis {
                let weight = n.dot(u);
                n -= weight * u;
                rhs -= weight * d;
            }
            let length = n.length();
            if length < EPS {
                if rhs.abs() > EPS {
                    return None;
                }
            } else {
                basis.push((n / length, rhs / length));
            }
        }
        let point = start + basis.iter().map(|(n, d)| *n * *d).sum::<DVec3>();
        if owners.iter().any(|&i| {
            (DVec3::from_array(self.panels[i].plane_normal).dot(point) + self.panels[i].plane_d)
                .abs()
                > EPS
        }) {
            None
        } else {
            Some(point)
        }
    }
    fn merge(&mut self, a: usize, b: usize, tolerance: f64) -> bool {
        let a = self.root(a);
        let b = self.root(b);
        if a == b || !self.owners[a].is_disjoint(&self.owners[b]) {
            return false;
        }
        let owners: BTreeSet<_> = self.owners[a].union(&self.owners[b]).copied().collect();
        let Some(point) = self.project((self.points[a] + self.points[b]) * 0.5, &owners) else {
            return false;
        };
        if !self.bounded(a, point, tolerance) || !self.bounded(b, point, tolerance) {
            return false;
        }
        let old = self.points[a];
        self.points[a] = point;
        self.parent[b] = a;
        if !self.valid(&owners) {
            self.points[a] = old;
            self.parent[b] = b;
            return false;
        }
        self.owners[a] = owners;
        let origins = self.origins[b].clone();
        self.origins[a].extend(origins);
        true
    }
    fn edges(&self) -> Vec<(usize, usize, usize)> {
        self.rings
            .iter()
            .enumerate()
            .flat_map(|(p, rings)| {
                rings.iter().flat_map(move |r| {
                    (0..r.len()).map(move |i| (p, self.root(r[i]), self.root(r[(i + 1) % r.len()])))
                })
            })
            .collect()
    }
    fn on_segment(&self, p: DVec3, a: DVec3, b: DVec3) -> (DVec3, f64) {
        let d = b - a;
        if d.length_squared() < EPS * EPS {
            return (a, 0.0);
        }
        let t = ((p - a).dot(d) / d.length_squared()).clamp(0.0, 1.0);
        (a + t * d, t)
    }
    fn snap_to_edge(
        &mut self,
        vertex: usize,
        a: usize,
        b: usize,
        target: usize,
        tolerance: f64,
        min_edge: f64,
    ) -> bool {
        let vertex = self.root(vertex);
        let a = self.root(a);
        let b = self.root(b);
        if vertex == a || vertex == b || self.owners[vertex].contains(&target) {
            return false;
        }
        let (near, t) = self.on_segment(self.points[vertex], self.points[a], self.points[b]);
        if t <= 0.0 || t >= 1.0 || near.distance(self.points[vertex]) > tolerance {
            return false;
        }
        let mut proposals = BTreeMap::new();
        let mut required = self.owners[vertex].clone();
        required.insert(target);
        // First try moving the vertex along the existing edge to all incident planes.
        let delta = self.points[b] - self.points[a];
        let mut exact_t = None;
        let mut possible = true;
        for &owner in &required {
            let n = DVec3::from_array(self.panels[owner].plane_normal);
            let distance = n.dot(self.points[a]) + self.panels[owner].plane_d;
            let denominator = n.dot(delta);
            if denominator.abs() < EPS {
                if distance.abs() > EPS {
                    possible = false;
                    break;
                }
            } else {
                let value = -distance / denominator;
                if exact_t.is_some_and(|old: f64| (value - old).abs() * delta.length() > EPS) {
                    possible = false;
                    break;
                }
                exact_t = Some(value);
            }
        }
        let point = if possible {
            let t = exact_t.unwrap_or(t);
            if t <= 0.0 || t >= 1.0 {
                return false;
            }
            self.points[a] + t * delta
        } else {
            // Shift both endpoints to the common plane intersection if their existing
            // incident planes and the displacement budget permit that change.
            for endpoint in [a, b] {
                let owners = self.owners[endpoint].union(&required).copied().collect();
                let Some(p) = self.project(self.points[endpoint], &owners) else {
                    return false;
                };
                if !self.bounded(endpoint, p, tolerance) {
                    return false;
                }
                proposals.insert(endpoint, p);
            }
            let (p, t) = self.on_segment(self.points[vertex], proposals[&a], proposals[&b]);
            if t <= 0.0 || t >= 1.0 {
                return false;
            }
            p
        };
        if !self.bounded(vertex, point, tolerance) {
            return false;
        }
        let pa = proposals.get(&a).copied().unwrap_or(self.points[a]);
        let pb = proposals.get(&b).copied().unwrap_or(self.points[b]);
        if point.distance(pa) < min_edge || point.distance(pb) < min_edge {
            return false;
        }
        if required.iter().any(|&i| {
            (DVec3::from_array(self.panels[i].plane_normal).dot(point) + self.panels[i].plane_d)
                .abs()
                > EPS
        }) {
            return false;
        }
        proposals.insert(vertex, point);
        let affected: BTreeSet<_> = proposals
            .keys()
            .flat_map(|&id| self.owners[id].iter().copied())
            .collect();
        let old: Vec<_> = proposals.keys().map(|&id| (id, self.points[id])).collect();
        for (&id, &p) in &proposals {
            self.points[id] = p;
        }
        if !self.valid(&affected) {
            for (id, p) in old {
                self.points[id] = p;
            }
            return false;
        }
        self.owners[vertex].insert(target);
        true
    }
}

pub fn conform_panels(panels: &mut [MacroPanel], config: &ReconstructionConfig) -> TopologySummary {
    let mut graph = Graph::new(panels, config.min_edge);
    let mut stats = TopologySummary::default();
    let tolerance = config.joint_tol;
    let mut pairs = Vec::new();
    for i in 0..graph.points.len() {
        for j in i + 1..graph.points.len() {
            let distance = graph.points[i].distance(graph.points[j]);
            if distance <= tolerance && graph.owners[i] != graph.owners[j] {
                pairs.push((distance, i, j));
            }
        }
    }
    pairs.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)).then(a.2.cmp(&b.2)));
    for (_, a, b) in pairs {
        if graph.merge(a, b, tolerance) {
            stats.merged_vertices += 1;
        }
    }
    let edges = graph.edges();
    let mut candidates = Vec::new();
    for id in 0..graph.points.len() {
        if graph.root(id) != id {
            continue;
        }
        for &(owner, a, b) in &edges {
            if graph.owners[id].contains(&owner) || id == a || id == b {
                continue;
            }
            let (point, t) = graph.on_segment(graph.points[id], graph.points[a], graph.points[b]);
            let distance = point.distance(graph.points[id]);
            if t > 0.0 && t < 1.0 && distance <= tolerance {
                candidates.push((distance, id, a, b, owner));
            }
        }
    }
    candidates.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)).then(a.4.cmp(&b.4)));
    for (_, id, a, b, owner) in candidates {
        if graph.snap_to_edge(id, a, b, owner, tolerance, config.min_edge) {
            stats.snapped_junctions += 1;
        }
    }
    // Node every coincident boundary, including opposite edge directions and T contacts.
    let mut output = graph.rings.clone();
    for (owner, rings) in graph.rings.iter().enumerate() {
        for (ri, ring) in rings.iter().enumerate() {
            let mut result = Vec::new();
            for i in 0..ring.len() {
                let a = graph.root(ring[i]);
                let b = graph.root(ring[(i + 1) % ring.len()]);
                result.push(a);
                let mut cuts = Vec::new();
                for id in 0..graph.points.len() {
                    if graph.root(id) != id || id == a || id == b {
                        continue;
                    }
                    let (p, t) =
                        graph.on_segment(graph.points[id], graph.points[a], graph.points[b]);
                    if t > 0.0 && t < 1.0 && p.distance(graph.points[id]) < EPS {
                        cuts.push((t, id));
                    }
                }
                cuts.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
                let mut previous = a;
                for (_, id) in cuts {
                    if graph.points[id].distance(graph.points[previous]) >= config.min_edge
                        && graph.points[id].distance(graph.points[b]) >= config.min_edge
                    {
                        result.push(id);
                        previous = id;
                        stats.inserted_vertices += 1;
                    }
                }
            }
            output[owner][ri] = result;
        }
    }
    graph.rings = output;
    let mut shared: BTreeMap<(usize, usize), BTreeSet<usize>> = BTreeMap::new();
    for (owner, a, b) in graph.edges() {
        shared
            .entry((a.min(b), a.max(b)))
            .or_default()
            .insert(owner);
    }
    let mut connections = vec![BTreeSet::new(); panels.len()];
    for owners in shared.values() {
        for &a in owners {
            for &b in owners {
                if a != b {
                    connections[a].insert(b);
                }
            }
        }
    }
    stats.connected_panel_pairs = connections.iter().map(|c| c.len()).sum::<usize>() / 2;
    let edges = graph.edges();
    for id in 0..graph.points.len() {
        if graph.root(id) != id {
            continue;
        }
        for &origin in &graph.origins[id] {
            stats.max_displacement = stats
                .max_displacement
                .max(origin.distance(graph.points[id]));
        }
        for &(owner, a, b) in &edges {
            if graph.owners[id].contains(&owner) || id == a || id == b {
                continue;
            }
            let (p, t) = graph.on_segment(graph.points[id], graph.points[a], graph.points[b]);
            let distance = p.distance(graph.points[id]);
            if t > 0.0 && t < 1.0 && distance > EPS && distance < tolerance {
                stats.unresolved_near_junctions += 1;
            }
        }
    }
    let coordinates: Vec<_> = (0..panels.len()).map(|i| graph.coords(i)).collect();
    drop(graph);
    let ids: Vec<_> = panels.iter().map(|p| p.id).collect();
    for (i, panel) in panels.iter_mut().enumerate() {
        panel.polygons = coordinates[i].clone();
        panel.connected_panel_ids = connections[i].iter().map(|&j| ids[j]).collect();
        panel.connected_panel_ids.sort_unstable();
    }
    stats
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::PanelType;
    fn panel(id: u32, normal: [f64; 3], d: f64, points: Vec<[f64; 3]>) -> MacroPanel {
        MacroPanel {
            id,
            panel_type: if normal[2] == 1.0 {
                PanelType::Slab
            } else {
                PanelType::Wall
            },
            stiffness_id: 1,
            plane_normal: normal,
            plane_d: d,
            polygons: vec![points],
            filled_holes: 0,
            filled_hole_area: 0.0,
            fe_count: 1,
            source_element_ids: vec![id],
            connected_panel_ids: vec![],
            constraint_points: vec![],
        }
    }
    fn slab(id: u32, z: f64) -> MacroPanel {
        panel(
            id,
            [0., 0., 1.],
            -z,
            vec![[0., 0., z], [1., 0., z], [1., 1., z], [0., 1., z]],
        )
    }
    #[test]
    fn near_t_junction_has_shared_edge_and_preserves_planes() {
        let mut panels = vec![
            slab(1, 0.),
            panel(
                2,
                [0., 1., 0.],
                -0.002,
                vec![
                    [0.25, 0.002, 0.003],
                    [0.75, 0.002, 0.003],
                    [0.75, 0.002, 1.],
                    [0.25, 0.002, 1.],
                ],
            ),
        ];
        let config = ReconstructionConfig::default();
        let stats = conform_panels(&mut panels, &config);
        assert!(stats.snapped_junctions >= 2);
        assert!(stats.inserted_vertices >= 2);
        assert_eq!(panels[0].connected_panel_ids, vec![2]);
        assert_eq!(panels[1].connected_panel_ids, vec![1]);
        assert!(stats.max_displacement <= config.joint_tol + EPS);
        for p in &panels {
            for point in p.polygons.iter().flatten() {
                assert!(
                    (DVec3::from_array(p.plane_normal).dot(DVec3::from_array(*point)) + p.plane_d)
                        .abs()
                        < EPS
                );
            }
        }
        let second = conform_panels(&mut panels, &config);
        assert_eq!(second.inserted_vertices, 0);
        assert_eq!(second.connected_panel_pairs, 1);
    }
    #[test]
    fn conflicting_parallel_planes_are_not_welded() {
        let mut panels = vec![slab(1, 0.), slab(2, 0.005)];
        let before = serde_json::to_value(&panels).unwrap();
        let stats = conform_panels(&mut panels, &ReconstructionConfig::default());
        assert_eq!(stats.merged_vertices, 0);
        assert_eq!(stats.connected_panel_pairs, 0);
        assert_eq!(before, serde_json::to_value(&panels).unwrap());
    }
    #[test]
    fn coplanar_gap_closes_without_creating_overlap() {
        let mut panels = vec![
            slab(1, 0.),
            panel(
                2,
                [0., 0., 1.],
                0.,
                vec![[1.005, 0., 0.], [2., 0., 0.], [2., 1., 0.], [1.005, 1., 0.]],
            ),
        ];
        let stats = conform_panels(&mut panels, &ReconstructionConfig::default());
        assert_eq!(stats.merged_vertices, 2);
        assert_eq!(stats.connected_panel_pairs, 1);
        let graph = Graph::new(&panels, 0.03);
        assert!(graph.overlap(0, 1) < 1e-10);
    }
    #[test]
    fn shared_corner_is_not_reported_as_shared_edge() {
        let mut panels = vec![
            slab(1, 0.),
            panel(
                2,
                [0., 0., 1.],
                0.,
                vec![[1., 1., 0.], [2., 1., 0.], [2., 2., 0.], [1., 2., 0.]],
            ),
        ];
        let stats = conform_panels(&mut panels, &ReconstructionConfig::default());
        assert_eq!(stats.merged_vertices, 1);
        assert_eq!(stats.connected_panel_pairs, 0);
    }
    #[test]
    fn validation_rejects_new_short_edge() {
        let panels = vec![slab(1, 0.)];
        let mut graph = Graph::new(&panels, 0.03);
        graph.points[1] = DVec3::new(0.02, 0., 0.);
        assert!(!graph.valid(&BTreeSet::from([0])));
    }
    #[test]
    fn validation_rejects_inverted_slender_face() {
        let panels = vec![panel(
            1,
            [0., 0., 1.],
            0.,
            vec![[0., 0., 0.], [1., 0., 0.], [0.5, 0.001, 0.]],
        )];
        let mut graph = Graph::new(&panels, 0.03);
        graph.points[2].y = -0.001;
        assert!(!graph.valid(&BTreeSet::from([0])));
    }
}
