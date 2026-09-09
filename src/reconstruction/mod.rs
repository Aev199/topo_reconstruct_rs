//! Planar surfaces on a shared topological graph. No implicit mechanical ties.
pub mod axes;
pub mod frame;
pub mod planes;
pub mod recognize;

use geo::{Contains, Intersects, Line, LineString, Polygon};
use glam::DVec3;
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, PartialEq)]
pub enum Error {
    InvalidPrecision,
    InvalidPlane,
    InvalidVertex,
    InvalidRing,
    NonPlanar,
    ShortEdge,
}

#[derive(Debug, Clone, Serialize)]
pub struct PlaneFrame {
    origin: [f64; 3],
    normal: [f64; 3],
    u: [f64; 3],
    v: [f64; 3],
}
impl PlaneFrame {
    pub fn new(origin: [f64; 3], normal: [f64; 3]) -> Result<Self, Error> {
        let o = DVec3::from_array(origin);
        let n = DVec3::from_array(normal);
        if !o.is_finite() || !n.is_finite() || !n.length().is_finite() || n.length() == 0.0 {
            return Err(Error::InvalidPlane);
        }
        let n = n.normalize();
        let helper = if n.z.abs() < 0.9 { DVec3::Z } else { DVec3::X };
        let u = helper.cross(n).normalize();
        Ok(Self {
            origin,
            normal: n.to_array(),
            u: u.to_array(),
            v: n.cross(u).to_array(),
        })
    }
    pub fn distance(&self, point: [f64; 3]) -> f64 {
        (DVec3::from_array(point) - DVec3::from_array(self.origin))
            .dot(DVec3::from_array(self.normal))
    }
    pub fn project(&self, point: [f64; 3]) -> [f64; 2] {
        let d = DVec3::from_array(point) - DVec3::from_array(self.origin);
        [
            d.dot(DVec3::from_array(self.u)),
            d.dot(DVec3::from_array(self.v)),
        ]
    }
    pub fn lift(&self, uv: [f64; 2]) -> [f64; 3] {
        (DVec3::from_array(self.origin)
            + uv[0] * DVec3::from_array(self.u)
            + uv[1] * DVec3::from_array(self.v))
        .to_array()
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct EdgeUse {
    pub edge: usize,
    pub reversed: bool,
}
#[derive(Debug, Clone, Serialize)]
pub struct Surface {
    pub plane: usize,
    /// Exterior first, followed by holes. Geometry lives in plane coordinates.
    pub contours: Vec<Vec<[f64; 2]>>,
    pub boundaries: Vec<Vec<EdgeUse>>,
    pub source_elements: Vec<u32>,
}
#[derive(Debug, Clone, Serialize)]
pub struct Model {
    precision: f64,
    minimum_edge: f64,
    planes: Vec<PlaneFrame>,
    vertices: Vec<[f64; 3]>,
    edges: Vec<[usize; 2]>,
    surfaces: Vec<Surface>,
    #[serde(skip)]
    edge_index: BTreeMap<[usize; 2], usize>,
}
impl Model {
    pub fn new(precision: f64, minimum_edge: f64) -> Result<Self, Error> {
        if !precision.is_finite()
            || !minimum_edge.is_finite()
            || precision <= 0.0
            || minimum_edge <= precision
        {
            return Err(Error::InvalidPrecision);
        }
        Ok(Self {
            precision,
            minimum_edge,
            planes: vec![],
            vertices: vec![],
            edges: vec![],
            surfaces: vec![],
            edge_index: BTreeMap::new(),
        })
    }
    pub fn add_plane(&mut self, plane: PlaneFrame) -> usize {
        let id = self.planes.len();
        self.planes.push(plane);
        id
    }
    /// Identity is explicit: coordinate proximity alone never creates a connection.
    pub fn add_vertex(&mut self, point: [f64; 3]) -> Result<usize, Error> {
        if !DVec3::from_array(point).is_finite() {
            return Err(Error::InvalidVertex);
        }
        let id = self.vertices.len();
        self.vertices.push(point);
        Ok(id)
    }
    pub fn surfaces(&self) -> &[Surface] {
        &self.surfaces
    }
    pub fn edges(&self) -> &[[usize; 2]] {
        &self.edges
    }
    pub fn planes(&self) -> &[PlaneFrame] {
        &self.planes
    }

    /// Transactional insertion. Reject warped input rather than projecting a
    /// shared vertex independently onto each owner's plane.
    pub fn add_surface(
        &mut self,
        plane: usize,
        rings: Vec<Vec<usize>>,
        mut source_elements: Vec<u32>,
    ) -> Result<usize, Error> {
        let frame = self.planes.get(plane).ok_or(Error::InvalidPlane)?;
        if rings.is_empty() {
            return Err(Error::InvalidRing);
        }
        let mut contours = vec![];
        for ring in &rings {
            if ring.len() < 3 || ring.iter().copied().collect::<BTreeSet<_>>().len() != ring.len() {
                return Err(Error::InvalidRing);
            }
            let mut uv = vec![];
            for (i, &id) in ring.iter().enumerate() {
                let p = *self.vertices.get(id).ok_or(Error::InvalidVertex)?;
                let q = *self
                    .vertices
                    .get(ring[(i + 1) % ring.len()])
                    .ok_or(Error::InvalidVertex)?;
                if frame.distance(p).abs() > self.precision {
                    return Err(Error::NonPlanar);
                }
                let a = frame.project(p);
                let b = frame.project(q);
                let length = (a[0] - b[0]).hypot(a[1] - b[1]);
                if !length.is_finite() || length < self.minimum_edge {
                    return Err(Error::ShortEdge);
                }
                uv.push(a);
            }
            validate_ring(&uv, self.precision)?;
            contours.push(uv);
        }
        let polygon = |r: &Vec<[f64; 2]>| {
            let points: Vec<_> = r.iter().chain(r.first()).map(|p| (p[0], p[1])).collect();
            Polygon::new(LineString::from(points), vec![])
        };
        let outer = polygon(&contours[0]);
        for i in 1..contours.len() {
            let hole = polygon(&contours[i]);
            if !outer.contains(&hole) || outer.exterior().intersects(hole.exterior()) {
                return Err(Error::InvalidRing);
            }
            for previous in &contours[1..i] {
                if polygon(previous).intersects(&hole) {
                    return Err(Error::InvalidRing);
                }
            }
        }
        // All validation precedes mutation, including edge interning.
        let boundaries = rings
            .iter()
            .map(|ring| {
                (0..ring.len())
                    .map(|i| {
                        let a = ring[i];
                        let b = ring[(i + 1) % ring.len()];
                        let key = [a.min(b), a.max(b)];
                        let edge = *self.edge_index.entry(key).or_insert_with(|| {
                            let id = self.edges.len();
                            self.edges.push(key);
                            id
                        });
                        EdgeUse {
                            edge,
                            reversed: a > b,
                        }
                    })
                    .collect()
            })
            .collect();
        source_elements.sort_unstable();
        source_elements.dedup();
        let id = self.surfaces.len();
        self.surfaces.push(Surface {
            plane,
            contours,
            boundaries,
            source_elements,
        });
        Ok(id)
    }
}

fn validate_ring(points: &[[f64; 2]], epsilon: f64) -> Result<(), Error> {
    let n = points.len();
    let line = |i: usize| {
        Line::new(
            (points[i][0], points[i][1]),
            (points[(i + 1) % n][0], points[(i + 1) % n][1]),
        )
    };
    let o = points[0];
    let area: f64 = (0..n)
        .map(|i| {
            let a = points[i];
            let b = points[(i + 1) % n];
            (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
        })
        .sum();
    if !area.is_finite() || area.abs() <= epsilon * epsilon {
        return Err(Error::InvalidRing);
    }
    for i in 0..n {
        for j in i + 1..n {
            if j == i + 1 || (i == 0 && j == n - 1) {
                continue;
            }
            if line(i).intersects(&line(j)) {
                return Err(Error::InvalidRing);
            }
        }
        let a = points[(i + n - 1) % n];
        let b = points[i];
        let c = points[(i + 1) % n];
        let cross = (a[0] - b[0]) * (c[1] - b[1]) - (a[1] - b[1]) * (c[0] - b[0]);
        let dot = (a[0] - b[0]) * (c[0] - b[0]) + (a[1] - b[1]) * (c[1] - b[1]);
        if cross.abs() <= epsilon * epsilon && dot > 0.0 {
            return Err(Error::InvalidRing);
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    fn square() -> (Model, usize, Vec<usize>) {
        let mut model = Model::new(1e-8, 0.03).unwrap();
        let plane = model.add_plane(PlaneFrame::new([0.; 3], [0., 0., 1.]).unwrap());
        let ring = [[0., 0., 0.], [2., 0., 0.], [2., 2., 0.], [0., 2., 0.]]
            .map(|p| model.add_vertex(p).unwrap())
            .to_vec();
        (model, plane, ring)
    }
    #[test]
    fn perpendicular_surfaces_reuse_one_boundary() {
        let (mut m, slab, r) = square();
        m.add_surface(slab, vec![r.clone()], vec![7]).unwrap();
        let wall = m.add_plane(PlaneFrame::new([0.; 3], [0., 1., 0.]).unwrap());
        let a = m.add_vertex([0., 0., 3.]).unwrap();
        let b = m.add_vertex([2., 0., 3.]).unwrap();
        m.add_surface(wall, vec![vec![r[1], r[0], a, b]], vec![8])
            .unwrap();
        assert_eq!(m.edges().len(), 7);
        let e = &m.surfaces()[0].boundaries[0][0];
        let f = &m.surfaces()[1].boundaries[0][0];
        assert_eq!(e.edge, f.edge);
        assert_ne!(e.reversed, f.reversed);
        for surface in m.surfaces() {
            let frame = &m.planes()[surface.plane];
            for &uv in surface.contours.iter().flatten() {
                assert!(frame.distance(frame.lift(uv)).abs() < 1e-12);
            }
        }
    }
    #[test]
    fn warped_surface_is_rejected_transactionally() {
        let (mut m, p, mut r) = square();
        r[2] = m.add_vertex([2., 2., 0.001]).unwrap();
        let before = serde_json::to_string(&m).unwrap();
        assert_eq!(m.add_surface(p, vec![r], vec![1]), Err(Error::NonPlanar));
        assert_eq!(before, serde_json::to_string(&m).unwrap());
    }
    #[test]
    fn rejects_self_intersections_and_short_edges() {
        let (mut m, p, r) = square();
        assert_eq!(
            m.add_surface(p, vec![vec![r[0], r[2], r[1], r[3]]], vec![]),
            Err(Error::InvalidRing)
        );
        let near = m.add_vertex([0.001, 0., 0.]).unwrap();
        assert_eq!(
            m.add_surface(p, vec![vec![r[0], near, r[2], r[3]]], vec![]),
            Err(Error::ShortEdge)
        );
        assert!(m.edges().is_empty());
    }
    #[test]
    fn holes_must_be_disjoint_and_strictly_inside() {
        let (mut m, p, r) = square();
        let hole = [
            [0.5, 0.5, 0.],
            [1.5, 0.5, 0.],
            [1.5, 1.5, 0.],
            [0.5, 1.5, 0.],
        ]
        .map(|q| m.add_vertex(q).unwrap())
        .to_vec();
        assert!(m
            .add_surface(p, vec![r.clone(), hole.clone()], vec![9, 8, 9])
            .is_ok());
        assert_eq!(m.surfaces()[0].source_elements, vec![8, 9]);
        let before = m.edges().len();
        assert_eq!(
            m.add_surface(p, vec![r.clone(), hole.clone(), hole], vec![]),
            Err(Error::InvalidRing)
        );
        let outside = [[3., 0., 0.], [4., 0., 0.], [4., 1., 0.], [3., 1., 0.]]
            .map(|q| m.add_vertex(q).unwrap())
            .to_vec();
        assert_eq!(
            m.add_surface(p, vec![r, outside], vec![]),
            Err(Error::InvalidRing)
        );
        assert_eq!(m.edges().len(), before);
    }
    #[test]
    fn rigid_transform_preserves_shared_topology() {
        let (base, _, _) = square();
        let rotation = glam::DQuat::from_axis_angle(DVec3::new(1., 2., 3.).normalize(), 0.73);
        let shift = DVec3::new(100., -250., 73.);
        let mut m = Model::new(1e-8, 0.03).unwrap();
        let p = m.add_plane(
            PlaneFrame::new(shift.to_array(), (rotation * DVec3::Z).to_array()).unwrap(),
        );
        let ids = base
            .vertices
            .iter()
            .map(|v| {
                m.add_vertex((rotation * DVec3::from_array(*v) + shift).to_array())
                    .unwrap()
            })
            .collect();
        m.add_surface(p, vec![ids], vec![1]).unwrap();
        assert_eq!(m.edges().len(), 4);
        for &uv in &m.surfaces()[0].contours[0] {
            assert!(m.planes()[p].distance(m.planes()[p].lift(uv)).abs() < 1e-10);
        }
    }
    #[test]
    fn invalid_numeric_inputs_are_rejected() {
        assert!(PlaneFrame::new([0.; 3], [0.; 3]).is_err());
        assert!(PlaneFrame::new([f64::NAN, 0., 0.], [0., 0., 1.]).is_err());
        assert!(Model::new(0., 0.03).is_err());
        assert!(Model::new(1e-8, f64::INFINITY).is_err());
    }
}
