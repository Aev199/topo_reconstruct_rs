//! Global synchronization of mesh constraints before local surface meshing.
//!
//! Every interval endpoint is materialized once on its axis.  Boundary edge
//! chains are then built from those same vertices and are shared by every
//! surface owning the edge.  This keeps surface constraints and bar chains
//! conforming without welding vertices by proximity.

use super::{parameter, point, sorted, subdivide, Policy};
use crate::reconstruction::{
    assembly::bars::{Axis, Contact},
    Model,
};
use glam::DVec3;
use serde::Serialize;
use std::collections::BTreeSet;

#[derive(Debug, Serialize)]
pub struct Report {
    pub shared_edge_count: usize,
    pub interval_contact_count: usize,
    pub synchronized_interval_endpoints: usize,
    pub endpoint_bindings: Vec<EndpointBinding>,
    pub axis_node_count: usize,
    pub edge_node_count: usize,
}

#[derive(Debug, Serialize)]
pub struct EndpointBinding {
    pub axis: usize,
    pub surface: usize,
    pub contact: usize,
    pub role: EndpointRole,
    pub parameter: f64,
    pub vertex: usize,
    pub generated: bool,
}

#[derive(Debug, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum EndpointRole {
    Start,
    End,
}

pub(super) struct Synchronized {
    pub axis_nodes: Vec<Vec<(f64, usize)>>,
    pub edge_nodes: Vec<Vec<(f64, usize)>>,
    pub report: Report,
}

fn axis_geometry(axis: &Axis, vertices: &[[f64; 3]]) -> Result<(DVec3, DVec3, f64), &'static str> {
    let [a, b] = axis.endpoints;
    let (Some(a), Some(b)) = (vertices.get(a), vertices.get(b)) else {
        return Err("invalid axis endpoint vertex");
    };
    let a = point(a);
    let b = point(b);
    let length = a.distance(b);
    if !length.is_finite() || length == 0. {
        return Err("degenerate axis constraint");
    }
    Ok((a, b, length))
}

fn normalized_parameter(t: f64, tolerance: f64) -> Result<f64, &'static str> {
    if !t.is_finite() || t < -tolerance || t > 1. + tolerance {
        return Err("axis interval parameter outside endpoints");
    }
    Ok(t.clamp(0., 1.))
}

fn insert_axis_node(
    chain: &mut Vec<(f64, usize)>,
    t: f64,
    vertex: usize,
    vertices: &[[f64; 3]],
    tolerance: f64,
    axis_start: DVec3,
    axis_end: DVec3,
    precision: f64,
) -> Result<bool, &'static str> {
    let axis = axis_end - axis_start;
    let p = point(
        vertices
            .get(vertex)
            .ok_or("invalid axis constraint vertex")?,
    );
    if p.distance(axis_start + axis * t) > precision {
        return Err("axis constraint vertex is off axis");
    }
    if let Some((u, _existing)) = chain.iter().find(|(_, n)| *n == vertex) {
        if (*u - t).abs() > tolerance {
            return Err("axis constraint vertex has inconsistent parameter");
        }
        return Ok(false);
    }
    if let Some((_, existing)) = chain.iter().find(|(u, _)| (*u - t).abs() <= tolerance) {
        if *existing != vertex {
            return Err("distinct axis constraint vertices share a parameter");
        }
        return Ok(false);
    }
    chain.push((t, vertex));
    Ok(true)
}

fn validate_surface(surface: usize, model: &Model) -> Result<(), &'static str> {
    model
        .surfaces
        .get(surface)
        .map(|_| ())
        .ok_or("invalid contact surface")
}

pub(super) fn synchronize(
    model: &Model,
    axes: &[Axis],
    contacts: &[Contact],
    vertices: &mut Vec<[f64; 3]>,
    policy: &Policy,
    precision: f64,
) -> std::result::Result<Synchronized, &'static str> {
    let mut axis_nodes = Vec::with_capacity(axes.len());
    let mut geometries = Vec::with_capacity(axes.len());
    for axis in axes {
        let geometry = axis_geometry(axis, vertices)?;
        let mut chain: Vec<_> = axis.anchors.iter().map(|a| (a.t, a.vertex)).collect();
        sorted(&mut chain, vertices, precision)?;
        if chain.len() < 2 {
            return Err("empty axis constraint");
        }
        geometries.push(geometry);
        axis_nodes.push(chain);
    }

    let mut endpoint_bindings = Vec::new();
    let mut interval_contact_count = 0;
    for (contact, item) in contacts.iter().enumerate() {
        match item {
            Contact::Point {
                axis,
                surface,
                vertex,
                t,
                ..
            } => {
                validate_surface(*surface, model)?;
                let (start, end, length) =
                    *geometries.get(*axis).ok_or("invalid point contact axis")?;
                let tolerance = precision / length;
                let t = normalized_parameter(*t, tolerance)?;
                insert_axis_node(
                    &mut axis_nodes[*axis],
                    t,
                    *vertex,
                    vertices,
                    tolerance,
                    start,
                    end,
                    precision,
                )?;
            }
            Contact::Interval {
                axis,
                surface,
                start_t,
                end_t,
                ..
            } => {
                validate_surface(*surface, model)?;
                let (start, end, length) = *geometries
                    .get(*axis)
                    .ok_or("invalid interval contact axis")?;
                let tolerance = precision / length;
                let start_t = normalized_parameter(*start_t, tolerance)?;
                let end_t = normalized_parameter(*end_t, tolerance)?;
                if end_t - start_t <= tolerance {
                    return Err("degenerate interval contact");
                }
                interval_contact_count += 1;
                for (role, t) in [(EndpointRole::Start, start_t), (EndpointRole::End, end_t)] {
                    let existing = axis_nodes[*axis]
                        .iter()
                        .find(|(u, _)| (*u - t).abs() <= tolerance)
                        .map(|(_, vertex)| *vertex);
                    let (vertex, generated) = if let Some(vertex) = existing {
                        (vertex, false)
                    } else {
                        let vertex = vertices.len();
                        vertices.push(start.lerp(end, t).to_array());
                        (vertex, true)
                    };
                    insert_axis_node(
                        &mut axis_nodes[*axis],
                        t,
                        vertex,
                        vertices,
                        tolerance,
                        start,
                        end,
                        precision,
                    )?;
                    endpoint_bindings.push(EndpointBinding {
                        axis: *axis,
                        surface: *surface,
                        contact,
                        role,
                        parameter: t,
                        vertex,
                        generated,
                    });
                }
            }
        }
    }

    for chain in &mut axis_nodes {
        sorted(chain, vertices, precision)?;
        // A synchronized interval endpoint is part of the global chain before
        // spacing is applied, so every owner receives the same split vertex.
        let expanded = subdivide(
            chain,
            vertices,
            policy.boundary_spacing,
            policy.maximum_added_vertices_per_surface,
            precision,
        )?;
        *chain = expanded;
    }

    let mut owners = vec![BTreeSet::new(); model.edges.len()];
    for (surface, item) in model.surfaces.iter().enumerate() {
        for edge_id in item
            .boundaries
            .iter()
            .flatten()
            .map(|edge| edge.edge)
            .chain(item.junctions.iter().copied())
        {
            let owner = owners
                .get_mut(edge_id)
                .ok_or("invalid surface edge reference")?;
            owner.insert(surface);
        }
    }
    let shared_edge_count = owners.iter().filter(|owner| owner.len() > 1).count();
    let mut edge_nodes = Vec::with_capacity(model.edges.len());
    for (edge_id, &[a, b]) in model.edges.iter().enumerate() {
        let pa = point(vertices.get(a).ok_or("invalid model edge vertex")?);
        let pb = point(vertices.get(b).ok_or("invalid model edge vertex")?);
        let mut chain = vec![(0., a), (1., b)];
        for item in contacts {
            let (surface, candidates) = match *item {
                Contact::Point {
                    vertex, surface, ..
                } => (surface, vec![vertex]),
                Contact::Interval {
                    axis,
                    surface,
                    start_t,
                    end_t,
                    ..
                } => (
                    surface,
                    axis_nodes[axis]
                        .iter()
                        .filter(|(t, _)| *t >= start_t && *t <= end_t)
                        .map(|(_, vertex)| *vertex)
                        .collect(),
                ),
            };
            if !owners[edge_id].contains(&surface) {
                continue;
            }
            for vertex in candidates {
                if let Some(t) = parameter(point(&vertices[vertex]), pa, pb, precision) {
                    chain.push((t, vertex));
                }
            }
        }
        sorted(&mut chain, vertices, precision)?;
        edge_nodes.push(subdivide(
            &chain,
            vertices,
            policy.boundary_spacing,
            policy.maximum_added_vertices_per_surface,
            precision,
        )?);
    }

    // A subdivision introduced on a common boundary is also a subdivision of
    // every bar interval lying on that boundary.  Insert it into the one
    // global axis chain, not into a surface-local copy.
    for item in contacts {
        let Contact::Interval {
            axis,
            surface,
            start_t,
            end_t,
            ..
        } = *item
        else {
            continue;
        };
        let [a, b] = axes[axis].endpoints;
        for edge_id in model.surfaces[surface]
            .boundaries
            .iter()
            .flatten()
            .map(|edge| edge.edge)
            .chain(model.surfaces[surface].junctions.iter().copied())
        {
            for &(_, vertex) in &edge_nodes[edge_id] {
                if let Some(t) = parameter(
                    point(&vertices[vertex]),
                    point(&vertices[a]),
                    point(&vertices[b]),
                    precision,
                ) {
                    let length = point(&vertices[a]).distance(point(&vertices[b]));
                    let tolerance = precision / length;
                    if t >= start_t - tolerance && t <= end_t + tolerance {
                        insert_axis_node(
                            &mut axis_nodes[axis],
                            t,
                            vertex,
                            vertices,
                            tolerance,
                            point(&vertices[a]),
                            point(&vertices[b]),
                            precision,
                        )?;
                    }
                }
            }
        }
    }
    for chain in &mut axis_nodes {
        sorted(chain, vertices, precision)?;
    }

    let axis_node_count = axis_nodes.iter().map(Vec::len).sum();
    let edge_node_count = edge_nodes.iter().map(Vec::len).sum();
    Ok(Synchronized {
        axis_nodes,
        edge_nodes,
        report: Report {
            shared_edge_count,
            interval_contact_count,
            synchronized_interval_endpoints: endpoint_bindings.len(),
            endpoint_bindings,
            axis_node_count,
            edge_node_count,
        },
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reconstruction::{
        assembly::bars::{Anchor, Location},
        PlaneFrame,
    };

    fn shared_edge_case() -> (Model, Vec<[f64; 3]>, Vec<Axis>, Vec<Contact>) {
        let mut model = Model::new(1e-8, 0.01).unwrap();
        let plane = model.add_plane(PlaneFrame::new([0., 0., 0.], [0., 0., 1.]).unwrap());
        let points = vec![
            [0., 0., 0.],
            [4., 0., 0.],
            [4., 2., 0.],
            [0., 2., 0.],
            [4., -2., 0.],
            [0., -2., 0.],
        ];
        let vertices: Vec<_> = points
            .iter()
            .map(|&p| model.add_vertex(p).unwrap())
            .collect();
        model
            .add_surface(plane, vec![vec![0, 1, 2, 3]], vec![1])
            .unwrap();
        model
            .add_surface(plane, vec![vec![1, 0, 5, 4]], vec![2])
            .unwrap();
        let axis = Axis {
            source_axis: 0,
            endpoints: [vertices[0], vertices[1]],
            anchors: vec![
                Anchor {
                    source_node: 10,
                    vertex: vertices[0],
                    t: 0.,
                },
                Anchor {
                    source_node: 11,
                    vertex: vertices[1],
                    t: 1.,
                },
            ],
            spans: vec![],
        };
        let contacts = vec![
            Contact::Interval {
                axis: 0,
                surface: 0,
                start_t: 0.25,
                end_t: 0.75,
                location: Location::Boundary,
            },
            Contact::Interval {
                axis: 0,
                surface: 1,
                start_t: 0.25,
                end_t: 0.75,
                location: Location::Boundary,
            },
        ];
        (model, points, vec![axis], contacts)
    }

    fn policy() -> Policy {
        Policy {
            boundary_spacing: 10.,
            maximum_area: 1.,
            minimum_angle_degrees: 20.,
            maximum_added_vertices_per_surface: 32,
        }
    }

    #[test]
    fn interval_endpoints_are_created_once_and_shared_by_surfaces_and_bar() {
        let (model, mut vertices, axes, contacts) = shared_edge_case();
        let synced = synchronize(&model, &axes, &contacts, &mut vertices, &policy(), 1e-8).unwrap();

        assert_eq!(synced.report.shared_edge_count, 1);
        assert_eq!(synced.report.interval_contact_count, 2);
        assert_eq!(synced.report.synchronized_interval_endpoints, 4);
        assert_eq!(synced.report.endpoint_bindings.len(), 4);
        assert_eq!(
            synced
                .report
                .endpoint_bindings
                .iter()
                .filter(|binding| binding.generated)
                .count(),
            2
        );
        let first = &synced.report.endpoint_bindings[0..2];
        let second = &synced.report.endpoint_bindings[2..4];
        assert_eq!(first[0].vertex, second[0].vertex);
        assert_eq!(first[1].vertex, second[1].vertex);
        assert_eq!(
            synced.axis_nodes[0]
                .iter()
                .map(|(_, vertex)| *vertex)
                .collect::<Vec<_>>(),
            vec![0, first[0].vertex, first[1].vertex, 1]
        );
        let common_edge = synced.edge_nodes[0]
            .iter()
            .map(|(_, vertex)| *vertex)
            .collect::<Vec<_>>();
        assert_eq!(common_edge, vec![0, first[0].vertex, first[1].vertex, 1]);
        assert_eq!(synced.report.axis_node_count, 4);
        assert_eq!(synced.report.edge_node_count, 16);
    }

    #[test]
    fn malformed_interval_is_rejected_before_any_mesh_constraint_is_added() {
        let (model, mut vertices, axes, mut contacts) = shared_edge_case();
        contacts[0] = Contact::Interval {
            axis: 0,
            surface: 0,
            start_t: 0.75,
            end_t: 0.75,
            location: Location::Interior,
        };
        assert!(matches!(
            synchronize(&model, &axes, &contacts, &mut vertices, &policy(), 1e-8),
            Err("degenerate interval contact")
        ));
        assert_eq!(vertices.len(), 6);
    }
}
