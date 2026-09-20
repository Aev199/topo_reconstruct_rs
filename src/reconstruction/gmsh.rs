//! Stable interchange contract for an external Gmsh/OpenCASCADE backend.
//!
//! Rust owns engineering reconstruction and provenance.  The backend receives
//! explicit 3D surface rings, axis geometry and source/property ownership; it
//! must not infer structural meaning from coordinate proximity.

use super::assembly;
use glam::DVec3;
use serde::Serialize;

pub const FORMAT: &str = "topo-reconstruct-gmsh-v1";

#[derive(Debug, Clone, Serialize)]
pub struct Policy {
    /// Numerical precision of the reconstructed topology.
    pub precision: f64,
    /// Explicit engineering lower bound used by deterministic post-mesh
    /// micro-edge cleanup. This is not a Boolean fuzzy tolerance.
    pub minimum_edge: f64,
    /// Target Gmsh surface-mesh size in model length units.
    pub target_mesh_size: f64,
}

impl Policy {
    pub fn validate(&self) -> Result<(), &'static str> {
        if !self.precision.is_finite()
            || self.precision <= 0.0
            || !self.minimum_edge.is_finite()
            || self.minimum_edge <= self.precision
            || !self.target_mesh_size.is_finite()
            || self.target_mesh_size < self.minimum_edge
        {
            return Err("invalid Gmsh interchange policy");
        }
        Ok(())
    }
}

#[derive(Debug, Serialize)]
pub struct Interchange {
    pub format: &'static str,
    /// Coordinates are deliberately unit-agnostic. The consumer must preserve
    /// the model unit; no hidden conversion is allowed.
    pub length_unit: &'static str,
    pub policy: Policy,
    pub source_coverage_complete: bool,
    pub surfaces: Vec<Surface>,
    pub axes: Vec<Axis>,
    pub contacts: Vec<Contact>,
    pub blockers: Vec<String>,
}

#[derive(Debug, Serialize)]
pub struct Surface {
    pub source_surface: usize,
    pub source_patch: usize,
    pub stiffness: u32,
    pub source_elements: Vec<u32>,
    /// Exterior ring first, then holes, all in global 3D coordinates.
    pub rings: Vec<Vec<[f64; 3]>>,
}

#[derive(Debug, Serialize)]
pub struct Axis {
    pub source_axis: usize,
    pub endpoints: [[f64; 3]; 2],
    pub property_spans: Vec<PropertySpan>,
    pub anchors: Vec<AxisAnchor>,
}

#[derive(Debug, Serialize)]
pub struct PropertySpan {
    pub source_element: u32,
    pub stiffness: u32,
    pub start_t: f64,
    pub end_t: f64,
}

#[derive(Debug, Serialize)]
pub struct AxisAnchor {
    pub source_node: u32,
    pub t: f64,
    pub point: [f64; 3],
}

#[derive(Debug, Clone, Copy, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum ContactLocation {
    Boundary,
    Interior,
}

impl From<assembly::bars::Location> for ContactLocation {
    fn from(value: assembly::bars::Location) -> Self {
        match value {
            assembly::bars::Location::Boundary => Self::Boundary,
            assembly::bars::Location::Interior => Self::Interior,
        }
    }
}

#[derive(Debug, Serialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum Contact {
    Point {
        axis: usize,
        surface: usize,
        t: f64,
        point: [f64; 3],
        location: ContactLocation,
    },
    Interval {
        axis: usize,
        surface: usize,
        start_t: f64,
        end_t: f64,
        endpoints: [[f64; 3]; 2],
        location: ContactLocation,
    },
}

fn point_on_axis(endpoints: [[f64; 3]; 2], t: f64) -> Result<[f64; 3], &'static str> {
    if !t.is_finite() || t < 0.0 || t > 1.0 {
        return Err("invalid axis contact parameter");
    }
    let a = DVec3::from_array(endpoints[0]);
    let b = DVec3::from_array(endpoints[1]);
    Ok(a.lerp(b, t).to_array())
}

/// Build a versioned, backend-neutral geometry package from a completed v2
/// assembly. No topology is repaired here.
pub fn from_assembly(
    source: &assembly::Report,
    target_mesh_size: f64,
) -> Result<Interchange, &'static str> {
    let policy = Policy {
        precision: source.policy.precision,
        minimum_edge: source.policy.minimum_edge,
        target_mesh_size,
    };
    policy.validate()?;

    let model = &source.preview;
    if model.surfaces.len() != source.surface_source_patches.len()
        || model.surfaces.len() != source.surface_stiffness.len()
    {
        return Err("surface provenance arrays do not match reconstructed surfaces");
    }

    let mut surfaces = Vec::with_capacity(model.surfaces.len());
    for (index, surface) in model.surfaces.iter().enumerate() {
        let plane = model
            .planes
            .get(surface.plane)
            .ok_or("invalid surface plane reference")?;
        let mut rings = Vec::with_capacity(surface.contours.len());
        for contour in &surface.contours {
            if contour.len() < 3 {
                return Err("invalid Gmsh surface ring");
            }
            rings.push(contour.iter().map(|&uv| plane.lift(uv)).collect());
        }
        surfaces.push(Surface {
            source_surface: index,
            source_patch: source.surface_source_patches[index],
            stiffness: source.surface_stiffness[index],
            source_elements: surface.source_elements.clone(),
            rings,
        });
    }

    let mut axes = Vec::with_capacity(source.axis_assembly.axes.len());
    for (axis_index, axis) in source.axis_assembly.axes.iter().enumerate() {
        if axis.source_axis != axis_index {
            // source_axis is a semantic provenance ID and need not be dense,
            // so this is deliberately not an error. Keep the explicit ID.
        }
        let endpoints = axis
            .endpoints
            .map(|vertex| *model.vertices.get(vertex).ok_or("invalid axis endpoint vertex"))
            .into_iter()
            .collect::<Result<Vec<_>, _>>()?;
        let endpoints: [[f64; 3]; 2] = endpoints
            .try_into()
            .map_err(|_| "invalid axis endpoint count")?;
        if DVec3::from_array(endpoints[0]).distance(DVec3::from_array(endpoints[1]))
            <= source.policy.precision
        {
            return Err("degenerate Gmsh axis");
        }

        let mut anchors = Vec::with_capacity(axis.anchors.len());
        for anchor in &axis.anchors {
            let point = *model
                .vertices
                .get(anchor.vertex)
                .ok_or("invalid axis anchor vertex")?;
            anchors.push(AxisAnchor {
                source_node: anchor.source_node,
                t: anchor.t,
                point,
            });
        }
        axes.push(Axis {
            source_axis: axis.source_axis,
            endpoints,
            property_spans: axis
                .spans
                .iter()
                .map(|span| PropertySpan {
                    source_element: span.element,
                    stiffness: span.stiffness,
                    start_t: span.start_t,
                    end_t: span.end_t,
                })
                .collect(),
            anchors,
        });
    }

    let mut contacts = Vec::with_capacity(source.axis_assembly.contacts.len());
    for item in &source.axis_assembly.contacts {
        match *item {
            assembly::bars::Contact::Point {
                axis,
                surface,
                t,
                location,
                ..
            } => {
                let geometry = axes.get(axis).ok_or("invalid contact axis")?;
                if surface >= surfaces.len() {
                    return Err("invalid contact surface");
                }
                contacts.push(Contact::Point {
                    axis,
                    surface,
                    t,
                    point: point_on_axis(geometry.endpoints, t)?,
                    location: location.into(),
                });
            }
            assembly::bars::Contact::Interval {
                axis,
                surface,
                start_t,
                end_t,
                location,
            } => {
                let geometry = axes.get(axis).ok_or("invalid contact axis")?;
                if surface >= surfaces.len() || end_t <= start_t {
                    return Err("invalid interval contact");
                }
                contacts.push(Contact::Interval {
                    axis,
                    surface,
                    start_t,
                    end_t,
                    endpoints: [
                        point_on_axis(geometry.endpoints, start_t)?,
                        point_on_axis(geometry.endpoints, end_t)?,
                    ],
                    location: location.into(),
                });
            }
        }
    }

    let mut blockers = Vec::new();
    if !source.all_surface_patches_built || !source.issues.is_empty() {
        blockers.push("unresolved_surface_assembly".into());
    }
    if !source.axis_assembly.all_axes_built || !source.axis_assembly.issues.is_empty() {
        blockers.push("unresolved_axis_assembly".into());
    }

    Ok(Interchange {
        format: FORMAT,
        length_unit: "model_unit",
        policy,
        source_coverage_complete: blockers.is_empty(),
        surfaces,
        axes,
        contacts,
        blockers,
    })
}
