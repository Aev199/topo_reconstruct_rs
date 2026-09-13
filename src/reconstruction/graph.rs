//! Dependency graph for joint geometric solving. Indices refer to the frame;
//! coordinates, properties and source provenance remain owned by that frame.
//! Only explicit source incidences and accepted support families connect objects.
use super::frame;
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize)]
#[serde(tag = "kind", content = "index", rename_all = "snake_case")]
pub enum Entity {
    Vertex(usize),
    Axis(usize),
    Plane(usize),
}

#[derive(Debug, Serialize)]
pub struct Component {
    pub vertices: Vec<usize>,
    pub axes: Vec<usize>,
    pub planes: Vec<usize>,
}

#[derive(Debug, Serialize)]
pub struct Graph {
    /// Structural dependencies, not proximity-based or mechanical ties.
    pub relations: Vec<[Entity; 2]>,
    pub components: Vec<Component>,
}

impl Graph {
    pub fn from_frame(frame: &frame::Report) -> Self {
        let entities = (0..frame.node_ids.len())
            .map(Entity::Vertex)
            .chain((0..frame.axes.len()).map(Entity::Axis))
            .chain((0..frame.surfaces.len()).map(Entity::Plane));
        let mut relations = Vec::new();
        for (i, axis) in frame.axes.iter().enumerate() {
            for node in axis
                .endpoints
                .into_iter()
                .chain(axis.anchors.iter().map(|a| a.node))
            {
                relations.push([Entity::Axis(i), Entity::Vertex(node)]);
            }
        }
        for (i, surface) in frame.surfaces.iter().enumerate() {
            for &node in &surface.nodes {
                relations.push([Entity::Plane(i), Entity::Vertex(node)]);
            }
        }
        for family in &frame.plane_families {
            for pair in family.windows(2) {
                relations.push([Entity::Plane(pair[0]), Entity::Plane(pair[1])]);
            }
        }
        Self::build(entities, relations)
    }

    fn build(entities: impl Iterator<Item = Entity>, relations: Vec<[Entity; 2]>) -> Self {
        let mut neighbors: BTreeMap<Entity, BTreeSet<Entity>> =
            entities.map(|e| (e, BTreeSet::new())).collect();
        let relations: BTreeSet<_> = relations
            .into_iter()
            .map(|[a, b]| if a <= b { [a, b] } else { [b, a] })
            .collect();
        for &[a, b] in &relations {
            neighbors.entry(a).or_default().insert(b);
            neighbors.entry(b).or_default().insert(a);
        }
        let mut seen = BTreeSet::new();
        let mut components = Vec::new();
        for &seed in neighbors.keys() {
            if seen.contains(&seed) {
                continue;
            }
            let mut pending = vec![seed];
            let mut members = BTreeSet::new();
            while let Some(entity) = pending.pop() {
                if seen.insert(entity) {
                    members.insert(entity);
                    pending.extend(&neighbors[&entity]);
                }
            }
            let mut component = Component {
                vertices: vec![],
                axes: vec![],
                planes: vec![],
            };
            for member in members {
                match member {
                    Entity::Vertex(i) => component.vertices.push(i),
                    Entity::Axis(i) => component.axes.push(i),
                    Entity::Plane(i) => component.planes.push(i),
                }
            }
            components.push(component);
        }
        Self {
            relations: relations.into_iter().collect(),
            components,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn joint_couples_whole_axes_and_all_incident_planes() {
        use Entity::*;
        let entities = vec![
            Vertex(0),
            Vertex(1),
            Vertex(2),
            Axis(0),
            Axis(1),
            Plane(0),
            Plane(1),
        ];
        let graph = Graph::build(
            entities.into_iter(),
            vec![
                [Axis(0), Vertex(0)],
                [Axis(0), Vertex(1)],
                [Axis(1), Vertex(1)],
                [Axis(1), Vertex(2)],
                [Plane(0), Vertex(0)],
                [Plane(1), Vertex(2)],
            ],
        );
        assert_eq!(graph.components.len(), 1);
        assert_eq!(graph.components[0].axes, vec![0, 1]);
        assert_eq!(graph.components[0].planes, vec![0, 1]);
        assert_eq!(graph.components[0].vertices, vec![0, 1, 2]);
    }

    #[test]
    fn no_proximity_connection_and_order_independent() {
        use Entity::*;
        let entities = vec![Vertex(0), Vertex(1), Plane(0), Plane(1), Axis(0)];
        let relations = vec![[Plane(0), Vertex(0)], [Plane(1), Vertex(1)]];
        let graph = Graph::build(entities.clone().into_iter(), relations.clone());
        assert_eq!(graph.components.len(), 3);
        let reversed = Graph::build(
            entities.into_iter().rev(),
            relations.into_iter().rev().collect(),
        );
        assert_eq!(
            serde_json::to_value(&graph).unwrap(),
            serde_json::to_value(&reversed).unwrap()
        );
    }
}
