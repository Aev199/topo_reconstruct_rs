//! Shared source FE data; independent of reconstruction versions.
use glam::DVec3;
use hashbrown::HashMap;

/// Сырые данные расчетной схемы КЭ
#[derive(Debug, Clone, Default)]
pub struct MeshData {
    /// Узлы: {node_id: DVec3(x, y, z)}
    pub nodes: HashMap<u32, DVec3>,
    /// Конечные элементы
    pub elements: Vec<ElementData>,
}

/// Описание конечного элемента
#[derive(Debug, Clone)]
pub struct ElementData {
    pub id: u32,
    pub elem_type: u32,
    pub stiff_id: u32,
    pub nodes: Vec<u32>,
}

/// Explicit supported LIRA element families; unknown types are retained, not guessed.
impl ElementData {
    pub fn is_shell(&self) -> bool {
        matches!((self.elem_type, self.nodes.len()), (41 | 44, 4) | (42, 3))
    }

    pub fn is_bar(&self) -> bool {
        self.elem_type == 10 && self.nodes.len() == 2
    }
}

/// Node pairs a non-shell element ties together (bars, rigid links,
/// two-node elastic links such as LIRA type 55), plus `groups` of nodes
/// moving together (absolutely rigid bodies: master node first). Structures
/// tied this way are joined, never left apart as at an expansion joint.
pub fn node_links(mesh: &MeshData, groups: &[Vec<u32>]) -> Vec<[u32; 2]> {
    let mut links: Vec<[u32; 2]> = mesh
        .elements
        .iter()
        .filter(|e| !e.is_shell() && e.nodes.len() >= 2)
        .flat_map(|e| e.nodes[1..].iter().map(move |&n| [e.nodes[0], n]))
        .chain(
            groups
                .iter()
                .filter(|g| g.len() >= 2)
                .flat_map(|g| g[1..].iter().map(move |&n| [g[0], n])),
        )
        .filter(|[a, b]| a != b)
        .map(|[a, b]| [a.min(b), a.max(b)])
        .collect();
    links.sort_unstable();
    links.dedup();
    links
}

/// Bar elements whose numerically given stiffness (EF, EIy, EIz, ...) no
/// real section has: the radius of gyration sqrt(EI / EF) exceeds both the
/// element length and 1 m (model units, metres; a solid section 3 m deep
/// has 0.87 m, building members a few centimetres to decimetres). These
/// are rigid links of the analysis model (LIRA "1000 200000 200000
/// 200000": 14 m for links about 1 m long), not structures. A short real
/// bar (5 cm, 8 cm radius) is not a link.
pub fn rigid_links(mesh: &MeshData, stiffness: &HashMap<u32, Vec<f64>>) -> Vec<u32> {
    mesh.elements
        .iter()
        .filter(|e| e.is_bar())
        .filter(|e| {
            let Some(values) = stiffness.get(&e.stiff_id) else {
                return false;
            };
            let (Some(&area), Some(&iy), Some(&iz)) =
                (values.first(), values.get(1), values.get(2))
            else {
                return false;
            };
            let (Some(a), Some(b)) = (mesh.nodes.get(&e.nodes[0]), mesh.nodes.get(&e.nodes[1]))
            else {
                return false;
            };
            let radius = (iy.max(iz) / area).sqrt();
            area > 0. && radius > 1. && radius > a.distance(*b)
        })
        .map(|e| e.id)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn links_tie_bars_springs_and_rigid_bodies_but_not_shells() {
        let element = |id, elem_type, nodes: Vec<u32>| ElementData {
            id,
            elem_type,
            stiff_id: 1,
            nodes,
        };
        let mesh = MeshData {
            nodes: Default::default(),
            elements: vec![
                element(1, 44, vec![1, 2, 3, 4]),
                element(2, 10, vec![5, 4]),
                element(3, 55, vec![6, 7]),
                element(4, 56, vec![8]),
            ],
        };
        let links = node_links(&mesh, &[vec![9, 10, 11]]);
        assert_eq!(links, vec![[4, 5], [6, 7], [9, 10], [9, 11]]);
    }

    #[test]
    fn bar_with_a_gyration_radius_beyond_its_length_is_a_rigid_link() {
        let mut mesh = MeshData::default();
        mesh.nodes.insert(1, DVec3::ZERO);
        mesh.nodes.insert(2, DVec3::new(1., 0., 0.));
        mesh.nodes.insert(3, DVec3::new(0., 0., 3.));
        mesh.nodes.insert(4, DVec3::new(0.05, 0., 0.));
        for (id, stiff_id, nodes) in [
            (1, 1, [1, 2]),
            (2, 2, [1, 3]),
            (3, 3, [1, 2]),
            (4, 2, [1, 4]),
        ] {
            mesh.elements.push(ElementData {
                id,
                elem_type: 10,
                stiff_id,
                nodes: nodes.to_vec(),
            });
        }
        let stiffness = HashMap::from([
            // A link: sqrt(200000 / 1000) = 14 m over 1 m.
            (1, vec![1000., 200000., 200000., 200000., 0., 0.]),
            // A column 3 m long: 8 cm.
            (2, vec![750000., 4687.5, 4687.5, 1118.11, 0., 0.]),
        ]);
        // Element 3 has a section type (not numeric), element 4 is a 5 cm
        // piece of the column section (8 cm radius): both kept.
        assert_eq!(rigid_links(&mesh, &stiffness), vec![1]);
    }
}
