use glam::DVec3;
use std::collections::BTreeMap;
use topo_reconstruct_rs::input::{ElementData, MeshData};

/// Exact FE fragment: slab with an opening and two properties, two walls,
/// a beam with a property change, and a column at an interior beam node.
pub fn source() -> MeshData {
    let mut mesh = MeshData::default();
    let mut nodes = BTreeMap::new();
    let mut add = |coordinates: Vec<(i32, i32, i32)>, kind, stiffness| {
        let ids = coordinates
            .into_iter()
            .map(|p| {
                *nodes.entry(p).or_insert_with(|| {
                    let id = mesh.nodes.len() as u32 + 1;
                    mesh.nodes
                        .insert(id, DVec3::new(p.0 as f64, p.1 as f64, p.2 as f64));
                    id
                })
            })
            .collect();
        mesh.elements.push(ElementData {
            id: mesh.elements.len() as u32 + 1,
            elem_type: kind,
            stiff_id: stiffness,
            nodes: ids,
        });
    };
    for x in 0..6 {
        for y in 0..4 {
            if (x, y) == (4, 2) {
                continue;
            }
            add(
                vec![(x, y, 0), (x + 1, y, 0), (x + 1, y + 1, 0), (x, y + 1, 0)],
                44,
                if x < 3 { 10 } else { 20 },
            );
        }
    }
    for x in 0..6 {
        for z in 0..2 {
            add(
                vec![(x, 0, z), (x + 1, 0, z), (x + 1, 0, z + 1), (x, 0, z + 1)],
                44,
                30,
            );
        }
    }
    for y in 0..4 {
        for z in 0..2 {
            add(
                vec![(0, y, z), (0, y + 1, z), (0, y + 1, z + 1), (0, y, z + 1)],
                44,
                40,
            );
        }
    }
    for x in 0..6 {
        add(
            vec![(x, 1, 0), (x + 1, 1, 0)],
            10,
            if x < 3 { 50 } else { 60 },
        );
    }
    for z in -2..0 {
        add(vec![(2, 1, z), (2, 1, z + 1)], 10, 70);
    }
    mesh
}
