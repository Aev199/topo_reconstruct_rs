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

/// Two source nodes merged by `MeshData::weld_unconnected`.
#[derive(Debug, Clone, serde::Serialize)]
pub struct Weld {
    pub dropped: u32,
    pub kept: u32,
    pub distance: f64,
}

impl MeshData {
    /// Merge nodes closer than `tolerance` that share no element: cracks and
    /// duplicated nodes left by mesh conversion. Nodes of one element are
    /// never merged, so no element collapses, and every node of a merged
    /// group stays within `tolerance` of the kept node. Returns every merge.
    pub fn weld_unconnected(&mut self, tolerance: f64) -> Vec<Weld> {
        use std::collections::{BTreeMap, BTreeSet};
        if !(tolerance > 0.) || !tolerance.is_finite() {
            return vec![];
        }
        let mut elements_of = BTreeMap::<u32, BTreeSet<usize>>::new();
        for (k, e) in self.elements.iter().enumerate() {
            for &n in &e.nodes {
                elements_of.entry(n).or_default().insert(k);
            }
        }
        let used: Vec<u32> = elements_of
            .keys()
            .copied()
            .filter(|n| self.nodes.contains_key(n))
            .collect();
        let cell = |p: DVec3| [p.x, p.y, p.z].map(|x| (x / tolerance).floor() as i64);
        let mut grid = BTreeMap::<[i64; 3], Vec<u32>>::new();
        for &n in &used {
            grid.entry(cell(self.nodes[&n])).or_default().push(n);
        }
        let mut pairs = vec![];
        for &a in &used {
            let pa = self.nodes[&a];
            let c = cell(pa);
            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        for &b in grid
                            .get(&[c[0] + dx, c[1] + dy, c[2] + dz])
                            .into_iter()
                            .flatten()
                        {
                            if a >= b {
                                continue;
                            }
                            let d = pa.distance(self.nodes[&b]);
                            if d <= tolerance && elements_of[&a].is_disjoint(&elements_of[&b]) {
                                pairs.push((d, a, b));
                            }
                        }
                    }
                }
            }
        }
        pairs.sort_by(|x, y| x.0.total_cmp(&y.0).then(x.1.cmp(&y.1)).then(x.2.cmp(&y.2)));
        // Groups: representative -> members and their elements.
        let mut group = BTreeMap::<u32, u32>::new();
        let mut members = BTreeMap::<u32, Vec<u32>>::new();
        let mut group_elements = BTreeMap::<u32, BTreeSet<usize>>::new();
        let find = |group: &BTreeMap<u32, u32>, n: u32| *group.get(&n).unwrap_or(&n);
        let mut welds = vec![];
        for (_, a, b) in pairs {
            let (ra, rb) = (find(&group, a), find(&group, b));
            if ra == rb {
                continue;
            }
            let elements = |r: u32| {
                group_elements
                    .get(&r)
                    .cloned()
                    .unwrap_or_else(|| elements_of[&r].clone())
            };
            let (ea, eb) = (elements(ra), elements(rb));
            if !ea.is_disjoint(&eb) {
                continue;
            }
            // Keep the representative with more elements (then lower id).
            let (keep, drop) =
                if (eb.len(), std::cmp::Reverse(rb)) > (ea.len(), std::cmp::Reverse(ra)) {
                    (rb, ra)
                } else {
                    (ra, rb)
                };
            let all: Vec<u32> = [ra, rb]
                .iter()
                .flat_map(|&r| members.get(&r).cloned().unwrap_or_else(|| vec![r]))
                .collect();
            let pk = self.nodes[&keep];
            if all.iter().any(|n| self.nodes[n].distance(pk) > tolerance) {
                continue;
            }
            for &n in &all {
                if n != keep {
                    group.insert(n, keep);
                }
            }
            members.remove(&drop);
            members.insert(keep, all.clone());
            let mut union = ea;
            union.extend(eb);
            group_elements.remove(&drop);
            group_elements.insert(keep, union);
            welds.push(Weld {
                dropped: drop,
                kept: keep,
                distance: self.nodes[&drop].distance(pk),
            });
        }
        // Resolve chains: every member points at its final representative.
        let resolve: Vec<(u32, u32)> = group
            .keys()
            .map(|&n| {
                let mut r = n;
                while let Some(&next) = group.get(&r) {
                    if next == r {
                        break;
                    }
                    r = next;
                }
                (n, r)
            })
            .collect();
        let map: BTreeMap<u32, u32> = resolve.into_iter().collect();
        for e in &mut self.elements {
            for n in &mut e.nodes {
                if let Some(&k) = map.get(n) {
                    *n = k;
                }
            }
        }
        for n in map.keys() {
            self.nodes.remove(n);
        }
        welds
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn quad(id: u32, nodes: [u32; 4]) -> ElementData {
        ElementData {
            id,
            elem_type: 44,
            stiff_id: 1,
            nodes: nodes.to_vec(),
        }
    }

    #[test]
    fn crack_between_elements_is_welded_but_elements_never_collapse() {
        let mut mesh = MeshData::default();
        // Two quads meeting at x = 1 with duplicated, slightly offset nodes,
        // and a third quad with two of its own nodes 1 mm apart.
        for (id, p) in [
            (1, [0., 0., 0.]),
            (2, [1., 0., 0.]),
            (3, [1., 1., 0.]),
            (4, [0., 1., 0.]),
            (5, [1.0015, 0., 0.]),
            (6, [2., 0., 0.]),
            (7, [2., 1., 0.]),
            (8, [1., 1.001, 0.]),
            (9, [3., 0., 0.]),
            (10, [3.001, 0., 0.]),
            (11, [3., 1., 0.]),
        ] {
            mesh.nodes.insert(id, DVec3::from_array(p));
        }
        mesh.elements = vec![
            quad(1, [1, 2, 3, 4]),
            quad(2, [5, 6, 7, 8]),
            quad(3, [9, 10, 11, 7]),
        ];
        let welds = mesh.weld_unconnected(0.005);
        assert_eq!(welds.len(), 2, "{welds:?}");
        assert!(welds.iter().all(|w| w.distance <= 0.005));
        // The two quads now share an edge; the third keeps nodes 9 and 10.
        let shared: Vec<_> = mesh.elements[0]
            .nodes
            .iter()
            .filter(|n| mesh.elements[1].nodes.contains(n))
            .collect();
        assert_eq!(shared.len(), 2);
        assert!(mesh.elements[2].nodes.contains(&9) && mesh.elements[2].nodes.contains(&10));
        for e in &mesh.elements {
            let unique: std::collections::BTreeSet<_> = e.nodes.iter().collect();
            assert_eq!(unique.len(), e.nodes.len());
            assert!(e.nodes.iter().all(|n| mesh.nodes.contains_key(n)));
        }
        // Idempotent.
        assert!(mesh.weld_unconnected(0.005).is_empty());
    }
}
