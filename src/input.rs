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
