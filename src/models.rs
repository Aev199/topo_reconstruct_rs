#![allow(dead_code)]

use glam::DVec3;
use hashbrown::HashMap;
use serde::{Deserialize, Serialize};

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

/// Тип макроэлемента панели
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum PanelType {
    Slab,
    Wall,
    InclinedPanel,
}

/// Восстановленная макропанель (плита, стена)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MacroPanel {
    pub id: u32,
    pub panel_type: PanelType,
    pub stiffness_id: u32,
    pub plane_normal: [f64; 3],
    pub plane_d: f64,
    /// Список 3D-контуров: [0] — внешний периметр, [1..] — внутренние проемы
    pub polygons: Vec<Vec<[f64; 3]>>,
    pub fe_count: usize,
    pub connected_panel_ids: Vec<u32>,
}

/// Тип стержневого макроэлемента
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum BarType {
    Column,
    Beam,
    Brace,
}

/// Восстановленный стержень
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MacroBar {
    pub bar_type: BarType,
    pub stiffness_id: u32,
    pub start_point: [f64; 3],
    pub end_point: [f64; 3],
    pub length: f64,
}

/// Итоговый сводный отчет
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReconstructionReport {
    pub slabs_count: usize,
    pub walls_count: usize,
    pub inclined_panels_count: usize,
    pub columns_count: usize,
    pub beams_count: usize,
    pub braces_count: usize,
    pub panels: Vec<MacroPanel>,
    pub bars: Vec<MacroBar>,
}