#![allow(dead_code)]

use glam::DVec3;
use hashbrown::HashMap;
use serde::{Deserialize, Serialize};

pub use crate::input::{ElementData, MeshData};

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
    pub filled_holes: usize,
    pub filled_hole_area: f64,
    pub fe_count: usize,
    pub source_element_ids: Vec<u32>,
    pub connected_panel_ids: Vec<u32>,
    pub constraint_points: Vec<[f64; 3]>,
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
    pub stiffness_id: Option<u32>,
    pub start_point: [f64; 3],
    pub end_point: [f64; 3],
    pub length: f64,
    pub start_node_id: u32,
    pub end_node_id: u32,
    pub source_element_ids: Vec<u32>,
    pub source_node_ids: Vec<u32>,
    pub start_panel_ids: Vec<u32>,
    pub end_panel_ids: Vec<u32>,
    pub constraints: Vec<BarConstraint>,
    pub property_spans: Vec<BarPropertySpan>,
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
    pub diagnostics: Vec<String>,
    pub bar_contacts: crate::geometry::bar_contacts::BarContactSummary,
    pub topology: crate::geometry::topology::TopologySummary,
    pub panels: Vec<MacroPanel>,
    pub bars: Vec<MacroBar>,
}

/// Meshing/attachment metadata on a two-endpoint geometric line; t is in [0,1].
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BarConstraint {
    pub t: f64,
    pub source_node_id: Option<u32>,
    pub kinds: Vec<String>,
    pub panel_ids: Vec<u32>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BarPropertySpan {
    pub start_t: f64,
    pub end_t: f64,
    pub stiffness_id: u32,
    pub source_element_ids: Vec<u32>,
}
