use crate::config::ReconstructionConfig;
use crate::geometry::utils::canonicalize_nodes;
use crate::models::{MeshData, PanelType, BarType, ReconstructionReport};
use super::bars::BarReconstructor;
use super::panels::PanelReconstructor;

pub struct TopologyPipeline<'a> {
    mesh_data: &'a MeshData,
    config: &'a ReconstructionConfig,
}

impl<'a> TopologyPipeline<'a> {
    pub fn new(mesh_data: &'a MeshData, config: &'a ReconstructionConfig) -> Self {
        Self { mesh_data, config }
    }

    pub fn run(&self) -> ReconstructionReport {
        // 1. Канонизация узлов
        let canonical_nodes = canonicalize_nodes(&self.mesh_data.nodes, self.config.canonical_precision);

        // 2. Восстановление плит и стен (параллельно через Rayon)
        let panels = PanelReconstructor::reconstruct(self.mesh_data, &canonical_nodes, self.config);

        // 3. Определение высотных отметок
        let slab_elevations = PanelReconstructor::extract_slab_elevations(&panels, self.config.tol_dist);

        // 4. Восстановление стержней
        let bars = BarReconstructor::reconstruct(self.mesh_data, &canonical_nodes, &slab_elevations, self.config);

        ReconstructionReport {
            slabs_count: panels.iter().filter(|p| p.panel_type == PanelType::Slab).count(),
            walls_count: panels.iter().filter(|p| p.panel_type == PanelType::Wall).count(),
            inclined_panels_count: panels.iter().filter(|p| p.panel_type == PanelType::InclinedPanel).count(),
            columns_count: bars.iter().filter(|b| b.bar_type == BarType::Column).count(),
            beams_count: bars.iter().filter(|b| b.bar_type == BarType::Beam).count(),
            braces_count: bars.iter().filter(|b| b.bar_type == BarType::Brace).count(),
            panels,
            bars,
        }
    }
}