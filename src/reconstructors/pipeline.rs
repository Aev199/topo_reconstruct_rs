use super::bars::BarReconstructor;
use super::panels::PanelReconstructor;
use crate::config::ReconstructionConfig;
use crate::geometry::utils::canonicalize_nodes;
use crate::models::{BarType, MeshData, PanelType, ReconstructionReport};

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
        let canonical_nodes = canonicalize_nodes(&self.mesh_data.nodes, self.config.weld_tol);

        // 2. Восстановление плит и стен (параллельно через Rayon)
        let panels = PanelReconstructor::reconstruct(self.mesh_data, &canonical_nodes, self.config);

        // 3. Определение высотных отметок
        let slab_elevations =
            PanelReconstructor::extract_slab_elevations(&panels, self.config.tol_dist);

        // 4. Восстановление стержней
        let bars = BarReconstructor::reconstruct(
            self.mesh_data,
            &canonical_nodes,
            &slab_elevations,
            self.config,
        );

        let represented: std::collections::HashSet<u32> = panels
            .iter()
            .flat_map(|p| p.source_element_ids.iter().copied())
            .collect();
        let missing: Vec<u32> = self
            .mesh_data
            .elements
            .iter()
            .filter(|e| e.is_shell() && !represented.contains(&e.id))
            .map(|e| e.id)
            .collect();
        let mut diagnostics = Vec::new();
        let mut excluded = std::collections::BTreeMap::<u32, Vec<u32>>::new();
        for el in &self.mesh_data.elements {
            if !el.is_shell() && !el.is_bar() {
                excluded.entry(el.elem_type).or_default().push(el.id);
            }
        }
        for (kind, ids) in excluded {
            diagnostics.push(format!("Тип КЭ {}: {} элементов сохранены во входных данных, не реконструируются как поверхности/стержни; ID: {:?}", kind, ids.len(), ids));
        }
        if !missing.is_empty() {
            diagnostics.push(format!(
                "Поддерживаемые оболочки без восстановленной поверхности: {:?}",
                missing
            ));
        }
        for panel in &panels {
            if panel.filled_holes > 0 {
                diagnostics.push(format!("Панель {}: закрыто узких проемов {}; площадь {}; ширина не более {} в единицах модели", panel.id,panel.filled_holes,panel.filled_hole_area,self.config.simplify_tol));
            }
            let short = panel.polygons.iter().any(|ring| {
                (0..ring.len()).any(|i| {
                    glam::DVec3::from_array(ring[i])
                        .distance(glam::DVec3::from_array(ring[(i + 1) % ring.len()]))
                        < self.config.min_edge
                })
            });
            if short {
                diagnostics.push(format!(
                    "Панель {}: остались короткие ребра; требуется проверка перед перебивкой",
                    panel.id
                ));
            }
        }
        let welded = canonical_nodes
            .iter()
            .filter(|(id, canonical)| id != canonical)
            .count();
        diagnostics.push(format!(
            "Объединено узлов: {}; допуск: {} в единицах исходной модели",
            welded, self.config.weld_tol
        ));
        ReconstructionReport {
            diagnostics,
            slabs_count: panels
                .iter()
                .filter(|p| p.panel_type == PanelType::Slab)
                .count(),
            walls_count: panels
                .iter()
                .filter(|p| p.panel_type == PanelType::Wall)
                .count(),
            inclined_panels_count: panels
                .iter()
                .filter(|p| p.panel_type == PanelType::InclinedPanel)
                .count(),
            columns_count: bars
                .iter()
                .filter(|b| b.bar_type == BarType::Column)
                .count(),
            beams_count: bars.iter().filter(|b| b.bar_type == BarType::Beam).count(),
            braces_count: bars.iter().filter(|b| b.bar_type == BarType::Brace).count(),
            panels,
            bars,
        }
    }
}
