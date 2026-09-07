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
        let mut panels =
            PanelReconstructor::reconstruct(self.mesh_data, &canonical_nodes, self.config);

        let topology = crate::geometry::topology::conform_panels(&mut panels, self.config);

        // 3. Определение высотных отметок
        let slab_elevations =
            PanelReconstructor::extract_slab_elevations(&panels, self.config.tol_dist);

        // 4. Восстановление стержней
        let mut bars = BarReconstructor::reconstruct(
            self.mesh_data,
            &canonical_nodes,
            &slab_elevations,
            self.config,
        );

        let bar_contacts = crate::geometry::bar_contacts::attach_bar_endpoints(
            &mut bars,
            &mut panels,
            self.mesh_data,
            &canonical_nodes,
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
        let mut diagnostics = vec![format!("Согласование границ: объединено вершин {}, вставлено {}, связанных пар панелей {}; максимальное перемещение {}", topology.merged_vertices, topology.inserted_vertices, topology.connected_panel_pairs, topology.max_displacement)];
        diagnostics.push(format!("Примыкания стержней: привязано групп {}, перемещено {}, не согласовано {}; максимум перемещения {}",bar_contacts.attached_groups,bar_contacts.moved_groups,bar_contacts.rejected_groups,bar_contacts.max_displacement));
        if !bar_contacts.joint_solve_converged {
            diagnostics.push("Совместное согласование осей не сошлось в заданном допуске; предварительные независимые смещения отменены. Требуется проверка оставшихся примыканий.".into());
        }
        if !bar_contacts.unresolved_bar_junctions.is_empty() {
            diagnostics.push(format!(
                "Прямые оси с несогласованными исходными примыканиями: {:?}",
                bar_contacts.unresolved_bar_junctions
            ));
        }
        let represented_bars: std::collections::HashSet<_> = bars
            .iter()
            .flat_map(|b| b.source_element_ids.iter().copied())
            .collect();
        let missing_bars: Vec<_> = self
            .mesh_data
            .elements
            .iter()
            .filter(|e| e.is_bar() && !represented_bars.contains(&e.id))
            .map(|e| e.id)
            .collect();
        if !missing_bars.is_empty() {
            diagnostics.push(format!(
                "Стержневые КЭ без восстановленной оси: {:?}",
                missing_bars
            ));
        }
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
            topology,
            bar_contacts,
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
