#![allow(dead_code)]

use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReconstructionConfig {
    /// Допуск по высоте перекрытий и расстоянию (м)
    pub tol_dist: f64,
    /// Допуск по углам нормалей (радианы)
    pub tol_angle: f64,
    /// Допуск упрощения контуров полигонов (м)
    pub simplify_tol: f64,
    /// Радиус объединения узлов в единицах исходной модели.
    pub weld_tol: f64,
    /// Минимальная длина ребра при очистке контура.
    pub min_edge: f64,
    /// Допуск согласования границ восстановленных панелей.
    pub joint_tol: f64,
    /// Флаг: делить плиты по стенам
    pub split_slabs_by_walls: bool,
    /// Флаг: делить плиты по балкам
    pub split_slabs_by_beams: bool,
}

impl Default for ReconstructionConfig {
    fn default() -> Self {
        Self {
            tol_dist: 0.15,
            tol_angle: 0.08,
            simplify_tol: 0.01,
            weld_tol: 0.001,
            min_edge: 0.03,
            joint_tol: 0.01,
            split_slabs_by_walls: true,
            split_slabs_by_beams: true,
        }
    }
}
