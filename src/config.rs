use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReconstructionConfig {
    /// Допуск по высоте перекрытий и расстоянию (м)
    pub tol_dist: f64,
    /// Допуск по углам нормалей (радианы)
    pub tol_angle: f64,
    /// Допуск упрощения контуров полигонов (м)
    pub simplify_tol: f64,
    /// Округление координат для сшивки дубликатов узлов (количество знаков)
    pub canonical_precision: u32,
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
            canonical_precision: 3,
            split_slabs_by_walls: true,
            split_slabs_by_beams: true,
        }
    }
}