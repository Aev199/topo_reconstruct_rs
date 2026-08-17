use crate::models::{BarType, PanelType, ReconstructionReport};
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

pub struct DxfExporter;

impl DxfExporter {
    /// Высокоскоростной прямой генератор CAD DXF R2010 (3D Polylines / Lines)
    pub fn export<P: AsRef<Path>>(report: &ReconstructionReport, filepath: P) -> io::Result<()> {
        let file = File::create(filepath)?;
        let mut w = BufWriter::new(file);

        // Заголовок и таблица слоев
        writeln!(w, "0\nSECTION\n2\nHEADER\n0\nENDSEC")?;
        writeln!(w, "0\nSECTION\n2\nTABLES\n0\nTABLE\n2\nLAYER\n70\n6")?;
        
        Self::write_layer(&mut w, "SLABS", 1)?;     // Красный
        Self::write_layer(&mut w, "WALLS", 3)?;     // Зеленый
        Self::write_layer(&mut w, "INCLINED", 2)?;  // Желтый
        Self::write_layer(&mut w, "COLUMNS", 5)?;   // Синий
        Self::write_layer(&mut w, "BEAMS", 4)?;     // Голубой
        Self::write_layer(&mut w, "BRACES", 6)?;    // Пурпурный

        writeln!(w, "0\nENDTAB\n0\nENDSEC")?;

        // Секция сущностей ENTITIES
        writeln!(w, "0\nSECTION\n2\nENTITIES")?;

        // 1. Отрисовка стержней (LINE)
        for b in &report.bars {
            let layer = match b.bar_type {
                BarType::Column => "COLUMNS",
                BarType::Beam => "BEAMS",
                BarType::Brace => "BRACES",
            };
            writeln!(w, "0\nLINE\n8\n{}\n10\n{}\n20\n{}\n30\n{}\n11\n{}\n21\n{}\n31\n{}",
                layer,
                b.start_point[0], b.start_point[1], b.start_point[2],
                b.end_point[0], b.end_point[1], b.end_point[2]
            )?;
        }

        // 2. Отрисовка чистых замкнутых 3D-полилиний (POLYLINE3D)
        for p in &report.panels {
            let layer = match p.panel_type {
                PanelType::Slab => "SLABS",
                PanelType::Wall => "WALLS",
                PanelType::InclinedPanel => "INCLINED",
            };

            for poly in &p.polygons {
                if poly.len() < 3 {
                    continue;
                }

                // Заголовок 3D-полилинии (70 -> 9: 8 (3D) + 1 (Closed))
                writeln!(w, "0\nPOLYLINE\n8\n{}\n66\n1\n70\n9", layer)?;

                for pt in poly {
                    // Вершина 3D-полилинии (70 -> 32: 3D polyline vertex)
                    writeln!(w, "0\nVERTEX\n8\n{}\n10\n{}\n20\n{}\n30\n{}\n70\n32", layer, pt[0], pt[1], pt[2])?;
                }

                writeln!(w, "0\nSEQEND")?;
            }
        }

        writeln!(w, "0\nENDSEC\n0\nEOF")?;
        w.flush()?;
        Ok(())
    }

    fn write_layer<W: Write>(w: &mut W, name: &str, color: i16) -> io::Result<()> {
        writeln!(w, "0\nLAYER\n2\n{}\n70\n0\n62\n{}\n6\nCONTINUOUS", name, color)
    }
}