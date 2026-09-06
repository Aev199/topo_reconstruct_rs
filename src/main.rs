#![allow(dead_code, unused_imports, unused_variables)]

mod config;
mod exporters;
mod geometry;
mod models;
mod parsers;
mod reconstructors;

use clap::Parser;
use config::ReconstructionConfig;
use exporters::{DxfExporter, JsonExporter};
use parsers::LiraParser;
use reconstructors::TopologyPipeline;
use std::time::Instant;

#[derive(Parser, Debug)]
#[command(
    author,
    version,
    about = "Реконструкция BIM-топологии из КЭ-сеток ПК ЛИРА (High-Performance Rust Core)"
)]
struct Args {
    /// Путь к текстовому файлу расчетной схемы ЛИРЫ (.txt)
    #[arg(default_value = "скала3.txt")]
    input: String,

    /// Радиус сшивки узлов (в единицах исходной модели)
    #[arg(long, default_value_t = 0.001)]
    weld_tol: f64,

    /// Допуск объединения плоскостей (в единицах исходной модели)
    #[arg(long, default_value_t = 0.15)]
    plane_tol: f64,

    /// Допуск дотягивания и упрощения контуров (в единицах модели)
    #[arg(long, default_value_t = 0.01)]
    simplify_tol: f64,

    /// Порог коротких ребер для попытки упрощения
    #[arg(long, default_value_t = 0.03)]
    min_edge: f64,

    /// Путь для сохранения JSON отчета
    #[arg(short, long, default_value = "building_topology.json")]
    json: String,

    /// Путь для сохранения DXF файла
    #[arg(short, long, default_value = "building_topology.dxf")]
    dxf: String,
}

fn main() {
    let args = Args::parse();
    let mut config = ReconstructionConfig::default();
    for value in [
        args.weld_tol,
        args.plane_tol,
        args.simplify_tol,
        args.min_edge,
    ] {
        if !value.is_finite() || value <= 0.0 {
            eprintln!("Допуски должны быть конечными положительными числами.");
            std::process::exit(2);
        }
    }
    config.weld_tol = args.weld_tol;
    config.tol_dist = args.plane_tol;
    config.simplify_tol = args.simplify_tol;
    config.min_edge = args.min_edge;

    println!("=== RECONSTRUCT TOPOLOGY CORE (RUST) ===");
    println!("Входной файл: {}", args.input);

    let total_start = Instant::now();

    // 1. Чтение и парсинг
    println!("\n1. Чтение и параллельный парсинг расчетной схемы...");
    let parse_start = Instant::now();
    let mesh_data = match LiraParser::parse(&args.input) {
        Ok(data) => {
            println!(
                "   [OK] Загружено за {:.2?}: {} узлов, {} КЭ",
                parse_start.elapsed(),
                data.nodes.len(),
                data.elements.len()
            );
            data
        }
        Err(e) => {
            eprintln!("   [ERROR] Ошибка при чтении файла: {}", e);
            std::process::exit(1);
        }
    };

    // 2. Реконструкция топологии
    println!("\n2. Параллельная реконструкция макроэлементов (Rayon)...");
    let recon_start = Instant::now();
    let pipeline = TopologyPipeline::new(&mesh_data, &config);
    let report = pipeline.run();
    println!(
        "   [OK] Реконструкция завершена за {:.2?}",
        recon_start.elapsed()
    );

    println!("\n--- РЕЗУЛЬТАТЫ РЕКОНСТРУКЦИИ ---");
    println!("   Плит перекрытий (Slabs):   {}", report.slabs_count);
    println!("   Стен / пилонов (Walls):     {}", report.walls_count);
    println!(
        "   Наклонных панелей:          {}",
        report.inclined_panels_count
    );
    println!("   Колонн (Columns):           {}", report.columns_count);
    println!("   Балок (Beams):              {}", report.beams_count);
    println!("   Связей / Раскосов (Braces): {}", report.braces_count);

    for message in &report.diagnostics {
        eprintln!("   [DIAGNOSTIC] {}", message);
    }

    // 3. Экспорт
    println!("\n3. Экспорт результатов...");
    let exp_start = Instant::now();
    if let Err(e) = JsonExporter::export(&report, &args.json) {
        eprintln!("   [ERROR] Ошибка экспорта JSON: {}", e);
    } else {
        println!("   [OK] JSON сохранен: {}", args.json);
    }

    if let Err(e) = DxfExporter::export(&report, &args.dxf) {
        eprintln!("   [ERROR] Ошибка экспорта DXF: {}", e);
    } else {
        println!("   [OK] CAD DXF сохранен: {}", args.dxf);
    }
    println!("   Экспорт выполнен за {:.2?}", exp_start.elapsed());

    println!("\n==========================================");
    println!("ИТОГОВОЕ ВРЕМЯ ВЫПОЛНЕНИЯ: {:.2?}", total_start.elapsed());
    println!("==========================================");
}
