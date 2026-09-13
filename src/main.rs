#![allow(dead_code, unused_imports, unused_variables)]

mod config;
mod exporters;
mod geometry;
mod input;
mod models;
mod parsers;
mod reconstructors;

use clap::Parser;
use config::ReconstructionConfig;
use exporters::{DxfExporter, JsonExporter};
use parsers::LiraParser;
use reconstructors::TopologyPipeline;
use std::{fs::File, time::Instant};

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

    /// Допуск согласования геометрии: 0.05 при координатах в метрах
    #[arg(long, default_value_t = 0.05)]
    joint_tol: f64,

    /// Верхний предел адаптивного ремонта контуров: 0.15 при координатах в метрах
    #[arg(long, default_value_t = 0.15)]
    repair_max_tol: f64,

    /// Путь для сохранения JSON отчета
    #[arg(short, long, default_value = "building_topology.json")]
    json: String,

    /// Путь для сохранения DXF файла
    #[arg(short, long, default_value = "building_topology.dxf")]
    dxf: String,

    /// Экспериментальный v2: записать recognize/frame/assembly JSON и не запускать legacy pipeline
    #[arg(long, value_name = "PATH")]
    v2_preview_json: Option<String>,

    /// Число итераций совместного v2 frame solver
    #[arg(long, default_value_t = 1000)]
    v2_iterations: usize,
}

fn run_v2_preview(input: &str, output: &str, iterations: usize) -> Result<(), Box<dyn std::error::Error>> {
    use topo_reconstruct_rs::{
        parsers::LiraParser as V2LiraParser,
        reconstruction::{assembly, frame, graph, planes, recognize},
    };

    if iterations == 0 {
        return Err("--v2-iterations must be positive".into());
    }
    let mesh = V2LiraParser::parse(input)?;
    let axes = recognize::recognize(
        &mesh,
        &recognize::Policy {
            angle: 0.02,
            line_tolerance: 0.01,
            numerical_precision: 1e-8,
        },
    )?;
    let plane_report = planes::recognize(
        &mesh,
        &planes::Policy {
            angle: 0.02,
            distance: 0.01,
            precision: 1e-8,
        },
    )?;
    let result = frame::solve(
        &mesh,
        &axes,
        &plane_report,
        &frame::Policy {
            up: [0., 0., 1.],
            angle: 0.02,
            maximum_movement: 0.15,
            relative_movement: 0.05,
            minimum_length: 0.03,
            residual_tolerance: 1e-7,
            iterations,
        },
    )?;
    let topology = assembly::assemble(
        &mesh,
        &result,
        &assembly::Policy {
            closure_tolerance: 0.001,
            junction_movement_limit: 0.05,
            precision: 1e-7,
            minimum_edge: 0.001,
        },
    )?;
    let report = serde_json::json!({
        "constraint_graph": graph::Graph::from_frame(&result),
        "frame": result,
        "topology": topology,
        "axis_recognition": axes,
        "plane_recognition": plane_report,
    });
    if output == "-" {
        serde_json::to_writer_pretty(std::io::stdout().lock(), &report)?;
        println!();
    } else {
        serde_json::to_writer_pretty(File::create(output)?, &report)?;
    }
    Ok(())
}

fn main() {
    let args = Args::parse();
    let mut config = ReconstructionConfig::default();
    for value in [
        args.weld_tol,
        args.plane_tol,
        args.simplify_tol,
        args.min_edge,
        args.joint_tol,
        args.repair_max_tol,
    ] {
        if !value.is_finite() || value <= 0.0 {
            eprintln!("Допуски должны быть конечными положительными числами.");
            std::process::exit(2);
        }
    }

    if let Some(path) = &args.v2_preview_json {
        if let Err(error) = run_v2_preview(&args.input, path, args.v2_iterations) {
            eprintln!("[V2 PREVIEW ERROR] {error}");
            std::process::exit(1);
        }
        eprintln!("[V2 PREVIEW] JSON сохранен: {path}");
        return;
    }

    config.weld_tol = args.weld_tol;
    config.tol_dist = args.plane_tol;
    config.simplify_tol = args.simplify_tol;
    config.min_edge = args.min_edge;
    config.joint_tol = args.joint_tol;
    config.repair_max_tol = args.repair_max_tol;

    println!("=== RECONSTRUCT TOPOLOGY CORE (RUST) ===");
    println!("Входной файл: {}", args.input);

    let total_start = Instant::now();

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
