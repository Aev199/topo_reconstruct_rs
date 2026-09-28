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

    /// Экспериментальный v2: записать mesh-preview, включая частично собранную модель
    #[arg(long, value_name = "PATH")]
    v2_mesh_preview_json: Option<String>,

    /// Сохранить даже вырожденные отверстия вместо геотехнического упрощения
    #[arg(long)]
    v2_preserve_details: bool,

    /// v2, геотехнический режим: допуск сведения стены на ось нижней несущей
    /// стены (смещение осей стен по этажам), в единицах модели. 0 — отключить.
    /// Предел перемещения вершин при замыкании поднимается до этого значения.
    #[arg(long, default_value_t = 0.05)]
    v2_stack_offset: f64,

    /// v2, геотехнический режим: допуск притяжки торца стены к оси другой
    /// стены и порог удаления лишних коллинеарных вершин у коротких рёбер,
    /// в единицах модели. 0 — отключить.
    #[arg(long, default_value_t = 0.05)]
    v2_wall_end_snap: f64,

    /// v2, геотехнический режим: максимальная ширина обрезаемой консоли за
    /// линией стыка, в единицах модели. 0 — не обрезать.
    #[arg(long, default_value_t = 0.25)]
    v2_console_width: f64,

    /// v2, геотехнический режим: максимальная ширина трещины конвертированной
    /// сетки внутри одной плоской конструкции (несвязанные узлы границы по
    /// разные стороны пустоты). Контур конструкции пересобирается без
    /// трещины; исходная сетка не сшивается. В единицах модели, 0 — отключить.
    #[arg(long, default_value_t = 0.01)]
    v2_crack_width: f64,

    /// v2, геотехнический режим: рёбра короче этого значения между двумя
    /// нужными углами (например, торцы нижней и верхней стен в нескольких
    /// миллиметрах друг от друга) схлопываются в одну вершину. По умолчанию
    /// 0,05 = 1/10 элемента PLAXIS 0,5 м. В единицах модели, 0 — отключить.
    #[arg(long, default_value_t = 0.05)]
    v2_edge_collapse: f64,

    /// v2, геотехнический режим (под PLAXIS): зазоры уже этого значения между
    /// вершиной одной конструкции и другой конструкцией закрываются (вершина
    /// переносится на плоскость и контур другой конструкции). По умолчанию
    /// 0,05 = 1/10 элемента PLAXIS 0,5 м. В единицах модели, 0 — сохранить
    /// все зазоры (например, настоящие деформационные швы).
    #[arg(long, default_value_t = 0.05)]
    v2_gap_closure: f64,
}

/// Geotechnical simplification tolerances of the v2 pipeline (model units).
struct V2Tolerances {
    stack_offset: f64,
    wall_end_snap: f64,
    console_width: f64,
    crack_width: f64,
    edge_collapse: f64,
    gap_closure: f64,
}

fn run_v2_preview(
    input: &str,
    output: &str,
    iterations: usize,
    include_mesh: bool,
    preserve_details: bool,
    tolerances: &V2Tolerances,
) -> Result<(), Box<dyn std::error::Error>> {
    use topo_reconstruct_rs::{
        parsers::LiraParser as V2LiraParser,
        reconstruction::{assembly, frame, graph, mesh, planes, recognize, reconcile},
    };

    if iterations == 0 {
        return Err("--v2-iterations must be positive".into());
    }
    for (name, value) in [
        ("--v2-stack-offset", tolerances.stack_offset),
        ("--v2-wall-end-snap", tolerances.wall_end_snap),
        ("--v2-console-width", tolerances.console_width),
        ("--v2-crack-width", tolerances.crack_width),
        ("--v2-edge-collapse", tolerances.edge_collapse),
        ("--v2-gap-closure", tolerances.gap_closure),
    ] {
        if !value.is_finite() || value < 0. {
            return Err(format!("{name} must be a finite non-negative length").into());
        }
    }
    // Stage timing on stderr when TOPO_TIMING is set (diagnostics only).
    let timing = std::env::var_os("TOPO_TIMING").is_some();
    let mut clock = Instant::now();
    let mut lap = |stage: &str| {
        if timing {
            eprintln!("[timing] {stage}: {:.1}s", clock.elapsed().as_secs_f64());
        }
        clock = Instant::now();
    };
    let mesh = V2LiraParser::parse(input)?;
    lap("parse");
    let axes = recognize::recognize(
        &mesh,
        &recognize::Policy {
            angle: 0.02,
            line_tolerance: 0.01,
            numerical_precision: 1e-8,
        },
    )?;
    lap("axis_recognition");
    let plane_report = planes::recognize(
        &mesh,
        &planes::Policy {
            angle: 0.02,
            distance: 0.01,
            precision: 1e-8,
        },
    )?;
    lap("plane_recognition");
    // A finite residual with no movement/axis failure is retried with a
    // doubled LSQR budget. The tolerance and all geometric budgets stay fixed.
    let result = frame::solve_with_retry(
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
        3,
    )?;
    lap("frame");
    let assembly_policy = assembly::Policy {
        closure_tolerance: 0.001,
        // A stacked wall moves by its offset when closing onto the lower axis.
        junction_movement_limit: tolerances.stack_offset.max(0.05),
        precision: 1e-7,
        minimum_edge: 0.001,
    };
    let topology = if preserve_details {
        assembly::assemble(&mesh, &result, &assembly_policy)?
    } else {
        assembly::assemble_geotechnical(
            &mesh,
            &result,
            &assembly_policy,
            &assembly::FeaturePolicy {
                maximum_console_width: tolerances.console_width,
                maximum_stack_offset: tolerances.stack_offset,
                maximum_wall_end_snap: tolerances.wall_end_snap,
                maximum_crack_width: tolerances.crack_width,
                maximum_collapsed_edge: tolerances.edge_collapse,
                maximum_gap: tolerances.gap_closure,
                ..Default::default()
            },
        )?
    };
    lap("assembly");
    let reconciliation =
        reconcile::solve(&mesh, &result, &topology, &reconcile::Policy::default())?;
    lap("reconciliation");
    let (mesh_report, mesh_error) = if include_mesh {
        match mesh::build_partial(
            &topology,
            &mesh::Policy {
                boundary_spacing: 0.5,
                maximum_area: 0.5,
                minimum_angle_degrees: 20.,
                maximum_added_vertices_per_surface: 10000,
            },
        ) {
            Ok(mesh) => (Some(mesh), None),
            Err(error) => (None, Some(error.to_string())),
        }
    } else {
        (None, None)
    };
    lap("mesh");
    let mut report = serde_json::json!({
        "constraint_graph": graph::Graph::from_frame(&result),
        "frame": result,
        "topology": topology,
        "reconciliation": reconciliation,
        "axis_recognition": axes,
        "plane_recognition": plane_report,
    });
    if include_mesh {
        report["mesh"] = serde_json::to_value(mesh_report)?;
        report["mesh_error"] = serde_json::to_value(mesh_error)?;
    }
    lap("report");
    // Buffered: the pretty report of a large model has hundreds of MB.
    use std::io::Write as _;
    if output == "-" {
        let mut out = std::io::BufWriter::new(std::io::stdout().lock());
        serde_json::to_writer_pretty(&mut out, &report)?;
        writeln!(out)?;
        out.flush()?;
    } else {
        let mut out = std::io::BufWriter::new(File::create(output)?);
        serde_json::to_writer_pretty(&mut out, &report)?;
        out.flush()?;
    }
    Ok(())
}

fn main() {
    let args = Args::parse();
    let tolerances = V2Tolerances {
        stack_offset: args.v2_stack_offset,
        wall_end_snap: args.v2_wall_end_snap,
        console_width: args.v2_console_width,
        crack_width: args.v2_crack_width,
        edge_collapse: args.v2_edge_collapse,
        gap_closure: args.v2_gap_closure,
    };
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
        if args.v2_mesh_preview_json.is_some() {
            eprintln!("Нельзя одновременно задавать --v2-preview-json и --v2-mesh-preview-json.");
            std::process::exit(2);
        }
        if let Err(error) = run_v2_preview(
            &args.input,
            path,
            args.v2_iterations,
            false,
            args.v2_preserve_details,
            &tolerances,
        ) {
            eprintln!("[V2 PREVIEW ERROR] {error}");
            std::process::exit(1);
        }
        eprintln!("[V2 PREVIEW] JSON сохранен: {path}");
        return;
    }
    if let Some(path) = &args.v2_mesh_preview_json {
        if let Err(error) = run_v2_preview(
            &args.input,
            path,
            args.v2_iterations,
            true,
            args.v2_preserve_details,
            &tolerances,
        ) {
            eprintln!("[V2 MESH PREVIEW ERROR] {error}");
            std::process::exit(1);
        }
        eprintln!("[V2 MESH PREVIEW] JSON сохранен: {path}");
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
