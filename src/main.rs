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

    /// v2, геотехнический режим: не закрывать зазоры поперёк плоскости
    /// (например, верх стены ниже плиты). По умолчанию они закрываются
    /// (стена достраивается в своей плоскости); параллельные конструкции
    /// (плиты на разных уровнях) не сводятся никогда.
    #[arg(long)]
    v2_keep_gap_offsets: bool,

    /// v2, геотехнический режим: верхняя граница допуска упрощения
    /// поверхности. Допуск поверхности — половина её толщины (из жёсткости
    /// GEI), но не меньше допуска зазоров и не больше этого значения.
    /// В единицах модели.
    #[arg(long, default_value_t = 0.2)]
    v2_simplification_cap: f64,

    /// v2, геотехнический режим: отверстия уже этого значения (меньшая
    /// сторона описанного прямоугольника) заделываются, если к ним ничего
    /// не примыкает и через них ничего не проходит. В единицах модели,
    /// 0 — сохранить все отверстия.
    #[arg(long, default_value_t = 1.0)]
    v2_min_opening: f64,

    /// v2, для разработки: файл кэша решённого каркаса. Если файл есть и
    /// записан для того же входного файла и тех же параметров каркаса,
    /// каркас берётся из него (сборка и сетка пересчитываются), иначе
    /// каркас решается и записывается. После изменения кода каркаса кэш
    /// нужно удалить.
    #[arg(long, value_name = "PATH")]
    v2_frame_cache: Option<String>,
}

/// Geotechnical simplification tolerances of the v2 pipeline (model units).
struct V2Tolerances {
    stack_offset: f64,
    wall_end_snap: f64,
    console_width: f64,
    crack_width: f64,
    edge_collapse: f64,
    gap_closure: f64,
    gap_offsets: bool,
    simplification_cap: f64,
    min_opening: f64,
    /// Development cache of the solved frame (not a tolerance).
    frame_cache: Option<String>,
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
        ("--v2-simplification-cap", tolerances.simplification_cap),
        ("--v2-min-opening", tolerances.min_opening),
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
    let mut mesh = V2LiraParser::parse(input)?;
    // Ties between structures in the analysis model (before rigid links are
    // dropped): structures tied together are joined, not kept apart.
    let connections = if preserve_details {
        vec![]
    } else {
        topo_reconstruct_rs::input::node_links(&mesh, &V2LiraParser::parse_rigid_bodies(input)?)
    };
    // Rigid links of the analysis model (fans spreading a column into a
    // slab, offsets) are not structures: geotechnical mode leaves them out.
    let rigid_links = if preserve_details {
        vec![]
    } else {
        let links =
            topo_reconstruct_rs::input::rigid_links(&mesh, &V2LiraParser::parse_stiffness(input)?);
        let removed: std::collections::HashSet<u32> = links.iter().copied().collect();
        mesh.elements.retain(|e| !removed.contains(&e.id));
        links
    };
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
    // Geotechnical mode also closes source gaps across planes (a wall end a
    // few centimetres from the walls it abuts) by moving planes in the solve.
    let gap_tolerance = if preserve_details || !tolerances.gap_offsets {
        0.
    } else {
        tolerances.gap_closure
    };
    let frame_policy = |over_constrained_panels: bool| frame::Policy {
        up: [0., 0., 1.],
        angle: 0.02,
        maximum_movement: 0.15,
        relative_movement: 0.05,
        minimum_length: 0.03,
        residual_tolerance: 1e-7,
        iterations,
        // Nearly parallel walls through common nodes may be joined into
        // one panel within the tolerance for aligning walls in one line.
        panel_tolerance: if preserve_details {
            0.
        } else {
            tolerances.stack_offset
        },
        // A short bar lying in a wall or slab follows its plane.
        geotechnical: !preserve_details,
        over_constrained_panels,
    };
    let solve = |policy: &frame::Policy| {
        frame::solve_closing_gaps(
            &mesh,
            &axes,
            &plane_report,
            policy,
            // Doubling iteration budgets: 1, 2, 4 and 8 times the base (a
            // curved-wall model needed about 4800 iterations).
            4,
            gap_tolerance,
            0.001,
        )
    };
    // Relaxation ladder: a frame the default rules cannot satisfy is solved
    // again with panels merged at over-constrained nodes (geotechnical
    // only). The step used is reported.
    let solve_ladder = || -> Result<(frame::Report, usize), Box<dyn std::error::Error>> {
        let result = solve(&frame_policy(false))?;
        if !result.accepted && !preserve_details {
            let relaxed = solve(&frame_policy(true))?;
            if relaxed.accepted {
                return Ok((relaxed, 1));
            }
        }
        Ok((result, 0))
    };
    let (result, relaxation) = match &tolerances.frame_cache {
        None => solve_ladder()?,
        Some(path) => {
            // The key covers what the frame depends on besides the code: the
            // input text, the policies and the gap tolerance.
            use std::hash::{Hash, Hasher};
            let mut hasher = std::collections::hash_map::DefaultHasher::new();
            std::fs::read(input)?.hash(&mut hasher);
            serde_json::to_string(&frame_policy(false))?.hash(&mut hasher);
            gap_tolerance.to_bits().hash(&mut hasher);
            preserve_details.hash(&mut hasher);
            // Only when links were removed: other caches stay valid.
            if !rigid_links.is_empty() {
                rigid_links.hash(&mut hasher);
            }
            let key = format!("{:016x}", hasher.finish());
            let cached = std::fs::read(path).ok().and_then(|bytes| {
                let value: serde_json::Value = serde_json::from_slice(&bytes).ok()?;
                (value["key"] == key.as_str()).then_some(())?;
                let frame = serde_json::from_value(value["frame"].clone()).ok()?;
                Some((frame, value["relaxation"].as_u64()? as usize))
            });
            match cached {
                Some(hit) => {
                    eprintln!("[V2] каркас взят из кэша {path}");
                    hit
                }
                None => {
                    let (frame, relaxation) = solve_ladder()?;
                    let value = serde_json::json!({
                        "key": key,
                        "relaxation": relaxation,
                        "frame": frame,
                    });
                    std::fs::write(path, serde_json::to_vec(&value)?)?;
                    (frame, relaxation)
                }
            }
        }
    };
    lap("frame");
    let assembly_policy = assembly::Policy {
        closure_tolerance: 0.001,
        // A stacked wall moves by its offset when closing onto the lower axis.
        // A vertex may move by the largest rule tolerance plus the closure
        // tolerance (a node within it counts as on a plane): an aligned wall
        // is offset by at most the tolerance at its junction, and its slight
        // non-parallelism elsewhere is measured against the closure tolerance.
        junction_movement_limit: tolerances.stack_offset.max(0.05) + 0.001,
        precision: assembly::PRECISION,
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
                close_offset_gaps: tolerances.gap_offsets,
                surface_thickness: V2LiraParser::parse_sections(input)?
                    .into_iter()
                    .filter_map(|(id, section)| match section {
                        topo_reconstruct_rs::parsers::Section::Plate { thickness } => {
                            Some((id, thickness))
                        }
                        _ => None,
                    })
                    .collect(),
                maximum_simplification: tolerances.simplification_cap,
                minimum_opening_width: tolerances.min_opening,
                connections,
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
        "relaxation": { "frame": relaxation },
        "rigid_links": { "removed": rigid_links.len(), "elements": rigid_links },
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
        gap_offsets: !args.v2_keep_gap_offsets,
        simplification_cap: args.v2_simplification_cap,
        min_opening: args.v2_min_opening,
        frame_cache: args.v2_frame_cache.clone(),
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
