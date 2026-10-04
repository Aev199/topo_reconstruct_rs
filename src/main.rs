use clap::Parser;
use std::{fs::File, path::PathBuf, time::Instant};
use topo_reconstruct_rs::pipeline::{self, Options, Profile};

/// Восстановление геотехнической геометрии (поверхности и оси) из
/// текстовой расчётной модели ПК ЛИРА. Допуски по умолчанию — профиль
/// PLAXIS (элемент 0,5 м, особенности меньше h/10 закрываются). Все длины —
/// в единицах модели (метры).
#[derive(Parser, Debug)]
#[command(author, version, about)]
struct Args {
    /// Текстовый файл расчётной схемы ЛИРА (.txt).
    input: PathBuf,

    /// Файл JSON с отчётом всех этапов ("-" — stdout).
    #[arg(short, long)]
    output: String,

    /// Построить пробную сетку (проверка пригодности к перебивке).
    #[arg(long)]
    mesh: bool,

    /// Файл кэша решённого каркаса: если записан для того же входного файла
    /// и тех же параметров каркаса, каркас берётся из него. После изменения
    /// кода каркаса кэш нужно удалить.
    #[arg(long, value_name = "PATH")]
    frame_cache: Option<PathBuf>,

    /// Целевой размер элемента сетки.
    #[arg(long, default_value_t = Profile::plaxis().element_size)]
    element_size: f64,

    /// Сведение стены на ось нижней несущей стены (и стен в одну линию).
    #[arg(long, default_value_t = Profile::plaxis().stack_offset)]
    stack_offset: f64,

    /// Притяжка торца стены к оси другой стены; порог лишних вершин у коротких рёбер.
    #[arg(long, default_value_t = Profile::plaxis().wall_end_snap)]
    wall_end_snap: f64,

    /// Максимальная ширина обрезаемой консоли за линией стыка.
    #[arg(long, default_value_t = Profile::plaxis().console_width)]
    console_width: f64,

    /// Максимальная ширина трещины конвертированной сетки внутри конструкции.
    #[arg(long, default_value_t = Profile::plaxis().crack_width)]
    crack_width: f64,

    /// Рёбра и стержни короче этого значения схлопываются.
    #[arg(long, default_value_t = Profile::plaxis().edge_collapse)]
    edge_collapse: f64,

    /// Зазоры уже этого значения между конструкциями закрываются (0 — все сохранить).
    #[arg(long, default_value_t = Profile::plaxis().gap_closure)]
    gap_closure: f64,

    /// Не закрывать зазоры поперёк плоскости (верх стены ниже плиты).
    #[arg(long)]
    keep_gap_offsets: bool,

    /// Верхняя граница допуска упрощения (половина толщины пластины).
    #[arg(long, default_value_t = Profile::plaxis().simplification_cap)]
    simplification_cap: f64,

    /// Свободные отверстия уже этого значения заделываются (0 — все сохранить).
    #[arg(long, default_value_t = Profile::plaxis().min_opening)]
    min_opening: f64,

    /// Отверстия длиннее этого значения сохраняются при любой ширине.
    #[arg(long, default_value_t = Profile::plaxis().max_opening_length)]
    max_opening_length: f64,

    /// Базовое число итераций решателя каркаса.
    #[arg(long, default_value_t = Profile::plaxis().iterations)]
    iterations: usize,
}

fn main() {
    let args = Args::parse();
    let profile = Profile {
        element_size: args.element_size,
        stack_offset: args.stack_offset,
        wall_end_snap: args.wall_end_snap,
        console_width: args.console_width,
        crack_width: args.crack_width,
        edge_collapse: args.edge_collapse,
        gap_closure: args.gap_closure,
        close_offset_gaps: !args.keep_gap_offsets,
        simplification_cap: args.simplification_cap,
        min_opening: args.min_opening,
        max_opening_length: args.max_opening_length,
        iterations: args.iterations,
    };
    let options = Options {
        mesh: args.mesh,
        frame_cache: args.frame_cache.clone(),
    };
    // Stage timing on stderr when TOPO_TIMING is set (diagnostics only).
    let timing = std::env::var_os("TOPO_TIMING").is_some();
    let mut clock = Instant::now();
    let mut stage = |name: &str| {
        if timing {
            eprintln!("[timing] {name}: {:.1}s", clock.elapsed().as_secs_f64());
        }
        clock = Instant::now();
    };
    let result = pipeline::run(&args.input, &profile, &options, &mut stage).and_then(|output| {
        // Buffered: the pretty report of a large model has hundreds of MB.
        let write = |out: &mut dyn std::io::Write| -> Result<(), pipeline::Error> {
            serde_json::to_writer_pretty(&mut *out, &output)?;
            writeln!(out)?;
            out.flush()?;
            Ok(())
        };
        if args.output == "-" {
            write(&mut std::io::BufWriter::new(std::io::stdout().lock()))
        } else {
            write(&mut std::io::BufWriter::new(File::create(&args.output)?))
        }
    });
    match result {
        Ok(()) => eprintln!("JSON сохранён: {}", args.output),
        Err(error) => {
            eprintln!("Ошибка: {error}");
            std::process::exit(1);
        }
    }
}
