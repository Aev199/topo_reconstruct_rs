#![allow(dead_code)]

mod config;
mod geometry;
mod models;
mod parsers;

use clap::Parser;
use config::ReconstructionConfig;
use parsers::LiraParser;
use std::time::Instant;

#[derive(Parser, Debug)]
#[command(author, version, about = "Реконструкция BIM-топологии из КЭ-сеток ПК ЛИРА (High-Performance Rust Core)")]
struct Args {
    /// Путь к текстовому файлу расчетной схемы ЛИРЫ (.txt)
    #[arg(default_value = "скала3.txt")]
    input: String,

    /// Путь для сохранения JSON отчета
    #[arg(short, long, default_value = "building_topology.json")]
    json: String,

    /// Путь для сохранения DXF файла
    #[arg(short, long, default_value = "building_topology.dxf")]
    dxf: String,
}

fn main() {
    let args = Args::parse();
    let config = ReconstructionConfig::default();

    println!("=== RECONSTRUCT TOPOLOGY CORE (RUST) ===");
    println!("Входной файл: {}", args.input);
    println!(
        "Допуск по высоте: {} м, Допуск упрощения: {} м",
        config.tol_dist, config.simplify_tol
    );

    // 1. Замер скорости параллельного парсинга
    println!("\n1. Чтение и парсинг расчетной схемы...");
    let parse_start = Instant::now();

    match LiraParser::parse(&args.input) {
        Ok(mesh_data) => {
            let parse_elapsed = parse_start.elapsed();
            println!(
                "   [OK] Схема успешно загружена за {:.2?}:",
                parse_elapsed
            );
            println!("   -> Узлов: {}", mesh_data.nodes.len());
            println!("   -> Конечных элементов: {}", mesh_data.elements.len());

            // 2. Тест канонизации узлов
            let canon_start = Instant::now();
            let canonical_map = geometry::canonicalize_nodes(&mesh_data.nodes, config.canonical_precision);
            println!(
                "   [OK] Канонизация узлов выполнена за {:.2?}: сшито {} узлов",
                canon_start.elapsed(),
                canonical_map.len()
            );
        }
        Err(e) => {
            eprintln!("   [ERROR] Ошибка при чтении файла: {}", e);
            std::process::exit(1);
        }
    }
}