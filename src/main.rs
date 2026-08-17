#![allow(dead_code)]

mod config;
mod models;

use clap::Parser;
use config::ReconstructionConfig;
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
    println!("Допуск по высоте: {} м, Допуск упрощения: {} м", config.tol_dist, config.simplify_tol);

    let start_time = Instant::now();

    println!("Инициализация завершена за {:?}", start_time.elapsed());
}