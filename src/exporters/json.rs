use crate::models::ReconstructionReport;
use std::fs::File;
use std::io::{self, BufWriter};
use std::path::Path;

pub struct JsonExporter;

impl JsonExporter {
    pub fn export<P: AsRef<Path>>(report: &ReconstructionReport, filepath: P) -> io::Result<()> {
        let file = File::create(filepath)?;
        let writer = BufWriter::new(file);
        serde_json::to_writer_pretty(writer, report)?;
        Ok(())
    }
}