//! Parse the loads of a LIRA file and print a summary.
fn main() -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let path = std::env::args().nth(1).ok_or("usage: loads_probe MODEL.txt")?;
    let bytes = std::fs::read(&path)?;
    let t = std::time::Instant::now();
    let set = topo_reconstruct_rs::parsers::loads::parse(&bytes);
    let mut by: std::collections::BTreeMap<(u16, u8), usize> = Default::default();
    for r in &set.rows {
        *by.entry((r.code, r.direction)).or_default() += 1;
    }
    println!("{} rows, {} cases, {} node axes, {} element axes, {:?}", set.rows.len(), set.cases.len(), set.node_axes.len(), set.element_axes.len(), t.elapsed());
    println!("cases: {:?}", set.cases.iter().take(4).collect::<Vec<_>>());
    println!("codes: {:?}", by);
    Ok(())
}
