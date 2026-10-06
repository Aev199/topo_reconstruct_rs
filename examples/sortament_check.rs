//! Compare the simplified shape (area, moments of inertia) of every catalogue profile with the catalogue's own
//! A, Ix, Iy: where a profile gets the wrong shape (a round tube as a square one) the ratio is far from 1.
use std::collections::BTreeMap;
use topo_reconstruct_rs::sections::{profile_shape, Profile};
use topo_reconstruct_rs::sortament;

fn main() {
    let mut worst: BTreeMap<String, (f64, String)> = BTreeMap::new();
    let mut counts: BTreeMap<String, usize> = BTreeMap::new();
    for e in sortament::entries() {
        // The block-13 section and designation that name this entry.
        let section = if e.kind == "Труба круглая" { "Pipe" } else if e.kind.contains("квадрат") || e.kind.contains("прямоугол") { "Tubing" } else if sortament::is_i_kind(&e.kind) { "DoubleT" } else { continue };
        let p = Profile { section: section.into(), shape: e.name.clone(), builtup: false };
        let Some(shape) = profile_shape(&p) else { continue };
        let (Some(a), Some(ix), Some(iy)) = (e.area, e.ix, e.iy) else { continue };
        let (sa, (sx, sy)) = (shape.area().unwrap_or(0.) * 1e4, shape.inertia().map_or((0., 0.), |(x, y)| (x * 1e8, y * 1e8)));
        for (what, ratio) in [("area", sa / a), ("ix", sx / ix), ("iy", sy / iy.max(1e-9))] {
            let key = format!("{} {}", e.kind, what);
            *counts.entry(key.clone()).or_default() += 1;
            let dev = (ratio - 1.).abs();
            if worst.get(&key).is_none_or(|w| dev > (w.0 - 1.).abs()) {
                worst.insert(key, (ratio, e.name.clone()));
            }
        }
    }
    for (k, (ratio, name)) in worst {
        println!("{k:55} n={:5} worst ratio {ratio:.3} at {name}", counts[&k]);
    }
}
