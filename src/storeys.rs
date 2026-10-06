//! The stiffness of the storeys cut off (`cutoff`), carried by the floor
//! slab at the level of the cut. The walls and columns that stood on the
//! level have a bending stiffness about the horizontal axes (their plan
//! sections, the walls' own and by the parallel-axis terms). A cantilever's
//! stiffness against rotation of its base falls with its height, so the
//! removed part counts as one storey: `EI' = EI * h / H` (h the storey
//! height below the level, H the height removed). A plate of the plan size
//! of the cap with the same bending stiffness, `E t^3 B / 12 = EI'` for each
//! horizontal axis, gives the thickness; the geometric mean of the two is
//! used and the user's factor scales EI'. The cap keeps its own weight (the
//! weight of what was cut off is carried by the loads). Only the bending
//! stiffness of the cap grows (the plate factor PLKE), its membrane stiffness
//! and weight stay.
use crate::parsers::lira::Material;
use crate::reconstruction::assembly::cutoff::{floors, live_supports, Cut};
use crate::reconstruction::assembly::edit::State;
use hashbrown::HashMap;
use serde::Serialize;

/// First number of the stiffness types made for cap slabs.
pub const CAP_STIFFNESS_BASE: u32 = 2_000_000;

#[derive(Debug, Clone, Serialize)]
pub struct CapReport {
    /// Cap surfaces whose stiffness was replaced.
    pub surfaces: usize,
    /// Equivalent thickness (m) before the comparison with their own.
    pub equivalent_thickness: f64,
    /// Bending stiffness about X and Y (m4, in units of the cap's E).
    pub ix: f64,
    pub iy: f64,
    /// Storey height and height removed (m).
    pub storey: f64,
    pub removed_height: f64,
    pub notes: Vec<String>,
}

/// Surfaces that close the lower part: horizontal at the level of the cut.
fn cap_surfaces(state: &State, cut: &Cut) -> Vec<usize> {
    let model = &state.model;
    (0..model.surfaces().len())
        .filter(|&s| {
            model.planes()[model.surfaces()[s].plane].normal()[2].abs() > 0.98
                && model.surface_edges(s).flat_map(|e| model.edges()[e]).all(|v| (model.vertices()[v][2] - cut.z).abs() < 0.05)
        })
        .collect()
}

fn plate(materials: &HashMap<u32, Material>, id: u32) -> Option<(f64, f64, f64)> {
    match materials.get(&id)? {
        Material::Plate { e, thickness, density, .. } => Some((*e, *thickness, density.unwrap_or(0.))),
        _ => None,
    }
}

/// The state and materials with the cap slabs' stiffness replaced by an
/// equivalent one, or `None` when there is nothing to do (no cut, no cap, the
/// slab is already stiffer).
pub fn with_cap(state: &State, materials: &HashMap<u32, Material>, factor: f64) -> Option<(State, HashMap<u32, Material>, CapReport)> {
    let cut = state.cut.as_ref()?;
    let caps = cap_surfaces(state, cut);
    let first = *caps.first()?;
    let (e_cap, _, _) = plate(materials, state.stiffness[first])?;
    // Plan sections of what stood on the level: (area, centre, own Ix, own Iy), by the E of the cap.
    let mut parts: Vec<(f64, [f64; 2], f64, f64)> = vec![];
    let mut notes = vec![];
    let (walls, columns) = live_supports(state);
    for (wa, wb, stiffness) in walls {
        let Some((e, t, _)) = plate(materials, stiffness) else { continue };
        let n = e / e_cap;
        let (dx, dy) = (wb.x - wa.x, wb.y - wa.y);
        let l = dx.hypot(dy);
        if l < 1e-9 {
            continue;
        }
        let (c, s) = (dx / l, dy / l);
        parts.push((n * t * l, [(wa.x + wb.x) / 2., (wa.y + wb.y) / 2.], n * (t * l.powi(3) * s * s + l * t.powi(3) * c * c) / 12., n * (t * l.powi(3) * c * c + l * t.powi(3) * s * s) / 12.));
    }
    for (p, stiffness) in columns {
        // E (t/m2), area and own moments of the column (a rectangle: a square of that area).
        let (e, a, own) = match materials.get(&stiffness) {
            Some(Material::Bar { e, width, height, .. }) => (*e, width * height, (width * height).powi(2) / 12.),
            Some(Material::Section { .. }) => {
                // S1..S6 sections and piles: E, steel or concrete, by the rules of the MIDAS export.
                let Some(spec) = crate::midas_stiffness::bar_spec(stiffness, materials, &Default::default()) else { continue };
                let (Some(a), Some((ix, iy))) = (spec.area, spec.shape.inertia()) else { continue };
                (spec.young_kn / 9.80665, a, ix.min(iy))
            }
            _ => continue,
        };
        // The modular ratio scales the area and the own moment alike (EI is linear in E).
        let n = e / e_cap;
        parts.push((n * a, [p.x, p.y], n * own, n * own));
    }
    let area: f64 = parts.iter().map(|p| p.0).sum();
    if parts.is_empty() || area <= 0. {
        return None;
    }
    let centre = [parts.iter().map(|p| p.0 * p.1[0]).sum::<f64>() / area, parts.iter().map(|p| p.0 * p.1[1]).sum::<f64>() / area];
    let ix: f64 = parts.iter().map(|p| p.3 * 0. + p.2 + p.0 * (p.1[1] - centre[1]).powi(2)).sum();
    let iy: f64 = parts.iter().map(|p| p.3 + p.0 * (p.1[0] - centre[0]).powi(2)).sum();
    // Plan size of the cap.
    let model = &state.model;
    let (mut lo, mut hi) = ([f64::MAX; 2], [f64::MIN; 2]);
    for &s in &caps {
        for v in model.surface_edges(s).flat_map(|e| model.edges()[e]) {
            let p = model.vertices()[v];
            for k in 0..2 {
                lo[k] = lo[k].min(p[k]);
                hi[k] = hi[k].max(p[k]);
            }
        }
    }
    let (bx, by) = ((hi[0] - lo[0]).max(1.), (hi[1] - lo[1]).max(1.));
    // One storey of the removed part.
    let levels = floors(state);
    // A storey is between floors that are real floors (a good part of the cap's area), not landings.
    let below = levels.iter().rev().find(|f| f.z < cut.z - 0.15 && f.major).map(|f| cut.z - f.z);
    let storey = below.unwrap_or_else(|| {
        notes.push("высота этажа под уровнем отсечения не определена: принята 3 м".into());
        3.
    });
    let removed = (cut.top - cut.z).max(storey);
    let scale = (storey / removed) * factor.max(0.);
    let t_x = (12. * ix * scale / bx).cbrt();
    let t_y = (12. * iy * scale / by).cbrt();
    let t_eq = (t_x * t_y).sqrt();
    let mut new_state = state.clone();
    let mut new_materials = materials.clone();
    let mut changed = 0;
    for &s in &caps {
        let base = state.stiffness[s];
        let Some((_, t, _)) = plate(materials, base) else { continue };
        if t_eq <= t * 1.0001 {
            continue;
        }
        let id = CAP_STIFFNESS_BASE + base;
        if let Some(Material::Plate { e, nu, thickness, density, membrane, bending, shear }) = materials.get(&base).cloned() {
            // Only the bending stiffness grows (PLKE x (t'/t)^3): the membrane
            // stiffness and the weight of the slab stay as they were.
            let kb = bending.unwrap_or(1.) * (t_eq / t).powi(3);
            new_materials.insert(id, Material::Plate { e, nu, thickness, density, membrane, bending: Some(kb), shear });
            new_state.stiffness[s] = id;
            changed += 1;
        }
    }
    if changed == 0 {
        return None;
    }
    notes.push(format!("жёсткость отброшенной части: EI' = EI·h/H·k, h = {storey:.2} м, H = {removed:.2} м, k = {factor}"));
    Some((new_state, new_materials, CapReport { surfaces: changed, equivalent_thickness: t_eq, ix, iy, storey, removed_height: removed, notes }))
}
