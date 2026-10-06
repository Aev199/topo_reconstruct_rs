//! Materials, sections and thicknesses of a MIDAS Civil `.mxt` file from the
//! stiffness types of a LIRA model, by the rules of the Lira_Midas-converter:
//!
//! * a material per distinct (E, nu, unit weight, plate or bar): types that
//!   agree share it; named `Plate_Beton_p12_h0.2`, `Beam_Steel_p3`,
//!   `Beam_Material_p1` after the unit weight (7850 kg/m3 +-20 % steel,
//!   2500 kg/m3 +-20 % concrete) and the lowest LIRA type of the group;
//! * Poisson's ratio from the file (`Mu`, GEI), else 0.3 for steel, 0.2 for
//!   concrete, 0.3 otherwise;
//! * steel sections (S1, S2, S3, S5, S6, S8, profiles of the block 13) get
//!   E = 206 GPa, unit weight 76.98 kN/m3 and nu = 0.3, unless the file's E
//!   is a concrete one (1e6 .. 6e6 t/m2: the converter makes every S6 steel,
//!   which turns round concrete columns into steel);
//! * bar sections are `DBUSER` shapes (SB, P, SR, H, T, B) of the dimensions,
//!   `VALUE` blocks (A, Iyy, Izz) for types with only EF, EIy, EIz; equal
//!   shapes share one section; rigid rods (E > 1e8 t/m2, 1 x 1 m) get `_AGT`
//!   in the name;
//! * a thickness per distinct plate thickness (to 1 mm), `Plate_0.2`.
//!
//! Differences from the converter, which reads some rows wrongly: the S0
//! dimensions and E are taken from the S0 row (the converter takes EIy of the
//! numeric row as the width when the type has one), the thickness of a plate
//! above 1 m is not read as centimetres, and E of a type with an EF row is not
//! taken as EF.
use crate::parsers::lira::Material;
use crate::sections::{profile_properties, profile_shape, Profile, Shape};
use hashbrown::HashMap;
use std::collections::{BTreeMap, BTreeSet};

/// Tonne-force to kN.
const G: f64 = 9.80665;
/// Steel: E (kN/m2), unit weight (kN/m3), nu.
const STEEL_E: f64 = 206_000_000.;
const STEEL_GAMMA: f64 = 7.85 * G;
const STEEL_NU: f64 = 0.3;
const CONCRETE_GAMMA: f64 = 2.5 * G;
const CONCRETE_NU: f64 = 0.2;
const DEFAULT_NU: f64 = 0.3;
/// Where the Young's modulus (t/m2) of a type is a concrete one.
const CONCRETE_E: std::ops::RangeInclusive<f64> = 1.0e6..=6.0e6;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mode {
    /// Nominal sections and thicknesses, the converter's file.
    Converter,
    /// The stiffness LIRA analysed with: EF, EIy, EIz of a type as a `VALUE` section,
    /// WLKE/PLKE of a plate as an equivalent thickness and E.
    Lira,
}

#[derive(Debug, Clone, Copy)]
pub struct Options {
    pub mode: Mode,
    /// Unit weights are multiplied by this (1: as LIRA; 0: the program adds no self-weight).
    pub density_multiplier: f64,
}

impl Default for Options {
    fn default() -> Self {
        Options { mode: Mode::Converter, density_multiplier: 1. }
    }
}

/// The sections of the file and which element uses which.
#[derive(Debug, Clone, Default)]
pub struct Stiffness {
    pub material_lines: Vec<String>,
    pub section_lines: Vec<String>,
    pub thickness_lines: Vec<String>,
    /// (is a plate, LIRA type) -> material number.
    pub material_of: BTreeMap<(bool, u32), usize>,
    /// LIRA type of a bar -> section number.
    pub section_of: BTreeMap<u32, usize>,
    /// LIRA type of a plate -> thickness number.
    pub thickness_of: BTreeMap<u32, usize>,
    /// What was changed or could not be transferred.
    pub notes: Vec<String>,
    /// LIRA types of the used elements without a material or section.
    pub missing: BTreeSet<u32>,
    pub sections: usize,
}

fn class(gamma: f64) -> &'static str {
    // The epsilon of the converter removes rounding differences.
    let eps = 0.02;
    if (gamma - STEEL_GAMMA).abs() <= STEEL_GAMMA * 0.2 + eps {
        "Steel"
    } else if (gamma - CONCRETE_GAMMA).abs() <= CONCRETE_GAMMA * 0.2 + eps {
        "Beton"
    } else {
        "Material"
    }
}

fn trim(v: f64, digits: usize) -> String {
    format!("{v:.digits$}").trim_end_matches('0').trim_end_matches('.').to_string()
}

/// Torsion constant of a rectangle b x h (m4).
fn rectangle_torsion(b: f64, h: f64) -> f64 {
    let (b, h) = (b.min(h), b.max(h));
    h * b.powi(3) * (1. / 3. - 0.21 * (b / h) * (1. - b.powi(4) / (12. * h.powi(4))))
}

struct MaterialSpec {
    plate: bool,
    young: f64,
    nu: f64,
    gamma: f64,
    /// Thickness of a plate, for the name.
    thickness: f64,
    /// A rigid rod: a weightless material of its own.
    rigid: bool,
}

/// The section of one bar type: the line without the number and a name.
struct SectionSpec {
    name: String,
    /// `DBUSER` (one body line) or `VALUE` (the body line and three of values).
    dbuser: bool,
    /// The lines after `YES, ` (the first) and the data lines of a VALUE block.
    lines: Vec<String>,
}

/// What a stiffness type of a bar is made of: the shape, E, nu and unit weight.
#[derive(Debug, Clone)]
pub struct BarSpec {
    pub shape: Shape,
    /// The numeric row EF, EIy, EIz, GIk of the type (tf, tf m2), when it has one.
    pub numeric: Option<[f64; 4]>,
    /// E (kN/m2), nu, unit weight (kN/m3; 0 for a rigid rod).
    pub young_kn: f64,
    pub nu: f64,
    pub gamma_kn: f64,
    /// Area (m2) of the section: the catalogue's for a profile, else of the shape.
    pub area: Option<f64>,
    pub steel: bool,
    pub rigid: bool,
    /// E was taken from EF / A (the E of the section row disagrees with the row EF).
    pub from_row: bool,
    /// The type has no E at all (only EF, EIy, EIz, GIk and no dimensions): steel is assumed.
    pub assumed_e: bool,
    /// A profile of the block 13 that is not in the catalogue.
    pub from_designation: bool,
}

/// The shape, E, nu and unit weight of a bar type by the converter's rules (see the module).
pub fn bar_spec(id: u32, materials: &HashMap<u32, Material>, profiles: &HashMap<u32, Profile>) -> Option<BarSpec> {
    let profile = profiles.get(&id);
    let material = materials.get(&id);
    // The shape, E (t/m2), nu, linear weight RO (t/m) and the numeric EF, EIy, EIz, GIk of the type.
    let (mut shape, file_e, file_nu, ro, numeric) = match material {
        Some(&Material::Bar { e, nu, width, height, density, stiffness }) => (Shape::Rect { h: height, b: width }, Some(e), nu, density, stiffness),
        Some(Material::Section { e, nu, density, shape, stiffness, .. }) => (shape.clone(), *e, *nu, *density, *stiffness),
        _ if profile.is_some() => (Shape::Explicit, None, None, None, None),
        _ => return None,
    };
    let mut from_designation = false;
    if let Some(p) = profile {
        match profile_shape(p) {
            Some(found) => shape = found,
            None => from_designation = true,
        }
    }
    let area = match profile.and_then(profile_properties) {
        Some(p) => Some(p.area),
        None => shape.area(),
    };
    // E of the type (t/m2): LIRA analyses with the numeric row EF, EIy, EIz, GIk when there is one,
    // so the type's E is EF / A of its section; the E written on the S line is kept when it agrees
    // (S0 rows) and replaced when it does not (the piles of S6 rows with E = 720000 whose EF belongs
    // to E of 4e6).
    let (mut e_type, mut from_row) = (file_e, false);
    if let (Some([ef, ..]), Some(a)) = (numeric, area) {
        if a > 0. && ef > 0. {
            let by_row = ef / a;
            if file_e.is_none_or(|e| (by_row / e - 1.).abs() > 0.02) {
                from_row = file_e.is_some();
                e_type = Some(by_row);
            }
        }
    }
    // Steel by the kind of the section (S1, S2, S3, S5, S6, S8) or a catalogue profile, unless E of
    // the type is a concrete one; a rectangular S0 bar is never steel by its kind. A type without E
    // (only EF, EIy, EIz, GIk and no dimensions) is taken as a steel one.
    let steel_kind = match material {
        Some(Material::Bar { .. }) => false,
        Some(Material::Section { section, .. }) => section != "explicit",
        _ => true,
    };
    let steel = (steel_kind || profile.is_some() || e_type.is_none()) && !e_type.is_some_and(|e| CONCRETE_E.contains(&e));
    // Rigid rods (АЖТ): E > 1e8 t/m2 and 1 x 1 m; a weightless material of their own.
    let rigid = matches!(material, Some(Material::Bar { e, width, height, .. }) if *e > 1e8 && (width - 1.).abs() < 1e-12 && (height - 1.).abs() < 1e-12);
    let (young_kn, nu, gamma_kn);
    if steel {
        (young_kn, nu, gamma_kn) = (STEEL_E, STEEL_NU, STEEL_GAMMA);
    } else {
        young_kn = e_type.map_or(STEEL_E, |e| e * G);
        gamma_kn = match (ro, area) {
            (Some(ro), Some(a)) if a > 0. => ro * G / a,
            _ => CONCRETE_GAMMA,
        };
        nu = file_nu.unwrap_or_else(|| match class(gamma_kn) {
            "Steel" => STEEL_NU,
            "Beton" => CONCRETE_NU,
            _ => DEFAULT_NU,
        });
    }
    Some(BarSpec { shape, numeric, young_kn, nu, gamma_kn: if rigid { 0. } else { gamma_kn }, area, steel, rigid, from_row, assumed_e: e_type.is_none(), from_designation })
}

/// Builds the materials, sections and thicknesses for the stiffness types
/// of the plates and of the bars of the mesh.
pub fn build(
    materials: &HashMap<u32, Material>,
    profiles: &HashMap<u32, Profile>,
    plates: &BTreeSet<u32>,
    bars: &BTreeSet<u32>,
    options: Options,
) -> Stiffness {
    let mut out = Stiffness::default();
    let multiplier = options.density_multiplier;
    let effective_bars = options.mode == Mode::Lira;
    let mut specs: BTreeMap<(bool, u32), MaterialSpec> = BTreeMap::new();
    let mut section_specs: BTreeMap<u32, SectionSpec> = BTreeMap::new();
    let mut plate_thickness: BTreeMap<u32, f64> = BTreeMap::new();
    let (mut steel_forced, mut reduced_ignored, mut from_designation, mut from_row_ef, mut assumed_e) = (vec![], vec![], vec![], vec![], vec![]);
    let mut equivalent: Vec<u32> = vec![];

    for &id in plates {
        let Some(&Material::Plate { e, nu, thickness, density, membrane, bending, .. }) = materials.get(&id) else {
            out.missing.insert(id);
            continue;
        };
        let (km, kb) = (membrane.unwrap_or(1.), bending.unwrap_or(1.));
        let mut gamma = density.unwrap_or(2.5) * G;
        let (mut young, mut t) = (e * G, thickness);
        let reduced = (km - 1.).abs() > 1e-12 || (kb - 1.).abs() > 1e-12;
        // The cap slab of a cut tower is the reason for a bending factor: always applied.
        let cap = id >= crate::storeys::CAP_STIFFNESS_BASE;
        if reduced && km > 0. && kb > 0. {
            if options.mode == Mode::Lira || cap {
                // E'd' = E d km, E'd'^3 = E d^3 kb and gamma' d' = gamma d keep membrane, bending and weight.
                let d2 = thickness * (kb / km).sqrt();
                young = young * km * thickness / d2;
                gamma = gamma * thickness / d2;
                t = d2;
                equivalent.push(id);
            } else {
                reduced_ignored.push(id);
            }
        }
        specs.insert((true, id), MaterialSpec { plate: true, young, nu, gamma: (gamma * multiplier).max(0.01), thickness: t, rigid: false });
        plate_thickness.insert(id, t);
    }

    for &id in bars {
        let Some(spec) = bar_spec(id, materials, profiles) else {
            out.missing.insert(id);
            continue;
        };
        let profile = profiles.get(&id);
        let (shape, numeric, young_kn, nu, rigid) = (spec.shape.clone(), spec.numeric, spec.young_kn, spec.nu, spec.rigid);
        if spec.from_designation {
            from_designation.push(id);
        }
        if spec.from_row {
            from_row_ef.push(id);
        }
        if spec.assumed_e {
            assumed_e.push(id);
        }
        if spec.steel {
            steel_forced.push(id);
        }
        let gamma = if rigid { 0. } else { (spec.gamma_kn * multiplier).max(0.01) };
        specs.insert((false, id), MaterialSpec { plate: false, young: young_kn, nu, gamma, thickness: 0., rigid });

        // The section.
        let mut name = crate::midas::safe_name(&shape.name(id, profile), &format!("Section-{id}"));
        if rigid {
            name += "_AGT";
        }
        let g = young_kn / (2. * (1. + nu));
        let value_block = |a: f64, iy: f64, iz: f64, j: f64| -> Vec<String> {
            vec![
                "H, BUILT, 0, 0, 0, 0, 0, 0, 0, 0, 0".to_string(),
                format!("       {a:.8}, 0, 0, {j:.8}, {iy:.8}, {iz:.8}, 0, 0, 0, 0"),
                "       0, 0, 0, 0, 0, 0, 0, 0, 0, 0".into(),
                "       0, 0, 0, 0, 0, 0, 0, 0, 0, 0".into(),
            ]
        };
        let (dbuser, lines) = if let (true, Some([ef, eiy, eiz, gik])) = (effective_bars, numeric) {
            // The stiffness LIRA analysed with: A = EF/E, I = EI/E (kN, m).
            let (ef, eiy, eiz, gik) = (ef * G, eiy * G, eiz * G, gik * G);
            let j = if gik > 0. { gik / g } else { match shape { Shape::Rect { h, b } => rectangle_torsion(b, h), _ => 0. } };
            (false, value_block(ef / young_kn, (eiy / young_kn).max(1e-12), (eiz / young_kn).max(1e-12), j))
        } else if let Some(body) = shape.dbuser_body() {
            if let (Some([ef, eiy, eiz, _]), Shape::Rect { h, b }) = (numeric, &shape) {
                // Nominal geometry: a reduction of the type in LIRA is not carried.
                let (a, e_t) = (b * h, young_kn / G);
                let nominal = [e_t * a, e_t * b * h.powi(3) / 12., e_t * h * b.powi(3) / 12.];
                if (ef / nominal[0] - 1.).abs() > 0.02 || (eiy / nominal[1] - 1.).abs() > 0.02 || (eiz / nominal[2] - 1.).abs() > 0.02 {
                    reduced_ignored.push(id);
                }
            }
            (true, vec![body])
        } else if let Some([ef, eiy, eiz, gik]) = numeric {
            // Only EF, EIy, EIz, GIk and no dimensions (a rolled profile not found in the catalogue).
            let (ef, eiy, eiz, gik) = (ef * G, eiy * G, eiz * G, gik * G);
            let j = if gik > 0. { gik / g } else { 0. };
            (false, value_block(ef / young_kn, (eiy / young_kn).max(1e-12), (eiz / young_kn).max(1e-12), j))
        } else {
            out.missing.insert(id);
            continue;
        };
        section_specs.insert(id, SectionSpec { name, dbuser, lines });
    }

    // Materials: one per signature, the number follows the lowest LIRA type of the group.
    let mut groups: BTreeMap<(i64, i64, i64, bool, bool), Vec<(bool, u32)>> = BTreeMap::new();
    for (&key, spec) in &specs {
        let sig = ((spec.young * 100.).round() as i64, (spec.nu * 1e6).round() as i64, (spec.gamma * 1e4).round() as i64, spec.plate, spec.rigid);
        groups.entry(sig).or_default().push(key);
    }
    let mut ordered: Vec<(u32, Vec<(bool, u32)>)> = groups.into_values().map(|g| (g.iter().map(|k| k.1).min().unwrap_or(0), g)).collect();
    ordered.sort_by_key(|g| (g.0, !g.1[0].0));
    let mut names: BTreeSet<String> = BTreeSet::new();
    for (number, (first, group)) in ordered.iter().enumerate() {
        let id = number + 1;
        let spec = &specs[&group[0]];
        let mut name = if spec.rigid { format!("AGT_p{first}") } else { format!("{}_{}_p{first}", if spec.plate { "Plate" } else { "Beam" }, class(spec.gamma)) };
        if spec.plate {
            // The thickness in the name when the group has one thickness.
            let ts: BTreeSet<i64> = group.iter().map(|k| (specs[k].thickness * 1e4).round() as i64).collect();
            if ts.len() == 1 {
                name += &format!("_h{}", trim(spec.thickness, 3));
            }
        }
        let base = name.clone();
        let mut k = 1;
        while !names.insert(name.to_lowercase()) {
            k += 1;
            name = format!("{base}_{k}");
        }
        out.material_lines.push(format!("{id:>5}, USER, {name}, 0, 0, , C, NO, 0, 2, {:.2}, {}, 0, {:.2}, 0", spec.young, spec.nu, spec.gamma));
        for &key in group {
            out.material_of.insert(key, id);
        }
    }

    // Sections: equal geometry, one section.
    let mut by_geometry: BTreeMap<(bool, Vec<String>), usize> = BTreeMap::new();
    let mut section_names: BTreeSet<String> = BTreeSet::new();
    for (&id, spec) in &section_specs {
        let number = match by_geometry.get(&(spec.dbuser, spec.lines.clone())) {
            Some(&n) => n,
            None => {
                let n = by_geometry.len() + 1;
                by_geometry.insert((spec.dbuser, spec.lines.clone()), n);
                let mut name = spec.name.clone();
                let mut k = 1;
                while !section_names.insert(name.to_lowercase()) {
                    k += 1;
                    name = format!("{}_{k}", spec.name);
                }
                if spec.dbuser {
                    out.section_lines.push(format!("{n:>5}, DBUSER, {name}, CC, 0, 0, 0, 0, 0, 0, YES, {}", spec.lines[0]));
                } else {
                    out.section_lines.push(format!("{n:>5}, VALUE, {name}, CC, 0, 0, 0, 0, 0, 0, YES, {}", spec.lines[0]));
                    out.section_lines.extend(spec.lines[1..].iter().cloned());
                }
                n
            }
        };
        out.section_of.insert(id, number);
    }
    out.sections = by_geometry.len();

    // Thicknesses continue the numbers of the sections.
    let mut by_thickness: BTreeMap<i64, usize> = BTreeMap::new();
    for (&id, &t) in &plate_thickness {
        let key = (t * 1000.).round() as i64;
        let number = *by_thickness.entry(key).or_insert_with(|| {
            let n = out.sections + out.thickness_lines.len() + 1;
            out.thickness_lines.push(format!("{n:>5}, VALUE, YES, {t:.3}, 0, NO, 0, 0, Plate_{}", trim(t, 3)));
            n
        });
        out.thickness_of.insert(id, number);
    }

    let list = |v: &[u32]| v.iter().map(u32::to_string).collect::<Vec<_>>().join(", ");
    if !steel_forced.is_empty() {
        out.notes.push(format!("сталь (E = 206 ГПа, 76,98 кН/м3, ν = 0,3 по правилам конвертера) у типов жёсткости: {}", list(&steel_forced)));
    }
    if !reduced_ignored.is_empty() {
        out.notes.push(format!(
            "приведённая жёсткость LIRA (EI, WLKE/PLKE) не применена, номинальные сечения и толщины, типы: {}",
            list(&reduced_ignored)
        ));
    }
    if !equivalent.is_empty() {
        out.notes.push(format!("плиты с WLKE/PLKE заменены эквивалентной толщиной и E (мембрана, изгиб и вес сохранены), типы: {}", list(&equivalent)));
    }
    if !from_row_ef.is_empty() {
        out.notes.push(format!("E типа взят как EF / A по строке жёсткости (E в строке сечения не согласуется с EF), типы: {}", list(&from_row_ef)));
    }
    if !assumed_e.is_empty() {
        out.notes.push(format!("у типа нет E и размеров (только EF, EIy, EIz, GIk): принят модуль стали, типы: {}", list(&assumed_e)));
    }
    if !from_designation.is_empty() {
        out.notes.push(format!("профиль не найден в сортаменте, размеры по обозначению или номинальные, типы: {}", list(&from_designation)));
    }
    if options.density_multiplier != 1. {
        out.notes.push(format!("удельные веса умножены на {}", options.density_multiplier));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parsers::lira::LiraParser;

    /// The stiffness block of the converter's test5.txt (LIRA 2022), with its block 13.
    const TEST5: &str = "( 3/\n1 30000 7.50001 7.50001 73804.1 0 0 /\n 0 RO 0.025/\n 0 S0 3.00001e+006 10 10/\n 0 Mu 0.2/\n2 GEI 3e+006 0.2 0.18 RO 2.5 /\n 0 WLKE 1 WLKG 1 PLKE 0.3 PLKG 0.3 /\n3 49154.4 270.979 17.3511 0.39419 0 0 /\n0 Y 0.00784425 0.00784425 Z 0.0612536 0.0612536 /\n0 RO 0.0184/\n0 S8 13 9 18 0 0/\n4 S3 0.305915 3 15 15 2 10 2/\n 0 RO 0.0220055/\n5 37306.9 53.6854 53.6854 35.1978 0 0 /\n0 Y 0.0287804 0.0287804 Z 0.0287804 0.0287804 /\n0 RO 0.01394/\n0 S8 13 10 10 0 0/\n6 342106 942.161 942.161 712.404 0 0 /\n0 Y 0.0306 0.0306 Z 0.0306 0.0306 /\n0 RO 0.127793/\n0 S8 13 18 18 0 0/\n7 134439 1282.56 77.3727 5.70476 0 0 /\n0 Y 0.0115104 0.0115104 Z 0.0781974 0.0781974 /\n0 RO 0.05024/\n0 S8 13 10 24.4 0 0/\n8 S5 3.05915e+007 20 20 2 3/\n 0 RO 0.0466622/\n9 S6 3.05915e+007 20 18/\n 0 RO 0.0158255/\n\n )\n\n{13/\n[3]\nSection = DoubleT \nMatId = STL \nFile  = |DSTU_8768_2018_short_UK.srt| \nShape = |18| \n\n[5]\nSection = Tubing \nMatId = STL \nFile  = |gn-kv_1.profiles.srt| \nShape = |100 x 5| \n\n[6]\nSection = Pipe \nMatId = STL \nFile  = |TRUBA2.profiles.srt| \nShape = |180 x 36| \n\n[7]\nSection = DoubleTBuiltup \nMatId = STL \nSteel  = |as flange has| \n\nFlange = Sheet \nFlangeShape = |100 x 22| \n\nFlangeRotation = 3\tFlangeBluePoint = 4 \tFlangeYellowPoint = 6\n\nWall = Sheet \nWallShape = |200 x 10| \n\n}\n";

    fn stiffness(mode: Mode) -> Stiffness {
        let materials = LiraParser::materials_from(TEST5.as_bytes());
        let profiles = LiraParser::profiles_from(TEST5.as_bytes());
        let ids: BTreeSet<u32> = (1..=9).collect();
        build(&materials, &profiles, &ids, &ids, Options { mode, density_multiplier: 1. })
    }

    #[test]
    fn the_converters_sections_of_its_test_model() {
        let s = stiffness(Mode::Converter);
        let text = s.section_lines.join("\n");
        // Lines of the converter's output for test5.txt (its first section is a different one: it reads the S0 row wrongly).
        for expected in [
            "DBUSER, DoubleT_18, CC, 0, 0, 0, 0, 0, 0, YES, H, 2, 0.1800, 0.0900, 0.0051, 0.0081, 0.0900, 0.0081, 0, 0, 0, 0",
            "DBUSER, I-Shape_0.150_X_0.15, CC, 0, 0, 0, 0, 0, 0, YES, H, 2, 0.1500, 0.1500, 0.0300, 0.0200, 0.1000, 0.0200, 0, 0, 0, 0",
            "DBUSER, Tube_100x5, CC, 0, 0, 0, 0, 0, 0, YES, B, 2, 0.1000, 0.1000, 0.0050, 0.0050, 0.0050, 0.0050, 0, 0, 0, 0",
            "DBUSER, Pipe_180x36, CC, 0, 0, 0, 0, 0, 0, YES, P, 2, 0.1800, 0.0360, 0, 0, 0, 0, 0, 0, 0, 0",
            "DBUSER, Flange100_x_22_Wall200_x_10, CC, 0, 0, 0, 0, 0, 0, YES, H, 2, 0.2440, 0.1000, 0.0100, 0.0220, 0.1000, 0.0220, 0, 0, 0, 0",
            "DBUSER, Box_0.200_X_0.2, CC, 0, 0, 0, 0, 0, 0, YES, B, 2, 0.2000, 0.2000, 0.0200, 0.0200, 0.0300, 0.0300, 0, 0, 0, 0",
            "DBUSER, Pipe_0.200_x_0.01, CC, 0, 0, 0, 0, 0, 0, YES, P, 2, 0.2000, 0.0100, 0, 0, 0, 0, 0, 0, 0, 0",
        ] {
            assert!(text.contains(expected), "{expected}\n{text}");
        }
        // The S0 row 10 x 10 cm: its own dimensions, not EIy of the numeric row.
        assert!(text.contains("DBUSER, Rect_0.100_X_0.1, CC, 0, 0, 0, 0, 0, 0, YES, SB, 2, 0.1000, 0.1000"), "{text}");
        // Steel sections of the converter's materials; the plate 3e6 t/m2 is concrete (Plate_Beton_p2_h0.18).
        let materials = s.material_lines.join("\n");
        assert!(materials.contains("Plate_Beton_p2_h0.18, 0, 0, , C, NO, 0, 2, 29419950.00, 0.2, 0, 24.52, 0"), "{materials}");
        assert!(materials.contains("Beam_Steel_p3, 0, 0, , C, NO, 0, 2, 206000000.00, 0.3, 0, 76.98, 0"), "{materials}");
        // One thickness line, after the sections.
        assert_eq!(s.thickness_lines.len(), 1);
        assert!(s.thickness_lines[0].ends_with("VALUE, YES, 0.180, 0, NO, 0, 0, Plate_0.18"), "{:?}", s.thickness_lines);
    }

    #[test]
    fn plates_keep_their_lira_factors_only_in_the_lira_mode_or_on_the_cap() {
        // PLKE 0.3: the converter's file keeps the nominal 0.18 m; the Lira mode an equivalent thickness.
        assert!(stiffness(Mode::Converter).thickness_lines[0].contains("0.180"));
        let lira = stiffness(Mode::Lira);
        let t = lira.thickness_lines[0].split(", ").nth(3).unwrap().parse::<f64>().unwrap();
        // d' = d sqrt(kb / km) = 0.18 sqrt(0.3).
        assert!((t - 0.18 * 0.3f64.sqrt()).abs() < 1e-3, "{t}");
        assert!(stiffness(Mode::Converter).notes.iter().any(|n| n.contains("не применена")));
    }

    #[test]
    fn concrete_round_columns_are_not_turned_into_steel() {
        // S6 of the island model: a solid round concrete column (E = 2.37e6 t/m2).
        let text = "( 3/\n6 S6 2.37e+006 28 0/\n 0 RO 0.169332/\n 0 Mu 0.2/\n)\n";
        let materials = LiraParser::materials_from(text.as_bytes());
        let ids: BTreeSet<u32> = [6].into();
        let s = build(&materials, &HashMap::new(), &BTreeSet::new(), &ids, Options::default());
        let m = s.material_lines.join("\n");
        assert!(m.contains("Beam_Beton_p6"), "{m}");
        assert!(m.contains("23240") || m.contains("2324"), "E = 2.37e6 t/m2 x 9.80665 = 23 242 000 kN/m2: {m}");
        assert!(s.section_lines[0].contains("SR, 2, 0.2800"), "{:?}", s.section_lines);
        assert!(s.notes.iter().all(|n| !n.contains("сталь")), "{:?}", s.notes);
    }

    #[test]
    fn equal_sections_share_one_and_zero_types_are_reported() {
        let text = "( 3/\n1 S0 3e+006 40 90/\n 0 RO 0.9/\n2 S0 3e+006 40 90/\n 0 RO 0.9/\n)\n";
        let materials = LiraParser::materials_from(text.as_bytes());
        let ids: BTreeSet<u32> = [1, 2, 7].into();
        let s = build(&materials, &HashMap::new(), &BTreeSet::new(), &ids, Options::default());
        assert_eq!(s.sections, 1);
        assert_eq!(s.section_of[&1], s.section_of[&2]);
        assert_eq!(s.material_of[&(false, 1)], s.material_of[&(false, 2)]);
        assert!(s.missing.contains(&7));
    }

    #[test]
    fn a_reduced_type_is_a_value_section_in_the_lira_mode_and_nominal_in_the_converters() {
        // test1.txt of the converter: EF, EIy, EIz, GIk of the row and S0 90 x 170 cm; EIy, EIz are 0.3 of the nominal.
        let text = "( 3/\n1 4.59e+006 331628 92947.5 73804.1 0 0 /\n 0 RO 3.825/\n 0 S0 3e+006 90 170/\n 0 Mu 0.2/\n)\n";
        let materials = LiraParser::materials_from(text.as_bytes());
        let ids: BTreeSet<u32> = [1].into();
        let converter = build(&materials, &HashMap::new(), &BTreeSet::new(), &ids, Options::default());
        // The dimensions of the S0 row (the converter reads EIy of the numeric row as the width here).
        assert!(converter.section_lines[0].contains("DBUSER, Rect_0.900_X_1.7, CC, 0, 0, 0, 0, 0, 0, YES, SB, 2, 1.7000, 0.9000"), "{:?}", converter.section_lines);
        assert!(converter.notes.iter().any(|n| n.contains("не применена") && n.contains('1')), "{:?}", converter.notes);
        let lira = build(&materials, &HashMap::new(), &BTreeSet::new(), &ids, Options { mode: Mode::Lira, density_multiplier: 1. });
        assert!(lira.section_lines[0].starts_with("    1, VALUE, Rect_0.900_X_1.7"), "{:?}", lira.section_lines);
        // A = EF / E = 1.53 m2, Iyy = EIy / E = 0.110543 m4, Izz = 0.030983 m4.
        assert!(lira.section_lines[1].contains("1.53000000, 0, 0,") && lira.section_lines[1].contains("0.11054267, 0.03098250"), "{:?}", lira.section_lines);
        // The unit weight of the bar: RO 3.825 t/m over 1.53 m2 = 2.5 t/m3.
        assert!(lira.material_lines[0].contains(", 24.52, 0"), "{:?}", lira.material_lines);
    }

    #[test]
    fn a_pile_with_a_stiffness_row_takes_e_from_ef_over_a_and_stays_concrete() {
        // «Для testa», type 392: EF = 2.77735e6 t, S6 line "720000 92 0" (D = 92 cm). E of the row is EF / A = 4.18e6 t/m2.
        let text = "( 3/\n392 2.77735e+006 141350 141350 191790 0 0 /\n 0 RO 2.22695/\n 0 S6 720000 92 0/\n 0 Mu 0.2/\n)\n";
        let materials = LiraParser::materials_from(text.as_bytes());
        let ids: BTreeSet<u32> = [392].into();
        let s = build(&materials, &HashMap::new(), &BTreeSet::new(), &ids, Options::default());
        // The solid round section of 92 cm, not a 1413 m rod.
        assert!(s.section_lines[0].contains("SR, 2, 0.9200"), "{:?}", s.section_lines);
        let a = std::f64::consts::PI / 4. * 0.92 * 0.92;
        let e = 2.77735e6 / a * 9.80665;
        assert!(s.material_lines[0].contains(&format!("{e:.2}")) && s.material_lines[0].contains("Beam_"), "{:?} (E = {e:.2})", s.material_lines);
        assert!(s.notes.iter().all(|n| !n.contains("сталь")) && s.notes.iter().any(|n| n.contains("EF / A")), "{:?}", s.notes);
    }

    #[test]
    fn the_s0_row_of_the_converters_test5_is_10_by_10_cm_with_its_e() {
        // "1 30000 7.5 7.5 73804 / S0 3.00001e+006 10 10": the converter writes 0.075 and E = EF x g.
        let text = "( 3/\n1 30000 7.50001 7.50001 73804.1 0 0 /\n 0 RO 0.025/\n 0 S0 3.00001e+006 10 10/\n 0 Mu 0.2/\n)\n";
        let materials = LiraParser::materials_from(text.as_bytes());
        let ids: BTreeSet<u32> = [1].into();
        let s = build(&materials, &HashMap::new(), &BTreeSet::new(), &ids, Options::default());
        assert!(s.section_lines[0].contains("SB, 2, 0.1000, 0.1000"), "{:?}", s.section_lines);
        assert!(s.material_lines[0].contains("29420009.80") || s.material_lines[0].contains("294200"), "{:?}", s.material_lines);
        assert!(s.material_lines[0].contains(", 24.52, 0"), "{:?}", s.material_lines);
    }

    #[test]
    fn a_rigid_rod_is_named_agt_and_the_density_multiplier_applies() {
        let text = "( 3/\n1 S0 1e+009 100 100/\n 0 RO 2.5/\n)\n";
        let materials = LiraParser::materials_from(text.as_bytes());
        let ids: BTreeSet<u32> = [1].into();
        let s = build(&materials, &HashMap::new(), &BTreeSet::new(), &ids, Options { mode: Mode::Converter, density_multiplier: 0. });
        assert!(s.section_lines[0].contains("_AGT"), "{:?}", s.section_lines);
        // A rigid rod: a weightless material of its own, as in the converter's corrected writer.
        assert!(s.material_lines[0].contains("AGT_p1") && s.material_lines[0].contains(", 0.00, 0"), "{:?}", s.material_lines);
    }
}
