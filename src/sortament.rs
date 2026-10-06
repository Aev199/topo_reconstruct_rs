//! The LIRA catalogue of rolled profiles (`data/sortament.tsv`, the one of the
//! Lira_Midas-converter, 5720 entries: I sections, square, rectangular and
//! round tubes, angles) with the dimensions (mm) and the section properties of
//! each entry, and the lookup of a profile by the section name and designation
//! of the block `{13/` of a LIRA file (`DoubleT 20Ш1`, `Tubing 100 x 5`,
//! `Pipe 273 x 20`).
use std::collections::HashMap;
use std::sync::OnceLock;

/// An entry of the catalogue (cm, cm2, cm4 as in the catalogue; mm for the dimensions).
#[derive(Debug, Clone, PartialEq)]
pub struct Entry {
    pub name: String,
    /// `нормальный (балочный)`, `квадратная`, `Труба круглая`, ...
    pub kind: String,
    pub h: Option<f64>,
    pub b: Option<f64>,
    /// Web thickness (mm); absent for tubes.
    pub s: Option<f64>,
    /// Flange thickness, or the wall of a tube (mm).
    pub t: Option<f64>,
    /// Fillet radius (mm).
    pub r: Option<f64>,
    /// Area (cm2) and mass (kg/m).
    pub area: Option<f64>,
    pub mass: Option<f64>,
    /// Second moments about the strong axis X (parallel to the flanges) and the weak axis Y (cm4).
    pub ix: Option<f64>,
    pub iy: Option<f64>,
}

/// Kinds of the catalogue that are I sections (the converter's `_SORTAMENT_TO_SECTION_TYPE`).
pub fn is_i_kind(kind: &str) -> bool {
    matches!(
        kind,
        "нормальный (балочный)"
            | "балочный нормальный"
            | "балочный широкополочный"
            | "широкополочный"
            | "с уклоном полок"
            | "с параллельными гранями полок"
            | "колонный"
            | "дополнительной серии балочный"
            | "дополнительной серии колонный"
            | "экономичный с параллельными гранями полок"
            | "легкой серии с параллельными гранями полок"
            | "специальный"
            | "свайный"
            | "дополнительной серии"
    )
}

struct Catalogue {
    entries: Vec<Entry>,
    by_name: HashMap<String, usize>,
    by_prefix: HashMap<String, Vec<usize>>,
}

fn catalogue() -> &'static Catalogue {
    static CATALOGUE: OnceLock<Catalogue> = OnceLock::new();
    CATALOGUE.get_or_init(|| {
        let mut entries = vec![];
        // The first line names the columns.
        for line in include_str!("../data/sortament.tsv").lines().skip(1) {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 14 {
                continue;
            }
            let num = |k: usize| f[k].trim().parse::<f64>().ok();
            entries.push(Entry {
                name: f[0].to_string(),
                kind: f[1].to_string(),
                h: num(2),
                b: num(3),
                s: num(4),
                t: num(5),
                r: num(6),
                area: num(7),
                mass: num(8),
                ix: num(9),
                iy: num(13),
            });
        }
        let mut by_name = HashMap::new();
        let mut by_prefix: HashMap<String, Vec<usize>> = HashMap::new();
        for (i, e) in entries.iter().enumerate() {
            by_name.insert(e.name.clone(), i);
            let prefix = numeric_prefix(&e.name);
            if !prefix.is_empty() {
                by_prefix.entry(prefix).or_default().push(i);
            }
        }
        Catalogue { entries, by_name, by_prefix }
    })
}

/// `18Б1` -> `18`, `100x5` -> `100`.
fn numeric_prefix(name: &str) -> String {
    name.chars().take_while(|c| c.is_ascii_digit() || *c == '.').collect()
}

/// A designation with the cyrillic x and the multiplication sign as latin `x`.
pub fn ascii_x(s: &str) -> String {
    s.replace(['х', 'Х', '×', 'X'], "x")
}

impl Catalogue {
    fn exact(&self, name: &str) -> Option<usize> {
        let variants = [
            name.to_string(),
            name.replace(" x ", "x"),
            name.replace(" x ", "х"),
            name.replace('х', "x"),
            name.replace('x', "х"),
        ];
        variants.iter().find_map(|v| self.by_name.get(v).copied())
    }

    /// A tube from `100 x 5` (round, or square `100x100x5`) or `120 x 80 x 5` (rectangular).
    fn tube(&self, shape: &str, square: bool) -> Option<usize> {
        let numbers: Vec<f64> = ascii_x(shape).split('x').filter_map(|p| p.trim().parse().ok()).collect();
        if numbers.len() == 3 {
            let name = |v: f64| if (v - v.round()).abs() < 1e-9 { format!("{}", v.round() as i64) } else { format!("{v}") };
            return self.exact(&format!("{}x{}x{}", name(numbers[0]), name(numbers[1]), name(numbers[2])));
        }
        if numbers.len() < 2 {
            return None;
        }
        let (a, t) = (numbers[0], numbers[1]);
        let search = if square { format!("{}x{}x{}", a as i64, a as i64, t as i64) } else { format!("{}x{}", a as i64, t as i64) };
        if let Some(&i) = self.by_name.get(&search) {
            return Some(i);
        }
        let kind = if square { "квадрат" } else { "труба круглая" };
        self.by_prefix.get(&format!("{}", a as i64))?.iter().copied().find(|&i| {
            let e = &self.entries[i];
            e.kind.to_lowercase().contains(kind) && e.t.is_some_and(|c| (c - t).abs() < 0.5)
        })
    }

    /// The `lookup_by_section_info` of the converter.
    fn lookup(&self, section: &str, shape: &str) -> Option<usize> {
        let (section, shape) = (section.to_lowercase(), shape.trim());
        if let Some(i) = self.exact(shape) {
            return Some(i);
        }
        if section.contains("doublet") || section.contains("двут") {
            if let Some(&i) = self.by_name.get(&format!("{shape}.0")) {
                if is_i_kind(&self.entries[i].kind) {
                    return Some(i);
                }
            }
            let candidates = self.by_prefix.get(&numeric_prefix(shape))?;
            if let Some(&i) = candidates.iter().find(|&&i| self.entries[i].name.replace(".0", "") == shape) {
                return Some(i);
            }
            return candidates.iter().copied().find(|&i| is_i_kind(&self.entries[i].kind)).or(candidates.first().copied());
        }
        if section.contains("tubing") || section.contains("tube") {
            return self.tube(shape, true);
        }
        if section.contains("pipe") {
            return self.tube(shape, false);
        }
        None
    }
}

/// The entry of a profile of the block 13: `section` is `DoubleT`, `Tubing`, `Pipe`, `shape` its designation.
pub fn lookup(section: &str, shape: &str) -> Option<&'static Entry> {
    let cat = catalogue();
    cat.lookup(section, shape).map(|i| &cat.entries[i])
}

/// The entry named exactly so (`20Б1`, `100x100x5`).
pub fn find(name: &str) -> Option<&'static Entry> {
    let cat = catalogue();
    cat.exact(name).map(|i| &cat.entries[i])
}

/// Number of entries.
pub fn len() -> usize {
    catalogue().entries.len()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_catalogue_is_complete_and_has_the_section_properties() {
        assert_eq!(len(), 5720);
        let i18 = find("18Б1").expect("I 18B1");
        // 18Б1 (GOST 26020): A = 19.58 cm2, Ix = 1063 cm4, Iy = 81.9 cm4, 15.4 kg/m.
        assert!((i18.area.unwrap() - 19.58).abs() < 1e-9 && (i18.ix.unwrap() - 1063.).abs() < 1e-9 && (i18.iy.unwrap() - 81.9).abs() < 1e-9 && (i18.mass.unwrap() - 15.4).abs() < 1e-9, "{i18:?}");
        let tube = find("100x100x5").expect("tube");
        // 100 x 100 x 5: A = 18.4 cm2.
        assert!((tube.area.unwrap() - 18.4).abs() < 0.5 && tube.r.is_some(), "{tube:?}");
        assert!(lookup("Tubing", "120 x 80 x 5").is_some());
        // A pipe 273 x 20 is not in the catalogue (the converter takes the designation then).
        assert!(lookup("Pipe", "273 x 20").is_none());
    }
}
