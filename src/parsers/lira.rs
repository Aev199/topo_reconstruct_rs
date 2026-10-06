use crate::input::{ElementData, MeshData};
use fast_float::parse as fast_parse_f64;
use glam::DVec3;
use hashbrown::HashMap;
use memmap2::Mmap;
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, Error, ErrorKind};
use std::path::Path;

pub struct LiraParser;

/// Габариты сечения из блока жёсткостей ЛИРА, в метрах.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Section {
    Plate { thickness: f64 },
    Bar { width: f64, height: f64 },
}

/// Material of a stiffness type of the LIRA stiffness block (3/), in model
/// units (no units block is read: by default t and m, E in t/m2).
#[derive(Debug, Clone, PartialEq, serde::Serialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum Material {
    /// GEI: shell of thickness `thickness` (m); `density` is per volume.
    /// `membrane` (WLKE) and `bending` (PLKE) scale E for in-plane and
    /// out-of-plane work; WLKG/PLKG scale G (not representable in PLAXIS).
    Plate {
        e: f64,
        nu: f64,
        thickness: f64,
        density: Option<f64>,
        membrane: Option<f64>,
        bending: Option<f64>,
        shear: [Option<f64>; 2],
    },
    /// S0: rectangular bar `width` (along local Y1) x `height` (along Z1),
    /// m (given in cm); `density` (RO) is per length. `stiffness`: the
    /// numeric EF, EIy, EIz, GIk of the type when given (they may differ
    /// from E and the S0 dimensions: LIRA allows independent stiffness).
    Bar {
        e: f64,
        nu: Option<f64>,
        width: f64,
        height: f64,
        density: Option<f64>,
        stiffness: Option<[f64; 4]>,
    },
    /// A bar of another section of the stiffness block: S1, S2, S3, S5, S6
    /// (parametric, the dimensions in `shape`), a named profile of the
    /// block 13 or a type with only the numeric EF, EIy, EIz, GIk
    /// (`Shape::Explicit`; `section` is `S8` for a rolled profile, `explicit` else).
    /// `e` is the Young's modulus of the type when it has one; `hint` the
    /// width and height (m) of an S8 row.
    Section {
        section: String,
        e: Option<f64>,
        nu: Option<f64>,
        density: Option<f64>,
        shape: crate::sections::Shape,
        stiffness: Option<[f64; 4]>,
        hint: Option<[f64; 2]>,
    },
}

impl LiraParser {
    /// Materials of the stiffness block: a row `id ...` with its
    /// continuation rows `0 KEY values` (RO, Mu, S0).
    pub fn parse_materials<P: AsRef<Path>>(filepath: P) -> io::Result<HashMap<u32, Material>> {
        let file = File::open(filepath)?;
        let mmap = unsafe { Mmap::map(&file)? };
        Ok(Self::materials_bytes(&mmap))
    }

    /// The named profiles of the block `{13/` (by the number of the stiffness type).
    pub fn profiles_from(content: &[u8]) -> HashMap<u32, crate::sections::Profile> {
        let text = String::from_utf8_lossy(content);
        let mut out = HashMap::new();
        let Some(start) = ["{13/", "{ 13/", "{13 /", "{ 13 /"].iter().filter_map(|m| text.find(m).map(|i| i + m.len())).min() else {
            return out;
        };
        let end = text[start..].find('}').map_or(text.len(), |e| start + e);
        let body = &text[start..end];
        let value = |block: &str, key: &str| -> Option<String> {
            let at = block.find(key)?;
            let rest = &block[at + key.len()..];
            let rest = rest.trim_start().strip_prefix('=')?.trim_start();
            if let Some(quoted) = rest.strip_prefix('|') {
                return Some(quoted.split('|').next()?.trim().to_string());
            }
            Some(rest.split_whitespace().next()?.to_string())
        };
        for block in body.split('[').skip(1) {
            let Some((head, rest)) = block.split_once(']') else { continue };
            let Ok(id) = head.trim().parse::<u32>() else { continue };
            let section = value(rest, "Section").unwrap_or_default();
            let builtup = rest.contains("Flange") || rest.contains("Wall");
            let shape = if builtup {
                match (value(rest, "FlangeShape"), value(rest, "WallShape")) {
                    (Some(f), Some(w)) => format!("Flange={f}; Wall={w}"),
                    _ => value(rest, "Shape").unwrap_or_default(),
                }
            } else {
                value(rest, "Shape").unwrap_or_default()
            };
            out.insert(id, crate::sections::Profile { section, shape, builtup });
        }
        out
    }

    /// Materials of input bytes already read (the opened snapshot).
    pub fn materials_from(content: &[u8]) -> HashMap<u32, Material> {
        Self::materials_bytes(content)
    }

    fn materials_bytes(content: &[u8]) -> HashMap<u32, Material> {
        let Some(block) = Self::extract_block(content, b"3") else {
            return HashMap::new();
        };
        // Words of each stiffness type, continuation rows appended.
        let mut types: Vec<(u32, Vec<String>)> = vec![];
        for row in block.split(|&b| b == b'/') {
            let mut words = Self::split_ascii_whitespace_bytes(row)
                .filter_map(|w| std::str::from_utf8(w).ok().map(str::to_string));
            let Some(Ok(id)) = words.next().map(|w| w.parse::<u32>()) else {
                continue;
            };
            if id == 0 {
                if let Some(last) = types.last_mut() {
                    last.1.extend(words);
                }
            } else {
                types.push((id, words.collect()));
            }
        }
        let number = |w: Option<&String>| {
            w.and_then(|w| fast_parse_f64::<f64, _>(w.as_bytes()).ok())
                .filter(|v| v.is_finite())
        };
        let after = |words: &[String], key: &str, k: usize| {
            words
                .iter()
                .position(|w| w.eq_ignore_ascii_case(key))
                .and_then(|i| number(words.get(i + k)))
        };
        types
            .into_iter()
            .filter_map(|(id, words)| {
                let density = after(&words, "RO", 1);
                let material = if let Some(i) = words.iter().position(|w| w == "GEI") {
                    Material::Plate {
                        e: number(words.get(i + 1)).filter(|v| *v > 0.)?,
                        nu: number(words.get(i + 2))?,
                        thickness: number(words.get(i + 3)).filter(|v| *v > 0.)?,
                        density,
                        membrane: after(&words, "WLKE", 1),
                        bending: after(&words, "PLKE", 1),
                        shear: [after(&words, "WLKG", 1), after(&words, "PLKG", 1)],
                    }
                } else if words.iter().any(|w| w == "S0") {
                    // A row starting with numbers: EF EIy EIz GIk ...
                    let numeric: Vec<f64> = words.iter().map_while(|w| number(Some(w))).collect();
                    Material::Bar {
                        e: after(&words, "S0", 1).filter(|v| *v > 0.)?,
                        width: after(&words, "S0", 2).filter(|v| *v > 0.)? / 100.,
                        height: after(&words, "S0", 3).filter(|v| *v > 0.)? / 100.,
                        nu: after(&words, "Mu", 1),
                        density,
                        // EF > 0; a bending stiffness of 0 is a real value (no bending in that plane).
                        stiffness: (numeric.len() >= 4 && numeric[0] > 0. && numeric[1] >= 0. && numeric[2] >= 0.)
                            .then(|| [numeric[0], numeric[1], numeric[2], numeric[3]]),
                    }
                } else if let Some((kind, i)) = ["S1", "S2", "S3", "S5", "S6"]
                    .iter()
                    .find_map(|k| words.iter().position(|w| w == k).map(|i| (*k, i)))
                {
                    use crate::sections::{cm, s6, Shape};
                    // E and the dimensions (cm, or m when given below 1) after the marker.
                    let v: Vec<f64> = words[i + 1..].iter().map_while(|w| number(Some(w))).collect();
                    let shape = match (kind, v.len()) {
                        ("S1", n) if n >= 5 => Shape::ISym { h: cm(v[2]), b: cm(v[1]), tw: cm(v[3]), tf: cm(v[4]) },
                        ("S2", n) if n >= 5 => Shape::T { h: cm(v[2]), b: cm(v[1]), tw: cm(v[3]), tf: cm(v[4]) },
                        ("S3", n) if n >= 7 => Shape::I { h: v[2] / 100., bt: v[3] / 100., tw: v[1] / 100., tft: v[4] / 100., bb: v[5] / 100., tfb: v[6] / 100. },
                        ("S5", n) if n >= 5 => Shape::Box { h: cm(v[1]), b: cm(v[2]), t1: cm(v[3]), t2: v[4] / 100. },
                        ("S6", n) if n >= 3 => s6(v[1], v[2]),
                        _ => return None,
                    };
                    Material::Section {
                        section: kind.to_string(),
                        e: v.first().copied().filter(|e| *e > 0.),
                        nu: after(&words, "Mu", 1).or_else(|| after(&words, "NU", 1)),
                        density,
                        shape,
                        stiffness: None,
                        hint: None,
                    }
                } else {
                    // Only the numeric row EF EIy EIz GIk (a rolled profile of the block 13 has S8 after it).
                    let numeric: Vec<f64> = words.iter().map_while(|w| number(Some(w))).collect();
                    if numeric.len() < 4 || numeric[0] <= 0. || numeric[1] < 0. || numeric[2] < 0. {
                        return None;
                    }
                    let s8 = words.iter().position(|w| w == "S8");
                    let hint = s8.and_then(|i| Some([number(words.get(i + 2))? / 100., number(words.get(i + 3))? / 100.]));
                    Material::Section {
                        section: if s8.is_some() { "S8".into() } else { "explicit".into() },
                        e: None,
                        nu: after(&words, "Mu", 1),
                        density,
                        shape: crate::sections::Shape::Explicit,
                        stiffness: Some([numeric[0], numeric[1], numeric[2], numeric[3]]),
                        hint,
                    }
                };
                Some((id, material))
            })
            .collect()
    }

    /// Потоковый параллельный парсинг текстового файла ЛИРА (.txt)
    pub fn parse<P: AsRef<Path>>(filepath: P) -> io::Result<MeshData> {
        let file = File::open(filepath)?;
        let mmap = unsafe { Mmap::map(&file)? };
        Self::parse_bytes(&mmap)
    }

    /// The mesh of input bytes already read (the opened snapshot).
    pub fn mesh_from(content: &[u8]) -> io::Result<MeshData> {
        Self::parse_bytes(content)
    }

    fn parse_bytes(content: &[u8]) -> io::Result<MeshData> {
        let nodes = Self::extract_block(content, b"4").ok_or_else(|| {
            Error::new(
                ErrorKind::InvalidData,
                "Блок координат узлов (4/) не найден",
            )
        })?;
        let elements = Self::extract_block(content, b"1")
            .ok_or_else(|| Error::new(ErrorKind::InvalidData, "Блок элементов (1/) не найден"))?;
        let node_rows: Vec<_> = nodes
            .split(|&b| b == b'/')
            .filter(|r| !Self::is_empty_or_ws(r))
            .collect();
        let parsed_nodes: io::Result<Vec<(u32, DVec3)>> = node_rows
            .par_iter()
            .enumerate()
            .map(|(i, row)| {
                let words: Vec<_> = Self::split_ascii_whitespace_bytes(row).collect();
                let error = || {
                    Error::new(
                        ErrorKind::InvalidData,
                        format!("Узел {}: ожидаются три конечные координаты", i + 1),
                    )
                };
                if words.len() != 3 {
                    return Err(error());
                }
                let mut xyz = [0.0; 3];
                for j in 0..3 {
                    xyz[j] = fast_parse_f64::<f64, _>(words[j]).map_err(|_| error())?;
                    if !xyz[j].is_finite() {
                        return Err(error());
                    }
                }
                Ok(((i + 1) as u32, DVec3::from_array(xyz)))
            })
            .collect();
        let nodes: HashMap<u32, DVec3> = parsed_nodes?.into_iter().collect();
        let element_rows: Vec<_> = elements
            .split(|&b| b == b'/')
            .filter(|r| !Self::is_empty_or_ws(r))
            .collect();
        let parsed_elements: io::Result<Vec<ElementData>> = element_rows
            .par_iter()
            .enumerate()
            .map(|(i, row)| {
                let error = || {
                    Error::new(
                        ErrorKind::InvalidData,
                        format!("КЭ {}: неверная запись типа, жесткости или узлов", i + 1),
                    )
                };
                let ints: io::Result<Vec<u32>> = Self::split_ascii_whitespace_bytes(row)
                    .map(|w| {
                        std::str::from_utf8(w)
                            .map_err(|_| error())?
                            .parse::<u32>()
                            .map_err(|_| error())
                    })
                    .collect();
                let ints = ints?;
                if ints.len() < 3 {
                    return Err(error());
                }
                if ints[2..].iter().any(|id| !nodes.contains_key(id)) {
                    return Err(Error::new(
                        ErrorKind::InvalidData,
                        format!("КЭ {}: ссылка на отсутствующий узел", i + 1),
                    ));
                }
                let element = ElementData {
                    id: (i + 1) as u32,
                    elem_type: ints[0],
                    stiff_id: ints[1],
                    nodes: ints[2..].to_vec(),
                };
                let expected = match element.elem_type {
                    10 => Some(2),
                    41 | 44 => Some(4),
                    42 => Some(3),
                    56 => Some(1),
                    _ => None,
                };
                if expected.is_some_and(|n| n != element.nodes.len()) {
                    return Err(error());
                }
                Ok(element)
            })
            .collect();
        Ok(MeshData {
            nodes,
            elements: parsed_elements?,
        })
    }

    /// Численно заданные жёсткости (блок 3/): номер -> значения (EF, EIy,
    /// EIz, GIk, ...). Жёсткости, заданные типом сечения (S0, GEI, ...),
    /// пропускаются; файл без блока даёт пустой набор.
    pub fn parse_stiffness<P: AsRef<Path>>(filepath: P) -> io::Result<HashMap<u32, Vec<f64>>> {
        let file = File::open(filepath)?;
        let mmap = unsafe { Mmap::map(&file)? };
        Ok(Self::stiffness_bytes(&mmap))
    }

    fn stiffness_bytes(content: &[u8]) -> HashMap<u32, Vec<f64>> {
        let Some(block) = Self::extract_block(content, b"3") else {
            return HashMap::new();
        };
        block
            .split(|&b| b == b'/')
            .filter_map(|row| {
                let mut words = Self::split_ascii_whitespace_bytes(row);
                let id = std::str::from_utf8(words.next()?)
                    .ok()?
                    .parse::<u32>()
                    .ok()?;
                let values: Option<Vec<f64>> = words
                    .map(|w| fast_parse_f64::<f64, _>(w).ok().filter(|v| v.is_finite()))
                    .collect();
                values.filter(|v| !v.is_empty()).map(|v| (id, v))
            })
            .collect()
    }

    /// Абсолютно жёсткие тела (блок 25/): группы узлов, первый — ведущий.
    /// Файл без блока даёт пустой список.
    pub fn parse_rigid_bodies<P: AsRef<Path>>(filepath: P) -> io::Result<Vec<Vec<u32>>> {
        let file = File::open(filepath)?;
        let mmap = unsafe { Mmap::map(&file)? };
        Ok(Self::rigid_bodies_bytes(&mmap))
    }

    fn rigid_bodies_bytes(content: &[u8]) -> Vec<Vec<u32>> {
        let Some(block) = Self::extract_block(content, b"25") else {
            return vec![];
        };
        block
            .split(|&b| b == b'/')
            .filter_map(|row| {
                let nodes: Option<Vec<u32>> = Self::split_ascii_whitespace_bytes(row)
                    .map(|w| std::str::from_utf8(w).ok()?.parse::<u32>().ok())
                    .collect();
                nodes.filter(|n| n.len() >= 2)
            })
            .collect()
    }

    /// Толщины пластин и габариты сечений стержней (блок 3/), в метрах:
    /// `GEI E nu H` — пластина толщиной H (м), `S0 E b h` — прямоугольное
    /// сечение b x h (см). Прочие жёсткости пропускаются.
    pub fn parse_sections<P: AsRef<Path>>(filepath: P) -> io::Result<HashMap<u32, Section>> {
        let file = File::open(filepath)?;
        let mmap = unsafe { Mmap::map(&file)? };
        Ok(Self::sections_bytes(&mmap))
    }

    fn sections_bytes(content: &[u8]) -> HashMap<u32, Section> {
        let Some(block) = Self::extract_block(content, b"3") else {
            return HashMap::new();
        };
        block
            .split(|&b| b == b'/')
            .filter_map(|row| {
                let mut words = Self::split_ascii_whitespace_bytes(row);
                let id = std::str::from_utf8(words.next()?)
                    .ok()?
                    .parse::<u32>()
                    .ok()?;
                let kind = words.next()?;
                let mut number = || {
                    words
                        .next()
                        .and_then(|w| fast_parse_f64::<f64, _>(w).ok())
                        .filter(|v| v.is_finite() && *v > 0.)
                };
                let section = match kind {
                    b"GEI" => {
                        let (_, _, h) = (number()?, number(), number()?);
                        Section::Plate { thickness: h }
                    }
                    b"S0" => {
                        let (_, b, h) = (number()?, number()?, number()?);
                        Section::Bar {
                            width: b / 100.,
                            height: h / 100.,
                        }
                    }
                    _ => return None,
                };
                (id != 0).then_some((id, section))
            })
            .collect()
    }

    /// Быстрый поиск содержимого блока `( <id>/ ... )` без аллокаций строк
    pub(crate) fn extract_block<'a>(content: &'a [u8], block_id: &[u8]) -> Option<&'a [u8]> {
        let mut i = 0;
        let len = content.len();

        while i < len {
            if content[i] == b'(' {
                let mut j = i + 1;
                // Пропускаем пробелы после '('
                while j < len
                    && (content[j] == b' '
                        || content[j] == b'\t'
                        || content[j] == b'\r'
                        || content[j] == b'\n')
                {
                    j += 1;
                }

                // Проверяем ID блока
                if j + block_id.len() < len && &content[j..j + block_id.len()] == block_id {
                    let mut k = j + block_id.len();
                    // Пропускаем пробелы перед '/'
                    while k < len
                        && (content[k] == b' '
                            || content[k] == b'\t'
                            || content[k] == b'\r'
                            || content[k] == b'\n')
                    {
                        k += 1;
                    }

                    if k < len && content[k] == b'/' {
                        let start = k + 1;
                        // Ищем закрывающую ')'
                        let mut depth = 1;
                        let mut end = start;
                        while end < len && depth > 0 {
                            if content[end] == b'(' {
                                depth += 1;
                            } else if content[end] == b')' {
                                depth -= 1;
                                if depth == 0 {
                                    return Some(&content[start..end]);
                                }
                            }
                            end += 1;
                        }
                    }
                }
            }
            i += 1;
        }
        None
    }

    /// Проверка, пустой ли слайс байт
    pub(crate) fn is_empty_or_ws(slice: &[u8]) -> bool {
        slice
            .iter()
            .all(|&b| b == b' ' || b == b'\t' || b == b'\r' || b == b'\n')
    }

    /// Итератор по непустым словам (whitespace-separated bytes)
    pub(crate) fn split_ascii_whitespace_bytes(slice: &[u8]) -> impl Iterator<Item = &[u8]> {
        let mut i = 0;
        let len = slice.len();

        std::iter::from_fn(move || {
            // Пропуск начальных пробелов
            while i < len
                && (slice[i] == b' ' || slice[i] == b'\t' || slice[i] == b'\r' || slice[i] == b'\n')
            {
                i += 1;
            }
            if i >= len {
                return None;
            }
            let start = i;
            // Поиск конца слова
            while i < len
                && !(slice[i] == b' '
                    || slice[i] == b'\t'
                    || slice[i] == b'\r'
                    || slice[i] == b'\n')
            {
                i += 1;
            }
            Some(&slice[start..i])
        })
    }
}
#[cfg(test)]
mod material_tests {
    use super::*;

    #[test]
    fn plates_and_bars_with_continuation_rows() {
        let text = b"( 3/\n1 1.2e+006 38400 15000 22572.8 0 0 /\n 0 RO 1/\n 0 S0 3e+006 50 80/\n 0 Mu 0.2/\n7 GEI 3e+006 0.2 0.2 RO 2.5 /\n 0 WLKE 1 WLKG 1 PLKE 0.6 PLKG 0.6 /\n2 S0 0.305915 20 20/\n 0 RO 0.0101972/\n)";
        let m = LiraParser::materials_bytes(text);
        assert_eq!(
            m[&1],
            Material::Bar {
                e: 3e6,
                nu: Some(0.2),
                width: 0.5,
                height: 0.8,
                density: Some(1.),
                stiffness: Some([1.2e6, 38400., 15000., 22572.8]),
            }
        );
        assert_eq!(
            m[&7],
            Material::Plate {
                e: 3e6,
                nu: 0.2,
                thickness: 0.2,
                density: Some(2.5),
                membrane: Some(1.),
                bending: Some(0.6),
                shear: [Some(1.), Some(0.6)],
            }
        );
        assert_eq!(
            m[&2],
            Material::Bar {
                e: 0.305915,
                nu: None,
                width: 0.2,
                height: 0.2,
                density: Some(0.0101972),
                stiffness: None,
            }
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn preserves_single_node_and_unknown_elements_without_renumbering() {
        let mesh = LiraParser::parse_bytes(
            b"(4/0 0 0/1 0 0/0 1 0/)(1/56 1 1/200 1 1 2 3/42 2 1 2 3/10 3 1 2/)",
        )
        .unwrap();
        assert_eq!(mesh.elements.len(), 4);
        assert_eq!(mesh.elements[2].id, 3);
        assert!(!mesh.elements[1].is_shell());
        assert!(!mesh.elements[1].is_bar());
        assert!(mesh.elements[2].is_shell());
    }
    #[test]
    fn rigid_bodies_are_node_groups_of_block_25() {
        let groups =
            LiraParser::rigid_bodies_bytes(b"(4/0 0 0/)( 25/ 1 9198 9199/ 7741 27 13445/ )");
        assert_eq!(groups, vec![vec![1, 9198, 9199], vec![7741, 27, 13445]]);
    }
    #[test]
    fn plate_thickness_and_bar_sections_are_read_in_metres() {
        let sections = LiraParser::sections_bytes(
            b"( 3/ 1 GEI 2.34e+006 0.2 0.3 RO 2.5 / 2 S0 3e+006 40 90/ 0 RO 0.9/ \
              0 Mu 0.2/ 3 1000 200000 200000 200000 0 0 / )",
        );
        assert_eq!(sections.len(), 2);
        assert_eq!(sections[&1], Section::Plate { thickness: 0.3 });
        assert_eq!(
            sections[&2],
            Section::Bar {
                width: 0.4,
                height: 0.9
            }
        );
    }
    #[test]
    fn a_zero_bending_stiffness_of_a_bar_is_read_as_a_value() {
        let materials = LiraParser::materials_from(b"( 3/ 2 3e+006 100 0 50 / 0 S0 3e+006 40 90/ )");
        let Some(Material::Bar { stiffness, .. }) = materials.get(&2) else { panic!("{materials:?}") };
        assert_eq!(*stiffness, Some([3e6, 100., 0., 50.]));
    }
    #[test]
    fn malformed_tokens_and_references_fail_instead_of_shifting_geometry() {
        for content in [
            "(4/0 junk 0 0/)(1/56 1 1/)",
            "(4/NaN 0 0/)(1/56 1 1/)",
            "(4/0 0 0/)(1/56 junk 1 1/)",
            "(4/0 0 0/)(1/56 1 2/)",
            "(4/0 0 0/)(1/42 1 1/)",
        ] {
            assert!(
                LiraParser::parse_bytes(content.as_bytes()).is_err(),
                "{content}"
            );
        }
    }
}
