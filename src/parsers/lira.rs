use crate::models::{ElementData, MeshData};
use fast_float::parse as fast_parse_f64;
use glam::DVec3;
use hashbrown::HashMap;
use memmap2::Mmap;
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, Error, ErrorKind};
use std::path::Path;

pub struct LiraParser;

impl LiraParser {
    /// Потоковый параллельный парсинг текстового файла ЛИРА (.txt)
    pub fn parse<P: AsRef<Path>>(filepath: P) -> io::Result<MeshData> {
        let file = File::open(filepath)?;
        let mmap = unsafe { Mmap::map(&file)? };
        Self::parse_bytes(&mmap)
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

    /// Быстрый поиск содержимого блока `( <id>/ ... )` без аллокаций строк
    fn extract_block<'a>(content: &'a [u8], block_id: &[u8]) -> Option<&'a [u8]> {
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
    fn is_empty_or_ws(slice: &[u8]) -> bool {
        slice
            .iter()
            .all(|&b| b == b' ' || b == b'\t' || b == b'\r' || b == b'\n')
    }

    /// Итератор по непустым словам (whitespace-separated bytes)
    fn split_ascii_whitespace_bytes(slice: &[u8]) -> impl Iterator<Item = &[u8]> {
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
