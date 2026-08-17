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
        let content = &mmap[..];

        // 1. Поиск блоков ( 4/ ... ) для узлов и ( 1/ ... ) для КЭ
        let block_4_data = Self::extract_block(content, b"4")
            .ok_or_else(|| Error::new(ErrorKind::InvalidData, "Блок координат узлов (4/) не найден"))?;

        let block_1_data = Self::extract_block(content, b"1")
            .ok_or_else(|| Error::new(ErrorKind::InvalidData, "Блок элементов (1/) не найден"))?;

        // 2. Параллельный парсинг узлов (Rayon)
        let raw_node_chunks: Vec<&[u8]> = block_4_data
            .split(|&b| b == b'/')
            .filter(|chunk| !Self::is_empty_or_ws(chunk))
            .collect();

        let parsed_nodes: Vec<(u32, DVec3)> = raw_node_chunks
            .par_iter()
            .enumerate()
            .filter_map(|(idx, chunk)| {
                let mut coords = [0.0f64; 3];
                let mut count = 0;

                // Читаем 3 координаты (x, y, z)
                for word in Self::split_ascii_whitespace_bytes(chunk) {
                    if count < 3 {
                        if let Ok(val) = fast_parse_f64::<f64, _>(word) {
                            coords[count] = val;
                            count += 1;
                        }
                    } else {
                        break;
                    }
                }

                if count >= 3 {
                    let node_id = (idx + 1) as u32;
                    Some((node_id, DVec3::new(coords[0], coords[1], coords[2])))
                } else {
                    None
                }
            })
            .collect();

        let mut nodes_map = HashMap::with_capacity(parsed_nodes.len());
        for (id, pt) in parsed_nodes {
            nodes_map.insert(id, pt);
        }

        // 3. Параллельный парсинг элементов (Rayon)
        let raw_elem_chunks: Vec<&[u8]> = block_1_data
            .split(|&b| b == b'/')
            .filter(|chunk| !Self::is_empty_or_ws(chunk))
            .collect();

        let elements: Vec<ElementData> = raw_elem_chunks
            .par_iter()
            .enumerate()
            .filter_map(|(idx, chunk)| {
                let mut ints = Vec::with_capacity(8);

                for word in Self::split_ascii_whitespace_bytes(chunk) {
                    // Быстрый парсинг u32
                    if let Ok(s) = std::str::from_utf8(word) {
                        if let Ok(val) = s.parse::<u32>() {
                            ints.push(val);
                        }
                    }
                }

                if ints.len() >= 4 {
                    let elem_id = (idx + 1) as u32;
                    let elem_type = ints[0];
                    let stiff_id = ints[1];
                    let elem_nodes = ints[2..].to_vec();

                    Some(ElementData {
                        id: elem_id,
                        elem_type,
                        stiff_id,
                        nodes: elem_nodes,
                    })
                } else {
                    None
                }
            })
            .collect();

        Ok(MeshData {
            nodes: nodes_map,
            elements,
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
                while j < len && (content[j] == b' ' || content[j] == b'\t' || content[j] == b'\r' || content[j] == b'\n') {
                    j += 1;
                }

                // Проверяем ID блока
                if j + block_id.len() < len && &content[j..j + block_id.len()] == block_id {
                    let mut k = j + block_id.len();
                    // Пропускаем пробелы перед '/'
                    while k < len && (content[k] == b' ' || content[k] == b'\t' || content[k] == b'\r' || content[k] == b'\n') {
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
        slice.iter().all(|&b| b == b' ' || b == b'\t' || b == b'\r' || b == b'\n')
    }

    /// Итератор по непустым словам (whitespace-separated bytes)
    fn split_ascii_whitespace_bytes(slice: &[u8]) -> impl Iterator<Item = &[u8]> {
        let mut i = 0;
        let len = slice.len();

        std::iter::from_fn(move || {
            // Пропуск начальных пробелов
            while i < len && (slice[i] == b' ' || slice[i] == b'\t' || slice[i] == b'\r' || slice[i] == b'\n') {
                i += 1;
            }
            if i >= len {
                return None;
            }
            let start = i;
            // Поиск конца слова
            while i < len && !(slice[i] == b' ' || slice[i] == b'\t' || slice[i] == b'\r' || slice[i] == b'\n') {
                i += 1;
            }
            Some(&slice[start..i])
        })
    }
}