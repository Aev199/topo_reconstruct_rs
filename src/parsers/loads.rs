//! Loads of a LIRA text file: documents 6 (load rows), 7 (their
//! parameters), 12 (local axes of nodes), 17 (local axes of bars) and the
//! load case names of statement 39 of document 0. The format is that of
//! the LIRA solver input description (appendix A, tables A.14-A.18).
use super::lira::LiraParser;
use fast_float::parse as fast_parse_f64;
use glam::DVec3;
use hashbrown::HashMap;
use rayon::prelude::*;

/// One row of document 6: a load on an element or a node.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LoadRow {
    /// Element or node number (by the load code).
    pub target: u32,
    pub code: u16,
    pub direction: u8,
    /// Row of document 7 holding the parameters.
    pub parameters: u32,
    pub case: u32,
}

#[derive(Debug, Clone, Default)]
pub struct LoadSet {
    /// Load case numbers and names.
    pub cases: Vec<(u32, String)>,
    pub rows: Vec<LoadRow>,
    /// Document 7: row number -> (start, length) in `values`.
    parameters: HashMap<u32, (u32, u32)>,
    values: Vec<f64>,
    /// Document 12: node -> local axes X and Y (unit vectors).
    pub node_axes: HashMap<u32, [DVec3; 2]>,
    /// Document 17: element -> its numbers (the Y vector of a bar, the
    /// angle of a plate).
    pub element_axes: HashMap<u32, Vec<f64>>,
}

impl LoadSet {
    /// Parameters of a row of document 7.
    pub fn parameters(&self, row: u32) -> &[f64] {
        self.parameters
            .get(&row)
            .map_or(&[][..], |&(start, len)| &self.values[start as usize..(start + len) as usize])
    }
}

/// Text of a file's bytes: UTF-8 when it is valid, else Windows-1251 (the
/// LIRA default).
pub fn decode(bytes: &[u8]) -> String {
    if let Ok(text) = std::str::from_utf8(bytes) {
        return text.trim_start_matches('\u{feff}').to_string();
    }
    bytes
        .iter()
        .map(|&b| match b {
            0..=0x7f => b as char,
            0xa8 => '\u{401}',
            0xb8 => '\u{451}',
            0xc0..=0xff => char::from_u32(0x410 + (b as u32 - 0xc0)).unwrap_or('?'),
            _ => '?',
        })
        .collect()
}

fn number(word: &[u8]) -> Option<f64> {
    fast_parse_f64::<f64, _>(word).ok().filter(|v| v.is_finite())
}

fn integer(word: &[u8]) -> Option<u32> {
    std::str::from_utf8(word).ok()?.parse().ok()
}

/// Load case names: statement 39 of document 0, `number: name ;` entries.
fn case_names(content: &[u8]) -> Vec<(u32, String)> {
    let Some(block) = LiraParser::extract_block(content, b"0") else {
        return vec![];
    };
    let text = decode(block);
    let Some(start) = text.find("39;") else {
        return vec![];
    };
    // Up to the next statement ("40;" ...): a number followed by ';' at a
    // line start.
    let rest = &text[start + 3..];
    let mut end = rest.len();
    let mut offset = 0;
    for line in rest.split_inclusive('\n') {
        let t = line.trim_start();
        let digits = t.chars().take_while(char::is_ascii_digit).count();
        if digits > 0 && t[digits..].starts_with(';') && offset > 0 {
            end = offset;
            break;
        }
        offset += line.len();
    }
    rest[..end]
        .split(';')
        .filter_map(|entry| {
            let (number, name) = entry.split_once(':')?;
            Some((number.trim().parse().ok()?, name.trim().trim_end_matches('/').trim().to_string()))
        })
        .collect()
}

/// Parse the loads of LIRA text bytes.
pub fn parse(content: &[u8]) -> LoadSet {
    let mut set = LoadSet {
        cases: case_names(content),
        ..Default::default()
    };
    if let Some(block) = LiraParser::extract_block(content, b"6") {
        let rows: Vec<&[u8]> = block.split(|&b| b == b'/').collect();
        set.rows = rows
            .par_iter()
            .filter_map(|row| {
                let mut words = LiraParser::split_ascii_whitespace_bytes(row);
                let (target, code, direction, parameters, case) = (
                    integer(words.next()?)?,
                    integer(words.next()?)?,
                    integer(words.next()?)?,
                    integer(words.next()?)?,
                    integer(words.next()?)?,
                );
                words.next().is_none().then_some(LoadRow {
                    target,
                    code: u16::try_from(code).ok()?,
                    direction: u8::try_from(direction).ok()?,
                    parameters,
                    case,
                })
            })
            .collect();
    }
    if let Some(block) = LiraParser::extract_block(content, b"7") {
        let rows: Vec<&[u8]> = block.split(|&b| b == b'/').collect();
        let parsed: Vec<(u32, Vec<f64>)> = rows
            .par_iter()
            .filter_map(|row| {
                let mut words = LiraParser::split_ascii_whitespace_bytes(row);
                let id = integer(words.next()?)?;
                let values: Option<Vec<f64>> = words.map(number).collect();
                Some((id, values?))
            })
            .collect();
        for (id, values) in parsed {
            let start = set.values.len() as u32;
            set.parameters.insert(id, (start, values.len() as u32));
            set.values.extend(values);
        }
    }
    let numbered = |id: &[u8]| -> Vec<(u32, Vec<f64>)> {
        let Some(block) = LiraParser::extract_block(content, id) else {
            return vec![];
        };
        block
            .split(|&b| b == b'/')
            .filter_map(|row| {
                let mut words = LiraParser::split_ascii_whitespace_bytes(row);
                let first = integer(words.next()?)?;
                let values: Option<Vec<f64>> = words.map(number).collect();
                Some((first, values?))
            })
            .collect()
    };
    for (node, v) in numbered(b"12") {
        if v.len() >= 6 {
            let x = DVec3::new(v[0], v[1], v[2]);
            let y = DVec3::new(v[3], v[4], v[5]);
            if x.length() > 0. && y.length() > 0. {
                set.node_axes.insert(node, [x.normalize(), y.normalize()]);
            }
        }
    }
    for (element, v) in numbered(b"17") {
        set.element_axes.insert(element, v);
    }
    set
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rows_parameters_cases_and_axes() {
        let text = "( 0/ 1; MODEL/ 2; 5/\n39;\n1: СВОЙ ;\n2: СНЕГ 1.4 ;\n40;\n)\n\
            ( 6/\n * comment /\n 3 16 3 1 1 / 7 6 2 2 2 / 12 0 1 3 2 /\n)\n\
            ( 7/\n 1 -0.5 / 2 0.2 0 0 0 / 3 4.5 /\n)\n\
            ( 12/ 5 1 0 0 0 1 0 / )\n( 17/ 7 0 0 1 / )\n";
        let s = parse(text.as_bytes());
        assert_eq!(s.cases, vec![(1, "СВОЙ".to_string()), (2, "СНЕГ 1.4".to_string())]);
        assert_eq!(s.rows.len(), 3);
        assert_eq!(
            s.rows[0],
            LoadRow { target: 3, code: 16, direction: 3, parameters: 1, case: 1 }
        );
        assert_eq!(s.parameters(1), &[-0.5]);
        assert_eq!(s.parameters(2), &[0.2, 0., 0., 0.]);
        assert_eq!(s.parameters(9), &[] as &[f64]);
        assert_eq!(s.node_axes[&5], [DVec3::X, DVec3::Y]);
        assert_eq!(s.element_axes[&7], vec![0., 0., 1.]);
    }

    #[test]
    fn windows_1251_names_are_decoded() {
        // "Снег" in Windows-1251.
        assert_eq!(decode(&[0xd1, 0xed, 0xe5, 0xe3]), "Снег");
        assert_eq!(decode("Снег".as_bytes()), "Снег");
    }
}
