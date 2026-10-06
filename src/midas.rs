//! The meshed model as a MIDAS Civil `.mxt` text file, in the layout of the
//! Lira_Midas-converter: *NODE, *ELEMENT, *MATERIAL, *SECTION, *THICKNESS,
//! *STLDCASE, *USE-STLD with *CONLOAD, *BEAMLOAD, *PRESSURE. kN and m.
//! Every LIRA load case (except those left out: the self-weight, stages) is
//! a load case of its own. Beam local axes are the program defaults (beta
//! angle 0), as in the converter.
use crate::mesh_loads::MeshLoads;
use crate::meshing::Mesh;
use crate::midas_stiffness::Stiffness;
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Default, Serialize)]
pub struct Report {
    pub nodes: usize,
    pub bars: usize,
    pub plates: usize,
    pub materials: usize,
    pub sections: usize,
    pub thicknesses: usize,
    pub load_cases: usize,
    pub warnings: Vec<String>,
}

/// Cyrillic to Latin as in the converter (MIDAS names are ASCII).
fn transliterate(s: &str) -> String {
    let table = [
        ('а', "a"), ('б', "b"), ('в', "v"), ('г', "g"), ('д', "d"), ('е', "e"), ('ё', "yo"), ('ж', "zh"),
        ('з', "z"), ('и', "i"), ('й', "y"), ('к', "k"), ('л', "l"), ('м', "m"), ('н', "n"), ('о', "o"),
        ('п', "p"), ('р', "r"), ('с', "s"), ('т', "t"), ('у', "u"), ('ф', "f"), ('х', "h"), ('ц', "ts"),
        ('ч', "ch"), ('ш', "sh"), ('щ', "sch"), ('ъ', ""), ('ы', "y"), ('ь', ""), ('э', "e"), ('ю', "yu"), ('я', "ya"),
    ];
    s.chars()
        .map(|c| {
            let lower = c.to_lowercase().next().unwrap_or(c);
            match table.iter().find(|(k, _)| *k == lower) {
                Some((_, t)) if c.is_uppercase() => t.to_uppercase(),
                Some((_, t)) => t.to_string(),
                None => c.to_string(),
            }
        })
        .collect()
}

/// A name MIDAS accepts: only `[0-9A-Za-z_.=<>+-]`, the rest becomes `_`.
pub fn safe_name(raw: &str, fallback: &str) -> String {
    let mut out = String::new();
    for c in transliterate(raw).chars() {
        let ok = c.is_ascii_alphanumeric() || "_.=<>+-".contains(c);
        let c = if ok { c } else { '_' };
        if !(c == '_' && out.ends_with('_')) {
            out.push(c);
        }
    }
    let out = out.trim_matches(|c| c == '_' || c == '.' || c == '-').to_string();
    if out.is_empty() { fallback.to_string() } else { out }
}

/// The file text. `cases` are the load cases to write (number, name); the
/// materials, sections and thicknesses come from `stiffness`
/// (`midas_stiffness::build`) by the stiffness number of the elements.
pub fn write_mxt(mesh: &Mesh, loads: &MeshLoads, stiffness: &Stiffness, cases: &[(u32, String)]) -> (String, Report) {
    let mut report = Report { nodes: mesh.nodes.len(), ..Default::default() };
    let mut out: Vec<String> = vec![];
    out.push("*VERSION\n   8.1.5\n".into());
    out.push("*UNIT    ; Unit System\n   KN   , M, BTU, F\n".into());
    out.push("*PROJINFO    ; Project Information\n   USER=TOPO\n   ADDRESS=MIDAS\n".into());

    // Materials, sections and thicknesses: `midas_stiffness`.
    let (material_of, section_of, thickness_of) = (&stiffness.material_of, &stiffness.section_of, &stiffness.thickness_of);
    let (material_lines, section_lines, thickness_lines) = (&stiffness.material_lines, &stiffness.section_lines, &stiffness.thickness_lines);
    report.materials = material_lines.len();
    report.sections = stiffness.sections;
    report.thicknesses = thickness_lines.len();

    // Nodes.
    let mut text = String::from("*NODE    ; Nodes\n");
    for (i, p) in mesh.nodes.iter().enumerate() {
        text += &format!("{:>6}, {:>12.3}, {:>12.3}, {:>12.3}\n", i + 1, p[0], p[1], p[2]);
    }
    out.push(text);

    // Elements: bars first, then plates.
    let mut text = String::from("*ELEMENT    ; Elements\n");
    let mut bar_id: Vec<Option<usize>> = vec![None; mesh.bars.len()];
    let mut next = 0;
    let mut skipped_bars = 0;
    for (i, b) in mesh.bars.iter().enumerate() {
        let (Some(mat), Some(sect)) = (material_of.get(&(false, b.stiffness)), section_of.get(&b.stiffness)) else {
            skipped_bars += 1;
            continue;
        };
        next += 1;
        bar_id[i] = Some(next);
        text += &format!("{next:>6}, BEAM  , {mat:>4}, {sect:>4}, {:>6}, {:>6}, 0\n", b.nodes[0] + 1, b.nodes[1] + 1);
    }
    report.bars = next;
    if skipped_bars > 0 {
        report.warnings.push(format!("{skipped_bars} bar elements without a material of the LIRA type were not written"));
    }
    let mut plate_id: Vec<Option<usize>> = vec![None; mesh.shells.len()];
    let mut skipped_plates = 0;
    for (i, s) in mesh.shells.iter().enumerate() {
        let (Some(mat), Some(thick)) = (material_of.get(&(true, s.stiffness)), thickness_of.get(&s.stiffness)) else {
            skipped_plates += 1;
            continue;
        };
        next += 1;
        plate_id[i] = Some(next);
        let n: Vec<usize> = s.nodes.iter().map(|n| n + 1).collect();
        if n.len() == 4 {
            text += &format!("{next:>6}, PLATE , {mat:>4}, {thick:>4}, {:>6}, {:>6}, {:>6}, {:>6}, 0\n", n[0], n[1], n[2], n[3]);
        } else {
            text += &format!("{next:>6}, PLATE , {mat:>4}, {thick:>4}, {:>6}, {:>6}, {:>6}, 0, 0\n", n[0], n[1], n[2]);
        }
    }
    report.plates = next - report.bars;
    if skipped_plates > 0 {
        report.warnings.push(format!("{skipped_plates} shell elements without a material of the LIRA type were not written"));
    }
    out.push(text);
    if !material_lines.is_empty() {
        out.push(format!("*MATERIAL    ; Material\n{}\n", material_lines.join("\n")));
    }
    if !section_lines.is_empty() {
        out.push(format!("*SECTION    ; Section\n{}\n", section_lines.join("\n")));
    }
    if !thickness_lines.is_empty() {
        out.push(format!("*THICKNESS    ; Thickness\n{}\n", thickness_lines.join("\n")));
    }

    // Load cases.
    let mut case_names: BTreeMap<u32, String> = BTreeMap::new();
    let mut used: BTreeSet<String> = BTreeSet::new();
    for (case, name) in cases {
        let base = safe_name(name, &format!("Case-{case}"));
        let mut unique_name = base.clone();
        let mut k = 1;
        while !used.insert(unique_name.to_lowercase()) {
            k += 1;
            unique_name = format!("{base}_{k}");
        }
        case_names.insert(*case, unique_name);
    }
    let written: Vec<u32> = case_names
        .keys()
        .copied()
        .filter(|c| {
            loads.pressures.iter().any(|p| p.0 == *c) || loads.bar_loads.iter().any(|b| b.case == *c) || loads.nodal.keys().any(|k| k.0 == *c)
        })
        .collect();
    report.load_cases = written.len();
    if !written.is_empty() {
        let mut text = String::from("*STLDCASE    ; Static Load Cases\n");
        for c in &written {
            text += &format!("   {}, D,\n", case_names[c]);
        }
        out.push(text);
    }
    const GLOBAL: [&str; 3] = ["GX", "GY", "GZ"];
    for case in written {
        let mut text = format!("*USE-STLD, {}\n\n", case_names[&case]);
        let mut nodal: Vec<(usize, [f64; 6])> = loads
            .nodal
            .iter()
            .filter(|(k, f)| k.0 == case && f.iter().any(|x| x.abs() > 1e-12))
            .map(|(k, f)| (k.1, *f))
            .collect();
        nodal.sort_by_key(|x| x.0);
        if !nodal.is_empty() {
            text += "*CONLOAD    ; Nodal Loads\n";
            for (n, f) in nodal {
                text += &format!("{:>6}, {:.3}, {:.3}, {:.3}, {:.3}, {:.3}, {:.3},\n", n + 1, f[0], f[1], f[2], f[3], f[4], f[5]);
            }
            text += "\n";
        }
        let mut beam_lines: Vec<(usize, usize, [f64; 2], [f64; 2])> = vec![];
        for l in &loads.bar_loads {
            if l.case != case {
                continue;
            }
            let Some(id) = bar_id[l.bar] else { continue };
            for k in 0..3 {
                if l.q[0][k].abs() > 1e-12 || l.q[1][k].abs() > 1e-12 {
                    beam_lines.push((id, k, l.t, [l.q[0][k], l.q[1][k]]));
                }
            }
        }
        beam_lines.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)).then(a.2[0].total_cmp(&b.2[0])));
        if !beam_lines.is_empty() {
            text += "*BEAMLOAD    ; Element Beam Loads\n";
            for (id, k, t, q) in beam_lines {
                text += &format!("{id:>6}, BEAM   , UNILOAD, {}, NO , NO, aDir[1], , , , {:.4}, {:.4}, {:.4}, {:.4}, 0, 0, 0, 0, , NO, 0, 0, NO,\n", GLOBAL[k], t[0], q[0], t[1], q[1]);
            }
            text += "\n";
        }
        let mut pressure_lines: Vec<(usize, usize, f64)> = vec![];
        for &(c, shell, p) in &loads.pressures {
            if c != case {
                continue;
            }
            let Some(id) = plate_id[shell] else { continue };
            for (k, v) in p.iter().enumerate() {
                if v.abs() > 1e-12 {
                    pressure_lines.push((id, k, *v));
                }
            }
        }
        pressure_lines.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        if !pressure_lines.is_empty() {
            text += "*PRESSURE    ; Pressure Loads\n";
            for (id, k, v) in pressure_lines {
                text += &format!("{id:>6}, PRES, PLATE, FACE, {}, 0, 0, 0, NO, {v:.4}, 0, 0, 0, 0,\n", GLOBAL[k]);
            }
            text += "\n";
        }
        out.push(text);
    }
    (out.join("\n"), report)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::meshing::{BarPiece, Shell};
    
    fn stiffness() -> Stiffness {
        let text = "( 3/\n1 S0 2.9e+007 50 80/\n 0 RO 1/\n 0 Mu 0.2/\n7 GEI 2.9e+007 0.2 0.2 RO 2.5 /\n)\n";
        let materials = crate::parsers::lira::LiraParser::materials_from(text.as_bytes());
        crate::midas_stiffness::build(&materials, &Default::default(), &[7].into(), &[1].into(), Default::default())
    }

    #[test]
    fn names_are_ascii_and_safe() {
        assert_eq!(safe_name("СНЕГ_1.4|0.5", "x"), "SNEG_1.4_0.5");
        assert_eq!(safe_name("ВЕТЕР_1.4_Y+", "x"), "VETER_1.4_Y+");
        assert_eq!(safe_name("ПОЛЕЗНАЯ НАГРУЗКА, K=1,2", "x"), "POLEZNAYA_NAGRUZKA_K=1_2");
        assert_eq!(safe_name("  ", "Case-3"), "Case-3");
    }

    #[test]
    fn mxt_has_the_converters_sections_and_loads_per_case() {
        let mesh = Mesh {
            nodes: vec![[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.], [0., 0., -3.]],
            vertex_nodes: vec![],
            shells: vec![Shell { nodes: vec![0, 1, 2, 3], surface: 0, stiffness: 7 }],
            bars: vec![BarPiece { nodes: [4, 0], axis: 0, stiffness: 1, t: [0., 1.] }],
        };
        let mut loads = MeshLoads::default();
        loads.pressures.push((1, 0, [0., 0., -5.]));
        loads.bar_loads.push(crate::mesh_loads::BarLoad { case: 2, bar: 0, t: [0., 1.], q: [[0., 2., 0.], [0., 2., 0.]] });
        loads.nodal.insert((2, 1), [1., 0., -10., 0., 0., 0.]);
        let cases = vec![(1, "ПОЛЫ".to_string()), (2, "СНЕГ".to_string()), (3, "ПУСТО".to_string())];
        let (text, report) = write_mxt(&mesh, &loads, &stiffness(), &cases);
        assert_eq!((report.nodes, report.bars, report.plates, report.load_cases), (5, 1, 1, 2));
        assert!(text.contains("*ELEMENT    ; Elements\n     1, BEAM  ,    1,    1,      5,      1, 0\n     2, PLATE ,    2,    2,      1,      2,      3,      4, 0\n"), "{text}");
        assert!(text.contains("*STLDCASE    ; Static Load Cases\n   POLY, D,\n   SNEG, D,\n"), "{text}");
        assert!(text.contains("*USE-STLD, POLY\n\n*PRESSURE    ; Pressure Loads\n     2, PRES, PLATE, FACE, GZ, 0, 0, 0, NO, -5.0000, 0, 0, 0, 0,\n"), "{text}");
        assert!(text.contains("*CONLOAD    ; Nodal Loads\n     2, 1.000, 0.000, -10.000, 0.000, 0.000, 0.000,\n"), "{text}");
        assert!(text.contains("     1, BEAM   , UNILOAD, GY, NO , NO, aDir[1], , , , 0.0000, 2.0000, 1.0000, 2.0000, 0, 0, 0, 0, , NO, 0, 0, NO,"), "{text}");
        assert!(text.contains("*THICKNESS    ; Thickness\n    2, VALUE, YES, 0.200, 0, NO, 0, 0, Plate_0.2"), "{text}");
        assert!(!text.contains("PUSTO"));
    }
}
