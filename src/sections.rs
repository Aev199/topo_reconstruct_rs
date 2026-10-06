//! Cross sections of the bars of a LIRA model as MIDAS Civil `DBUSER` shapes,
//! after the Lira_Midas-converter: parametric LIRA sections (S1, S2, S3, S5,
//! S6), named profiles of the block `{13/` looked up in the LIRA catalogue
//! (`data/sortament.tsv`: name, kind, h, b, s, t in mm) and builtup I sections.
//! Dimensions in metres.
use serde::Serialize;

/// Shape of a bar section.
#[derive(Debug, Clone, PartialEq, Serialize)]
#[serde(tag = "shape", rename_all = "snake_case")]
pub enum Shape {
    /// Solid rectangle (MIDAS SB): height `h`, width `b`.
    Rect { h: f64, b: f64 },
    /// Round tube (P): outer diameter `d`, wall `t`.
    Pipe { d: f64, t: f64 },
    /// Solid round bar (SR).
    Round { d: f64 },
    /// Symmetric I section of LIRA S1 (H): the flanges are equal.
    ISym { h: f64, b: f64, tw: f64, tf: f64 },
    /// T section (T).
    T { h: f64, b: f64, tw: f64, tf: f64 },
    /// I section (H): height, top flange width, web, top flange, bottom flange width and thickness.
    I { h: f64, bt: f64, tw: f64, tft: f64, bb: f64, tfb: f64 },
    /// Box of LIRA S5 (B): sides `t1`, flanges `t2`.
    Box { h: f64, b: f64, t1: f64, t2: f64 },
    /// Rectangular or square tube with equal walls (B).
    Tube { h: f64, b: f64, t: f64 },
    /// No dimensions: only the stiffness values (EF, EIy, EIz, GIk) of the type.
    Explicit,
}

/// A named profile of the block `{13/`.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Profile {
    /// `DoubleT`, `Tubing`, `Pipe`, `Round`, `DoubleTBuiltup`, ...
    pub section: String,
    /// The designation: `18`, `100 x 5`, `Flange=100 x 22; Wall=200 x 10`.
    pub shape: String,
    pub builtup: bool,
}

/// `Pipe` or `Round` for a LIRA S6 row: outer and inner diameters (cm, as the
/// dimensions of every LIRA section row; the converter takes values up to 1 as metres).
pub fn s6(d_outer: f64, d_inner: f64) -> Shape {
    let (d, inner) = (d_outer / 100., d_inner / 100.);
    let tw = (d - inner) / 2.;
    if inner > 0. && tw > 1e-4 && d > inner {
        Shape::Pipe { d, t: tw }
    } else {
        Shape::Round { d }
    }
}

impl Shape {
    /// Area (m2) of the nominal shape; `None` without dimensions.
    pub fn area(&self) -> Option<f64> {
        Some(match *self {
            Shape::Rect { h, b } => h * b,
            Shape::Pipe { d, t } => std::f64::consts::PI / 4. * (d * d - (d - 2. * t).powi(2)),
            Shape::Round { d } => std::f64::consts::PI / 4. * d * d,
            Shape::ISym { h, b, tw, tf } => 2. * b * tf + (h - 2. * tf) * tw,
            Shape::T { h, b, tw, tf } => b * tf + (h - tf) * tw,
            Shape::I { h, bt, tw, tft, bb, tfb } => bt * tft + bb * tfb + (h - tft - tfb) * tw,
            Shape::Box { h, b, t1, t2 } => h * b - (h - 2. * t2).max(0.) * (b - 2. * t1).max(0.),
            Shape::Tube { h, b, t } => h * b - (h - 2. * t).max(0.) * (b - 2. * t).max(0.),
            Shape::Explicit => return None,
        })
    }

    /// Second moments (m4) about the horizontal axis (the bending of the height `h`, the strong
    /// one) and the vertical axis, for the section standing as drawn; `None` without dimensions.
    pub fn inertia(&self) -> Option<(f64, f64)> {
        let pi = std::f64::consts::PI;
        Some(match *self {
            Shape::Rect { h, b } => (b * h.powi(3) / 12., h * b.powi(3) / 12.),
            Shape::Round { d } => (pi * d.powi(4) / 64., pi * d.powi(4) / 64.),
            Shape::Pipe { d, t } => (pi * (d.powi(4) - (d - 2. * t).powi(4)) / 64., pi * (d.powi(4) - (d - 2. * t).powi(4)) / 64.),
            Shape::Tube { h, b, t } => ((b * h.powi(3) - (b - 2. * t) * (h - 2. * t).powi(3)) / 12., (h * b.powi(3) - (h - 2. * t) * (b - 2. * t).powi(3)) / 12.),
            Shape::Box { h, b, t1, t2 } => ((b * h.powi(3) - (b - 2. * t1) * (h - 2. * t2).powi(3)) / 12., (h * b.powi(3) - (h - 2. * t2) * (b - 2. * t1).powi(3)) / 12.),
            Shape::ISym { h, b, tw, tf } => return Shape::I { h, bt: b, tw, tft: tf, bb: b, tfb: tf }.inertia(),
            Shape::I { h, bt, tw, tft, bb, tfb } => {
                // Parts from the bottom: bottom flange, web, top flange.
                let hw = h - tft - tfb;
                let parts = [(bb * tfb, tfb / 2., bb * tfb.powi(3) / 12.), (tw * hw, tfb + hw / 2., tw * hw.powi(3) / 12.), (bt * tft, h - tft / 2., bt * tft.powi(3) / 12.)];
                let a: f64 = parts.iter().map(|p| p.0).sum();
                let yc = parts.iter().map(|p| p.0 * p.1).sum::<f64>() / a;
                let ix = parts.iter().map(|p| p.2 + p.0 * (p.1 - yc).powi(2)).sum();
                (ix, tft * bt.powi(3) / 12. + tfb * bb.powi(3) / 12. + hw * tw.powi(3) / 12.)
            }
            Shape::T { h, b, tw, tf } => {
                let hw = h - tf;
                let parts = [(tw * hw, hw / 2., tw * hw.powi(3) / 12.), (b * tf, h - tf / 2., b * tf.powi(3) / 12.)];
                let a: f64 = parts.iter().map(|p| p.0).sum();
                let yc = parts.iter().map(|p| p.0 * p.1).sum::<f64>() / a;
                (parts.iter().map(|p| p.2 + p.0 * (p.1 - yc).powi(2)).sum(), tf * b.powi(3) / 12. + hw * tw.powi(3) / 12.)
            }
            Shape::Explicit => return None,
        })
    }

    /// Torsion constant (m4): exact for the solid round, the ring and the rectangle (series), Bredt's for the
    /// closed tubes and boxes (mean line), the thin-wall sum for the open I and T sections.
    pub fn torsion(&self) -> Option<f64> {
        let pi = std::f64::consts::PI;
        Some(match *self {
            Shape::Rect { h, b } => {
                let (b, h) = (b.min(h), b.max(h));
                h * b.powi(3) * (1. / 3. - 0.21 * (b / h) * (1. - b.powi(4) / (12. * h.powi(4))))
            }
            Shape::Round { d } => pi * d.powi(4) / 32.,
            Shape::Pipe { d, t } => pi * (d.powi(4) - (d - 2. * t).powi(4)) / 32.,
            Shape::Tube { h, b, t } => {
                let (bm, hm) = (b - t, h - t);
                4. * (bm * hm).powi(2) * t / (2. * (bm + hm))
            }
            Shape::Box { h, b, t1, t2 } => {
                // Webs t1 (height h), flanges t2 (width b), along the mean line.
                let (bm, hm) = (b - t1, h - t2);
                4. * (bm * hm).powi(2) / (2. * bm / t2 + 2. * hm / t1)
            }
            Shape::ISym { h, b, tw, tf } => (2. * b * tf.powi(3) + (h - 2. * tf) * tw.powi(3)) / 3.,
            Shape::I { h, bt, tw, tft, bb, tfb } => (bt * tft.powi(3) + bb * tfb.powi(3) + (h - tft - tfb) * tw.powi(3)) / 3.,
            Shape::T { h, b, tw, tf } => (b * tf.powi(3) + (h - tf) * tw.powi(3)) / 3.,
            Shape::Explicit => return None,
        })
    }

    /// Height and width (m) of the bounding box of the section.
    pub fn extent(&self) -> Option<(f64, f64)> {
        Some(match *self {
            Shape::Rect { h, b } | Shape::ISym { h, b, .. } | Shape::T { h, b, .. } | Shape::Box { h, b, .. } | Shape::Tube { h, b, .. } => (h, b),
            Shape::Round { d } | Shape::Pipe { d, .. } => (d, d),
            Shape::I { h, bt, bb, .. } => (h, bt.max(bb)),
            Shape::Explicit => return None,
        })
    }

    /// The part of the converter's section name that depends on the shape
    /// (`Rect_0.2_X_0.3`, `Pipe_0.2_x_0.01`, `Box_0.2_X_0.2`, `I-Shape_0.15_X_0.15`).
    fn plain_name(&self) -> String {
        let name = match *self {
            Shape::Rect { h, b } => format!("Rect_{b:.3}_X_{h:.3}"),
            Shape::Pipe { d, t } => format!("Pipe_{d:.3}_x_{t:.3}"),
            Shape::Round { d } => format!("Rod_D{d:.3}"),
            Shape::ISym { h, b, .. } | Shape::I { h, bt: b, .. } => format!("I-Shape_{h:.3}_X_{b:.3}"),
            Shape::T { h, b, .. } => format!("T-Shape_{h:.3}_X_{b:.3}"),
            Shape::Box { h, b, .. } => format!("Box_{:.3}_X_{:.3}", h.max(b), h.min(b)),
            Shape::Tube { h, b, .. } => format!("Tube_{:.3}_X_{:.3}", h.max(b), h.min(b)),
            Shape::Explicit => return String::new(),
        };
        name.trim_end_matches('0').trim_end_matches('.').to_string()
    }

    /// The part of the `*SECTION` line of a `DBUSER` section after `YES, `
    /// (the converter's layout): `SB, 2, 0.4000, 0.2000, 0, ...`.
    pub fn dbuser_body(&self) -> Option<String> {
        Some(match *self {
            Shape::Rect { h, b } => format!("SB, 2, {h:.4}, {b:.4}, 0, 0, 0, 0, 0, 0, 0, 0"),
            Shape::Pipe { d, t } => format!("P, 2, {d:.4}, {t:.4}, 0, 0, 0, 0, 0, 0, 0, 0"),
            // "SR, 2" as the converter's corrected writer has it (its base writer wrote "SR, 1").
            Shape::Round { d } => format!("SR, 2, {d:.4}, 0, 0, 0, 0, 0, 0, 0, 0, 0"),
            Shape::ISym { h, b, tw, tf } => format!("H, 2, {h:.4}, {b:.4}, {tw:.4}, {tf:.4}, {b:.4}, {tf:.4}, 0, 0, 0, 0"),
            Shape::T { h, b, tw, tf } => format!("T, 2, {h:.4}, {b:.4}, {tw:.4}, {tf:.4}, 0, 0, 0, 0, 0, 0"),
            Shape::I { h, bt, tw, tft, bb, tfb } => format!("H, 2, {h:.4}, {bt:.4}, {tw:.4}, {tft:.4}, {bb:.4}, {tfb:.4}, 0, 0, 0, 0"),
            Shape::Box { h, b, t1, t2 } => format!("B, 2, {:.4}, {:.4}, {t1:.4}, {t1:.4}, {t2:.4}, {t2:.4}, 0, 0, 0, 0", h.max(b), h.min(b)),
            Shape::Tube { h, b, t } => format!("B, 2, {h:.4}, {b:.4}, {t:.4}, {t:.4}, {t:.4}, {t:.4}, 0, 0, 0, 0"),
            Shape::Explicit => return None,
        })
    }

    /// The whole `*SECTION` line of a `DBUSER` section.
    pub fn dbuser(&self, id: usize, name: &str) -> Option<String> {
        Some(format!("{id:>5}, DBUSER, {name}, CC, 0, 0, 0, 0, 0, 0, YES, {}", self.dbuser_body()?))
    }

    /// The name of a section of this shape from a LIRA stiffness type
    /// (a profile of the block 13 names it by its designation).
    pub fn name(&self, id: u32, profile: Option<&Profile>) -> String {
        if let Some(p) = profile {
            let shape = p.shape.replace(" x ", "x");
            return match p.section.as_str() {
                "DoubleT" => format!("DoubleT_{}", p.shape),
                "Tubing" => format!("Tube_{shape}"),
                "Pipe" => format!("Pipe_{shape}"),
                "DoubleTBuiltup" => p.shape.replace("; ", "_").replace('=', ""),
                other => format!("{other}_{}", p.shape),
            };
        }
        let name = self.plain_name();
        if name.is_empty() { format!("Section-{id}") } else { name }
    }
}

/// Dimensions in mm of `a x b`: `100 x 5` -> `[100, 5]`.
fn dimensions(shape: &str) -> Option<Vec<f64>> {
    crate::sortament::ascii_x(shape).split('x').map(|p| p.trim().parse::<f64>().ok()).collect()
}

/// The shape of a named profile of the block 13; `None` when it cannot be
/// found in the catalogue or read from its designation.
pub fn profile_shape(profile: &Profile) -> Option<Shape> {
    if profile.builtup || profile.section == "DoubleTBuiltup" {
        // "Flange=100 x 22; Wall=200 x 10": the flange width x thickness, the wall height x thickness.
        let part = |key: &str| {
            let at = profile.shape.find(key)? + key.len();
            let text = profile.shape[at..].split(';').next()?;
            dimensions(text)
        };
        let (flange, wall) = (part("Flange=")?, part("Wall=")?);
        if flange.len() < 2 || wall.len() < 2 {
            return None;
        }
        let (bf, tf, hw, tw) = (flange[0], flange[1], wall[0], wall[1]);
        return Some(Shape::I { h: ((hw + 2. * tf) / 1000.), bt: bf / 1000., tw: tw / 1000., tft: tf / 1000., bb: bf / 1000., tfb: tf / 1000. });
    }
    let section = profile.section.to_lowercase();
    if section == "round" {
        return dimensions(&profile.shape).and_then(|d| d.first().copied()).map(|d| Shape::Round { d: d / 1000. });
    }
    if let Some(e) = crate::sortament::lookup(&profile.section, &profile.shape) {
        let h = e.h? / 1000.;
        if e.kind == "Труба круглая" || section.contains("pipe") {
            let t = e.t.or(e.s).unwrap_or(5.) / 1000.;
            return Some(Shape::Pipe { d: h, t });
        }
        if e.kind.contains("квадрат") || e.kind.contains("прямоугол") || section.contains("tub") {
            let t = e.t.or(e.s).unwrap_or(5.) / 1000.;
            return Some(Shape::Tube { h, b: e.b.map_or(h, |b| b / 1000.), t });
        }
        if let Some(b) = e.b {
            let tw = e.s.map_or(h * 0.03, |s| s / 1000.);
            let tf = e.t.map_or(h * 0.06, |t| t / 1000.);
            return Some(Shape::I { h, bt: b / 1000., tw, tft: tf, bb: b / 1000., tfb: tf });
        }
    }
    // Not in the catalogue: the dimensions of the designation.
    let d = dimensions(&profile.shape)?;
    if section.contains("tub") && d.len() == 2 {
        return Some(Shape::Tube { h: d[0] / 1000., b: d[0] / 1000., t: d[1] / 1000. });
    }
    if section.contains("tub") && d.len() == 3 {
        return Some(Shape::Tube { h: d[0] / 1000., b: d[1] / 1000., t: d[2] / 1000. });
    }
    if section.contains("pipe") && d.len() == 2 {
        return Some(Shape::Pipe { d: d[0] / 1000., t: d[1] / 1000. });
    }
    None
}

/// Section properties of a profile from the catalogue (m2, m4, kg/m).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Properties {
    pub area: f64,
    /// Second moment about the strong axis (parallel to the flanges) and the weak one.
    pub ix: f64,
    pub iy: f64,
    pub mass: f64,
}

/// The catalogue properties of a named profile: exact area and moments of inertia
/// (with the fillets), not those of the simplified shape of `profile_shape`.
pub fn profile_properties(profile: &Profile) -> Option<Properties> {
    if profile.builtup || profile.section == "DoubleTBuiltup" || profile.section.eq_ignore_ascii_case("round") {
        return None;
    }
    let e = crate::sortament::lookup(&profile.section, &profile.shape)?;
    Some(Properties { area: e.area? / 1e4, ix: e.ix? / 1e8, iy: e.iy? / 1e8, mass: e.mass.unwrap_or(0.) })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn profile(section: &str, shape: &str) -> Profile {
        Profile { section: section.into(), shape: shape.into(), builtup: shape.contains("Flange=") }
    }

    #[test]
    fn the_sections_of_the_converters_test_model() {
        // test5.txt of the converter: S3, S5, S6 rows and the profiles of the block 13.
        let s3 = Shape::I { h: 0.15, bt: 0.15, tw: 0.03, tft: 0.02, bb: 0.10, tfb: 0.02 };
        assert_eq!(s3.dbuser(3, "I-Shape_0.150_X_0.15").unwrap(), "    3, DBUSER, I-Shape_0.150_X_0.15, CC, 0, 0, 0, 0, 0, 0, YES, H, 2, 0.1500, 0.1500, 0.0300, 0.0200, 0.1000, 0.0200, 0, 0, 0, 0");
        assert_eq!(s3.name(4, None), "I-Shape_0.150_X_0.15");
        let s5 = Shape::Box { h: 0.2, b: 0.2, t1: 0.02, t2: 0.03 };
        assert_eq!(s5.dbuser(7, "Box").unwrap(), "    7, DBUSER, Box, CC, 0, 0, 0, 0, 0, 0, YES, B, 2, 0.2000, 0.2000, 0.0200, 0.0200, 0.0300, 0.0300, 0, 0, 0, 0");
        let pipe = s6(20., 18.);
        let Shape::Pipe { d, t } = pipe else { panic!("{pipe:?}") };
        assert!((d - 0.2).abs() < 1e-12 && (t - 0.01).abs() < 1e-12);
        assert_eq!(pipe.dbuser(8, "P").unwrap(), "    8, DBUSER, P, CC, 0, 0, 0, 0, 0, 0, YES, P, 2, 0.2000, 0.0100, 0, 0, 0, 0, 0, 0, 0, 0");
        assert_eq!(pipe.name(9, None), "Pipe_0.200_x_0.01");
        // A solid round bar (diameter 28 cm, no inner one).
        assert_eq!(s6(28., 0.), Shape::Round { d: 0.28 });
        // DoubleT 18: H 180, B 90, tw 5.1, tf 8.1 (converter output of test5).
        let t = profile_shape(&profile("DoubleT", "18")).unwrap();
        let Shape::I { h, bt, tw, tft, .. } = t else { panic!("{t:?}") };
        assert!((h - 0.18).abs() < 1e-9 && (bt - 0.09).abs() < 1e-9 && (tw - 0.0051).abs() < 1e-9 && (tft - 0.0081).abs() < 1e-9, "{t:?}");
        // Tubing 100 x 5 -> box 100 x 100 x 5; Pipe 180 x 36 is not in the catalogue: its designation.
        assert_eq!(profile_shape(&profile("Tubing", "100 x 5")), Some(Shape::Tube { h: 0.1, b: 0.1, t: 0.005 }));
        assert_eq!(profile_shape(&profile("Pipe", "180 x 36")), Some(Shape::Pipe { d: 0.18, t: 0.036 }));
        // Builtup: flange 100 x 22, wall 200 x 10 -> H = 200 + 2 x 22.
        let b = profile_shape(&profile("DoubleTBuiltup", "Flange=100 x 22; Wall=200 x 10")).unwrap();
        assert_eq!(b, Shape::I { h: 0.244, bt: 0.1, tw: 0.01, tft: 0.022, bb: 0.1, tfb: 0.022 });
        assert_eq!(b.name(7, Some(&profile("DoubleTBuiltup", "Flange=100 x 22; Wall=200 x 10"))), "Flange100 x 22_Wall200 x 10");
    }

    #[test]
    fn profiles_of_the_island_model_are_found() {
        // Cyrillic designations of GOST profiles and a round bar.
        assert!(matches!(profile_shape(&profile("DoubleT", "20Ш1")), Some(Shape::I { .. })));
        assert!(matches!(profile_shape(&profile("DoubleT", "20К2")), Some(Shape::I { .. })));
        assert_eq!(profile_shape(&profile("Round", "100")), Some(Shape::Round { d: 0.1 }));
        assert!(matches!(profile_shape(&profile("Pipe", "273 x 20")), Some(Shape::Pipe { .. })));
        assert!(matches!(profile_shape(&profile("Pipe", "244.5х8")), Some(Shape::Pipe { .. })));
        assert!(matches!(profile_shape(&profile("Tubing", "120 x 80 x 5")), Some(Shape::Tube { .. })));
    }

    #[test]
    fn moments_of_inertia_of_the_shapes() {
        // I 18Б1 by its simplified dimensions (no fillets): within 8 % of the catalogue's 1063 cm4 and 81.9 cm4.
        let i = Shape::I { h: 0.177, bt: 0.091, tw: 0.0043, tft: 0.0065, bb: 0.091, tfb: 0.0065 };
        let (ix, iy) = i.inertia().unwrap();
        assert!((ix * 1e8 / 1063. - 1.).abs() < 0.08 && (iy * 1e8 / 81.9 - 1.).abs() < 0.12, "{} {}", ix * 1e8, iy * 1e8);
        let (rx, ry) = Shape::Rect { h: 0.5, b: 0.3 }.inertia().unwrap();
        assert!((rx - 0.3 * 0.125 / 12.).abs() < 1e-12 && (ry - 0.5 * 0.027 / 12.).abs() < 1e-12);
        let (t, _) = Shape::Tube { h: 0.1, b: 0.1, t: 0.005 }.inertia().unwrap();
        // 100 x 100 x 5 without fillets: (0.1^4 - 0.09^4) / 12 = 286.6 cm4.
        assert!((t * 1e8 - 286.58).abs() < 0.01, "{}", t * 1e8);
    }

    #[test]
    fn areas_of_the_shapes() {
        let a = Shape::Rect { h: 0.2, b: 0.3 }.area().unwrap();
        assert!((a - 0.06).abs() < 1e-12);
        let p = Shape::Pipe { d: 0.2, t: 0.01 }.area().unwrap();
        assert!((p - std::f64::consts::PI / 4. * (0.04 - 0.0324)).abs() < 1e-12);
        assert!(Shape::Explicit.area().is_none());
    }
}
