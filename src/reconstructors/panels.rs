use crate::config::ReconstructionConfig;
use crate::geometry::utils::{clean_polygon_coords_3d, get_plane_basis};
use crate::models::{ElementData, MacroPanel, MeshData, PanelType};
use geo::{BooleanOps, Coord, LineString, MultiPolygon, Polygon, Simplify};
use glam::DVec3;
use hashbrown::HashMap;
use rayon::prelude::*;

pub struct PanelReconstructor;

struct ParsedShell {
    id: u32,
    stiff_id: u32,
    normal: DVec3,
    d: f64,
    centroid: DVec3,
    pts: Vec<DVec3>,
}

struct PlaneCluster {
    is_horiz: bool,
    is_vert: bool,
    z: Option<f64>,
    normal: DVec3,
    d: f64,
    elements: Vec<ParsedShell>,
}

impl PanelReconstructor {
    /// Параллельное восстановление плит и стен
    pub fn reconstruct(
        mesh_data: &MeshData,
        canonical_nodes: &HashMap<u32, u32>,
        config: &ReconstructionConfig,
    ) -> Vec<MacroPanel> {
        let shell_elems: Vec<&ElementData> = mesh_data
            .elements
            .iter()
            .filter(|e| e.nodes.len() >= 3)
            .collect();

        let parsed_shells: Vec<ParsedShell> = shell_elems
            .into_iter()
            .filter_map(|el| {
                let pts: Vec<DVec3> = el
                    .nodes
                    .iter()
                    .filter_map(|nid| canonical_nodes.get(nid).and_then(|cid| mesh_data.nodes.get(cid).copied()))
                    .collect();

                if pts.len() < 3 {
                    return None;
                }

                let centroid = pts.iter().copied().sum::<DVec3>() / (pts.len() as f64);
                let v1 = pts[1] - pts[0];
                let v2 = pts[2] - pts[0];
                let mut norm = v1.cross(v2);
                let mut n_len = norm.length();

                if n_len < 1e-7 && pts.len() >= 4 {
                    let v2_alt = pts[3] - pts[0];
                    norm = v1.cross(v2_alt);
                    n_len = norm.length();
                }

                if n_len < 1e-7 {
                    return None;
                }

                norm /= n_len;

                // Каноническая ориентация нормали
                if norm.x.abs() > 1e-4 {
                    if norm.x < 0.0 { norm = -norm; }
                } else if norm.y.abs() > 1e-4 {
                    if norm.y < 0.0 { norm = -norm; }
                } else if norm.z < 0.0 {
                    norm = -norm;
                }

                let d = -norm.dot(centroid);

                Some(ParsedShell {
                    id: el.id,
                    stiff_id: el.stiff_id,
                    normal: norm,
                    d,
                    centroid,
                    pts,
                })
            })
            .collect();

        // Кластеризация элементов по плоскостям
        let mut clusters: Vec<PlaneCluster> = Vec::new();
        for el in parsed_shells {
            let nz = el.normal.z;
            let is_horiz = nz.abs() > 0.85;
            let is_vert = nz.abs() < 0.15;

            let mut matched = false;
            for cl in &mut clusters {
                if is_horiz && cl.is_horiz {
                    if let Some(cl_z) = cl.z {
                        if (el.centroid.z - cl_z).abs() < config.tol_dist {
                            cl.elements.push(el);
                            matched = true;
                            break;
                        }
                    }
                } else if is_vert && cl.is_vert {
                    let n2d_el = el.normal.truncate().normalize();
                    let n2d_cl = cl.normal.truncate().normalize();
                    if n2d_el.dot(n2d_cl) > 0.98 && (el.d - cl.d).abs() < config.tol_dist {
                        cl.elements.push(el);
                        matched = true;
                        break;
                    }
                } else if !is_horiz && !is_vert && !cl.is_horiz && !cl.is_vert {
                    if el.normal.dot(cl.normal) > 0.98 && (el.d - cl.d).abs() < config.tol_dist {
                        cl.elements.push(el);
                        matched = true;
                        break;
                    }
                }
            }

            if !matched {
                clusters.push(PlaneCluster {
                    is_horiz,
                    is_vert,
                    z: if is_horiz { Some(el.centroid.z) } else { None },
                    normal: el.normal,
                    d: el.d,
                    elements: vec![el],
                });
            }
        }

        // Параллельное булево объединение через Rayon
        let panels_results: Vec<Vec<MacroPanel>> = clusters
            .into_par_iter()
            .enumerate()
            .map(|(cl_idx, cl)| {
                let norm = cl.normal;
                let orig = cl.elements[0].centroid;
                let (u_axis, v_axis) = get_plane_basis(norm);

                let mut geo_polygons = Vec::with_capacity(cl.elements.len());

                for el in &cl.elements {
                    let mut pts_2d: Vec<(f64, f64, f64)> = el
                        .pts
                        .iter()
                        .map(|&p| {
                            let rel = p - orig;
                            let u = rel.dot(u_axis);
                            let v = rel.dot(v_axis);
                            let angle = v.atan2(u);
                            (u, v, angle)
                        })
                        .collect();

                    // Сортировка по углу (CCW) против диагональных крестов
                    pts_2d.sort_by(|a, b| a.2.total_cmp(&b.2));

                    let mut coords: Vec<Coord<f64>> = pts_2d
                        .iter()
                        .map(|&(u, v, _)| Coord { x: u, y: v })
                        .collect();

                    if coords.len() >= 3 {
                        coords.push(coords[0]); // Замыкание
                        let poly = Polygon::new(LineString::new(coords), vec![]);
                        geo_polygons.push(poly);
                    }
                }

                if geo_polygons.is_empty() {
                    return Vec::new();
                }

                // Каскадное булево объединение
                let mut merged: MultiPolygon<f64> = MultiPolygon::new(vec![]);
                for poly in geo_polygons {
                    merged = merged.union(&poly);
                }

                let simplified = if config.simplify_tol > 0.0 {
                    merged.simplify(&config.simplify_tol)
                } else {
                    merged
                };

                let panel_type = if cl.is_horiz {
                    PanelType::Slab
                } else if cl.is_vert {
                    PanelType::Wall
                } else {
                    PanelType::InclinedPanel
                };

                let mut local_panels = Vec::new();
                let mut sub_id = 1;

                for geom in simplified.0 {
                    let mut panel_polygons_3d = Vec::new();

                    // Внешний контур
                    let ext_coords: Vec<Coord<f64>> = geom.exterior().coords().copied().collect();
                    let raw_ext: Vec<DVec3> = ext_coords
                        .iter()
                        .take(ext_coords.len().saturating_sub(1))
                        .map(|c| orig + c.x * u_axis + c.y * v_axis)
                        .collect();

                    let clean_ext = clean_polygon_coords_3d(&raw_ext, 0.9995);
                    if clean_ext.len() >= 3 {
                        panel_polygons_3d.push(clean_ext.into_iter().map(|p| [p.x, p.y, p.z]).collect());
                    }

                    // Внутренние проемы
                    for hole in geom.interiors() {
                        let hole_coords: Vec<Coord<f64>> = hole.coords().copied().collect();
                        let raw_hole: Vec<DVec3> = hole_coords
                            .iter()
                            .take(hole_coords.len().saturating_sub(1))
                            .map(|c| orig + c.x * u_axis + c.y * v_axis)
                            .collect();

                        let clean_hole = clean_polygon_coords_3d(&raw_hole, 0.9995);
                        if clean_hole.len() >= 3 {
                            panel_polygons_3d.push(clean_hole.into_iter().map(|p| [p.x, p.y, p.z]).collect());
                        }
                    }

                    if !panel_polygons_3d.is_empty() {
                        local_panels.push(MacroPanel {
                            id: (cl_idx * 1000 + sub_id) as u32,
                            panel_type,
                            stiffness_id: cl.elements[0].stiff_id,
                            plane_normal: [norm.x, norm.y, norm.z],
                            plane_d: cl.d,
                            polygons: panel_polygons_3d,
                            fe_count: cl.elements.len(),
                            connected_panel_ids: vec![],
                        });
                        sub_id += 1;
                    }
                }

                local_panels
            })
            .collect();

        panels_results.into_iter().flatten().collect()
    }

    /// Извлечение уникальных высотных отметок плит перекрытий
    pub fn extract_slab_elevations(panels: &[MacroPanel], tol_dist: f64) -> Vec<f64> {
        let mut slab_z: Vec<f64> = Vec::new();
        for p in panels {
            if p.panel_type == PanelType::Slab {
                for poly in &p.polygons {
                    for pt in poly {
                        slab_z.push(pt[2]);
                    }
                }
            }
        }

        slab_z.sort_by(|a, b| a.total_cmp(b));
        let mut unique_levels: Vec<f64> = Vec::new();
        for z in slab_z {
            if unique_levels.is_empty()
                || (z - unique_levels[unique_levels.len() - 1]).abs() > tol_dist
            {
                unique_levels.push((z * 1000.0).round() / 1000.0);
            }
        }
        unique_levels
    }
}