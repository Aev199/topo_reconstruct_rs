# Prior art: geometry from FE meshes

Reviewed 2026-09-28 to decide what not to reinvent. Sources are listed at the
end; statements about third-party tools come from their documentation.

## PLAXIS-Structure Interaction (PSI, NIP-Informatika)

User manual read in full (30 pages). PSI exports a structural model (LIRA,
SCAD, SAP2000, ETABS, Robot, CSiBridge) into PLAXIS 3D and couples the
results back (prescribed displacements or updated spring stiffnesses, or
an iterative spring update). Its geometry handling is minimal:

- plate elements with the same stiffness, in one plane and sharing 2 nodes
  (an edge) are merged into one plate; beams with the same stiffness and
  local axes are merged into one beam;
- variable thickness is averaged (arithmetic mean of nodal thicknesses);
- coordinates are rounded to a configurable number of digits, and nodes are
  compared with a configurable coordinate precision;
- plate holes are activated per phase; the soil "trace" (след) is the set of
  plates/points carrying springs;
- nothing repairs cracks, junctions, near misses or small features: the
  model is expected to be clean, PLAXIS intersects everything, and the user
  is told to inspect the export.

Our region rule (plane + stiffness + edge connectivity) is the same as the
PSI plate merge; everything after it (cracks, junctions, merges, consoles,
stacked walls, bar contacts, audits) has no counterpart in PSI.

## Mesh-to-geometry tools

- Altair HyperMesh "Surfaces from FE" / FE geometry: fits surfaces to shell
  elements with a feature angle and a deviation tolerance.
- Coreform/Sandia Cubit mesh-based geometry (MBG): builds a facetted model
  from an Exodus mesh or facets.
- Gmsh `classifySurfaces` + `createGeometry`: splits a discrete mesh by a
  feature angle and reparametrizes patches for remeshing; fails on
  non-manifold input ("wrong topology of boundary mesh").

All assume one conforming, clean mesh. None reconnects separately meshed
structures (unshared nodes, cracks, near misses) or applies structural
rules.

## Geometry kernels (building blocks)

- Open CASCADE: `BRepBuilderAPI_Sewing` (edges of different faces within a
  tolerance become shared: gap closure), General Fuse / Boolean operations
  with a fuzzy value (imprints all faces against each other, creating shared
  topology: our junction insertion), `ShapeUpgrade_UnifySameDomain` (merges
  faces on one surface: our coplanar merge), ShapeFix/BRepCheck.
- CGAL: region-growing plane detection and Variational Shape Approximation
  on meshes, `stitch_borders` (exact duplicate borders only),
  `autorefine_triangle_soup` with snap rounding (robust resolution of
  intersections).

They provide the operations but not the decisions (which gap is a crack,
which wall end to close, what never moves, provenance, audit).

## Structural software practice

- Revit analytical model: auto-detect adjustments align walls, floors,
  columns and beams within tolerances (our stacked-wall and wall-end rules
  follow the same idea).
- ETABS/SAP2000: automatic edge constraints "zip" non-matching meshes by
  displacement interpolation instead of making the geometry conform;
  mortar / tie constraints / rigid links in general FE practice. This is an
  option for source gaps that are not joints (mesh-level connection without
  moving geometry), usable for a MIDAS mesh but not for PLAXIS geometry.

## Research

- Attene, Campen, Kobbelt 2013, mesh repair survey: defect classes
  (degenerate, self-intersections, holes, gaps) and repair families.
- Barequet, Sharir 1995: filling gaps in polyhedron boundaries by partial
  curve matching; the optimal matching is NP-hard, so practical methods are
  heuristic (as our crack rules are).
- Cohen-Steiner, Alliez, Desbrun 2004 (VSA): planar proxy segmentation.
- Nan, Wonka 2017 (PolyFit): watertight polygonal models by selecting faces
  of a plane arrangement with binary optimization; designed for closed
  buildings from point clouds, not open shell-and-bar structures.

No source reconstructs a connected, simplified geotechnical B-rep from an
imperfect building FE model; the building rules remain ours.

## Consequences for this project

1. Keep the reconstruction rules, provenance and audits.
2. Do not maintain a production mesher: use Gmsh (surfaces with embedded
   junction curves) or the target program for the final mesh; the Spade
   mesh stays a verification tool.
3. Use a kernel (Open CASCADE) for export (STEP for PLAXIS) and as an
   independent validity check (BRepCheck), not for the decisions.
4. For source gaps that are not joints, a mesh-level tie/rigid-link option
   (MIDAS) is an alternative to geometric closure.

## Sources

- PSI user manual, "PLAXIS Structure Interaction" (NIP-Informatika).
- https://2022.help.altair.com/2022.3/hwdesktop/hwx/topics/pre_processing/geometry/surfaces_create_from_fe_t.htm
- https://coreform.com/cubit_help/geometry/model_definitions/mesh_based_geometry.htm
- https://gmsh.info/doc/texinfo/gmsh.html
- https://arxiv.org/pdf/2001.02542
- https://dev.opencascade.org/doc/refman/html/class_b_rep_builder_a_p_i___sewing.html
- https://dev.opencascade.org/doc/overview/html/specification__boolean_operations.html
- https://dev.opencascade.org/doc/refman/html/class_shape_upgrade___unify_same_domain.html
- https://doc.cgal.org/latest/Shape_detection/index.html
- https://www.cgal.org/2019/01/29/VSA/
- https://doc.cgal.org/latest/Polygon_mesh_processing/index.html
- https://www.cgal.org/2025/06/13/autorefine-and-snap/
- https://dl.acm.org/doi/10.1145/2431211.2431214
- https://doi.org/10.1016/0167-8396(94)00011-G
- https://dl.acm.org/doi/10.1145/1015706.1015817
- https://github.com/LiangliangNan/PolyFit
- https://help.autodesk.com/cloudhelp/2018/ENU/Revit-Analyze/files/GUID-D2F0E7B1-68A9-4BB8-88C1-EE53E662313D.htm
- https://structuralacademy.com/article/en/edge-constraints-csi
