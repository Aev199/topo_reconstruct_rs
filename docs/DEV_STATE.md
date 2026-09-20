# Development state

Updated: 2026-09-20

This file is intentionally short. Detailed reconstruction rules remain in
`docs/GEOTECHNICAL_GEOMETRY.md` and `docs/RECONSTRUCTION_V2.md`. The verified
Gmsh experiment is documented in `docs/GMSH_OCC_TRIAL.md`.

## Product objective

Produce robust structural geometry from imperfect FE input for PLAXIS geometry
and a quality conforming MIDAS mesh. Exact reproduction of the source FE
tessellation is secondary.

## Verified Rust baseline

Before the Gmsh/OpenCASCADE experiment the private full model had a usable local
surface mesh, approximately:

- 22,170 shell triangles;
- maximum triangle area about 0.400 m²;
- minimum triangle angle about 2.39°;
- 163 triangles below 20°.

The remaining blocker was global surface-junction conformity.

## Verified Gmsh/OpenCASCADE experiment

On an existing v2 preview of the private full model, the independent input audit
found 1,258 finite surface-contact segments and 524 unrepresented junction
segments.

Passing the reconstructed planar surfaces through OpenCASCADE General Fuse via
Gmsh `occ.fragment()` fixed **524 / 524** previously unrepresented segments.
A strict shared-mesh-edge audit found:

- 0 remaining unrepresented tested surface contacts;
- 0 regressions among previously conforming contacts;
- all tested T-junction and interior-crossing classes conforming.

The exact Boolean operation produced 444 output faces from 387 input surfaces.
The returned input-to-output map is sufficient to carry source-surface
provenance through fragmentation.

A global fuzzy Boolean tolerance is rejected: even 0.1 mm changed the
fragmentation too aggressively and lost expected junction representation.

## Selected experimental mesh policy

Use exact General Fuse, then:

- explicit target size, currently 0.75 m for the full-model trial;
- `MeshSizeFromPoints=0`;
- `MeshSizeFromCurvature=0`;
- `MeshSizeExtendFromBoundary=0`;
- 10 smoothing passes;
- `Relocate2D` optimization;
- local minimum-edge cleanup only for connected sub-resolution mesh edges.

The full-model General Fuse produced one pathological CAD micro-edge of about
0.036 mm. Its local mesh component is safely handled by keeping the existing
multi-surface junction node fixed and collapsing only the adjacent parasitic
nodes. No nearest-neighbor weld is allowed.

With a 0.75 m target size after this cleanup:

- 20,108 triangles;
- maximum triangle area about 0.372 m²;
- minimum triangle angle about 8.13°;
- 8 triangles below 20°;
- 0 tested surface-contact defects;
- only two parasitic nodes merged;
- maximum node movement about 0.084 mm;
- three degenerate sliver triangles removed.

This is currently better than the custom Rust meshing baseline on both tested
junction conformity and triangle-quality outliers while keeping a comparable
mesh size.

## Architecture direction

Do **not** continue building a custom general-purpose surface Boolean/junction
kernel as the primary path.

Preferred split:

1. Rust reconstructs engineering meaning: axes, planes, contours/openings,
   properties, movement bounds and provenance.
2. Gmsh/OpenCASCADE performs exact surface fragmentation and conformal CAD
   topology.
3. Gmsh generates the surface mesh for MIDAS.
4. A small deterministic post-process removes only sub-resolution micro-edge
   artifacts while preserving proven junction nodes and logging movement.
5. Independent project audits verify the result.
6. PLAXIS can consume the healed conformal geometry and remesh it itself.

The custom surface-junction branch/PR remains a fallback/reference until the
Gmsh adapter is integrated, but should not be merged as the primary solution.

## Next coherent development batch

Integrate the experimental backend without disturbing the existing working
pipeline:

- define a stable Rust -> Gmsh interchange representation;
- preserve stiffness/property/source-element ownership through
  `occ.fragment()` output mapping;
- implement deterministic minimum-edge mesh cleanup with movement/provenance
  diagnostics;
- carry rods/axes and their surface contacts into the fragmented model;
- add synthetic regression fixtures for T-junctions, crossings, openings,
  property boundaries and sub-resolution edge cleanup;
- keep the old mesher available behind a backend choice until the new path
  passes MIDAS import and PLAXIS geometry verification.

Do not optimize the old custom surface-junction implementation further unless
the Gmsh path exposes a class it cannot represent.

## Still not globally complete

The current successful audit covers finite surface-surface contact segments. It
does not yet prove:

- isolated point contacts;
- non-parallel near misses;
- bar-bar intersections;
- complete loads/properties transfer;
- actual MIDAS import;
- actual PLAXIS geometry import.

Final readiness flags must remain false until those checks exist.

## Efficient execution / Actions

Use `AGENTS.md` as the development contract.

Normal push/PR CI on `main` is one Ubuntu `cargo test --locked` job.
Cross-platform release builds are manual only. The Gmsh runtime probe in this
experiment branch is also manual-only; the private full model must never be put
in Actions.

Use local/synthetic regressions as the inner loop and the private full model only
as an integration gate.
