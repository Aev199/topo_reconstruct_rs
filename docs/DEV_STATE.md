# Development state

Updated: 2026-09-28
Repository baseline reviewed through the surface-junction batch (see git log).

This file is intentionally short. It is the entry point for the next development
session; detailed rationale belongs in `docs/GEOTECHNICAL_GEOMETRY.md` and
`docs/RECONSTRUCTION_V2.md`.

## Product objective

Produce valid, connected structural geometry from imperfect FE input for:

- reliable remeshing / geometry import in PLAXIS;
- a quality conforming mesh for MIDAS;
- later transfer of structural properties, loads and provenance.

Exact reproduction of the source FE tessellation is not the objective.

## Private fixtures

Three private LIRA text models are used as the Tier C gate (never committed):
`скала1` (right section, 388 surfaces / 581 axes), `типовая секция`
(23 / 291) and `тест 5` (108 / 8).

## Current verified state

Surface-surface junctions are explicit shared topology
(`assembly::junctions`, surface `embedded_edges`), and both meshes use one
subdivision of every junction edge. Verified on all three private models:

| Global surface audit | скала1 | типовая секция | тест 5 |
|---|---:|---:|---:|
| Unrepresented junction segments (before → after) | 520 → 0 | 400 → 0 | 2194 → 0 |
| Segments lacking shared mesh edges | 455 → 0 | 292 → 0 | 2194 → 0 |
| Audit passed | yes | yes | yes |

Everywhere: 0 invalid surfaces, 0 coplanar overlaps, 0 unresolved source
elements, trial topology and external-mesher gate pass, surface/axis counts
unchanged, 0 reconciliation problems. Debug and release outputs are identical.

Geotechnical assembly trims thin consoles beyond junction lines
(`assembly::consoles`; тест 5: 26.61 m² of 0.10–0.19 m slab consoles) and
aligns stacked walls to the axis of the wall below within 50 mm
(`assembly::stacking`; тест 5: 8 walls, 25 mm), closes wall ends and removes
redundant collinear vertices at short edges within 50 mm
(`assembly::cleanup`), audit still passing. Tolerances are CLI options:
`--v2-stack-offset`, `--v2-wall-end-snap`, `--v2-console-width`.
Duplicated vertices and wall ends next to contour corners are merged,
carrying bar axes along (`assembly::cleanup`, `Model::merge_vertices`).

Mesh quality (below 20°): скала1 163 → 73 (min 2.39°), типовая секция 0,
тест 5 28 (min 6.98°, none below 5°). The тест 5 baseline (53) was not comparable: its
slab meshes ignored 2194 junction segments.

## Current blockers

All four private fixtures pass the strict extended audit, including the
mesh-level connectivity and missing-property checks. Reconstruction targets
PLAXIS (`docs/GEOTECHNICAL_GEOMETRY.md`, "PLAXIS profile"):
`scripts/check_plaxis_profile.py` is the readiness gate at a target element
size of 0.5 m. типовая секция and тест 5 pass it. Remaining:

1. тест 6: 49 gaps below 50 mm left by the user's rules (21 between
   parallel structures, never joined; 28 needing a surface to leave its
   plane), 4 short edges;
2. скала1: one contour corner of 8 degrees;
3. trial mesh quality on тест 6 (minimum 3.33 degrees); the final mesh
   should come from Gmsh or the target program (`docs/PRIOR_ART.md`);
4. two junction diagnostics on тест 6 (a 19.6 mm near touch, one refused
   crossing split).

Fixed after тест 6: junction ends micrometres from existing vertices or
surface edges, near-planar wall tops, crack mouths shared across a patch,
columns through slabs without shared nodes, cracks of the converted mesh (contour rebuild, not
node welding: `--v2-node-weld` is removed), duplicate edge keys after splits (released edges stay
indexed), a stacked-alignment identification that ignored other supports of
the upper node, unbuffered report output, slow console trimming.

## Next coherent development batch

### Goal

Actual solver handoff: an exporter for the trial mesh/geometry and an import
of one fixture into MIDAS GTS NX and/or PLAXIS 3D by the user.

### Required approach

- Agree the exchange format with the user (for example NASTRAN bulk data for a
  MIDAS mesh, DXF/STEP faces for PLAXIS geometry).
- Export shells per stiffness and bars per span with source provenance;
  shared nodes must stay shared in the file.
- Keep `export_ready` false until the user confirms a successful import.

### Also pending

- Load transfer (loads are not in the reconstruction input).
- Actual MIDAS/PLAXIS import verification.
- `panic = "abort"` in release makes the Spade `catch_unwind` ineffective.

## Explicitly not complete yet

Even where the surface audit passes, global readiness is not yet
proven. Remaining audit scope includes:

- non-parallel near misses (review items only);
- load transfer verification;
- actual MIDAS/PLAXIS import verification.

Do not set `export_ready` or equivalent final-readiness flags merely because
the current surface audit becomes green.

## Efficient execution policy

Use `AGENTS.md` for the full development workflow.

In particular:

- use small synthetic regressions as the inner loop;
- use the private full model only as an integration gate;
- perform one coherent implementation/test/review batch before asking the user
  for another "continue";
- do not add model-specific exceptions;
- update this file only after verified state or priorities change.

## GitHub Actions budget

The repository is public, so standard hosted-runner compute is currently free
for public-repository Actions. Still avoid waste because artifact/cache storage,
queue time and developer attention remain finite.

Current `.github/workflows/build.yml` runs `cargo test --locked` on pushes and
PRs to `main`/`master`; Windows/Linux release binaries are built only by a
manual dispatch with `build_release`. Feature branches do not trigger CI.
During geometry iteration:

- use `[skip ci]` when a cross-platform release build adds no information;
- do not use Actions for the private full-model loop;
- avoid unnecessary reruns;
- keep artifact-producing runs for points where a binary is actually useful.

Revisit the workflow itself only if artifact storage or unnecessary
cross-platform builds become a practical issue; do not redesign CI merely for
the sake of redesigning it.
