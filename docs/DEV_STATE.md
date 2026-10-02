# Development state

Updated: 2026-10-02
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

Private LIRA text models used as the Tier C gate (never committed):
`скала1` (right section, 59 surfaces / 525 axes since planar panel
growing), `типовая секция` (23 / 291), `тест 5` (108 / 8), `тест 6`, and
since 2026-09-30 eight more: a test slab, two АЖТ tests, скала with seismic
loads, Багратион concrete bedding, "для testa", ЖК Остров (piles, curved
walls) and Багратион v4 (160 MB).

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

Scoreboard 2026-10-02 (12 fixtures; assembly checker, strict audit, trial
mesh valid / external-mesher ready, PLAXIS profile): скала1, типовая секция,
тест 5, тест 6, the test slab, both АЖТ tests, скала seismic and "для testa"
pass as before (byte-identical where no rule applies). Open:

- ЖК Остров: frame accepted (default step, 11 min), all surfaces and bars
  built, assembly checker passes, trial mesh valid and mesher-ready.
  Strict audit: 26 coplanar overlaps of 1-3 mm2, 23 unshared point
  contacts (0.2 um - 0.9 mm), 4 unrepresented intersections;
- Багратион bedding: mesh valid, 18 bars in a slab plane without contact,
  1 overlapping bar, trial mesh angles (mesher gate false). Closed crack
  voids keep only shared nodes (2026-10-02): triangles under 1 degree
  195 -> 99, sites 154 -> 63;
- Багратион v4: frame accepted on the relaxation step, assembly checker
  passes, every surface built, 1 of 40 600 bars rejected
  (axis_direction_conflict), trial mesh topologically valid (2.1 M
  triangles). Micrometre vertex pairs removed (interval ends snap to
  contour vertices, micro-overlaps of coplanar slabs close): no triangle
  edge under 0.1 mm, 51 triangles under 1 degree (minimum 0.037), all at
  millimetre features (0.45 mm - 1 cm), as for ЖК Остров and Багратион
  bedding.

Rules added 2026-10-01/02 (geotechnical; deviation from the source is
allowed, user decision 2026-10-01): over-constrained panel merges only as a
frame relaxation step (`relaxation.frame` in the report); bar-bar crossings
share a generated vertex; bar ends within 50 mm of another bar's span join
it (`bar_tees`), merging into a bar node they project onto; bars lying on
each other exchange nodes and same-stiffness contained bars are removed as
duplicates; points within ten precisions are one point (interval ends, bar
nodes); short bars may shorten within their end budgets; bar nodes slide
onto an adjacent crossing.

Rules added 2026-10-02 (Багратион v4 mesh): a bar end merging into the
bar node it projects onto collapses the bar piece joining them; a node of
one bar next to a bar crossing slides onto it (bars ending there follow);
a crossing next to another generated crossing is still shared; an invalid
contour with a hanging node of the source mesh closes its seam (crack
width, filled-area budget; reported as a `seam` simplified hole); mesh
promotion of Steiner vertices and edge/bar chain synchronization use the
point slack (ten precisions) like the assembly.

Rigid links (user decision 2026-10-02): geotechnical mode removes bars
whose numeric stiffness gives a radius of gyration sqrt(EI/EF) beyond
their length and 1 m (LIRA "1000 200000 200000 200000"; Багратион v4
36 011 elements, Багратион bedding 2143, none elsewhere). Point slack is
10 um. Other rules for the reviewed defect classes use the 5 cm gap
tolerance (variant A).

State after this batch:
- Багратион v4: frame 25 min without links (relaxation step), no
  rejected bar or surface, assembly checker passes, trial mesh valid and
  mesher-ready; 22 triangles under 1 degree at 4 bars along slab edges
  with a 1-2 mm kink (one contour neighbour on the bar).
- Багратион bedding: mesher-ready, no bar audit issue; 12 sites, all in
  slab contours; 5 surfaces flagged by the global audit's ring check
  (pre-existing, to examine).
- ЖК Остров unchanged (mesher-ready; strict audit: 25 coplanar overlaps,
  21 unshared point contacts, 4 unrepresented intersections).

Speed: frame retries continue one LSQR run (ЖК Остров frame 34 -> 11 min),
`--v2-frame-cache PATH` reuses a solved frame (development), trial mesh
contacts indexed by surface (Багратион v4 mesh 9 min -> under 1 min),
`TOPO_DIAG=1` prints mesh failure details.

Previous fixtures:

1. тест 6: PLAXIS profile passes (compound alignment corners closed);
2. скала1: one contour corner of 8 degrees (a real wedge-shaped stiffness
   zone, kept by user decision);
3. trial mesh quality on тест 6 (minimum 7.09 degrees); the final mesh
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
