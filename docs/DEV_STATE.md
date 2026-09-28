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
| Unrepresented junction segments (before → after) | 520 → 7 | 400 → 0 | 2194 → 0 |
| Segments lacking shared mesh edges | 455 → 7 | 292 → 0 | 2194 → 0 |
| Audit passed | no | yes | yes |

Everywhere: 0 invalid surfaces, 0 coplanar overlaps, 0 unresolved source
elements, trial topology and external-mesher gate pass, surface/axis counts
unchanged, 0 reconciliation problems. Debug and release outputs are identical.

Geotechnical assembly trims thin consoles beyond junction lines
(`assembly::consoles`): тест 5 loses 26.61 m² of 0.10–0.19 m slab consoles,
audit still passing.

Mesh quality (below 20°): скала1 163 → 73 (min 2.39°), типовая секция 0,
тест 5 241 (min 1.93°). The тест 5 baseline (53) was not comparable: its
slab meshes ignored 2194 junction segments.

## Current blockers

1. скала1 residual (7 segments, 3 pairs), reported, not repaired:
   duplicated source nodes at identical positions joining coplanar panels
   (`coincident_distinct_vertices`) and a wall corner 0.08 mm from a slab
   vertex that anchors a bar (`junction_vertex_near_edge_end`). Both need an
   explicit vertex-identity rule (merge with provenance, or a declared seam),
   never coordinate welding.
2. тест 5: walls above and below one slab with axes 25–30 mm apart (aligned
   outer faces, different thicknesses), often short piers near the edge,
   produce two junction lines millimetres apart and the remaining small and
   acute elements. Needs an explicit stacked-wall axis alignment rule
   (moves a wall plane; awaiting the user's decision on tolerance).

## Next coherent development batch

### Goal

Stacked-wall axis alignment and vertex identity with provenance, so small
offsets stop forcing tiny elements and the global audit passes on all three
fixtures.

### Required approach

- Classify coincident distinct vertices: same support planes and no source
  release/hinge information → merge identity with provenance; otherwise keep
  as a declared seam that the auditor accepts explicitly.
- Near-miss vertex/edge pairs inside a surface within an adaptive tolerance
  (local feature size, not a global constant): move onto the edge only when
  every plane of the vertex is kept and the move is recorded; bar anchors move
  only with their axis.
- Extend the independent audit to near misses (the current audit ignores them).
- Regression cases before/with the fix, with transforms and idempotence.

### Also pending

- Isolated point contacts, bar-bar intersections, load/property transfer.
- Actual MIDAS/PLAXIS import verification.
- `panic = "abort"` in release makes the Spade `catch_unwind` ineffective.

## Explicitly not complete yet

Even where the surface audit passes, global readiness is not yet
proven. Remaining audit scope includes:

- isolated point contacts;
- non-parallel near misses;
- bar-bar intersections;
- full load/property transfer verification;
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
