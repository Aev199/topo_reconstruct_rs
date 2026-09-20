# Development state

Updated: 2026-09-20
Repository baseline reviewed through commit `2e1a9a2`.

This file is intentionally short. It is the entry point for the next development
session; detailed rationale belongs in `docs/GEOTECHNICAL_GEOMETRY.md` and
`docs/RECONSTRUCTION_V2.md`.

## Product objective

Produce valid, connected structural geometry from imperfect FE input for:

- reliable remeshing / geometry import in PLAXIS;
- a quality conforming mesh for MIDAS;
- later transfer of structural properties, loads and provenance.

Exact reproduction of the source FE tessellation is not the objective.

## Current verified state

The current geotechnical assembly preserves supported source coverage on the
private full `скала1` fixture:

- 388 assembled surfaces;
- 581 assembled axes;
- 0 unresolved supported surface source elements in the mesh;
- 0 unresolved supported axis source elements in the mesh;
- trial topology gate passes;
- external-mesher trial gate passes;
- no invalid individual surfaces in the current global surface audit;
- maximum reported boundary planarity error: 3.795e-15 m;
- no positive-area coplanar overlaps detected;
- no parallel projected-overlap near-face findings within the implemented
  50 mm audit check.

Interior mesh refinement was repaired. On the last verified full-model run:

- maximum triangle area: 0.399988 m²;
- minimum triangle angle: 2.39098°;
- triangles below 20°: 163;
- 22,170 shell triangles;
- 4,482 bar segments;
- refinement cap is not reached.

The user accepted the present local triangle quality as adequate for the next
MIDAS-oriented step. Do not spend the next batch chasing the 20° target unless a
new solver/import failure shows it is necessary.

## Current blocker

The global surface-junction audit does **not** pass.

Last verified audit:

- 1,426 candidate surface pairs;
- 830 conforming boundary-contact segments;
- 520 unrepresented intersection segments;
- 109 distinct affected surface pairs;
- 498 T-junction segments requiring conformity;
- 20 interior-crossing segments requiring conformity;
- 2 boundary-junction segments requiring conformity;
- 65 of those segments already conforming in the mesh;
- 455 still lacking shared mesh-edge coverage.

These are structural junctions to represent, not surfaces to delete.

## Next coherent development batch

### Goal

Make surface-surface junctions topologically explicit and mesh-conforming for
the intersection classes already detected by the global auditor.

### Required approach

- Derive intersection lines geometrically from the participating surfaces.
- Insert the junction into the shared topology, not independently into two
  coincident meshes.
- Split affected surface regions as needed.
- Synchronize mesh vertices/edges on both sides of the junction.
- Preserve property/material regions and source provenance through the split.
- Treat T-junctions, interior crossings and boundary contacts by geometry class,
  not by source IDs.
- Keep the global auditor independent/read-only.

### Regression cases before/with implementation

At minimum cover:

- a T-junction terminating in another panel interior;
- two panels with an interior crossing;
- a junction ending on an existing boundary vertex;
- a junction that crosses a property-region boundary;
- reversed surface orientation/normals;
- translated/rotated geometry;
- a nearby but non-intersecting surface pair that must remain separate.

Where practical, verify idempotence: assembling an already conforming junction
must not continue splitting or moving it.

### Batch acceptance

Tier A:

```sh
cargo test --offline
```

Relevant global-audit Python unit tests must also pass.

Before the private full-model run, the synthetic cases must demonstrate shared
model-edge identity and shared mesh-edge coverage along the whole expected
junction.

On the private full model, compare at least:

- invalid surfaces;
- coplanar overlaps;
- unrepresented intersection segments;
- affected surface pairs;
- mesh-conforming intersection segments;
- unresolved source-element coverage;
- surface/axis counts;
- triangle/bar-segment counts;
- provenance/reconciliation failures.

The batch must not trade fewer junction defects for new invalid faces,
overlaps, lost provenance or missing source-element coverage.

A useful target is to eliminate the currently implemented classes of
unrepresented junctions. If some remain, classify the residual geometry and add
a regression case instead of relaxing the audit.

## Explicitly not complete yet

Even if the next surface-junction batch passes, global readiness is not yet
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

Current `.github/workflows/build.yml` builds Windows and Ubuntu release
binaries on pushes/PRs and can also be dispatched manually. During geometry
iteration:

- use `[skip ci]` when a cross-platform release build adds no information;
- do not use Actions for the private full-model loop;
- avoid unnecessary reruns;
- keep artifact-producing runs for points where a binary is actually useful.

Revisit the workflow itself only if artifact storage or unnecessary
cross-platform builds become a practical issue; do not redesign CI merely for
the sake of redesigning it.
