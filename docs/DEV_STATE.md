# Development state

Updated: 2026-09-20

This is the short handoff for the next development session. Detailed geometry
rules remain in `docs/GEOTECHNICAL_GEOMETRY.md` and
`docs/RECONSTRUCTION_V2.md`; the Gmsh evidence is in
`docs/GMSH_OCC_TRIAL.md`.

## Product objective

Produce robust structural geometry from imperfect FE input for PLAXIS geometry
and a quality conforming MIDAS mesh. Exact reproduction of the source FE
tessellation is secondary.

## Current architecture

Do **not** build a custom general-purpose Boolean/junction kernel as the primary
path.

1. Rust reconstructs engineering meaning: axes, planes, contours/openings,
   stiffness/property regions, movement budgets and source provenance.
2. Rust emits the versioned `topo-reconstruct-gmsh-v1` interchange.
3. Gmsh/OpenCASCADE performs exact mixed-dimensional General Fuse.
4. Only Rust-declared rod/surface contacts are explicitly mesh-embedded.
5. Gmsh creates the conforming 2D/1D mesh.
6. A conservative local cleanup handles only sub-resolution connected
   micro-edge components and logs every node move.
7. Independent audits decide whether the result satisfies the geometry contract.

The old custom Rust mesher/junction implementation remains a fallback/reference
until solver import is verified.

## Verified surface result

On the private full model, before Gmsh the independent surface audit found:

- 387 reconstructed surfaces;
- 1,258 finite surface-contact segments;
- 524 unrepresented segments, including 503 T-junction and 19
  interior-crossing segments.

The surface-only exact OpenCASCADE General Fuse fixed **524 / 524** with no
regression among already conforming contacts.

Do not use a global Boolean fuzzy tolerance for healing. Even 0.1 mm changed the
fragmentation too aggressively on the private model.

## Verified mixed-dimensional full-model result

The current private v2 preview contains 387 reconstructed surfaces and 564
assembled axes. With surfaces and rods passed through the same exact General
Fuse:

- 536 output CAD faces;
- 2,601 input rod curve pieces and 2,601 mapped output curve pieces;
- 0 surface ownership conflicts;
- 0 curve ownership conflicts;
- 0 unmapped output surfaces;
- 4,354 declared rod/surface contacts.

General Fuse alone does not force every coplanar interior rod curve to become a
2D mesh constraint. Therefore the backend explicitly calls Gmsh mesh embedding
only for contacts already established by Rust semantics. On the private model:

- 66 rod curves required explicit surface embedding;
- 5 point contacts required explicit point embedding;
- 3,106 point contacts were already covered by declared interval contacts;
- the remaining contacts were already represented by OCC boundary identity.

Strict audit **before** micro-edge cleanup: **4,354 / 4,354** rod/surface
contacts conforming.

Strict audit **after** cleanup: 4,352 / 4,354. The two coordinate failures are
the same physical point and are caused solely by one logged cleanup operation.
The final healing-aware audit accepts them only because all of the following are
true:

- the contact was strictly conforming before cleanup;
- it remains a shared mesh node after cleanup;
- its original expected point is the exact logged `from` node;
- its current node is the exact logged `to` representative;
- the logged movement is below the explicit `minimum_edge` policy.

Final result: **4,354 / 4,354 contacts conforming**, with no global tolerance
relaxation and `backend_ready=true`.

## Current mesh metrics

At target size 0.75 model units:

- 22,220 shell triangles;
- 3,189 bar segments;
- maximum triangle area about 0.372234 m²;
- minimum triangle angle about 0.80084°;
- 34 triangles below 20°;
- 8 triangles below 5°.

The mixed-dimensional constraints worsen a few local angles compared with the
surface-only Gmsh trial, but the outlier population is still small. Do not tune
the whole model merely to improve the single minimum-angle number; first
localize any element that blocks solver import.

Micro-edge cleanup on this full run:

- 2 short mesh edges in one connected component;
- 2 parasitic nodes merged;
- maximum movement about 0.08365 mm;
- 3 degenerate shell triangles removed;
- 0 bar elements removed;
- 0 unresolved cleanup components.

## Semantic shared-node gate

Geometric coincidence produced by OCC is not automatically treated as a
structural connection.

Before healing, every mesh node with more than one surface/axis owner is checked
against a semantic owner graph. Allowed direct relations are:

- surface/surface: a finite shared mesh edge or a common source boundary node;
- axis/axis: a common source anchor node;
- axis/surface: an explicit Rust-declared point or interval contact.

Transitive connectivity through valid relations is allowed. Coordinate
proximity alone is never a relation.

Private full-model result:

- 3,340 shared mesh nodes checked;
- 0 unintended shared nodes;
- 1,335 direct axis/axis shared-node relations;
- 1 axis/axis relation connected transitively through valid surface semantics;
- 0 unsupported axis/axis shared-node relations.

This audit is a blocker: any disconnected semantic owner components sharing one
mesh node make the backend not ready.

## Frozen solver-mesh contract

The Gmsh backend now produces a separate
`topo-reconstruct-solver-mesh-v1` package. Solver adapters depend only on this
package, not on OCC tags or internal reconstruction structures.

It contains:

- compact global vertices;
- triangular shell elements;
- two-node bar elements;
- surface regions with stiffness/source provenance;
- bar regions with axis/span/stiffness/source provenance.

The standalone package is written only after a clean coverage audit.

Independent private-model coverage verification:

- 387 / 387 surface regions represented;
- 2,601 / 2,601 bar regions represented;
- 13,301 / 13,301 source shell elements covered;
- 2,601 / 2,601 source bar elements covered;
- 0 invalid or degenerate compact shell elements;
- 0 invalid or degenerate compact bar elements;
- 0 shell stiffness mismatches;
- 0 bar stiffness mismatches.

## MIDAS adapter state

`scripts/midas_fpn.py` now validates the solver mesh and produces a deterministic
`topo-reconstruct-midas-gts-plan-v1`:

- 1-based MIDAS node numbering;
- globally unique shell + bar element numbering;
- planned TRIA/LINE connectivity;
- dimensionality-separated stiffness mesh sets;
- retained region/provenance mapping.

Direct FPN serialization remains intentionally disabled until the exact modern
GTS NX record layout is calibrated on one tiny target-version FPN export.
See `docs/MIDAS_FPN_ADAPTER.md`.

## Current readiness

The implemented surface/surface, rod/surface and shared-node semantic topology
classes are verified on the private full model. The final compact solver mesh
also passes property/source coverage checks.

The project is **not** final solver readiness yet. Still required:

- calibrate the exact target-version GTS NX FPN serialization;
- import a tiny generated fixture into MIDAS;
- import the private full solver mesh into MIDAS;
- PLAXIS geometry export/import verification;
- loads, materials, stages and other analysis semantics later.

Final `export_ready` flags must remain false until target-program import is
verified.

## Next coherent development batch

1. Keep `topo-reconstruct-solver-mesh-v1` as the stable adapter boundary.
2. Obtain/inspect one minimal target-version GTS NX FPN reference file.
3. Implement the calibrated FPN writer from the already deterministic MIDAS
   import plan.
4. Validate a tiny generated file in MIDAS before the private full model.
5. Only after actual solver rejection consider further local mesh-quality
   repairs.

## Efficient execution / Actions

Use `AGENTS.md` as the development contract.

Normal CI is the one Ubuntu Rust test job. Cross-platform release builds are
manual. The Gmsh probe is manual-only; private full-model runs stay local and
are never uploaded to Actions.

Use synthetic regressions as the inner loop and the private full model only as
an integration gate.
