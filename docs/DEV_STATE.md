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

## Current readiness

The geometry backend is now verified for the implemented surface/surface and
rod/surface contact classes on the private full model. It is **not** final
solver readiness yet.

Still required:

- independent bar/bar intersection audit;
- property/source coverage audit on the final compact mesh package;
- choose and implement the MIDAS GTS NX import path;
- actual MIDAS import verification;
- PLAXIS geometry export/import verification;
- loads, materials, stages and other analysis semantics later.

Final `export_ready` flags must remain false until target-program import is
verified.

## Next coherent development batch

1. Add a bar/bar intersection audit and prove that General Fuse gives shared
   mesh-node identity at every declared/actual rod intersection.
2. Freeze the compact solver-mesh result schema: global vertices, shell
   triangles, bar segments, stiffness IDs and source provenance.
3. Build the first MIDAS adapter around that schema instead of around internal
   Rust/Gmsh structures.
4. Use a tiny synthetic import fixture before trying the private full model.
5. Only after a successful target-program import consider local mesh-quality
   repairs for any solver-rejected elements.

## Efficient execution / Actions

Use `AGENTS.md` as the development contract.

Normal CI is the one Ubuntu Rust test job. Cross-platform release builds are
manual. The Gmsh probe is manual-only; private full-model runs stay local and
are never uploaded to Actions.

Use synthetic regressions as the inner loop and the private full model only as
an integration gate.
