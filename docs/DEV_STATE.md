# Development state

Updated: 2026-09-22

## Finite junction coverage correction

The earlier `backend_ready=true` gate only compared finite surface-pair
identities before and after local repair. It did not check source geometric
intersections that never acquired shared mesh edges. The private `test5` model
exposed a wall/floor line 10.85 m long with zero shared mesh edges, despite
complete region and source-element coverage.

The OCC adapter now normalizes outer/hole wire orientation for its own
convention and fragments all input faces and curves together. The independent
finite-intersection auditor checks every source surface intersection against
actual shared mesh edge IDs, accepting a changed endpoint only with the exact
logged bounded repair move. `nonconforming_surface_junction` blocks export and
solver-mesh writing. `numpy` and `shapely` are now backend requirements for
this audit; the manual Gmsh workflow installs them.

Fresh full-model results at mesh size 0.75: `skala1` 534/534, `skala2`
272/272, `test5` 354/354 finite intersection segments conforming (one `test5`
endpoint covered by a logged 15.09 mm movement). All three retain complete
source/property/contact coverage, no unintended shared nodes and zero backend
blockers. `test5` now has 40,551 shell triangles and 39 below 5 degrees
(minimum 0.00088 degrees), versus 37,183 and 4 in the incomplete prior mesh.
These acute elements remain a solver-mesh quality risk: this is evidence of
shared topology, not evidence of solver import or acceptable target mesh
quality. Older counters below describe the superseded unjoined mesh.

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
6. Conservative local repairs handle only proven sub-resolution micro-edges
   and sliver-producing near-vertex junctions inside the explicit engineering
   movement budget; every node move is logged.
7. Unsupported geometric point-only contacts are separated by node identity
   without moving their coordinates.
8. Independent audits decide whether the result satisfies the geometry contract.

Geotechnical policy: a sub-resolution micro-opening may be closed when the
transactional repair remains valid, bounded and provenance-logged. Tiny
openings do not require preservation solely because they exist in the source;
meaningful openings and structural/material ownership still require the normal
geometry gates.

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

## Multi-model geometry regression

Geometry is now exercised on three private real FE fixtures with one common
policy. The detailed matrix is in `docs/GEOMETRY_REGRESSION_MATRIX.md`.

All three currently satisfy `backend_ready=true`, complete source/property
coverage, zero duplicate/degenerate compact elements, zero failed declared
bar/surface contacts, zero unintended final shared nodes, zero finite
surface-pair changes caused by repair and numerical surface planarity.

The new fixtures established two general rules:

- near-vertex junction regularization is allowed only for sliver-producing
  edges, inside the 50 mm engineering movement budget, never for a real source
  boundary edge, with explicit move provenance;
- unsupported isolated point-only CAD contacts are separated by node identity
  at unchanged coordinates.

These rules are inactive on the clean `skala2` control fixture and the
near-vertex rule is inactive on `skala1`, providing a non-regression check.

## MIDAS adapter state

The neutral solver-mesh boundary and the preliminary MIDAS import-plan code
remain in the branch, but FPN work is **paused**. Current development priority
is broader geometry validation.

## Current readiness

The geometry backend is multi-model regression-tested, but additional diverse
real FE fixtures are still valuable before declaring the geometry problem
closed.

Current strict gates are:

- complete region/source-element/property coverage;
- no invalid, degenerate or duplicate compact elements;
- all declared mixed-dimensional contacts conforming;
- no unintended shared-node semantics;
- no repair-induced loss or creation of finite shared-surface pairs;
- numerical surface planarity;
- bounded, provenance-logged geometry movement.

Mesh-angle statistics remain diagnostics. A few local low-angle triangles are
not by themselves a topology failure if all geometry gates pass.

Final `export_ready` remains false; solver-specific export/import work is
deliberately deferred.

## Next coherent development batch

1. Keep the three current private fixtures as mandatory geometry regressions.
2. Add more structurally different real FE fixtures when available.
3. Improve geometry logic only when a new fixture exposes a general topology
   class, not to optimize isolated mesh-angle outliers.
4. Keep `topo-reconstruct-solver-mesh-v1` stable while geometry validation
   continues.
5. Do not resume FPN serialization until geometry coverage is broad enough.

## Efficient execution / Actions

Use `AGENTS.md` as the development contract.

Normal CI is the one Ubuntu Rust test job. Cross-platform release builds are
manual. The Gmsh probe is manual-only; private full-model runs stay local and
are never uploaded to Actions.

Use synthetic regressions as the inner loop and the private full model only as
an integration gate.
