# Gmsh / OpenCASCADE topology backend trial

Status: successful surface and mixed-dimensional architecture experiment on
2026-09-20.

This document contains aggregate metrics only. The private full-model fixture is
not committed.

## Architecture tested

The versioned Rust interchange contains explicit 3D surface rings, rod axes,
property spans, source-element provenance and declared rod/surface contacts.

The external backend then:

1. builds OpenCASCADE faces and property/contact-aware rod curve pieces;
2. runs exact Gmsh `model.occ.fragment()` / OpenCASCADE General Fuse;
3. preserves source/property ownership through the Boolean output map;
4. explicitly mesh-embeds only rod/surface contacts already declared by Rust;
5. creates one conforming mixed 2D/1D mesh;
6. collapses only connected sub-resolution micro-edge components with an
   existing proven shared-junction representative;
7. audits the result independently.

No model-specific source IDs or coordinates are used as repair rules.

## Surface-only evidence

Input private-model reconstruction:

- 387 surfaces;
- 1,258 finite surface-contact segments;
- 524 unrepresented intersection segments;
- 503 T-junction segments;
- 19 interior crossings;
- 2 boundary-junction segments;
- 0 invalid surfaces.

Exact General Fuse fixed **524 / 524** previously unrepresented segments and
introduced no regression among previously conforming contacts.

A surface-only 0.75 m mesh after conservative micro-edge repair had about
20,108 triangles, maximum area 0.372 m², minimum angle 8.13° and 8 triangles
below 20°.

## Why fuzzy Boolean tolerance was rejected

The exact operation produced one pathological CAD micro-edge about 0.036 mm
long. Raising the global OpenCASCADE Boolean tolerance is not a safe repair:
even 0.1 mm changed the full-model fragmentation drastically and lost expected
junction representation.

The accepted policy is exact General Fuse followed by an explicit local
post-mesh cleanup. General nearest-neighbor welding is forbidden.

## Rod integration

The current private v2 preview contains 564 assembled rod axes. Axes are split
before OCC at:

- property-span boundaries;
- reconstructed anchors;
- point-contact parameters;
- interval-contact endpoints.

The full mixed-dimensional General Fuse produced:

- 536 output surface faces from 387 input surfaces;
- 2,601 input rod curve pieces;
- 2,601 mapped output rod curve pieces;
- 0 surface ownership conflicts;
- 0 curve ownership conflicts;
- 0 unmapped output surfaces.

### Coplanar rods

General Fuse preserves coplanar rod curves as 1D OCC entities, but does not
guarantee that every such interior curve constrains the triangulation of the
owning face.

The correct repair is **not** geometric proximity inference. For each contact
already established by Rust, the backend explicitly calls Gmsh mesh embedding
only into fragments of that declared source surface.

On the private model this required:

- 66 embedded curve/face pairs;
- 5 embedded point/face pairs;
- 3,106 point contacts already covered by declared interval contacts.

## Rod/surface contact audit

There are 4,354 declared rod/surface contacts.

Before micro-edge cleanup the strict audit found:

- conforming: **4,354 / 4,354**;
- failed: **0**.

After cleanup, strict coordinate matching alone gives 4,352 / 4,354. Both
failures are the same physical point shifted by the one allowed micro-edge
collapse.

The production audit therefore remains strict in two stages. A post-healing
point is accepted only when it was strictly conforming before cleanup, remains
topologically shared afterwards, and the exact logged node move links the
original expected point to the current representative with movement below
`minimum_edge`.

For the private model the two accepted contacts both follow the logged move
`1037 -> 1038`, movement about **0.03568 mm**. No global geometric tolerance is
expanded.

Final result:

- conforming: **4,354 / 4,354**;
- accepted specifically through healing provenance: 2;
- unresolved: **0**;
- backend blockers: **0**;
- `backend_ready=true`.

## Mixed mesh metrics

With target size 0.75 model units:

| metric | mixed Gmsh backend |
| --- | ---: |
| shell triangles | 22,220 |
| bar segments | 3,189 |
| maximum triangle area | 0.372234 m² |
| minimum triangle angle | 0.80084° |
| triangles below 20° | 34 |
| triangles below 5° | 8 |

The additional rod constraints create several local low-angle shell triangles.
This is not currently treated as a reason for global remeshing; target-program
import should first identify whether any of these elements are actually
unacceptable.

## Micro-edge cleanup result

The full mixed model contained two sub-resolution mesh edges in one connected
component.

The cleanup:

- keeps the already shared multi-surface junction node fixed;
- merges two parasitic nodes;
- moves no node more than about **0.08365 mm**;
- removes three shell triangles made exactly degenerate;
- removes zero bar elements;
- leaves zero unresolved micro-edge components;
- records every `from -> to` node move and coordinate displacement.

## Architecture conclusion

The primary route is now:

```text
FE mesh
  -> Rust engineering reconstruction
  -> topo-reconstruct-gmsh-v1
  -> exact OpenCASCADE General Fuse
  -> explicit embedding of Rust-declared mixed-dimensional contacts
  -> Gmsh 2D/1D mesh
  -> provenance-preserving micro-edge cleanup
  -> independent audits
  -> solver adapter
```

For PLAXIS the pipeline can stop earlier and hand off conformal/healed CAD
geometry for remeshing.

The custom Rust surface Boolean/junction implementation should remain a
fallback/reference, not the primary path.

## Shared-node semantic audit

A bar/bar check based only on geometric intersection initially appeared to show
one unsupported relation. Inspecting the complete local topology showed why a
pairwise test is insufficient: two axes can legitimately share a node through
two surfaces whose own junction is valid.

The production gate therefore audits the complete owner graph at every shared
mesh node. Direct relations are allowed only through a finite shared surface
edge, a common source node, or an explicit Rust rod/surface contact; transitive
paths through those valid relations are accepted.

Private full-model result before healing:

- 3,340 shared nodes checked;
- 0 unintended shared nodes;
- 1,335 direct axis/axis relations;
- 1 surface-mediated axis/axis relation;
- 0 unsupported axis/axis relations.

The backend is blocked if a single mesh node contains more than one disconnected
semantic owner component.

Surface boundary source-node IDs are now exported explicitly in
`topo-reconstruct-gmsh-v1`; the audit never infers source identity from
coordinates.

## Solver-mesh handoff

The final adapter boundary is
`topo-reconstruct-solver-mesh-v1`, containing compact vertices, triangular
shells, two-node bars and separate surface/bar region registries.

The full-model coverage audit found:

- 387 / 387 surface regions represented;
- 2,601 / 2,601 bar regions represented;
- 13,301 / 13,301 source shell elements represented;
- 2,601 / 2,601 source bar elements represented;
- no invalid or degenerate compact elements;
- no stiffness mismatches.

This means MIDAS/PLAXIS adapters no longer need to understand Gmsh entity tags
or OpenCASCADE fragmentation maps.

## Remaining validation

The geometry/backend phase is sufficiently closed to move to the target-program
boundary.

Next work is:

- calibrate the exact modern GTS NX FPN element/group record layout;
- generate and import a tiny MIDAS fixture;
- then import the private full solver mesh;
- verify PLAXIS geometry handoff separately;
- add loads/materials/stages only after geometry import is stable.

The Gmsh GitHub Actions probe is manual-only. Private full-model data is never
put in Actions.
