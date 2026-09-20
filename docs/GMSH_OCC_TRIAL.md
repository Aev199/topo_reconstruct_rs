# Gmsh / OpenCASCADE topology backend trial

Status: successful architecture experiment on 2026-09-20.

This document contains only aggregate metrics. The private full-model fixture is
not committed.

## Question

Can the low-level problem of intersecting reconstructed structural surfaces,
splitting them into conformal fragments, and producing one shared mesh topology
be delegated to an established CAD/meshing kernel instead of implementing a
custom surface-junction kernel in Rust?

## Prototype

The trial in `scripts/gmsh_occ_trial.py` converts the already reconstructed
planar surfaces to OpenCASCADE plane surfaces and calls Gmsh
`model.occ.fragment()` (OpenCASCADE General Fuse / BooleanFragments).

Responsibilities remain separated:

- Rust: FE recognition, engineering interpretation, plane/axis restoration,
  contours/openings, properties and provenance.
- OpenCASCADE: exact surface-surface fragmentation and conformal CAD topology.
- Gmsh: 2D meshing of the fragmented topology.
- A small explicit post-process may remove sub-resolution mesh micro-edges; it
  must never weld unrelated geometry by proximity.

The experiment does not use source element IDs or coordinates in any rule.

## Full private-model result

Input reconstruction used the existing v2 preview of the private `скала1`
fixture.

Independent pre-Gmsh surface audit:

- 387 reconstructed surfaces;
- 1,258 finite surface-contact segments;
- 524 unrepresented intersection segments in the reconstructed topology;
- 503 T-junction segments among those defects;
- 19 interior-crossing segments among those defects;
- 2 boundary-junction segments among those defects;
- 0 invalid surfaces;
- maximum surface planarity error about 4.82e-15 model units.

After exact OpenCASCADE General Fuse:

- 444 output CAD surfaces;
- 22 input surfaces were split into more than one output face;
- maximum fragments for one input surface: 9;
- every output face retained one input-surface owner through the Gmsh
  BooleanFragments output map.

A strict independent check then compared every original finite contact segment
against global mesh-edge identity on both owning surfaces:

- previously unrepresented: 524;
- fixed: **524 / 524**;
- still unrepresented: **0**;
- previously conforming contacts regressed: **0**;
- total unrepresented across all 1,258 checked contacts: **0**.

This includes all tested T-junctions and interior crossings.

## Mesh-quality finding

A naive `fragment -> mesh` run initially had a misleadingly bad minimum angle:
about 0.004 degrees. The bulk mesh was already good; only three triangles were
below 5 degrees.

The cause was isolated and deterministic: General Fuse produced one CAD edge of
about 3.57e-5 m (0.036 mm), with a neighboring mesh edge of about 0.084 mm.
These lengths are far below the reconstruction's engineering resolution. Global
Boolean fuzzy tolerances are **not** suitable for removing them: even a 0.1 mm
Boolean tolerance changed the fragmentation drastically and lost expected
junction representation.

The safe repair is local and topological:

1. detect mesh edges below the explicit minimum-edge policy;
2. form only connected micro-edge components;
3. when a component contains an existing multi-surface junction node, keep that
   junction node fixed as the representative;
4. map the parasitic neighboring nodes to that representative;
5. discard only triangles made degenerate by that exact collapse;
6. rerun the independent contact-coverage audit.

On the full model this affected one micro-edge component:

- two parasitic nodes merged;
- maximum movement: about **0.084 mm**;
- three degenerate sliver triangles removed;
- **0** surface-contact regressions after healing.

No general nearest-neighbor weld is allowed.

## Selected mesh scale

With the explicit uniform size sources configured as recommended by Gmsh
(`MeshSizeFromPoints=0`, `MeshSizeFromCurvature=0`,
`MeshSizeExtendFromBoundary=0`), 10 smoothing passes and `Relocate2D`, the
following full-model results were measured after the local micro-edge repair:

| target size, m | triangles | max area, m² | min angle | triangles < 20° | contact defects |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.50 | 36,335 | 0.160 | 8.13° | 4 | 0 |
| 0.65 | 22,666 | 0.273 | 8.13° | 10 | 0 |
| **0.75** | **20,108** | **0.372** | **8.13°** | **8** | **0** |
| 0.90 | 16,988 | 0.515 | 8.13° | 8 | 0 |
| 1.00 | 15,458 | 0.689 | 8.13° | 8 | 0 |

For comparison, the previously verified custom Rust mesh had approximately:

- 22,170 shell triangles;
- maximum triangle area about 0.400 m²;
- minimum triangle angle about 2.39°;
- 163 triangles below 20°.

The 0.75 m Gmsh trial therefore has a comparable element count and maximum
area, substantially better worst-angle behavior, and fully conforming tested
surface junctions.

## Architecture conclusion

The custom Rust surface-junction implementation should not be the primary path.

Preferred architecture:

```text
FE mesh
  -> Rust engineering recognition/restoration
  -> valid planar surfaces + openings + semantic provenance
  -> OpenCASCADE General Fuse (Gmsh occ.fragment)
  -> conformal fragmented CAD topology
  -> Gmsh surface mesh
  -> explicit sub-resolution micro-edge cleanup
  -> MIDAS adapter

For PLAXIS:
  ... -> conformal/healed CAD geometry -> PLAXIS remeshing
```

The independent Rust/Python audits remain valuable and should stay separate from
Gmsh: Gmsh produces geometry/mesh; our code decides whether the result satisfies
the geotechnical geometry contract.

## What remains before replacing the existing path

- Implement the micro-edge cleanup as a production transformation with
  provenance and movement diagnostics, not only as a trial script.
- Transfer source stiffness/property ownership through the
  `occ.fragment()` input-to-output map.
- Integrate rods/axes and their intersections with the fragmented surface
  topology.
- Export a MIDAS-consumable mesh while preserving properties and loads.
- Export or hand off CAD surfaces suitable for PLAXIS.
- Add independent regression fixtures for the Gmsh adapter.
- Verify isolated point contacts and bar-bar intersections, which are outside
  the current full-model surface audit.

The temporary Gmsh GitHub Actions probe is only for obtaining/verifying the
runtime in the current development environment. It is not intended to become a
normal per-commit CI job.
