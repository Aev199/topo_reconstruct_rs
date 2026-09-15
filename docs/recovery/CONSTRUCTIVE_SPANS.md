# Recovery checkpoint: constructive spans (2026-09-14)

## Status

This checkpoint was reconstructed from the conversation's recorded edits and
tool outputs after the execution environment became unavailable. It is NOT a
byte-for-byte backup of the lost worktree. The companion
[constructive-spans.patch](constructive-spans.patch) restores the core behavior
against commit `babdb88b9fda379d93dc852b646d5aff823e109c`.
The patch has now been applied to the build sources and extended with the
required test rewrites, a transverse-floor regression, and a mesh-domain
barrier fix. The implementation was verified locally on 2026-09-15.

## User-approved requirement

An axis needs to be straight between directly connected constructions, not over
an entire multi-storey or multi-span chain. A beam may change direction at a
column/beam/wall connection; a column at a beam or floor connection. Segments
retain a shared vertex. Geometric division must not introduce a hinge/release.
Source FE provenance and property intervals must remain intact.

## Recovered implementation

- Split recognized axes BEFORE frame equations are constructed.
- Cut at nodes belonging to another recognized axis.
- Cut at transverse recognized shell supports.
- For continuous coplanar shell support, cut at its entry/exit rather than at
  every shell FE node. A property change alone is not a cut criterion.
- Renormalize anchors and property spans to each segment's [0,1] interval.
- Mark divided chain segments with `constructive_segment`.
- Do not impose exact vertical/horizontal direction on these segments; retain
  the existing angle, movement, length, straightness and incidence checks.
- In bar assembly use the segment's own direction.
- Before dividing a plane residual by a tiny line-plane slope, check whether
  the point already satisfies the plane within existing numerical precision.
- Update the manually constructed frame::Axis test fixture with the new field.
- When a boundary-to-boundary construction is split into connected segments,
  retain its material-domain barrier across the union of explicit constraint
  edges. True dangling branches remain open constraints.

These cuts are based on source-node incidence, not on proximity welding.
Connections without a shared source node are not newly inferred by this patch.
Continuous coplanar supports currently use recognized patch identity; review
transitions between compatible patches rather than assuming all such cuts are
necessarily structural. Short-feature protection remains active.

## Recorded real-fragment result (before environment loss)

Command used the existing real_fragment example, source skala1.txt and the
full boundary-repair report; selected preview surfaces 267, 32 and 9:

```sh
cargo run --example real_fragment -- INPUT.txt --report FULL_REPORT.json \
  --surface 267 --surface 32 --surface 9 > fragment.json
python3 scripts/check_v2_assembly.py fragment.json
python3 scripts/check_fragment_mesh.py fragment.json
```

Observed: 3 accepted surfaces, 65/65 accepted constructive spans, 0 axis
rejections, 651 shell source FEs, 395 bar source FEs, 271 contacts.
Assembly checker passed.

Trial mesh: 689 triangles and 708 bar elements; topology_valid=true.
quality_passed=false; blockers minimum_angle_not_met and maximum_area_not_met.
Minimum angle: 2.088400229501563 degrees (required 20).
Maximum area: 5.910512460709631 m² (required 0.5).
Independent mesh checker reached and failed its minimum-angle assertion,
after checking topology/provenance/contact consistency. It did NOT pass.
Export remained false. These are recorded earlier observations, not a rerun
of this recovery patch.

The cut removes external constraints: 783 omitted FEs touch selected nodes.
Success of this fragment is not acceptance of the complete building.

## Test work completed

Existing assertions based on one straight chain are intentionally obsolete:
1. assembly/bars: interior_joint_shares_one_vertex_without_splitting_axis_or_properties:
   expect 3 axes instead of 2; two beam segments retain stiffnesses 10/20,
   both share the middle vertex with the column; verify source FE coverage,
   contacts, straightness, scale/rotation/renumbering.
2. assembly/bars: fixed_middle_joint_rejects_bent_beam_without_moving_surface_vertices:
   change to acceptance of 3 axes at an offset constructive joint;
   preserve the joint coordinate (1,1.002,0), all FEs [5,6,7], no issues.
3. frame: joint_slides_along_whole_beam_under_scale_rotation_and_renumbering:
   expect 3 one-span axes; retain transform, plane-incidence and provenance
   checks; a shared endpoint no longer needs an internal sliding parameter.
4. frame: plane_and_axes_share_interior_anchor_after_regularization:
   frame and graph have 3 axes, retain shared-node and plane checks.
5. tests/mesh_pipeline.rs: expected assembled axes changes from 2 to 3;
   do not weaken mesh quality, coverage or property assertions.
6. Add a column-through-floors test: five bar nodes at z=0..4 with the
   middle node x=0.003; four bar FEs with distinct stiffnesses; transverse
   triangle shell patches through nodes z=0,2,4. Recognition returns one
   chain; segmentation returns two spans with two FEs and three anchors each.
   Both share z=2. Repeat under rotation, translation, reversed traversal and
   renumbering. Intermediate mesh nodes and stiffness changes must not cut.
7. Keep the coplanar interior FE sliding/property test (no other axis):
   it must still return one axis rather than splitting at each shell node.

The obsolete assertions were updated and the floor-through-column scenario was
added. The complete suite now passes: 59 library tests, 49 binary tests and
four integration tests (112 total). The coplanar interior FE/property test
still remains one axis, while the explicit beam/column joint is represented by
two constructive beam segments sharing one vertex.

## Verification completed

The full Rust suite and changed-file rustfmt checks pass. The built-in skewed
fragment also passes the independent mesh checker: four surfaces, three axes,
238 triangles, 16 bars, minimum angle 20.460094 degrees and maximum triangle
area 0.46875. The full-model real-fragment input is not present in this
worktree, so the historical `skala1` result above remains a recorded result,
not a new verification claim. Export and full-model equivalence remain out of
scope; `export_ready` stays false.

No workflows were manually dispatched for this checkpoint.
