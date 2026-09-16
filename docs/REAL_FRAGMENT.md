# Auditable real-model fragments

`examples/real_fragment.rs` selects source shells from surface provenance in a
full-model preview, includes whole recognized bar axes touching those shells,
and reruns recognition, frame assembly, topology and the mesh gate. Surface
indices are command-line inputs; no building name, node ID or coordinate is
special-cased by the algorithm.

```sh
cargo run --example real_fragment -- INPUT.txt --report FULL_REPORT.json \
  --surface 267 --surface 32 --surface 9 > fragment.json
python3 scripts/check_v2_assembly.py fragment.json
```

The reference coordinates must match the input. The report preserves raw source
elements, all omitted elements sharing selected nodes (`cut_connections`), and
selected axes rejected in the global report. Cutting away external constraints
changes the problem: local acceptance is never full-model acceptance.

## Recorded skala1 fragment

Using the full preview from boundary repair (commit `3af2093`), surfaces 267, 32
and 9 select a slab around elevation 27.53 m and two adjacent walls. The input
contains 877 selected nodes, 651 shell elements and 395 bar elements. The local
assembly produces three surfaces and 31 of 33 axes, with 222 reported contacts.
Two axes fail `shared_anchors_not_collinear`. The mesh gate consequently returns
`mesh: null` and an explicit error. The local frame is accepted; its maximum
movement is approximately 1.394 mm. There are 783 omitted elements touching
selected nodes; 21 selected axes overlap failures in the global report.

`display_triangles` only fills accepted surface contours for visualization.
It is not a new finite-element mesh and carries no mesh-quality acceptance.
`export_ready` remains false. Loads, supports and releases are not transferred.

Next: reconcile the frame and shared-anchor numerical acceptance criteria,
then rerun the mesh gate with both rejected axes retained. Do not silently
remove those axes or interpret a cut-fragment result as global success.

## Latest rerun (2026-09-16)

The uploaded `скала1.txt` was parsed and rerun with the current v2 pipeline.
The full preview recognizes 348 axes with no recognition rejections. The
selected fragment (surfaces 267, 32 and 9) contains 877 nodes, 1046 source
elements (651 shells and 395 bars), 65/65 accepted axes, 271 contacts, and
783 omitted elements touching selected nodes. Frame and assembly are accepted;
`export_ready` remains false.

The mesh gate now reaches the emitted FE-candidate mesh. It contains 1065
triangles, 708 bar elements and 4336 vertices; `topology_valid=true`, while
`quality_passed=false` because the minimum angle is 12.17047° against the
20° profile threshold. The maximum triangle area is 0.492107 m² against the
0.5 m² limit, and the maximum edge ratio is 4.62125. The independent handoff
gate reports `external_mesher_ready=true` with no
`external_mesher_blockers`: this local angle warning does not by itself block
passing the checked geometry to an external mesher. `export_ready` remains
false. The run also exercises a guard for the upstream Spade 2.15.1 panic
(`Failed to locate position`); a conservative retry prevents that library
failure from aborting the report.

The refinement uses a relative minimum-area hint (0.1% of the configured
maximum area) to stop an acute constrained fan from generating microscopic
triangles. This is not a quality waiver: every emitted triangle is still
measured, and the fragment remains rejected until the 20° criterion is met.
`display_triangles` contains 152 visualization triangles; it is not an
accepted FE mesh. No source elements, releases, loads or supports are
invented or transferred.
