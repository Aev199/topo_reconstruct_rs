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
