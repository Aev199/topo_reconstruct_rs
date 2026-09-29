# Geotechnical geometry acceptance

The primary objective is valid, connected geometry suitable for PLAXIS geometry
import or a quality MIDAS mesh. Exact reproduction of the source FE discretization
is secondary. Source element IDs remain provenance; they are not a requirement for
one-to-one reproduction of nodes or triangles. Material regions, structural
connections, meaningful openings and load paths must survive simplification.

## Implemented: collapsed opening repair

The v2 CLI now uses `assembly::assemble_geotechnical` by default. The existing
`assembly::assemble` API remains conservative; `--v2-preserve-details` selects it
from the CLI for comparison. Legacy reconstruction is unaffected.

An interior opening can be closed only when all checks pass:

- Its source ring and exterior are valid, and the opening is strictly inside.
- The source vertices lie within 0.05 model units of the longest hole chord
  (50 mm for models expressed in metres).
- After shared-plane closure, the maximum distance to that chord is at most
  numerical precision, and the area is at most precision times hole extent.
- Cumulative filled source area is at most 0.1% of the exterior source area in
  the property region. A large opening collapsed by a bad solve is not accepted.
- The remaining surface passes the normal transactional contour validator.

Ordinary narrow openings that remain numerically valid are not closed. No
model-specific IDs are used. Every applied change records source nodes/elements,
source and candidate dimensions, and thresholds in `topology.simplified_holes`.
No coordinates or material regions are deleted. Hole nodes remain mandatory
interior mesh points; existing shared edges remain constraints. Axis contacts
are assembled against the repaired surfaces. Missing points/constraints fail
the trial topology check. The Python checker also verifies that all retained
hole nodes occur in triangles of the repaired surface.

## Full-model verification, 2026-09-19

Input: the private full `скала1` fixture, not a synthetic fragment. Baseline:
upstream `de0e380`. Source fixture and full reports are not committed.

| Measurement | Baseline | Geotechnical repair |
|---|---:|---:|
| Assembled surfaces | 387 | 388 |
| Assembled axes | 565 | 581 |
| Unresolved surface source elements in mesh | 1811 | 0 |
| Unresolved axis source elements in mesh | 140 | 0 |
| Triangles | 15203 | 17348 |
| Bar segments | 4273 | 4482 |
| Trial topology valid | false | true |
| External mesher trial gate | false | true |

Two collapsed openings were closed. Their source widths were approximately
0.0969 mm and 0.01235 mm, and their source areas were 2.40636e-5 m² and
2.31470e-6 m². All 13370 recognized shell elements and 2742 recognized bar
elements retain provenance. These counts concern supported structural elements,
not every record in the source file. Reconciliation reports zero remaining
problems for this run. All 125 Rust tests and the independent assembly/mesh
checker pass, including rigid-transform and source-junction retention cases.

## Remaining acceptance work

`external_mesher_ready` is the existing trial-mesh handoff gate, not proof of
successful import or a complete global surface-intersection audit. `export_ready`
remains false. The local 20° / 0.5 m² mesh-quality profile still fails: minimum
angle 0.309°, maximum area 3.611 m², and one surface reaches its refinement cap.
This mesh must not be presented as a finished MIDAS analysis mesh.

Next work: global intersections/near-coincident faces and conformity at all
surface junctions; targeted geometric simplification of constraints that force
poor triangles; actual solver import and load/property transfer verification.
Do not relax geometric validity or relabel an incomplete import as ready merely
because a source-deviation metric is acceptable.

Reproduction:

```sh
cargo test --offline
cargo run --offline -- --v2-mesh-preview-json geotechnical.json model.txt
python3 scripts/check_v2_assembly.py geotechnical.json
cargo run --offline -- --v2-preserve-details --v2-mesh-preview-json strict.json model.txt
```

## Interior refinement correction, 2026-09-20

The earlier missing interior points were a meshing defect, not evidence that
their surfaces needed geometric simplification. Internal beam constraints could
be mistaken for winding boundaries by the triangulator. Disabling exterior
exclusion avoided that problem but spent the shared surface budget outside the
material; subsequent regions could then receive no refinement.

Size control now seeds the entire material domain before local angle refinement.
Only actual surface boundaries toggle inside/outside membership; internal bars
remain constraints. Exterior refinement is disabled. A material-only angle pass
revisits faces the library may have excluded because of internal constraints,
protects fixed edges from encroachment, and recalculates circumcenters after
each insertion. Both stages remain bounded by the configured vertex budget.
Coordinates, contours, source-property provenance and shared boundaries are
unchanged; the additional nodes are mesh nodes, not restored source FE nodes.

| Full скала1 measurement | Before | After |
|---|---:|---:|
| Maximum triangle area, m² | 3.61125 | 0.399988 |
| Minimum triangle angle | 0.30886° | 2.39098° |
| Triangles below 20° | 462 | 163 |
| Triangles | 17348 | 22170 |
| Bar segments | 4482 | 4482 |
| Refinement cap reached | yes | no |

All 388 surfaces and 581 axes remain represented, with complete supported
source-element coverage. Trial topology and the external-mesher gate pass.
The only remaining local quality blocker is `minimum_angle_not_met`; export
readiness is still false. Do not attribute the remaining acute triangles to
bad source geometry without a targeted constraint audit.

Validation: 126 Rust tests pass, including a new scale-varied case with a hole
and a dangling beam verifying material-only seeding, preserved constraints,
area coverage and the zero-budget case. The independent assembly checker passes;
an independent calculation from all triangle coordinates confirms maximum area,
minimum angle, positive triangle areas and the 163 remaining angle violations.

## Global surface-junction audit, 2026-09-20

The user accepted the current triangle quality for a trial MIDAS import; further
20-degree refinement is not the priority. The independent global audit of the
`b18e2d4` full-model report **does not pass**. The earlier trial handoff flag must
not be interpreted as global geometric validity.

`scripts/check_v2_global_geometry.py` checks planar polygon validity, positive-area
coplanar overlap, and positive-length surface intersections. It distinguishes
boundary contacts, T-junctions and interior crossings. A junction is represented
only when shared model edge IDs cover its entire length. When a mesh is supplied,
shared triangle edge IDs must independently cover that length. Merely coincident
coordinates do not establish a shared connection. Each issue includes surface IDs,
endpoints, length and both conformity flags. Counts refer to intersection segments,
which may be subdivided by source boundary vertices, not independent defects.

| Full-model check | Result |
|---|---:|
| Surfaces checked | 388 |
| Candidate surface pairs | 1426 |
| Invalid individual surfaces | 0 |
| Maximum boundary planarity error, m | 3.795e-15 |
| Positive-area coplanar overlaps detected | 0 |
| Parallel face pairs within 50 mm with projected overlap | 0 |
| Conforming boundary contact segments | 830 |
| Unrepresented intersection segments | 520 |
| Distinct affected surface pairs | 109 |
| T-junction segments requiring geometry conformity | 498 |
| Interior crossing segments requiring geometry conformity | 20 |
| Boundary junction segments requiring conformity | 2 |
| Above segments already conforming in the mesh | 65 |
| Above segments also lacking shared mesh edges | 455 |

These are meaningful structural junctions, not surfaces to delete. For example,
surfaces 0/34 meet along an 8.5 m T-junction and surfaces 35/270 have an interior
crossing about 5.1134 m long. The next repair must insert their intersection lines
as shared constraints and split affected surface regions, preserving property
provenance and synchronizing mesh vertices on both sides. The auditor itself is
read-only and does not change the accepted mesh or geometry.

The audit is deliberately incomplete: isolated point contacts, nonparallel
near-misses, bar-bar intersections and load/property transfer are not checked.
`audit_complete` and `solver_import_verified` remain false even if the implemented
surface checks pass. Near-face findings are review items, not automatic merges.
A successful PLAXIS or MIDAS import has not been demonstrated.

Reproduction (optional Python audit dependencies):

```sh
python3 -m pip install -r scripts/requirements-audit.txt
python3 -m unittest discover -s scripts -p test_v2_global_geometry.py
python3 scripts/check_v2_global_geometry.py refined-final.json --output global-audit.json --strict
```

Ten analytic regression tests cover shared boundaries, T-junctions, crossings,
coplanar overlaps, near faces, disjoint surfaces, holes, duplicated mesh node IDs,
and translated/reversed-normal geometry. The full-model strict run exits 1 as
expected and writes the complete diagnostic report before exiting. Private source
coordinates and the full report are not committed.

## Explicit surface junctions, 2026-09-28

Surface-surface junctions are now shared topology. After axis assembly (the
last stage that may move boundary anchors) `assembly::junctions::insert`:

- intersects every pair of non-parallel surfaces (plane-plane line clipped by
  both closed contours) and every pair of coplanar surfaces with overlapping
  collinear boundary edges;
- builds one ordered chain of vertices along each junction segment from the
  existing vertices of both surfaces and their mandatory interior nodes (bar
  contacts, retained nodes of closed openings);
- creates a vertex only where a junction ends inside an edge, at the exact
  intersection of that edge with the other plane, and splits the edge
  globally, for every surface using it;
- keeps a junction edge on a contour as a boundary edge and records it inside
  the material as a surface `embedded_edges` entry with the same edge id.

No property region is divided and no source element changes owner. Embedded
edges are mesh constraints: the mesh of both surfaces uses one global
subdivision of each junction edge. The pass is idempotent.

Sub-millimetre source misalignments are closed rather than turned into
parasitic edges. When a required junction vertex falls within the minimum
edge length (1 mm by default) of an existing unlocked vertex, that vertex is
moved onto the junction; likewise an embedded line ending within that
distance of an edge of the same surface is extended onto it. A move is
accepted only if the vertex stays on every plane it lies on and all affected
contours revalidate; each move is recorded in `junctions.snapped_vertices`.
Bar anchors and mandatory interior nodes never move. Distinct vertices at
one location are reported (`coincident_distinct_vertices`) and never merged,
since a duplicated source node may be an intentional seam or hinge.

Mesh changes made necessary by embedded junction lines:

- Spade classifies refinement faces by constraint parity. A dangling internal
  line (a wall ending inside a slab) made it treat entire material regions as
  exterior, so their angle refinement was skipped. Regions are now refined with
  their closed boundary only; internal constraints are then restored after
  removing Steiner points that encroach them, and a material-only pass
  (circumcenter, off-center or centroid; independent points inserted in
  batches) repairs angle and area near them.
- Shared constraint edges are subdivided with density `1/size`, where size is
  the distance to the nearest non-adjacent edge of any owning surface, capped
  by the nominal spacing. Uniform subdivision is reproduced exactly when no
  feature is nearby. Constraint chains stay global, so conformity is unchanged.

JSON schema additions: `topology.junctions` (report), surface
`embedded_edges`, `preview.orphaned_edges`. Vertices with index
`>= vertex_source_nodes.len()` are generated junction vertices listed in
`topology.junctions.generated_vertices`; they have no source node. The Python
auditors accept embedded edges only after checking that each lies inside its
surface and on its plane, and that none duplicates a contour edge.

Private full-model verification (fixtures and reports not committed;
unrepresented = audit segments lacking shared model edges):

| Measurement | скала1 before | after | типовая секция before | after | тест 5 before | after |
|---|---:|---:|---:|---:|---:|---:|
| Surfaces / axes | 388 / 581 | 388 / 581 | 23 / 291 | 23 / 291 | 108 / 8 | 108 / 8 |
| Invalid surfaces, coplanar overlaps | 0, 0 | 0, 0 | 0, 0 | 0, 0 | 0, 0 | 0, 0 |
| Unrepresented junction segments | 520 | 7 | 400 | 0 | 2194 | 0 |
| Affected surface pairs | 109 | 3 | 40 | 0 | 229 | 0 |
| Segments lacking shared mesh edges | 455 | 7 | 292 | 0 | 2194 | 0 |
| Global surface audit passed | no | no | no | yes | no | yes |
| Unresolved surface / axis source elements | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| Trial topology valid, external-mesher gate | yes, yes | yes, yes | yes, yes | yes, yes | yes, yes | yes, yes |
| Triangles | 22170 | 23196 | 14514 | 14346 | 34065 | 49485 |
| Bar segments | 4482 | 4537 | 3032 | 3034 | 79 | 79 |
| Minimum triangle angle | 2.391° | 2.391° | 20.17° | 20.21° | 9.67° | 1.93° |
| Triangles below 20° | 163 | 73 | 0 | 0 | 53 | 261 |
| Maximum triangle area, m² | 0.400 | 0.400 | 0.469 | 0.469 | 0.432 | 0.400 |
| Reconciliation problems | 0 | 0 | 0 | 0 | 0 | 0 |

The "before" meshes of тест 5 were locally better only because 2194 junction
segments were not in the slab meshes at all; that mesh was not conforming.
The remaining low angles there are driven by source features 4.6–25 mm apart
(for example wall lines offset by a few millimetres on one slab), above the
1 mm snap limit. They need an explicit near-miss rule with its own evidence,
not a larger silent tolerance.

The 7 скала1 residual segments are explained and reported: three coplanar
panels of one slab meet with duplicated source nodes at identical positions
(`coincident_distinct_vertices`), and one wall corner lies 0.08 mm from a slab
vertex that anchors a bar and therefore may not move
(`junction_vertex_near_edge_end`). Both need an explicit vertex-identity
decision, not coordinate merging.

Validation: 138 Rust tests pass, including 11 synthetic junction cases (T-junction,
interior crossing, junction leaving a panel, ending on an existing vertex,
crossing a property-region boundary, walls meeting on a slab, interior node on
a junction, micro-offset corner, wall ending short of a slab edge, coincident
vertices, nearby separate panels) under three scales, a rigid transform and
reversed normals, all with an idempotence check; and an FE-to-mesh test with a
T-shaped and a crossing wall verifying shared mesh edges on both sides and the
full 20° / area quality profile. 12 Python audit tests pass. Debug and release
outputs of all three private models were compared field by field.

Toolchain note: with `lto = true` and one codegen unit, rustc 1.94.1 produced a
wrong value for an `Option<DVec3>` returned by an inlined helper (reproduced only
in release; debug and non-LTO release were correct). The helper was restructured
and the three-model debug/release comparison is identical. Keep that comparison
in the full-model gate. Separately, `panic = "abort"` in the release profile
means the existing `catch_unwind` around Spade refinement cannot turn a library
panic into a diagnostic in release builds.

## Console trimming, 2026-09-28

A mid-surface FE model lets a slab reach the outer face of a wall whose
mid-plane lies inside it: a console of half the wall thickness (0.10–0.20 m on
тест 5, 171 m of free edge in total). It carries no structure but forces
elements far below the target size. Geotechnical assembly
(`assembly::consoles`, after junction insertion) now trims such consoles; the
conservative `--v2-preserve-details` path does not.

Each surface is divided by its contour and junction lines into faces
(separate components such as openings or floating junction lines belong to the
smallest face containing them). A connected group of faces is removed when:

- every face lies within `maximum_console_width` of its own junction lines
  (default 0.25 model units: half the default mesh spacing, so such a strip
  cannot hold target-size elements);
- each face touches the free contour and the group borders material that stays
  across a junction line (it is beyond a junction, not an isolated panel);
- no removed contour edge is used by another surface (a parapet or facade
  panel on the edge keeps the console);
- no opening, dangling line, retained node, bar anchor or bar interval lies in
  it;
- the remaining material is one valid region.

The junction line becomes the surface contour; its edge id is already shared
with the wall, so the connection stays topological. Passes repeat until no
change: trimming a slab console frees the contour of a wall stub beyond a
perpendicular wall, which is then trimmed; a stub standing on a wider plate
stays. Source elements remain assigned to their surface; each trim reports
width, area, free length, junction edges and released edges
(`topology.consoles`), and faces kept for a reason are counted by reason.

тест 5 (private): 6 surfaces trimmed, 26.61 m², maximum width 0.19 m, three
consoles kept because a parapet/facade surface shares their edge. The global
surface audit still passes; triangles 49485 → 41975; below 20° 261 → 241;
minimum angle unchanged at 1.93°. скала1 and типовая секция have no such
consoles and are unchanged. Debug and release outputs are identical.

The remaining small elements on тест 5 come from walls above and below one
slab whose axes are 25–30 mm apart (outer faces aligned, different
thicknesses), often short piers near the slab edge. Their two junction lines
cannot be meshed without millimetre elements. Aligning stacked walls to a
common axis moves a wall plane and is left as an explicit next decision.

Tests: 5 synthetic cases (T-junction and crossing wall, too wide, parapet on
the edge, retained node, opening plus floating line, perimeter annulus,
corner consoles with wall stubs) under three scales, rigid transforms and
reversed normals, and an FE-to-mesh test verifying trim width/area, full mesh
quality and shared slab/wall mesh edges.

## Stacked wall alignment, 2026-09-28

User decision: a wall standing on a slab adopts the axis of the wall carrying
it (the lower, bearing wall), with a 50 mm tolerance. Geotechnical assembly
only (`assembly::stacking`, `FeaturePolicy::maximum_stack_offset`, default
0.05 model units = the default junction movement limit).

A pair qualifies when both patches are vertical and parallel (frame angle
tolerance), lie on opposite sides of one horizontal slab, are in line contact
with it and overlap along the contact, and the upper patch lies within the
tolerance of the lower plane. Contact is geometric: a wall resting on a slab
often shares no source node with it; a node counts when it lies on the slab
plane within the closure tolerance and within three median node spacings of a
slab node. Stacks are processed bottom-up, so every wall adopts the plane of
the lowest wall of its stack, and the offset to that plane is checked again.

The upper patch receives the lower support (`support_representatives`); every
vertex then closes through the ordinary support intersection. Before
committing a pair, every affected vertex is checked against the junction
movement limit and its cumulative movement budget; a failing pair is kept and
reported. After alignment, upper and lower junction nodes lie on one line;
one-to-one pairs less than the minimum edge length apart along it share one
vertex (`stacked_walls.identified`, bar anchors excluded, movement rechecked).
This identity is justified by the alignment itself; coincident source nodes
elsewhere are still only reported.

тест 5 (private): 8 walls aligned (offset 25 mm each), 56 node pairs
identified (at most 18 µm apart), global surface audit passes, triangles
41975 → 36780, below 20° 241 → 131, below 5° 33 → 25, minimum angle 1.59°.
скала1 and типовая секция have no stacked offsets and are unchanged. Debug and
release outputs are identical.

Remaining acute triangles on тест 5 come from wall ends 5–30 mm short of, or
beyond, the axis of a perpendicular wall (a 12 mm node spacing on the junction
chain next to a 0.4 m element). Closing them needs a larger, evidence-based
snap tolerance for wall ends than the 1 mm used now.

Tests: an FE-to-mesh case (upper wall 25 mm off, resting on the slab without
shared nodes) under three scales and rigid transforms checks the alignment,
seven identified nodes, full mesh quality and one shared junction line for
both walls and the slab; conservative assembly leaves the walls unchanged.

## Wall ends, redundant vertices and CLI tolerances, 2026-09-28

User decision: close wall-end offsets up to 50 mm, and make the tolerances
configurable. The remaining acute triangles on тест 5 had two causes:

- A source node left on a straight contour or junction line a few
  millimetres from a wall end (typically the slab node where a wall stood
  before stacked alignment). `assembly::cleanup` removes a vertex joining
  exactly two edges with identical users when it bounds an edge shorter than
  the wall-end tolerance, is not a bar anchor or retained node, and lies within
  1/1000 of the tolerance of the straight line of its neighbours (50 µm at
  50 mm). Contours change by at most that deviation, which is reported per
  vertex; real corners are never removed.
- The free end of a junction line (a wall end inside a slab) short of another
  line of the same surface. The in-surface near-touch rule now moves such an
  end up to the wall-end tolerance (other vertices still less than the minimum
  edge). The target is the intersection of every plane of the vertex with the
  planes of the target edge, so the wall keeps its plane; an oblique target
  requiring a longer move is rejected and reported (`embedded_near_touch`).

Command-line options (model units; 0 disables the rule):

| Option | Default | Controls |
|---|---:|---|
| `--v2-stack-offset` | 0.05 | stacked-wall alignment tolerance; the vertex closure movement limit is raised to at least this value |
| `--v2-wall-end-snap` | 0.05 | wall-end closure and redundant-vertex short-edge threshold |
| `--v2-console-width` | 0.25 | maximum trimmed console width |

They apply to geotechnical assembly only; `--v2-preserve-details` ignores
them.

тест 5 (private), after stacked alignment → after this batch: redundant
vertices removed 16 (maximum deviation 10.8 µm), wall ends closed 8 (largest
move 0.2 mm; the rest were already within the minimum edge), triangles 36780 →
36513, below 20° 131 → 37, below 5° 25 → 0, minimum angle 1.59° → 6.98°.
The global surface audit still passes. скала1 and типовая секция are unchanged.
Debug and release outputs are identical.

Six wall ends on тест 5 (three walls on two levels) stop 43–49 mm from a slab
contour corner that itself lies 11 mm off the wall axis. Closing them requires
merging two vertices and moving the corner, i.e. a vertex-identity decision;
they are reported, not repaired. Five short edges remain.

Tests: a wall end 30 mm short of a perpendicular wall axis closes with a
50 mm tolerance and stays open without it; a stray collinear contour node is
removed while a locked one and a real 12 mm notch are kept; all under three
scales, rigid transforms and reversed normals.

## Vertex merging, 2026-09-28

User decision: duplicated nodes are merged; the geotechnical model is meant
to be simplified. Geotechnical assembly now merges (`assembly::cleanup`,
`Model::merge_vertices`, transactional, every affected surface revalidated):

- before junction insertion, distinct vertices closer than the minimum edge
  length (duplicated source nodes), keeping a bar anchor if one is involved
  (`topology.coincident_vertices`);
- after junction insertion, the free end of a junction line (a wall end) into
  the nearest vertex of the same surface within `--v2-wall-end-snap`, the kept
  vertex moving onto the planes of both (`topology.wall_ends`).

The kept vertex moves to the intersection of every plane of both vertices,
within the tolerance; a bar anchor never moves and two bar anchors are never
merged (their axes are not re-solved). Each merge records both vertices, their
source nodes, the distance and the movement; each rejection records a reason.
Near-touch diagnostics of the junction pass that a later merge resolves are
dropped from `topology.junctions.issues`.

Private fixtures:

| | скала1 | типовая секция | тест 5 |
|---|---:|---:|---:|
| Duplicates merged / rejected | 1 / 2 | 0 / 0 | 0 / 0 |
| Wall ends merged into corners | 0 | 0 | 6 (corner moved 11–19 mm) |
| Unrepresented junction segments | 7 → 5 (1 pair) | 0 | 0 |
| Triangles below 20° | 73 | 0 | 37 → 28 |
| Minimum angle | 2.39° | 20.21° | 6.98° |

скала1 rejections: two duplicated nodes that are both bar anchors (0.2 µm
apart), and a bar anchor 36 µm off the plane of a wall whose corner lies next
to it. The last remaining unrepresented pair comes from that anchor. Closing
it needs moving a bar anchor with its axis, i.e. re-solving the axis after a
merge. Debug and release outputs are identical.

Tests: duplicated nodes merge into one shared edge and are idempotent;
duplicated bar anchors are kept; a wall end 30 mm short of a slab corner that
is 10 mm off the wall axis merges into the moved corner; all under three
scales, rigid transforms and reversed normals.

## Axis-aware vertex merging, 2026-09-28

Merges may now involve bar anchors. A bar anchor moves only as the end of its
axes: the other end stays, and every other anchor is re-placed at its own
parameter on the new straight line, keeping all planes of its surfaces
(`Model::move_vertex` validates each). An interior anchor, an anchor shared
with another axis, a move beyond the tolerance, or a merge of two vertices of
one axis is rejected with its reason. After a merge the axes' endpoints and
anchors and the point contacts refer to the kept vertex; parameters and spans
are unchanged. Retained nodes of closed openings never move. The assembly
checker treats a dropped source node as the node it was merged into.

Private fixtures, all with the global surface audit **passing**:

| | скала1 | типовая секция | тест 5 |
|---|---:|---:|---:|
| Coincident merges (rejected) | 3 (0) | 0 | 0 |
| Wall-end merges | 0 | 0 | 6 |
| Unrepresented junction segments | 0 | 0 | 0 |
| Junction diagnostics | 0 | 0 | 0 |
| Triangles below 20° / minimum angle | 73 / 2.39° | 0 / 20.21° | 28 / 6.98° |

скала1 merges: two duplicated nodes (0 and 0.2 µm, the latter both bar ends)
and a bar end 91 µm from a wall corner, moved 36 µm onto the wall plane with
its bar. Surface/axis counts, source coverage, trial topology and the
external-mesher gate are unchanged; reconciliation reports no problem. Debug
and release outputs are identical.

The implemented surface audit is complete for its scope; bar-bar
intersections, isolated point contacts, non-parallel near misses, load and
property transfer, and an actual MIDAS/PLAXIS import remain unverified.

Tests: duplicated bar ends merge and both bars stay straight with the contact
following the kept vertex; a bar passing through a duplicate is not bent.

## Extended independent audit, 2026-09-28

`scripts/check_v2_global_geometry.py` now also checks, read-only:

- bar-bar: crossings or touches without a shared node (`unshared_bar_intersection`)
  and overlapping collinear bars (`overlapping_bars`);
- bar-surface: a bar piercing or touching a panel without a point contact or a
  shared vertex (`unshared_bar_surface_intersection`), and a bar lying in a
  panel without an interval contact (`bar_in_surface_without_contact`);
- surface point contacts: a vertex of one panel on another panel that is not
  one of its vertices (`unshared_point_contact`);
- property transfer: every triangle carries its surface stiffness and every
  bar segment a stiffness of its axis spans.

Review items (reported, never failing): surface-surface and bar-surface
near misses within `--near-distance` (default 0.05), bar near misses between
bars not already joined through a neighbouring bar, and bars shorter than that
distance. `global_checks_passed` combines all failing classes; `--strict`
exits 1 on it. Loads are not part of the reconstruction input and remain
unverified; `solver_import_verified` stays false.

Geotechnical assembly additionally merges the ends of two different bars
within `--v2-wall-end-snap` (`topology.bar_ends`), carrying both bars; a merge
that would collapse a bar joining them is rejected.

Private fixtures, strict audit passing on all three:

| | скала1 | типовая секция | тест 5 |
|---|---:|---:|---:|
| Failing issues (all classes) | 0 | 0 | 0 |
| Triangles / bars with wrong stiffness | 0 / 0 | 0 / 0 | 0 / 0 |
| Review items | 2 short bars (6.9 and 28.6 mm) | 0 | 5 surface near misses (30–36 mm) |
| Bar-end merges | 0 (2 rejected: joined by short bars) | 0 | 0 |

The short bars of скала1 are source elements joining two beam pieces; they are
kept (removing them would drop a source element). The near misses of тест 5
are a wall corner 30 mm from a perpendicular wall and a stepped common edge of
two coplanar walls 36 mm from a neighbouring vertex. Debug and release outputs
are identical.

Tests: six analytic auditor cases (bar crossings with and without a node, a
column through a slab with and without a contact, a bar in a slab with and
without an interval contact, a gap versus a short joining bar, a tilted panel
touching a slab at a corner, a stiffness mismatch) and a Rust case joining a
beam split by a 20 mm gap while keeping a 20 mm joining bar.

## Mesh-level connectivity audit, 2026-09-28

An external review found three gaps in the extended audit: any bar-surface
contact record passed a whole in-plane intersection; bar checks read contact
records, not the mesh; missing stiffness data counted as zero errors. Now:

- a bar lying in a surface must have its whole in-surface length covered by
  surface edges or interval contacts (`bar_in_surface_without_contact`), and,
  with a mesh, by bar segments that are also triangle edges of that surface
  (`bar_in_surface_not_shared_in_mesh`);
- a piercing/touching bar, two crossing bars and a vertex shared by two
  surfaces must use one mesh node (`bar_surface_point_not_shared_in_mesh`,
  `bar_intersection_not_shared_in_mesh`, `shared_vertex_not_shared_in_mesh`);
  every axis must be covered by its mesh bars (`bar_not_covered_by_mesh`);
- a triangle without stiffness or without a `surface_stiffness` record, a
  mesh bar without a valid axis/spans, and a surface without triangles fail.

The new checks found a real defect on скала1: two beams in a slab plane were
connected to the slab only at their ends. Causes and general fixes:

- contacts were derived before later vertex moves (a coincident merge moved a
  beam end by 0.036 mm off the slab plane numerically, so the in-plane
  interval was never recorded). `bars::refresh_contacts` recomputes all
  bar-surface point and interval contacts from the final geometry;
- a beam node stayed 0.038 mm from the slab edge vertex it belongs to, two
  parallel lines side by side. `cleanup::merge_bar_anchors` identifies a bar
  node with a surface vertex within `minimum_edge`, moving the surface vertex
  onto the bar (a bar never bends; reported in `topology.bar_anchors`).

After the fix all three fixtures pass the strict audit with mesh-level
checks; mesh metrics unchanged (скала1 +2 triangles). Debug and release
outputs are identical. Tests: 10 new analytic auditor cases and 2 Rust cases
(beam node beside a slab edge vertex, stale contact without its interval),
each over rotations, translations, reversed normals and scales.

## Cracks of converted meshes: region contours rebuilt, 2026-09-28

The large private fixture тест 6 (converted mesh, 490k elements) had 20
planar patches that did not assemble. Causes, all inside one planar region:

- a wedge crack: two neighbouring elements end at distinct nodes 0.4-3 mm
  apart and share the node at the crack tip;
- a zero-width crack: overlapping collinear element edges end at distinct
  nodes (up to 50 mm apart), for example along a line where two meshes with
  non-matching nodes meet (patch 57: nine real 0.15 x 0.2 m openings joined
  by such zero-width channels into one self-intersecting contour);
- an edge collapsed by a stacked-wall alignment: an upper perpendicular wall
  moved 50 mm onto the wall below identifies its node with the lower node,
  and the 50 mm step edge of the wall containing both collapses.

Source nodes are no longer welded (`--v2-node-weld` is removed: it regressed
скала1 and тест 5). Instead the contour of the affected region is rebuilt
(`assembly::cracks`, `--v2-crack-width`, default 0.01 m):

1. every contour edge is split at contour nodes lying on it (within ten
   times the precision); a piece traversed twice has material on both sides
   and is dropped, the remaining pieces are relinked into simple contours;
2. an excursion of a contour between two unconnected nodes within the crack
   width (path longer than twice their distance, enclosing void of mean
   width within the crack width, chord crossing no element) is cut and its
   mouth becomes one contour point; a mouth wider than the crack width is
   allowed only for a zero-width crack, and a thin fin of material (its loop
   has the material orientation) is never cut;
3. consecutive contour vertices identified by a stacked-wall alignment
   collapse to one vertex.

Enclosed narrow voids of positive width stay with the collapsed-opening
policy, which keeps their nodes. Every rebuilt region is reported in
`topology.cracks` (pieces removed, cuts, identified mouth nodes, removed
nodes, widest mouth).

тест 6 before/after: plane patches not assembled 20 -> 0, unresolved source
elements 12461 -> 0, trial mesh topology valid and external-mesher gate pass
(before: blocked). 23 regions rebuilt: 162 overlapping pieces, 13 cuts with
mouths up to 3.0 mm. скала1, типовая секция, тест 5 have no crack and are
unchanged. The whole тест 6 run takes about 80 s (`TOPO_TIMING=1`).

Tests: a wedge crack closed and a 20 mm slit kept, a non-matching interface
removed while a real notch stays, a 5 mm fin of material never cut, each
over rotations, translations, scales and renumbering.

## тест 6 residuals closed, 2026-09-28

After the crack rebuild тест 6 assembled completely but failed the strict
audit (58 unrepresented junction segments on 5 pairs, 1 coplanar overlap,
55 unshared point contacts, 2 unshared bar-surface contacts). General fixes,
each with a regression test that fails without it:

- a vertex already on a junction line within the minimum edge of the
  required junction end becomes that end (it moves onto it) instead of a
  second vertex 2.9 um away that rejected the whole junction (56 segments
  and 54 point contacts along one wall/slab line);
- a junction chain end just outside a surface (a wall corner held by two
  walls 30 um beyond a slab edge) bends the nearest boundary edge through
  it (`Model::split_edge_within`, less than the minimum edge, contours
  revalidated); boundary splits during junction insertion accept the
  detection tolerance;
- a crack mouth identified in one region applies to every region of the
  same patch: a part split off at the crack otherwise kept the old node and
  overlapped the main region by a 0.6 cm2 sliver;
- a contour edge within the minimum edge of another surface's plane and on
  its material, but beyond the detection tolerance (a wall top 0.7 um below
  a slab it was not connected to in the source), is settled onto that plane
  before junction detection (`junctions.snapped_vertices`);
- a bar node on the plane and material of a surface that does not use it (a
  column through a slab without a shared source node) is a point contact of
  that surface.

тест 6 strict audit: all failing classes 0 (was 58 + 1 + 55 + 2); 203 surface
near misses remain review items; triangles below 1 degree 41 -> 10, minimum
angle 0.00 -> 0.24 degrees. Two junction diagnostics remain (a 19.6 mm near
touch left open, one crossing split refused as a short edge); the audit
confirms the geometry they concern is connected. скала1, типовая секция,
тест 5: output byte-identical; debug and release outputs identical.

## Short edges between needed corners, 2026-09-28

All ten triangles below 1 degree on тест 6 sat at model edges of 1-9 mm
between two needed corners (for example the ends of a lower and an upper
wall 3 mm apart on one slab line), which no earlier rule touched: redundant
collinear vertices are removed and wall ends snap to contour corners, but
two corners stay. `cleanup::collapse_short_edges` (`--v2-edge-collapse`,
default 0.01 m) merges the ends of a surface edge shorter than the
tolerance into the vertex of more surfaces, shortest first, keeping every
plane; a bar node on such an edge is not merged (`bar_node_on_collapsed_edge`).
Merges are reported in `topology.short_edge_merges`.

тест 6: 29 edges collapsed (up to 9.4 mm); minimum angle 0.24 -> 3.33
degrees, triangles below 1/5/20 degrees 10/36/326 -> 0/6/271; strict audit
still passing. тест 5: 2 edges (6.7 mm), mesh metrics unchanged. скала1 and
типовая секция: no such edge, unchanged. A 20 mm step between wall ends is
kept (test).

## Review fixes: plane check of refreshed contacts, whole crack chords, 2026-09-28

An external code review found two gaps, neither observed on the fixtures
(all four outputs byte-identical after the fix):

- `bars::refresh_contacts` accepted a point contact from an earlier record
  or ownership without checking that the node still lies on the surface
  plane; the mesh would then have attached an off-plane node to the panel.
  The plane check is now mandatory for every point contact.
- the crack chord test only checked its midpoint; a chord could cross an
  element near one end. The whole chord is now clipped against every
  element (convex pieces; a nonconvex quadrilateral as the two triangles of
  its inner diagonal) and rejected if any part lies strictly inside, more
  than the precision from the boundary; running along element edges (a
  zero-width crack) is not crossing.

Tests: a node moved 1 mm off the slab loses its stale contact; a chord with
its midpoint in the void but crossing an element corner near one end is
rejected, a chord along an edge or through a corner is not, and the notch
of a nonconvex quadrilateral is void.

## PLAXIS profile, 2026-09-28

Reconstruction now targets PLAXIS 3D's sensitivity. PLAXIS intersects the
imported geometry itself (parametric geometry only since 2016; its snap
tolerance defaults to 1 mm); small gaps or overlaps break the intersection
or force a very fine mesh, triangles below its tolerance are dropped, and it
warns about edges many times smaller than the target element. PSI/LIRA
practice recommends elements of at least 0.5 m.

`scripts/check_plaxis_profile.py` measures the geometry against a target
element size h (default 0.5 m): short surface edges and bar pieces (< h/10),
sharp contour corners (< 10 degrees), narrow faces (a contour vertex within
h/10 of a non-neighbouring contour edge of its surface) and gaps (a surface
vertex or bar node within h/10 of another surface without being one of its
vertices). Five analytic tests.

New rules, both at h/10 = 0.05 m by default:

- `assembly::gaps` (`--v2-gap-closure`, 0 keeps every gap, e.g. real
  joints): a vertex within the tolerance of another surface moves onto its
  plane keeping its own planes; if still outside it merges into a contour
  vertex, slides along its planes onto a contour edge (split there), or the
  edge bends through it. Sub-millimetre touches stay with junction
  insertion; coincident vertices are merged again afterwards
  (`topology.gaps`);
- `--v2-edge-collapse` default raised from 0.01 to 0.05 m.

| PLAXIS profile | скала1 | типовая секция | тест 5 | тест 6 |
|---|---:|---:|---:|---:|
| Short edges before -> after | 0 -> 0 | 0 -> 0 | 3 -> 0 | 20 -> 4 |
| Gaps < 50 mm before -> after | 0 -> 0 | 0 -> 0 | 2 -> 0 | 117 -> 49 |
| Short bar pieces | 7 | 0 | 0 | 0 |
| Sharp corners | 1 (8 degrees) | 0 | 0 | 0 |

The strict global audit passes on all four fixtures; скала1 and типовая
секция geometry and mesh are unchanged; debug and release outputs are
identical. The 49 remaining gaps on тест 6 need a whole surface to move or
rotate (e.g. an oblique wall end lying in the plane of another wall 26 mm
beyond its end, which is held by a third wall; slabs 25 mm apart in height),
or are near-parallel surfaces; they are reported in `topology.gaps.rejected`.
The short bar pieces of скала1 are source bars (7-28 mm) joining beams.

## PLAXIS decisions, 2026-09-29

User decisions after the PLAXIS profile batch:

1. Surfaces stay in their planes: no plane-level alignment to close gaps
   (vertices only move within their own planes, as before).
2. Bars shorter than the collapse tolerance collapse
   (`cleanup::collapse_short_bars`, `--v2-edge-collapse`, 0.05 m): a whole
   short bar disappears and its ends merge; a short piece between two nodes
   of one axis collapses into the node of more connections; the axis stays
   straight, its source elements are reported in `topology.short_bars`
   and the assembly checker accounts for them.
3. Gaps close only where the two structures lie in one plane in the source
   (the vertex within the minimum edge of the other surface's plane).
   Offsets across a plane (a wall top below a slab) are kept unless
   `--v2-gap-offsets` is given. Deformation joints are modelled wider than
   the 50 mm tolerance and are never closed.

Short edge collapse now also tries the other end when one merge is
rejected.

| PLAXIS profile after | скала1 | типовая секция | тест 5 | тест 6 |
|---|---:|---:|---:|---:|
| Short bar pieces | 7 -> 0 | 0 | 0 | 0 |
| Short edges | 0 | 0 | 1 | 4 |
| Gaps < 50 mm (offsets kept) | 0 | 0 | 0 | 63 |
| Sharp corners | 1 (8 degrees) | 0 | 0 | 0 |

скала1 trial mesh: minimum angle 2.39 -> 8.13 degrees, triangles below 20
degrees 73 -> 11 (the needles came from the short bars). The strict audit
passes on all four fixtures, типовая секция is unchanged, debug and release
outputs are identical. The remaining short edge of тест 5 (12 mm, a slab
tongue touching a wall end) is rejected by contour validation in both merge
directions.

Update 2026-09-29: offsets across a plane (a wall top below a slab) are
closed by default (the wall extends within its own plane);
`--v2-keep-gap-offsets` keeps them. Parallel structures (slabs at
different levels) are never brought together: a vertex of a surface
parallel to the other one (within 0.02 rad) is not a gap candidate. Result:
тест 5 passes the PLAXIS profile (16 gaps closed, including the 12 mm
edge); тест 6 keeps 49 gaps, 21 of them between parallel structures (kept
by rule) and 28 that would need a surface to leave its plane (corners held
by several planes, edges shared with a third structure).

## Source gaps closed in the frame solve, 2026-09-29

User decision (after the corner cases of тест 6, a wall end 26 mm beyond
the corner of two other walls): such gaps are source defects and are closed
by the frame solve ("variant 3"), not by moving vertices out of their
planes afterwards.

`frame::gaps::candidates` finds virtual incidences: a contour node of its
own patches within the gap tolerance (`--v2-gap-closure`, 50 mm) of the
material of a non-parallel patch it does not belong to, off its plane by
at least 1 mm. Excluded, with the reason in the rule:

- nodes of a patch parallel to the target (within 0.02 rad): slabs at
  different levels and parallel walls are never joined;
- interior nodes of a patch (a gap is at the edge of a structure);
- nodes near two mutually parallel candidate planes (between the walls of
  a stacked offset pair: wall alignment's job);
- steps: a node of the same structure already on the target plane within
  the tolerance of the projection (the offset is a contour step or crack
  mouth to that node; pulling onto the plane would collapse it to zero).

`frame::solve_closing_gaps` adds `normal . node - offset = 0` for each
incidence to the least-squares frame solve (normals fixed: planes only
translate, they never rotate), with the usual movement budgets. Up to four
rounds drop incidences named among the largest constraint failures or at
over-budget nodes; a failed solve falls back to the base solve. The
applied and dropped incidences are reported in `frame.virtual_incidences`
(constraint origin `VirtualIncidence`). `--v2-keep-gap-offsets` and
`--preserve-details` disable it.

Assembly gap closure also handles point touches below the minimum edge (a
corner micrometres from another surface's contour edge, lower bound
10 x precision): junction insertion imprints only intersection lines.

| After | скала1 | типовая секция | тест 5 | тест 6 |
|---|---:|---:|---:|---:|
| Virtual incidences applied / dropped | 0 | 0 | 7 / 0 | 115 / 0 |
| Heights closed (mm) | - | - | 11-19 | 1.5-40 |
| Maximum frame movement (mm) | 1.5 | 0 | 0.8 -> 18.9 | 5.8 -> 39.6 |
| PLAXIS profile items | 1 | 0 | 0 | 53 -> 26 |
| Trial mesh minimum angle (degrees) | 8.13 | 20.21 | 7.15 | 3.33 -> 7.09 |

All four pass the strict extended audit and the assembly checker; debug and
release outputs are identical (тест 5). The corner of тест 6 (wall ending
short of two walls, 26 items) is closed. тест 6 residuals: 10 gaps of
49.97 mm between parallel walls (kept by rule, at the tolerance), 4 of
49.6 mm on a slab at a parallel slab edge, a 19 mm offset between two
nearly coplanar walls (0.2 degrees) with its short edges, a 25 mm step at a
wall corner and one 47 mm gap at a 53-degree wall.
