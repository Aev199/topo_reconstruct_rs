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
