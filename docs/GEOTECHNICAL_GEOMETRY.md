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
