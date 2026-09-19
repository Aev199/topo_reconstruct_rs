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
