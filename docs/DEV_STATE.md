# Development state

Updated: 2026-10-04
Repository baseline: geotechnical mode only, editor application (see git log).

This file is intentionally short. It is the entry point for the next development
session; detailed rationale belongs in `docs/GEOTECHNICAL_GEOMETRY.md` and
`docs/RECONSTRUCTION_V2.md`.

## Product objective

Produce valid, connected structural geometry from imperfect FE input for:

- reliable remeshing / geometry import in PLAXIS;
- a quality conforming mesh for MIDAS;
- later transfer of structural properties, loads and provenance.

Exact reproduction of the source FE tessellation is not the objective.

## Private fixtures

Private LIRA text models used as the Tier C gate (never committed):
`скала1` (right section, 59 surfaces / 525 axes since planar panel
growing), `типовая секция` (23 / 291), `тест 5` (108 / 8), `тест 6`, and
since 2026-09-30 eight more: a test slab, two АЖТ tests, скала with seismic
loads, Багратион concrete bedding, "для testa", ЖК Остров (piles, curved
walls) and Багратион v4 (160 MB).

## Current verified state

Surface-surface junctions are explicit shared topology
(`assembly::junctions`, surface `embedded_edges`), and both meshes use one
subdivision of every junction edge. Verified on all three private models:

| Global surface audit | скала1 | типовая секция | тест 5 |
|---|---:|---:|---:|
| Unrepresented junction segments (before → after) | 520 → 0 | 400 → 0 | 2194 → 0 |
| Segments lacking shared mesh edges | 455 → 0 | 292 → 0 | 2194 → 0 |
| Audit passed | yes | yes | yes |

Everywhere: 0 invalid surfaces, 0 coplanar overlaps, 0 unresolved source
elements, trial topology and external-mesher gate pass, surface/axis counts
unchanged, 0 reconciliation problems. Debug and release outputs are identical.

Geotechnical assembly trims thin consoles beyond junction lines
(`assembly::consoles`; тест 5: 26.61 m² of 0.10–0.19 m slab consoles) and
aligns stacked walls to the axis of the wall below within 50 mm
(`assembly::stacking`; тест 5: 8 walls, 25 mm), closes wall ends and removes
redundant collinear vertices at short edges within 50 mm
(`assembly::cleanup`), audit still passing. Tolerances are CLI options:
`--stack-offset`, `--wall-end-snap`, `--console-width`.
Duplicated vertices and wall ends next to contour corners are merged,
carrying bar axes along (`assembly::cleanup`, `Model::merge_vertices`).

Mesh quality (below 20°): скала1 163 → 73 (min 2.39°), типовая секция 0,
тест 5 28 (min 6.98°, none below 5°). The тест 5 baseline (53) was not comparable: its
slab meshes ignored 2194 junction segments.

## Current blockers

All four private fixtures pass the strict extended audit, including the
mesh-level connectivity and missing-property checks. Reconstruction targets
PLAXIS (`docs/GEOTECHNICAL_GEOMETRY.md`, "PLAXIS profile"):
`scripts/check_plaxis_profile.py` is the readiness gate at a target element
size of 0.5 m. типовая секция and тест 5 pass it. Remaining:

Scoreboard 2026-10-02 (12 fixtures; assembly checker, strict audit, trial
mesh valid / external-mesher ready, PLAXIS profile): скала1, типовая секция,
тест 5, тест 6, the test slab, both АЖТ tests, скала seismic and "для testa"
pass as before (byte-identical where no rule applies). Open:

- ЖК Остров: frame accepted (default step, 11 min), all surfaces and bars
  built, assembly checker passes, trial mesh valid and mesher-ready.
  Strict audit: 26 coplanar overlaps of 1-3 mm2, 23 unshared point
  contacts (0.2 um - 0.9 mm), 4 unrepresented intersections;
- Багратион bedding: mesh valid, 18 bars in a slab plane without contact,
  1 overlapping bar, trial mesh angles (mesher gate false). Closed crack
  voids keep only shared nodes (2026-10-02): triangles under 1 degree
  195 -> 99, sites 154 -> 63;
- Багратион v4: frame accepted on the relaxation step, assembly checker
  passes, every surface built, 1 of 40 600 bars rejected
  (axis_direction_conflict), trial mesh topologically valid (2.1 M
  triangles). Micrometre vertex pairs removed (interval ends snap to
  contour vertices, micro-overlaps of coplanar slabs close): no triangle
  edge under 0.1 mm, 51 triangles under 1 degree (minimum 0.037), all at
  millimetre features (0.45 mm - 1 cm), as for ЖК Остров and Багратион
  bedding.

Rules added 2026-10-01/02 (geotechnical; deviation from the source is
allowed, user decision 2026-10-01): over-constrained panel merges only as a
frame relaxation step (`relaxation.frame` in the report); bar-bar crossings
share a generated vertex; bar ends within 50 mm of another bar's span join
it (`bar_tees`), merging into a bar node they project onto; bars lying on
each other exchange nodes and same-stiffness contained bars are removed as
duplicates; points within ten precisions are one point (interval ends, bar
nodes); short bars may shorten within their end budgets; bar nodes slide
onto an adjacent crossing.

Rules added 2026-10-02 (Багратион v4 mesh): a bar end merging into the
bar node it projects onto collapses the bar piece joining them; a node of
one bar next to a bar crossing slides onto it (bars ending there follow);
a crossing next to another generated crossing is still shared; an invalid
contour with a hanging node of the source mesh closes its seam (crack
width, filled-area budget; reported as a `seam` simplified hole); mesh
promotion of Steiner vertices and edge/bar chain synchronization use the
point slack (ten precisions) like the assembly.

Rigid links (user decision 2026-10-02): geotechnical mode removes bars
whose numeric stiffness gives a radius of gyration sqrt(EI/EF) beyond
their length and 1 m (LIRA "1000 200000 200000 200000"; Багратион v4
36 011 elements, Багратион bedding 2143, none elsewhere). Point slack is
10 um. Other rules for the reviewed defect classes use the 5 cm gap
tolerance (variant A).

State (2026-10-02, latest):
- All 12 fixtures: trial mesh valid and mesher-ready; strict global audit
  passes on all but ЖК Остров.
- Багратион v4: 7 sharp-triangle sites left (2 wall/slab corners 2-7 mm
  off a beam: moving the corner would leave the wall plane).
- Багратион bedding: strict audit passes (contour pinch from an exact
  edge split fixed); 12 sharp-triangle sites inside slab contours.
- ЖК Остров: 25 coplanar overlaps (single elements of warped wall corners
  folded onto each other by the frame and gap closures; containment shows
  only after gap closure, so removing them needs surface removal after
  gaps), 19 unshared point contacts and 4 unrepresented intersections in
  one knot of walls (surfaces 566-569, 595, 596).

Automation closed (2026-10-04). Final rules of the batch: numerical
planarity 1 um (`assembly::Policy.precision`), point slack ten precisions
(10 um), contacts on a plane within five precisions (as the audit);
surfaces covered by, or narrower than twice their tolerance and partly
folded onto, a surface of their plane and stiffness are removed
(`removed_slivers` reasons `covered` / `absorbed`); junctions are
inserted again after generalization; hole rings collinear in the source
are seams, holes the frame closed to no area under the minimum opening
are filled.

Final scoreboard (strict global audit; <1/<5/<20 degree trial triangles):
all 12 fixtures mesher-ready; strict passes on 11 of 12 — ЖК Остров now
passes (8/53/6830; was 25 overlaps, 19 point contacts, 4 intersections),
Для testa (0/4/260), Багратион bedding (40/112/333), m1, m2, m5, m6,
test slab, АЖТ x2, скала seismic. Багратион v4 (11/58/1263) keeps 10
residual items: 2 coplanar overlaps, 1 point contact, 1 bar intersection,
6 bars in a surface without contact, and one 25 mm unlinked joint (kept
by rule) — candidates for the manual editor.

Tolerance ladder (user decisions 2026-10-03): section sizes change by
decimetres in the target model, so the simplification tolerance of a
surface is half its plate thickness (LIRA `GEI E nu H`), at least the 5 cm
gap tolerance and at most 0.2 m (`--simplification-cap`); a shared
chain takes the tolerance of its thinnest surface. Free openings narrower
than 1 m (smaller side of the minimum bounding rectangle,
`--min-opening`) are filled when no surface or bar is attached and
nothing passes through (`topology.filled_openings`). Expansion joints are
never closed: an in-plane gap wider than the crack width with constant
width along twenty widths of contour stays open, unless the analysis
model ties the two surfaces (bars, rigid links, two-node links such as
LIRA type 55, rigid bodies of block 25): such structures are joined
(Для testa: 2 cm slab/balcony gaps on type-55 links and 2 cm bars).
Generalization keeps clear of bars and junction lines: no new wedge under
20 degrees, no approach nearer than before without touching. Consoles are
trimmed again after generalization.

Ladder results (<1/<5/<20 degree trial triangles, vs 1d57ec3):
Багратион bedding 71/259/2913 -> 40/113/334, strict passes, 4 openings
filled; Багратион v4 10/54/1289 -> 11/56/1270, overlaps 8 -> 2, point
contacts 29 -> 1, 138 openings filled, one 25 mm unlinked gap kept as a
joint; Для testa 0/4/195 -> 0/4/259, strict passes, 96 openings filled,
balconies joined; m6 0/0/204 -> 0/0/507 (min 7.1 -> 7.0); m5 0/0/25 ->
0/0/79; ЖК Остров 11/185/5701 -> 12/98/6970, 660 openings filled,
PLAXIS gaps 95 -> 61, 25 overlaps unchanged; m1, m2, n1-n4 unchanged.

Geometry generalization (2026-10-03, `cleanup::generalize_contours`, last
assembly step, after bar imprinting): every chain of contour/junction edges
between fixed vertices (branch points, changes of the user-surface set,
bar anchors and ends, contacts, retained hole nodes) is simplified once by
Douglas-Peucker within the 5 cm gap tolerance, so all surfaces sharing it
stay conforming. Vertices are only removed (`generalized_contours.removed`
with source node and deviation; the assembly checker verifies them). A
removal is skipped when the shortcut would sweep over another edge, bar or
vertex, or close a wedge under 20 degrees with a bar (at a kept vertex, or
a bar within 5 cm) narrower than before.

Generalization A/B (1d57ec3 vs this batch, frame caches, <1/<5/<20 degree
trial triangles; strict audit unchanged unless noted):

| Model | before | after | notes |
|---|---|---|---|
| Багратион bedding | 71/259/2913 | 34/109/343 | PLAXIS narrow faces 12 -> 1, gaps 51 -> 15 |
| ЖК Остров | 11/185/5701 | 11/87/7822 | PLAXIS gaps 95 -> 62, narrow 5 -> 0; 25 overlaps unchanged |
| Багратион v4 | 10/54/1289 | 11/54/1279 | strict fails on both: overlaps 8 -> 2, point contacts 29 -> 1, bar-in-surface without contact 14 -> 6 |
| для testa | 0/4/195 | 0/4/258 | |
| скала seismic; m5; m1 | 0/0/11; 0/0/25; 0/0/11 | 0/0/11; 0/0/84; 0/0/13 | |
| m6 | 0/0/204 | 0/3/525 | min angle 7.1 -> 3.0 |
| m2, test slab, АЖТ x2 | unchanged | unchanged | |

Багратион v4 at 1d57ec3 re-run with the current frame cache (relaxation
step 1) fails the strict audit, unlike the earlier record below.
The growth of 5-20 degree triangles is a trial-mesh effect, not invalid
geometry: long generalized edges are subdivided by the refinement of the
neighbouring surface (a Steiner point 26 mm from a wall corner, because a
slab vertex lies 5 cm beyond the wall end). Graded subdivision of shared
edges in the trial mesh, or the external mesher, addresses it.

Next for generalization: snap a contour chain lying within tolerance and
parallel to a bar or another contour onto it (instead of leaving a 4-5 cm
strip); bars as straight axes between real joints; panels absorbing warped
corner elements (ЖК Остров overlaps; covered-surface removal is still in
`git stash` as an alternative).

Rules of this batch: slab edges kinked off a beam (one neighbour on it,
within the recognition angle) straightened; a surface vertex next to a
beam node merges into it; an edge split never pinches a contour; a vertex
touching another surface's contour edge splits it.

Speed: frame retries continue one LSQR run (ЖК Остров frame 34 -> 11 min),
`--frame-cache PATH` reuses a solved frame (development), trial mesh
contacts indexed by surface (Багратион v4 mesh 9 min -> under 1 min),
`TOPO_DIAG=1` prints mesh failure details.

Previous fixtures:

1. тест 6: PLAXIS profile passes (compound alignment corners closed);
2. скала1: one contour corner of 8 degrees (a real wedge-shaped stiffness
   zone, kept by user decision);
3. trial mesh quality on тест 6 (minimum 7.09 degrees); the final mesh
   should come from Gmsh or the target program (`docs/PRIOR_ART.md`);
4. two junction diagnostics on тест 6 (a 19.6 mm near touch, one refused
   crossing split).

Fixed after тест 6: junction ends micrometres from existing vertices or
surface edges, near-planar wall tops, crack mouths shared across a patch,
columns through slabs without shared nodes, cracks of the converted mesh (contour rebuild, not
node welding: `--v2-node-weld` is removed), duplicate edge keys after splits (released edges stay
indexed), a stacked-alignment identification that ignored other supports of
the upper node, unbuffered report output, slow console trimming.

## Review of 9203e05..8c0f264 (2026-10-04)

Fixed in 762f635: generalization could leave a removed vertex up to
twice its tolerance from the final contour (and reported the distance to
the Douglas-Peucker line); `covered` was tested by vertices only (a
surface bridging a notch of a nonconvex one was deleted); joint length
counted a contour wholly within the band twice; reports before the
removal of covered/absorbed surfaces keep the old surface numbering
(`topology.surface_renumbering`); synthetic tests ran at 0.1 um while
production uses 1 um (`assembly::PRECISION`). Tier C: no metric change.

User decisions 2026-10-04 (18d3ec5, 63efb20, 0f63556):

- `absorbed` only on an overlap the global audit would count (more than
  five precisions along the shorter contour), with the lost part within
  the piece's simplification tolerance of the absorbing surface; looked
  for again after generalization; chains allowed. Records give the
  absorbing surface, overlap, lost area and distance.
- Generalization never nears an edge or bar without touching it, except
  an edge strictly inside the material of the chain's own surface (an
  existing overlap, which nearing shrinks).
- The audits cap the precision taken from the report at 1 um.
- Openings longer than 3 m are kept (`--max-opening-length`).
- Geotechnical mode only: V1 and `--v2-preserve-details` removed;
  `pipeline::run` with `Profile::plaxis()`; CLI `model.txt -o out.json
  [--mesh] [--frame-cache PATH]`.

Tier C (0f63556 vs 8c0f264; <1/<5/<20 degree trial triangles, failing
audit items, PLAXIS items):

| Model | before | after |
|---|---|---|
| Багратион v4 | 11/58/1263, 10, 59 | 11/55/1257, 10 (same sites), 57 |
| ЖК Остров | 8/53/6830, 0, 86 | 6/50/7017, 0, 93 (gaps 50 -> 57, 27 long openings kept) |
| Багратион bedding | 40/112/333, 0, 63 | 33/103/335, 0, 63 |
| скала1, скала seismic | 0/0/13, 0, 1 | 0/0/11, 0, 1 |
| тест 5, тест 6, Для testa, скала2, test slab, АЖТ x2 | | unchanged |

Strict audit and assembly checker pass on 11 of 12, as before.

The in-application audit (`audit`, Rust) agrees with the Python
auditors class by class on all 12 fixtures, except 6
`bar_in_surface_without_contact` items on Багратион v4: bars leaving a
slab at 0.13 degrees, where the Python audit's 5 um buffer stretches the
in-surface part 2.2 mm beyond the contact interval; the Rust audit clips
exactly and does not count them. Open for a user decision: treat them as
a grazing artefact of the Python audit (then v4 has 4 real residual
items) or keep them.

## Editor (2026-10-04)

Windows application (user decision): Rust core in process, three.js UI
in WebView2 (Tauri 2, `app/`). `service::Service::dispatch` is the single
command layer (the Tauri command `call` and `examples/editor_server.rs`
for a browser). Edits (`assembly::edit`, `editor::Session`): move/merge a
vertex, delete/join surfaces (one plane and stiffness), split an edge,
gap as joint / close it; transactional, journaled with a replay check,
re-audited at once; undo replays the journal; a project stores input,
hash, profile and journal; the saved result carries `user_edits` and
`edits`, accepted by the assembly checker. UI checks:
`app/ui/tests/e2e.mjs`, `gap.mjs` (Playwright through the bridge).
Windows installer: manual `build_app` job (one windows-latest job, NSIS,
artifact 7 days), first build green (run 37190104074).

Portable build (user decision 2026-10-04): one `topo-editor.exe`, data
(frame cache, WebView2 profile) in `topo-editor-data` beside it; needs the
WebView2 runtime (part of Windows 11 and updated Windows 10).

Speed (2026-10-04, Багратион v4, 4 cores): frame 25 min -> about 1 min
(LSQR over active unknowns only, CSR/CSC, parallel products, blocked
norms instead of a sequential hypot; 2/3 of the old time went to the
default and snapped solves that never converge before the relaxed one),
assembly 8 -> about 3 min (vertex removal searches only the chain's
surfaces and validates only the edited ring). Whole run without cache
33 -> about 5 min. The frame of ill-conditioned models depends on where
LSQR stops (ЖК Остров: nodes up to 6 cm apart between two valid
solutions); a final touch pass makes the assembly robust to it. Tier C
after the change: strict audit 11/12 as before, v4 the same 10 sites.

External audit of 825f963 (fixed): opening a model/project with unsaved
edits asks save / discard / cancel (`summary.dirty` against the journal
last saved or fully replayed); a project stores the hash of the bytes
actually reconstructed and saving reports a changed file on disk; the
Tier C runner takes the PLAXIS verdict from `plaxis.json.passed`; the
Python profile lists `user_edits.accepted_joints` separately
(`accepted_joints`, `accepted`), not as gaps; the UI shows two verdicts
(geometry/connectivity, PLAXIS profile) and accepted joints; empty scene
has finite bounds; Russian journal labels; materials disposed.

Editor sufficiency (2026-10-04, checked on Багратион v4 with
`examples/edit_probe.rs`, which applies journal edits to a real model and
prints the audit change): all 4 failures close in 4 clicks — the two
coplanar overlaps (a stiffer zone whose contour skips a vertex of the
slab hole, or detours 2 cm) and the point contact by `merge_vertices`, the
bar intersection (bar bent 9 um at a shared node) by `connect_bars`; the
edited report passes the assembly checker and the strict audit except the
6 grazing `bar_in_surface_without_contact` items (open question). New
edits: `connect_bars` (shared node; a node slides along its bars, a bar
end is drawn onto the other bar, never bent), `connect_bar_to_surfaces`,
`connect_surfaces` (junction of a pair as shared edges), `delete_bar`,
`merge_vertices` of consecutive bar nodes collapses the piece (provenance
in `user_edits.removed_bars`), `move_vertex` slides an interior bar node
along its bars; edited bars are listed in `user_edits.edited_bars`. Not
yet closable in the editor: gaps between near-parallel planes (3.7 cm
wall foot / slab bottom on v4: needs a "move surface to plane" edit),
short edges where a merge would create another short edge.

PLAXIS 3D export (2026-10-04, user decision: directly through the Python
API, no MIDAS yet): `plaxis::exchange` writes `topo-plaxis-1` (plates on
planar hole-free polygons, beams per stiffness run of a bar, elastic
materials from LIRA block 3: GEI E/nu/H/RO, S0 b x h cm with RO per length;
forces x9.80665 t -> kN, configurable; stiffness types without a material
are listed). Holed surfaces: cuts along v (or u, the better variant) from
the lowest/highest points of each hole's extreme columns to the first
contour, cuts within 5 mm of a vertex go to it when clear, cuts between
pieces sharing only the cut are removed; area check falls back to CDT
triangles. Багратион v4: 4041 plates -> 5623 polygons (226 holed), all
valid and planar, areas exact, 9 short edges / 7 sharp corners added by
cuts; 6464 beams; 3 stiffness types without material. Loader
`scripts/plaxis_export.py` (plxscripting: gotostructures, platemat/beammat
+ setproperties with fallback property sets 2022+/older, surface + plate,
line + beam, `.Material`), recorder tests; the editor runs it with the
PLAXIS Python (`run_plaxis`). Bentley documentation hosts are blocked by the
network policy: command names come from PLAXIS command logs quoted in
search results and plxscripting 1.0.4 from PyPI. Not verified against a
real PLAXIS 3D 2022 yet (Windows + licence needed).

Re-audit of 9419b92 (2026-10-05), fixed:
- R1 effective stiffness: the parser keeps the numeric EF/EIy/EIz/GIk of
  a bar type and WLKE/PLKE/WLKG/PLKG of a shell; export mode `effective`
  (default) gives A = EF/E, I3 = EIy/E, I2 = EIz/E and, for shells, an
  equivalent thickness and E keeping membrane and bending stiffness and
  the weight; `nominal` uses b x h and the GEI values. What cannot be
  transferred (GIk, WLKG/PLKG) is listed per material.
- I2/I3 and orientation (PLAXIS Reference Manual: section height along
  local axis 2, I3 = width x height^3 / 12): local axis 2 = LIRA Z1 by the
  LIRA default rule (upward in the bar's vertical plane, global X for a
  vertical bar); the loader sets the line's AxisFunction Manual / Axis2
  for rectangular sections and reports when PLAXIS refuses. LIRA rotation
  angles of sections are not read (warning in the file).
- R3 property sets from Bentley's table: plates Identification/d/
  Isotropic/StructNu12/Gamma (V22.02+), D3d (V22.00), MaterialName/d/
  IsIsotropic/Nu12/w (V21); beams CrossSectionType "User-defined", A/I2/I3/
  E/Gamma (V21: BeamType, Iyy/Izz, w). `setmaterial` as PLAXIS logs it. A
  failed run deletes the objects it created.
- R4 materials come from the opened input bytes, not the file on disk.
- R2/R5/R6 bars: `bars::axis_defects` (node off the straight bar, nodes not
  strictly ordered or closer than the minimum edge, ends not nodes at 0/1,
  empty spans) is a failure in the Rust audit (`broken_bar`) and the
  session refuses any edit adding one; a bar end slid along its line keeps
  its spans and refuses passing an interior node; connect-to-surfaces works
  on the chosen bar with every bar guarding its nodes and also splits
  contours a bar lying in a surface's plane crosses; connect-bars slides an
  interior node onto the crossing. The assembly checker checks node order
  of edited bars. Zero-length beam pieces are skipped with a warning.

First real PLAXIS run (user, 2026-10-05, скала3 and тест 5): "Cannot
intersect while the geometry contains invalid objects" for polygons whose
first three points were collinear (straight contour runs through junction
vertices): PLAXIS fits the plane through the first points ("Define plane:
First points") and called them not coplanar (0.45 m). Fix:
`plaxis::plaxis_polygon` drops straight-line vertices (PLAXIS intersects
and recreates junction points) and starts each ring at the corner with the
largest first triangle. Simulating PLAXIS's check on 9 fixtures: 0
degenerate starts, deviation from the first-points plane at most 0.7 um
(тест 5 had 56 of 127 polygons with a collinear start before).

Connectivity diagnostics (2026-10-05, `audit/connectivity.rs`, Rust audit
only): `floating_group` (a connected group of bars/surfaces apart from the
main structure; PLAXIS class for bars only, review otherwise), `free_bar_end`
(review; legitimate for pile tips and cantilevers), `lost_bar_link` and
`lost_surface_link` (the source model shared a node between two bars or a
bar and a shell patch, the result has no shared vertex/contact: failure when
the two sit in different connected groups, PLAXIS item otherwise).
Source links come from the frame (`assembly::Report.source_links`).
Findings carry repair proposals (`fixes`, editor JSON): connect bars,
merge a bar end into a vertex, move a bar end onto a surface, connect a
bar to surfaces, delete a floating group (`delete_bars`, one journal entry).
New editing behaviour: `connect_bars` seats an end on a parallel bar.
Tier C: no lost links on the 13 fixtures (the floating groups of
«Для testa» — 1, 40 and 40 bars — are separate in the source too).

Real PLAXIS round 2 (user, 2026-10-05): тест 5 imports and meshes. скала3
(new private fixture, 13 now): a floor slab (z 7.33, 2956 elements) and 6
more regions (4797 elements) were missing while every audit was green.
Cause: the frame did not close (residual 0.3 mm): a merged wall family
(26+27, nodes 25 mm off its fit, accepted by the panel tolerance) leaning
1.6e-4 was not snapped vertical because the snap rule compared with the
1 cm plane distance; three walls then met at one point instead of their
common vertical edge, the edge nodes would move metres, and the regions
were dropped (`support_intersection_or_movement_budget`). Fixes: a family
snaps when snapping worsens its fit by at most the plane distance
(regression test fails without it); unbuilt regions and bars are audit
failures (`surface_not_built`, `bar_not_built`, Rust audit and session)
and fail `check_v2_assembly.py` (it used to count their elements as
accounted); gaps are closed once more after the later stages (generalized
contours, junctions). The Python global audit recomputes a coplanar
overlap on a grid of tolerance/1000: GEOS returned a 6.4 m2 phantom
overlap between two walls touching along a line (скала3 surfaces 6/7).
Tier C, 13 fixtures (fresh frames): скала3 85 surfaces, all regions built,
strict audit and PLAXIS profile pass; ЖК Остров PLAXIS items gap 60 -> 44,
short edge 27 -> 24, sharp corner 4 -> 5; all other fixtures unchanged
(strict 12/13, v4 the same residuals).

## Load transfer (2026-10-05)

Static loads are parsed from the source LIRA file (`parsers/loads.rs`) and
mapped onto the reconstructed geometry (`loads.rs`): node forces → point
loads, bar and plate-edge loads → line loads along the bar axes, plate
pressures → surface loads on hole-free polygons (whole surface when all its
elements carry the value, otherwise the region of the elements clipped to the
surface contour; at most 40 value groups per surface). Plate stamps (code
5/15 point forces) are clustered by case, direction, level and force and
spread evenly over their elements. Cross-checked against the user's
Lira_Midas-converter: the sign of forces is reversed (positive LIRA value acts
against the axis; verified by the wind case names X±/Y± and gravity cases);
rows of doc 6 without a case number belong to the previous case (one fixture
has 1.5 M such rows). Per-case resultants of source and exported loads agree
(1e-5 on скала1/2/3, up to 4 % on тест 5, тест 6, ЖК Остров where clipping to the
contours loses small areas); loads on plate elements not in the geometry are
reported separately. Skipped (static settlement task): thermal, dynamic,
stage codes 8/88, prescribed displacements, plate moments, arbitrary
trapezoids on plates. Loader: `pointload/lineload/surfload`, phases per case —
not verified against a real PLAXIS yet (command and property names from the
documentation).

## PLAXIS load combination (2026-10-05)

`Settings.combination` (`loads.rs`): the cases of the source become ONE
case (`COMBINATION`) with a factor per case chosen by the user (cases
without a factor, the self-weight case first, are dropped; PLAXIS applies
the self-weight itself). With `Simplify` the combined loads are reduced so
that the PLAXIS geometry gets no new contours and the resultant stays:
plate pressure → one uniform load over the whole plate when the loaded part
is at least `min_fraction` (30 %) of it and the centre of pressure is within
`center_tolerance` (15 %) of its size from the plate centre, else point
loads at the centres of pressure of the connected loaded parts; bar loads →
uniform over the whole bar or point loads at the segment centres; plate-edge
lines stay lines between model vertices, otherwise point loads; node loads
by node. Check on the fixtures: resultant of the combination equals the
weighted source (exact within 1e-5); many partly loaded plates become point
loads (e.g. тест 5: 98 of 99 plates) — the tolerances are in the dialog.
UI: dialog table of cases (use / factor; self-weight and seismic cases off by
default), mode «сочетание» / «по загружениям, исходные контуры». Not yet
verified in a real PLAXIS.

## MIDAS export (2026-10-05)

Format decision (user): MIDAS Civil `.mxt` as written by Lira_Midas-converter
(that converter has no FPN output); gmsh embedded in the exe.
Chain: `meshing::mesh_state` (Gmsh C API via libloading, `gmsh.rs`; shell
elements per surface with shared boundary curves, bar pieces between
anchors, vertices on an edge's inside split it, unused nodes dropped) →
`loads::transfer` with `combination: None` and `cases` = chosen cases →
`mesh_loads::transfer` (pressure by the exact area fraction of the contour in
each element, bar lines on bar pieces, plate-edge lines to boundary nodes,
point loads to the bar/nearest nodes) → `midas::write_mxt` (*NODE, *ELEMENT
bars then plates, *MATERIAL, *SECTION VALUE with torsion J of the rectangle,
*THICKNESS, *STLDCASE, *USE-STLD with *CONLOAD/*BEAMLOAD/*PRESSURE, names
transliterated and made safe as in the converter). Service command
`export_midas`, button «Экспорт в MIDAS…» (size, quads, cases table; self-weight,
stages and dynamic cases off by default; construction stages are modelled as
load cases in LIRA and are not transferred to PLAXIS or MIDAS).
Checked on the fixtures: no duplicate nodes, resultants on the mesh match the
geometry loads (≤ 1.2 %), скала2 → 14.3 k nodes, 25 k plates, 3 k bars,
13 cases, 5.8 MB. Gmsh library: `TOPO_GMSH_LIB`, beside the exe, or embedded
(`GMSH_DLL_PATH` at build time; the CI app job fetches the SDK's
gmsh-4.15.dll, 89 MB, and embeds it; extracted to the temp folder on first
use). Gmsh is GPL: the embedded build is for internal use or must be released
under GPL-compatible terms. Not verified in MIDAS: beam local axes (beta
angle 0 as in the converter), supports (not exported), quad ordering.

## Cutting off upper storeys (2026-10-05)

`reconstruction/assembly/cutoff.rs`: `floors(state)` (horizontal surfaces
grouped within 0.15 m; `major` = at least 20 % of the largest floor, landings
and pits are not listed), `State::cut_above(z)` = journal edit `cut_above`
(undo, project replay, editing before and after work as for any edit).
Vertices within the shortest edge of the level are on it; an edge crossing
farther gets one new vertex shared by everything that uses it; surfaces
crossing the level are clipped (faces of the planar graph of the edges below
plus chords along the level inside the material; holes below the level go to
their face, embedded edges stay) and re-added with `Model::add_surface`, so
edge identity with the neighbours is kept; bars are trimmed (new anchor or
an existing one near the crossing); everything above is removed with
provenance (`removed`, `removed_bars` reason `cut`). `State.cut` keeps the
lines of the walls and points of the columns at the level (supports) and
the top elevation. Tried on all large fixtures at every major floor: no new
audit failures (v4: 4 before, 3 after); Для testa 11, Остров 31, скала3 4,
тест 5 4 levels cut without errors.

Export (`loads.rs`, `storeys.rs`): every load on an element above the
level (plate pressures, edge lines, bar loads, node loads) and the weight of
the elements above (pseudo-case `CUT_WEIGHT_CASE`, from density x thickness
or RO x length, tf) go to the support nearest in plan (line load along a
wall, point load on a column); the overturning moment of the items about
the supports is kept by a couple of vertical forces on the supports (least
squares; the vertical-axis moment is not kept). Resultants of every case
equal those of the uncut model (checked on Для testa at floor 8). The cap
slab (horizontal surfaces at the level) gets a bending stiffness for the
plan sections of the removed walls and columns: `EI' = EI h / H k` (h storey
height below, H height removed, k the user's factor), plate of the cap's plan
size `t = (12 EI' / B)^(1/3)` per axis, geometric mean; only the plate factor
PLKE grows, membrane stiffness and weight stay (typically a few metres of
equivalent thickness: tune k). Export options `cut: {loads, cap, factor}`
in `export_plaxis` and `export_midas`; UI button «Этажи…» and options block
in both export dialogs. Not verified in PLAXIS or MIDAS; the equivalence is
an engineering approximation to be calibrated by the user.

### Audit of 073700f: fixes (2026-10-06)

- Loads: bar edge loads follow the element's local node numbers (L1); a
  combination of pressures of both signs is kept as force + couple (L2,
  `sign_groups`); crossing wall/column loads and weight are clipped by the cut
  and supports are vertex ids, so they follow edits (C2, C4, `live_supports`);
  the moment of removed loads is kept by a least-squares couple about the
  supports (C3); modular ratio is linear in the column inertia (C5).
- Cut: new vertices through `Model::split_edge` and shared-vertex reuse, no
  coincident unshared vertices (C1, test `a_cut_leaves_no_coincident_vertices`);
  vertex renumbering after cuts keeps `connected` consistent (G1, unit test
  `renumber_pairs`; no end-to-end repro available).
- Mesh loads: exact area fraction instead of centre tests; a patch well off the
  centre of a large element goes to the nodes by its centroid (keeps moment);
  partial bar loads clip the piece on the line of the load (N1); point loads
  use barycentric weights (N3); tests in `tests/cutoff.rs`.
- EI = 0 of a bar type is read as a value and floored with a note (M1).
- The new check found a regression at once: `load_resultant` for surface loads summed |fan triangles|, overcounting concave contours 4-5x (скала3, тест 5); now signed (unit test with an L-shape). After the fix: Для testa/скала3 forces within 0.3 %, тест 5 within 2.4 % (СВ 0.6 %, полезная 1.2 2.4 % are flagged by `problems`).
- Honest check: every case reports force and moment (source → geometry → mesh,
  `Report::problems`, tolerances 2 % force / 5 % moment of the moment scale);
  `export_plaxis` and `export_midas` return `load_problems`, the UI shows an
  error status «НЕ ГОТОВО» instead of success when any is non-empty.
- Still unverified: real import into PLAXIS 3D 2022 and MIDAS Civil, comparison
  with an independent reference model (user side).

### Re-audit of e8e4693: fixes R1–R7 (2026-10-06)

- R1: centres of source loads, removed weight and stamps are centres of area
  (`area_centroid`, signed fan), not the vertex average; for an element across
  the cut level the source tally uses the centroid of the kept part
  (`kept_part`), so the reference no longer shares the error of the export.
- R2: a wall in `live_supports` needs a non-horizontal surface with an edge on its
  line (the top edge of a deleted wall stays on the slab).
- R3/R5: a line load along plate edges is carried by every edge on its line (grid
  index of distinct shell edges) as consistent nodal forces of the linear load
  over the covered part: partial and triangular loads keep force and moment;
  what lies over no edge goes to `lost`.
- R4: the whole-coverage shortcut of pressure patches applies to convex contours
  only; concave ones use the exact intersection.
- R6: a covered part more than 1 % of the diagonal off the element centre goes to
  the nodes by its centroid (moment kept). Smaller offsets remain a constant
  pressure (moment error ≤ 1 % of the diagonal × force).
- R7: plates outside the geometry are tallied with force, moment, sum of |F| and
  sum of |F||r| (`not_in_geometry_*`): a couple that sums to zero is a problem.
- Known approximation: trapezoidal pressure on a plate replaced by its mean
  (counted in `approximated`; its source tally is the same mean).

### First real PLAXIS import (Остров, cut) (2026-10-06)

User report: very slow transfer and `Lines overlap: Line_112..435` at gotostages/gotomesh,
a dense field of point loads on the first slab above the foundation.
- Cause 1: line loads of different cases (and of the combination, which only relabels
  the cases of the loads made while reading) lay on the same lines as separate objects.
  Now `loads::consolidate` (end of `transfer`): lines of all cases are split at each
  other's ends where they overlap on a line and added per case and segment; equal
  points and equal surface outlines of one case are added. Force and moment are
  unchanged (unit test). `plaxis_export.py` makes one object per place, whatever the
  number of cases, and sets each case's values in its phase with
  `g_i.set(obj.prop, phase, value)` (unverified against a real PLAXIS: a failure is
  reported as `unset_phase_values` and a warning, the first case's values stay).
- Cause 2: in the combination a plate whose centre of pressure is outside it (courtyards,
  openings) got one point load per element: 32 089 points on Остров. Now the points of a
  plate are merged by position (bisection, force and first moment kept) to at most
  `max_points` (default 12, UI field) per plate and sign: 830 points, 22 uniform surface
  loads, force and moment identical to the source.

- Cause 3 (second PLAXIS import): "Point is mesh-independent" — PLAXIS deletes a point
  load that is on no plate and no beam, so its load was lost. The loads are placed at
  source positions (centres of edge loads, pressure centres); the reconstructed geometry
  differs from them by up to the repair tolerance (33 mm on Остров). 53 of 234 and 301 of
  830 points were off. `loads::attach_points` (called by `add_loads` with the exported
  plate polygons and beams) moves each point to the nearest point of a polygon or a beam
  and keeps its moment about any point by a couple F × shift; counted in `approximated`.
  Line loads checked the same way: none off the structure.

### MIDAS stiffness after Lira_Midas-converter (2026-10-06)

`midas_stiffness::build` (materials, sections, thicknesses of the `.mxt`), `sections` (shapes,
catalogue `data/sortament.tsv` from the converter's `material.json`: name, kind, h, b, s, t), the
LIRA parser now reads S1, S2, S3, S5, S6, rows with only EF/EIy/EIz/GIk (S8) and the profiles of
`{13/` (`LiraParser::profiles_from`). Rules taken over: material per (E, nu, unit weight,
plate/bar) named `Plate_Beton_p12_h0.2` / `Beam_Steel_p3` after the unit weight class; nu from
the file, else steel 0.3, concrete 0.2; steel sections E = 206 GPa, 76.98 kN/m3, nu 0.3;
`DBUSER` sections SB, P, SR, H, T, B by dimensions, profiles of the block 13 from the catalogue
(DoubleT, Tubing, Pipe, Round, builtup I; designation as the fallback), `_AGT` for rigid rods,
equal shapes one section, thickness per mm, density multiplier. Two modes in the MIDAS dialog:
`converter` (nominal, the converter's file) and `lira` (EF/EIy/EIz as a `VALUE` section,
WLKE/PLKE as an equivalent thickness and E). The cap slab of a cut is always equivalent.
Where the converter is wrong and was not copied (checked by running it on its test1/test5):
S0 dimensions are taken from the S0 row (it reads EIy of the numeric row as the width), E is
not EF, a plate over 1 m thick is not read as cm (Остров: slabs 1.0, 1.2, 1.4 m), S6 and S3
with a concrete E (1e6..6e6 t/m2) stay concrete (it makes every S6 steel: round concrete
columns of «Для testa» become steel), a 3-number tube designation is looked up by all three.
Types with E = 0 (S0 rods of «Для testa», 1279 bars) are reported as without material and not
written. LIRA reductions (EI x 0.3, PLKE 0.3) are listed in `stiffness_notes` in converter mode.
PLAXIS still reads S0 and GEI only (bars of other sections: `missing_materials`).

Next (user's plan): MIDAS pipeline = every load except the
self-weight as its own load case with the source contours (`combination:
None` path of `transfer` gives that); (done: see above).

## Next coherent development batch

- Use the editor on the Багратион v4 residuals and record what the edits
  need (missing operations, picking on large models).
- Exporter (MIDAS/PLAXIS exchange format) — postponed by the user.
- Progress and cancellation for long reconstructions in the UI; frame
  cache next to the project.

### Also pending

- Rigid-body/coupled-displacement docs other than LIRA block 25 (the
  "объединение перемещений" block is not identified yet).
- The LIRA units block (33/) is not read: lengths and plate thickness are
  taken in metres.

## Explicitly not complete yet

Even where the surface audit passes, global readiness is not yet
proven. Remaining audit scope includes:

- non-parallel near misses (review items only);
- load transfer verification;
- actual MIDAS/PLAXIS import verification.

Do not set `export_ready` or equivalent final-readiness flags merely because
the current surface audit becomes green.

## Efficient execution policy

Use `AGENTS.md` for the full development workflow.

In particular:

- use small synthetic regressions as the inner loop;
- use the private full model only as an integration gate;
- perform one coherent implementation/test/review batch before asking the user
  for another "continue";
- do not add model-specific exceptions;
- update this file only after verified state or priorities change.

## GitHub Actions budget

The repository is public, so standard hosted-runner compute is currently free
for public-repository Actions. Still avoid waste because artifact/cache storage,
queue time and developer attention remain finite.

Current `.github/workflows/build.yml` runs `cargo test --locked` on pushes and
PRs to `main`/`master`; Windows/Linux release binaries are built only by a
manual dispatch with `build_release`, the Windows editor installer by a
manual dispatch with `build_app` (one job). Feature branches do not
trigger CI.
During geometry iteration:

- use `[skip ci]` when a cross-platform release build adds no information;
- do not use Actions for the private full-model loop;
- avoid unnecessary reruns;
- keep artifact-producing runs for points where a binary is actually useful.

Revisit the workflow itself only if artifact storage or unnecessary
cross-platform builds become a practical issue; do not redesign CI merely for
the sake of redesigning it.
