# Mesh quality: source review, 2026-09-14

Status: static review only; execution environment still unavailable.
This is an implementation handoff, not a completed mesh-quality fix.

## Evidence in current source

In src/reconstruction/mesh.rs, axis_nodes and edge_nodes are subdivided before
per-surface meshing. Shared boundary nodes are also incorporated into bar chains.
Per-surface constraint sets then use these shared IDs. Final topology checks
expect each constrained edge itself to exist, not an independently refined chain.

In src/reconstruction/mesh/domain.rs, local material regions call refinement
with keep_constraint_edges(). Their generated interior nodes are appended
independently. There is no feedback stage that requests new subdivisions of a
global constraint and propagates them to every owning surface and bar.

The recovered fragment report had topology_valid=true and quality_passed=false,
with minimum_angle_not_met and maximum_area_not_met, but no refinement_limit
blocker. Therefore simply increasing the resource budget is not an
evidence-supported fix.

## Interpretation, not yet a demonstrated root cause

Frozen constraint discretization is a plausible obstacle to improving triangles
near boundaries or bar contacts. The exact cause of the worst 2.088-degree and
5.911-m² triangles cannot be established without the original report/runtime.
Other candidates include a genuinely acute domain corner, close mandatory
anchors, or faces excluded by refinement. Do not label all of them as a
single confirmed defect.

## Next implementation

1. Locate failing triangles by surface, coordinates, neighboring constrained
   edges and nearby source-node IDs. Classify angle and area failures separately.
2. For failures requiring constraint subdivision, create requests against a
   canonical shared constraint and its parameter interval.
3. Apply each subdivision once globally, reusing the new vertex ID on every
   owning surface and applicable bar chain. Preserve source FE/property span
   boundaries. New geometric nodes must not introduce mechanical releases.
4. Rebuild only affected surfaces; update area and contact checks to accept
   the canonical subdivided edge chain rather than the old unsplit edge.
5. Keep strict resource/round limits and detect lack of progress. A truly acute
   input corner may prevent the requested minimum angle without a geometry
   change; report it, do not silently relax quality or delete the feature.

Do not merely remove keep_constraint_edges(): independently refined surfaces
could disagree on a common edge and disconnect the mesh. Geometry may stay
fixed while its shared discretization changes.

Acceptance cases: wall/slab common edge; coplanar bar interval; property-boundary
intersection; dangling bar; actual acute corner; close but distinct anchors;
rotation, translation and renumbering; full source/property coverage.
The same real fragment must be rerun after applying and validating the recovered
constructive-span patch. Full-building acceptance remains separate.

No runtime metrics were newly measured for this review. No build workflow
was manually dispatched.
