# Geometry regression matrix

Updated: 2026-09-20

This matrix is the current private integration gate for the Gmsh/OpenCASCADE
backend. Only aggregate results are committed; source fixtures remain private.

All models use the same reconstruction and backend policy:

- precision: 1e-7 model units;
- minimum edge: 0.001 model units;
- junction movement limit: 0.05 model units;
- target Gmsh surface size: 0.75 model units;
- exact OCC General Fuse; no global Boolean fuzzy tolerance.

## Current result

| metric | skala1 | skala2 | test5 |
| --- | ---: | ---: | ---: |
| backend ready | yes | yes | yes |
| reconstructed surfaces | 387 | 23 | 108 |
| reconstructed axes | 564 | 291 | 8 |
| shell triangles | 22,220 | 16,872 | 37,183 |
| bar segments | 3,189 | 2,522 | 56 |
| min triangle angle | 0.80084° | 21.95616° | 2.07156° |
| triangles < 20° | 34 | 0 | 413 |
| triangles < 5° | 8 | 0 | 4 |
| declared bar/surface contacts | 4,354 | 2,787 | 12 |
| conforming bar/surface contacts | 4,354 | 2,787 | 12 |
| final unintended shared nodes | 0 | 0 | 0 |
| duplicate shell/bar elements | 0 / 0 | 0 / 0 | 0 / 0 |
| finite surface-pair changes after repair | 0 | 0 | 0 |
| maximum surface planarity error | 1.65e-8 | 1.78e-13 | 2.55e-9 |

Coverage is complete on all three fixtures.

### skala1 coverage

- surface regions: 387 / 387;
- bar regions: 2,601 / 2,601;
- source shell elements: 13,301 / 13,301;
- source bar elements: 2,601 / 2,601;
- no invalid/degenerate compact elements;
- no stiffness mismatches.

Near-vertex regularization is inactive on this model. Existing sub-millimetre
micro-edge healing remains the only geometry movement:

- 2 parasitic nodes merged;
- maximum movement about 0.08365 mm;
- 3 exactly degenerate shell triangles removed;
- 2 point contacts retained through exact logged movement provenance.

This is an explicit non-regression check for the original full-model case.

### skala2 coverage

- surface regions: 23 / 23;
- bar regions: 2,522 / 2,522;
- source shell elements: 9,964 / 9,964;
- source bar elements: 2,522 / 2,522;
- no repair or healing was required;
- no invalid/degenerate elements;
- no stiffness mismatches.

This is the clean control fixture: the new repair logic performs zero geometry
operations when the exact OCC/Gmsh result is already good.

### test5 coverage

- surface regions: 108 / 108;
- bar regions: 8 / 8;
- source shell elements: 36,866 / 36,866;
- source bar elements: 8 / 8;
- no invalid/degenerate compact elements;
- no stiffness mismatches.

This fixture exposed two geometry classes not present in the original full
model.

## Near-vertex junction regularization

The raw exact-OCC mesh of test5 contained 30 triangles below 5 degrees, with
the worst minimum angle about 0.0004 degrees. The source surface contours
themselves were not degenerate; the slivers were created by surface
intersections passing very close to existing structural vertices.

The accepted rule is intentionally narrower than proximity welding:

1. only the shortest edge of a triangle below 5 degrees is considered;
2. the edge must be shorter than the explicit 50 mm engineering junction
   movement budget;
3. an original source boundary edge is never collapsed;
4. a weak source vertex may move to an already exact multi-surface OCC junction;
5. two weak source vertices already lying on the same exact shared OCC edge may
   move to the midpoint of that edge;
6. every movement is logged with from/to node, coordinates and movement limit;
7. surface planarity, finite surface-pair identity, duplicates, coverage and
   semantic topology are re-audited afterwards.

On test5:

- 14 candidate edges;
- 11 accepted repairs;
- 15 node mappings;
- maximum movement about 30.10 mm;
- 22 triangles made exactly degenerate by the collapse and removed;
- triangles below 5 degrees reduced from 30 to 4;
- minimum angle improved to 2.07156 degrees;
- no finite surface pair was lost or created;
- maximum final planarity error remained about 2.55e-9 model units.

The four remaining sub-5-degree triangles are not treated as a topology
failure. Two are associated with a short generated shared-junction segment and
two with generated interior mesh edges. Chasing these further would turn the
geometry backend into a general mesh optimizer; leave them for the established
mesher/target-solver quality gate unless a real import rejects them.

## Isolated point-contact semantics

Exact CAD geometry can make two surfaces share a mesh node even when the
reconstructed engineering semantics do not justify a connection.

The backend now audits the complete owner graph of every shared mesh node. If
one geometric node contains multiple disconnected semantic components, the
components receive duplicate node IDs at the exact same coordinates.

This is a topology-only operation:

- coordinates do not move;
- finite shared surface edges are unchanged;
- declared rod/surface contacts are preserved;
- only unsupported point-only identity is separated.

On test5:

- 16 unsupported point-only shared nodes detected;
- 16 duplicate IDs created;
- final unintended shared nodes: 0.

On skala1 and skala2 this operation is inactive.

## Required regression gates

A future geometry change is not acceptable merely because the mesh looks
better. All private fixtures must continue to satisfy:

- backend_ready=true;
- complete surface and bar region coverage;
- complete source-element coverage;
- zero invalid or degenerate compact elements;
- zero duplicate compact shell/bar elements;
- zero stiffness mismatches;
- zero failed declared bar/surface contacts;
- zero unintended final shared nodes;
- zero lost or newly created finite shared-surface pairs due to repair;
- surface planarity within numerical tolerance;
- every geometry movement bounded and provenance logged.

Mesh-angle statistics are diagnostics, not a substitute for these topology
gates.

## Current priority

Continue geometry validation on diverse real FE fixtures. Do not resume FPN
serialization work until geometry regression coverage is considered broad
enough.
