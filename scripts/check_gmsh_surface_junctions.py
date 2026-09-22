#!/usr/bin/env python3
"""Independently compare finite source-surface intersections with shared mesh edges.

Read-only audit of the versioned Gmsh interchange and final backend report.
Unlike the repair-delta audit, it checks every intersecting source-surface pair.
Requires the optional numpy/shapely geometry-audit dependencies.
"""

import argparse
import json
from pathlib import Path

import numpy as np

from check_v2_global_geometry import audit, covers


def repaired_segment_covered(item: dict, result: dict, mesh_edges: dict,
                             precision: float, surfaces: list[dict]) -> bool:
    """Allow only an exact, logged endpoint move, then require shared edge IDs."""
    endpoints = [np.asarray(point) for point in item["endpoints"]]
    moved = False
    for index, point in enumerate(endpoints):
        candidates = []
        for stage, key in (("junction_regularization", "accepted"), ("healing", "components")):
            for entry in result.get(stage, {}).get(key, []):
                if stage == "healing" and entry.get("status") != "collapsed":
                    continue
                for move in entry.get("moves", []):
                    before = np.asarray(move["from_coordinate"], dtype=float)
                    after = np.asarray(move["to_coordinate"], dtype=float)
                    values = (float(move["movement"]), float(move["movement_limit"]))
                    if not (np.isfinite(before).all() and np.isfinite(after).all()
                            and all(np.isfinite(value) for value in values)):
                        continue
                    displacement = np.linalg.norm(after - before)
                    if (values[1] <= 0 or values[0] <= 0 or displacement >= values[1]
                            or abs(displacement-values[0]) > precision * 5
                            or np.linalg.norm(point-before) > precision * 5):
                        continue
                    # An unrelated repair at the same coordinate does not
                    # justify this pair: the source endpoint must belong to
                    # the source boundary of one of the intersecting faces.
                    if not any(np.linalg.norm(np.asarray(vertex)-before) <= precision * 5
                               for owner in item["surfaces"]
                               for ring in surfaces[owner]["rings"] for vertex in ring):
                        continue
                    candidates.append(after)
        # Ambiguous provenance is not a licence to pick the most convenient
        # move. Only a single exact move may alter this expected endpoint.
        if len(candidates) > 1:
            return False
        if candidates:
            endpoints[index] = candidates[0]
            moved = True
    if not moved:
        return False
    start, end = endpoints
    length = np.linalg.norm(end-start)
    if length <= precision:
        return False
    direction = (end-start)/length
    left, right = item["surfaces"]
    vertices = np.asarray(result["mesh"]["vertices"])
    spans = []
    for n1, n2 in mesh_edges[left] & mesh_edges[right]:
        vector = vertices[[n1, n2]] - start
        t = vector @ direction
        if np.max(np.linalg.norm(vector-t[:, None]*direction, axis=1)) <= precision * 5:
            spans.append((float(t.min()), float(t.max())))
    return covers(spans, 0, length, precision * 5)


def check(input_data: dict, result: dict) -> dict:
    surfaces = input_data["surfaces"]
    model = {"vertices": [], "edges": [], "planes": [], "surfaces": []}
    for source, surface in enumerate(surfaces):
        rings = [np.asarray(ring, dtype=float) for ring in surface["rings"]]
        points = np.concatenate(rings)
        origin = points.mean(axis=0)
        _, _, basis = np.linalg.svd(points - origin, full_matrices=False)
        u, v, normal = basis
        plane = len(model["planes"])
        model["planes"].append(dict(origin=origin.tolist(), normal=normal.tolist(),
                                    u=u.tolist(), v=v.tolist()))
        boundaries, contours = [], []
        for ring in rings:
            start = len(model["vertices"])
            model["vertices"].extend(ring.tolist())
            boundaries.append([])
            for index in range(len(ring)):
                boundaries[-1].append({"edge": len(model["edges"])})
                model["edges"].append([start + index, start + (index + 1) % len(ring)])
            contours.append(np.column_stack(((ring-origin) @ u, (ring-origin) @ v)).tolist())
        model["surfaces"].append({"plane": plane, "contours": contours,
                                  "boundaries": boundaries,
                                  "source_elements": surface["source_elements"]})
    package = {"topology": {"preview": model, "policy": {
        "precision": float(input_data["policy"]["precision"])}},
        "mesh": {"vertices": result["mesh"]["vertices"], "triangles": [
            {"vertices": tri["vertices"], "surface": tri["source_surface"]}
            for tri in result["mesh"]["triangles"]]}}
    raw = audit(package)
    edges = {source: set() for source in range(len(surfaces))}
    for triangle in result["mesh"]["triangles"]:
        nodes = triangle["vertices"]
        owner = triangle["source_surface"]
        edges[owner].update(tuple(sorted((nodes[k], nodes[(k+1) % 3]))) for k in range(3))
    missing = [item for item in raw["contacts"] if item["mesh_conforming"] is False]
    precision = float(input_data["policy"]["precision"])
    accepted = [item for item in missing if repaired_segment_covered(item, result, edges, precision, surfaces)]
    missing = [item for item in missing if item not in accepted]
    return {"checked_segment_count": len(raw["contacts"]),
            "missing_segment_count": len(missing),
            "missing_surface_pair_count": len({tuple(x["surfaces"]) for x in missing}),
            "accepted_by_exact_repair_provenance": len(accepted),
            "missing_segments": missing, "invalid_source_surfaces": raw["invalid_surfaces"],
            "clean": not missing and not raw["invalid_surfaces"]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("result", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = check(json.loads(args.input.read_text()), json.loads(args.result.read_text()))
    if args.output:
        args.output.write_text(json.dumps(report, indent=2))
    print(f"finite surface junctions: {report['checked_segment_count']}, "
          f"unrepresented: {report['missing_segment_count']}")
    raise SystemExit(0 if report["clean"] else 1)
