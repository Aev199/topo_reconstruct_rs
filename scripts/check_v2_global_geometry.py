"""Pairwise audit of assembled v2 surfaces and their actual trial-mesh junctions.

Requires numpy and shapely>=2. Coordinates/near_distance use model length units.
This is a read-only geometric audit, not proof of solver import or load transfer.
"""
import argparse
from collections import Counter
import json
from pathlib import Path

import numpy as np
from shapely.geometry import LineString, Polygon
from shapely.validation import explain_validity


def intervals(geometry, origin, direction, lift):
    result = []
    if geometry.is_empty:
        return result
    if geometry.geom_type in ("LineString", "LinearRing"):
        xyz = lift(np.asarray(geometry.coords))
        t = (xyz - origin) @ direction
        result.append((float(t.min()), float(t.max())))
    elif hasattr(geometry, "geoms"):
        for item in geometry.geoms:
            result.extend(intervals(item, origin, direction, lift))
    return result


def covers(parts, start, end, eps):
    cursor = start
    for a, b in sorted(parts):
        if b < cursor - eps:
            continue
        if a > cursor + eps:
            return False
        cursor = max(cursor, b)
        if cursor >= end - eps:
            return True
    return cursor >= end - eps


class Surface:
    def __init__(self, index, record, model, vertices):
        self.index = index
        self.source_elements = record["source_elements"]
        plane = model["planes"][record["plane"]]
        self.o, self.n, self.u, self.v = (
            np.asarray(plane[k], dtype=float) for k in ("origin", "normal", "u", "v")
        )
        self.shape = Polygon(record["contours"][0], record["contours"][1:])
        self.rings = [self.lift(np.asarray(r)) for r in record["contours"]]
        self.points = np.concatenate(self.rings)
        self.low, self.high = self.points.min(axis=0), self.points.max(axis=0)
        self.edge_ids = {e["edge"] for ring in record["boundaries"] for e in ring}
        self.edge_ids.update(record.get("junctions", []))
        ids = {v for e in self.edge_ids for v in model["edges"][e]}
        self.planarity = float(np.max(np.abs((vertices[list(ids)] - self.o) @ self.n)))

    def project(self, points):
        d = points - self.o
        return np.column_stack((d @ self.u, d @ self.v))

    def lift(self, points):
        return self.o + points[:, :1] * self.u + points[:, 1:] * self.v


def audit(data, near_distance=0.05):
    topology = data["topology"]
    model = topology["preview"]
    eps = float(topology["policy"]["precision"]) * 5
    if not np.isfinite(near_distance) or near_distance < eps:
        raise ValueError("near distance must be finite and at least audit precision")
    vertices = np.asarray(model["vertices"], dtype=float)
    if not np.isfinite(vertices).all():
        raise ValueError("nonfinite geometry")
    surfaces = [Surface(i, s, model, vertices) for i, s in enumerate(model["surfaces"])]
    invalid = []
    for s in surfaces:
        if not s.shape.is_valid or s.shape.area <= eps * eps or s.planarity > eps:
            invalid.append(dict(surface=s.index, reason=explain_validity(s.shape),
                                planarity=s.planarity, area=s.shape.area))
    invalid_ids = {s["surface"] for s in invalid}
    mesh = data.get("mesh")
    mesh_vertices = np.asarray(mesh["vertices"], dtype=float) if mesh else None
    mesh_edges = [set() for _ in surfaces]
    if mesh:
        for t in mesh["triangles"]:
            ids = t["vertices"]
            mesh_edges[t["surface"]].update(tuple(sorted((ids[k], ids[(k+1) % 3]))) for k in range(3))
    issues, contacts, near = [], [], []
    candidates = 0
    for i, a in enumerate(surfaces):
        if i in invalid_ids:
            continue
        for b in surfaces[i+1:]:
            if b.index in invalid_ids or np.any(a.high + near_distance < b.low) or np.any(b.high + near_distance < a.low):
                continue
            candidates += 1
            pair = [i, b.index]
            direction = np.cross(a.n, b.n)
            sine = np.linalg.norm(direction)
            if sine < 1e-8:
                gap = abs(float((b.o-a.o) @ a.n))
                if gap > near_distance:
                    continue
                rings = [a.project(r) for r in b.rings]
                other = Polygon(rings[0], rings[1:])
                if not other.is_valid:
                    issues.append(dict(surfaces=pair, kind="invalid_projected_contour"))
                    continue
                common = a.shape.intersection(other)
                area_threshold = eps * max(min(a.shape.length, other.length), eps)
                if common.area > area_threshold:
                    entry = dict(surfaces=pair, area=float(common.area), gap=gap)
                    if gap <= eps:
                        issues.append(dict(entry, kind="coplanar_overlap"))
                    else:
                        near.append(dict(entry, kind="near_parallel_faces"))
                if gap > eps or common.length <= eps or common.area > area_threshold:
                    continue
                # Coplanar touching boundaries may have multiple disjoint segments.
                lines = common.geoms if hasattr(common, "geoms") else [common]
                segments = []
                for line in lines:
                    if line.geom_type != "LineString":
                        continue
                    points = a.lift(np.asarray(line.coords))
                    segments.extend(zip(points[:-1], points[1:]))
            else:
                direction /= sine
                # Compute the intersection relative to A's origin for stability.
                offset = float((b.o-a.o) @ b.n)
                origin = a.o + np.cross(direction, a.n) * (offset / sine)
                span = max(np.linalg.norm(p-origin) for p in np.concatenate((a.points, b.points))) + 1.
                line = np.array([origin-span*direction, origin+span*direction])
                ia = intervals(a.shape.intersection(LineString(a.project(line))), origin, direction, a.lift)
                ib = intervals(b.shape.intersection(LineString(b.project(line))), origin, direction, b.lift)
                segments = []
                for x0, x1 in ia:
                    for y0, y1 in ib:
                        start, end = max(x0, y0), min(x1, y1)
                        if end-start > eps:
                            segments.append((origin+start*direction, origin+end*direction))
            for start, end in segments:
                length = float(np.linalg.norm(end-start))
                if length <= eps:
                    continue
                d = (end-start)/length
                boundary_a = a.shape.boundary.buffer(eps).covers(LineString(a.project(np.array([start,end]))))
                boundary_b = b.shape.boundary.buffer(eps).covers(LineString(b.project(np.array([start,end]))))
                kind = "boundary_junction" if boundary_a and boundary_b else "t_junction" if boundary_a or boundary_b else "crossing"
                shared_geometry = a.edge_ids & b.edge_ids
                def edge_intervals(edges, points):
                    out = []
                    for edge in edges:
                        p = points[list(edge)] - start
                        t = p @ d
                        if np.max(np.linalg.norm(p-t[:,None]*d,axis=1)) <= eps:
                            out.append((float(t.min()),float(t.max())))
                    return out
                geometry_conforming = covers(edge_intervals([model["edges"][e] for e in shared_geometry],vertices),0,length,eps)
                mesh_conforming = None if mesh is None else covers(edge_intervals(mesh_edges[i] & mesh_edges[b.index],mesh_vertices),0,length,eps)
                entry = dict(surfaces=pair, kind=kind, length=length,
                             endpoints=[start.tolist(),end.tolist()],
                             geometry_conforming=geometry_conforming, mesh_conforming=mesh_conforming)
                contacts.append(entry)
                if not geometry_conforming or mesh_conforming is False:
                    issues.append(dict(entry, kind="unrepresented_intersection", contact_kind=kind))
    return dict(
        scope="surface contours, coplanar overlaps, intersection lines and shared edge identities",
        audit_complete=False,
        precision=eps, near_distance=near_distance, surfaces=len(surfaces),
        candidate_pairs=candidates, invalid_surfaces=invalid, issues=issues,
        issue_counts=dict(Counter(i["kind"] for i in issues)),
        issue_pair_count=len({tuple(i["surfaces"]) for i in issues}),
        contact_counts=dict(Counter(c["kind"] for c in contacts)),
        geometry_unrepresented_segments=sum(not c["geometry_conforming"] for c in contacts),
        mesh_unrepresented_segments=None if mesh is None else sum(c["mesh_conforming"] is False for c in contacts),
        near_face_pair_count=len(near), contacts=contacts,
        near_faces=near, maximum_planarity_error=max((s.planarity for s in surfaces),default=0.),
        global_surface_checks_passed=not invalid and not issues,
        solver_import_verified=False,
        limitations=["isolated point contacts and nonparallel near-misses are not audited",
                     "near parallel faces require interpretation; no automatic merging",
                     "bar-bar intersections and load transfer are not audited"],
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("report", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--near-distance", type=float, default=0.05)
    parser.add_argument("--strict", action="store_true",
                        help="exit 1 on detected defects (passing does not certify solver readiness)")
    args = parser.parse_args()
    result = audit(json.loads(args.report.read_text()), args.near_distance)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2))
    print(json.dumps({k:v for k,v in result.items() if k not in ("issues","contacts","near_faces")},ensure_ascii=False,indent=2))

    if args.strict and not result["global_surface_checks_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
