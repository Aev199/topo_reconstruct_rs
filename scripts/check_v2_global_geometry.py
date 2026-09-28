"""Pairwise audit of assembled v2 surfaces, bars and their trial-mesh junctions.

Requires numpy and shapely>=2. Coordinates/near_distance use model length units.
Checked: surface validity, coplanar overlaps, surface-surface intersection lines,
isolated surface point contacts, bar-bar and bar-surface intersections, and
property transfer to the trial mesh. Near misses are review items, not failures.
This is a read-only geometric audit, not proof of solver import or load transfer.
"""
import argparse
from collections import Counter
import json
from pathlib import Path

import numpy as np
import shapely
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
        boundary = {e["edge"] for ring in record["boundaries"] for e in ring}
        # Embedded edges are junction lines inside the material. They count as
        # surface topology only if they are checked to lie in the surface.
        self.embedded = set(record.get("embedded_edges", []))
        self.edge_ids = boundary | self.embedded
        ids = {v for e in self.edge_ids for v in model["edges"][e]}
        self.planarity = float(np.max(np.abs((vertices[list(ids)] - self.o) @ self.n)))
        self.embedded_on_boundary = sorted(self.embedded & boundary)
        self.embedded_segments = [(e, LineString(self.project(vertices[list(model["edges"][e])])))
                                  for e in sorted(self.embedded - boundary)]

    def project(self, points):
        d = points - self.o
        return np.column_stack((d @ self.u, d @ self.v))

    def lift(self, points):
        return self.o + points[:, :1] * self.u + points[:, 1:] * self.v


def closest_points(p0, p1, q0, q1):
    """Parameters and distance of the closest points of two segments."""
    d1, d2, r = p1 - p0, q1 - q0, p0 - q0
    a, e, f = d1 @ d1, d2 @ d2, d2 @ r
    c, b = d1 @ r, d1 @ d2
    denom = a * e - b * b
    s = np.clip((b * f - c * e) / denom, 0., 1.) if denom > 1e-18 * a * e else 0.
    t = (b * s + f) / e
    if t < 0.:
        t, s = 0., np.clip(-c / a, 0., 1.)
    elif t > 1.:
        t, s = 1., np.clip((b - c) / a, 0., 1.)
    x, y = p0 + s * d1, q0 + t * d2
    return s, t, float(np.linalg.norm(x - y)), x


def surface_distances(surface, points):
    """Distance of 3D points to a surface polygon, and their plane heights."""
    h = (points - surface.o) @ surface.n
    xy = surface.project(points)
    d2 = shapely.distance(surface.shape, shapely.points(xy))
    return np.hypot(h, d2), h


def audit_point_contacts(surfaces, invalid_ids, vertices, model, eps, near_distance):
    """Surface vertices touching another surface without a shared vertex."""
    vertex_sets = [{v for e in s.edge_ids for v in model["edges"][e]} for s in surfaces]
    issues, review = [], []
    for a in surfaces:
        if a.index in invalid_ids:
            continue
        ids = np.array(sorted(vertex_sets[a.index]))
        points = vertices[ids]
        for b in surfaces:
            if b.index == a.index or b.index in invalid_ids:
                continue
            if np.any(a.high + near_distance < b.low) or np.any(b.high + near_distance < a.low):
                continue
            dist, _ = surface_distances(b, points)
            for v, d in zip(ids[dist <= near_distance], dist[dist <= near_distance]):
                v = int(v)
                if v in vertex_sets[b.index]:
                    continue
                entry = dict(vertex=v, surfaces=[a.index, b.index], distance=float(d),
                             point=vertices[v].tolist())
                if d <= eps:
                    # On b but not one of its vertices: unshared unless it lies
                    # on one of b's edges (a split would be required) - both
                    # cases mean the contact has no shared identity.
                    issues.append(dict(entry, kind="unshared_point_contact"))
                else:
                    review.append(dict(entry, kind="surface_near_miss"))
    return issues, review


def audit_bars(data, surfaces, invalid_ids, vertices, eps, near_distance):
    bars = data["topology"].get("axis_assembly", {})
    axes = bars.get("axes", [])
    contacts = bars.get("contacts", [])
    segments = []
    for i, axis in enumerate(axes):
        p0, p1 = (vertices[v] for v in axis["endpoints"])
        anchors = {a["vertex"] for a in axis["anchors"]} | set(axis["endpoints"])
        segments.append((i, p0, p1, anchors, np.minimum(p0, p1), np.maximum(p0, p1)))
    issues, review = [], []
    # Bars sharing a vertex are connected; a near miss between bars joined
    # through one intermediate bar is a short bar, not a gap.
    by_vertex = {}
    for i, _, _, anchors, _, _ in segments:
        for v in anchors:
            by_vertex.setdefault(v, set()).add(i)
    neighbours = [set().union(*(by_vertex[v] for v in anchors)) - {i}
                  for i, _, _, anchors, _, _ in segments]
    for i, p0, p1, _, _, _ in segments:
        length = float(np.linalg.norm(p1 - p0))
        if length < near_distance:
            review.append(dict(bars=[i], kind="short_bar", length=length, point=p0.tolist()))
    for k, (i, p0, p1, anchors_i, lo_i, hi_i) in enumerate(segments):
        for j, q0, q1, anchors_j, lo_j, hi_j in segments[k + 1:]:
            if np.any(hi_i + near_distance < lo_j) or np.any(hi_j + near_distance < lo_i):
                continue
            _, _, dist, x = closest_points(p0, p1, q0, q1)
            if dist > near_distance:
                continue
            d1, d2 = p1 - p0, q1 - q0
            parallel = np.linalg.norm(np.cross(d1, d2)) <= 1e-9 * np.linalg.norm(d1) * np.linalg.norm(d2)
            entry = dict(bars=[i, j], distance=dist, point=x.tolist())
            if dist <= eps:
                if parallel:
                    u = d1 / np.linalg.norm(d1)
                    ts = sorted([0., float(d1 @ u)])
                    tq = sorted([float((q0 - p0) @ u), float((q1 - p0) @ u)])
                    overlap = min(ts[1], tq[1]) - max(ts[0], tq[0])
                    if overlap > eps:
                        issues.append(dict(entry, kind="overlapping_bars", length=overlap))
                        continue
                shared = [v for v in anchors_i & anchors_j
                          if np.linalg.norm(vertices[v] - x) <= eps]
                if not shared:
                    issues.append(dict(entry, kind="unshared_bar_intersection"))
            elif j not in neighbours[i] and not neighbours[i] & neighbours[j]:
                review.append(dict(entry, kind="bar_near_miss"))
    surface_vertices = [{v for e in s.edge_ids for v in data["topology"]["preview"]["edges"][e]}
                        for s in surfaces]
    contact_pairs = {(c["axis"], c["surface"]) for c in contacts}
    point_contacts = {}
    for c in contacts:
        if c["kind"] == "point":
            point_contacts.setdefault((c["axis"], c["surface"]), []).append(c["vertex"])
    for i, p0, p1, anchors, lo, hi in segments:
        for s in surfaces:
            if s.index in invalid_ids:
                continue
            if np.any(hi + near_distance < s.low) or np.any(s.high + near_distance < lo):
                continue
            h0, h1 = ((p - s.o) @ s.n for p in (p0, p1))
            entry = dict(bar=i, surface=s.index)
            if abs(h0) <= eps and abs(h1) <= eps:
                inside = s.shape.buffer(eps).intersection(LineString(s.project(np.array([p0, p1]))))
                # A touch at one point is checked like a piercing point below.
                if inside.length > 10 * eps and (i, s.index) not in contact_pairs:
                    issues.append(dict(entry, kind="bar_in_surface_without_contact",
                                       length=float(inside.length)))
                if inside.length > 10 * eps:
                    continue
                # Touching at an end: that end must be connected.
                d, _ = surface_distances(s, np.array([p0, p1]))
                ends = data["topology"]["axis_assembly"]["axes"][i]["endpoints"]
                for v, dist in zip(ends, d):
                    if dist <= eps and v not in surface_vertices[s.index] \
                            and v not in point_contacts.get((i, s.index), []):
                        issues.append(dict(entry, kind="unshared_bar_surface_intersection",
                                           point=vertices[v].tolist()))
                if np.min(d) > eps and np.min(d) <= near_distance \
                        and (i, s.index) not in contact_pairs:
                    review.append(dict(entry, kind="bar_surface_near_miss", distance=float(np.min(d))))
                continue
            if h0 * h1 < 0 or abs(h0) <= eps or abs(h1) <= eps:
                t = h0 / (h0 - h1) if abs(h0 - h1) > 0 else 0.
                x = p0 + np.clip(t, 0., 1.) * (p1 - p0)
                d, _ = surface_distances(s, x[None, :])
                if d[0] <= eps:
                    # Connected through a point contact or a vertex of the surface.
                    shared = [v for v in set(point_contacts.get((i, s.index), [])) | anchors
                              if np.linalg.norm(vertices[v] - x) <= eps
                              and (v in surface_vertices[s.index]
                                   or v in point_contacts.get((i, s.index), []))]
                    if not shared:
                        issues.append(dict(entry, kind="unshared_bar_surface_intersection",
                                           point=x.tolist()))
                    continue
            if (i, s.index) in contact_pairs:
                continue
            d, _ = surface_distances(s, np.array([p0, p1]))
            if np.min(d) <= near_distance:
                review.append(dict(entry, kind="bar_surface_near_miss", distance=float(np.min(d))))
    return issues, review


def audit_properties(data):
    topology, mesh = data["topology"], data.get("mesh")
    result = dict(loads="not available: the reconstruction input carries geometry and stiffness only")
    if not mesh:
        return result
    stiffness = topology.get("surface_stiffness", [])
    wrong_triangles = sum(1 for t in mesh["triangles"]
                          if "stiffness" in t and t["surface"] < len(stiffness)
                          and t["stiffness"] != stiffness[t["surface"]])
    axes = topology.get("axis_assembly", {}).get("axes", [])
    wrong_bars = sum(1 for b in mesh.get("bars", [])
                     if (b["source_element"], b["stiffness"]) not in
                     {(sp["element"], sp["stiffness"]) for sp in axes[b["axis"]]["spans"]})
    result.update(triangles_with_wrong_stiffness=wrong_triangles,
                  bars_with_wrong_stiffness=wrong_bars)
    return result


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
            continue
        region = s.shape.buffer(eps)
        outside = [e for e, line in s.embedded_segments
                   if line.length <= eps or not region.covers(line)]
        if outside or s.embedded_on_boundary:
            invalid.append(dict(surface=s.index, reason="invalid embedded edges",
                                embedded_outside=outside,
                                embedded_on_boundary=s.embedded_on_boundary,
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
    point_issues, point_review = audit_point_contacts(
        surfaces, invalid_ids, vertices, model, eps, near_distance)
    bar_issues, bar_review = audit_bars(data, surfaces, invalid_ids, vertices, eps, near_distance)
    properties = audit_properties(data)
    extra_issues = point_issues + bar_issues
    property_ok = not properties.get("triangles_with_wrong_stiffness") and \
        not properties.get("bars_with_wrong_stiffness")
    surface_passed = not invalid and not issues
    return dict(
        scope="surfaces, intersection lines, point contacts, bars, property transfer",
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
        global_surface_checks_passed=surface_passed,
        point_and_bar_issues=extra_issues,
        point_and_bar_issue_counts=dict(Counter(i["kind"] for i in extra_issues)),
        review_items=point_review + bar_review,
        review_counts=dict(Counter(i["kind"] for i in point_review + bar_review)),
        properties=properties,
        global_checks_passed=surface_passed and not extra_issues and property_ok,
        solver_import_verified=False,
        limitations=["near misses and near parallel faces are review items; no automatic merging",
                     "loads are not part of the reconstruction input and are not audited",
                     "a passing audit does not prove a successful MIDAS/PLAXIS import"],
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
    hidden = ("issues", "contacts", "near_faces", "point_and_bar_issues", "review_items")
    print(json.dumps({k:v for k,v in result.items() if k not in hidden},ensure_ascii=False,indent=2))

    if args.strict and not result["global_checks_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
