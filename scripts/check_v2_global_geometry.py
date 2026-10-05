"""Pairwise audit of assembled v2 surfaces, bars and their trial-mesh junctions.

Requires numpy and shapely>=2. Coordinates/near_distance use model length units.
Checked: surface validity, coplanar overlaps, surface-surface intersection lines,
isolated surface point contacts, bar-bar and bar-surface intersections (in the
topology and as shared nodes/edges of the trial mesh), and property transfer to
the trial mesh, where missing data is a failure. Near misses are review items, not failures.
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


def audit_point_contacts(surfaces, invalid_ids, vertices, model, eps, near_distance,
                         mesh_index=None):
    """Surface vertices touching another surface without a shared vertex.

    With a mesh, a shared vertex must also be one mesh node of both surfaces.
    """
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
                    if a.index < b.index and mesh_index is not None and not any(
                            m in mesh_index.surface_nodes[a.index]
                            and m in mesh_index.surface_nodes[b.index]
                            for m in mesh_index.near(vertices[v], eps)):
                        issues.append(dict(vertex=v, surfaces=[a.index, b.index], distance=float(d),
                                           point=vertices[v].tolist(),
                                           kind="shared_vertex_not_shared_in_mesh"))
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


class MeshIndex:
    """Mesh node identity: which surfaces/axes use a node, and a spatial hash."""

    def __init__(self, mesh, surface_count, axis_count, cell):
        self.vertices = np.asarray(mesh["vertices"], dtype=float)
        self.cell = cell
        self.surface_edges = [set() for _ in range(surface_count)]
        self.surface_nodes = [set() for _ in range(surface_count)]
        self.axis_edges = [set() for _ in range(axis_count)]
        self.axis_nodes = [set() for _ in range(axis_count)]
        for t in mesh["triangles"]:
            ids = t["vertices"]
            if 0 <= t["surface"] < surface_count:
                self.surface_nodes[t["surface"]].update(ids)
                self.surface_edges[t["surface"]].update(
                    tuple(sorted((ids[k], ids[(k + 1) % 3]))) for k in range(3))
        for b in mesh.get("bars", []):
            if 0 <= b.get("axis", -1) < axis_count:
                self.axis_nodes[b["axis"]].update(b["vertices"])
                self.axis_edges[b["axis"]].add(tuple(sorted(b["vertices"])))
        self.grid = {}
        used = set().union(*self.surface_nodes, *self.axis_nodes)
        for m in used:
            self.grid.setdefault(self.key(self.vertices[m]), []).append(m)

    def key(self, p):
        return tuple(np.floor(p / self.cell).astype(np.int64))

    def near(self, x, radius):
        cx = self.key(x)
        found = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for m in self.grid.get((cx[0] + dx, cx[1] + dy, cx[2] + dz), ()):
                        if np.linalg.norm(self.vertices[m] - x) <= radius:
                            found.append(m)
        return found


def segment_intervals(pairs, points, start, direction, eps):
    """Parameters along start+t*direction of the segments lying on that line."""
    out = []
    for pair in pairs:
        p = points[list(pair)] - start
        t = p @ direction
        if np.max(np.linalg.norm(p - t[:, None] * direction, axis=1)) <= eps:
            out.append((float(t.min()), float(t.max())))
    return out


def audit_bars(data, surfaces, invalid_ids, vertices, eps, near_distance, mesh_index=None):
    """Bar connectivity in the topology and, when a mesh exists, in the mesh.

    Topology: every crossing needs a shared vertex; a bar lying in a surface
    needs its whole in-surface length covered by surface edges or interval
    contacts. Mesh: the same crossings need a shared mesh node, and the
    in-surface length must be covered by bar segments that are also triangle
    edges of that surface. Contact records alone never prove a mesh junction.
    """
    bars = data["topology"].get("axis_assembly", {})
    axes = bars.get("axes", [])
    contacts = bars.get("contacts", [])
    model = data["topology"]["preview"]
    segments = []
    for i, axis in enumerate(axes):
        p0, p1 = (vertices[v] for v in axis["endpoints"])
        anchors = {a["vertex"] for a in axis["anchors"]} | set(axis["endpoints"])
        segments.append((i, p0, p1, anchors, np.minimum(p0, p1), np.maximum(p0, p1)))
    issues, review = [], []

    def mesh_shared(x, *groups):
        """A mesh node at x used by every group (sets of mesh node ids)."""
        return any(all(m in g for g in groups) for m in mesh_index.near(x, eps))

    if mesh_index is not None:
        for i, p0, p1, _, _, _ in segments:
            length = float(np.linalg.norm(p1 - p0))
            if length <= eps:
                continue
            u = (p1 - p0) / length
            parts = segment_intervals(mesh_index.axis_edges[i], mesh_index.vertices, p0, u, eps)
            if not covers(parts, 0., length, eps):
                issues.append(dict(bar=i, kind="bar_not_covered_by_mesh", length=length,
                                   point=p0.tolist()))
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
                elif mesh_index is not None and not mesh_shared(
                        x, mesh_index.axis_nodes[i], mesh_index.axis_nodes[j]):
                    issues.append(dict(entry, kind="bar_intersection_not_shared_in_mesh"))
            elif j not in neighbours[i] and not neighbours[i] & neighbours[j]:
                review.append(dict(entry, kind="bar_near_miss"))
    surface_vertices = [{v for e in s.edge_ids for v in model["edges"][e]} for s in surfaces]
    point_contacts, interval_contacts = {}, {}
    for c in contacts:
        if c["kind"] == "point":
            point_contacts.setdefault((c["axis"], c["surface"]), []).append(c["vertex"])
        elif c["kind"] == "interval":
            interval_contacts.setdefault((c["axis"], c["surface"]), []).append(
                (float(c["start_t"]), float(c["end_t"])))

    def check_point(i, s, x, anchors, entry, radius=None):
        """A bar meeting a surface at one point x, located within radius."""
        radius = eps if radius is None else max(radius, eps)
        records = point_contacts.get((i, s.index), [])
        shared = [v for v in set(records) | anchors
                  if np.linalg.norm(vertices[v] - x) <= radius
                  and (v in surface_vertices[s.index] or v in records)]
        if not shared:
            issues.append(dict(entry, kind="unshared_bar_surface_intersection", point=x.tolist()))
        elif mesh_index is not None and not any(
                m in mesh_index.axis_nodes[i] and m in mesh_index.surface_nodes[s.index]
                for m in mesh_index.near(vertices[shared[0]], eps)):
            issues.append(dict(entry, kind="bar_surface_point_not_shared_in_mesh", point=x.tolist()))

    for i, p0, p1, anchors, lo, hi in segments:
        length = float(np.linalg.norm(p1 - p0))
        for s in surfaces:
            if s.index in invalid_ids:
                continue
            if np.any(hi + near_distance < s.low) or np.any(s.high + near_distance < lo):
                continue
            h0, h1 = ((p - s.o) @ s.n for p in (p0, p1))
            entry = dict(bar=i, surface=s.index)
            if abs(h0) <= eps and abs(h1) <= eps and length > eps:
                u = (p1 - p0) / length
                inside = s.shape.buffer(eps).intersection(LineString(s.project(np.array([p0, p1]))))
                parts = [(max(a, 0.), min(b, length))
                         for a, b in intervals(inside, p0, u, s.lift)]
                lines = [(a, b) for a, b in parts if b - a > 10 * eps]
                # The buffered region extends each part by about eps where the
                # bar leaves the surface obliquely.
                slack = 4 * eps
                for a, b in lines:
                    item = dict(entry, start=a, end=b, length=b - a,
                                point=(p0 + (a + b) / 2 * u).tolist())
                    edges = [model["edges"][e] for e in s.edge_ids]
                    topo = segment_intervals(edges, vertices, p0, u, eps) + [
                        (t0 * length, t1 * length)
                        for t0, t1 in interval_contacts.get((i, s.index), [])]
                    if not covers(topo, a + slack, b - slack, slack):
                        issues.append(dict(item, kind="bar_in_surface_without_contact"))
                    elif mesh_index is not None:
                        shared = mesh_index.axis_edges[i] & mesh_index.surface_edges[s.index]
                        mesh_parts = segment_intervals(shared, mesh_index.vertices, p0, u, eps)
                        if not covers(mesh_parts, a + slack, b - slack, slack):
                            issues.append(dict(item, kind="bar_in_surface_not_shared_in_mesh"))
                if lines:
                    continue
                # Touching at points only (typically an end): each touch point
                # must be connected.
                # A touch found in the eps-buffered region is up to its own
                # length away from the exact contact point.
                for a, b in parts:
                    check_point(i, s, p0 + (a + b) / 2 * u, anchors, entry, (b - a) / 2 + eps)
                if not parts and (i, s.index) not in point_contacts:
                    d, _ = surface_distances(s, np.array([p0, p1]))
                    if np.min(d) <= near_distance:
                        review.append(dict(entry, kind="bar_surface_near_miss",
                                           distance=float(np.min(d))))
                continue
            if h0 * h1 < 0 or abs(h0) <= eps or abs(h1) <= eps:
                t = h0 / (h0 - h1) if abs(h0 - h1) > 0 else 0.
                x = p0 + np.clip(t, 0., 1.) * (p1 - p0)
                d, _ = surface_distances(s, x[None, :])
                if d[0] <= eps:
                    check_point(i, s, x, anchors, entry)
                    continue
            if (i, s.index) in point_contacts or (i, s.index) in interval_contacts:
                continue
            d, _ = surface_distances(s, np.array([p0, p1]))
            if np.min(d) <= near_distance:
                review.append(dict(entry, kind="bar_surface_near_miss", distance=float(np.min(d))))
    return issues, review


def audit_properties(data, surface_count):
    """Property transfer to the trial mesh; missing data counts as a failure."""
    topology, mesh = data["topology"], data.get("mesh")
    result = dict(loads="not available: the reconstruction input carries geometry and stiffness only")
    if not mesh:
        return result
    stiffness = topology.get("surface_stiffness")
    axes = topology.get("axis_assembly", {}).get("axes", [])
    missing_triangles = wrong_triangles = 0
    meshed = set()
    for t in mesh["triangles"]:
        s = t.get("surface")
        if isinstance(s, int) and 0 <= s < surface_count:
            meshed.add(s)
        if "stiffness" not in t or stiffness is None or not isinstance(s, int) \
                or not 0 <= s < len(stiffness) or stiffness[s] is None:
            missing_triangles += 1
        elif t["stiffness"] != stiffness[s]:
            wrong_triangles += 1
    missing_bars = wrong_bars = 0
    for b in mesh.get("bars", []):
        a = b.get("axis")
        if not isinstance(a, int) or not 0 <= a < len(axes) or not axes[a].get("spans") \
                or "stiffness" not in b or "source_element" not in b:
            missing_bars += 1
        elif (b["source_element"], b["stiffness"]) not in \
                {(sp["element"], sp["stiffness"]) for sp in axes[a]["spans"]}:
            wrong_bars += 1
    result.update(triangles_with_wrong_stiffness=wrong_triangles,
                  triangles_with_missing_stiffness=missing_triangles,
                  bars_with_wrong_stiffness=wrong_bars,
                  bars_with_missing_stiffness=missing_bars,
                  surfaces_without_triangles=sorted(set(range(surface_count)) - meshed))
    return result


# Coarsest numerical precision the audit accepts (model units, 1 um): the
# audited report cannot loosen the audit by declaring a coarser precision.
MAXIMUM_PRECISION = 1e-6


def audit(data, near_distance=0.05, maximum_precision=MAXIMUM_PRECISION):
    topology = data["topology"]
    model = topology["preview"]
    eps = min(float(topology["policy"]["precision"]), maximum_precision) * 5
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
    axis_count = len(topology.get("axis_assembly", {}).get("axes", []))
    mesh_index = MeshIndex(mesh, len(surfaces), axis_count, max(eps, 1e-12)) if mesh else None
    mesh_vertices = mesh_index.vertices if mesh else None
    mesh_edges = mesh_index.surface_edges if mesh else [set() for _ in surfaces]
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
                if common.area > 0:
                    # GEOS overlay can return area outside both inputs along
                    # nearly coincident collinear edges (a 6 m2 phantom
                    # between two walls touching along a line): the overlap
                    # is recomputed on a grid a thousandth of the tolerance.
                    grid = eps / 1000.
                    common = shapely.set_precision(a.shape, grid).intersection(
                        shapely.set_precision(other, grid))
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
                    return segment_intervals(edges, points, start, d, eps)
                geometry_conforming = covers(edge_intervals([model["edges"][e] for e in shared_geometry],vertices),0,length,eps)
                mesh_conforming = None if mesh is None else covers(edge_intervals(mesh_edges[i] & mesh_edges[b.index],mesh_vertices),0,length,eps)
                entry = dict(surfaces=pair, kind=kind, length=length,
                             endpoints=[start.tolist(),end.tolist()],
                             geometry_conforming=geometry_conforming, mesh_conforming=mesh_conforming)
                contacts.append(entry)
                if not geometry_conforming or mesh_conforming is False:
                    issues.append(dict(entry, kind="unrepresented_intersection", contact_kind=kind))
    point_issues, point_review = audit_point_contacts(
        surfaces, invalid_ids, vertices, model, eps, near_distance, mesh_index)
    bar_issues, bar_review = audit_bars(data, surfaces, invalid_ids, vertices, eps,
                                        near_distance, mesh_index)
    properties = audit_properties(data, len(surfaces))
    extra_issues = point_issues + bar_issues
    property_ok = not any(properties.get(k) for k in (
        "triangles_with_wrong_stiffness", "triangles_with_missing_stiffness",
        "bars_with_wrong_stiffness", "bars_with_missing_stiffness", "surfaces_without_triangles"))
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
    parser.add_argument("--maximum-precision", type=float, default=MAXIMUM_PRECISION,
                        help="coarsest precision accepted from the report (model units)")
    parser.add_argument("--strict", action="store_true",
                        help="exit 1 on detected defects (passing does not certify solver readiness)")
    args = parser.parse_args()
    result = audit(json.loads(args.report.read_text()), args.near_distance,
                   args.maximum_precision)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2))
    hidden = ("issues", "contacts", "near_faces", "point_and_bar_issues", "review_items")
    print(json.dumps({k:v for k,v in result.items() if k not in hidden},ensure_ascii=False,indent=2))

    if args.strict and not result["global_checks_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
