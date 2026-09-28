"""PLAXIS readiness profile of an assembled v2 geometry (read-only).

PLAXIS 3D intersects all imported geometry itself and fails, or meshes very
finely, where features are far below the target element size: small gaps or
overlaps between objects, tiny edges, narrow faces and needle corners (its
snap tolerance defaults to 1 mm; it warns about edges many times smaller than
the target element size). This checker measures those features relative to
a target element size h:

- short edges: surface edges and bar pieces shorter than `edge` (h/10);
- sharp corners: contour corners below `angle` degrees;
- narrow faces: a contour vertex closer than `width` (h/10) to a contour
  edge of the same surface that does not touch its neighbourhood;
- gaps: a surface vertex or bar node closer than `gap` (h/10) to another
  surface without being one of its vertices (0 excluded: that is an
  unshared contact, a failure of the global audit).

Requires numpy and shapely>=2 (same as check_v2_global_geometry.py).
"""
import argparse
from collections import Counter
import json
from pathlib import Path

import numpy as np
import shapely

from check_v2_global_geometry import Surface, surface_distances


def point_segment(p, a, b):
    d = b - a
    t = np.clip(((p - a) @ d) / max(d @ d, 1e-300), 0., 1.)
    return float(np.linalg.norm(p - (a + t * d)))


def rings_of(record, model):
    """Vertex rings of a surface from its oriented boundary edge uses."""
    out = []
    for ring in record["boundaries"]:
        ids = []
        for use in ring:
            a, b = model["edges"][use["edge"]]
            ids.append(b if use["reversed"] else a)
        out.append(ids)
    return out


def profile(data, element_size=0.5, edge=None, width=None, gap=None, angle=10.):
    edge = element_size / 10 if edge is None else edge
    width = element_size / 10 if width is None else width
    gap = element_size / 10 if gap is None else gap
    topology = data["topology"]
    model = topology["preview"]
    eps = float(topology["policy"]["precision"]) * 5
    vertices = np.asarray(model["vertices"], dtype=float)
    surfaces = [Surface(i, s, model, vertices) for i, s in enumerate(model["surfaces"])]
    items = []

    # Short surface edges (boundary and embedded), once per model edge.
    used = sorted({e for s in surfaces for e in s.edge_ids})
    for e in used:
        a, b = model["edges"][e]
        length = float(np.linalg.norm(vertices[a] - vertices[b]))
        if length < edge:
            items.append(dict(kind="short_edge", edge=e, length=length,
                              surfaces=[s.index for s in surfaces if e in s.edge_ids],
                              point=vertices[a].tolist()))
    # Short bar pieces between consecutive anchors.
    axes = topology.get("axis_assembly", {}).get("axes", [])
    bar_nodes = set()
    for i, axis in enumerate(axes):
        anchors = sorted(axis["anchors"], key=lambda a: a["t"])
        ids = [a["vertex"] for a in anchors]
        bar_nodes.update(ids)
        for u, w in zip(ids, ids[1:]):
            length = float(np.linalg.norm(vertices[u] - vertices[w]))
            if 0 < length < edge:
                items.append(dict(kind="short_bar_piece", bar=i, length=length,
                                  point=vertices[u].tolist()))

    for s, record in zip(surfaces, model["surfaces"]):
        rings = rings_of(record, model)
        # Sharp corners.
        for ring in rings:
            n = len(ring)
            for k in range(n):
                p, q, r = (vertices[ring[(k + j) % n]] for j in (-1, 0, 1))
                u, v = p - q, r - q
                cos = u @ v / max(np.linalg.norm(u) * np.linalg.norm(v), 1e-300)
                a = float(np.degrees(np.arccos(np.clip(cos, -1., 1.))))
                if a < angle:
                    items.append(dict(kind="sharp_corner", surface=s.index, vertex=ring[k],
                                      angle=a, point=q.tolist()))
        # Narrow parts: vertex to a non-neighbouring contour edge.
        segments = [(ring[k], ring[(k + 1) % len(ring)]) for ring in rings for k in range(len(ring))]
        if not segments:
            continue
        pa = np.array([vertices[a] for a, _ in segments])
        pb = np.array([vertices[b] for _, b in segments])
        neighbours = {}
        for a, b in segments:
            neighbours.setdefault(a, set()).update((a, b))
            neighbours.setdefault(b, set()).update((a, b))
        for ring in rings:
            for v in ring:
                near = set().union(*(neighbours[w] for w in neighbours[v]))
                p = vertices[v]
                d = pb - pa
                t = np.clip(np.einsum("ij,ij->i", p - pa, d) /
                            np.maximum(np.einsum("ij,ij->i", d, d), 1e-300), 0., 1.)
                dist = np.linalg.norm(p - (pa + t[:, None] * d), axis=1)
                for k in np.nonzero(dist < width)[0]:
                    a, b = segments[k]
                    if a in near or b in near:
                        continue
                    items.append(dict(kind="narrow_face", surface=s.index, vertex=v,
                                      edge=[a, b], width=float(dist[k]), point=p.tolist()))

    # Gaps between a surface vertex / bar node and another surface.
    vertex_sets = [{v for e in s.edge_ids for v in model["edges"][e]} for s in surfaces]
    owners = {}
    for s in surfaces:
        for v in vertex_sets[s.index]:
            owners.setdefault(v, []).append(s)
    candidates = sorted(set(owners) | bar_nodes)
    points = vertices[candidates]
    for b in surfaces:
        low, high = b.low - gap, b.high + gap
        inside = np.all((points >= low) & (points <= high), axis=1)
        ids = np.array(candidates)[inside]
        if not len(ids):
            continue
        dist, _ = surface_distances(b, vertices[ids])
        for v, d in zip(ids, dist):
            v = int(v)
            if eps < d < gap and v not in vertex_sets[b.index]:
                items.append(dict(kind="gap", vertex=v, surface=b.index, distance=float(d),
                                  of=[s.index for s in owners.get(v, [])],
                                  bar_node=v in bar_nodes, point=vertices[v].tolist()))

    counts = Counter(i["kind"] for i in items)
    return dict(
        scope="PLAXIS readiness: short edges, sharp corners, narrow faces, gaps",
        element_size=element_size,
        thresholds=dict(edge=edge, width=width, gap=gap, angle_degrees=angle),
        counts=dict(counts),
        passed=not items,
        items=items,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("report", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--element-size", type=float, default=0.5)
    parser.add_argument("--edge", type=float)
    parser.add_argument("--width", type=float)
    parser.add_argument("--gap", type=float)
    parser.add_argument("--angle", type=float, default=10.)
    parser.add_argument("--strict", action="store_true", help="exit 1 unless every check passes")
    args = parser.parse_args()
    result = profile(json.loads(args.report.read_text()), args.element_size,
                     args.edge, args.width, args.gap, args.angle)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2))
    print(json.dumps({k: v for k, v in result.items() if k != "items"}, ensure_ascii=False, indent=2))
    if args.strict and not result["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
