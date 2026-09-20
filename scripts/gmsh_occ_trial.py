#!/usr/bin/env python3
"""Trial OpenCASCADE/Gmsh backend for reconstructed v2 surfaces.

This is deliberately an external experiment: Rust remains responsible for
engineering recognition/reconstruction; Gmsh/OpenCASCADE is asked to fragment
the reconstructed planar surfaces into conformal CAD topology and mesh it.

Input is a --v2-preview-json report produced by topo_reconstruct_rs.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Iterable

import gmsh


def lift(plane: dict, uv: Iterable[float]) -> list[float]:
    u, v = uv
    return [
        plane["origin"][k] + u * plane["u"][k] + v * plane["v"][k]
        for k in range(3)
    ]


def add_ring(points3d: list[list[float]]) -> int:
    point_tags = [gmsh.model.occ.addPoint(*p) for p in points3d]
    curves = [
        gmsh.model.occ.addLine(point_tags[i], point_tags[(i + 1) % len(point_tags)])
        for i in range(len(point_tags))
    ]
    return gmsh.model.occ.addWire(curves)


def add_surface(surface: dict, plane: dict) -> int:
    wires = []
    for ring in surface["contours"]:
        if len(ring) < 3:
            raise ValueError("surface ring has fewer than 3 points")
        wires.append(add_ring([lift(plane, uv) for uv in ring]))
    return gmsh.model.occ.addPlaneSurface(wires)


def triangle_edges_for_surface(tag: int) -> tuple[set[tuple[int, int]], int, int]:
    """Return global first-order triangle edges, triangle count, other 2D count."""
    result: set[tuple[int, int]] = set()
    triangle_count = 0
    other_count = 0
    types, _element_tags, node_tags = gmsh.model.mesh.getElements(2, tag)
    for element_type, nodes in zip(types, node_tags):
        props = gmsh.model.mesh.getElementProperties(element_type)
        name = props[0]
        nodes_per_element = int(props[3])
        primary_nodes = int(props[5])
        if name.lower().startswith("triangle") and primary_nodes >= 3:
            for i in range(0, len(nodes), nodes_per_element):
                tri = [int(x) for x in nodes[i : i + 3]]
                triangle_count += 1
                result.update(
                    {
                        tuple(sorted((tri[0], tri[1]))),
                        tuple(sorted((tri[1], tri[2]))),
                        tuple(sorted((tri[2], tri[0]))),
                    }
                )
        else:
            other_count += len(nodes) // max(nodes_per_element, 1)
    return result, triangle_count, other_count


def triangle_quality() -> dict:
    node_tags, xyz, _ = gmsh.model.mesh.getNodes()
    coords = {
        int(tag): (xyz[i], xyz[i + 1], xyz[i + 2])
        for i, tag in zip(range(0, len(xyz), 3), node_tags)
    }
    minimum_angle = 180.0
    maximum_area = 0.0
    triangles = 0
    for _dim, surface_tag in gmsh.model.getEntities(2):
        types, _tags, node_lists = gmsh.model.mesh.getElements(2, surface_tag)
        for element_type, nodes in zip(types, node_lists):
            props = gmsh.model.mesh.getElementProperties(element_type)
            if not props[0].lower().startswith("triangle"):
                continue
            nodes_per_element = int(props[3])
            for i in range(0, len(nodes), nodes_per_element):
                a, b, c = (coords[int(n)] for n in nodes[i : i + 3])
                ab = tuple(b[k] - a[k] for k in range(3))
                ac = tuple(c[k] - a[k] for k in range(3))
                bc = tuple(c[k] - b[k] for k in range(3))
                lab = math.dist(a, b)
                lac = math.dist(a, c)
                lbc = math.dist(b, c)
                cross = (
                    ab[1] * ac[2] - ab[2] * ac[1],
                    ab[2] * ac[0] - ab[0] * ac[2],
                    ab[0] * ac[1] - ab[1] * ac[0],
                )
                area = 0.5 * math.sqrt(sum(x * x for x in cross))
                maximum_area = max(maximum_area, area)
                if min(lab, lac, lbc) > 0:
                    angles = []
                    for opposite, left, right in (
                        (lbc, lab, lac),
                        (lac, lab, lbc),
                        (lab, lac, lbc),
                    ):
                        cosine = (left * left + right * right - opposite * opposite) / (
                            2 * left * right
                        )
                        angles.append(math.degrees(math.acos(max(-1.0, min(1.0, cosine)))))
                    minimum_angle = min(minimum_angle, *angles)
                triangles += 1
    return {
        "triangles": triangles,
        "minimum_triangle_angle_degrees": minimum_angle if triangles else None,
        "maximum_triangle_area": maximum_area if triangles else None,
    }


def fragment_and_mesh(data: dict, mesh_size: float | None = None, write_msh: Path | None = None) -> dict:
    topology = data["topology"]
    model = topology["preview"]
    surfaces = model["surfaces"]
    planes = model["planes"]

    gmsh.clear()
    gmsh.model.add("topo_reconstruct_occ_trial")
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("Geometry.OCCBooleanPreserveNumbering", 1)
    gmsh.option.setNumber("Mesh.ElementOrder", 1)

    input_entities: list[tuple[int, int]] = []
    for surface in surfaces:
        tag = add_surface(surface, planes[surface["plane"]])
        input_entities.append((2, tag))

    if len(input_entities) > 1:
        # General Fuse all surfaces in one OCC operation. The output map is
        # ordered as objects followed by tools, so source lineage is retained.
        objects = [input_entities[0]]
        tools = input_entities[1:]
        output, output_map = gmsh.model.occ.fragment(
            objects, tools, removeObject=True, removeTool=True
        )
        source_to_output = [
            sorted({tag for dim, tag in mapped if dim == 2}) for mapped in output_map
        ]
        output_surfaces = sorted({tag for dim, tag in output if dim == 2})
    else:
        source_to_output = [[input_entities[0][1]]] if input_entities else []
        output_surfaces = [input_entities[0][1]] if input_entities else []

    gmsh.model.occ.synchronize()

    source_curves: list[set[int]] = []
    for mapped in source_to_output:
        curves: set[int] = set()
        for surface_tag in mapped:
            curves.update(
                tag
                for dim, tag in gmsh.model.getBoundary(
                    [(2, surface_tag)], combined=False, oriented=False, recursive=False
                )
                if dim == 1
            )
        source_curves.append(curves)

    shared_curve_pairs = []
    for i in range(len(source_curves)):
        for j in range(i + 1, len(source_curves)):
            shared = sorted(source_curves[i] & source_curves[j])
            if shared:
                shared_curve_pairs.append(
                    {"surfaces": [i, j], "curve_count": len(shared), "curves": shared}
                )

    if mesh_size is not None:
        if not math.isfinite(mesh_size) or mesh_size <= 0:
            raise ValueError("mesh size must be finite and positive")
        gmsh.option.setNumber("Mesh.MeshSizeMin", mesh_size)
        gmsh.option.setNumber("Mesh.MeshSizeMax", mesh_size)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

    gmsh.model.mesh.generate(2)

    source_mesh_edges: list[set[tuple[int, int]]] = []
    per_source_triangles = []
    non_triangle_elements = 0
    for mapped in source_to_output:
        edges: set[tuple[int, int]] = set()
        triangles = 0
        for surface_tag in mapped:
            e, count, other = triangle_edges_for_surface(surface_tag)
            edges.update(e)
            triangles += count
            non_triangle_elements += other
        source_mesh_edges.append(edges)
        per_source_triangles.append(triangles)

    conforming_mesh_pairs = []
    for pair in shared_curve_pairs:
        i, j = pair["surfaces"]
        shared_edges = source_mesh_edges[i] & source_mesh_edges[j]
        conforming_mesh_pairs.append(
            {
                "surfaces": [i, j],
                "shared_curve_count": pair["curve_count"],
                "shared_mesh_edge_count": len(shared_edges),
                "mesh_conforming": bool(shared_edges),
            }
        )

    all_nodes, _xyz, _ = gmsh.model.mesh.getNodes()
    quality = triangle_quality()
    mapped_output_owners: dict[int, list[int]] = {}
    for source, mapped in enumerate(source_to_output):
        for output_tag in mapped:
            mapped_output_owners.setdefault(output_tag, []).append(source)

    if write_msh:
        gmsh.write(str(write_msh))

    return {
        "backend": "gmsh_opencascade_fragment",
        "gmsh_version": gmsh.option.getString("General.Version"),
        "input_surfaces": len(surfaces),
        "output_surfaces": len(output_surfaces),
        "source_to_output_surface_count": [len(x) for x in source_to_output],
        "output_surfaces_with_multiple_source_owners": sum(
            len(owners) > 1 for owners in mapped_output_owners.values()
        ),
        "shared_curve_pairs": len(shared_curve_pairs),
        "shared_curves": sum(x["curve_count"] for x in shared_curve_pairs),
        "mesh_conforming_shared_curve_pairs": sum(
            x["mesh_conforming"] for x in conforming_mesh_pairs
        ),
        "mesh_nonconforming_shared_curve_pairs": sum(
            not x["mesh_conforming"] for x in conforming_mesh_pairs
        ),
        "mesh_nodes": len(all_nodes),
        "non_triangle_2d_elements": non_triangle_elements,
        "per_source_triangle_count": per_source_triangles,
        **quality,
        "pair_details": conforming_mesh_pairs,
    }


def synthetic_report() -> dict:
    # Two surfaces use independent vertices and cross in the interior. This is
    # the exact class that currently forces custom shared-junction logic.
    planes = [
        {
            "origin": [0.0, 0.0, 0.0],
            "normal": [0.0, 0.0, 1.0],
            "u": [1.0, 0.0, 0.0],
            "v": [0.0, 1.0, 0.0],
        },
        {
            "origin": [0.0, 0.0, 0.0],
            "normal": [1.0, 0.0, 0.0],
            "u": [0.0, 1.0, 0.0],
            "v": [0.0, 0.0, 1.0],
        },
    ]
    surfaces = [
        {
            "plane": 0,
            "contours": [[[-2.0, -1.0], [2.0, -1.0], [2.0, 1.0], [-2.0, 1.0]]],
            "source_elements": [1],
        },
        {
            "plane": 1,
            "contours": [[[-2.0, -1.0], [2.0, -1.0], [2.0, 1.0], [-2.0, 1.0]]],
            "source_elements": [2],
        },
    ]
    return {"topology": {"preview": {"planes": planes, "surfaces": surfaces}}}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", nargs="?", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--msh", type=Path)
    parser.add_argument("--mesh-size", type=float, default=0.5)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()

    gmsh.initialize()
    try:
        data = synthetic_report() if args.self_test else json.loads(args.input.read_text())
        result = fragment_and_mesh(data, args.mesh_size, args.msh)
        if args.self_test:
            assert result["input_surfaces"] == 2
            assert result["shared_curve_pairs"] >= 1
            assert result["mesh_nonconforming_shared_curve_pairs"] == 0
            assert result["triangles"] > 0
        text = json.dumps(result, indent=2)
        if args.output:
            args.output.write_text(text)
        print(text)
        return 0
    finally:
        gmsh.finalize()


if __name__ == "__main__":
    raise SystemExit(main())
