#!/usr/bin/env python3
"""Gmsh/OpenCASCADE backend prototype for topo_reconstruct_rs.

The preferred input is the versioned `topo-reconstruct-gmsh-v1` interchange
emitted by Rust. Legacy v2 preview JSON is accepted only to keep the private
full-model regression reproducible during migration.

Engineering recognition stays in Rust. This backend:
1. creates OpenCASCADE planar faces from explicit 3D rings;
2. runs exact General Fuse (occ.fragment) with no global fuzzy tolerance;
3. carries source/property ownership through the Boolean output map;
4. builds one conformal 2D mesh;
5. performs conservative sub-resolution micro-edge cleanup;
6. emits a versioned mesh/provenance result for later MIDAS/PLAXIS adapters.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import gmsh

INPUT_FORMAT = "topo-reconstruct-gmsh-v1"
RESULT_FORMAT = "topo-reconstruct-gmsh-result-v1"


def distance(a: Iterable[float], b: Iterable[float]) -> float:
    return math.dist(tuple(a), tuple(b))


def triangle_area(a: Iterable[float], b: Iterable[float], c: Iterable[float]) -> float:
    a = tuple(a)
    b = tuple(b)
    c = tuple(c)
    ab = tuple(b[i] - a[i] for i in range(3))
    ac = tuple(c[i] - a[i] for i in range(3))
    cross = (
        ab[1] * ac[2] - ab[2] * ac[1],
        ab[2] * ac[0] - ab[0] * ac[2],
        ab[0] * ac[1] - ab[1] * ac[0],
    )
    return 0.5 * math.sqrt(sum(x * x for x in cross))


def triangle_angles(a: Iterable[float], b: Iterable[float], c: Iterable[float]) -> list[float]:
    a = tuple(a)
    b = tuple(b)
    c = tuple(c)
    ab = distance(a, b)
    ac = distance(a, c)
    bc = distance(b, c)
    if min(ab, ac, bc) <= 0:
        return [0.0, 0.0, 180.0]
    result = []
    for opposite, left, right in ((bc, ab, ac), (ac, ab, bc), (ab, ac, bc)):
        cosine = (left * left + right * right - opposite * opposite) / (2 * left * right)
        result.append(math.degrees(math.acos(max(-1.0, min(1.0, cosine)))))
    return result


def lift(plane: dict, uv: Iterable[float]) -> list[float]:
    u, v = uv
    return [
        plane["origin"][k] + u * plane["u"][k] + v * plane["v"][k]
        for k in range(3)
    ]


def normalize_input(data: dict, mesh_size_override: float | None) -> dict:
    if data.get("format") == INPUT_FORMAT:
        normalized = data
    else:
        # Migration path from the previous full v2 preview report.
        topology = data["topology"]
        model = topology["preview"]
        surfaces = []
        stiffness = topology.get("surface_stiffness", [])
        patches = topology.get("surface_source_patches", [])
        for index, surface in enumerate(model["surfaces"]):
            plane = model["planes"][surface["plane"]]
            surfaces.append(
                {
                    "source_surface": index,
                    "source_patch": patches[index] if index < len(patches) else index,
                    "stiffness": stiffness[index] if index < len(stiffness) else 0,
                    "source_elements": surface.get("source_elements", []),
                    "rings": [
                        [lift(plane, uv) for uv in ring]
                        for ring in surface["contours"]
                    ],
                }
            )
        policy = topology.get("policy", {})
        normalized = {
            "format": INPUT_FORMAT,
            "length_unit": "model_unit",
            "policy": {
                "precision": float(policy.get("precision", 1e-7)),
                "minimum_edge": float(policy.get("minimum_edge", 1e-3)),
                "target_mesh_size": 0.75,
            },
            "source_coverage_complete": bool(
                topology.get("all_surface_patches_built", True)
                and topology.get("axis_assembly", {}).get("all_axes_built", True)
            ),
            "surfaces": surfaces,
            "axes": [],
            "contacts": [],
            "blockers": [],
        }

    policy = normalized["policy"]
    if mesh_size_override is not None:
        policy = dict(policy)
        policy["target_mesh_size"] = mesh_size_override
        normalized = dict(normalized)
        normalized["policy"] = policy

    precision = float(policy["precision"])
    minimum_edge = float(policy["minimum_edge"])
    target = float(policy["target_mesh_size"])
    if (
        not math.isfinite(precision)
        or not math.isfinite(minimum_edge)
        or not math.isfinite(target)
        or precision <= 0
        or minimum_edge <= precision
        or target < minimum_edge
    ):
        raise ValueError("invalid backend policy")
    return normalized


def add_ring(points3d: list[list[float]]) -> int:
    if len(points3d) < 3:
        raise ValueError("surface ring has fewer than 3 points")
    point_tags = [gmsh.model.occ.addPoint(*p) for p in points3d]
    curves = [
        gmsh.model.occ.addLine(point_tags[i], point_tags[(i + 1) % len(point_tags)])
        for i in range(len(point_tags))
    ]
    return gmsh.model.occ.addWire(curves)


def add_surface(surface: dict) -> int:
    wires = [add_ring(ring) for ring in surface["rings"]]
    return gmsh.model.occ.addPlaneSurface(wires)


def element_triangles(surface_tag: int) -> list[tuple[int, int, int]]:
    triangles = []
    types, _element_tags, node_tags = gmsh.model.mesh.getElements(2, surface_tag)
    for element_type, nodes in zip(types, node_tags):
        name, _dim, _order, nodes_per_element, _local, primary_nodes = (
            gmsh.model.mesh.getElementProperties(element_type)
        )
        if not name.lower().startswith("triangle") or int(primary_nodes) < 3:
            raise RuntimeError(
                f"unsupported 2D element on surface {surface_tag}: {name}"
            )
        nodes_per_element = int(nodes_per_element)
        for i in range(0, len(nodes), nodes_per_element):
            triangles.append(tuple(int(x) for x in nodes[i : i + 3]))
    return triangles


def all_node_coordinates() -> dict[int, tuple[float, float, float]]:
    tags, xyz, _ = gmsh.model.mesh.getNodes()
    return {
        int(tag): (float(xyz[i]), float(xyz[i + 1]), float(xyz[i + 2]))
        for i, tag in zip(range(0, len(xyz), 3), tags)
    }


def quality(coords: dict[int, tuple[float, float, float]], triangles: list[dict]) -> dict:
    minimum_angle = 180.0
    maximum_area = 0.0
    below_20 = 0
    below_5 = 0
    for triangle in triangles:
        a, b, c = (coords[n] for n in triangle["nodes"])
        area = triangle_area(a, b, c)
        angles = triangle_angles(a, b, c)
        local_min = min(angles)
        minimum_angle = min(minimum_angle, local_min)
        maximum_area = max(maximum_area, area)
        below_20 += local_min < 20.0
        below_5 += local_min < 5.0
    return {
        "triangle_count": len(triangles),
        "minimum_triangle_angle_degrees": minimum_angle if triangles else None,
        "maximum_triangle_area": maximum_area if triangles else None,
        "triangles_below_20_degrees": below_20,
        "triangles_below_5_degrees": below_5,
    }


class Dsu:
    def __init__(self) -> None:
        self.parent: dict[int, int] = {}

    def find(self, x: int) -> int:
        parent = self.parent.setdefault(x, x)
        if parent != x:
            self.parent[x] = self.find(parent)
        return self.parent[x]

    def union(self, a: int, b: int) -> None:
        a = self.find(a)
        b = self.find(b)
        if a != b:
            self.parent[max(a, b)] = min(a, b)


def mesh_edges(triangles: list[dict]) -> set[tuple[int, int]]:
    result = set()
    for triangle in triangles:
        a, b, c = triangle["nodes"]
        result.update(
            {
                tuple(sorted((a, b))),
                tuple(sorted((b, c))),
                tuple(sorted((c, a))),
            }
        )
    return result


def node_source_owners(triangles: list[dict]) -> dict[int, set[int]]:
    owners: dict[int, set[int]] = defaultdict(set)
    for triangle in triangles:
        for node in triangle["nodes"]:
            owners[node].add(triangle["source_surface"])
    return owners


def heal_micro_edges(
    coords: dict[int, tuple[float, float, float]],
    triangles: list[dict],
    minimum_edge: float,
    precision: float,
) -> tuple[list[dict], dict]:
    """Collapse only proven junction-local sub-resolution edge components.

    A component is eligible only when it contains an already shared node owned
    by at least two reconstructed source surfaces. That node remains fixed.
    This deliberately refuses general proximity welding.
    """

    edges = mesh_edges(triangles)
    short_edges = [
        edge
        for edge in edges
        if precision < distance(coords[edge[0]], coords[edge[1]]) < minimum_edge
    ]
    dsu = Dsu()
    for a, b in short_edges:
        dsu.union(a, b)

    components: dict[int, set[int]] = defaultdict(set)
    for a, b in short_edges:
        components[dsu.find(a)].update((a, b))

    owners = node_source_owners(triangles)
    mapping: dict[int, int] = {}
    diagnostics = []
    maximum_movement = 0.0
    unresolved = 0

    for nodes in sorted(components.values(), key=lambda x: (min(x), len(x))):
        candidates = [
            node for node in nodes if len(owners.get(node, set())) >= 2
        ]
        if not candidates:
            diagnostics.append(
                {
                    "nodes": sorted(nodes),
                    "status": "unresolved_no_shared_junction_node",
                }
            )
            unresolved += 1
            continue
        maximum_owner_count = max(len(owners[node]) for node in candidates)
        representatives = sorted(
            node for node in candidates if len(owners[node]) == maximum_owner_count
        )
        # Two distinct already-shared junction nodes are structural evidence,
        # not numerical noise. Never choose between them by an arbitrary tag.
        if len(representatives) != 1:
            diagnostics.append(
                {
                    "nodes": sorted(nodes),
                    "candidate_representatives": representatives,
                    "status": "unresolved_multiple_shared_junction_nodes",
                }
            )
            unresolved += 1
            continue

        representative = representatives[0]
        movements = {node: distance(coords[node], coords[representative]) for node in nodes}
        component_max = max(movements.values(), default=0.0)
        # A transitive chain of individually short edges can span much farther
        # than the engineering minimum-edge threshold. Refuse such a collapse.
        if component_max >= minimum_edge:
            diagnostics.append(
                {
                    "nodes": sorted(nodes),
                    "representative": representative,
                    "status": "unresolved_component_extent",
                    "maximum_movement": component_max,
                }
            )
            unresolved += 1
            continue

        for node in nodes:
            if node != representative:
                mapping[node] = representative
        maximum_movement = max(maximum_movement, component_max)
        diagnostics.append(
            {
                "nodes": sorted(nodes),
                "representative": representative,
                "representative_source_owner_count": len(owners[representative]),
                "status": "collapsed",
                "maximum_movement": component_max,
            }
        )

    healed = []
    removed_degenerate = 0
    near_degenerate = 0
    for triangle in triangles:
        nodes = tuple(mapping.get(node, node) for node in triangle["nodes"])
        if len(set(nodes)) < 3:
            removed_degenerate += 1
            continue
        if triangle_area(*(coords[node] for node in nodes)) <= precision * precision:
            near_degenerate += 1
            # Do not silently discard a non-topologically-degenerate triangle.
            healed.append({**triangle, "nodes": nodes})
            continue
        healed.append({**triangle, "nodes": nodes})

    return healed, {
        "minimum_edge": minimum_edge,
        "short_edge_count_before": len(short_edges),
        "component_count": len(components),
        "collapsed_component_count": sum(x["status"] == "collapsed" for x in diagnostics),
        "unresolved_component_count": unresolved,
        "merged_node_count": len(mapping),
        "maximum_node_movement": maximum_movement,
        "removed_degenerate_triangles": removed_degenerate,
        "near_degenerate_triangles_after": near_degenerate,
        "components": diagnostics,
    }


def shared_mesh_edges_by_source(triangles: list[dict]) -> dict[tuple[int, int], int]:
    edge_owners: dict[tuple[int, int], set[int]] = defaultdict(set)
    for triangle in triangles:
        a, b, c = triangle["nodes"]
        for edge in (
            tuple(sorted((a, b))),
            tuple(sorted((b, c))),
            tuple(sorted((c, a))),
        ):
            edge_owners[edge].add(triangle["source_surface"])

    result: dict[tuple[int, int], int] = defaultdict(int)
    for owners in edge_owners.values():
        owners = sorted(owners)
        for i in range(len(owners)):
            for j in range(i + 1, len(owners)):
                result[(owners[i], owners[j])] += 1
    return dict(result)


@dataclass
class Fragmentation:
    source_to_output: list[list[int]]
    output_owners: dict[int, list[int]]
    output_surfaces: list[int]


def fragment(input_surfaces: list[dict]) -> Fragmentation:
    input_entities = [(2, add_surface(surface)) for surface in input_surfaces]
    if not input_entities:
        return Fragmentation([], {}, [])

    if len(input_entities) == 1:
        source_to_output = [[input_entities[0][1]]]
        output_surfaces = [input_entities[0][1]]
    else:
        output, output_map = gmsh.model.occ.fragment(
            [input_entities[0]],
            input_entities[1:],
            removeObject=True,
            removeTool=True,
        )
        source_to_output = [
            sorted({tag for dim, tag in mapped if dim == 2}) for mapped in output_map
        ]
        output_surfaces = sorted({tag for dim, tag in output if dim == 2})

    gmsh.model.occ.synchronize()
    output_owners: dict[int, list[int]] = defaultdict(list)
    for source, tags in enumerate(source_to_output):
        for tag in tags:
            output_owners[tag].append(source)
    return Fragmentation(
        source_to_output=source_to_output,
        output_owners={tag: sorted(set(owners)) for tag, owners in output_owners.items()},
        output_surfaces=output_surfaces,
    )


def physical_groups(fragmentation: Fragmentation, surfaces: list[dict]) -> list[dict]:
    by_stiffness: dict[int, list[int]] = defaultdict(list)
    for tag in fragmentation.output_surfaces:
        owners = fragmentation.output_owners.get(tag, [])
        if len(owners) != 1:
            continue
        by_stiffness[int(surfaces[owners[0]]["stiffness"])].append(tag)

    result = []
    for stiffness, tags in sorted(by_stiffness.items()):
        group = gmsh.model.addPhysicalGroup(2, sorted(set(tags)))
        name = f"stiffness_{stiffness}"
        gmsh.model.setPhysicalName(2, group, name)
        result.append(
            {
                "dimension": 2,
                "physical_tag": group,
                "name": name,
                "stiffness": stiffness,
                "surface_tags": sorted(set(tags)),
            }
        )
    return result


def run_backend(
    data: dict,
    mesh_size_override: float | None = None,
    write_msh: Path | None = None,
) -> dict:
    data = normalize_input(data, mesh_size_override)
    surfaces = data["surfaces"]
    policy = data["policy"]
    precision = float(policy["precision"])
    minimum_edge = float(policy["minimum_edge"])
    mesh_size = float(policy["target_mesh_size"])

    gmsh.clear()
    gmsh.model.add("topo_reconstruct_occ")
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("Geometry.OCCBooleanPreserveNumbering", 1)
    gmsh.option.setNumber("Mesh.ElementOrder", 1)
    # Exact General Fuse: do not use Geometry.ToleranceBoolean as a repair
    # budget. The full-model trial showed that even 0.1 mm changed topology too
    # aggressively.
    gmsh.option.setNumber("Geometry.ToleranceBoolean", 0.0)

    fragmentation = fragment(surfaces)

    ownership_conflicts = []
    fragmented_surfaces = []
    unmapped_output_surfaces = []
    for tag in fragmentation.output_surfaces:
        owners = fragmentation.output_owners.get(tag, [])
        if not owners:
            unmapped_output_surfaces.append(tag)
            continue
        if len(owners) > 1:
            stiffnesses = sorted({int(surfaces[owner]["stiffness"]) for owner in owners})
            ownership_conflicts.append(
                {
                    "output_surface": tag,
                    "source_surfaces": owners,
                    "stiffnesses": stiffnesses,
                }
            )
            continue
        source = owners[0]
        item = surfaces[source]
        fragmented_surfaces.append(
            {
                "output_surface": tag,
                "source_surface": source,
                "source_patch": item["source_patch"],
                "stiffness": item["stiffness"],
                "source_elements": item["source_elements"],
            }
        )

    groups = physical_groups(fragmentation, surfaces)

    gmsh.option.setNumber("Mesh.MeshSizeMin", mesh_size)
    gmsh.option.setNumber("Mesh.MeshSizeMax", mesh_size)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.Smoothing", 10)

    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.optimize("Relocate2D", force=True, niter=10)

    coords = all_node_coordinates()
    triangles = []
    non_triangle_blockers = []
    for item in fragmented_surfaces:
        tag = item["output_surface"]
        try:
            local = element_triangles(tag)
        except RuntimeError as error:
            non_triangle_blockers.append(str(error))
            continue
        for nodes in local:
            triangles.append(
                {
                    "nodes": nodes,
                    "output_surface": tag,
                    "source_surface": item["source_surface"],
                    "stiffness": item["stiffness"],
                }
            )

    raw_quality = quality(coords, triangles)
    raw_shared = shared_mesh_edges_by_source(triangles)
    healed_triangles, healing = heal_micro_edges(
        coords, triangles, minimum_edge, precision
    )
    healed_quality = quality(coords, healed_triangles)
    healed_shared = shared_mesh_edges_by_source(healed_triangles)

    used_nodes = sorted({node for triangle in healed_triangles for node in triangle["nodes"]})
    compact_index = {tag: index for index, tag in enumerate(used_nodes)}
    compact_vertices = [coords[tag] for tag in used_nodes]
    compact_triangles = [
        {
            "vertices": [compact_index[node] for node in triangle["nodes"]],
            "output_surface": triangle["output_surface"],
            "source_surface": triangle["source_surface"],
            "stiffness": triangle["stiffness"],
        }
        for triangle in healed_triangles
    ]

    blockers = list(data.get("blockers", []))
    if ownership_conflicts:
        blockers.append("ambiguous_fragment_ownership")
    if unmapped_output_surfaces:
        blockers.append("unmapped_fragment_surface")
    if non_triangle_blockers:
        blockers.append("unsupported_2d_elements")
    if healing["unresolved_component_count"]:
        blockers.append("unresolved_micro_edge_component")
    if healing["near_degenerate_triangles_after"]:
        blockers.append("near_degenerate_triangle_after_healing")

    if write_msh:
        gmsh.write(str(write_msh))

    return {
        "format": RESULT_FORMAT,
        "backend": "gmsh_opencascade_fragment",
        "gmsh_version": gmsh.option.getString("General.Version"),
        "length_unit": data.get("length_unit", "model_unit"),
        "policy": policy,
        "source_coverage_complete": bool(data.get("source_coverage_complete", False)),
        "backend_ready": not blockers,
        "blockers": sorted(set(blockers)),
        "input_surface_count": len(surfaces),
        "output_surface_count": len(fragmentation.output_surfaces),
        "source_to_output_surface_count": [
            len(tags) for tags in fragmentation.source_to_output
        ],
        "fragmented_surfaces": fragmented_surfaces,
        "ownership_conflicts": ownership_conflicts,
        "unmapped_output_surfaces": unmapped_output_surfaces,
        "physical_groups": groups,
        "raw_mesh": {
            **raw_quality,
            "shared_surface_pair_count": len(raw_shared),
            "shared_mesh_edge_count": sum(raw_shared.values()),
        },
        "healing": healing,
        "mesh": {
            **healed_quality,
            "shared_surface_pair_count": len(healed_shared),
            "shared_mesh_edge_count": sum(healed_shared.values()),
            "vertices": compact_vertices,
            "triangles": compact_triangles,
        },
        # Carry reconstructed line semantics forward unchanged. The next backend
        # stage will embed/synchronize them against the fragmented surfaces.
        "axes": data.get("axes", []),
        "contacts": data.get("contacts", []),
    }


def synthetic_report() -> dict:
    surfaces = [
        {
            "source_surface": 0,
            "source_patch": 0,
            "stiffness": 10,
            "source_elements": [1],
            "rings": [[
                [-2.0, -1.0, 0.0],
                [2.0, -1.0, 0.0],
                [2.0, 1.0, 0.0],
                [-2.0, 1.0, 0.0],
            ]],
        },
        {
            "source_surface": 1,
            "source_patch": 1,
            "stiffness": 20,
            "source_elements": [2],
            "rings": [[
                [0.0, -2.0, -1.0],
                [0.0, 2.0, -1.0],
                [0.0, 2.0, 1.0],
                [0.0, -2.0, 1.0],
            ]],
        },
    ]
    return {
        "format": INPUT_FORMAT,
        "length_unit": "model_unit",
        "policy": {
            "precision": 1e-8,
            "minimum_edge": 1e-3,
            "target_mesh_size": 0.5,
        },
        "source_coverage_complete": True,
        "surfaces": surfaces,
        "axes": [],
        "contacts": [],
        "blockers": [],
    }


def self_test() -> dict:
    result = run_backend(synthetic_report())
    assert result["input_surface_count"] == 2
    assert result["output_surface_count"] >= 3
    assert result["mesh"]["triangle_count"] > 0
    assert result["mesh"]["shared_surface_pair_count"] >= 1
    assert result["mesh"]["shared_mesh_edge_count"] >= 1
    assert not result["ownership_conflicts"]
    assert result["backend_ready"], result["blockers"]

    # Pure healing regression: only the pre-existing multi-surface junction node
    # is allowed to absorb the two adjacent parasitic nodes.
    coords = {
        1: (0.0, 0.0, 0.0),
        2: (0.00004, 0.0, 0.0),
        3: (0.0, 0.00008, 0.0),
        4: (1.0, 0.0, 0.0),
        5: (0.0, 1.0, 0.0),
        6: (0.0, 0.0, 1.0),
        7: (0.0, 1.0, 1.0),
    }
    triangles = [
        {"nodes": (1, 2, 4), "source_surface": 0, "output_surface": 1, "stiffness": 10},
        {"nodes": (1, 4, 5), "source_surface": 0, "output_surface": 1, "stiffness": 10},
        {"nodes": (1, 3, 6), "source_surface": 1, "output_surface": 2, "stiffness": 20},
        {"nodes": (1, 6, 7), "source_surface": 1, "output_surface": 2, "stiffness": 20},
    ]
    healed, report = heal_micro_edges(coords, triangles, 0.001, 1e-8)
    assert report["collapsed_component_count"] == 1
    assert report["unresolved_component_count"] == 0
    assert report["maximum_node_movement"] < 0.001
    assert all(1 in triangle["nodes"] for triangle in healed)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", nargs="?", type=Path)
    parser.add_argument("--output", type=Path, help="summary/full backend JSON")
    parser.add_argument("--msh", type=Path)
    parser.add_argument("--mesh-size", type=float)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()

    gmsh.initialize()
    try:
        if args.self_test:
            result = self_test()
        else:
            if args.input is None:
                parser.error("input is required unless --self-test is used")
            result = run_backend(
                json.loads(args.input.read_text()),
                mesh_size_override=args.mesh_size,
                write_msh=args.msh,
            )
        text = json.dumps(result, indent=2)
        if args.output:
            args.output.write_text(text)
        print(text)
        return 0
    finally:
        gmsh.finalize()


if __name__ == "__main__":
    raise SystemExit(main())
