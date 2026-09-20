#!/usr/bin/env python3
"""Gmsh/OpenCASCADE backend prototype for topo_reconstruct_rs.

The JSON result keeps scalar quality counters separate from mesh arrays.
Mixed-dimensional regression covers both embedded rod intervals and point crossings.

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
SOLVER_MESH_FORMAT = "topo-reconstruct-solver-mesh-v1"


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
                    "ring_source_nodes": [
                        [
                            int(topology["vertex_source_nodes"][
                                (
                                    model["edges"][edge_use["edge"]][1]
                                    if edge_use["reversed"]
                                    else model["edges"][edge_use["edge"]][0]
                                )
                            ])
                            for edge_use in boundary
                        ]
                        for boundary in surface["boundaries"]
                    ],
                }
            )
        axis_report = topology.get("axis_assembly", {})
        axes = []
        for axis in axis_report.get("axes", []):
            endpoints = [model["vertices"][int(vertex)] for vertex in axis["endpoints"]]
            axes.append(
                {
                    "source_axis": int(axis["source_axis"]),
                    "endpoints": endpoints,
                    "property_spans": [
                        {
                            "source_element": int(span["element"]),
                            "stiffness": int(span["stiffness"]),
                            "start_t": float(span["start_t"]),
                            "end_t": float(span["end_t"]),
                        }
                        for span in axis.get("spans", [])
                    ],
                    "anchors": [
                        {
                            "source_node": int(anchor["source_node"]),
                            "t": float(anchor["t"]),
                            "point": model["vertices"][int(anchor["vertex"])],
                        }
                        for anchor in axis.get("anchors", [])
                    ],
                }
            )

        contacts = []
        for contact in axis_report.get("contacts", []):
            axis_index = int(contact["axis"])
            item = {
                "kind": contact["kind"],
                "axis": axis_index,
                "surface": int(contact["surface"]),
                "location": contact["location"],
            }
            if contact["kind"] == "point":
                t = float(contact["t"])
                item.update(
                    {
                        "t": t,
                        "point": list(axis_point(axes[axis_index], t)),
                    }
                )
            else:
                start_t = float(contact["start_t"])
                end_t = float(contact["end_t"])
                item.update(
                    {
                        "start_t": start_t,
                        "end_t": end_t,
                        "endpoints": [
                            list(axis_point(axes[axis_index], start_t)),
                            list(axis_point(axes[axis_index], end_t)),
                        ],
                    }
                )
            contacts.append(item)

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
                and axis_report.get("all_axes_built", True)
            ),
            "surfaces": surfaces,
            "axes": axes,
            "contacts": contacts,
            "blockers": [],
        }

    policy = normalized["policy"]
    if mesh_size_override is not None:
        policy = dict(policy)
        policy["target_mesh_size"] = mesh_size_override
        normalized = dict(normalized)
        normalized["policy"] = policy

    for surface in normalized.get("surfaces", []):
        rings = surface.get("rings", [])
        source_rings = surface.get("ring_source_nodes")
        if source_rings is None or len(rings) != len(source_rings):
            raise ValueError("missing surface source-node provenance")
        if any(len(ring) != len(nodes) for ring, nodes in zip(rings, source_rings)):
            raise ValueError("surface source-node provenance length mismatch")

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


def axis_point(axis: dict, t: float) -> tuple[float, float, float]:
    a = axis["endpoints"][0]
    b = axis["endpoints"][1]
    return tuple(float(a[i]) + t * (float(b[i]) - float(a[i])) for i in range(3))


def axis_length(axis: dict) -> float:
    return distance(axis["endpoints"][0], axis["endpoints"][1])


def unique_parameters(values: Iterable[float], length: float, precision: float) -> list[float]:
    ordered = sorted(float(value) for value in values)
    result = []
    for value in ordered:
        if not math.isfinite(value) or value < -precision / max(length, precision) or value > 1.0 + precision / max(length, precision):
            raise ValueError("axis split parameter outside endpoints")
        value = min(1.0, max(0.0, value))
        if not result or abs(value - result[-1]) * length > precision:
            result.append(value)
    return result


def build_axis_curve_inputs(
    axes: list[dict],
    contacts: list[dict],
    precision: float,
) -> tuple[list[tuple[int, int]], list[dict], list[str]]:
    """Create property/contact-aware OCC curve pieces for all reconstructed axes."""

    entities: list[tuple[int, int]] = []
    curve_inputs: list[dict] = []
    blockers: list[str] = []
    contacts_by_axis: dict[int, list[dict]] = defaultdict(list)
    for contact in contacts:
        contacts_by_axis[int(contact["axis"])].append(contact)

    for axis_index, axis in enumerate(axes):
        length = axis_length(axis)
        if not math.isfinite(length) or length <= precision:
            blockers.append(f"degenerate_axis:{axis_index}")
            continue

        cuts = [0.0, 1.0]
        for span in axis.get("property_spans", []):
            cuts.extend((float(span["start_t"]), float(span["end_t"])))
        for anchor in axis.get("anchors", []):
            cuts.append(float(anchor["t"]))
        for contact in contacts_by_axis.get(axis_index, []):
            if contact["kind"] == "point":
                cuts.append(float(contact["t"]))
            elif contact["kind"] == "interval":
                cuts.extend((float(contact["start_t"]), float(contact["end_t"])))

        try:
            cuts = unique_parameters(cuts, length, precision)
        except ValueError:
            blockers.append(f"invalid_axis_parameter:{axis_index}")
            continue

        point_tags = {
            t: gmsh.model.occ.addPoint(*axis_point(axis, t))
            for t in cuts
        }
        tolerance = precision / length
        spans = axis.get("property_spans", [])
        for start_t, end_t in zip(cuts, cuts[1:]):
            if end_t - start_t <= tolerance:
                continue
            midpoint = 0.5 * (start_t + end_t)
            active = [
                span
                for span in spans
                if float(span["start_t"]) - tolerance <= midpoint <= float(span["end_t"]) + tolerance
            ]
            if not active:
                blockers.append(
                    f"axis_interval_without_property:axis={axis_index},start={start_t},end={end_t}"
                )
                continue
            stiffnesses = {int(span["stiffness"]) for span in active}
            if len(stiffnesses) != 1:
                blockers.append(
                    f"axis_interval_property_conflict:axis={axis_index},start={start_t},end={end_t}"
                )
                continue

            tag = gmsh.model.occ.addLine(point_tags[start_t], point_tags[end_t])
            source_elements = sorted(
                {int(span["source_element"]) for span in active}
            )
            curve_inputs.append(
                {
                    "input_curve": len(curve_inputs),
                    "axis": axis_index,
                    "source_axis": int(axis["source_axis"]),
                    "start_t": start_t,
                    "end_t": end_t,
                    "stiffness": next(iter(stiffnesses)),
                    "source_elements": source_elements,
                }
            )
            entities.append((1, tag))

    return entities, curve_inputs, blockers


def element_lines(curve_tag: int) -> list[tuple[int, int]]:
    lines = []
    types, _element_tags, node_tags = gmsh.model.mesh.getElements(1, curve_tag)
    for element_type, nodes in zip(types, node_tags):
        name, _dim, _order, nodes_per_element, _local, primary_nodes = (
            gmsh.model.mesh.getElementProperties(element_type)
        )
        if not name.lower().startswith("line") or int(primary_nodes) < 2:
            raise RuntimeError(
                f"unsupported 1D element on curve {curve_tag}: {name}"
            )
        nodes_per_element = int(nodes_per_element)
        for i in range(0, len(nodes), nodes_per_element):
            lines.append(tuple(int(x) for x in nodes[i : i + 2]))
    return lines


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
) -> tuple[list[dict], dict[int, int], dict]:
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

        moves = []
        for node in nodes:
            if node != representative:
                mapping[node] = representative
                moves.append(
                    {
                        "from": node,
                        "to": representative,
                        "from_coordinate": coords[node],
                        "to_coordinate": coords[representative],
                        "movement": movements[node],
                    }
                )
        maximum_movement = max(maximum_movement, component_max)
        diagnostics.append(
            {
                "nodes": sorted(nodes),
                "representative": representative,
                "representative_source_owner_count": len(owners[representative]),
                "status": "collapsed",
                "maximum_movement": component_max,
                "moves": moves,
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

    return healed, mapping, {
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


def surface_edges_by_source(triangles: list[dict]) -> dict[int, set[tuple[int, int]]]:
    result: dict[int, set[tuple[int, int]]] = defaultdict(set)
    for triangle in triangles:
        a, b, c = triangle["nodes"]
        result[triangle["source_surface"]].update(
            {
                tuple(sorted((a, b))),
                tuple(sorted((b, c))),
                tuple(sorted((c, a))),
            }
        )
    return result


def surface_nodes_by_source(triangles: list[dict]) -> dict[int, set[int]]:
    result: dict[int, set[int]] = defaultdict(set)
    for triangle in triangles:
        result[triangle["source_surface"]].update(triangle["nodes"])
    return result


def axis_parameter(point: Iterable[float], axis: dict) -> tuple[float, float]:
    a = tuple(float(x) for x in axis["endpoints"][0])
    b = tuple(float(x) for x in axis["endpoints"][1])
    p = tuple(float(x) for x in point)
    d = tuple(b[i] - a[i] for i in range(3))
    length2 = sum(x * x for x in d)
    if length2 == 0:
        return 0.0, float("inf")
    t = sum((p[i] - a[i]) * d[i] for i in range(3)) / length2
    q = tuple(a[i] + t * d[i] for i in range(3))
    return t, distance(p, q)


def audit_contacts(
    data: dict,
    coords: dict[int, tuple[float, float, float]],
    triangles: list[dict],
    bars: list[dict],
    precision: float,
) -> dict:
    """Verify that declared reconstructed bar/surface contacts became mesh identity."""

    surface_edges = surface_edges_by_source(triangles)
    surface_nodes = surface_nodes_by_source(triangles)
    axis_nodes: dict[int, set[int]] = defaultdict(set)
    bars_by_axis: dict[int, list[dict]] = defaultdict(list)
    for bar in bars:
        axis_nodes[bar["axis"]].update(bar["nodes"])
        bars_by_axis[int(bar["axis"])].append(bar)

    details = []
    failed = 0
    tolerance = precision * 10.0
    for index, contact in enumerate(data.get("contacts", [])):
        axis_index = int(contact["axis"])
        surface = int(contact["surface"])
        axis = data["axes"][axis_index]

        if contact["kind"] == "point":
            expected = tuple(float(x) for x in contact["point"])
            shared = axis_nodes.get(axis_index, set()) & surface_nodes.get(surface, set())
            nearest = min(
                (distance(coords[node], expected), node) for node in shared
            ) if shared else None
            conforming = nearest is not None and nearest[0] <= tolerance
            if not conforming:
                failed += 1
            details.append(
                {
                    "contact": index,
                    "kind": "point",
                    "axis": axis_index,
                    "surface": surface,
                    "mesh_conforming": conforming,
                    "topologically_shared": nearest is not None,
                    "shared_node": nearest[1] if nearest else None,
                    "distance": nearest[0] if nearest else None,
                }
            )
            continue

        start_t = float(contact["start_t"])
        end_t = float(contact["end_t"])
        expected_length = axis_length(axis) * (end_t - start_t)
        axis_tol = tolerance / max(axis_length(axis), tolerance)
        total = 0.0
        shared_total = 0.0
        for bar in bars_by_axis.get(axis_index, []):
            p0 = coords[bar["nodes"][0]]
            p1 = coords[bar["nodes"][1]]
            midpoint = tuple((p0[i] + p1[i]) * 0.5 for i in range(3))
            t, residual = axis_parameter(midpoint, axis)
            if residual > tolerance:
                continue
            if start_t - axis_tol <= t <= end_t + axis_tol:
                length = distance(p0, p1)
                total += length
                if tuple(sorted(bar["nodes"])) in surface_edges.get(surface, set()):
                    shared_total += length

        length_tolerance = max(tolerance, expected_length * 1e-10)
        conforming = (
            abs(total - expected_length) <= length_tolerance
            and abs(shared_total - expected_length) <= length_tolerance
        )
        if not conforming:
            failed += 1
        details.append(
            {
                "contact": index,
                "kind": "interval",
                "axis": axis_index,
                "surface": surface,
                "mesh_conforming": conforming,
                "expected_length": expected_length,
                "bar_length": total,
                "shared_mesh_edge_length": shared_total,
            }
        )

    return {
        "contact_count": len(details),
        "conforming_contact_count": len(details) - failed,
        "failed_contact_count": failed,
        "details": details,
    }


def summarize_contact_audit(report: dict) -> dict:
    return {
        "contact_count": report["contact_count"],
        "conforming_contact_count": report["conforming_contact_count"],
        "failed_contact_count": report["failed_contact_count"],
    }


def reconcile_healed_contact_audit(
    data: dict,
    raw: dict,
    strict_healed: dict,
    healing: dict,
    precision: float,
    minimum_edge: float,
) -> dict:
    """Accept a moved point contact only through an explicit healing move.

    The contact must have been strictly conforming before healing, must remain
    topologically shared afterwards, and its expected point must coincide with
    the `from` node of the exact logged move whose `to` node is now shared.
    No global tolerance is relaxed.
    """
    tolerance = precision * 10.0
    moves = [
        move
        for component in healing.get("components", [])
        if component.get("status") == "collapsed"
        for move in component.get("moves", [])
    ]
    result = []
    failed = 0
    accepted_by_healing = 0
    if len(raw["details"]) != len(strict_healed["details"]):
        raise RuntimeError("contact audit length changed across healing")

    for raw_item, healed_item in zip(raw["details"], strict_healed["details"]):
        item = dict(healed_item)
        if healed_item["mesh_conforming"]:
            item["conformity"] = "strict"
            result.append(item)
            continue

        accepted = False
        if (
            healed_item["kind"] == "point"
            and raw_item["mesh_conforming"]
            and healed_item.get("topologically_shared")
            and healed_item.get("shared_node") is not None
        ):
            contact = data["contacts"][healed_item["contact"]]
            expected = tuple(float(x) for x in contact["point"])
            for move in moves:
                if int(move["to"]) != int(healed_item["shared_node"]):
                    continue
                if distance(move["from_coordinate"], expected) > tolerance:
                    continue
                movement = float(move["movement"])
                if movement >= minimum_edge:
                    continue
                accepted = True
                item["mesh_conforming"] = True
                item["conformity"] = "healed_shared_node"
                item["healing_movement"] = movement
                item["healing_from_node"] = int(move["from"])
                item["healing_to_node"] = int(move["to"])
                accepted_by_healing += 1
                break

        if not accepted:
            item["conformity"] = "failed"
            failed += 1
        result.append(item)

    return {
        "contact_count": len(result),
        "conforming_contact_count": len(result) - failed,
        "failed_contact_count": failed,
        "accepted_by_healing_count": accepted_by_healing,
        "details": result,
    }


def audit_semantic_shared_nodes(
    data: dict,
    coords: dict[int, tuple[float, float, float]],
    triangles: list[dict],
    bars: list[dict],
    precision: float,
) -> dict:
    """Reject shared mesh nodes that have no reconstructed semantic path.

    Surface/surface identity is justified by a finite shared mesh edge or a
    common source boundary node. Axis/axis identity requires a common source
    anchor node. Axis/surface identity requires an explicit Rust contact.
    Transitive paths through those relations are valid; coordinate proximity
    alone is never a relation.
    """
    tolerance = precision * 10.0
    node_surfaces: dict[int, set[int]] = defaultdict(set)
    node_axes: dict[int, set[int]] = defaultdict(set)
    edge_surfaces: dict[tuple[int, int], set[int]] = defaultdict(set)

    for triangle in triangles:
        surface = int(triangle["source_surface"])
        a, b, c = triangle["nodes"]
        node_surfaces[a].add(surface)
        node_surfaces[b].add(surface)
        node_surfaces[c].add(surface)
        for edge in (
            tuple(sorted((a, b))),
            tuple(sorted((b, c))),
            tuple(sorted((c, a))),
        ):
            edge_surfaces[edge].add(surface)

    for bar in bars:
        axis = int(bar["axis"])
        for node in bar["nodes"]:
            node_axes[node].add(axis)

    finite_surface_pairs_at_node: dict[int, set[tuple[int, int]]] = defaultdict(set)
    for edge, owners in edge_surfaces.items():
        owners = sorted(owners)
        for i in range(len(owners)):
            for j in range(i + 1, len(owners)):
                pair = (owners[i], owners[j])
                finite_surface_pairs_at_node[edge[0]].add(pair)
                finite_surface_pairs_at_node[edge[1]].add(pair)

    surface_source_points: list[dict[int, list[tuple[float, float, float]]]] = []
    for surface in data.get("surfaces", []):
        mapping: dict[int, list[tuple[float, float, float]]] = defaultdict(list)
        for ring, source_nodes in zip(
            surface.get("rings", []), surface.get("ring_source_nodes", [])
        ):
            for point, source_node in zip(ring, source_nodes):
                mapping[int(source_node)].append(tuple(float(x) for x in point))
        surface_source_points.append(mapping)

    axis_source_points: list[dict[int, tuple[float, float, float]]] = []
    for axis in data.get("axes", []):
        axis_source_points.append(
            {
                int(anchor["source_node"]): tuple(float(x) for x in anchor["point"])
                for anchor in axis.get("anchors", [])
            }
        )

    contacts_by_pair: dict[tuple[int, int], list[dict]] = defaultdict(list)
    for contact in data.get("contacts", []):
        contacts_by_pair[(int(contact["axis"]), int(contact["surface"]))].append(contact)

    def relation(
        node: int,
        left: tuple[str, int],
        right: tuple[str, int],
    ) -> tuple[bool, str | None]:
        point = coords[node]
        if left[0] == "surface" and right[0] == "surface":
            s1, s2 = left[1], right[1]
            pair = tuple(sorted((s1, s2)))
            if pair in finite_surface_pairs_at_node.get(node, set()):
                return True, "shared_surface_edge"
            common = set(surface_source_points[s1]) & set(surface_source_points[s2])
            for source_node in common:
                if any(
                    distance(candidate, point) <= tolerance
                    for candidate in surface_source_points[s1][source_node]
                ) and any(
                    distance(candidate, point) <= tolerance
                    for candidate in surface_source_points[s2][source_node]
                ):
                    return True, f"shared_surface_source_node:{source_node}"
            return False, None

        if left[0] == "axis" and right[0] == "axis":
            a1, a2 = left[1], right[1]
            common = set(axis_source_points[a1]) & set(axis_source_points[a2])
            for source_node in common:
                if (
                    distance(axis_source_points[a1][source_node], point) <= tolerance
                    and distance(axis_source_points[a2][source_node], point) <= tolerance
                ):
                    return True, f"shared_axis_source_node:{source_node}"
            return False, None

        axis = left[1] if left[0] == "axis" else right[1]
        surface = left[1] if left[0] == "surface" else right[1]
        geometry = data["axes"][axis]
        parameter, residual = axis_parameter(point, geometry)
        parameter_tolerance = tolerance / max(axis_length(geometry), tolerance)
        for contact in contacts_by_pair.get((axis, surface), []):
            if contact["kind"] == "point":
                if distance(contact["point"], point) <= tolerance:
                    return True, "declared_point_contact"
            elif (
                residual <= tolerance
                and float(contact["start_t"]) - parameter_tolerance
                <= parameter
                <= float(contact["end_t"]) + parameter_tolerance
            ):
                return True, "declared_interval_contact"
        return False, None

    shared_node_count = 0
    unintended = []
    direct_axis_axis = 0
    mediated_axis_axis = 0
    unsupported_axis_axis = 0
    mediated_examples = []

    for node in sorted(set(node_surfaces) | set(node_axes)):
        owners = [
            ("surface", surface) for surface in sorted(node_surfaces[node])
        ] + [("axis", axis) for axis in sorted(node_axes[node])]
        if len(owners) <= 1:
            continue
        shared_node_count += 1

        adjacency = {owner: set() for owner in owners}
        direct_relations: dict[tuple[tuple[str, int], tuple[str, int]], str] = {}
        for i, left in enumerate(owners):
            for right in owners[i + 1 :]:
                connected, reason = relation(node, left, right)
                if connected:
                    adjacency[left].add(right)
                    adjacency[right].add(left)
                    direct_relations[(left, right)] = reason or "semantic_relation"

        components = []
        owner_component = {}
        seen = set()
        for owner in owners:
            if owner in seen:
                continue
            stack = [owner]
            seen.add(owner)
            component = []
            while stack:
                current = stack.pop()
                component.append(current)
                for other in adjacency[current]:
                    if other not in seen:
                        seen.add(other)
                        stack.append(other)
            component_id = len(components)
            for item in component:
                owner_component[item] = component_id
            components.append(component)

        axes_here = sorted(node_axes[node])
        for i, a1 in enumerate(axes_here):
            for a2 in axes_here[i + 1 :]:
                left = ("axis", a1)
                right = ("axis", a2)
                key = (left, right) if (left, right) in direct_relations else (right, left)
                if key in direct_relations:
                    direct_axis_axis += 1
                elif owner_component[left] == owner_component[right]:
                    mediated_axis_axis += 1
                    if len(mediated_examples) < 20:
                        mediated_examples.append(
                            {"node": node, "axes": [a1, a2], "coordinate": coords[node]}
                        )
                else:
                    unsupported_axis_axis += 1

        if len(components) > 1:
            unintended.append(
                {
                    "node": node,
                    "coordinate": coords[node],
                    "components": components,
                }
            )

    return {
        "shared_node_count": shared_node_count,
        "unintended_shared_node_count": len(unintended),
        "direct_axis_axis_shared_node_count": direct_axis_axis,
        "surface_mediated_axis_axis_shared_node_count": mediated_axis_axis,
        "unsupported_axis_axis_shared_node_count": unsupported_axis_axis,
        "surface_mediated_axis_axis_examples": mediated_examples,
        "unintended_shared_nodes": unintended,
    }


@dataclass
class Fragmentation:
    surface_to_output: list[list[int]]
    surface_output_owners: dict[int, list[int]]
    output_surfaces: list[int]
    curve_to_output: list[list[int]]
    curve_output_owners: dict[int, list[int]]
    output_curves: list[int]
    curve_inputs: list[dict]
    axis_build_blockers: list[str]


def fragment(
    input_surfaces: list[dict],
    axes: list[dict],
    contacts: list[dict],
    precision: float,
) -> Fragmentation:
    surface_entities = [(2, add_surface(surface)) for surface in input_surfaces]
    curve_entities, curve_inputs, axis_build_blockers = build_axis_curve_inputs(
        axes, contacts, precision
    )
    input_entities = surface_entities + curve_entities

    if not input_entities:
        return Fragmentation([], {}, [], [], {}, [], curve_inputs, axis_build_blockers)

    if len(input_entities) == 1:
        output_map = [[input_entities[0]]]
        output = [input_entities[0]]
    else:
        output, output_map = gmsh.model.occ.fragment(
            [input_entities[0]],
            input_entities[1:],
            removeObject=True,
            removeTool=True,
        )

    surface_count = len(surface_entities)
    surface_map_raw = output_map[:surface_count]
    curve_map_raw = output_map[surface_count:]

    surface_to_output = [
        sorted({tag for dim, tag in mapped if dim == 2}) for mapped in surface_map_raw
    ]
    curve_to_output = [
        sorted({tag for dim, tag in mapped if dim == 1}) for mapped in curve_map_raw
    ]
    output_surfaces = sorted(
        {tag for mapped in surface_to_output for tag in mapped}
        | {tag for dim, tag in output if dim == 2}
    )
    output_curves = sorted({tag for mapped in curve_to_output for tag in mapped})

    gmsh.model.occ.synchronize()

    surface_output_owners: dict[int, list[int]] = defaultdict(list)
    for source, tags in enumerate(surface_to_output):
        for tag in tags:
            surface_output_owners[tag].append(source)

    curve_output_owners: dict[int, list[int]] = defaultdict(list)
    for source, tags in enumerate(curve_to_output):
        for tag in tags:
            curve_output_owners[tag].append(source)

    return Fragmentation(
        surface_to_output=surface_to_output,
        surface_output_owners={
            tag: sorted(set(owners)) for tag, owners in surface_output_owners.items()
        },
        output_surfaces=output_surfaces,
        curve_to_output=curve_to_output,
        curve_output_owners={
            tag: sorted(set(owners)) for tag, owners in curve_output_owners.items()
        },
        output_curves=output_curves,
        curve_inputs=curve_inputs,
        axis_build_blockers=axis_build_blockers,
    )


def entity_point(tag: int) -> tuple[float, float, float]:
    value = gmsh.model.getValue(0, tag, [])
    return (float(value[0]), float(value[1]), float(value[2]))


def curve_midpoint(tag: int) -> tuple[float, float, float]:
    low_x, low_y, low_z, high_x, high_y, high_z = gmsh.model.getBoundingBox(1, tag)
    return (
        0.5 * (low_x + high_x),
        0.5 * (low_y + high_y),
        0.5 * (low_z + high_z),
    )


def embed_declared_contacts(
    fragmentation: Fragmentation,
    data: dict,
    precision: float,
) -> dict:
    """Mesh-embed only contacts already established by Rust semantics."""
    tolerance = precision * 10.0
    curve_inputs_by_axis: dict[int, list[int]] = defaultdict(list)
    for index, item in enumerate(fragmentation.curve_inputs):
        curve_inputs_by_axis[int(item["axis"])].append(index)

    # Point contacts already lying inside a declared interval contact need no
    # separate 0D embed: embedding the interval curve makes every split point
    # on that curve a surface mesh node by identity.
    interval_coverage: dict[tuple[int, int], list[tuple[float, float]]] = defaultdict(list)
    for contact in data.get("contacts", []):
        if contact["kind"] == "interval":
            interval_coverage[(int(contact["axis"]), int(contact["surface"]))].append(
                (float(contact["start_t"]), float(contact["end_t"]))
            )

    embedded_curves: set[tuple[int, int]] = set()
    embedded_points: set[tuple[int, int]] = set()
    existing_boundary_curves = 0
    existing_boundary_points = 0
    covered_point_contacts = 0
    blockers: list[str] = []

    boundary_curves: dict[int, set[int]] = {}
    boundary_points: dict[int, set[int]] = {}
    for source_faces in fragmentation.surface_to_output:
        for face in source_faces:
            if face in boundary_curves:
                continue
            boundary_curves[face] = {
                int(tag)
                for dim, tag in gmsh.model.getBoundary(
                    [(2, face)], combined=False, oriented=False, recursive=False
                )
                if dim == 1
            }
            points = set()
            for curve in boundary_curves[face]:
                points.update(
                    int(tag)
                    for dim, tag in gmsh.model.getBoundary(
                        [(1, curve)], combined=False, oriented=False, recursive=False
                    )
                    if dim == 0
                )
            boundary_points[face] = points

    # Interval contacts first. The curve pieces are already split at property,
    # anchor and contact parameters; General Fuse may split them further at
    # surface intersections, so each output curve can be assigned by midpoint.
    for contact_index, contact in enumerate(data.get("contacts", [])):
        if contact["kind"] != "interval":
            continue
        axis = int(contact["axis"])
        surface = int(contact["surface"])
        faces = fragmentation.surface_to_output[surface]
        axis_length_value = axis_length(data["axes"][axis])
        param_tol = tolerance / max(axis_length_value, tolerance)
        start_t = float(contact["start_t"])
        end_t = float(contact["end_t"])
        selected = []
        for input_curve in curve_inputs_by_axis.get(axis, []):
            item = fragmentation.curve_inputs[input_curve]
            if (
                float(item["start_t"]) >= start_t - param_tol
                and float(item["end_t"]) <= end_t + param_tol
            ):
                selected.extend(fragmentation.curve_to_output[input_curve])
        if not selected:
            blockers.append(f"interval_contact_without_curve:{contact_index}")
            continue

        for curve in sorted(set(selected)):
            midpoint = curve_midpoint(curve)
            owners = []
            for face in faces:
                if curve in boundary_curves[face]:
                    existing_boundary_curves += 1
                    owners.append(face)
                    continue
                # Rust already established this curve/surface contact. isInside
                # selects only the OCC fragment containing this exact piece; no
                # nearest-surface or fuzzy geometric inference is performed.
                if gmsh.model.isInside(2, face, midpoint) > 0:
                    gmsh.model.mesh.embed(1, [curve], 2, face)
                    embedded_curves.add((curve, face))
                    owners.append(face)
            if not owners:
                blockers.append(
                    f"interval_curve_not_on_contact_surface:contact={contact_index},curve={curve}"
                )

    # Point-only contacts. Use the actual OCC endpoint created by the split
    # axis curve and embed it only in the Rust-declared source surface.
    for contact_index, contact in enumerate(data.get("contacts", [])):
        if contact["kind"] != "point":
            continue
        axis = int(contact["axis"])
        surface = int(contact["surface"])
        t = float(contact["t"])
        axis_length_value = axis_length(data["axes"][axis])
        param_tol = tolerance / max(axis_length_value, tolerance)
        if any(
            start - param_tol <= t <= end + param_tol
            for start, end in interval_coverage.get((axis, surface), [])
        ):
            covered_point_contacts += 1
            continue

        faces = fragmentation.surface_to_output[surface]
        expected = tuple(float(x) for x in contact["point"])
        candidates = set()
        for input_curve in curve_inputs_by_axis.get(axis, []):
            item = fragmentation.curve_inputs[input_curve]
            if (
                abs(float(item["start_t"]) - t) > param_tol
                and abs(float(item["end_t"]) - t) > param_tol
            ):
                continue
            for curve in fragmentation.curve_to_output[input_curve]:
                candidates.update(
                    int(tag)
                    for dim, tag in gmsh.model.getBoundary(
                        [(1, curve)], combined=False, oriented=False, recursive=False
                    )
                    if dim == 0 and distance(entity_point(int(tag)), expected) <= tolerance
                )
        if not candidates:
            blockers.append(f"point_contact_without_axis_point:{contact_index}")
            continue

        placed = False
        for point in sorted(candidates):
            xyz = entity_point(point)
            for face in faces:
                if point in boundary_points[face]:
                    existing_boundary_points += 1
                    placed = True
                    continue
                if gmsh.model.isInside(2, face, xyz) > 0:
                    gmsh.model.mesh.embed(0, [point], 2, face)
                    embedded_points.add((point, face))
                    placed = True
        if not placed:
            blockers.append(f"point_contact_not_on_contact_surface:{contact_index}")

    return {
        "embedded_curve_count": len(embedded_curves),
        "embedded_point_count": len(embedded_points),
        "covered_point_contact_count": covered_point_contacts,
        "existing_boundary_curve_count": existing_boundary_curves,
        "existing_boundary_point_count": existing_boundary_points,
        "blockers": blockers,
    }


def physical_groups(
    fragmentation: Fragmentation,
    surfaces: list[dict],
) -> list[dict]:
    surface_by_stiffness: dict[int, list[int]] = defaultdict(list)
    for tag in fragmentation.output_surfaces:
        owners = fragmentation.surface_output_owners.get(tag, [])
        if len(owners) != 1:
            continue
        surface_by_stiffness[int(surfaces[owners[0]]["stiffness"])].append(tag)

    curve_by_stiffness: dict[int, list[int]] = defaultdict(list)
    for tag in fragmentation.output_curves:
        owners = fragmentation.curve_output_owners.get(tag, [])
        if len(owners) != 1:
            continue
        curve_by_stiffness[
            int(fragmentation.curve_inputs[owners[0]]["stiffness"])
        ].append(tag)

    result = []
    for stiffness, tags in sorted(surface_by_stiffness.items()):
        group = gmsh.model.addPhysicalGroup(2, sorted(set(tags)))
        name = f"surface_stiffness_{stiffness}"
        gmsh.model.setPhysicalName(2, group, name)
        result.append(
            {
                "dimension": 2,
                "physical_tag": group,
                "name": name,
                "stiffness": stiffness,
                "entity_tags": sorted(set(tags)),
            }
        )

    for stiffness, tags in sorted(curve_by_stiffness.items()):
        group = gmsh.model.addPhysicalGroup(1, sorted(set(tags)))
        name = f"bar_stiffness_{stiffness}"
        gmsh.model.setPhysicalName(1, group, name)
        result.append(
            {
                "dimension": 1,
                "physical_tag": group,
                "name": name,
                "stiffness": stiffness,
                "entity_tags": sorted(set(tags)),
            }
        )

    return result


def build_solver_mesh_package(
    data: dict,
    fragmentation: Fragmentation,
    vertices: list[tuple[float, float, float]],
    triangles: list[dict],
    bars: list[dict],
    precision: float,
) -> tuple[dict, dict]:
    """Build and audit the backend-neutral mesh consumed by solver adapters."""

    surface_regions = [
        {
            "region": int(surface["source_surface"]),
            "source_patch": int(surface["source_patch"]),
            "stiffness": int(surface["stiffness"]),
            "source_elements": sorted({int(x) for x in surface["source_elements"]}),
        }
        for surface in data["surfaces"]
    ]
    bar_regions = [
        {
            "region": int(item["input_curve"]),
            "axis": int(item["axis"]),
            "source_axis": int(item["source_axis"]),
            "stiffness": int(item["stiffness"]),
            "source_elements": sorted({int(x) for x in item["source_elements"]}),
            "start_t": float(item["start_t"]),
            "end_t": float(item["end_t"]),
        }
        for item in fragmentation.curve_inputs
    ]

    shell_elements = [
        {
            "vertices": [int(x) for x in triangle["vertices"]],
            "region": int(triangle["source_surface"]),
        }
        for triangle in triangles
    ]
    bar_elements = [
        {
            "vertices": [int(x) for x in bar["vertices"]],
            "region": int(bar["input_curve"]),
        }
        for bar in bars
    ]

    surface_region_by_id = {item["region"]: item for item in surface_regions}
    bar_region_by_id = {item["region"]: item for item in bar_regions}
    vertex_count = len(vertices)

    invalid_shell_indices = []
    degenerate_shell_elements = []
    shell_stiffness_mismatches = []
    represented_surface_regions = set()
    for index, (source, element) in enumerate(zip(triangles, shell_elements)):
        ids = element["vertices"]
        if len(ids) != 3 or any(node < 0 or node >= vertex_count for node in ids):
            invalid_shell_indices.append(index)
            continue
        if len(set(ids)) != 3 or triangle_area(*(vertices[node] for node in ids)) <= precision * precision:
            degenerate_shell_elements.append(index)
        region = element["region"]
        represented_surface_regions.add(region)
        expected = surface_region_by_id.get(region)
        if expected is None or int(source["stiffness"]) != expected["stiffness"]:
            shell_stiffness_mismatches.append(index)

    invalid_bar_indices = []
    degenerate_bar_elements = []
    bar_stiffness_mismatches = []
    represented_bar_regions = set()
    for index, (source, element) in enumerate(zip(bars, bar_elements)):
        ids = element["vertices"]
        if len(ids) != 2 or any(node < 0 or node >= vertex_count for node in ids):
            invalid_bar_indices.append(index)
            continue
        if len(set(ids)) != 2 or distance(vertices[ids[0]], vertices[ids[1]]) <= precision:
            degenerate_bar_elements.append(index)
        region = element["region"]
        represented_bar_regions.add(region)
        expected = bar_region_by_id.get(region)
        if expected is None or int(source["stiffness"]) != expected["stiffness"]:
            bar_stiffness_mismatches.append(index)

    expected_surface_regions = set(surface_region_by_id)
    expected_bar_regions = set(bar_region_by_id)
    missing_surface_regions = sorted(expected_surface_regions - represented_surface_regions)
    unexpected_surface_regions = sorted(represented_surface_regions - expected_surface_regions)
    missing_bar_regions = sorted(expected_bar_regions - represented_bar_regions)
    unexpected_bar_regions = sorted(represented_bar_regions - expected_bar_regions)

    expected_surface_elements = {
        source_element
        for region in surface_regions
        for source_element in region["source_elements"]
    }
    represented_surface_elements = {
        source_element
        for region_id in represented_surface_regions
        if region_id in surface_region_by_id
        for source_element in surface_region_by_id[region_id]["source_elements"]
    }
    expected_bar_elements = {
        source_element
        for region in bar_regions
        for source_element in region["source_elements"]
    }
    represented_bar_elements = {
        source_element
        for region_id in represented_bar_regions
        if region_id in bar_region_by_id
        for source_element in bar_region_by_id[region_id]["source_elements"]
    }

    audit = {
        "surface_region_count": len(surface_regions),
        "represented_surface_region_count": len(represented_surface_regions),
        "bar_region_count": len(bar_regions),
        "represented_bar_region_count": len(represented_bar_regions),
        "expected_surface_source_element_count": len(expected_surface_elements),
        "represented_surface_source_element_count": len(represented_surface_elements),
        "expected_bar_source_element_count": len(expected_bar_elements),
        "represented_bar_source_element_count": len(represented_bar_elements),
        "missing_surface_regions": missing_surface_regions,
        "unexpected_surface_regions": unexpected_surface_regions,
        "missing_bar_regions": missing_bar_regions,
        "unexpected_bar_regions": unexpected_bar_regions,
        "missing_surface_source_elements": sorted(
            expected_surface_elements - represented_surface_elements
        ),
        "unexpected_surface_source_elements": sorted(
            represented_surface_elements - expected_surface_elements
        ),
        "missing_bar_source_elements": sorted(
            expected_bar_elements - represented_bar_elements
        ),
        "unexpected_bar_source_elements": sorted(
            represented_bar_elements - expected_bar_elements
        ),
        "invalid_shell_element_indices": invalid_shell_indices,
        "invalid_bar_element_indices": invalid_bar_indices,
        "degenerate_shell_element_indices": degenerate_shell_elements,
        "degenerate_bar_element_indices": degenerate_bar_elements,
        "shell_stiffness_mismatch_indices": shell_stiffness_mismatches,
        "bar_stiffness_mismatch_indices": bar_stiffness_mismatches,
    }
    audit["clean"] = not any(
        value
        for key, value in audit.items()
        if key not in {
            "clean",
            "surface_region_count",
            "represented_surface_region_count",
            "bar_region_count",
            "represented_bar_region_count",
            "expected_surface_source_element_count",
            "represented_surface_source_element_count",
            "expected_bar_source_element_count",
            "represented_bar_source_element_count",
        }
    )

    package = {
        "format": SOLVER_MESH_FORMAT,
        "length_unit": data.get("length_unit", "model_unit"),
        "index_base": 0,
        "vertices": vertices,
        "surface_regions": surface_regions,
        "bar_regions": bar_regions,
        "shell_elements": shell_elements,
        "bar_elements": bar_elements,
    }
    return package, audit


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

    fragmentation = fragment(
        surfaces,
        data.get("axes", []),
        data.get("contacts", []),
        precision,
    )

    surface_ownership_conflicts = []
    fragmented_surfaces = []
    unmapped_output_surfaces = []
    for tag in fragmentation.output_surfaces:
        owners = fragmentation.surface_output_owners.get(tag, [])
        if not owners:
            unmapped_output_surfaces.append(tag)
            continue
        if len(owners) > 1:
            stiffnesses = sorted({int(surfaces[owner]["stiffness"]) for owner in owners})
            surface_ownership_conflicts.append(
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

    curve_ownership_conflicts = []
    fragmented_curves = []
    for tag in fragmentation.output_curves:
        owners = fragmentation.curve_output_owners.get(tag, [])
        if len(owners) != 1:
            curve_ownership_conflicts.append(
                {
                    "output_curve": tag,
                    "input_curves": owners,
                    "reason": "unmapped" if not owners else "multiple_input_curves",
                }
            )
            continue
        source = owners[0]
        item = fragmentation.curve_inputs[source]
        fragmented_curves.append(
            {
                "output_curve": tag,
                **item,
            }
        )

    contact_embedding = embed_declared_contacts(fragmentation, data, precision)
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
    bars = []
    non_triangle_blockers = []
    non_line_blockers = []
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

    for item in fragmented_curves:
        tag = item["output_curve"]
        try:
            local = element_lines(tag)
        except RuntimeError as error:
            non_line_blockers.append(str(error))
            continue
        for nodes in local:
            bars.append(
                {
                    "nodes": nodes,
                    "output_curve": tag,
                    "input_curve": item["input_curve"],
                    "axis": item["axis"],
                    "source_axis": item["source_axis"],
                    "start_t": item["start_t"],
                    "end_t": item["end_t"],
                    "stiffness": item["stiffness"],
                    "source_elements": item["source_elements"],
                }
            )

    raw_quality = quality(coords, triangles)
    raw_shared = shared_mesh_edges_by_source(triangles)
    raw_contact_audit = audit_contacts(data, coords, triangles, bars, precision)
    semantic_node_audit = audit_semantic_shared_nodes(
        data, coords, triangles, bars, precision
    )
    healed_triangles, node_mapping, healing = heal_micro_edges(
        coords, triangles, minimum_edge, precision
    )

    healed_bars = []
    removed_degenerate_bars = 0
    for bar in bars:
        nodes = tuple(node_mapping.get(node, node) for node in bar["nodes"])
        if nodes[0] == nodes[1]:
            removed_degenerate_bars += 1
            continue
        healed_bars.append({**bar, "nodes": nodes})
    healing["removed_degenerate_bars"] = removed_degenerate_bars

    healed_quality = quality(coords, healed_triangles)
    healed_shared = shared_mesh_edges_by_source(healed_triangles)
    strict_healed_contact_audit = audit_contacts(
        data, coords, healed_triangles, healed_bars, precision
    )
    contact_audit = reconcile_healed_contact_audit(
        data,
        raw_contact_audit,
        strict_healed_contact_audit,
        healing,
        precision,
        minimum_edge,
    )

    used_nodes = sorted(
        {node for triangle in healed_triangles for node in triangle["nodes"]}
        | {node for bar in healed_bars for node in bar["nodes"]}
    )
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
    compact_bars = [
        {
            "vertices": [compact_index[node] for node in bar["nodes"]],
            "output_curve": bar["output_curve"],
            "input_curve": bar["input_curve"],
            "axis": bar["axis"],
            "source_axis": bar["source_axis"],
            "stiffness": bar["stiffness"],
            "source_elements": bar["source_elements"],
        }
        for bar in healed_bars
    ]

    solver_mesh, solver_mesh_audit = build_solver_mesh_package(
        data,
        fragmentation,
        compact_vertices,
        compact_triangles,
        compact_bars,
        precision,
    )

    blockers = list(data.get("blockers", []))
    blockers.extend(fragmentation.axis_build_blockers)
    blockers.extend(contact_embedding["blockers"])
    if surface_ownership_conflicts:
        blockers.append("ambiguous_surface_fragment_ownership")
    if curve_ownership_conflicts:
        blockers.append("ambiguous_curve_fragment_ownership")
    if unmapped_output_surfaces:
        blockers.append("unmapped_fragment_surface")
    if non_triangle_blockers:
        blockers.append("unsupported_2d_elements")
    if non_line_blockers:
        blockers.append("unsupported_1d_elements")
    if healing["unresolved_component_count"]:
        blockers.append("unresolved_micro_edge_component")
    if healing["near_degenerate_triangles_after"]:
        blockers.append("near_degenerate_triangle_after_healing")
    if contact_audit["failed_contact_count"]:
        blockers.append("nonconforming_bar_surface_contact")
    if semantic_node_audit["unintended_shared_node_count"]:
        blockers.append("unintended_shared_mesh_node")
    if not solver_mesh_audit["clean"]:
        blockers.append("solver_mesh_coverage_failure")

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
            len(tags) for tags in fragmentation.surface_to_output
        ],
        "input_axis_count": len(data.get("axes", [])),
        "input_curve_piece_count": len(fragmentation.curve_inputs),
        "output_curve_piece_count": len(fragmentation.output_curves),
        "source_to_output_curve_count": [
            len(tags) for tags in fragmentation.curve_to_output
        ],
        "fragmented_surfaces": fragmented_surfaces,
        "fragmented_curves": fragmented_curves,
        "surface_ownership_conflicts": surface_ownership_conflicts,
        "curve_ownership_conflicts": curve_ownership_conflicts,
        "unmapped_output_surfaces": unmapped_output_surfaces,
        "physical_groups": groups,
        "contact_embedding": contact_embedding,
        "bar_surface_contacts_before_healing": summarize_contact_audit(raw_contact_audit),
        "semantic_shared_nodes_before_healing": semantic_node_audit,
        "bar_surface_contacts_strict_after_healing": summarize_contact_audit(
            strict_healed_contact_audit
        ),
        "raw_mesh": {
            **raw_quality,
            "shared_surface_pair_count": len(raw_shared),
            "shared_mesh_edge_count": sum(raw_shared.values()),
        },
        "healing": healing,
        "bar_surface_contacts": contact_audit,
        "solver_mesh_audit": solver_mesh_audit,
        "solver_mesh": solver_mesh,
        "mesh": {
            **healed_quality,
            "bar_count": len(compact_bars),
            "shared_surface_pair_count": len(healed_shared),
            "shared_mesh_edge_count": sum(healed_shared.values()),
            "vertices": compact_vertices,
            "triangles": compact_triangles,
            "bars": compact_bars,
        },
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
            "ring_source_nodes": [[1, 2, 3, 4]],
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
            "ring_source_nodes": [[5, 6, 7, 8]],
        },
    ]
    axes = [
        {
            "source_axis": 0,
            "endpoints": [[-1.5, 0.0, 0.0], [1.5, 0.0, 0.0]],
            "property_spans": [
                {
                    "source_element": 100,
                    "stiffness": 50,
                    "start_t": 0.0,
                    "end_t": 1.0,
                }
            ],
            "anchors": [
                {"source_node": 101, "t": 0.0, "point": [-1.5, 0.0, 0.0]},
                {"source_node": 102, "t": 1.0, "point": [1.5, 0.0, 0.0]},
            ],
        }
    ]
    contacts = [
        {
            "kind": "interval",
            "axis": 0,
            "surface": 0,
            "start_t": 0.0,
            "end_t": 1.0,
            "endpoints": [[-1.5, 0.0, 0.0], [1.5, 0.0, 0.0]],
            "location": "interior",
        },
        {
            "kind": "point",
            "axis": 0,
            "surface": 1,
            "t": 0.5,
            "point": [0.0, 0.0, 0.0],
            "location": "interior",
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
        "axes": axes,
        "contacts": contacts,
        "blockers": [],
    }


def self_test() -> dict:
    result = run_backend(synthetic_report())
    assert result["input_surface_count"] == 2
    assert result["output_surface_count"] >= 3
    assert result["mesh"]["triangle_count"] > 0
    assert result["mesh"]["shared_surface_pair_count"] >= 1
    assert result["mesh"]["shared_mesh_edge_count"] >= 1
    assert result["mesh"]["bar_count"] > 0
    assert result["bar_surface_contacts"]["contact_count"] == 2
    assert result["bar_surface_contacts"]["failed_contact_count"] == 0
    assert result["semantic_shared_nodes_before_healing"]["unintended_shared_node_count"] == 0
    assert result["solver_mesh"]["format"] == SOLVER_MESH_FORMAT
    assert result["solver_mesh_audit"]["clean"], result["solver_mesh_audit"]
    assert len(result["solver_mesh"]["surface_regions"]) == 2
    assert len(result["solver_mesh"]["bar_regions"]) >= 1
    assert len(result["solver_mesh"]["shell_elements"]) == result["mesh"]["triangle_count"]
    assert len(result["solver_mesh"]["bar_elements"]) == result["mesh"]["bar_count"]
    assert not result["surface_ownership_conflicts"]
    assert not result["curve_ownership_conflicts"]
    assert any(
        group["dimension"] == 1 and group["stiffness"] == 50
        for group in result["physical_groups"]
    )
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
    healed, mapping, report = heal_micro_edges(coords, triangles, 0.001, 1e-8)
    assert report["collapsed_component_count"] == 1
    assert report["unresolved_component_count"] == 0
    assert report["maximum_node_movement"] < 0.001
    assert mapping == {2: 1, 3: 1}
    assert all(1 in triangle["nodes"] for triangle in healed)

    fake_data = {"contacts": [{"kind": "point", "point": list(coords[2])}]}
    raw_audit = {
        "details": [
            {
                "contact": 0,
                "kind": "point",
                "mesh_conforming": True,
                "topologically_shared": True,
                "shared_node": 2,
                "distance": 0.0,
            }
        ]
    }
    strict_audit = {
        "details": [
            {
                "contact": 0,
                "kind": "point",
                "mesh_conforming": False,
                "topologically_shared": True,
                "shared_node": 1,
                "distance": distance(coords[1], coords[2]),
            }
        ]
    }
    reconciled = reconcile_healed_contact_audit(
        fake_data, raw_audit, strict_audit, report, 1e-8, 0.001
    )
    assert reconciled["failed_contact_count"] == 0
    assert reconciled["accepted_by_healing_count"] == 1
    assert reconciled["details"][0]["conformity"] == "healed_shared_node"

    unsupported = audit_semantic_shared_nodes(
        {
            "surfaces": [],
            "axes": [
                {
                    "endpoints": [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
                    "anchors": [{"source_node": 1, "point": [0.0, 0.0, 0.0]}],
                },
                {
                    "endpoints": [[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
                    "anchors": [{"source_node": 2, "point": [0.0, 0.0, 0.0]}],
                },
            ],
            "contacts": [],
        },
        {1: (0.0, 0.0, 0.0), 2: (1.0, 0.0, 0.0), 3: (0.0, 1.0, 0.0)},
        [],
        [
            {"nodes": (1, 2), "axis": 0},
            {"nodes": (1, 3), "axis": 1},
        ],
        1e-8,
    )
    assert unsupported["unintended_shared_node_count"] == 1
    assert unsupported["unsupported_axis_axis_shared_node_count"] == 1
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", nargs="?", type=Path)
    parser.add_argument("--output", type=Path, help="summary/full backend JSON")
    parser.add_argument(
        "--solver-mesh-output",
        type=Path,
        help="write only the backend-neutral topo-reconstruct-solver-mesh-v1 package",
    )
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
        if args.solver_mesh_output:
            solver_mesh = result.get("solver_mesh")
            audit = result.get("solver_mesh_audit", {})
            if solver_mesh is None:
                raise RuntimeError("backend result does not contain solver_mesh")
            if not audit.get("clean", False):
                raise RuntimeError(
                    "refusing to write solver mesh because solver_mesh_audit is not clean"
                )
            args.solver_mesh_output.write_text(
                json.dumps(solver_mesh, indent=2)
            )
        print(text)
        return 0
    finally:
        gmsh.finalize()


if __name__ == "__main__":
    raise SystemExit(main())
