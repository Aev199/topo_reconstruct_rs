#!/usr/bin/env python3
"""Prepare a MIDAS GTS NX import plan from topo-reconstruct solver mesh.

This adapter intentionally stops one layer before writing FPN records. Public
GTS NX material confirms that FPN is the native text neutral exchange format,
and public parsers confirm NODE/LINE/TRIA plus MSET/MSETE/MSETN records, but the
exact current LINE/TRIA field layout is not sufficiently documented to treat a
guessed writer as production-safe.

The deterministic plan produced here freezes everything that is independent of
that last serialization detail: 1-based node IDs, global element IDs, element
connectivity, stiffness grouping and source-region provenance. One small FPN
export from the target GTS NX version is enough to calibrate the final writer.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from collections import defaultdict
from pathlib import Path
from typing import Iterable

SOLVER_MESH_FORMAT = "topo-reconstruct-solver-mesh-v1"
PLAN_FORMAT = "topo-reconstruct-midas-gts-plan-v1"
INSPECTION_FORMAT = "topo-reconstruct-fpn-inspection-v1"

# These facts are deliberately separated from unknown serialization details.
OBSERVED_GTS_NX = {
    "neutral_extension": ".fpn",
    "text_format": True,
    "record_families": ["NODE", "LINE", "TRIA", "MSET", "MSETE", "MSETN"],
    "observed_field_positions": {
        "NODE": {"id": 1, "x": 2, "y": 3, "z": 4},
        "LINE": {"node_fields": [3, 4]},
        "TRIA": {"node_fields": [3, 4, 5]},
        "RECT": {"node_fields": [3, 4, 5, 6]},
        "MSET": {"name": 2},
        "MSETE": {"count": 2},
        "MSETN": {"count": 2},
    },
    "writer_verified": False,
}


class ValidationError(ValueError):
    pass


def _finite_triplet(value: object, label: str) -> tuple[float, float, float]:
    if not isinstance(value, list) or len(value) != 3:
        raise ValidationError(f"{label} must contain exactly 3 coordinates")
    result = tuple(float(x) for x in value)
    if not all(math.isfinite(x) for x in result):
        raise ValidationError(f"{label} contains a non-finite coordinate")
    return result


def _integer(value: object, label: str, minimum: int | None = None) -> int:
    if isinstance(value, bool):
        raise ValidationError(f"{label} must be an integer")
    try:
        result = int(value)
    except (TypeError, ValueError) as error:
        raise ValidationError(f"{label} must be an integer") from error
    if result != value:
        raise ValidationError(f"{label} must be an integer")
    if minimum is not None and result < minimum:
        raise ValidationError(f"{label} must be >= {minimum}")
    return result


def load_solver_mesh(path: Path) -> dict:
    data = json.loads(path.read_text())
    if data.get("format") == SOLVER_MESH_FORMAT:
        return data
    nested = data.get("solver_mesh")
    if isinstance(nested, dict) and nested.get("format") == SOLVER_MESH_FORMAT:
        return nested
    raise ValidationError(
        f"expected {SOLVER_MESH_FORMAT} or a backend result containing solver_mesh"
    )


def validate_solver_mesh(data: dict) -> dict:
    if data.get("format") != SOLVER_MESH_FORMAT:
        raise ValidationError(f"unsupported solver mesh format: {data.get('format')!r}")
    if data.get("index_base") != 0:
        raise ValidationError("solver mesh must use zero-based connectivity")

    vertices = data.get("vertices")
    if not isinstance(vertices, list) or not vertices:
        raise ValidationError("solver mesh has no vertices")
    normalized_vertices = [
        _finite_triplet(vertex, f"vertices[{index}]")
        for index, vertex in enumerate(vertices)
    ]

    def normalize_regions(name: str, required: tuple[str, ...]) -> tuple[list[dict], dict[int, dict]]:
        regions = data.get(name)
        if not isinstance(regions, list):
            raise ValidationError(f"{name} must be an array")
        by_id: dict[int, dict] = {}
        normalized = []
        for index, region in enumerate(regions):
            if not isinstance(region, dict):
                raise ValidationError(f"{name}[{index}] must be an object")
            region_id = _integer(region.get("region"), f"{name}[{index}].region", 0)
            if region_id in by_id:
                raise ValidationError(f"duplicate {name} region id {region_id}")
            item = dict(region)
            item["region"] = region_id
            item["stiffness"] = _integer(
                region.get("stiffness"), f"{name}[{index}].stiffness", 0
            )
            for field in required:
                if field not in region:
                    raise ValidationError(f"{name}[{index}] missing {field}")
            source_elements = region.get("source_elements", [])
            if not isinstance(source_elements, list) or not source_elements:
                raise ValidationError(f"{name}[{index}] has no source_elements")
            item["source_elements"] = sorted(
                {
                    _integer(
                        value,
                        f"{name}[{index}].source_elements",
                        1,
                    )
                    for value in source_elements
                }
            )
            by_id[region_id] = item
            normalized.append(item)
        return normalized, by_id

    surface_regions, surface_by_id = normalize_regions(
        "surface_regions", ("source_patch",)
    )
    bar_regions, bar_by_id = normalize_regions(
        "bar_regions", ("axis", "source_axis", "start_t", "end_t")
    )

    def normalize_elements(
        name: str,
        expected_nodes: int,
        regions: dict[int, dict],
    ) -> list[dict]:
        elements = data.get(name)
        if not isinstance(elements, list):
            raise ValidationError(f"{name} must be an array")
        result = []
        for index, element in enumerate(elements):
            if not isinstance(element, dict):
                raise ValidationError(f"{name}[{index}] must be an object")
            nodes = element.get("vertices")
            if not isinstance(nodes, list) or len(nodes) != expected_nodes:
                raise ValidationError(
                    f"{name}[{index}] must have {expected_nodes} vertex indices"
                )
            normalized_nodes = [
                _integer(node, f"{name}[{index}].vertices", 0) for node in nodes
            ]
            if len(set(normalized_nodes)) != expected_nodes:
                raise ValidationError(f"{name}[{index}] is topologically degenerate")
            if any(node >= len(normalized_vertices) for node in normalized_nodes):
                raise ValidationError(f"{name}[{index}] references a missing vertex")
            region = _integer(element.get("region"), f"{name}[{index}].region", 0)
            if region not in regions:
                raise ValidationError(f"{name}[{index}] references missing region {region}")
            result.append({"vertices": normalized_nodes, "region": region})
        return result

    shells = normalize_elements("shell_elements", 3, surface_by_id)
    bars = normalize_elements("bar_elements", 2, bar_by_id)
    if not shells and not bars:
        raise ValidationError("solver mesh has no elements")

    represented_surfaces = {element["region"] for element in shells}
    represented_bars = {element["region"] for element in bars}
    missing_surfaces = sorted(set(surface_by_id) - represented_surfaces)
    missing_bars = sorted(set(bar_by_id) - represented_bars)
    if missing_surfaces or missing_bars:
        raise ValidationError(
            f"unrepresented regions: surfaces={missing_surfaces}, bars={missing_bars}"
        )

    used_vertices = {
        node
        for element in shells + bars
        for node in element["vertices"]
    }
    return {
        "length_unit": data.get("length_unit", "model_unit"),
        "vertices": normalized_vertices,
        "surface_regions": surface_regions,
        "surface_by_id": surface_by_id,
        "bar_regions": bar_regions,
        "bar_by_id": bar_by_id,
        "shell_elements": shells,
        "bar_elements": bars,
        "unused_vertex_indices": sorted(set(range(len(vertices))) - used_vertices),
    }


def _safe_group_name(prefix: str, stiffness: int) -> str:
    raw = f"{prefix}_K{stiffness}"
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", raw)


def build_import_plan(data: dict) -> dict:
    model = validate_solver_mesh(data)
    shell_count = len(model["shell_elements"])

    shell_records = []
    for index, element in enumerate(model["shell_elements"]):
        region = model["surface_by_id"][element["region"]]
        shell_records.append(
            {
                "element_id": index + 1,
                "kind": "TRIA",
                "node_ids": [node + 1 for node in element["vertices"]],
                "region": element["region"],
                "stiffness": region["stiffness"],
            }
        )

    bar_records = []
    for index, element in enumerate(model["bar_elements"]):
        region = model["bar_by_id"][element["region"]]
        bar_records.append(
            {
                "element_id": shell_count + index + 1,
                "kind": "LINE",
                "node_ids": [node + 1 for node in element["vertices"]],
                "region": element["region"],
                "stiffness": region["stiffness"],
            }
        )

    grouped: dict[tuple[str, int], dict] = {}
    for record in shell_records + bar_records:
        dimension = "SHELL" if record["kind"] == "TRIA" else "BAR"
        key = (dimension, record["stiffness"])
        group = grouped.setdefault(
            key,
            {
                "name": _safe_group_name(dimension, record["stiffness"]),
                "dimension": 2 if dimension == "SHELL" else 1,
                "stiffness": record["stiffness"],
                "element_ids": [],
                "node_ids": set(),
                "regions": set(),
            },
        )
        group["element_ids"].append(record["element_id"])
        group["node_ids"].update(record["node_ids"])
        group["regions"].add(record["region"])

    mesh_sets = []
    for key in sorted(grouped):
        group = grouped[key]
        mesh_sets.append(
            {
                "mesh_set_id": len(mesh_sets) + 1,
                "name": group["name"],
                "dimension": group["dimension"],
                "stiffness": group["stiffness"],
                "element_ids": group["element_ids"],
                "node_ids": sorted(group["node_ids"]),
                "regions": sorted(group["regions"]),
            }
        )

    property_registry = []
    for dimension, stiffness in sorted(
        {(item["dimension"], item["stiffness"]) for item in mesh_sets}
    ):
        property_registry.append(
            {
                "property_slot": len(property_registry) + 1,
                "dimension": dimension,
                "source_stiffness": stiffness,
                "assignment": "unresolved_until_fpn_profile_calibration",
            }
        )

    return {
        "format": PLAN_FORMAT,
        "target": "MIDAS GTS NX",
        "target_format": "FPN",
        "format_profile": {
            **OBSERVED_GTS_NX,
            "status": "serialization_not_yet_calibrated",
        },
        "length_unit": model["length_unit"],
        "id_policy": {
            "node_ids": "solver vertex index + 1",
            "shell_element_ids": [1, shell_count],
            "bar_element_ids": [
                shell_count + 1,
                shell_count + len(bar_records),
            ],
            "global_element_ids_unique": True,
        },
        "counts": {
            "nodes": len(model["vertices"]),
            "shell_elements": len(shell_records),
            "bar_elements": len(bar_records),
            "surface_regions": len(model["surface_regions"]),
            "bar_regions": len(model["bar_regions"]),
            "mesh_sets": len(mesh_sets),
            "unused_vertices": len(model["unused_vertex_indices"]),
        },
        "nodes": [
            {"node_id": index + 1, "xyz": list(point)}
            for index, point in enumerate(model["vertices"])
        ],
        "shell_elements": shell_records,
        "bar_elements": bar_records,
        "surface_regions": model["surface_regions"],
        "bar_regions": model["bar_regions"],
        "mesh_sets": mesh_sets,
        "property_registry": property_registry,
        "calibration_required": [
            "target GTS NX FPN version/header records",
            "exact modern TRIA field layout and property field semantics",
            "exact modern LINE field layout and property/orientation field semantics",
            "exact MSET/MSETE/MSETN record syntax and continuation formatting",
            "text encoding emitted by the target GTS NX version",
        ],
    }


def decode_fpn(raw: bytes) -> tuple[str, str]:
    for encoding in ("utf-8-sig", "gbk", "cp1252"):
        try:
            return raw.decode(encoding), encoding
        except UnicodeDecodeError:
            pass
    raise ValidationError("could not decode FPN as utf-8-sig, GBK or cp1252")


def inspect_fpn_text(text: str, encoding: str = "unknown") -> dict:
    section = None
    sections = []
    records: dict[str, dict] = {}
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith("$$"):
            section = line[2:].strip()
            sections.append(section)
            continue
        token = line.split(",", 1)[0].strip().upper()
        if token not in {"NODE", "LINE", "TRIA", "RECT", "MSET", "MSETE", "MSETN"}:
            continue
        fields = [part.strip() for part in line.rstrip(",").split(",")]
        item = records.setdefault(
            token,
            {
                "count": 0,
                "minimum_field_count": len(fields),
                "maximum_field_count": len(fields),
                "sections": set(),
                "examples": [],
            },
        )
        item["count"] += 1
        item["minimum_field_count"] = min(item["minimum_field_count"], len(fields))
        item["maximum_field_count"] = max(item["maximum_field_count"], len(fields))
        if section:
            item["sections"].add(section)
        if len(item["examples"]) < 3:
            item["examples"].append(fields)

    for item in records.values():
        item["sections"] = sorted(item["sections"])

    return {
        "format": INSPECTION_FORMAT,
        "encoding": encoding,
        "sections": sections,
        "records": records,
        "observed_profile": OBSERVED_GTS_NX,
        "calibration_ready": all(
            token in records for token in ("NODE", "TRIA", "LINE", "MSET", "MSETE", "MSETN")
        ),
    }


def inspect_fpn(path: Path) -> dict:
    text, encoding = decode_fpn(path.read_bytes())
    return inspect_fpn_text(text, encoding)


def self_test() -> None:
    mesh = {
        "format": SOLVER_MESH_FORMAT,
        "length_unit": "m",
        "index_base": 0,
        "vertices": [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ],
        "surface_regions": [
            {
                "region": 0,
                "source_patch": 3,
                "stiffness": 10,
                "source_elements": [11, 12],
            }
        ],
        "bar_regions": [
            {
                "region": 7,
                "axis": 0,
                "source_axis": 5,
                "stiffness": 20,
                "source_elements": [13],
                "start_t": 0.0,
                "end_t": 1.0,
            }
        ],
        "shell_elements": [{"vertices": [0, 1, 2], "region": 0}],
        "bar_elements": [{"vertices": [0, 3], "region": 7}],
    }
    plan = build_import_plan(mesh)
    assert plan["counts"]["nodes"] == 4
    assert plan["shell_elements"][0]["element_id"] == 1
    assert plan["bar_elements"][0]["element_id"] == 2
    assert plan["shell_elements"][0]["node_ids"] == [1, 2, 3]
    assert {group["name"] for group in plan["mesh_sets"]} == {
        "SHELL_K10",
        "BAR_K20",
    }
    assert plan["format_profile"]["writer_verified"] is False

    sample = """$$      Node
NODE,1,0.,0.,0.
NODE,2,1.,0.,0.
$$      Element
TRIA,1,10,1,2,3
LINE,2,20,1,2,0.,0.
$$      Mesh Set
MSET,1,SHELL_K10
MSETE,1,1
1,
MSETN,1,3
1,2,3,
"""
    inspection = inspect_fpn_text(sample, "synthetic")
    assert inspection["records"]["NODE"]["minimum_field_count"] == 5
    assert inspection["records"]["TRIA"]["maximum_field_count"] == 6
    assert inspection["records"]["LINE"]["maximum_field_count"] == 7
    assert inspection["records"]["MSET"]["count"] == 1


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", nargs="?", type=Path, help="solver mesh or backend JSON")
    parser.add_argument("--plan-output", type=Path)
    parser.add_argument("--inspect-fpn", type=Path)
    parser.add_argument("--inspection-output", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()

    if args.self_test:
        self_test()
        print("MIDAS FPN adapter self-test: OK")
        return 0

    if args.inspect_fpn:
        result = inspect_fpn(args.inspect_fpn)
        text = json.dumps(result, indent=2, ensure_ascii=False)
        if args.inspection_output:
            args.inspection_output.write_text(text)
        print(text)
        return 0

    if args.input is None:
        parser.error("input is required unless --inspect-fpn or --self-test is used")

    solver_mesh = load_solver_mesh(args.input)
    plan = build_import_plan(solver_mesh)
    text = json.dumps(plan, indent=2, ensure_ascii=False)
    output = args.plan_output
    if output:
        output.write_text(text)
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
