"""Independent, stdlib-only checks of recognize_frame JSON (not a mesh gate)."""
import argparse
import collections
import json
import math


def _id(value, size):
    assert isinstance(value, int) and not isinstance(value, bool)
    assert 0 <= value < size


def _finite_vector(value, length):
    assert isinstance(value, list) and len(value) == length
    assert all(isinstance(x, (int, float)) and math.isfinite(x) for x in value)


def _sorted_unique_ids(value):
    assert isinstance(value, list)
    assert all(isinstance(x, int) and not isinstance(x, bool) for x in value)
    assert value == sorted(set(value))


def check_mesh(mesh, model, bars):
    """Check optional diagnostic mesh data without treating it as export proof."""
    assert isinstance(mesh, dict)
    vertices = mesh["vertices"]
    surface_source_elements = mesh["surface_source_elements"]
    surface_count = len(model["surfaces"])
    axis_count = len(bars["axes"])
    assert vertices
    for vertex in vertices:
        _finite_vector(vertex, 3)

    expected_surface_sources = [s["source_elements"] for s in model["surfaces"]]
    assert surface_source_elements == expected_surface_sources
    assert len(surface_source_elements) == surface_count

    for triangle in mesh["triangles"]:
        ids = triangle["vertices"]
        assert isinstance(ids, list) and len(ids) == 3
        assert len(set(ids)) == 3
        for vertex in ids:
            _id(vertex, len(vertices))
        _id(triangle["surface"], surface_count)
        assert isinstance(triangle["stiffness"], int)

    for bar in mesh["bars"]:
        ids = bar["vertices"]
        assert isinstance(ids, list) and len(ids) == 2
        assert ids[0] != ids[1]
        for vertex in ids:
            _id(vertex, len(vertices))
        _id(bar["axis"], axis_count)
        source = (bar["source_element"], bar["stiffness"])
        spans = {
            (span["element"], span["stiffness"])
            for span in bars["axes"][bar["axis"]]["spans"]
        }
        assert source in spans

    triangle_surfaces = {triangle["surface"] for triangle in mesh["triangles"]}
    diagnostic_surfaces = set()
    for diagnostic in mesh["mesh_surface_errors"]:
        _id(diagnostic["surface"], surface_count)
        assert diagnostic["surface"] not in diagnostic_surfaces
        diagnostic_surfaces.add(diagnostic["surface"])
        assert isinstance(diagnostic["reason"], str) and diagnostic["reason"]
        assert set(diagnostic["source_elements"]) <= set(
            surface_source_elements[diagnostic["surface"]]
        )
        assert diagnostic["surface"] not in triangle_surfaces

    for diagnostic in mesh["quality_diagnostics"]:
        _id(diagnostic["surface"], surface_count)
        assert set(diagnostic["source_elements"]) <= set(
            surface_source_elements[diagnostic["surface"]]
        )
        assert isinstance(diagnostic["vertices"], list)
        assert len(diagnostic["vertices"]) == 3
        for vertex in diagnostic["vertices"]:
            _id(vertex, len(vertices))
        assert len(diagnostic["source_nodes"]) == 3
        for node in diagnostic["source_nodes"]:
            if node is not None:
                assert isinstance(node, int) and not isinstance(node, bool)
        for edge in diagnostic["constrained_edges"]:
            assert isinstance(edge, list) and len(edge) == 2
            for vertex in edge:
                _id(vertex, len(vertices))
        assert diagnostic["reasons"]
        assert all(isinstance(reason, str) for reason in diagnostic["reasons"])
        for field in ("area", "minimum_angle_degrees", "edge_ratio"):
            assert math.isfinite(diagnostic[field])

    unresolved_surface = mesh["unresolved_surface_source_elements"]
    unresolved_axis = mesh["unresolved_axis_source_elements"]
    _sorted_unique_ids(unresolved_surface)
    _sorted_unique_ids(unresolved_axis)
    for diagnostic in mesh["mesh_surface_errors"]:
        assert set(diagnostic["source_elements"]) <= set(unresolved_surface)

    assert isinstance(mesh["source_coverage_complete"], bool)
    assert isinstance(mesh["topology_valid"], bool)
    assert isinstance(mesh["quality_passed"], bool)
    assert isinstance(mesh["export_ready"], bool)
    assert not mesh["export_ready"]
    assert all(isinstance(blocker, str) for blocker in mesh["blockers"])
    assert all(
        isinstance(blocker, str) for blocker in mesh["external_mesher_blockers"]
    )
    assert mesh["external_mesher_ready"] == (not mesh["external_mesher_blockers"])
    if mesh["source_coverage_complete"]:
        assert not unresolved_surface and not unresolved_axis
        assert not mesh["mesh_surface_errors"]
    else:
        assert not mesh["external_mesher_ready"]
    if mesh["mesh_surface_errors"]:
        assert not mesh["source_coverage_complete"]
        assert not mesh["topology_valid"]

    synchronization = mesh["constraint_synchronization"]
    bindings = synchronization["endpoint_bindings"]
    assert synchronization["synchronized_interval_endpoints"] == len(bindings)
    assert synchronization["axis_node_count"] <= len(vertices)
    assert synchronization["edge_node_count"] <= len(vertices)
    assert synchronization["interval_contact_count"] <= len(bars["contacts"])
    for binding in bindings:
        _id(binding["axis"], axis_count)
        _id(binding["surface"], surface_count)
        _id(binding["contact"], len(bars["contacts"]))
        _id(binding["vertex"], len(vertices))
        assert binding["role"] in ("start", "end")
        assert math.isfinite(binding["parameter"])
        assert 0.0 <= binding["parameter"] <= 1.0

    for field in ("minimum_angle_degrees", "maximum_triangle_area", "maximum_edge_ratio"):
        assert math.isfinite(mesh[field])
    return {
        "vertices": len(vertices),
        "triangles": len(mesh["triangles"]),
        "bars": len(mesh["bars"]),
        "source_coverage_complete": mesh["source_coverage_complete"],
        "external_mesher_ready": mesh["external_mesher_ready"],
    }


def check(data, baseline=None):
    frame, topology = data["frame"], data["topology"]
    model, bars = topology["preview"], topology["axis_assembly"]
    epsilon = topology["policy"]["precision"]
    vertices = model["vertices"]
    source_nodes = topology["vertex_source_nodes"]
    assert len(vertices) == len(source_nodes) == len(set(source_nodes))
    lookup = {n: i for i, n in enumerate(frame["node_ids"])}
    budget = {n: frame["policy"]["maximum_movement"] for n in source_nodes}
    for axis in frame["axes"]:
        limit = frame["policy"]["relative_movement"] * math.dist(
            *(frame["reference_points"][i] for i in axis["endpoints"])
        )
        for anchor in axis["anchors"]:
            n = frame["node_ids"][anchor["node"]]
            if n in budget:
                budget[n] = min(budget[n], limit)
    for n, point in zip(source_nodes, vertices):
        i = lookup[n]
        assert math.dist(point, frame["reference_points"][i]) <= budget[n] + epsilon
        assert math.dist(point, frame["candidate_points"][i]) <= topology["policy"]["junction_movement_limit"] + epsilon

    expected_shells = [e for s in frame["surfaces"] for e in s["source_elements"]]
    actual_shells = [e for s in model["surfaces"] for e in s["source_elements"]]
    actual_shells += [e for issue in topology["issues"] for e in issue["source_elements"]]
    assert collections.Counter(expected_shells) == collections.Counter(actual_shells)
    assert len(actual_shells) == len(set(actual_shells))
    assert len(model["edges"]) == len(set(tuple(e) for e in model["edges"]))

    def distance(surface, point):
        plane = model["planes"][surface["plane"]]
        return sum((a - b) * n for a, b, n in zip(point, plane["origin"], plane["normal"]))

    for surface, patch, stiffness in zip(model["surfaces"], topology["surface_source_patches"], topology["surface_stiffness"]):
        assert set(surface["source_elements"]) <= set(frame["surfaces"][patch]["stiffness_regions"][str(stiffness)])
        for ring, contour in zip(surface["boundaries"], surface["contours"]):
            assert len(ring) == len(contour)
            plane = model["planes"][surface["plane"]]
            for edge, uv in zip(ring, contour):
                start = model["edges"][edge["edge"]][1 if edge["reversed"] else 0]
                lifted = [plane["origin"][k] + uv[0] * plane["u"][k] + uv[1] * plane["v"][k] for k in range(3)]
                assert math.dist(lifted, vertices[start]) <= epsilon
                for v in model["edges"][edge["edge"]]:
                    assert abs(distance(surface, vertices[v])) <= epsilon

    expected_bars = [s["element"] for a in frame["axes"] for s in a["spans"]]
    actual_bars = [s["element"] for a in bars["axes"] for s in a["spans"]]
    actual_bars += [e for i in bars["issues"] for e in i["source_elements"]]
    assert collections.Counter(expected_bars) == collections.Counter(actual_bars)
    assert len(actual_bars) == len(set(actual_bars))
    represented = [a["source_axis"] for a in bars["axes"]] + [i["source_axis"] for i in bars["issues"]]
    assert sorted(represented) == list(range(len(frame["axes"])))
    for axis in bars["axes"]:
        original = frame["axes"][axis["source_axis"]]
        assert len(axis["endpoints"]) == 2
        a, b = (vertices[v] for v in axis["endpoints"])
        assert math.dist(a, b) > epsilon
        assert [source_nodes[v] for v in axis["endpoints"]] == [frame["node_ids"][i] for i in original["endpoints"]]
        parameters = {}
        for anchor in axis["anchors"]:
            assert source_nodes[anchor["vertex"]] == anchor["source_node"]
            assert 0 <= anchor["t"] <= 1
            interpolated = [x + (y - x) * anchor["t"] for x, y in zip(a, b)]
            assert math.dist(interpolated, vertices[anchor["vertex"]]) <= epsilon
            parameters[anchor["source_node"]] = anchor["t"]
        original_parameters = {x["t"]: frame["node_ids"][x["node"]] for x in original["anchors"]}
        assert len(axis["spans"]) == len(original["spans"])
        for span, source in zip(axis["spans"], original["spans"]):
            assert (span["element"], span["stiffness"]) == (source["element"], source["stiffness"])
            assert span["start_t"] == parameters[original_parameters[source["start_t"]]]
            assert span["end_t"] == parameters[original_parameters[source["end_t"]]]
            assert 0 <= span["start_t"] < span["end_t"] <= 1

    for contact in bars["contacts"]:
        axis, surface = bars["axes"][contact["axis"]], model["surfaces"][contact["surface"]]
        if contact["kind"] == "point":
            assert any(x["vertex"] == contact["vertex"] and x["t"] == contact["t"] for x in axis["anchors"])
            assert abs(distance(surface, vertices[contact["vertex"]])) <= epsilon
        else:
            assert 0 <= contact["start_t"] < contact["end_t"] <= 1
            assert all(abs(distance(surface, vertices[v])) <= epsilon for v in axis["endpoints"])
    if baseline:
        before = baseline["topology"]["preview"]
        old_axes = baseline["topology"].get("axis_assembly", {}).get("axes", [])
        assert {a["source_axis"] for a in old_axes} <= {a["source_axis"] for a in bars["axes"]}
        for field in ("edges", "planes"):
            assert model[field] == before[field], field
        assert len(model["surfaces"]) == len(before["surfaces"])
        for old, new in zip(before["surfaces"], model["surfaces"]):
            for field in ("plane", "boundaries", "source_elements"):
                assert old[field] == new[field]
        logged = {}
        for repair in bars.get("boundary_repairs", []):
            if not repair["accepted"]:
                assert not repair["changes"]
            for change in repair["changes"]:
                assert change["source_node"] not in logged
                logged[change["source_node"]] = change
        current_by_node = dict(zip(source_nodes, vertices))
        boundary_vertices = {v for edge in before["edges"] for v in edge}
        for v in boundary_vertices:
            n = baseline["topology"]["vertex_source_nodes"][v]
            old = before["vertices"][v]
            if n in logged:
                assert logged[n]["before"] == old
                assert logged[n]["after"] == current_by_node[n]
            else:
                assert current_by_node[n] == old
    assert not topology["export_ready"]
    assert not bars["mesh_constraints_complete"]
    result = {"surfaces": len(model["surfaces"]), "axes": len(bars["axes"]),
              "rejected_axes": len(bars["issues"]), "shell_elements": len(actual_shells),
              "bar_elements": len(actual_bars), "contacts": len(bars["contacts"])}
    if "mesh" in data:
        mesh = data["mesh"]
        if mesh is None:
            assert isinstance(data.get("mesh_error"), str) and data["mesh_error"]
        else:
            assert data.get("mesh_error") is None
            result["mesh"] = check_mesh(mesh, model, bars)
            # Removed material boundaries must not erase shared source joints.
            for hole in topology.get("simplified_holes", []):
                surfaces = [i for i, surface in enumerate(model["surfaces"])
                            if surface["source_elements"] == hole["source_elements"]]
                assert len(surfaces) == 1
                used = {v for tri in mesh["triangles"] if tri["surface"] == surfaces[0]
                        for v in tri["vertices"]}
                if mesh["topology_valid"]:
                    assert all(source_nodes.index(n) in used for n in hole["source_nodes"])

    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("report")
    parser.add_argument("--baseline", help="Check retained topology, logged boundary moves, and previously accepted axes")
    args = parser.parse_args()
    with open(args.report, encoding="utf-8") as stream:
        data = json.load(stream)
    baseline = None
    if args.baseline:
        with open(args.baseline, encoding="utf-8") as stream:
            baseline = json.load(stream)
    print(json.dumps(check(data, baseline), ensure_ascii=False, indent=2))

