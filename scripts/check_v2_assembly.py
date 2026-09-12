"""Independent, stdlib-only checks of recognize_frame JSON (not a mesh gate)."""
import argparse
import collections
import json
import math


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
        for ring in surface["boundaries"]:
            for edge in ring:
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
        for field in ("surfaces", "edges", "planes"):
            assert model[field] == before[field], field
        assert vertices[:len(before["vertices"])] == before["vertices"]
    assert not topology["export_ready"]
    assert not bars["mesh_constraints_complete"]
    return {"surfaces": len(model["surfaces"]), "axes": len(bars["axes"]),
            "rejected_axes": len(bars["issues"]), "shell_elements": len(actual_shells),
            "bar_elements": len(actual_bars), "contacts": len(bars["contacts"])}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("report")
    parser.add_argument("--baseline", help="Optional surface-only report; check it remains unchanged")
    args = parser.parse_args()
    with open(args.report, encoding="utf-8") as stream:
        data = json.load(stream)
    baseline = None
    if args.baseline:
        with open(args.baseline, encoding="utf-8") as stream:
            baseline = json.load(stream)
    print(json.dumps(check(data, baseline), ensure_ascii=False, indent=2))
