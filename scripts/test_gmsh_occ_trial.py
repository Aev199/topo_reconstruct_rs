import contextlib
import importlib
import io
import json
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock


try:
    import gmsh  # noqa: F401
except BaseException:
    sys.modules["gmsh"] = types.ModuleType("gmsh")

sys.path.insert(0, str(Path(__file__).parent))
gmsh_occ_trial = importlib.import_module("gmsh_occ_trial")


class GmshOccTrialTests(unittest.TestCase):
    def test_occ_fragments_tool_tool_wall_slab_and_hole(self):
        import gmsh

        def face(rings, region):
            return {"rings": rings, "source_surface": region, "source_patch": region,
                    "stiffness": region + 1, "source_elements": [region + 1]}

        unrelated = face([[[10, 0, 0], [11, 0, 0], [11, 1, 0], [10, 1, 0]]], 0)
        slab = face([
            [[0, 0, 0], [4, 0, 0], [4, 4, 0], [0, 4, 0]],
            [[1, 1, 0], [1, 2, 0], [2, 2, 0], [2, 1, 0]],
        ], 1)
        wall = face([[[3, 0, -1], [3, 4, -1], [3, 4, 1], [3, 0, 1]]], 2)
        gmsh.initialize()
        try:
            gmsh.option.setNumber("General.Terminal", 0)
            for shift, reverse in ((0.0, False), (100000.0, True)):
                gmsh.clear()
                gmsh.model.add("tool_tool_junction")
                items = json.loads(json.dumps([unrelated, slab, wall]))
                for item in items:
                    for ring in item["rings"]:
                        if reverse:
                            ring.reverse()
                        for point in ring:
                            point[0] += shift
                            point[1] -= shift
                fragmented = gmsh_occ_trial.fragment(items, [], [], 1e-8)
                slab_edges = {tag for face_tag in fragmented.surface_to_output[1]
                              for dim, tag in gmsh.model.getBoundary([(2, face_tag)], oriented=False)
                              if dim == 1}
                wall_edges = {tag for face_tag in fragmented.surface_to_output[2]
                              for dim, tag in gmsh.model.getBoundary([(2, face_tag)], oriented=False)
                              if dim == 1}
                self.assertTrue(slab_edges & wall_edges)
                slab_area = sum(gmsh.model.occ.getMass(2, tag)
                                for tag in fragmented.surface_to_output[1])
                self.assertAlmostEqual(slab_area, 15.0, places=7)
        finally:
            gmsh.finalize()

    def test_near_vertex_inversion_is_transactionally_rejected(self):
        coords = {
            1: (0.0, 0.0, 0.0), 2: (0.02, 0.0, 0.0), 3: (0.02, 1.0, 0.0),
            4: (0.2, 0.0, 0.0), 5: (-1.0, -1.0, 0.0), 6: (0.01, -1.0, 0.0),
            7: (0.014, -0.5, 0.0), 8: (0.02, 1.0, 1.0), 9: (0.02, 0.0, 1.0),
            10: (0.02, -1.0, 1.0), 11: (0.02, -1.0, 0.0),
        }
        rings = [[1, 4, 3, 5, 6, 7], [11, 3, 8, 10]]
        data = {
            "surfaces": [{
                "source_surface": surface, "rings": [[coords[node] for node in ring]],
                "ring_source_nodes": [[surface * 100 + node for node in ring]],
            } for surface, ring in enumerate(rings)],
            "axes": [], "contacts": [],
        }
        triangles = [
            {"nodes": nodes, "source_surface": surface, "output_surface": surface + 1,
             "stiffness": surface + 10}
            for surface, nodes in [
                (0, (1, 2, 3)), (0, (2, 4, 3)), (0, (1, 3, 5)),
                (0, (1, 5, 6)), (0, (1, 6, 7)), (1, (2, 3, 8)),
                (1, (2, 8, 9)), (1, (2, 9, 10)), (1, (2, 10, 11)),
            ]
        ]
        repaired, repaired_triangles, _bars, report = gmsh_occ_trial.regularize_near_vertex_junctions(
            data, coords, triangles, [], 0.05, 0.001, 1e-8
        )
        self.assertEqual(repaired[1], coords[1])
        self.assertEqual(report["accepted_count"], 0)
        self.assertTrue(any(item.get("reason") == "inverted_surviving_face" for item in report["skipped"]))
        self.assertEqual(len(repaired_triangles), len(triangles))

    def test_safe_micro_edge_collapse_still_applies(self):
        coords = {
            1: (0.0, 0.0, 0.0), 2: (0.00004, 0.0, 0.0), 3: (0.0, 0.00008, 0.0),
            4: (1.0, 0.0, 0.0), 5: (0.0, 1.0, 0.0), 6: (0.0, 0.0, 1.0),
            7: (0.0, 1.0, 1.0),
        }
        triangles = [
            {"nodes": (1, 2, 4), "source_surface": 0},
            {"nodes": (1, 4, 5), "source_surface": 0},
            {"nodes": (1, 3, 6), "source_surface": 1},
            {"nodes": (1, 6, 7), "source_surface": 1},
        ]
        healed, mapping, report = gmsh_occ_trial.heal_micro_edges(coords, triangles, 0.001, 1e-8)
        self.assertEqual(mapping, {2: 1, 3: 1})
        self.assertEqual(report["collapsed_component_count"], 1)
        self.assertEqual(report["unresolved_component_count"], 0)
        self.assertEqual(len(healed), 2)

    def test_subresolution_hole_edge_may_close_with_logged_bounded_move(self):
        # A plate has a 0.5 mm triangular hole edge.  Node 5 is also the
        # declared wall junction, so the existing junction-local policy may
        # absorb node 6 and simplify this micro-opening.
        coords = {
            1: (-1.0, -1.0, 0.0), 2: (1.0, -1.0, 0.0), 3: (1.0, 2.0, 0.0),
            4: (-1.0, 2.0, 0.0), 5: (0.0, 0.0, 0.0), 6: (0.0005, 0.0, 0.0),
            7: (0.0, 1.0, 0.0), 8: (1.0, 0.0, 0.0), 9: (1.0, 1.0, 0.0),
            10: (0.2, 1.5, 0.0), 11: (-1.0, 1.5, 0.0),
            12: (0.0, 0.0, 1.0), 13: (0.0, 0.2, 2.0),
        }
        triangles = [
            {"nodes": nodes, "source_surface": 0} for nodes in [
                (5, 6, 8), (5, 8, 9), (5, 9, 7), (5, 7, 10),
                (5, 10, 11), (5, 11, 4),
            ]
        ] + [{"nodes": (5, 12, 13), "source_surface": 1}]
        healed, mapping, report = gmsh_occ_trial.heal_micro_edges(coords, triangles, 0.001, 1e-8)
        self.assertEqual(mapping, {6: 5})
        self.assertEqual(report["collapsed_component_count"], 1)
        self.assertEqual(report["removed_degenerate_triangles"], 1)
        move = report["components"][0]["moves"][0]
        self.assertEqual((move["from"], move["to"]), (6, 5))
        self.assertLess(move["movement"], 0.001)
        self.assertTrue(all(
            gmsh_occ_trial.triangle_area(*(coords[node] for node in triangle["nodes"])) > 1e-16
            for triangle in healed
        ))

    def test_overlap_guard_is_stable_for_translated_scaled_coordinates(self):
        for scale, shift in ((1.0, 0.0), (1e-3, 1e9)):
            left = [(shift + scale * x, shift + scale * y, 0.0) for x, y in ((0, 0), (2, 0), (0, 2))]
            right = [(shift + scale * x, shift + scale * y, 0.0) for x, y in ((0.5, 0.5), (2.5, 0.5), (0.5, 2.5))]
            self.assertTrue(gmsh_occ_trial._coplanar_triangle_overlap(left, right, max(1e-12, scale * 1e-8)))

    def test_incomplete_coverage_is_a_blocker_for_both_formats(self):
        versioned = {
            "format": gmsh_occ_trial.INPUT_FORMAT,
            "policy": {"precision": 1e-8, "minimum_edge": 1e-3, "target_mesh_size": 1.0},
            "surfaces": [], "axes": [], "contacts": [],
            "source_coverage_complete": False,
        }
        normalized = gmsh_occ_trial.normalize_input(versioned, None)
        self.assertIn("incomplete_source_coverage", normalized["blockers"])

        legacy = {
            "topology": {
                "preview": {"surfaces": [], "planes": [], "edges": [], "vertices": []},
                "vertex_source_nodes": [],
                "surface_stiffness": [], "surface_source_patches": [],
                "policy": {},
            }
        }
        normalized = gmsh_occ_trial.normalize_input(legacy, None)
        self.assertFalse(normalized["source_coverage_complete"])
        self.assertIn("incomplete_source_coverage", normalized["blockers"])

    def test_msh_writer_uses_final_solver_vertices_and_provenance_groups(self):
        package = {
            "vertices": [[0, 0, 0], [0, 0, 0], [1, 0, 0], [0, 1, 0]],
            "surface_regions": [{"region": 4, "stiffness": 12}],
            "bar_regions": [{"region": 9, "stiffness": 30}],
            "shell_elements": [{"vertices": [0, 2, 3], "region": 4}],
            "bar_elements": [{"vertices": [0, 1], "region": 9}],
        }
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "final.msh"
            gmsh_occ_trial.write_solver_mesh_msh(path, package)
            text = path.read_text()
        self.assertIn("$PhysicalNames", text)
        self.assertIn('2 1 "surface_region_4_k12"', text)
        self.assertIn('1 2 "bar_region_9_k30"', text)
        self.assertIn("$Nodes\n4\n", text)
        self.assertIn("$Elements\n2\n", text)

    def test_cli_writes_diagnostic_but_blocks_solver_output(self):
        result = {
            "backend_ready": False,
            "blockers": ["incomplete_source_coverage"],
            "solver_mesh": {"format": gmsh_occ_trial.SOLVER_MESH_FORMAT},
            "solver_mesh_audit": {"clean": True},
        }
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "diagnostic.json"
            solver = Path(directory) / "solver.json"
            input_path = Path(directory) / "input.json"
            input_path.write_text("{}")
            argv = ["gmsh_occ_trial.py", str(input_path), "--output", str(output), "--solver-mesh-output", str(solver)]
            with mock.patch.object(sys, "argv", argv), mock.patch.object(gmsh_occ_trial.gmsh, "initialize", create=True), mock.patch.object(gmsh_occ_trial.gmsh, "finalize", create=True), mock.patch.object(gmsh_occ_trial, "run_backend", return_value=result), contextlib.redirect_stdout(io.StringIO()):
                status = gmsh_occ_trial.main()
            self.assertEqual(status, 2)
            self.assertTrue(output.exists())
            self.assertFalse(solver.exists())
            self.assertEqual(json.loads(output.read_text())["blockers"], ["incomplete_source_coverage"])


if __name__ == "__main__":
    unittest.main()
