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
