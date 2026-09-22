import unittest

from check_gmsh_surface_junctions import check, repaired_segment_covered


class FiniteJunctionTests(unittest.TestCase):
    @staticmethod
    def model(shared: bool):
        surfaces = [
            {"rings": [[[0, 0, 0], [2, 0, 0], [2, 2, 0], [0, 2, 0]]],
             "source_elements": [1]},
            {"rings": [[[1, 0, -1], [1, 2, -1], [1, 2, 1], [1, 0, 1]]],
             "source_elements": [2]},
        ]
        vertices = [[0, 0, 0], [1, 0, 0], [1, 2, 0], [0, 2, 0], [2, 0, 0], [2, 2, 0],
                    [1, 0, -1], [1, 2, -1], [1, 2, 1], [1, 0, 1]]
        if not shared:
            vertices += [[1, 0, 0], [1, 2, 0]]
        tris = [([0, 1, 2], 0), ([0, 2, 3], 0), ([1, 4, 5], 0), ([1, 5, 2], 0),
                ([6, 7, 2 if shared else 11], 1),
                ([6, 2 if shared else 11, 1 if shared else 10], 1),
                ([1 if shared else 10, 2 if shared else 11, 8], 1),
                ([1 if shared else 10, 8, 9], 1)]
        return {"surfaces": surfaces, "policy": {"precision": 1e-8}}, {
            "mesh": {"vertices": vertices, "triangles": [
                {"vertices": nodes, "source_surface": source} for nodes, source in tris]}}

    def test_shared_wall_slab_line_is_required_not_only_area_coverage(self):
        inp, disconnected = self.model(False)
        self.assertFalse(check(inp, disconnected)["clean"])
        self.assertGreater(check(inp, disconnected)["missing_segment_count"], 0)
        _, conforming = self.model(True)
        self.assertTrue(check(inp, conforming)["clean"])

    def test_only_logged_bounded_move_can_justify_shortened_line(self):
        item = {"surfaces": [0, 1], "endpoints": [[0, 0, 0], [1, 0, 0]]}
        result = {"mesh": {"vertices": [[.01, 0, 0], [1, 0, 0]]},
                  "junction_regularization": {"accepted": []}}
        edges = {0: {(0, 1)}, 1: {(0, 1)}}
        surfaces = [{"rings": [[[0, 0, 0]]]}, {"rings": [[[2, 0, 0]]]}]
        self.assertFalse(repaired_segment_covered(item, result, edges, 1e-8, surfaces))
        result["junction_regularization"]["accepted"] = [{"moves": [{
            "from_coordinate": [0, 0, 0], "to_coordinate": [.01, 0, 0],
            "movement": .01, "movement_limit": .05}]}]
        self.assertTrue(repaired_segment_covered(item, result, edges, 1e-8, surfaces))
        result["junction_regularization"]["accepted"][0]["moves"][0]["movement_limit"] = .001
        self.assertFalse(repaired_segment_covered(item, result, edges, 1e-8, surfaces))

    def test_repair_provenance_rejects_forged_unrelated_or_ambiguous_moves(self):
        item = {"surfaces": [0, 1], "endpoints": [[0, 0, 0], [1, 0, 0]]}
        result = {"mesh": {"vertices": [[.01, 0, 0], [1, 0, 0]]},
                  "junction_regularization": {"accepted": [{"moves": [{
                      "from_coordinate": [0, 0, 0], "to_coordinate": [.01, 0, 0],
                      "movement": .001, "movement_limit": .05}]}]}}
        edges = {0: {(0, 1)}, 1: {(0, 1)}}
        surfaces = [{"rings": [[[0, 0, 0]]]}, {"rings": [[[2, 0, 0]]]}]
        self.assertFalse(repaired_segment_covered(item, result, edges, 1e-8, surfaces))
        move = result["junction_regularization"]["accepted"][0]["moves"][0]
        move["movement"] = .01
        self.assertTrue(repaired_segment_covered(item, result, edges, 1e-8, surfaces))
        surfaces[0]["rings"] = [[[3, 0, 0]]]
        self.assertFalse(repaired_segment_covered(item, result, edges, 1e-8, surfaces))
        surfaces[0]["rings"] = [[[0, 0, 0]]]
        result["junction_regularization"]["accepted"].append({"moves": [dict(move)]})
        self.assertFalse(repaired_segment_covered(item, result, edges, 1e-8, surfaces))


if __name__ == "__main__":
    unittest.main()
