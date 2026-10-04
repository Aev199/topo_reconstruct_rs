"""The PLAXIS loader against a recording stand-in for PLAXIS Input."""
import unittest
from plaxis_export import Recorder, build

DATA = {
    "format": "topo-plaxis-1",
    "plate_materials": [{"name": "GEI_7_h200", "stiffness": 7, "e": 2.94e7, "nu": 0.2, "d": 0.2, "gamma": 24.5}],
    "beam_materials": [{"name": "S0_1_50x80", "stiffness": 1, "e": 2.94e7, "nu": 0.2, "width": 0.5,
                        "height": 0.8, "a": 0.4, "i2": 0.02133, "i3": 0.00833, "gamma": 24.5}],
    "plates": [
        {"surface": 0, "stiffness": 7, "material": "GEI_7_h200",
         "polygons": [[[0, 0, 0], [4, 0, 0], [4, 3, 0]], [[0, 0, 0], [4, 3, 0], [0, 3, 0]]]},
        {"surface": 1, "stiffness": 9, "material": None, "polygons": [[[0, 0, 0], [0, 3, 0], [0, 3, 3]]]},
    ],
    "beams": [{"bar": 0, "stiffness": 1, "material": "S0_1_50x80", "start": [2, 1, 0], "end": [2, 1, -3]}],
    "missing_materials": [9],
}


class LoaderTests(unittest.TestCase):
    def test_objects_and_materials(self):
        g = Recorder()
        r = build(g, DATA, progress=lambda _: None)
        self.assertEqual((r["plates"], r["beams"], r["plate_materials"], r["beam_materials"]), (3, 1, 1, 1))
        self.assertEqual(g.commands[0], "gotostructures")
        self.assertIn('surface (0.0 0.0 0.0) (4.0 0.0 0.0) (4.0 3.0 0.0)', g.commands)
        self.assertIn("set Plate_1.Material PlateMat_1", g.commands)
        self.assertNotIn("set Plate_3.Material PlateMat_1", g.commands)
        self.assertIn("set Beam_1.Material BeamMat_1", g.commands)
        plate = next(c for c in g.commands if c.startswith("setproperties PlateMat_1"))
        self.assertIn('"D3d" 0.2', plate)
        self.assertIn('"StructNu12" 0.2', plate)
        self.assertEqual(r["property_sets"], ["beam set 0", "plate set 0"])

    def test_older_property_names_are_tried_next(self):
        g = Recorder(reject={"D3d", "CrossSectionType"})
        r = build(g, DATA, progress=lambda _: None)
        self.assertEqual(r["property_sets"], ["beam set 1", "plate set 1"])
        self.assertEqual(len(r["rejected_property_sets"]), 2)
        self.assertIn("delete PlateMat_1", g.commands)
        self.assertIn("set Plate_1.Material PlateMat_2", g.commands)

    def test_shift(self):
        g = Recorder()
        build(g, DATA, shift=(2.0, 1.0, 0.0), progress=lambda _: None)
        self.assertIn("line (0.0 0.0 0.0) (0.0 0.0 -3.0)", g.commands)


if __name__ == "__main__":
    unittest.main()
