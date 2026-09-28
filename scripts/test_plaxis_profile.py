"""Analytic fixtures for the PLAXIS readiness profile."""
import unittest
from check_plaxis_profile import profile
from test_v2_global_geometry import model, xy, wall, with_bars


def slab(ring, z=0):
    return ([0, 0, z], [1, 0, 0], [0, 1, 0], [0, 0, 1], [ring])


class ProfileTests(unittest.TestCase):
    def test_clean_junction_passes(self):
        r = profile(model([xy(), wall()]), element_size=0.5)
        self.assertTrue(r['passed'], r['counts'])

    def test_short_edge_and_sharp_corner(self):
        # A 20 mm step in the outline and a 5 degree needle.
        r = profile(model([slab([[0, 0], [1, 0], [1, 0.5], [1.02, 0.5], [1.02, 1], [0, 1]])]))
        self.assertEqual(r['counts'].get('short_edge'), 1)
        r = profile(model([slab([[0, 0], [1, 0], [0, 0.0875]])]))
        self.assertIn('sharp_corner', r['counts'])

    def test_narrow_face(self):
        # A slot 20 mm wide cut 0.5 m into a slab.
        ring = [[0, 0], [1, 0], [1, 1], [0.51, 1], [0.51, 0.5], [0.49, 0.5], [0.49, 1], [0, 1]]
        r = profile(model([slab(ring)]))
        self.assertIn('narrow_face', r['counts'])
        # 0.2 m wide: fine for 0.5 m elements.
        ring = [[0, 0], [1, 0], [1, 1], [0.6, 1], [0.6, 0.5], [0.4, 0.5], [0.4, 1], [0, 1]]
        self.assertNotIn('narrow_face', profile(model([slab(ring)]))['counts'])

    def test_gap_between_surfaces_and_bar_node(self):
        # A wall top 25 mm below a slab; a column ending 30 mm under it.
        r = profile(model([xy(), wall(y=0.5, z0=-1, z1=-0.025)]))
        self.assertIn('gap', r['counts'])
        data = with_bars(model([xy()]), [[(0.5, 0.5, -1), (0.5, 0.5, -0.03)]])
        self.assertEqual(profile(data)['counts'], {'gap': 1})
        # 0.1 m apart: not a gap at this element size.
        self.assertTrue(profile(model([xy(), wall(y=0.5, z0=-1, z1=-0.1)]))['passed'])

    def test_scale_follows_element_size(self):
        data = model([xy(), wall(y=0.5, z0=-1, z1=-0.025)])
        self.assertTrue(profile(data, element_size=0.1)['passed'])


if __name__ == '__main__':
    unittest.main()
