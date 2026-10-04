"""Independent analytic fixtures for the global junction auditor."""
import unittest
from check_v2_global_geometry import audit


def model(polygons):
    result = dict(vertices=[], edges=[], planes=[], surfaces=[])
    vertex_ids, edge_ids = {}, {}
    for origin, u, v, normal, rings in polygons:
        plane = len(result['planes'])
        result['planes'].append(dict(origin=origin, u=u, v=v, normal=normal))
        boundaries = []
        for ring in rings:
            ids = []
            for x, y in ring:
                point = tuple(origin[k] + x*u[k] + y*v[k] for k in range(3))
                if point not in vertex_ids:
                    vertex_ids[point] = len(result['vertices'])
                    result['vertices'].append(point)
                ids.append(vertex_ids[point])
            boundary = []
            for a, b in zip(ids, ids[1:] + ids[:1]):
                edge = tuple(sorted((a, b)))
                if edge not in edge_ids:
                    edge_ids[edge] = len(result['edges'])
                    result['edges'].append(edge)
                boundary.append(dict(edge=edge_ids[edge], reversed=a>b))
            boundaries.append(boundary)
        result['surfaces'].append(dict(plane=plane, contours=rings,
                                      boundaries=boundaries, source_elements=[plane]))
    return dict(topology=dict(preview=result, policy=dict(precision=1e-7)))


def xy(x0=0, x1=1, z=0):
    return ([0,0,z], [1,0,0], [0,1,0], [0,0,1],
            [[[x0,0],[x1,0],[x1,1],[x0,1]]])


def wall(y=0, z0=0, z1=1):
    return ([0,y,0], [1,0,0], [0,0,1], [0,-1,0],
            [[[0,z0],[1,z0],[1,z1],[0,z1]]])


class GlobalAuditTests(unittest.TestCase):
    def test_report_precision_cannot_loosen_the_audit(self):
        data = model([xy(), wall(y=.5)])
        data['topology']['policy']['precision'] = 1e-3
        r = audit(data)
        self.assertAlmostEqual(r['precision'], 5e-6)
        self.assertFalse(r['global_surface_checks_passed'])
        self.assertAlmostEqual(audit(data, maximum_precision=1e-3)['precision'], 5e-3)

    def test_shared_boundary(self):
        r = audit(model([xy(), wall()]))
        self.assertTrue(r['global_surface_checks_passed'])
        self.assertEqual(len(r['contacts']), 1)
        self.assertTrue(r['contacts'][0]['geometry_conforming'])

    def test_t_junction(self):
        r = audit(model([xy(), wall(y=.5)]))
        self.assertFalse(r['global_surface_checks_passed'])
        self.assertEqual(r['issues'][0]['contact_kind'], 't_junction')
        self.assertAlmostEqual(r['issues'][0]['length'], 1)

    def test_crossing(self):
        r = audit(model([xy(), wall(y=.5, z0=-1)]))
        self.assertEqual(r['issues'][0]['contact_kind'], 'crossing')

    def test_coplanar_shared_boundary(self):
        r = audit(model([xy(), xy(1,2)]))
        self.assertTrue(r['global_surface_checks_passed'])
        self.assertEqual(len(r['contacts']), 1)

    def test_overlap(self):
        r = audit(model([xy(), xy(.5,1.5)]))
        self.assertEqual(r['issues'][0]['kind'], 'coplanar_overlap')
        self.assertAlmostEqual(r['issues'][0]['area'], .5)

    def test_near_faces(self):
        r = audit(model([xy(), xy(z=.01)]))
        self.assertEqual(len(r['near_faces']), 1)
        self.assertFalse(r['issues'])

    def test_disjoint(self):
        r = audit(model([xy(), xy(2,3), xy(z=1)]))
        self.assertFalse(r['contacts'])
        self.assertFalse(r['issues'])

    def test_hole_splits_contact(self):
        p = xy()
        p[-1].append([[.25,.25],[.25,.75],[.75,.75],[.75,.25]])
        r = audit(model([p, wall(y=.5, z0=-1)]))
        self.assertEqual(len(r['contacts']), 2)
        self.assertAlmostEqual(sum(c['length'] for c in r['contacts']), .5)

    def test_mesh_requires_shared_identifiers(self):
        data = model([xy(), wall()])
        vertices = data['topology']['preview']['vertices']
        # Surface 0 and 1 share model edge (0,1); mesh duplicates its endpoints.
        data['mesh'] = dict(vertices=vertices + [vertices[0], vertices[1]],
                            triangles=[dict(surface=0, vertices=[0,1,2]),
                                       dict(surface=1, vertices=[6,7,4])])
        r = audit(data)
        self.assertTrue(r['contacts'][0]['geometry_conforming'])
        self.assertFalse(r['contacts'][0]['mesh_conforming'])
        data['mesh']['triangles'][1]['vertices'] = [0,1,4]
        self.assertTrue(audit(data)['global_surface_checks_passed'])

    def test_translation_and_reversed_normal(self):
        data = model([xy(), wall(y=.5,z0=-1)])
        for p in data['topology']['preview']['planes']:
            p['origin'] = [q+1e6 for q in p['origin']]
            p['normal'] = [-q for q in p['normal']]
        data['topology']['preview']['vertices'] = [
            [q+1e6 for q in p] for p in data['topology']['preview']['vertices']]
        r = audit(data)
        self.assertEqual(r['issues'][0]['contact_kind'], 'crossing')
        self.assertAlmostEqual(r['issues'][0]['length'], 1)

    @staticmethod
    def wall_edge(data, z):
        preview = data['topology']['preview']
        ids = [i for i, p in enumerate(preview['vertices'])
               if abs(p[1]-.5) < 1e-12 and abs(p[2]-z) < 1e-12]
        return preview['edges'].index(tuple(sorted(ids)))

    def test_embedded_edge_represents_t_junction(self):
        data = model([xy(), wall(y=.5)])
        preview = data['topology']['preview']
        preview['surfaces'][0]['embedded_edges'] = [self.wall_edge(data, 0)]
        r = audit(data)
        self.assertTrue(r['global_surface_checks_passed'], r['issues'])
        self.assertEqual(r['contacts'][0]['kind'], 't_junction')
        self.assertTrue(r['contacts'][0]['geometry_conforming'])

    def test_embedded_edge_must_lie_inside_surface(self):
        data = model([xy(), wall(y=.5)])
        preview = data['topology']['preview']
        # The wall top edge is off the slab plane and outside its material.
        preview['surfaces'][0]['embedded_edges'] = [self.wall_edge(data, 1)]
        r = audit(data)
        self.assertFalse(r['global_surface_checks_passed'])
        self.assertTrue(r['invalid_surfaces'])
        # A boundary edge listed as embedded is also rejected.
        preview['surfaces'][0]['embedded_edges'] = [preview['surfaces'][0]['boundaries'][0][0]['edge']]
        self.assertEqual(audit(data)['invalid_surfaces'][0]['reason'], 'invalid embedded edges')


def with_bars(data, bars, contacts=()):
    """Append bar axes given by coordinates; equal points share one vertex."""
    preview = data['topology']['preview']
    vertices = preview['vertices']
    index = {tuple(p): i for i, p in enumerate(vertices)}
    def vid(p):
        if tuple(p) not in index:
            index[tuple(p)] = len(vertices)
            vertices.append(tuple(p))
        return index[tuple(p)]
    axes = []
    for points in bars:
        ids = [vid(p) for p in points]
        axes.append(dict(endpoints=[ids[0], ids[-1]], spans=[],
                         anchors=[dict(vertex=v, t=k / (len(ids) - 1)) for k, v in enumerate(ids)]))
    data['topology']['axis_assembly'] = dict(
        axes=axes, contacts=[dict(c, vertex=vid(c['vertex'])) if 'vertex' in c else c
                             for c in contacts])
    return data


class BarAndPointTests(unittest.TestCase):
    def test_crossing_bars_need_a_shared_node(self):
        data = with_bars(model([xy(z=5)]), [[(0, 0, 1), (2, 2, 1)], [(0, 2, 1), (2, 0, 1)]])
        r = audit(data)
        self.assertEqual(r['point_and_bar_issue_counts'], {'unshared_bar_intersection': 1})
        self.assertFalse(r['global_checks_passed'])
        data = with_bars(model([xy(z=5)]), [[(0, 0, 1), (1, 1, 1), (2, 2, 1)],
                                           [(0, 2, 1), (1, 1, 1), (2, 0, 1)]])
        self.assertTrue(audit(data)['global_checks_passed'])

    def test_bar_piercing_a_panel_needs_a_contact(self):
        column = [(0.5, 0.5, -1), (0.5, 0.5, 0), (0.5, 0.5, 1)]
        data = with_bars(model([xy()]), [column])
        r = audit(data)
        self.assertEqual(r['point_and_bar_issue_counts'], {'unshared_bar_surface_intersection': 1})
        data = with_bars(model([xy()]), [column],
                         [dict(kind='point', axis=0, surface=0, vertex=(0.5, 0.5, 0))])
        self.assertTrue(audit(data)['global_checks_passed'])

    def test_bar_in_panel_needs_an_interval_contact(self):
        data = with_bars(model([xy()]), [[(0.1, 0.5, 0), (0.9, 0.5, 0)]])
        r = audit(data)
        self.assertEqual(r['point_and_bar_issue_counts'], {'bar_in_surface_without_contact': 1})
        data['topology']['axis_assembly']['contacts'] = [
            dict(kind='interval', axis=0, surface=0, start_t=0, end_t=1)]
        self.assertTrue(audit(data)['global_checks_passed'])

    def test_short_bars_and_gaps_are_review_items(self):
        data = with_bars(model([xy(z=5)]), [[(0, 0, 1), (1, 0, 1)], [(1.02, 0, 1), (2, 0, 1)],
                                           [(0, 1, 1), (1, 1, 1)], [(1, 1, 1), (1.01, 1, 1)],
                                           [(1.01, 1, 1), (2, 1, 1)]])
        r = audit(data)
        self.assertTrue(r['global_checks_passed'])
        self.assertEqual(r['review_counts'], {'bar_near_miss': 1, 'short_bar': 1})

    def test_corner_touching_a_panel_needs_a_shared_vertex(self):
        # A tilted panel whose corner touches the slab interior at one point.
        tilted = ([0.5, 0.5, 0], [1, 0, 0], [0, 0.6, 0.8], [0, -0.8, 0.6],
                  [[[0, 0], [0.3, 0.2], [0, 1]]])
        r = audit(model([xy(), tilted]))
        self.assertIn('unshared_point_contact', r['point_and_bar_issue_counts'])
        self.assertFalse(r['global_checks_passed'])

    def test_property_transfer_mismatch_fails(self):
        data = model([xy()])
        data['topology']['surface_stiffness'] = [7]
        data['mesh'] = dict(vertices=data['topology']['preview']['vertices'], bars=[],
                            triangles=[dict(surface=0, vertices=[0, 1, 2], stiffness=8)])
        r = audit(data)
        self.assertEqual(r['properties']['triangles_with_wrong_stiffness'], 1)
        self.assertFalse(r['global_checks_passed'])

    def test_partial_interval_contact_does_not_cover_the_line(self):
        data = with_bars(model([xy()]), [[(0.1, 0.5, 0), (0.9, 0.5, 0)]],
                         [dict(kind='interval', axis=0, surface=0, start_t=0, end_t=0.2)])
        r = audit(data)
        self.assertEqual(r['point_and_bar_issue_counts'], {'bar_in_surface_without_contact': 1})
        data['topology']['axis_assembly']['contacts'].append(
            dict(kind='interval', axis=0, surface=0, start_t=0.2, end_t=1))
        self.assertTrue(audit(data)['global_checks_passed'])

    def test_bar_crossing_a_panel_edge_needs_full_coverage(self):
        # In the slab plane, the bar leaves the slab: only the inside part is
        # required, and a contact covering it is sufficient.
        data = with_bars(model([xy()]), [[(0.5, 0.5, 0), (1.5, 0.5, 0)]],
                         [dict(kind='interval', axis=0, surface=0, start_t=0, end_t=0.25)])
        self.assertEqual(audit(data)['point_and_bar_issue_counts'],
                         {'bar_in_surface_without_contact': 1})
        data['topology']['axis_assembly']['contacts'][0]['end_t'] = 0.5
        self.assertTrue(audit(data)['global_checks_passed'])


def meshed(data, triangles, bars, extra=()):
    """Trial mesh over the preview vertices plus extra (duplicate) nodes."""
    preview = data['topology']['preview']
    data['topology']['surface_stiffness'] = [1] * len(preview['surfaces'])
    for axis in data['topology'].get('axis_assembly', {}).get('axes', []):
        axis['spans'] = [dict(element=1, stiffness=2, start_t=0, end_t=1)]
    data['mesh'] = dict(vertices=list(preview['vertices']) + list(extra),
                        triangles=[dict(surface=s, vertices=v, stiffness=1) for s, v in triangles],
                        bars=[dict(axis=a, vertices=v, source_element=1, stiffness=2)
                              for a, v in bars])
    return data


class MeshConnectivityTests(unittest.TestCase):
    """Contact records never prove that the mesh is connected."""

    def column(self):
        data = with_bars(model([xy()]), [[(0.5, 0.5, -1), (0.5, 0.5, 0), (0.5, 0.5, 1)]],
                         [dict(kind='point', axis=0, surface=0, vertex=(0.5, 0.5, 0))])
        # Preview vertices: 0-3 slab corners, 4-6 column nodes (5 = pierce point).
        return data

    def test_piercing_bar_must_share_a_mesh_node(self):
        fan = [(0, [0, 1, 5]), (0, [1, 2, 5]), (0, [2, 3, 5]), (0, [3, 0, 5])]
        data = meshed(self.column(), fan, [(0, [4, 5]), (0, [5, 6])])
        self.assertTrue(audit(data)['global_checks_passed'])
        # Same geometry, but the bar uses a duplicate node at the pierce point.
        data = meshed(self.column(), fan, [(0, [4, 7]), (0, [7, 6])], extra=[(0.5, 0.5, 0)])
        r = audit(data)
        self.assertEqual(r['point_and_bar_issue_counts'],
                         {'bar_surface_point_not_shared_in_mesh': 1})
        self.assertFalse(r['global_checks_passed'])

    def test_bar_in_panel_must_follow_triangle_edges(self):
        data = with_bars(model([xy()]), [[(0, 0.5, 0), (1, 0.5, 0)]],
                         [dict(kind='interval', axis=0, surface=0, start_t=0, end_t=1)])
        # Preview vertices 4 and 5 are the bar ends on the slab edges.
        split = [(0, [0, 1, 5]), (0, [0, 5, 4]), (0, [4, 5, 2]), (0, [4, 2, 3])]
        self.assertTrue(audit(meshed(data, split, [(0, [4, 5])]))['global_checks_passed'])
        diagonal = [(0, [0, 1, 2]), (0, [0, 2, 3])]
        data = with_bars(model([xy()]), [[(0, 0.5, 0), (1, 0.5, 0)]],
                         [dict(kind='interval', axis=0, surface=0, start_t=0, end_t=1)])
        r = audit(meshed(data, diagonal, [(0, [4, 5])]))
        self.assertIn('bar_in_surface_not_shared_in_mesh', r['point_and_bar_issue_counts'])

    def test_unmeshed_bar_fails(self):
        fan = [(0, [0, 1, 5]), (0, [1, 2, 5]), (0, [2, 3, 5]), (0, [3, 0, 5])]
        r = audit(meshed(self.column(), fan, [(0, [4, 5])]))
        self.assertIn('bar_not_covered_by_mesh', r['point_and_bar_issue_counts'])

    def test_crossing_bars_must_share_a_mesh_node(self):
        data = with_bars(model([xy(z=5)]), [[(0, 0, 1), (1, 1, 1), (2, 2, 1)],
                                           [(0, 2, 1), (1, 1, 1), (2, 0, 1)]])
        slab = [(0, [0, 1, 2]), (0, [0, 2, 3])]
        # Vertices 4-6 and 7, 5, 8 form the two bars through shared node 5.
        good = meshed(data, slab, [(0, [4, 5]), (0, [5, 6]), (1, [7, 5]), (1, [5, 8])])
        self.assertTrue(audit(good)['global_checks_passed'])
        bad = meshed(data, slab, [(0, [4, 5]), (0, [5, 6]), (1, [7, 9]), (1, [9, 8])],
                     extra=[(1, 1, 1)])
        self.assertEqual(audit(bad)['point_and_bar_issue_counts'],
                         {'bar_intersection_not_shared_in_mesh': 1})

    def test_shared_surface_vertex_must_be_one_mesh_node(self):
        data = model([xy(), wall()])
        good = meshed(data, [(0, [0, 1, 2]), (0, [0, 2, 3]), (1, [0, 1, 4]), (1, [0, 4, 5])], [])
        self.assertTrue(audit(good)['global_checks_passed'])
        data = model([xy(), wall()])
        bad = meshed(data, [(0, [0, 1, 2]), (0, [0, 2, 3]), (1, [6, 7, 4]), (1, [6, 4, 5])], [],
                     extra=data['topology']['preview']['vertices'][:2])
        r = audit(bad)
        self.assertEqual(r['point_and_bar_issue_counts'], {'shared_vertex_not_shared_in_mesh': 2})


class PropertyTests(unittest.TestCase):
    def slab(self, triangle):
        data = model([xy(), xy(1, 2)])
        data['topology']['surface_stiffness'] = [7, 7]
        data['mesh'] = dict(vertices=data['topology']['preview']['vertices'], bars=[],
                            triangles=[dict(surface=0, vertices=[0, 1, 2], stiffness=7), triangle])
        return data

    def test_missing_stiffness_field_fails(self):
        self.assertTrue(audit(self.slab(dict(surface=1, vertices=[1, 5, 2], stiffness=7)))
                        ['global_checks_passed'])
        r = audit(self.slab(dict(surface=1, vertices=[1, 5, 2])))
        self.assertEqual(r['properties']['triangles_with_missing_stiffness'], 1)
        self.assertFalse(r['global_checks_passed'])

    def test_missing_surface_record_fails(self):
        data = self.slab(dict(surface=1, vertices=[1, 5, 2], stiffness=7))
        data['topology']['surface_stiffness'] = [7]
        r = audit(data)
        self.assertEqual(r['properties']['triangles_with_missing_stiffness'], 1)
        self.assertFalse(r['global_checks_passed'])
        del data['topology']['surface_stiffness']
        self.assertEqual(audit(data)['properties']['triangles_with_missing_stiffness'], 2)

    def test_unmeshed_surface_and_bar_without_spans_fail(self):
        data = self.slab(dict(surface=0, vertices=[0, 2, 3], stiffness=7))
        r = audit(data)
        self.assertEqual(r['properties']['surfaces_without_triangles'], [1])
        self.assertFalse(r['global_checks_passed'])
        data = with_bars(self.slab(dict(surface=1, vertices=[1, 5, 2], stiffness=7)),
                         [[(0, 0, 3), (1, 0, 3)]])
        n = len(data['topology']['preview']['vertices'])
        data['mesh']['vertices'] = data['topology']['preview']['vertices']
        data['mesh']['bars'] = [dict(axis=0, vertices=[n - 2, n - 1], source_element=1, stiffness=2)]
        r = audit(data)
        self.assertEqual(r['properties']['bars_with_missing_stiffness'], 1)
        data['mesh']['bars'][0]['axis'] = 5
        self.assertEqual(audit(data)['properties']['bars_with_missing_stiffness'], 1)


if __name__ == '__main__':
    unittest.main()
