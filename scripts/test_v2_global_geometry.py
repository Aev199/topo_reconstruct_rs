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


if __name__ == '__main__':
    unittest.main()
