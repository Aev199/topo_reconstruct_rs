"""Independent checks of mesh_fragment JSON. No solver or mesher dependency."""
import argparse
import collections
import json
import math


def sub(a, b):
    return [x - y for x, y in zip(a, b)]


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def cross(a, b):
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]


def ring_area(ring):
    return abs(sum(a[0]*b[1]-a[1]*b[0] for a, b in zip(ring, ring[1:]+ring[:1])))/2


def inside(p, ring):
    result = False
    for a, b in zip(ring, ring[1:]+ring[:1]):
        if (a[1] > p[1]) != (b[1] > p[1]) and p[0] < a[0]+(p[1]-a[1])*(b[0]-a[0])/(b[1]-a[1]):
            result = not result
    return result


def check(data):
    top, mesh = data['topology'], data['mesh']
    model, axes = top['preview'], top['axis_assembly']['axes']
    vertices = mesh['vertices']
    eps = top['policy']['precision']
    assert vertices[:len(model['vertices'])] == model['vertices']
    assert mesh['surface_source_elements'] == [s['source_elements'] for s in model['surfaces']]
    assert top['all_surface_patches_built'] and top['axis_assembly']['all_axes_built']
    assert not top['issues'] and not top['axis_assembly']['issues']
    areas = collections.Counter()
    edges = [collections.Counter() for _ in model['surfaces']]
    used = [set() for _ in model['surfaces']]
    minimum_angle, maximum_area = 180., 0.
    for triangle in mesh['triangles']:
        ids, s = triangle['vertices'], triangle['surface']
        assert len(set(ids)) == 3
        assert triangle['stiffness'] == top['surface_stiffness'][s]
        points = [vertices[n] for n in ids]
        normal = cross(sub(points[1], points[0]), sub(points[2], points[0]))
        area = math.hypot(*normal)/2
        assert area > 0
        plane = model['planes'][model['surfaces'][s]['plane']]
        assert dot(normal, plane['normal']) > 0
        assert all(abs(dot(sub(p, plane['origin']), plane['normal'])) <= eps*2 for p in points)
        center = [sum(p[k] for p in points)/3 for k in range(3)]
        uv = [dot(sub(center, plane['origin']), plane[k]) for k in ('u', 'v')]
        rings = model['surfaces'][s]['contours']
        assert inside(uv, rings[0]) and not any(inside(uv, ring) for ring in rings[1:])
        areas[s] += area
        maximum_area = max(maximum_area, area)
        used[s].update(ids)
        for i in range(3):
            a, b = sub(points[(i+1)%3], points[i]), sub(points[(i+2)%3], points[i])
            angle = math.degrees(math.acos(max(-1., min(1., dot(a, b)/math.hypot(*a)/math.hypot(*b)))))
            minimum_angle = min(minimum_angle, angle)
            edges[s][tuple(sorted((ids[i], ids[(i+1)%3])))] += 1
    shared_chains = {}
    for s, surface in enumerate(model['surfaces']):
        target = ring_area(surface['contours'][0])-sum(map(ring_area, surface['contours'][1:]))
        assert abs(areas[s]-target) <= 10*eps*math.sqrt(target)
        expected_boundary = set()
        for use in sum(surface['boundaries'], []):
            e = use['edge']
            a, b = model['edges'][e]
            start, end = vertices[a], vertices[b]
            direction = sub(end, start)
            chain = []
            for n in used[s]:
                t = dot(sub(vertices[n], start), direction)/dot(direction, direction)
                q = [p+t*d for p, d in zip(start, direction)]
                if -1e-9 <= t <= 1+1e-9 and math.dist(q, vertices[n]) <= eps:
                    chain.append((t, n))
            chain = [n for _, n in sorted(chain)]
            assert chain[0] == a and chain[-1] == b
            if e in shared_chains:
                assert shared_chains[e] == chain
            shared_chains[e] = chain
            expected_boundary.update(tuple(sorted(pair)) for pair in zip(chain, chain[1:]))
        assert {edge for edge, count in edges[s].items() if count == 1} == expected_boundary
        assert all(count == (1 if edge in expected_boundary else 2) for edge, count in edges[s].items())
    lengths = collections.Counter()
    for bar in mesh['bars']:
        axis = axes[bar['axis']]
        spans = [s for s in axis['spans'] if s['element'] == bar['source_element']]
        assert len(spans) == 1 and spans[0]['stiffness'] == bar['stiffness']
        start, end = [vertices[n] for n in axis['endpoints']]
        direction = sub(end, start)
        ts = []
        for n in bar['vertices']:
            t = dot(sub(vertices[n], start), direction)/dot(direction, direction)
            assert math.dist(vertices[n], [p+t*d for p, d in zip(start, direction)]) <= eps
            assert spans[0]['start_t']-1e-9 <= t <= spans[0]['end_t']+1e-9
            ts.append(t)
        assert ts[1] > ts[0]
        lengths[bar['source_element']] += math.dist(*(vertices[n] for n in bar['vertices']))
        for contact in top['axis_assembly']['contacts']:
            if contact['kind'] == 'interval' and contact['axis'] == bar['axis']:
                if contact['start_t'] <= sum(ts)/2 <= contact['end_t']:
                    edge = tuple(sorted(bar['vertices']))
                    assert edges[contact['surface']][edge] == (1 if contact['location'] == 'boundary' else 2)
    for axis in axes:
        length = math.dist(*(vertices[n] for n in axis['endpoints']))
        for span in axis['spans']:
            assert abs(lengths[span['element']]-length*(span['end_t']-span['start_t'])) <= eps
    for c in top['axis_assembly']['contacts']:
        if c['kind'] == 'point':
            assert c['vertex'] in used[c['surface']]
            assert any(b['axis'] == c['axis'] and c['vertex'] in b['vertices'] for b in mesh['bars'])
    assert minimum_angle+1e-7 >= mesh['policy']['minimum_angle_degrees']
    assert maximum_area <= mesh['policy']['maximum_area']*(1+1e-7)
    assert math.isclose(minimum_angle, mesh['minimum_angle_degrees'], abs_tol=1e-7)
    assert math.isclose(maximum_area, mesh['maximum_triangle_area'], abs_tol=eps*eps)
    assert mesh['topology_valid'] and mesh['quality_passed'] and not mesh['blockers']
    assert not mesh['export_ready']
    return {'surfaces': len(areas), 'axes': len(axes), 'triangles': len(mesh['triangles']),
            'bars': len(mesh['bars']), 'minimum_angle_degrees': minimum_angle,
            'maximum_triangle_area': maximum_area}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('report')
    args = parser.parse_args()
    with open(args.report, encoding='utf-8') as file:
        print(json.dumps(check(json.load(file)), indent=2))
