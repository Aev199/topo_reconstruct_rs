"""Check nonlinear proposal geometry directly from JSON, without solver matrices."""
import argparse
import json
import math


def sub(a, b):
    return [x - y for x, y in zip(a, b)]


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def unit(a):
    length = math.hypot(*a)
    return [x / length for x in a]


def cross(a, b):
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]


def check(data):
    r = data['frame']
    assert data['proposal_only']
    points, original = r['candidate_points'], r['reference_points']
    policy = r['policy']
    tolerance = policy['residual_tolerance']
    up = unit(policy['up'])
    u = unit(cross([0, 0, 1] if abs(up[2]) < .9 else [1, 0, 0], up))
    v = cross(up, u)
    errors = []
    budgets = [policy['maximum_movement']] * len(points)
    valid_parameters = True
    moved_parameters = 0
    invalid_axes = set()
    assert len(r['axes']) == len(r['sliding_parameters']) == len(data['axis_recognition']['axes'])
    for index, (axis, ts, source) in enumerate(zip(r['axes'], r['sliding_parameters'], data['axis_recognition']['axes'])):
        assert axis['spans'] == source['spans']
        assert len(ts) == len(axis['anchors'])
        a, b = axis['endpoints']
        direction = sub(points[b], points[a])
        reference = sub(original[b], original[a])
        length = math.hypot(*reference)
        new_length = math.hypot(*direction)
        if (new_length + tolerance < min(policy['minimum_length'], length)
                or dot(direction, reference) <= 0
                or new_length == 0
                or dot(unit(direction), unit(reference)) < math.cos(policy['angle'])):
            invalid_axes.add(index)
        cap = min(policy['maximum_movement'], policy['relative_movement'] * length)
        ordered = sorted(zip(axis['anchors'], ts), key=lambda pair: pair[0]['t'])
        valid_parameters &= all(math.isfinite(t) and 0 <= t <= 1 for t in ts)
        valid_parameters &= all(x[1] < y[1] for x, y in zip(ordered, ordered[1:]))
        for anchor, t, source_anchor in zip(axis['anchors'], ts, source['anchors']):
            assert anchor['t'] == source_anchor['t']
            n = anchor['node']
            assert r['node_ids'][n] == source_anchor['node']
            budgets[n] = min(budgets[n], cap)
            errors.extend(abs(points[n][k] - points[a][k] - t*direction[k]) for k in range(3))
            moved_parameters += abs(t - anchor['t']) > 1e-10
        if length < policy['minimum_length']:
            errors.extend(abs(x) for x in sub(direction, reference))
        elif abs(dot(unit(reference), up)) >= math.cos(policy['angle']):
            errors.extend([abs(dot(direction, u)), abs(dot(direction, v))])
        elif abs(dot(unit(reference), up)) <= math.sin(policy['angle']):
            errors.append(abs(dot(direction, up)))
    assert len(r['surfaces']) == len(r['candidate_planes'])
    for surface, plane in zip(r['surfaces'], r['candidate_planes']):
        errors.extend(abs(dot(sub(points[n], plane['origin']), plane['normal'])) for n in surface['nodes'])
    maximum = max(errors, default=0)
    assert math.isclose(maximum, r['candidate_max_residual'], abs_tol=tolerance * .01, rel_tol=1e-7)
    assert (maximum <= tolerance) == r['candidate_constraints_satisfied']
    assert valid_parameters == r['candidate_parameters_valid']
    movement = [math.dist(p, q) for p, q in zip(points, original)]
    failures = {r['node_ids'][i] for i, (d, cap) in enumerate(zip(movement, budgets)) if d > cap + 1e-10}
    assert failures == {f['node_id'] for f in r['movement_failures']}
    assert invalid_axes == {f['axis'] for f in r['axis_failures']}
    if r['accepted']:
        assert r['points'] == points and not failures and valid_parameters and maximum <= tolerance
        assert not r['axis_failures']
    else:
        assert r['points'] == original and r['maximum_movement'] == 0
    return {'accepted': r['accepted'], 'maximum_residual': maximum,
            'maximum_movement': max(movement, default=0),
            'changed_parameters': moved_parameters, 'over_budget_nodes': len(failures),
            'invalid_axes': len(invalid_axes)}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('report')
    args = parser.parse_args()
    with open(args.report, encoding='utf-8') as file:
        print(json.dumps(check(json.load(file)), indent=2))
