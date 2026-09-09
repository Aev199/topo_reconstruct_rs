"""Independent provenance and transformation checks for the v2 plane example.
Usage: python scripts/check_v2_axes.py model.txt target/debug/examples/recognize_planes
Only temporary geometry copies are created. Standard library only.
"""
import json
import math
from pathlib import Path
import re
import subprocess
import sys
import tempfile


def check(model, executable):
    text = Path(model).read_text(encoding='utf-8-sig')
    def block(number):
        start = re.search(r'\(\s*' + str(number) + r'\s*/', text).end()
        return [r.split() for r in text[start:text.index(')', start)].split('/') if r.strip()]
    nodes = {i: list(map(float, r)) for i, r in enumerate(block(4), 1)}
    elements = {i: list(map(int, r)) for i, r in enumerate(block(1), 1)}
    expected = {i for i, r in elements.items() if (r[0] in (41, 44) and len(r) == 6) or (r[0] == 42 and len(r) == 5)}
    binary = str(Path(executable).resolve())
    def run(path):
        return json.loads(subprocess.check_output([binary, str(path)]))
    current_nodes = nodes
    def memberships(report, mapping):
        represented = [mapping[e] for p in report['patches'] for e in p['source_elements']]
        rejected = [mapping[e['element']] for e in report['rejected']]
        assert len(represented + rejected) == len(set(represented + rejected))
        assert set(represented + rejected) == expected
        for patch in report['patches']:
            assert patch['maximum_deviation'] <= report['policy']['distance'] + 1e-9
            normal = patch['plane']['normal']
            assert abs(sum(v*v for v in normal)-1) < 1e-9
            origin = patch['plane']['origin']
            measured = max(abs(sum((current_nodes[n][i]-origin[i])*normal[i] for i in range(3))) for n in patch['source_nodes'])
            assert abs(measured-patch['maximum_deviation']) < 1e-8
            assert measured <= report['policy']['distance'] + 1e-9
            property_ids = []
            for stiffness, members in patch['stiffness_regions'].items():
                for element in members:
                    assert elements[mapping[element]][1] == int(stiffness)
                property_ids.extend(members)
            assert sorted(property_ids) == sorted(patch['source_elements'])
        return {tuple(sorted(mapping[e] for e in p['source_elements'])) for p in report['patches']}
    baseline = run(Path(model).resolve())
    assert baseline == run(Path(model).resolve()), 'Non-deterministic output'
    groups = memberships(baseline, {i: i for i in elements})
    angle = math.radians(37)
    with tempfile.TemporaryDirectory(prefix='v2-axes-') as directory:
        cases = [
            ('translated', lambda p: [p[0] + 100, p[1] - 200, p[2] + 30], False),
            ('rotated', lambda p: [math.cos(angle)*p[0]-math.sin(angle)*p[1], math.sin(angle)*p[0]+math.cos(angle)*p[1], p[2]], False),
            ('renumbered', lambda p: p, True),
        ]
        for name, transform, reverse in cases:
            ns = list(nodes)[::-1] if reverse else list(nodes)
            es = list(elements)[::-1] if reverse else list(elements)
            node_ids = {old: new for new, old in enumerate(ns, 1)}
            mapping = {new: old for new, old in enumerate(es, 1)}
            coordinates = '/'.join(' '.join(format(v, '.15g') for v in transform(nodes[i])) for i in ns)
            rows = '/'.join(' '.join(map(str, elements[i][:2] + [node_ids[n] for n in elements[i][2:]])) for i in es)
            path = Path(directory) / (name + '.txt')
            path.write_text('(4/' + coordinates + '/)\n(1/' + rows + '/)', encoding='utf-8')
            current_nodes = {node_ids[n]: transform(nodes[n]) for n in ns}
            assert memberships(run(path), mapping) == groups, name
            print(name + ': plane membership, FE coverage and property regions preserved')
    print(f"Verified {len(expected)} source shells, {len(groups)} plane patches; repeat identical")


if __name__ == '__main__':
    check(*sys.argv[1:])
