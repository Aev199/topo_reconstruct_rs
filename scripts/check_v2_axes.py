"""Independent provenance and transformation checks for the v2 axis example.
Usage: python scripts/check_v2_axes.py model.txt target/debug/examples/recognize_lira
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
    expected = {i for i, r in elements.items() if r[0] == 10 and len(r) == 4}
    binary = str(Path(executable).resolve())
    def run(path):
        return json.loads(subprocess.check_output([binary, str(path)]))
    def memberships(report, mapping):
        represented = [mapping[s['element']] for a in report['axes'] for s in a['spans']]
        rejected = [mapping[e['element']] for e in report['rejected']]
        assert len(represented + rejected) == len(set(represented + rejected))
        assert set(represented + rejected) == expected
        for axis in report['axes']:
            assert len(axis['endpoints']) == 2
            spans = axis['spans']
            assert abs(spans[0]['start_t']) < 1e-9 and abs(spans[-1]['end_t'] - 1) < 1e-9
            for span in spans:
                assert span['start_t'] < span['end_t']
                assert span['stiffness'] == elements[mapping[span['element']]][1]
            assert all(abs(a['end_t'] - b['start_t']) < 1e-9 for a, b in zip(spans, spans[1:]))
            for anchor in axis['anchors']:
                assert anchor['distance_to_axis'] <= 0.01 + 1e-9
        return {tuple(sorted(mapping[s['element']] for s in a['spans'])) for a in report['axes']}
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
            assert memberships(run(path), mapping) == groups, name
            print(name + ': axis membership, FE coverage and property spans preserved')
    print(f"Verified {len(expected)} source bars, {len(groups)} axes; repeat identical")


if __name__ == '__main__':
    check(*sys.argv[1:])
