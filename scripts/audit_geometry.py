"""Independent geometric audit. Usage: python scripts/audit_geometry.py model.txt report.json.
Requires numpy and shapely. Does not modify source files.
"""
import argparse
import json
import re
from collections import Counter
from pathlib import Path
import numpy as np
from shapely.geometry import Polygon
from shapely.ops import unary_union


def read_model(path):
    text = Path(path).read_bytes().decode('utf-8-sig', errors='replace')
    def block(number):
        match = re.search(r'\(\s*' + str(number) + r'\s*/', text)
        return [r.split() for r in text[match.end():text.index(')', match.end())].split('/') if r.strip()]
    nodes = {i: np.array(list(map(float, r))) for i, r in enumerate(block(4), 1)}
    elements = {i: list(map(int, r)) for i, r in enumerate(block(1), 1)}
    return nodes, elements


def audit(model, report):
    nodes, elements = read_model(model)
    panels = report['panels']
    represented = [i for p in panels for i in p['source_element_ids']]
    shells = {i for i, e in elements.items() if e[0] in (41, 42, 44)}
    deviations, invalid, shapes = [], [], []
    narrow = []
    vertices = []
    edges = []
    max_plane_error = 0.0
    for p in panels:
        normal = np.array(p['plane_normal'])
        helper = np.array([0.,0.,1.]) if abs(normal[2]) < .9 else np.array([1.,0.,0.])
        u = np.cross(normal, helper); u /= np.linalg.norm(u)
        v = np.cross(normal, u)
        def project(pts):
            pts=np.array(pts); return np.column_stack((pts@u, pts@v))
        rings = [project(r) for r in p['polygons']]
        shape = Polygon(rings[0], rings[1:]); shapes.append((p,shape,u,v))
        if not shape.is_valid: invalid.append(p['id'])
        if shape.minimum_clearance < 0.03:
            narrow.append(dict(panel=p['id'],minimum_clearance=shape.minimum_clearance))
        for ring in p['polygons']:
            ring=np.array(ring)
            max_plane_error=max(max_plane_error,float(np.max(np.abs(ring@normal+p['plane_d']))))
            for start,end in zip(ring,np.roll(ring,-1,axis=0)):
                vertices.append((p['id'],start)); edges.append((p['id'],start,end))
        original=[]
        for i in p['source_element_ids']:
            pts=project([nodes[n] for n in elements[i][2:]])
            mid=pts.mean(axis=0); delta=pts-mid
            pts=pts[np.argsort(np.arctan2(delta[:,1],delta[:,0]))]
            polygon=Polygon(pts)
            if polygon.is_valid and polygon.area>0: original.append(polygon)
        if original and shape.is_valid:
            source=unary_union(original)
            difference=shape.symmetric_difference(source).area
            deviations.append(dict(panel=p['id'],source_area=source.area,output_area=shape.area,
                                   symmetric_difference=difference,relative=difference/source.area))
    overlaps=[]
    for i,(a,pa,u,v) in enumerate(shapes):
        for b,pb,_,_ in shapes[i+1:]:
            if np.linalg.norm(np.array(a['plane_normal'])-np.array(b['plane_normal']))>1e-8: continue
            if abs(a['plane_d']-b['plane_d'])>1e-6: continue
            rings=[np.column_stack((np.array(r)@u,np.array(r)@v)) for r in b['polygons']]
            pb=Polygon(rings[0],rings[1:])
            if pa.is_valid and pb.is_valid:
                area=pa.intersection(pb).area
                if area>1e-8: overlaps.append(dict(panels=[a['id'],b['id']],area=area))
    near_junctions=0
    if edges:
        starts=np.array([e[1] for e in edges]); ends=np.array([e[2] for e in edges])
        owners=np.array([e[0] for e in edges]); direction=ends-starts
        lengths=np.sum(direction*direction,axis=1)
        for owner,point in vertices:
            t=np.divide(np.sum((point-starts)*direction,axis=1),lengths,out=np.zeros_like(lengths),where=lengths>0)
            distances=np.linalg.norm(point-(starts+np.clip(t,0,1)[:,None]*direction),axis=1)
            near_junctions+=int(np.sum((owners!=owner)&(t>1e-6)&(t<1-1e-6)&(distances>=1e-7)&(distances<0.01)))
    return dict(narrow_feature_threshold=0.03,junction_distance_threshold=0.01,
        narrow_features=narrow,near_vertex_edge_candidates=near_junctions,
        max_plane_error=max_plane_error,
        nodes=len(nodes),elements=len(elements),types=dict(Counter(e[0] for e in elements.values())),
        source_shells=len(shells),represented_shells=len(shells.intersection(represented)),
        missing_shells=sorted(shells-set(represented)),non_shells_in_panels=sorted(set(represented)-shells),
        duplicate_source_ids=[i for i,c in Counter(represented).items() if c>1],
        panels=len(panels),invalid_polygons=invalid,coplanar_overlaps=overlaps,
        total_projected_source_area=sum(d['source_area'] for d in deviations),
        total_output_area=sum(d['output_area'] for d in deviations),
        total_symmetric_difference=sum(d['symmetric_difference'] for d in deviations),
        largest_changes=sorted(deviations,key=lambda d:d['relative'],reverse=True)[:10])

if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('model'); parser.add_argument('report')
    args=parser.parse_args()
    print(json.dumps(audit(args.model,json.loads(Path(args.report).read_text())),ensure_ascii=False,indent=2))
