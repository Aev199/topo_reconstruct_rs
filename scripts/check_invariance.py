"""Model-independent checks under repeated execution, translation and renumbering.
Usage: python scripts/check_invariance.py model.txt path/to/binary
Requires numpy (through audit_geometry). Temporary input copies contain geometry only.
"""
import argparse
import json
import subprocess
import tempfile
from pathlib import Path
import numpy as np
from audit_geometry import read_model


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('model'); parser.add_argument('binary')
    args=parser.parse_args(); binary=str(Path(args.binary).resolve())
    nodes,elements=read_model(args.model)
    keys=['slabs_count','walls_count','inclined_panels_count','columns_count','beams_count','braces_count']
    with tempfile.TemporaryDirectory(prefix='topo-invariance-') as directory:
        root=Path(directory)
        def run(path):
            subprocess.run([binary,str(path),'--json',str(root/'out.json'),'--dxf',str(root/'out.dxf')],check=True,capture_output=True)
            return json.loads((root/'out.json').read_text())
        baseline=run(Path(args.model).resolve())
        repeated=run(Path(args.model).resolve())
        assert baseline==repeated,'Output is not deterministic between runs'
        print('Repeated run: identical JSON')
        def write(name,coordinate_transform,node_order,element_order):
            idmap={old:new for new,old in enumerate(node_order,1)}
            block4='(4/'+''.join(' '.join(format(x,'.15g') for x in coordinate_transform(nodes[n]))+'/' for n in node_order)+')'
            block1='(1/'+''.join(' '.join(map(str,[*elements[e][:2],*[idmap[n] for n in elements[e][2:]]]))+'/' for e in element_order)+')'
            path=root/name;path.write_text(block4+'\n'+block1);return path
        angle=np.deg2rad(37.0)
        rotation=np.array([[np.cos(angle),-np.sin(angle),0.],[np.sin(angle),np.cos(angle),0.],[0.,0.,1.]])
        cases=[('translated',lambda p:p+np.array([1000.,-500.,10.]),list(nodes),list(elements)),
               ('rotated',lambda p:np.array([-p[1],p[0],p[2]]),list(nodes),list(elements)),
               ('rotated_37',lambda p:rotation@p,list(nodes),list(elements)),
               ('renumbered',lambda p:p,list(reversed(nodes)),list(reversed(elements)))]
        for name,transform,norder,eorder in cases:
            report=run(write(name+'.txt',transform,norder,eorder))
            summary={k:report[k] for k in keys}
            assert summary=={k:baseline[k] for k in keys},(name,summary,{k:baseline[k] for k in keys})
            original_groups=sorted(tuple(sorted(p['source_element_ids'])) for p in baseline['panels'])
            restored_groups=sorted(tuple(sorted(eorder[i-1] for i in p['source_element_ids'])) for p in report['panels'])
            assert restored_groups==original_groups,(name,'panel source membership changed')
            assert sum(p['filled_holes'] for p in report['panels'])==sum(p['filled_holes'] for p in baseline['panels'])
            print(name+': counts, panel membership and filled holes preserved')

if __name__=='__main__': main()
