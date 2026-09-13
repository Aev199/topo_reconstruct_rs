#!/usr/bin/env python3
"""Early M4 meshability gate for reconstruction-v2 JSON; never repairs geometry."""
import argparse, json, math
from collections import Counter, defaultdict, deque
from pathlib import Path


def sub(a,b): return tuple(x-y for x,y in zip(a,b))
def dot(a,b): return sum(x*y for x,y in zip(a,b))
def cross(a,b): return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def norm(a): return math.sqrt(dot(a,a))
def turn(a,b,c): return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
def area(poly): return .5*sum(poly[i][0]*poly[(i+1)%len(poly)][1]-poly[(i+1)%len(poly)][0]*poly[i][1] for i in range(len(poly)))
def inside(p,a,b,c,e):
    q=(turn(a,b,p),turn(b,c,p),turn(c,a,p)); return not(min(q)<-e and max(q)>e)

def ring(surface,edges):
    return [edges[u['edge']][1 if u['reversed'] else 0] for u in surface['boundaries'][0]]

def earclip(ids,uv,e):
    ids=list(ids); a=area([uv[v] for v in ids])
    if len(ids)<3 or len(set(ids))!=len(ids) or abs(a)<=e*e: raise ValueError('invalid_ring')
    ccw=a>0; out=[]; guard=0
    while len(ids)>3:
        guard+=1
        if guard>len(ids)*len(ids)+100: raise ValueError('ear_clipping_stalled')
        found=False
        for i,b in enumerate(ids):
            x,z=ids[i-1],ids[(i+1)%len(ids)]; t=turn(uv[x],uv[b],uv[z])
            if (ccw and t<=e) or ((not ccw) and t>=-e): continue
            if any(inside(uv[p],uv[x],uv[b],uv[z],e) for p in ids if p not in (x,b,z)): continue
            out.append((x,b,z) if ccw else (x,z,b)); del ids[i]; found=True; break
        if not found: raise ValueError('no_ear_found')
    x,y,z=ids; out.append((x,y,z) if ccw else (x,z,y)); return out

def project(pl,p):
    d=sub(p,pl['origin']); return (dot(d,pl['u']),dot(d,pl['v']))

def insert_point(tris,v,uv,e):
    hits=[i for i,t in enumerate(tris) if inside(uv[v],*(uv[x] for x in t),e)]
    if not hits: raise ValueError('point_contact_outside_mesh')
    out=[]
    for i,t in enumerate(tris):
        if i not in hits: out.append(t); continue
        a,b,c=t; p=uv[v]
        on=None
        for x,y in ((a,b),(b,c),(c,a)):
            if abs(turn(uv[x],uv[y],p))<=e: on=(x,y,next(q for q in t if q not in (x,y))); break
        if on:
            x,y,o=on
            if v not in (x,y): out.extend([(x,v,o),(v,y,o)])
        else: out.extend([(a,b,v),(b,c,v),(c,a,v)])
    return [t for t in out if len(set(t))==3]

def tri_metrics(t,verts):
    p=[verts[i] for i in t]; ls=[norm(sub(p[(i+1)%3],p[i])) for i in range(3)]
    ar=.5*norm(cross(sub(p[1],p[0]),sub(p[2],p[0])))
    angles=[]
    for i in range(3):
        u,v=sub(p[i-1],p[i]),sub(p[(i+1)%3],p[i]); c=max(-1,min(1,dot(u,v)/(norm(u)*norm(v))))
        angles.append(math.degrees(math.acos(c)))
    return ar,min(angles),max(ls)/min(ls)

def blockers(data):
    top=data['topology']; m=top['preview']; bars=top.get('axis_assembly',{}); src=top['vertex_source_nodes']; out=defaultdict(set)
    boundary={s:{src[v] for v in ring(x,m['edges'])} for s,x in enumerate(m['surfaces'])}
    for s,x in enumerate(m['surfaces']):
        if len(x['contours'])!=1: out[s].add('holes_require_constrained_mesher')
    for c in bars.get('contacts',[]):
        if c['kind']=='interval' and c['location']=='interior': out[c['surface']].add('interior_axis_interval_requires_constrained_edge')
    bad=[set(x.get('source_nodes',[])) for x in bars.get('issues',[])]
    for s,nodes in boundary.items():
        if any(nodes & b for b in bad): out[s].add('rejected_axis_touches_surface_boundary')
    return {k:sorted(v) for k,v in out.items()}

def components(ids,m):
    ids=set(ids); own=defaultdict(list)
    for s in ids:
        for r in m['surfaces'][s]['boundaries']:
            for u in r: own[u['edge']].append(s)
    g={s:set() for s in ids}
    for o in own.values():
        for a in o: g[a].update(x for x in o if x!=a)
    out=[]
    while ids:
        seed=min(ids); ids.remove(seed); q=deque([seed]); c=[]
        while q:
            a=q.popleft(); c.append(a)
            for b in sorted(g[a]):
                if b in ids: ids.remove(b); q.append(b)
        out.append(sorted(c))
    return out

def choose(data,n):
    m=data['topology']['preview']; bl=blockers(data); ok=[i for i in range(len(m['surfaces'])) if i not in bl]; cs=components(ok,m)
    if not cs: return [],bl
    def score(c):
        cnt=Counter(u['edge'] for s in c for r in m['surfaces'][s]['boundaries'] for u in r)
        return (sum(v>1 for v in cnt.values()),len(c),-min(c))
    c=max(cs,key=score)
    if len(c)<=n:return c,bl
    own=defaultdict(list)
    for s in c:
        for r in m['surfaces'][s]['boundaries']:
            for u in r: own[u['edge']].append(s)
    g={s:set() for s in c}
    for o in own.values():
        for a in o:g[a].update(x for x in o if x!=a)
    seed=max(c,key=lambda s:(len(g[s]),-s)); q=deque([seed]); seen={seed}; out=[]
    while q and len(out)<n:
        a=q.popleft();out.append(a)
        for b in sorted(g[a],key=lambda x:(-len(g[x]),x)):
            if b not in seen:seen.add(b);q.append(b)
    return sorted(out),bl

def build(data,max_surfaces=8):
    top=data['topology'];m=top['preview'];bars=top.get('axis_assembly',{});e=top['policy']['precision']; sel,bl=choose(data,max_surfaces)
    r={'format':'topo-reconstruct-m4-probe-v1','meshable':False,'selected_surfaces':sel,'selected_axes':[],'triangles':[],'metrics':{},'blockers':[],'global_surface_blockers':{str(k):v for k,v in sorted(bl.items())}}
    if not sel:r['blockers'].append('no_eligible_connected_fragment');return r
    by=defaultdict(list)
    for c in bars.get('contacts',[]):by[c['surface']].append(c)
    alltris=[]; axes=set()
    try:
        for s in sel:
            sf=m['surfaces'][s];pl=m['planes'][sf['plane']];ids=ring(sf,m['edges']);uv={v:tuple(sf['contours'][0][i]) for i,v in enumerate(ids)};ts=earclip(ids,uv,e)
            for c in by.get(s,[]):
                axes.add(c['axis'])
                if c['kind']=='point' and c['location']=='interior':
                    v=c['vertex'];uv.setdefault(v,project(pl,m['vertices'][v]));ts=insert_point(ts,v,uv,e)
            for a,b,c in ts:
                n=cross(sub(m['vertices'][b],m['vertices'][a]),sub(m['vertices'][c],m['vertices'][a]))
                if dot(n,pl['normal'])<0:b,c=c,b
                alltris.append({'surface':s,'vertices':[a,b,c]})
    except (KeyError,ValueError,ZeroDivisionError) as x:r['blockers'].append('triangulation_failed:'+str(x));return r
    for s in sel:
        cnt=Counter(tuple(sorted((t['vertices'][i],t['vertices'][(i+1)%3]))) for t in alltris if t['surface']==s for i in range(3))
        boundary={tuple(sorted(m['edges'][u['edge']])) for u in m['surfaces'][s]['boundaries'][0]}
        if any(cnt[x]!=1 for x in boundary):r['blockers'].append(f'boundary_edge_multiplicity:{s}')
        if any(v!=2 for k,v in cnt.items() if k not in boundary):r['blockers'].append(f'internal_edge_multiplicity:{s}')
    qs=[tri_metrics(tuple(t['vertices']),m['vertices']) for t in alltris]
    r['triangles']=alltris;r['selected_axes']=sorted(axes);r['metrics']={'surface_count':len(sel),'axis_count':len(axes),'vertex_count':len({v for t in alltris for v in t['vertices']}),'triangle_count':len(alltris),'minimum_triangle_area':min((x[0] for x in qs),default=0),'minimum_triangle_angle_deg':min((x[1] for x in qs),default=0),'maximum_triangle_edge_ratio':max((x[2] for x in qs),default=0)}
    r['meshable']=bool(alltris) and not r['blockers'];return r

def write_obj(path,data,r):
    m=data['topology']['preview']; used=sorted({v for t in r['triangles'] for v in t['vertices']}); mp={v:i+1 for i,v in enumerate(used)};lines=['# M4 probe']
    lines += ['v %.12g %.12g %.12g'%tuple(m['vertices'][v]) for v in used]
    cur=None
    for t in r['triangles']:
        if t['surface']!=cur:cur=t['surface'];lines.append('g surface_%s'%cur)
        lines.append('f '+' '.join(str(mp[v]) for v in t['vertices']))
    Path(path).write_text('\n'.join(lines)+'\n',encoding='utf-8')

def selftest():
    verts=[[0,0,0],[1,0,0],[1,1,0],[0,1,0],[2,0,0],[2,1,0],[.5,.5,0]];edges=[[0,1],[1,2],[2,3],[0,3],[1,4],[4,5],[2,5]]
    sf=[{'plane':0,'contours':[[[0,0],[1,0],[1,1],[0,1]]],'boundaries':[[{'edge':0,'reversed':False},{'edge':1,'reversed':False},{'edge':2,'reversed':False},{'edge':3,'reversed':True}]]},{'plane':0,'contours':[[[1,0],[2,0],[2,1],[1,1]]],'boundaries':[[{'edge':4,'reversed':False},{'edge':5,'reversed':False},{'edge':6,'reversed':True},{'edge':1,'reversed':True}]]}]
    d={'topology':{'policy':{'precision':1e-8},'vertex_source_nodes':[1,2,3,4,5,6,7],'preview':{'vertices':verts,'edges':edges,'planes':[{'origin':[0,0,0],'normal':[0,0,1],'u':[1,0,0],'v':[0,1,0]}],'surfaces':sf},'axis_assembly':{'issues':[],'contacts':[{'kind':'point','axis':0,'surface':0,'vertex':6,'location':'interior'}]}}}
    r=build(d);assert r['meshable'] and r['metrics']['surface_count']==2 and r['metrics']['triangle_count']==6 and any(6 in t['vertices'] for t in r['triangles'])

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('report',nargs='?');p.add_argument('--max-surfaces',type=int,default=8);p.add_argument('--output');p.add_argument('--obj');p.add_argument('--self-test',action='store_true');a=p.parse_args()
    if a.self_test:selftest();print('m4_mesh_probe self-test: OK');return 0
    if not a.report:p.error('report required unless --self-test')
    data=json.load(open(a.report,encoding='utf-8'));r=build(data,a.max_surfaces);text=json.dumps(r,ensure_ascii=False,indent=2)
    Path(a.output).write_text(text+'\n',encoding='utf-8') if a.output else print(text)
    if a.obj and r['triangles']:write_obj(a.obj,data,r)
    return 0 if r['meshable'] else 2
if __name__=='__main__':raise SystemExit(main())
