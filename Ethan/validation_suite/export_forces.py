"""Adapter for native project solvers. Does not claim to evaluate Zeno forces."""
import argparse,json,sys
from pathlib import Path
import numpy as np
from core import load,law,digest

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('trajectory',type=Path)
    p.add_argument('--repo-root',type=Path,default=Path('.'));p.add_argument('--solver',choices=['direct','tree'],required=True)
    p.add_argument('--frame',type=int,default=0);p.add_argument('--theta',type=float,default=.5);p.add_argument('--samples',type=int,default=96)
    p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    if a.out.exists():p.error('Output already exists')
    r=load(a.trajectory,snapshot=True);eps,G,kernel=law(r)
    if not 0<=a.frame<len(r.t) or a.samples<1 or a.theta<=0:p.error('Invalid frame/sample count/theta')
    sys.path.insert(0,str(a.repo_root.resolve()))
    if a.solver=='direct':
        from nbody.gravity.parallel_direct import acceleration
        import nbody.gravity.parallel_direct as module
        acc=acceleration(r.x[a.frame],r.m,eps=eps,G=G)
    else:
        from nbody.gravity.fast_barnes_hut import acceleration
        import nbody.gravity.fast_barnes_hut as module
        acc=acceleration(r.x[a.frame],r.m,eps=eps,theta=a.theta,G=G)
    # Radial stratification per component/galaxy plus random coverage.
    rng=np.random.default_rng(1467);ids=[]
    groups=[np.flatnonzero((r.g==g)&(r.c==c)) for g in np.unique(r.g) for c in np.unique(r.c)]
    for group in groups:
        if not len(group):continue
        x=r.x[a.frame,group];radius=np.linalg.norm(x-np.average(x,axis=0,weights=r.m[group]),axis=1)
        order=group[np.argsort(radius)];n=max(3,a.samples//max(1,len(groups)))
        ids.extend(order[np.linspace(0,len(order)-1,min(n,len(order))).astype(int)].tolist())
    ids=np.unique(np.r_[ids,rng.choice(len(r.m),min(a.samples,len(r.m)),replace=False)]).astype(int)
    meta=dict(trajectory_sha256=digest(a.trajectory),frame=a.frame,epsilon=eps,G=G,kernel=kernel,solver=a.solver,theta=a.theta if a.solver=='tree' else None,solver_source_sha256=digest(module.__file__))
    np.savez_compressed(a.out,indices=ids,a=acc[ids],metadata=json.dumps(meta));print('Saved',a.out,'targets',len(ids))

if __name__=='__main__':main()
