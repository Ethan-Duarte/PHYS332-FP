"""Prepare matched low-resolution homologous size trials; no evolution launched."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import numpy as np
from nbody.io import load,save
from nbody.particles import Particles


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('source',type=Path)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--scales',nargs='+',type=float,default=[1.,1.25,1.5])
    ap.add_argument('--seed',type=int,default=1234)
    ap.add_argument('--disk',type=int,default=6144)
    ap.add_argument('--bulge',type=int,default=2048)
    ap.add_argument('--halo',type=int,default=8192)
    args=ap.parse_args()
    if not all(np.isfinite(s) and s>0 for s in args.scales):ap.error('Positive finite scales required')
    p,original=load(args.source)
    if original.get('model_id')!='zeno_reference_bdh_v1':ap.error('Use the original full-resolution jeans_candidate.npz')
    counts={0:args.disk,1:args.bulge,2:args.halo}
    if min(counts.values())<1:ap.error('Component counts must be positive')
    if len(p.m)!=65536 or any(np.sum(p.component==k)<n for k,n in counts.items()):ap.error('Expected the full 65536-particle reference galaxy')
    if not np.isfinite(p.x).all() or not np.isfinite(p.v).all() or not np.isfinite(p.m).all() or np.any(p.m<=0):ap.error('Invalid source')
    rng=np.random.default_rng(args.seed);ids=[];mass=[]
    for k,n in counts.items():
        pool=np.flatnonzero(p.component==k)
        chosen=rng.choice(pool,n,replace=False)
        ids.extend(chosen);selected=p.m[chosen];mass.extend(selected*p.m[pool].sum()/selected.sum())
    ids=np.array(ids);mass=np.array(mass)
    disk=p.component==0
    origin=np.average(p.x[disk],axis=0,weights=p.m[disk]);bulk=np.average(p.v[disk],axis=0,weights=p.m[disk])
    x=p.x[ids]-origin;v=p.v[ids]-bulk
    args.out.mkdir(parents=True,exist_ok=False)
    source_hash=hashlib.sha256(args.source.read_bytes()).hexdigest()
    manifest=[]
    for scale in args.scales:
        md=deepcopy(original)
        md['parent_model_id']=original['model_id'];md['model_id']='zeno_reference_bdh_size_screen_v1'
        md['epsilon']=float(original['epsilon']*scale)
        md['construction_force_model']=dict(kernel='plummer',epsilon=md['epsilon'],G=1.)
        md['equilibrium_status']='unvalidated reduced-particle geometry screen'
        md['size_screen']=dict(scale=scale,seed=args.seed,source_sha256=source_hash,
                               source_model=deepcopy(original),counts={str(k):n for k,n in counts.items()},
                               transformation='All components: x -> s*x; v -> v/sqrt(s); epsilon -> s*epsilon; component masses unchanged.',
                               internal_time_scale_factor=float(scale**1.5),
                               sampling='Same selected particle IDs at every scale; reduced counts increase finite-particle noise.',
                               interpretation='Whole-galaxy size relative to fixed encounter orbit; changes halo and bulge as well as disk. Not a disk-only experiment or preserved-equilibrium certification.')
        name='candidate_s'+f'{scale:.2f}'.replace('.','p')+'.npz'
        save(args.out/name,Particles(x*scale,v/np.sqrt(scale),mass,p.component[ids],np.zeros(len(ids),dtype=int)),md)
        manifest.append(dict(scale=scale,file=name,particles_per_galaxy=len(ids),pair_particles=2*len(ids),epsilon=md['epsilon']))
        print(f'Prepared {name}: {len(ids)} particles/galaxy, size x{scale:g}, epsilon {md["epsilon"]:g}',flush=True)
    np.save(args.out/'selected_source_indices.npy',ids)
    (args.out/'manifest.json').write_text(json.dumps(manifest,indent=2))
    print('No evolution launched. Compare a matched-resolution control before interpreting parameter changes.',flush=True)

if __name__=='__main__':main()
