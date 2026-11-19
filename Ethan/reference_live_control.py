"""Live reference disk in its smooth fixed spheroid field, without AGAMA.

This removes BOTH spheroid graininess and response. It is an isolation
diagnostic, not a merger model or a unique measurement of halo scattering.
"""
import argparse
import csv
import hashlib
import json
import time
from functools import partial
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicSpline
from nbody.io import load, save_trajectory
from nbody.particles import Particles
from nbody.integrate import evolve
from nbody.diagnostics import structure
from nbody.galaxies.provenance import check_evolution_softening


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('candidate',type=Path)
    parser.add_argument('spheroid_field',type=Path)
    parser.add_argument('--solver',choices=['direct','tree'],default='direct')
    parser.add_argument('--theta',type=float,default=.5)
    parser.add_argument('--threads',type=int,default=4)
    parser.add_argument('--dt',type=float,default=1/512)
    parser.add_argument('--steps',type=int,default=64)
    parser.add_argument('--save-every',type=int,default=16)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    if min(args.steps,args.threads,args.save_every)<1 or not np.isfinite(args.dt) or args.dt<=0 or not np.isfinite(args.theta) or args.theta<=0:
        parser.error('steps, threads, save-every, dt and theta must be positive')
    particles,meta=load(args.candidate)
    if not np.isfinite(particles.x).all() or not np.isfinite(particles.v).all() or not np.isfinite(particles.m).all() or np.any(particles.m<=0):
        raise SystemExit('Invalid particle data')
    with np.load(args.spheroid_field) as field:
        eps=float(field['epsilon']);r=field['radius'];phi=field['potential']
        if float(field['G'])!=1:raise SystemExit('Expected G=1')
        if field['model_id'].item()!=meta.get('model_id'):
            raise SystemExit('Candidate and smooth field belong to different models')
    check_evolution_softening(meta,eps)
    if meta.get('backend')!='Zeno reference, experimental moment adjustment':
        raise SystemExit('Use the reference candidate produced by zeno_smooth_control.py --vertical-jeans')
    disk=particles.component==0
    p=Particles(particles.x[disk].copy(),particles.v[disk].copy(),particles.m[disk].copy(),
                particles.component[disk].copy(),particles.galaxy[disk].copy())
    potential=CubicSpline(r,phi,bc_type=((1,0),'not-a-knot'))
    if args.solver=='direct':
        import numba
        numba.set_num_threads(min(args.threads,numba.config.NUMBA_NUM_THREADS))
        from nbody.gravity.parallel_direct import acceleration
        self_force=partial(acceleration,eps=eps)
    else:
        from nbody.gravity.fast_barnes_hut import acceleration
        self_force=partial(acceleration,eps=eps,theta=args.theta)

    def force(x,m):
        radius=np.linalg.norm(x,axis=1)
        if np.any(radius>r[-1]):raise RuntimeError('Particle outside external field; extrapolation forbidden')
        radial=-potential(radius,1)
        external=x*np.divide(radial,radius,out=np.zeros_like(radius),where=radius>0)[:,None]
        return self_force(x,m)+external

    print('Warming up the selected solver...',flush=True)
    force(p.x,p.m)
    started=time.perf_counter();force(p.x,p.m)
    seconds=time.perf_counter()-started
    print(f'Complete force evaluation: {seconds:.3f} s; approximate integration cost: {seconds*(args.steps+1)/60:.1f} minutes, excluding output.',flush=True)
    args.out.mkdir(parents=True,exist_ok=False)
    metadata=dict(setup='isolated live reference disk; smooth fixed spheroids',
        source=str(args.candidate),source_sha256=hashlib.sha256(args.candidate.read_bytes()).hexdigest(),
        epsilon=eps,dt=args.dt,steps=args.steps,solver=args.solver,theta=args.theta,
        equilibrium_status='unvalidated',caveat=__doc__)
    (args.out/'config.json').write_text(json.dumps(metadata,indent=2))
    rows=[];frames=[]
    for t,state in evolve(p,force,args.dt,args.steps,args.save_every):
        rows.append(dict(time=float(t),**structure(state)))
        frames.append((t,state))
        temporary=args.out/'checkpoint.tmp.npz'
        save_trajectory(temporary,frames,metadata)
        temporary.replace(args.out/'disk_trajectory.npz')
        with (args.out/'metrics.csv').open('w',newline='') as handle:
            writer=csv.DictWriter(handle,fieldnames=list(rows[0]))
            writer.writeheader();writer.writerows(rows)
        print(f't={t:.5f}/{args.steps*args.dt:.5f}',flush=True)
    summary=dict(relative_changes={k:rows[-1][k]/rows[0][k]-1 for k in ['z_rms','sigma_z','r50']},
                 duration=rows[-1]['time'],caveat=__doc__)
    (args.out/'summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':main()
