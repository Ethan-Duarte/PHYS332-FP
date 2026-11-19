"""Diagnostic only: live disk in the candidate's fixed smooth bulge+halo field.

Removes spheroid discreteness AND its dynamical response; this is not a
production merger model or a unique measurement of halo collisional heating.
"""
import argparse
import json
from datetime import datetime
from pathlib import Path
import numpy as np
import numba
import agama
from agama_softened_control import metrics
from nbody.io import load,save_trajectory
from nbody.particles import Particles
from nbody.integrate import evolve
from nbody.gravity.parallel_direct import acceleration


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source',type=Path)
    parser.add_argument('--steps',type=int,default=512)
    parser.add_argument('--dt',type=float,default=1/1024)
    parser.add_argument('--threads',type=int,default=4)
    args = parser.parse_args()
    if args.steps<1 or args.dt<=0 or args.threads<1:parser.error('positive steps,dt,threads required')
    p,meta = load(args.source)
    if 'numerical_force_checks' not in meta:raise SystemExit('Use the softened candidate')
    numba.set_num_threads(min(args.threads,numba.config.NUMBA_NUM_THREADS))
    agama.setNumThreads(args.threads)
    potential = agama.Potential(str(args.source.parent/'smooth_potential.ini'))
    if len(potential)!=2:raise SystemExit('Expected spheroid and disk potential components')
    spheroid = potential[0]
    disk = p.component==0
    disk = Particles(p.x[disk].copy(),p.v[disk].copy(),p.m[disk].copy(),p.component[disk].copy(),p.galaxy[disk].copy())
    epsilon=float(meta['epsilon'])
    def force(x,m):return acceleration(x,m,epsilon)+spheroid.force(x)
    out=args.source.parent/('rigid_spheroids_'+datetime.now().strftime('%Y%m%d_%H%M%S'))
    out.mkdir()
    rows,frames=[],[]
    every=max(8,8*int(np.ceil(args.steps/256)))
    for t,state in evolve(disk,force,args.dt,args.steps,every):
        rows.append(dict(time=float(t),**metrics(state)))
        frames.append((t,state))
        temp=out/'checkpoint.tmp.npz'
        save_trajectory(temp,frames,dict(status='diagnostic only; static spheroids',epsilon=epsilon,dt=args.dt,source=str(args.source)))
        temp.replace(out/'disk_trajectory.npz')
        (out/'metrics.json').write_text(json.dumps(rows,indent=2))
        print(f't={t:.6f}/{args.dt*args.steps:.6f}',flush=True)
    summary=dict(relative_changes={k:rows[-1][k]/rows[0][k]-1 for k in ('height','sigma_z','r50')},duration=args.dt*args.steps,
                 caveat='Removes both spheroid granularity and response; cannot distinguish them alone; not a merger model')
    (out/'summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2),flush=True)
    print(f'Results: {out}',flush=True)


if __name__=='__main__':main()
