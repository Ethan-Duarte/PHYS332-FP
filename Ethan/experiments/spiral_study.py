"""Run: python -m experiments.spiral_study --steps 400 --particles 1600

Or reanalyze an isolated disk trajectory without another gravity run:
python -m experiments.spiral_study --trajectory results/runs/NAME/run_trajectory.npz
"""
import argparse
import json
from functools import partial
from pathlib import Path
import numpy as np

from nbody.galaxies.generate import galaxy
from nbody.gravity.fast_direct import acceleration
from nbody.integrate import evolve
from nbody.io import load_trajectory, save_trajectory, save
from nbody.particles import Particles
from nbody.spirals import spiral_study


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--trajectory',type=Path)
    parser.add_argument('--particles',type=int,default=1200,help='Total live particles (disk, bulge, halo)')
    parser.add_argument('--steps',type=int,default=400)
    parser.add_argument('--dt',type=float,default=.003)
    parser.add_argument('--eps',type=float,default=.0125)
    parser.add_argument('--save-every',type=int,default=10)
    parser.add_argument('--seed',type=int,default=42)
    parser.add_argument('--out',type=Path,default=Path('results/spiral_study'))
    a=parser.parse_args()
    if a.trajectory:
        t,x,v,m,comp,ids,meta=load_trajectory(a.trajectory)
        if meta.get('setup')!='Isolated galaxy' or not np.any(comp==0):
            raise ValueError('This study requires an isolated disk–bulge–halo trajectory.')
        snapshots=[(float(time),Particles(px,pv,m,comp,ids)) for time,px,pv in zip(t,x,v)]
        rd=meta['config']['disk']['scale']
        eps=meta['eps']
        out=a.out if a.out!=Path('results/spiral_study') else a.trajectory.parent
    else:
        config=json.loads((Path(__file__).resolve().parents[1]/'presets/mice.json').read_text())
        config['seed']=a.seed
        nd=int(a.particles*.6);nb=int(a.particles*.1);nh=a.particles-nd-nb
        config['particles']=dict(disk=nd,bulge=nb,halo=nh)
        p=galaxy(config,'preview')
        a.out.mkdir(parents=True,exist_ok=True)
        save(a.out/'initial_galaxy.npz',p,dict(backend='preview',config=config))
        snapshots=list(evolve(p,partial(acceleration,eps=a.eps),a.dt,a.steps,a.save_every))
        save_trajectory(a.out/'isolation_trajectory.npz',snapshots,
                        dict(config=config,eps=a.eps,dt=a.dt,setup='Isolated galaxy',backend='preview',
                             caveat='Native preview disk velocities NOT in equilibrium'))
        rd=config['disk']['scale'];eps=a.eps;out=a.out
    report,plot=spiral_study(snapshots,rd,eps,out)
    print(json.dumps(report,indent=2))
    print('Saved',plot or out/'spiral_report.json')

if __name__=='__main__':
    main()
