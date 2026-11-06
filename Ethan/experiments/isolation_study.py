"""Run: python -m experiments.isolation_study --model plummer [--backend agama]."""
import argparse
import json
from pathlib import Path
from functools import partial
from nbody.galaxies.sources.native import plummer
from nbody.galaxies.generate import galaxy
from nbody.gravity.fast_direct import acceleration
from nbody.integrate import evolve
from nbody.io import save,save_trajectory
from nbody.diagnostics import inspect_isolation


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--model',choices=('plummer','disk'),default='plummer')
    ap.add_argument('--backend',choices=('preview','agama'),default='preview')
    ap.add_argument('--particles',type=int,default=500)
    ap.add_argument('--steps',type=int,default=100)
    ap.add_argument('--dt',type=float,default=.002)
    ap.add_argument('--eps',type=float,default=.0125)
    ap.add_argument('--seed',type=int,default=42)
    ap.add_argument('--out',default='results/isolation_demo')
    a=ap.parse_args();out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    if a.model=='plummer':
        p=plummer(a.particles,mass=1.25,a=.2,seed=a.seed);backend='analytic Plummer DF'
    else:
        cfg=json.loads((Path(__file__).resolve().parents[1]/'presets/mice.json').read_text())
        cfg['seed']=a.seed
        cfg['particles']={k:max(10,round(a.particles*w)) for k,w in [('disk',.55),('bulge',.10),('halo',.35)]}
        p=galaxy(cfg,a.backend);backend=a.backend
    save(out/'initial_galaxy.npz',p,dict(backend=backend,model=a.model))
    force=partial(acceleration,eps=a.eps)
    snaps=list(evolve(p,force,a.dt,a.steps,max(1,a.steps//20)))
    save_trajectory(out/'isolation_trajectory.npz',snaps,dict(backend=backend,eps=a.eps,dt=a.dt))
    report,_=inspect_isolation(snaps,a.eps,out,backend)
    print(json.dumps(report,indent=2))
    print('Saved to',out.resolve())

if __name__=='__main__':main()
