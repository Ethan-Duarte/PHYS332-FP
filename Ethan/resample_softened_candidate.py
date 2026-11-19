"""Change particle resolution while holding the underlying softened DF fixed.

Reconstructs the DF reference field; does not repeat the SCM iterations.
Unchanged components are copied exactly from the supplied source snapshot.
"""
import argparse
import json
import shutil
from datetime import datetime
from pathlib import Path
import numpy as np
import agama
from agama_softened_control import make_distribution_functions, run_screen
from nbody.galaxies.softening import soften_expansion
from nbody.particles import Particles, combine
from nbody.io import load, save


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('--disk', type=int)
    parser.add_argument('--bulge', type=int)
    parser.add_argument('--halo', type=int)
    parser.add_argument('--seed', type=int)
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--steps', type=int, default=128)
    parser.add_argument('--dt', type=float, default=1/1024)
    parser.add_argument('--build-only', action='store_true')
    args = parser.parse_args()
    source, meta = load(args.source)
    if 'numerical_force_checks' not in meta or 'quadrature_order' not in meta:
        raise SystemExit('Supply a candidate from agama_softened_control.py')
    cfg = meta['config']
    for index,name in enumerate(('disk','bulge','halo')):
        if getattr(args,name) is None:
            setattr(args,name,int(np.sum(source.component==index)))
    if args.seed is None:
        args.seed = cfg['seed']
    if min(args.disk,args.bulge,args.halo) < 1000 or args.steps<1 or args.dt<=0 or args.threads<1:
        parser.error('>=1000 particles/component; positive steps, dt, threads required')
    d,b,h = [cfg[k] for k in ('disk','bulge','halo')]
    rd,zd = float(d['scale']),float(d['height'])
    args.epsilon = float(meta['epsilon'])
    agama.setNumThreads(args.threads)
    sg = dict(rminSph=min(b['scale'],h['scale'],rd)/40,
              rmaxSph=max(2*h.get('cutoff',0),30*h['scale'],15*b['scale'],15*rd),
              sizeRadialSph=50,lmaxAngularSph=6)
    dg = dict(RminCyl=rd/30,RmaxCyl=12*rd,sizeRadialCyl=30,
              zminCyl=zd/4,zmaxCyl=10*rd,sizeVerticalCyl=26)
    density = [agama.Density(type='Disk',mass=d['mass'],scaleRadius=rd,scaleHeight=-zd/2),
               agama.Density(type='Spheroid',mass=b['mass'],scaleRadius=b['scale'],alpha=1,beta=4,gamma=1),
               agama.Density(type='Spheroid',mass=h['mass'],scaleRadius=h['scale'],alpha=1,beta=4,gamma=0)]
    initial = agama.SelfConsistentModel(**sg,**dg)
    initial.components = [agama.Component(density=rho,disklike=i==0) for i,rho in enumerate(density)]
    print('Reconstructing the fixed DF reference field (no SCM loop)...',flush=True)
    initial.iterate()
    reference = soften_expansion(initial.potential,args.epsilon,sg,dg,order=tuple(meta['quadrature_order']))
    dfs = make_distribution_functions(agama,reference,density,cfg)
    potential_file = args.source.parent/'smooth_potential.ini'
    final = agama.Potential(str(potential_file))
    af = agama.ActionFinder(final)
    components = []
    reused = []
    for index,count in enumerate((args.disk,args.bulge,args.halo)):
        mask = source.component==index
        mass = float(source.m[mask].sum())
        if abs(float(dfs[index].totalMass())/mass-1) > 1e-5:
            raise SystemExit('Reconstructed DF mass differs from source; stop rather than silently changing model')
        if count==int(mask.sum()) and args.seed==cfg['seed']:
            print(f'Component {index}: copying {count} original particles exactly',flush=True)
            components.append(Particles(source.x[mask].copy(),source.v[mask].copy(),source.m[mask].copy(),source.component[mask].copy(),source.galaxy[mask].copy()))
            reused.append(index)
        else:
            print(f'Component {index}: sampling {count} particles',flush=True)
            seed = args.seed if index==0 else args.seed+1 if index==1 else args.seed+100+count
            agama.setRandomSeed(seed)
            xv,m = agama.GalaxyModel(potential=final,df=dfs[index],af=af).sample(count)
            m = np.asarray(m,float);m *= mass/m.sum()
            components.append(Particles(np.array(xv[:,:3]),np.array(xv[:,3:]),m,np.full(count,index),np.zeros(count,dtype=int)))
    output = args.source.parent.parent/('agama_resampled_'+datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
    output.mkdir()
    shutil.copyfile(potential_file,output/'smooth_potential.ini')
    updated = dict(meta, config=dict(cfg,seed=args.seed,particles=dict(disk=args.disk,bulge=args.bulge,halo=args.halo)),
                   source_candidate=str(args.source),unchanged_components_copied=reused,
                   status='Resolution test candidate; live stability unverified')
    save(output/'candidate.npz',combine(*components),updated)
    print(f'Saved {output}/candidate.npz',flush=True)
    if not args.build_only:
        run_screen(output/'candidate.npz',args)


if __name__=='__main__':
    main()
