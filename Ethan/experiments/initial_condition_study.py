"""Native vs globally virialized disk candidates: a repeatable V1.6 demonstration.

Run: python -m experiments.initial_condition_study
"""
import json
from pathlib import Path
from nbody.galaxies.generate import galaxy
from nbody.io import save
from nbody.initial_conditions import run_relax_restart

root=Path(__file__).resolve().parents[1]
cfg=json.loads((root/'presets/mice.json').read_text())
cfg['particles']=dict(disk=240,bulge=40,halo=160)
cfg['seed']=42
p=galaxy(cfg,'preview')
out=root/'results'/'v16_initial_conditions_demo'
out.mkdir(parents=True,exist_ok=True)
save(out/'native_initial.npz',p,dict(backend='native preview',config=cfg))

for label,target in [('native',None),('globally_virialized',1.0)]:
    result=run_relax_restart(p,out/label,eps=.0125,dt=.004,
                             relax_steps=60,restart_steps=60,
                             save_every=10,target_virial=target,
                             source='same initial native disk-bulge-halo snapshot')
    print(label,'status:',result['status'],
          'Qvir:',result['virial_ratio'],
          'z_rms:',{k:round(v.get('z_rms',0),5)
                   for k,v in result['structure'].items()})
print('Results:',out)
print('These are exploratory native models, not validated Mice progenitors.')
