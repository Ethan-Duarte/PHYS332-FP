"""Controlled Zeno reruns from an existing saved encounter initial frame.
Place beside zeno_mice_run.py in the project's Ethan directory.
"""
import argparse,json,hashlib,sys,time,subprocess,os
from pathlib import Path
import numpy as np
from zeno_mice_run import linux_path,linux_run,convert_frames
import zeno_mice_run
from nbody.particles import Particles

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('trajectory',type=Path,nargs='?');p.add_argument('--find-results',type=Path)
    p.add_argument('--mode',choices=['isolate','half-dt','tight-theta','baseline'],required=True)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--duration',type=float)
    p.add_argument('--dt',type=float,help='Optional explicit timestep, e.g. for isolation refinement')
    p.add_argument('--initial-particles',type=Path,help='Optional original double-precision initial_particles.npz')
    p.add_argument('--zeno-root',default='/home/whisky/nbody-tools/zeno-reference')
    a=p.parse_args()
    if a.trajectory is None:
        if a.find_results is None:p.error('Give trajectory or --find-results')
        matches=[]
        for f in a.find_results.rglob('run_trajectory.npz'):
            try:
                with np.load(f,allow_pickle=False) as d:
                    z=json.loads(d['metadata'].item())
                    if (len(d['m'])==16384 and z.get('setup')=='Mice geometry trial'
                        and z.get('source_sha256')=='a259ae858b7385a05eeaa2ffb83b5bdc55c95950180e502a0059695f3eeeafce'
                        and np.isclose(z.get('nominal_pericenter',-1),.2)
                        and np.isclose(z.get('dt',-1),1/512) and np.isclose(z.get('theta',-1),.7)):
                        matches.append(f)
            except (OSError,ValueError,KeyError):continue
        if not matches:p.error('No matching completed 16k size1.25 rp0.20 original trial found')
        a.trajectory=max(matches,key=lambda f:f.stat().st_mtime)
    source=a.trajectory.resolve()
    with np.load(source,allow_pickle=False) as d:
        old=json.loads(d['metadata'].item());original_dtype=str(d['x'].dtype);x=d['x'][0].astype(float);v=d['v'][0].astype(float)
        m=d['m'].copy();c=d['component'].copy();g=d['galaxy'].copy()
    original=a.initial_particles
    if original is None and source.name=='run_trajectory.npz':
        original=source.with_name('initial_particles.npz')
    if a.initial_particles is not None and not original.exists():p.error('Explicit initial-particles file does not exist')
    precision=original_dtype+' input initial frame loaded as float64'
    if original is not None and original.exists():
        with np.load(original,allow_pickle=False) as d:
            for k,expected in [('m',m),('component',c),('galaxy',g)]:
                if not np.array_equal(d[k],expected):p.error('Original snapshot labels/masses do not match replay')
            if not np.allclose(d['x'],x,rtol=2e-6,atol=2e-6) or not np.allclose(d['v'],v,rtol=2e-6,atol=2e-6):p.error('Initial snapshot does not match replay')
            x=d['x'].astype(float);v=d['v'].astype(float)
        precision='original initial_particles.npz; checked against replay'
    if a.mode=='isolate':
        keep=g==0;x=x[keep];v=v[keep];m=m[keep];c=c[keep];g=g[keep]
        x-=np.average(x,axis=0,weights=m);v-=np.average(v,axis=0,weights=m)
    dt=float(old['dt'])/(2 if a.mode=='half-dt' else 1)
    if a.dt is not None:
        if a.dt<=0:p.error('dt must be positive')
        dt=a.dt
    theta=.4 if a.mode=='tight-theta' else float(old['theta'])
    duration=float(a.duration if a.duration is not None else old['duration'])
    if duration<=0 or not np.isclose(duration/dt,round(duration/dt)):p.error('Duration must be positive integer multiple of dt')
    out=a.out.resolve();out.mkdir(parents=True,exist_ok=False)
    initial=Particles(x,v,m,c,g)
    meta=dict(old,validation_mode=a.mode,parent_trajectory_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
              input_precision=precision,dt=dt,theta=theta,duration=duration,
              note='Controlled numerical/isolated test; no automatic validation verdict.')
    with (out/'initial.txt').open('w') as f:
        f.write(str(len(m))+'\n');np.savetxt(f,np.column_stack([m,x,v]),fmt='%.17g')
    root=a.zeno_root.rstrip('/');work=linux_path(out)
    bridge=work+'/snapshot_bridge';bridge_source=linux_path(Path(zeno_mice_run.__file__).with_name('zeno_snapshot_bridge.c'))
    linux_run(['gcc','-std=gnu99','-DLINUX','-DDOUBLEPREC','-include','unistd.h',
        '-Wno-error=implicit-function-declaration','-I'+root+'/inc','-o',bridge,bridge_source,
        '-L'+root+'/lib','-lNBody','-lClib','-lgsl','-lgslcblas','-lm'])
    linux_run([bridge,'mode=import','in='+work+'/initial.txt','out='+work+'/initial.snap'])
    cmd=[root+'/bin/treecode_dp','in='+work+'/initial.snap','out='+work+'/evolution.snap',
        'log='+work+'/tree.log','save='+work+'/state_%d','eps='+str(old['eps']),
        'theta='+str(theta),'usequad='+str(old.get('usequad',True)).lower(),
        'dtime='+str(dt),'tstop='+str(duration),'dtout='+str(float(old['dt'])*old.get('save_every',8)),
        'outputs=Mass,Position,Velocity']
    sha_cmd=(['wsl','-d','Ubuntu','--'] if os.name=='nt' else [])+['sha256sum',root+'/bin/treecode_dp']
    meta['parent_solver_sha256']=old.get('solver_sha256');meta['solver_sha256']=subprocess.check_output(sha_cmd,text=True).split()[0]
    meta['tree_command']=cmd;(out/'config.json').write_text(json.dumps(meta,indent=2))
    print('Source:',source,'\nMode:',a.mode,'N:',len(m),'dt:',dt,'theta:',theta,'\nInput:',precision,flush=True)
    linux_run(cmd,progress=out/'tree.log')
    linux_run([bridge,'mode=export','in='+work+'/evolution.snap','out='+work+'/frames.bin'])
    times=convert_frames(out/'frames.bin',initial,meta,out/'run_trajectory.npz')
    if not np.isclose(times[-1],duration,atol=dt/8):raise RuntimeError('Incomplete run')
    if a.mode=='isolate':
        from nbody.diagnostics import structure
        with np.load(out/'run_trajectory.npz') as d:
            rows=[dict(time=float(t),**structure(Particles(xx,vv,m,c,g))) for t,xx,vv in zip(d['t'],d['x'],d['v'])]
        import csv
        with (out/'isolation_history.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
        changes={k:rows[-1][k]/rows[0][k]-1 for k in ['r50','z_rms','sigma_z','sigma_R']}
        result=dict(duration=duration,relative_changes=changes,verdict='Measured isolation evolution; no automatic pass',
            input_precision=precision,parent_solver_sha256=meta['parent_solver_sha256'],solver_sha256=meta['solver_sha256'])
        (out/'isolation_summary.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
    print('Completed:',out,flush=True)

if __name__=='__main__':main()
