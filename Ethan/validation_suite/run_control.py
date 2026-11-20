"""Native direct/tree controls using the project's own force and integrator code."""
import argparse,hashlib,json,sys,time
from pathlib import Path
import numpy as np
from core import load,law,digest,write_json

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('input',type=Path);p.add_argument('--repo-root',type=Path,default=Path('.'))
    p.add_argument('--mode',choices=['isolate','replay-initial'],required=True)
    p.add_argument('--solver',choices=['direct','tree'],default='direct');p.add_argument('--galaxy',type=int,default=0)
    p.add_argument('--duration',type=float);p.add_argument('--dt',type=float);p.add_argument('--theta',type=float,default=.5)
    p.add_argument('--save-interval',type=float,default=1/64);p.add_argument('--threads',type=int,default=4)
    p.add_argument('--epsilon',type=float);p.add_argument('--G',type=float);p.add_argument('--kernel',choices=['plummer'])
    p.add_argument('--prepare-only',action='store_true');p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();r=load(a.input,snapshot=True);eps,G,kernel=law(r,a.epsilon,a.G,a.kernel)
    dt=float(a.dt if a.dt is not None else r.meta.get('dt',1/512));duration=float(a.duration if a.duration is not None else r.meta.get('duration',1.25))
    if dt<=0 or duration<=0 or a.save_interval<=0 or a.threads<1 or a.theta<=0:p.error('Positive settings required')
    steps=round(duration/dt);stride=round(a.save_interval/dt)
    if steps<1 or stride<1 or not np.isclose(steps*dt,duration,rtol=0,atol=1e-10) or not np.isclose(stride*dt,a.save_interval,rtol=0,atol=1e-10):p.error('Duration and save interval must be integer multiples of dt')
    x=r.x[0].copy();v=r.v[0].copy();m=r.m.copy();c=r.c.copy();g=r.g.copy();precision='replay initial frame promoted to float64'
    initial_path=a.input.with_name('initial_particles.npz')
    if initial_path.exists() and initial_path.resolve()!=a.input.resolve():
        initial=load(initial_path,snapshot=True)
        if not np.array_equal(initial.m,m) or not np.array_equal(initial.c,c) or not np.array_equal(initial.g,g) or not np.allclose(initial.x[0],x,rtol=2e-6,atol=2e-6) or not np.allclose(initial.v[0],v,rtol=2e-6,atol=2e-6):p.error('Adjacent initial snapshot does not match trajectory')
        x=initial.x[0].copy();v=initial.v[0].copy();precision='adjacent original initial_particles.npz, correspondence checked'
    if a.mode=='isolate':
        if not r.labels:p.error('Isolation extraction requires galaxy labels')
        keep=g==a.galaxy
        if not keep.any():p.error('Galaxy label not found')
        x=x[keep];v=v[keep];m=m[keep];c=c[keep];g=g[keep]
        x-=np.average(x,axis=0,weights=m);v-=np.average(v,axis=0,weights=m)
    a.out.mkdir(parents=True,exist_ok=False)
    sys.path.insert(0,str(a.repo_root.resolve()))
    from nbody.particles import Particles
    from nbody.integrate import evolve
    import nbody.integrate as integration_module
    from numba import set_num_threads
    set_num_threads(a.threads)
    if a.solver=='direct':
        import nbody.gravity.parallel_direct as force_module
        force=lambda xx,mm:force_module.acceleration(xx,mm,eps=eps,G=G)
    else:
        import nbody.gravity.fast_barnes_hut as force_module
        force=lambda xx,mm:force_module.acceleration(xx,mm,eps=eps,G=G,theta=a.theta)
    solver_files=[Path(force_module.__file__),Path(integration_module.__file__)]
    if a.solver=='tree':
        import nbody.gravity.barnes_hut as tree_builder
        solver_files.append(Path(tree_builder.__file__))
    combined=hashlib.sha256(b''.join(p.read_bytes() for p in solver_files)).hexdigest()
    source=r.meta.get('source',{});parent=source.get('size_screen',{}).get('source_sha256',r.meta.get('validation_parent_model_id'))
    meta=dict(r.meta,backend='native '+a.solver,dynamics_model='self_gravitating_collisionless',eps=eps,G=G,evolution_kernel=kernel,dt=dt,theta=None if a.solver=='direct' else a.theta,
        duration=duration,validation_mode=a.mode,validation_parent_model_id=parent,input_precision=precision,
        parent_trajectory_sha256=digest(a.input),parent_solver_sha256=r.meta.get('solver_sha256'),solver_sha256=combined,
        solver_identity_kind='SHA256(force + integrator + tree builder when applicable); not compiled binary',solver_source_files={str(p):digest(p) for p in solver_files},force_source_path=str(force_module.__file__),integration_source_path=str(integration_module.__file__))
    np.savez_compressed(a.out/'initial_particles.npz',x=x,v=v,m=m,component=c,galaxy=g,metadata=json.dumps(meta))
    write_json(a.out/'config.json',meta)
    if a.prepare_only:print('Prepared initial state only:',a.out);return
    # A real full-size evaluation includes any current tree-building overhead. Warm JIT on tiny input first.
    force(x[:min(8,len(x))],m[:min(8,len(x))]);start=time.perf_counter();force(x,m);seconds=time.perf_counter()-start
    print(f'N={len(m)}; force={seconds:.3f}s; forecast {(steps+1)*seconds/60:.2f} minutes plus diagnostics/I/O',flush=True)
    snapshots=[];started=time.perf_counter()
    for t,state in evolve(Particles(x,v,m,c,g),force,dt,steps,save_every=stride):
        snapshots.append((t,state))
        np.savez_compressed(a.out/'checkpoint.npz',time=t,x=state.x,v=state.v,m=m,component=c,galaxy=g,metadata=json.dumps(meta))
        print(f't={t:.6f}/{duration:g}; elapsed={(time.perf_counter()-started)/60:.2f} min',flush=True)
    np.savez_compressed(a.out/'run_trajectory.npz',t=np.array([z[0] for z in snapshots]),x=np.stack([z[1].x for z in snapshots]),v=np.stack([z[1].v for z in snapshots]),m=m,component=c,galaxy=g,metadata=json.dumps(meta))
    print('Saved:',a.out/'run_trajectory.npz')

if __name__=='__main__':main()
