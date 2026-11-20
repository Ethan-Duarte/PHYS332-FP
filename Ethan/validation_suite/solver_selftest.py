"""Test the project's actual force and leapfrog implementations, not a duplicate."""
import hashlib,sys
from pathlib import Path
import numpy as np
from core import reference_force

def run(a,finish,output_dir):
    out=output_dir(a.out);root=a.repo_root.resolve();sys.path.insert(0,str(root))
    from nbody.gravity.parallel_direct import acceleration
    from nbody.integrate import evolve
    from nbody.particles import Particles
    from scipy.integrate import solve_ivp
    import nbody.gravity.parallel_direct as force_module
    import nbody.integrate as integration_module
    rng=np.random.default_rng(7231);x=rng.normal(size=(32,3));m=rng.uniform(.1,2,32);eps=.05;G=1.
    expected=reference_force(x,m,np.arange(len(m)),eps,G);actual=acceleration(x,m,eps=eps,G=G)
    err=float(np.max(np.linalg.norm(actual-expected,axis=1)/np.maximum(np.linalg.norm(expected,axis=1),1e-30)))
    # Independent scipy DOP853 integration of the softened two-body ODE.
    masses=np.array([.7,1.3]);M=masses.sum();separation=1.;omega=np.sqrt(G*M/(separation**2+eps**2)**1.5);period=2*np.pi/omega
    xx=np.array([[-masses[1]/M,0,0],[masses[0]/M,0,0]])
    vv=np.array([[0,-omega*masses[1]/M,0],[0,omega*masses[0]/M,0]])
    def rhs(t,y):
        pos=y[:6].reshape(2,3);vel=y[6:].reshape(2,3)
        return np.r_[vel.ravel(),reference_force(pos,masses,[0,1],eps,G).ravel()]
    final=2*period;y0=np.r_[xx.ravel(),vv.ravel()]
    ref=solve_ivp(rhs,[0,final],y0,method='DOP853',rtol=1e-12,atol=1e-13)
    if not ref.success:raise RuntimeError('Independent orbit integration failed')
    errors=[]
    for steps in [256,512]:
        particle=Particles(xx.copy(),vv.copy(),masses,np.zeros(2,dtype=int),np.zeros(2,dtype=int))
        snapshots=list(evolve(particle,lambda x,m:acceleration(x,m,eps=eps,G=G),final/steps,steps,save_every=steps))
        state=snapshots[-1][1];errors.append(float(np.sqrt(np.mean((state.x-ref.y[:6,-1].reshape(2,3))**2))/separation))
    ratio=errors[0]/errors[1]
    files=[Path(force_module.__file__),Path(integration_module.__file__)]
    checks=[dict(name='Direct force vs independent vectorized sum',status='PASS' if err<1e-10 else 'FAIL',value=err,limit=1e-10),
        dict(name='Softened circular orbit: endpoint RMS / separation',status='PASS' if errors[1]<.01 else 'FAIL',value=errors[1],limit=.01),
        dict(name='Leapfrog second-order convergence on fixture',status='PASS' if 3<ratio<5 else 'FAIL',value=ratio,expected='3 < coarse/fine error < 5')]
    return finish(out,'selftest',files,checks,dict(force_max_relative_error=err,orbit_coarse_error=errors[0],orbit_fine_error=errors[1],orbit_error_ratio=ratio,period=period),
        ['Acceptance bounds are declared for this numerical fixture, not galaxy or observational tolerances.',
         'This tests the native project direct solver and leapfrog. It does not test Zeno or GADGET binaries.',
         'Two orbital periods; independent DOP853 rtol=1e-12. Passing this fixture is not full galaxy validation.'])
