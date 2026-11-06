"""Matched-initial-conditions direct vs Barnes-Hut full-trajectory experiment."""
import json
from pathlib import Path
from time import perf_counter
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from functools import partial
from .gravity.fast_direct import acceleration as direct
from .gravity.fast_barnes_hut import acceleration as tree
from .gravity.direct import acceleration as numpy_direct
from .integrate import evolve
from .diagnostics import total_energy
from .analysis import core_separation
from .io import save_trajectory


def compare_trajectories(initial, folder, dt=.006, steps=30, save_every=5, eps=.0125, theta=.5):
    """Warm JIT, integrate identical initial particles, save curves and measured runtimes."""
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
    solvers={'Fast direct (Numba)':partial(direct,eps=eps),
             'Barnes-Hut (Numba)':partial(tree,eps=eps,theta=theta)}
    if len(initial.m)<=800:
        solvers['Direct (NumPy)']=partial(numpy_direct,eps=eps)
    # Small warmups remove one-time JIT compilation from timing.
    direct(initial.x[:5],initial.m[:5],eps=eps)
    tree(initial.x[:5],initial.m[:5],eps=eps,theta=theta)
    reference=direct(initial.x,initial.m,eps=eps)
    results={}
    for name,force in solvers.items():
        a=force(initial.x,initial.m)
        err=np.linalg.norm(a-reference,axis=1)
        rel=err/np.maximum(np.linalg.norm(reference,axis=1),1e-12)
        start=perf_counter()
        snapshots=list(evolve(initial,force,dt,steps,save_every))
        elapsed=perf_counter()-start
        labels=name.lower().replace(' ','_').replace('(','').replace(')','')
        save_trajectory(folder/f'{labels}.npz',snapshots,dict(solver=name,eps=eps,theta=theta,dt=dt,steps=steps))
        times=np.array([t for t,_ in snapshots])
        sep=np.array([core_separation(s) for _,s in snapshots])
        E0=total_energy(snapshots[0][1],eps)
        E1=total_energy(snapshots[-1][1],eps)
        results[name]=dict(seconds=float(elapsed),force_median_error=float(np.median(rel)),
                           force_p95_error=float(np.percentile(rel,95)),
                           energy_relative_change=float((E1-E0)/abs(E0)),
                           times=times.tolist(),separation=sep.tolist(),
                           states=snapshots)
    ref_final=results['Fast direct (Numba)']['states'][-1][1].x.copy()
    for name,entry in results.items():
        entry['rms_position_error_final']=float(np.sqrt(np.mean(np.sum((entry['states'][-1][1].x-ref_final)**2,axis=1))))
        del entry['states']
    (folder/'solver_comparison.json').write_text(json.dumps(dict(n=len(initial.m),eps=eps,dt=dt,theta=theta,solvers=results),indent=2))
    fig,ax=plt.subplots(1,2,figsize=(10,3.5))
    for name,row in results.items(): ax[0].plot(row['times'],row['separation'],label=name)
    ax[0].set(xlabel='Time',ylabel='Bulge-centroid separation',title='Matched trajectories');ax[0].legend(fontsize=8)
    names=list(results);ax[1].barh(names,[results[k]['seconds'] for k in names]);ax[1].set(xlabel='Seconds (excluding initial JIT)',title='Evolution runtime')
    for a in ax:a.grid(alpha=.2)
    fig.tight_layout();fig.savefig(folder/'solver_comparison.png',dpi=140);plt.close(fig)
    return results, folder/'solver_comparison.png'
