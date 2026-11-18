"""Repeatable isolated stability + timestep-refinement *screen* for one fixed galaxy.

This screens one finite-N realization. Multiple particle resolutions and random
seeds are still required before describing a disk as scientifically validated.
"""
from functools import partial
import json
from pathlib import Path
import numpy as np
from .integrate import evolve
from .diagnostics import inspect_isolation, structure
from .io import save_trajectory


def validate_initial_state(initial, folder, dt=0.003, steps=80, eps=0.0125,
                           save_every=10, backend='saved snapshot', solver='Fast direct'):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    if solver == 'Barnes–Hut':
        from .gravity.fast_barnes_hut import acceleration
        force = partial(acceleration, eps=eps, theta=0.5)
    else:
        from .gravity.fast_direct import acceleration
        force = partial(acceleration, eps=eps)
    trials = {}
    for name, factor in [('base', 1), ('half_dt', 2)]:
        dt_trial = dt / factor
        frames = list(evolve(initial, force, dt_trial, steps*factor, save_every*factor))
        run = folder / name
        run.mkdir(exist_ok=True)
        report, _ = inspect_isolation(frames, eps, run, backend)
        save_trajectory(run / 'trajectory.npz', frames,
                        dict(setup='Isolated galaxy', solver=solver, backend=backend,
                             eps=eps, dt=dt_trial, steps=steps*factor, save_every=save_every*factor))
        final = frames[-1][1]
        trials[name] = dict(report=report, final=structure(final), x=final.x)
    a, b = trials['base'], trials['half_dt']
    refinements = {
        'r50_relative': abs(a['final']['r50']-b['final']['r50'])/b['final']['r50'],
        'position_rms': float(np.sqrt(np.mean(np.sum((a['x']-b['x'])**2, axis=1))))
    }
    if 'z_rms' in a['final']:
        refinements['z_rms_relative'] = abs(a['final']['z_rms']-b['final']['z_rms'])/b['final']['z_rms']
    results = dict(
        title='Single-realization isolation + timestep refinement screen',
        backend=backend, solver=solver, n_particles=len(initial.m), eps=float(eps),
        total_time=float(dt*steps), timestep_base=float(dt), timestep_refined=float(dt/2),
        base=a['report'], half_dt=b['report'], refined_comparison=refinements,
        verdict='screening comparison only — not a validation of the galaxy or Mice encounter',
        required_next='Repeat at multiple particle resolutions, softening lengths and independent seeds; verify target density/velocity profiles and galaxy orbital periods.')
    (folder/'validation_comparison.json').write_text(json.dumps(results, indent=2))
    return results
