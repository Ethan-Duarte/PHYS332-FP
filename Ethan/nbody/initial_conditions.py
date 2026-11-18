"""Unconstrained relaxation and restart: compare, don't assume, equilibrium."""
import json
from pathlib import Path
from functools import partial

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from .integrate import evolve
from .io import save, save_trajectory
from .diagnostics import structure, inspect_isolation
from .virial import virial_state, rescale_velocities


def force_for(solver, eps, theta):
    if solver == 'Barnes–Hut':
        from .gravity.fast_barnes_hut import acceleration
        return partial(acceleration, eps=eps, theta=theta)
    from .gravity.fast_direct import acceleration
    return partial(acceleration, eps=eps)


def run_relax_restart(initial, folder, *, eps=.0125, dt=.003,
                      relax_steps=100, restart_steps=100, save_every=10,
                      solver='Fast direct', theta=.5, target_virial=None,
                      cancelled=None, source='saved initial state'):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    force = force_for(solver, eps, theta)
    p = initial
    preparation = None
    if target_virial is not None:
        p, preparation = rescale_velocities(p, force, target_virial)
    save(folder/'prepared_initial.npz', p,
         dict(method='virial+relaxation' if preparation else 'relaxation', source=source,
              eps=eps, dt=dt, solver=solver, preparation=preparation))

    relax = list(evolve(p, force, dt, int(relax_steps), int(save_every), cancelled=cancelled))
    meta = dict(setup='Isolated galaxy', method='unconstrained relaxation',
                source=source, eps=eps, dt=dt, solver=solver, theta=theta)
    save_trajectory(folder/'relaxation_trajectory.npz', relax, meta)
    if relax[-1][0] < relax_steps*dt - 1e-10:
        report = dict(status='cancelled', elapsed=relax[-1][0],
                      caution='Incomplete relaxation; no final galaxy candidate saved.')
        (folder/'initial_condition_report.json').write_text(json.dumps(report, indent=2))
        return report

    relaxed = relax[-1][1]
    save(folder/'relaxed_candidate.npz', relaxed,
         dict(method='unconstrained relaxation', source=source, eps=eps, dt=dt,
              relax_steps=relax_steps, solver=solver, target_virial=target_virial,
              validation='candidate; restart and convergence screening required'))
    restart = list(evolve(relaxed, force, dt, int(restart_steps), int(save_every), cancelled=cancelled))
    save_trajectory(folder/'restart_trajectory.npz', restart,
                    {**meta, 'method':'restart from relaxed candidate'})

    pre = structure(p)
    post = structure(relaxed)
    final = structure(restart[-1][1])
    virials = [virial_state(state, force)['ratio'] for state in
               [p, relaxed, restart[-1][1]]]
    reports = {}
    for name, frames in [('relaxation', relax), ('restart', restart)]:
        if len(frames) > 1:
            stage_dir = folder / name
            _, _ = inspect_isolation(frames, eps, stage_dir, source)
            reports[name] = json.loads((stage_dir/'isolation_report.json').read_text())

    fig, axs = plt.subplots(1, 3, figsize=(13, 3.5))
    for name, frames, shift in [('Relax', relax, 0.), ('Restart', restart, relax[-1][0])]:
        ts = [t+shift for t, _ in frames]
        r = [structure(s)['r50'] for _, s in frames]
        axs[0].plot(ts, r, marker='.', label=name)
        if 'z_rms' in pre:
            z = [structure(s)['z_rms'] for _, s in frames]
            sig = [structure(s)['sigma_R'] for _, s in frames]
            axs[1].plot(ts, z, marker='.', label=name)
            axs[2].plot(ts, sig, marker='.', label=name)
    axs[0].set(ylabel='Half-mass radius', xlabel='Time (code units)')
    axs[1].set(ylabel='Disk RMS thickness' if 'z_rms' in pre else 'Not a disk',
               xlabel='Time (code units)')
    axs[2].set(ylabel='Disk radial velocity dispersion' if 'sigma_R' in pre else 'Not a disk',
               xlabel='Time (code units)')
    for ax in axs:
        ax.axvline(relax[-1][0], color='gray', ls='--', lw=1)
        if ax.get_legend_handles_labels()[0]: ax.legend()
        ax.grid(alpha=.2)
    fig.tight_layout()
    fig.savefig(folder/'relax_restart.png', dpi=145)
    plt.close(fig)

    report = dict(status='complete' if restart[-1][0] >= dt*restart_steps-1e-10 else 'restart cancelled',
                  n_particles=len(p.m), source=source, solver=solver, eps=eps, theta=theta,
                  dt=dt, relax_steps=relax_steps, restart_steps=restart_steps,
                  preparation=preparation,
                  virial_ratio=dict(prepared=virials[0], after_relaxation=virials[1],
                                    after_restart=virials[2]),
                  structure=dict(prepared=pre, after_relaxation=post, after_restart=final),
                  isolation_reports=reports,
                  caution='Unconstrained relaxation changes the physical density/velocity distributions. Virial ratio near 1 and a short restart are not proof of equilibrium.')
    (folder/'initial_condition_report.json').write_text(json.dumps(report, indent=2))
    return report
