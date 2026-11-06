"""Measurements for isolation screens and solver convergence. G=1 code units."""
import csv
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numba import njit


@njit(cache=True)
def _potential(x, m, eps):
    result = 0.
    for i in range(len(m)):
        for j in range(i):
            dx = x[i, 0]-x[j, 0]
            dy = x[i, 1]-x[j, 1]
            dz = x[i, 2]-x[j, 2]
            result -= m[i]*m[j]/np.sqrt(dx*dx+dy*dy+dz*dz+eps*eps)
    return result


def total_energy(p, eps):
    kinetic = 0.5*np.sum(p.m*np.sum(p.v*p.v, axis=1))
    return float(kinetic+_potential(p.x, p.m, float(eps)))


def weighted_half_radius(radius, mass):
    order = np.argsort(radius)
    i = np.searchsorted(np.cumsum(mass[order]), 0.5*np.sum(mass))
    return float(radius[order[min(i,len(order)-1)]])


def structure(p):
    """Internal disk structural statistics, or spherical half-mass radius for Plummer."""
    disk = p.component == 0
    if np.any(disk):
        x, v, m = p.x[disk], p.v[disk], p.m[disk]
        # Isolated disk aligned with xy plane; recenter on disk centroid.
        cm = np.average(x, axis=0, weights=m)
        vm = np.average(v, axis=0, weights=m)
        x, v = x-cm, v-vm
        R = np.linalg.norm(x[:, :2], axis=1)
        vr = (x[:, 0]*v[:, 0]+x[:, 1]*v[:, 1])/np.maximum(R,1e-12)
        return dict(r50=weighted_half_radius(R,m),
                    z_rms=float(np.sqrt(np.average(x[:, 2]**2, weights=m))),
                    sigma_R=float(np.sqrt(np.average((vr-np.average(vr,weights=m))**2, weights=m))))
    cm = np.average(p.x, axis=0, weights=p.m)
    r = np.linalg.norm(p.x-cm, axis=1)
    return dict(r50=weighted_half_radius(r,p.m))


def inspect_isolation(snapshots, eps, folder, backend='preview'):
    """Creates CSV, PNG, JSON report. Thresholds are screening heuristics, not validation."""
    folder.mkdir(parents=True, exist_ok=True)
    samples = []
    stride=max(1,(len(snapshots)-1)//10)
    for i,(t,p) in enumerate(snapshots):
        row = dict(time=float(t), **structure(p))
        if (i % stride == 0 or i==len(snapshots)-1) and len(p.m)<=6000:
            row['energy'] = total_energy(p,eps)
        samples.append(row)
    start,end = samples[0],samples[-1]
    changes = {key: float((end[key]-start[key])/start[key]) for key in start if key in ('r50','z_rms','sigma_R') and start[key]}
    energy_change = ((end['energy']-start['energy'])/abs(start['energy'])) if 'energy' in end and start['energy'] else None
    max_energy_drift=(max(abs((row['energy']-start['energy'])/start['energy']) for row in samples if 'energy' in row)
                      if 'energy' in start and start['energy'] else None)
    screen = abs(changes['r50'])<=.05 and abs(changes.get('z_rms',0))<=.10 and (max_energy_drift is None or max_energy_drift<=.01)
    report = dict(backend=backend, verdict='screening pass' if screen else 'screening fail',
                  validation='NOT validated: test multiple resolutions, seeds and several orbital periods',
                  relative_changes=changes, relative_energy_change=energy_change,
                  max_abs_energy_drift=max_energy_drift,
                  thresholds=dict(r50=.05, z_rms=.10, energy=.01),
                  n_particles=len(snapshots[0][1].m), duration=float(end['time']))
    keys=['time','r50','z_rms','sigma_R','energy']
    with (folder/'isolation_metrics.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys,extrasaction='ignore');w.writeheader();w.writerows(samples)
    (folder/'isolation_report.json').write_text(json.dumps(report,indent=2))
    fig,axs=plt.subplots(1,2,figsize=(10,3.3))
    times=[s['time'] for s in samples]
    axs[0].plot(times,[s['r50']/start['r50'] for s in samples],label='r50 / initial')
    if 'z_rms' in start:
        axs[0].plot(times,[s['z_rms']/start['z_rms'] for s in samples],label='thickness / initial')
    axs[0].axhline(1,color='gray',ls=':',lw=1)
    axs[0].set(xlabel='Time (code units)',ylabel='Relative structural metric',title='Isolated structure')
    axs[0].legend()
    e = [(s['energy']-start['energy'])/abs(start['energy']) for s in samples if 'energy' in s] if 'energy' in start else []
    if e:
        axs[1].plot([s['time'] for s in samples if 'energy' in s],e,'o-',ms=3)
    axs[1].set(xlabel='Time (code units)',ylabel='(E - E0) / |E0|',title='Sampled energy drift')
    for a in axs:a.grid(alpha=.2)
    fig.tight_layout();fig.savefig(folder/'isolation_diagnostics.png',dpi=140);plt.close(fig)
    return report, fig if False else folder/'isolation_diagnostics.png'


def encounter_report(snapshots, eps, folder):
    """Energy and sampled bulge centers. One interval is not a convergence proof."""
    folder.mkdir(parents=True,exist_ok=True)
    from .analysis import core_separation
    times=np.array([t for t,_ in snapshots])
    sep=np.array([core_separation(p) for _,p in snapshots])
    e0=total_energy(snapshots[0][1],eps) if len(snapshots[0][1].m)<=10000 else None
    e1=total_energy(snapshots[-1][1],eps) if e0 is not None else None
    drift=(e1-e0)/abs(e0) if e0 else None
    report=dict(min_bulge_separation=float(np.min(sep)),time_min_separation=float(times[np.argmin(sep)]),
                energy_relative_change=drift, n_particles=len(snapshots[0][1].m),
                note='Bulge-centroid minima and endpoint energy drift are diagnostics, not Mice validation.')
    (folder/'encounter_report.json').write_text(json.dumps(report,indent=2))
    with (folder/'encounter_metrics.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['time','bulge_separation']);w.writerows(zip(times,sep))
    return report


def profile_checks(p, config, folder, backend='preview'):
    """CDF sampling checks; discrepancies are diagnostics, not confidence statements."""
    from .galaxies.profiles import dehnen_mass, disk_mass
    folder.mkdir(parents=True,exist_ok=True)
    fig, ax = plt.subplots(1,3,figsize=(12,3.1))
    results={}
    for k,(index,label) in enumerate([(0,'Disk R'),(1,'Bulge r'),(2,'Halo r')]):
        mask=p.component==index
        if not mask.any():continue
        part=p.x[mask]
        cm=np.average(part,axis=0,weights=p.m[mask]); part=part-cm
        r=np.linalg.norm(part[:,:2],axis=1) if index==0 else np.linalg.norm(part,axis=1)
        grid=np.linspace(0,np.percentile(r,99.5)*1.15,120)
        observed=np.searchsorted(np.sort(r),grid,side='right')/len(r)
        spec=config[['disk','bulge','halo'][index]]
        if index==0:
            reference=disk_mass(grid,spec['mass'],spec['scale'])/spec['mass']
        else:
            reference=dehnen_mass(grid,spec['mass'],spec['scale'],1 if index==1 else 0)/spec['mass']
        if backend=='preview' and spec.get('cutoff'):
            rmax=spec['cutoff']
            total=(disk_mass(rmax,spec['mass'],spec['scale']) if index==0 else
                   dehnen_mass(rmax,spec['mass'],spec['scale'],1 if index==1 else 0))
            reference=np.minimum(reference*spec['mass']/total,1)
        ks=float(np.max(np.abs(observed-reference)))
        results[label]=dict(max_cdf_discrepancy=ks,n=int(len(r)))
        ax[k].plot(grid,observed,label='Realization')
        ax[k].plot(grid,reference,'--',label='Target analytic')
        ax[k].set(xlabel='R (code units)' if index==0 else 'r (code units)',
                  ylabel='Enclosed mass fraction',title=f'{label} · max gap {ks:.3f}')
        ax[k].grid(alpha=.15)
    ax[0].legend(fontsize=8);fig.tight_layout()
    fig.savefig(folder/'initial_profiles.png',dpi=140);plt.close(fig)
    (folder/'initial_profiles.json').write_text(json.dumps(results,indent=2))
    return results
