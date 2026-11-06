"""Spiral modes and disk-heating diagnostics for an isolated, initially XY-aligned disk.

Modes measure the stellar disk (component 0). Gravity for approximate Toomre Q
includes disk, bulge, and halo; G=1 and the run's Plummer softening are used.
"""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from numba import njit


@njit(cache=True)
def circular_speed_squared(x, m, center, radii, eps, angles=8):
    """Azimuthal average of -R a_R at z=center_z, using all live particles."""
    result = np.zeros(len(radii))
    for k in range(len(radii)):
        r = radii[k]
        for a in range(angles):
            phi = 2 * np.pi * a / angles
            c, s = np.cos(phi), np.sin(phi)
            px, py, pz = center[0] + r*c, center[1] + r*s, center[2]
            radial_acc = 0.
            for j in range(len(m)):
                dx, dy, dz = x[j, 0]-px, x[j, 1]-py, x[j, 2]-pz
                inv = 1 / ((dx*dx + dy*dy + dz*dz + eps*eps)**1.5)
                radial_acc += m[j] * (dx*c + dy*s) * inv
            result[k] -= r * radial_acc / angles
    return result


def disk_modes(p, rd, eps, edges=None, min_count=20, compute_q=True):
    """Return annular Fourier modes, kinematics and an approximate axisymmetric Q."""
    if edges is None:
        edges = rd * np.array([.5, 1, 1.5, 2, 2.5, 3, 4, 5])
    disk = p.component == 0
    if not np.any(disk):
        return []
    x, v, m = p.x[disk], p.v[disk], p.m[disk]
    center = np.average(x, axis=0, weights=m)
    bulk = np.average(v, axis=0, weights=m)
    x, v = x-center, v-bulk
    R = np.hypot(x[:, 0], x[:, 1])
    phi = np.arctan2(x[:, 1], x[:, 0])
    vr = (x[:, 0]*v[:, 0] + x[:, 1]*v[:, 1]) / np.maximum(R, 1e-12)
    vphi = (x[:, 0]*v[:, 1] - x[:, 1]*v[:, 0]) / np.maximum(R, 1e-12)
    mid = .5*(edges[1:]+edges[:-1])
    vc2 = circular_speed_squared(p.x, p.m, center, mid, float(eps)) if compute_q else np.full(len(mid), np.nan)
    omega2 = vc2 / mid**2
    kappa2 = np.gradient(mid**4 * omega2, mid, edge_order=1) / mid**3
    rows=[]
    for i, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        choose = (R >= lo) & (R < hi)
        n = int(choose.sum())
        if n < min_count:
            continue
        mass = m[choose]
        total = mass.sum()
        neff = total**2 / (mass @ mass)
        density = total / (np.pi*(hi**2-lo**2))
        def dispersion(values):
            subset = values[choose]
            mu = np.average(subset, weights=mass)
            return float(np.sqrt(np.average((subset-mu)**2, weights=mass)))
        sigR = dispersion(vr)
        q = sigR*np.sqrt(kappa2[i])/(3.36*density) if kappa2[i] > 0 and density > 0 else np.nan
        modes = [(mass @ np.exp(1j*k*phi[choose])) / total for k in range(1,5)]
        rows.append(dict(bin=i, r_min=float(lo), r_max=float(hi), radius=float(mid[i]),
                         count=n, mass=float(total), neff=float(neff), noise_rms=float(1/np.sqrt(neff)),
                         A1=float(abs(modes[0])), A2=float(abs(modes[1])),
                         A3=float(abs(modes[2])), A4=float(abs(modes[3])),
                         phase2=float(np.angle(modes[1])/2),
                         sigma_R=sigR, sigma_phi=dispersion(vphi),
                         z_rms=float(np.sqrt(np.average(x[choose,2]**2, weights=mass))),
                         surface_density=float(density),
                         vcirc=float(np.sqrt(vc2[i])) if vc2[i] > 0 else np.nan,
                         kappa=float(np.sqrt(kappa2[i])) if kappa2[i] > 0 else np.nan,
                         Q_star=float(q)))
    return rows


def spiral_study(snapshots, rd, eps, folder):
    """Save radial/time-series CSV, 2 figures, and a cautious summary JSON."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    annuli, summary = [], []
    for time, p in snapshots:
        rows = disk_modes(p, rd, eps)
        for row in rows:
            row['time'] = float(time)
        annuli.extend(rows)
        if not rows:
            continue
        weights = np.array([r['mass'] for r in rows])
        summary.append(dict(time=float(time),
                            **{f'A{k}': float(np.average([r[f'A{k}'] for r in rows], weights=weights)) for k in range(1,5)},
                            noise_rms=float(np.average([r['noise_rms'] for r in rows], weights=weights)),
                            sigma_R=float(np.average([r['sigma_R'] for r in rows], weights=weights)),
                            z_rms=float(np.average([r['z_rms'] for r in rows], weights=weights))))
    if not summary:
        report = {'status':'insufficient disk resolution',
                  'caution':'No annulus contained 20 disk particles; increase stellar particle count.'}
        (folder/'spiral_report.json').write_text(json.dumps(report, indent=2))
        return report, None
    def csv_file(name, data):
        with (folder/name).open('w', newline='') as f:
            w=csv.DictWriter(f, fieldnames=list(data[0])); w.writeheader(); w.writerows(data)
    csv_file('spiral_modes.csv',summary)
    csv_file('spiral_annuli.csv',annuli)
    ts = np.array([r['time'] for r in summary])
    a2 = np.array([r['A2'] for r in summary])
    peak = int(np.argmax(a2))
    report = dict(
        model='isolated disk; components 0 only for Fourier and stellar kinematics',
        A2_initial=float(a2[0]), A2_peak=float(a2[peak]),
        time_of_A2_peak=float(ts[peak]), A2_final=float(a2[-1]),
        noise_rms_at_peak=float(summary[peak]['noise_rms']),
        A2_to_noise_at_peak=float(a2[peak]/summary[peak]['noise_rms']),
        sigma_R_change=float(summary[-1]['sigma_R']/summary[0]['sigma_R']-1) if summary[0]['sigma_R'] else None,
        thickness_change=float(summary[-1]['z_rms']/summary[0]['z_rms']-1) if summary[0]['z_rms'] else None,
        definition='A_m is mass-weighted mean of abs(annular complex Fourier amplitudes), over fixed radial bins',
        q_definition='Qstar=sigma_R*kappa/(3.36*G*Sigma); kappa from azimuthally averaged softened live-particle radial force, G=1',
        caution='Finite-N noise, self-gravity, relaxation, bars and spirals can contribute. No equilibrium/spiral-mode validation implied.')
    (folder/'spiral_report.json').write_text(json.dumps(report, indent=2))

    fig, axes = plt.subplots(2, 2, figsize=(11,7.5))
    for k in range(1,5):
        axes[0,0].plot(ts,[row[f'A{k}'] for row in summary],label=f'm={k}')
    axes[0,0].plot(ts,[row['noise_rms'] for row in summary],':',color='black',label='shot-noise RMS')
    axes[0,0].set(ylabel='Mean annular |Aₘ|',title='Stellar Fourier modes')
    axes[0,0].legend(ncol=3,fontsize=8)
    sr0 = summary[0]['sigma_R'] or np.nan
    z0 = summary[0]['z_rms'] or np.nan
    axes[0,1].plot(ts,[r['sigma_R']/sr0 for r in summary],label='σR / initial')
    axes[0,1].plot(ts,[r['z_rms']/z0 for r in summary],label='height / initial')
    axes[0,1].set(ylabel='Relative to initial',title='Disk heating / thickening')
    axes[0,1].legend(fontsize=8)
    selected_times = sorted(set([0,peak,len(summary)-1]))
    for j in selected_times:
        rows = [r for r in annuli if abs(r['time']-ts[j])<1e-10]
        axes[1,0].plot([r['radius']/rd for r in rows],[r['A2'] for r in rows], 'o-',ms=3,label=f't={ts[j]:.3f}')
        axes[1,1].plot([r['radius']/rd for r in rows],[r['Q_star'] for r in rows], 'o-',ms=3,label=f't={ts[j]:.3f}')
    axes[1,0].set(xlabel='R/Rd',ylabel='|A₂(R)|',title='Two-armed amplitude by radius')
    axes[1,1].set(xlabel='R/Rd',ylabel='Q* (approx.)',title='Local axisymmetric Q estimate')
    axes[1,0].legend(fontsize=8)
    axes[1,1].legend(fontsize=8)
    for ax in axes.flat:
        ax.grid(alpha=.2)
    for ax in axes[0]: ax.set_xlabel('Time (code units)')
    fig.tight_layout()
    fig.savefig(folder/'spiral_diagnostics.png',dpi=145)
    plt.close(fig)

    fig, (ax1,ax2) = plt.subplots(1,2,figsize=(10,3.5))
    for j in selected_times:
        rows = [r for r in annuli if abs(r['time']-ts[j])<1e-10]
        if not rows: continue
        rr=np.array([r['radius']/rd for r in rows]);phase=np.array([r['phase2'] for r in rows])
        ax1.plot(rr,[r['A2'] for r in rows],'o-',label=f't={ts[j]:.3f}')
        ax2.plot(rr,np.unwrap(2*phase)/2,'o-',label=f't={ts[j]:.3f}')
    ax1.set(xlabel='R/Rd',ylabel='|A₂|',title='Mode strength')
    ax2.set(xlabel='R/Rd',ylabel='Phase φ₂ (rad; unwrapped)',title='Mode phase across radius')
    for ax in (ax1,ax2): ax.legend(fontsize=8);ax.grid(alpha=.2)
    fig.tight_layout();fig.savefig(folder/'spiral_phase.png',dpi=145);plt.close(fig)
    return report, folder/'spiral_diagnostics.png'
