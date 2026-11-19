"""Run the moment-adjusted galaxy pair with installed Zeno (Windows/WSL or Linux).
Exploratory collisionless geometry trial; not an exact Barnes/Privon reproduction.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import struct
import subprocess
import sys
import time

import numpy as np
from nbody.io import load, save
from nbody.particles import Particles
from nbody.encounters import parabolic_pair
from nbody.galaxies.provenance import check_evolution_softening


def linux_path(path):
    path = Path(path).resolve()
    if os.name != 'nt':
        return str(path)
    return subprocess.check_output(['wsl', '-d', 'Ubuntu', '--', 'wslpath', '-a', str(path)], text=True).strip()


def linux_run(arguments, progress=None):
    command = (['wsl', '-d', 'Ubuntu', '--'] if os.name == 'nt' else []) + arguments
    if progress is None:
        subprocess.run(command, check=True)
        return
    process = subprocess.Popen(command)
    last = ''
    while process.poll() is None:
        if Path(progress).exists():
            lines = Path(progress).read_text(errors='replace').splitlines()
            # Zeno's diagnostic rows start with time and contain eight numbers.
            rows = [line.strip() for line in lines if len(line.split()) == 8
                    and all(_number(word) for word in line.split())]
            if rows and rows[-1] != last:
                last = rows[-1]
                print('Zeno: time=' + last.split()[0] + '; CPU minutes=' + last.split()[-1], flush=True)
        time.sleep(5)
    if process.returncode:
        raise subprocess.CalledProcessError(process.returncode, command)


def _number(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


def convert_frames(binary, initial, metadata, destination):
    times, positions, velocities = [], [], []
    with Path(binary).open('rb') as stream:
        while header := stream.read(12):
            if len(header) != 12:
                raise ValueError('Incomplete frame header')
            n, time_value = struct.unpack('<id', header)
            if n != len(initial.m):
                raise ValueError('Particle count changed unexpectedly')
            rows = np.fromfile(stream, dtype='<f8', count=n*7)
            if rows.size != n*7:
                raise ValueError('Incomplete output frame')
            rows = rows.reshape(n, 7)
            if not np.allclose(rows[:, 0], initial.m, rtol=1e-12, atol=0):
                raise ValueError('Output masses/order differ from input')
            if not np.isfinite(rows).all():
                raise ValueError('Nonfinite output')
            times.append(time_value)
            positions.append(rows[:, 1:4].astype(np.float32))
            velocities.append(rows[:, 4:7].astype(np.float32))
    if not times or times[0] != 0 or np.any(np.diff(times) <= 0):
        raise ValueError('Missing initial frame or invalid output times')
    if not np.allclose(positions[0], initial.x, rtol=2e-6, atol=2e-6):
        raise ValueError('Initial snapshot round trip changed positions')
    if not np.allclose(velocities[0], initial.v, rtol=2e-6, atol=2e-6):
        raise ValueError('Initial snapshot round trip changed velocities')
    # float32 storage is for replay; solver itself runs in double precision.
    np.savez_compressed(destination, t=np.asarray(times), x=np.stack(positions),
                        v=np.stack(velocities), m=initial.m,
                        component=initial.component, galaxy=initial.galaxy,
                        metadata=json.dumps(metadata))
    return times


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('candidate', type=Path)
    ap.add_argument('--out', required=True, type=Path)
    ap.add_argument('--zeno-root', default='/home/whisky/nbody-tools/zeno-reference')
    ap.add_argument('--duration', type=float, default=1.5)
    ap.add_argument('--dt', type=float, default=1/512)
    ap.add_argument('--save-every', type=int, default=8)
    ap.add_argument('--theta', type=float, default=.7)
    ap.add_argument('--separation', type=float, default=1.)
    ap.add_argument('--pericenter', type=float, default=.25)
    ap.add_argument('--prepare-only', action='store_true')
    args = ap.parse_args()
    if min(args.duration, args.dt, args.theta) <= 0 or args.save_every < 1:
        ap.error('Duration, dt, theta and save-every must be positive')
    if not 0 < args.pericenter <= args.separation:
        ap.error('Require 0 < pericenter <= separation')
    if args.duration < args.dt*args.save_every:
        ap.error('Duration must allow at least one saved evolved frame; reduce save-every')
    if not np.isclose(args.duration/args.dt, round(args.duration/args.dt), rtol=0, atol=1e-8):
        ap.error('Duration must be an integer multiple of dt')
    p, source = load(args.candidate)
    if source.get('model_id') not in ('zeno_reference_bdh_v1', 'zeno_reference_bdh_size_screen_v1'):
        ap.error('Use jeans_candidate.npz from the corrected reference control')
    eps = float(source['epsilon'])
    check_evolution_softening(source, eps)
    if not all(np.isfinite(a).all() for a in (p.x,p.v,p.m)) or np.any(p.m <= 0):
        ap.error('Invalid particles')
    disk = p.component == 0
    if not disk.any() or not np.any(p.component == 1):
        ap.error('Missing disk or bulge')
    center = np.average(p.x[disk], axis=0, weights=p.m[disk])
    velocity = np.average(p.v[disk], axis=0, weights=p.m[disk])
    galaxy = Particles(p.x-center, p.v-velocity, p.m, p.component, p.galaxy)
    phi = np.degrees(np.arctan2(2*np.sqrt(args.pericenter*(args.separation-args.pericenter))/args.separation,
                               2*args.pericenter/args.separation-1))
    orientations = [[25., -30.+phi], [40., 60.+phi]]
    initial = parabolic_pair(galaxy, galaxy, args.separation, args.pericenter, orientations)
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    metadata = dict(setup='Mice geometry trial', backend='Zeno treecode_dp',
                    source=source, source_sha256=hashlib.sha256(args.candidate.read_bytes()).hexdigest(),
                    sampling='Same particle realization rotated into both galaxies; not independent samples',
                    eps=eps, dt=args.dt, theta=args.theta, usequad=True, duration=args.duration,
                    save_every=args.save_every, nominal_separation=args.separation,
                    nominal_pericenter=args.pericenter, orientations=orientations,
                    reference_pericenter_angle_deg=float(phi),
                    removed_disk_position=center.tolist(), removed_disk_velocity=velocity.tolist(),
                    caveat='Moment-adjusted Identikit-style progenitors with Barnes-inspired orbit. Not an exact published fit or validated live equilibrium.',
                    references=['https://arxiv.org/abs/astro-ph/0402248',
                                'https://arxiv.org/abs/0811.3039',
                                'https://github.com/joshuabarnes/zeno'])
    save(out/'initial_particles.npz', initial, metadata)
    (out/'config.json').write_text(json.dumps(metadata,indent=2),encoding='utf-8')
    with (out/'initial.txt').open('w') as stream:
        stream.write(f'{len(initial.m)}\n')
        np.savetxt(stream, np.column_stack([initial.m,initial.x,initial.v]), fmt='%.17g')
    # Pass individual arguments, preserving Windows paths containing spaces.
    work = linux_path(out)
    bridge_source = linux_path(Path(__file__).with_name('zeno_snapshot_bridge.c'))
    bridge = work+'/snapshot_bridge'
    root = args.zeno_root.rstrip('/')
    linux_run(['gcc','-std=gnu99','-DLINUX','-DDOUBLEPREC','-include','unistd.h',
               '-Wno-error=implicit-function-declaration','-I'+root+'/inc',
               '-o',bridge,bridge_source,'-L'+root+'/lib','-lNBody','-lClib','-lgsl','-lgslcblas','-lm'])
    linux_run([bridge,'mode=import','in='+work+'/initial.txt','out='+work+'/initial.snap'])
    print(f'Prepared {len(initial.m):,} particles. Zeno double-precision tree; epsilon={eps}; theta={args.theta}.',flush=True)
    if args.prepare_only:
        linux_run([bridge,'mode=export','in='+work+'/initial.snap','out='+work+'/frames.bin'])
        convert_frames(out/'frames.bin',initial,metadata,out/'prepared_trajectory.npz')
        print('Snapshot round trip verified; no evolution launched.',flush=True)
        return
    print('Running collisionless encounter. Progress goes to tree.log; saved snapshots survive interruption.',flush=True)
    started=time.perf_counter()
    tree_command = [root+'/bin/treecode_dp','in='+work+'/initial.snap','out='+work+'/encounter.snap',
               'save='+work+'/state_%d','log='+work+'/tree.log',
               'eps='+str(eps),'theta='+str(args.theta),'usequad=true',
               'dtime='+str(args.dt),'tstop='+str(args.duration),
               'dtout='+str(args.dt*args.save_every),'outputs=Mass,Position,Velocity']
    metadata['tree_command'] = tree_command
    sha_command = (['wsl','-d','Ubuntu','--'] if os.name == 'nt' else []) + ['sha256sum', root+'/bin/treecode_dp']
    metadata['solver_sha256'] = subprocess.check_output(sha_command,text=True).split()[0]
    (out/'config.json').write_text(json.dumps(metadata,indent=2),encoding='utf-8')
    linux_run(tree_command, progress=out/'tree.log')
    linux_run([bridge,'mode=export','in='+work+'/encounter.snap','out='+work+'/frames.bin'])
    times=convert_frames(out/'frames.bin',initial,metadata,out/'run_trajectory.npz')
    if not np.isclose(times[-1],args.duration,rtol=0,atol=args.dt/8):
        raise RuntimeError('Run did not reach requested end time')
    print(f'Finished: {len(times)} frames in {(time.perf_counter()-started)/60:.1f} minutes.\n{out / "run_trajectory.npz"}',flush=True)

if __name__ == '__main__':
    main()
