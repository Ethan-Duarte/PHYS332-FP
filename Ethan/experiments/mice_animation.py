"""High-particle-count encounter animation. Run with: python -m experiments.mice_animation"""
import argparse
import json
from functools import partial
from pathlib import Path

from nbody.galaxies.generate import galaxy
from nbody.encounters import parabolic_pair
from nbody.integrate import evolve
from nbody.io import save_trajectory, load_trajectory
from nbody.animation import render

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description='Simulate and animate a Mice-style encounter')
    parser.add_argument('--disk', type=int, default=800, help='stellar disk particles per galaxy')
    parser.add_argument('--bulge', type=int, default=100, help='bulge particles per galaxy')
    parser.add_argument('--halo', type=int, default=350, help='live halo particles per galaxy')
    parser.add_argument('--solver', choices=['fast', 'direct', 'tree'], default='fast')
    parser.add_argument('--backend', choices=['preview', 'agama'], default='preview')
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--dt', type=float, default=.006)
    parser.add_argument('--steps', type=int, default=240)
    parser.add_argument('--save-every', type=int, default=4)
    parser.add_argument('--fps', type=int, default=15)
    parser.add_argument('--lim', type=float, default=1.5)
    parser.add_argument('--view', choices=['xy', 'xz', 'both'], default='both')
    parser.add_argument('--format', choices=['gif', 'mp4', 'both'], default='gif')
    parser.add_argument('--name', default='mice_animation', help='output filename prefix')
    parser.add_argument('--show-halo', action='store_true')
    parser.add_argument('--load', type=Path, help='render existing trajectory without rerunning physics')
    args = parser.parse_args()
    out = ROOT / 'results'
    out.mkdir(exist_ok=True)

    if args.load:
        t, x, v, m, comp, ids, metadata = load_trajectory(args.load)
        base = args.load.stem.removesuffix('_trajectory') if args.name == 'mice_animation' else args.name
    else:
        config = json.loads((ROOT / 'presets/mice.json').read_text())
        config['seed'] = args.seed
        config['particles'] = dict(disk=args.disk, bulge=args.bulge, halo=args.halo)
        g1 = galaxy(config, args.backend)
        g2 = galaxy(dict(config, seed=args.seed+1), args.backend)
        e = config['encounter']
        initial = parabolic_pair(g1, g2, e['separation'], e['pericenter'], e['orientations'])
        if args.solver == 'fast':
            from nbody.gravity.fast_direct import acceleration
        elif args.solver == 'tree':
            from nbody.gravity.barnes_hut import acceleration
        else:
            from nbody.gravity.direct import acceleration
        force = partial(acceleration, eps=.01, **({'theta': .5} if args.solver == 'tree' else {}))
        print(f'Running {len(initial.m)} live particles for {args.steps} steps using {args.solver}...')
        snapshots = list(evolve(initial, force, args.dt, args.steps, args.save_every))
        metadata = dict(config=config, backend=args.backend, solver=args.solver,
                        eps=.01, dt=args.dt, steps=args.steps, save_every=args.save_every,
                        caveat='Preview velocities are not in verified equilibrium.')
        base = args.name
        trajectory = out / f'{base}_trajectory.npz'
        save_trajectory(trajectory, snapshots, metadata)
        print(f'Saved {len(snapshots)} frames to {trajectory}')
        t, x, v, m, comp, ids, _ = load_trajectory(trajectory)

    for fmt in ('gif', 'mp4') if args.format == 'both' else (args.format,):
        path = out / f'{base}.{fmt}'
        render(t, x, comp, ids, path, fps=args.fps, lim=args.lim,
               view=args.view, show_halo=args.show_halo)
        print(f'Saved animation: {path}')


if __name__ == '__main__':
    main()
