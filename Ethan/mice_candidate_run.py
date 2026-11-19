"""Run cached AGAMA galaxies without importing AGAMA.

Place this file in PHYS332-FP/Ethan.
Use --a and --b to supply corrected snapshots. This is an exploratory setup.
Original galaxy files are never overwritten.
"""
import argparse
import csv
import hashlib
import json
import time as clock
from functools import partial
from pathlib import Path

import numpy as np

from nbody.encounters import parabolic_pair
from nbody.integrate import evolve
from nbody.io import load, save, save_trajectory
from nbody.particles import Particles


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["encounter"])
    parser.add_argument("--a", type=Path, required=True)
    parser.add_argument("--b", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--checkpoint-every",type=int,default=8,help="Save checkpoint every N stored frames")
    parser.add_argument("--steps", type=int, default=128)
    parser.add_argument("--dt", type=float, default=1 / 1024)
    parser.add_argument("--eps", type=float, default=None)
    parser.add_argument("--theta", type=float, default=0.5)
    parser.add_argument("--solver", choices=["direct", "tree"],
                        default="direct")
    parser.add_argument("--save-every", type=int, default=16)
    parser.add_argument("--separation", type=float, default=1.0)
    parser.add_argument("--pericenter", type=float, default=0.25)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.threads < 1 or args.checkpoint_every < 1: parser.error('threads must be positive')
    if args.a.resolve()==args.b.resolve():
        parser.error('Use independently sampled progenitors, not the same snapshot twice')

    if (args.steps < 1 or args.dt <= 0 or (args.eps is not None and args.eps <= 0)
            or args.save_every < 1):
        parser.error("steps, dt, eps and save-every must be positive")
    if not 0 < args.pericenter <= args.separation:
        parser.error("Require 0 < pericenter <= separation")

    root = Path(__file__).resolve().parent
    out = args.out if args.out.is_absolute() else root / args.out
    out.mkdir(parents=True, exist_ok=False)
    sources = []

    def read_galaxy(path):
        if not path.is_absolute(): path = root/path
        particles, metadata = load(path)

        if (
            not all(
                np.isfinite(array).all()
                for array in (
                    particles.x, particles.v, particles.m
                )
            )
            or np.any(particles.m <= 0)
        ):
            raise ValueError(f"Invalid particles in {path}")

        expected = metadata.get('epsilon')
        if expected is None:
            raise ValueError('Use a corrected candidate with explicit epsilon metadata')
        if args.eps is None: args.eps = float(expected)
        if not np.isclose(args.eps,float(expected),rtol=0,atol=1e-15):
            raise ValueError('Encounter epsilon must match BOTH candidate equilibria')
        disk = particles.component == 0
        if not disk.any():
            raise ValueError(f"No stellar disk particles in {path}")

        # Locate the visible galaxy using its stellar disk.
        # Distant halo outliers remain in the simulation.
        origin = np.average(
            particles.x[disk], axis=0, weights=particles.m[disk]
        )
        bulk_velocity = np.average(
            particles.v[disk], axis=0, weights=particles.m[disk]
        )

        sources.append({
            "file": str(path.resolve()),
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "metadata": metadata,
            "removed_disk_position": origin.tolist(),
            "removed_disk_velocity": bulk_velocity.tolist(),
        })

        return Particles(
            particles.x - origin,
            particles.v - bulk_velocity,
            particles.m,
            particles.component,
            particles.galaxy,
        )

    galaxy_a = read_galaxy(args.a)

    galaxy_b = read_galaxy(args.b)

    # The initial separation points along +x. Calculate the
    # point-mass reference pericenter direction, then express
    # the literature-inspired disk orientations in that frame.
    phi = np.arctan2(
        2 * np.sqrt(
            args.pericenter * (args.separation - args.pericenter)
        ) / args.separation,
        2 * args.pericenter / args.separation - 1,
    )
    pericenter_angle = float(np.rad2deg(phi))
    orientations = [
        [25.0, -30.0 + pericenter_angle],
        [40.0, 60.0 + pericenter_angle],
    ]

    initial = parabolic_pair(
        galaxy_a,
        galaxy_b,
        args.separation,
        args.pericenter,
        orientations,
    )

    if args.solver == "tree":
        from nbody.gravity.fast_barnes_hut import acceleration

        force = partial(
            acceleration, eps=args.eps, theta=args.theta
        )
    else:
        import numba
        numba.set_num_threads(min(args.threads,numba.config.NUMBA_NUM_THREADS))
        from nbody.gravity.parallel_direct import acceleration

        force = partial(acceleration, eps=args.eps)

    metadata = {
        "setup": "Mice encounter",
        "backend": "stored AGAMA candidate",
        "model": "Disk\u2013bulge\u2013halo",
        "config": sources[0]["metadata"].get("config"),
        "sources": sources,
        "solver": args.solver,
        "eps": args.eps,
        "theta": args.theta,
        "dt": args.dt,
        "steps": args.steps,
        "save_every": args.save_every,
        "nominal_separation": args.separation,
        "nominal_pericenter": args.pericenter,
        "reference_pericenter_angle_deg": pericenter_angle,
        "orientations": orientations,
        "centering": (
            "mass-weighted stellar disk; halo outliers retained"
        ),
        "caveat": (
            "Literature-inspired trial, not a reproduced observational "
            "fit. Pericenter is a point-mass reference, not the "
            "measured closest approach of the live galaxies."
        ),
    }

    save(out / "initial_particles.npz", initial, metadata)
    (out / "config.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )

    force(initial.x[:32],initial.m[:32])
    started=clock.perf_counter()
    force(initial.x,initial.m)
    seconds=clock.perf_counter()-started
    print(f'Initial full force evaluation: {seconds:.3f} s; estimated force-only run: {seconds*(args.steps+1)/60:.1f} min',flush=True)
    if args.prepare_only:
        print('Prepared and benchmarked; no evolution requested.',flush=True)
        return
    frames = []
    try:
        for time, particles in evolve(initial,force,args.dt,args.steps,args.save_every):
            frames.append((time,particles))
            if len(frames)==1 or len(frames)%args.checkpoint_every==0 or time >= args.steps*args.dt-1e-12:
                temporary=out/'checkpoint.tmp.npz'
                save_trajectory(temporary,frames,metadata)
                temporary.replace(out/'run_trajectory.npz')
            print(f't={time:.4f}/{args.dt*args.steps:.4f}; {len(frames)} frames',flush=True)
    except KeyboardInterrupt:
        if frames:
            temporary=out/'checkpoint.tmp.npz'
            save_trajectory(temporary,frames,metadata)
            temporary.replace(out/'run_trajectory.npz')
        print('Interrupted; complete saved checkpoints retained.',flush=True)
        return

    # Median bulge positions resist distant-particle contamination.
    # This is a location diagnostic, not a bound-core finder.
    rows = []
    for time, particles in frames:
        centers = [
            np.median(
                particles.x[
                    (particles.component == 1)
                    & (particles.galaxy == galaxy_id)
                ],
                axis=0,
            )
            for galaxy_id in (0, 1)
        ]
        separation = float(np.linalg.norm(centers[1] - centers[0]))
        rows.append([time, separation])

    with (out / "core_separation.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "time_code",
            "bulge_median_separation_code",
        ])
        writer.writerows(rows)

    print("Saved:", out / "run_trajectory.npz")
    print("Open it in the GUI's Replay a saved run tab.")


if __name__ == "__main__":
    main()