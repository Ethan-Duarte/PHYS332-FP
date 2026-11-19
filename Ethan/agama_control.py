import csv
import json
from datetime import datetime
from pathlib import Path

import numpy as np

from nbody.integrate import evolve
from nbody.io import save, save_trajectory
from nbody.particles import Particles, combine


def metrics(p):
    x = p.x - np.average(p.x, axis=0, weights=p.m)
    v = p.v - np.average(p.v, axis=0, weights=p.m)
    cov = (x * p.m[:, None]).T @ x / p.m.sum()
    values, vectors = np.linalg.eigh(cov)
    normal = vectors[:, 0]
    z = x @ normal
    radius = np.sqrt(
        np.maximum(np.sum(x*x, axis=1) - z*z, 0)
    )
    order = np.argsort(radius)
    r50 = np.interp(
        0.5, np.cumsum(p.m[order]) / p.m.sum(), radius[order]
    )
    return {
        "r50": float(r50),
        "height": float(np.sqrt(max(values[0], 0))),
        "sigma_z": float(
            np.sqrt(np.average((v @ normal)**2, weights=p.m))
        ),
    }


def main():
    try:
        import agama
    except ImportError as exc:
        raise SystemExit(
            "AGAMA must import in this Python environment."
        ) from exc

    root = Path(__file__).resolve().parent
    out = root / "results" / (
        "agama_control_" + datetime.now().strftime("%Y%m%d_%H%M%S")
    )
    out.mkdir(parents=True, exist_ok=False)

    cfg = json.loads((root / "presets/mice.json").read_text())
    cfg.update(Q=1.5, seed=42)
    d, b, h = [cfg[k] for k in ("disk", "bulge", "halo")]
    rd, zd = float(d["scale"]), float(d["height"])

    # Default dimensionless units, G=1. Do not call setUnits().
    agama.setRandomSeed(42)
    densities = [
        agama.Density(
            type="Disk", mass=d["mass"],
            scaleRadius=rd, scaleHeight=-zd/2
        ),
        agama.Density(
            type="Spheroid", mass=b["mass"],
            scaleRadius=b["scale"], alpha=1, beta=4, gamma=1
        ),
        agama.Density(
            type="Spheroid", mass=h["mass"],
            scaleRadius=h["scale"], alpha=1, beta=4, gamma=0
        ),
    ]

    sg = dict(
        rminSph=min(b["scale"], h["scale"], rd)/40,
        rmaxSph=max(
            2*h.get("cutoff", 0), 30*h["scale"],
            15*b["scale"], 15*rd
        ),
        sizeRadialSph=50, lmaxAngularSph=6
    )
    dg = dict(
        RminCyl=rd/30, RmaxCyl=12*rd,
        sizeRadialCyl=30, zminCyl=zd/4,
        zmaxCyl=10*rd, sizeVerticalCyl=26
    )

    model = agama.SelfConsistentModel(**sg, **dg)
    model.components = [
        agama.Component(density=rho, disklike=(i == 0))
        for i, rho in enumerate(densities)
    ]
    print("Constructing the initial potential...", flush=True)
    model.iterate()
    pot0 = model.potential

    def omega2(r):
        return -pot0.force(r, 0., 0.)[0] / r

    step = rd*1e-3
    kappa2 = (
        rd*(omega2(rd+step)-omega2(rd-step))/(2*step)
        + 4*omega2(rd)
    )
    if not np.isfinite(kappa2) or kappa2 <= 0:
        raise RuntimeError("Invalid radial epicycle frequency")

    sigma = (
        1.5*3.36*d["mass"]
        / (2*np.pi*rd**2*np.e)
        / np.sqrt(kappa2)
    )
    params = dict(
        type="QuasiIsothermal", potential=pot0,
        Sigma0=d["mass"]/(2*np.pi*rd**2),
        Rdisk=rd, Hdisk=zd/2,
        sigmar0=sigma*np.exp(.5),
        Rsigmar=2*rd, sigmamin=.01*sigma
    )
    df0 = agama.DistributionFunction(**params)
    params["Sigma0"] *= d["mass"] / df0.totalMass()

    dfs = [agama.DistributionFunction(**params)] + [
        agama.DistributionFunction(
            type="QuasiSpherical",
            potential=pot0, density=rho
        )
        for rho in densities[1:]
    ]
    model.components = [
        agama.Component(
            df=df, disklike=(i == 0),
            **(dg if i == 0 else sg)
        )
        for i, df in enumerate(dfs)
    ]

    # Compare radial/vertical forces and component densities.
    # Convergence here applies to these sampled disk-region points.
    probes = np.array([
        [r, 0., z]
        for r in np.geomspace(.25*rd, 6*rd, 20)
        for z in (0., .5*zd, zd, 2*zd)
    ])

    previous = None
    history = []
    consecutive = 0

    for iteration in range(1, 21):
        model.iterate()
        force = np.asarray(model.potential.force(probes))
        fields = [force[:, 0], force[:, 2]] + [
            np.asarray(comp.density.density(probes))
            for comp in model.components
        ]

        errors = []
        if previous is not None:
            for new, old in zip(fields, previous):
                floor = max(
                    float(np.max(np.abs(old)))*.01, 1e-12
                )
                errors.append(float(np.max(
                    np.abs(new-old)
                    / np.maximum(np.abs(old), floor)
                )))

        change = max(errors) if errors else None
        history.append({
            "iteration": iteration,
            "max_fractional_change": change,
            "field_changes": errors
        })
        (out / "convergence.json").write_text(
            json.dumps(history, indent=2)
        )
        print(
            f"Iteration {iteration}: "
            f"maximum sampled change = {change}",
            flush=True
        )

        consecutive = (
            consecutive+1
            if change is not None and change < .005
            else 0
        )
        previous = [field.copy() for field in fields]

        if iteration >= 6 and consecutive >= 2:
            break
    else:
        raise SystemExit(
            "Convergence threshold not met. Send this file: "
            + str(out / "convergence.json")
        )

    model.potential.export(str(out / "smooth_potential.ini"))

    def sample(component, count, seed):
        agama.setRandomSeed(seed)
        xv, masses = agama.GalaxyModel(
            potential=model.potential,
            df=dfs[component], af=model.af
        ).sample(count)
        xv = np.asarray(xv, float)
        masses = np.asarray(masses, float)

        # Fix physical masses across the controlled pair.
        masses *= (
            [d["mass"], b["mass"], h["mass"]][component]
            / masses.sum()
        )
        return Particles(
            xv[:, :3].copy(), xv[:, 3:].copy(), masses,
            np.full(count, component, dtype=int),
            np.zeros(count, dtype=int)
        )

    print("Sampling the shared disk and bulge...", flush=True)
    disk = sample(0, 4800, 42)
    bulge = sample(1, 1200, 43)

    metadata = {
        "backend": "agama",
        "config": cfg,
        "convergence": history,
        "status": "candidate; stability remains unverified",
        "frame": "smooth model origin, without COM recentering",
        "mass_normalization": (
            "sampled masses normalized to target component masses"
        )
    }

    for count in (6000, 24000):
        print(
            f"Sampling halo with {count} particles...",
            flush=True
        )
        halo = sample(2, count, 100+count)
        config = {
            **cfg,
            "particles": dict(
                disk=4800, bulge=1200, halo=count
            )
        }
        save(
            out / f"candidate_halo{count}.npz",
            combine(disk, bulge, halo),
            {**metadata, "config": config}
        )

    print(
        "Evolving disk tracers in the fixed smooth potential...",
        flush=True
    )

    def fixed_force(x, masses):
        return np.asarray(model.potential.force(x))

    # The smooth field includes disk, bulge and halo.
    # Tracers add no live particle self-gravity in this control.
    frames, rows = [], []
    for t, p in evolve(disk, fixed_force, 1/1024, 512, 32):
        frames.append((t, p))
        rows.append(dict(time=float(t), **metrics(p)))
        print(f"Smooth control t={t:.5f}/0.5", flush=True)

    save_trajectory(
        out / "smooth_control_trajectory.npz",
        frames,
        dict(
            setup="Isolated galaxy",
            force="fixed smooth potential",
            disk_tracers_only=True,
            dt=1/1024, steps=512,
            caveat=(
                "Diagnostic control, "
                "not a live self-gravitating galaxy"
            )
        )
    )

    changes = {
        key: rows[-1][key]/rows[0][key]-1
        for key in ("r50", "height", "sigma_z")
    }
    summary = {
        "relative_changes": changes,
        "time_span": .5,
        "caveat": (
            "One realization and timestep. "
            "No automatic validation verdict."
        )
    }
    (out / "smooth_control_summary.json").write_text(
        json.dumps(summary, indent=2)
    )
    with (out / "smooth_metrics.csv").open(
        "w", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0])
        )
        writer.writeheader()
        writer.writerows(rows)

    print(json.dumps(summary, indent=2))
    print("Results:", out)


if __name__ == "__main__":
    main()