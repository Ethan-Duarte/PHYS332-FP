import argparse
import time
import csv
import json
from datetime import datetime
from pathlib import Path

import numpy as np

from nbody.galaxies.softening import PlummerConvolution, soften_expansion
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
    result = {
        "r50": float(r50),
        "height": float(np.sqrt(max(values[0], 0))),
        "sigma_z": float(
            np.sqrt(np.average((v @ normal)**2, weights=p.m))
        ),
    }
    # Remove a separately fitted plane in each annulus to distinguish local
    # thickening from a coherent tilt/warp.  Sparse annuli are reported as null.
    for inner, outer in ((0.,1.),(1.,2.),(2.,3.),(3.,5.)):
        mask = (radius >= inner/12) & (radius < outer/12)
        key = f'height_Rd_{inner:g}_{outer:g}'
        if mask.sum() < 30:
            result[key] = None
            continue
        local_x = x[mask]-np.average(x[mask],axis=0,weights=p.m[mask])
        covariance = (local_x*p.m[mask,None]).T@local_x/p.m[mask].sum()
        result[key] = float(np.sqrt(max(np.linalg.eigvalsh(covariance)[0],0)))
    angle = np.arctan2(x[:,1],x[:,0])
    inner_disk = (radius > .5/12) & (radius < 3/12)
    result['A2_inner_disk'] = float(abs(np.average(np.exp(2j*angle[inner_disk]),weights=p.m[inner_disk])))
    result['tilt_degrees'] = float(np.degrees(np.arccos(np.clip(abs(normal[2]),0,1))))
    return result


def make_distribution_functions(agama, pot0, densities, cfg):
    d = cfg['disk']
    rd, zd = float(d['scale']), float(d['height'])
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
        cfg["Q"]*3.36*d["mass"]
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
    return dfs


def main():
    parser = argparse.ArgumentParser(description="Experimental Plummer-compatible AGAMA equilibrium; no stability claim")
    parser.add_argument('--epsilon', type=float, default=.00625)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--steps', type=int, default=32)
    parser.add_argument('--dt', type=float, default=1/1024)
    parser.add_argument('--halo', type=int, default=6000)
    parser.add_argument('--disk', type=int, default=4800)
    parser.add_argument('--bulge', type=int, default=1200)
    parser.add_argument('--q', type=float, default=1.5)
    parser.add_argument('--iterations', type=int, default=20)
    parser.add_argument('--build-only', action='store_true')
    parser.add_argument('--interpolated-actions', action='store_true', help='Faster SCM action approximation; exact-AF density check still required')
    parser.add_argument('--screen', type=Path, help='Previously built candidate snapshot; skip reconstruction')
    args = parser.parse_args()
    if args.epsilon <= 0 or args.dt <= 0 or args.steps < 1 or args.halo < 1000 or args.disk < 1000 or args.bulge < 1000 or args.threads < 1:
        parser.error('positive epsilon/dt/steps and >=1000 halo particles required')
    try:
        import agama
    except ImportError as exc:
        raise SystemExit(
            "AGAMA must import in this Python environment."
        ) from exc

    agama.setNumThreads(args.threads)
    if args.screen is not None:
        return run_screen(args.screen, args)

    root = Path(__file__).resolve().parent
    out = root / "results" / (
        "agama_softened_" + datetime.now().strftime("%Y%m%d_%H%M%S")
    )
    out.mkdir(parents=True, exist_ok=False)

    cfg = json.loads((root / "presets/mice.json").read_text())
    cfg.update(Q=args.q, seed=args.seed, softening=args.epsilon)
    d, b, h = [cfg[k] for k in ("disk", "bulge", "halo")]
    rd, zd = float(d["scale"]), float(d["height"])

    # Default dimensionless units, G=1. Do not call setUnits().
    agama.setRandomSeed(args.seed)
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

    model = agama.SelfConsistentModel(useActionInterpolation=args.interpolated_actions, **sg, **dg)
    model.components = [
        agama.Component(density=rho, disklike=(i == 0))
        for i, rho in enumerate(densities)
    ]
    print("Constructing the initial potential...", flush=True)
    model.iterate()
    print("Building the initial Plummer-smoothed field...", flush=True)
    initial_newtonian = model.potential
    early_points = np.array([[rr*rd,0.,zz*zd] for rr in (.5,1.,2.,3.) for zz in (.5,1.,2.)])
    fine_initial = PlummerConvolution(initial_newtonian,args.epsilon,(48,24,48)).force(early_points)
    quad_order = None
    for candidate_order in ((16,8,16),(24,12,24),(36,18,36)):
        trial = PlummerConvolution(initial_newtonian,args.epsilon,candidate_order).force(early_points)
        mismatch = max(float(np.max(np.abs(trial[:,j]-fine_initial[:,j])/np.maximum(np.abs(fine_initial[:,j]),.02*np.max(np.abs(fine_initial[:,j]))))) for j in (0,2))
        if mismatch < .01:
            quad_order = candidate_order
            break
    if quad_order is None:
        raise SystemExit('Initial quadrature refinement did not converge below 1%; candidate NOT constructed')
    print(f'Selected quadrature {quad_order}; initial force refinement error {mismatch:.4%}',flush=True)
    pot0 = soften_expansion(model.potential, args.epsilon, sg, dg, order=quad_order)
    model.potential = pot0

    dfs = make_distribution_functions(agama,pot0,densities,cfg)
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

    for iteration in range(1, args.iterations+1):
        started = time.perf_counter()
        model.iterate()
        newtonian = model.potential
        model.potential = soften_expansion(newtonian, args.epsilon, sg, dg, order=quad_order)
        print(f"Softened iteration {iteration}: {time.perf_counter()-started:.1f} seconds", flush=True)
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
    # Independent checks of quadrature AND its potential expansion on disk probes.
    # Exclude z=0 for vertical relative errors and use a floor near weak fields.
    print("Checking quadrature and potential expansion...", flush=True)
    coarse = PlummerConvolution(newtonian, args.epsilon, quad_order)
    fine = PlummerConvolution(newtonian, args.epsilon, (48, 24, 48))
    check_points = np.array([[rr*rd, 0., zz*zd] for rr in (.5, 1., 2., 3.) for zz in (.5, 1., 2.)])
    expected = fine.force(check_points)
    approximate = coarse.force(check_points)
    expanded = model.potential.force(check_points)
    errors = {}
    for label, measured in [('quadrature', approximate), ('expansion', expanded)]:
        for axis, index in [('radial', 0), ('vertical', 2)]:
            denominator = np.maximum(np.abs(expected[:, index]), .02*np.max(np.abs(expected[:, index])))
            errors[f'{label}_{axis}_max_fraction'] = float(np.max(np.abs(measured[:, index]-expected[:, index])/denominator))
    (out/'numerical_force_checks.json').write_text(json.dumps(errors, indent=2))
    print(json.dumps(errors, indent=2), flush=True)
    if max(errors.values()) > .02:
        raise SystemExit(f'Numerical force error exceeds provisional 2% gate. Candidate NOT accepted. Send {out}/numerical_force_checks.json')
    final_af = agama.ActionFinder(model.potential)
    print('Checking DF densities using a fresh non-interpolated action finder...',flush=True)
    exact_density_errors = []
    for index,df in enumerate(dfs):
        density_exact = np.asarray(agama.GalaxyModel(potential=model.potential,df=df,af=final_af).moments(check_points,dens=True,vel=False,vel2=False))
        density_grid = np.asarray(model.components[index].density.density(check_points))
        denominator = np.maximum(np.abs(density_exact), .02*np.max(np.abs(density_exact)))
        exact_density_errors.append(float(np.max(np.abs(density_grid-density_exact)/denominator)))
    (out/'density_closure_checks.json').write_text(json.dumps(exact_density_errors,indent=2))
    print('Exact-AF density discrepancies:',exact_density_errors,flush=True)
    if max(exact_density_errors) > .03:
        raise SystemExit('DF/grid density mismatch >3%; candidate NOT accepted')
    df_masses = [float(df.totalMass()) for df in dfs]
    density_masses = [float(comp.density.totalMass()) for comp in model.components]
    (out/'mass_checks.json').write_text(json.dumps(dict(df_masses=df_masses, density_masses=density_masses), indent=2))
    if max(abs(a/b-1) for a,b in zip(density_masses,df_masses)) > .02:
        raise SystemExit('DF and grid-density masses disagree by >2%; candidate NOT accepted')

    def sample(component, count, seed):
        agama.setRandomSeed(seed)
        xv, masses = agama.GalaxyModel(
            potential=model.potential,
            df=dfs[component], af=final_af
        ).sample(count)
        xv = np.asarray(xv, float)
        masses = np.asarray(masses, float)

        # Fix physical masses across the controlled pair.
        masses *= (
            df_masses[component]
            / masses.sum()
        )
        return Particles(
            xv[:, :3].copy(), xv[:, 3:].copy(), masses,
            np.full(count, component, dtype=int),
            np.zeros(count, dtype=int)
        )

    print("Sampling the shared disk and bulge...", flush=True)
    disk = sample(0, args.disk, args.seed)
    bulge = sample(1, args.bulge, args.seed+1)

    metadata = {
        "backend": "agama",
        "config": cfg,
        "convergence": history,
        "status": "softening-compatible candidate; live stability remains unverified",
        "epsilon": args.epsilon,
        "scm_interpolated_actions": args.interpolated_actions,
        "density_closure_max_errors": exact_density_errors,
        "quadrature_order": list(quad_order),
        "numerical_force_checks": errors,
        "frame": "smooth model origin, without COM recentering",
        "mass_normalization": (
            "sampled masses normalized to DF masses, not analytic target masses"
        )
    }

    for count in (args.halo,):
        print(
            f"Sampling halo with {count} particles...",
            flush=True
        )
        halo = sample(2, count, args.seed+100+count)
        config = {
            **cfg,
            "particles": dict(
                disk=args.disk, bulge=args.bulge, halo=count
            )
        }
        save(
            out / f"candidate_halo{count}.npz",
            combine(disk, bulge, halo),
            {**metadata, "config": config}
        )

    print(f"Candidate saved: {out}", flush=True)
    if not args.build_only:
        run_screen(out / f"candidate_halo{args.halo}.npz", args)


def run_screen(snapshot, args):
    import agama
    from nbody.io import load
    import numba
    if args.threads > 1:
        from nbody.gravity.parallel_direct import acceleration
        numba.set_num_threads(min(args.threads,numba.config.NUMBA_NUM_THREADS))
    else:
        from nbody.gravity.fast_direct import acceleration
    from numba import njit

    @njit(cache=False)
    def total_energy(x,v,m,eps):
        kinetic = .5*np.sum(m[:,None]*v*v)
        potential_energy = 0.
        for i in range(len(m)):
            for j in range(i):
                delta = x[j]-x[i]
                potential_energy -= m[i]*m[j]/np.sqrt(np.sum(delta*delta)+eps*eps)
        return kinetic+potential_energy
    source, meta = load(snapshot)
    if meta.get('backend') != 'agama' or 'numerical_force_checks' not in meta:
        raise SystemExit('Use a candidate produced by this softened-equilibrium script')
    eps = float(meta['epsilon'])
    if not np.isclose(eps, args.epsilon, rtol=0, atol=1e-15):
        raise SystemExit(f"Snapshot epsilon={eps}; specify --epsilon {eps} (no silent force changes)")
    potential = agama.Potential(str(snapshot.parent/'smooth_potential.ini'))
    output = snapshot.parent / ('screen_' + datetime.now().strftime('%Y%m%d_%H%M%S'))
    output.mkdir()
    force_audit = []
    for rr in (.5,1.,2.,3.):
        for zz in (.5,1.,2.):
            radius, height = rr/12, zz*.005
            probe = np.array([[radius,0.,height],[radius,0.,-height]])
            smooth = potential.force(probe)
            restoring = -(smooth[0,2]-smooth[1,2])/2
            samples = []
            for phi in np.arange(48)*2*np.pi/48:
                planar = (source.x[:,0]-radius*np.cos(phi))**2+(source.x[:,1]-radius*np.sin(phi))**2
                up, down = source.x[:,2]-height, source.x[:,2]+height
                component_force = .5*source.m*(down/(planar+down*down+eps*eps)**1.5-up/(planar+up*up+eps*eps)**1.5)
                samples.append([component_force[source.component==i].sum()/restoring for i in range(3)])
            samples = np.asarray(samples)
            force_audit.append(dict(R_over_Rd=rr,z_over_zd=zz,mean_support_ratio=float(samples.sum(axis=1).mean()),azimuth_std_over_smooth=float(samples.sum(axis=1).std()),component_std_over_smooth=samples.std(axis=0).tolist()))
    (output/'initial_force_audit.json').write_text(json.dumps(force_audit,indent=2))
    mask = source.component == 0
    disk = Particles(source.x[mask].copy(), source.v[mask].copy(), source.m[mask].copy(), source.component[mask], source.galaxy[mask])
    phases = [('smooth', disk, lambda x,m: np.asarray(potential.force(x))),
              ('live', source, lambda x,m: acceleration(x,m,eps=eps,G=1.))]
    summaries = {}
    for label, initial, force in phases:
        print(f"Warming up {label} force...", flush=True)
        force(initial.x[:16], initial.m[:16])
        start = time.perf_counter()
        force(initial.x,initial.m)
        force_seconds = time.perf_counter()-start
        save_every = max(8,8*int(np.ceil(args.steps/256)))
        print(f"{label}: estimated force time {(args.steps+1)*force_seconds:.1f} s; checkpoints every {save_every} steps", flush=True)
        initial_energy = total_energy(initial.x,initial.v,initial.m,eps) if label == 'live' else None
        frames, rows = [], []
        selected = initial.component == 0
        start = time.perf_counter()
        try:
            for t, state in evolve(initial, force, args.dt, args.steps, save_every):
                frames.append((t,state))
                tracer = Particles(state.x[selected],state.v[selected],state.m[selected])
                row = dict(time=float(t), **metrics(tracer))
                rows.append(row)
                # Replace only after a complete compressed archive is written.
                temporary = output / f'{label}_checkpoint.tmp.npz'
                save_trajectory(temporary,frames,dict(force=label,epsilon=eps,dt=args.dt,requested_steps=args.steps,source=str(snapshot),status='screening only'))
                temporary.replace(output/f'{label}_trajectory.npz')
                (output/f'{label}_metrics.json').write_text(json.dumps(rows,indent=2))
                print(f"{label}: t={t:.6f}/{args.steps*args.dt:.6f}, elapsed {time.perf_counter()-start:.1f}s",flush=True)
        except KeyboardInterrupt:
            print(f"Interrupted; completed checkpoints retained in {output}",flush=True)
            return
        summaries[label] = dict(relative_changes={k:rows[-1][k]/rows[0][k]-1 for k in ('r50','height','sigma_z')}, duration=args.steps*args.dt,steps=args.steps,dt=args.dt,wall_seconds=time.perf_counter()-start)
    summaries['live']['relative_energy_change'] = float((total_energy(state.x,state.v,state.m,eps)-initial_energy)/abs(initial_energy))
    summaries['status'] = 'Diagnostic screening only; no automatic stability verdict'
    (output/'summary.json').write_text(json.dumps(summaries,indent=2))
    print(json.dumps(summaries,indent=2),flush=True)
    print(f'Results: {output}',flush=True)


if __name__ == "__main__":
    main()