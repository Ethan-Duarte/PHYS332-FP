"""Benchmark and accuracy comparison for four N-body gravity solvers.

Solvers compared (all expose ``acceleration(x, m, eps=..., [theta=...], G=...)``):

    direct            vectorized NumPy O(N^2)          (reference for small N)
    fast_direct       Numba O(N^2)                     (reference for larger N)
    barnes_hut        pure-Python octree + traversal
    fast_barnes_hut   Python tree build, Numba traversal

Three experiments, each written to CSV in the output directory:

    force_error.csv   per-solver relative force error vs. exact direct summation
    timing.csv        wall-clock time vs. N (min / median of repeats)
    tree_split.csv    Barnes-Hut time split into tree build vs. traversal
    conservation.csv  energy / momentum / angular momentum / virial ratio vs. time

plus ``sanity.json`` (cross-checks between solver pairs) and ``meta.json``
(machine, library versions, seed, settings).

Usage (the solver files live in an importable package, so pass its dotted name):

    python compare.py --package myproject.gravity --out results
    python compare.py --package myproject.gravity --quick

or from a notebook:

    import compare
    compare.run_all("myproject.gravity", out="results", quick=True)
"""
from __future__ import annotations

import argparse
import csv
import importlib
import json
import os
import platform
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

# --------------------------------------------------------------------------
# Solver loading
# --------------------------------------------------------------------------
SOLVER_MODULES = {
    "direct": "direct",
    "fast_direct": "fast_direct",
    "barnes_hut": "barnes_hut",
    "fast_barnes_hut": "fast_barnes_hut",
}
TREE_SOLVERS = {"barnes_hut", "fast_barnes_hut"}

# Largest N each solver is allowed to run in the timing sweep. Pure-Python
# Barnes-Hut is very slow; NumPy direct needs ~ N^2 * 3 * 8 bytes per temporary.
DEFAULT_MAX_N = {
    "direct": 6_000,
    "fast_direct": 100_000,
    "barnes_hut": 4_000,
    "fast_barnes_hut": 100_000,
}


def load_solvers(package: str, only=None) -> dict:
    """Import whichever solver modules are available. Missing ones (e.g. no numba)
    are skipped with a warning rather than aborting the whole benchmark."""
    modules = {}
    for name, mod in SOLVER_MODULES.items():
        if only and name not in only:
            continue
        try:
            modules[name] = importlib.import_module(f"{package}.{mod}")
        except Exception as exc:  # ImportError, numba missing, ...
            print(f"[warn] could not import {package}.{mod}: {type(exc).__name__}: {exc}",
                  file=sys.stderr)
    if not modules:
        raise ImportError(f"No solver modules found in package {package!r}")
    return modules


def make_callable(name: str, module, eps: float, theta: float, G: float = 1.0):
    """Return f(x, m) -> accelerations with eps/theta/G bound."""
    fn = module.acceleration
    if name in TREE_SOLVERS:
        return lambda x, m: fn(x, m, eps=eps, theta=theta, G=G)
    return lambda x, m: fn(x, m, eps=eps, G=G)


# --------------------------------------------------------------------------
# Initial conditions
# --------------------------------------------------------------------------
def plummer_sphere(n: int, seed: int = 0, a: float = 1.0, total_mass: float = 1.0, G: float = 1.0):
    """Plummer model in approximate virial equilibrium (Aarseth, Henon & Wielen 1974).
    Returned in the centre-of-mass frame. Equal-mass particles."""
    rng = np.random.default_rng(seed)
    # radii from the inverse cumulative mass; clip the extreme tail
    u = rng.uniform(1e-6, 1 - 1e-6, n)
    r = a / np.sqrt(u ** (-2.0 / 3.0) - 1.0)

    def isotropic(size, scale):
        cos_t = rng.uniform(-1, 1, size)
        phi = rng.uniform(0, 2 * np.pi, size)
        sin_t = np.sqrt(1 - cos_t**2)
        return scale[:, None] * np.stack([sin_t * np.cos(phi), sin_t * np.sin(phi), cos_t], axis=1)

    x = isotropic(n, r)

    # speed: rejection sampling of g(q) = q^2 (1 - q^2)^3.5
    q = np.empty(n)
    filled = 0
    while filled < n:
        qs = rng.uniform(0, 1, 2 * (n - filled) + 16)
        ys = rng.uniform(0, 0.1, qs.size)
        ok = qs[ys < qs**2 * (1 - qs**2) ** 3.5]
        take = ok[: n - filled]
        q[filled: filled + take.size] = take
        filled += take.size
    v_esc = np.sqrt(2 * G * total_mass) * (r**2 + a**2) ** -0.25
    v = isotropic(n, q * v_esc)

    m = np.full(n, total_mass / n)
    x -= np.average(x, weights=m, axis=0)
    v -= np.average(v, weights=m, axis=0)
    return x, v, m


# --------------------------------------------------------------------------
# Diagnostics
# --------------------------------------------------------------------------
def kinetic_energy(v, m):
    return 0.5 * float(np.sum(m * np.sum(v * v, axis=1)))


def potential_energy(x, m, eps, G=1.0, max_block_elems=20_000_000):
    """Softened potential U = -G sum_{i<j} m_i m_j / sqrt(r_ij^2 + eps^2).

    Computed in row blocks so memory stays O(block * N) instead of O(N^2 * 3).
    This is the potential whose gradient is the Plummer-softened force used by
    all four solvers. It is always evaluated exactly (direct), independent of
    the solver under test, and is NOT included in solver timings."""
    n = len(m)
    block = max(1, max_block_elems // n)
    total = 0.0
    for i0 in range(0, n, block):
        i1 = min(n, i0 + block)
        d2 = np.zeros((i1 - i0, n))
        for k in range(3):
            dk = x[i0:i1, None, k] - x[None, :, k]
            d2 += dk * dk
        d2 += eps**2
        inv = 1.0 / np.sqrt(d2)
        idx = np.arange(i0, i1)
        inv[idx - i0, idx] = 0.0  # remove self-interaction
        total += float(np.sum(m[i0:i1, None] * m[None, :] * inv))
    return -0.5 * G * total  # 0.5 because each pair counted twice


def diagnostics(x, v, m, eps, G=1.0):
    ke = kinetic_energy(v, m)
    pe = potential_energy(x, m, eps, G)
    p = np.sum(m[:, None] * v, axis=0)
    L = np.sum(m[:, None] * np.cross(x, v), axis=0)
    com = np.average(x, weights=m, axis=0)
    r = np.linalg.norm(x - com, axis=1)
    return {
        "ke": ke,
        "pe": pe,
        "energy": ke + pe,
        "virial": 2 * ke / abs(pe) if pe != 0 else float("nan"),
        "p_norm": float(np.linalg.norm(p)),
        "L_norm": float(np.linalg.norm(L)),
        "r_half": float(np.median(r)),  # half-mass radius for equal masses
    }


# --------------------------------------------------------------------------
# Force error
# --------------------------------------------------------------------------
def relative_force_error(a, a_ref):
    """Per-particle |a - a_ref| / |a_ref|."""
    denom = np.linalg.norm(a_ref, axis=1)
    denom = np.where(denom > 0, denom, np.finfo(float).tiny)
    return np.linalg.norm(a - a_ref, axis=1) / denom


def error_stats(err):
    return {
        "median": float(np.median(err)),
        "mean": float(np.mean(err)),
        "p99": float(np.percentile(err, 99)),
        "max": float(np.max(err)),
    }


def pick_reference(modules):
    """Exact reference: Numba direct if available, else NumPy direct."""
    return "fast_direct" if "fast_direct" in modules else "direct"


# --------------------------------------------------------------------------
# Timing
# --------------------------------------------------------------------------
def warm_up(modules, eps, theta, n=256, seed=123):
    """Trigger JIT compilation / cache loading outside the timed region."""
    x, _, m = plummer_sphere(n, seed)
    for name, mod in modules.items():
        make_callable(name, mod, eps, theta)(x, m)


def time_call(fn, x, m, repeats=5, slow_s=5.0):
    """Best/median of `repeats` calls. If a single call takes longer than `slow_s`
    seconds, stop early (repeating an 80 s call five times adds nothing)."""
    times = []
    out = None
    for _ in range(repeats):
        t0 = time.perf_counter()
        out = fn(x, m)
        times.append(time.perf_counter() - t0)
        if times[-1] > slow_s:
            break
    return {"min": min(times), "median": float(np.median(times)), "repeats": len(times)}, out


def estimate_crossover(timing_rows, slow="fast_direct", fast="fast_barnes_hut", n_fit=4):
    """Fit t = c * N^p to the largest `n_fit` points of each solver (log-log) and
    solve for the N where the curves intersect. Works even when the crossover lies
    beyond the measured range (then `extrapolated` is True). Returns None if the
    fitted slopes do not cross (the tree is not scaling better than direct)."""
    import math
    fits, nmax = {}, {}
    for name in (slow, fast):
        pts = sorted((r["n"], r["time_min_s"]) for r in timing_rows if r["solver"] == name)
        if len(pts) < 2:
            return None
        pts = pts[-n_fit:]
        lx = np.log([p[0] for p in pts]); ly = np.log([p[1] for p in pts])
        slope, icpt = np.polyfit(lx, ly, 1)
        fits[name] = (slope, icpt)
        nmax[name] = pts[-1][0]
    (p1, c1), (p2, c2) = fits[slow], fits[fast]
    if p1 <= p2:
        return None
    n_cross = math.exp((c2 - c1) / (p1 - p2))
    return {"n_cross": n_cross, "slopes": {slow: p1, fast: p2},
            "extrapolated": n_cross > min(nmax.values())}


# --------------------------------------------------------------------------
# Time integration (kick-drift-kick leapfrog, one force evaluation per step)
# --------------------------------------------------------------------------
def leapfrog_run(accel, x0, v0, m, eps, dt, steps, diag_every, G=1.0):
    x, v = x0.copy(), v0.copy()
    a = accel(x, m)
    rows = []
    force_time = 0.0

    def record(step):
        d = diagnostics(x, v, m, eps, G)
        d["step"] = step
        d["t"] = step * dt
        rows.append(d)

    record(0)
    wall0 = time.perf_counter()
    for step in range(1, steps + 1):
        v += 0.5 * dt * a
        x += dt * v
        t0 = time.perf_counter()
        a = accel(x, m)
        force_time += time.perf_counter() - t0
        v += 0.5 * dt * a
        if step % diag_every == 0 or step == steps:
            record(step)
    wall = time.perf_counter() - wall0

    e0 = rows[0]["energy"]
    p0 = rows[0]["p_norm"]
    L0 = rows[0]["L_norm"]
    r0 = rows[0]["r_half"]
    for r in rows:
        r["dE_rel"] = (r["energy"] - e0) / abs(e0)
        r["p_drift"] = r["p_norm"] - p0  # initial momentum ~ 0 in COM frame, so absolute
        r["dL_rel"] = (r["L_norm"] - L0) / L0 if L0 > 0 else float("nan")
        r["dr_half_rel"] = (r["r_half"] - r0) / r0
    return rows, {"wall_total": wall, "wall_force": force_time}


# --------------------------------------------------------------------------
# Experiments
# --------------------------------------------------------------------------
@dataclass
class Config:
    package: str
    out: str = "results"
    seed: int = 42
    eps: float = 0.05
    G: float = 1.0
    thetas: tuple = (0.3, 0.5, 0.7, 1.0)
    theta_default: float = 0.5
    # force-error experiment
    n_force: int = 2000
    # timing sweep
    n_sweep: tuple = (250, 500, 1000, 2000, 4000, 8000, 16000, 32000)
    repeats: int = 5
    solvers: tuple | None = None   # None = all four
    max_n: dict = field(default_factory=lambda: dict(DEFAULT_MAX_N))
    # conservation run
    n_cons: int = 500
    dt: float = 0.01
    steps: int = 300
    diag_every: int = 10

    @classmethod
    def quick(cls, package, out="results"):
        return cls(package=package, out=out, thetas=(0.5, 1.0), n_force=500,
                   n_sweep=(250, 500, 1000, 2000), repeats=3,
                   n_cons=200, steps=100, diag_every=10)


def _write_csv(path, rows):
    if not rows:
        return
    keys = list(rows[0].keys())
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def experiment_force_error(modules, cfg: Config):
    ref_name = pick_reference(modules)
    x, _, m = plummer_sphere(cfg.n_force, cfg.seed)
    a_ref = make_callable(ref_name, modules[ref_name], cfg.eps, cfg.theta_default, cfg.G)(x, m)
    rows = []
    for name, mod in modules.items():
        thetas = cfg.thetas if name in TREE_SOLVERS else (float("nan"),)
        for theta in thetas:
            fn = make_callable(name, mod, cfg.eps, theta if theta == theta else cfg.theta_default, cfg.G)
            t, a = time_call(fn, x, m, repeats=max(1, min(3, cfg.repeats)))
            row = {"solver": name, "theta": theta, "n": cfg.n_force, "reference": ref_name,
                   "time_min_s": t["min"]}
            row.update({f"err_{k}": v for k, v in error_stats(relative_force_error(a, a_ref)).items()})
            rows.append(row)
            print(f"  force error  {name:16s} theta={theta:<4}  median={row['err_median']:.2e}  "
                  f"p99={row['err_p99']:.2e}  t={t['min']:.3f}s")
    return rows


def experiment_sanity(modules, cfg: Config, n=500):
    """Pairs that should agree to roundoff."""
    x, _, m = plummer_sphere(n, cfg.seed + 1)
    out = {"n": n}

    def acc(name, theta=None):
        return make_callable(name, modules[name], cfg.eps, theta or cfg.theta_default, cfg.G)(x, m)

    if "direct" in modules and "fast_direct" in modules:
        out["direct_vs_fast_direct_max_rel"] = float(
            relative_force_error(acc("fast_direct"), acc("direct")).max())
    if "barnes_hut" in modules and "fast_barnes_hut" in modules:
        out["bh_vs_fast_bh_max_rel"] = float(
            relative_force_error(acc("fast_barnes_hut"), acc("barnes_hut")).max())
    # tiny-theta limit: tree code must converge to direct summation
    ref = acc(pick_reference(modules))
    for name in TREE_SOLVERS & modules.keys():
        out[f"{name}_theta0.05_vs_direct_max_rel"] = float(
            relative_force_error(acc(name, 0.05), ref).max())
    return out


def experiment_timing(modules, cfg: Config, path=None):
    """Rows are written to `path` after every measurement, so a crash or kill
    (e.g. out-of-memory) does not lose the data gathered so far."""
    rows = []
    for n in cfg.n_sweep:
        x, _, m = plummer_sphere(n, cfg.seed)
        for name, mod in modules.items():
            if n > cfg.max_n.get(name, np.inf):
                continue
            fn = make_callable(name, mod, cfg.eps, cfg.theta_default, cfg.G)
            try:
                t, _ = time_call(fn, x, m, repeats=cfg.repeats)
            except MemoryError:
                print(f"  timing  {name:16s} N={n:>6}  MemoryError")
                continue
            rows.append({"solver": name, "n": n, "theta": cfg.theta_default if name in TREE_SOLVERS else float("nan"),
                         "time_min_s": t["min"], "time_median_s": t["median"], "repeats": t["repeats"]})
            if path is not None:
                _write_csv(path, rows)
            print(f"  timing  {name:16s} N={n:>6}  min={t['min']:.4f}s  median={t['median']:.4f}s")
    return rows


def experiment_tree_split(modules, cfg: Config):
    """Where does Barnes-Hut spend its time? build (Python) vs. traversal."""
    rows = []
    for n in cfg.n_sweep:
        x, _, m = plummer_sphere(n, cfg.seed)
        first_new = len(rows)

        if "barnes_hut" in modules and n <= cfg.max_n.get("barnes_hut", np.inf):
            bh = modules["barnes_hut"]
            tb = []
            for _ in range(cfg.repeats):
                t0 = time.perf_counter()
                bh.build_tree(x, m)
                tb.append(time.perf_counter() - t0)
            total, _ = time_call(make_callable("barnes_hut", bh, cfg.eps, cfg.theta_default, cfg.G), x, m,
                                 repeats=max(1, min(2, cfg.repeats)))
            rows.append({"solver": "barnes_hut", "n": n, "build_s": min(tb),
                         "pack_s": 0.0, "traverse_s": max(total["min"] - min(tb), 0.0),
                         "total_s": total["min"]})

        if "fast_barnes_hut" in modules and n <= cfg.max_n.get("fast_barnes_hut", np.inf):
            fbh = modules["fast_barnes_hut"]
            tb, tp, tt = [], [], []
            for _ in range(cfg.repeats):
                t0 = time.perf_counter()
                root = fbh.build_tree(x, m)  # imported into fast module namespace
                t1 = time.perf_counter()
                packed = fbh.packed_tree(x, m)  # includes its own build_tree call
                t2 = time.perf_counter()
                fbh._traverse(x, m, *packed, cfg.eps, cfg.theta_default, cfg.G)
                t3 = time.perf_counter()
                tb.append(t1 - t0)
                tp.append(max((t2 - t1) - (t1 - t0), 0.0))  # packing = packed_tree minus its build
                tt.append(t3 - t2)
            rows.append({"solver": "fast_barnes_hut", "n": n, "build_s": min(tb),
                         "pack_s": min(tp), "traverse_s": min(tt),
                         "total_s": min(tb) + min(tp) + min(tt)})
        for r in rows[first_new:]:
            print(f"  tree split  {r['solver']:16s} N={n:>6}  build={r['build_s']:.3f}s  "
                  f"pack={r['pack_s']:.3f}s  traverse={r['traverse_s']:.3f}s")
    return rows


def experiment_conservation(modules, cfg: Config):
    x0, v0, m = plummer_sphere(cfg.n_cons, cfg.seed)
    rows = []
    for name, mod in modules.items():
        accel = make_callable(name, mod, cfg.eps, cfg.theta_default, cfg.G)
        print(f"  conservation  {name:16s} N={cfg.n_cons} steps={cfg.steps} dt={cfg.dt} ...", flush=True)
        series, timing = leapfrog_run(accel, x0, v0, m, cfg.eps, cfg.dt, cfg.steps, cfg.diag_every, cfg.G)
        for r in series:
            r.update({"solver": name, "wall_total_s": timing["wall_total"],
                      "wall_force_s": timing["wall_force"]})
        rows.extend(series)
        print(f"      max |dE/E| = {max(abs(r['dE_rel']) for r in series):.2e}   "
              f"final p drift = {series[-1]['p_drift']:.2e}   wall = {timing['wall_total']:.1f}s")
    return rows


# --------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------
def machine_info():
    info = {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "processor": platform.processor(),
        "cpu_count": os.cpu_count(),
        "numpy": np.__version__,
        "env_threads": {k: os.environ.get(k) for k in
                        ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS")},
    }
    try:
        import numba
        info["numba"] = numba.__version__
        info["numba_threads"] = numba.get_num_threads()
    except Exception:
        info["numba"] = None
    return info


def run_all(package: str, out: str = "results", quick: bool = False, cfg: Config | None = None,
            experiments=("sanity", "force_error", "timing", "tree_split", "conservation")):
    cfg = cfg or (Config.quick(package, out) if quick else Config(package=package, out=out))
    out_dir = Path(cfg.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    modules = load_solvers(cfg.package, cfg.solvers)
    print(f"Solvers loaded: {', '.join(modules)}")
    print("Warming up (JIT compile / cache load)...")
    warm_up(modules, cfg.eps, cfg.theta_default)

    meta = {"config": {k: (list(v) if isinstance(v, tuple) else v) for k, v in cfg.__dict__.items()},
            "solvers": list(modules), "machine": machine_info()}
    (out_dir / "meta.json").write_text(json.dumps(meta, indent=2))

    if "sanity" in experiments:
        print("[1/5] Sanity checks")
        sanity = experiment_sanity(modules, cfg)
        (out_dir / "sanity.json").write_text(json.dumps(sanity, indent=2))
        for k, v in sanity.items():
            print(f"  {k}: {v}")
    if "force_error" in experiments:
        print("[2/5] Force error vs. exact direct summation")
        _write_csv(out_dir / "force_error.csv", experiment_force_error(modules, cfg))
    if "timing" in experiments:
        print("[3/5] Timing sweep")
        _write_csv(out_dir / "timing.csv", experiment_timing(modules, cfg, path=out_dir / "timing.csv"))
    if "tree_split" in experiments:
        print("[4/5] Barnes-Hut build vs. traversal")
        _write_csv(out_dir / "tree_split.csv", experiment_tree_split(modules, cfg))
    if "conservation" in experiments:
        print("[5/5] Conservation over time (leapfrog)")
        _write_csv(out_dir / "conservation.csv", experiment_conservation(modules, cfg))

    print(f"Done. Results in {out_dir.resolve()}")
    return out_dir


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--package", required=True, help="dotted package containing direct.py, barnes_hut.py, ...")
    p.add_argument("--out", default="results")
    p.add_argument("--quick", action="store_true", help="small, fast settings for a smoke test")
    p.add_argument("--seed", type=int, default=None)
    p.add_argument("--eps", type=float, default=None)
    p.add_argument("--theta", type=float, default=None, help="default opening angle for timing/conservation")
    p.add_argument("--n-cons", type=int, default=None)
    p.add_argument("--steps", type=int, default=None)
    p.add_argument("--dt", type=float, default=None)
    p.add_argument("--n-sweep", nargs="+", type=int, help="particle counts for the timing/tree-split sweep")
    p.add_argument("--max-n", nargs="+", metavar="SOLVER=N",
                   help="per-solver N cap, e.g. direct=4000 barnes_hut=2000 fast_direct=60000")
    p.add_argument("--solvers", nargs="+", choices=list(SOLVER_MODULES),
                   help="run only these solvers (e.g. the two fast ones for large N)")
    p.add_argument("--repeats", type=int, default=None)
    p.add_argument("--only", nargs="+", choices=["sanity", "force_error", "timing", "tree_split", "conservation"])
    args = p.parse_args(argv)

    cfg = Config.quick(args.package, args.out) if args.quick else Config(package=args.package, out=args.out)
    for attr, val in (("seed", args.seed), ("eps", args.eps), ("theta_default", args.theta),
                      ("n_cons", args.n_cons), ("steps", args.steps), ("dt", args.dt)):
        if val is not None:
            setattr(cfg, attr, val)
    if args.n_sweep:
        cfg.n_sweep = tuple(args.n_sweep)
    if args.max_n:
        for item in args.max_n:
            k, v = item.split("=")
            cfg.max_n[k] = int(v)
    if args.solvers:
        cfg.solvers = tuple(args.solvers)
    if args.repeats:
        cfg.repeats = args.repeats
    kwargs = {"experiments": tuple(args.only)} if args.only else {}
    run_all(args.package, cfg=cfg, **kwargs)


if __name__ == "__main__":
    main()
