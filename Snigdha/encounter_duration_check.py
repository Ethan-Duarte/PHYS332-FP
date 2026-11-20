"""Check whether an encounter has evolved far enough after first pericenter.

Save this file beside revised_nbody_analysis.py, then run it from the same
Python environment/repository setup as the main analysis script. It performs
ONE encounter integration (not a parameter sweep), and records diagnostics at
several checkpoints along that same trajectory.

This is a duration diagnostic, not a validated production result. The preview
initial conditions and the disk-energy/unbinding diagnostic have limitations;
see the main analysis script and project notes.
"""

import copy
import csv
from functools import partial
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Reuse the exact galaxy setup and diagnostic definitions from the main script.
import snigdha_galaxy_analysis as analysis
from nbody.integrate import evolve
from nbody.analysis import core_separation
from nbody.gravity.fast_barnes_hut import acceleration as fast_tree

# ---------------------------------------------------------------------
# Test settings: edit these if you want a different representative case.
# ---------------------------------------------------------------------
Q = 1.0
INCLINATION_DEG = 0
DT = analysis.DT
THETA = analysis.THETA
EPS = analysis.EPS
N_PRIMARY = analysis.N_PRIMARY
SAVE_EVERY = 10
CHECKPOINT_STEPS = [200, 400, 800, 1000]
MAX_STEPS = max(CHECKPOINT_STEPS)

# Keep output separate from the full sweep so it cannot overwrite its files.
OUT = analysis.ROOT / "Ethan" / "results" / "encounter_duration_check"
OUT.mkdir(parents=True, exist_ok=True)


def main():
    print("Single-encounter duration check")
    print(f"q={Q}, inclination={INCLINATION_DEG} deg, N_primary={N_PRIMARY}")
    print(f"dt={DT}, max steps={MAX_STEPS}, max time={DT * MAX_STEPS:.3f}")
    print(f"Checkpoints (steps): {CHECKPOINT_STEPS}")

    initial, cfg, counts_a, counts_b = analysis.make_initial_pair(
        Q, INCLINATION_DEG, N_PRIMARY
    )
    rd = analysis.disk_scale_length(cfg)
    force = partial(fast_tree, eps=EPS, theta=THETA)

    # Capture the actual t=0 state and then save only requested checkpoints.
    p0 = copy.deepcopy(initial)
    checkpoint_states = {}
    separation_curve = []
    wanted = set(CHECKPOINT_STEPS)

    for t, p in evolve(copy.deepcopy(initial), force, DT, MAX_STEPS, SAVE_EVERY):
        step = int(round(float(t) / DT))
        separation_curve.append((float(t), float(core_separation(p))))
        if step in wanted:
            checkpoint_states[step] = copy.deepcopy(p)
            print(f"Captured step {step} at t={t:.4f}")

    missing = sorted(wanted - checkpoint_states.keys())
    if missing:
        raise RuntimeError(
            f"Integrator did not yield requested checkpoint steps {missing}. "
            "Check evolve() save/yield behavior and align checkpoints with SAVE_EVERY."
        )

    # Summarize how the endpoint diagnostics change as the same run is extended.
    rows = []
    for step in CHECKPOINT_STEPS:
        state = checkpoint_states[step]
        diagnostics = analysis.measure_encounter_physics(p0, state, rd, eps=EPS)
        delta_e_orb, outside_pct, mean_delta_eps, eps90, positive_pct = diagnostics
        nearest_sep = min(separation_curve, key=lambda pair: abs(pair[0] - step * DT))[1]
        row = {
            "steps": step,
            "time": step * DT,
            "core_separation": nearest_sep,
            "delta_E_orb_proxy": delta_e_orb,
            "final_fraction_3d_radius_gt_3p5Rd_percent": outside_pct,
            "mean_delta_disk_specific_energy": mean_delta_eps,
            "final_disk_eps90": eps90,
            "positive_host_specific_energy_percent": positive_pct,
        }
        rows.append(row)
        print(
            f"step={step:4d}, t={step * DT:6.3f}, sep={nearest_sep:.4f}, "
            f"outside(3.5Rd)={outside_pct:6.2f}%, "
            f"mean_delta_eps={mean_delta_eps:+.5f}, "
            f"positive_host_eps={positive_pct:.2f}%"
        )

    csv_path = OUT / "duration_check.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    # Figure 1: core separation over the entire run, with endpoint checkpoints.
    times = np.array([item[0] for item in separation_curve])
    seps = np.array([item[1] for item in separation_curve])
    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.plot(times, seps, color="tab:blue", lw=1.8, label="Core separation")
    for step in CHECKPOINT_STEPS:
        t = step * DT
        ax.axvline(t, color="0.5", ls="--", lw=0.8)
        ax.scatter([t], [np.interp(t, times, seps)], s=24, zorder=3)
    ax.set_xlabel("Time (code units)")
    ax.set_ylabel("Bulge-core separation (code units)")
    ax.set_title(f"Encounter duration check (q={Q:g}, i={INCLINATION_DEG}°)")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "duration_separation.png", dpi=160)
    plt.close(fig)

    # Figure 2: endpoint diagnostics versus elapsed time.
    checkpoint_times = np.array([row["time"] for row in rows])
    outside_vals = np.array([
        row["final_fraction_3d_radius_gt_3p5Rd_percent"] for row in rows
    ])
    energy_vals = np.array([row["mean_delta_disk_specific_energy"] for row in rows])
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    axes[0].plot(checkpoint_times, outside_vals, "o-")
    axes[0].set_xlabel("Endpoint time (code units)")
    axes[0].set_ylabel("Final fraction at 3D radius > 3.5 Rd (%)")
    axes[0].set_title("Spatial diagnostic vs. endpoint")
    axes[0].grid(alpha=0.25)
    axes[1].plot(checkpoint_times, energy_vals, "o-", color="tab:orange")
    axes[1].set_xlabel("Endpoint time (code units)")
    axes[1].set_ylabel("Mean disk specific-energy change")
    axes[1].set_title("Energy diagnostic vs. endpoint")
    axes[1].grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(OUT / "duration_diagnostics.png", dpi=160)
    plt.close(fig)

    print(f"\nSaved CSV and figures to: {OUT}")
    print("Interpretation: compare successive checkpoint values; no universal "
          "settling tolerance is assumed by this script.")


if __name__ == "__main__":
    main()
