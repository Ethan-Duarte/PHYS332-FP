# # """
# # Parameter-sweep and energy analysis script (no modifications to nbody/).
# # Evaluates:
# #   1. Numerical energy conservation across time and Barnes-Hut opening angle (theta)
# #      based on Barnes & Hut (1986).
# #   2. Two-variable physical encounter outcomes (inclination vs. mass ratio)
# #      based on Holmberg (1941) and Toomre & Toomre (1972).
# # """
# # import json
# # from pathlib import Path
# # from functools import partial
# # from time import perf_counter
# # import numpy as np
# # import matplotlib

# # # Use non-interactive 'Agg' backend so figures render directly to PNG files
# # # without opening GUI windows during batch loops.
# # matplotlib.use('Agg')
# # import matplotlib.pyplot as plt

# # # Import existing simulation engine, galaxy builders, and diagnostics directly
# # from nbody.galaxies.generate import galaxy
# # from nbody.encounters import parabolic_pair, place, rotation
# # from nbody.particles import combine
# # from nbody.integrate import evolve
# # from nbody.diagnostics import total_energy
# # from nbody.analysis import core_separation
# # from nbody.gravity.fast_direct import acceleration as fast_direct
# # from nbody.gravity.fast_barnes_hut import acceleration as fast_tree

# # #parallel processing for simulations
# # from concurrent.futures import ProcessPoolExecutor, as_completed
# # import multiprocessing as mp
# # import copy


# # # Resolve project root directory and create a dedicated output folder for plots
# # ROOT = Path(__file__).resolve().parent.parent
# # OUT = ROOT / 'Ethan' / 'results' / 'my_analysis'
# # OUT.mkdir(parents=True, exist_ok=True)


# # def make_galaxy(mass_scale=1.0, n_scale=1.0, seed=42, backend='preview'):
# #     """
# #     Construct a multi-component (disk + bulge + dark matter halo) galaxy
# #     from the shared Mice preset, scaling mass and particle resolution.
# #     """
# #     cfg = json.loads((ROOT / 'Ethan' / 'presets' / 'mice.json').read_text())
# #     cfg['seed'] = seed

# #     # Scale component masses and particle counts proportionally so particle mass
# #     # resolution stays consistent across different galaxy mass ratios.
# #     for comp in ('disk', 'bulge', 'halo'):
# #         cfg[comp]['mass'] *= mass_scale
# #         cfg['particles'][comp] = max(10, int(cfg['particles'][comp] * n_scale))

# #     return galaxy(cfg, backend)


# # def specific_binding_energy(p,galaxy_id=0,eps=0.0125,G=1.0,chunk_size=512):
# #     """
# #     Calculate approximate internal specific energy of disk particles
# #     relative to their own galaxy's bulk velocity.

# #     Includes:
# #       1. Relative kinetic energy.
# #       2. Gravitational potential from the same galaxy.
# #     """

# #     # Identify galaxy particles and bulge core
# #     galaxy_mask = p.galaxy == galaxy_id
# #     core_mask = galaxy_mask & (p.component == 1)
# #     disk_mask = galaxy_mask & (p.component == 0)

# #     # Galaxy bulk velocity, estimated using its bulge
# #     core_velocity = np.average(
# #         p.v[core_mask],
# #         axis=0,
# #         weights=p.m[core_mask]
# #     )

# #     # Relative kinetic energy
# #     relative_velocity = p.v[disk_mask] - core_velocity

# #     kinetic = 0.5 * np.sum(
# #         relative_velocity**2,
# #         axis=1
# #     )

# #     # Galaxy-only gravitational potential
# #     galaxy_positions = p.x[galaxy_mask]
# #     galaxy_masses = p.m[galaxy_mask]

# #     disk_positions = p.x[disk_mask]

# #     N_disk = len(disk_positions)
# #     potential = np.zeros(N_disk)

# #     galaxy_indices = np.flatnonzero(galaxy_mask)
# #     disk_indices = np.flatnonzero(disk_mask)

# #     # Chunked calculation to limit memory consumption
# #     for start in range(0, N_disk, chunk_size):

# #         end = min(start + chunk_size, N_disk)

# #         dx = (
# #             disk_positions[start:end, None, :]
# #             - galaxy_positions[None, :, :]
# #         )

# #         r2 = np.sum(dx**2, axis=2) + eps**2

# #         # Exclude self-interaction
# #         rows = np.arange(end - start)
# #         cols = np.searchsorted(
# #             galaxy_indices,
# #             disk_indices[start:end]
# #         )

# #         r2[rows, cols] = np.inf

# #         potential[start:end] = -G * np.sum(
# #             galaxy_masses / np.sqrt(r2),
# #             axis=1
# #         )

# #     return kinetic + potential

# # def run_encounter(initial, force, dt, steps, save_every,
# #                   track_separation=False):

# #     initial_state = None
# #     final_state = None
# #     separation_curve = []

# #     for t, p in evolve(initial, force, dt, steps, save_every):

# #         if initial_state is None:
# #             initial_state = copy.deepcopy(p)

# #         final_state = copy.deepcopy(p)

# #         if track_separation:
# #             separation_curve.append((t, core_separation(p)))

# #     return initial_state, final_state, separation_curve

# # # ==============================================================================
# # # PART 1: Energy Conservation Over Time & Across Opening Angle Theta (Barnes-Hut)
# # # ==============================================================================
# # def study_energy_vs_theta(eps=0.0125, dt=0.006, steps=200, save_every=10):
# #     """
# #     Tests numerical accuracy and computational cost of the Barnes-Hut tree code
# #     as a function of time and opening angle theta (Barnes & Hut 1986).
# #     """
# #     print("Running Study 1: Energy conservation vs. time & opening angle theta...")

# #     # Build two identical baseline galaxies on a parabolic encounter trajectory
# #     g1 = make_galaxy(mass_scale=1.0, n_scale=1, seed=42)
# #     g2 = make_galaxy(mass_scale=1.0, n_scale=1, seed=43)
# #     initial = parabolic_pair(
# #         g1, g2, separation=1.0, pericenter=0.25,
# #         orientations=((25, -30), (40, 60))
# #     )

# #     # Run a tiny dummy force calculation first to trigger Numba JIT compilation;
# #     # this prevents one-time compiler overhead from skewing our runtime benchmarks.
# #     fast_direct(initial.x[:5], initial.m[:5], eps=eps)
# #     fast_tree(initial.x[:5], initial.m[:5], eps=eps, theta=0.5)

# #     # Sweep theta from 0.0 (exact O(N^2) direct summation) up to 1.0 (aggressive tree grouping)
# #     thetas = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
# #     fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))

# #     max_drifts = []
# #     runtimes = []

# #     for theta in thetas:
# #         # Select either exact direct summation or Barnes-Hut tree with opening angle theta
# #         if theta == 0.0:
# #             force = partial(fast_direct, eps=eps)
# #             label = 'Direct (θ = 0)'
# #         else:
# #             force = partial(fast_tree, eps=eps, theta=theta)
# #             label = f'BH θ = {theta:.1f}'

# #         # Time the full leapfrog integration and collect saved snapshots
# #         t0 = perf_counter()
# #         snapshots = list(
# #             evolve(copy.deepcopy(initial), force, dt, steps, save_every)
# #         )   
# #         elapsed = perf_counter() - t0

# #         # Compute total energy E(t) = Kinetic + Softened Potential at every snapshot
# #         times = np.array([t for t, _ in snapshots])
# #         energies = np.array([total_energy(p, eps) for _, p in snapshots])

# #         # Fractional energy error |E(t) - E(0)| / |E(0)|
# #         rel_drift = np.abs(energies - energies[0]) / abs(energies[0])

# #         # Plot error evolution over time on left panel
# #         ax1.plot(times, rel_drift * 100, lw=1.8, label=label)
# #         max_drifts.append(np.max(rel_drift) * 100)
# #         runtimes.append(elapsed)
# #         print(f"  {label:15s} | max |ΔE/E0| = {np.max(rel_drift)*100:.4f}% | time = {elapsed:.2f}s")

# #     ax1.set(
# #         xlabel='Time (code units)',
# #         ylabel='Fractional Energy Drift |ΔE / E₀| (%)',
# #         title='Energy Conservation Over Time'
# #     )
# #     ax1.set_yscale('log')  # Log scale highlights orders-of-magnitude differences across theta
# #     ax1.grid(alpha=0.25)
# #     ax1.legend(fontsize=8)

# #     # Right panel: Dual-axis plot comparing max energy error (left axis) vs runtime (right axis)
# #     color1, color2 = 'tab:blue', 'tab:red'
# #     ax2.plot(thetas, max_drifts, 'o-', color=color1, lw=2, label='Max Energy Error (%)')
# #     ax2.set_xlabel('Barnes–Hut Opening Angle θ (0 = Direct)')
# #     ax2.set_ylabel('Max |ΔE / E₀| (%)', color=color1)
# #     ax2.tick_params(axis='y', labelcolor=color1)
# #     ax2.grid(alpha=0.25)

# #     ax2_r = ax2.twinx()
# #     ax2_r.plot(thetas, runtimes, 's--', color=color2, lw=2, label='Runtime (s)')
# #     ax2_r.set_ylabel('Wall-Clock Runtime (s)', color=color2)
# #     ax2_r.tick_params(axis='y', labelcolor=color2)
# #     ax2.set_title('Accuracy vs. Computational Cost')

# #     fig.tight_layout()
# #     out_path = OUT / 'energy_vs_theta.png'
# #     fig.savefig(out_path, dpi=150)
# #     plt.close(fig)
# #     print(f"Saved -> {out_path}\n")


# # # ==============================================================================
# # # PART 2: Two-Variable Encounter Sweep (Inclination vs. Mass Ratio)
# # # ==============================================================================

# # def measure_encounter_physics(p0, pf, eps=0.0125):
# #     """
# #     Quantifies encounter outcomes:
# #       1. ΔE_orb: Change in relative orbital energy between the two galaxy cores.
# #       2. stripped_fraction: Percentage of Galaxy A disk particles that get 
# #          displaced beyond 3.5x their initial scale length due to tidal disruption.
# #     """
# #     def orbital_energy(p):
# #         pos, vel, mass = [], [], []
# #         for gid in (0, 1):
# #             core = (p.galaxy == gid) & (p.component == 1)
# #             all_g = (p.galaxy == gid)
# #             pos.append(np.average(p.x[core], axis=0, weights=p.m[core]))
# #             vel.append(np.average(p.v[core], axis=0, weights=p.m[core]))
# #             mass.append(np.sum(p.m[all_g]))
# #         mu = (mass[0] * mass[1]) / (mass[0] + mass[1])
# #         v_rel = np.linalg.norm(vel[1] - vel[0])
# #         r_sep = np.sqrt(np.sum((pos[1] - pos[0])**2) + eps**2)
# #         return 0.5 * mu * v_rel**2 - (mass[0] * mass[1]) / r_sep

# #     # 1. Orbital energy change
# #     e_orb_init = orbital_energy(p0)
# #     e_orb_final = orbital_energy(pf)
# #     delta_e_orb = e_orb_final - e_orb_init

# #     # 2. Robust spatial tidal stripping measurement

# #     # Identify Galaxy A's disk particles
# #     disk_a_mask = (
# #         (p0.galaxy == 0) &
# #         (p0.component == 0)
# #     )

# #     # Identify Galaxy A's bulge core
# #     core_a0_mask = (
# #         (p0.galaxy == 0) &
# #         (p0.component == 1)
# #     )

# #     core_af_mask = (
# #         (pf.galaxy == 0) &
# #         (pf.component == 1)
# #     )

# #     # Initial and final Galaxy A core positions
# #     core_a0 = np.average(
# #         p0.x[core_a0_mask],
# #         axis=0,
# #         weights=p0.m[core_a0_mask]
# #     )

# #     core_af = np.average(
# #         pf.x[core_af_mask],
# #         axis=0,
# #         weights=pf.m[core_af_mask]
# #     )

# #     # Initial disk particle distances from galaxy center
# #     r0 = np.linalg.norm(
# #         p0.x[disk_a_mask] - core_a0,
# #         axis=1
# #     )

# #     # Final disk particle distances from galaxy center
# #     rf = np.linalg.norm(
# #         pf.x[disk_a_mask] - core_af,
# #         axis=1
# #     )

# #     # Radial displacement diagnostic
# #     # Use the same threshold consistently across encounters
# #     r_threshold = 3.5 * np.median(r0)

# #     stripped_fraction = float(
# #         np.mean(rf > r_threshold) * 100
# #     )

# #     # Internal specific energy at initial and final times
# #     eps_init = specific_binding_energy(
# #         p0, galaxy_id=0, eps=eps
# #     )

# #     eps_final = specific_binding_energy(
# #         pf, galaxy_id=0, eps=eps
# #     )

# #     # Change in mean internal disk energy
# #     mean_delta_eps = float(
# #         np.mean(eps_final - eps_init)
# #     )

# #     # 90th percentile of final disk specific energy
# #     eps_90_final = float(
# #         np.percentile(eps_final, 90)
# #     )

# #     # Approximate positive-energy fraction
# #     unbound_fraction = float(
# #         np.mean(eps_final > 0) * 100
# #     )

# #     print(
# #         f"  Radial displacement: {stripped_fraction:.2f}% | "
# #         f"Positive energy: {unbound_fraction:.2f}% | "
# #         f"Mean Δε: {mean_delta_eps:.4f}"
# #     )

# #     return (
# #         delta_e_orb,
# #         stripped_fraction,
# #         mean_delta_eps,
# #         eps_90_final,
# #         unbound_fraction
# #     )



# # def run_single_encounter(args):
# #     """
# #     Run one independent galaxy encounter in a separate process.
# #     Returns the physical diagnostics and optional separation curve.
# #     """
# #     q, inc, eps, theta, dt, steps, save_every = args

# #     # Each worker creates its own galaxies to avoid sharing mutable particle data
# #     g1 = make_galaxy(mass_scale=1.0, n_scale=1, seed=42)
# #     g2 = make_galaxy(
# #         mass_scale=q,
# #         n_scale=max(0.25, 0.5 * q),
# #         seed=43
# #     )

# #     # Place galaxies on a parabolic trajectory with Galaxy A tilted by inclination `inc`
# #     initial = parabolic_pair(
# #         g1, g2,
# #         separation=1.0,
# #         pericenter=0.25,
# #         orientations=((inc, 0), (0, 0))
# #     )

# #     # Select Barnes-Hut force calculation
# #     force = partial(fast_tree, eps=eps, theta=theta)

# #     # Run the encounter
# #     p0, pf, separation_curve = run_encounter(
# #         initial, force, dt, steps, save_every, 
# #         track_separation=(q == 1.0 and inc in (0, 90, 180))
# #     )

# #     # Extract quantitative encounter diagnostics
# #     dE_orb, stripped, mean_deps, eps_90, unbound = (measure_encounter_physics(p0, pf, eps=eps))

# #     # Only retain separation curves for equal-mass prograde, polar, and retrograde cases
# #     if q == 1.0 and inc in (0, 90, 180):
# #         saved_curve = separation_curve
# #     else:
# #         saved_curve = None

# #     return (q, inc, dE_orb, stripped,mean_deps,eps_90,unbound, saved_curve)


# # def study_encounter_grid(eps=0.0125,theta=0.6,dt=0.006,steps=200,save_every=10,max_workers=4):
# #     """
# #     Runs a 2D parameter sweep across orbital inclination (i) and companion mass ratio
# #     (q = M_B / M_A) to map orbital energy loss, tidal stripping, and orbital decay.
# #     """
# #     print("Running Study 2: 2D Sweep (Inclination vs. Mass Ratio)...")

# #     # Parameter grid: 0 deg = prograde/direct spin-orbit alignment; 180 deg = retrograde
# #     inclinations = list(range(0, 181, 10))
# #     mass_ratios = [0.25, 0.5, 1.0]

# #     # Allocate 2D arrays to store results for heatmaps
# #     dE_grid = np.zeros((len(mass_ratios), len(inclinations)))
# #     strip_grid = np.zeros((len(mass_ratios), len(inclinations)))
# #     sep_curves = {}

# #     deps_grid = np.zeros((len(mass_ratios), len(inclinations)))
# #     eps90_grid = np.zeros((len(mass_ratios), len(inclinations)))

# #     unbound_grid = np.zeros((len(mass_ratios), len(inclinations)))

# #     # Create all 57 independent encounter jobs
# #     jobs = [
# #         (q, inc, eps, theta, dt, steps, save_every)
# #         for q in mass_ratios
# #         for inc in inclinations
# #     ]

# #     # Map mass ratios and inclinations to array indices
# #     q_indices = {q: i for i, q in enumerate(mass_ratios)}
# #     inc_indices = {inc: i for i, inc in enumerate(inclinations)}

# #     print(f"Submitting {len(jobs)} encounters using {max_workers} workers...")

# #     # Run independent encounters in parallel
# #     with ProcessPoolExecutor(
# #         max_workers=max_workers,
# #         mp_context=mp.get_context("spawn")
# #     ) as executor:

# #         futures = [
# #             executor.submit(run_single_encounter, job)
# #             for job in jobs
# #         ]

# #         for completed, future in enumerate(as_completed(futures), start=1):

# #             # Retrieve results from completed worker
# #             (
# #                 q,
# #                 inc,
# #                 dE_orb,
# #                 stripped,
# #                 mean_deps,
# #                 eps_90,
# #                 unbound,
# #                 separation_curve
# #             ) = future.result()

# #             r_idx = q_indices[q]
# #             c_idx = inc_indices[inc]

# #             # Store numerical results
# #             dE_grid[r_idx, c_idx] = dE_orb
# #             strip_grid[r_idx, c_idx] = stripped
# #             deps_grid[r_idx, c_idx] = mean_deps
# #             eps90_grid[r_idx, c_idx] = eps_90
# #             unbound_grid[r_idx, c_idx] = unbound

# #             # Save separation vs. time curves for equal-mass prograde, polar, and retrograde cases
# #             if separation_curve is not None:
# #                 times = [t for t, _ in separation_curve]
# #                 seps = [sep for _, sep in separation_curve]
# #                 sep_curves[f'i = {inc}°'] = (times, seps)

# #             print(
# #                 f"[{completed}/{len(jobs)}] "
# #                 f"q={q:.2f}, i={inc:3d}° | "
# #                 f"ΔE_orb = {dE_orb:+.4f} | "
# #                 f"Stripped disk = {stripped:.1f}%"
# #             )

# #     print("All 57 encounters completed.")

# #     # Create a 3-panel summary figure
# #     fig, axs = plt.subplots(3, 2, figsize=(12, 8))
# #     axs = axs.ravel()

# #     def _plot_by_q(ax, grid, title, ylabel):
# #         for r_idx, q in enumerate(mass_ratios):
# #             ax.plot(inclinations, grid[r_idx], marker='o', lw=2, label=f'q = {q:.2f}')
# #         ax.set_xticks(inclinations, [f'{i}°' for i in inclinations])
# #         ax.set(xlabel='Inclination i', ylabel=ylabel, title=title)
# #         ax.grid(alpha=0.25)
# #         ax.legend()

# #     # Plot 1: Orbital energy change vs. inclination for each mass ratio q
# #     _plot_by_q(axs[0], dE_grid, 'Orbital Energy Change ΔE_orb', 'ΔE_orb (code units)')

# #     # Plot 2: Stripped disk fraction vs. inclination for each mass ratio q
# #     _plot_by_q(axs[1], strip_grid, 'Stripped Disk Fraction (%)', '% Disk Particles Stripped')

# #     # Plot 3: Bulge-to-bulge separation over time (equal mass, q = 1.0)
# #     for label, (times, seps) in sep_curves.items():
# #         axs[2].plot(times, seps, lw=2, label=label)
# #     axs[2].set(
# #         xlabel='Time (code units)',
# #         ylabel='Bulge Separation',
# #         title='Orbital Decay (Equal Mass, q = 1.0)'
# #     )
# #     axs[2].grid(alpha=0.25)
# #     axs[2].legend()

# #     # Plot 4: Mean disk specific energy gain (or swap in eps90_grid for ε_90)
# #     _plot_by_q(axs[3], deps_grid, 'Mean Disk Specific Energy Gain Δε', 'Mean Δε (code units)')

# #     _plot_by_q(axs[4],unbound_grid,'Unbound Disk Fraction','% Unbound Disk Particles')
# #     axs[5].axis('off')

# #     fig.tight_layout()
# #     out_path = OUT / 'encounter_2d_sweep.png'
# #     fig.savefig(out_path, dpi=150)
# #     plt.close(fig)
# #     print(f"Saved -> {out_path}\n")


# # if __name__ == '__main__':
# #     study_energy_vs_theta()
# #     study_encounter_grid()

# """
# Revised analysis script for the PHYS332-FP project.

# Key changes from the submitted version:
# - Keeps Barnes-Hut as the force solver for all production runs.
# - Uses 200 integration steps (not 200 particles).
# - Uses a configurable particle resolution, with 2048 particles in the primary
#   galaxy as an initial PERFORMANCE-TEST setting, not a proven crossover.
# - Scales companion particle counts approximately with q to keep particle
#   masses similar, subject to integer rounding.
# - Uses the Mice preset's exponential disk scale length for the 3.5 Rd
#   radial-displacement diagnostic.
# - Saves 200-step encounter results and makes the post-encounter interpretation
#   explicit: 200 steps may not cover the whole encounter; inspect separation.
# - Keeps the energy-vs-theta experiment separate from the inclination/mass-ratio
#   sweep.

# IMPORTANT:
# There is no universal "optimal N" for Barnes-Hut. The crossover depends on the
# actual tree builder, hardware, theta, particle distribution and direct solver.
# The repository's tree construction includes Python work, so do not assume BH
# beats direct summation at a particular N. This script deliberately does NOT
# switch to direct summation. Use the included timing output to characterize the
# Barnes-Hut runtime as N increases. If your project needs to claim BH is faster
# than direct, a separately controlled head-to-head benchmark is required.

# The native preview galaxy generator uses approximate velocities; it is not a
# validated equilibrium model. Treat encounter outcomes as exploratory until
# isolated-galaxy stability has been checked.
# """

# import copy
# import json
# from functools import partial
# from pathlib import Path
# from time import perf_counter
# from concurrent.futures import ProcessPoolExecutor, as_completed
# import multiprocessing as mp

# import numpy as np
# import matplotlib
# matplotlib.use("Agg")
# import matplotlib.pyplot as plt

# from nbody.galaxies.generate import galaxy
# from nbody.encounters import parabolic_pair
# from nbody.integrate import evolve
# from nbody.diagnostics import total_energy
# from nbody.analysis import core_separation
# from nbody.gravity.fast_barnes_hut import acceleration as fast_tree
# import csv


# # ---------------------------------------------------------------------
# # Configuration
# # ---------------------------------------------------------------------

# ROOT = Path(__file__).resolve().parent.parent
# PRESET_PATH = ROOT / "Ethan" / "presets" / "mice.json"
# OUT = ROOT / "Ethan" / "results" / "my_analysis_revised"
# OUT.mkdir(parents=True, exist_ok=True)

# EPS = 0.0125
# G = 1.0
# THETA = 0.6

# DT = 0.006
# STEPS = 200                  # User-requested integration length.
# SAVE_EVERY = 10              # Save every 10 steps: 21 snapshots incl. t=0.
# MAX_WORKERS = 4

# # Initial performance-test resolution. This is NOT a guaranteed BH crossover.
# # It is particles in the primary galaxy; companion count scales approximately q.
# N_PRIMARY = 2048

# # Coarser grid for a quick 200-step development run.
# INCLINATIONS = [0, 30, 60, 90, 120, 150, 180]
# MASS_RATIOS = [0.25, 0.5, 1.0]

# # Exponential disk scale length from the supplied mice.json preset.
# # The preset's disk["scale"] is read at runtime below; this is only fallback.
# FALLBACK_RD = 1.0 / 12.0


# # ---------------------------------------------------------------------
# # Galaxy generation and particle-count control
# # ---------------------------------------------------------------------

# def load_preset():
#     return json.loads(PRESET_PATH.read_text())


# def component_particle_counts(total, preset_particles):
#     """Allocate total particles in the preset's disk:bulge:halo proportions."""
#     names = ("disk", "bulge", "halo")
#     weights = np.array([preset_particles[k] for k in names], dtype=float)
#     raw = total * weights / weights.sum()
#     counts = np.floor(raw).astype(int)

#     # Guarantee at least one particle per component.
#     counts = np.maximum(counts, 1)
#     while counts.sum() > total:
#         candidates = np.where(counts > 1)[0]
#         if len(candidates) == 0:
#             raise ValueError("Total N too small for three components.")
#         j = candidates[np.argmax(counts[candidates] - raw[candidates])]
#         counts[j] -= 1
#     while counts.sum() < total:
#         j = int(np.argmax(raw - counts))
#         counts[j] += 1
#         raw[j] = counts[j]  # prevent repeatedly selecting same component

#     return dict(zip(names, counts.tolist()))


# def make_galaxy(mass_scale=1.0, seed=42, n_particles=None, backend="preview"):
#     cfg = load_preset()
#     cfg["seed"] = int(seed)

#     for comp in ("disk", "bulge", "halo"):
#         cfg[comp]["mass"] *= mass_scale

#     if n_particles is not None:
#         cfg["particles"].update(n_particles)

#     p = galaxy(cfg, backend)
#     return p, cfg


# def make_initial_pair(q, inclination, n_primary=N_PRIMARY):
#     """
#     Companion count approximately q*N_primary. Matching component allocations
#     approximately preserves per-particle mass between corresponding components.
#     """
#     preset = load_preset()
#     counts_a = component_particle_counts(n_primary, preset["particles"])
#     n_companion = max(3, int(round(q * n_primary)))
#     counts_b = component_particle_counts(n_companion, preset["particles"])

#     g1, cfg = make_galaxy(
#         mass_scale=1.0, seed=42, n_particles=counts_a
#     )
#     g2, _ = make_galaxy(
#         mass_scale=q, seed=43, n_particles=counts_b
#     )

#     initial = parabolic_pair(
#         g1, g2,
#         separation=1.0,
#         pericenter=0.25,
#         orientations=((inclination, 0), (0, 0))
#     )
#     return initial, cfg, counts_a, counts_b


# def disk_scale_length(cfg):
#     return float(cfg["disk"].get("scale", FALLBACK_RD))


# # ---------------------------------------------------------------------
# # Diagnostics
# # ---------------------------------------------------------------------

# def specific_binding_energy(p, galaxy_id=0, eps=EPS, G=G, chunk_size=512):
#     """Approximate disk specific energy in the galaxy-only potential."""
#     galaxy_mask = p.galaxy == galaxy_id
#     core_mask = galaxy_mask & (p.component == 1)
#     disk_mask = galaxy_mask & (p.component == 0)

#     if not np.any(core_mask) or not np.any(disk_mask):
#         raise ValueError("Galaxy needs bulge-core and disk particles.")

#     core_velocity = np.average(
#         p.v[core_mask], axis=0, weights=p.m[core_mask]
#     )
#     relative_velocity = p.v[disk_mask] - core_velocity
#     kinetic = 0.5 * np.sum(relative_velocity**2, axis=1)

#     galaxy_positions = p.x[galaxy_mask]
#     galaxy_masses = p.m[galaxy_mask]
#     disk_positions = p.x[disk_mask]
#     galaxy_indices = np.flatnonzero(galaxy_mask)
#     disk_indices = np.flatnonzero(disk_mask)

#     potential = np.zeros(len(disk_positions))
#     for start in range(0, len(disk_positions), chunk_size):
#         end = min(start + chunk_size, len(disk_positions))
#         dx = disk_positions[start:end, None, :] - galaxy_positions[None, :, :]
#         r2 = np.sum(dx * dx, axis=2) + eps**2
#         rows = np.arange(end - start)
#         cols = np.searchsorted(galaxy_indices, disk_indices[start:end])
#         r2[rows, cols] = np.inf
#         potential[start:end] = -G * np.sum(
#             galaxy_masses / np.sqrt(r2), axis=1
#         )

#     return kinetic + potential


# def run_encounter(initial, force, dt, steps, save_every):
#     """Run fixed-length integration and retain first/final states + separation."""
#     p0 = None
#     pf = None
#     separation_curve = []

#     for t, p in evolve(initial, force, dt, steps, save_every):
#         if p0 is None:
#             p0 = copy.deepcopy(p)
#         pf = copy.deepcopy(p)
#         separation_curve.append((float(t), float(core_separation(p))))

#     return p0, pf, separation_curve


# def measure_encounter_physics(p0, pf, rd, eps=EPS):
#     """Return orbital-energy change and clearly labelled disk diagnostics."""

#     def orbital_energy(p):
#         positions, velocities, masses = [], [], []
#         for gid in (0, 1):
#             core = (p.galaxy == gid) & (p.component == 1)
#             all_g = p.galaxy == gid
#             positions.append(np.average(p.x[core], axis=0, weights=p.m[core]))
#             velocities.append(np.average(p.v[core], axis=0, weights=p.m[core]))
#             masses.append(np.sum(p.m[all_g]))
#         mu = masses[0] * masses[1] / sum(masses)
#         vrel = np.linalg.norm(velocities[1] - velocities[0])
#         rsep = np.sqrt(np.sum((positions[1] - positions[0])**2) + eps**2)
#         return 0.5 * mu * vrel**2 - G * masses[0] * masses[1] / rsep

#     delta_e_orb = orbital_energy(pf) - orbital_energy(p0)

#     disk0 = (p0.galaxy == 0) & (p0.component == 0)
#     core0 = (p0.galaxy == 0) & (p0.component == 1)
#     coref = (pf.galaxy == 0) & (pf.component == 1)

#     center0 = np.average(p0.x[core0], axis=0, weights=p0.m[core0])
#     centerf = np.average(pf.x[coref], axis=0, weights=pf.m[coref])

#     # This is 3D displacement from the core, not cylindrical disk-plane radius.
#     r_final_3d = np.linalg.norm(pf.x[disk0] - centerf, axis=1)
#     outside_fraction = float(np.mean(r_final_3d > 3.5 * rd) * 100.0)

#     e0 = specific_binding_energy(p0, galaxy_id=0, eps=eps)
#     ef = specific_binding_energy(pf, galaxy_id=0, eps=eps)
#     mean_delta_eps = float(np.mean(ef - e0))
#     eps90 = float(np.percentile(ef, 90))
#     positive_energy_fraction = float(np.mean(ef > 0) * 100.0)

#     return (
#         float(delta_e_orb), outside_fraction, mean_delta_eps,
#         eps90, positive_energy_fraction
#     )


# # ---------------------------------------------------------------------
# # One worker per encounter
# # ---------------------------------------------------------------------

# def run_single_encounter(args):
#     q, inc, eps, theta, dt, steps, save_every, n_primary = args
#     initial, cfg, counts_a, counts_b = make_initial_pair(q, inc, n_primary)
#     disk_a_count = np.count_nonzero(
#     (initial.galaxy == 0) & (initial.component == 0))
#     print(
#         f"DEBUG q={q}, i={inc}: "
#         f"primary total={counts_a['disk'] + counts_a['bulge'] + counts_a['halo']}, "
#         f"primary disk={counts_a['disk']}, "
#         f"actual disk mask={disk_a_count}, "
#         f"pair total={len(initial.m)}"
#     )

#     rd = disk_scale_length(cfg)

#     force = partial(fast_tree, eps=eps, theta=theta)
#     p0, pf, separation_curve = run_encounter(
#         initial, force, dt, steps, save_every
#     )

#     diagnostics = measure_encounter_physics(p0, pf, rd, eps=eps)
#     return q, inc, *diagnostics, separation_curve, len(initial.m), rd


# # ---------------------------------------------------------------------
# # Study 1: energy conservation vs theta
# # ---------------------------------------------------------------------

# def study_energy_vs_theta(eps=EPS, dt=DT, steps=STEPS, save_every=SAVE_EVERY,
#                           n_primary=N_PRIMARY):
#     print("Study 1: Barnes-Hut energy conservation vs opening angle")
#     initial, _, _, _ = make_initial_pair(1.0, 0, n_primary)

#     # Warm up JIT compilation before timing; the force remains Barnes-Hut.
#     fast_tree(initial.x[:8], initial.m[:8], eps=eps, theta=0.5)

#     thetas = [0.2, 0.4, 0.6, 0.8, 1.0]
#     fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
#     max_drifts, runtimes = [], []

#     for theta in thetas:
#         force = partial(fast_tree, eps=eps, theta=theta)
#         t0 = perf_counter()
#         snapshots = list(evolve(
#             copy.deepcopy(initial), force, dt, steps, save_every
#         ))
#         elapsed = perf_counter() - t0

#         times = np.array([t for t, _ in snapshots])
#         energies = np.array([total_energy(p, eps) for _, p in snapshots])
#         drift = np.abs(energies - energies[0]) / max(abs(energies[0]), 1e-30)

#         ax1.plot(times, np.maximum(drift * 100, 1e-12),
#                  lw=1.7, label=f"BH theta={theta:.1f}")
#         max_drifts.append(float(np.max(drift) * 100))
#         runtimes.append(elapsed)
#         print(
#             f"theta={theta:.1f}: max |dE/E0|={max_drifts[-1]:.5g}% "
#             f"| runtime={elapsed:.3f}s"
#         )

#     ax1.set_yscale("log")
#     ax1.set_xlabel("Time (code units)")
#     ax1.set_ylabel("Absolute fractional energy drift (%)")
#     ax1.set_title("Energy conservation over time")
#     ax1.grid(alpha=0.25)
#     ax1.legend(fontsize=8)

#     ax2.plot(thetas, max_drifts, "o-", label="Maximum energy drift (%)")
#     ax2.set_xlabel("Barnes-Hut opening angle theta")
#     ax2.set_ylabel("Maximum energy drift (%)")
#     ax2.grid(alpha=0.25)
#     ax2r = ax2.twinx()
#     ax2r.plot(thetas, runtimes, "s--", color="tab:red", label="Runtime")
#     ax2r.set_ylabel("Runtime (s)")
#     ax2.set_title("Accuracy/runtime trade-off")

#     fig.tight_layout()
#     fig.savefig(OUT / "energy_vs_theta.png", dpi=150)
#     plt.close(fig)


# # ---------------------------------------------------------------------
# # Study 2: inclination vs mass ratio
# # ---------------------------------------------------------------------

# def study_encounter_grid(eps=EPS, theta=THETA, dt=DT, steps=STEPS,
#                          save_every=SAVE_EVERY, max_workers=MAX_WORKERS,
#                          n_primary=N_PRIMARY):
#     print(
#         f"Study 2: inclination/mass-ratio sweep | "
#         f"{steps} steps | primary N={n_primary}"
#     )

#     dE = np.zeros((len(MASS_RATIOS), len(INCLINATIONS)))
#     outside = np.zeros_like(dE)
#     deps = np.zeros_like(dE)
#     eps90 = np.zeros_like(dE)
#     positive = np.zeros_like(dE)
#     sep_curves = {}

#     jobs = [
#         (q, inc, eps, theta, dt, steps, save_every, n_primary)
#         for q in MASS_RATIOS for inc in INCLINATIONS
#     ]
#     qi = {q: i for i, q in enumerate(MASS_RATIOS)}
#     ii = {inc: i for i, inc in enumerate(INCLINATIONS)}

#     with ProcessPoolExecutor(
#         max_workers=max_workers,
#         mp_context=mp.get_context("spawn")
#     ) as executor:
#         futures = [executor.submit(run_single_encounter, job) for job in jobs]

#         for k, future in enumerate(as_completed(futures), start=1):
#             (q, inc, de, outfrac, md, e90, posfrac,
#              curve, n_total, rd) = future.result()

#             r, c = qi[q], ii[inc]
#             dE[r, c] = de
#             outside[r, c] = outfrac
#             deps[r, c] = md
#             eps90[r, c] = e90
#             positive[r, c] = posfrac

#             if q == 1.0 and inc in (0, 90, 180):
#                 sep_curves[f"i={inc} deg"] = curve

#             print(
#                 f"[{k}/{len(jobs)}] q={q:.2f}, i={inc:3d} deg, "
#                 f"N={n_total}, Rd={rd:.5f}, "
#                 f"deltaEorb={de:+.4g}, outside(3.5Rd)={outfrac:.1f}%"
#             )

#     # Save raw arrays in CSV for reproducibility.
#     csv_path = OUT / "encounter_grid.csv"
#     with csv_path.open("w", newline="") as f:
#         writer = csv.writer(f)
#         writer.writerow([
#             "q", "inclination_deg", "delta_E_orb",
#             "outside_3p5Rd_percent", "mean_delta_disk_specific_energy",
#             "final_disk_eps90", "positive_disk_energy_percent"
#         ])
#         for r, q in enumerate(MASS_RATIOS):
#             for c, inc in enumerate(INCLINATIONS):
#                 writer.writerow([
#                     q, inc, dE[r, c], outside[r, c], deps[r, c],
#                     eps90[r, c], positive[r, c]
#                 ])

#     fig, axs = plt.subplots(3, 2, figsize=(12, 9))
#     axs = axs.ravel()

#     def plot_by_q(ax, grid, title, ylabel):
#         for r, q in enumerate(MASS_RATIOS):
#             ax.plot(INCLINATIONS, grid[r], marker="o", lw=1.7,
#                     label=f"q={q:.2f}")
#         ax.set_xlabel("Inclination (degrees)")
#         ax.set_ylabel(ylabel)
#         ax.set_title(title)
#         ax.grid(alpha=0.25)
#         ax.legend(fontsize=8)

#     plot_by_q(axs[0], dE, "Change in approximate orbital energy",
#               "Delta E_orb (code units)")
#     plot_by_q(axs[1], outside, "Disk particles outside 3.5 Rd",
#               "Particles outside threshold (%)")

#     for label, curve in sep_curves.items():
#         times = [x[0] for x in curve]
#         seps = [x[1] for x in curve]
#         axs[2].plot(times, seps, label=label)
#     axs[2].set_xlabel("Time (code units)")
#     axs[2].set_ylabel("Bulge-core separation")
#     axs[2].set_title("Equal-mass core separation")
#     axs[2].grid(alpha=0.25)
#     axs[2].legend(fontsize=8)

#     plot_by_q(axs[3], deps, "Mean disk specific-energy change",
#               "Mean Delta epsilon")
#     plot_by_q(axs[4], positive, "Positive disk-specific-energy fraction",
#               "Positive energy (%)")
#     axs[5].axis("off")

#     fig.tight_layout()
#     fig.savefig(OUT / "encounter_2d_sweep.png", dpi=150)
#     plt.close(fig)

#     print(f"Saved plots and data in {OUT}")


# # ---------------------------------------------------------------------
# # Simple Barnes-Hut scaling benchmark (BH only)
# # ---------------------------------------------------------------------

# def benchmark_barnes_hut_scaling():
#     """
#     Times the project's Barnes-Hut force call for increasing N.
#     This measures this implementation's scaling, NOT the crossover against
#     direct summation. Run locally on the target machine.
#     """
#     rng = np.random.default_rng(1234)
#     sizes = [256, 512, 1024, 2048, 4096, 8192]
#     rows = []

#     # Compile/warm up before collecting timings.
#     x0 = rng.normal(size=(32, 3))
#     m0 = np.full(32, 1.0 / 32)
#     fast_tree(x0, m0, eps=EPS, theta=THETA)

#     for n in sizes:
#         x = rng.normal(size=(n, 3))
#         m = np.full(n, 1.0 / n)

#         # One warmup call, then three measured calls.
#         fast_tree(x, m, eps=EPS, theta=THETA)
#         samples = []
#         for _ in range(3):
#             t0 = perf_counter()
#             fast_tree(x, m, eps=EPS, theta=THETA)
#             samples.append(perf_counter() - t0)

#         median = float(np.median(samples))
#         rows.append((n, median))
#         print(f"Barnes-Hut N={n:5d}: median force-call time={median:.5f}s")

#     with (OUT / "barnes_hut_scaling.csv").open("w", newline="") as f:
#         writer = csv.writer(f)
#         writer.writerow(["N", "median_force_call_seconds"])
#         writer.writerows(rows)

#     ns, ts = zip(*rows)
#     plt.figure(figsize=(6.5, 4.2))
#     plt.plot(ns, ts, "o-")
#     plt.xlabel("Particle count N")
#     plt.ylabel("Median Barnes-Hut force-call time (s)")
#     plt.title("Barnes-Hut runtime scaling on this machine")
#     plt.grid(alpha=0.25)
#     plt.tight_layout()
#     plt.savefig(OUT / "barnes_hut_scaling.png", dpi=150)
#     plt.close()


# if __name__ == "__main__":
#     print(f"Preset: {PRESET_PATH}")
#     print(f"dt={DT}, steps={STEPS}, simulated duration={DT * STEPS:.3f}")
#     print(f"Primary particle target N={N_PRIMARY}; q={MASS_RATIOS}")
#     print(f"Barnes-Hut theta={THETA}; workers={MAX_WORKERS}")
#     benchmark_barnes_hut_scaling()
#     study_energy_vs_theta()
#     study_encounter_grid()

"""
Revised analysis script for the PHYS332-FP project.

Key changes from the submitted version:
- Keeps Barnes-Hut as the force solver for all production runs.
- Uses 200 integration steps (not 200 particles).
- Uses a configurable particle resolution, with 2048 particles in the primary
  galaxy as an initial PERFORMANCE-TEST setting, not a proven crossover.
- Scales companion particle counts approximately with q to keep particle
  masses similar, subject to integer rounding.
- Uses the Mice preset's exponential disk scale length for the 3.5 Rd
  radial-displacement diagnostic.
- Saves 200-step encounter results and makes the post-encounter interpretation
  explicit: 200 steps may not cover the whole encounter; inspect separation.
- Keeps the energy-vs-theta experiment separate from the inclination/mass-ratio
  sweep.

IMPORTANT:
There is no universal "optimal N" for Barnes-Hut. The crossover depends on the
actual tree builder, hardware, theta, particle distribution and direct solver.
The repository's tree construction includes Python work, so do not assume BH
beats direct summation at a particular N. This script deliberately does NOT
switch to direct summation. Use the included timing output to characterize the
Barnes-Hut runtime as N increases. If your project needs to claim BH is faster
than direct, a separately controlled head-to-head benchmark is required.

The native preview galaxy generator uses approximate velocities; it is not a
validated equilibrium model. Treat encounter outcomes as exploratory until
isolated-galaxy stability has been checked.
"""

import copy
import csv
import json
from functools import partial
from pathlib import Path
from time import perf_counter
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from nbody.galaxies.generate import galaxy
from nbody.encounters import parabolic_pair
from nbody.integrate import evolve
from nbody.diagnostics import total_energy
from nbody.analysis import core_separation
from nbody.gravity.fast_barnes_hut import acceleration as fast_tree


# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------

ROOT = Path(__file__).resolve().parent.parent
PRESET_PATH = ROOT / "Ethan" / "presets" / "mice.json"
OUT = ROOT / "Ethan" / "results" / "my_analysis_revised"
OUT.mkdir(parents=True, exist_ok=True)

EPS = 0.0125
G = 1.0
THETA = 0.6

DT = 0.006
STEPS = 200                  # User-requested integration length.
SAVE_EVERY = 10              # Save every 10 steps: 21 snapshots incl. t=0.
MAX_WORKERS = 4

# Initial performance-test resolution. This is NOT a guaranteed BH crossover.
# It is particles in the primary galaxy; companion count scales approximately q.
N_PRIMARY = 2048

# Coarser grid for a quick 200-step development run.
INCLINATIONS = [0, 30, 60, 90, 120, 150, 180]
MASS_RATIOS = [0.25, 0.5, 1.0]

# Exponential disk scale length from the supplied mice.json preset.
# The preset's disk["scale"] is read at runtime below; this is only fallback.
FALLBACK_RD = 1.0 / 12.0


# ---------------------------------------------------------------------
# Galaxy generation and particle-count control
# ---------------------------------------------------------------------

def load_preset():
    return json.loads(PRESET_PATH.read_text())


def component_particle_counts(total, preset_particles):
    """Allocate total particles in the preset's disk:bulge:halo proportions."""
    names = ("disk", "bulge", "halo")
    weights = np.array([preset_particles[k] for k in names], dtype=float)
    raw = total * weights / weights.sum()
    counts = np.floor(raw).astype(int)

    # Guarantee at least one particle per component.
    counts = np.maximum(counts, 1)
    while counts.sum() > total:
        candidates = np.where(counts > 1)[0]
        if len(candidates) == 0:
            raise ValueError("Total N too small for three components.")
        j = candidates[np.argmax(counts[candidates] - raw[candidates])]
        counts[j] -= 1
    while counts.sum() < total:
        j = int(np.argmax(raw - counts))
        counts[j] += 1
        raw[j] = counts[j]  # prevent repeatedly selecting same component

    return dict(zip(names, counts.tolist()))


def make_galaxy(mass_scale=1.0, seed=42, n_particles=None, backend="preview"):
    cfg = load_preset()
    cfg["seed"] = int(seed)

    for comp in ("disk", "bulge", "halo"):
        cfg[comp]["mass"] *= mass_scale

    if n_particles is not None:
        cfg["particles"].update(n_particles)

    p = galaxy(cfg, backend)
    return p, cfg


def make_initial_pair(q, inclination, n_primary=N_PRIMARY):
    """
    Companion count approximately q*N_primary. Matching component allocations
    approximately preserves per-particle mass between corresponding components.
    """
    preset = load_preset()
    counts_a = component_particle_counts(n_primary, preset["particles"])
    n_companion = max(3, int(round(q * n_primary)))
    counts_b = component_particle_counts(n_companion, preset["particles"])

    g1, cfg = make_galaxy(
        mass_scale=1.0, seed=42, n_particles=counts_a
    )
    g2, _ = make_galaxy(
        mass_scale=q, seed=43, n_particles=counts_b
    )

    initial = parabolic_pair(
        g1, g2,
        separation=1.0,
        pericenter=0.25,
        orientations=((inclination, 0), (0, 0))
    )
    return initial, cfg, counts_a, counts_b


def disk_scale_length(cfg):
    return float(cfg["disk"].get("scale", FALLBACK_RD))


# ---------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------

def specific_binding_energy(p, galaxy_id=0, eps=EPS, G=G, chunk_size=512):
    """Approximate disk specific energy in the galaxy-only potential."""
    galaxy_mask = p.galaxy == galaxy_id
    core_mask = galaxy_mask & (p.component == 1)
    disk_mask = galaxy_mask & (p.component == 0)

    if not np.any(core_mask) or not np.any(disk_mask):
        raise ValueError("Galaxy needs bulge-core and disk particles.")

    core_velocity = np.average(
        p.v[core_mask], axis=0, weights=p.m[core_mask]
    )
    relative_velocity = p.v[disk_mask] - core_velocity
    kinetic = 0.5 * np.sum(relative_velocity**2, axis=1)

    galaxy_positions = p.x[galaxy_mask]
    galaxy_masses = p.m[galaxy_mask]
    disk_positions = p.x[disk_mask]
    galaxy_indices = np.flatnonzero(galaxy_mask)
    disk_indices = np.flatnonzero(disk_mask)

    potential = np.zeros(len(disk_positions))
    for start in range(0, len(disk_positions), chunk_size):
        end = min(start + chunk_size, len(disk_positions))
        dx = disk_positions[start:end, None, :] - galaxy_positions[None, :, :]
        r2 = np.sum(dx * dx, axis=2) + eps**2
        rows = np.arange(end - start)
        cols = np.searchsorted(galaxy_indices, disk_indices[start:end])
        r2[rows, cols] = np.inf
        potential[start:end] = -G * np.sum(
            galaxy_masses / np.sqrt(r2), axis=1
        )

    return kinetic + potential


def run_encounter(initial, force, dt, steps, save_every):
    """Run fixed-length integration and retain first/final states + separation."""
    p0 = None
    pf = None
    separation_curve = []

    for t, p in evolve(initial, force, dt, steps, save_every):
        if p0 is None:
            p0 = copy.deepcopy(p)
        pf = copy.deepcopy(p)
        separation_curve.append((float(t), float(core_separation(p))))

    return p0, pf, separation_curve

def measure_displaced_disk_particles(p0, pf, multiplier=1.5, radius_percentile=95.0):
    """
    Identifies Galaxy A disk particles displaced beyond `multiplier` times
    the initial disk's characteristic/outer radius relative to Galaxy A's core.

    Returns:
        displaced_percent (float): Percentage of Galaxy A disk particles displaced.
        displaced_mask (np.ndarray): Boolean mask (length = N_disk_a) of displaced disk stars.
        r_threshold (float): Radial cutoff used for displacement.
        r0 (np.ndarray): Initial radii of Galaxy A disk particles from initial core.
        rf (np.ndarray): Final radii of those same disk particles from final core.
    """
    disk0 = (p0.galaxy == 0) & (p0.component == 0)
    core0 = (p0.galaxy == 0) & (p0.component == 1)
    coref = (pf.galaxy == 0) & (pf.component == 1)

    center0 = np.average(p0.x[core0], axis=0, weights=p0.m[core0])
    centerf = np.average(pf.x[coref], axis=0, weights=pf.m[coref])

    # Initial radii of Galaxy A's disk particles from its initial core
    r0 = np.linalg.norm(p0.x[disk0] - center0, axis=1)

    # Use e.g. 95th percentile (or np.max(r0)) as the initial disk extent,
    # then flag particles pushed past `multiplier` * initial disk radius.
    initial_disk_extent = float(np.percentile(r0, radius_percentile))
    r_threshold = float(multiplier) * initial_disk_extent

    # Final 3D radii of the same disk particles from Galaxy A's final core
    rf = np.linalg.norm(pf.x[disk0] - centerf, axis=1)

    displaced_mask = rf > r_threshold
    displaced_percent = float(np.mean(displaced_mask) * 100.0)
    return displaced_percent, displaced_mask, r_threshold, r0, rf

def measure_encounter_physics(p0, pf, rd, eps=EPS, displace_mult=2.0, radius_percentile=95.0):
    """Return orbital-energy change and clearly labelled disk diagnostics."""

    def orbital_energy(p):
        positions, velocities, masses = [], [], []
        for gid in (0, 1):
            core = (p.galaxy == gid) & (p.component == 1)
            all_g = p.galaxy == gid
            positions.append(np.average(p.x[core], axis=0, weights=p.m[core]))
            velocities.append(np.average(p.v[core], axis=0, weights=p.m[core]))
            masses.append(np.sum(p.m[all_g]))
        mu = masses[0] * masses[1] / sum(masses)
        vrel = np.linalg.norm(velocities[1] - velocities[0])
        rsep = np.sqrt(np.sum((positions[1] - positions[0])**2) + eps**2)
        return 0.5 * mu * vrel**2 - G * masses[0] * masses[1] / rsep

    delta_e_orb = orbital_energy(pf) - orbital_energy(p0)

    disk0 = (p0.galaxy == 0) & (p0.component == 0)
    core0 = (p0.galaxy == 0) & (p0.component == 1)
    coref = (pf.galaxy == 0) & (pf.component == 1)

    center0 = np.average(p0.x[core0], axis=0, weights=p0.m[core0])
    centerf = np.average(pf.x[coref], axis=0, weights=pf.m[coref])

    # Initial radii of Galaxy A's disk particles relative to initial core center
    r0 = np.linalg.norm(p0.x[disk0] - center0, axis=1)
    r_disk_95 = float(np.percentile(r0, radius_percentile))
    r_threshold = float(displace_mult) * r_disk_95

    # Final 3D displacement from final core center (not cylindrical disk-plane radius)
    r_final_3d = np.linalg.norm(pf.x[disk0] - centerf, axis=1)
    displaced_mask = r_final_3d > r_threshold
    displaced_fraction = float(np.mean(displaced_mask) * 100.0)

    # Indices of displaced particles in the full snapshot arrays (useful if you want to inspect them)
    disk0_indices = np.flatnonzero(disk0)
    displaced_particle_indices = disk0_indices[displaced_mask]

    e0 = specific_binding_energy(p0, galaxy_id=0, eps=eps)
    ef = specific_binding_energy(pf, galaxy_id=0, eps=eps)
    mean_delta_eps = float(np.mean(ef - e0))
    eps90 = float(np.percentile(ef, 90))
    positive_energy_fraction = float(np.mean(ef > 0) * 100.0)

    return (
        float(delta_e_orb), displaced_fraction, mean_delta_eps,
        eps90, positive_energy_fraction
    )


# ---------------------------------------------------------------------
# One worker per encounter
# ---------------------------------------------------------------------

def run_single_encounter(args):
    q, inc, eps, theta, dt, steps, save_every, n_primary = args
    initial, cfg, counts_a, counts_b = make_initial_pair(q, inc, n_primary)
    rd = disk_scale_length(cfg)

    force = partial(fast_tree, eps=eps, theta=theta)
    p0, pf, separation_curve = run_encounter(
        initial, force, dt, steps, save_every
    )

    diagnostics = measure_encounter_physics(p0, pf, rd, eps=eps)
    return q, inc, *diagnostics, separation_curve, len(initial.m), rd


# ---------------------------------------------------------------------
# Study 1: energy conservation vs theta
# ---------------------------------------------------------------------

def study_energy_vs_theta(eps=EPS, dt=DT, steps=STEPS, save_every=SAVE_EVERY,
                          n_primary=N_PRIMARY):
    print("Study 1: Barnes-Hut energy conservation vs opening angle")
    initial, _, _, _ = make_initial_pair(1.0, 0, n_primary)

    # Warm up JIT compilation before timing; the force remains Barnes-Hut.
    fast_tree(initial.x[:8], initial.m[:8], eps=eps, theta=0.5)

    thetas = [0.2, 0.4, 0.6, 0.8, 1.0]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
    max_drifts, runtimes = [], []

    for theta in thetas:
        force = partial(fast_tree, eps=eps, theta=theta)
        t0 = perf_counter()
        snapshots = list(evolve(
            copy.deepcopy(initial), force, dt, steps, save_every
        ))
        elapsed = perf_counter() - t0

        times = np.array([t for t, _ in snapshots])
        energies = np.array([total_energy(p, eps) for _, p in snapshots])
        drift = np.abs(energies - energies[0]) / max(abs(energies[0]), 1e-30)

        ax1.plot(times, np.maximum(drift * 100, 1e-12),
                 lw=1.7, label=f"BH theta={theta:.1f}")
        max_drifts.append(float(np.max(drift) * 100))
        runtimes.append(elapsed)
        print(
            f"theta={theta:.1f}: max |dE/E0|={max_drifts[-1]:.5g}% "
            f"| runtime={elapsed:.3f}s"
        )

    ax1.set_yscale("log")
    ax1.set_xlabel("Time (code units)")
    ax1.set_ylabel("Absolute fractional energy drift (%)")
    ax1.set_title("Energy conservation over time")
    ax1.grid(alpha=0.25)
    ax1.legend(fontsize=8)

    ax2.plot(thetas, max_drifts, "o-", label="Maximum energy drift (%)")
    ax2.set_xlabel("Barnes-Hut opening angle theta")
    ax2.set_ylabel("Maximum energy drift (%)")
    ax2.grid(alpha=0.25)
    ax2r = ax2.twinx()
    ax2r.plot(thetas, runtimes, "s--", color="tab:red", label="Runtime")
    ax2r.set_ylabel("Runtime (s)")
    ax2.set_title("Accuracy/runtime trade-off")

    fig.tight_layout()
    fig.savefig(OUT / "energy_vs_theta.png", dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------
# Study 2: inclination vs mass ratio
# ---------------------------------------------------------------------

def study_encounter_grid(eps=EPS, theta=THETA, dt=DT, steps=STEPS,
                         save_every=SAVE_EVERY, max_workers=MAX_WORKERS,
                         n_primary=N_PRIMARY):
    print(
        f"Study 2: inclination/mass-ratio sweep | "
        f"{steps} steps | primary N={n_primary}"
    )

    dE = np.zeros((len(MASS_RATIOS), len(INCLINATIONS)))
    outside = np.zeros_like(dE)
    deps = np.zeros_like(dE)
    eps90 = np.zeros_like(dE)
    positive = np.zeros_like(dE)
    sep_curves = {}

    jobs = [
        (q, inc, eps, theta, dt, steps, save_every, n_primary)
        for q in MASS_RATIOS for inc in INCLINATIONS
    ]
    qi = {q: i for i, q in enumerate(MASS_RATIOS)}
    ii = {inc: i for i, inc in enumerate(INCLINATIONS)}

    with ProcessPoolExecutor(
        max_workers=max_workers,
        mp_context=mp.get_context("spawn")
    ) as executor:
        futures = [executor.submit(run_single_encounter, job) for job in jobs]

        for k, future in enumerate(as_completed(futures), start=1):
            (q, inc, de, outfrac, md, e90, posfrac,
             curve, n_total, rd) = future.result()

            r, c = qi[q], ii[inc]
            dE[r, c] = de
            outside[r, c] = outfrac
            deps[r, c] = md
            eps90[r, c] = e90
            positive[r, c] = posfrac

            if q == 1.0 and inc in (0, 90, 180):
                sep_curves[f"i={inc} deg"] = curve

            print(
                f"[{k}/{len(jobs)}] q={q:.2f}, i={inc:3d} deg, "
                f"N={n_total}, Rd={rd:.5f}, "
                f"deltaEorb={de:+.4g}, outside(3.5Rd)={outfrac:.1f}%"
            )

    # Save raw arrays in CSV for reproducibility.
    csv_path = OUT / "encounter_grid.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "q", "inclination_deg", "delta_E_orb",
            "outside_3p5Rd_percent", "mean_delta_disk_specific_energy",
            "final_disk_eps90", "positive_disk_energy_percent"
        ])
        for r, q in enumerate(MASS_RATIOS):
            for c, inc in enumerate(INCLINATIONS):
                writer.writerow([
                    q, inc, dE[r, c], outside[r, c], deps[r, c],
                    eps90[r, c], positive[r, c]
                ])

    fig, axs = plt.subplots(3, 2, figsize=(12, 9))
    axs = axs.ravel()

    def plot_by_q(ax, grid, title, ylabel):
        for r, q in enumerate(MASS_RATIOS):
            ax.plot(INCLINATIONS, grid[r], marker="o", lw=1.7,
                    label=f"q={q:.2f}")
        ax.set_xlabel("Inclination (degrees)")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)

    plot_by_q(axs[0], dE, "Change in approximate orbital energy",
              "Delta E_orb (code units)")
    plot_by_q(axs[1], outside, "Disk particles outside 3.5 Rd",
              "Particles outside threshold (%)")

    for label, curve in sep_curves.items():
        times = [x[0] for x in curve]
        seps = [x[1] for x in curve]
        axs[2].plot(times, seps, label=label)
    axs[2].set_xlabel("Time (code units)")
    axs[2].set_ylabel("Bulge-core separation")
    axs[2].set_title("Equal-mass core separation")
    axs[2].grid(alpha=0.25)
    axs[2].legend(fontsize=8)

    plot_by_q(axs[3], deps, "Mean disk specific-energy change",
              "Mean Delta epsilon")
    plot_by_q(axs[4], positive, "Positive disk-specific-energy fraction",
              "Positive energy (%)")
    axs[5].axis("off")

    fig.tight_layout()
    fig.savefig(OUT / "encounter_2d_sweep.png", dpi=150)
    plt.close(fig)

    print(f"Saved plots and data in {OUT}")


# ---------------------------------------------------------------------
# Simple Barnes-Hut scaling benchmark (BH only)
# ---------------------------------------------------------------------

def benchmark_barnes_hut_scaling():
    """
    Times the project's Barnes-Hut force call for increasing N.
    This measures this implementation's scaling, NOT the crossover against
    direct summation. Run locally on the target machine.
    """
    rng = np.random.default_rng(1234)
    sizes = [256, 512, 1024, 2048, 4096, 8192]
    rows = []

    # Compile/warm up before collecting timings.
    x0 = rng.normal(size=(32, 3))
    m0 = np.full(32, 1.0 / 32)
    fast_tree(x0, m0, eps=EPS, theta=THETA)

    for n in sizes:
        x = rng.normal(size=(n, 3))
        m = np.full(n, 1.0 / n)

        # One warmup call, then three measured calls.
        fast_tree(x, m, eps=EPS, theta=THETA)
        samples = []
        for _ in range(3):
            t0 = perf_counter()
            fast_tree(x, m, eps=EPS, theta=THETA)
            samples.append(perf_counter() - t0)

        median = float(np.median(samples))
        rows.append((n, median))
        print(f"Barnes-Hut N={n:5d}: median force-call time={median:.5f}s")

    with (OUT / "barnes_hut_scaling.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["N", "median_force_call_seconds"])
        writer.writerows(rows)

    ns, ts = zip(*rows)
    plt.figure(figsize=(6.5, 4.2))
    plt.plot(ns, ts, "o-")
    plt.xlabel("Particle count N")
    plt.ylabel("Median Barnes-Hut force-call time (s)")
    plt.title("Barnes-Hut runtime scaling on this machine")
    plt.grid(alpha=0.25)
    plt.tight_layout()
    plt.savefig(OUT / "barnes_hut_scaling.png", dpi=150)
    plt.close()


if __name__ == "__main__":
    print(f"Preset: {PRESET_PATH}")
    print(f"dt={DT}, steps={STEPS}, simulated duration={DT * STEPS:.3f}")
    print(f"Primary particle target N={N_PRIMARY}; q={MASS_RATIOS}")
    print(f"Barnes-Hut theta={THETA}; workers={MAX_WORKERS}")
    benchmark_barnes_hut_scaling()
    study_energy_vs_theta()
    study_encounter_grid()
