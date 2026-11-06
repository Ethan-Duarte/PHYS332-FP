# Verification and validation roadmap

**Current status: not observationally validated.** Basic tests do not establish galaxy equilibrium or a successful reconstruction of NGC 4676.

## Verification of numerical implementation

1. Gravity: one-/two-body analytic acceleration, action–reaction and net force zero.
2. Barnes–Hut: match direct forces as `theta` shrinks on random, spherical and thin-disk configurations; record median, 95th-percentile and absolute force errors.
3. Integrator: test an equal-mass circular binary; halve `dt` twice and check second-order global trajectory convergence.
4. Conservation: compare energy under the *same softened force law*, plus linear and angular momentum.
5. Timing: benchmark direct and tree `N=50..5000` (or practical limit), reporting crossover and memory constraints.

## Galaxy generation

6. Draw many independent initial positions; compare cumulative `M(<r)` for bulge/halo and `M(<R)` for disk with analytic curves including chosen cutoffs.
7. Check velocities: disk rotation curve, `sigma_R`, `sigma_phi`, `sigma_z`, bulge/halo anisotropy, and Toomre Q where relevant.
8. Evolve each galaxy *in isolation* for multiple local orbital periods. Track mass profiles, disk half-mass radius, thickness, dispersions and bar growth.
9. Repeat with finer `dt`, different `eps`, larger *halo* and stellar particle counts, and multiple seeds; numerical effects must be smaller than the proposed physical signal.

## Encounter comparison

10. Reproduce Barnes reference encounter with suitable equilibrium galaxies (inclinations 25°,40°, line-of-node angles −30°,60°, nominal pericenter 0.25); actual pericenter must be measured from simulated centers.
11. Match projected galaxy centers and optical tidal-tail outlines to observed Mice; use held-out line-of-sight kinematics where suitable data permit. A viewing transformation must be fitted explicitly.
12. Predefine tidal-tail metrics and sample at equal times after *actual* periapsis; use seed variability and convergence bands to defend conclusions.

Thresholds, adopted source measurements, and code/initial-condition revisions must be recorded before presenting results as validated.

## Version 1.4: implemented diagnostics and remaining acceptance gates

- `nbody/diagnostics.py`: vectorized/compiled softened pairwise energy, isolated disk/sphere structural time series, truncated/native profile CDF comparison, and exported CSV/PNG/JSON. A **screening pass** means only that selected threshold checks passed over that particular duration. The test must be repeated with different seeds, timesteps, force settings and particle counts.
- `nbody/solver_study.py`: initializes all solvers from the **same** particle realization, measures acceleration errors relative to compiled direct summation, compares full-trajectory time, endpoint energy and structural galaxy separation. Long-time per-particle RMS separation is not itself a physical convergence metric for chaotic encounters; a future version must quantify convergence of macroscopic tidal metrics.
- AGAMA is code-integrated as a selectable generator and automatically caches a snapshot, but the compiled extension is unavailable in the build environment. Do not mark the disk–bulge–halo galaxy as physically validated until this backend has run, shown consistent component masses/velocity distributions, and passed isolated stability at converged resolution.
- The Mice preset adopts Barnes's published nominal pericenter/orientations; it is **not** yet fitted to the observed Mice's sky projection, line-of-sight velocity gradients, or tail geometry. No observational validation is claimed.

### Initial screening demonstration (2026-11-05)

| Case | Particle count | Duration | Radius drift | Thickness drift | Max sampled energy drift | Screening |
|---|---:|---:|---:|---:|---:|---|
| Isolated analytic Plummer | 500 | 1.00 | −3.63% | N/A | 0.00229% | Pass |
| Isolated native preview disk–bulge–halo | 600 | 0.16 | +17.08% | +99.59% | 0.0179% | Fail |

The second case explicitly demonstrates that nearly conserved total energy does not imply approximate collisionless equilibrium. The contrast depends on the stated resolution, timestep, softening and sampling; it is not a universal comparison of galaxy families.
