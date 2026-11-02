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
