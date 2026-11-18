# Validation protocol — V1.5

## Separate the questions

- **Library/API verification:** can our Python code use AGAMA's actual compiled classes?
- **Initial-condition verification:** does the particle sample reproduce masses and phase-space
  properties of the AGAMA DF model to sampling uncertainty?
- **Solver verification:** does our N-body engine integrate the selected softened equations
  consistently as `dt`, Barnes–Hut `theta`, and particle count are changed?
- **Model validation:** is the realized galaxy sufficiently stable to interpret a merger?
- **Encounter validation:** does the evolved model reproduce appropriate Mice observables?

These are **different** levels of evidence. V1.5 provides the workflow for the first four,
but does not claim any AGAMA-dependent level has been achieved in this build environment.

## First-pass isolated test

For one initial particle snapshot, run an isolation experiment at common physical duration `T`:

- Base: `dt`, `steps`.
- Refinement: `dt/2`, `2*steps`.
- Save at the same physical times, using `save_every` and `2*save_every`.
- Keep `eps`, force law, particle masses, positions and velocities identical.

Compare:

\[
\Delta_{r50} = \frac{r_{50}(T)-r_{50}(0)}{r_{50}(0)}, \quad
\Delta_z = \frac{z_{\rm rms}(T)-z_{\rm rms}(0)}{z_{\rm rms}(0)}.
\]

The code also records soft-potential energy drift and RMS particle-position difference
between the two integrations.

**Draft numerical screens**, not natural physical tolerances: |Δr50| <= 5%,
|Δz| <= 10%, sampled max relative energy drift <= 1%. A disk with a real bar or spiral can
legitimately evolve; the thresholds are prompts to investigate, not a proof of failure.

## Baseline scientific validation grid (future work)

| Variation | Purpose | Proposed values |
|---|---|---|
| Stellar/halo resolution | Discern two-body artificial heating | at least 3 total-N levels, holding mass fractions |
| Independent random seeds | Measure realization scatter | at least 3 per chosen resolution |
| Time step | Check integrator accuracy | Δt, Δt/2, Δt/4 |
| Softening | Check regularization dependence | ε/2, ε, 2ε |
| BH opening angle | Tree approximation | 0.3, 0.5, 0.7, reference direct |
| Time span | Structural secular changes | at least several orbital periods at R=3Rd |

For each configuration report `(r50,z_rms,sigma_R)` and annular `A_m(t)` along with energy
and momentum diagnostics. Record actual orbital period and actual steps, not just code time.

The key question is not whether the model is perfectly stationary. It is whether
changes caused by unphysical relaxation are small compared with the *physical effects we
intend to measure*. In a stellar disk, finite-N halo particles can heat the disk.

## Why a tiny timestep error does not prove physical equilibrium

If both coarse and fine integrations agree that `z_rms` doubles, then time integration
is likely not the principal source of that change — but the starting galaxy is still
inadequate for our desired encounter interpretation.

The old V1.4.1 native preview showed exactly this kind of behavior: large changes in
disk thickness while total energy remained well conserved.

## Key files

- `nbody/galaxies/sources/agama.py`: smooth DF model and sampling adapter.
- `nbody/validation.py`: isolated two-timestep study of one fixed snapshot.
- `nbody/diagnostics.py`: morphology/energy measurements and provisional screens.
- `nbody/spirals.py`: annular Fourier modes, radial kinematics, approximate Q.
- `experiments/agama_selftest.py`: actual AGAMA import/construction smoke test.
- `experiments/validate_galaxy.py`: non-GUI isolation and numerical refinement run.
- GUI tab **Test a saved equilibrium galaxy**: launch comparison, collect reports.

## Mice encounter eligibility

Only send a galaxy to the Mice reference comparison if:

1. Its source (AGAMA release/commit, configuration and random seeds) is documented.
2. Its mass/position/velocity distributions have passed sampling checks.
3. Its isolation tests have sufficiently small **resolution-dependent** drift.
4. Its force solver's error and energy behavior are numerically characterized.
5. The intended numerical resolution covers the tidal structures that we plan to measure.

This is an acceptance protocol, **not** a blanket assertion that the current preset passes.
