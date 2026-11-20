# Reproducible N-body validation suite

Extract the `validation_suite` folder into `PHYS332-FP/Ethan`. No existing project files need replacing. Commands below are PowerShell commands run from `Ethan`.

This suite measures evidence, checks controlled comparisons, and applies explicitly chosen tolerances. It does not turn a pretty image, energy conservation, or a single stability screen into a claim that the Mice model is physically validated.

## Start here: actual native solver verification

```powershell
.\.venv\Scripts\python.exe -m pip install -r .\validation_suite\requirements.txt
$out = "results\solver_selftest_" + (Get-Date -Format "yyyyMMdd_HHmmss")
.\.venv\Scripts\python.exe .\validation_suite\validate.py selftest --repo-root . --out $out
```

This imports your real `nbody.gravity.parallel_direct` and `nbody.integrate` modules. It checks direct forces against an independent vectorized sum, then evolves a softened circular two-body orbit at two timesteps against SciPy DOP853. The force error, orbit errors and convergence ratio are saved. These are numerical test fixtures with stated tolerances; they do **not** test Zeno or GADGET executables.

## Inspect a trajectory

Select the file explicitly; the newest file is not necessarily the selected encounter.

```powershell
$trajectory = (Read-Host "Paste the full path to run_trajectory.npz").Trim('"')
$out = "results\trajectory_audit_" + (Get-Date -Format "yyyyMMdd_HHmmss")
.\.venv\Scripts\python.exe .\validation_suite\validate.py inspect $trajectory --context encounter --energy-samples 9 --self-gravity-only --out $out
```

Use `--context isolated` **only** for a trajectory evolved with no companion. Use `unspecified` when unsure. Without `--energy-samples`, the fast inspection skips expensive pairwise energy calculations and labels energy evidence missing. Full energy is O(N^2) per sampled frame; nine saved-frame values do not bound all-step drift.

`--self-gravity-only` is an explicit assertion about the run being analysed. Do not use it on disk tracers in a fixed smooth field, gas runs, externally forced systems or cosmological trajectories. Native controls generated here record this dynamics model automatically. Other models require their own energy terms; no pair-only energy pass is silently assigned to them.

Outputs include `report.json`, readable `report.md`, plots, conservation history, sampled energy, disk structural history and radial diagnostics. Disk/component labels follow the current project: 0=disk, 1=bulge, 2=halo.

Required trajectory NPZ fields:

- `t`: [frames], strictly increasing code times.
- `x`, `v`: [frames,N,3], positions and velocities at the **same** times.
- `m`: [N], fixed positive particle masses; particle ordering must remain constant between frames.
- `component`, `galaxy`: [N] integer labels for component diagnostics.
- `metadata`: a scalar JSON string recording force law, softening, G, dt and solver identity when available.

Snapshots with `x,v` [N,3] and no `t` are accepted by `initial`, `forces` and `run_control.py`.

`m` cannot vary across time in this version. Gas conversion/star formation, changing particle counts, rotating coordinates and cosmological/comoving data are outside the supported format. Do not feed a GADGET gas run into this collisionless energy calculation without a separate adapter and appropriate thermodynamic terms.

## Generate a direct-summation isolation control

```powershell
$out = "results\isolation_direct_" + (Get-Date -Format "yyyyMMdd_HHmmss")
.\.venv\Scripts\python.exe .\validation_suite\run_control.py $trajectory --repo-root . --mode isolate --solver direct --duration 1.25 --dt 0.001953125 --threads 4 --out $out
.\.venv\Scripts\python.exe .\validation_suite\validate.py inspect "$out\run_trajectory.npz" --context isolated --energy-samples 9 --out "${out}_analysis"
```

The runner extracts galaxy 0 (change `--galaxy` if needed), removes uniform bulk motion, and preserves masses, labels, softening and individual internal velocities. It prefers an adjacent original `initial_particles.npz` after checking it against replay frame 0. Otherwise it explicitly records use of the replay's initial precision.

The native direct solver and native leapfrog are used. A warmed-up full-size force timing and runtime forecast are printed before evolution. Checkpoints survive interruption, but this runner does not resume checkpoints automatically. It writes float64 trajectories. Use `--prepare-only` to inspect the initial state without evolving anything.

For timestep refinement, repeat the **same input** with `--dt 0.0009765625`; keep duration and save interval unchanged. For a native tree comparison use `--solver tree --theta 0.4`. Use `--mode replay-initial` to retain the full encounter instead of extracting one galaxy. No additional encounter geometry is regenerated.

## Compare two runs

```powershell
$a = (Read-Host "Path to baseline trajectory").Trim('"')
$b = (Read-Host "Path to refined trajectory").Trim('"')
$out = "results\timestep_comparison_" + (Get-Date -Format "yyyyMMdd_HHmmss")
.\.venv\Scripts\python.exe .\validation_suite\validate.py compare $a $b --kind timestep --out $out
```

Kinds: `timestep`, `solver`, `resolution`, `seed`.

- Timestep/solver comparisons require matching initial particles, masses, labels and force law, allowing only a uniform position/velocity offset and replay precision. Timestep comparisons require the same recorded solver identity and theta, with different recorded dt.
- Solver comparisons require the same recorded dt. Different tree criteria/solvers are permitted.
- Resolution/seed comparisons check component masses and flag absent parent-construction provenance. Common construction is not established by similar-looking snapshots. The optional `validation_parent_model_id` is an explicit provenance declaration, not a hash proof by itself.
- Only common saved times are compared. No temporal interpolation, orbital-phase shifting or morphology alignment is done.
- Scalar structural differences are reported. This does not calculate image-tail residuals or observational scores.

## Check force errors on any snapshot

```powershell
.\.venv\Scripts\python.exe .\validation_suite\export_forces.py $trajectory --repo-root . --solver direct --frame 0 --samples 96 --out .\results\native_direct_forces.npz
.\.venv\Scripts\python.exe .\validation_suite\validate.py forces $trajectory .\results\native_direct_forces.npz --frame 0 --out .\results\native_force_check
```

The exporter uses native project code and selects radial-stratified plus random targets. It evaluates the full source distribution. For the tree, use `--solver tree --theta 0.5` and a different output name.

For an external engine, supply an NPZ containing `indices` [S], `a` [S,3], and scalar JSON `metadata` with `trajectory_sha256`, `frame`, `epsilon`, `G`, `kernel`, plus solver identity. Accelerations must be actual same-frame engine accelerations, not velocity differences across widely spaced outputs. No external-engine force adapter is silently substituted.

## Initial-condition target checks

```powershell
.\.venv\Scripts\python.exe .\validation_suite\validate.py initial .\path\initial_particles.npz .\validation_suite\examples\targets.template.json --out .\results\initial_target_check
```

Edit the target JSON for the **actual** construction. The supplied component masses are a labelled Mice-style example, not universal galaxy parameters. Add `cdf_csv` per component to compare a target cumulative radial profile. CSV columns: `radius,cumulative_fraction`; values must include your taper/truncation and normalization. Coordinate is `disk_radius` or `spherical_radius`. The maximum discrepancy is on the table grid, not an exact KS statistic or p-value.

Optional `moments_csv` per component compares initial velocity moments in radial bins. Columns:

```text
R_lo,R_hi,rotation_abs,sigma_R,sigma_z
```

Leave individual target moment cells empty to omit that check. Threshold keys are `initial_rotation_abs_relative_error`, `initial_sigma_R_relative_error` and `initial_sigma_z_relative_error`. Rotation is absolute because fitted-plane eigenvectors have arbitrary signs; test physical spin orientations separately from the encounter manifest. Neither spatial nor moment agreement establishes a stationary full distribution function.

## Compare observed spatial regions and velocities

```powershell
.\.venv\Scripts\python.exe .\validation_suite\validate.py observations $trajectory .\path\observation_setup.json --out .\results\observational_comparison
```

Start with the **empty format templates** in `examples`. No real observational velocity measurements or calibrated scales are invented or bundled.

The setup specifies one time, camera, roll, length scale (kpc/code-length), velocity scale (km/s/code-velocity), systemic velocity and sky offset. Origin is the whole system's instantaneous mass COM, with system bulk velocity removed. The camera convention is X right, Y down, positive LOS velocity away from the observer; the roll rotates both spatial axes. Transform measured sky coordinates to this convention before comparing. The command does not fit a camera or adopt the older Hubble nucleus registration implicitly.

CSV columns: `id,role,kind,galaxy,x,y,radius,vlos,sigma_v`.

- Each region is a circular sky aperture in kpc, associated with one labelled galaxy (or -1 for either).
- `kind=tail` uses disk particles; `kind=nucleus` uses bulge particles.
- `role=fit` means used to fit the model. `role=holdout` must actually have been withheld.
- Empty `vlos,sigma_v` supplies morphology-only coverage. A supplied velocity needs a positive measurement uncertainty.
- Reported quantities: population, stellar mass fraction, weighted median velocity, 16th/84th velocity percentiles, residual, and residual divided by observational sigma.

Mass fraction is not optical intensity or HI emissivity. Median residual/sigma is a descriptive diagnostic, not a calibrated likelihood. Broad model velocity distributions are reported separately and must not be treated as a free improvement in fit. Parameters and observation provenance must be documented externally; do not designate fit points as independent holdouts after selecting a model with them.

Time and mass scales are derived consistently from L and V. This release does not digitize published figures, identify tails automatically, fit observational parameters, calculate a calibrated chi-square, or certify unique encounter parameters. Those tasks require actual measurements and a declared fitting protocol.

## Decisions and tolerance policy

Statuses: `PASS`, `FAIL`, `MEASURED`, `MISSING`, `INCONCLUSIVE`.

- Hard format errors return exit code 3 and do not produce a pass report.
- A failed declared check returns exit code 2.
- Exit code 0 means successful analysis, **not physical validation**. Read `status` and individual checks.
- Without a policy, physical quantities remain measured. Missing evidence remains missing.
- `examples/policy.example.json` contains explicitly **unreviewed illustrative engineering limits**. It deliberately does not activate pass/fail gates. Review the scope and scientific justification, edit the numbers/rationale, and set `reviewed:true` only after choosing a policy. Then supply `--policy your_policy.json`.
- Report-level `PASS_FOR_SPECIFIED_CHECKS` is restricted to that command. There is no whole-project green badge. A threshold crossing is not automatically proof that the underlying physical model is wrong.

For isolation, the configurable gates use maximum departure during the tested interval, not just the endpoint. They apply to both whole-disk and fixed inner subsets. Collective evolution may be physical; investigate it rather than treating every threshold failure as a numerical bug. The radial diagnostics include mean z/vz and within-cell residual heights, but cannot uniquely separate collisional heating from collective modes.

All output directories must be new. Existing output is never overwritten. Reports record input paths, SHA256 hashes, metadata, settings, measured values, and applied policies. Sampling/observer target CSV files are hashed as inputs too. Negative results are retained.

## Batch execution

Create a JSON file with jobs such as:

```json
{"jobs":[
 {"args":["inspect","results/isolated/run_trajectory.npz","--context","isolated"]},
 {"args":["compare","results/isolated/run_trajectory.npz","results/halfdt/run_trajectory.npz","--kind","timestep"]}
]}
```

Paths in job arguments are relative to the manifest directory. Omit `--out` in each job; batch assigns separate folders. Only analysis commands and `selftest` are allowed, not shell commands or automatic simulations.

```powershell
.\.venv\Scripts\python.exe .\validation_suite\validate.py batch .\validation_jobs.json --out .\results\validation_batch
```

## Tested release evidence

`sample_reports` contains actual analyses of the supplied selected encounter and matched isolation timesteps, the native solver selftest, and an independent force comparison. These outputs are examples from the supplied data, not thresholds for future galaxies. Absolute paths in them refer to the analysis environment.

Run package regression tests:

```powershell
.\.venv\Scripts\python.exe -m unittest discover -s .\validation_suite\tests -v
```

## Scientific scope and references

The suite supports isolated, collisionless, fixed-mass Plummer-softened G=specified simulations in an inertial frame. It checks the softened model actually evolved; it does not imply unsoftened gravity is sufficiently resolved on every relevant scale. Use a resolution/softening study for that claim, regenerating consistent initial conditions when appropriate.

- Barnes (2012), softened fields and equilibrium: https://arxiv.org/abs/1205.2729
- Barnes & Hibbard (2009), Identikit construction and model matching: https://arxiv.org/abs/0811.3039
- Privon et al. (2013), Mice sky and position-velocity comparison: https://arxiv.org/abs/1305.5865
- Mortazavi et al., observational region scoring and tracer limitations: https://arxiv.org/abs/1701.03819
- Hibbard & van Gorkom observations: https://www.cv.nrao.edu/~jhibbard/MergeSeq/

Passing numerical tests, construction targets, isolated controls and convergence checks establishes evidence about the calculation. Observational consistency additionally requires measured morphology and velocities, their uncertainties, one consistent physical calibration, and explicit model limitations. The suite keeps those claims separate.
