# Equilibrium repair checkpoint — 19 November 2026

Base repository: e9705496bfb6215988481d21dbdbcde488153d60.
These are local proposed changes, not a published release or a certification
that the full live galaxy has reached equilibrium.

## What changes

1. `nbody.diagnostics.structure` measures disk thickness and radial dispersion
   in its fitted plane. `z_rms_fixed_axis` retains the former axis-based height;
   `sigma_z` and `tilt_degrees` are now recorded in isolation CSVs.
2. `zeno_smooth_control.py` integrates the existing reference disk in a verified
   fixed softened continuum potential. It has no dependency on compiled AGAMA.
3. Its opt-in `--vertical-jeans` builds a separate experimental candidate with
   sigma_z^2(R,z) = integral_z^infinity rho(R,z') g_z(R,z') dz' / rho(R,z).
   This assumes zero radial-vertical stress in the initial target. It changes
   only disk v_z; positions, in-plane velocities and spheroid particles remain
   unchanged. It is a local second-moment correction, not an exact 3D DF.
4. `reference_live_control.py` tests that candidate's live disk in fixed smooth
   bulge and halo fields, using this project's own direct or tree solver. It
   prints a warmed-up runtime forecast and writes checkpoints. This control
   removes both spheroid graininess and response; it is not a merger model.
5. Saved candidates with an explicit epsilon cannot be silently run at another
   epsilon in the GUI or encounter runner. Legacy annotations remain labelled
   as unrecorded construction force laws. Newtonian AGAMA construction remains
   explicitly distinct from softened construction and live stability.
6. The Zeno generator exposes particle counts and seeds, defaults to 65,536
   total particles, and writes source/patch hashes and a model manifest.
7. Encounter inputs are validated before creating the run directory, so an
   invalid input does not consume an output name.

## Reproduced result

Using the supplied full-resolution original initial.txt (24,576 disk, 8,192
bulge, 32,768 halo particles), at t=1.5:

| Experiment | RMS fitted-plane thickness change |
|---|---:|
| Existing fully live Zeno run | +14.6458% |
| Original particles, fixed smooth field, dt=1/1024 | +3.9544% |
| Jeans-moment-adjusted candidate, fixed field, dt=1/1024 | +1.9789% |

The adjusted candidate gives +1.9970% at dt=1/512; endpoint thickness differs
by 0.0178% relative between timesteps. Maximum sampled vertical force error
against refined quadrature is approximately 0.00505%. This demonstrates a
partial repair of the smooth-field transient. It does not establish live
stability or assign the live growth uniquely to particle scattering.

No original cached galaxy is overwritten, and the reference generator's
published sheet recipe is retained as a reproducible baseline. The adjustment
is explicitly opt-in. Larger particle counts are not described as a universal
fix, and no damping or unconstrained pre-relaxation is applied automatically.

## Windows PowerShell: reproduce the candidate

From `PHYS332-FP/Ethan`, install optional science dependencies:

```powershell
.\.venv\Scripts\python.exe -m pip install -r requirements-science.txt
```

Run the correction and both smooth-field timestep checks. The input is the
original full-resolution initial.txt in WSL; the output folder must be new:

```powershell
.\.venv\Scripts\python.exe .\zeno_smooth_control.py `
  --initial "\\wsl.localhost\Ubuntu\home\whisky\nbody-tools\zeno-reference\full-resolution-20261119-032100\initial.txt" `
  --vertical-jeans `
  --out results\reference_jeans_candidate
```

Outputs include `jeans_candidate.npz`, `smooth_spheroid_field.npz`,
`field_checks.json`, `base_history.csv`, `half_dt_history.csv`, and `summary.json`.
The input geometry and component masses must correspond to the documented
reference model; the script deliberately does not accept arbitrary galaxies.

## Next live-disk control

The default 64 steps are a short screen. The runtime forecast is measured
locally before evolution; increasing to 768 steps matches t=1.5 at dt=1/512.
Do not compare the short screen's endpoint directly with the long live run.

```powershell
.\.venv\Scripts\python.exe .\reference_live_control.py `
  results\reference_jeans_candidate\jeans_candidate.npz `
  results\reference_jeans_candidate\smooth_spheroid_field.npz `
  --solver direct --threads 4 --steps 64 `
  --out results\reference_live_disk_screen
```

The fully live corrected candidate must also be tested from this exact saved
snapshot, with matched epsilon and a recorded solver. The GUI can load the
candidate after it is copied into `results/galaxies`, but smooth-field stability
does not certify a fully live encounter. Do not use the existing AGAMA-only
`rigid_spheroid_control.py` on this Zeno candidate.

## Remaining work

- Matched live-disk and fully live tests of the corrected candidate.
- Separate halo/disk resolution and independent-seed checks over the required
  encounter duration; study collective amplitudes alongside dispersion.
- A fully stationary DF may be preferable if the remaining smooth adjustment
  is unacceptable. The local Jeans correction alone cannot guarantee it.
- No change to leapfrog or the pair-force kernel is supported by the completed
  short timestep and opening-angle checks.

## Sources

- Barnes & Hibbard (2009), Appendix B:
  https://doi.org/10.1088/0004-6256/137/2/3071
- Barnes (2012), softened fields and equilibrium:
  https://arxiv.org/abs/1205.2729
- McMillan & Dehnen (2007), staged disk tests:
  https://arxiv.org/abs/astro-ph/0703534
- Rodionov & Sotnikova (2006), limits of moment-based construction:
  https://arxiv.org/abs/astro-ph/0611176
