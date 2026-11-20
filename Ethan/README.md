# Final planned isolated stellar-resolution comparison

Extract the files directly into your existing `Ethan` directory, alongside `zeno_mice_run.py`, `zeno_snapshot_bridge.c` and the `nbody` package. Replace the two supplied helper scripts if asked. Do not run the preparation audit script: the initial state is already supplied.

From PowerShell in Ethan:

```powershell
$finalOut = "results\final_stellar_resolution_" + (Get-Date -Format "yyyyMMdd_HHmmss")

.\.venv\Scripts\python.exe .\run_final_stellar_test.py `
  --out $finalOut
```

The driver checks the prepared state, runs one isolated Zeno simulation, analyzes the saved trajectory, then creates `analysis_to_share.zip` inside the run folder. Send that ZIP. The full trajectory stays available in the same folder for rendering. No second simulation is launched by the analysis step.

## Fixed experimental settings

| Quantity | Choice | Purpose |
|---|---|---|
| Disk particles | 12,288 | Four times the previous sampling at fixed disk mass |
| Bulge particles | 4,096 | Four times the previous sampling at fixed bulge mass |
| Halo particles | 32,768 | Exact previous halo particles, velocities and masses retained before common recentering |
| Total particles | 49,152 | Targeted stellar refinement rather than increasing all components |
| Component masses | 0.1875, 0.0625, 1.0000047683731454 | Unchanged |
| Plummer epsilon | 0.009375 | Same construction/evolution setting as the baseline |
| Timestep | 1/512 = 0.001953125 | Same timestep; halving it previously did not cure growth |
| Tree theta | 0.7 | Same criterion; prior direct/tree structural comparison agreed closely |
| Quadrupoles | On | Same solver settings |
| Duration | 1.25 code units, 640 steps | Same endpoint as the completed baseline, including late mode growth |
| Saved interval | Every 8 steps, 81 scheduled frames including t=0 | Same cadence, with positions and velocities for analysis and later rendering |

The prepared state retains the original 3,072 disk and 1,024 bulge coordinates and velocities, and adds distinct particles from the verified same parent realization. Stellar particle masses decrease fourfold; they are not cloned particles. All component masses are preserved. Sampling uses a recorded fixed seed. The runner subtracts a common overall center and bulk velocity, which preserves internal distances and relative velocities.

This is one joint disk/bulge sampling comparison, not separate identification of their individual contributions or a seed ensemble. It can show whether the measured structural growth is sensitive to stellar sampling. A quieter result does not by itself prove an exact equilibrium distribution function.

No velocity damping, Q changes, force-kernel changes or new softening are introduced. Prograde versus retrograde orientation has no orbital meaning for this isolated test; encounter spin correction will be applied when constructing the animation's collision setup.

## What was checked before delivery

- Parent candidate identity and unique matching of retained stars to the parent realization.
- Exact halo-state retention before common recentering.
- Component total masses and distinct added parent IDs.
- Initial-state checksum, particle counts and finite coordinates/velocities.
- Driver preflight and Python compilation. The long final run was not executed here.
- The analysis helper was previously tested against the original isolation diagnostics and synthetic bending fields.

`prepare_final_stellar_control.py` and `final_stellar_initial.preparation.json` preserve the preparation procedure and source hashes for auditing. You do not need the large parent trajectory to run the supplied prepared NPZ.

After this comparison we move to animation development. Its measured limitations will be recorded rather than triggering an automatic sequence of more validation simulations.
