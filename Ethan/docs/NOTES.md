# Architecture and decision register, Version 1

| Decision | Reason | Status |
|---|---|---|
| Work from fresh source layout, reuse recovered gravity algorithms | Minimize coupling to old demos | Implemented |
| Use arrays `(x,v,m)` and persistent `(component,galaxy)` labels | Any force solver or analyzer can consume one structure | Implemented |
| Force function `(x,m) -> accelerations` | Leapfrog doesn't know tree internals | Implemented |
| Separate profile math, sampling, encounter geometry | Extensibility and mathematical validation | Implemented |
| Native exact Plummer DF as first equilibrium test | Truly consistent spherical phase space, simple probability distributions | Implemented, basic tests |
| Native Mice preview deliberately approximate | Allow full end-to-end demo without claiming realistic stellar equilibrium | Implemented, **not validated** |
| Optional AGAMA for disk-bulge-halo equilibrium | Distribution-function-based self-consistent modeling | Written, not runtime tested |
| `G=1`, explicit Barnes-style preset | Consistent with earlier solver and paper | Implemented |
| Seeded native sampling, exact snapshots saved | Reproducibility | Implemented |
| No front end or general plugin architecture yet | Avoid premature complexity | Deferred |
| Target system: the Mice (NGC 4676) | Published model and observational tidal tails | Selected |

Source links: Barnes (2004) https://doi.org/10.1111/j.1365-2966.2004.07725.x ; AGAMA https://github.com/GalacticDynamics-Oxford/Agama ; Vasiliev (2019) https://doi.org/10.1093/mnras/sty2672 ; Dehnen (1993) https://doi.org/10.1093/mnras/265.1.250.

## V1.1 animation / resolution decision

- Keep the original small `mice_preview.py` as a quick smoke test.
- `experiments/mice_animation.py` evolves high-count particles and saves a trajectory for repeatable playback. Dark matter remains live, but is hidden by default in the visual display to emphasize stellar tidal features.
- Add an **optional** Numba-accelerated direct-summation implementation to obtain high-resolution teaching animations without changing the original algorithms. The force law is identical to `gravity/direct.py`; benchmarks must still compare the named implementations honestly.
- An attractive animation is **not** a validation test. Native disk kinematics remain approximate, as described in `README.md` and `docs/VALIDATION.md`.

### V1.3 update, 2026-11-03

GUI enhancements completed: cancel event, compiled Barnes–Hut force traversal (compared with original), interactive 3D Plotly viewer with self-contained HTML export, arbitrary initial relative position/velocity vectors, offline trajectory replay, and per-snapshot force-comparison tab. AGAMA equilibrium initial conditions are **not** validated or enabled. A working animated collision is not yet evidence of an astrophysically validated Mice reconstruction. Full-resolution convergence and observational reconstruction remain the next core milestones.

### V1.4 — 2026-11-05

- Retain the V1.3 interactive viewer, cancellation, and all gravitational solvers. New experiment mode: isolated galaxy (Plummer analytical DF / disk–bulge–halo native preview / optional AGAMA candidate / saved snapshot).
- Add `diagnostics.py` for isolated structural drift, softened energy and sampled-vs-analytic density CDFs, and explicit **screening pass/fail**, never silently equating this to validation.
- Add matched solver evolution (`solver_study.py`) using identical initial particles. Distinguish force accuracy, numerical energy drift, full evolution time, and final particle-position discrepancies; early performance results favor compiled direct summation at low N.
- Freeze Mice reference **nominal** parameters: Barnes 2004 parabolic rp=0.25, disks (25°, −30°) and (40°, 60°), eps=0.0125, dt=1/512. Our initial separation and viewing frame are still development choices, so observational fit is incomplete.
- AGAMA code adapted toward documented DF API; extension unavailable here. Its outputs remain equilibrium **candidates** until masses, units, kinematics and long-term isolated stability are measured.
- Use saved model snapshots to avoid repeated expensive AGAMA generation. When reused for both progenitors the particle noise is correlated; future releases can support independent seeded samples from one converged DF.
- Defer GUI website deployment, other galaxy types and multiple-body encounters until a robust physical baseline exists.

### V1.4.1 — 2026-11-06

- Added annular, mass-weighted Fourier modes `A_m(R,t)`, `m=1..4`, plus phase and finite-N RMS noise floor. Average the **absolute annular amplitudes**, not the complex amplitudes over the whole disk, to avoid spiral-phase cancellation. Reuse the exact saved trajectory for reanalysis.
- Measure `sigma_R` and disk thickness alongside `A_m`. Estimate Toomre `Q` with κ inferred from azimuthally averaged softened live-particle force, but label it **approximate**, especially for non-equilibrium disks. Low particle counts may make annuli unresolvable.
- Added live AGAMA import status in the GUI, identifying the Python interpreter; the compiled AGAMA dependency remains absent in our environment and its adapter remains unverified. No fallback to preview models on import failure.
- **Validation unchanged:** native disk velocity fields remain heuristic and have demonstrated major relaxation; amplitude of a transient `m=2` feature does not prove a true spiral pattern. Research goal is still a stable AGAMA-derived reference disk for V1.5, then the Mice comparison.
