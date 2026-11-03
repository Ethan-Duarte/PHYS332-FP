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
