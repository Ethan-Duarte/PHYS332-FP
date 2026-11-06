# N-Baddies V1.4 — Galaxy Dynamics Lab

A student-built gravitational N-body simulator with a local browser GUI, controllable 3D playback, isolated-galaxy diagnostic experiments, a Barnes-inspired Mice encounter, and direct-versus-Barnes–Hut trajectory comparisons.

## Launching on Windows

Extract the **entire** folder and double-click `start_gui.bat`. The launcher creates `.venv`, installs dependencies, and opens `http://127.0.0.1:7860`. Keep the console window open. Any Python startup traceback is saved to `startup_error.log`. The Windows batch file has the same corrected structure as V1.3; it has not been executed on Windows in this testing environment.

For terminal use in the project directory:

```bash
pip install -r requirements.txt -r requirements-fast.txt -r requirements-gui.txt -r requirements-viewer.txt
python -m pytest -q
python -m experiments.isolation_study --model plummer --particles 500
python -m experiments.isolation_study --model disk --backend preview --particles 1000
python -m experiments.solver_study --particles-per-galaxy 200
python -m experiments.mice_animation --disk 400 --bulge 50 --halo 200 --steps 160
```

AGAMA is a **separate optional compiled dependency**. The GUI can select it; there is no silent substitution with preview particles. See `docs/AGAMA_SETUP.md` for the API and installation cautions.

## GUI modes

1. **Isolated galaxy**: either Plummer (exact isotropic DF in the continuum limit) or disk–bulge–halo (native preview, AGAMA candidate, or a saved particle snapshot). Output includes structure-vs-time plots, initial profile checks for known galaxy configurations, numerical energy drift, and a *screening pass/fail* verdict. A screening pass is not scientific validation.
2. **Mice encounter**: Barnes (2004)-inspired disk, bulge and cored halo parameters, approximate parabolic initial orbit (nominal pericenter 0.25 in G=1 units), disk orientations 25°/−30° and 40°/60°. Initial separation and sky viewing orientation are not an observationally fitted reconstruction. The initial halo/stellar configurations are live, but there is no gas.
3. **Custom encounter**: initial relative position and velocity vectors, masses and orientations, with a saved trajectory and 3D playback.
4. **Compare gravitational solvers**: *one-force-snapshot* acceleration errors and runtimes.
5. **Full-trajectory solver study**: identical initial conditions evolved with compiled direct, compiled Barnes–Hut and NumPy direct (small N). Reports force errors, full evolution time, endpoint energy drift, galaxy-core tracks and final RMS per-particle position differences. The error metrics do not establish a converged merger by themselves.
6. **Replay saved run**: choose an NPZ trajectory and change camera bounds, point count and halo visibility without recomputing forces.

The 3D viewer supports Play/Pause, timeline scrubbing, individual-frame stepping and camera rotation. Cancellation is cooperative and is checked between full force evaluations; compilation and a single force evaluation cannot be interrupted instantly. Partial snapshots are retained when possible.

The **Apply Barnes (2004) softening and timestep** button sets `epsilon = 0.0125` and `dt = 1/512`. It does **not** automatically make all other settings identical to the published experiment.

## Reproducibility and results

Each GUI run creates `results/runs/<name_timestamp>/` containing:

- `initial_particles.npz`: exact initial state, including a reference for fair solver comparisons.
- `run_trajectory.npz`: all saved particle positions, velocities, component IDs, progenitor IDs, and settings.
- `config.json`: physical and numerical provenance, solver and generator labels.
- `viewer_3d.html`: self-contained interactive 3D playback.
- `isolation_metrics.csv`, `isolation_report.json`, `isolation_diagnostics.png` for isolated runs; plus `initial_profiles.png` / `.json` when the generator configuration is available.
- `encounter_metrics.csv`, `encounter_report.json` for collision runs.
- Optional MP4 and GIF exports.

AGAMA-generated particle realizations are also cached to `results/galaxies/`; choose **Saved snapshot** to reuse one. Saved run folders are intentionally separate from the generated galaxy library. For encounters, non-preview sources currently reuse the *same sampled realization* for both progenitors (with different orientations and orbital states); this avoids rerunning a costly DF builder but creates correlated finite-N noise. The default native preview generates two different seeds.

The results are small enough for teaching experiments but can grow rapidly for high particle counts and many saved frames. Large runs should use the compiled force solvers.

## Important scientific boundaries

- The *native disk–bulge–halo* galaxy generator samples motivated **positions** but assigns approximate, not DF-equilibrium, velocities. It is for demos only.
- The *Plummer* generator samples an analytic isotropic equilibrium DF, but finite-particle simulations with softened forces still need resolution and stability testing. Our short demonstration is only a screen.
- The *AGAMA* adapter follows the official three-component self-consistent DF workflow, including component DFs, iterative potential updates, and DF particle sampling. **AGAMA is not installed here; runtime compatibility, gravity units, model masses, and isolated stability remain unverified**. We intentionally do not label those results validated.
- The observed Mice (NGC 4676) contain gas and stellar populations absent from our collisionless model; a qualitative tidal picture cannot establish observational agreement. Actual orbital/pericenter measurements require tracking live galaxy centers, not only the nominal point-mass initialization.
- The energy diagnostic is the total *softened* N-body energy sampled during the run (up to roughly a dozen stored frames for up to 6000 particles). Drift between sampled frames is not measured; for higher N, the expensive pairwise energy diagnostic is skipped.
- **Do not treat the GUI's isolated screening thresholds (5% half-mass drift, 10% disk thickness drift, 1% endpoint energy) as universal scientific standards.** Longer isolation runs, several particle resolutions, multiple seeds, and comparisons against the published Mice are required.

## Layout

- `nbody/gravity/`: NumPy and Numba direct force, Python and Numba Barnes–Hut octree.
- `nbody/galaxies/`: analytic profiles, seeded native samplers, AGAMA DF adapter.
- `nbody/particles.py`, `integrate.py`, `encounters.py`: shared particle state, leapfrog integration, live encounters.
- `nbody/diagnostics.py`: isolated structure, energy, analytic profile comparisons, and encounter summaries.
- `nbody/solver_study.py`: end-to-end matched-initial-conditions force solver comparison.
- `nbody/io.py`: saved particle snapshots and trajectories.
- `nbody/viewer.py` / `animation.py`: interactive 3D and optional 2D video.
- `experiments/`: reproducible standalone scripts; `app.py` is the Gradio GUI.
- `docs/`: modeling choices, scientific status, usage, and cumulative bibliography.

## Scientific references

Barnes, J. E. (2004), *Shock-induced star formation in a model of the Mice*, MNRAS 350, 798, https://doi.org/10.1111/j.1365-2966.2004.07725.x ; Vasiliev, E. (2019), *AGAMA: Action-based galaxy modelling architecture*, https://doi.org/10.1093/mnras/sty2672 ; AGAMA official worked example, https://github.com/GalacticDynamics-Oxford/Agama/blob/master/py/example_self_consistent_model3.py . See `docs/BIBLIOGRAPHY.md` for the running inventory.
