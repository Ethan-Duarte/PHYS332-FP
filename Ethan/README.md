# N-Baddies: N-body Galaxy Simulator — V1.3

A small, fresh Python project. The direct force calculation, leapfrog, and octree descend from our earlier student implementations, but the file layout is new.

## Browser interface (V1.3) — interactive 3D + cancellation

**Extract the entire zip** to a normal folder and double-click `start_gui.bat`.
The corrected launcher creates `.venv` if needed, installs missing dependencies, and
starts the local Gradio server. **The window stays open and shows errors** if startup fails.

If the browser doesn't open, manually visit **http://127.0.0.1:7860** while the
launcher is still running. If Python fails during startup, its full traceback is
written to `startup_error.log` in the project folder. Use that file to diagnose
what's missing. Avoid launching `app.py` by double-clicking it.

If you installed dependencies earlier, keep the existing `.venv` folder when
replacing the project; they shouldn't need to be installed again. If dependency
errors persist, delete `.venv` once and rerun the launcher with Python 3.11 or 3.12.

**Run an encounter** supports the Mice-inspired parabolic setup or custom relative position/velocity vectors; the main result is a draggable, time-animated 3D viewer. **Replay a saved run** changes camera bounds and particle visibility without rerunning gravity. Results are saved under `results/runs/`.

For manual startup in a Windows Command Prompt opened in this folder:

```bat
.venv\Scripts\python.exe -u launch_gui.py
```

See `docs/GUI_V13.md` for the current interface and `docs/GUI_GUIDE.md` for the earlier interface. The interface runs only on
localhost. It does not make the physical model scientifically validated.

## Status (important)

- **Working/tested:** softened direct gravity, Barnes–Hut tree gravity, leapfrog, particle structure, saved snapshots, exact isotropic Plummer sphere, analytic density profiles, illustrative disk–bulge–halo sampling, simple parabolic encounter, basic diagnostics.
- **Illustrative only:** `mice_preview()` generates the published-style *spatial* components, but uses heuristic velocities. Its galaxy is **not** a validated collisionless equilibrium; resulting tidal shapes cannot yet be presented as a Mice reconstruction.
- **Unverified optional backend:** `galaxies/sources/agama.py` follows AGAMA's documented self-consistent DF process, but AGAMA is not installed in the test environment. Its API behavior, equilibrium quality, sampling reproducibility, and unit consistency require testing on a local AGAMA installation.
- **Still needed for real validation:** isolated galaxy stability for disk–bulge–halo, adequate resolution, velocity profiles, convergence, physically matched orbital/viewing parameters, observational comparisons.

## Quick start

Run from this folder:

```bash
pip install -r requirements.txt
python -m pytest -q
python -m experiments.two_body
python -m experiments.plummer_validation
python -m experiments.mice_preview
python -m experiments.force_benchmark
```

Mice preview writes `results/mice_preview.png` and `results/mice_preview.npz`.

Optional AGAMA is independently installed (`pip install agama`, which compiles a C++ extension and may require platform-specific dependencies); use:

```python
import json
from nbody.galaxies.generate import galaxy
cfg = json.load(open('presets/mice.json'))
model = galaxy(cfg, backend='agama')
```

**AGAMA seed warning:** the `seed` parameter governs native NumPy sampling only. AGAMA's internal sampler isn't established here to use that seed; save the resulting `.npz` for exact replay.

## Directory map

- `nbody/particles.py`: common particle arrays `(x,v,m,component,galaxy)`, combine and recenter.
- `nbody/gravity/`: direct and octree-based Barnes–Hut acceleration.
- `nbody/integrate.py`: kick–drift–kick timestep and periodic snapshots.
- `nbody/galaxies/profiles.py`: Hernquist/Dehnen and exponential disk analytic profiles.
- `nbody/galaxies/sources/native.py`: Plummer exact equilibrium and approximate Mice preview.
- `nbody/galaxies/sources/agama.py`: separate, optional multicomponent equilibrium backend.
- `nbody/galaxies/generate.py`: generation entry point.
- `nbody/encounters.py`: orientations, offsets, Keplerian parabolic initial orbit.
- `nbody/analysis.py`: energy, momentum, radius and core separation.
- `nbody/io.py`: portable snapshot writing/loading.
- `presets/mice.json`: Barnes-inspired mass and spatial scales, developer-level particle counts.
- `experiments/`: independently executable scripts.
- `tests/`: automated mathematical and plumbing checks.

## Units and model details

Internal gravitational constant is **G=1**. Barnes (2004) provides a physical mapping for his reference Mice fit: 1 length unit ≈35.7 kpc, 1 speed unit ≈180 km/s, 1 time unit ≈194 Myr, 1 mass unit ≈2.69×10¹¹ solar masses. These are *model-matching scales*, not fundamental conversion constants.

The `mice.json` reference uses Barnes's idealized collisionless density shapes: stellar disk ∝ exp(-R/Rd) sech²(z/zd), Hernquist bulge (Dehnen gamma=1), cored Dehnen halo (gamma=0). It does **not** reproduce Barnes's SPH gas component. Truncated spatial sampling is conditionally normalized to the specified component mass, so this is not an exact untruncated Barnes realization.

`encounters.parabolic_pair` places the galaxy COMs on an initial *point-mass* parabolic Kepler orbit. The live, extended halos alter the subsequent orbit; measure the actual pericenter from the evolution. Disk `rotation` adopts the explicit convention `Rz(node_angle) @ Rx(inclination)`; it is not a published sky-viewing transformation.

**Softening**: both gravity solvers calculate acceleration proportional to `(r² + eps²)^(-3/2)` with the same `eps`. The diagnostic energy uses the matching potential `-Gmimj/sqrt(r² + eps²)`.

**Barnes–Hut performance**: this clean Python implementation is useful for learning and force comparison. It may be slower than NumPy-vectorized direct summation for modest N; it needs a separate performance/optimization study before large simulations.

See `docs/WALKTHROUGH.md` for a conceptual lesson, `docs/VALIDATION.md` for the tests required before making scientific claims, and `docs/BIBLIOGRAPHY.md` for the running references.

## References

- Barnes, J. E. (2004), *Shock-induced star formation in a model of the Mice*: https://doi.org/10.1111/j.1365-2966.2004.07725.x
- Hernquist, L. (1990), analytic spherical galaxy: https://doi.org/10.1086/168845
- Dehnen, W. (1993), spherical gamma-family profiles: https://doi.org/10.1093/mnras/265.1.250
- Barnes & Hut (1986), tree algorithm: https://doi.org/10.1038/324446a0
- Vasiliev (2019), AGAMA: https://doi.org/10.1093/mnras/sty2672
- AGAMA three-component worked example: https://github.com/GalacticDynamics-Oxford/Agama/blob/master/py/example_self_consistent_model_bdh.py

## Higher-resolution galaxy animation (V1.1)

The older `mice_preview.py` is intentionally a tiny static demonstration. The new `mice_animation.py` evolves a larger **live** particle system, saves all selected frames to a replayable `.npz`, and exports a GIF or MP4. It shows disk and bulge particles but hides halo particles by default; the halo still contributes gravity.

Install the optional compiled acceleration dependency:

```bash
pip install -r requirements.txt -r requirements-fast.txt
```

Run a 2,500-particle encounter (1,600 stellar disk, 200 bulge, 700 halo):

```bash
python -m experiments.mice_animation --format both
```

For a higher-resolution run like the included example (3,500 live particles, including 2,600 disk stars):

```bash
python -m experiments.mice_animation --disk 1300 --bulge 100 --halo 350 --steps 240 --dt 0.006 --save-every 4 --format both
```

All particle counts are **per galaxy**, so the total count is twice their sum. The optional `--solver fast` uses Numba-compiled **direct summation** (the same softened force law as `gravity/direct.py`, with O(N) temporary memory). You can select `--solver direct` or `--solver tree` to run the original two implementations, although these are significantly slower at high counts. The first Numba call includes one-time compilation overhead.

Outputs:

- `results/mice_animation_trajectory.npz`: saved times, full positions and velocities, masses, particle IDs, component IDs, and metadata
- `results/mice_animation.gif`: automatically looping animation (Pillow writer)
- `results/mice_animation.mp4`: compressed video (requires FFmpeg available to Matplotlib)

To **rerender the same simulation without recalculating gravity**, use:

```bash
python -m experiments.mice_animation --load results/mice_animation_trajectory.npz --view xy --lim 1.8 --fps 24 --format mp4
```

Useful options: `--disk`, `--bulge`, `--halo`, `--seed`, `--steps`, `--dt`, `--save-every`, `--fps`, `--view xy|xz|both`, `--lim`, `--format gif|mp4|both`, `--show-halo`, `--solver fast|direct|tree`, `--name my_run`. Use a unique `--name` for parameter comparisons so files are not overwritten.

**Scientific limitation:** More particles help visualize fine structure but don't fix the approximate velocity distributions in `native.mice_preview`. The animations remain illustrative until the disk–bulge–halo equilibrium model, timestep/softening convergence, and encounter configuration are verified and validated. A bigger N is not equivalent to a converged scientific result. See `docs/VALIDATION.md`.

---

## V1.3 — GUI, live cancellation, interactive 3D, custom encounter vectors

Run `start_gui.bat` on Windows or `python -m pip install -r requirements.txt -r requirements-fast.txt -r requirements-gui.txt -r requirements-viewer.txt` followed by `python app.py`.

- The main viewer is an **interactive WebGL 3D plot**: drag to rotate, scroll to zoom, play/pause or select a saved time. A self-contained HTML copy can be downloaded.
- **Cancel run** requests termination at the next timestep/force-check boundary. Partial completed snapshots are retained.
- **Mice-inspired** and **Custom relative vectors** modes select nominal parabolic initial conditions or arbitrary relative position/velocity vectors.
- **Barnes–Hut (Numba)** accelerates the existing octree's traversal, retaining the original Python implementation as a test reference. A Python Barnes–Hut option is not offered for long GUI runs.
- **Compare gravitational solvers** benchmarks acceleration error and speed against accelerated direct summation on one shared particle distribution.
- The GUI hides the untested AGAMA backend until its API and isolated equilibrium are validated.

See `docs/GUI_V13.md` for parameters, limitations, and the scientific roadmap. The simulator remains an **exploratory preview** of a disk-galaxy encounter, not a validated Mice reconstruction.


## V1.3.1 — 3D playback fix (2026-11-03)

The embedded Gradio Plotly control did not reliably apply animation frames, even though
the saved trajectories contained changing particle positions. The GUI now embeds an
independent offline HTML viewer with its own **Play / Pause**, step buttons,
**time slider**, and interactive 3D camera. These controls directly update Plotly's
particle coordinates, rather than depending on Gradio's Plotly frame renderer.

Use **Replay a saved run** to open existing `.npz` trajectories; no new simulation
is necessary. Every replay also saves `viewer_3d.html`, which can be opened
independently in your browser.

The embedded preview may limit playback to 180 evenly spaced snapshots for browser
responsiveness. All simulation snapshots are still preserved in the `.npz` file.

If the embedded viewer does not appear, download the **Standalone interactive 3D viewer**
file in the GUI and open the HTML directly in your browser.
