# N-Baddies browser UI — user guide

## Start without a terminal (Windows)

1. Extract the entire zip archive to a writable folder.
2. Double-click **start_gui.bat** in the extracted folder. The first launch creates a `.venv` Python environment and downloads dependencies. Use Python **3.11 or 3.12** on Windows. You need an Internet connection for that first setup.
3. The launcher will print the local URL and try to open your browser. If it does not open, go to http://127.0.0.1:7860 yourself while the window is running. If startup fails, check startup_error.log. Your browser opens the local interface. It runs on your own computer at `http://127.0.0.1:7860` (or a nearby available port).
4. Keep the small launcher console open while using the GUI; closing it stops the local server.

The GUI does not upload your galaxy particles to a cloud service. The Gradio server binds to localhost (`127.0.0.1`) and `share=False`.

If your browser does not open automatically, copy the local URL printed in the launcher window into your browser.

For macOS/Linux or users comfortable with Python:

```bash
python -m pip install -r requirements.txt -r requirements-fast.txt -r requirements-gui.txt
python app.py
```

## Run an encounter

1. In **Run an encounter**, choose a run name (e.g., `tail_test_01`).
2. Start with **Native preview**, **Fast direct (Numba)**, and the default galaxy profile. The original v1.1 code handles the force evaluation and integration.
3. Choose the counts **per galaxy**: `disk=800`, `bulge=100`, `halo=350` means **2,500 live particles** in the two-galaxy simulation.
4. Change the encounter separation, nominal pericenter, and disk angles in **Encounter geometry**.
5. Choose number of steps, dt, snapshot spacing, view, FPS, and GIF/MP4 format.
6. Click **Run simulation and animate**. The GIF/video, a plot of sampled bulge-centroid separation, a trajectory (`.npz`) and run settings (`.json`) appear on the right.

The first Numba run includes a compilation delay. Barnes–Hut is currently implemented with a pure-Python octree and can be slower than accelerated direct summation at a few thousand particles. MP4 export requires FFmpeg visible to Matplotlib. GIF export works using Pillow.

## Replay an existing simulation

1. Click **Replay a saved run**.
2. Select a trajectory and click **Render saved trajectory**. The included v1.1 `mice_animation_trajectory.npz` is already in the list.
3. Change the view (`xy`, `xz`, or both), zoom, frame rate, or whether halo particles are drawn.
4. Render again. **The gravitational simulation is not rerun.** Only the saved frames are plotted.
5. After running a new simulation, click **Refresh list** to find the new trajectory.

Every new run writes a unique subfolder under `results/runs/` containing the exact saved particles and a JSON settings record. The history can be replayed or analyzed in Python later. The existing terminal experiment scripts still work.

## Physics notes

- The native generator's disk/bulge/halo **spatial distributions** are Barnes-inspired. Its **velocities are approximate** and have not passed isolated equilibrium tests. The animations are development demos, not validated astronomical reconstructions.
- The angle inputs use the code's `rotation = Rz(node) @ Rx(inclination)` convention, **not** published sky-viewing angles.
- `pericenter` defines a nominal *point-mass* initial orbit. It is not the guaranteed actual closest approach for live galaxies.
- Changing particle count, dt, softening, or the Barnes–Hut opening angle changes numerical accuracy. Increasing particles improves visual sampling but does not guarantee convergence.
- The center-separation plot uses **sampled bulge centroids** and is noisy when bulge particle counts are small.
- The AGAMA backend is marked experimental and requires separate installation and validation. Its seeding behavior is **not verified**. Save its generated `.npz` if you need exact replay.

## Technical layout

- `app.py`: UI controls and orchestration. Creates a configuration dictionary and calls the existing library.
- `nbody/galaxies/`: galaxy initial conditions.
- `nbody/encounters.py`: rotations, center offsets, and initial bulk orbital velocities.
- `nbody/gravity/`: interchangeable force solvers (NumPy direct, Numba direct, Barnes–Hut).
- `nbody/integrate.py`: leapfrog integration; yields selected frames.
- `nbody/animation.py`: GIF/MP4 from already computed frames.
- `nbody/io.py`: saved trajectories and metadata.
- `results/runs/`: user-generated trajectories and animations.

When constructing a new experiment, keep physical models in `nbody/` and use `experiments/` for numerical/astrophysical analysis. The GUI is only one interface to the same code.
