# N-Baddies V1.3 GUI guide

## Launch (Windows)

Keep your existing virtual environment if upgrading from V1.2: copy the new files over, or extract the new zip into a fresh directory. The launcher `start_gui.bat` checks and installs the new Plotly dependency, then opens `http://127.0.0.1:7860`.

**Suggested first run:** Leave the Mice-inspired parabolic preset, select **Fast direct (Numba)**, click Run simulation. The main output is a **3D Plotly viewer**: rotate by dragging, zoom using the wheel, play/pause or move the time slider. Halo particles are simulated regardless of whether they're drawn.

## Custom initial conditions

Select **Custom relative vectors** and supply:

- `Δx, Δy, Δz`: initial position of galaxy B's center relative to galaxy A.
- `Δvx, Δvy, Δvz`: initial velocity of galaxy B relative to galaxy A.
- Inclination and node: internal disk orientation for each galaxy.

The two galaxies are offset into the center-of-mass frame. The custom vectors are *initial state vectors*, not orbital elements; a nominal pericenter is not enforced. All quantities are in the simulation's G=1 code units.

## Cancel and replay

Use **Cancel run** while a calculation is ongoing. The application checks cancellation between timesteps and (for compiled Barnes–Hut) after each force evaluation. Compilation on the first use cannot be interrupted immediately. The last completed snapshots are saved as a partial trajectory where available. A cancelled run skips optional GIF/MP4 rendering.

In **Replay a saved run**, you can change zoom, particle display cap, and halo visibility without repeating the gravitational integration. Download the self-contained **viewer_3d.html** file to inspect the encounter offline. **3D play controls animate the saved particle states; they do not re-integrate the system.**

## Force comparisons

In the **Compare gravitational solvers** tab, compare a single particle distribution using a shared Plummer-softened force law. The reference is compiled O(N²) direct summation. The benchmark reports runtime, median/95th percentile relative acceleration error, and maximum absolute acceleration error. It excludes compilation time. The NumPy version is only included at feasible particle counts. This is *force verification*, not proof of a physically validated galactic merger.

## Known limitations

1. Mice initial positions and masses use Barnes-inspired parameters, but the native velocities are **heuristic** and the initial disk–bulge–halo is **not in proven collisionless equilibrium**. We have not reproduced the 2004 simulation's observed morphology.
2. The AGAMA adapter is still an optional, unverified backend and is deliberately not exposed in this GUI yet.
3. The 3D Plotly playback may slow at very high displayed particle counts; reduce **Maximum displayed particles** while retaining the full-resolution gravitational evolution and saved trajectory.
4. The measured bulge-centroid separation is noisy when each bulge contains only a small number of particles; this is not a robust dynamical-core tracker.
5. Barnes–Hut forces are not pairwise symmetric, so momentum conservation can differ from direct summation; quantify it in convergence studies.

## Next scientific milestones

1. Fix and test an AGAMA-generated (or equivalent) equilibrium disk–bulge–halo model in isolation.
2. Reproduce a Barnes-style Mice encounter with a documented orbital orientation and viewing direction, then compare observed tidal-tail geometry.
3. Extend the force comparison to full-trajectory convergence of tail observables, energy, momentum, timesteps, particle count, and opening angle.
