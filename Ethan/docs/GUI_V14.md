# V1.4 Guide: one trusted galaxy before many galaxy types

## A first experiment

Open the GUI, choose **Isolated galaxy**, set **Plummer (analytic DF)**, 500–1000 total particles (the counts in the three particle sliders are added together), a timestep of 0.002 and a softened force of 0.0125, then run. The isolated model starts from a known analytical phase-space distribution rather than guessed velocities. Watch the 3D timeline and the graph showing radius changes. Inspect `isolation_report.json` and `isolation_diagnostics.png`.

Try a longer run. An energy change small at a very short duration cannot establish stability for several orbital periods, and a radius change may arise from finite sampling and softening. Repeat with different seeds and resolutions.

## Galaxy building

Choose **Disk–bulge–halo**. The **Native preview** generator creates correct intended spatial distributions but not equilibrium velocities. Use that for GUI experimentation, not for astrophysical claims.

Choose **AGAMA** only if installed and compiled for the Python version on your machine. A successful sample is an **equilibrium candidate** until tested. The GUI saves the sampled galaxy to `results/galaxies/`. It also exposes a target Toomre Q parameter around one disk scale; actual Q(R) and velocity moments must still be measured from the DF model.

Choose **Saved snapshot** to reuse an existing galaxy without regeneration. Snapshot provenance remains in its `.npz`; if not accompanied by a compatible source configuration, theoretical density checks are skipped rather than inferred.

## The Mice model

Use **Mice encounter**. The starting masses and scales are Barnes-inspired and the nominal parabolic pericenter is 0.25. The published inclination/node pairs are (25°, −30°) and (40°, 60°). The program's initial COM separation is a development choice, not an observable. Use the **Barnes numerical settings** button for `eps = .0125` and `dt = 1/512`, but do not assume observational matching.

Run the same initial condition at increasingly fine timesteps and with different force solvers. Save every run. The diagnostics track sampled bulge separation and endpoint softened energy drift; they do not identify tidal tails or fit observed sky images.

## Solver study

Use **Compare gravitational solvers** for a fast force comparison. The `p95` is the 95th percentile of relative acceleration error against fast direct summation.

Use **Full-trajectory solver study** to integrate the same initial particles with several solvers. The study writes complete NPZ trajectories, a plot, a JSON report, energy drift and RMS final position discrepancy. A few thousand-body accuracy tests can be expensive; start with 200–400 total particles and 20–30 steps.

## Interpretation

A *screening pass* says that chosen numerical/structural thresholds were met over one run. **A validated galaxy** additionally requires phase-space equilibrium diagnostics, longer isolated stability, repeated seeds, spatial/temporal resolution convergence and an observationally appropriate comparison. A *visually appealing interaction* alone does not count.
