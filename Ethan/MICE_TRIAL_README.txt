Extract into Ethan, after the Equilibrium Diagnostics patch. This adds a trial runner and a second independently sampled 12k corrected progenitor. It does not overwrite the original GUI or cached models.

PREPARATION ONLY (no evolution):
.\.venv\Scripts\python.exe .\mice_candidate_run.py encounter --a results\softened_review\baseline\candidate.npz --b results\softened_review\baseline_seed43\candidate.npz --prepare-only --out results\mice_prepare

SHORT SMOKE TEST (t=0.125; not enough for the Mice tails):
.\.venv\Scripts\python.exe .\mice_candidate_run.py encounter --a results\softened_review\baseline\candidate.npz --b results\softened_review\baseline_seed43\candidate.npz --steps 128 --out results\mice_smoke

Each output directory must be new. The full-force timing forecasts runtime before integration. Ctrl+C retains completed checkpoints.

This trial uses the corrected LOW-RESOLUTION progenitors, whose isolated thickening is known. It is for encounter plumbing and exploratory geometry; it is not a validated Mice fit. Default dt=1/1024 and epsilon comes from both snapshots (0.00625). Do not change epsilon without regenerating the initial conditions.

Nominal separation=1 and pericenter=0.25; disk inclinations 25/40 degrees and literature-relative nodes -30/+60, converted to the point-mass reference periapsis direction. These orientation conventions and the live closest approach must be checked before a published-model reproduction is claimed. A point-mass nominal pericenter is not the measured live pericenter.

The runner exports initial_particles.npz, config.json, run_trajectory.npz and, on completion, core_separation.csv. The center diagnostic is the coordinate-wise median bulge position, not a bound-core identification algorithm. It uses stellar disk centering to avoid remote halo COM outliers; the associated changes in the relative component frame are recorded in source metadata.

Future longer trials can use --steps 1536 (t=1.5), but inspect the forecast first. Time zero is the initial configuration, NOT first pericenter. Identify measured first passage and compare snapshots by time since that event. Do not assume the final frame is the observed Mice epoch. Rendering the trial does not establish agreement with observational morphology or kinematics.

The tested preparation used the actual 12k seed42/43 snapshots. A full collision was not reproduced in this delivery. Use the higher-resolution candidates only after isolation checks, and generate an independent matching seed43 with resample_softened_candidate.py. No gas or star formation is included.

Reference: Barnes 2004 https://arxiv.org/abs/astro-ph/0402248
