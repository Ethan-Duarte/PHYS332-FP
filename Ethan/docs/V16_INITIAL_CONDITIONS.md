# V1.6 — Global virialization, unconstrained relaxation, and AGAMA

## Three distinct experiments

1. **Native disk–bulge–halo**: sample specified positions and approximate velocities. Its isolated disk previously heated and thickened dramatically; treat this as a deliberately unvalidated control.
2. **Globally virialized native**: multiply all internal velocities by the common factor needed to reach a specified global virial ratio. This preserves particle positions and the bulk velocity, but *changes disk rotation, angular momentum, radial and vertical dispersions*. It does not solve the collisionless Boltzmann equation. Compare this candidate with (1), never silently replace the original reference galaxy.
3. **Unconstrained relaxed native**: evolve a saved single galaxy in isolation, save its late snapshot, reset simulation time to zero, and evolve again. Look for structural drifts in both phases. The late galaxy may be thicker, hotter or less concentrated than the intended original density profile.
4. **AGAMA-generated candidate (when importable)**: independently generate positions and velocities from a self-consistent potential/distribution-function model using `experiments.agama_selftest`, then apply exactly the same isolated tests. The AGAMA source is bundled, but the compiled extension is a separate installation.

## Why the virial ratio is not the energy ratio

In the center-of-mass frame, define the system moment of inertia

\[I=\sum_i m_i|\mathbf r_i-\mathbf R_{\rm CM}|^2.\]

Differentiating twice under Newton's second law gives

\[\tfrac12\ddot I=2K+W,\quad
K=\tfrac12\sum_i m_i|\mathbf v_i-\mathbf V_{\rm CM}|^2,\quad
W=\sum_i m_i(\mathbf r_i-\mathbf R_{\rm CM})\cdot\mathbf a_i.\]

For a stationary isolated system (or a suitable long-term average), the right-hand side approaches zero. We measure `Q_vir = 2K / |W|` and make a *global* velocity-rescaling candidate with

\[\mathbf v'_i=\mathbf V_{\rm CM}+\sqrt{Q_{\rm target}/Q_{\rm measured}}\,(\mathbf v_i-\mathbf V_{\rm CM}).\]

With Plummer-softened pair forces, the virial term **is not equal to the potential energy**. The code therefore computes `W` from particle positions and actual accelerations. `nbody/virial.py` contains these two functions.

Even Q=1 at the initial instant does **not** imply an equilibrium phase-space distribution. Stars can still be on unsuitable individual orbits; a disk can still develop bars or heat.

## GUI walkthrough

- Start with **Run a simulation → Isolated galaxy → Disk–bulge–halo → Native preview**. Save a run with a manageable number of particles.
- Open **Virialize & relax galaxies**, click **Refresh sources**. Choose **Saved isolated trajectory**, then select your isolated run.
- Select **First frame**. Set softening to the value used in the run; choose `Fast direct` and target global virial ratio `1.0`.
- Click **Make virialized candidate**. A new snapshot appears under `results/galaxies/`. You can select it in **Test a saved equilibrium galaxy** after clicking Refresh.
- Alternatively, click **Relax, save and independently restart**. Defaults: `dt=0.004`, 100 relaxation steps, 100 restart steps. Longer runs are necessary before making stability claims.
- Compare `results/runs/initial_conditions_.../relax_restart.png` and `initial_condition_report.json`: structural profiles and the virial ratio before and after relaxation. Check that the original target density profile has not drifted unacceptably.
- Enable **Globally virialize velocities before relaxation** to combine approaches, but treat that as a **separate fourth candidate**.
- The button **Cancel relaxation** requests cancellation between timesteps; it cannot interrupt an in-progress compiled force evaluation. Partial relaxation trajectories are saved and not promoted to galaxy candidates.

For a fully scripted demonstration, run `python -m experiments.initial_condition_study` from the project root. It uses the **same initial particle realization** for native and globally virialized candidates.

## Output files

`results/runs/initial_conditions_.../` contains `prepared_initial.npz`, `relaxation_trajectory.npz`, `relaxed_candidate.npz`, `restart_trajectory.npz`, both isolation reports, `relax_restart.png`, and `initial_condition_report.json`. Completed candidates are also copied into `results/galaxies/` for direct reuse in the GUI.

If the source is an existing trajectory, you may choose its first or final frame. Only *isolated* trajectories are accepted. Merger snapshots aren't valid initial single-galaxy candidates.

## Presenting this study scientifically

For each method, plot half-mass radius, vertical disk thickness, radial velocity dispersion and global virial ratio. Prefer matching particle counts, seeds, softening, integration times and force solvers, and vary resolution later. A longer stationary-looking interval is necessary but not sufficient evidence of equilibrium. Explain why unconstrained relaxation can change the target galaxy and why AGAMA addresses a different problem by constructing self-consistent DFs.

## Decision register

- Keep G=1 code units, the native preview as a diagnostic control, and AGAMA as an unverified candidate until built and isolation-tested.
- Use the *force-consistent* virial term W, not the Newtonian unsoftened potential energy as a shortcut.
- Preserve every pre- and post-relaxation particle snapshot with source settings. A restarted candidate is not relabeled as the original Barnes galaxy if its parameters have drifted.
- A short integration timestep-refinement check addresses numerical accuracy, not the full equilibrium question.

## Bibliography additions / retained baseline

- Rodionov & Sotnikova (2006), iterative N-body galaxy modeling, https://arxiv.org/abs/astro-ph/0611176
- Yurin & Springel (2014), GALIC initial conditions, https://doi.org/10.1093/mnras/stu1421
- Vasiliev (2019), AGAMA and self-consistent modeling, https://doi.org/10.1093/mnras/sty2672
- Hernquist (1993), composite galaxy realizations, https://doi.org/10.1086/191784
- Barnes (2004), Mice numerical encounter, https://doi.org/10.1111/j.1365-2966.2004.07725.x

See `docs/BIBLIOGRAPHY.md` for the running reference register.
