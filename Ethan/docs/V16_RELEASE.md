# N-Baddies V1.6 release notes

**Built from V1.5.1 Mac/Windows release.**

- New `nbody/virial.py`: consistent global virial measurement `Q = 2K/|W|` using the actual softened gravitational acceleration; optional scalar internal-velocity rescaling.
- New `nbody/initial_conditions.py`: unconstrained isolated relaxation, saved late snapshot, independent restart with time reset, structural and virial history report and comparative plot.
- New **Virialize & relax galaxies** GUI tab: select a galaxy snapshot or isolated trajectory, generate virialized candidate, relax/restart, save output for future comparisons, and cancel between timesteps.
- New `scripts/install_agama_windows.ps1`: find installed x64 MSVC build tools and install supplied AGAMA source into the GUI `.venv`; step-by-step guide in `docs/WINDOWS_AGAMA_POWERSHELL.md`.
- New `experiments/initial_condition_study.py`: reproducible matching seed comparison with native and globally virialized candidate.

**Verified:** 41 Python tests pass; local GUI serves HTTP 200; the 440-particle native and virialized experiments finish with saved snapshots, trajectories, and plots.

**Not verified:** AGAMA native extension compilation or self-consistent model creation on the user's Windows machine; Mac/Windows launcher behavior with these new files. No galaxy has been demonstrated as scientifically equilibrium-validated. AGAMA remains optional and we never silently substitute the native generator.
