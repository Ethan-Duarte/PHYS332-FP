# N-Baddies: N-body galaxy simulations

Direct and Barnes–Hut gravity, galaxy construction, isolated controls and exploratory Mice encounters. The current priority is presentation plots and animations; the planned isolated validation sequence is complete for now.

## Start here

- **Open the GUI:** double-click `start_gui.bat` on Windows, or run `bash start_gui.sh` on macOS/Linux.
- **Current model and limitations:** [docs/CURRENT_MODEL.md](docs/CURRENT_MODEL.md).
- **CLI commands and tool index:** [docs/WORKFLOW.md](docs/WORKFLOW.md). Run commands from this `Ethan` directory.
- **Sources:** [docs/BIBLIOGRAPHY.md](docs/BIBLIOGRAPHY.md).
- **Cleanup details:** [docs/CLEANUP.md](docs/CLEANUP.md).

## Folder guide

| Folder | Contents |
|---|---|
| `nbody/` | Gravity solvers, integrator, galaxy models, diagnostics and rendering |
| `tools/` | Command-line generation, analysis, reference and encounter tools |
| `experiments/` | Existing module-based demonstrations and benchmark studies |
| `validation_suite/` | Reusable measurements and reference reports |
| `assets/initial_conditions/` | Prepared final and halo-refined states, with provenance |
| `assets/reference/` | Original reference particle export |
| `assets/galaxies/legacy_agama/` | Earlier cached candidates for reproducibility, not the selected model |
| `presets/` | Model configurations |
| `results/` | Your local trajectories, generated galaxies, plots and exports; ignored by Git |
| `docs/` | Guides, scientific notes, bibliography and historical checkpoint instructions |
| `scripts/` | AGAMA installation helpers |
| `tests/` | Software regression checks |
| `vendor/` | Supplied AGAMA source archive; generated local builds remain local |

## Environments

Keep using the existing `.venv` in this directory. The cleanup does not require rebuilding AGAMA. For a fresh environment:

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements-science.txt
```

The GUI launcher installs its GUI requirements if missing. AGAMA is optional for generators that need it; replay and Zeno reference tools do not require importing AGAMA. Existing WSL Zeno/GADGET installations are separate and are not moved.

## Tests

```powershell
.\.venv\Scripts\python.exe -m pip install -r requirements-dev.txt
.\.venv\Scripts\python.exe -m pytest tests -q
.\.venv\Scripts\python.exe -m unittest discover -s validation_suite/tests -q
```

The core solvers and physical parameters were not changed by the cleanup. Software tests do not certify a stationary galaxy or an observational fit.
