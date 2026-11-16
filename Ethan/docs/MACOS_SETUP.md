# Running N-Baddies V1.5.1 on macOS

The **same Python simulation engine, GUI, 3D viewer, presets and analyses** are used on Windows and macOS. The difference is the launcher and optional compiled AGAMA installation.

## 1. Start the GUI (AGAMA not required)

1. Install **Python 3.12** from https://www.python.org/downloads/macos/ . Python 3.11 also works. Both **Apple Silicon (M-series)** and **Intel** Macs are supported by the general Python stack, subject to compatible packages for your chosen Python and architecture.
2. Unzip the release and keep the whole `nbody_v1_5` folder together.
3. In Finder, double-click `start_gui.command` in that folder. The first run creates `.venv` and installs NumPy, Numba, Gradio, Plotly, Matplotlib and the other required Python packages. Allow an internet connection for this first installation.
4. The app should open in the default browser. If it doesn't, visit **http://127.0.0.1:7860** while its Terminal window is open. Closing that window stops the local server.

If macOS refuses to open a downloaded `.command` file, open **Terminal**, navigate to the extracted folder, and run:

```bash
cd /path/to/nbody_v1_5
bash start_gui.sh
```

This does not require changing macOS security settings. If Finder reports the file is not executable, run `chmod +x start_gui.command start_gui.sh` once (or use the Terminal command above).

The GUI only listens on `127.0.0.1` and does not publish a server to the internet. Running the simulator does not require AGAMA.

## 2. Verify installation

In Terminal from the extracted project directory:

```bash
.venv/bin/python -m pytest -q
.venv/bin/python -m experiments.two_body
.venv/bin/python -m experiments.isolation_study --model plummer --particles 500
```

The source-code test suite currently passes on the development environment. **A real macOS device has not been available for runtime verification.** If you see an error, copy the Terminal traceback or `startup_error.log`. Do not repeatedly reinstall Python unless the message identifies a missing or incompatible Python runtime.

## 3. Optional: install AGAMA for disk equilibrium initial conditions

AGAMA is a **compiled C++/Python package**, and its source archive is in `vendor/Agama-master.zip`. The source archive alone does not make `import agama` work.

1. Run the GUI launcher once, creating `.venv`, then close the application.
2. Install Apple's developer command line tools with `xcode-select --install` (if not already installed).
3. Install [Homebrew](https://brew.sh/) if you don't have it, then run:

```bash
brew install gsl eigen libomp
```

4. From the project root:

```bash
bash scripts/install_agama_macos.sh
.venv/bin/python -m experiments.agama_selftest
```

The installer runs pip **inside this project's virtual environment**, adds Homebrew's GSL/Eigen/libomp locations to compiler flags, and verifies that `import agama` succeeds. It does **not** automatically install Homebrew, Xcode tools, or those libraries. The compiler may prompt for additional optional dependencies. See [AGAMA upstream INSTALL](https://github.com/GalacticDynamics-Oxford/Agama/blob/master/INSTALL).

**Apple Silicon note:** Use arm64 Python on arm64 hardware, and match it to your Homebrew libraries under `/opt/homebrew`. A Rosetta x86_64 Python paired with arm64 libraries will cause linker/import errors. On an Intel Mac, Homebrew typically lives at `/usr/local`.

AGAMA has been integrated at the API and mock-test level, but a compiled AGAMA run and isolated-disk validation remain outstanding. A successful import is not scientific validation.

## 4. Carry over existing runs

You can move the old `results/runs/` and `results/galaxies/` folders from Windows to Mac. The snapshots and trajectories are NumPy `.npz` files, not compiled executables. Keep the scientific configuration/metadata with each run.

**Don't copy `.venv` between Windows and Mac.** It includes OS-specific executables and compiled extensions. Each machine should create its own `.venv`.

## 5. Performance and known limits

- Numba compiles force kernels for the local architecture on first use, so the first force evaluation may be slower. Subsequent calls reuse compiled code/caches.
- The compiled direct solver still scales as O(N²); the tree code uses an approximate force and has a different runtime/error tradeoff.
- A laptop can heat up substantially during high-particle-count runs. Start at modest N and stop long calculations with the GUI's Cancel button (checked between force evaluations).
- The native disk generator remains a *preview*, not a verified equilibrium model. Select Plummer for analytic DF tests or install AGAMA and validate the generated disk in isolation.
