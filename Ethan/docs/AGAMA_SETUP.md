> **Windows update (V1.6):** use the direct PowerShell instructions in [WINDOWS_AGAMA_POWERSHELL.md](WINDOWS_AGAMA_POWERSHELL.md). The new `scripts/install_agama_windows.ps1` locates Visual Studio and runs the bundled source build in the GUI's `.venv`. Build and real-model runtime remain unverified on Windows.

# AGAMA — installing the supplied master source, and what "implemented" means

N-Baddies V1.5 ships with `vendor/Agama-master.zip`, copied from the supplied AGAMA-master archive.
The official AGAMA source is https://github.com/GalacticDynamics-Oxford/Agama .
Documentation in that archive: `INSTALL`, `py/example_self_consistent_model3.py`,
`py/example_self_consistent_model_bdh.py`, and `doc/reference.tex`.

## Source versus importable extension

AGAMA is a Python interface to a C++ library. Merely having `Agama-master.zip` on disk does **not**
make `import agama` work; C++ code must be compiled to a `.pyd` on Windows or `.so`/`.dylib`-linked extension on macOS/Linux.

**Test environment:** The V1.5 builder had Python, NumPy, SciPy, Numba and Gradio, but no AGAMA
extension and no GSL *development headers*. Accordingly, we could not execute the real AGAMA
physics backend. The included mock integration tests only check the calls that our adapter makes.
Do not infer that AGAMA-generated galaxies have already been checked for equilibrium.

## macOS installation (Intel and Apple Silicon)

The same AGAMA archive can be compiled on Mac. See **`docs/MACOS_SETUP.md`** for complete launcher instructions. The optional Mac build steps are:

```bash
xcode-select --install  # only if the Apple command-line tools are missing
brew install gsl eigen libomp
bash scripts/install_agama_macos.sh
.venv/bin/python -m experiments.agama_selftest
```

AGAMA's own `setup.py` detects macOS and specifies the host architecture. Use an arm64 Python with arm64 dependencies on Apple Silicon; do not mix Intel/Rosetta and ARM libraries. The script passes Homebrew include/library paths through `CFLAGS`/`LDFLAGS`, as permitted by the supplied AGAMA setup.py. No real macOS build has been exercised by this release. An imported AGAMA extension still requires model mass, kinematic, and isolated stability tests before any physical validation claim.

## Windows installation (x64 Python 3.11 or 3.12)

1. Install **Visual Studio 2022 Build Tools** with the Desktop development with C++ workload.
2. Open **x64 Native Tools Command Prompt for VS 2022** (sets up `cl.exe`, `nmake`, etc.).
3. In that prompt, `cd` to the extracted N-Baddies V1.5 root.
4. Run `start_gui.bat` at least once to make `.venv` (you can close the GUI afterwards).
5. Run `scripts\install_agama_windows.bat` from the *developer prompt*.
6. If GSL is missing, AGAMA's build system may prompt to obtain/build the GSL port. This step
   may require internet access and may fail if the upstream mirror is unavailable. Read errors.
7. Verify the same interpreter used by our GUI:

```bat
.venv\Scripts\python.exe -c "import agama; print(agama.__file__, agama.G)"
.venv\Scripts\python.exe -m experiments.agama_selftest
```

8. Refresh the GUI's **Check AGAMA availability**, or restart the app.

The source installer **does not** silently substitute our approximate disk generator on failure.
It also does not overwrite saved galaxies or experimental runs.

**Alternative:** From the same developer prompt and in the V1.5 root, run:

```bat
.venv\Scripts\python.exe -m pip install -U pip setuptools wheel numpy
.venv\Scripts\python.exe -m pip install -v --no-build-isolation vendor\Agama-master.zip
```

The `--no-build-isolation` flag follows AGAMA's instructions for finding the current NumPy.
The source archive is bundled, but any missing C++ dependencies can still require downloads.
Inspect the bundled upstream `setup.py` and `LICENSE` if desired.

**Troubleshooting**

- `ModuleNotFoundError: No module named 'agama'`: the extension was not installed *in `.venv`*.
- `cl is not recognized`: use the Visual Studio **x64 Native Tools** developer prompt.
- C++ compile/link error mentioning GSL: GSL headers or libraries are missing or wrong-architecture.
- A C++ compilation succeeds but `import agama` fails: confirm 64-bit architecture, the Python
  version used for building, and the presence of its MSVC runtime dependencies.
- GUI says AGAMA unavailable after installation: check its displayed `sys.executable`, then
  run the two `-c "import agama"` tests using that **same** executable.

Consult the precise upstream build instructions:
https://github.com/GalacticDynamics-Oxford/Agama/blob/master/INSTALL

## What the adapter actually does

`nbody/galaxies/sources/agama.py`:

1. Builds a negative-scale-height AGAMA `Disk` density plus Hernquist bulge and cored Dehnen halo.
2. Constructs an initial potential from these density components, with spherical and cylindrical
   grids extending to the halo's radius range (much larger than the old V1.4 grid).
3. Uses `Q(Rd)` to estimate a radial velocity dispersion and builds a quasi-isothermal disk DF.
4. Corrects its normalization using `df_disk.totalMass()`.
5. Builds `QuasiSpherical` bulge and halo DFs from the initial potential and component densities.
6. Iterates the self-consistent model's potential and density, using AGAMA's own solver.
7. Samples full 6D phase-space particles from each DF; when making a binary, samples **two independent
   realizations** from the one model using different seeds, instead of duplicating the same particles.
8. Recenters positions and mean velocities and saves `agama_generation_report.json` and cached
   `.npz` initial conditions.

**Physics caution:** source `config['halo']['cutoff']` is only a native-preview sampling cutoff.
The AGAMA profile in this release is a finite-mass `gamma=0, beta=4` Dehnen halo **without a hard
cutoff**. The spherical grid reaches beyond the reference cutoff but does not change the halo law.
The final DF-based density can differ from the initial analytic profiles. This is why the report
records both requested masses and actual DF/sample masses instead of silently forcing a match.

**Units:** All modeled coordinates are dimensionless `G=1` code units. **Do not call**
`agama.setUnits(mass=...,length=...,velocity=...)` with arbitrary physical scales and then
pass the same dimensionless parameter values. AGAMA explicitly documents that `setUnits` changes
numerical `G`. Converting to physical kpc, km/s, solar masses requires a consistent scale transform.

## Scientific verification after first successful import

1. Run `python -m experiments.agama_selftest` to test actual construction and sampling.
2. Open `results/galaxies/agama_smoke_900.json` and compare requested, DF, and sampled masses.
3. In the GUI: **Test a saved equilibrium galaxy** → choose `agama_smoke_900.npz` → Run with
   `dt` and `dt/2`. This is only a **numerical screen**.
4. Inspect `r50`, disk `z_rms`, `sigma_R`, `energy` and disk Fourier amplitudes.
5. Rerun with more particles and independent seeds; vary force softening, solver and time step.
6. Only then permit using a configuration for quantitative Mice claims.

A self-consistent DF in AGAMA's smooth continuum potential does **not** guarantee equilibrium
in our finite-N, Plummer-softened N-body realization.
