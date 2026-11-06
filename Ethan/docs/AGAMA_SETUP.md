# AGAMA: dependency status and installation

The file `nbody/galaxies/sources/agama.py` is an **adapter** that constructs a three-component self-consistent equilibrium *candidate* using AGAMA's documented API. It does **not** contain or substitute for the AGAMA C++/Python library, and the backend has **not** been run-verified in the release environment.

## Check your environment first

Open the GUI's **Initial-condition model** accordion: it shows whether `import agama` succeeds **in the interpreter running the GUI**. You can click **Check AGAMA availability** again after installation. The message identifies the exact Python executable being used.

You can independently check with:

```powershell
.\.venv\Scripts\python.exe -c "import sys; print(sys.executable); import agama; print(agama.__file__)"
```

To attempt installation in that same environment:

```powershell
.\.venv\Scripts\python.exe -m pip install -U pip setuptools wheel numpy
.\.venv\Scripts\python.exe -m pip install agama
```

**Windows build note:** AGAMA contains compiled C++ code. Installation may need Microsoft C++ Build Tools and a correctly initialized MSVC environment (e.g., an x64 Developer Command Prompt / `vcvarsall.bat`); dependencies such as GSL may also need compiling. Do not assume a standard pip installation will work on every Windows machine. See the **current** official installation guide before changing compiler options:

https://github.com/GalacticDynamics-Oxford/Agama/blob/master/INSTALL

Official three-component reference example:
https://github.com/GalacticDynamics-Oxford/Agama/blob/master/py/example_self_consistent_model3.py

## What to verify once AGAMA imports

1. Generate a candidate; check that disk/bulge/halo sampled masses and density profiles match the requested parameters.
2. Compare orbit speeds and velocity dispersions against the equilibrium DF solution.
3. Evolve in isolation for several orbital periods; track structural and kinematic drift with V1.4.1 diagnostics.
4. Repeat at different particle counts, time steps, softening lengths and random seeds.
5. Save an actual particle snapshot with metadata so every gravitational solver uses exactly the same initial realization.

AGAMA import success **does not imply** the adapter, masses, units, or sampled galaxy have been validated. The native preview is still explicitly non-equilibrium, and the GUI will **not** silently switch to it when AGAMA fails.
