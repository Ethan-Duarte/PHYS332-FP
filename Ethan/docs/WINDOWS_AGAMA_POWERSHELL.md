# Install AGAMA using Windows PowerShell — N-Baddies V1.6

The AGAMA source zip is bundled at `vendor/Agama-master.zip`. It is **not** a precompiled module. AGAMA's `setup.py` compiles C++ with MSVC and requires GSL. Having the zip in your project folder is not sufficient for `import agama`.

## Prerequisites (one-time)

- Windows 10/11, x64; Python 3.12 x64 is recommended for consistency with the project.
- **Visual Studio Build Tools 2022**, from https://visualstudio.microsoft.com/visual-cpp-build-tools/ . In Visual Studio Installer select **Desktop development with C++**. Include MSVC v143 x64/x86 tools and a Windows SDK; restart PowerShell after installation.
- Network access for Python dependencies and potentially GSL. AGAMA's Windows source setup may attempt to download an external GSL port from an upstream server. If that URL is unavailable, the source build can fail and will need manual GSL provisioning.

## Commands to copy into PowerShell

**1.** Open Windows PowerShell in the extracted N-Baddies project folder (`nbody_v1_5`; it has `app.py`, `start_gui.bat` and `scripts` in it). You can type `cd ` and drag the folder into the console. For example:

```powershell
Set-Location 'C:\Users\YOU\Documents\nbody_v1_5'
```

Replace the sample path with your actual folder. Check:

```powershell
Get-ChildItem .\start_gui.bat, .\vendor\Agama-master.zip
```

**2.** Create the virtual environment if it does not already exist. If `start_gui.bat` worked, skip this step.

```powershell
py -3.12 -m venv .venv
```

**3.** Install project dependencies if the GUI has not already done so:

```powershell
.\.venv\Scripts\python.exe -m pip install -r requirements.txt -r requirements-fast.txt -r requirements-gui.txt -r requirements-viewer.txt
```

**4.** Run the new PowerShell installer, which finds Visual Studio's x64 compiler, initializes its environment and invokes AGAMA's `setup.py` through pip:

```powershell
powershell.exe -ExecutionPolicy Bypass -File .\scripts\install_agama_windows.ps1
```

The `-ExecutionPolicy Bypass` applies to this child PowerShell process; it does not change your persistent system policy. The installer builds **from the bundled zip** and into **this project's `.venv`**, not a different Python installation.

**5.** Verify the module import:

```powershell
.\.venv\Scripts\python.exe -c "import agama; print(agama.__file__); print('G=', agama.G)"
```

**6.** Exercise the actual galaxy-generation API (may take some time):

```powershell
.\.venv\Scripts\python.exe -m experiments.agama_selftest
```

**7.** Restart `start_gui.bat`, select the AGAMA generator, and click **Check AGAMA availability**. A successfully generated galaxy is **an equilibrium candidate**, not yet an isolation-validated galaxy. Save it and run **Test a saved equilibrium galaxy**, then the new **Virialize & relax galaxies** comparisons.

## Troubleshooting

- `py -3.12` not found: install Python 3.12 x64, or use your already working `.venv` and skip environment creation.
- `Visual Studio Build Tools not found` or `MSVC C++ tools not found`: install the C++ workload in Visual Studio Installer, then rerun step 4.
- `GSL library not found`: the AGAMA source may ask permission to download and compile GSL. Agree if you trust the upstream source. Its optional download currently depends on an external host and may fail; use the installer output to diagnose the precise issue.
- Compiler error: attach the last 60–100 lines of the PowerShell build output. Do not install AGAMA in another unrelated Python interpreter to fix this.
- `ModuleNotFoundError: agama` after install: run step 5 using the **same** `.venv` interpreter. If the build failed, the extension was never installed.
- `agama_selftest` fails *after a successful import*: the source build is working, but the model-generation adapter or configuration still needs debugging. Share that traceback; do not switch to native preview and call it equilibrium.

## Source authority

The installation flow follows the bundled `Agama-master/INSTALL` and https://github.com/GalacticDynamics-Oxford/Agama/blob/master/INSTALL . The upstream guide describes `pip install .`, `--no-build-isolation`, the `vcvarsall`/MSVC environment requirement, and GSL on Windows. This project automates activation of MSVC only; **the Windows build has not been tested on a real Windows PC in the build environment**.
