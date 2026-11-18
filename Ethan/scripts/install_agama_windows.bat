@echo off
setlocal EnableExtensions
cd /d "%~dp0\.."
echo === AGAMA source build using N-Baddies .venv ===
if not exist ".venv\Scripts\python.exe" (
  echo ERROR: Run start_gui.bat first to create .venv.
  exit /b 1
)
where cl >nul 2>nul
if errorlevel 1 (
  echo ERROR: MSVC cl.exe not found in the current command prompt.
  echo Use scripts\install_agama_windows.ps1 in PowerShell, or launch x64 Native Tools Command Prompt.
  exit /b 1
)
where nmake >nul 2>nul
if errorlevel 1 (
  echo ERROR: nmake.exe not found. Install Desktop development with C++ in Visual Studio Build Tools.
  exit /b 1
)
echo Installing source into the SAME environment as the GUI.
".venv\Scripts\python.exe" -m pip install -U pip setuptools wheel numpy
if errorlevel 1 exit /b 1
if "%~1"=="" (
  set "AGAMA_SOURCE=vendor\Agama-master.zip"
) else (
  set "AGAMA_SOURCE=%~1"
)
".venv\Scripts\python.exe" -m pip install -v --no-build-isolation "%AGAMA_SOURCE%"
if errorlevel 1 exit /b 1
".venv\Scripts\python.exe" -c "import agama; print('AGAMA imported:', agama.__file__); print('G =', agama.G)"
if errorlevel 1 exit /b 1
echo AGAMA imported successfully. Next: .venv\Scripts\python.exe -m experiments.agama_selftest
exit /b 0
