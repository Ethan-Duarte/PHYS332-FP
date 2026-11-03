@echo off
setlocal EnableExtensions
title N-Baddies Galaxy Simulator
cd /d "%~dp0"

echo ==========================================
echo N-Baddies - Local Galaxy Simulator
echo ==========================================
echo Project folder: %CD%
echo.

set "PYTHON=.venv\Scripts\python.exe"
if exist "%PYTHON%" goto :check_deps

echo Creating a local Python environment...
where py >nul 2>nul
if errorlevel 1 goto :try_python

py -3.12 -m venv .venv
if exist "%PYTHON%" goto :check_deps
py -3.11 -m venv .venv
if exist "%PYTHON%" goto :check_deps

:try_python
python -m venv .venv
if not exist "%PYTHON%" goto :python_error

:check_deps
echo Using:
"%PYTHON%" --version
"%PYTHON%" -c "import gradio, numba, numpy, matplotlib, PIL" >nul 2>nul
if not errorlevel 1 goto :start_server

echo.
echo Installing missing dependencies. This may take a while on first launch.
"%PYTHON%" -m pip install -r requirements.txt -r requirements-fast.txt -r requirements-gui.txt
if errorlevel 1 goto :install_error

:start_server
echo.
echo Starting N-Baddies...
echo If the browser does not open, visit http://127.0.0.1:7860
echo Leave this window open while using the interface.
echo.
"%PYTHON%" -u launch_gui.py
if errorlevel 1 goto :server_error

echo.
echo The server exited normally.
goto :finish

:python_error
echo.
echo ERROR: Could not create the Python environment.
echo Install Python 3.11 or 3.12, then run this launcher again.
goto :finish

:install_error
echo.
echo ERROR: Dependency installation failed. Review the messages above.
echo Try deleting the .venv folder, then start again if needed.
goto :finish

:server_error
echo.
echo ERROR: The GUI did not start. Review the message above.
echo A traceback may also be saved to startup_error.log.
goto :finish

:finish
echo.
echo Press any key to close this window.
pause >nul
endlocal
