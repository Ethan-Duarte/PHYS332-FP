#!/bin/bash
# N-Baddies cross-platform Unix launcher (macOS / Linux).
set -u
cd "$(dirname "$0")" || exit 1

pause_and_exit() {
    printf '\n%s\n' "$1"
    if [[ -t 0 ]]; then
        read -r -p 'Press Return to close this window...' _
    fi
    exit 1
}

printf '\n=== N-Baddies Galaxy Simulator ===\nFolder: %s\n' "$PWD"

PYTHON='.venv/bin/python'
if [[ ! -x "$PYTHON" ]]; then
    SYSTEM_PYTHON=''
    for candidate in python3.12 python3.11 python3.13 python3; do
        if command -v "$candidate" >/dev/null 2>&1 &&
           "$candidate" -c 'import sys; sys.exit(not ((3, 11) <= sys.version_info[:2] <= (3, 13)))' 2>/dev/null; then
            SYSTEM_PYTHON="$candidate"
            break
        fi
    done
    [[ -n "$SYSTEM_PYTHON" ]] || pause_and_exit 'Python 3.11–3.13 not found. Install Python 3.12 from python.org, then retry.'
    printf 'Creating Python environment with %s...\n' "$SYSTEM_PYTHON"
    "$SYSTEM_PYTHON" -m venv .venv || pause_and_exit 'Could not create .venv. Check the Python installation.'
fi

printf 'Python: '
"$PYTHON" --version || pause_and_exit 'The Python environment is invalid. Delete .venv and retry.'

if ! "$PYTHON" -c 'import gradio, numba, numpy, matplotlib, PIL, plotly' >/dev/null 2>&1; then
    printf '\nInstalling Python dependencies into .venv...\n'
    "$PYTHON" -m pip install -r requirements.txt -r requirements-fast.txt -r requirements-gui.txt -r requirements-viewer.txt ||
        pause_and_exit 'Dependency installation failed. Read the error above; see docs/MACOS_SETUP.md.'
fi

printf '\nLaunching local interface...\n'
printf 'If the browser does not open, visit http://127.0.0.1:7860\n'
printf 'Keep this Terminal window open while using the simulator.\n\n'
"$PYTHON" -u launch_gui.py || pause_and_exit 'GUI startup failed. See startup_error.log for the traceback.'

printf '\nGUI stopped.\n'
if [[ -t 0 ]]; then
    read -r -p 'Press Return to close this window...' _
fi
