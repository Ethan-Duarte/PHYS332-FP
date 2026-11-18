#!/bin/bash
# Optional: compile the supplied AGAMA source into the GUI's Python environment.
set -u
cd "$(dirname "$0")/.." || exit 1

if [[ "$(uname -s)" != 'Darwin' ]]; then
    echo 'This AGAMA installation script is for macOS only.'
    exit 1
fi

PYTHON='.venv/bin/python'
if [[ ! -x "$PYTHON" ]]; then
    echo 'Launch N-Baddies once with start_gui.command to create .venv first.'
    exit 1
fi
if ! xcrun --find clang++ >/dev/null 2>&1; then
    echo 'C++ compiler missing. Install Xcode Command Line Tools with: xcode-select --install'
    exit 1
fi

printf '\n=== AGAMA source build (macOS) ===\n'
"$PYTHON" -c 'import platform,sys; print("Python:",sys.executable); print("Architecture:",platform.machine())'
echo "Host architecture: $(uname -m)"
printf 'The Python and host architectures should match (arm64 with arm64, or x86_64 with x86_64).\n'

# Apple's compiler does not automatically search Homebrew's headers/libraries.
if command -v brew >/dev/null 2>&1; then
    for pkg in gsl eigen libomp; do
        prefix="$(brew --prefix "$pkg" 2>/dev/null || true)"
        if [[ -d "$prefix/include" ]]; then
            export CFLAGS="${CFLAGS:-} -I$prefix/include"
            if [[ "$pkg" == 'eigen' && -d "$prefix/include/eigen3" ]]; then
                export CFLAGS="$CFLAGS -I$prefix/include/eigen3"
            fi
        fi
        if [[ -d "$prefix/lib" ]]; then
            export LDFLAGS="${LDFLAGS:-} -L$prefix/lib"
        fi
    done
else
    echo 'Homebrew not found; AGAMA may need to download/build its own GSL dependency.'
    echo 'Recommended prerequisite: install Homebrew, then: brew install gsl eigen libomp'
fi

"$PYTHON" -m pip install -U pip setuptools wheel numpy || exit 1
"$PYTHON" -m pip install -v --no-build-isolation vendor/Agama-master.zip || exit 1
"$PYTHON" -c 'import agama; print("AGAMA module:", agama.__file__); print("G:", agama.G)' || exit 1
printf '\nAGAMA imports successfully. Now test its real modeling API:\n'
printf '  %s -m experiments.agama_selftest\n' "$PYTHON"
printf 'A successful import is not evidence of a stable disk galaxy.\n'
