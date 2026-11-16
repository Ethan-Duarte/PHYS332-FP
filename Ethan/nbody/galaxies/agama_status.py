"""Inspect the active Python environment; a source zip is not a compiled AGAMA module."""
import sys
from pathlib import Path

SOURCE = Path(__file__).resolve().parents[2] / 'vendor' / 'Agama-master.zip'


def status():
    try:
        import agama
        return {'available': True, 'python': sys.executable,
                'module': getattr(agama, '__file__', 'unknown'),
                'source_included': SOURCE.exists(),
                'message': 'AGAMA imported; adapter and isolated stability still require verification.'}
    except Exception as exc:
        return {'available': False, 'python': sys.executable, 'module': None,
                'source_included': SOURCE.exists(),
                'message': f'{type(exc).__name__}: {exc}'}


def message():
    s = status()
    if s['available']:
        return (f'**AGAMA available:** `{s["module"]}`. This does **not** validate a particle model. '
                'Generate a candidate, inspect masses, and run isolated stability tests.')
    source_note = ('The AGAMA-master source archive is included under `vendor/`, but it must be **compiled**.'
                   if s['source_included'] else 'The source archive was not found under `vendor/`.')
    if sys.platform == 'darwin':
        install = ('On macOS, run `bash scripts/install_agama_macos.sh` from Terminal '
                   'after installing Xcode Command Line Tools and Homebrew `gsl`, `eigen`, `libomp`.')
    elif sys.platform == 'win32':
        install = ('On Windows, run `scripts/install_agama_windows.bat` from an '
                   '**x64 Native Tools** developer prompt.')
    else:
        install = 'On Linux, build AGAMA in the active `.venv` using the upstream INSTALL guide.'
    return (f'**AGAMA unavailable:** `{s["message"]}`\n\n'
            f'GUI Python: `{s["python"]}`\n\n'
            f'{source_note} {install} AGAMA requires a C++ toolchain and GSL; '
            'see `docs/AGAMA_SETUP.md`. We will **not** silently substitute native preview particles.')
