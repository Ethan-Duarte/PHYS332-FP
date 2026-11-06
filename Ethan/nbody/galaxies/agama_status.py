"""Report whether AGAMA imports in the interpreter used by the GUI."""
import sys


def status():
    try:
        import agama
        return {'available': True, 'python': sys.executable,
                'module': getattr(agama, '__file__', 'unknown'),
                'message': 'AGAMA imported successfully. Galaxy model/runtime validation is still required.'}
    except Exception as exc:
        return {'available': False, 'python': sys.executable,
                'module': None, 'message': f'{type(exc).__name__}: {exc}'}


def message():
    s = status()
    if s['available']:
        return f'**AGAMA available:** `{s["module"]}`. This does **not** validate the equilibrium adapter.'
    return (f'**AGAMA unavailable:** `{s["message"]}`\n\n'
            f'GUI Python: `{s["python"]}`\n\n'
            'To install into this environment, use that Python executable with `-m pip install agama`. '
            'Windows installation may require MSVC C++ build tools; see `docs/AGAMA_SETUP.md`. '
            'No automatic fallback to preview will occur.')
