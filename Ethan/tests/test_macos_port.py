"""Packaging and routing checks; not substitutes for testing on macOS."""
from pathlib import Path
import subprocess
import sys
import os
import shutil
from nbody.galaxies import agama_status

ROOT = Path(__file__).resolve().parents[1]


def test_unix_launchers_are_well_formed():
    paths = [ROOT / 'start_gui.sh', ROOT / 'start_gui.command',
             ROOT / 'scripts/install_agama_macos.sh']
    for path in paths:
        assert path.exists() and path.read_text().startswith('#!/bin/bash')
        if os.name != 'nt':
            assert path.stat().st_mode & 0o111
        if shutil.which('bash'):
            subprocess.run(['bash', '-n', str(path)], check=True)


def test_macos_launcher_installs_into_local_environment():
    s = (ROOT / 'start_gui.sh').read_text()
    assert '.venv/bin/python' in s
    assert 'launch_gui.py' in s
    assert 'requirements-fast.txt' in s
    assert 'requirements-gui.txt' in s
    assert 'requirements-viewer.txt' in s


def test_agama_status_platform_guidance(monkeypatch):
    monkeypatch.setattr(agama_status, 'status', lambda: {
        'available': False, 'python': '/tmp/.venv/bin/python',
        'message': 'ModuleNotFoundError: No module named agama', 'source_included': True})
    monkeypatch.setattr(sys, 'platform', 'darwin')
    assert 'install_agama_macos.sh' in agama_status.message()
    monkeypatch.setattr(sys, 'platform', 'win32')
    assert 'install_agama_windows.bat' in agama_status.message()
