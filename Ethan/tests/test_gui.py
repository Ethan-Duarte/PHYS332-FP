"""Tiny integration checks: UI uses the existing simulation and replay modules."""
from pathlib import Path
import pytest

pytest.importorskip('gradio')

import app


def test_gui_run_and_replay(tmp_path, monkeypatch):
    monkeypatch.setattr(app, 'RESULTS', tmp_path)
    monkeypatch.setattr(app, 'RUNS', tmp_path / 'runs')
    options = dict(zip(app.FIELD_NAMES, [c.value for c in app.FIELD_CONTROLS]))
    options.update(name='gui_test', disk=50, bulge=10, halo=50,
                   solver='Direct (NumPy)', steps=20, save_every=10,
                   export='GIF', view='xy')
    progress = lambda *args, **kwargs: None
    gif, mp4, trajectory, settings, figure, status = app.simulate(
        *(options[k] for k in app.FIELD_NAMES), progress=progress)
    assert Path(gif).exists() and Path(trajectory).exists() and Path(settings).exists()
    assert mp4 is None and '220 live particles' in status
    selected = str(Path(trajectory).relative_to(tmp_path))
    assert selected in app.available_runs()
    replay = app.replay_saved(selected, 'xz', 2.0, 12, False, 'GIF', progress=progress)
    assert Path(replay[0]).exists()
    assert '220 live particles' in replay[4]


def test_gui_rejects_invalid_orbit():
    options = dict(zip(app.FIELD_NAMES, [c.value for c in app.FIELD_CONTROLS]))
    options['pericenter'] = options['separation']
    with pytest.raises(Exception, match='Pericenter'):
        app.simulate(*(options[k] for k in app.FIELD_NAMES), progress=lambda *args, **kwargs: None)
