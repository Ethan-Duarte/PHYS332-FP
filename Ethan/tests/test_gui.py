"""GUI end-to-end and cancellation tests with small particle counts."""
from pathlib import Path
import pytest
pytest.importorskip('gradio')
import app


def params(**updates):
    p = dict(zip(app.FIELD_NAMES, [c.value for c in app.FIELD_CONTROLS]))
    p.update(dict(name='test_gui', disk=50, bulge=10, halo=50, steps=20,
             save_every=5, solver='Barnes–Hut (Numba)', export='None')); p.update(updates)
    return [p[k] for k in app.FIELD_NAMES]


def test_run_and_replay(tmp_path, monkeypatch):
    monkeypatch.setattr(app,'RESULTS',tmp_path)
    monkeypatch.setattr(app,'RUNS',tmp_path/'runs')
    final = list(app.simulate(*params(),progress=lambda *a,**k:None))[-1]
    assert "<iframe" in final[0] and "viewer_3d.html" in final[0]
    assert Path(final[4]).exists() and Path(final[5]).exists() and Path(final[6]).exists()
    assert final[2] is None and final[3] is None
    assert 'Simulation complete' in final[7]
    selected = str(Path(final[4]).relative_to(tmp_path))
    replay = app.replay_saved(selected,1.5,False,1000,15,'None',progress=lambda *a,**k:None)
    assert "<iframe" in replay[0] and Path(replay[4]).exists()


def test_custom_vectors_preserve_relative_state():
    from nbody.particles import Particles
    from nbody.encounters import relative_pair
    import numpy as np
    one = Particles(np.zeros((1,3)),np.zeros((1,3)),np.array([1.]),np.array([1]),np.array([0]))
    two = Particles(np.zeros((1,3)),np.zeros((1,3)),np.array([2.]),np.array([1]),np.array([0]))
    p = relative_pair(one,two,[1,2,3],[-1,.5,.2])
    np.testing.assert_allclose(p.x[1]-p.x[0], [1,2,3])
    np.testing.assert_allclose(p.v[1]-p.v[0], [-1,.5,.2])
    np.testing.assert_allclose(p.m @ p.x, [0,0,0],atol=1e-14)
    np.testing.assert_allclose(p.m @ p.v, [0,0,0],atol=1e-14)


def test_cancel_saves_partial_run(tmp_path, monkeypatch):
    monkeypatch.setattr(app,'RESULTS',tmp_path)
    monkeypatch.setattr(app,'RUNS',tmp_path/'runs')
    results = app.simulate(*params(steps=100),progress=lambda *a,**k:None)
    for item in results:
        if 'frames saved' in str(item[7]):
            app.cancel_active()
            break
    final = list(results)[-1]
    assert 'Cancelled' in final[7]
    assert Path(final[4]).exists()
    from nbody.io import load_trajectory
    assert load_trajectory(final[4])[-1]['cancelled'] is True
