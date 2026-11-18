"""V1.5 isolated refinement and GUI regression tests."""
import json
import numpy as np
from nbody.galaxies.sources.native import plummer
from nbody.validation import validate_initial_state
from nbody.io import save


def test_isolation_refinement_saves_identical_end_times(tmp_path):
    p=plummer(60, mass=1.0, a=1.0, seed=8)
    report=validate_initial_state(p,tmp_path,dt=.005,steps=12,eps=.025,save_every=4)
    assert (tmp_path/'validation_comparison.json').exists()
    assert (tmp_path/'base'/'trajectory.npz').exists()
    assert (tmp_path/'half_dt'/'trajectory.npz').exists()
    assert np.isclose(report['base']['duration'],report['half_dt']['duration'])
    assert report['refined_comparison']['position_rms'] >= 0
    assert 'not a validation' in report['verdict']


def test_gui_saved_galaxy_screen(tmp_path, monkeypatch):
    import app
    folder=tmp_path/'galaxies'
    folder.mkdir()
    save(folder/'test.npz',plummer(60,seed=13),{'backend':'analytic Plummer DF'})
    monkeypatch.setattr(app,'GALAXIES',folder)
    monkeypatch.setattr(app,'RUNS',tmp_path/'runs')
    txt,report,a,b=app.validate_saved_galaxy('test.npz',.005,20,.03,'Fast direct',progress=lambda *args,**kwargs:None)
    assert 'Single-run' in txt
    assert report.endswith('.json')
    assert a.endswith('.png') and b.endswith('.png')
