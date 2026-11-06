"""V1.4.1 Fourier-mode, saved replay, and AGAMA-status regression tests."""
import json
from pathlib import Path
from types import SimpleNamespace
import numpy as np

from nbody.particles import Particles
from nbody.spirals import disk_modes, spiral_study
from nbody.galaxies import agama_status


def ring(n=120, two_arms=False):
    a = np.linspace(0, 2*np.pi, n, endpoint=False)
    if two_arms:
        a = np.where(np.arange(n) % 2 == 0, 0.0, np.pi)
    r = np.full(n, 1.25)
    x = np.column_stack([r*np.cos(a), r*np.sin(a), np.zeros(n)])
    v = np.column_stack([-np.sin(a), np.cos(a), .03*np.sin(3*a)])
    v[:,0] += .03*np.sin(np.arange(n))
    return Particles(x, v, np.full(n,1/n), np.zeros(n,dtype=int), np.zeros(n,dtype=int))


def test_annular_modes_distinguish_arms_from_uniform_ring():
    even = disk_modes(ring(two_arms=False), rd=1., eps=.05, compute_q=False)
    arms = disk_modes(ring(two_arms=True), rd=1., eps=.05, compute_q=False)
    assert len(even) == 1 and len(arms) == 1
    assert even[0]['A2'] < 1e-10
    assert arms[0]['A2'] > .99
    assert arms[0]['A1'] < 1e-10
    np.testing.assert_allclose(arms[0]['noise_rms'], 1/np.sqrt(120), rtol=1e-12)


def test_spiral_study_outputs_peak_and_phase(tmp_path):
    states=[(0.,ring()),(.5,ring(two_arms=True)),(1.,ring())]
    report, plot=spiral_study(states, rd=1, eps=.05, folder=tmp_path)
    assert plot.exists() and report['time_of_A2_peak']==.5
    assert report['A2_peak'] > .99
    for name in ('spiral_modes.csv','spiral_annuli.csv','spiral_phase.png','spiral_report.json'):
        assert (tmp_path/name).exists()
    data=json.loads((tmp_path/'spiral_report.json').read_text())
    assert data['A2_peak']==report['A2_peak']


def test_agama_environment_status(monkeypatch):
    import sys
    monkeypatch.setitem(sys.modules,'agama',None)
    assert agama_status.status()['available'] is False
    assert 'unavailable' in agama_status.message()
    monkeypatch.setitem(sys.modules,'agama',SimpleNamespace(__file__='mock-agama.pyd'))
    assert agama_status.status()['available'] is True


def test_gui_isolated_disk_and_saved_analysis(tmp_path,monkeypatch):
    import app
    monkeypatch.setattr(app,'RESULTS',tmp_path)
    monkeypatch.setattr(app,'RUNS',tmp_path/'runs')
    monkeypatch.setattr(app,'GALAXIES',tmp_path/'galaxies')
    app.GALAXIES.mkdir()
    fields=dict(zip(app.FIELD_NAMES,[field.value for field in app.FIELD_CONTROLS]))
    fields.update(name='spiral_test',setup='Isolated galaxy',model='Disk–bulge–halo',
                  generator='Native preview (not equilibrium)',disk=500,bulge=20,halo=80,
                  solver='Fast direct (Numba)',steps=10,save_every=5,export='None')
    final=list(app.simulate(*[fields[k] for k in app.FIELD_NAMES],
                            progress=lambda *a,**kw: None))[-1]
    assert 'Simulation complete' in final[7],final[7]
    assert len(final)==10 and final[9] is not None
    assert Path(final[9]).exists()
    selected=str(Path(final[4]).relative_to(tmp_path))
    plot,phase,mode_csv,annulus_csv,message=app.analyze_saved_disk(selected)
    assert Path(plot).exists() and Path(phase).exists()
    assert Path(mode_csv).exists() and Path(annulus_csv).exists()
    assert 'A₂ peak' in message
