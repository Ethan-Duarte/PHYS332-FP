"""V1.4 physics and GUI workflow regression tests."""
from functools import partial
import json
import numpy as np
from nbody.galaxies.sources.native import plummer, mice_preview
from nbody.particles import combine
from nbody.galaxies.generate import galaxy
from nbody.diagnostics import total_energy, inspect_isolation, profile_checks, structure
from nbody.gravity.fast_direct import acceleration
from nbody.integrate import evolve
from nbody.solver_study import compare_trajectories
from nbody.encounters import parabolic_pair
from nbody.viewer import write_viewer


def test_energy_matches_pairwise_potential():
    p=plummer(40,mass=1.,a=.2,seed=7)
    kinetic=.5*np.sum(p.m*np.sum(p.v*p.v,axis=1))
    potential=0.
    for i in range(len(p.m)):
        for j in range(i):
            d=p.x[i]-p.x[j]
            potential-=p.m[i]*p.m[j]/np.sqrt(np.dot(d,d)+.01**2)
    np.testing.assert_allclose(total_energy(p,.01),kinetic+potential,rtol=1e-12)


def test_isolation_outputs_and_3d_plummer(tmp_path):
    p=plummer(100,mass=1.25,a=.2,seed=2)
    snaps=list(evolve(p,partial(acceleration,eps=.015),dt=.002,steps=10,save_every=5))
    report,img=inspect_isolation(snaps,.015,tmp_path,'analytic Plummer DF')
    assert report['n_particles']==100
    assert 'NOT validated' in report['validation']
    assert img.exists() and (tmp_path/'isolation_metrics.csv').exists()
    path=write_viewer(tmp_path/'viewer.html',np.array([t for t,_ in snaps]),
                      np.array([s.x for _,s in snaps]),p.component,p.galaxy)
    html=path.read_text()
    assert '100 displayed' in html and 'Plotly.restyle' in html


def test_profile_checks_realization(tmp_path):
    from pathlib import Path
    config=json.loads((Path(__file__).resolve().parents[1]/'presets/mice.json').read_text())
    config['particles']=dict(disk=700,bulge=300,halo=600)
    p=mice_preview(config)
    report=profile_checks(p,config,tmp_path,'preview')
    assert set(report)=={'Disk R','Bulge r','Halo r'}
    assert all(v['max_cdf_discrepancy']<.16 for v in report.values()),report
    assert (tmp_path/'initial_profiles.png').exists()


def test_matched_solver_study(tmp_path):
    from pathlib import Path
    config=json.loads((Path(__file__).resolve().parents[1]/'presets/mice.json').read_text())
    config['particles']=dict(disk=30,bulge=10,halo=30)
    one=galaxy(config,'preview');two=galaxy({**config,'seed':43},'preview')
    initial=parabolic_pair(one,two,1.,.25,config['encounter']['orientations'])
    report,plot=compare_trajectories(initial,tmp_path,dt=.003,steps=4,save_every=2)
    assert set(report)=={'Fast direct (Numba)','Barnes-Hut (Numba)','Direct (NumPy)'}
    assert report['Fast direct (Numba)']['rms_position_error_final']==0
    assert report['Direct (NumPy)']['force_p95_error']<1e-10
    assert report['Barnes-Hut (Numba)']['force_p95_error']>=0
    assert plot.exists() and (tmp_path/'solver_comparison.json').exists()


def test_gui_isolated_and_cached(tmp_path,monkeypatch):
    import app
    from nbody.io import save
    monkeypatch.setattr(app,'RESULTS',tmp_path)
    monkeypatch.setattr(app,'RUNS',tmp_path/'runs')
    monkeypatch.setattr(app,'GALAXIES',tmp_path/'galaxies')
    app.GALAXIES.mkdir()
    def run(**overrides):
        values=dict(zip(app.FIELD_NAMES,[x.value for x in app.FIELD_CONTROLS]))
        values.update(dict(name='study',setup='Isolated galaxy',model='Plummer (analytic DF)',
                           disk=50,bulge=10,halo=50,steps=10,save_every=5,export='None',
                           solver='Fast direct (Numba)'))
        values.update(overrides)
        return list(app.simulate(*[values[k] for k in app.FIELD_NAMES],
                                 progress=lambda *a,**kw:None))[-1]
    first=run()
    assert 'Simulation complete' in first[7] and 'screening' in first[7]
    assert (next(app.RUNS.glob('study*'))/'isolation_report.json').exists()
    from nbody.particles import Particles
    p=plummer(110,mass=1.25,a=.2,seed=8)
    # Cached snapshots are explicit data; no claim of an equilibrium disk is inferred.
    save(app.GALAXIES/'example.npz',p,dict(backend='external snapshot'))
    second=run(model='Disk–bulge–halo',generator='Saved snapshot',snapshot='example.npz')
    assert 'Simulation complete' in second[7] and 'saved: external snapshot' in second[7]
