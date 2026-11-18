import json
import numpy as np

from nbody.particles import Particles
from nbody.gravity.fast_direct import acceleration
from nbody.galaxies.sources.native import plummer
from nbody.virial import virial_state, rescale_velocities
from nbody.initial_conditions import run_relax_restart, force_for
from nbody.io import load, load_trajectory


def test_virial_term_is_softened_force_not_potential_energy():
    p = Particles(np.array([[-1.,0,0],[1.,0,0]]),
                  np.array([[0,.2,0],[0,-.2,0]]), np.ones(2),
                  np.full(2, 3), np.zeros(2, dtype=int))
    eps = .2
    v = virial_state(p,lambda x,m: acceleration(x,m,eps))
    expected_W = -4/(4+eps**2)**1.5
    assert np.isclose(v['virial_term'],expected_W,rtol=1e-12)
    assert v['ratio'] > 0


def test_virial_rescaling_preserves_geometry_and_bulk_velocity():
    p=plummer(140, mass=1., a=1., seed=13)
    force=force_for('Fast direct',.05,.5)
    q,report=rescale_velocities(p,force,target=1.)
    assert np.array_equal(p.x,q.x)
    assert np.isclose(report['after']['ratio'],1.,rtol=1e-12)
    assert np.allclose(report['before']['mean_velocity'],report['after']['mean_velocity'])
    assert not np.array_equal(p.v,q.v)


def test_relax_then_restart_and_save(tmp_path):
    p=plummer(120,mass=1.,a=1.,seed=7)
    r=run_relax_restart(p,tmp_path,eps=.04,dt=.01,relax_steps=8,
                        restart_steps=6,save_every=3,solver='Fast direct',
                        target_virial=1.)
    assert r['status']=='complete'
    assert abs(r['virial_ratio']['prepared']-1.)<1e-10
    assert (tmp_path/'relaxed_candidate.npz').exists()
    assert (tmp_path/'relax_restart.png').exists()
    q,meta=load(tmp_path/'relaxed_candidate.npz')
    assert meta['validation'].startswith('candidate')
    t,x,v,m,comp,ids,meta=load_trajectory(tmp_path/'restart_trajectory.npz')
    assert t[0]==0 and np.isclose(t[-1],.06)
    assert np.allclose(x[0],q.x)
    assert json.loads((tmp_path/'initial_condition_report.json').read_text())['status']=='complete'


def test_cancel_relax_does_not_save_unfinished_candidate(tmp_path):
    from threading import Event
    event=Event();event.set()
    p=plummer(80,mass=1.,a=1.,seed=1)
    r=run_relax_restart(p,tmp_path,dt=.01,relax_steps=20,restart_steps=20,
                        save_every=5,cancelled=event)
    assert r['status']=='cancelled'
    assert not (tmp_path/'relaxed_candidate.npz').exists()


def test_gui_source_and_candidate_creation(tmp_path,monkeypatch):
    import app
    p=plummer(100, mass=1, a=1, seed=2)
    monkeypatch.setattr(app,'GALAXIES',tmp_path)
    from nbody.io import save
    save(tmp_path/'base.npz',p,{'backend':'analytic plummer'})
    msg,npz,report=app.virialize_candidate('Saved galaxy','base.npz',None,'First frame',
                                            1.,.04,'Fast direct',.5)
    assert 'Saved' in msg
    assert np.isclose(json.loads(open(report).read())['after']['ratio'],1.)
