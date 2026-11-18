"""Structural API tests of the adapter; NOT a real compiled-AGAMA physics test."""
import json
import sys
from pathlib import Path
from types import SimpleNamespace
import numpy as np

from nbody.galaxies.sources.agama import generate

ROOT = Path(__file__).resolve().parents[1]


class FakeDensity:
    def __init__(self, **kw):
        self.params = kw
    def totalMass(self):
        return float(self.params.get('mass', 0))


class FakeComponent:
    def __init__(self, density=None, df=None, **kw):
        self.density = density or FakeDensity(mass=df.totalMass())
        self.df = df
        self.params = kw


class FakePotential:
    def force(self, R, y, z):
        return np.array([-2. / R, 0, 0])  # omega^2 = 2/R^2; kappa^2=4/R^2 > 0


class FakeSCM:
    def __init__(self, **kw):
        self.parameters = kw
        self.components = []
        self.potential = FakePotential()
        self.af = object()
        self.iters = 0
    def iterate(self):
        self.iters += 1


class FakeDF:
    def __init__(self, **kw):
        self.params = kw
    def totalMass(self):
        if self.params['type'] == 'QuasiIsothermal':
            p = self.params
            return float(p['Sigma0'] * 2 * np.pi * p['Rdisk']**2 * 1.1)
        return self.params['density'].totalMass()


class FakeGalaxy:
    def __init__(self, potential, df, af):
        self.df = df
    def sample(self, n):
        rng = FakeAgama.rng
        x = rng.normal(size=(n, 6))
        m = np.full(n, self.df.totalMass()/n)
        return x, m


class FakeAgama:
    rng = np.random.default_rng(42)
    models = []
    densities = []
    @staticmethod
    def setRandomSeed(seed):
        FakeAgama.rng = np.random.default_rng(seed)
    @staticmethod
    def Density(**kw):
        d = FakeDensity(**kw)
        FakeAgama.densities.append(d)
        return d
    @staticmethod
    def SelfConsistentModel(**kw):
        m = FakeSCM(**kw)
        FakeAgama.models.append(m)
        return m
    @staticmethod
    def Component(**kw):
        return FakeComponent(**kw)
    @staticmethod
    def DistributionFunction(**kw):
        return FakeDF(**kw)
    GalaxyModel = FakeGalaxy


def test_supplied_agama_source_archive_is_packaged():
    assert (ROOT / 'vendor' / 'Agama-master.zip').exists()


def test_agama_adapter_uses_source_api_correctly(monkeypatch):
    FakeAgama.models, FakeAgama.densities = [], []
    monkeypatch.setitem(sys.modules, 'agama', FakeAgama)
    cfg = json.loads((ROOT/'presets/mice.json').read_text())
    cfg['particles'] = {'disk': 40, 'bulge': 20, 'halo': 60}
    (p1, p2), report = generate(cfg, with_report=True, realizations=2)
    assert len(p1.m) == len(p2.m) == 120
    assert report['n_realizations'] == 2
    assert report['iterations'] == 4
    assert report['model_radius_max'] >= 2 * cfg['halo']['cutoff']
    assert FakeAgama.densities[0].params['scaleHeight'] == -cfg['disk']['height']/2
    assert np.isclose(sum(p1.m), sum(cfg[k]['mass'] for k in ('disk','bulge','halo')))
    assert not np.allclose(p1.x, p2.x)
    assert np.allclose(np.average(p1.x,axis=0,weights=p1.m),0,atol=1e-12)
    # Two distinct seeds for independently sampled progenitors.
    assert report['sampled_masses'][0] == report['sampled_masses'][1]
    assert report['status'].startswith('DF candidate')


def test_single_galaxy_agama_interface(monkeypatch):
    monkeypatch.setitem(sys.modules, 'agama', FakeAgama)
    cfg = json.loads((ROOT/'presets/mice.json').read_text())
    p = generate(cfg)
    assert len(p.m) == sum(cfg['particles'].values())


def test_mock_agama_gui_encounter_is_reproducible(tmp_path,monkeypatch):
    """Test integration plumbing only; FakeAgama does not represent real equilibrium physics."""
    import app
    from pathlib import Path
    monkeypatch.setitem(sys.modules,'agama',FakeAgama)
    monkeypatch.setattr(app,'RESULTS',tmp_path)
    monkeypatch.setattr(app,'RUNS',tmp_path/'runs')
    monkeypatch.setattr(app,'GALAXIES',tmp_path/'galaxies')
    app.GALAXIES.mkdir()
    values=dict(zip(app.FIELD_NAMES,[field.value for field in app.FIELD_CONTROLS]))
    values.update(name='mock_agama', setup='Mice encounter',
                  model='Disk–bulge–halo',generator='AGAMA (equilibrium candidate)',
                  disk=50,bulge=10,halo=50,steps=8,save_every=4,export='None',
                  solver='Fast direct (Numba)')
    last=list(app.simulate(*[values[k] for k in app.FIELD_NAMES],progress=lambda *args,**kwargs:None))[-1]
    assert 'Simulation complete' in last[7],last[7]
    assert 'AGAMA DF model generated' in last[7]
    folder=next(app.RUNS.glob('mock_agama*'))
    meta=json.loads((folder/'config.json').read_text())
    assert meta['generation_report']['n_realizations'] == 2
    assert meta['correlated_realizations'] is False
    assert (folder/'agama_generation_report.json').exists()
    assert len(list(app.GALAXIES.glob('agama_*_?.npz'))) == 2
