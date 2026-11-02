import json
from pathlib import Path
from functools import partial
import numpy as np
from nbody import Particles, evolve
from nbody.gravity import direct, barnes_hut
from nbody.analysis import energy, momentum, half_mass_radius, core_separation
from nbody.galaxies.sources.native import plummer, mice_preview
from nbody.galaxies.profiles import dehnen_mass, disk_mass
from nbody.encounters import parabolic_pair
from nbody.io import save, load


ROOT = Path(__file__).resolve().parents[1]


def binary():
    return Particles(np.array([[-1., 0, 0], [1., 0, 0]]),
                     np.array([[0., -.5, 0], [0., .5, 0]]), np.ones(2),
                     np.zeros(2, dtype=int), np.zeros(2, dtype=int))


def test_direct_binary_force():
    p = binary()
    assert np.allclose(direct.acceleration(p.x, p.m, 0), [[0.25, 0, 0], [-0.25, 0, 0]])


def test_direct_momentum_balance():
    rng = np.random.default_rng(17)
    x, m = rng.normal(size=(20,3)), rng.random(20)+.1
    assert np.allclose((direct.acceleration(x, m, .02)*m[:, None]).sum(axis=0), 0, atol=1e-13)


def test_tree_accuracy_and_self_interaction():
    rng = np.random.default_rng(11)
    x, m = rng.normal(size=(45, 3)), rng.uniform(.2, 2., size=45)
    exact = direct.acceleration(x, m, .01)
    assert np.allclose(barnes_hut.acceleration(x, m, .01, theta=0), exact, atol=1e-12)
    approx = barnes_hut.acceleration(x, m, .01, theta=.3)
    assert np.linalg.norm(approx-exact)/np.linalg.norm(exact) < .02


def test_two_body_orbit_and_energy():
    p = binary()
    snapshots = list(evolve(p, partial(direct.acceleration, eps=0),
                            4*np.pi/4000, 4000, save_every=4000))
    last = snapshots[-1][1]
    assert np.max(np.abs(last.x-p.x)) < 1e-4
    assert abs(sum(energy(last, 0))-sum(energy(p, 0))) < 1e-8


def test_plummer_spatial_cdf_and_velocities():
    p = plummer(8000, seed=10)
    radii = np.linalg.norm(p.x, axis=1)
    # Centroid subtraction introduces a small finite-sample translation.
    empirical = (radii < 1).mean()
    theoretical = 1/2**1.5
    assert abs(empirical-theoretical) < .025
    speeds = np.linalg.norm(p.v, axis=1)
    escape = np.sqrt(2/np.sqrt(radii**2+1))
    assert (speeds <= escape*1.04).all()


def test_plummer_seed():
    assert np.array_equal(plummer(80, seed=2).x, plummer(80, seed=2).x)
    assert not np.array_equal(plummer(80, seed=2).x, plummer(80, seed=3).x)


def test_profiles():
    assert np.isclose(dehnen_mass(1, 1, 1, gamma=1), 0.25)
    assert np.isclose(disk_mass(1, 1, 1), 1-2/np.e)


def test_preview_mass_and_encounter():
    cfg = json.loads((ROOT/'presets/mice.json').read_text())
    p1 = mice_preview(cfg)
    p2 = mice_preview(dict(cfg, seed=43))
    assert len(p1.m) == sum(cfg['particles'].values())
    assert np.isclose(sum(p1.m), 1.25)
    encounter = parabolic_pair(p1, p2, 1, .25, ((25, -30), (40, 60)))
    # Full progenitor centers have the requested separation; a 20-particle bulge fluctuates.
    for gid, expected in ((0, -.5), (1, .5)):
        mask = encounter.galaxy == gid
        assert np.isclose(np.average(encounter.x[mask, 0], weights=encounter.m[mask]), expected)
    assert np.isfinite(core_separation(encounter))
    assert np.linalg.norm(momentum(encounter)) < 1e-9


def test_io(tmp_path):
    p = plummer(12, seed=9)
    path = tmp_path/'example.npz'
    save(path, p, {'seed':9})
    loaded, meta = load(path)
    assert np.array_equal(loaded.x, p.x)
    assert np.array_equal(loaded.component, p.component)
    assert meta == {'seed':9}
