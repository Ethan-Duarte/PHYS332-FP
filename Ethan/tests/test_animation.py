"""Tests for the optional accelerated solver and repeatable animation output."""
import numpy as np
import pytest
from nbody.gravity.direct import acceleration as direct
from nbody.io import save_trajectory, load_trajectory
from nbody.animation import render
from nbody.galaxies.sources.native import plummer


def test_fast_direct():
    pytest.importorskip('numba')
    from nbody.gravity.fast_direct import acceleration as fast
    rng = np.random.default_rng(8)
    x = rng.normal(size=(75, 3))
    m = rng.uniform(.1, 2, 75)
    assert np.allclose(fast(x, m, .02), direct(x, m, .02), atol=1e-12)


def test_trajectory_and_gif(tmp_path):
    p = plummer(16, seed=2)
    states = [(0., p), (.1, type(p)(p.x + .1*p.v, p.v, p.m, p.component,
                                     np.concatenate([np.zeros(8, int), np.ones(8, int)])))]
    path = tmp_path / 'run.npz'
    save_trajectory(path, states, {'note': 'test'})
    t, x, v, m, comp, ids, meta = load_trajectory(path)
    assert x.shape == (2, 16, 3)
    assert np.allclose(t, [0, .1])
    assert meta == {'note': 'test'}
    gif = tmp_path / 'test.gif'
    render(t, x, comp, ids, gif, fps=4, lim=4, view='xy')
    assert gif.stat().st_size > 0
