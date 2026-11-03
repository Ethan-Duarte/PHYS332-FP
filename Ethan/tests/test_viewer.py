import numpy as np
from nbody.viewer import figure3d
from nbody.gravity.barnes_hut import acceleration as python_tree
from nbody.gravity.fast_barnes_hut import acceleration as fast_tree
from nbody.benchmarks import compare


def test_3d_viewer_has_animated_frames():
    rng=np.random.default_rng(4)
    x=rng.normal(size=(4,40,3))
    comp=np.r_[np.zeros(20),np.ones(20)]
    ids=np.r_[np.zeros(20),np.ones(20)]
    fig=figure3d(np.arange(4),x,comp,ids)
    assert len(fig.frames)==4
    assert fig.data[0].type=='scatter3d'
    assert len(fig.layout.sliders[0].steps)==4


def test_accelerated_barnes_hut_matches_original():
    rng=np.random.default_rng(3)
    x=rng.normal(size=(80,3)); m=rng.uniform(.01,1,size=80)
    for theta in [.25,.5,.9]:
        np.testing.assert_allclose(fast_tree(x,m,.025,theta),
                                   python_tree(x,m,.025,theta),rtol=1e-12,atol=1e-12)


def test_benchmark_has_direct_reference():
    rng=np.random.default_rng(9)
    x=rng.normal(size=(24,3));m=rng.random(24)
    rows=compare(x,m,.03,.5)
    assert rows[0][0]=='Fast direct (Numba)'
    assert rows[0][2]==0
    assert 'Barnes' in rows[1][0]
