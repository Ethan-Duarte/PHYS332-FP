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


def test_standalone_viewer_contains_distinct_positions_and_controls(tmp_path):
    import json
    import re
    from nbody.viewer import write_viewer, viewer_iframe
    times = np.arange(3, dtype=float)
    positions = np.zeros((3, 4, 3))
    positions[:, 0, 0] = [0, 2, 5]
    positions[:, 2, 2] = [0, -1, -3]
    components = np.array([0, 0, 1, 1])
    origins = np.array([0, 0, 1, 1])
    path = write_viewer(tmp_path / 'viewer_3d.html', times, positions, components, origins)
    html = path.read_text()
    assert 'Plotly.restyle' in html and 'id="timeline"' in html and 'id="play"' in html
    payload = json.loads(re.search(r'const payload=(\{.*?\});\s*const frames=', html, re.S).group(1))
    assert len(payload['frames']) == 3
    assert payload['frames'][0] != payload['frames'][-1]
    assert '<iframe' in viewer_iframe(path)
