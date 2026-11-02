"""Run: python -m experiments.force_benchmark"""
import time
import numpy as np
from nbody.gravity import direct, barnes_hut

rng = np.random.default_rng(42)
for n in (32, 64, 128, 256):
    x = rng.normal(size=(n, 3))
    m = rng.random(n)/n
    t0 = time.perf_counter()
    reference = direct.acceleration(x, m, eps=0.02)
    direct_ms = 1e3*(time.perf_counter()-t0)
    for theta in (0.3, 0.5, 0.8):
        t0 = time.perf_counter()
        approx = barnes_hut.acceleration(x, m, eps=0.02, theta=theta)
        tree_ms = 1e3*(time.perf_counter()-t0)
        frac = np.linalg.norm(approx-reference, axis=1)/np.maximum(np.linalg.norm(reference, axis=1), 1e-15)
        print(f'N={n:4d}, theta={theta:.1f}: direct {direct_ms:7.2f}ms, '
              f'tree {tree_ms:7.2f}ms, median error {np.median(frac):.2%}, '
              f'95th percentile {np.percentile(frac,95):.2%}')
