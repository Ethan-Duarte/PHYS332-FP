"""One-snapshot force comparison. Benchmark physics and timing independently."""
from time import perf_counter
import numpy as np
from .gravity.direct import acceleration as numpy_direct
from .gravity.fast_direct import acceleration as fast_direct
from .gravity.fast_barnes_hut import acceleration as fast_tree
from .gravity.barnes_hut import acceleration as python_tree


def compare(x, m, eps=.01, theta=.5):
    # Warm JIT with tiny arrays to separate compilation time from timing.
    fast_direct(x[:4], m[:4], eps=eps)
    fast_tree(x[:4], m[:4], eps=eps, theta=theta)
    reference = fast_direct(x, m, eps=eps)
    solvers = [('Fast direct (Numba)',fast_direct),('Barnes–Hut (Numba)',fast_tree)]
    if len(m) <= 1200:
        solvers.append(('Direct (NumPy)',numpy_direct))
    if len(m) <= 150:
        solvers.append(('Barnes–Hut (Python)',python_tree))
    rows = []
    for name, fn in solvers:
        t0 = perf_counter()
        a = fn(x, m, eps=eps, **({'theta':theta} if 'Barnes' in name else {}))
        elapsed = perf_counter()-t0
        abs_error = np.linalg.norm(a-reference, axis=1)
        rel_error = abs_error/np.maximum(np.linalg.norm(reference,axis=1),1e-12)
        rows.append([name, round(elapsed,4), round(float(np.median(rel_error))*100,3),
                     round(float(np.percentile(rel_error,95))*100,3),
                     round(float(np.max(abs_error)),5)])
    return rows
