"""The acceleration must be identical when only target-loop parallelism changes."""
import numba
import numpy as np
from nbody.gravity.fast_direct import acceleration as serial
from nbody.gravity.parallel_direct import acceleration as parallel


def test_parallel_matches_serial_bitwise():
    rng = np.random.default_rng(3)
    x = rng.normal(size=(256,3))
    m = rng.uniform(.1,1,256)
    previous = numba.get_num_threads()
    try:
        for epsilon in (.0015625,.00625):
            expected = serial(x,m,epsilon)
            for count in (1,2,4):
                numba.set_num_threads(min(count,numba.config.NUMBA_NUM_THREADS))
                assert np.array_equal(parallel(x,m,epsilon),expected)
    finally:
        numba.set_num_threads(previous)
