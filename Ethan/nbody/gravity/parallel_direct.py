"""Direct Plummer gravity parallelized only over independent target particles.

No tree, changed kernel, fastmath, or reduction across threads. Each target
retains the serial implementation's source-particle summation order.
"""
import numpy as np
from numba import njit, prange


@njit(cache=True, parallel=True)
def acceleration(x, m, eps=0.01, G=1.0):
    n = len(m)
    a = np.zeros((n, 3))
    for i in prange(n):
        ax, ay, az = 0., 0., 0.
        for j in range(n):
            if i == j:
                continue
            dx = x[j, 0] - x[i, 0]
            dy = x[j, 1] - x[i, 1]
            dz = x[j, 2] - x[i, 2]
            r2 = dx*dx + dy*dy + dz*dz + eps*eps
            factor = G*m[j]/(r2*np.sqrt(r2))
            ax += dx*factor
            ay += dy*factor
            az += dz*factor
        a[i, 0], a[i, 1], a[i, 2] = ax, ay, az
    return a
