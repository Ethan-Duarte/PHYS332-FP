"""Optional fast O(N²) direct gravity; O(N) temporary memory. Requires numba."""
import numpy as np
from numba import njit


@njit(cache=True)
def acceleration(x, m, eps=0.01, G=1.0):
    n = len(m)
    a = np.zeros((n, 3))
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            dx = x[j, 0] - x[i, 0]
            dy = x[j, 1] - x[i, 1]
            dz = x[j, 2] - x[i, 2]
            r2 = dx*dx + dy*dy + dz*dz + eps*eps
            factor = G * m[j] / (r2 * np.sqrt(r2))
            a[i, 0] += dx * factor
            a[i, 1] += dy * factor
            a[i, 2] += dz * factor
    return a
