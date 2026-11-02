"""Vectorized direct O(N^2) Plummer-softened gravity. G=1 by default."""
import numpy as np


def acceleration(x, m, eps=0.01, G=1.0):
    d = x[None, :, :] - x[:, None, :]
    r2 = np.sum(d*d, axis=2) + eps**2
    np.fill_diagonal(r2, np.inf)
    return G * np.sum(m[None, :, None] * d / r2[:, :, None]**1.5, axis=1)
