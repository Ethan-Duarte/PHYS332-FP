"""Simple diagnostics; O(N^2) potential energy is for small validation runs."""
import numpy as np
from .particles import center


def energy(p, eps=0.01, G=1):
    kinetic = np.sum(p.m * np.sum(p.v**2, axis=1)) / 2
    i, j = np.triu_indices(len(p.m), k=1)
    r2 = np.sum((p.x[j] - p.x[i])**2, axis=1) + eps**2
    potential = -G*np.sum(p.m[i]*p.m[j]/np.sqrt(r2))
    return kinetic, potential


def momentum(p):
    return np.sum(p.m[:, None]*p.v, axis=0)


def angular_momentum(p):
    return np.sum(p.m[:, None]*np.cross(p.x, p.v), axis=0)


def half_mass_radius(p, component=None):
    mask = np.ones(len(p.m), dtype=bool) if component is None else p.component == component
    x, m = p.x[mask], p.m[mask]
    r = np.linalg.norm(x-np.average(x, axis=0, weights=m), axis=1)
    order = np.argsort(r)
    return r[order][np.searchsorted(np.cumsum(m[order]), 0.5*m.sum())]


def core_separation(p):
    # Bulge centers are less contaminated by tidal tails than all-particle COMs.
    cores = [np.average(p.x[(p.galaxy == i) & (p.component == 1)], axis=0,
                        weights=p.m[(p.galaxy == i) & (p.component == 1)])
             for i in (0, 1)]
    return np.linalg.norm(cores[1]-cores[0])
