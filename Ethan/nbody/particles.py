"""Minimal particle representation shared by every module."""
from dataclasses import dataclass
import numpy as np


@dataclass
class Particles:
    x: np.ndarray                 # (N, 3) positions
    v: np.ndarray                 # (N, 3) velocities
    m: np.ndarray                 # (N,) masses
    component: np.ndarray = None  # 0 disk, 1 bulge, 2 halo, 3 test sphere
    galaxy: np.ndarray = None     # original galaxy ID, retained after encounters


def combine(*systems):
    return Particles(*(np.concatenate([getattr(p, key) for p in systems])
                       for key in ('x', 'v', 'm', 'component', 'galaxy')))


def center(p):
    return np.average(p.x, axis=0, weights=p.m)


def recenter(p):
    x0 = center(p)
    v0 = np.average(p.v, axis=0, weights=p.m)
    return Particles(p.x - x0, p.v - v0, p.m, p.component, p.galaxy)
