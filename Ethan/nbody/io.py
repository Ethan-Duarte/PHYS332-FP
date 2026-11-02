"""Portable NumPy snapshots with provenance metadata."""
import json
import numpy as np
from .particles import Particles


def save(path, p, metadata=None):
    np.savez_compressed(path, x=p.x, v=p.v, m=p.m,
                        component=p.component, galaxy=p.galaxy,
                        metadata=json.dumps(metadata or {}))


def load(path):
    with np.load(path) as data:
        p = Particles(*(data[k].copy() for k in ('x', 'v', 'm', 'component', 'galaxy')))
        return p, json.loads(data['metadata'].item())
