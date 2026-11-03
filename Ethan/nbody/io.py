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


def save_trajectory(path, snapshots, metadata=None):
    """Save replayable frames and all per-particle labels; keep velocity histories for analysis."""
    times = np.array([t for t, _ in snapshots])
    states = [p for _, p in snapshots]
    np.savez_compressed(path, t=times,
                        x=np.stack([p.x for p in states]),
                        v=np.stack([p.v for p in states]),
                        m=states[0].m, component=states[0].component,
                        galaxy=states[0].galaxy,
                        metadata=json.dumps(metadata or {}))


def load_trajectory(path):
    with np.load(path) as d:
        return (d['t'].copy(), d['x'].copy(), d['v'].copy(),
                d['m'].copy(), d['component'].copy(), d['galaxy'].copy(),
                json.loads(d['metadata'].item()))
