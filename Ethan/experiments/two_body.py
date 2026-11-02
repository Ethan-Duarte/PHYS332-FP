"""Run from project root: python -m experiments.two_body"""
import numpy as np
from functools import partial
from nbody import Particles, evolve
from nbody.gravity import direct, barnes_hut
from nbody.analysis import energy


p = Particles(np.array([[-1.,0,0], [1.,0,0]]),
              np.array([[0.,-0.5,0], [0.,0.5,0]]), np.ones(2),
              np.zeros(2, dtype=int), np.zeros(2, dtype=int))
for label, solver in [('direct', direct.acceleration), ('tree', barnes_hut.acceleration)]:
    force = partial(solver, eps=0.0)
    snapshots = list(evolve(p, force, dt=2*np.pi/2000, steps=4000, save_every=4000))
    last = snapshots[-1][1]
    print(f'{label}: one orbit position error = {np.linalg.norm(last.x-p.x):.3e}, '
          f'energy drift = {energy(last, 0)[0]+energy(last, 0)[1] - sum(energy(p, 0)):.3e}')
