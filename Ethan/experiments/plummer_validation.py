"""Quick finite-N equilibrium sanity check, NOT comprehensive validation."""
from functools import partial
import numpy as np
from nbody.galaxies.sources.native import plummer
from nbody.gravity.direct import acceleration
from nbody.integrate import evolve
from nbody.analysis import energy, half_mass_radius

p = plummer(300, seed=42)
E0 = sum(energy(p, eps=0.05))
r0 = half_mass_radius(p)
final = list(evolve(p, partial(acceleration, eps=0.05), 0.01, 20, 20))[-1][1]
E1 = sum(energy(final, eps=0.05))
print(f'Plummer half-mass radius: {r0:.4f} -> {half_mass_radius(final):.4f}')
print(f'Energy drift (short test): {(E1-E0)/abs(E0):+.3e}')
