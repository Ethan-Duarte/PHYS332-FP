"""Global scalar virial diagnostics for isolated, softened N-body systems.

With Plummer force softening, W=sum(m_i r_i·a_i), not the potential energy U.
A global virial ratio near 1 is necessary, not sufficient, for equilibrium.
"""
import numpy as np
from .particles import Particles


def virial_state(p, force):
    m = p.m
    center = np.average(p.x, axis=0, weights=m)
    mean_velocity = np.average(p.v, axis=0, weights=m)
    v = p.v - mean_velocity
    a = force(p.x, m)
    kinetic = 0.5 * np.sum(m * np.sum(v*v, axis=1))
    virial_term = np.sum(m * np.sum((p.x-center)*a, axis=1))
    return dict(kinetic=float(kinetic), virial_term=float(virial_term),
                ratio=float(2*kinetic / abs(virial_term)),
                center=center.tolist(), mean_velocity=mean_velocity.tolist())


def rescale_velocities(p, force, target=1.0):
    """Set *global* 2K/|W| to target without modifying the input particles."""
    before = virial_state(p, force)
    if before['kinetic'] <= 0 or before['virial_term'] >= 0 or target <= 0:
        raise ValueError('Velocity rescaling needs a bound, moving system and positive target.')
    vbar = np.asarray(before['mean_velocity'])
    scale = np.sqrt(target / before['ratio'])
    new_p = Particles(p.x.copy(), vbar + (p.v-vbar)*scale, p.m.copy(),
                      p.component.copy(), p.galaxy.copy())
    after = virial_state(new_p, force)
    return new_p, dict(before=before, after=after, speed_factor=float(scale),
                       caution='Global virial balance does not guarantee local Jeans equilibrium, a stable disk, or the original rotation curve.')
