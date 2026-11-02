"""Kick-drift-kick leapfrog; yields snapshots rather than storing every step."""
from .particles import Particles


def evolve(p, force, dt, steps, save_every=10):
    x, v, m = p.x.copy(), p.v.copy(), p.m
    a = force(x, m)
    yield 0.0, Particles(x.copy(), v.copy(), m, p.component, p.galaxy)
    for step in range(1, steps + 1):
        v += 0.5 * dt * a                # kick
        x += dt * v                      # drift
        a = force(x, m)
        v += 0.5 * dt * a                # kick
        if step % save_every == 0 or step == steps:
            yield step*dt, Particles(x.copy(), v.copy(), m, p.component, p.galaxy)
