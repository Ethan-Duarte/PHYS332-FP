"""Kick–drift–kick leapfrog, with optional cancellation between timesteps."""
from .particles import Particles


def evolve(p, force, dt, steps, save_every=10, cancelled=None):
    x, v, m = p.x.copy(), p.v.copy(), p.m
    a = force(x, m)
    yield 0.0, Particles(x.copy(), v.copy(), m, p.component, p.galaxy)
    for step in range(1, steps + 1):
        if cancelled is not None and cancelled.is_set():
            return
        v += 0.5 * dt * a
        x += dt * v
        a = force(x, m)
        v += 0.5 * dt * a
        if step % save_every == 0 or step == steps:
            yield step*dt, Particles(x.copy(), v.copy(), m, p.component, p.galaxy)
