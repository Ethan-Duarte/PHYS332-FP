"""Native initial conditions. Plummer is an equilibrium DF; Mice preview is NOT."""
import numpy as np
from ...particles import Particles, combine, recenter
from ..profiles import dehnen_mass, disk_mass, disk_surface_density


def directions(n, rng):
    q = rng.normal(size=(n, 3))
    return q / np.linalg.norm(q, axis=1)[:, None]


def plummer(n, mass=1., a=1., seed=None, G=1.):
    """Exact isotropic Plummer DF: spatial CDF plus beta-distributed speeds."""
    rng = np.random.default_rng(seed)
    u = rng.random(n)
    r = a / np.sqrt(u**(-2/3) - 1)
    x = r[:, None] * directions(n, rng)
    vesc = np.sqrt(2*G*mass / np.sqrt(r*r + a*a))
    v = vesc[:, None] * np.sqrt(rng.beta(1.5, 4.5, n))[:, None] * directions(n, rng)
    return recenter(Particles(x, v, np.full(n, mass/n),
                              np.full(n, 3), np.zeros(n, dtype=int)))


def sample_dehnen(n, mass, a, gamma, rmax, rng):
    u_max = (rmax / (rmax + a))**(3-gamma)
    u = (rng.random(n) * u_max)**(1/(3-gamma))
    r = a*u/(1-u)
    return r[:, None] * directions(n, rng)


def sample_disk(n, rd, zd, rmax, rng):
    # Gamma(k=2) radial law is p(R) proportional to R exp(-R/Rd).
    R = rng.gamma(shape=2, scale=rd, size=n)
    bad = R > rmax
    while np.any(bad):
        R[bad] = rng.gamma(shape=2, scale=rd, size=bad.sum())
        bad = R > rmax
    phi = rng.uniform(0, 2*np.pi, n)
    z = zd * np.arctanh(2*rng.random(n)-1)  # sech^2(z/zd)
    return np.column_stack((R*np.cos(phi), R*np.sin(phi), z))


def mice_preview(config):
    """Illustrative model: correct requested *spatial* profiles, approximate velocities.

    NOT a self-consistent equilibrium. Use for plumbing/visualization only.
    """
    rng = np.random.default_rng(config['seed'])
    counts = config['particles']
    d, b, h = (config[k] for k in ('disk', 'bulge', 'halo'))
    nd, nb, nh = (counts[k] for k in ('disk', 'bulge', 'halo'))
    xd = sample_disk(nd, d['scale'], d['height'], d['cutoff'], rng)
    xb = sample_dehnen(nb, b['mass'], b['scale'], 1, b['cutoff'], rng)
    xh = sample_dehnen(nh, h['mass'], h['scale'], 0, h['cutoff'], rng)

    # Preliminary circular velocity: spherical approximation for all components.
    def vc(r):
        enclosed = (disk_mass(r, d['mass'], d['scale']) +
                    dehnen_mass(r, b['mass'], b['scale'], 1) +
                    dehnen_mass(r, h['mass'], h['scale'], 0))
        return np.sqrt(np.maximum(enclosed / np.maximum(r, 1e-8), 0))

    R = np.linalg.norm(xd[:, :2], axis=1)
    vd = np.column_stack((-xd[:, 1]/R, xd[:, 0]/R, np.zeros(nd))) * vc(R)[:, None]
    # Warm disk velocities (heuristic; NOT a solved Jeans/DF model).
    vd += rng.normal(size=(nd, 3)) * (0.05*vc(R))[:, None]
    rb, rh = np.linalg.norm(xb, axis=1), np.linalg.norm(xh, axis=1)
    vb = rng.normal(size=(nb, 3)) * (vc(rb)/np.sqrt(2))[:, None]
    vh = rng.normal(size=(nh, 3)) * (vc(rh)/np.sqrt(2))[:, None]

    systems = [Particles(x, v, np.full(n, spec['mass']/n), np.full(n, idx),
                         np.zeros(n, dtype=int))
               for idx, (x, v, n, spec) in enumerate(((xd, vd, nd, d),
                                                      (xb, vb, nb, b),
                                                      (xh, vh, nh, h)))]
    return recenter(combine(*systems))
