"""Galactic encounter setup: barycentric custom state or nominal parabolic orbit."""
import numpy as np
from .particles import Particles, combine


def rotation(inclination=0, node_angle=0):
    i, w = np.deg2rad([inclination, node_angle])
    ci, si, cw, sw = np.cos(i), np.sin(i), np.cos(w), np.sin(w)
    rz = np.array([[cw, -sw, 0], [sw, cw, 0], [0, 0, 1]])
    rx = np.array([[1, 0, 0], [0, ci, -si], [0, si, ci]])
    return rz @ rx


def place(p, rot=None, position=(0, 0, 0), velocity=(0, 0, 0), galaxy_id=0):
    rot = np.eye(3) if rot is None else rot
    return Particles(p.x @ rot.T + position, p.v @ rot.T + velocity,
                     p.m.copy(), p.component.copy(), np.full(len(p.m), galaxy_id))


def relative_pair(g1, g2, displacement, relative_velocity, orientations=((0, 0), (0, 0))):
    """Set r_B-r_A and v_B-v_A; place both systems at zero total COM position/momentum."""
    m1, m2 = g1.m.sum(), g2.m.sum()
    r, v = np.asarray(displacement, float), np.asarray(relative_velocity, float)
    a = place(g1, rotation(*orientations[0]), -r*m2/(m1+m2), -v*m2/(m1+m2), 0)
    b = place(g2, rotation(*orientations[1]), r*m1/(m1+m2), v*m1/(m1+m2), 1)
    return combine(a, b)


def parabolic_pair(g1, g2, separation, pericenter, orientations=((0, 0), (0, 0)), G=1):
    """Point-mass parabolic reference; extended live galaxies need not reach this pericenter."""
    m1, m2 = g1.m.sum(), g2.m.sum()
    mu = G*(m1+m2)
    vr = -np.sqrt(2*mu*(separation-pericenter)/separation**2)
    vt = np.sqrt(2*mu*pericenter)/separation
    return relative_pair(g1, g2, (separation, 0, 0), (vr, vt, 0), orientations)
