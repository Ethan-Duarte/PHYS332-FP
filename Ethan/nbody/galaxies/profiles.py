"""Analytic galaxy mass profiles; total masses assume untruncated models."""
import numpy as np


def dehnen_mass(r, mass, a, gamma=0):
    return mass * (r / (r + a))**(3-gamma)


def dehnen_density(r, mass, a, gamma=0):
    return (3-gamma)*mass*a / (4*np.pi*r**gamma*(r+a)**(4-gamma))


def disk_mass(R, mass, rd):
    x = R / rd
    return mass * (1 - (1+x)*np.exp(-x))


def disk_surface_density(R, mass, rd):
    return mass / (2*np.pi*rd**2) * np.exp(-R/rd)
