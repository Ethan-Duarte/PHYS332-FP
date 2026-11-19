"""Independent spherical integral checks; no AGAMA installation needed."""
import numpy as np
from scipy.integrate import quad
from nbody.galaxies.softening import PlummerConvolution


class SmoothSphere:
    scale = .02

    def potential(self, x):
        return -1/np.sqrt(np.sum(x*x, axis=1)+self.scale**2)

    def force(self, x):
        return -x/(np.sum(x*x, axis=1)+self.scale**2)[:, None]**1.5


def independent_spherical_integral(radius, epsilon):
    # Integrate the angular average analytically and the radial integral
    # adaptively.  This does not share the 3D implementation's quadrature.
    def integrand(t):
        s = epsilon*t/np.sqrt(1-t*t)
        angular = -2/(np.sqrt((radius+s)**2+SmoothSphere.scale**2)
                      + np.sqrt((radius-s)**2+SmoothSphere.scale**2))
        return 3*t*t*angular
    return quad(integrand, 0, 1, epsabs=1e-9, epsrel=1e-10)[0]


def test_smoothing_matches_independent_integral():
    convolution = PlummerConvolution(SmoothSphere(), .00625, (36,18,36))
    for radius in (0., .003125, .00625, .02, .08, .25):
        expected = independent_spherical_integral(radius, .00625)
        assert np.isclose(convolution([radius,0,0]), expected, rtol=5e-4)


def test_force_is_gradient_and_preserves_symmetry():
    convolution = PlummerConvolution(SmoothSphere(), .00625, (36,18,36))
    points = np.array([[.01,0,.003],[-.01,0,-.003]])
    force = convolution.force(points)
    assert np.allclose(force[0], -force[1], atol=1e-10)
    assert np.allclose(force[:,1], 0, atol=1e-10)
    for dim in (0,2):
        delta = np.zeros(3); delta[dim] = 1e-6
        derivative = -(convolution(points+delta)-convolution(points-delta))/(2e-6)
        assert np.allclose(derivative, force[:,dim], rtol=1e-6)


def test_small_softening_limit():
    source = SmoothSphere()
    convolution = PlummerConvolution(source, 1e-6)
    points = np.array([[.02,0,.005],[.08,0,.01]])
    assert np.allclose(convolution(points), source.potential(points), rtol=1e-6)
