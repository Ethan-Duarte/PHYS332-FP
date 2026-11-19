"""Plummer smoothing of a smooth potential, not a velocity rescaling.

For a mass density rho, Phi_e = W_e * Phi_Newton, with
W_e(r) = 3 e^2 / [4 pi (r^2+e^2)^(5/2)].  Particles must still be
sampled from rho, NOT from W_e*rho.  Reference: Barnes (2012),
https://arxiv.org/abs/1205.2729 .

The finite quadrature and potential expansion require convergence checks.
This module does not assert that any generated live model is stable.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss


class PlummerConvolution:
    def __init__(self, potential, epsilon, order=(16, 8, 16), max_points=250000):
        if not np.isfinite(epsilon) or epsilon <= 0:
            raise ValueError('epsilon must be positive and finite')
        nr, nm, nphi = map(int, order)
        if min(nr, nm, nphi) < 4 or nphi % 2:
            raise ValueError('quadrature orders >=4; azimuth order must be even')
        # t=s/sqrt(s*s+epsilon*epsilon); radial probability is 3*t*t dt.
        t, wt = leggauss(nr)
        t, wt = (t+1)/2, wt/2
        mu, wm = leggauss(nm)
        phi = 2*np.pi*np.arange(nphi)/nphi
        tt, mm, pp = np.meshgrid(t, mu, phi, indexing='ij')
        radius = epsilon*tt/np.sqrt(1-tt*tt)
        planar = radius*np.sqrt(1-mm*mm)
        self.offset = np.stack((planar*np.cos(pp), planar*np.sin(pp), radius*mm), axis=-1).reshape(-1, 3)
        self.weight = np.broadcast_to(
            (3*t*t*wt)[:, None, None]*wm[None, :, None]/(2*nphi),
            tt.shape).reshape(-1)
        self.potential = potential
        self.order = (nr, nm, nphi)
        self.epsilon = float(epsilon)
        self.chunk = max(1, int(max_points)//len(self.weight))
        if abs(self.weight.sum()-1) > 1e-12:
            raise RuntimeError('quadrature weights do not conserve mass')

    def _evaluate(self, xyz, force=False):
        xyz = np.asarray(xyz, dtype=float)
        scalar = xyz.ndim == 1
        xyz = np.atleast_2d(xyz)
        result = np.empty((len(xyz), 3) if force else len(xyz))
        for start in range(0, len(xyz), self.chunk):
            pts = xyz[start:start+self.chunk, None, :] + self.offset[None, :, :]
            if force:
                val = self.potential.force(pts.reshape(-1, 3)).reshape(len(pts), -1, 3)
                result[start:start+len(pts)] = np.einsum('ijk,j->ik', val, self.weight)
            else:
                val = self.potential.potential(pts.reshape(-1, 3)).reshape(len(pts), -1)
                result[start:start+len(pts)] = val @ self.weight
        return result[0] if scalar else result

    def __call__(self, xyz):
        return self._evaluate(xyz)

    def force(self, xyz):
        return self._evaluate(xyz, force=True)


def soften_expansion(base, epsilon, sph_grid, disk_grid, order=(16, 8, 16), refinement=1):
    """Smooth SCM's spheroid/disk potential components separately.

    A fresh ActionFinder must be made for the returned potential.  Assigning
    scm.potential invalidates scm.af in AGAMA's supported Python interface.
    """
    import agama
    if len(base) != 2:
        raise ValueError('Expected the SCM composite with spheroid and disk expansions')
    parts = []
    for index, part in enumerate(base):
        callback = PlummerConvolution(part, epsilon, order)
        if index == 0:
            params = dict(type='Multipole', lmax=sph_grid['lmaxAngularSph'],
                          mmax=0, gridSizeR=int(60*refinement),
                          rmin=sph_grid['rminSph'], rmax=sph_grid['rmaxSph'])
        else:
            params = dict(type='CylSpline', mmax=0, gridSizeR=int(40*refinement),
                          gridSizeZ=int(40*refinement),
                          Rmin=disk_grid['RminCyl'], Rmax=disk_grid['RmaxCyl'],
                          zmin=disk_grid['zminCyl'], zmax=disk_grid['zmaxCyl'])
        # Python callbacks serialize on the GIL; avoid OpenMP thread switching.
        with agama.setNumThreads(1):
            parts.append(agama.Potential(potential=callback, symmetry='a', **params))
    return agama.Potential(*parts)
