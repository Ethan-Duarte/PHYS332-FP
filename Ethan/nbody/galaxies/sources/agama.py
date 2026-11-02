"""Optional AGAMA equilibrium backend, adapted from official self-consistent examples.

Requires `agama` (compiled extension). This adapter is not exercised by our local tests;
verify the generated model with isolated runs before scientific use.
"""
import numpy as np
from ...particles import Particles, combine, recenter


def generate(config):
    import agama
    d, b, h = (config[k] for k in ('disk', 'bulge', 'halo'))
    nd, nb, nh = (config['particles'][k] for k in ('disk', 'bulge', 'halo'))
    rd, zd = d['scale'], d['height']

    # AGAMA's negative scaleHeight means sech^2(z/(2h)); Barnes uses sech^2(z/zd).
    disk = agama.Density(type='Disk', surfaceDensity=d['mass']/(2*np.pi*rd**2),
                         scaleRadius=rd, scaleHeight=-zd/2)
    bulge = agama.Density(type='Dehnen', mass=b['mass'], gamma=1, scaleRadius=b['scale'])
    halo = agama.Density(type='Dehnen', mass=h['mass'], gamma=0, scaleRadius=h['scale'])

    # First solve the potential from the requested density components.
    model = agama.SelfConsistentModel(
        rminSph=b['scale']/40, rmaxSph=10*h['cutoff'], sizeRadialSph=40,
        lmaxAngularSph=6, RminCyl=rd/40, RmaxCyl=12*rd,
        sizeRadialCyl=25, zminCyl=zd/5, zmaxCyl=8*rd, sizeVerticalCyl=25)
    model.components = [agama.Component(density=disk, disklike=True),
                        agama.Component(density=bulge, disklike=False),
                        agama.Component(density=halo, disklike=False)]
    model.iterate()

    # Estimate the disk DF's radial dispersion from a target Q at R=Rd.
    potential = model.potential
    R = rd
    def omega2(radius):
        return -potential.force(radius, 0, 0)[0] / radius
    delta = R*1e-3
    kappa = np.sqrt(R*(omega2(R+delta)-omega2(R-delta))/(2*delta) + 4*omega2(R))
    sigma_r_at_rd = config.get('Q', 1.5)*3.36*disk.projectedDensity(R, 0)/kappa
    sigma_r0 = sigma_r_at_rd*np.exp(0.5)

    df_disk = agama.DistributionFunction(type='QuasiIsothermal', potential=potential,
        Sigma0=d['mass']/(2*np.pi*rd**2), Rdisk=rd, Hdisk=zd/2,
        sigmar0=sigma_r0, Rsigmar=2*rd, sigmamin=0.01*sigma_r0)
    df_bulge = agama.DistributionFunction(type='QuasiSpherical',
                                            potential=potential, density=bulge)
    df_halo = agama.DistributionFunction(type='QuasiSpherical',
                                           potential=potential, density=halo)
    model.components = [
        agama.Component(df=df_disk, disklike=True, rminCyl=rd/40,
                        rmaxCyl=10*rd, sizeRadialCyl=20, zminCyl=zd/5,
                        zmaxCyl=8*zd, sizeVerticalCyl=15),
        agama.Component(df=df_bulge, disklike=False,
                        rminSph=b['scale']/40, rmaxSph=20*b['scale'], sizeRadialSph=30),
        agama.Component(df=df_halo, disklike=False,
                        rminSph=h['scale']/100, rmaxSph=10*h['cutoff'], sizeRadialSph=35)]
    for _ in range(4):
        model.iterate()

    def sample(df, n, component):
        xv, masses = agama.GalaxyModel(potential=model.potential,
                                       df=df, af=model.af).sample(n)
        return Particles(xv[:, :3], xv[:, 3:], masses,
                         np.full(n, component), np.zeros(n, dtype=int))

    # AGAMA's sampler random state is separate from NumPy's seed; save the snapshot.
    return recenter(combine(sample(df_disk, nd, 0), sample(df_bulge, nb, 1),
                            sample(df_halo, nh, 2)))
