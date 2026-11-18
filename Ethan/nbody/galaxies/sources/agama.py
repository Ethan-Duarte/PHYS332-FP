"""AGAMA self-consistent disk+bulge+halo equilibrium *candidate*.

Ported against the supplied AGAMA-master source, especially
py/example_self_consistent_model3.py and doc/reference.tex.
G=1 throughout: do not call agama.setUnits().
"""
import numpy as np
from ...particles import Particles, combine, recenter


def generate(config, with_report=False, realizations=1):
    import agama
    d, b, h = (config[k] for k in ('disk', 'bulge', 'halo'))
    nd, nb, nh = (int(config['particles'][k]) for k in ('disk', 'bulge', 'halo'))
    rd, zd = float(d['scale']), float(d['height'])
    seed, Q = int(config['seed']), float(config.get('Q', 1.5))
    agama.setRandomSeed(seed)

    # In AGAMA the negative scaleHeight indicates sech²(z/(2h)); h = zd/2.
    disk = agama.Density(type='Disk', mass=d['mass'], scaleRadius=rd, scaleHeight=-zd/2)
    bulge = agama.Density(type='Spheroid', mass=b['mass'], alpha=1, beta=4,
                          gamma=1, scaleRadius=b['scale'])
    halo = agama.Density(type='Spheroid', mass=h['mass'], alpha=1, beta=4,
                         gamma=0, scaleRadius=h['scale'])

    # The old adapter used rmaxSph=10*a_h=1.0 and failed to encompass
    # the extended halo. This grid now extends beyond the prescribed halo
    # outer scale/cutoff (the cutoff itself is preview-only, not an AGAMA taper).
    rmax = max(2*float(h.get('cutoff', 0)), 30*float(h['scale']),
               15*float(b['scale']), 15*rd)
    rmin = min(float(b['scale']), float(h['scale']), rd) / 40
    sph_grid = dict(rminSph=rmin, rmaxSph=rmax,
                    sizeRadialSph=50, lmaxAngularSph=6)
    disk_grid = dict(RminCyl=rd/30, RmaxCyl=12*rd,
                     sizeRadialCyl=30, zminCyl=zd/4,
                     zmaxCyl=10*rd, sizeVerticalCyl=26)
    model = agama.SelfConsistentModel(**sph_grid, **disk_grid)
    model.components = [agama.Component(density=disk, disklike=True),
                        agama.Component(density=bulge, disklike=False),
                        agama.Component(density=halo, disklike=False)]
    model.iterate()
    potential = model.potential

    # Q = kappa sigma_R / (3.36 G Sigma), at R=Rd, with G=1.
    def omega2(R):
        return -potential.force(R, 0., 0.)[0] / R
    step = rd * 1e-3
    kappa2 = rd * (omega2(rd+step) - omega2(rd-step)) / (2*step) + 4*omega2(rd)
    kappa = float(np.sqrt(kappa2))
    surface_density_rd = float(d['mass']) / (2*np.pi*rd**2*np.e)
    sigmaR_rd = Q * 3.36 * surface_density_rd / kappa
    df_par = dict(type='QuasiIsothermal', potential=potential,
                  Sigma0=d['mass']/(2*np.pi*rd**2), Rdisk=rd, Hdisk=zd/2,
                  sigmar0=sigmaR_rd*np.exp(.5), Rsigmar=2*rd,
                  sigmamin=0.01*sigmaR_rd)
    disk_df = agama.DistributionFunction(**df_par)
    df_par['Sigma0'] *= float(d['mass']) / disk_df.totalMass()
    disk_df = agama.DistributionFunction(**df_par)
    bulge_df = agama.DistributionFunction(type='QuasiSpherical', potential=potential, density=bulge)
    halo_df = agama.DistributionFunction(type='QuasiSpherical', potential=potential, density=halo)
    dfs = [disk_df, bulge_df, halo_df]

    # Iterate in the combined, live potential rather than sampling the first guess.
    model.components = [
        agama.Component(df=disk_df, disklike=True, **disk_grid),
        agama.Component(df=bulge_df, disklike=False, **sph_grid),
        agama.Component(df=halo_df, disklike=False, **sph_grid)]
    history = []
    for _ in range(int(config.get('agama_iterations', 4))):
        model.iterate()
        history.append([float(c.density.totalMass()) for c in model.components])

    def sample(df, n, component):
        xv, masses = agama.GalaxyModel(potential=model.potential, df=df, af=model.af).sample(n)
        return Particles(np.asarray(xv[:, :3], float), np.asarray(xv[:, 3:], float),
                         np.asarray(masses, float), np.full(n, component, dtype=int),
                         np.zeros(n, dtype=int))
    realizations_out, sampled_masses = [], []
    for realization in range(realizations):
        agama.setRandomSeed(seed + realization)
        components = [sample(df, n, i) for i, (df, n) in enumerate(zip(dfs, (nd, nb, nh)))]
        sampled_masses.append([float(x.m.sum()) for x in components])
        realizations_out.append(recenter(combine(*components)))
    report = dict(
        generator='AGAMA QuasiIsothermal disk + QuasiSpherical bulge/halo',
        source='AGAMA-master/py/example_self_consistent_model3.py',
        status='DF candidate: isolated evolution and convergence NOT YET VALIDATED',
        units='dimensionless, G=1', seed=seed, Q_target_at_Rd=Q,
        kappa_at_Rd=kappa, sigmaR_target_at_Rd=float(sigmaR_rd),
        model_radius_max=float(rmax), iterations=len(history),
        requested_masses=[float(d['mass']), float(b['mass']), float(h['mass'])],
        initial_density_masses=[float(x.totalMass()) for x in (disk,bulge,halo)],
        distribution_function_masses=[float(x.totalMass()) for x in dfs],
        iteration_density_masses=history,
        sampled_masses=sampled_masses,
        n_particles=[nd, nb, nh], n_realizations=realizations,
        notes='AGAMA result may differ from initial analytic density, and native Plummer softening changes equilibrium.')
    output = realizations_out[0] if realizations == 1 else tuple(realizations_out)
    return (output, report) if with_report else output
