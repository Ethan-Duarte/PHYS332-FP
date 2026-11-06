"""AGAMA equilibrium-candidate disk + bulge + halo; optional compiled dependency.

Adapted from official Agama/py/example_self_consistent_model3.py. Save particle
snapshots and test in isolation: DF self-consistency does not validate finite-N runs.
"""
import numpy as np
from ...particles import Particles, combine, recenter


def generate(config):
    import agama
    d, b, h = (config[k] for k in ('disk','bulge','halo'))
    nd, nb, nh = (int(config['particles'][k]) for k in ('disk','bulge','halo'))
    rd, zd = d['scale'], d['height']
    # AGAMA's native default dimensionless units use G=1; do not call setUnits()
    # with kpc/km/s/Msun numbers here or we'll change the gravitational constant.
    agama.setRandomSeed(int(config['seed']))
    disk = agama.Density(type='Disk', mass=d['mass'], scaleRadius=rd, scaleHeight=-zd/2)
    bulge = agama.Density(type='Spheroid', mass=b['mass'], alpha=1, beta=4,
                          gamma=1, scaleRadius=b['scale'])
    halo = agama.Density(type='Spheroid', mass=h['mass'], alpha=1, beta=4,
                         gamma=0, scaleRadius=h['scale'])
    model = agama.SelfConsistentModel(
        rminSph=min(b['scale'],h['scale'])/40, rmaxSph=10*h['scale'], sizeRadialSph=40,
        lmaxAngularSph=6, RminCyl=rd/40, RmaxCyl=12*rd,
        sizeRadialCyl=25, zminCyl=zd/5, zmaxCyl=8*rd, sizeVerticalCyl=25)
    model.components = [agama.Component(density=disk, disklike=True),
                        agama.Component(density=bulge, disklike=False),
                        agama.Component(density=halo, disklike=False)]
    model.iterate()
    potential=model.potential
    def omega2(R):
        return -potential.force(R,0,0)[0]/R
    delta=rd*1e-3
    kappa=np.sqrt(rd*(omega2(rd+delta)-omega2(rd-delta))/(2*delta)+4*omega2(rd))
    Sigma=d['mass']/(2*np.pi*rd**2)*np.exp(-1)
    sigma_r_at_rd=config.get('Q',1.5)*3.36*Sigma/kappa
    sigma_r0=sigma_r_at_rd*np.exp(.5)
    disk_df_params=dict(type='QuasiIsothermal', potential=potential,
        Sigma0=d['mass']/(2*np.pi*rd**2), Rdisk=rd, Hdisk=zd/2,
        sigmar0=sigma_r0, Rsigmar=2*rd, sigmamin=.01*sigma_r0)
    df_disk=agama.DistributionFunction(**disk_df_params)
    # A quasi-isothermal DF's integrated mass needn't exactly equal Sigma0*2*pi*Rd².
    disk_df_params['Sigma0']*=d['mass']/df_disk.totalMass()
    df_disk=agama.DistributionFunction(**disk_df_params)
    df_bulge=agama.DistributionFunction(type='QuasiSpherical',potential=potential,density=bulge)
    df_halo=agama.DistributionFunction(type='QuasiSpherical',potential=potential,density=halo)
    model.components = [
        agama.Component(df=df_disk,disklike=True,rminCyl=rd/40,rmaxCyl=12*rd,
                        sizeRadialCyl=25,zminCyl=zd/5,zmaxCyl=8*rd,sizeVerticalCyl=25),
        agama.Component(df=df_bulge,disklike=False,rminSph=b['scale']/40,
                        rmaxSph=10*h['scale'],sizeRadialSph=40),
        agama.Component(df=df_halo,disklike=False,rminSph=h['scale']/40,
                        rmaxSph=10*h['scale'],sizeRadialSph=40)]
    for _ in range(4): model.iterate()
    def sample(df,n,i):
        xv,m=agama.GalaxyModel(potential=model.potential,df=df,af=model.af).sample(n)
        return Particles(xv[:,:3],xv[:,3:],m,np.full(n,i),np.zeros(n,dtype=int))
    return recenter(combine(sample(df_disk,nd,0),sample(df_bulge,nb,1),sample(df_halo,nh,2)))
