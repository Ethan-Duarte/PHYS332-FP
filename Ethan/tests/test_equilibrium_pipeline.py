import numpy as np
import pytest
from nbody.particles import Particles
from nbody.diagnostics import structure
from nbody.galaxies.provenance import check_evolution_softening
from zeno_smooth_control import vertical_jeans_variance


def test_measured_thickness_is_invariant_under_rigid_tilt():
    rng=np.random.default_rng(17)
    x=rng.normal(size=(1200,3))*[1,.7,.02]
    v=rng.normal(size=(1200,3))*[.1,.1,.03]
    p=Particles(x,v,np.ones(len(x)),np.zeros(len(x),int),np.zeros(len(x),int))
    angle=.7
    rotation=np.array([[1,0,0],[0,np.cos(angle),-np.sin(angle)],[0,np.sin(angle),np.cos(angle)]])
    q=Particles(x@rotation.T,v@rotation.T,p.m,p.component,p.galaxy)
    a,b=structure(p),structure(q)
    for key in ['z_rms','sigma_R','sigma_z','r50']:
        assert np.isclose(a[key],b[key],rtol=1e-11)
    assert b['z_rms_fixed_axis']>10*a['z_rms']


def test_vertical_jeans_integral_recovers_exact_isothermal_sheet():
    z=np.linspace(0,14,12001)
    sigma2=.7
    restoring=2*sigma2*np.tanh(z)
    result=vertical_jeans_variance(z,restoring,1.)
    assert np.allclose(result[z<7],sigma2,rtol=3e-6)


def test_vertical_jeans_matches_combined_sheet_and_external_force():
    z=np.linspace(0,14,12001)
    # An external compressive field must add support, not be ignored.
    base=vertical_jeans_variance(z,1.4*np.tanh(z),1.)
    total=vertical_jeans_variance(z,1.4*np.tanh(z)+.2*z,1.)
    assert np.all(total[z<7]>base[z<7])


def test_force_metadata_does_not_turn_legacy_annotation_into_validation():
    note=check_evolution_softening({'epsilon':.0075},.0075)
    assert 'unrecorded' in note
    with pytest.raises(ValueError,match='declares epsilon'):
        check_evolution_softening({'epsilon':.0075},.0125)
    assert 'mismatch' in check_evolution_softening({'construction_force_model':{'kernel':'newtonian'}},.0075)
