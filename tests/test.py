import pytest
import numpy as np
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from grale.core import GWEvent, LensingCalculator

def test_chirp_mass_correctness():
    event = GWEvent(m1=1.4, m2=1.2, M0=1.186)
    cm = event.chirp_mass()
    expected = (1.4 * 1.2)**(3/5) / (1.4 + 1.2)**(1/5)
    assert np.isclose(cm, expected)

def test_true_redshift_single():
    event = GWEvent(m1=1.4, m2=1.2, M0=1.186)
    z = event.true_redshift_single(1.0)
    assert np.isclose(z, 0.186)

def test_invalid_inputs():
    with pytest.raises(ValueError):
        GWEvent(m1=-1.0, m2=1.2, M0=1.186)
    with pytest.raises(ValueError):
        GWEvent(m1=1.4, m2=-1.2, M0=1.186)
    with pytest.raises(ValueError):
        GWEvent(m1=1.4, m2=1.2, M0=-1.0)
    with pytest.raises(ValueError):
        GWEvent(m1=2.0, m2=1.0, M0=1.186)

def test_redshift_range_output():
    event = GWEvent(m1=1.4, m2=1.2, M0=1.186)
    results = event.redshift_range()
    assert "chirp_masses" in results
    assert results["chirp_masses"].shape[0] > 0
    assert results["plausible_redshifts"].min() >= 0

def test_magnification_positive():
    cosmo = FlatLambdaCDM(H0=67.9, Om0=0.3065)
    lens = LensingCalculator(cosmo, D_mu1=40 * u.Mpc, sigma=160, theta_offset=10.07)
    D_true = lens.luminosity_distance(0.1)
    mu = lens.magnification(D_true)
    assert mu > 1.0

def test_einstein_radius_units():
    cosmo = FlatLambdaCDM(H0=67.9, Om0=0.3065)
    lens = LensingCalculator(cosmo, D_mu1=40 * u.Mpc, sigma=160, theta_offset=10.07)
    rE = lens.einstein_radius(80 * u.Mpc, 100 * u.Mpc)
    assert rE.unit == u.arcsec
    assert rE.value > 0

def test_compute_over_redshift_range_structure():
    cosmo = FlatLambdaCDM(H0=67.9, Om0=0.3065)
    lens = LensingCalculator(cosmo, D_mu1=40 * u.Mpc, sigma=160, theta_offset=10.07)
    z_array = np.linspace(0.01, 0.1, 5)
    results, summary = lens.compute_over_redshift_range(z_array, z_lens=0.0098)
    assert "mu_geo" in results
    assert isinstance(summary, str)

def test_reverse_calc_output_format():
    cosmo = FlatLambdaCDM(H0=67.9, Om0=0.3065)
    lens = LensingCalculator(cosmo, D_mu1=40 * u.Mpc, sigma=160, theta_offset=10.07)
    magn_range = np.linspace(1.0, 10.0, 5)
    z, D, summary = lens.reverse_calc(magn_range)
    assert len(z) == len(D)
    assert isinstance(summary, str)
