import numpy as np
from astropy import units as u
from astropy.cosmology import z_at_value

class GWEvent:
    def __init__(self, m1, m2, M0):
        self.m1 = m1
        self.m2 = m2
        self.M0 = M0

    def chirp_mass(self):
        return (self.m1 * self.m2)**(3/5) / (self.m1 + self.m2)**(1/5)

    def true_redshift(self):
        M_true = self.chirp_mass()
        return (self.M0 / M_true) - 1


class LensingCalculator:
    def __init__(self, cosmo, D_mu1, sigma, theta_offset):
        self.cosmo = cosmo
        self.D_mu1 = D_mu1
        self.sigma = sigma
        self.theta_offset = theta_offset  # in arcseconds
        self.c = 3e5  # speed of light in km/s

    def luminosity_distance(self, z):
        return self.cosmo.luminosity_distance(z)

    def comoving_distance(self, z):
        return self.cosmo.comoving_distance(z)

    def angular_diameter_distance(self, z):
        return self.cosmo.angular_diameter_distance(z)

    def angular_diameter_distance_z1z2(self, z1, z2):
        return self.cosmo.angular_diameter_distance_z1z2(z1, z2)

    def comoving_distance_diff(self, z_source, z_lens):
        return (self.comoving_distance(z_source) - self.comoving_distance(z_lens)) / (1 + z_source)

    def magnification(self, D_true):
        return (D_true / self.D_mu1)**2

    def einstein_radius(self, DLS, DS):
        rE = 4 * np.pi * (self.sigma / self.c)**2 * DLS / DS
        rE = rE.to(u.rad).to(u.arcsec)
        return rE

    def magnifying_power(self, einstein_radius):
        theta = self.theta_offset * u.arcsec
        return theta / (theta - einstein_radius)

    def reverse_calc(self, magn_range):
        """
        Given a range of magnifications, compute the redshift values that would 
        correspond to those magnifications assuming observed D_mu1.
        """
        reverse_distances = self.D_mu1 * np.sqrt(magn_range)
        reverse_redshifts = [z_at_value(self.cosmo.luminosity_distance, d) for d in reverse_distances]
        return reverse_redshifts, reverse_distances
