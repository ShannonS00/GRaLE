import numpy as np
from astropy import units as u
from astropy.cosmology import z_at_value

class GWEvent:
    def __init__(self, m1, m2, M0):
        if m1 <= 0 or m2 <= 0 or M0 <= 0:
            raise ValueError("Masses and observed chirp mass must be positive.")
        if m1 > m2:
            raise ValueError("m1 must be smaller or equal than m2 (m1 < m2).")
        self.m1 = m1
        self.m2 = m2
        self.M0 = M0

    def chirp_mass(self, m1=None, m2=None):
        """
        Calculate chirp mass for given m1 and m2 (or use self.m1 and self.m2).
        """
        m1 = m1 if m1 is not None else self.m1
        m2 = m2 if m2 is not None else self.m2
        return (m1 * m2)**(3/5) / (m1 + m2)**(1/5)
    
    def true_redshift_single(self, M):
        """
        Calculate redshift for a single chirp mass value.
        """
        return (self.M0 / M) - 1

    
    def true_redshift(self):
        """
        Calculate true redshift from the default chirp mass.
        """
        M_true = self.chirp_mass()
        return self.true_redshift_single(M_true)
    
    def redshift_range(self, 
                       delta=0.4, 
                       step=0.01, 
                       m1_range=None, 
                       m2_range=None, 
                       z_lens=None):
        """
        Compute a grid of redshifts over m1/m2 combinations.
        
        You can:
        - Specify m1_range and m2_range manually (numpy arrays)
        - OR vary around self.m1/self.m2 ± delta with custom step
        - Optionally provide z_lens to filter lensed scenarios
        """
        # If no mass ranges are given, construct them from delta and step
        if m1_range is None:
            m1_range = np.arange(self.m1 - delta, self.m1 + delta + step, step)
        if m2_range is None:
            m2_range = np.arange(self.m2 - delta, self.m2 + delta + step, step)

        # Grid calculation
        chirp_masses = np.array([
            [self.chirp_mass(m1, m2) for m1 in m1_range]
            for m2 in m2_range
        ])
        redshifts = np.array([
            [self.true_redshift_single(M) for M in row]
            for row in chirp_masses
        ])

        # Filters
        plausible = redshifts[redshifts >= 0]
        plausible_lensed = redshifts[redshifts >= z_lens] if z_lens is not None else []

        # Info output
        print(f"m1 ∈ [{m1_range.min():.2f}, {m1_range.max():.2f}]")
        print(f"m2 ∈ [{m2_range.min():.2f}, {m2_range.max():.2f}]")
        print(f"Chirp mass range: [{chirp_masses.min():.2f}, {chirp_masses.max():.2f}]")
        print(f"Plausible redshift range: [{plausible.min():.4f}, {plausible.max():.4f}]")

        if z_lens is not None and len(plausible_lensed) > 0:
            print(f"Lensed redshift range (z ≥ {z_lens}): "
                  f"[{plausible_lensed.min():.4f}, {plausible_lensed.max():.4f}]")

        return {
            "m1_range": m1_range,
            "m2_range": m2_range,
            "chirp_masses": chirp_masses,
            "redshifts": redshifts,
            "plausible_redshifts": plausible,
            "plausible_redshifts_lensed": plausible_lensed
        }



class LensingCalculator:
    def __init__(self, cosmo, D_mu1, sigma, theta_offset):
        if D_mu1 <= 0 * u.Mpc or sigma <= 0 or theta_offset <= 0:
            raise ValueError("D_mu1, sigma, and theta_offset must be positive.")
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
        if D_true <= 0 * u.Mpc:
            raise ValueError("True distance must be positive.")
        return (D_true / self.D_mu1)**2

    def einstein_radius(self, DLS, DS):
        """
        Compute Einstein radius in arcseconds.
        DLS and DS must be Quantities with units (e.g., Mpc).
        """
        rE_rad = 4 * np.pi * (self.sigma / self.c)**2 * (DLS / DS)  # dimensionless
        rE = (rE_rad * u.rad).to(u.arcsec)  # manually add rad unit, then convert
        return rE


    def magnifying_power(self, einstein_radius):
        theta = self.theta_offset * u.arcsec
        if np.isclose(theta.value, einstein_radius.value):
            raise ValueError("Attention: unphysical lensing scenario (infinite magnification).")
        return (theta / (theta - einstein_radius)).decompose().value

    def compute_over_redshift_range(self, z_array, z_lens, return_summary=True):
        """
        Compute DS, DLS, Einstein radius, and magnifying power over redshift array.

        Parameters:
            z_array: array-like of source redshifts
            z_lens: float, redshift of the lens
            return_summary: if True, also returns a string summary for printing

        Returns:
            results: dict of lensing quantities
            summary (optional): string with ranges
        """
        DS = np.array([self.angular_diameter_distance(z).value for z in z_array])
        DLS = np.array([self.angular_diameter_distance_z1z2(z_lens, z).value for z in z_array])
        r_E = np.array([
            self.einstein_radius(dls * u.Mpc, ds * u.Mpc).value
            for dls, ds in zip(DLS, DS)
        ])
        mu_geo = np.array([self.magnifying_power(re * u.arcsec) for re in r_E])

        z_min, z_max = np.min(z_array), np.max(z_array)
        DS_min, DS_max = np.min(DS), np.max(DS)
        DLS_min, DLS_max = np.min(DLS), np.max(DLS)
        rE_min, rE_max = np.min(r_E), np.max(r_E)
        mu_min, mu_max = np.min(mu_geo), np.max(mu_geo)

        summary = (
            f"Lensed redshift range (z ≥ {z_lens}): [{z_min:.4f}, {z_max:.4f}]\n"
            f"DS range: [{DS_min:.2f}, {DS_max:.2f}] Mpc\n"
            f"DLS range: [{DLS_min:.2f}, {DLS_max:.2f}] Mpc\n"
            f"Einstein radius range: [{rE_min:.3f}, {rE_max:.3f}] arcsec\n"
            f"Magnification range: [{mu_min:.2f}, {mu_max:.2f}]"
        )

        results = {
            "z": z_array,
            "DS": DS,
            "DLS": DLS,
            "r_E": r_E,
            "mu_geo": mu_geo
        }

        return (results, summary) if return_summary else results

    
    def reverse_calc(self, magn_range):
        """
        Given a range of magnifications, compute the redshift values and 
        corresponding true luminosity distances assuming observed D_mu1.

        Always returns a summary string describing the result ranges.

        Parameters:
            magn_range: array-like of magnification values (μ)

        Returns:
            redshifts: list of redshifts
            distances: list of luminosity distances
            summary: formatted string with min/max of inputs and outputs
        """
        reverse_distances = self.D_mu1 * np.sqrt(magn_range)
        reverse_redshifts = [z_at_value(self.cosmo.luminosity_distance, d) for d in reverse_distances]

        mu_min, mu_max = np.min(magn_range), np.max(magn_range)
        z_min, z_max = np.min(reverse_redshifts), np.max(reverse_redshifts)
        d_min, d_max = np.min(reverse_distances).value, np.max(reverse_distances).value

        summary = (
            f"Magnification range: [{mu_min:.2f}, {mu_max:.2f}]\n"
            f"Reverse calc redshift range: [{z_min:.4f}, {z_max:.4f}]\n"
            f"Reverse calc distance range: [{d_min:.2f}, {d_max:.2f}] Mpc"
        )

        return reverse_redshifts, reverse_distances, summary
    