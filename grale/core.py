import numpy as np
from astropy import units as u
from astropy.cosmology import z_at_value

class GWEvent:
    def __init__(self, m1, m2):
        """
        Initialize a gravitational wave event with component masses m1 and m2.

        Parameters:
            m1 (float): Mass of the lighter compact object (must be > 0).
            m2 (float): Mass of the heavier compact object (must be > 0 and m2 ≥ m1).

        Raises:
            ValueError: If masses are non-positive or if m1 > m2.
        """
        if m1 <= 0 or m2 <= 0:
            raise ValueError("Masses must be positive.")
        if m1 > m2:
            raise ValueError("m1 must be smaller or equal to m2 (m1 ≤ m2).")
        self.m1 = m1
        self.m2 = m2

    def chirp_mass(self, m1=None, m2=None):
        """
        Calculate the chirp mass M₀ = (m₁ * m₂)^{3/5} / (m₁ + m₂)^{1/5}, 
        which determines the amplitude and frequency evolution of a GW signal.

        Parameters:
            m1 (float, optional): First mass. Defaults to self.m1.
            m2 (float, optional): Second mass. Defaults to self.m2.

        Returns:
            float: Chirp mass value.

        Notes:
            The chirp mass is stored in self.M0 as a side effect.
        """
        m1 = m1 if m1 is not None else self.m1
        m2 = m2 if m2 is not None else self.m2
        self.M0 = (m1 * m2)**(3/5) / (m1 + m2)**(1/5)
        return (m1 * m2)**(3/5) / (m1 + m2)**(1/5)
    
    def true_redshift(self, M):
        """
        Compute the true redshift of a gravitational wave source,
        given the observed (redshifted) chirp mass and the intrinsic chirp mass.

        Parameters:
            M (float): Observed redshifted chirp mass.

        Returns:
            float: Estimated redshift (z = M₀ / M - 1).
        """
        return (self.M0 / M) - 1

    
    def redshift_range(self, 
                       delta=0.4, 
                       step=0.01, 
                       m1_range=None, 
                       m2_range=None, 
                       z_lens=None):
        """
        Compute a grid of chirp masses and their corresponding redshifts over a 
        parameter space of m1 and m2. Optionally apply a redshift threshold to model lensing.

        Parameters:
            delta (float): Variation range ± around self.m1 and self.m2 if no ranges provided.
            step (float): Step size for mass grid.
            m1_range (np.ndarray, optional): Custom range of m1 values.
            m2_range (np.ndarray, optional): Custom range of m2 values.
            z_lens (float, optional): Minimum redshift to qualify as lensed.

        Returns:
            dict: Contains mass ranges, chirp mass grid, redshift grid, and filters:
                - 'm1_range', 'm2_range'
                - 'chirp_masses'
                - 'redshifts'
                - 'plausible_redshifts' (z ≥ 0)
                - 'plausible_redshifts_lensed' (z ≥ z_lens, if given)
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
            [self.true_redshift(M) for M in row]
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
        """
        Initialize lensing calculator with observational and lens model parameters.

        Parameters:
            cosmo (Cosmology): Astropy cosmology instance.
            D_mu1 (Quantity): Observed luminosity distance (with units).
            sigma (float): Velocity dispersion of lens (in km/s).
            theta_offset (float): Angular offset between image and lens center (in arcseconds).

        Raises:
            ValueError: If any parameter is non-positive.
        """
        if D_mu1 <= 0 * u.Mpc or sigma <= 0 or theta_offset <= 0:
            raise ValueError("D_mu1, sigma, and theta_offset must be positive.")
        self.cosmo = cosmo
        self.D_mu1 = D_mu1
        self.sigma = sigma
        self.theta_offset = theta_offset  # in arcseconds
        self.c = 3e5  # speed of light in km/s

    def luminosity_distance(self, z):
        """Return luminosity distance (with units) for a given redshift z."""

        return self.cosmo.luminosity_distance(z)

    def comoving_distance(self, z):
        """Return comoving distance (with units) for a given redshift z."""

        return self.cosmo.comoving_distance(z)

    def angular_diameter_distance(self, z):
        """Return angular diameter distance (with units) for a given redshift z."""
        return self.cosmo.angular_diameter_distance(z)

    def angular_diameter_distance_z1z2(self, z1, z2):
        """Return angular diameter distance between redshift z1 and z2."""
        return self.cosmo.angular_diameter_distance_z1z2(z1, z2)

    def comoving_distance_diff(self, z_source, z_lens):
        """
        Compute (D_C(z_source) - D_C(z_lens)) / (1 + z_source).

        Parameters:
            z_source (float): Source redshift.
            z_lens (float): Lens redshift.

        Returns:
            Quantity: Distance value in Mpc.
        """
        return (self.comoving_distance(z_source) - self.comoving_distance(z_lens)) / (1 + z_source)

    def magnification(self, D_true):
        """
        Compute lensing magnification μ = (D_true / D_mu1)^2.

        Parameters:
            D_true (Quantity): True luminosity distance (must be positive).

        Returns:
            float: Lensing magnification factor.
        """
        if D_true <= 0 * u.Mpc:
            raise ValueError("True distance must be positive.")
        return (D_true / self.D_mu1)**2

    def einstein_radius(self, DLS, DS):
        """
        Compute Einstein radius in arcseconds.

        Parameters:
            DLS (Quantity): Angular diameter distance between lens and source.
            DS (Quantity): Angular diameter distance to source.

        Returns:
            Quantity: Einstein radius in arcseconds.
        """
        rE_rad = 4 * np.pi * (self.sigma / self.c)**2 * (DLS / DS)  # dimensionless
        rE = (rE_rad * u.rad).to(u.arcsec)  # manually add rad unit, then convert
        return rE


    def magnifying_power(self, einstein_radius):
        """
        Compute geometric magnification μ_geo = θ / (θ - θ_E).

        Parameters:
            einstein_radius (Quantity): Einstein radius in arcseconds.

        Returns:
            float: Geometric magnification factor.

        Raises:
            ValueError: If θ ≈ θ_E (unphysical case of infinite magnification).
        """
        theta = self.theta_offset * u.arcsec
        if np.isclose(theta.value, einstein_radius.value):
            raise ValueError("Attention: unphysical lensing scenario (infinite magnification).")
        return (theta / (theta - einstein_radius)).decompose().value

    def compute_over_redshift_range(self, z_array, z_lens):
        """
        Compute lensing quantities (D_S, D_LS, θ_E, μ_geo) over a source redshift array.

        Parameters:
            z_array (array-like): Source redshift values.
            z_lens (float): Lens redshift.

        Returns:
            dict: Computed lensing properties for each redshift:
                - 'z', 'DS', 'DLS', 'r_E', 'mu_geo'
                - 'plausible_magnifications': μ_geo for z ≥ z_lens
        """
        DS = np.array([self.angular_diameter_distance(z).value for z in z_array])
        DLS = np.array([self.angular_diameter_distance_z1z2(z_lens, z).value for z in z_array])
        r_E = np.array([
            self.einstein_radius(dls * u.Mpc, ds * u.Mpc).value
            for dls, ds in zip(DLS, DS)
        ])
        mu_geo = np.array([self.magnifying_power(re * u.arcsec) for re in r_E])

        # Redshift filters
        plausible_z = z_array[z_array >= z_lens] if z_lens is not None else z_array
        plausible_mu = mu_geo[z_array >= z_lens] if z_lens is not None else mu_geo
      
        # Info output (matching redshift_range style)
        print(f"z ∈ [{z_array.min():.4f}, {z_array.max():.4f}]")
        print(f"DS range: [{DS.min():.2f}, {DS.max():.2f}] Mpc")
        print(f"DLS range: [{DLS.min():.2f}, {DLS.max():.2f}] Mpc")
        print(f"Einstein radius range: [{r_E.min():.3f}, {r_E.max():.3f}] arcsec")
        print(f"Magnification range: [{mu_geo.min():.2f}, {mu_geo.max():.2f}]")

        return {
            "z": z_array,
            "DS": DS,
            "DLS": DLS,
            "r_E": r_E,
            "mu_geo": mu_geo,
            "plausible_magnifications": plausible_mu
        }

    
    def reverse_calc(self, magn_range):
        """
        Estimate the redshifts and true luminosity distances corresponding to a 
        given range of magnifications, assuming the observed distance is D_mu1.

        Parameters:
            magn_range (array-like): Array of magnification values μ > 0.

        Returns:
            tuple:
                - reverse_redshifts (list): Redshifts corresponding to true distances.
                - reverse_distances (list): True luminosity distances.
        """
        reverse_distances = self.D_mu1 * np.sqrt(magn_range)
        reverse_redshifts = [z_at_value(self.cosmo.luminosity_distance, d) for d in reverse_distances]

        mu_min, mu_max = np.min(magn_range), np.max(magn_range)
        z_min, z_max = np.min(reverse_redshifts), np.max(reverse_redshifts)
        d_min, d_max = np.min(reverse_distances).value, np.max(reverse_distances).value

        # Info output 
        print(f"Magnification range: [{mu_min:.2f}, {mu_max:.2f}]")
        print(f"Reverse calc redshift range: [{z_min:.4f}, {z_max:.4f}]")
        print(f"Reverse calc distance range: [{d_min:.2f}, {d_max:.2f}] Mpc")

        return reverse_redshifts, reverse_distances