from grale import GWEvent, calculate_distances, LensingCalculator
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np

# SETUP
cosmo = FlatLambdaCDM(H0=67.9, Om0=0.3065)
z_lens = 0.0098
distance_mu1 = 40 * u.Mpc

# STEP 1: Create GWEvent with one m1, m2 pair (required)
holy = GWEvent(m1=1.1, m2=1.5)
print("GWEvent Results:")
# Calculate chirp mass
chirp_mass = holy.chirp_mass()
print(f"Chirp Mass: {chirp_mass:.3f} M_sun")

# Get M0 also by:
holy.M0 


# Get redshift results over those mass ranges 
# Delta is the mass difference 
redshift_results = holy.redshift_range(delta=0.5, step=0.01, z_lens=z_lens) 
print("-----------------------------")
print("Distance Calculator Results:")
# Setup the distance calculator
dist = calculate_distances(cosmo = cosmo)

# Calculate luminosity distance
dlum = dist.luminosity_distance(z_lens)
print(f"Luminosity Distance: {dlum:.3f} Mpc")


print("-----------------------------")
print("Lensing Calculator Results:")

# Setup lensing calculator
z_lensed = redshift_results["plausible_redshifts_lensed"] #get an array of redshifts for later use 


# STEP 2: Compute lensing quantities

# Iniitialize the LensingCalculator with the cosmology, distance, and lens parameters

biblicly_accurate = LensingCalculator(
    cosmo=cosmo,
    D_mu1=distance_mu1,
    sigma=160,
    theta_offset=10.07
)


magnification = biblicly_accurate.magnification(dlum)
print(f"Magnification of source: {magnification:.3f}")


lens_results = biblicly_accurate.compute_over_redshift_range(z_lensed, z_lens=z_lens)

print("-----------------------------")


# STEP 3: Reverse-calculate redshift range from magnification

mu_geo = lens_results["mu_geo"]
magn_range = np.linspace(mu_geo.min(), mu_geo.max(), 100)

print("Reverse Calculator Results:")
reverse_redshifts, reverse_distances = biblicly_accurate.reverse_calc(magn_range)


# PLOT: Redshift vs Magnification
plt.figure(figsize=(8, 5))
plt.scatter(lens_results["z"], lens_results["mu_geo"], color='purple')
plt.xlabel("Redshift (z)")
plt.ylabel("Magnification (μ)")
plt.title("Magnification vs Redshift")
plt.grid(True)
plt.tight_layout()
plt.show()
