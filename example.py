from grale import GWEvent, LensingCalculator
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np



# SETUP
# -------------------------
cosmo = FlatLambdaCDM(H0=67.9, Om0=0.3065)
z_lens = 0.0098
distance_mu1 = 40 * u.Mpc

# -------------------------
# STEP 1: Create GWEvent with one m1, m2 pair (required)
# -------------------------
event = GWEvent(m1=1.1, m2=2.3, M0=1.186)

# Get redshift results over those mass ranges 
# Delta is the mass difference 
redshift_results = event.redshift_range(delta=0.5, step=0.01, z_lens=z_lens) 

z_lensed = redshift_results["plausible_redshifts_lensed"]

# -------------------------
# STEP 2: Compute lensing quantities
# -------------------------
lens = LensingCalculator(
    cosmo=cosmo,
    D_mu1=distance_mu1,
    sigma=160,
    theta_offset=10.07
)

lens_results, lens_summary = lens.compute_over_redshift_range(z_lensed, z_lens=z_lens)

# -------------------------
# STEP 3: Reverse-calculate redshift range from magnification
# -------------------------
mu_geo = lens_results["mu_geo"]
magn_range = np.linspace(mu_geo.min(), mu_geo.max(), 100)

reverse_z, reverse_D, reverse_summary = lens.reverse_calc(magn_range)

# -------------------------
# Print Result summaries
# -------------------------
print("=== LENSING RESULTS ===")
print(lens_summary)

print("\n=== REVERSE CALCULATION RESULTS ===")
print(reverse_summary)

# -------------------------
# EXAMPLES OF SINGLE PARAMETERS
# -------------------------
print("\n=== EXAMPLES ===")
print("First lensed redshift:", z_lensed[0])
print("First Einstein radius:", lens_results["r_E"][0])
print("First reverse redshift:", reverse_z[0])
print("First reverse distance:", reverse_D[0])


# -------------------------
# PLOT: Redshift vs Magnification
# -------------------------
plt.figure(figsize=(8, 5))
plt.scatter(lens_results["z"], lens_results["mu_geo"], color='purple')
plt.xlabel("Redshift (z)")
plt.ylabel("Magnification (μ)")
plt.title("Magnification vs Redshift")
plt.grid(True)
plt.tight_layout()
plt.show()
