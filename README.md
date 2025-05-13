# GRaLE Gravitational-wave Lensing Explorer
This package investigates whether gravitational wave events, such as GW170817 may have been gravitationally lensed by a foreground galaxy, by computing lensing-related quantities.

---

##  Requirements

- Python 3.8+
- numpy
- matplotlib
- astropy

---

## What Problem the Project Aims to Solve

Gravitational wave signals, such as GW170817, may have been affected by gravitational lensing due to a nearby foreground galaxy (e.g., NGC 4993). Lensing could alter the apparent luminosity distance and energy signature of the event. This package aims to quantify and explore such effects and determine if the observed properties of GW170817 are consistent with a lensing scenario.

---

##  Expected Functionality

- Calculate chirp mass from neutron star masses
- Estimate true redshift from observed chirp mass
- Compute luminosity distance from redshift using cosmology
- Estimate magnification from true and observed distance
- Compute Einstein radius from lens parameters
- Estimate magnifying power of a galaxy
- Provide visualizations for key lensing diagnostics

---

## Why Each Piece of Functionality Is Needed

| Function                     | Purpose                                                                        |
|-----------------------------|---------------------------------------------------------------------------------|
| `chirp_mass()`              | Calculates intrinsic property of the binary system                              |
| `true_redshift()`           | Infers redshift assuming lensing magnified the signal                           |
| `luminosity_distance()`     | Converts redshift to distance using cosmology                                   |
| `magnification()`           | Determines how much brighter the event appears due to lensing                   |
| `einstein_radius()`         | Evaluates lens strength from velocity dispersion and geometry                   |
| `magnifying_power()`        | Calculates total lensing effect from galaxy-lens configuration                  |

---

## Expected Interfaces / Dependencies

GRaLE is a pure Python library that works in:

- Jupyter notebooks
- Python scripts
- Scientific workflows or simulations

**Dependencies:**

- `numpy`
- `matplotlib`
- `astropy`

---

## 🧑‍🚀 User Stories

###  Main Functionality

**Story 1:**  
_As an astrophysics student, I want to calculate the true redshift of GW170817 using different combinations of neutron star masses so that I can explore whether lensing might explain the observed chirp mass._

**Story 2:**  
_As a researcher, I want to compute the Einstein radius and magnifying power of a potential lensing galaxy, so I can determine whether lensing is physically plausible in the case of GW170817._

---

### ⚠ Edge Cases

**Story 3:**  
_As a user, I want the package to raise warnings or errors when unphysical values (e.g. negative redshifts or infinite magnifications) occur, so that I can avoid incorrect conclusions._

**Story 4:**  
_As a user, I want the lensing calculator to handle very small angular offsets and low velocity dispersions without crashing, so I can model weak lensing scenarios reliably._

---

## Pseudocode Usage Examples

```python
from grale import GWEvent, LensingCalculator
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

# SETUP
cosmo = FlatLambdaCDM(H0=67.9, Om0=0.3065)
z_lens = 0.0098
distance_mu1 = 40 * u.Mpc

# STEP 1: Create GWEvent with one m1, m2 pair (required)
event = GWEvent(m1=1.1, m2=1.5)

# Calculate chirp mass
chirp_mass = event.chirp_mass()
print(f"Chirp Mass: {chirp_mass:.3f} M_sun")

# Get M0 :
event.M0 


# Get redshift results over those mass ranges 
# Delta is the mass difference 
redshift_results = event.redshift_range(delta=0.5, step=0.01, z_lens=z_lens) 

z_lensed = redshift_results["plausible_redshifts_lensed"] #get an array of redshifts for later use
```

---
## Link to ReadTheDocs: 
https://grale.readthedocs.io/en/main/index.html
---

##  License

This project is licensed under the MIT License. See the `LICENSE` file for more information.

---

## Author

Shannon Schroeder  
