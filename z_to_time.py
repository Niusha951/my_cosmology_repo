import numpy as np
from scipy.integrate import quad
from scipy.optimize import fsolve

# Define cosmological parameters
H0 = 70.0  # Hubble constant in km/s/Mpc
Om0 = 0.3  # Present-day matter density parameter
Ol0 = 0.7  # Present-day dark energy density parameter
Ob0 = 0.05  # Present-day baryon density parameter

# Define unit conversion constants
km_s = 1e3  # m/s
Mpc = 3.08e16  # m
Gyr = 365.25 * 24 * 60 * 60 * 1e9  # s
unit_conv = Gyr * km_s / Mpc  # Conversion factor for time units

# Calculate Hubble parameter at present time in units of 1/Gyr
H0_inv = 1 / (H0 / unit_conv)

# Define the function for the integrand
def integrand(z):
    """
    Calculate the integrand for cosmic time calculation.

    Parameters:
    - z (float): Redshift value.

    Returns:
    - float: Value of the integrand at the given redshift.
    """
    E = np.sqrt(Om0 * (1 + z)**3 + Ob0 * (1 + z)**3 + Ol0)
    return 1 / (E * (1 + z))

# Function to calculate cosmic time given a redshift
def redshift_to_time(z):
    """
    Convert redshift to cosmic time.

    Parameters:
    - z (float): Redshift value.

    Returns:
    - float: Cosmic time in Gyr corresponding to the given redshift.
    """
    integral, _ = quad(integrand, z, np.inf)
    t = H0_inv * integral  # Cosmic time in Gyr
    return t

# Function to calculate redshift given cosmic time (inverse problem)
def time_to_redshift(t):
    """
    Convert cosmic time to redshift (inverse problem).

    Parameters:
    - t (float): Cosmic time in Gyr.

    Returns:
    - float: Redshift corresponding to the given cosmic time.
    """
    if t == 0:
        return float('inf')
    
    def solve_z(z):
        integral, _ = quad(integrand, z, np.inf)
        t_guess = H0_inv * integral
        return t - t_guess
    
    # Calculate the redshift corresponding to the cosmic time
    z = fsolve(solve_z, 0)[0]
    return z

# Function to calculate lookback time given cosmic time
def look_back_time(t):
    """
    Calculate the lookback time given cosmic time.

    Parameters:
    - t (float): Cosmic time in Gyr.

    Returns:
    - float: Lookback time in Gyr corresponding to the given cosmic time.
    """
    age_Univ = redshift_to_time(0)
    look_back = age_Univ - t
    return look_back

# Example usage and demonstration
if __name__ == "__main__":
    # Generate a range of redshifts and calculate cosmic times and lookback times
    redshifts = np.linspace(0, 10, 100)
    cosmic_times = np.array([redshift_to_time(z) for z in redshifts])
    lookback_times = np.array([look_back_time(t) for t in cosmic_times])
    
    # Print results for redshift to time conversion
    print("Redshift to Cosmic Time Conversion:")
    for z, t, t_l in zip(redshifts, cosmic_times, lookback_times):
        print(f"Redshift: {z:.2f}, Cosmic Time: {t:.2f} Gyr, Lookback Time: {t_l:.2f} Gyr")

    # Generate a range of cosmic times and calculate corresponding redshifts and lookback times
    cosmic_times_range = np.linspace(0, 10, 100)
    redshifts_calculated = np.array([time_to_redshift(t) for t in cosmic_times_range])
    lookback_times = np.array([look_back_time(t) for t in cosmic_times_range])
    
    # Print results for time to redshift conversion
    print("\nCosmic Time to Redshift Conversion:")
    for t, z, t_l in zip(cosmic_times_range, redshifts_calculated, lookback_times):
        print(f"Cosmic Time: {t:.2f} Gyr, Redshift: {z:.2f}, Lookback Time: {t_l:.2f} Gyr")
