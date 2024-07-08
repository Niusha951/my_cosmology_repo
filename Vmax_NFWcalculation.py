import numpy as np

# Define constants
G = 4.30091e-6  # Gravitational constant in (kpc/M_sun) * (km/s)^2
h = 0.7

def calculate_rho_s(M_vir, r_s, r_vir):
    """
    Calculate the characteristic density rho_s for an NFW profile.

    Parameters:
    - M_vir (float): Virial mass in solar masses.
    - r_s (float): Scale radius in kpc.
    - r_vir (float): Virial radius in kpc.

    Returns:
    - rho_s (float): Characteristic density in solar masses per cubic kpc.
    """
    x_vir = r_vir / r_s
    rho_s = M_vir / (4 * np.pi * r_s**3 * (np.log(1 + x_vir) - x_vir / (1 + x_vir)))
    return rho_s

def calculate_V_max(rho_s, r_s):
    """
    Calculate the maximum circular velocity V_max assuming an NFW profile.

    Parameters:
    - rho_s (float): Characteristic density in solar masses per cubic kpc.
    - r_s (float): Scale radius in kpc.

    Returns:
    - V_max (float): Maximum circular velocity in km/s.
    """
    V_max = 0.465 * np.sqrt(4 * np.pi * G * rho_s * r_s**2)
    return V_max

# Example parameters
M_vir = 1e12  # Example virial mass in solar masses
R_s = 10.0  # Example scale radius in kpc
R_vir = 100.0  # Example virial radius in kpc

# Calculate rho_s and V_max
rho_s = calculate_rho_s(M_vir, R_s, R_vir)
V_max = calculate_V_max(rho_s, R_s)

print(f"V_max calculated assuming NFW: {V_max:.2f} km/s")

