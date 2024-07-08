import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

# Constants
G = 4.30091e-6  # Gravitational constant in (kpc/M_sun) * (km/s)^2
H0 = 70.0  # Hubble constant in km/s/Mpc
h = 0.7
Omega_m0 = 0.2815
Omega_lambda0 = 0.7185

# Bryan & Norman (1998) virial overdensity
def delta_vir(Omega_m):
    return 18 * np.pi**2 + 82 * (Omega_m - 1) - 39 * (Omega_m - 1)**2

# Hubble parameter at redshift z
def hubble_parameter(H0, Omega_m0, Omega_lambda0, z):
    return H0 * np.sqrt(Omega_m0 * (1 + z)**3 + Omega_lambda0)

# Critical density of the Universe at redshift z
def critical_density(H0, Omega_m0, Omega_lambda0, z):
    H_z = hubble_parameter(H0, Omega_m0, Omega_lambda0, z) / 3.086e19  # Convert H_z to 1/s
    return 3 * H_z**2 / (8 * np.pi * 6.67430e-11) * ((3.086e19)**3 / (1.9884e30))  # Convert to M_sun/kpc^3

# Solve for M_vir
def solve_M_vir(M_vir_guess, r_s, rho_s, Delta_vir, rho_crit_z):
    r_vir = (3 * M_vir_guess / (4 * np.pi * Delta_vir * rho_crit_z))**(1/3)
    x_vir = r_vir / r_s
    return M_vir_guess - (4 * np.pi * rho_s * r_s**3 * (np.log(1 + x_vir) - x_vir / (1 + x_vir)))

# Load data and process
def load_and_process_data(filename):
    data = np.loadtxt(filename, skiprows=1)
    scale_factor = data[:,0]
    R_max = data[:,7] * scale_factor / h
    V_max = data[:,6]
    M_vir_initial_guess = data[:,1] / h
    return scale_factor, R_max, V_max, M_vir_initial_guess

def calculate_virial_mass(scale_factor, R_max, V_max, M_vir_initial_guess):
    z = (1 / scale_factor) - 1
    Omega_m = Omega_m0 * (1 + z)**3 / (Omega_m0 * (1 + z)**3 + Omega_lambda0)
    Delta_vir = delta_vir(Omega_m)
    rho_crit_z = critical_density(H0, Omega_m0, Omega_lambda0, z)
    r_s = R_max / 2.163
    rho_s = (V_max**2) / (0.465**2 * 4 * np.pi * G * r_s**2)
    
    M_vir_solution = np.zeros_like(M_vir_initial_guess)
    R_vir_solution = np.zeros_like(M_vir_initial_guess)
    
    for i in range(len(M_vir_initial_guess)):
        M_vir_solution[i] = fsolve(solve_M_vir, M_vir_initial_guess[i], args=(r_s[i], rho_s[i], Delta_vir[i], rho_crit_z[i]))
        R_vir_solution[i] = (3 * M_vir_solution[i] / (4 * np.pi * Delta_vir[i] * rho_crit_z[i]))**(1/3)
    
    return M_vir_solution, R_vir_solution

def save_results(filename, scale_factor, M_vir_solution, R_vir_solution, R_max):
    results = np.stack((scale_factor, M_vir_solution * h, R_vir_solution * h / scale_factor, R_max * h / scale_factor), axis=1)
    column_names = "#scale M_vir R_vir R_max"
    with open(filename, 'w') as file:
        file.write(column_names + '\n')
        np.savetxt(file, results, fmt='%.6e')

def plot_results(scale_factor, M_vir_solution, R_vir_solution, data):
    plt.figure(figsize=(10, 6))
    plt.plot(scale_factor, R_vir_solution, 'o-', alpha=0.7, label='Python Solution')
    plt.plot(scale_factor, data[:,4] * data[:,0] / h, 'o-', alpha=0.5, label='Test')
    plt.xlabel('Scale Factor')
    plt.ylabel(r'R$_{\rm vir}$ (kpc)')
    plt.title('CDM test')
    plt.legend()
    plt.grid(True)
    
    plt.figure(figsize=(10, 6))
    plt.plot(scale_factor, M_vir_solution, 'o-', alpha=0.7, label='Python Solution')
    plt.plot(scale_factor, data[:,1] / h, 'o-', alpha=0.5, label='Test')
    plt.xlabel('Scale Factor')
    plt.ylabel(r'M$_{\rm vir}$ (M$_{\rm \odot}$)')
    plt.title('CDM test')
    plt.legend()
    plt.grid(True)
    plt.show()

def main():
    input_filename = #####
    output_filename = #####
    
    scale_factor, R_max, V_max, M_vir_initial_guess = load_and_process_data(input_filename)
    M_vir_solution, R_vir_solution = calculate_virial_mass(scale_factor, R_max, V_max, M_vir_initial_guess)
    save_results(output_filename, scale_factor, M_vir_solution, R_vir_solution, R_max)
    data = np.loadtxt(input_filename, skiprows=1)
    plot_results(scale_factor, M_vir_solution, R_vir_solution, data)

if __name__ == '__main__':
    main()

