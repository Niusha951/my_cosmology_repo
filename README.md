# My Cosmology Repository

This repository contains various Python scripts for cosmological calculations. Below is a brief overview of each script and how to use them.

## Scripts

### `Mvir_Rvir_NFWcalculation.py`

This script calculates the virial mass of halos based on the NFW profile. It reads `V_max` and `R_max` from an input file, performs calculations to solve `M_vir` and `R_vir`, and saves the results to an output file. 

**Usage:**
```bash
python halo_mass_calculation.py
```

**Customization**:
- **Input Structure**: Adjustments may be required depending on the structure of your input file.
- **Cosmology Parameters**: You can easily modify the cosmology parameters by changing the constants defined at the beginning of the script.
- **Virial Overdensity**: The relation for virial overdensity can be changed by modifying the corresponding function.

### `calculate_V_max.py`

This Python script calculates the `V_max` assuming an NFW profile for dark matter halos.

**Customization**:
- **Input Parameters**: Adjust the values of `M_vir`, `R_s`, and `R_vir` according to your specific halo parameters.
- **Cosmology Constants**: Modify `G` or `h` as needed for different cosmological models or simulations.

### `z_to_time.py`

This Python script contains functions for converting between redshift (z) and cosmic time (t) in a cosmological context, based on specified cosmological parameters. It includes utilities to calculate lookback time given cosmic time.

**Customization**:
- **Input Parameters**: Ensure the cosmological parameters (`H0`, `Om0`, `Ol0`, `Ob0`) are set correctly according to your specific cosmology.

