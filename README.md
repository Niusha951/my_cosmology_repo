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


