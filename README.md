# PFD Dispersion Calculator

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![NumPy](https://img.shields.io/badge/numpy-%23013243.svg?logo=numpy&logoColor=white)](https://numpy.org/)
[![Tests](https://img.shields.io/badge/tests-passing-brightgreen.svg)]()

A Python implementation of the Petersson-Frisch Dispersion (PFD) model for computing two-body and three-body dispersion energies with hydrogen bond corrections, as described in the PFD-3B exchange-correlation functional.

## Features

- **Two-body (C₆) dispersion** with element-specific damping functions
- **Three-body (Axilrod-Teller-Muto)** dispersion corrections
- **Hydrogen bond corrections** with angular switching functions
- **Analytical gradients** for geometry optimizations
- **Analytical Hessians** for frequency calculations
- **Numerical derivatives** for validation
- **Gaussian output parsing** for integration with existing workflows
- Pure Python with NumPy — no compilation required

## Installation

```bash
git clone https://github.com/BLZ11/pfd_dispersion.git
cd pfd_dispersion
pip install numpy
```

**Requirements:** Python 3.8+ and NumPy

## Quick Start

### Command Line

```bash
# Full calculation (energy + gradient + Hessian)
python pfd_calculator.py molecule.out results.pfd

# Energy only (fastest)
python pfd_calculator.py molecule.out results.pfd --energy-only

# Energy + gradient (no Hessian)
python pfd_calculator.py molecule.out results.pfd --grad-only

# Validate analytical gradients
python pfd_calculator.py molecule.out results.pfd --numerical-grad

# Validate analytical Hessians
python pfd_calculator.py molecule.out results.pfd --numerical-hessian
```

### Python API

```python
from pfd import compute_pfd, parse_gaussian, get_atomic_data, HARTREE_TO_MH

# Parse Gaussian output
title, n, scf_energy, z, coords = parse_gaussian('molecule.out')
ehomo, alpha = get_atomic_data(z)

# Compute PFD correction
result = compute_pfd(n, z, coords, ehomo, alpha)

# Print energies in mEₕ
print(f"2-body:  {result['energy_2body'] * HARTREE_TO_MH:12.6f} mEₕ")
print(f"3-body:  {result['energy_3body'] * HARTREE_TO_MH:12.6f} mEₕ")
print(f"H-bond:  {result['energy_hbond'] * HARTREE_TO_MH:12.6f} mEₕ")
print(f"Total:   {result['energy_total'] * HARTREE_TO_MH:12.6f} mEₕ")

# Access forces and Hessians
forces = result['forces_total']      # Shape: (3, n), units: Eₕ/Bohr
hessian = result['hessian_total']    # Shape: (3, 3, n), units: Eₕ/Bohr²
```

### Skip Hessians for Faster Calculations

```python
# For geometry optimizations (no Hessians needed)
result = compute_pfd(n, z, coords, ehomo, alpha, compute_hessians=False)
```

## Output Dictionary

The `compute_pfd()` function returns:

| Key | Shape | Units | Description |
|-----|-------|-------|-------------|
| `energy_2body` | scalar | Eₕ | Two-body dispersion |
| `energy_3body` | scalar | Eₕ | Three-body dispersion |
| `energy_hbond` | scalar | Eₕ | H-bond correction |
| `energy_total` | scalar | Eₕ | Total correction |
| `forces_2body` | (3, n) | Eₕ/Bohr | Two-body forces |
| `forces_3body` | (3, n) | Eₕ/Bohr | Three-body forces |
| `forces_hbond` | (3, n) | Eₕ/Bohr | H-bond forces |
| `forces_total` | (3, n) | Eₕ/Bohr | Total forces |
| `hessian_2body` | (3, 3, n) | Eₕ/Bohr² | Two-body Hessian blocks |
| `hessian_3body` | (3, 3, n) | Eₕ/Bohr² | Three-body Hessian blocks |
| `hessian_hbond` | (3, 3, n) | Eₕ/Bohr² | H-bond Hessian blocks |
| `hessian_total` | (3, 3, n) | Eₕ/Bohr² | Total Hessian blocks |
| `pair_energies` | (n, n) | Eₕ | Pairwise energy matrix |

## Numerical Validation

Analytical derivatives can be validated against numerical finite differences:

```python
# Numerical gradients (finite differences of energy)
result = compute_pfd(n, z, coords, ehomo, alpha, use_numerical_gradients=True)

# Numerical Hessians (finite differences of analytical gradients)
result = compute_pfd(n, z, coords, ehomo, alpha, use_numerical_hessians=True)
```

Numerical Hessians are computed from analytical gradients (not from double finite differences of energy), requiring only one level of finite differencing for improved accuracy.

## Package Structure

```
pfd/
├── __init__.py       # Main API: compute_pfd()
├── parameters.py     # Physical constants, atomic data (HOMO, α)
├── core.py           # Geometry, damping functions, C₆/C₉ coefficients
├── energy.py         # Energy calculations
├── gradients.py      # Analytical first derivatives
├── hessian.py        # Analytical second derivatives
├── numerical.py      # Numerical derivatives for validation
└── io.py             # Gaussian parsing, output formatting

pfd_calculator.py     # Command-line interface
```

## Supported Elements

Parameterized for H–Rn (Z = 1–86) with element-specific HOMO energies and polarizabilities.

## Validation

Analytical gradients and Hessians have been validated against numerical finite differences with agreement to < 10⁻⁹ Eₕ/Bohr.

## Citation

If you use this code, please cite:

> A. Austin, G. A. Petersson, M. J. Frisch, F. J. Dobek, G. Scalmani, and K. Throssell,
> "A Density Functional with Spherical Atom Dispersion Terms,"
> *J. Chem. Theory Comput.* **8**, 4989–5007 (2012).
> [DOI: 10.1021/ct300778e](https://doi.org/10.1021/ct300778e)

> G. A. Petersson, M. J. Frisch, F. Dobek, and B. Zulueta,
> "Three-Body Dispersion Corrections to the Spherical Atom Model: The PFD-3B Density Functional,"
> *J. Phys. Chem. A* **124**, 10296–10311 (2020).
> [DOI: 10.1021/acs.jpca.0c08976](https://doi.org/10.1021/acs.jpca.0c08976)

## License

MIT License

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
