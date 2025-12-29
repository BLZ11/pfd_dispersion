"""
PFD - Petersson-Frisch Dispersion Model
========================================

A Python implementation of the PFD empirical dispersion correction for
quantum chemical calculations.

Overview
--------
PFD corrects for long-range electron correlation (dispersion/van der Waals)
interactions not captured by most density functional theory (DFT) methods.
The correction consists of three components:

1. **Two-body dispersion**: Classical C6/R⁶ pairwise interactions with
   damping to avoid double-counting at short range.

2. **Three-body dispersion**: Axilrod-Teller-Muto (ATM) correction for
   non-additive dispersion effects in molecular clusters.

3. **Hydrogen bond correction**: Angular-dependent modification of
   dispersion in hydrogen-bonded systems.

Quick Start
-----------
>>> from pfd import compute_pfd, parse_gaussian, get_atomic_data
>>>
>>> # Parse Gaussian output for geometry
>>> title, n, scf, z, coords = parse_gaussian("molecule.out")
>>> ehomo, alpha = get_atomic_data(z)
>>>
>>> # Compute dispersion energy and forces
>>> results = compute_pfd(n, z, coords, ehomo, alpha)
>>> print(f"Dispersion energy: {results['energy_total']:.6f} Hartree")

Module Structure
----------------
- parameters.py: Model parameters and atomic data
- core.py: Geometry, damping, and coefficient calculations
- energy.py: Energy calculations
- gradients.py: Analytical first derivatives (forces)
- hessian.py: Analytical second derivatives
- numerical.py: Numerical derivatives for validation
- io.py: Input/output handling
- pfd_calculator.py: Command-line interface

Units
-----
- Energies: Hartree
- Distances: Bohr (1 Bohr = 0.529177 Å)
- Forces: Hartree/Bohr
- Polarizabilities: Bohr³

References
----------
Petersson et al., J. Chem. Phys. (full citation in documentation)
"""

from .parameters import PFDParams, PARAMS, get_atomic_data, HARTREE_TO_MH
from .io import parse_gaussian, write_output, PFDResult
from .energy import compute_pfd_energy, compute_2body_energy, compute_3body_energy, compute_hbond_energy
from .gradients import compute_all_gradients
from .hessian import compute_all_hessians

__all__ = ['compute_pfd', 'parse_gaussian', 'write_output', 'get_atomic_data',
           'PFDParams', 'PARAMS', 'PFDResult', 'HARTREE_TO_MH']


def compute_pfd(num_atoms: int, z, coords, ehomo, alpha,
                params: PFDParams = PARAMS,
                compute_gradients: bool = True,
                compute_hessians: bool = False,
                use_numerical_gradients: bool = False,
                use_numerical_hessians: bool = False) -> dict:
    """
    Main interface for PFD dispersion calculations.
    
    Parameters
    ----------
    num_atoms : int
        Number of atoms in the molecule.
    z : array-like
        Atomic numbers, shape (num_atoms,).
    coords : array-like
        Cartesian coordinates in Bohr, shape (3, num_atoms).
    ehomo : array-like
        HOMO energies in Hartree, shape (num_atoms,).
    alpha : array-like
        Polarizabilities in Bohr³, shape (num_atoms,).
    params : PFDParams, optional
        Model parameters. Default uses standard PFD parameters.
    compute_gradients : bool, optional
        Whether to compute forces. Default True.
    compute_hessians : bool, optional
        Whether to compute Hessians. Default False.
    use_numerical_gradients : bool, optional
        Use numerical instead of analytical gradients (for testing).
    use_numerical_hessians : bool, optional
        Use numerical instead of analytical Hessians (for testing).
        
    Returns
    -------
    dict
        Dictionary containing:
        - energy_2body, energy_3body, energy_hbond, energy_total: Energies (Hartree)
        - forces_2body, forces_3body, forces_hbond, forces_total: Forces (Hartree/Bohr)
        - hessian_2body, hessian_3body, hessian_hbond, hessian_total: Hessians
        - pair_energies: Matrix of pairwise energies
        
    Examples
    --------
    Energy only:
    >>> results = compute_pfd(n, z, coords, ehomo, alpha, compute_gradients=False)
    
    Energy + forces:
    >>> results = compute_pfd(n, z, coords, ehomo, alpha)
    
    Full calculation with Hessians:
    >>> results = compute_pfd(n, z, coords, ehomo, alpha, compute_hessians=True)
    """
    import numpy as np
    z = np.asarray(z, dtype=np.int32)
    coords = np.asarray(coords, dtype=np.float64)
    ehomo = np.asarray(ehomo, dtype=np.float64)
    alpha = np.asarray(alpha, dtype=np.float64)
    
    if use_numerical_gradients or use_numerical_hessians:
        from .numerical import compute_all_numerical
        return compute_all_numerical(num_atoms, z, coords, ehomo, alpha, params,
                                     compute_hess=use_numerical_hessians)
    
    if compute_hessians:
        h_2b, h_3b, h_hb, f_2b, f_3b, f_hb, e_2b, e_3b, e_hb, pe = compute_all_hessians(
            num_atoms, z, coords, ehomo, alpha, params)
    elif compute_gradients:
        f_2b, f_3b, f_hb, e_2b, e_3b, e_hb, pe = compute_all_gradients(
            num_atoms, z, coords, ehomo, alpha, params)
        h_2b = h_3b = h_hb = np.zeros((3, 3, num_atoms))
    else:
        e_2b, e_3b, e_hb, pe = compute_pfd_energy(num_atoms, z, coords, ehomo, alpha, params)
        f_2b = f_3b = f_hb = np.zeros((3, num_atoms))
        h_2b = h_3b = h_hb = np.zeros((3, 3, num_atoms))
    
    return {
        'energy_2body': e_2b, 'energy_3body': e_3b, 'energy_hbond': e_hb,
        'energy_total': e_2b + e_3b + e_hb,
        'forces_2body': f_2b, 'forces_3body': f_3b, 'forces_hbond': f_hb,
        'forces_total': f_2b + f_3b + f_hb,
        'hessian_2body': h_2b, 'hessian_3body': h_3b, 'hessian_hbond': h_hb,
        'hessian_total': h_2b + h_3b + h_hb,
        'pair_energies': pe,
    }
