"""
PFD Numerical Module
====================

This module provides finite difference calculations for validating
analytical gradients and Hessians.

Finite difference formulas use central differences for accuracy:
- Gradient: f'(x) ≈ -[f(x+h) - f(x-h)] / (2h)
- Hessian: f''(x) ≈ -[f'(x+h) - f'(x-h)] / (2h)

The Hessian is computed by finite differences of the analytical gradient,
not by double finite differences of the energy.

The default step size (1e-5 Bohr) balances truncation error (wants larger h)
and round-off error (wants smaller h).

Usage
-----
Numerical derivatives are primarily used for testing:

>>> from pfd import compute_pfd
>>> results_analytical = compute_pfd(n, z, coords, ehomo, alpha)
>>> results_numerical = compute_pfd(n, z, coords, ehomo, alpha,
...                                  use_numerical_gradients=True)
>>> # Compare forces_total between the two

Functions
---------
- numerical_gradient: Generic gradient via finite differences
- numerical_hessian: Generic Hessian via finite differences of analytical gradient
- compute_*_numerical: Numerical versions of specific energy components
- compute_all_numerical: Full numerical calculation
"""

import numpy as np
from typing import Tuple, Callable
from .parameters import PFDParams, PARAMS

# Default finite difference step size in Bohr
STEP = np.float64(1e-5)


def numerical_gradient(energy_func: Callable, coords: np.ndarray, num_atoms: int, h: float = STEP) -> np.ndarray:
    """
    Compute gradient numerically via central finite differences.
    
    Parameters
    ----------
    energy_func : callable
        Function that takes coords and returns energy.
    coords : np.ndarray
        Atomic coordinates, shape (3, num_atoms).
    num_atoms : int
        Number of atoms.
    h : float, optional
        Step size in Bohr.
        
    Returns
    -------
    np.ndarray
        Forces = -∇E, shape (3, num_atoms).
    """
    forces = np.zeros((3, num_atoms), dtype=np.float64)
    for i in range(num_atoms):
        for a in range(3):
            cp, cm = coords.copy(), coords.copy()
            cp[a, i] += h
            cm[a, i] -= h
            forces[a, i] = -(energy_func(cp) - energy_func(cm)) / (2.0 * h)
    return forces


def numerical_hessian(grad_func: Callable, coords: np.ndarray, num_atoms: int, h: float = STEP) -> np.ndarray:
    """
    Compute Hessian numerically via central finite differences of the analytical gradient.
    
    Uses the formula:
        ∂²E/∂x_a∂x_b ≈ -[∂E/∂x_a(x_b+h) - ∂E/∂x_a(x_b-h)] / (2h)
    
    Parameters
    ----------
    grad_func : callable
        Function that takes coords and returns analytical forces array, shape (3, num_atoms).
    coords : np.ndarray
        Atomic coordinates, shape (3, num_atoms).
    num_atoms : int
        Number of atoms.
    h : float, optional
        Step size in Bohr.
        
    Returns
    -------
    np.ndarray
        Diagonal Hessian blocks, shape (3, 3, num_atoms).
    """
    hess = np.zeros((3, 3, num_atoms), dtype=np.float64)
    for i in range(num_atoms):
        for b in range(3):
            cp, cm = coords.copy(), coords.copy()
            cp[b, i] += h
            cm[b, i] -= h
            fp, fm = grad_func(cp), grad_func(cm)
            for a in range(3):
                hess[a, b, i] = -(fp[a, i] - fm[a, i]) / (2.0 * h)
    return hess


def compute_2body_gradient_numerical(num_atoms: int, z: np.ndarray, coords: np.ndarray,
                                     ehomo: np.ndarray, alpha: np.ndarray,
                                     p: PFDParams = PARAMS) -> Tuple[np.ndarray, float, np.ndarray]:
    """
    Compute two-body gradient numerically.
    
    Parameters
    ----------
    num_atoms : int
        Number of atoms.
    z : np.ndarray
        Atomic numbers.
    coords : np.ndarray
        Coordinates in Bohr, shape (3, num_atoms).
    ehomo : np.ndarray
        HOMO energies in Hartree.
    alpha : np.ndarray
        Polarizabilities in Bohr³.
    p : PFDParams, optional
        Model parameters.
        
    Returns
    -------
    forces : np.ndarray
        Numerical forces, shape (3, num_atoms).
    energy : float
        Two-body energy in Hartree.
    pair_energies : np.ndarray
        Matrix of pair energies.
    """
    from .energy import compute_2body_energy
    e, pe = compute_2body_energy(num_atoms, z, coords, ehomo, alpha, p)
    f = numerical_gradient(lambda c: compute_2body_energy(num_atoms, z, c, ehomo, alpha, p)[0], coords, num_atoms)
    return f, e, pe


def compute_3body_gradient_numerical(num_atoms: int, z: np.ndarray, coords: np.ndarray,
                                     ehomo: np.ndarray, alpha: np.ndarray,
                                     p: PFDParams = PARAMS) -> Tuple[np.ndarray, float]:
    """
    Compute three-body gradient numerically.
    
    Parameters
    ----------
    num_atoms : int
        Number of atoms.
    z : np.ndarray
        Atomic numbers.
    coords : np.ndarray
        Coordinates in Bohr, shape (3, num_atoms).
    ehomo : np.ndarray
        HOMO energies in Hartree.
    alpha : np.ndarray
        Polarizabilities in Bohr³.
    p : PFDParams, optional
        Model parameters.
        
    Returns
    -------
    forces : np.ndarray
        Numerical forces, shape (3, num_atoms).
    energy : float
        Three-body energy in Hartree.
    """
    from .energy import compute_3body_energy
    e = compute_3body_energy(num_atoms, z, coords, ehomo, alpha, p)
    f = numerical_gradient(lambda c: compute_3body_energy(num_atoms, z, c, ehomo, alpha, p), coords, num_atoms)
    return f, e


def compute_hbond_gradient_numerical(num_atoms: int, z: np.ndarray, coords: np.ndarray,
                                     ehomo: np.ndarray, alpha: np.ndarray, pair_e: np.ndarray,
                                     p: PFDParams = PARAMS) -> Tuple[np.ndarray, float]:
    """
    Compute H-bond gradient numerically.
    
    Note: The H-bond energy depends on pair energies which change with geometry,
    so pair energies are recomputed at each displaced geometry.
    
    Parameters
    ----------
    num_atoms : int
        Number of atoms.
    z : np.ndarray
        Atomic numbers.
    coords : np.ndarray
        Coordinates in Bohr, shape (3, num_atoms).
    ehomo : np.ndarray
        HOMO energies in Hartree.
    alpha : np.ndarray
        Polarizabilities in Bohr³.
    pair_e : np.ndarray
        Two-body pair energies (used for energy at reference geometry).
    p : PFDParams, optional
        Model parameters.
        
    Returns
    -------
    forces : np.ndarray
        Numerical forces, shape (3, num_atoms).
    energy : float
        H-bond energy in Hartree.
    """
    from .energy import compute_2body_energy, compute_hbond_energy
    e = compute_hbond_energy(num_atoms, z, coords, ehomo, pair_e, p)
    def efunc(c):
        _, pe = compute_2body_energy(num_atoms, z, c, ehomo, alpha, p)
        return compute_hbond_energy(num_atoms, z, c, ehomo, pe, p)
    return numerical_gradient(efunc, coords, num_atoms), e


def compute_2body_hessian_numerical(num_atoms: int, z: np.ndarray, coords: np.ndarray,
                                    ehomo: np.ndarray, alpha: np.ndarray,
                                    p: PFDParams = PARAMS) -> Tuple[np.ndarray, np.ndarray, float, np.ndarray]:
    """
    Compute two-body Hessian numerically.
    
    Parameters
    ----------
    num_atoms : int
        Number of atoms.
    z : np.ndarray
        Atomic numbers.
    coords : np.ndarray
        Coordinates in Bohr, shape (3, num_atoms).
    ehomo : np.ndarray
        HOMO energies in Hartree.
    alpha : np.ndarray
        Polarizabilities in Bohr³.
    p : PFDParams, optional
        Model parameters.
        
    Returns
    -------
    hessian : np.ndarray
        Numerical Hessian, shape (3, 3, num_atoms).
    forces : np.ndarray
        Analytical forces (from gradient function), shape (3, num_atoms).
    energy : float
        Two-body energy in Hartree.
    pair_energies : np.ndarray
        Matrix of pair energies.
    """
    from .gradients import compute_2body_gradient
    def gfunc(c):
        return compute_2body_gradient(num_atoms, z, c, ehomo, alpha, p)[0]
    h = numerical_hessian(gfunc, coords, num_atoms)
    f, e, pe, _ = compute_2body_gradient(num_atoms, z, coords, ehomo, alpha, p)
    return h, f, e, pe


def compute_3body_hessian_numerical(num_atoms: int, z: np.ndarray, coords: np.ndarray,
                                    ehomo: np.ndarray, alpha: np.ndarray,
                                    p: PFDParams = PARAMS) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Compute three-body Hessian numerically.
    
    Parameters
    ----------
    num_atoms : int
        Number of atoms.
    z : np.ndarray
        Atomic numbers.
    coords : np.ndarray
        Coordinates in Bohr, shape (3, num_atoms).
    ehomo : np.ndarray
        HOMO energies in Hartree.
    alpha : np.ndarray
        Polarizabilities in Bohr³.
    p : PFDParams, optional
        Model parameters.
        
    Returns
    -------
    hessian : np.ndarray
        Numerical Hessian, shape (3, 3, num_atoms).
    forces : np.ndarray
        Analytical forces, shape (3, num_atoms).
    energy : float
        Three-body energy in Hartree.
    """
    from .gradients import compute_3body_gradient
    def gfunc(c):
        return compute_3body_gradient(num_atoms, z, c, ehomo, alpha, p)[0]
    h = numerical_hessian(gfunc, coords, num_atoms)
    f, e = compute_3body_gradient(num_atoms, z, coords, ehomo, alpha, p)
    return h, f, e


def compute_hbond_hessian_numerical(num_atoms: int, z: np.ndarray, coords: np.ndarray,
                                    ehomo: np.ndarray, alpha: np.ndarray, pair_e: np.ndarray,
                                    p: PFDParams = PARAMS) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Compute H-bond Hessian numerically.
    
    Note: Pair energies and force factors are recomputed at each displaced
    geometry to ensure correct numerical derivatives.
    
    Parameters
    ----------
    num_atoms : int
        Number of atoms.
    z : np.ndarray
        Atomic numbers.
    coords : np.ndarray
        Coordinates in Bohr, shape (3, num_atoms).
    ehomo : np.ndarray
        HOMO energies in Hartree.
    alpha : np.ndarray
        Polarizabilities in Bohr³.
    pair_e : np.ndarray
        Two-body pair energies (used for energy/forces at reference geometry).
    p : PFDParams, optional
        Model parameters.
        
    Returns
    -------
    hessian : np.ndarray
        Numerical Hessian, shape (3, 3, num_atoms).
    forces : np.ndarray
        Analytical forces, shape (3, num_atoms).
    energy : float
        H-bond energy in Hartree.
    """
    from .gradients import compute_2body_gradient, compute_hbond_gradient
    def gfunc(c):
        _, _, pe, pff = compute_2body_gradient(num_atoms, z, c, ehomo, alpha, p)
        return compute_hbond_gradient(num_atoms, z, c, ehomo, alpha, pe, pff, p)[0]
    h = numerical_hessian(gfunc, coords, num_atoms)
    _, _, pe, pff = compute_2body_gradient(num_atoms, z, coords, ehomo, alpha, p)
    f, e = compute_hbond_gradient(num_atoms, z, coords, ehomo, alpha, pe, pff, p)
    return h, f, e


def compute_all_numerical(num_atoms: int, z: np.ndarray, coords: np.ndarray,
                          ehomo: np.ndarray, alpha: np.ndarray,
                          p: PFDParams = PARAMS, compute_hess: bool = True) -> dict:
    """
    Compute all PFD terms using numerical derivatives.
    
    This function provides the same interface as compute_pfd() but uses
    finite differences instead of analytical derivatives. Primarily used
    for validating the analytical implementations.
    
    Parameters
    ----------
    num_atoms : int
        Number of atoms.
    z : np.ndarray
        Atomic numbers.
    coords : np.ndarray
        Coordinates in Bohr, shape (3, num_atoms).
    ehomo : np.ndarray
        HOMO energies in Hartree.
    alpha : np.ndarray
        Polarizabilities in Bohr³.
    p : PFDParams, optional
        Model parameters.
    compute_hess : bool, optional
        Whether to compute Hessians. Default True.
        
    Returns
    -------
    dict
        Dictionary with keys:
        - energy_2body, energy_3body, energy_hbond, energy_total
        - forces_2body, forces_3body, forces_hbond, forces_total
        - hessian_2body, hessian_3body, hessian_hbond, hessian_total
        - pair_energies
    """
    f_2b, e_2b, pe = compute_2body_gradient_numerical(num_atoms, z, coords, ehomo, alpha, p)
    f_3b, e_3b = compute_3body_gradient_numerical(num_atoms, z, coords, ehomo, alpha, p)
    f_hb, e_hb = compute_hbond_gradient_numerical(num_atoms, z, coords, ehomo, alpha, pe, p)
    
    if compute_hess:
        h_2b, _, _, _ = compute_2body_hessian_numerical(num_atoms, z, coords, ehomo, alpha, p)
        h_3b, _, _ = compute_3body_hessian_numerical(num_atoms, z, coords, ehomo, alpha, p)
        h_hb, _, _ = compute_hbond_hessian_numerical(num_atoms, z, coords, ehomo, alpha, pe, p)
    else:
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
