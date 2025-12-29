"""
PFD Energy Module
=================

This module computes the three components of the PFD dispersion energy:

1. **Two-body (pairwise) dispersion**: The dominant C6/R⁶ interaction between atom pairs
2. **Three-body (ATM) dispersion**: Axilrod-Teller-Muto correction for non-additivity
3. **Hydrogen bond correction**: Angular modification of acceptor-acceptor dispersion

Physical Background
-------------------

**Two-Body Dispersion**

The pairwise dispersion energy between atoms i and j is:
    E_ij = -C6_ij · f²(r_ij) / (r²_ij - Rs²)³

where:
- C6_ij is the dispersion coefficient from the London formula
- f(r) is the damping function that attenuates at short range
- Rs² is the singularity shift that modifies the R⁻⁶ form

**Three-Body Dispersion**

The Axilrod-Teller-Muto (ATM) triple-dipole interaction captures
the non-additive component of dispersion:
    E_ijk = -C9 · (3·cos_i·cos_j·cos_k + 1) · f² / (r_ij·r_ik·r_jk)³

The angular factor (3·cos_i·cos_j·cos_k + 1) is negative for linear
arrangements and positive for equilateral triangles.

**Hydrogen Bond Correction**

In H-bonded systems (H-X···Y), the dispersion between acceptors X and Y
is modified by an angular factor depending on the X-H-Y angle:
    E_hb = E_XY · g(θ_XHY)

This accounts for the directional nature of hydrogen bonding.

Usage
-----
>>> from pfd import parse_gaussian, get_atomic_data
>>> from pfd.energy import compute_pfd_energy
>>>
>>> title, n, scf, z, coords = parse_gaussian("molecule.out")
>>> ehomo, alpha = get_atomic_data(z)
>>> e_2b, e_3b, e_hb, pair_e = compute_pfd_energy(n, z, coords, ehomo, alpha)
>>> print(f"Total dispersion: {e_2b + e_3b + e_hb:.6f} Hartree")
"""

import numpy as np
from typing import Tuple
from .parameters import PFDParams, PARAMS, HBOND_ATOMS
from .core import (compute_pair_params, compute_damping_radius, compute_c9, damping_function,
                   compute_triangle_geometry, hbond_angular_factor, identify_hbond_pattern,
                   iter_pairs, iter_triplets)


def compute_2body_energy(num_atoms: int, z: np.ndarray, coords: np.ndarray,
                         ehomo: np.ndarray, alpha: np.ndarray,
                         p: PFDParams = PARAMS) -> Tuple[float, np.ndarray]:
    """
    Compute the two-body (pairwise) dispersion energy.
    
    This is the dominant component of dispersion, summing over all atom pairs:
        E_2body = Σ_{i<j} -C6_ij · f²(r_ij) / (r²_ij - Rs²)³
    
    Parameters
    ----------
    num_atoms : int
        Number of atoms.
    z : np.ndarray
        Atomic numbers, shape (num_atoms,).
    coords : np.ndarray
        Atomic coordinates in Bohr, shape (3, num_atoms).
    ehomo : np.ndarray
        HOMO energies in Hartree, shape (num_atoms,).
    alpha : np.ndarray
        Polarizabilities in Bohr³, shape (num_atoms,).
    p : PFDParams, optional
        Model parameters.
        
    Returns
    -------
    tuple of (float, np.ndarray)
        - energy: Total two-body energy in Hartree
        - pair_energies: Matrix of pair energies, pair_e[i,j] for i<j
    """
    energy = 0.0
    pair_e = np.zeros((num_atoms, num_atoms), dtype=np.float64)
    
    for i, j in iter_pairs(num_atoms, z):
        delta = coords[:, i] - coords[:, j]
        r2 = delta @ delta
        r = np.sqrt(r2)
        if r < 0.1:
            continue
        
        pp = compute_pair_params(ehomo[i], ehomo[j], alpha[i], alpha[j], z[i], z[j], p)
        if r <= 2.0 * pp.rd:
            continue
        
        damp = damping_function(r, pp.rd)
        denom = r2 - pp.rs2
        if abs(denom) < 1e-10:
            continue
        
        e_pair = -pp.c6 * damp.f**2 / denom**3
        energy += e_pair
        pair_e[i, j] = e_pair
    
    return energy, pair_e


def compute_3body_energy(num_atoms: int, z: np.ndarray, coords: np.ndarray,
                         ehomo: np.ndarray, alpha: np.ndarray,
                         p: PFDParams = PARAMS) -> float:
    """
    Compute the three-body (Axilrod-Teller-Muto) dispersion energy.
    
    The ATM term accounts for non-additive dispersion interactions:
        E_3body = Σ_{i<j<k} -C9 · (3·cos_i·cos_j·cos_k + 1) · f² / (r_ij·r_ik·r_jk)³
    
    The angular factor makes the energy:
    - Negative (attractive) for equilateral/acute triangles
    - Positive (repulsive) for linear arrangements
    
    Parameters
    ----------
    num_atoms : int
        Number of atoms.
    z : np.ndarray
        Atomic numbers.
    coords : np.ndarray
        Atomic coordinates in Bohr.
    ehomo : np.ndarray
        HOMO energies in Hartree.
    alpha : np.ndarray
        Polarizabilities in Bohr³.
    p : PFDParams, optional
        Model parameters.
        
    Returns
    -------
    float
        Three-body energy in Hartree.
    """
    energy = 0.0
    ds = p.threebody_damping_scale
    
    for i, j, k in iter_triplets(num_atoms, z):
        tri = compute_triangle_geometry(coords, i, j, k)
        c9 = compute_c9(ehomo[i], ehomo[j], ehomo[k], alpha[i], alpha[j], alpha[k])
        
        rd_ij = ds * compute_damping_radius(ehomo[i], ehomo[j], z[i], z[j], p)
        rd_ik = ds * compute_damping_radius(ehomo[i], ehomo[k], z[i], z[k], p)
        rd_jk = ds * compute_damping_radius(ehomo[j], ehomo[k], z[j], z[k], p)
        
        f_ij = damping_function(tri.r_ij, rd_ij).f
        f_ik = damping_function(tri.r_ik, rd_ik).f
        f_jk = damping_function(tri.r_jk, rd_jk).f
        
        angular = 3.0 * tri.cos_i * tri.cos_j * tri.cos_k + 1.0
        damp = (f_ij * f_ik * f_jk) ** 2
        geom = (tri.r_ij_inv * tri.r_ik_inv * tri.r_jk_inv) ** 3
        energy -= c9 * angular * damp * geom
    
    return energy


def compute_hbond_energy(num_atoms: int, z: np.ndarray, coords: np.ndarray,
                         ehomo: np.ndarray, pair_e: np.ndarray,
                         p: PFDParams = PARAMS) -> float:
    """
    Compute the hydrogen bond dispersion correction.
    
    For H-X···Y systems where X and Y are electronegative atoms (N, O, F, etc.),
    the X-Y dispersion is scaled by an angular factor depending on the X-H-Y angle:
        E_hb = Σ E_XY · g(θ_XHY)
    
    The angular factor g(θ) reduces the X-Y dispersion when H is between them,
    modeling the directional character of hydrogen bonds.
    
    Parameters
    ----------
    num_atoms : int
        Number of atoms.
    z : np.ndarray
        Atomic numbers.
    coords : np.ndarray
        Atomic coordinates in Bohr.
    ehomo : np.ndarray
        HOMO energies (not used but kept for API consistency).
    pair_e : np.ndarray
        Two-body pair energies from compute_2body_energy.
    p : PFDParams, optional
        Model parameters.
        
    Returns
    -------
    float
        H-bond correction energy in Hartree.
    """
    energy = 0.0
    
    for i, j, k in iter_triplets(num_atoms, z):
        h_vertex, _, acc1_rel, acc2_rel = identify_hbond_pattern(z[i], z[j], z[k], HBOND_ATOMS)
        if h_vertex is None:
            continue
        
        idx_map = {0: i, 1: j, 2: k}
        acc1, acc2 = idx_map[acc1_rel], idx_map[acc2_rel]
        tri = compute_triangle_geometry(coords, i, j, k)
        cos_h = {'i': tri.cos_i, 'j': tri.cos_j, 'k': tri.cos_k}[h_vertex]
        
        lo, hi = (acc1, acc2) if acc1 < acc2 else (acc2, acc1)
        hb_f, _ = hbond_angular_factor(cos_h, p)
        if z[acc1] == 34 or z[acc2] == 34:  # Se reduction
            hb_f /= 2.0
        energy += pair_e[lo, hi] * hb_f
    
    return energy


def compute_pfd_energy(num_atoms: int, z: np.ndarray, coords: np.ndarray,
                       ehomo: np.ndarray, alpha: np.ndarray,
                       p: PFDParams = PARAMS) -> Tuple[float, float, float, np.ndarray]:
    """
    Compute all PFD dispersion energies.
    
    This is the main interface for energy calculations.
    
    Parameters
    ----------
    num_atoms : int
        Number of atoms.
    z : np.ndarray
        Atomic numbers.
    coords : np.ndarray
        Atomic coordinates in Bohr.
    ehomo : np.ndarray
        HOMO energies in Hartree.
    alpha : np.ndarray
        Polarizabilities in Bohr³.
    p : PFDParams, optional
        Model parameters.
        
    Returns
    -------
    tuple of (float, float, float, np.ndarray)
        - e_2body: Two-body dispersion energy
        - e_3body: Three-body dispersion energy
        - e_hbond: Hydrogen bond correction
        - pair_energies: Matrix of pair energies
        
    Example
    -------
    >>> e_2b, e_3b, e_hb, _ = compute_pfd_energy(n, z, coords, ehomo, alpha)
    >>> e_total = e_2b + e_3b + e_hb
    """
    e_2b, pair_e = compute_2body_energy(num_atoms, z, coords, ehomo, alpha, p)
    e_3b = compute_3body_energy(num_atoms, z, coords, ehomo, alpha, p)
    e_hb = compute_hbond_energy(num_atoms, z, coords, ehomo, pair_e, p)
    return e_2b, e_3b, e_hb, pair_e
