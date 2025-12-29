"""
PFD Gradients Module
====================

This module computes analytical gradients (forces) of the PFD dispersion energy.

Overview for Non-Specialists
----------------------------
In molecular simulations, we often need to know not just the energy of a system,
but also how the energy changes when atoms move. This "rate of change" is the
gradient, and its negative gives the force on each atom:

    Force = -∇E = -(∂E/∂x, ∂E/∂y, ∂E/∂z)

Forces tell us which direction atoms would move to lower the energy, which is
essential for geometry optimization and molecular dynamics simulations.

Mathematical Background
-----------------------

**Two-Body Gradient**

For a pair of atoms i and j separated by distance r, the two-body dispersion
energy is:
    E_ij = -C6 · f²(r) / (r² - Rs²)³

The force requires computing dE/dr, then projecting onto Cartesian coordinates.
Using the chain rule:

    dE/dr = ∂E/∂(denom) · d(denom)/dr + ∂E/∂(f²) · d(f²)/dr

where denom = r² - Rs². This gives:

    dE/dr = 6r · C6 · f² / (r² - Rs²)⁴  -  2 · C6 · f · (df/dr) / (r² - Rs²)³

The Cartesian force on atom i from atom j is then:
    F_i = -(dE/dr) · (r_ij / |r_ij|)

where r_ij = r_i - r_j is the vector from j to i.

**Three-Body Gradient**

The three-body energy depends on three distances (r_ij, r_ik, r_jk) and three
angles (θ_i, θ_j, θ_k). The gradient requires:

1. ∂E/∂r_ij, ∂E/∂r_ik, ∂E/∂r_jk (derivatives w.r.t. distances)
2. ∂E/∂cos_i, ∂E/∂cos_j, ∂E/∂cos_k (derivatives w.r.t. angles)
3. ∂cos/∂r (how angles depend on distances, from law of cosines)

These are combined via the chain rule:
    dE/dr_ij = ∂E/∂r_ij + Σ_c (∂E/∂cos_c) · (∂cos_c/∂r_ij)

Then transformed to Cartesian coordinates using the unit vectors along each edge.

**H-Bond Gradient**

The H-bond correction modifies the acceptor-acceptor dispersion by an angular
factor g(θ) that depends on the X-H-Y angle. The gradient has two contributions:
1. The original pair force scaled by g(θ)
2. A new term from dg/dθ as the angle changes with atomic positions

Implementation Notes
--------------------
- Forces are computed simultaneously with energies to reuse intermediate values
- pair_ff (pair force factors) are saved for use in H-bond gradient
- Unit vectors (r_ij/|r_ij|) transform scalar dE/dr to vector forces

Functions
---------
- compute_2body_gradient: Pairwise dispersion forces
- compute_3body_gradient: Three-body (ATM) forces  
- compute_hbond_gradient: H-bond correction forces
- compute_all_gradients: Combined interface returning all components

Output Format
-------------
Forces are returned as numpy arrays with shape (3, num_atoms):
- forces[0, i] = x-component of force on atom i
- forces[1, i] = y-component of force on atom i
- forces[2, i] = z-component of force on atom i

Units are Hartree/Bohr (atomic units).
"""

import numpy as np
from typing import Tuple
from .parameters import PFDParams, PARAMS, HBOND_ATOMS
from .core import (compute_pair_params, compute_damping_radius, compute_c9, damping_function,
                   compute_triangle_geometry, compute_cosine_derivatives,
                   hbond_angular_factor, identify_hbond_pattern, iter_pairs, iter_triplets)


def compute_2body_gradient(num_atoms: int, z: np.ndarray, coords: np.ndarray,
                           ehomo: np.ndarray, alpha: np.ndarray,
                           p: PFDParams = PARAMS) -> Tuple[np.ndarray, float, np.ndarray, np.ndarray]:
    """
    Compute two-body dispersion forces.
    
    For each atom pair, computes the force arising from the pairwise dispersion
    energy E_ij = -C6·f²/(r²-Rs²)³.
    
    Parameters
    ----------
    num_atoms : int
        Number of atoms in the molecule.
    z : np.ndarray
        Atomic numbers, shape (num_atoms,).
    coords : np.ndarray
        Atomic coordinates in Bohr, shape (3, num_atoms).
    ehomo : np.ndarray
        HOMO energies in Hartree, shape (num_atoms,).
    alpha : np.ndarray
        Polarizabilities in Bohr³, shape (num_atoms,).
    p : PFDParams, optional
        Model parameters (default: standard PFD parameters).
        
    Returns
    -------
    forces : np.ndarray
        Forces on each atom, shape (3, num_atoms), units Hartree/Bohr.
    energy : float
        Total two-body dispersion energy in Hartree.
    pair_energies : np.ndarray
        Matrix of pair energies, pair_e[i,j] for i < j.
    pair_force_factors : np.ndarray
        Force factors ff = -dE/dr/r for each pair, used by H-bond gradient.
        
    Notes
    -----
    The force factor ff is defined such that F_i = -ff · (r_i - r_j).
    This is stored for reuse in the H-bond gradient calculation.
    """
    forces = np.zeros((3, num_atoms), dtype=np.float64)
    energy = 0.0
    pair_e = np.zeros((num_atoms, num_atoms), dtype=np.float64)
    pair_ff = np.zeros((num_atoms, num_atoms), dtype=np.float64)
    
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
        
        # Energy: E = -C6 · f² / denom³
        disp = pp.c6 / denom**3
        e_pair = -disp * damp.f**2
        energy += e_pair
        pair_e[i, j] = e_pair
        
        # Gradient: dE/dr via chain rule
        # d(1/denom³)/dr = -3/denom⁴ · d(denom)/dr = -3/denom⁴ · 2r = -6r/denom⁴
        ddisp_dr = -6.0 * r * disp / denom
        # d(f²)/dr = 2f · df/dr
        dE_dr = ddisp_dr * damp.f**2 + 2.0 * disp * damp.f * damp.df
        
        # Force factor: F = -dE/dr · (delta/r) = ff · delta
        ff = -dE_dr / r
        pair_ff[i, j] = ff
        
        # Apply Newton's 3rd law: F_i = -F_j
        forces[:, i] -= delta * ff
        forces[:, j] += delta * ff
    
    return forces, energy, pair_e, pair_ff


def compute_3body_gradient(num_atoms: int, z: np.ndarray, coords: np.ndarray,
                           ehomo: np.ndarray, alpha: np.ndarray,
                           p: PFDParams = PARAMS) -> Tuple[np.ndarray, float]:
    """
    Compute three-body (Axilrod-Teller-Muto) dispersion forces.
    
    The three-body energy for triplet i-j-k is:
        E_ijk = -C9 · (3·cos_i·cos_j·cos_k + 1) · f² / (r_ij·r_ik·r_jk)³
    
    The gradient requires derivatives with respect to all three distances,
    accounting for both the 1/r³ geometric factor and the angular factor.
    
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
        Forces, shape (3, num_atoms), Hartree/Bohr.
    energy : float
        Three-body energy in Hartree.
        
    Notes
    -----
    The gradient computation involves:
    1. ∂E/∂r_m = -3·E/r_m (from the 1/r³ factors)
    2. ∂E/∂f_m² = E/f_m² (from the damping)  
    3. ∂E/∂cos_c = C9·prefac·3·(product of other two cosines)
    4. ∂cos_c/∂r_m from the law of cosines derivatives
    
    These are combined and transformed to Cartesian forces using unit vectors.
    """
    forces = np.zeros((3, num_atoms), dtype=np.float64)
    energy = 0.0
    ds = p.threebody_damping_scale
    
    for i, j, k in iter_triplets(num_atoms, z):
        tri = compute_triangle_geometry(coords, i, j, k)
        dcos = compute_cosine_derivatives(tri)
        
        c9 = compute_c9(ehomo[i], ehomo[j], ehomo[k], alpha[i], alpha[j], alpha[k])
        rd_ij = ds * compute_damping_radius(ehomo[i], ehomo[j], z[i], z[j], p)
        rd_ik = ds * compute_damping_radius(ehomo[i], ehomo[k], z[i], z[k], p)
        rd_jk = ds * compute_damping_radius(ehomo[j], ehomo[k], z[j], z[k], p)
        
        d_ij, d_ik, d_jk = damping_function(tri.r_ij, rd_ij), damping_function(tri.r_ik, rd_ik), damping_function(tri.r_jk, rd_jk)
        f_ij, f_ik, f_jk = d_ij.f, d_ik.f, d_jk.f
        
        if abs(f_ij) < 1e-15 or abs(f_ik) < 1e-15 or abs(f_jk) < 1e-15:
            continue
        
        # Energy: E = -C9 · angular · (f_ij·f_ik·f_jk)² · (1/r_ij·r_ik·r_jk)³
        prefac = (f_ij * f_ik * f_jk)**2 * (tri.r_ij_inv * tri.r_ik_inv * tri.r_jk_inv)**3
        angular = 3.0 * tri.cos_i * tri.cos_j * tri.cos_k + 1.0
        disp = c9 * angular * prefac
        energy -= disp
        
        # Partial derivatives of E with respect to intermediate quantities
        # ∂E/∂r_m = -3·E/r_m (from geometric factor)
        dE_dr_ij, dE_dr_ik, dE_dr_jk = -3.0 * disp * tri.r_ij_inv, -3.0 * disp * tri.r_ik_inv, -3.0 * disp * tri.r_jk_inv
        
        # ∂E/∂cos_c = C9 · prefac · 3 · (product of other cosines)
        dE_dcos_i = c9 * prefac * 3.0 * tri.cos_j * tri.cos_k
        dE_dcos_j = c9 * prefac * 3.0 * tri.cos_i * tri.cos_k
        dE_dcos_k = c9 * prefac * 3.0 * tri.cos_i * tri.cos_j
        
        # ∂E/∂(f²) = E/f² (from damping factors)
        dE_df_ij, dE_df_ik, dE_df_jk = disp / f_ij**2, disp / f_ik**2, disp / f_jk**2
        
        # d(f²)/dr for each edge (chain rule through damping function)
        x_ij, x_ik, x_jk = (tri.r_ij - 2*rd_ij)/rd_ij, (tri.r_ik - 2*rd_ik)/rd_ik, (tri.r_jk - 2*rd_jk)/rd_jk
        df2_ij = 2.0 * f_ij * (1.0 + 2.0*x_ij) * np.exp(-2.0*x_ij) / rd_ij
        df2_ik = 2.0 * f_ik * (1.0 + 2.0*x_ik) * np.exp(-2.0*x_ik) / rd_ik
        df2_jk = 2.0 * f_jk * (1.0 + 2.0*x_jk) * np.exp(-2.0*x_jk) / rd_jk
        
        # Total dE/dr_m combining all contributions via chain rule
        dE_rij = (dE_dr_ij + dE_df_ij*df2_ij + dE_dcos_k*dcos.dcos_k_drij +
                  dE_dcos_j*dcos.dcos_j_drij + dE_dcos_i*dcos.dcos_i_drij)
        dE_rik = (dE_dr_ik + dE_df_ik*df2_ik + dE_dcos_k*dcos.dcos_k_drik +
                  dE_dcos_j*dcos.dcos_j_drik + dE_dcos_i*dcos.dcos_i_drik)
        dE_rjk = (dE_dr_jk + dE_df_jk*df2_jk + dE_dcos_k*dcos.dcos_k_drjk +
                  dE_dcos_j*dcos.dcos_j_drjk + dE_dcos_i*dcos.dcos_i_drjk)
        
        # Transform to Cartesian forces: F = -∇E
        # For atom i: connected to j via r_ij and to k via r_ik
        forces[:, i] += tri.vec_ij * tri.r_ij_inv * dE_rij + tri.vec_ik * tri.r_ik_inv * dE_rik
        # For atom j: connected to i via -r_ij and to k via r_jk
        forces[:, j] += -tri.vec_ij * tri.r_ij_inv * dE_rij + tri.vec_jk * tri.r_jk_inv * dE_rjk
        # For atom k: connected to i via -r_ik and to j via -r_jk
        forces[:, k] += -tri.vec_ik * tri.r_ik_inv * dE_rik - tri.vec_jk * tri.r_jk_inv * dE_rjk
    
    return forces, energy


def compute_hbond_gradient(num_atoms: int, z: np.ndarray, coords: np.ndarray,
                           ehomo: np.ndarray, alpha: np.ndarray, pair_e: np.ndarray,
                           pair_ff: np.ndarray, p: PFDParams = PARAMS) -> Tuple[np.ndarray, float]:
    """
    Compute hydrogen bond correction forces.
    
    For H-X···Y systems (H bonded to electronegative X, interacting with Y),
    the X-Y dispersion is modified by an angular factor g(θ) where θ is the
    X-H-Y angle:
        E_hb = E_XY · g(θ)
    
    The gradient has two terms:
    1. The original X-Y force scaled by g(θ)
    2. A torque-like term from dg/dθ as atoms move
    
    Parameters
    ----------
    num_atoms : int
        Number of atoms.
    z : np.ndarray
        Atomic numbers.
    coords : np.ndarray
        Coordinates in Bohr.
    ehomo : np.ndarray
        HOMO energies (not used, kept for API consistency).
    alpha : np.ndarray
        Polarizabilities (not used, kept for API consistency).
    pair_e : np.ndarray
        Two-body pair energies from compute_2body_gradient.
    pair_ff : np.ndarray
        Pair force factors from compute_2body_gradient.
    p : PFDParams, optional
        Model parameters.
        
    Returns
    -------
    forces : np.ndarray
        Forces, shape (3, num_atoms), Hartree/Bohr.
    energy : float
        H-bond correction energy in Hartree.
        
    Notes
    -----
    The angular factor g(θ) reduces to zero for linear X-H-Y (θ = 180°)
    and is maximal for bent geometries. The derivative dg/d(cos θ) 
    contributes forces that resist changes in the H-bond angle.
    """
    forces = np.zeros((3, num_atoms), dtype=np.float64)
    energy = 0.0
    
    for i, j, k in iter_triplets(num_atoms, z):
        h_vertex, _, acc1_rel, acc2_rel = identify_hbond_pattern(z[i], z[j], z[k], HBOND_ATOMS)
        if h_vertex is None:
            continue
        
        idx_map = {0: i, 1: j, 2: k}
        acc1, acc2 = idx_map[acc1_rel], idx_map[acc2_rel]
        tri = compute_triangle_geometry(coords, i, j, k)
        dcos = compute_cosine_derivatives(tri)
        
        # Get angle at hydrogen and identify acceptor-acceptor vector
        if h_vertex == 'i':
            cos_h, dcos_r = tri.cos_i, (dcos.dcos_i_drij, dcos.dcos_i_drik, dcos.dcos_i_drjk)
            acc_vec, acc_r_inv = tri.vec_jk, tri.r_jk_inv
        elif h_vertex == 'j':
            cos_h, dcos_r = tri.cos_j, (dcos.dcos_j_drij, dcos.dcos_j_drik, dcos.dcos_j_drjk)
            acc_vec, acc_r_inv = tri.vec_ik, tri.r_ik_inv
        else:
            cos_h, dcos_r = tri.cos_k, (dcos.dcos_k_drij, dcos.dcos_k_drik, dcos.dcos_k_drjk)
            acc_vec, acc_r_inv = tri.vec_ij, tri.r_ij_inv
        
        # Get pair energy and force factor for acceptor-acceptor interaction
        lo, hi = (acc1, acc2) if acc1 < acc2 else (acc2, acc1)
        pe, ff = pair_e[lo, hi], pair_ff[lo, hi]
        
        # Angular correction factor and its derivative
        hb_f, hb_df = hbond_angular_factor(cos_h, p)
        if z[acc1] == 34 or z[acc2] == 34:  # Selenium reduction
            hb_f, hb_df = hb_f / 2.0, hb_df / 2.0
        
        energy += pe * hb_f
        
        # Cartesian derivatives of cos(θ_H) for each atom
        # d(cos)/dx_i = Σ_m (∂cos/∂r_m) · (∂r_m/∂x_i)
        dcos_xi = tri.vec_ij*tri.r_ij_inv*dcos_r[0] + tri.vec_ik*tri.r_ik_inv*dcos_r[1]
        dcos_xj = -tri.vec_ij*tri.r_ij_inv*dcos_r[0] + tri.vec_jk*tri.r_jk_inv*dcos_r[2]
        dcos_xk = -tri.vec_ik*tri.r_ik_inv*dcos_r[1] - tri.vec_jk*tri.r_jk_inv*dcos_r[2]
        
        # Force contributions:
        # 1. Scaled pair force: g(θ) · F_pair on acceptors
        # 2. Angular term: E_pair · dg/d(cos) · d(cos)/dx on all atoms
        if h_vertex == 'i':
            forces[:, j] -= hb_f*ff*acc_vec + pe*hb_df*dcos_xj
            forces[:, k] -= -hb_f*ff*acc_vec + pe*hb_df*dcos_xk
            forces[:, i] -= pe*hb_df*dcos_xi
        elif h_vertex == 'j':
            forces[:, i] -= hb_f*ff*acc_vec + pe*hb_df*dcos_xi
            forces[:, k] -= -hb_f*ff*acc_vec + pe*hb_df*dcos_xk
            forces[:, j] -= pe*hb_df*dcos_xj
        else:
            forces[:, i] -= hb_f*ff*acc_vec + pe*hb_df*dcos_xi
            forces[:, j] -= -hb_f*ff*acc_vec + pe*hb_df*dcos_xj
            forces[:, k] -= pe*hb_df*dcos_xk
    
    return forces, energy


def compute_all_gradients(num_atoms: int, z: np.ndarray, coords: np.ndarray,
                          ehomo: np.ndarray, alpha: np.ndarray,
                          p: PFDParams = PARAMS) -> Tuple[np.ndarray, np.ndarray, np.ndarray, float, float, float, np.ndarray]:
    """
    Compute all PFD gradient components.
    
    This is the main interface for gradient calculations, returning forces
    and energies for all three dispersion components.
    
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
    f_2body : np.ndarray
        Two-body forces, shape (3, num_atoms).
    f_3body : np.ndarray
        Three-body forces, shape (3, num_atoms).
    f_hbond : np.ndarray
        H-bond correction forces, shape (3, num_atoms).
    e_2body : float
        Two-body energy in Hartree.
    e_3body : float
        Three-body energy in Hartree.
    e_hbond : float
        H-bond correction energy in Hartree.
    pair_energies : np.ndarray
        Matrix of pairwise energies.
        
    Example
    -------
    >>> f_2b, f_3b, f_hb, e_2b, e_3b, e_hb, pe = compute_all_gradients(n, z, coords, ehomo, alpha)
    >>> f_total = f_2b + f_3b + f_hb
    >>> e_total = e_2b + e_3b + e_hb
    """
    f_2b, e_2b, pair_e, pair_ff = compute_2body_gradient(num_atoms, z, coords, ehomo, alpha, p)
    f_3b, e_3b = compute_3body_gradient(num_atoms, z, coords, ehomo, alpha, p)
    f_hb, e_hb = compute_hbond_gradient(num_atoms, z, coords, ehomo, alpha, pair_e, pair_ff, p)
    return f_2b, f_3b, f_hb, e_2b, e_3b, e_hb, pair_e
