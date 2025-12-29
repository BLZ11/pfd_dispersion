"""
PFD Core Module
===============

This module provides the fundamental building blocks for PFD dispersion calculations:
geometry computations, damping functions, dispersion coefficients, and iterators.

Mathematical Background
-----------------------

**Triangle Geometry (for 3-body terms)**

For three atoms i, j, k, the triangle is defined by:
- Edge vectors: vec_ij = r_i - r_j, vec_ik = r_i - r_k, vec_jk = r_j - r_k  
- Edge lengths: r_ij, r_ik, r_jk
- Interior angles via law of cosines: cos(θ_i) = (r_ij² + r_ik² - r_jk²) / (2·r_ij·r_ik)

**Damping Function**

The damping function f(r) smoothly turns off dispersion at short range to avoid
divergence and double-counting with the underlying electronic structure method:

    f(r) = 1 - (1 + x)·exp(-2x)    where x = (r - 2·rd) / rd

This gives f → 0 as r → 0 and f → 1 as r → ∞.

**Dispersion Coefficients**

C6 coefficients are computed from atomic polarizabilities (α) and HOMO energies (ε)
using a modified London formula:

    C6_ij = (3/2) · α_i · α_j · ε_i · ε_j / (ε_i + ε_j)

C9 three-body coefficients use the triple-dipole formula:

    C9_ijk = (3/2) · ε_i·ε_j·ε_k · (ε_i+ε_j+ε_k) · α_i·α_j·α_k / [(ε_i+ε_j)(ε_i+ε_k)(ε_j+ε_k)]

**Hydrogen Bond Correction**

For H-X···Y systems (X, Y = electronegative atoms), the dispersion between
acceptors X and Y is modified by an angular factor depending on the X-H-Y angle.

Data Structures
---------------
- TriangleGeometry: Complete geometric description of an atom triplet
- CosineDerivatives: ∂cos/∂r derivatives for gradient/Hessian calculations  
- DampingResult: f, df/dr, d²f/dr² for energy, gradient, Hessian
- PairParams: Precomputed pair parameters (rd, C6, rs²)

Key Functions
-------------
- compute_triangle_geometry: Set up triplet geometry
- damping_function: Evaluate damping and derivatives
- compute_c6, compute_c9: Dispersion coefficients
- compute_3body_hessian_terms: Analytical second derivatives for triplets
- iter_pairs, iter_triplets: Loop over valid atom combinations
"""

import numpy as np
from typing import Tuple, NamedTuple
from .parameters import PFDParams, PARAMS, PI


# =============================================================================
# Data Structures
# =============================================================================

class TriangleGeometry(NamedTuple):
    """
    Complete geometry for a triangle formed by atoms i, j, k.
    
    The triangle is characterized by three edge vectors, their lengths,
    and the interior angles at each vertex (as cosines).
    
    Attributes
    ----------
    vec_ij, vec_ik, vec_jk : np.ndarray
        Edge vectors (3D Cartesian), pointing from second to first index.
    r_ij, r_ik, r_jk : float
        Edge lengths in Bohr.
    r_ij2, r_ik2, r_jk2 : float  
        Squared edge lengths (precomputed for efficiency).
    r_ij_inv, r_ik_inv, r_jk_inv : float
        Inverse edge lengths (1/r).
    cos_i, cos_j, cos_k : float
        Cosines of interior angles at vertices i, j, k.
        
    Notes
    -----
    The angle at vertex i is the angle ∠jik, opposite to edge jk.
    Cosines are clipped to [-1, 1] to handle numerical edge cases.
    """
    vec_ij: np.ndarray; vec_ik: np.ndarray; vec_jk: np.ndarray
    r_ij: float; r_ik: float; r_jk: float
    r_ij2: float; r_ik2: float; r_jk2: float
    r_ij_inv: float; r_ik_inv: float; r_jk_inv: float
    cos_i: float; cos_j: float; cos_k: float


class CosineDerivatives(NamedTuple):
    """
    First derivatives of triangle cosines with respect to edge lengths.
    
    These are needed for computing gradients and Hessians of the three-body
    dispersion energy, which depends on the angles via the Axilrod-Teller term.
    
    Attributes
    ----------
    dcos_i_drij, dcos_i_drik, dcos_i_drjk : float
        Derivatives of cos(θ_i) with respect to r_ij, r_ik, r_jk.
    dcos_j_drij, dcos_j_drik, dcos_j_drjk : float
        Derivatives of cos(θ_j) with respect to r_ij, r_ik, r_jk.
    dcos_k_drij, dcos_k_drik, dcos_k_drjk : float
        Derivatives of cos(θ_k) with respect to r_ij, r_ik, r_jk.
        
    Notes
    -----
    Derived from the law of cosines by differentiation. For example:
    cos_i = (r_ij² + r_ik² - r_jk²) / (2·r_ij·r_ik)
    """
    dcos_i_drij: float; dcos_i_drik: float; dcos_i_drjk: float
    dcos_j_drij: float; dcos_j_drik: float; dcos_j_drjk: float
    dcos_k_drij: float; dcos_k_drik: float; dcos_k_drjk: float


class DampingResult(NamedTuple):
    """
    Damping function value and its derivatives.
    
    The damping function smoothly attenuates dispersion at short range.
    
    Attributes
    ----------
    f : float
        Damping function value, f ∈ [0, 1].
    df : float
        First derivative df/dr (for gradient calculations).
    d2f : float
        Second derivative d²f/dr² (for Hessian calculations).
    """
    f: float; df: float; d2f: float


class PairParams(NamedTuple):
    """
    Precomputed parameters for an atom pair.
    
    Attributes
    ----------
    rd : float
        Damping radius in Bohr. Dispersion is attenuated for r < 2·rd.
    c6 : float
        C6 dispersion coefficient in Hartree·Bohr⁶.
    rs2 : float
        Squared singularity shift parameter. Used in the modified dispersion
        formula: E = -C6·f² / (r² - rs²)³ to avoid the R⁻⁶ singularity.
    """
    rd: float; c6: float; rs2: float


# =============================================================================
# Geometry Functions
# =============================================================================

def compute_triangle_geometry(coords: np.ndarray, i: int, j: int, k: int) -> TriangleGeometry:
    """
    Compute complete triangle geometry for three atoms.
    
    Parameters
    ----------
    coords : np.ndarray
        Atomic coordinates, shape (3, n_atoms), in Bohr.
    i, j, k : int
        Atom indices (0-based).
        
    Returns
    -------
    TriangleGeometry
        Complete geometric description of the i-j-k triangle.
    """
    vec_ij, vec_ik, vec_jk = coords[:, i] - coords[:, j], coords[:, i] - coords[:, k], coords[:, j] - coords[:, k]
    r_ij2, r_ik2, r_jk2 = vec_ij @ vec_ij, vec_ik @ vec_ik, vec_jk @ vec_jk
    r_ij, r_ik, r_jk = np.sqrt(r_ij2), np.sqrt(r_ik2), np.sqrt(r_jk2)
    r_ij_inv, r_ik_inv, r_jk_inv = 1.0 / r_ij, 1.0 / r_ik, 1.0 / r_jk
    cos_i = np.clip((r_ij2 + r_ik2 - r_jk2) / (2.0 * r_ij * r_ik), -1.0, 1.0)
    cos_j = np.clip((r_ij2 + r_jk2 - r_ik2) / (2.0 * r_ij * r_jk), -1.0, 1.0)
    cos_k = np.clip((r_ik2 + r_jk2 - r_ij2) / (2.0 * r_ik * r_jk), -1.0, 1.0)
    return TriangleGeometry(vec_ij, vec_ik, vec_jk, r_ij, r_ik, r_jk, r_ij2, r_ik2, r_jk2,
                           r_ij_inv, r_ik_inv, r_jk_inv, cos_i, cos_j, cos_k)


def compute_cosine_derivatives(t: TriangleGeometry) -> CosineDerivatives:
    """
    Compute derivatives of cosines with respect to edge lengths.
    
    Uses analytical derivatives of the law of cosines formula.
    
    Parameters
    ----------
    t : TriangleGeometry
        Precomputed triangle geometry.
        
    Returns
    -------
    CosineDerivatives
        All nine ∂cos/∂r derivatives.
    """
    return CosineDerivatives(
        (t.r_ij2 - t.r_ik2 + t.r_jk2) * t.r_ij_inv**2 * t.r_ik_inv / 2.0,
        (t.r_ik2 - t.r_ij2 + t.r_jk2) * t.r_ik_inv**2 * t.r_ij_inv / 2.0,
        -t.r_jk * t.r_ij_inv * t.r_ik_inv,
        (t.r_ij2 - t.r_jk2 + t.r_ik2) * t.r_ij_inv**2 * t.r_jk_inv / 2.0,
        -t.r_ik * t.r_ij_inv * t.r_jk_inv,
        (t.r_jk2 - t.r_ij2 + t.r_ik2) * t.r_jk_inv**2 * t.r_ij_inv / 2.0,
        -t.r_ij * t.r_ik_inv * t.r_jk_inv,
        (t.r_ik2 - t.r_jk2 + t.r_ij2) * t.r_ik_inv**2 * t.r_jk_inv / 2.0,
        (t.r_jk2 - t.r_ik2 + t.r_ij2) * t.r_jk_inv**2 * t.r_ik_inv / 2.0,
    )


def compute_cosine_second_derivatives(t: TriangleGeometry, dc: CosineDerivatives) -> dict:
    """
    Compute second derivatives of cosines with respect to edge lengths.
    
    These are needed for Hessian calculations of the three-body energy.
    
    Parameters
    ----------
    t : TriangleGeometry
        Precomputed triangle geometry.
    dc : CosineDerivatives
        First derivatives (used for chain rule, but not directly here).
        
    Returns
    -------
    dict
        Dictionary with keys like 'd2cos_i_drij_drjk' for ∂²cos_i/(∂r_ij·∂r_jk).
        All 27 combinations of (i,j,k) × (ij,ik,jk) × (ij,ik,jk) are provided.
    """
    r_ij, r_ik, r_jk = t.r_ij, t.r_ik, t.r_jk
    r_ij2, r_ik2, r_jk2 = t.r_ij2, t.r_ik2, t.r_jk2
    sum_r2 = r_ij2 + r_ik2 + r_jk2
    return {
        'd2cos_i_drij_drij': (r_ik2 - r_jk2) / (r_ij**3 * r_ik),
        'd2cos_i_drij_drik': -sum_r2 / (2.0 * r_ij2 * r_ik2),
        'd2cos_i_drij_drjk': r_jk / (r_ij2 * r_ik),
        'd2cos_i_drik_drij': -sum_r2 / (2.0 * r_ij2 * r_ik2),
        'd2cos_i_drik_drik': (r_ij2 - r_jk2) / (r_ij * r_ik**3),
        'd2cos_i_drik_drjk': r_jk / (r_ij * r_ik2),
        'd2cos_i_drjk_drij': r_jk / (r_ij2 * r_ik),
        'd2cos_i_drjk_drik': r_jk / (r_ij * r_ik2),
        'd2cos_i_drjk_drjk': -1.0 / (r_ij * r_ik),
        'd2cos_j_drij_drij': (r_jk2 - r_ik2) / (r_ij**3 * r_jk),
        'd2cos_j_drij_drik': r_ik / (r_ij2 * r_jk),
        'd2cos_j_drij_drjk': -sum_r2 / (2.0 * r_ij2 * r_jk2),
        'd2cos_j_drik_drij': r_ik / (r_ij2 * r_jk),
        'd2cos_j_drik_drik': -1.0 / (r_ij * r_jk),
        'd2cos_j_drik_drjk': r_ik / (r_ij * r_jk2),
        'd2cos_j_drjk_drij': -sum_r2 / (2.0 * r_ij2 * r_jk2),
        'd2cos_j_drjk_drik': r_ik / (r_ij * r_jk2),
        'd2cos_j_drjk_drjk': (r_ij2 - r_ik2) / (r_ij * r_jk**3),
        'd2cos_k_drij_drij': -1.0 / (r_ik * r_jk),
        'd2cos_k_drij_drik': r_ij / (r_ik2 * r_jk),
        'd2cos_k_drij_drjk': r_ij / (r_ik * r_jk2),
        'd2cos_k_drik_drij': r_ij / (r_ik2 * r_jk),
        'd2cos_k_drik_drik': (r_jk2 - r_ij2) / (r_ik**3 * r_jk),
        'd2cos_k_drik_drjk': -sum_r2 / (2.0 * r_ik2 * r_jk2),
        'd2cos_k_drjk_drij': r_ij / (r_ik * r_jk2),
        'd2cos_k_drjk_drik': -sum_r2 / (2.0 * r_ik2 * r_jk2),
        'd2cos_k_drjk_drjk': (r_ik2 - r_ij2) / (r_ik * r_jk**3),
    }


# =============================================================================
# Damping and Dispersion Coefficients
# =============================================================================

def damping_function(r: float, rd: float) -> DampingResult:
    """
    Evaluate the PFD damping function and its derivatives.
    
    The damping function smoothly transitions from 0 (short range) to 1 (long range):
        f(r) = 1 - (1 + x)·exp(-2x),    where x = (r - 2·rd) / rd
    
    Parameters
    ----------
    r : float
        Interatomic distance in Bohr.
    rd : float
        Damping radius in Bohr.
        
    Returns
    -------
    DampingResult
        Named tuple with (f, df/dr, d²f/dr²).
    """
    x = (r - 2.0 * rd) / rd
    exp_x = np.exp(-2.0 * x)
    return DampingResult(1.0 - (1.0 + x) * exp_x,
                        (1.0 + 2.0 * x) * exp_x / rd,
                        -4.0 * x * exp_x / rd**2)


def compute_damping_radius(ehomo_i: float, ehomo_j: float, z_i: int, z_j: int, p: PFDParams) -> float:
    """
    Compute the damping radius for an atom pair.
    
    The damping radius determines where the damping function begins to
    attenuate the dispersion interaction.
    
    Parameters
    ----------
    ehomo_i, ehomo_j : float
        HOMO energies of atoms i and j (Hartree, negative values).
    z_i, z_j : int
        Atomic numbers of atoms i and j.
    p : PFDParams
        Model parameters.
        
    Returns
    -------
    float
        Damping radius in Bohr.
        
    Notes
    -----
    Special scaling is applied for hydrogen-containing pairs.
    """
    rd = p.damping_scale / np.sqrt(-ehomo_i - ehomo_j) + p.damping_offset
    sqrt_he = np.sqrt(p.helium_factor)
    if z_i == 1 or z_j == 1:
        rd *= sqrt_he
    if (z_i == 1 and z_j == 6) or (z_i == 6 and z_j == 1):
        rd /= p.hydrogen_scale * sqrt_he
    return rd


def compute_c6(ehomo_i: float, ehomo_j: float, alpha_i: float, alpha_j: float, p: PFDParams) -> float:
    """
    Compute the C6 dispersion coefficient for an atom pair.
    
    Uses a modified London formula based on atomic polarizabilities and HOMO energies:
        C6 = scale · (3/2) · α_i · α_j · ε_i · ε_j / (ε_i + ε_j)
    
    Parameters
    ----------
    ehomo_i, ehomo_j : float
        HOMO energies (Hartree, negative).
    alpha_i, alpha_j : float
        Polarizabilities (Bohr³).
    p : PFDParams
        Model parameters.
        
    Returns
    -------
    float
        C6 coefficient in Hartree·Bohr⁶.
    """
    return -p.c6_scale * 1.5 * ehomo_i * ehomo_j / (ehomo_i + ehomo_j) * alpha_i * alpha_j


def compute_rs2(ehomo_i: float, ehomo_j: float, z_i: int, z_j: int, p: PFDParams) -> float:
    """
    Compute the singularity-shift parameter Rs².
    
    This parameter modifies the R⁻⁶ term to R⁻⁶ → (R² - Rs²)⁻³, avoiding
    the singularity at R = 0.
    
    Parameters
    ----------
    ehomo_i, ehomo_j : float
        HOMO energies (Hartree, negative).
    z_i, z_j : int
        Atomic numbers.
    p : PFDParams
        Model parameters.
        
    Returns
    -------
    float
        Rs² parameter in Bohr².
    """
    rs2 = p.singularity_shift / (ehomo_i + ehomo_j)
    if z_i == 1 or z_j == 1:
        rs2 /= p.hydrogen_scale ** 2
    if (z_i == 1 and z_j == 6) or (z_i == 6 and z_j == 1):
        rs2 /= 4.0
    if z_i == 1 and z_j == 1:
        rs2 *= p.helium_factor * p.hydrogen_scale ** 2
    if z_i == 2 or z_j == 2:
        rs2 *= p.helium_factor
    return rs2


def compute_pair_params(ehomo_i: float, ehomo_j: float, alpha_i: float, alpha_j: float,
                        z_i: int, z_j: int, p: PFDParams = PARAMS) -> PairParams:
    """
    Compute all interaction parameters for an atom pair.
    
    This is a convenience function that bundles rd, C6, and Rs².
    
    Parameters
    ----------
    ehomo_i, ehomo_j : float
        HOMO energies (Hartree).
    alpha_i, alpha_j : float
        Polarizabilities (Bohr³).
    z_i, z_j : int
        Atomic numbers.
    p : PFDParams
        Model parameters.
        
    Returns
    -------
    PairParams
        Named tuple with (rd, c6, rs2).
    """
    return PairParams(compute_damping_radius(ehomo_i, ehomo_j, z_i, z_j, p),
                     compute_c6(ehomo_i, ehomo_j, alpha_i, alpha_j, p),
                     compute_rs2(ehomo_i, ehomo_j, z_i, z_j, p))


def compute_c9(ehomo_i: float, ehomo_j: float, ehomo_k: float,
               alpha_i: float, alpha_j: float, alpha_k: float) -> float:
    """
    Compute the C9 three-body dispersion coefficient.
    
    Uses the triple-dipole formula derived from third-order perturbation theory:
        C9 = (3/2) · ε_i·ε_j·ε_k · (ε_i+ε_j+ε_k) · α_i·α_j·α_k / [(ε_i+ε_j)(ε_i+ε_k)(ε_j+ε_k)]
    
    where ε are HOMO energies and α are polarizabilities.
    
    Parameters
    ----------
    ehomo_i, ehomo_j, ehomo_k : float
        HOMO energies of the three atoms (Hartree, negative).
    alpha_i, alpha_j, alpha_k : float
        Polarizabilities (Bohr³).
        
    Returns
    -------
    float
        C9 coefficient in Hartree·Bohr⁹.
    """
    num = 1.5 * ehomo_i * ehomo_j * ehomo_k * (ehomo_i + ehomo_j + ehomo_k)
    den = (ehomo_i + ehomo_j) * (ehomo_i + ehomo_k) * (ehomo_j + ehomo_k)
    return num * alpha_i * alpha_j * alpha_k / den


# =============================================================================
# Hydrogen Bond Functions
# =============================================================================

def hbond_angular_factor(cos_angle: float, p: PFDParams) -> Tuple[float, float]:
    """
    Compute the hydrogen bond angular correction factor.
    
    In H-X···Y systems, the dispersion between acceptors X and Y is modified
    by an angular-dependent factor based on the X-H-Y angle.
    
    The switching function is:
        f(θ) = strength · (1 + tanh(ρ)) / 2
    where ρ = sharpness · (θ - θ₀) / π
    
    Parameters
    ----------
    cos_angle : float
        Cosine of the X-H-Y angle.
    p : PFDParams
        Model parameters.
        
    Returns
    -------
    tuple of (float, float)
        (f, df/d(cos_angle)) - the angular factor and its derivative.
    """
    if abs(cos_angle) >= 1.0 - 1e-10:
        return 0.0, 0.0
    theta = np.arccos(cos_angle)
    rho = p.hbond_sharpness * (theta - p.hbond_angle_param * PI / 2.0) / PI
    tanh_rho = np.tanh(rho)
    f = p.hbond_strength * (1.0 + tanh_rho) / 2.0
    sin_theta = np.sqrt(1.0 - cos_angle**2)
    df = p.hbond_strength * p.hbond_sharpness / (2.0 * PI) * (1.0 - tanh_rho**2) * (-1.0 / sin_theta)
    return f, df


def identify_hbond_pattern(z_i: int, z_j: int, z_k: int, hbond_atoms: set) -> Tuple[str, int, int, int]:
    """
    Identify if a triplet forms an H-bond pattern H-X···Y.
    
    An H-bond pattern requires:
    - One hydrogen atom (the donor H)
    - Two electronegative atoms (the acceptors X and Y)
    
    Parameters
    ----------
    z_i, z_j, z_k : int
        Atomic numbers of the three atoms.
    hbond_atoms : set
        Set of atomic numbers that can act as H-bond acceptors.
        
    Returns
    -------
    tuple of (str or None, int, int, int)
        (h_vertex, h_idx, acc1_idx, acc2_idx) where:
        - h_vertex: 'i', 'j', or 'k' indicating which atom is H, or None if no H-bond
        - h_idx, acc1_idx, acc2_idx: indices (0, 1, 2) of H and acceptor atoms
    """
    if z_i == 1 and z_j in hbond_atoms and z_k in hbond_atoms:
        return 'i', 0, 1, 2
    if z_j == 1 and z_i in hbond_atoms and z_k in hbond_atoms:
        return 'j', 1, 0, 2
    if z_k == 1 and z_i in hbond_atoms and z_j in hbond_atoms:
        return 'k', 2, 0, 1
    return None, -1, -1, -1


# =============================================================================
# Iteration Utilities
# =============================================================================

def iter_pairs(num_atoms: int, z: np.ndarray, max_z: int = 115):
    """
    Iterate over valid atom pairs (i < j).
    
    Skips atoms with atomic number > max_z (e.g., dummy atoms).
    
    Parameters
    ----------
    num_atoms : int
        Total number of atoms.
    z : np.ndarray
        Atomic numbers.
    max_z : int, optional
        Maximum valid atomic number (default 115).
        
    Yields
    ------
    tuple of (int, int)
        Atom indices (i, j) with i < j.
    """
    for i in range(num_atoms - 1):
        if z[i] <= max_z:
            for j in range(i + 1, num_atoms):
                if z[j] <= max_z:
                    yield i, j


def iter_triplets(num_atoms: int, z: np.ndarray, max_z: int = 115):
    """
    Iterate over valid atom triplets (i < j < k).
    
    Skips atoms with atomic number > max_z.
    
    Parameters
    ----------
    num_atoms : int
        Total number of atoms.
    z : np.ndarray
        Atomic numbers.
    max_z : int, optional
        Maximum valid atomic number (default 115).
        
    Yields
    ------
    tuple of (int, int, int)
        Atom indices (i, j, k) with i < j < k.
    """
    for i in range(num_atoms - 2):
        if z[i] <= max_z:
            for j in range(i + 1, num_atoms - 1):
                if z[j] <= max_z:
                    for k in range(j + 1, num_atoms):
                        if z[k] <= max_z:
                            yield i, j, k


def compute_3body_hessian_terms(tri, dc, d2cos, c9, ep, disp3b, fe, fd, exp_fd, fs, f_dr,
                                 ed_r, ded_cos, edf, edt_r, dcos):
    """
    Compute the Hessian (second derivative) contributions from a single atom triplet.
    
    Mathematical Background
    -----------------------
    The three-body energy depends on distances and angles:
        E = -C9 · (3·cos_i·cos_j·cos_k + 1) · f²_ij·f²_ik·f²_jk / (r_ij·r_ik·r_jk)³
    
    The Hessian requires computing ∂²E/∂x_a∂x_b for all coordinate pairs.
    This is done in two stages:
    
    1. **Internal coordinate Hessian**: Compute the 3×3 matrix d²E/dr_m·dr_n
       where m,n ∈ {ij, ik, jk}. This requires second derivatives of:
       - The 1/r³ geometric factors
       - The f² damping factors  
       - The angular factor (3·cos_i·cos_j·cos_k + 1)
    
    2. **Cartesian transformation**: Convert to atomic coordinates using:
       ∂²E/∂x_a∂x_b = Σ_mn (∂²E/∂r_m∂r_n)·(∂r_m/∂x_a)·(∂r_n/∂x_b)
                    + Σ_m (∂E/∂r_m)·(∂²r_m/∂x_a∂x_b)
    
    The second term involves the "projection operator" P = (I - u⊗u)/r which
    accounts for how the distance changes when atoms move perpendicular to
    the bond direction.
    
    Implementation Details
    ----------------------
    The function computes three types of derivative matrices (all 3×3):
    
    - d_edf[m,n]: ∂(edf_m)/∂r_n where edf = E/(f²) is the damping contribution
    - d_ed_r[m,n]: ∂(ed_r_m)/∂r_n where ed_r = -3E/r is the geometric contribution
    - d_ded_cos[c,n]: ∂(ded_cos_c)/∂r_n where ded_cos is the angular contribution
    
    These are combined into d²E/dr_m·dr_n, then transformed to Cartesian Hessians.
    
    Parameters
    ----------
    tri : TriangleGeometry
        Precomputed triangle geometry (distances, angles, vectors).
    dc : CosineDerivatives
        First derivatives ∂cos/∂r.
    d2cos : dict
        Second derivatives ∂²cos/∂r_a∂r_b from compute_cosine_second_derivatives.
    c9 : float
        Three-body dispersion coefficient C9.
    ep : float
        Combined prefactor: (f_ij·f_ik·f_jk)² · (r_ij·r_ik·r_jk)⁻³
    disp3b : float
        Three-body dispersion energy magnitude: C9 · angular · ep
    fe : list of float
        Damping function values [f_ij, f_ik, f_jk].
    fd : list of float
        Scaled distance arguments [(r-2rd)/rd] for each edge.
    exp_fd : list of float
        Exponential terms exp(-2·fd) for each edge.
    fs : list of float
        Inverse damping radii [1/rd_ij, 1/rd_ik, 1/rd_jk].
    f_dr : list of float
        First derivatives d(f²)/dr for each edge.
    ed_r : list of float
        Energy derivatives from 1/r³ factors: -3·disp3b/r for each edge.
    ded_cos : list of float
        Energy derivatives from angular factor: ∂E/∂cos_c for c = i,j,k.
    edf : list of float
        Energy derivatives from damping: disp3b/f² for each edge.
    edt_r : list of float
        Total first derivatives dE/dr for each edge (sum of all contributions).
    dcos : list of list
        Matrix of ∂cos_c/∂r_m, indexed as dcos[c][m].
        
    Returns
    -------
    h_i, h_j, h_k : np.ndarray
        Diagonal Hessian blocks for atoms i, j, k, each shape (3, 3).
        These are the ∂²E/∂x_atom² contributions to be added to the full Hessian.
        
    Notes
    -----
    - The returned Hessians are negated because E_3body = -disp3b
    - Explicit variable names (d_ed_r00, d_ed_r01, ...) are used instead of
      loops for performance, as this function is called for every triplet
    - The projection operator P_m = (I - u_m⊗u_m)/r_m handles the geometric
      contribution from ∂²r/∂x_a∂x_b
    """
    r_inv0, r_inv1, r_inv2 = tri.r_ij_inv, tri.r_ik_inv, tri.r_jk_inv
    cos0, cos1, cos2 = tri.cos_i, tri.cos_j, tri.cos_k
    fe2_0, fe2_1, fe2_2 = fe[0]**2, fe[1]**2, fe[2]**2
    
    # d²(f²)/dr² for each edge (second derivatives of squared damping)
    t0 = 2.0*fs[0]*(1.0+2.0*fd[0])*exp_fd[0]
    t1 = 2.0*fs[1]*(1.0+2.0*fd[1])*exp_fd[1]
    t2 = 2.0*fs[2]*(1.0+2.0*fd[2])*exp_fd[2]
    d2f0 = t0*t0/2.0 - 8.0*fs[0]**2*fe[0]*fd[0]*exp_fd[0]
    d2f1 = t1*t1/2.0 - 8.0*fs[1]**2*fe[1]*fd[1]*exp_fd[1]
    d2f2 = t2*t2/2.0 - 8.0*fs[2]**2*fe[2]*fd[2]*exp_fd[2]
    
    # Precompute common terms
    f_dr_fe2_0, f_dr_fe2_1, f_dr_fe2_2 = f_dr[0]/fe2_0, f_dr[1]/fe2_1, f_dr[2]/fe2_2
    dcos_ded0 = ded_cos[0]*dcos[0][0] + ded_cos[1]*dcos[1][0] + ded_cos[2]*dcos[2][0]
    dcos_ded1 = ded_cos[0]*dcos[0][1] + ded_cos[1]*dcos[1][1] + ded_cos[2]*dcos[2][1]
    dcos_ded2 = ded_cos[0]*dcos[0][2] + ded_cos[1]*dcos[1][2] + ded_cos[2]*dcos[2][2]
    
    # d(edf)/dr - 3x3 matrix of derivatives of damping energy terms
    base00, base01, base02 = (ed_r[0]+dcos_ded0)/fe2_0, (ed_r[1]+dcos_ded1)/fe2_0, (ed_r[2]+dcos_ded2)/fe2_0
    base10, base11, base12 = (ed_r[0]+dcos_ded0)/fe2_1, (ed_r[1]+dcos_ded1)/fe2_1, (ed_r[2]+dcos_ded2)/fe2_1
    base20, base21, base22 = (ed_r[0]+dcos_ded0)/fe2_2, (ed_r[1]+dcos_ded1)/fe2_2, (ed_r[2]+dcos_ded2)/fe2_2
    d_edf00, d_edf01, d_edf02 = base00, base01 + edf[0]*f_dr_fe2_1, base02 + edf[0]*f_dr_fe2_2
    d_edf10, d_edf11, d_edf12 = base10 + edf[1]*f_dr_fe2_0, base11, base12 + edf[1]*f_dr_fe2_2
    d_edf20, d_edf21, d_edf22 = base20 + edf[2]*f_dr_fe2_0, base21 + edf[2]*f_dr_fe2_1, base22
    
    # d(ed_r)/dr - 3x3 matrix of derivatives of 1/r³ energy terms
    d_ed_r00 = -4.0*ed_r[0]*r_inv0 - 3.0*r_inv0*dcos_ded0 + ed_r[0]*f_dr_fe2_0
    d_ed_r01 = -3.0*ed_r[0]*r_inv1 - 3.0*r_inv0*dcos_ded1 + ed_r[0]*f_dr_fe2_1
    d_ed_r02 = -3.0*ed_r[0]*r_inv2 - 3.0*r_inv0*dcos_ded2 + ed_r[0]*f_dr_fe2_2
    d_ed_r10 = -3.0*ed_r[1]*r_inv0 - 3.0*r_inv1*dcos_ded0 + ed_r[1]*f_dr_fe2_0
    d_ed_r11 = -4.0*ed_r[1]*r_inv1 - 3.0*r_inv1*dcos_ded1 + ed_r[1]*f_dr_fe2_1
    d_ed_r12 = -3.0*ed_r[1]*r_inv2 - 3.0*r_inv1*dcos_ded2 + ed_r[1]*f_dr_fe2_2
    d_ed_r20 = -3.0*ed_r[2]*r_inv0 - 3.0*r_inv2*dcos_ded0 + ed_r[2]*f_dr_fe2_0
    d_ed_r21 = -3.0*ed_r[2]*r_inv1 - 3.0*r_inv2*dcos_ded1 + ed_r[2]*f_dr_fe2_1
    d_ed_r22 = -4.0*ed_r[2]*r_inv2 - 3.0*r_inv2*dcos_ded2 + ed_r[2]*f_dr_fe2_2
    
    # d(ded_cos)/dr - 3x3 matrix
    c9ep3 = 3.0*ep*c9
    # cos_cross: row c has products of other two cosines
    cc00, cc01, cc02 = 0, cos2, cos1
    cc10, cc11, cc12 = cos2, 0, cos0
    cc20, cc21, cc22 = cos1, cos0, 0
    
    d_ded_cos00 = -3.0*ded_cos[0]*r_inv0 + c9ep3*(cc01*dcos[1][0]+cc02*dcos[2][0]) + ded_cos[0]*f_dr_fe2_0
    d_ded_cos01 = -3.0*ded_cos[0]*r_inv1 + c9ep3*(cc01*dcos[1][1]+cc02*dcos[2][1]) + ded_cos[0]*f_dr_fe2_1
    d_ded_cos02 = -3.0*ded_cos[0]*r_inv2 + c9ep3*(cc01*dcos[1][2]+cc02*dcos[2][2]) + ded_cos[0]*f_dr_fe2_2
    d_ded_cos10 = -3.0*ded_cos[1]*r_inv0 + c9ep3*(cc10*dcos[0][0]+cc12*dcos[2][0]) + ded_cos[1]*f_dr_fe2_0
    d_ded_cos11 = -3.0*ded_cos[1]*r_inv1 + c9ep3*(cc10*dcos[0][1]+cc12*dcos[2][1]) + ded_cos[1]*f_dr_fe2_1
    d_ded_cos12 = -3.0*ded_cos[1]*r_inv2 + c9ep3*(cc10*dcos[0][2]+cc12*dcos[2][2]) + ded_cos[1]*f_dr_fe2_2
    d_ded_cos20 = -3.0*ded_cos[2]*r_inv0 + c9ep3*(cc20*dcos[0][0]+cc21*dcos[1][0]) + ded_cos[2]*f_dr_fe2_0
    d_ded_cos21 = -3.0*ded_cos[2]*r_inv1 + c9ep3*(cc20*dcos[0][1]+cc21*dcos[1][1]) + ded_cos[2]*f_dr_fe2_1
    d_ded_cos22 = -3.0*ded_cos[2]*r_inv2 + c9ep3*(cc20*dcos[0][2]+cc21*dcos[1][2]) + ded_cos[2]*f_dr_fe2_2
    
    # Build d²E/dr² with direct d2cos access
    d2E00 = (d_ed_r00 + f_dr[0]*d_edf00 + edf[0]*d2f0 +
             dcos[0][0]*d_ded_cos00 + dcos[1][0]*d_ded_cos10 + dcos[2][0]*d_ded_cos20 +
             ded_cos[0]*d2cos.get('d2cos_i_drij_drij',0) + ded_cos[1]*d2cos.get('d2cos_j_drij_drij',0) + ded_cos[2]*d2cos.get('d2cos_k_drij_drij',0))
    d2E01 = (d_ed_r01 + f_dr[0]*d_edf01 +
             dcos[0][0]*d_ded_cos01 + dcos[1][0]*d_ded_cos11 + dcos[2][0]*d_ded_cos21 +
             ded_cos[0]*d2cos.get('d2cos_i_drij_drik',0) + ded_cos[1]*d2cos.get('d2cos_j_drij_drik',0) + ded_cos[2]*d2cos.get('d2cos_k_drij_drik',0))
    d2E02 = (d_ed_r02 + f_dr[0]*d_edf02 +
             dcos[0][0]*d_ded_cos02 + dcos[1][0]*d_ded_cos12 + dcos[2][0]*d_ded_cos22 +
             ded_cos[0]*d2cos.get('d2cos_i_drij_drjk',0) + ded_cos[1]*d2cos.get('d2cos_j_drij_drjk',0) + ded_cos[2]*d2cos.get('d2cos_k_drij_drjk',0))
    d2E10 = (d_ed_r10 + f_dr[1]*d_edf10 +
             dcos[0][1]*d_ded_cos00 + dcos[1][1]*d_ded_cos10 + dcos[2][1]*d_ded_cos20 +
             ded_cos[0]*d2cos.get('d2cos_i_drik_drij',0) + ded_cos[1]*d2cos.get('d2cos_j_drik_drij',0) + ded_cos[2]*d2cos.get('d2cos_k_drik_drij',0))
    d2E11 = (d_ed_r11 + f_dr[1]*d_edf11 + edf[1]*d2f1 +
             dcos[0][1]*d_ded_cos01 + dcos[1][1]*d_ded_cos11 + dcos[2][1]*d_ded_cos21 +
             ded_cos[0]*d2cos.get('d2cos_i_drik_drik',0) + ded_cos[1]*d2cos.get('d2cos_j_drik_drik',0) + ded_cos[2]*d2cos.get('d2cos_k_drik_drik',0))
    d2E12 = (d_ed_r12 + f_dr[1]*d_edf12 +
             dcos[0][1]*d_ded_cos02 + dcos[1][1]*d_ded_cos12 + dcos[2][1]*d_ded_cos22 +
             ded_cos[0]*d2cos.get('d2cos_i_drik_drjk',0) + ded_cos[1]*d2cos.get('d2cos_j_drik_drjk',0) + ded_cos[2]*d2cos.get('d2cos_k_drik_drjk',0))
    d2E20 = (d_ed_r20 + f_dr[2]*d_edf20 +
             dcos[0][2]*d_ded_cos00 + dcos[1][2]*d_ded_cos10 + dcos[2][2]*d_ded_cos20 +
             ded_cos[0]*d2cos.get('d2cos_i_drjk_drij',0) + ded_cos[1]*d2cos.get('d2cos_j_drjk_drij',0) + ded_cos[2]*d2cos.get('d2cos_k_drjk_drij',0))
    d2E21 = (d_ed_r21 + f_dr[2]*d_edf21 +
             dcos[0][2]*d_ded_cos01 + dcos[1][2]*d_ded_cos11 + dcos[2][2]*d_ded_cos21 +
             ded_cos[0]*d2cos.get('d2cos_i_drjk_drik',0) + ded_cos[1]*d2cos.get('d2cos_j_drjk_drik',0) + ded_cos[2]*d2cos.get('d2cos_k_drjk_drik',0))
    d2E22 = (d_ed_r22 + f_dr[2]*d_edf22 + edf[2]*d2f2 +
             dcos[0][2]*d_ded_cos02 + dcos[1][2]*d_ded_cos12 + dcos[2][2]*d_ded_cos22 +
             ded_cos[0]*d2cos.get('d2cos_i_drjk_drjk',0) + ded_cos[1]*d2cos.get('d2cos_j_drjk_drjk',0) + ded_cos[2]*d2cos.get('d2cos_k_drjk_drjk',0))
    
    # Cartesian Hessian
    u0, u1, u2 = tri.vec_ij*r_inv0, tri.vec_ik*r_inv1, tri.vec_jk*r_inv2
    I3 = np.eye(3)
    P0, P1, P2 = (I3-np.outer(u0,u0))*r_inv0, (I3-np.outer(u1,u1))*r_inv1, (I3-np.outer(u2,u2))*r_inv2
    
    h_i = np.outer(u0,u0)*d2E00 + P0*edt_r[0] + np.outer(u1,u1)*d2E11 + P1*edt_r[1] + np.outer(u0,u1)*d2E01 + np.outer(u1,u0)*d2E10
    h_j = np.outer(u0,u0)*d2E00 + P0*edt_r[0] + np.outer(u2,u2)*d2E22 + P2*edt_r[2] - np.outer(u0,u2)*d2E02 - np.outer(u2,u0)*d2E20
    h_k = np.outer(u1,u1)*d2E11 + P1*edt_r[1] + np.outer(u2,u2)*d2E22 + P2*edt_r[2] + np.outer(u1,u2)*d2E12 + np.outer(u2,u1)*d2E21
    return -h_i, -h_j, -h_k
