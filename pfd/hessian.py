"""
PFD Hessian Module
==================

This module computes analytical second derivatives (Hessians) of PFD dispersion.

The Hessian matrix ∂²E/∂x_i∂x_j is needed for:
- Vibrational frequency calculations
- Geometry optimization (Newton-Raphson)
- Thermodynamic property calculations

Implementation Notes
--------------------

The Hessian is stored in a compressed diagonal format: hess[a, b, i] represents
the (a,b) element of atom i's diagonal Hessian block. This captures the dominant
contributions efficiently.

For two-body terms:
    ∂²E/∂x_i∂x_j = (∂²E/∂r²) · (∂r/∂x_i)(∂r/∂x_j) + (∂E/∂r) · (∂²r/∂x_i∂x_j)

For three-body terms, the chain rule becomes significantly more complex,
involving derivatives with respect to all three distances and three angles.
The heavy lifting is done by compute_3body_hessian_terms in core.py.

Functions
---------
- compute_2body_hessian: Pairwise dispersion Hessian
- compute_3body_hessian: Three-body (ATM) Hessian
- compute_hbond_hessian: H-bond correction Hessian
- compute_all_hessians: Combined interface
"""

import numpy as np
from .parameters import PFDParams, PARAMS, HBOND_ATOMS
from .core import (compute_pair_params, compute_damping_radius, compute_c9, damping_function,
                   compute_triangle_geometry, compute_cosine_derivatives, compute_cosine_second_derivatives,
                   hbond_angular_factor, identify_hbond_pattern, iter_pairs, iter_triplets,
                   compute_3body_hessian_terms)


def compute_2body_hessian(n, z, coords, ehomo, alpha, p=PARAMS):
    """
    Compute two-body dispersion Hessian.
    
    Returns
    -------
    tuple of (hess, forces, energy, pair_e, pair_d1, pair_d2)
        - hess: Diagonal Hessian blocks, shape (3, 3, n_atoms)
        - forces: Gradient, shape (3, n_atoms)
        - energy: Two-body energy (Hartree)
        - pair_e: Pair energies matrix
        - pair_d1: First derivative factors (for H-bond Hessian)
        - pair_d2: Second derivative factors (for H-bond Hessian)
    """
    forces, hess, energy = np.zeros((3, n)), np.zeros((3, 3, n)), 0.0
    pair_e, pair_d1, pair_d2 = [np.zeros((n, n)) for _ in range(3)]
    
    for i, j in iter_pairs(n, z):
        d = coords[:, i] - coords[:, j]
        r = np.sqrt(d @ d)
        if r < 0.1: continue
        pp = compute_pair_params(ehomo[i], ehomo[j], alpha[i], alpha[j], z[i], z[j], p)
        if r <= 2.0 * pp.rd: continue
        df = damping_function(r, pp.rd)
        denom = r*r - pp.rs2
        if abs(denom) < 1e-10: continue
        
        r_inv, disp = 1.0/r, pp.c6/denom**3
        pair_e[i, j] = e_ij = -disp * df.f**2
        energy += e_ij
        
        ddisp = -6.0*r*disp/denom
        dE = ddisp*df.f**2 + 2.0*disp*df.f*df.df
        pair_d1[i, j] = ff = -r_inv*dE
        forces[:, i] -= d*ff; forces[:, j] += d*ff
        
        d2disp = (r_inv - 8.0*r/denom)*ddisp
        d2E = d2disp*df.f**2 + 4.0*df.f*ddisp*df.df + 2.0*disp*df.df**2 + 2.0*disp*df.f*df.d2f
        pair_d2[i, j] = hf = d2E*r_inv**2 - dE/r**3
        
        for a in range(3):
            for b in range(3):
                h_ab = -d[a]*d[b]*hf - (r_inv*dE if a == b else 0.0)
                hess[a,b,i] += h_ab; hess[a,b,j] += h_ab
    return hess, forces, energy, pair_e, pair_d1, pair_d2


def compute_3body_hessian(n, z, coords, ehomo, alpha, p=PARAMS):
    """
    Compute the three-body (Axilrod-Teller-Muto) dispersion Hessian.
    
    The three-body energy for triplet i-j-k is:
        E_ijk = -C9 · (3·cos_i·cos_j·cos_k + 1) · f²_ij·f²_ik·f²_jk / (r_ij·r_ik·r_jk)³
    
    The Hessian computation is complex because E depends on:
    - Three distances: r_ij, r_ik, r_jk
    - Three angles: θ_i, θ_j, θ_k (via their cosines)
    - Three damping functions: f_ij, f_ik, f_jk
    
    The chain rule gives:
        ∂²E/∂x_a∂x_b = Σ_m,n (∂²E/∂r_m∂r_n)·(∂r_m/∂x_a)·(∂r_n/∂x_b)
                     + Σ_m (∂E/∂r_m)·(∂²r_m/∂x_a∂x_b)
    
    The first term involves the 3×3 internal Hessian d²E/dr_m·dr_n.
    The second term involves projection operators for each edge.
    
    Parameters
    ----------
    n : int
        Number of atoms.
    z : np.ndarray
        Atomic numbers.
    coords : np.ndarray
        Coordinates in Bohr, shape (3, n).
    ehomo : np.ndarray
        HOMO energies in Hartree.
    alpha : np.ndarray
        Polarizabilities in Bohr³.
    p : PFDParams, optional
        Model parameters.
        
    Returns
    -------
    hess : np.ndarray
        Diagonal Hessian blocks, shape (3, 3, n).
    forces : np.ndarray
        Forces (gradient), shape (3, n).
    energy : float
        Three-body energy in Hartree.
        
    Notes
    -----
    The heavy computation is delegated to compute_3body_hessian_terms() in core.py,
    which handles the explicit chain-rule derivatives for a single triplet.
    This function loops over all triplets and accumulates their contributions.
    
    The intermediate quantities passed to compute_3body_hessian_terms are:
    - ed_r: ∂E/∂r contributions from the 1/r³ geometric factor
    - ded_cos: ∂E/∂cos contributions from the angular factor
    - edf: ∂E/∂(f²) contributions from damping
    - f_dr: d(f²)/dr damping derivatives
    - dcos: ∂cos/∂r matrix from law of cosines
    """
    forces, hess, energy = np.zeros((3, n)), np.zeros((3, 3, n)), 0.0
    scale = p.threebody_damping_scale
    
    for i, j, k in iter_triplets(n, z):
        tri = compute_triangle_geometry(coords, i, j, k)
        dc = compute_cosine_derivatives(tri)
        
        rd = [scale*compute_damping_radius(ehomo[a], ehomo[b], z[a], z[b], p) for a,b in [(i,j),(i,k),(j,k)]]
        damp = [damping_function(r, rd_) for r, rd_ in zip([tri.r_ij, tri.r_ik, tri.r_jk], rd)]
        fe = [d.f for d in damp]
        if any(abs(f) < 1e-15 for f in fe): continue
        
        c9 = compute_c9(ehomo[i], ehomo[j], ehomo[k], alpha[i], alpha[j], alpha[k])
        ep = (fe[0]*fe[1]*fe[2])**2 * (tri.r_ij_inv*tri.r_ik_inv*tri.r_jk_inv)**3
        disp3b = c9 * (3.0*tri.cos_i*tri.cos_j*tri.cos_k + 1.0) * ep
        energy -= disp3b
        
        ed_r = [-3.0*disp3b*r for r in [tri.r_ij_inv, tri.r_ik_inv, tri.r_jk_inv]]
        ded_cos = [c9*ep*3.0*c for c in [tri.cos_j*tri.cos_k, tri.cos_i*tri.cos_k, tri.cos_i*tri.cos_j]]
        edf = [disp3b/f**2 for f in fe]
        
        fs = [1.0/r for r in rd]
        fd = [(tri.r_ij-2*rd[0])/rd[0], (tri.r_ik-2*rd[1])/rd[1], (tri.r_jk-2*rd[2])/rd[2]]
        exp_fd = [np.exp(-2.0*x) for x in fd]
        f_dr = [2.0*fs[m]*fe[m]*(1.0+2.0*fd[m])*exp_fd[m] for m in range(3)]
        
        dcos = [[dc.dcos_i_drij, dc.dcos_i_drik, dc.dcos_i_drjk],
                [dc.dcos_j_drij, dc.dcos_j_drik, dc.dcos_j_drjk],
                [dc.dcos_k_drij, dc.dcos_k_drik, dc.dcos_k_drjk]]
        
        edt_r = [ed_r[m] + edf[m]*f_dr[m] + sum(ded_cos[c]*dcos[c][m] for c in range(3)) for m in range(3)]
        
        forces[:, i] += tri.vec_ij*tri.r_ij_inv*edt_r[0] + tri.vec_ik*tri.r_ik_inv*edt_r[1]
        forces[:, j] += -tri.vec_ij*tri.r_ij_inv*edt_r[0] + tri.vec_jk*tri.r_jk_inv*edt_r[2]
        forces[:, k] += -tri.vec_ik*tri.r_ik_inv*edt_r[1] - tri.vec_jk*tri.r_jk_inv*edt_r[2]
        
        d2cos = compute_cosine_second_derivatives(tri, dc)
        h_i, h_j, h_k = compute_3body_hessian_terms(tri, dc, d2cos, c9, ep, disp3b, fe, fd, exp_fd, fs, f_dr, ed_r, ded_cos, edf, edt_r, dcos)
        hess[:,:,i] += h_i; hess[:,:,j] += h_j; hess[:,:,k] += h_k
    return hess, forces, energy


def compute_hbond_hessian(n, z, coords, ehomo, alpha, pair_e, pair_d1, pair_d2, p=PARAMS):
    """
    Compute the hydrogen bond correction Hessian.
    
    The H-bond energy correction is:
        E_hb = E_XY · g(θ)
    
    where E_XY is the X-Y pair dispersion and g(θ) is an angular switching function.
    
    The Hessian has four contributions:
    1. g · ∂²E_XY/∂x²  (scaled pair Hessian)
    2. E_XY · (∂²g/∂cos²) · (∂cos/∂x)⊗(∂cos/∂x)  (angular 2nd derivative)
    3. E_XY · (∂g/∂cos) · (∂²cos/∂x²)  (cosine 2nd derivative)  
    4. Cross terms between pair force and angular gradient
    
    Returns
    -------
    hess : np.ndarray
        Diagonal Hessian blocks, shape (3, 3, n).
    forces : np.ndarray
        Forces, shape (3, n).
    energy : float
        H-bond correction energy in Hartree.
    """
    forces, hess, energy = np.zeros((3, n)), np.zeros((3, 3, n)), 0.0
    
    for i, j, k in iter_triplets(n, z):
        h_vtx, _, a1_rel, a2_rel = identify_hbond_pattern(z[i], z[j], z[k], HBOND_ATOMS)
        if h_vtx is None:
            continue
        
        idx = {0: i, 1: j, 2: k}
        a1, a2 = idx[a1_rel], idx[a2_rel]
        
        tri = compute_triangle_geometry(coords, i, j, k)
        dc = compute_cosine_derivatives(tri)
        d2cos = compute_cosine_second_derivatives(tri, dc)
        
        # Identify angle at hydrogen and acceptor-acceptor geometry
        if h_vtx == 'i':
            cos_h = tri.cos_i
            dcos_drij, dcos_drik, dcos_drjk = dc.dcos_i_drij, dc.dcos_i_drik, dc.dcos_i_drjk
            acc_vec, acc_r_inv = tri.vec_jk, tri.r_jk_inv
            acc_pair = (j, k)
            d2c_rij_rij = d2cos['d2cos_i_drij_drij']
            d2c_rij_rik = d2cos['d2cos_i_drij_drik']
            d2c_rij_rjk = d2cos['d2cos_i_drij_drjk']
            d2c_rik_rik = d2cos['d2cos_i_drik_drik']
            d2c_rik_rjk = d2cos['d2cos_i_drik_drjk']
            d2c_rjk_rjk = d2cos['d2cos_i_drjk_drjk']
        elif h_vtx == 'j':
            cos_h = tri.cos_j
            dcos_drij, dcos_drik, dcos_drjk = dc.dcos_j_drij, dc.dcos_j_drik, dc.dcos_j_drjk
            acc_vec, acc_r_inv = tri.vec_ik, tri.r_ik_inv
            acc_pair = (i, k)
            d2c_rij_rij = d2cos['d2cos_j_drij_drij']
            d2c_rij_rik = d2cos['d2cos_j_drij_drik']
            d2c_rij_rjk = d2cos['d2cos_j_drij_drjk']
            d2c_rik_rik = d2cos['d2cos_j_drik_drik']
            d2c_rik_rjk = d2cos['d2cos_j_drik_drjk']
            d2c_rjk_rjk = d2cos['d2cos_j_drjk_drjk']
        else:
            cos_h = tri.cos_k
            dcos_drij, dcos_drik, dcos_drjk = dc.dcos_k_drij, dc.dcos_k_drik, dc.dcos_k_drjk
            acc_vec, acc_r_inv = tri.vec_ij, tri.r_ij_inv
            acc_pair = (i, j)
            d2c_rij_rij = d2cos['d2cos_k_drij_drij']
            d2c_rij_rik = d2cos['d2cos_k_drij_drik']
            d2c_rij_rjk = d2cos['d2cos_k_drij_drjk']
            d2c_rik_rik = d2cos['d2cos_k_drik_drik']
            d2c_rik_rjk = d2cos['d2cos_k_drik_drjk']
            d2c_rjk_rjk = d2cos['d2cos_k_drjk_drjk']
        
        lo, hi = (a1, a2) if a1 < a2 else (a2, a1)
        pe, ff, hf = pair_e[lo, hi], pair_d1[lo, hi], pair_d2[lo, hi]
        
        # Get angular factor with BOTH first and second derivatives
        hb_g, hb_dg, hb_d2g = hbond_angular_factor(cos_h, p)
        
        if z[a1] == 34 or z[a2] == 34:
            hb_g, hb_dg, hb_d2g = hb_g/2.0, hb_dg/2.0, hb_d2g/2.0
        
        energy += pe * hb_g
        
        # Cartesian derivatives of cos(θ_H)
        dcos_xi = tri.vec_ij*tri.r_ij_inv*dcos_drij + tri.vec_ik*tri.r_ik_inv*dcos_drik
        dcos_xj = -tri.vec_ij*tri.r_ij_inv*dcos_drij + tri.vec_jk*tri.r_jk_inv*dcos_drjk
        dcos_xk = -tri.vec_ik*tri.r_ik_inv*dcos_drik - tri.vec_jk*tri.r_jk_inv*dcos_drjk
        
        # Forces
        if h_vtx == 'i':
            forces[:, j] -= hb_g*ff*acc_vec + pe*hb_dg*dcos_xj
            forces[:, k] -= -hb_g*ff*acc_vec + pe*hb_dg*dcos_xk
            forces[:, i] -= pe*hb_dg*dcos_xi
        elif h_vtx == 'j':
            forces[:, i] -= hb_g*ff*acc_vec + pe*hb_dg*dcos_xi
            forces[:, k] -= -hb_g*ff*acc_vec + pe*hb_dg*dcos_xk
            forces[:, j] -= pe*hb_dg*dcos_xj
        else:
            forces[:, i] -= hb_g*ff*acc_vec + pe*hb_dg*dcos_xi
            forces[:, j] -= -hb_g*ff*acc_vec + pe*hb_dg*dcos_xj
            forces[:, k] -= pe*hb_dg*dcos_xk
        
        # ===== HESSIAN TERMS =====
        
        # Term 1: g · ∂²E_XY/∂x² (scaled pair Hessian on acceptors)
        # Note: ff = -dE/r, and diagonal term in 2-body is -(r_inv*dE) = +ff
        for a in range(3):
            for b in range(3):
                h_ab = -acc_vec[a]*acc_vec[b]*hf + (ff if a == b else 0.0)
                hess[a, b, acc_pair[0]] += hb_g * h_ab
                hess[a, b, acc_pair[1]] += hb_g * h_ab
        
        # Term 2: E_XY · (∂²g/∂cos²) · (∂cos/∂x)⊗(∂cos/∂x)
        # NOTE: Using hb_d2g (second derivative), not hb_dg!
        for m, dcos_xm in [(i, dcos_xi), (j, dcos_xj), (k, dcos_xk)]:
            hess[:, :, m] += pe * hb_d2g * np.outer(dcos_xm, dcos_xm)
        
        # Term 3: E_XY · (∂g/∂cos) · (∂²cos/∂x²)
        u_ij = tri.vec_ij * tri.r_ij_inv
        u_ik = tri.vec_ik * tri.r_ik_inv
        u_jk = tri.vec_jk * tri.r_jk_inv
        I3 = np.eye(3)
        P_ij = (I3 - np.outer(u_ij, u_ij)) * tri.r_ij_inv
        P_ik = (I3 - np.outer(u_ik, u_ik)) * tri.r_ik_inv
        P_jk = (I3 - np.outer(u_jk, u_jk)) * tri.r_jk_inv
        
        # ∂²cos/∂x_i² (atom i connects to j via r_ij, to k via r_ik)
        d2cos_xi_xi = (np.outer(u_ij, u_ij) * d2c_rij_rij +
                       np.outer(u_ik, u_ik) * d2c_rik_rik +
                       (np.outer(u_ij, u_ik) + np.outer(u_ik, u_ij)) * d2c_rij_rik +
                       P_ij * dcos_drij + P_ik * dcos_drik)
        
        # ∂²cos/∂x_j² (atom j connects to i via -r_ij, to k via r_jk)
        d2cos_xj_xj = (np.outer(u_ij, u_ij) * d2c_rij_rij +
                       np.outer(u_jk, u_jk) * d2c_rjk_rjk -
                       (np.outer(u_ij, u_jk) + np.outer(u_jk, u_ij)) * d2c_rij_rjk +
                       P_ij * dcos_drij + P_jk * dcos_drjk)
        
        # ∂²cos/∂x_k² (atom k connects to i via -r_ik, to j via -r_jk)
        d2cos_xk_xk = (np.outer(u_ik, u_ik) * d2c_rik_rik +
                       np.outer(u_jk, u_jk) * d2c_rjk_rjk +
                       (np.outer(u_ik, u_jk) + np.outer(u_jk, u_ik)) * d2c_rik_rjk +
                       P_ik * dcos_drik + P_jk * dcos_drjk)
        
        hess[:, :, i] += pe * hb_dg * d2cos_xi_xi
        hess[:, :, j] += pe * hb_dg * d2cos_xj_xj
        hess[:, :, k] += pe * hb_dg * d2cos_xk_xk
        
        # Term 4: Cross terms 2·(∂g/∂cos)·(∂cos/∂x)·(∂E_XY/∂x) on acceptors
        dE_acc0 = ff * acc_vec
        dE_acc1 = -ff * acc_vec
        
        if h_vtx == 'i':
            dcos_acc0, dcos_acc1 = dcos_xj, dcos_xk
        elif h_vtx == 'j':
            dcos_acc0, dcos_acc1 = dcos_xi, dcos_xk
        else:
            dcos_acc0, dcos_acc1 = dcos_xi, dcos_xj
        
        hess[:, :, acc_pair[0]] += hb_dg * (np.outer(dcos_acc0, dE_acc0) + np.outer(dE_acc0, dcos_acc0))
        hess[:, :, acc_pair[1]] += hb_dg * (np.outer(dcos_acc1, dE_acc1) + np.outer(dE_acc1, dcos_acc1))
    
    return hess, forces, energy


def compute_all_hessians(n, z, coords, ehomo, alpha, p=PARAMS):
    """
    Compute all PFD Hessian components.
    
    This is the main interface for Hessian calculations, computing second
    derivatives for two-body, three-body, and H-bond dispersion terms.
    
    Parameters
    ----------
    n : int
        Number of atoms.
    z : np.ndarray
        Atomic numbers, shape (n,).
    coords : np.ndarray
        Coordinates in Bohr, shape (3, n).
    ehomo : np.ndarray
        HOMO energies in Hartree, shape (n,).
    alpha : np.ndarray
        Polarizabilities in Bohr³, shape (n,).
    p : PFDParams, optional
        Model parameters.
        
    Returns
    -------
    h_2body : np.ndarray
        Two-body Hessian, shape (3, 3, n).
    h_3body : np.ndarray
        Three-body Hessian, shape (3, 3, n).
    h_hbond : np.ndarray
        H-bond Hessian, shape (3, 3, n).
    f_2body : np.ndarray
        Two-body forces, shape (3, n).
    f_3body : np.ndarray
        Three-body forces, shape (3, n).
    f_hbond : np.ndarray
        H-bond forces, shape (3, n).
    e_2body : float
        Two-body energy in Hartree.
    e_3body : float
        Three-body energy in Hartree.
    e_hbond : float
        H-bond energy in Hartree.
    pair_energies : np.ndarray
        Matrix of pairwise energies.
        
    Example
    -------
    >>> h_2b, h_3b, h_hb, f_2b, f_3b, f_hb, e_2b, e_3b, e_hb, pe = compute_all_hessians(n, z, coords, ehomo, alpha)
    >>> h_total = h_2b + h_3b + h_hb
    >>> f_total = f_2b + f_3b + f_hb
    >>> e_total = e_2b + e_3b + e_hb
    
    Notes
    -----
    The two-body calculation is performed first because it generates pair_d1
    and pair_d2 arrays that are reused by the H-bond Hessian calculation.
    This avoids redundant computation of pair derivatives.
    """
    h_2b, f_2b, e_2b, pair_e, pair_d1, pair_d2 = compute_2body_hessian(n, z, coords, ehomo, alpha, p)
    h_3b, f_3b, e_3b = compute_3body_hessian(n, z, coords, ehomo, alpha, p)
    h_hb, f_hb, e_hb = compute_hbond_hessian(n, z, coords, ehomo, alpha, pair_e, pair_d1, pair_d2, p)
    return h_2b, h_3b, h_hb, f_2b, f_3b, f_hb, e_2b, e_3b, e_hb, pair_e
