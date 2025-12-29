#!/usr/bin/env python3
"""
PFD Dispersion Calculator - Command Line Interface
===================================================

This script provides a command-line interface for computing PFD dispersion
corrections from Gaussian output files.

Usage
-----
Basic calculation (energy + forces):
    python pfd_calculator.py molecule.out

Write output to file:
    python pfd_calculator.py molecule.out results.pfd

Energy only (faster):
    python pfd_calculator.py molecule.out --energy-only

Validate analytical gradients against numerical:
    python pfd_calculator.py molecule.out --numerical-grad

Full calculation with Hessians:
    python pfd_calculator.py molecule.out --full

Options
-------
--energy-only, -e     Compute only energies (no gradients)
--grad-only, -g       Compute energies and gradients (no Hessian)
--numerical-grad, -n  Compare analytical and numerical gradients
--numerical-hessian   Compare analytical and numerical Hessians
-q, --quiet           Suppress output (only write to file)

Output
------
Results are printed to stdout and optionally written to a file.
Energies are displayed in mEₕ.
Forces are in mEₕ/Bohr for two-body and three-body,
and nEₕ/Bohr for H-bond corrections.
Hessians are in mEₕ/Bohr².

Examples
--------
# Quick energy check
python pfd_calculator.py water.out -e

# Validate gradients for debugging
python pfd_calculator.py water.out --numerical-grad

# Production calculation with output file
python pfd_calculator.py protein.out protein.pfd
"""

import sys
import os
import argparse
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pfd import parse_gaussian, write_output, get_atomic_data, HARTREE_TO_MH
from pfd.energy import compute_pfd_energy
from pfd.gradients import compute_all_gradients
from pfd.hessian import compute_all_hessians


def _print_gradient_comparison(label: str, n: int, anal: np.ndarray, num: np.ndarray):
    """
    Print element-by-element gradient comparison table.
    
    Output is displayed in mEₕ/Bohr.
    
    Parameters
    ----------
    label : str
        Component label (e.g., "TWO-BODY").
    n : int
        Number of atoms.
    anal : np.ndarray
        Analytical forces in Eₕ/Bohr, shape (3, n).
    num : np.ndarray
        Numerical forces in Eₕ/Bohr, shape (3, n).
    """
    print(f"\n{label} GRADIENTS")
    print("-" * 70)
    print(f"{'Atom':>4} {'Comp':>4} {'Analytical':>15} {'Numerical':>15} {'Difference':>12}")
    print("-" * 70)
    for atom in range(n):
        for a, xyz in enumerate(['x', 'y', 'z']):
            av = anal[a, atom] * HARTREE_TO_MH
            nv = num[a, atom] * HARTREE_TO_MH
            print(f"{atom+1:>4} {xyz:>4} {av:15.9f} {nv:15.9f} {abs(av-nv):12.2e}")
    print(f"Max difference: {np.max(np.abs(anal - num)) * HARTREE_TO_MH:.2e} mH/Bohr")


def _print_hessian_comparison(label: str, n: int, anal: np.ndarray, num: np.ndarray):
    """
    Print Hessian comparison table for a single component.
    
    Output is displayed in mEₕ/Bohr².
    
    Parameters
    ----------
    label : str
        Component label (e.g., "TWO-BODY").
    n : int
        Number of atoms.
    anal : np.ndarray
        Analytical Hessian in Eₕ/Bohr², shape (3, 3, n).
    num : np.ndarray
        Numerical Hessian in Eₕ/Bohr², shape (3, 3, n).
    """
    print(f"\n{label} HESSIANS")
    print("-" * 70)
    print(f"{'Atom':>4} {'a,b':>4} {'Analytical':>15} {'Numerical':>15} {'Difference':>12}")
    print("-" * 70)
    for atom in range(n):
        for a in range(3):
            for b in range(3):
                av = anal[a, b, atom] * HARTREE_TO_MH
                nv = num[a, b, atom] * HARTREE_TO_MH
                ab = f"{'xyz'[a]}{'xyz'[b]}"
                print(f"{atom+1:>4} {ab:>4} {av:15.9f} {nv:15.9f} {abs(av-nv):12.2e}")
    print(f"Max difference: {np.max(np.abs(anal - num)) * HARTREE_TO_MH:.2e} mH/Bohr²")


def main():
    """
    Main entry point for the PFD command-line interface.
    
    Parses command-line arguments, reads molecular geometry from a Gaussian
    output file, computes PFD dispersion corrections at the requested level
    (energy, gradient, or Hessian), and writes results to an output file.
    
    Returns
    -------
    int
        Exit code (0 for success).
    """
    parser = argparse.ArgumentParser(description='PFD Dispersion Calculator')
    parser.add_argument('input', nargs='?', default='Gaussian.txt', help='Input Gaussian file')
    parser.add_argument('output', nargs='?', default='pfd.out', help='Output file')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--energy-only', action='store_true', help='Energy only')
    group.add_argument('--grad-only', action='store_true', help='Energy + gradients')
    parser.add_argument('--numerical-grad', action='store_true', help='Validate gradients')
    parser.add_argument('--numerical-hessian', action='store_true', help='Validate Hessians')
    parser.add_argument('-q', '--quiet', action='store_true', help='Suppress output')
    args = parser.parse_args()
    
    title, n, scf, z, coords = parse_gaussian(args.input)
    ehomo, alpha = get_atomic_data(z)
    
    if not args.quiet:
        print(f"PFD Calculator: {args.input} -> {args.output}")
    
    # Determine level and compute
    if args.energy_only:
        level = 'energy'
        e_2b, e_3b, e_hb, pe = compute_pfd_energy(n, z, coords, ehomo, alpha)
        f_2b = f_3b = f_hb = np.zeros((3, n))
        h_2b = h_3b = h_hb = np.zeros((3, 3, n))
    elif args.grad_only or args.numerical_grad:
        level = 'gradient'
        f_2b, f_3b, f_hb, e_2b, e_3b, e_hb, pe = compute_all_gradients(n, z, coords, ehomo, alpha)
        h_2b = h_3b = h_hb = np.zeros((3, 3, n))
    else:
        level = 'hessian'
        h_2b, h_3b, h_hb, f_2b, f_3b, f_hb, e_2b, e_3b, e_hb, pe = compute_all_hessians(n, z, coords, ehomo, alpha)
    
    result = {
        'energy_2body': e_2b, 'energy_3body': e_3b, 'energy_hbond': e_hb,
        'energy_total': e_2b + e_3b + e_hb,
        'forces_2body': f_2b, 'forces_3body': f_3b, 'forces_hbond': f_hb,
        'forces_total': f_2b + f_3b + f_hb,
        'hessian_2body': h_2b, 'hessian_3body': h_3b, 'hessian_hbond': h_hb,
        'hessian_total': h_2b + h_3b + h_hb,
        'pair_energies': pe,
    }
    
    if not args.quiet:
        print(f"  2-body: {e_2b*HARTREE_TO_MH:12.6f} mH")
        print(f"  3-body: {e_3b*HARTREE_TO_MH:12.6f} mH")
        print(f"  H-bond: {e_hb*HARTREE_TO_MH:12.6f} mH")
        print(f"  Total:  {(e_2b+e_3b+e_hb)*HARTREE_TO_MH:12.6f} mH")
    
    # Numerical gradient validation
    if args.numerical_grad:
        from pfd.numerical import (compute_2body_gradient_numerical, compute_3body_gradient_numerical,
                                   compute_hbond_gradient_numerical)
        print("\n" + "=" * 70)
        print("NUMERICAL GRADIENT VALIDATION")
        print("=" * 70)
        
        nf_2b, _, npe = compute_2body_gradient_numerical(n, z, coords, ehomo, alpha)
        nf_3b, _ = compute_3body_gradient_numerical(n, z, coords, ehomo, alpha)
        nf_hb, _ = compute_hbond_gradient_numerical(n, z, coords, ehomo, alpha, npe)
        
        _print_gradient_comparison("TWO-BODY", n, f_2b, nf_2b)
        _print_gradient_comparison("THREE-BODY", n, f_3b, nf_3b)
        _print_gradient_comparison("H-BOND", n, f_hb, nf_hb)
        
        nf_total = nf_2b + nf_3b + nf_hb
        print("\n" + "-" * 70)
        print("TOTAL GRADIENT SUMMARY")
        print("-" * 70)
        print(f"  2-body max diff:  {np.max(np.abs(f_2b - nf_2b)) * HARTREE_TO_MH:.2e} mH/Bohr")
        print(f"  3-body max diff:  {np.max(np.abs(f_3b - nf_3b)) * HARTREE_TO_MH:.2e} mH/Bohr")
        print(f"  H-bond max diff:  {np.max(np.abs(f_hb - nf_hb)) * HARTREE_TO_MH:.2e} mH/Bohr")
        print(f"  Total max diff:   {np.max(np.abs(result['forces_total'] - nf_total)) * HARTREE_TO_MH:.2e} mH/Bohr")
    
    # Numerical Hessian validation
    if args.numerical_hessian:
        from pfd.numerical import (compute_2body_hessian_numerical, compute_3body_hessian_numerical,
                                   compute_hbond_hessian_numerical)
        print("\n" + "=" * 70)
        print("NUMERICAL HESSIAN VALIDATION (from analytical gradients)")
        print("=" * 70)
        
        nh_2b, _, _, npe = compute_2body_hessian_numerical(n, z, coords, ehomo, alpha)
        nh_3b, _, _ = compute_3body_hessian_numerical(n, z, coords, ehomo, alpha)
        nh_hb, _, _ = compute_hbond_hessian_numerical(n, z, coords, ehomo, alpha, npe)
        
        _print_hessian_comparison("TWO-BODY", n, h_2b, nh_2b)
        _print_hessian_comparison("THREE-BODY", n, h_3b, nh_3b)
        _print_hessian_comparison("H-BOND", n, h_hb, nh_hb)
        
        nh_total = nh_2b + nh_3b + nh_hb
        print("\n" + "-" * 70)
        print("TOTAL HESSIAN SUMMARY")
        print("-" * 70)
        print(f"  2-body max diff:  {np.max(np.abs(h_2b - nh_2b)) * HARTREE_TO_MH:.2e} mH/Bohr²")
        print(f"  3-body max diff:  {np.max(np.abs(h_3b - nh_3b)) * HARTREE_TO_MH:.2e} mH/Bohr²")
        print(f"  H-bond max diff:  {np.max(np.abs(h_hb - nh_hb)) * HARTREE_TO_MH:.2e} mH/Bohr²")
        print(f"  Total max diff:   {np.max(np.abs(result['hessian_total'] - nh_total)) * HARTREE_TO_MH:.2e} mH/Bohr²")
    
    write_output(args.output, title, n, scf, z, coords, ehomo, alpha, result,
                 include_gradients=(level in ['gradient', 'hessian']),
                 include_hessians=(level == 'hessian'))
    
    if not args.quiet:
        print(f"\nOutput written to: {args.output}")
        if args.energy_only:
            print("  (energy only)")
        elif args.grad_only or args.numerical_grad:
            print("  (energy + gradients)")
        else:
            print("  (energy + gradients + Hessians)")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
