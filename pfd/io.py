"""
PFD I/O Module
==============

This module handles input parsing and output formatting for PFD calculations.

Input Parsing
-------------
The primary input format is Gaussian output files. The parse_gaussian function
extracts the molecular geometry and SCF energy from a Gaussian log file.

Output Formatting
-----------------
Results can be written to files in a format compatible with post-processing
tools. The output includes energies, forces, and Hessians with appropriate
unit conversions.

Data Structures
---------------
- PFDResult: Named tuple containing all calculation results
  - energies (2-body, 3-body, H-bond) in Hartree
  - forces as (3, n_atoms) arrays in Hartree/Bohr
  - Hessians as (3, 3, n_atoms) diagonal block arrays

Unit Conversions
----------------
- Coordinates: Ångströms (input) → Bohr (internal)
- Energies: Hartree → milliHartree (output display)
- Forces: Hartree/Bohr (internal)
"""

import numpy as np
from typing import Tuple, NamedTuple
import re
from .parameters import ANG_TO_BOHR, HARTREE_TO_MH, HARTREE_TO_NH


class PFDResult(NamedTuple):
    """
    Complete results from a PFD calculation.
    
    Attributes
    ----------
    energy_2body, energy_3body, energy_hbond : float
        Energy components in Hartree.
    forces_2body, forces_3body, forces_hbond : np.ndarray
        Force arrays, shape (3, n_atoms), in Hartree/Bohr.
    hessian_2body, hessian_3body, hessian_hbond : np.ndarray
        Hessian diagonal blocks, shape (3, 3, n_atoms).
    """
    energy_2body: float; energy_3body: float; energy_hbond: float
    forces_2body: np.ndarray; forces_3body: np.ndarray; forces_hbond: np.ndarray
    hessian_2body: np.ndarray; hessian_3body: np.ndarray; hessian_hbond: np.ndarray


def parse_gaussian(filename: str) -> Tuple[str, int, float, np.ndarray, np.ndarray]:
    """
    Parse molecular geometry from a Gaussian output file.
    
    Extracts:
    - Job title
    - Number of atoms
    - SCF energy
    - Atomic numbers
    - Cartesian coordinates (converted to Bohr)
    
    Parameters
    ----------
    filename : str
        Path to Gaussian output file.
        
    Returns
    -------
    tuple of (str, int, float, np.ndarray, np.ndarray)
        - title: Job title or "PFD Calculation"
        - num_atoms: Number of atoms
        - scf_energy: SCF energy in Hartree
        - atomic_numbers: Array of atomic numbers
        - coords: Coordinates in Bohr, shape (3, num_atoms)
        
    Notes
    -----
    Looks for "Standard orientation" or "Input orientation" geometry blocks.
    Coordinates are converted from Ångströms to Bohr.
    """
    with open(filename) as f:
        content = f.read()
    lines = content.split('\n')
    
    # Title
    title = "PFD Calculation"
    for i, line in enumerate(lines):
        if 'l101.exe' in line:
            for j in range(i + 1, min(i + 15, len(lines))):
                if '---' in lines[j] and len(lines[j].strip()) < 30:
                    title = lines[j + 1].strip() if j + 1 < len(lines) else title
                    break
            break
    
    # NAtoms
    num_atoms = 0
    for line in lines:
        m = re.search(r'NAtoms=\s*(\d+)', line)
        if m:
            num_atoms = int(m.group(1))
            break
    
    # Coordinates
    coords, atomic_nums = [], []
    for i, line in enumerate(lines):
        if 'Input orientation:' in line:
            for j in range(i + 5, i + 5 + num_atoms):
                p = lines[j].split()
                atomic_nums.append(int(p[1]))
                coords.append([float(p[3]), float(p[4]), float(p[5])])
            break
    
    # SCF energy
    scf_energy = 0.0
    for line in lines:
        if 'SCF Done:' in line:
            m = re.search(r'=\s*([-\d.]+)', line)
            if m:
                scf_energy = float(m.group(1))
                break
    
    return (title, num_atoms, scf_energy,
            np.array(atomic_nums, dtype=np.int32),
            np.array(coords).T * ANG_TO_BOHR)


def write_output(filename: str, title: str, num_atoms: int, scf_energy: float,
                 z: np.ndarray, coords: np.ndarray, ehomo: np.ndarray, alpha: np.ndarray,
                 result, include_gradients: bool = True, include_hessians: bool = True):
    """
    Write PFD calculation results to a formatted output file.
    
    The output format is compatible with Gaussian-style post-processing tools
    and includes molecular geometry, atomic parameters, energies, forces, and
    Hessians with appropriate unit conversions.
    
    Parameters
    ----------
    filename : str
        Path to output file.
    title : str
        Job title for the header.
    num_atoms : int
        Number of atoms.
    scf_energy : float
        SCF energy in Hartree (for reference).
    z : np.ndarray
        Atomic numbers, shape (num_atoms,).
    coords : np.ndarray
        Coordinates in Bohr, shape (3, num_atoms).
    ehomo : np.ndarray
        HOMO energies in Hartree, shape (num_atoms,).
    alpha : np.ndarray
        Polarizabilities in Bohr³, shape (num_atoms,).
    result : dict or PFDResult
        Calculation results containing energies, forces, and Hessians.
    include_gradients : bool, optional
        Whether to include force output. Default True.
    include_hessians : bool, optional
        Whether to include Hessian output. Default True.
        
    Notes
    -----
    Output units:
    - Coordinates: Ångströms
    - Energies: milliHartrees
    - Forces: milliHartrees/Bohr (or nanoHartrees/Bohr for H-bond)
    - Hessians: milliHartrees/Bohr²
    """
    # Extract from dict or PFDResult
    if isinstance(result, dict):
        e_2b, e_3b, e_hb = result['energy_2body'], result['energy_3body'], result['energy_hbond']
        f_2b, f_3b, f_hb = result['forces_2body'], result['forces_3body'], result['forces_hbond']
        h_2b, h_3b, h_hb = result['hessian_2body'], result['hessian_3body'], result['hessian_hbond']
    else:
        e_2b, e_3b, e_hb = result.energy_2body, result.energy_3body, result.energy_hbond
        f_2b, f_3b, f_hb = result.forces_2body, result.forces_3body, result.forces_hbond
        h_2b, h_3b, h_hb = result.hessian_2body, result.hessian_3body, result.hessian_hbond
    
    with open(filename, 'w') as f:
        # Header
        f.write(f" {title:<79}\n")
        f.write(f" NAtoms={num_atoms:12d}\n")
        f.write(f"   E(RAPFD) = {scf_energy:18.9f}\n")
        
        # Coordinates
        f.write("Atom Atomic Number    x           y           z          Ehomo       Alpha\n")
        coords_ang = coords / ANG_TO_BOHR
        for i in range(num_atoms):
            f.write(f"{i+1:3d}{z[i]:10d}{coords_ang[0,i]:12.6f}{coords_ang[1,i]:12.6f}"
                   f"{coords_ang[2,i]:12.6f}{ehomo[i]:12.8f}{alpha[i]:12.8f}\n")
        
        # Energies
        d2b, d3b, dhb = e_2b * HARTREE_TO_MH, e_3b * HARTREE_TO_MH, e_hb * HARTREE_TO_MH
        f.write(f"PFD-2B dispersion energy = {d2b:16.10f} milliHartrees.\n")
        f.write(f"PFD-3B dispersion energy = {d3b:16.10f} milliHartrees.\n")
        f.write(f"    Hydrogen Bond energy = {dhb:16.10f} milliHartrees.\n")
        f.write(f"Total correction to Epf  = {d2b+d3b+dhb:16.10f} milliHartrees.\n")
        
        if include_gradients:
            _write_forces(f, num_atoms, z, "Two-Body", f_2b, HARTREE_TO_MH, "milliHartrees/Bohr")
            _write_forces(f, num_atoms, z, "Three-Body", f_3b, HARTREE_TO_MH, "milliHartrees/Bohr")
            _write_forces(f, num_atoms, z, "H-Bonding", f_hb, HARTREE_TO_NH, "nanoHartrees/Bohr", exp=True)
        
        if include_hessians:
            _write_hessian(f, num_atoms, z, "Two-Body", h_2b, HARTREE_TO_MH)
            _write_hessian(f, num_atoms, z, "Three-Body", h_3b, HARTREE_TO_MH)
            _write_hessian(f, num_atoms, z, "H-Bonding", h_hb, HARTREE_TO_MH)


def _write_forces(f, n: int, z: np.ndarray, label: str, forces: np.ndarray, scale: float, unit: str, exp: bool = False):
    """Write force array to file in tabular format."""
    f.write(f"                   {label} Dispersion Forces\n")
    f.write("-" * 67 + "\n")
    f.write(f"  Center  Atomic               Forces ({unit})        \n")
    f.write("  Number  Number             X              Y              Z\n")
    f.write("-" * 67 + "\n")
    for i in range(n):
        if exp:
            f.write(f"{i+1:6d}{z[i]:8d}{forces[0,i]*scale:17.8E}{forces[1,i]*scale:17.8E}{forces[2,i]*scale:17.8E}\n")
        else:
            f.write(f"{i+1:6d}{z[i]:8d}{forces[0,i]*scale:15.9f}{forces[1,i]*scale:15.9f}{forces[2,i]*scale:15.9f}\n")
    f.write("-" * 67 + "\n")


def _write_hessian(f, n: int, z: np.ndarray, label: str, hess: np.ndarray, scale: float):
    """Write Hessian diagonal blocks to file in tabular format."""
    f.write(f"                   {label} Second Derivatives\n")
    f.write("-" * 67 + "\n")
    f.write("  Center  Atomic         Curvatures (milliHartrees/Bohr^2)        \n")
    f.write("  Number  Number            d/dx            d/dy          d/dz\n")
    f.write("-" * 67 + "\n")
    for i in range(n):
        f.write(f"{i+1:6d}{z[i]:9d}   d/dx{hess[0,0,i]*scale:15.9f}{hess[0,1,i]*scale:15.9f}{hess[0,2,i]*scale:15.9f}\n")
        f.write(f"                  d/dy{hess[1,0,i]*scale:15.9f}{hess[1,1,i]*scale:15.9f}{hess[1,2,i]*scale:15.9f}\n")
        f.write(f"                  d/dz{hess[2,0,i]*scale:15.9f}{hess[2,1,i]*scale:15.9f}{hess[2,2,i]*scale:15.9f}\n")
    f.write("-" * 67 + "\n")


def print_results(result: PFDResult):
    """
    Print PFD results summary to console.
    
    Parameters
    ----------
    result : PFDResult
        Calculation results.
        
    Notes
    -----
    Energies are displayed in milliHartrees (mH).
    """
    d2b = result.energy_2body * HARTREE_TO_MH
    d3b = result.energy_3body * HARTREE_TO_MH
    dhb = result.energy_hbond * HARTREE_TO_MH
    print(f"\nPFD Results:")
    print(f"  2-body: {d2b:16.10f} mH")
    print(f"  3-body: {d3b:16.10f} mH")
    print(f"  H-bond: {dhb:16.10f} mH")
    print(f"  Total:  {d2b+d3b+dhb:16.10f} mH")
