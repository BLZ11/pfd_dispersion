"""
PFD Parameters Module
=====================

This module defines the parameters, atomic data, and physical constants used in the
Petersson-Frisch Dispersion (PFD) correction for quantum chemical calculations.

Overview
--------
PFD is an empirical dispersion correction that accounts for long-range electron 
correlation effects (van der Waals/London dispersion forces) not captured by most
density functional theory (DFT) methods. The correction includes:

1. Two-body (pairwise) dispersion: C6/R^6 terms between atom pairs
2. Three-body (Axilrod-Teller-Muto) dispersion: C9 terms for atom triplets  
3. Hydrogen-bond corrections: Angular-dependent scaling for H-bonded systems

Key Parameters
--------------
The PFDParams dataclass contains fitted parameters that control:

- Damping function: Prevents divergence at short interatomic distances
  - damping_offset, damping_scale: Define the damping radius
  - singularity_shift: Shifts the R^6 denominator to avoid singularities

- Dispersion coefficients:
  - c6_scale: Global scaling of C6 coefficients
  - helium_factor: Special treatment for noble gas interactions
  - hydrogen_scale: Enhanced dispersion for hydrogen atoms

- Three-body terms:
  - threebody_damping_scale: Damping radius scaling for triplet interactions

- Hydrogen bonding:
  - hbond_strength, hbond_sharpness, hbond_angle_param: Control the angular
    dependence of dispersion in H-bonded systems

Atomic Data
-----------
- _EHOMO: HOMO (highest occupied molecular orbital) energies in Hartree
- _ALPHA: Static dipole polarizabilities in Bohr³

These atomic properties are used to compute C6 and C9 dispersion coefficients
via combination rules based on the London formula.

Units
-----
- Energies: Hartree (1 Hartree = 627.5 kcal/mol = 2625.5 kJ/mol)
- Distances: Bohr (1 Bohr = 0.529177 Å)
- Polarizabilities: Bohr³

References
----------
Petersson et al., J. Chem. Phys. (complete citation in main documentation)
"""

import numpy as np
from dataclasses import dataclass
from typing import FrozenSet

# =============================================================================
# Physical Constants
# =============================================================================

PI = np.pi
ANG_TO_BOHR = np.float64(1.8897261246257702)  # 1 Ångström = 1.8897... Bohr
BOHR_TO_ANG = 1.0 / ANG_TO_BOHR               # 1 Bohr = 0.5292... Ångström

# =============================================================================
# Unit Conversions for Output
# =============================================================================

HARTREE_TO_MH = np.float64(1000.0)      # Hartree to milliHartree
HARTREE_TO_NH = np.float64(1.0e6)       # Hartree to nanoHartree

# =============================================================================
# Hydrogen Bond Acceptor Atoms
# =============================================================================

# Atomic numbers of electronegative atoms that can act as H-bond acceptors
# N(7), O(8), F(9), P(15), S(16), Cl(17), As(33), Se(34), Br(35)
HBOND_ATOMS: FrozenSet[int] = frozenset({7, 8, 9, 15, 16, 17, 33, 34, 35})


# =============================================================================
# PFD Model Parameters
# =============================================================================

@dataclass(frozen=True)
class PFDParams:
    """
    Petersson-Frisch Dispersion (PFD) model parameters.
    
    These parameters were optimized to reproduce benchmark dispersion energies
    across a diverse set of molecular complexes.
    
    Attributes
    ----------
    damping_offset : float
        Offset in the damping radius calculation (Bohr). Controls where the
        damping function begins to attenuate the dispersion interaction.
        
    damping_scale : float
        Scaling factor for the damping radius. Larger values extend the
        damping region to longer distances.
        
    singularity_shift : float
        Shift applied to R² in the dispersion denominator to prevent
        singularities: E ~ C6/(R² - rs²)³ instead of E ~ C6/R⁶.
        
    c6_scale : float
        Global scaling factor for all C6 dispersion coefficients.
        
    helium_factor : float
        Special scaling for interactions involving helium atoms.
        
    hydrogen_scale : float
        Enhanced dispersion scaling for hydrogen atoms, which have
        characteristically different polarizability behavior.
        
    threebody_damping_scale : float
        Scaling factor for damping radii in three-body (ATM) terms.
        Typically 0.5 to account for the shorter-range nature of 3-body effects.
        
    hbond_strength : float
        Strength of the hydrogen bond angular correction. Negative value
        indicates reduction of acceptor-acceptor dispersion in H-bonds.
        
    hbond_sharpness : float
        Controls how rapidly the H-bond correction turns on/off with angle.
        
    hbond_angle_param : float
        Angular parameter for the H-bond correction function.
    """
    damping_offset: float = np.float64(-0.44382)
    damping_scale: float = np.float64(2.829285)
    singularity_shift: float = np.float64(-13.141740)
    c6_scale: float = np.float64(1.186044)
    helium_factor: float = np.float64(0.872259905)
    hydrogen_scale: float = np.float64(1.259456)
    threebody_damping_scale: float = np.float64(0.5)
    hbond_strength: float = -np.sqrt(3.0) / 2.0
    hbond_sharpness: float = np.float64(50.0)
    hbond_angle_param: float = np.float64(8.0 / 7.0)


# Default parameter set
PARAMS = PFDParams()

# =============================================================================
# Atomic HOMO Energies (Hartree)
# =============================================================================

# HOMO energies from atomic calculations, used in C6/C9 combination rules.
# These represent the ionization-related energy scale for each element.
_EHOMO = {
    1: -0.59277, 2: -0.91786, 3: -2.48668, 4: -0.83643, 5: -0.49781,
    6: -0.45579, 7: -0.57073, 8: -0.67794, 9: -0.43097, 10: -0.85114,
    11: -1.51908, 12: -0.31883, 13: -0.21790, 14: -0.27664, 15: -0.37275,
    16: -0.41589, 17: -0.25031, 18: -0.59132, 19: -0.72581, 20: -0.36045,
    21: -0.48732, 22: -0.47558, 23: -0.49997, 24: -0.52235, 25: -0.86849,
    26: -0.80746, 27: -0.78606, 28: -0.78761, 29: -0.81001, 30: -0.29256,
    31: -0.62872, 32: -0.29600, 33: -0.37030, 34: -0.41904, 35: -0.15934,
    36: -0.52431, 37: -0.81007, 38: -0.25846, 39: -0.43555, 40: -0.53274,
    41: -0.60845, 42: -0.61026, 43: -0.75892, 44: -0.89013, 45: -0.84578,
    46: -0.86934, 47: -1.38234, 48: -0.26485, 49: -0.19728, 50: -0.26504,
    51: -0.29447, 52: -0.31657, 53: -0.36359, 54: -0.45784, 55: -0.86183,
    56: -1.33004, 57: -0.68288, 58: -0.62554, 59: -0.63543, 60: -0.87969,
    61: -0.90269, 62: -0.95965, 63: -1.01869, 64: -0.68117, 65: -0.79650,
    66: -0.90495, 67: -0.83479, 68: -0.90156, 69: -0.97413, 70: -1.03326,
    71: -0.71549, 72: -0.78396, 73: -0.82838, 74: -0.56104, 75: -1.03841,
    76: -0.86250, 77: -0.95930, 78: -1.01755, 79: -0.76490, 80: -1.29467,
    81: -0.19274, 82: -0.25601, 83: -0.31982, 84: -0.19466, 85: -0.26009,
    86: -0.42799, 116: -1.29000, 117: -1.21000, 118: -1.21000,
}

# =============================================================================
# Atomic Polarizabilities (Bohr³)
# =============================================================================

# Static dipole polarizabilities, used in London formula for C6 coefficients.
# Larger polarizability = stronger dispersion interactions.
_ALPHA = {
    1: 1.9044, 2: 1.3643, 3: 0.1204, 4: 6.1686, 5: 7.4437,
    6: 8.7820, 7: 5.6129, 8: 4.8334, 9: 10.9428, 10: 2.7109,
    11: 0.4334, 12: 15.0282, 13: 14.6788, 14: 26.8194, 15: 23.8290,
    16: 19.8022, 17: 30.0468, 18: 10.9745, 19: 5.0994, 20: 15.1949,
    21: 15.7342, 22: 23.2528, 23: 16.1366, 24: 20.6438, 25: 24.5676,
    26: 15.0295, 27: 17.9818, 28: 16.3963, 29: 18.2705, 30: 19.2898,
    31: 23.4477, 32: 42.0000, 33: 28.3542, 34: 31.2402, 35: 84.5580,
    36: 16.9594, 37: 15.6479, 38: 18.6897, 39: 19.3530, 40: 34.3444,
    41: 19.8480, 42: 27.2140, 43: 30.2182, 44: 24.1099, 45: 22.1176,
    46: 21.2047, 47: 7.5669, 48: 23.7264, 49: 29.3879, 50: 30.1278,
    51: 41.9382, 52: 38.5203, 53: 33.9511, 54: 27.5774, 55: 15.646,
    56: 10.446, 57: 20.183, 58: 49.954, 59: 47.643, 60: 11.869,
    61: 10.867, 62: 9.961, 63: 9.187, 64: 36.798, 65: 13.407,
    66: 9.382, 67: 9.817, 68: 8.908, 69: 8.407, 70: 8.058,
    71: 12.960, 72: 26.792, 73: 24.930, 74: 34.563, 75: 11.186,
    76: 10.944, 77: 10.127, 78: 9.404, 79: 13.666, 80: 7.810,
    81: 69.360, 82: 58.663, 83: 47.969, 84: 45.274, 85: 39.229,
    86: 33.066, 116: -3.2000, 117: -2.1000, 118: -2.3000,
}


# =============================================================================
# Accessor Functions
# =============================================================================

def get_ehomo(z: int) -> float:
    """
    Get HOMO energy for an atom.
    
    Parameters
    ----------
    z : int
        Atomic number (1=H, 6=C, etc.)
        
    Returns
    -------
    float
        HOMO energy in Hartree. Returns -0.5 Ha for unknown elements.
    """
    return _EHOMO.get(z, -0.5)


def get_alpha(z: int) -> float:
    """
    Get static dipole polarizability for an atom.
    
    Parameters
    ----------
    z : int
        Atomic number (1=H, 6=C, etc.)
        
    Returns
    -------
    float
        Polarizability in Bohr³. Returns 10.0 Bohr³ for unknown elements.
    """
    return _ALPHA.get(z, 10.0)


def get_atomic_data(atomic_numbers: np.ndarray) -> tuple:
    """
    Get HOMO energies and polarizabilities for a set of atoms.
    
    This is the main interface for retrieving atomic data needed for
    PFD calculations.
    
    Parameters
    ----------
    atomic_numbers : np.ndarray
        Array of atomic numbers for all atoms in the molecule.
        
    Returns
    -------
    tuple of (np.ndarray, np.ndarray)
        - homo_energies: HOMO energies in Hartree, shape (n_atoms,)
        - polarizabilities: Polarizabilities in Bohr³, shape (n_atoms,)
        
    Example
    -------
    >>> z = np.array([6, 1, 1, 1, 1])  # Methane: C + 4H
    >>> ehomo, alpha = get_atomic_data(z)
    >>> print(ehomo)  # HOMO energies for C, H, H, H, H
    """
    return (np.array([get_ehomo(z) for z in atomic_numbers], dtype=np.float64),
            np.array([get_alpha(z) for z in atomic_numbers], dtype=np.float64))
