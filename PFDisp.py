# MIT License
# 
# Copyright (c) 2024, Barbaro Zulueta
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""Compute two-body and three-body Petersson-Frisch dispersion model and hydrogen bonding corrections"""
import numpy as np
from parameters import EhomoDat, AlphaDat, HbndDat,PFD3Bparam

def PFD_2B(atoms,atom_number, Distance_Matrix):
    """Compute the Petersson-Frisch two-body dispersion
       correction.
       
       Parameters
       ----------
       atoms: obj:'np.ndarray'
           the atomic elements
       atom_number: obj:'np.ndarray'
           the atomic number
       Distance_Matrix: obj:'np.ndarray'
           a matrix containing the atom-pair distance (in Angstroms)
       
       Returns
       -------
       E_disp_2B: obj:'np.ndarray'
           two-body dispersion energy (in Hartrees) 
    """
    
    size = atoms.shape[0]
    Distance_Matrix = Distance_Matrix / 0.52917721092 # convert distances from Angstroms to Bohr
    Disp2B = np.zeros((size,size)) # initialize array
    E_disp_2B = 0
    fR0, fRd, fRs, fC6, fHe, fZH, fRd6b = PFD3Bparam('two-body')
    
    for i in range(1,size):
        AlphaI = AlphaDat(atoms[i])
        EhomoI = EhomoDat(atoms[i])
        Ihal = 0
        if (atom_number[i] == '116') or (atom_number[i] == '117') or (atom_number[i] == '118'):
            Ihal = 1
        for j in range(i):
            Jhal = 0
            if (atom_number[j] == '116') or (atom_number[j] == '117') or (atom_number[j] == '118'):
                Jhal = 1
            if (atom_number[i] == '9') and (atom_number[j] == '115'):
                continue
            if (atom_number[i] == '17') and (atom_number[j] == '115'):
                continue
            if (atom_number[i] == '35') and (atom_number[j] == '115'):
                continue
            Rij = Distance_Matrix[i][j]
            if Rij < 0.1:
                continue
            RDij = fRd / (np.sqrt(-(EhomoDat(atoms[i]) + EhomoDat(atoms[j])))) + fR0
            if (Ihal == 1) and (Jhal == 1):
                RDij = 1.5 * RDij
            
            # Modify Rd for 1s orbitals
            fHe2 = np.sqrt(fHe)
            if (atom_number[i] == '1') or (atom_number[j] == '1'): 
                RDij = fHe2 * RDij
            if (atom_number[i] == '1') and (atom_number[j] == '6'): 
                RDij = RDij / (fZH * fHe2)
            if (atom_number[i] == '6') and (atom_number[j] == '1'): 
                RDij = RDij / (fZH * fHe2)
            
            EhIJ = EhomoDat(atoms[i]) * EhomoDat(atoms[j]) / (EhomoDat(atoms[i]) + EhomoDat(atoms[j]))
            C6IJ = -fC6 * (3 / 2) * EhIJ * AlphaDat(atoms[i]) * AlphaDat(atoms[j]) 
            RSIJ2 = fRs / (EhomoDat(atoms[i]) + EhomoDat(atoms[j]))
    
            if ((atom_number[i] == '9') and (atom_number[j] == '116')) or ((atom_number[i] == '17') and (atom_number[j] == '117')) or ((atom_number[i] == '35') and (atom_number[j] == '118')):
                RijInv = 0
            else:
                RijInv = 1 / Rij
            Fs1 = 1 / (fRd6b * RDij)
            
            # Modify Rs for 1s orbitals
            if (atom_number[i] == '1') or (atom_number[j] == '1'):
                RSIJ2 = RSIJ2 / (fZH * fZH) 
            if (atom_number[i] == '1') and (atom_number[j] == '6'):
                RSIJ2 = RSIJ2 / 4
            if (atom_number[i] == '6') and (atom_number[j] == '1'):
                RSIJ2 = RSIJ2 / 4
            if (atom_number[i] == '1') and (atom_number[j] == '1'):
                RSIJ2 = fHe * fZH * fZH * RSIJ2
            if (atom_number[i] == '2') or (atom_number[j] == '2'):
                RSIJ2 = fHe * RSIJ2  
            FD = (Rij - 2 * RDij) / RDij
            FS = C6IJ / (Rij**2 - RSIJ2)**3 
            FE = 1 - (1 + FD) * np.exp(-2 * FD)
                           
            # Compute PFD-2 Body
            # Note E_disp_2B = 0, if Rij < 2.0*RDIJ
            if (Rij > 2 * RDij):
                Disp2B[i][j] = -FS * (FE**2)
                E_disp_2B += Disp2B[i][j]
            else:
                continue
                           
    return (E_disp_2B, Disp2B) # in Hartrees
                           
def PFD_3B(atoms,atom_number, Distance_Matrix):
    """Compute the Petersson-Frisch three-body dispersion
       correction.
       
       Parameters
       ----------
       atoms: obj:'np.ndarray'
           the atomic elements
       atom_number: obj:'np.ndarray'
           the atomic number
       Distance_Matrix: obj:'np.ndarray'
           a matrix containing the atom-pair distance (in Ansgtroms)
       
       Returns
       -------
       E_disp_3B: obj:'np.ndarray'
           three-body dispersion energy (in Hartrees)
    """
    E_disp_3B = 0
    size = atoms.shape[0]
    Distance_Matrix = Distance_Matrix / 0.52917721092 # convert distances from Angstroms to Bohr
    fR0, fRd, fRs, fC6, fHe, fZH, fRd6b = PFD3Bparam('three-body')
    
    for i in range(size):
        if np.int(atom_number[i]) > 115:
            continue
        AlphaI = AlphaDat(atoms[i])
        EhomoI = EhomoDat(atoms[i])
        for j in range(i):
            if np.int(atom_number[j]) > 115:
                continue
            Rij = Distance_Matrix[i][j]
            RDij = fRd / np.sqrt(-(EhomoDat(atoms[i]) + EhomoDat(atoms[j]))) + fR0
            
            # Modify Rd for 1s orbitals
            fHe2 = np.sqrt(fHe)
            if (atom_number[i] == '1') or (atom_number[j] == '1'): 
                RDij = fHe2 * RDij
            if (atom_number[i] == '1') and (atom_number[j] == '6'): 
                RDij = RDij / (fZH * fHe2)
            if (atom_number[i] == '6') and (atom_number[j] == '1'): 
                RDij = RDij / (fZH * fHe2)
            
            EhIJ = EhomoDat(atoms[i]) * EhomoDat(atoms[j]) / (EhomoDat(atoms[i]) + EhomoDat(atoms[j]))
            C6IJ = -fC6 * (3 / 2) * EhIJ * AlphaDat(atoms[i]) * AlphaDat(atoms[j]) 
            RSIJ2 = fRs / (EhomoDat(atoms[i]) + EhomoDat(atoms[j]))
            RijInv = 1 / Rij
            Fs1 = 1 / (fRd6b * RDij)  
            
            # Modify Rs for 1s orbitals
            if (atom_number[i] == '1') or (atom_number[j] == '1'):
                RSIJ2 = RSIJ2 / (fZH * fZH) 
            if (atom_number[i] == '1') and (atom_number[j] == '6'):
                RSIJ2 = RSIJ2 / 4
            if (atom_number[i] == '6') and (atom_number[j] == '1'):
                RSIJ2 = RSIJ2 / 4
            if (atom_number[i] == '1') and (atom_number[j] == '1'):
                RSIJ2 = fHe * fZH * fZH * RSIJ2
            if (atom_number[i] == '2') or (atom_number[j] == '2'):
                RSIJ2 = fHe * RSIJ2  
            FD = (Rij - 2 * RDij) / RDij
            FS = C6IJ / (Rij**2 - RSIJ2)**3 
            FE = 1 - (1 + FD) * np.exp(-2 * FD)
            
            # Calculate 3-body dispersion terms of i and j with each unique third atom, k
            for k in range(j):
                if atom_number[k] == '115':
                    continue
                # Compute atom-pair distances for three-body interactions
                Rik = Distance_Matrix[i][k]
                Rjk = Distance_Matrix[j][k]
                # Calculate two-body damping radii
                RDik = fRd / np.sqrt(-(EhomoDat(atoms[i]) + EhomoDat(atoms[k]))) + fR0
                RDjk = fRd / np.sqrt(-(EhomoDat(atoms[j]) + EhomoDat(atoms[k]))) + fR0
                
                # Modify Rd for 1s orbitals
                ## Rik corrections
                if (atom_number[i] == '1') or (atom_number[k] == '1'): 
                    RDik = fHe2 * RDik
                if (atom_number[i] == '1') and (atom_number[k] == '6'): 
                    RDik = RDik / (fZH * fHe2)
                if (atom_number[i] == '6') and (atom_number[k] == '1'): 
                    RDik = RDik / (fZH * fHe2)
                EhIK = EhomoDat(atoms[i]) * EhomoDat(atoms[k]) / (EhomoDat(atoms[i]) + EhomoDat(atoms[k]))
                C6IK = -fC6 * (3 / 2) * EhIK * AlphaDat(atoms[i]) * AlphaDat(atoms[k]) 
                FD2 = (Rik - (2 * fRd6b * RDik)) / (fRd6b * RDik)
                Fs2 = 1 / (fRd6b * RDik)
                RikInv = 1 / Rik  
                
                ## Rjk corrections
                if (atom_number[j] == '1') or (atom_number[k] == '1'): 
                    RDjk = fHe2 * RDjk
                if (atom_number[j] == '1') and (atom_number[k] == '6'): 
                    RDjk = RDjk / (fZH * fHe2)
                if (atom_number[j] == '6') and (atom_number[k] == '1'): 
                    RDjk = RDjk / (fZH * fHe2)
                EhJK = EhomoDat(atoms[j]) * EhomoDat(atoms[k]) / (EhomoDat(atoms[j]) + EhomoDat(atoms[k]))
                C6JK = -fC6 * (3 / 2) * EhJK * AlphaDat(atoms[j]) * AlphaDat(atoms[k]) 
                FD3 = (Rjk - (2 * fRd6b * RDjk)) / (fRd6b * RDjk)
                Fs3 = 1 / (fRd6b * RDjk)
                RjkInv = 1 / Rjk
                
                # Calculate C9
                HOMO_prod = EhomoDat(atoms[i]) * EhomoDat(atoms[j]) * EhomoDat(atoms[k]) 
                HOMO_sum = (EhomoDat(atoms[i]) + EhomoDat(atoms[j]) + EhomoDat(atoms[k]))
                Alpha_prod =  AlphaDat(atoms[i]) * AlphaDat(atoms[j]) * AlphaDat(atoms[k])
                Eh = HOMO_prod * HOMO_sum * Alpha_prod
                C9ijk = (3/2) * Eh / ((EhomoDat(atoms[i]) + EhomoDat(atoms[j])) * (EhomoDat(atoms[i]) + EhomoDat(atoms[k])) * (EhomoDat(atoms[j]) + EhomoDat(atoms[k])))
                
                # Calculate the cosine values using the law of cosine equation
                Cosij = ((Rik**2 + Rjk**2) - (Rij**2))/(2 * Rik * Rjk) 
                if (Cosij > 1):
                    Cosij = 1
                if (Cosij < -1):
                    Cosij = -1
                Cosik = ((Rij**2 + Rjk**2) - (Rik**2))/(2 * Rij * Rjk) 
                if (Cosik > 1):
                    Cosik = 1
                if (Cosik < -1):
                    Cosik = -1                
                Cosjk = ((Rij**2 + Rik**2) - (Rjk**2))/(2 * Rij * Rik)
                if (Cosjk > 1):
                    Cosjk = 1
                if (Cosjk < -1):
                    Cosjk = -1 
                
                # Calculate the damp three-body dispersion energy
                ## Damping three-body radii corrections
                RD3Bij = fRd6b * RDij
                RD3Bik = fRd6b * RDik
                RD3Bjk = fRd6b * RDjk
                
                ## Damping functions for PFD-3B (Part I)
                FD1 = (Rij - 2 * RD3Bij) / RD3Bij 
                FD2 = (Rik - 2 * RD3Bik) / RD3Bik 
                FD3 = (Rjk - 2 * RD3Bjk) / RD3Bjk 
                
                ## Damping functions for PFD-3B (Part II)
                FE1 = 1 - (1 + FD1) * np.exp(-2 * FD1)
                FE2 = 1 - (1 + FD2) * np.exp(-2 * FD2)
                FE3 = 1 - (1 + FD3) * np.exp(-2 * FD3)
                
                ## Three-body dispersion energy
                Ep = ((FE1 * FE2 * FE3)**2) * ((RijInv * RikInv * RjkInv)**3)
                T = ((3 * Cosij * Cosik * Cosjk) + 1)
                Disp_ijk = -C9ijk * T * Ep
                E_disp_3B += Disp_ijk
           
                
    return E_disp_3B

def H_bonding(atoms,atom_number, Distance_Matrix, Disp2B):
    """Compute the Petersson-Frisch hydrogen bonding corrections.
       
       Parameters
       ----------
       atoms: obj:'np.ndarray'
           the atomic elements
       atom_number: obj:'np.ndarray'
           the atomic number
       Distance_Matrix: obj:'np.ndarray'
           a matrix containing the atom-pair distance (in Ansgtroms)
       Disp2B: obj:'np.ndarray'
           two-body dispersion matrix 
       
       Returns
       -------
       E_Hbond: obj:'np.float64'
           hydrogen bonding energy (in Hartrees)
    """
    E_Hbond = 0
    size = atoms.shape[0]
    Distance_Matrix = Distance_Matrix / 0.52917721092 # convert distances from Angstroms to Bohr
    H_bondmat = np.zeros((size,size)) # H-bonding matrix contirbutions
    fHB, HBexp, HBangle = PFD3Bparam('h-bonding')
    
    # Convert HBangle from degrees to radians
    HBaRad = HBangle * np.pi / 2
    
    for i in range(size):
        for j in range(i):
            for k in range(j):
                # Get distances
                Rij = Distance_Matrix[i][j]
                Rik = Distance_Matrix[i][k]
                Rjk = Distance_Matrix[j][k]
                
                # Compute Cosines
                Cosij = ((Rik**2 + Rjk**2) - Rij**2)/(2 * Rik * Rjk) 
                if (Cosij > 1):
                    Cosij = 1
                if (Cosij < -1):
                    Cosij = -1
                Cosik = ((Rij**2 + Rjk**2) - Rik**2)/(2 * Rij * Rjk) 
                if (Cosik > 1):
                    Cosik = 1
                if (Cosik < -1):
                    Cosik = -1                
                Cosjk = ((Rij**2 + Rik**2) - Rjk**2)/(2 * Rij * Rik)
                if (Cosjk > 1):
                    Cosjk = 1
                if (Cosjk < -1):
                    Cosjk = -1
                    
                # Determine the angle 
                Angleij = np.arccos(Cosij)
                Angleik = np.arccos(Cosik)
                Anglejk = np.arccos(Cosjk)
                
                # Calculate hydrogen bonding corrections
                if (atom_number[i] == '1') and (HbndDat(atoms[j]) == True) and (HbndDat(atoms[k]) == True):
                    DispHBjk = fHB * (1 + np.tanh(HBexp * (Anglejk - HBaRad) / np.pi)) / 2
                    if (atom_number[j] == '34') or (atom_number[k] == '34'):
                        DispHBjk = DispHBjk / 2
                    H_bondmat[j][k] = Disp2B[j][k] * DispHBjk
                    E_Hbond += H_bondmat[j][k] 
                    
                if (atom_number[j] == '1') and (HbndDat(atoms[i]) == True) and (HbndDat(atoms[k]) == True):
                    DispHBik = fHB * (1 + np.tanh(HBexp * (Angleik - HBaRad) / np.pi)) / 2
                    if (atom_number[i] == '34') or (atom_number[k] == '34'):
                        DispHBik = DispHBik / 2
                    H_bondmat[i][k] = Disp2B[i][k] * DispHBik
                    E_Hbond += H_bondmat[i][k]
                    
                if (atom_number[k] == '1') and (HbndDat(atoms[i]) == True) and (HbndDat(atoms[j]) == True):
                    DispHBij = fHB * (1 + np.tanh(HBexp * (Angleij - HBaRad) / np.pi)) / 2
                    if (atom_number[i] == '34') or (atom_number[j] == '34'):
                        DispHBij = DispHBij / 2
                    H_bondmat[i][j] = Disp2B[i][j] * DispHBij
                    E_Hbond += H_bondmat[i][j]
    
    return E_Hbond