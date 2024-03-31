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

"""Return the output file containing two-body, three-body, and hydrogen bonding corrections for the Petersson-Frisch dispersion model."""
import numpy as np
from read_output import read_entire_output
from PFDisp import PFD_2B, PFD_3B, H_bonding
from parameters import EhomoDat, AlphaDat

class PFD:
    def __init__(self,name):
        data = read_entire_output(name)
        Atoms,AtomN, XYZ, DistanceMatrix, APF, Disp2B = data
        self.Atoms = Atoms
        self.AtomN = AtomN
        self.XYZ = XYZ
        self.DistanceMatrix = DistanceMatrix
        self.APF = APF
        self.Disp2B = Disp2B
        
        return None
        
    def output(self):
        """Retrive the output of the result"""
        
        
        print(f' NAtoms={self.Atoms.shape[0]:>13}')
        print(f'   E(RAPFD) = {self.APF:>4}')
        print(f'Atom Atomic Number    x           y           z          Ehomo       Alpha')
        for l in range(self.Atoms.shape[0]):
            print(f'{l+1:>3} {self.AtomN[l]:>9} {self.XYZ[l][0]:>14.6f} {self.XYZ[l][1]:>11.6f} {self.XYZ[l][2]:>11.6f} {EhomoDat(self.Atoms[l]):>11.8f}  {AlphaDat(self.Atoms[l]):>11.8f}')
        # Two-body corrections
        two_body, Disp2B = PFD_2B(self.Atoms, self.AtomN, self.DistanceMatrix)
        print(f'PFD-2B dispersion energy =  {np.round(two_body * 10**3,10):0.10f} milliHartrees.')

        # Three-body corrections
        three_body = PFD_3B(self.Atoms, self.AtomN, self.DistanceMatrix)
        print(f'PFD-3B dispersion energy =   {np.round(three_body * 10**3,10):0.10f} milliHartrees.')

        # Hydrogen bonding corrections
        E_hbond = H_bonding(self.Atoms, self.AtomN, self.DistanceMatrix, Disp2B)
        print(f'    Hydrogen bond energy =   {np.round(E_hbond * 10**3,10):0.10f} milliHartrees.')

        # Two-body + Three-body + H-bonding
        print(f'Total correction to Epf  =  {np.round((E_hbond + three_body + two_body)* 10**3,10):0.10f} milliHartrees.')
        
        # Create inheritance over two-body, three-body, and H-bonding
        self.two_body = np.round(two_body * 10**3,10) # in milliHartrees
        self.three_body = np.round(three_body * 10**3,10) # in milliHartrees
        self.hbond = np.round(E_hbond * 10**3,10) # in milliHartrees
        self.total_Epf = np.round((two_body + three_body + E_hbond)*10**3,10) # in milliHartrees
        
        return None
        