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

"""Read the Gaussian output file"""

import numpy as np

def periodic_table(atom_number):
    """Convert the atom numbers into atom symbols
    
       Parameter 
       ----------
       atom_number: obj: 'np.int'
           the atom number 
        
       Return
       -------
       nAtoms: obj: 'np.ndarray'
           the atom symbol
       """
    
                                    
    Atoms = {'1':'H',  
             '2':'He',
             '3':'Li', 
             '4':'Be',  
             '5':'B',  
             '6':'C',   
             '7':'N', 
             '8':'O',  
             '9':'F',   
             '10':'Ne', 
             '11':'Na', 
             '12':'Mg',
             '13':'Al', 
             '14':'Si', 
             '15':'P', 
             '16':'S',
             '17':'Cl',
             '18':'Ar',
             '19':'K',
             '20':'Ca',
             '21':'Sc',
             '22':'Ti', 
             '23':'V',
             '24':'Cr',
             '25':'Mn',
             '26':'Fe',
             '27':'Co',
             '28':'Ni',
             '29':'Cu',
             '30':'Zn',
             '31':'Ga',
             '32':'Ge',
             '33':'As',
             '34':'Se',
             '35':'Br',
             '36':'Kr',
             '37':'Rb',
             '38':'Sr',
             '39':'Y',
             '40':'Zr',
             '41':'Nb',
             '42':'Mo',
             '43':'Tc',
             '44':'Ru',
             '45':'Rh',
             '46':'Pd',
             '47':'Ag',
             '48':'Cd',
             '49':'In',
             '50':'Sn',
             '51':'Sb',
             '52':'Te',
             '53':'I',
             '54':'Xe',
             '55':'Cs',
             '56':'Ba',
             '57':'La',
             '58':'Ce',
             '59':'Pr',
             '60':'Nd',
             '61':'Pm',
             '62':'Sm',
             '63':'Eu',
             '64':'Gd',
             '65':'Tb',
             '66':'Dy',
             '67':'Ho',
             '68':'Er',
             '69':'Tm',
             '70':'Yb',
             '71':'Lu',
             '72':'Hf',
             '73':'Ta',
             '74':'W',
             '75':'Re',
             '76':'Os',
             '77':'Ir',
             '78':'Pt',
             '79':'Au',
             '80':'Hg',
             '81':'Tl',
             '82':'Pb',
             '83':'Bi',
             '84':'Po',
             '85':'At',
             '86':'Rn',
             '87':'Fr',
             '88':'Ra',
             '89':'Ac',
             '90':'Th',
             '91':'Pa',
             '92':'U',
             '93':'Np',
             '94':'Pu',
             '95':'Am',
             '96':'Cm',
             '97':'Bk',
             '98':'Cf',
             '99':'Es',
             '100':'Fm',
             '101':'Md',
             '102':'No',
             '103':'Lr',
             '104':'Rf',
             '105':'Db',
             '106':'Sg',
             '107':'Bh',
             '108':'Hs',
             '109':'Mt',
             '110':'Ds',
             '111':'Uu',
             '112':'Ub',
             '113':'Ut',
             '114':'Uq', 
             '115':'Up'
            }
                
    return Atoms[atom_number]

def read_coordinates(vertical_position):
    """Read the 'Standard Orientation' xyz geometries
    
       Parameter
       ---------
       vertical_position: obj:'str' or 'list' of 'str' 
           The list of strings present in the ROHF file 
           
       Return
       ------
       XYZ: obj:'np.ndarray'
           The standard orientation xyz geometries
       AtomN: obj:'np.ndarray'
           Array of atom numbers as ordered by the atom array 
       """
    
    # Read the XYZ file
    
    n = 0
    for l in vertical_position:
        if l.startswith(' --------'):
            break          
        else:       
            X, Y, Z = l[34::].split()
            atom_number = l[9:21].replace(' ','')
            if n == 0:
                XYZ = np.array([np.array([np.float64(X),np.float64(Y),np.float64(Z)])])
                AtomN = np.array([atom_number])
                n = 1
            else:
                XYZ = np.concatenate((XYZ,np.array([np.array([np.float64(X),np.float64(Y),np.float64(Z)])])))
                AtomN = np.concatenate((AtomN,np.array([atom_number])))
    return (XYZ, AtomN)

def distances(XYZ):
    """Compute the atom-pair distances
    
       Parameter 
       ----------
       XYZ: obj: 'np.int'
           the xyz coordinates 
        
       Return
       -------
       DistanceMatrix: obj: 'np.ndarray'
           a lower diagional matrix containing atom-pair distances 
    """
    
    length = XYZ.shape[0] # row size of the distances 
    DistanceMatrix = np.zeros((length,length), dtype=np.float64) # set distance matrix
    for param1 in range(1,length):
        for param2 in range(param1):
            d = np.sqrt((XYZ[param2][0] - XYZ[param1][0])**2 + (XYZ[param2][1] - XYZ[param1][1])**2 + (XYZ[param2][2] - XYZ[param1][2])**2)
            DistanceMatrix[param1][param2] = d
    return DistanceMatrix
    
def read_entire_output(FILE):
    """Read the entire APF output file. 
       
       Parameters
       ----------
       FILE: obj:'str'
           Name of the file/path
       
       Returns
       -------
       nAtoms: obj:'np.ndarray'
           Atoms present in the molecule
       AtomN : obj:'np.ndarray'
           Atomic number of the atoms
       XYZ: obj:'np.ndarray'
           Standard orientation xyz geometries
       DistanceMatrix: obj:'np.ndarray'
           Atom pair distances
       APF: obj:'np.float64'
           Energy from APF energy
    """
    
    o = 0
    with open(FILE, mode='r') as reader:
        p = reader.readlines()[20:]
        for n in p:
            o += 1 
            if n.startswith('                          Input orientation:'): # get coordinates 
                XYZ, AtomN = read_coordinates(p[o+4::])
                DistanceMatrix = distances(XYZ)
            elif n.startswith(' SCF Done:'): # get APF energy (in Hartrees)
                APF = np.float64(n[22:42].replace(' ',''))
            elif n.startswith(' R6Disp:  Petersson-Frisch Dispersion energy'): # get PF two-body dispersion energy 
                line = n[49:].replace('Hartrees.\n','')
                Disp2B = np.float64(line.replace(' ','')) # get PF dispersion energy(in Hartrees)
    reader.close()
    Atoms = np.array([])
    for l in range(AtomN.shape[0]):
        Atoms = np.concatenate((Atoms,np.array([periodic_table(AtomN[l])])))
    return (Atoms,AtomN, XYZ, DistanceMatrix, APF, Disp2B)