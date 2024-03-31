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

"""Parameters for the PFD-3B dispersion model. 
   These include hydrogen bonding, polarizability, and HOMO energy"""

import numpy as np

def EhomoDat(atom):
    """HOMO energies usef in the empirical dispersion energy according 
       to the APF formula.
       
       Parameters
       ----------
       atom: obj:'str'
            A single atom
            
       Returns
       -------
       Ehomo: obj:'np.float64'
           The HOMO energy for the atomic species (in Hartrees)
       """
                                         # Source:
    Ehomo = {'TV': np.float64(0),        # dummy variables
             'X': np.float64(0),         # dummy variables
             'Bq': np.float64(0),        # dummy variables
             'H': np.float64(-0.59277),  # H2
             'He':np.float64(-0.91786),
             'Li': np.float64(-2.48668), # Li(1s)
             'Be':np.float64(-0.83643),  # Be
             'B':np.float64(-0.49781),   # B
             'C':np.float64(-0.45579),   # C2
             'N': np.float64(-0.57073),  # N
             'O':np.float64(-0.67794),   # O
             'F':np.float64(-0.43097),   # F-
             'Ne':np.float64(-0.85114),  # Ne
             'Na':np.float64(-1.51908),  # Na(2p)
             'Mg':np.float64(-0.31883),  # Mg
             'Al':np.float64(-0.2179),   # Al
             'Si':np.float64(-0.27664),  # Si(Si_2H_4)
             'P':np.float64(-0.37275),   # P
             'S':np.float64(-0.41589),   # S
             'Cl':np.float64(-0.25031),  # Cl-
             'Ar':np.float64(-0.59132),  # Ar
             'K':np.float64(-0.72581),   # K(3p,4s)
             'Ca':np.float64(-0.36045),  # Ca
             'Sc':np.float64(-0.48732),  # Sc+
             'Ti':np.float64(-0.47558),  # Ti+
             'V':np.float64(-0.49997),   # V+
             'Cr':np.float64(-0.52235),  # Cr+
             'Mn':np.float64(-0.86849),  # Mn+
             'Fe':np.float64(-0.80746),  # Fe+
             'Co':np.float64(-0.78606),  # Co
             'Ni':np.float64(-0.78761),  # Ni+
             'Cu':np.float64(-0.81001),  # Cu(3d)
             'Zn':np.float64(-0.29256),  # Zn
             'Ga':np.float64(-0.62872),  # Ga(GaF3)
             'Ge':np.float64(-0.29600),  # Ge(Ge2H4)
             'As':np.float64(-0.37030),  # As
             'Se':np.float64(-0.41904),  # Se(SeH2)
             'Br':np.float64(-0.15934),  # Br-
             'Kr':np.float64(-0.52431),  # Kr
             'Rb':np.float64(-0.81007),  # Rb(4p)
             'Sr':np.float64(-0.25846),  # Sr
             'Y':np.float64(-0.43555),   # Y+
             'Zr':np.float64(-0.53274),  # Zr+
             'Nb':np.float64(-0.60845),  # Nb+
             'Mo':np.float64(-0.61026),  # Mo+
             'Tc':np.float64(-0.75892),  # Tc+
             'Ru':np.float64(-0.89013),  # Ru+
             'Rh':np.float64(-0.84578),  # Rh+
             'Pd':np.float64(-0.86934),  # Pd+
             'Ag':np.float64(-1.38234),  # Ag(4d)
             'Cd':np.float64(-0.26485),  # Cd 
             'In':np.float64(-0.19728),  # In
             'Sn':np.float64(-0.26504),  # Sn
             'Sb':np.float64(-0.29447),  # Sb
             'Te':np.float64(-0.31657),  # Te
             'I':np.float64(-0.36359),   # I
             'Xe':np.float64(-0.45784),  # Xe
             'Cs':np.float64(-0.86183),  # Cs
             'Ba':np.float64(-1.33004),  # Ba
             'La':np.float64(-0.68288),  # La
             'Ce':np.float64(-0.62554),  # Ce
             'Pr':np.float64(-0.63543),  # Pr
             'Nd':np.float64(-0.87969),  # Nd
             'Pm':np.float64(-0.90269),  # Pm
             'Sm':np.float64(-0.95965),  # Sm
             'Eu':np.float64(-1.01869),  # Eu
             'Gd':np.float64(-0.68117),  # Gd
             'Tb':np.float64(-0.79650),  # Tb
             'Dy':np.float64(-0.90495),  # Dy
             'Ho':np.float64(-0.83479),  # Ho
             'Er':np.float64(-0.90156),  # Er
             'Tm':np.float64(-0.97413),  # Tm
             'Yb':np.float64(-1.03326),  # Yb
             'Lu':np.float64(-0.71549),  # Lu
             'Hf':np.float64(-0.78396),  # Hf
             'Ta':np.float64(-0.82838),  # Ta
             'W':np.float64(-0.56104),   # W
             'Re':np.float64(-1.03841),  # Re
             'Os':np.float64(-0.86250),  # Os
             'Ir':np.float64(-0.95930),  # Ir
             'Pt':np.float64(-1.01755),  # Pt
             'Au':np.float64(-0.76490),  # Au
             'Hg':np.float64(-1.29467),  # Hg
             'Tl':np.float64(-0.19274),  # Tl
             'Pb':np.float64(-0.25601),  # Pb
             'Bi':np.float64(-0.31982),  # Bi
             'Po':np.float64(-0.19466),  # Po
             'At':np.float64(-0.26009),  # At
             'Rn':np.float64(-0.42799),  # Rn
             'Fr':np.float64(0),         # Fr
             'Ra':np.float64(0),         # Ra
             'Ac':np.float64(0),         # Ac
             'Th':np.float64(0),         # Th
             'Pa':np.float64(0),         # Pa
             'U':np.float64(0),          # U
             'Np':np.float64(0),         # Np
             'Pu':np.float64(0),         # Pu
             'Am':np.float64(0),         # Am
             'Cm':np.float64(0),         # Cm
             'Bk':np.float64(0),         # Bk
             'Cf':np.float64(0),         # Cf
             'Es':np.float64(0),         # Es
             'Fm':np.float64(0),         # Fm
             'Md':np.float64(0),         # Md
             'No':np.float64(0),         # No
             'Lr':np.float64(0),         # Lr
             'Rf':np.float64(0),         # Rf
             'Db':np.float64(0),         # Db
             'Sg':np.float64(0),         # Sg
             'Bh':np.float64(0),         # Bh
             'Hs':np.float64(0),         # Hs
             'Mt':np.float64(0),         # Mt
             'Ds':np.float64(0),         # Ds
             'Uu':np.float64(0),         # Uu
             'Ub':np.float64(0),         # Ub
             'Ut':np.float64(0),         # Ut
             'Uq':np.float64(0),         # Uq
             'Up':np.float64(0),         # Up
             'PF0->F':np.float64(-1.29000), # PF0 -> F
             'PF0->Cl':np.float64(-1.21000), # PF0 -> Cl
             'PF0->Br':np.float64(-1.21000) # PF0 -> Br
            } 
    
    return Ehomo[atom]

def AlphaDat(atom):
    """Isotropic polarizabilities used in the empirical dispersion energy according to the PF0 formula.
       
       Parameters
       ----------
       atom: obj:'str'
            A single atom
            
       Returns
       -------
       alpha: obj:'np.float64'
           The isotropic polarizability for the atomic species (in A^3)
       """
                                         # Source:
    alpha = {'TV': np.float64(0),        # dummy variables
             'X': np.float64(0),         # dummy variables
             'Bq': np.float64(0),        # dummy variables
             'H': np.float64(1.9044),    # H
             'He':np.float64(1.3643),    # He
             'Li': np.float64(0.1204),   # Li(1s)
             'Be':np.float64(6.1686),    # Be
             'B':np.float64(7.4437),     # B
             'C':np.float64(8.7820),     # C
             'N': np.float64(5.6129),    # N
             'O':np.float64(4.8334),     # O
             'F':np.float64(10.9428),    # F
             'Ne':np.float64(2.7109),    # Ne
             'Na':np.float64(0.4334),    # Na2.0113
             'Mg':np.float64(15.0282),   # Mg
             'Al':np.float64(14.6788),   # Al
             'Si':np.float64(26.8194),   # Si(Si_2H_4)
             'P':np.float64(23.8290),    # P
             'S':np.float64(19.8022),    # S
             'Cl':np.float64(30.0468),   # Cl
             'Ar':np.float64(10.9745),   # Ar
             'K':np.float64(5.0994),     # K
             'Ca':np.float64(15.1949),   # Ca
             'Sc':np.float64(15.7342),   # Sc
             'Ti':np.float64(23.2528),   # Ti
             'V':np.float64(16.1366),    # V
             'Cr':np.float64(20.6438),   # Cr
             'Mn':np.float64(24.5676),   # Mn
             'Fe':np.float64(15.0295),   # Fe
             'Co':np.float64(17.9818),   # Co
             'Ni':np.float64(16.3963),   # Ni
             'Cu':np.float64(18.2705),   # Cu
             'Zn':np.float64(19.2898),   # Zn
             'Ga':np.float64(23.4477),   # Ga
             'Ge':np.float64(42.0000),   # Ge(Ge2H4)
             'As':np.float64(28.3542),   # As
             'Se':np.float64(31.2402),   # Se(SeH2)
             'Br':np.float64(84.5580),   # Br
             'Kr':np.float64(16.9594),   # Kr
             'Rb':np.float64(15.6479),   # Rb
             'Sr':np.float64(18.6897),   # Sr
             'Y':np.float64(19.3530),    # Y
             'Zr':np.float64(34.3444),   # Zr
             'Nb':np.float64(19.8480),   # Nb
             'Mo':np.float64(27.2140),   # Mo
             'Tc':np.float64(30.2182),   # Tc
             'Ru':np.float64(24.1099),   # Ru
             'Rh':np.float64(22.1176),   # Rh
             'Pd':np.float64(21.2047),   # Pd
             'Ag':np.float64(7.5669),    # Ag
             'Cd':np.float64(23.7264),   # Cd 
             'In':np.float64(29.3879),   # In
             'Sn':np.float64(30.1278),   # Sn
             'Sb':np.float64(41.9382),   # Sb
             'Te':np.float64(38.5203),   # Te
             'I':np.float64(33.9511),    # I
             'Xe':np.float64(27.5774),   # Xe
             'Cs':np.float64(15.646),    # Cs
             'Ba':np.float64(10.446),    # Ba
             'La':np.float64(20.183),    # La
             'Ce':np.float64(49.954),    # Ce
             'Pr':np.float64(47.643),    # Pr
             'Nd':np.float64(11.869),    # Nd
             'Pm':np.float64(10.867),    # Pm
             'Sm':np.float64(9.961),     # Sm
             'Eu':np.float64(9.187),     # Eu
             'Gd':np.float64(36.798),    # Gd
             'Tb':np.float64(13.407),    # Tb
             'Dy':np.float64(9.382),     # Dy
             'Ho':np.float64(9.817),     # Ho
             'Er':np.float64(8.908),     # Er
             'Tm':np.float64(8.407),     # Tm
             'Yb':np.float64(8.058),     # Yb
             'Lu':np.float64(12.960),    # Lu
             'Hf':np.float64(26.792),    # Hf
             'Ta':np.float64(24.930),    # Ta
             'W':np.float64(34.563),     # W
             'Re':np.float64(11.186),    # Re
             'Os':np.float64(10.944),    # Os
             'Ir':np.float64(10.127),    # Ir
             'Pt':np.float64(9.404),     # Pt
             'Au':np.float64(13.666),    # Au
             'Hg':np.float64(7.810),     # Hg
             'Tl':np.float64(69.360),    # Tl
             'Pb':np.float64(58.663),    # Pb
             'Bi':np.float64(47.969),    # Bi
             'Po':np.float64(45.274),    # Po
             'At':np.float64(39.229),    # At
             'Rn':np.float64(33.066),    # Rn
             'Fr':np.float64(0),         # Fr
             'Ra':np.float64(0),         # Ra
             'Ac':np.float64(0),         # Ac
             'Th':np.float64(0),         # Th
             'Pa':np.float64(0),         # Pa
             'U':np.float64(0),          # U
             'Np':np.float64(0),         # Np
             'Pu':np.float64(0),         # Pu
             'Am':np.float64(0),         # Am
             'Cm':np.float64(0),         # Cm
             'Bk':np.float64(0),         # Bk
             'Cf':np.float64(0),         # Cf
             'Es':np.float64(0),         # Es
             'Fm':np.float64(0),         # Fm
             'Md':np.float64(0),         # Md
             'No':np.float64(0),         # No
             'Lr':np.float64(0),         # Lr
             'Rf':np.float64(0),         # Rf
             'Db':np.float64(0),         # Db
             'Sg':np.float64(0),         # Sg
             'Bh':np.float64(0),         # Bh
             'Hs':np.float64(0),         # Hs
             'Mt':np.float64(0),         # Mt
             'Ds':np.float64(0),         # Ds
             'Uu':np.float64(0),         # Uu
             'Ub':np.float64(0),         # Ub
             'Ut':np.float64(0),         # Ut
             'Uq':np.float64(0),         # Uq
             'Up':np.float64(0),         # Up
             'PF0->F':np.float64(-3.2000), # PF0 -> F
             'PF0->Cl':np.float64(-2.1000), # PF0 -> Cl
             'PF0->Br':np.float64(-2.3000) # PF0 -> Br
            } 
    
    return alpha[atom]

def HbndDat(atom):
    """Table of electronegative atoms capable of forming hydrogen bonds. 
       
       Parameters
       ----------
       atom: obj:'str'
            A single atom
            
       Returns
       -------
       hbond: obj:'bool'
           true or false if the atoms are electronegative to form a hydrogen bond
       """
                          # Source:
    hbond = {'TV':False,  # dummy variables
             'X':False,   # dummy variables
             'Bq':False,  # dummy variables
             'H': False,  # H
             'He':False,  # He
             'Li':False,  # Li(1s)
             'Be':False,  # Be
             'B':False,   # B
             'C':False,   # C
             'N':True,    # N
             'O':True,    # O
             'F':True,    # F
             'Ne':False,  # Ne
             'Na':False,  # Na
             'Mg':False,  # Mg
             'Al':False,  # Al
             'Si':False,  # Si(Si_2H_4)
             'P':True,    # P
             'S':True,    # S
             'Cl':True,   # Cl
             'Ar':False,  # Ar
             'K':False,   # K
             'Ca':False,  # Ca
             'Sc':False,  # Sc
             'Ti':False,  # Ti
             'V':False,   # V
             'Cr':False,  # Cr
             'Mn':False,  # Mn
             'Fe':False,  # Fe
             'Co':False,  # Co
             'Ni':False,  # Ni
             'Cu':False,  # Cu
             'Zn':False,  # Zn
             'Ga':False,  # Ga
             'Ge':False,  # Ge(Ge2H4)
             'As':True,   # As
             'Se':True,   # Se(SeH2)
             'Br':True,   # Br
             'Kr':False,  # Kr
             'Rb':False,  # Rb
             'Sr':False,  # Sr
             'Y':False,   # Y
             'Zr':False,  # Zr
             'Nb':False,  # Nb
             'Mo':False,  # Mo
             'Tc':False,  # Tc
             'Ru':False,  # Ru
             'Rh':False,  # Rh
             'Pd':False,  # Pd
             'Ag':False,  # Ag
             'Cd':False,  # Cd 
             'In':False,  # In
             'Sn':False,  # Sn
             'Sb':False,  # Sb
             'Te':False,  # Te
             'I':False,   # I
             'Xe':False,  # Xe
             'Cs':False,  # Cs
             'Ba':False,  # Ba
             'La':False,  # La
             'Ce':False,  # Ce
             'Pr':False,  # Pr
             'Nd':False,  # Nd
             'Pm':False,  # Pm
             'Sm':False,  # Sm
             'Eu':False,  # Eu
             'Gd':False,  # Gd
             'Tb':False,  # Tb
             'Dy':False,  # Dy
             'Ho':False,  # Ho
             'Er':False,  # Er
             'Tm':False,  # Tm
             'Yb':False,  # Yb
             'Lu':False,  # Lu
             'Hf':False,  # Hf
             'Ta':False,  # Ta
             'W':False,   # W
             'Re':False,  # Re
             'Os':False,  # Os
             'Ir':False,  # Ir
             'Pt':False,  # Pt
             'Au':False,  # Au
             'Hg':False,  # Hg
             'Tl':False,  # Tl
             'Pb':False,  # Pb
             'Bi':False,  # Bi
             'Po':False,  # Po
             'At':False,  # At
             'Rn':False,  # Rn
             'Fr':False,  # Fr
             'Ra':False,  # Ra
             'Ac':False,  # Ac
             'Th':False,  # Th
             'Pa':False,  # Pa
             'U':False,   # U
             'Np':False,  # Np
             'Pu':False,  # Pu
             'Am':False,  # Am
             'Cm':False,  # Cm
             'Bk':False,  # Bk
             'Cf':False,  # Cf
             'Es':False,  # Es
             'Fm':False,  # Fm
             'Md':False,  # Md
             'No':False,  # No
             'Lr':False,  # Lr
             'Rf':False,  # Rf
             'Db':False,  # Db
             'Sg':False,  # Sg
             'Bh':False,  # Bh
             'Hs':False,  # Hs
             'Mt':False,  # Mt
             'Ds':False,  # Ds
             'Uu':False,  # Uu
             'Ub':False,  # Ub
             'Ut':False,  # Ut
             'Uq':False,  # Uq
             'Up':False,  # Up
             'Uh':False,  # Uh
             'Us':False,  # Us
             'Uo':False   # Uo
            } 
    
    return hbond[atom]

def PFD3Bparam(whatDisp):
    """Return parameters two-body and three-body dispersion parameters.
    
       Parameters
       ----------
       whatDisp: obj:'str'
            name of the correction (e.g., 'two-body', 'three-body', and 'h-bonding')
            
       Returns
       -------
       data: obj:'tuple'
          parameters for the dispersion models
    """
    
    fR0 =  -0.443820000
    fRd =   2.829285000
    fRs = -13.141740000
    fC6 =   1.186044000
    fHe =   0.872259905
    fZH =   1.259456000
    fRd6b = 0.500000000
    fHB =  -0.866025403784438
    HBexp = 50.00000000
    HBangle = 1.14285714285714
    
    if whatDisp == 'two-body':
        return (fR0, fRd, fRs, fC6, fHe, fZH, fRd6b)
    elif whatDisp == 'three-body':
        return (fR0, fRd, fRs, fC6,fHe, fZH, fRd6b)
    elif whatDisp == 'h-bonding':
        return (fHB, HBexp, HBangle)
    else:
        return (fR0, fRd, fRs, fC6,fHe, fZH, fRd6b, fHB, HBexp, HBangle)