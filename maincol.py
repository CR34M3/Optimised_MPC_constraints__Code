#!/usr/bin/env python
"""
Main file for M Project (Optimised MPC Constraints).
Testing file for distillation column
Author: Andre Campher
"""
# Dependencies: - SciPy
#               - convertfuns
#               - auxfuns
#               - conclasses

from conclasses import ConSet
from scipy import array, linalg
from auxfuns import mat2ab
from fitting import fitset
from convertfuns import con2pscon

#MAIN START =================================================================

# define AIS and DOS (equations : Ax<b)
# equations in the form Ax<b, matrix = [A s b] 
# with s the sign vector [1:>, -1:<]
AISA, AISs, AISb = mat2ab(array([[1., 0., 1.,  11],
                                 [1., 0., -1., 15],
                                 [0., 1.,  1., 25],
                                 [0., 1., -1., 90.]]))
AIS = ConSet(AISA, AISs, AISb) 

DOSA, DOSs, DOSb = mat2ab(array([[1., 0.,  1., 66.],
                                 [1., 0., -1., 68.],
                                 [0., 1.,  1., 78.],
                                 [0., 1., -1., 82.]]))
DOS = ConSet(DOSA, DOSs, DOSb)

# define G (steady-state model)
#              R     T10sp
G = array([[-0.0575, 0.96],    # T1
           [ -0.146, 0.518]])  # T8
Gi = linalg.inv(G)

# calc AOS (from G and AIS)
lss = array([[68., 78.]]) # nominal operating point (used for model generation)
AOSA, AOSs, AOSb = AIS.outconlin(G, lss)
AOS = ConSet(AOSA, AOSs, AOSb)
#TODO: check nominal op.point use

# Calc intersection of AOS|DOS
DOSi = ConSet(*mat2ab(array([[-0.47486499,  0.88005866, -1, 36.55612415],
                             [ 1.,         0.,         -1,  68.],
                             [ 0.,        -1.,         -1, -78.]])))


# Calc additional spaces
DIS = ConSet(*DOS.outconlin(Gi, AIS.cscent))
DOSn = fitset(DOSi, 'r', 'a')
DISn = ConSet(*DOSn.outconlin(Gi, AIS.cscent))
DISi = ConSet(*DOSi.outconlin(Gi, AIS.cscent))
DIS2 = ConSet(*AIS.intersect(DIS))
DOS2 = ConSet(*DIS2.outconlin(G, lss))

# Calc modified constraints and model for high/low limits
modA, mods, modb, modG = con2pscon(DOSi, G, 'o')