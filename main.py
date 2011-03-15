#!/usr/bin/env python
"""
Main file for M Project (Optimised MPC Constraints).
Author: Andre Campher
"""
# Dependencies: - SciPy
#               - convertfuns
#               - auxfuns
#               - conclasses

from conclasses import ConSet
from scipy import array
from auxfuns import mat2ab

#MAIN START =================================================================

# define AIS and DOS (equations : Ax<b)
# equations in the form Ax<b, matrix = [A s b] 
# with s the sign vector [1:>, -1:<]
AISA, AISs, AISb = mat2ab(array([[1., 0., 1.,  -0.0525],
                                 [1., 0., -1., 0.125],
                                 [0., 1.,  1., -10],
                                 [0., 1., -1., 10.]]))
AIS = ConSet(AISA, AISs, AISb) 

DOSA, DOSs, DOSb = mat2ab(array([[1., 0.,  1., -1.],
                                 [1., 0., -1.,  1.],
                                 [0., 1.,  1., -1.],
                                 [0., 1., -1.,  1.]]))
DOS = ConSet(DOSA, DOSs, DOSb)

# define G (steady-state model)
lss = array([[50., 50.]]) # nominal operating point (used for model generation)
G = array([[1, 0.0025],
           [2, 0.0025]])  # gain matrix - atm linear steady state matrix

# calc AOS (from G and AIS)
AOSA, AOSs, AOSb = AIS.outconlin(G, AIS.cscent, lss)  
AOS = ConSet(AOSA, AOSs, AOSb)
print AOS.vert

# Calc intersection of AOS|DOS
#print AOS.intersect(DOS)