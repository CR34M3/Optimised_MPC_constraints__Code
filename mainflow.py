#!/usr/bin/env python
"""
Main file for M Project (Optimised MPC Constraints).
Testing file for level and flow rig
Author: Andre Campher
"""
# Dependencies: - SciPy
#               - convertfuns
#               - auxfuns
#               - conclasses

from conclasses import ConSet
from scipy import array, linalg
from auxfuns import mat2ab
from fitting import fitset, fitmaxbox

#MAIN START =================================================================

# define AIS and DOS (equations : Ax<b)
# equations in the form Ax<b, matrix = [A s b] 
# with s the sign vector [1:>, -1:<]
AISA, AISs, AISb = mat2ab(array([[1., 0., 1.,  4],
                                 [1., 0., -1., 20],
                                 [0., 1.,  1., 4],
                                 [0., 1., -1., 20.]]))
AIS = ConSet(AISA, AISs, AISb) 

DOSA, DOSs, DOSb = mat2ab(array([[1., 0.,  1., 1.2],
                                 [1., 0., -1., 1.5],
                                 [0., 1.,  1., 15.],
                                 [0., 1., -1., 18.]]))
DOS = ConSet(DOSA, DOSs, DOSb)

POS = ConSet(*mat2ab(array([[1., 0.,  1., 0],
                            [1., 0., -1., 1.9],
                            [0., 1.,  1., 6.],
                            [0., 1., -1., 25.]])))

# define G (steady-state model)
G = array([[-0.0476, -0.0498],
           [ 0.0111, -0.0604]])

Gi = linalg.inv(G)

# calc AOS (from G and AIS)
lss = array([[1.322, 16.4048]]) # nominal operating point (used for model generation)
AOS = ConSet(*AIS.outconlin(G, AIS.cscent, lss))
mbox = ConSet(fitmaxbox(AOS, 0.2).intersect(POS))

# -- AOS/DOS intersection
DOSi = ConSet(AOS.intersect(DOS))
DIS = ConSet(*DOS.outconlin(Gi, lss, AIS.cscent))

# -- Fitted constraints
DOSn = fitset(DOSi, 'r', 'a')
DISn = ConSet(*DOSn.outconlin(Gi, lss, AIS.cscent))
DISi = ConSet(DIS.intersect(AIS))

print DISi