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
from pylab import plot, show, legend, xlabel, ylabel
from fitting import fitset

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
DIS = ConSet(*DOS.outconlin(Gi, AIS.cscent))
DOSn = fitset(DOSi, 'r', 'a')
DISn = ConSet(*DOSn.outconlin(Gi, AIS.cscent))
DISi = ConSet(*DOSi.outconlin(Gi, AIS.cscent))
DIS2 = ConSet(*AIS.intersect(DIS))
DOS2 = ConSet(*DIS2.outconlin(G, lss))

# Plotting
# -- AIS
pl = array([1, 0, 2, 3, 1])
pl2 = array([0, 1, 2, 0])
#plot(AIS.vert[pl, 0], AIS.vert[pl, 1], 'r', linewidth=3)
#xlabel("R (mA)")
#ylabel("T$_{10sp}$ ($^\circ$C)")
#show()

# -- AOS / DOS
plot(AOS.vert[pl, 0], AOS.vert[pl, 1], 'b', linewidth=1)
plot(DOS.vert[pl, 0], DOS.vert[pl, 1], 'k--', linewidth=1)
#plot(DOS2.vert[pl, 0], DOS2.vert[pl, 1], 'k--', linewidth=3)
xlabel("T$_{1}$ ($^\circ$C)")
ylabel("T$_{8}$ ($^\circ$C)")
#show()

#DOSn = fitset(DOSi, 'r', 'a')
#
plot(DOSn.vert[pl, 0], DOSn.vert[pl, 1], 'k--', linewidth=3)
plot(DOSi.vert[pl2, 0], DOSi.vert[pl2, 1], 'g', linewidth=3)
show()

# -- Input space constriants
#plot(AIS.vert[pl, 0], AIS.vert[pl, 1], 'r', linewidth=3)
#plot(DIS.vert[pl, 0], DIS.vert[pl, 1], 'b--')
#plot(DISn.vert[pl, 0], DISn.vert[pl, 1], 'b--', linewidth=3)
#plot(DISi.vert[pl2, 0], DISi.vert[pl2, 1], 'k-.', linewidth=3)
#xlabel("R (mA)")
#ylabel("T$_{10sp}$ ($^\circ$C)")
#show()
