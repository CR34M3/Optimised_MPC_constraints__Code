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
from pylab import plot, show, legend, xlabel, ylabel
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
           [ 0.0111, -0.0604]])  # gain matrix - atm linear steady state matrix

Gi = linalg.inv(G)

# calc AOS (from G and AIS)
lss = array([[1.322, 16.4048]]) # nominal operating point (used for model generation)
AOS = ConSet(*AIS.outconlin(G, lss))

pl = array([1, 0, 2, 3, 1])
# -- AIS
#plot(AIS.vert[pl, 0], AIS.vert[pl, 1], 'r', linewidth=3)
#xlabel("x$_1$ (mA)")
#ylabel("x$_2$ (mA)")
#show()
mbox = ConSet(*fitmaxbox(AOS, 0.2).intersect(POS))
# -- AOS / DOS
plot(mbox.vert[pl, 0], mbox.vert[pl, 1], 'r.-.', linewidth=3)
plot(POS.vert[pl, 0], POS.vert[pl, 1], 'k', linewidth=3)
plot(AOS.vert[pl, 0], AOS.vert[pl, 1], 'b', linewidth=1)
#plot(DOS.vert[pl, 0], DOS.vert[pl, 1], 'k--', linewidth=1)
xlabel("F (gpm)")
ylabel("L (cm)")
show()

# -- AOS/DOS intersection
#DOSi = ConSet(*AOS.intersect(DOS))
#DIS = ConSet(*DOS.outconlin(Gi, AIS.cscent))

# -- Fitted constraints
#DOSn = fitset(DOSi, 'r', 'a')

#plot(DOSn.vert[pl, 0], DOSn.vert[pl, 1], 'g', linewidth=3)
#show()
#DISn = ConSet(*DOSn.outconlin(Gi, AIS.cscent))
#DISi = ConSet(*DIS.intersect(AIS))
#
#print DISn.vol()/DISi.vol()

# -- Input space constriants
#plot(AIS.vert[pl, 0], AIS.vert[pl, 1], 'r', linewidth=3)
#plot(DIS.vert[pl, 0], DIS.vert[pl, 1], 'b--')
#plot(DISn.vert[pl, 0], DISn.vert[pl, 1], 'b--', linewidth=3)
#xlabel("x$_1$ (mA)")
#ylabel("x$_2$ (mA)")
#show()