#!/usr/bin/env python
"""
Main file for M Project (Optimised MPC Constraints).
Author: Andre Campher
"""
# Dependencies: - SciPy
#               - convertfuns
#               - auxfuns
#               - conclasses

from conclasses import conset
from scipy import mat

#MAIN START ===================================================================================================

# define AIS and DOS (equations : Ax<b)
# equations in the form Ax<b, matrix = [A s b] with s the sign vector [1:>, -1:<]
AIS = conset(mat('1 0  1 -1;\
                  1 0 -1  1;\
                  0 1  1 -1;\
                  0 1 -1  1')) 

DOS = conset(mat('1 0  1 -1;\
                  1 0 -1  1;\
                  0 1  1 -1;\
                  0 1 -1  1'))

# define G (steady-state model)
G = mat('1.2 0.5; 0.5 1.2') #gain matrix - atm linear steady state matrix

# calc AOS (from G and AIS)
AOS = conset(AIS.outconlin(G))

# Calc intersection of AOS|DOS
print AOS.intersect(DOS)

#TODO
# (Vectors defining inside/outside to calc)
# (Plotting for lower dimensional shapes)