# main.py
# Main file for M Project (Optimised MPC Constraints)
#
# Author: Andre Campher
# 2010-04-27
#
# Dependancies: NumPy

from numpy import matrix

# --skeleton--
# insert AIS (equations : Ax<b)
#	- convert to vertices on feasible region
#   - con2vert [in progress]
AIS = matrix('1 0 1 0; 1 0 -1 1; 0 1 1 0; 0 1 -1 1') #equations in the form Ax<b, matrix = [A s b] with s the sign vector
#AIS_s =
#AIS_A =
#AIS_b =  
AISvert = matrix('0 0; 0 1; 1 0; 1 1')
print AISvert

# insert DOS (equations : Ax<b)
#	- convert to vertices on feasible region (necessary?)

# insert G (steady-state model)

# calc AOS (from G and AIS)
#	- convert vertices to equations (Ax<b)

# Calc intersection of AOS|DOS

# (Vectors defining inside/outside to calc)
