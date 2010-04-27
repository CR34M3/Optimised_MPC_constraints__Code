# main.py
# Main file for M Project (Optimised MPC Constraints)
#
# Author: Andre Campher
# 2010-04-27
#
# Dependencies: SciPy

from scipy import *

# --skeleton--
# insert AIS (equations : Ax<b)
#	- convert to vertices on feasible region
#   - con2vert [translation in progress]
AIS = mat('1 0 1 -1; 1 0 -1 1; 0 1 1 -1; 0 1 -1 1') #equations in the form Ax<b, matrix = [A s b] with s the sign vector
#AIS_s =
#AIS_A =
#AIS_b =  
AISvert = mat('-1 -1; -1 1; 1 -1; 1 1') #manual vertex insertion while con2vert is being translated

# insert DOS (equations : Ax<b)
#	- convert to vertices on feasible region (necessary?)
DOS = mat('1 0 1 -1; 1 0 -1 1; 0 1 1 -1; 0 1 -1 1') #equations in the form Ax<b, matrix = [A s b] with s the sign vector
#DOS_s =
#DOS_A =
#DOS_b =  
DOSvert = mat('-1 -1; -1 1; 1 -1; 1 1') #manual vertex insertion while con2vert is being translated
tester = mat('1; 1')

# insert G (steady-state model)
G = mat('1.2 0.5; 0.5 1.2') #gain matrix - atm linear steady state matrix

# calc AOS (from G and AIS)
#	- convert vertices to equations (Ax<b)
#   - vert2con [translation in progress]
AOSverttemp = empty([1,2])
for vert in AISvert:
	x = G*vert.transpose()
	AOSverttemp = vstack((AOSverttemp,x.transpose()))

AOSvert = AOSverttemp[1:,:] #remove first line of junk data from AOSverttemp
print AOSvert
# Calc intersection of AOS|DOS

# (Vectors defining inside/outside to calc)
