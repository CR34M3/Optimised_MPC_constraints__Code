#!/usr/bin/env python

# main.py
# Main file for M Project (Optimised MPC Constraints)
#
# Author: Andre Campher
# 2010-04-27
#
# Dependencies: SciPy

from scipy import *
from con2vert import con2vert
from vert2con import vert2con

def mat2ab(Asbmat): 	# transform to standard Ax<b notation with sign vector
	A = Asbmat[:,:-2]
	s = Asbmat[:,-2]
	b = Asbmat[:,-1]
	A = multiply(A,-s)
	b = multiply(b,-s)
	return A,s,b

#===== MAIN START =====

# insert AIS (equations : Ax<b)
AIS = mat('1 0  1 -1;\
	   1 0 -1  1;\
	   0 1  1 -1;\
	   0 1 -1  1') #equations in the form Ax<b, matrix = [A s b] with s the sign vector [1:>, -1:<]
# insert DOS (equations : Ax<b)
DOS = mat('1 0  1 -1;\
	   1 0 -1  1;\
	   0 1  1 -1;\
	   0 1 -1  1') #equations in the form Ax<b, matrix = [A s b] with s the sign vector [1:>, -1:<]
# insert G (steady-state model)
G = mat('1.2 0.5; 0.5 1.2') #gain matrix - atm linear steady state matrix


# convert to A,b form
AIS_A,AIS_s,AIS_b = mat2ab(AIS)
DOS_A,DOS_s,DOS_b = mat2ab(DOS)

# convert to vertices on feasible region
AISvert = con2vert(AIS_A,AIS_b)
DOSvert = con2vert(DOS_A,DOS_b)

# calc AOS (from G and AIS)
AOSverttemp = empty([1,AISvert.shape[1]])
for vert in AISvert:
	x = G*vert.transpose()
	AOSverttemp = vstack((AOSverttemp,x.transpose()))
AOSvert = AOSverttemp[1:,:] #remove first line of junk data from AOSverttemp

# convert AOS vertices to equations (Ax<b)
AOS_A,AOS_s,AOS_b = vert2con(AOSvert)

# Calc intersection of AOS|DOS
#	- using con2vert with both AOS constraints and DOS constraints should result in the intersection
combA = vstack((AOS_A,DOS_A))
combb = vstack((AOS_b,DOS_b))
intAOSDOSvert = con2vert(combA,combb)
intAOSDOScon = vert2con(intAOSDOSvert)

print intAOSDOScon

#TODO
# (Vectors defining inside/outside to calc)
# (Plotting for lower dimensional shapes)
