#!/usr/bin/env python

# vert2con.py
# Converts sets of vertices to a list of constraints (of the feasible region)
# Output, A,b and s of the set;  Ax < b   (s is the sign-vector [to be used later])
#
# Author: Andre Campher
#
# Dependencies: * qhull (libqhull5, qhull-bin)
#               * scipy

from scipy import *
import string
import subprocess #to use qhull

V = mat('0 0 0; 0 0 2; 0 2 0; 0 2 2; 2 0 0; 2 0 2; 2 1 2; 1 2 2; 2 2 0; 2 2 1') #vertex points matrix to test

# ================ generate data file for qhull
#format:
# n //dimensions (columns of V)
# p //number of points (rows of V)
# v //points e.g. 1 3 12 \n 2 4 5

infile = open("qhullin",'w') #generate input file
infile.write(str(V.shape[1])+"\n"+str(V.shape[0])+"\n")
for row in range(0,V.shape[0]):
	stringmat = ""
	for col in range(0,V.shape[1]):
		stringmat = string.join([stringmat,str(V[row,col])]," ")
	infile.write(stringmat+"\n")
infile.close()

# subprocess.call(["function","arguments"]) or subprocess.Popen('function expression', shell=True)
# run qhull with: qhull < data or cat data | qhull

qhullp = subprocess.Popen('qhull n < qhullin', shell=True, stdout=subprocess.PIPE) #calc convex hull and get normals
Vc = qhullp.communicate()[0] #qhull output to Vc
ks = Vc.split('\n')
ks = string.join(ks[2:],';') #remove leading dimension output
k = mat(ks[:-1]) #convert to martrix with vertices

# k is a (n+1)x(p) matrix in the form [A b] (from qhull doc: Ax < -b is satisfied), thus;
A = k[:,:-1]
b = -k[:,-1]
s = -ones([k.shape[0],1])
