#!/usr/bin/env python

# vert2con.py
# Python implementation of vert2con.m by Michael Kleder (July 2005),
#  available: http://www.mathworks.co.jp/matlabcentral/fileexchange/7895-vert2con-vertices-to-constraints
# Converts sets of vertices to a list of constraints (of the feasible region)
#
# Author: Michael Kelder (Original)
#         Andre Campher (Python implementation)

# Dependencies: * qhull (libqhull5, qhull-bin)


# ============ MATLAB code ==================
# function [A,b] = vert2con(V)
# k = convhulln(V);
# c = mean(V(unique(k),:));
# V=V-repmat(c,[size(V,1) 1]);
# A  = NaN*zeros(size(k,1),size(V,2));
# rc=0;
# for ix = 1:size(k,1)
#     F = V(k(ix,:),:);
#     if rank(F,1e-5) == size(F,1)
#         rc=rc+1;
#         A(rc,:)=F\ones(size(F,1),1);
#     end
# end
# A=A(1:rc,:);
# b=ones(size(A,1),1);
# b=b+A*c';
# % eliminate dumplicate constraints:
# [null,I]=unique(num2str([A b],6),'rows');
# A=A(I,:); % rounding is NOT done for actual returned results
# b=b(I);
# return

from scipy import *
import string
import subprocess #to use qhull

V = mat('0 0 0; 0 1 1;0    0 1; 1 0    8; 7 8 9') #vertex points matrix to test

# ================ generate data file for qhull
#format:
# n //dimensions (columns of V)
# p //number of points (rows of V)
# v //points e.g. 1 3 12 \n 2 4 5

infile = open("qhullin",'w')
infile.write(str(V.shape[1])+"\n"+str(V.shape[0])+"\n")
for row in range(0,V.shape[0]):
	stringmat = ""
	for col in range(0,V.shape[1]):
		stringmat = string.join([stringmat,str(V[row,col])]," ")
	infile.write(stringmat+"\n")
infile.close()

# subprocess.call(["function","arguments"]) or subprocess.Popen('function expression', shell=True)
# run qhull with: qhull < data or cat data | qhull

#subprocess.call("qhull")
