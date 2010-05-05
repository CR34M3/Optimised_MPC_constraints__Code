#!/usr/bin/env python

# uniqmat.py
# Removes duplicate entries from a matrix
#	INPUT:  A - matrix (possibly containing duplicates)
#		t - tolerance (float) [positive]
#	OUTPUT: B - matrix (A with the duplicates removed)
#
# Author: Andre Campher

from scipy import vstack,zeros
from numpy import all as nall

def uniqm(A,t):
	B = zeros((1,A.shape[1]))
	for r1 in range(0,A.shape[0]):
		dup = False
		row = A[r1,:]
		for r2 in range(r1+1,A.shape[0]):
			rowt = A[r2,:]
			if nall(abs(row-rowt)<t):
				dup = True #is duplicate
				break
		if not dup:
			B = vstack((B,row))
	B = B[1:,:]
	return B
