#!/usr/bin/env python

# uniqmat.py
# Removes duplicate entries from a matrix
#	INPUT:  A - matrix (possibly containing duplicates)
#		t - tolerance (float) [positive]
#	OUTPUT: B - matrix (A with the duplicates removed)
#
# Author: Andre Campher

from scipy import zeros
import numpy

def uniqm(A,t):
	Nrows = A.shape[0]
	uniquerows = [r1 for r1 in range(Nrows)
		      if not all(numpy.all(abs(A[r1,:]-A[r2,:])<t)
				 for r2 in range(r1+1,Nrows))] 
	return A[uniquerows,:].copy()
