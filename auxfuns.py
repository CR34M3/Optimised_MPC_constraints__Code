#!/usr/bin/env python

from scipy import zeros
import numpy

# uniqm ====================================================================================================
# Removes duplicate entries from a matrix
#    INPUT:  A - matrix (possibly containing duplicates)
#        t - tolerance (float) [positive]
#    OUTPUT: matrix (A with the duplicates removed)
def uniqm(A,t):
    Nrows = A.shape[0]
    uniquerows = [r1 for r1 in range(Nrows)
              if not all(numpy.all(abs(A[r1,:]-A[r2,:])<t)
                 for r2 in range(r1+1,Nrows))] 
    return A[uniquerows,:].copy()


# mat2ab ====================================================================================================
# Transform to standard Ax<b notation using sign vector
#    INPUT: Asbmat - inequality matrix in the form Ax<b, matrix = [A s b] with s the sign vector [1:>, -1:<]
def mat2ab(Asbmat):
    A = Asbmat[:,:-2]
    s = Asbmat[:,-2]
    b = Asbmat[:,-1]
    A = numpy.multiply(A,-s)
    b = numpy.multiply(b,-s)
    return A,s,b
