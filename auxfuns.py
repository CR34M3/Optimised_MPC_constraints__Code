#!/usr/bin/env python
"""Auxiliary functions to manipulate data-types and and change data-formats."""
from scipy import zeros
import numpy

def uniqm(A,t):
    """
    Return input matrix with duplicate entries removed.
     A - [matrix] input matrix (possibly containing duplicates)
     t - [float]  tolerance (default=0.01)
    """
    Nrows = A.shape[0]
    uniquerows = [r1 for r1 in range(Nrows)
              if not any(numpy.all(abs(A[r1,:]-A[r2,:])<t)
             for r2 in range(r1+1,Nrows))]
    return A[uniquerows,:].copy()


def mat2ab(Asbmat):
    """
    Transform [A s b]-form matrix to standard Ax<b notation.
     Asbmat - [matrix] inequality matrix in the form Ax<b, matrix = [A s b] with s the sign vector [1:>, -1:<]
    """    
    A = Asbmat[:,:-2]
    s = Asbmat[:,-2]
    b = Asbmat[:,-1]
    A = numpy.multiply(A,-s)
    b = numpy.multiply(b,-s)
    return A,s,b
