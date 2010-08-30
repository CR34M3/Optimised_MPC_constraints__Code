#!/usr/bin/env python
"""Functions to calculate vertices from constraints and vice versa."""
from scipy import *
from numpy import linalg, matlib
import subprocess #to use qhull
from gendatafile import *
from auxfuns import *
from os import remove
from sys import exit

def vert2con(V):
    """
    Convert sets of vertices to a list of constraints (of the feasible region).
    Return A, b and s of the set;  Ax < b   (s is the sign-vector [to be used later])
    """
    # Dependencies: * qhull (libqhull5, qhull-bin)
    #               * scipy
    #               * gendatafile
    k = qhull(V,"n") #convert to martrix with vertices
    # k is a (n+1)x(p) matrix in the form [A b] (from qhull doc: Ax < -b is satisfied), thus;
    A = k[:,:-1]
    b = -k[:,-1]
    s = -ones([k.shape[0],1])
    return hstack([A,s,b])

def con2vert(A,b):
    """Convert sets of constraints to a list of vertices (of the feasible region)."""
    # Python implementation of con2vert.m by Michael Kleder (July 2005),
    #  available: http://www.mathworks.com/matlabcentral/fileexchange/7894-con2vert-constraints-to-vertices
    # Author: Michael Kelder (Original)
    #         Andre Campher (Python implementation)
    #
    # Dependencies : - Scipy
    #         - Numpy
    #         - gendatafile
    #         - uniqmat
    c = linalg.lstsq(A,b)[0]
    b = b-A*c
    D = A / matlib.repmat(b,1,A.shape[1])
    Dtest = vstack((D,zeros([1,D.shape[1]])))

#== Volume error check ==
    VolDt = qhull(Dtest,"FA") #get volume of D-hull
    VolD = qhull(D,"FA") #get volume of D-hull

    if VolDt > VolD:
        print 'error : Non-bounding constraints detected (consider box constraints on variables). Exiting...'
        exit(1)
#== ==

    fmatv = qhull(D,"Ft") #vertices on facets

    G  = zeros((fmatv.shape[0],D.shape[1]));
    for ix in range(0,fmatv.shape[0]):
        F = D[fmatv[ix,:],:].squeeze()
        G[ix,:] = linalg.lstsq(F,ones((F.shape[0],1)))[0].transpose()

    V = G + matlib.repmat(c.transpose(),G.shape[0],1)
    ux = uniqm(V,0.01)

    return ux

def con2pscon(A,s,b):
    """Convert a set of constraints to a set of pseudo constraints using only high/low limits."""
    #Take Asb matrix as input
    #Determine necessary conversions
    ###Check for single row entries (and ignore)
    checkmat = zeros((1,A.shape[1]))
    tempcvmat = zeros((0,A.shape[1]))
    tempA = zeros((0,A.shape[1]))
    tempsb = zeros((0,2))
    origA = zeros((0,A.shape[1]))
    origsb = zeros((0,2))
    nconv = 0
    for r in range(0,A.shape[0]):
        if sum(A[r,:]==checkmat) < A.shape[1]-1:  #if less than n-1 zeros, convert
            nconv = nconv + 1 #flag as conversion
            tempA = vstack((tempA,A[r,:]))
            tempsb = vstack((tempsb,hstack((s[r],b[r]/abs(b[r]))))) #normalise b
            tempcvmat = vstack((tempcvmat,A[r,:]/abs(b[r]))) #add to tempcvmat
        else:                                #if n-1 zeros, just build
            origA = vstack((origA,A[r,:])) #keep values and positions
            origsb = vstack((origsb,hstack((s[r],b[r]))))
    if nconv>0:
        padcol = zeros((origA.shape[0],nconv))
        padrow = zeros((nconv,origA.shape[1]))
        padvar = eye(nconv)
        psA = vstack((hstack((origA,padcol)),hstack((padrow,padvar))))
        pssb = vstack((origsb,tempsb))
        convertmat = vstack((ones((A.shape[0]-nconv,A.shape[1])),tempcvmat))
        return hstack([psA,pssb]),convertmat
    else:   #if no conversions were made
        return hstack([A,s,b]),ones(A.shape)
    #Express cset as combination of linear inequalities with high/low limits
    

if __name__ == "__main__":
    import doctest
    doctest.testfile("tests/convertfunstests.txt")
    
    
#TODO con2vert =====
# error-checking
#    - fix volume check (for redundant constraints)
#TODO general
# check that floating-point math is used