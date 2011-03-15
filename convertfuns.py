#!/usr/bin/env python
"""Functions to calculate vertices from constraints and vice versa."""
from scipy import ones, mat, vstack, zeros, hstack, eye, array, dot, tile, sum
from scipy import all as sciall
from numpy import linalg, matlib
from auxfuns import uniqm, qhull

def vert2con(V):
    """
    Convert sets of vertices to a list of constraints (of the feasible region).
    Return A, b and s of the set;  Ax < b   (s is the sign-vector [to be used later])
    vert2con always closes the shape and generates inequalities accordingly.
    """
    # Dependencies: * qhull (libqhull5, qhull-bin)
    #               * scipy
    k = qhull(V,"n") #convert to martrix with vertices
    # k is a (n+1)x(p) matrix in the form [A b] (from qhull doc: Ax < -b is
    # satisfied), thus;
    A = k[:, :-1]
    b = array([-k[:, -1]]).T
    s = -ones(b.shape)
    return A, s, b

def con2vert(A, b):
    """
    Convert sets of constraints to a list of vertices (of the feasible region).
    If the shape is open, con2vert returns False for the closed property.
    """
    # Python implementation of con2vert.m by Michael Kleder (July 2005),
    #  available: http://www.mathworks.com/matlabcentral/fileexchange/7894
    #  -con2vert-constraints-to-vertices
    # Author: Michael Kelder (Original)
    #         Andre Campher (Python implementation)
    c = linalg.lstsq(mat(A), mat(b))[0]
    btmp = mat(b)-mat(A)*c
    D = mat(A)/matlib.repmat(btmp, 1, A.shape[1])

    fmatv = qhull(D, "Ft") #vertices on facets

    G  = zeros((fmatv.shape[0], D.shape[1]))
    for ix in range(0, fmatv.shape[0]):
        F = D[fmatv[ix, :], :].squeeze()
        G[ix, :] = linalg.lstsq(F, ones((F.shape[0], 1)))[0].transpose()

    V = G + matlib.repmat(c.transpose(), G.shape[0], 1)
    ux = uniqm(V)

    eps = 1e-13
    Av = dot(A, ux.T)
    bv = tile(b, (1, ux.shape[0]))
    closed = sciall(Av - bv <= eps)

    return ux, closed

def con2pscon(cset, G, type):
    """
    Convert a constraint set to another constraint with only
    high/low limits.
        cset - constraint set to convert [ConSet]
        G - model of original set [array]
        type - type of model conversion to be done [string]
                'i'nput constraint set
                'o'utput constraint set
    """
    #Check for single variable entries
    checkmat = sum(zeros(cset.A.shape) == cset.A, axis=1)
    #Determine number of conversions
    nconv = sum(checkmat < cset.nd-1)
    if nconv:
        #keep original high/low limits (s remains unchanged)
        keepA = vstack([cset.A[x, :] for x in range(cset.A.shape[0]) 
                        if checkmat[x]])
        keepb = vstack([cset.b[x, :] for x in range(cset.b.shape[0]) 
                        if checkmat[x]])
        fixA = vstack([cset.A[x, :] for x in range(cset.A.shape[0]) 
                        if not checkmat[x]])
        fixb = vstack([cset.b[x, :] for x in range(cset.b.shape[0]) 
                        if not checkmat[x]])
        tempA = zeros((keepA.shape[0]+nconv, keepA.shape[1]+nconv))
        tempA[:keepA.shape[0], :keepA.shape[1]] = keepA
        tempA[-nconv:, -nconv:] = eye(nconv)
        tempb = vstack((keepb, fixb))
        if type in "iI":
            tempG = vstack((G, fixA))
        elif type in "oO":
            fixG = sum(G, axis=0)*fixA
            tempG = vstack((G, fixG))
        else:
            tempG = G
        return tempA, cset.s, tempb, tempG
    else:
        return cset.A, cset.s, cset.b, G

if __name__ == "__main__":
    import doctest
    doctest.testfile("tests/convertfunstests.txt")
    
#TODO con2vert =====
# error-checking
#    - fix volume check (for redundant constraints)
#TODO general
# check that floating-point math is used
# fix output of con2pscon