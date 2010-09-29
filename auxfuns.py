#!/usr/bin/env python
"""Auxiliary functions to manipulate data-types and and change data-formats."""
from scipy import mat, array, ones
import numpy
import subprocess #to use qhull

def uniqm(A, t=1e-13):
    """
    Return input matrix with duplicate entries removed.
     A - [matrix] input matrix (possibly containing duplicates)
     t - [float]  tolerance (default=1e-13)
    """
    Nrows = A.shape[0]
    uniquerows = [r1 for r1 in range(Nrows)
                  if not any(numpy.all(abs(A[r1, :] - A[r2, :]) < t)
                             for r2 in range(r1 + 1, Nrows))]
    return A[uniquerows, :].copy()


def mat2ab(Asbmat):
    """
    Transform [A s b]-form matrix to standard Ax<b notation.
     Asbmat - [array] inequality matrix in the form Ax<b, matrix = [A s b] 
                        with s the sign vector [1:>, -1:<]
    """    
    stmp = array([Asbmat[:, -2]])
    b = Asbmat[:, -1]*-stmp
    A = Asbmat[:, :-2]*-stmp.T
    s = -ones(stmp.shape)
    return A, s.T, b.T

def qhullstr(V):
    """ 
    generate string qhull input format.
    
    yields a newline separated stirng of format:
        dimensions (columns of V)
        number of points (rows of V)
        one string for each row of V
    """
    V = numpy.array(V)
    return "%i\n%i\n" % (V.shape[1], V.shape[0]) \
           + "\n".join(" ".join(str(e) for e in row) for row in V) 

def qhull(V, qstring):
    """
    Use qhull to determine convex hull / volume / normals.
     V - [matrix] vertices
     qstring - [string] arguments to pass to qhull
    """
    try:
        qhullp = subprocess.Popen(["qhull", qstring],
                              stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        Vc = qhullp.communicate(qhullstr(V))[0] #qhull output to Vc
        
        if qstring == "FS": #calc area and volume
            ks = Vc.split('\n')[-2]
            Vol = float(ks.split(' ')[-2]) #get volume of D-hull
            return Vol
        elif qstring == "Ft": #calc vertices and facets
            ks = Vc.split('\n')
            fms = int(ks[1].split(' ')[1]) #get size of facet matrix
            fmat = ks[-fms-1:-1]
            fmat = mat(';'.join(fmat)) #generate matrix
            fmatv = fmat[:, 1:] #vertices on facets
            return array(fmatv)
        elif qstring == "n": #calc convex hull and get normals
            ks = ';'.join(Vc.split('\n')[2:]) #remove leading dimension output
            k = mat(ks[:-1]) #convert to martrix with vertices
            return array(k)
        else:
            exit(1)
    except:
        raise NameError('QhullError')
    
if __name__ == "__main__":
    import doctest
    doctest.testfile("tests/auxfunstests.txt")   
    
    
#TODO - auxfuns
# qhull error handling   
# fix auxfuns tests