#!/usr/bin/env python
"""Auxiliary functions to manipulate data-types and and change data-formats."""
from scipy import zeros,mat,optimize,ceil,tile,multiply
from gendatafile import *
import numpy
from os import remove
import subprocess #to use qhull

def uniqm(A,t=0.01):
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

def qhull(V,qstring):
    """
    Use qhull to determine convex hull / volume / normals.
     V - [matrix] vertices
     qstring - [string] arguments to pass to qhull
    """
    # subprocess.call(["function","arguments"]) or subprocess.Popen('function expression', shell=True)
    # run qhull with: qhull < data or cat data | qhull
    filename = genfile(V)
    qstringfull = "qhull " + qstring + " < " + filename
    qhullp = subprocess.Popen(qstringfull, shell=True, stdout=subprocess.PIPE) #calc convex hull and get normals
    Vc = qhullp.communicate()[0] #qhull output to Vc
    remove(filename)
        
    if qstring == "FA": #calc summary and volume
        ks = Vc.split('\n')[-3]
        Vol = float(ks.split(' ')[-1]) #get volume of D-hull
        return Vol
    if qstring == "Ft": #calc vertices and facets
        ks = Vc.split('\n')
        fms = int(ks[1].split(' ')[1]) #get size of facet matrix
        fmat = ks[-fms-1:-1]
        fmat = mat(';'.join(fmat)) #generate matrix
        fmatn = fmat[:,0] #number of points on facets
        fmatv = fmat[:,1:] #vertices on facets
        return fmatv
    if qstring == "n": #calc convex hull and get normals
        ks = ';'.join(Vc.split('\n')[2:]) #remove leading dimension output
        k = mat(ks[:-1]) #convert to martrix with vertices
        return k
    else:
        exit(1)

def fitshape(cset,ncon):
    """
    Fit a constraint set (specified by the number of constraints) within an existing
    constraint set.
     cset - [conset] existing constraint set
     ncon - [int] number of constraints to fit
    """
#### State problem
   # ncongiven > ncon >= nD+1
   # (current dirty method - specify edges as x*dimensions)
   # subject to
    # Using the Euler characteristic for spherical polyhedra: V-E+F = 2 (check dimensions)
    # Within constraint set
#### Define parameters
    nv = 2 - ncon + (ceil(ncon/cset.nd)*cset.nd) #determine number of vertices (noting E > F)
    sp = cset.vert[0:nnv,:] #use the first nv vertices (of initial space) as starting points
#### Objective fn
    def objfn(V,*args):
        shapevol = args[0]
        return shapevol - qhull(V,"FA")
    def ieconsfn(V,*args):
        initcset = args[1]
        tmpV = zeros((V.shape[0]*initcset.A.shape[0],V.shape[1])) #temp V matrix for all inequalities
        #convert to standard form (all s = -1)
        initcset.A = multiply(initcset.A,-1*initcset.s)
        initcset.b = multiply(initcset.b,-1*initcset.s)
        tmpA = tile(initcset.A,(V.shape[0],1))
        tmpb = tile(initcset.b,(V.shape[0],1))
        #newvertcons = initcset*
        return newvertcons
#### Maximise volume
    sol = optimize.fmin_slsqp(objfn,sp,f_ieqcon=[])    

if __name__ == "__main__":
    import doctest
    doctest.testfile("tests/auxfunstests.txt")
    
#TODO - auxfuns
# qhull error handling   