#!/usr/bin/env python
"""Auxiliary functions to manipulate data-types and and change data-formats."""
from scipy import zeros,mat,optimize,ceil,tile,multiply,vstack
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
     Asbmat - [matrix] inequality matrix in the form Ax<b, matrix = [A s b] 
                        with s the sign vector [1:>, -1:<]
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
    from convertfuns import con2vert
#### State problem
   # ncongiven > ncon >= nD+1
   # subject to
    # Constraints are bounding
    # All vertices within constraint set
#### Define parameters
#    sp = #starting point - combined Ab matrix to optimise
#### Objective fn
    def objfn(Ab, *args):
        A = Ab[:,:-1] #split Ab into A and b (assuming all -1 for s)
        b = Ab[:,-1]
        initcset = args[0]
        try:
            V = con2vert(A,b) # get vertices for constraint set iteration
        except: #catch unbound shapes
            #create large box around original shape
            fixA = vstack((eye(initcset.nd),-eye(initcset.nd)))
            fixb = vstack((mat([numpy.max(xv) for xv in initcset.vert.transpose()]).transpose(),
                          -mat([numpy.min(xv) for xv in initcset.vert.transpose()]).transpose()))
            fA = vstack((fixA,A))
            fb = vstack((fixb,b))
            V = con2vert(fA,fb)
        return 1 - initcset.vol/qhull(V,"FA") #maximise fraction orig vol / new vol 
#### Constraints    
    def ieconsfn(Ab,*args):
        A = Ab[:,:-1] #split Ab into A and b (assuming all -1 for s)
        b = Ab[:,-1]
        V = con2vert(A,b).transpose() # get vertices for constraint set iteration
        initcset = args[0]
        #constraint checking for vertices
        return -(initcset.A*V - initcset.b) #negative as to make valid entries > 0
    
#### Maximise volume
    return optimize.fmin_slsqp(objfn,sp,f_ieqcon=ieconsfn,args=cset)

if __name__ == "__main__":
#    import doctest
#    doctest.testfile("tests/auxfunstests.txt")
        
    
#TODO - auxfuns
# qhull error handling   