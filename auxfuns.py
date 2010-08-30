#!/usr/bin/env python
"""Auxiliary functions to manipulate data-types and and change data-formats."""
from scipy import zeros,mat,optimize,ceil,tile,multiply,vstack,eye,array,hstack
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


def fitshape(cset,ncon,sp):
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
#### Helping fns
    def splitab(inAb,nd):
        nd1 = nd+1
        #split Ab into A and b (assuming all -1 for s)
        Ab = mat(inAb)
        tempAb = vstack([Ab[:,p*nd1:p*nd1+nd1] for p in range(0,Ab.shape[1]/nd1)])
        return tempAb[:,:-1], tempAb[:,-1] #return A and b
    def tryvol(A,b,cs):
        try:
            V = con2vert(A,b) # get vertices for constraint set iteration
        except: #catch unbound shapes
            #create large box around original shape
            fixA = vstack((eye(cs.nd),-eye(cs.nd)))
            fixb = vstack((mat([numpy.max(xv) for xv in cs.vert.transpose()]).transpose(),
                          -mat([numpy.min(xv) for xv in cs.vert.transpose()]).transpose()))
            fA = vstack((fixA,A))
            fb = vstack((fixb,b))
            V = con2vert(fA,fb)
        return qhull(V,"FA"),V #return vol (normal or fixed) and vertices (normal or fixed) 
#### Objective fn
    def objfn(Ab,*args):
        initcset = args[0]
        A,b = splitab(Ab,initcset.nd)
        return -tryvol(A,b,initcset)[0] #-volume to minimise 
#### Constraints    
    def ieconsfn(Ab,*args):
        initcset = args[0]
        A,b = splitab(Ab,initcset.nd)
        V = tryvol(A,b,initcset)[1].transpose() # get vertices for constraint set iteration
        #constraint checking for vertices
        conscheck = -(initcset.A*V - initcset.b) #negative as to make valid entries > 0
        return hstack([conscheck[x,:] for x in range(0,conscheck.shape[0])]) #compress to single line
        
#### Maximise volume
    return optimize.fmin_slsqp(objfn,sp,f_ieqcons=ieconsfn,args=[cset])

if __name__ == "__main__":
#    import doctest
#    doctest.testfile("tests/auxfunstests.txt")
    from conclasses import conset
    ineqs = mat('1 0 -1 1;-1 0 -1 0;0 1 -1 1;0 -1 -1 0')
    cs = conset(ineqs)
    print fitshape(cs,4,mat('1 0 0.8 -1 0 -0.2 0 1 0.8 0 -1 -0.2'))
    
        
    
#TODO - auxfuns
# qhull error handling   