#!/usr/bin/env python
"""Functions to optimally fit one 'shape' into another."""
from scipy import array, mat, optimize, vstack, eye, hstack, ones, tile, sqrt
from scipy import power, real
import numpy
import random
from conclasses import ConSet
from convertfuns import con2vert
from auxfuns import qhull

def tryvol(A, b, cs):
    """
    Try to determine volume of feasible region, otherwise impose box 
    constraint.
    """
    try:
        V = con2vert(A, b) # get vertices for constraint set iteration
    except: #catch unbound shapes
        #create large box around original shape
        fixA = vstack((eye(cs.nd), -eye(cs.nd)))
        fixb = vstack((mat([numpy.max(xv) 
                            for xv in cs.vert.transpose()]).transpose(),
                      -mat([numpy.min(xv) 
                            for xv in cs.vert.transpose()]).transpose()))
        fA = vstack((fixA, A))
        fb = vstack((fixb, b))
        V = con2vert(fA, fb)
    #return vol (normal or fixed) and vertices (normal or fixed) 
    return qhull(V, "FA"), V

def splitab(inAb, nd):
    """Split input Ab matrix into separate A and b."""
    nd1 = nd + 1
    #split Ab into A and b (assuming all -1 for s)
    Ab = mat(inAb)
    tempAb = vstack([Ab[:, p*nd1:p*nd1 + nd1]
                     for p in range(0, Ab.shape[1]/nd1)])
    return tempAb[:, :-1], tempAb[:, -1] #return A and b

def genstart(cs, ncon):
    """
    Generate a starting shape (for optimisation) with given number of faces.
    """
    #1. Determine centre of initial region, cscent
    c0 = ones((1, cs.vert.shape[1]))
    def cntrobjfn(cntr):
        """Objective function to determine center of shape."""
        cmat = tile(cntr, (cs.vert.shape[0], 1))
        #distance between vertices and centre
        dist = sqrt(power((cs.vert - cmat), 2)*mat(ones((cmat.shape[1], 1)))) 
        return real(sum(power(dist, 2)))
    cscent = optimize.fmin(cntrobjfn, c0)
    #2. Generate ncon-1 Gaussian vectors
    spherevecs = ones((ncon, cs.vert.shape[1]))
    for rows in range(spherevecs.shape[0] - 1):
        for cols in range(spherevecs.shape[1]):
            spherevecs[rows, cols] = random.gauss(0, 0.33)  # TODO: check sigma
    #3. Determine resultant of vectors, add last vector as mirror resultant
    spherevecs[-1, :] = -sum(spherevecs[:-1, :])
    srad = 1.0
    # points on sphere
    spherepts = spherevecs.T*(srad/sqrt(sum((spherevecs.T)**2)))  
    spherepts = spherepts.T + cscent  # move to center of constraint set
    #4. Generate tangent planes on sphere at points, convert to inequalities
    A = 2*spherepts
    b = array([2*sum(spherepts.T**2)]).T
    s = array(-ones(b.shape))
    #5. Check if vertices of new feasible region is within initial shape
    #6. Optimise sphere-radius, r, to have all points within initial shape
    spset = ConSet(A,s,b)
    return spherepts

def fitshape(cset, ncon, sp):
    """
    Fit a constraint set (specified by the number of constraints) within an existing
    constraint set.
     cset - [ConSet] existing constraint set
     ncon - [int] number of constraints to fit
    """
    #### State problem
    # ncongiven > ncon >= nD+1
    # Constraints are bounding
    # All vertices within constraint set
    #### Define parameters
    #    sp = #starting point - combined Ab matrix to optimise
    #### Helping fns
    #### Objective fn
    def objfn(Ab, *args):
        """Volume objective function."""
        initcs = args[0]
        A, b = splitab(Ab, initcs.nd)
        return -tryvol(A, b, initcs)[0] #-volume to minimise 
    #### Constraints    
    def ieconsfn(Ab, *args):
        """Optimiser inequality constraint function."""
        initcs = args[0]
        A, b = splitab(Ab, initcs.nd)
        # get vertices for constraint set iteration
        V = tryvol(A, b, initcs)[1].transpose()
        #constraint checking for vertices
        conscheck = -(initcs.A*V - initcs.b)
        #compress to single line  
        return hstack([conscheck[x, :] 
                       for x in range(0, conscheck.shape[0])])         
    #### Maximise volume
    return optimize.fmin_slsqp(objfn, sp, f_ieqcons=ieconsfn, args=[cset])

if __name__ == "__main__":
    v = array([[1, 0], [0, 1], [0, 0], [1, 1]])
    initcset = ConSet(v)
    print genstart(initcset, 4)