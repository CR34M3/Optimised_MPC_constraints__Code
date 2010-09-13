#!/usr/bin/env python
"""Functions to optimally fit one 'shape' into another."""
from scipy import array, mat, optimize, eye, c_, r_, ones, tile, sqrt
from scipy import power, real, size
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
    except NameError: #catch unbound shapes
        #create large box around original shape
        fixA = r_[eye(cs.nd), -eye(cs.nd)]
        fixb = r_[[numpy.max(cs.vert[:, x]) for x in range(cs.vert.shape[1])],
                  [-numpy.min(cs.vert[:, x]) for x in range(cs.vert.shape[1])]]
        fixb = array([fixb]).T
        fA = r_[fixA, A]
        fb = 2 * r_[fixb, b]  # double the box size
        V = con2vert(fA, fb)
    #return vol (normal or fixed) and vertices (normal or fixed) 
    return qhull(V, "FA"), V

def splitAb(inAb, nd):
    """Split input Ab matrix into separate A and b."""
    tmpAb = r_[[inAb[2*x:2*x + nd + 1] for x in range(inAb.shape[0]/(nd+1))]]
    A = tmpAb[:, :-1]
    b = array([tmpAb[:, -1]]).T
    return A, b

def genstart(ncon, cs):
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
    #3. Determine resultant of vectors, add last vector as mirror of resultant
    spherevecs[-1, :] = -sum(spherevecs[:-1, :])
    # points on sphere
    def spherefn(srad):
        """
        Function to create points on a sphere (of radius, srad), convert them to
        halfspaces and convert the halfspaces to a constraint set.
        """
        spherepts = spherevecs.T*(srad/sqrt(sum((spherevecs.T)**2)))  
        spherepts = spherepts.T + cscent  # move to center of constraint set
        #4. Generate tangent planes on sphere at points, convert to inequalities
        A = -(spherepts - cscent)
        b = array([(sum((A*spherepts).T))]).T
        s = array(-ones(b.shape))
        return ConSet(A, s, b)  # constraints around sphere
    
    #6. Optimise sphere-radius, r, to have all points within initial shape
    startrad = 1.0
    if spherefn(startrad).allinside(cs)[0]:
        while spherefn(startrad).allinside(cs)[0]:
            startrad = startrad * 1.1
        finrad = startrad / 1.1    
    elif not spherefn(startrad).allinside(cs)[0]:
        while not spherefn(startrad).allinside(cs)[0]:
            startrad = startrad / 1.1
        finrad = startrad
    return spherefn(finrad)

def fitshape(ncon, cset):
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
    spset = genstart(ncon, cset)  
    sp = c_[spset.A, spset.b] # starting point - combined Ab matrix to optimise
    #### Objective fn
    def objfn(Ab, *args):
        """Volume objective function."""
        initcs = args[0]
        A, b = splitAb(Ab, initcs.nd)
        return -tryvol(A, b, initcs)[0] #-volume to minimise 
    #### Constraints    
    def ieconsfn(Ab, *args):
        """Optimiser inequality constraint function."""
        initcs = args[0]
        A, b = splitAb(Ab, initcs.nd)
        # get vertices for constraint set iteration
        V = tryvol(A, b, initcs)[1]
        iterset = ConSet(V)
        #constraint checking for vertices
        ineqs = iterset.allinside(initcs)[1] 
        return ineqs.reshape(size(ineqs), )
    #### Maximise volume
    return optimize.fmin_slsqp(objfn, sp, f_ieqcons=ieconsfn, args=[cset])

if __name__ == "__main__":
    v = array([[10, 0], [0, 10], [0, 0], [10, 10]])
    initcset = ConSet(v)
    tA, tb = splitAb(fitshape(3, initcset), initcset.nd)
    ts = -ones(tb.shape)
    print ConSet(tA, ts, tb).vert