#!/usr/bin/env python
"""Functions to optimally fit one 'shape' into another."""
from scipy import array, optimize, eye, c_, r_, ones, sqrt
from scipy import vstack, linalg
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
        V = con2vert(A, b)[0] # get vertices for constraint set iteration
    except NameError: #catch unbound shapes
        #create large box around original shape
        fixA = r_[eye(cs.nd), -eye(cs.nd)]
        fixb = r_[[numpy.max(cs.vert[:, x]) for x in range(cs.vert.shape[1])],
                  [-numpy.min(cs.vert[:, x]) for x in range(cs.vert.shape[1])]]
        fixb = array([fixb]).T
        fA = r_[fixA, A]
        fb = 2. * r_[fixb, b]  # double the box size
        V = con2vert(fA, fb)[0]
    #return vol (normal or fixed) and vertices (normal or fixed) 
    return qhull(V, "FS"), V

def splitAb(inAb, nd):
    """Split input Ab array (1d) into separate A and b."""
    tmpAb = inAb.reshape(len(inAb)/(nd + 1), nd + 1)
    A = tmpAb[:, :-1]
    b = array([tmpAb[:, -1]]).T
    return A, b

def genstart(ncon, cs):
    """
    Generate a starting shape (for optimisation) with given number of faces.
    """
    #1. Determine centre of initial region, cscent
    cscent = sum(cs.vert)/len(cs.vert)  # assuming the region is convex
    #2. Generate ncon-1 Gaussian vectors
    spherevecs = ones((ncon, cs.vert.shape[1]))
    for rows in range(spherevecs.shape[0] - 1):
        for cols in range(spherevecs.shape[1]):
            spherevecs[rows, cols] = random.gauss(0., 0.33)  # TODO: check sigma
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
        s = array(ones(b.shape))
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
#    spset = ConSet(array([[8, 4], [8.2, 4.5], [8.5, 3.8]]))
    snorm = linalg.norm(spset.b)
    sp = c_[spset.A, spset.b]/snorm # starting point - combined Ab matrix to optimise
    vp2 = vstack([spset.vert, spset.vert[0, :]])
    plot(vp2[:, 0], vp2[:, 1], 'g--')
    
    #### Objective fn
    def objfn(Ab, *args):
        """Volume objective function."""
        initcs = args[0]
        A, b = splitAb(Ab, initcs.nd)
        vol, V = tryvol(A, b, initcs)
        Pv = 200.
        Pn = 100.
        #Penalties
        # large b norm
        bnorm = abs(linalg.norm(b) - 1)
        # points outside of init space
        iterset = ConSet(V)
        outnorm = linalg.norm(iterset.allinside(initcs)[1])
        # open shape
        cl = con2vert(A, b)[1]
        if cl:
            closed = 1
        else:
            closed = -1
        return (-vol*closed) + Pn*(bnorm**2) + Pv*(outnorm**3)
    #### Maximise volume
    optAb = optimize.fmin(objfn, sp, args=[cset], maxiter=20000, disp=True)
    return optAb

if __name__ == "__main__":
    from pylab import plot, show
    print "GO!"
    v = array([[0, 0], [10, 0], [10, 10], [0, 10]])
    initcset = ConSet(v)
    Abopt = fitshape(3, initcset)
    tA, tb = splitAb(Abopt, initcset.nd)
    ts = -ones(tb.shape)
    optset = ConSet(tA, ts, tb)
    vp2 = vstack([optset.vert, optset.vert[0, :]])
    vp = vstack([v, v[0, :]])
    vst = array([[1, 1], [1, 9], [9, 1], [1, 1]])
    plot(vp[:, 0], vp[:, 1], 'b', vp2[:, 0], vp2[:, 1], 'r')
    show()