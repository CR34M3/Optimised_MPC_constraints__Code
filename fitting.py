#!/usr/bin/env python
"""Functions to optimally fit one 'shape' into another."""
from scipy import array, mat, optimize, eye, c_, r_, ones, tile, sqrt
from scipy import power, real, size, vstack, linalg, average
from scipy.spatial.distance import pdist
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
    c0 = ones((1, cs.vert.shape[1]))
    def cntrobjfn(cntr):
        """Objective function to determine center of shape."""
        cmat = tile(cntr, (cs.vert.shape[0], 1))
        #distance between vertices and centre
        dist = sqrt(power((cs.vert - cmat), 2)*mat(ones((cmat.shape[1], 1)))) 
        return real(sum(power(dist, 2)))
    cscent = optimize.fmin(cntrobjfn, c0, disp = False)
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
#    spset = genstart(ncon, cset) 
    spset = ConSet(array([[1, 1], [1, 9], [9, 1]]))
    sp = c_[spset.A, spset.b] # starting point - combined Ab matrix to optimise
    
    #### Objective fn
    def objfn(Ab, *args):
        """Volume objective function."""
        initcs = args[0]
        A, b = splitAb(Ab, initcs.nd)
        vol, V = tryvol(A, b, initcs)
        dists = pdist(V)  # pairwise distances between vertices
        return -vol*(average(dists))
    #### Constraints
    def eqconsfn(Ab, *args):
        """Optimiser equality constraint function."""
        initcs = args[0]
        b = splitAb(Ab, initcs.nd)[1]
        return array([linalg.norm(b) - 100])
    def ieqconsfn(Ab, *args):
        """Optimiser inequality constraint function."""
        initcs = args[0]
        A, b = splitAb(Ab, initcs.nd)
        # get vertices for constraint set iteration
        V = tryvol(A, b, initcs)[1]
        iterset = ConSet(V)
        #constraint checking for vertices
        ineqs = iterset.allinside(initcs)[1][:,1]
        return ineqs.reshape(size(ineqs), )
    #### Maximise volume
    optAb = optimize.fmin_slsqp(objfn, sp, f_eqcons=eqconsfn, 
                               f_ieqcons=ieqconsfn, args=[cset], iprint=3)
    return optAb

if __name__ == "__main__":
    from pylab import plot, show
    v = array([[0, 0], [10, 0], [10, 10], [0, 10]])
    initcset = ConSet(v)
    Abopt = fitshape(3, initcset)
    print Abopt
    tA, tb = splitAb(Abopt, initcset.nd)
    ts = -ones(tb.shape)
    optset = ConSet(tA, ts, tb)
    vp2 = vstack([optset.vert, optset.vert[0, :]])
    vp = vstack([v, v[0, :]])
    vst = array([[1, 1], [1, 9], [9, 1], [1, 1]])
    plot(vp[:, 0], vp[:, 1], 'b', vp2[:, 0], vp2[:, 1], 'r', vst[:, 0], vst[:, 1], 'g')
    show()