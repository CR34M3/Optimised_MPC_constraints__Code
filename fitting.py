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
        fb = 2 * r_[fixb, b]  # double the box size
        V = con2vert(fA, fb)[0]
    #return vol (normal or fixed) and vertices (normal or fixed) 
    return qhull(V, "FA"), V

def splitAb(inAb, nd):
    """Split input Ab array (1d) into separate A and b."""
    tmpAb = inAb.reshape(len(inAb)/(nd + 1), nd + 1)
    A = tmpAb[:, :-1]
    b = array([tmpAb[:, -1]]).T
    return A, b

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
    
    #### Objective fn
    def objfn(Ab, *args):
        """Volume objective function."""
        initcs = args[0]
        A, b = splitAb(Ab, initcs.nd)
        try:
            cl = con2vert(A, b)[1]
        except NameError:
            cl = -1
        if cl:
            closed = 1
        else:
            closed = -1
        vol = tryvol(A, b, initcs)[0]
        if not cl:
            print vol
        return closed * -vol
    #### Constraints (generate ranges)
    Abranges = ((-1, 1), (-1, 1), (-1, 1), (-1, 1), (-1, 1), (-1, 1), (-1, 1), (-1, 1), (-1, 1))
#    for k in range(9-1):
#        Abranges = Abranges, (-1, 1)
    #### Maximise volume
    optAb = optimize.brute(objfn, Abranges, args=[cset], Ns=3, finish=None)
    return optAb

if __name__ == "__main__":
    from pylab import plot, show
# ===== TEST 1
    v = array([[0, 0], [10, 0], [10, 10], [0, 10]])
    initcset = ConSet(v)
    Abopt = fitshape(3, initcset)
    tA, tb = splitAb(Abopt, initcset.nd)
    ts = -ones(tb.shape)
    optset = ConSet(tA, ts, tb)
    print optset.vert
    print optset.allinside(initcset)
    vp2 = vstack([optset.vert, optset.vert[0, :]])
    vp = vstack([v, v[0, :]])
    vst = array([[1, 1], [1, 9], [9, 1], [1, 1]])
    plot(vp[:, 0], vp[:, 1], 'b', vp2[:, 0], vp2[:, 1], 'r')
    show()
# ======= TEST 2
#    from auxfuns import mat2ab
#    AISA, AISs, AISb = mat2ab(array([[1., 0., 1., 0.],
#                                     [0., 1., 1., 0.],
#                                     [0.5, 0.5, -1, 1.]]))
#    AIS = ConSet(AISA, AISs, AISb)
#    vp2 = vstack([AIS.vert, AIS.vert[0, :]])
#    plot(vp2[:, 0], vp2[:, 1], 'r')
#    show()