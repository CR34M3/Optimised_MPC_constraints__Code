#!/usr/bin/env python
"""Functions to optimally fit one 'shape' into another."""
from scipy import array, optimize, eye, c_, r_, ones, sqrt
from scipy import linalg, sum, size
import numpy
import random
from conclasses import ConSet
from convertfuns import con2vert
from auxfuns import qhull, splitAb

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
    #TODO: check box fitting
     
def genstart(shapetype, *args):
    """
    Generate a starting shape (for optimisation) with given number of faces.
        shapetype : [string] (r)ectangle, (a)rbitrary
        args -> (initial constraint set, number of faces)
    """
    if shapetype in 'rR':  # Rectangle
        cs = args[0]
        ncon = cs.nd*2
    elif shapetype in 'aA':  # Arbitrary shape
        cs = args[0]
        ncon = args[1]
    
    #1. Determine centre of initial region, cscent
    cscent = cs.cscent  # assuming the region is convex
    def spherefn(srad, cent, svecs):
        """
        Function to create points on a sphere (of radius, srad), convert them to
        halfspaces and convert the halfspaces to a constraint set.
        """
        spherepts = svecs.T*(srad/sqrt(sum((svecs.T)**2)))  
        spherepts = spherepts.T + cent  # move to center of constraint set
        #4. Generate tangent planes on sphere at points, convert to inequalities
        A = -(spherepts - cent)
        b = array([(sum((A*spherepts).T, axis=0))]).T
        s = array(ones(b.shape))
        return ConSet(A, s, b)  # constraints around sphere

    if shapetype in 'rR':  # Rectangle
        Astart = r_[eye(cs.nd), -eye(cs.nd)]
        sstart = -ones((cs.nd*2, 1))
        bstart = lambda k: array(r_[cscent.T + k, -(cscent.T - k)]).reshape(size(cscent)*2,1)
        kstart = 1.0
        if ConSet(Astart, sstart, bstart(kstart)).allinside(cs)[0]:
            while ConSet(Astart, sstart, bstart(kstart)).allinside(cs)[0]:
                kstart = kstart * 1.01
            kfin = kstart / 1.01
        elif not ConSet(Astart, sstart, bstart(kstart)).allinside(cs)[0]:
            while not ConSet(Astart, sstart, bstart(kstart)).allinside(cs)[0]:
                kstart = kstart / 1.01
            kfin = kstart        
        return ConSet(Astart, sstart, bstart(kfin))
    elif shapetype in 'aA':  # Arbitrary shape
        #2. Generate ncon-1 Gaussian vectors
        spherevecs = ones((ncon, cs.vert.shape[1]))
        for rows in range(spherevecs.shape[0] - 1):
            for cols in range(spherevecs.shape[1]):
                spherevecs[rows, cols] = random.gauss(0, 0.33)  # TODO: check sigma                
        #3. Determine resultant of vectors, add last vector as mirror of resultant
        spherevecs[-1, :] = -sum(spherevecs[:-1, :], axis=0)
        # points on sphere
        #6. Optimise sphere-radius, r, to have all points within initial shape
        startrad = cs.vol()**(1./cs.nd)
        if spherefn(startrad, cscent, spherevecs).allinside(cs)[0]:
            while spherefn(startrad, cscent, spherevecs).allinside(cs)[0]:
                startrad = startrad * 1.01
            finrad = startrad / 1.01    
        elif not spherefn(startrad, cscent, spherevecs).allinside(cs)[0]:
            while not spherefn(startrad, cscent, spherevecs).allinside(cs)[0]:
                startrad = startrad / 1.01
            finrad = startrad
        return spherefn(finrad, cscent, spherevecs)

def fitshape(cset, spset, solver):
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
    snorm = linalg.norm(spset.b)
    sp = c_[spset.A, spset.b]/snorm # starting point - combined Ab matrix to optimise
    
    #### Objective fn (FMIN)
    def objfn(Ab, *args):
        """Volume objective function for fmin (Simplex)."""
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
        #outnorm = iterset.allinside(initcs)[2]
        # open shape
        cl = con2vert(A, b)[1]
        if cl:
            closed = 1
        else:
            closed = -1
        vol = tryvol(A, b, initcs)[0]
        return (-vol*closed) + Pn*(bnorm**initcs.nd) + Pv*(outnorm**(initcs.nd+3))
       
    #### Objective fn (SLSQP)
    def objfn2(Ab, *args):
        """Volume objective function for SLSQP."""
        initcs = args[0]
        A, b = splitAb(Ab, initcs.nd)
        cl = con2vert(A, b)[1]
        if cl:
            closed = 1
        else:
            closed = -1
        vol = tryvol(A, b, initcs)[0]
        return closed * -vol
    def eqconsfn(Ab, *args):
        """Optimiser equality constraint function."""
        initcs = args[0]
        b = splitAb(Ab, initcs.nd)[1]
        # get vertices
        bn = linalg.norm(b) - 1
        return array([bn])
    def ieqconsfn(Ab, *args):
        """Optimiser inequality constraint function."""
        initcs = args[0]
        A, b = splitAb(Ab, initcs.nd)
        # get vertices
        V = tryvol(A, b, initcs)[1]
        iterset = ConSet(V)
        # constraint checking
        if iterset.allinside(initcs):
            return 0
        else:
            return -1
    
    #### Maximise volume
    if solver in 'aA':
        optAb = optimize.fmin_slsqp(objfn2, sp, f_eqcons=eqconsfn,
                                    f_ieqcons=ieqconsfn, args=[cset], iprint=0)
    elif solver in 'bB':
        optAb = optimize.fmin(objfn, sp, args=[cset], maxiter=20000, disp=False)
    elif solver in 'cC':
        optAb = optimize.fmin_cobyla(objfn2, sp, ieqconsfn, args=[cset], iprint=0)
        optAb = optAb.ravel()
    tA, tb = splitAb(optAb, cset.nd)
    ts = -ones(tb.shape)
    optsol = ConSet(tA, ts, tb)
    if solver in 'bB':
        optcent = sum(optsol.vert, axis=0)/len(optsol.vert)
        itersol = ConSet(optsol.vert)
        while not itersol.allinside(cset)[0]:
            vi = (itersol.vert - optcent)*0.9999 + optcent
            itersol = ConSet(vi)
        optsol = itersol
    return optsol

def fitcube(cset, spset, solver):
    """
    Fit a rectangle (high/low limits on outputs) within an existing
    constraint set.
     cset - [ConSet] existing constraint set
    """
    #### State problem
    # nvar = cset.nd
    # Constraints are bounding
    # All vertices within constraint set
    #### Define parameters
    sp = spset.b.ravel() # starting point - b matrix to optimise
    #### Objective fn (SLSQP)
    def objfn(Ab, *args):
        """Volume objective function for SLSQP."""
        initcs = args[0]
        A = r_[eye(initcs.nd), -eye(initcs.nd)]
        b = array([Ab]).T
        cl = con2vert(A, b)[1]
        if cl:
            closed = 1
        else:
            closed = -1
        vol = tryvol(A, b, initcs)[0]
        return closed * -vol
    #### Constraints
    def ieqconsfn(Ab, *args):
        """Optimiser inequality constraint function."""
        initcs = args[0]
        A = r_[eye(initcs.nd), -eye(initcs.nd)]
        b = array([Ab]).T
        # get vertices for constraint set iteration
        V = tryvol(A, b, initcs)[1]
        iterset = ConSet(V)
        #constraint checking for vertices
        ineqs = iterset.allinside(initcs)[1]
        return array([-linalg.norm(ineqs)])
    
    #### Objective fn (FMIN)
    def objfn2(Ab, *args):
        """Volume objective function for fmin (Simplex)."""
        initcs = args[0]
        A = r_[eye(initcs.nd), -eye(initcs.nd)]
        b = array([Ab]).T
        vol, V = tryvol(A, b, initcs)
        Pv = 200.
        #Penalties
        # points outside of init space
        iterset = ConSet(V)
        outnorm = linalg.norm(iterset.allinside(initcs)[1])
        # open shape
        cl = con2vert(A, b)[1]
        if cl:
            closed = 1
        else:
            closed = -1
        vol = tryvol(A, b, initcs)[0]
        return (-vol*closed) + Pv*(outnorm**3)
    
    #### Maximise volume
    if solver in 'aA':
        optAb = optimize.fmin_slsqp(objfn, sp, f_ieqcons=ieqconsfn, args=[cset], 
                                iprint=0)
    elif solver in 'bB':
        optAb = optimize.fmin(objfn2, sp, args=[cset], maxiter=50000, disp=False)
    if solver in 'cC':
        optAb = optimize.fmin_cobyla(objfn, sp, ieqconsfn, args=[cset], 
                                iprint=0)
    tA = r_[eye(cset.nd), -eye(cset.nd)]
    tb = array([optAb]).T
    ts = -ones(tb.shape)
    optsol = ConSet(tA, ts, tb)
    if solver in 'bBcC':
        optcent = sum(optsol.vert, axis=0)/len(optsol.vert)
        itersol = ConSet(optsol.vert)
        while not itersol.allinside(cset)[0]:
            vi = (itersol.vert - optcent)*0.9999 + optcent
            itersol = ConSet(vi)
        optsol = itersol
    return optsol

def fitset(cset, *args):
    """
    Return a fitted constraint set within cset.
      cset - [ConSet] to fit within
      *args (in this order) - (stype)[string] 'a'rbitrary / 'r'ectangular
                            - (solver)[string] 'a' SLSQP / 'b' fmin / 'c' Cobyla
                            - (ncon)[integer] number of constraints to fit
                                cset.nd + 1 < ncon < cset.A.shape[0]
    """
    #Get args
    stype = args[0]
    solver = args[1]
    if stype in 'aA':
        ncon = args[2]
        nostarts = 1
    else:
        ncon = 0
        nostarts = 1
    #Starting point
    refvol = 0
    for k in range(nostarts):
        spshape = genstart(stype, cset, ncon)
        if spshape.vol() > refvol:
            finalsp = spshape
            refvol = finalsp.vol()
    #Fit set
    if stype in 'aA':
        optsol = fitshape(cset, finalsp, solver)
    elif stype in 'rR':
        optsol = fitcube(cset, finalsp, solver)
    #Return
    return optsol

def fitmaxbox(cset, sf):
    """
    Return high/low constraint set around the given constraint set.
      cset - [ConSet] to fit over
      sf - [integer] safety factor (fraction 0-1)
    """
    #Get args
    dev = sf*((cset.vert.max(0))-cset.cscent)
    maxvals = cset.vert.max(0) + dev#-cset.cscent) + cset.cscent
    minvals = cset.vert.min(0) - dev#-cset.cscent) + cset.cscent
    Abox = r_[eye(cset.nd), -eye(cset.nd)]
    bbox = c_[maxvals, -minvals].T
    sbox = -ones((Abox.shape[0], 1))
    return ConSet(Abox, sbox, bbox)