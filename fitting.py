#!/usr/bin/env python
"""Functions to optimally fit a 'shape' into another."""
from scipy import mat,optimize,vstack,eye,hstack,ones,tile,sqrt,power,real
import numpy

def tryvol(A,b,cs):
    """Try to determine volume of feasible region, otherwise impose box constraint."""
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

def splitab(inAb,nd):
    """Split input Ab matrix into separate A and b."""
    nd1 = nd+1
    #split Ab into A and b (assuming all -1 for s)
    Ab = mat(inAb)
    tempAb = vstack([Ab[:,p*nd1:p*nd1+nd1] for p in range(0,Ab.shape[1]/nd1)])
    return tempAb[:,:-1], tempAb[:,-1] #return A and b

def genstart(cs,ncon):
    """Generate a starting shape (for optimisation) with given number of faces."""
    #1. Determine centre of initial region, cscent
    c0 = ones((1,cs.shape[1]))
    def cntrobjfn(cntr):
        cmat = tile(cntr,(cs.shape[0],1))
        dist = sqrt(power((cs-cmat),2)*mat(ones((cmat.shape[1],1)))) #return distance between vertices and centre
        return real(sum(power(dist,2)))
    cscent = optimize.fmin(cntrobjfn,c0)
    #2. Generate n-sphere, radius r and centre cscent, with ncon number of points on surface
    def nsphere(ptsn,cntr,rad):
        cmatn = tile(cntr[:,:-1],(ptsn.shape[0],1))
        ptsn1 = sqrt((rad**2)-power(ptsn-cmatn,2)*mat(ones((cmatn.shape[1],1)))) + cntr[:,-1]
        return real(ptsn1)
    #3. Equally space point on sphere's surface (maximise distances to all other points)
    #def sphereobjfn(pts):
    #output as fraction of radius
    #4. Generate tangent planes on sphere at points
    #5. Convert tangent planes to inequalities and generate feasible region (always closed?)
    #6. Check if vertices of new feasible region is within initial shape
    #7. Optimise sphere-radius, r, to have all points within initial shape
    #return cscent

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
    from conclasses import conset
    v = mat('1 0;0 1;0 0;1 1')
    print genstart(v,2)