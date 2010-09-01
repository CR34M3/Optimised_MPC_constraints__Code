#!/usr/bin/env python
"""Functions to optimally fit a 'shape' into another."""
from scipy import mat,optimize,vstack,eye,hstack
import numpy

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