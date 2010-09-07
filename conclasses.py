#!/usr/bin/env python
"""
Class definitions for constraint set
Author: Andre Campher
"""
# Dependencies: - convertfuns
#               - auxfuns
#               - SciPy

from auxfuns import qhull
from convertfuns import *
from scipy import *

class ConSet:
    def __init__(self,*inargs):
        if len(inargs) == 1:
            self.vert = inargs[0]
            self.A, self.s, self.b = vert2con(self.vert)
        elif len(inargs) == 3:
            self.A, self.s, self.b = inargs
            self.vert = con2vert(self.A,self.b)
        else:
            exit(1)  # TODO: Raise exception
        self.nd = self.A.shape[1]
        self.cons = hstack((self.A, self.s, self.b))  # (for historical reasons) TODO: remove
        self.vol = qhull(self.vert,"FA")
        
    def outconlin(self,model):
        """Convert constraints to output space using a linear model"""       
        # calc AOS (from G and AIS)
        outverttemp = empty([1,self.vert.shape[1]])
        for v in self.vert:
            x = model*v.transpose()
            outverttemp = vstack((outverttemp,x.transpose()))
        return vert2con(outverttemp[1:,:]) #remove first line of junk data from outverttemp and convert
    
    def intersect(self,conset2):
        """Determine intersection between current constraint set and another"""
        combA = vstack((self.A,conset2.A))
        combb = vstack((self.b,conset2.b))
        intcombvert = con2vert(combA,combb)
        return vert2con(intcombvert)