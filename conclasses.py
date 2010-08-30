#!/usr/bin/env python
"""
Class definitions for constraint set
Author: Andre Campher
"""
# Dependencies: - convertfuns
#               - auxfuns
#               - SciPy

from auxfuns import mat2ab, qhull
from convertfuns import *
from scipy import *

class conset:
    def __init__(self,conmatrix):
        self.A, self.s, self.b = mat2ab(conmatrix)
        self.nd = self.A.shape[1]
        self.vert = con2vert(self.A,self.b)
        self.cons = conmatrix
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
