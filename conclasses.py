#!/usr/bin/env python
"""
Class definitions for constraint set
Author: Andre Campher
"""
# Dependencies: - convertfuns
#               - auxfuns
#               - SciPy

from auxfuns import qhull
from convertfuns import vert2con, con2vert
from scipy import hstack, empty, vstack, dot

class ConSet:
    """
    Class for constraint sets. Generated by either a constraint set [A s b] or
    by a set of vertices [v].
    """
    def __init__(self, *inargs):
        if len(inargs) == 1:
            self.vert = inargs[0]
            self.A, self.s, self.b = vert2con(self.vert)
        elif len(inargs) == 3:
            self.A, self.s, self.b = inargs
            self.vert = con2vert(self.A, self.b)
        else:
            exit(1)  # TODO: Raise exception
        self.nd = self.A.shape[1]
        self.cons = hstack((self.A, self.s, self.b))  # TODO: remove this
        self.vol = qhull(self.vert,"FA")
        
    def outconlin(self, model):
        """Convert constraints to output space using a linear model"""       
        # calc AOS (from G and AIS)
        outverttemp = empty([1, self.vert.shape[1]])
        for v in self.vert:
            x = model * v.transpose()
            outverttemp = vstack((outverttemp, x.transpose()))
        #remove first line of junk data from outverttemp and convert
        return vert2con(outverttemp[1:, :]) 
    
    def intersect(self, conset2):
        """Determine intersection between current constraint set and another"""
        combA = vstack((self.A, conset2.A))
        combb = vstack((self.b, conset2.b))
        intcombvert = con2vert(combA, combb)
        return vert2con(intcombvert)
    
    def allinside(self, conset2):
        """
        Determine if all vertices of self is within conset2. allvinside merely
        returns True/False whereas insidenorm returns a measure of 'inside-ness'
        better suited for optimisers.
        """
        # Inside check
        print dot(conset2.A, self.vert.T)
        allvinside = True
        # Inside norm
        insidenorm = 1
        return allvinside, insidenorm