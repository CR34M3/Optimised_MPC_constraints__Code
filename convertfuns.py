#!/usr/bin/env python
"""Functions to calculate vertices from constraints and vice versa."""
from scipy import *
from numpy import linalg, matlib
import subprocess #to use qhull
from gendatafile import *
from auxfuns import *
from os import remove
from sys import exit

def vert2con(V):
    """
    Convert sets of vertices to a list of constraints (of the feasible region).
    Return A, b and s of the set;  Ax < b   (s is the sign-vector [to be used later])
    """
    # Dependencies: * qhull (libqhull5, qhull-bin)
    #               * scipy
    #               * gendatafile
    genfile(V)
# subprocess.call(["function","arguments"]) or subprocess.Popen('function expression', shell=True)
# run qhull with: qhull < data or cat data | qhull
    qhullp = subprocess.Popen('qhull n < qhullin', shell=True, stdout=subprocess.PIPE) #calc convex hull and get normals
    Vc = qhullp.communicate()[0] #qhull output to Vc
    ks = ';'.join(Vc.split('\n')[2:]) #remove leading dimension output
    k = mat(ks[:-1]) #convert to martrix with vertices
# k is a (n+1)x(p) matrix in the form [A b] (from qhull doc: Ax < -b is satisfied), thus;
    A = k[:,:-1]
    b = -k[:,-1]
    s = -ones([k.shape[0],1])

    remove('qhullin')
    return A,s,b

def con2vert(A,b):
    """Convert sets of constraints to a list of vertices (of the feasible region)."""
    # Python implementation of con2vert.m by Michael Kleder (July 2005),
    #  available: http://www.mathworks.com/matlabcentral/fileexchange/7894-con2vert-constraints-to-vertices
    # Author: Michael Kelder (Original)
    #         Andre Campher (Python implementation)
    #
    # Dependencies : - Scipy
    #         - Numpy
    #         - gendatafile
    #         - uniqmat
    c = linalg.lstsq(A,b)[0]
    b = b-A*c
    D = A / matlib.repmat(b,1,A.shape[1])
    Dtest = vstack((D,zeros([1,D.shape[1]])))

#== Volume error check ==
    genfile(Dtest)
    qhullp = subprocess.Popen('qhull FA < qhullin', shell=True, stdout=subprocess.PIPE) #calc summary and volume
    Vc = qhullp.communicate()[0] #qhull output to Vc
    ks = Vc.split('\n')[-3]
    VolDt = float(ks.split(' ')[-1]) #get volume of D-hull

    genfile(D)
    qhullp = subprocess.Popen('qhull FA < qhullin', shell=True, stdout=subprocess.PIPE) #calc summary and volume
    Vc = qhullp.communicate()[0] #qhull output to Vc
    ks = Vc.split('\n')[-3]
    VolD = float(ks.split(' ')[-1]) #get volume of D-hull

    if VolDt > VolD:
        print 'error : Non-bounding constraints detected (consider box constraints on variables). Exiting...'
        exit(1)
#== ==

    qhullp = subprocess.Popen('qhull Ft < qhullin', shell=True, stdout=subprocess.PIPE) #calc vertices and facets
    Vc = qhullp.communicate()[0] #qhull output to Vc
    ks = Vc.split('\n')
    fms = int(ks[1].split(' ')[1]) #get size of facet matrix
    fmat = ks[-fms-1:-1]
    fmat = mat(';'.join(fmat)) #generate matrix
    fmatn = fmat[:,0] #number of points on facets
    fmatv = fmat[:,1:] #vertices on facets

    G  = zeros((fmatv.shape[0],D.shape[1]));
    for ix in range(0,fmatv.shape[0]):
        F = D[fmatv[ix,:],:].squeeze()
        G[ix,:] = linalg.lstsq(F,ones((F.shape[0],1)))[0].transpose()

    V = G + matlib.repmat(c.transpose(),G.shape[0],1)
    ux = uniqm(V,0.01)

    remove('qhullin')
    return ux

#TODO con2vert =====
# error-checking
#    - fix volume check (for redundant constraints)
