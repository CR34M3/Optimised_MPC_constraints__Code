#!/usr/bin/env python

# con2vert.py
# Python implementation of con2vert.m by Michael Kleder (July 2005),
#  available: http://www.mathworks.com/matlabcentral/fileexchange/7894-con2vert-constraints-to-vertices
# Converts sets of constraints to a list of vertices (of the feasible region)
#
# Author: Michael Kelder (Original)
#         Andre Campher (Python implementation)

# ===== Original MATLAB Code =====
#c = A\b;
#if ~all(A*c < b);
#    [c,f,ef] = fminsearch(@obj,c,'params',{A,b});
#    if ef ~= 1
#        error('Unable to locate a point within the interior of a feasible region.')
#    end
#end
#b = b - A*c;
#D = A ./ repmat(b,[1 size(A,2)]);
#[k,v2] = convhulln([D;zeros(1,size(D,2))]);
#[k,v1] = convhulln(D);
#if v2 > v1
#    error('Non-bounding constraints detected. (Consider box constraints on variables.)')
#end
#nr = unique(k(:));
#G  = zeros(size(k,1),size(D,2));
#for ix = 1:size(k,1)
#    F = D(k(ix,:),:);
#    G(ix,:)=F\ones(size(F,1),1);
#end
#V = G + repmat(c',[size(G,1),1]);
#[null,I]=unique(num2str(V,6),'rows');
#V=V(I,:);
#return
#function d = obj(c,params)
#A=params{1};
#b=params{2};
#d = A*c-b;
#k=(d>=-1e-15);
#d(k)=d(k)+1;
#d = max([0;d]);
#return

# Dependencies : - Scipy
#		 - Numpy
#		 - gendatafile
#		 - uniqmat

from scipy import *
from numpy import linalg, matlib
from gendatafile import *
from uniqmat import *
from os import remove
from string import join
from sys import exit
from itertools import groupby
import subprocess

#unit cube in 3D for testing
#A = mat('-1 0 0; 0 -1 0; 0 0 -1; 1 0 0; 0 1 0; 0 0 1')
#b = mat('0;0;0;1;1;1')
#cut cube in 3D for testing
#A = mat('-1 0 0;0 -1 0;0 0 -1;1 0 0;0 1 0; 0 0 1; 1 1 1')
#b = mat('0;0;0;2;2;2;5')
#2D example - 5 equations
A = mat('1 0;0 1;-1 0;0 -1;-1 -1')
b = mat('200;170;-100;-80;-200')

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
	print 'error : Non-bounding constraints detected. (Consider box constraints on variables.)'
	exit(1)
#== ==

qhullp = subprocess.Popen('qhull Ft < qhullin', shell=True, stdout=subprocess.PIPE) #calc vertices and facets
Vc = qhullp.communicate()[0] #qhull output to Vc
ks = Vc.split('\n')
fms = int(ks[1].split(' ')[1]) #get size of facet matrix
fmat = ks[-fms-1:-1]
fmat = mat(join(fmat,';')) #generate matrix
fmatn = fmat[:,0] #number of points on facets
fmatv = fmat[:,1:] #vertices on facets

G  = zeros((fmatv.shape[0],D.shape[1]));
for ix in range(0,fmatv.shape[0]):
	F = D[fmatv[ix,:],:].squeeze()
	G[ix,:] = linalg.lstsq(F,ones((F.shape[0],1)))[0].transpose()

V = G + matlib.repmat(c.transpose(),G.shape[0],1)

ux = uniqm(V,0.01)

print ux.round(3)
remove('qhullin')

#TODO =====
# error-checking
#	- fix volume check (for redundant constraints)
