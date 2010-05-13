#!/usr/bin/env python

# gendatafile.py
# Generate data file for qhull process
#    INPUT:  V - matrix
#    OUTPUT: qhullin - plaintext file with matrix (A) data  [file is saved]
#
# Author: Andre Campher

def genfile(V):
# generate data file for qhull
# format:
# n //dimensions (columns of V)
# p //number of points (rows of V)
# v //points e.g. 1 3 12 \n 2 4 5

    infile = open("qhullin",'w') #generate input file
    infile.write(str(V.shape[1])+"\n"+str(V.shape[0])+"\n")
    for row in range(0,V.shape[0]):
        stringmat = ""
        for col in range(0,V.shape[1]):
            stringmat = " ".join([stringmat,str(V[row,col])])
        infile.write(stringmat+"\n")
    infile.close()
