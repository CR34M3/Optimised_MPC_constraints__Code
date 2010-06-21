#!/usr/bin/env python
"""Functions to generate and manipulate data files needed for qhull process."""
# Author: Andre Campher

def genfile(V):
    """
    Generate data file for qhull process and save as plaintext file.
     V - [matrix] vertices to be put in data file
    """
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



if __name__ == "__main__":
    import doctest
    doctest.testfile("tests/gendatafiletests.txt")
