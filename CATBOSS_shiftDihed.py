#!/usr/bin/env python

import argparse
parser = argparse.ArgumentParser(description='This script takes a dihedral trajectory and shift the values to reflect the dihedral change evolution considering circularity.')
parser.add_argument('dihed', default='lG6.dihed.txt', help='dihedral traj file, e.g. lG6.dihedral.txt')
args = parser.parse_args()

import numpy as np

inputfile=args.dihed
outputfile=inputfile[:-4]+"_shifted.dat"
        
def DihedShift(before,after):
    '''
    Dihedral angle in degree, this function calculates the absoulte dihedral angle changes considering the circulation
    '''
    delta=after-before
    return after-360*int(round(delta/360))

def format3f(x):
    if type(x)==str:
        return x
    else:
        return str(round(x,3))
    
with open(inputfile) as f:
    with open(outputfile,'w') as fo:
	#for discard in range(0,17):
	#	n = f.readline()
        for i, line in enumerate(f):
            if i == 0:
                oldline=line.strip().split()
            newline=line.strip().split() 
            for j in range(1,len(newline)):
                current=float(newline[j])
                previous=float(oldline[j])
                newline[j]=DihedShift(previous,current)
            oldline=newline
            newline+='\n'
            #print (newline)
            fo.write(''.join(['%11s' % (i,) for i in list(map(format3f,newline))]))
           

