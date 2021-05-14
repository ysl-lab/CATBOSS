import argparse
parser = argparse.ArgumentParser(description='This script takes a dihedral trajectory and detects change points using SIMPLE (simultaneous Penalized Likelihood Estimation, see Fan et al. P. Natl. Acad. Sci, 2015, 112, 7454-7459). Two parameters alpha and lambda are controlling the extent of simultaneous changes and total number of changes detected, respectively. alpha -> 1 means more simultaneous changes (0<alpha<1), and smaller lambda gives more changes.')
parser.add_argument('shifteddihed', default='shifted_dihedral.dat', help='input shifted dihedral file')
parser.add_argument('--alpha', type=float, default=0.7, help='extent of simultaneous changes, 0.7 by default suggested by the author if no prior ')
parser.add_argument('--lam', type=float, default=10, help='sensitivity of detecting changes, 10 by default')
args = parser.parse_args()

import numpy as np
from SIMPLEchangepoint import ComputeChanges
import collections

inputfile=args.shifteddihed
lam=args.lam
alpha=args.alpha

outputfile=inputfile[:-4]+".lam"+str(lam)+"alpha"+str(alpha)+".transitionProba.dat"
outputfile2=inputfile[:-4]+".lam"+str(lam)+"alpha"+str(alpha)+".transitionSummary.dat"

alldata=np.loadtxt(inputfile).T

time=alldata[0]
data=alldata[1:]
CPDresults = ComputeChanges(data,lam,alpha,lam_min=0,parallel=False)

def changeORnot(con_set,size):
    x=[0]*size
    for i in con_set:
        x[i] = 1
    return '   '.join(map(str,x))
    
od = collections.OrderedDict(sorted(CPDresults.items()))
with open(outputfile,'w') as fo:
    for t in range(len(time)):
        if t not in od.keys():
            fo.write(str(time[t])+'   '+'   '.join(map(str,[0]*len(data)))+'\n')
        else:
            fo.write(str(time[t])+'   '+changeORnot(od[t],len(data))+'\n')

def strplus1(x):
    return str(x+1)
with open(outputfile2,'w') as fo2:
    for k, v in od.iteritems():
        fo2.write('{:10.1f}  {:5d}  {:s}\n'.format(time[k],len(v),','.join(map(strplus1,v))))
