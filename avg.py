import numpy as np
import os
import sys
nfold = 19
fold  = sys.argv[1]
n     = len(np.loadtxt(fold+"/f0/fs.txt")[:,0])
col   = len(np.loadtxt(fold+"/f0/fs.txt")[0,:])
dat   = np.zeros((n,col)) 
for i in range(nfold):
    fl = np.loadtxt(fold+"/f%s/fs.txt"%(i))
    for c in range(col):
        dat[:,c] += fl[:,c]/nfold

np.savetxt("%s.txt"%(fold),dat)
    


