import numpy as np 
import os
traj  = 1000
wc    = 0.1/27.2114 
nfold = 20 
#os.system("rm -rf N-*")
for n in np.arange(nfold):
    os.system("rm -rf N-%s"%(n))
    os.mkdir("traj-%s"%(n)) 
    os.chdir("traj-%s"%(n)) 
    os.mkdir("log")
    os.system("cp ../vv.py ./")
    os.system("cp ../main.py ./")
    os.system("cp ../model.py ./")
    os.system("cp ../condor* ./")
    os.system("cp ../input.txt ./")
    #dat = np.array([traj, wc])

    os.system("condor_submit condor.sub")
 
    os.chdir("../") 