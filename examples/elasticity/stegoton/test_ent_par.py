#Run a sequence of entropy tests
import os
import setrun
#from setrun import setrun
from entropy import entropy
import pylab as pl
from pyclaw.runclaw import runclaw
import numpy as np
from multiprocessing import pool

reload(setrun)
#grids = [450,900,1800,3600]
#lims=[1,2,3,4]
Ks=np.array([1,2,4,8,16,32,64])
rhos=np.array([1,2,4,8,16,32,64])
ntrp = np.zeros([len(Ks),len(rhos)])

for i,K in enumerate(Ks):
  for j,rho in enumerate(rhos):
    ntrp[i,j]=ent[-1]
    pl.pcolor(Ks,rhos,ntrp)
    pl.colorbar()
    pl.draw()

#t=10*np.arange(8,51)
#pl.hold(False)

#for nt in ntrp:
#    pl.plot(t,nt)
#    pl.hold(True)

#pl.hold(False)
#pl.legend(grids)

pl.pcolor(Ks,rhos,ntrp)
pl.colorbar()
pl.draw()

def getent(K,rho):
  rundata=setrun.setrun()
  rundata.probdata.K_B=K
  rundata.probdata.rho_B=rho
  print K,rho
  rundata.write()
  runclaw(outdir='./_output')

  ent=entropy(rundata.clawdata.nout)
  return ent[-1]
