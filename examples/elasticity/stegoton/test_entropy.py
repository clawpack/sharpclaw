#Run a sequence of entropy tests
import os
import setrun
#from setrun import setrun
from entropy import entropy
import pylab as pl
from pyclaw.runclaw import runclaw
import numpy as np

reload(setrun)
grids = [1800,3600,7200,14400,28800]
lims=[1,2,3,4]
ntrp = np.zeros([len(lims),len(grids)])

for i,lim in enumerate(lims):
  for j,mx in enumerate(grids):
    print mx,lim

    rundata=setrun.setrun()

    rundata.probdata.K_B=4.0
    rundata.probdata.rho_B=4.0
    rundata.probdata.K_A=1.0
    rundata.probdata.rho_A=1.0

    rundata.clawdata.mthlim = [lim,lim]
    rundata.clawdata.mx = mx

    rundata.write()
    runclaw(outdir='./_output')

    ent=entropy(rundata.clawdata.nout)
    ntrp[i,j]=ent[-1]


linestyles=['-k','--k','-.k',':k']
limnames=['Minmod','Superbee','van Leer','MC']

pl.hold(False)
for i,lim in enumerate(lims):
  pl.loglog(np.array(grids)/150.,abs(ntrp[i,:]-1.),linestyles[i],linewidth=3)
  pl.hold(True)

pl.legend(limnames)

for i,lim in enumerate(lims):
  for j,mx in enumerate(grids):
    if ntrp[i,j]>1: pl.plot(mx/150.,abs(ntrp[i,j]-1.),'+k',markersize=10)

pl.axis([4,150,1.e-4,10])
pl.xlabel('Grid cells per layer',fontsize=20)
pl.ylabel('$\Delta \eta/\eta_0$',fontsize=20)
pl.title('Entropy Convergence',fontsize=20)
pl.hold(False)
