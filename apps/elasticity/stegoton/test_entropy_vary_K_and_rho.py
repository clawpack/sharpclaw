#Run a sequence of entropy tests
import os
import setrun
#from setrun import setrun
from entropy import entropy
import pylab as pl
from pyclaw.runclaw import runclaw
import numpy as np

reload(setrun)
#grids = [450,900,1800,3600]
#lims=[1,2,3,4]
Ks=np.array([1,1.5,2,3,4,6,8,16,32,64])
rhos=np.array([1,1.5,2,3,4,6,8,16,32,64])
ntrp = np.zeros([len(Ks),len(rhos)])

for i,K in enumerate(Ks):
  for j,rho in enumerate(rhos):
    rundata=setrun.setrun()
    rundata.probdata.K_B=K
    rundata.probdata.rho_B=rho
    rundata.probdata.K_A=1.
    rundata.probdata.rho_A=1.
    print K,rho

    #Rescale time and IC:
    cmean = 2*np.sqrt(K/((K+1.)*(rho+1.)))
    Kmean = 2*K/(K+1.)
    rhomean = 0.5*(rho+1.)
    rundata.probdata.a1 = rundata.probdata.a1/np.sqrt(rhomean)
    rundata.clawdata.tfinal = 500/(cmean/0.8)
    rundata.write()
    runclaw(outdir='./_output')

    ent=entropy(rundata.clawdata.nout)
    ntrp[i,j]=ent[-1]
    pl.pcolor(np.log2(Ks),np.log2(rhos),np.log10(abs(1.-ntrp)),cmap='gray')
    pl.colorbar(); pl.axis('image')
    pl.colorbar()
    pl.draw()

#t=10*np.arange(8,51)
#pl.hold(False)

#for nt in ntrp:
#    pl.plot(t,nt)
#    pl.hold(True)

#pl.hold(False)
#pl.legend(grids)

pl.pcolor(np.log2(Ks),np.log2(rhos),np.log10(abs(1.-ntrp)),cmap='gray')
pl.colorbar(); pl.axis('image')
pl.xlabel('$\log_2(K_B)$',fontsize=20)
pl.ylabel('$\log_2(p_B)$',fontsize=20)
pl.title('$\log_{10}|\eta-\eta_0|$',fontsize=25)
pl.draw()
