#Script to plot evolution of total entropy over time
#
#IMPORTANT: Assumes that output files contain velocity (u) and stress (sigma)
# and that sigma = exp(K*epsilon)
#
# In general, the entropy is given by
# eta = 0.5 * rho * u^2 + \int_0^epsilon sigma(s) ds.

import numpy as np
import pylab as pl
import time
from pyclaw.plotters.data import ClawPlotData

def entropy(nend,nstart=8,outdir='./_output',sigfun='exp'):
    aux = np.loadtxt('_output/fort.aux')
    rho = aux[:,0]
    K   = aux[:,1]

    plotdata = ClawPlotData()
    plotdata.outdir='./_output'
    ents=np.zeros([nend-nstart+1,1])

    for i in range(nstart,nend+1):
        dat = plotdata.getframe(i)
        u = dat.q[:,1]
        mom=rho*u**2/2.

        # Change the next line according to the functional form of sigma:
        # For exponential relation sigma(eps) = e^(K*eps) - 1:
        # MUST OUTPUT SIGMA!
        sigma = dat.q[:,0]
        sigint = (sigma-np.log(sigma+1.))/K 

        # For sigma(eps) = K*eps + 0.5*eps**2:
        # MUST OUTPUT EPS!
        #eps = dat.q[:,0]
        #sigint = K*eps/2. + eps**3/6.

        ent = (mom+sigint)
        #pl.clf()
        #pl.plot(nrg)
        #pl.axis([0,1800,0,0.015])
        #pl.draw()
        #time.sleep(2.1)

        ents[i-nstart]=np.sum(mom+sigint)*dat.grid.d[0]
        print ents[i-nstart]

    ents=ents/ents[0]

    pl.clf()
    pl.plot(ents)
    pl.draw()

    return ents
