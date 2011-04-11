"""Compute error in acoustics solution"""

import numpy as np
import matplotlib.pyplot as pl
from pyclaw.plotters.data import ClawPlotData
from sharpclawtest import compute_exact_solution
import setrun
from pyclaw.runclaw import runclaw
reload(setrun)

outdir='./_output'
frame=1
err=[]

#for mx in [200,400,800,1600,3200]:
for mx in [3200,6400]:
    rundata=setrun.setrun('sharpclaw')
    rundata.probdata.ic=9
    rundata.probdata.x0=-4.0
    rundata.probdata.a=4.0  # 1 for narrow pulse, 4 for wide pulse
    rundata.clawdata.mx=mx
    rundata.clawdata.nout=1
    rundata.clawdata.tfinal=8
    rundata.write()
    runclaw(xclawcmd='xsclaw',outdir=outdir)

    #Get the material parameters
    aux = np.loadtxt(outdir+'/fort.aux')
    rho = aux[:,0]; K   = aux[:,1]

    plotdata = ClawPlotData()
    plotdata.outdir=outdir

    #Read in the solution
    dat = plotdata.getframe(frame)
    eps = dat.q[:,0]
    p = -eps*K
    u = dat.q[:,1]

    exact = compute_exact_solution(outdir,frame)
    error = sum(abs(exact[1,:]-u))*dat.d[0]
    err.append(error)
    #    xc=dat.center[0]
    #    pl.clf()
    #    pl.plot(xc,exact[1,:],'-k',xc,u,'or')
    #    pl.draw()
    print error

conv=np.zeros(len(err))

for i in range(1,len(err)):
    conv[i]=-np.log2(err[i]/err[i-1])
print conv
