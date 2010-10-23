#Run a sequence of time-reversal tests
import os
import setrun
#from setrun import setrun
from entropy import entropy
import pylab as pl
from pyclaw.runclaw import runclaw

def test_tr(grids=[1800,3600,7200],frame=2):
    reload(setrun)
    errs = []

    for i,mx in enumerate(grids):
        print i,mx
        rundata=setrun.setrun()

        rundata.probdata.trtime=600.
        rundata.probdata.t1= 10.
        rundata.probdata.a1= 0.1
        rundata.probdata.tw1=10.

        rundata.clawdata.mx=mx
        rundata.clawdata.tfinal=1200.
        rundata.clawdata.nout=60
        rundata.write()
        runclaw(xclawcmd='xsclaw',outdir='./_output')

        exact,approx,xc = timereverse(frame=frame)

        Linf_err = max(abs(approx-exact))
        print mx, Linf_err

        errs.append(Linf_err)

    return grids,errs

def timereverse(frame=2,nframes=60):
    from pyclaw.plotters.data import ClawPlotData

    plotdata = ClawPlotData()
    plotdata.outdir='./_output'

    #"Exact" solution
    dat0 = plotdata.getframe(frame)
    u0 = dat0.q[:,0]

    #Time-reversed solution
    dat1 = plotdata.getframe(nframes-frame)
    u1 = dat1.q[:,0]

    xc = dat0.c_center

    return u0,u1,xc[0]

def plot_tr(frame=2,nframes=60):

    exact,approx,xc = timereverse(frame,nframes)

    pl.plot(xc,exact,'-k',linewidth=4)
    pl.hold(True)
    N=len(xc); skip=10
    sub=range(0,N,10)
    pl.plot(xc[sub],approx[sub],'b--',linewidth=4)
    pl.hold(False)
    pl.axis([10,40,-0.02,0.7])
    pl.xlabel('x',fontsize=20)
    pl.ylabel('$\sigma$',fontsize=20)
