import matplotlib as mpl
from pylab import *
from clawtools import *
import actools as ac

ap=ac.AcousticsProblem('setprob.data','sharpclaw.data')
print ap.interface,ap.rho,ap.c
add=0
sep=0.5

fig_size = [3,3]
params = {'axes.labelsize': 10,
          'text.fontsize':  30,
          'title.fontsize':  30,
          'font.size':  20,
          'legend.fontsize': 10,
          'xtick.labelsize': 8,
          'ytick.labelsize': 8,
          'text.usetex': False,
          'backend': 'ps',
          'figure.figsize': fig_size}

mpl.rcParams.update(params)
ion()
for t in linspace(3,5,30):
    xc = linspace(-1,1,100)

    qq=zeros([2,size(xc)])
    for ii in range(size(xc)):
        qq[:,ii] = ap.qchar(xc[ii],t)
    subplot(1,2,1)
    plot(xc,qq[0,:]+add,'k')
    title('Pressure at t='+str(t))
    hold(True)
    subplot(1,2,2)
    plot(xc,qq[1,:]+add,'k')
    title('Velocity at t='+str(t))
    hold(True)
    draw()
    add+=sep

add-=sep
plot([0.,0.],[-sep/2.,sep/2.],'-k',lw=2.)
axis([-1.,1.,-sep/2.,add+qq.max()])
title('Pressure')
axis('off')
hold(False)
subplot(1,2,1)
plot([0.,0.],[-sep/2.,sep/2.],'-k',lw=2.)
axis([-1.,1.,-sep/2.,add+qq.max()],'-r')
title('Velocity')
axis('off')
hold(False)
draw()
ioff()

mpl.rcParams.update(params)
savefig('fig1.eps')
