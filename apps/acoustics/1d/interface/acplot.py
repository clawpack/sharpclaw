#!/usr/bin/env python

from numpy import *
from clawtools import *
from actools import *

ac=AcousticsProblem('setprob.data','claw1ez.data')

tmin=0.
tmax=30.1
Nt=40
Nx=817
t=linspace(tmin,tmax,Nt)
x=linspace(ac.xlower,ac.xupper,Nx)
q=zeros([2,size(x)])

for i in range(len(t)):
    for j in range(len(x)):
        q[:,j]=ac.qchar(x[j],t[i])

    ion()
    clf()
    plot(x,q[0,:],'k')
    draw()
    time.sleep(0.1)
    ioff

