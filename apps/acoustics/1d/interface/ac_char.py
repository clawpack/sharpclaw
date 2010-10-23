#! /usr/bin/env python

#Solve 1D acoustics problems by characteristics

from numpy import *
from pylab import *
from actools import *
import time

#Space grid
xmin=-10.
xmax=10.
N=400
x=linspace(xmin,xmax,N)
dx=x[2]-x[1]

#Time grid
tmin=0.
tmax=10.1
Nt=40
t=linspace(tmin,tmax,Nt)

q=zeros([Nt,2,N],float) # q(i,1,:)= pressure; q(i,2,:)=velocity at time t[i]

#Set up material layers
#aux.interface=array([5])
#aux.c=array([1,0.5])
#aux.rho=array([1,4])
aux.interface=array([5,6])
aux.c=array([1,0.5,2])
aux.rho=array([1,4,3])
aux.z=aux.c*aux.rho
nint=size(aux.interface)
nlayer=size(aux.interface)+1
#Set left and right eigenvectors
aux.r1=zeros([nlayer,2])   
aux.r2=zeros([nlayer,2])
aux.l1=zeros([nlayer,2])
aux.l2=zeros([nlayer,2])
aux.int_r1=zeros([nlayer,2])
aux.int_r2=zeros([nlayer,2])
aux.int_l1=zeros([nlayer,2])
aux.int_l2=zeros([nlayer,2])
aux.CT_lr=zeros([nint])
aux.CR_lr=zeros([nint])
aux.CT_rl=zeros([nint])
aux.CR_rl=zeros([nint])
for ilayer in range(nlayer):
	aux.r1[ilayer,:]=array([-aux.z[ilayer],1])
	aux.r2[ilayer,:]=array([ aux.z[ilayer],1])
	aux.l1[ilayer,:]=1./(2*aux.z[ilayer])*array([-1,aux.z[ilayer]])
	aux.l2[ilayer,:]=1./(2*aux.z[ilayer])*array([ 1,aux.z[ilayer]])

for iint in range(nint):
	aux.int_r1[iint,:]=array([-aux.z[iint],1])
	aux.int_r2[iint,:]=array([ aux.z[iint+1],1])
	aux.int_l1[iint,:]=1./(aux.z[iint]+aux.z[iint+1])*array([-1,aux.z[iint+1]])
	aux.int_l2[iint,:]=1./(aux.z[ilayer]+aux.z[iint+1])*array([1,aux.z[iint]])

	#Transmission and reflection coefficients
	ziave=1./(aux.z[iint+1]+aux.z[iint])
	aux.CT_lr[iint]=2*aux.z[iint+1]*ziave
	aux.CR_lr[iint]=(aux.z[iint+1]-aux.z[iint])*ziave
	aux.CT_rl[iint]=2*aux.z[iint]*ziave
	aux.CR_rl[iint]=(aux.z[iint]-aux.z[iint+1])*ziave

ion()

#Begin time stepping
for i in range(Nt):
	for j in range(N):
		q[i,:,j]=qchar(x[j],t[i],aux)

	clf()
	plot(x,q[i,0,:],'b',x,q[i,1,:],'r')
	draw()
	time.sleep(0.1)
	print(t[i])
ioff()
