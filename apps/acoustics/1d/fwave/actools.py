from clawtools import *
import pylab as pl
from pylab import *


#======================
class AcousticsProblem:
#======================
    """
    Gives exact solution for 1D acoustics problem with piecewise constant
    medium with arbitrary number of interfaces
    """

    #==============================
    def __init__(self,*datafilelist):
    #==============================
        """
        Set up the medium, eigenvectors, and transmission/reflection 
	coefficients.  Typically, datafilelist='setprob.data','claw1ez.data'.
        """
        ProblemData = ClawData()
	ProblemData.datafiles=[]
	for datafile in datafilelist:
	    ProblemData.read(datafile)
        try:
	    #Is this really necessary, or can we just use ProblemData directly?
	    #say self=ProblemData?
            #Locations of material interfaces:
            self.interface=squeeze(pl.array([ProblemData.interface]))
            #Densities:
            self.rho=squeeze(pl.array([ProblemData.rho]))
            #Sound speeds:
            self.c=squeeze(pl.array([ProblemData.c]))
            self.ic=ProblemData.ic
            self.bcLeft=ProblemData.mthbc[0]
            self.bcRight=ProblemData.mthbc[1]
            self.xlower=ProblemData.xlower
            self.xupper=ProblemData.xupper
            self.a = ProblemData.a
            self.x0 = ProblemData.x0
        except:
            self.interface=squeeze(pl.array([ProblemData.interface]))
            #Densities:
            self.rho=squeeze(pl.array([ProblemData.rho]))
            #Sound speeds:
            self.c=squeeze(pl.array([ProblemData.c]))
            self.ic=ProblemData.ic
            self.bcLeft=ProblemData.mthbc_xlower
            self.bcRight=ProblemData.mthbc_xupper
            self.xlower=ProblemData.xlower
            self.xupper=ProblemData.xupper
            self.a = ProblemData.a
            self.x0 = ProblemData.x0
#        except:
#            print "ERROR: ProblemData lacks required fields"
#            print "Make sure you read in both setprob.data and claw1ez.data"
#            stop
        try: nint=len(self.interface)
        except: 
            nint=1
            self.interface=pl.array([self.interface])
        nlayer=nint+1
        #Impedances:
        self.z=self.c*self.rho
        #Set left and right eigenvectors
        self.r1=pl.zeros([nlayer,2])   
        self.r2=pl.zeros([nlayer,2])
        self.l1=pl.zeros([nlayer,2])
        self.l2=pl.zeros([nlayer,2])
        self.CT_lr=pl.zeros([nint+(self.bcLeft==2)])
        self.CR_lr=pl.zeros([nint+(self.bcLeft==2)])
        self.CT_rl=pl.zeros([nint+(self.bcLeft==2)])
        self.CR_rl=pl.zeros([nint+(self.bcLeft==2)])
        for ilayer in range(nlayer):
            #Right and left eigenvectors
            self.r1[ilayer,:]=pl.array([-self.z[ilayer],1])
            self.r2[ilayer,:]=pl.array([ self.z[ilayer],1])
            self.l1[ilayer,:]=1./(2*self.z[ilayer])*pl.array([-1,self.z[ilayer]])
            self.l2[ilayer,:]=1./(2*self.z[ilayer])*pl.array([ 1,self.z[ilayer]])
        for iint in range(nint):
            ziave=1./(self.z[iint+1]+self.z[iint])
            #Transmission and reflection coefficients
            self.CT_lr[iint]=2*self.z[iint+1]*ziave
            self.CR_lr[iint]=(self.z[iint+1]-self.z[iint])*ziave
            self.CT_rl[iint]=2*self.z[iint]*ziave
            self.CR_rl[iint]=(self.z[iint]-self.z[iint+1])*ziave
	if self.bcLeft==2: #Add the periodic boundary interface
            self.interface=insert(self.interface,nint,self.xupper)
            nint+=1
            #Transmission and reflection coefficients for periodic boundary
            ziave=1./(self.z[0]+self.z[nint-1])
            self.CT_lr[nint-1]=2*self.z[0]*ziave
            self.CR_lr[nint-1]=(self.z[0]-self.z[nint-1])*ziave
            self.CT_rl[nint-1]=2*self.z[nint-1]*ziave
            self.CR_rl[nint-1]=(self.z[nint-1]-self.z[0])*ziave
        self.nint=nint
        self.nlayer=nlayer


    #==============================
    def uvec(self,xarr,t):
    #==============================
        """Find solution to 1D acoustics by characteristics"""
        try: xarr.__iter__()
        except: xarr = [xarr]
        q=pl.zeros([pl.size(xarr)])
        i=0
        for x in xarr:
            imat=self.getmat(x)
            p1=self.w1x(x,t)
            p2=self.w2x(x,t)
            q[i] = (p2-p1)/self.z[imat] #Velocity
            i+=1
        return q


    #==============================
    def pvec(self,xarr,t):
    #==============================
        """Find solution to 1D acoustics by characteristics"""
        try: xarr.__iter__()
        except: xarr = [xarr]
        q=pl.zeros([pl.size(xarr)])
        i=0
        for x in xarr:
            imat=self.getmat(x)
            p1=self.w1x(x,t)
            p2=self.w2x(x,t)
            q[i] = p1+p2 #Pressure
            i+=1
        return q


    #==============================
    def qchar(self,x,t):
    #==============================
        """Find solution to 1D acoustics by characteristics"""
        imat=self.getmat(x)
        q=pl.array([0.,0.])
        p1=self.w1x(x,t)
        p2=self.w2x(x,t)
        q[0] = p1+p2 #Pressure
        q[1] = (p2-p1)/self.z[imat] #Velocity
        return q


    #==============================
    def w1x(self,x,t):
    #==============================
        """Evaluate 1-wave component at (x,t)"""
        #Returns the 1-wave component magnitude
        #First determine if there is an interface involved
        imat=self.getmat(x)
        iint=self.whichint(x,t,1) #Reverse 1-characteristics go to right
        if iint==-100: #No interface involved
            x0=x+self.c[imat]*t
            w1val=-pl.dot(self.l1[imat,:],self.q0(x0))*self.z[imat]
        else: #Characteristic crosses interface iint
            tnew=t+(x-self.interface[iint])/self.c[imat]
            w1val=self.w1int(iint,tnew)*self.CT_rl[iint] + self.w2int(iint,tnew)*self.CR_lr[iint]
        return w1val

    #==============================
    def w1int(self,iint,t):
    #==============================
        """Evaluate 1-wave component at (self.interface[iint],t)"""
        #Returns the 1-wave component magnitude
        #First determine if there is an interface involved
        if iint==self.nint-1:
            iint=-1
        x0=self.interface[iint]+self.c[iint+1]*t #Characteristic traceback
        if iint==-1: #At rightmost interface
            if self.bcLeft==2 and x0>self.xupper and \
              (mod(x0-self.xlower,self.xupper-self.xlower)+self.xlower)>self.interface[0]:#Crosses interface 0
                tnew=t-(self.interface[0]-self.interface[iint]+self.xupper-self.xlower)/self.c[0]
                w1val=self.w1int(0,tnew)*self.CT_rl[0] + self.w2int(0,tnew)*self.CR_lr[0]
            else: #Not periodic BCs or doesn't reach interface 0
                w1val=-pl.dot(self.l1[iint+1,:],self.q0(x0))*self.z[iint+1]
        elif x0<self.interface[iint+1]: #Doesn't reach next interface
            w1val=-pl.dot(self.l1[iint+1,:],self.q0(x0))*self.z[iint+1]
        else: #Characteristic crosses next interface to right
            tnew=t-(self.interface[iint+1]-self.interface[iint])/self.c[iint+1]
            w1val=self.w1int(iint+1,tnew)*self.CT_rl[iint+1] + self.w2int(iint+1,tnew)*self.CR_lr[iint+1]
        return w1val


    #==============================
    def w2x(self,x,t):
    #==============================
        """Evaluate 2-wave component at (x,t)"""
        #Returns the 2-wave component magnitude
        #First determine if there is an interface involved
        imat=self.getmat(x)
        iint=self.whichint(x,t,-1) #Reverse 2-characteristics go to left
        if iint==-100: #No interface involved
            x0=x-self.c[imat]*t
            w2val=pl.dot(self.l2[imat,:],self.q0(x0))*self.z[imat]
        else: #Characteristic crosses interface iint
            if iint==-1:
                tnew=t-(x-self.xlower)/self.c[imat]
	    else:
                tnew=t-(x-self.interface[iint])/self.c[imat]
            w2val=self.w1int(iint,tnew)*self.CR_rl[iint] + self.w2int(iint,tnew)*self.CT_lr[iint]
        return w2val

    #==============================
    def w2int(self,iint,t):
    #==============================
        """Evaluates 2-wave component at (self.interface[iint],t)"""
        #Returns the 2-wave component magnitude
        #First determine if there is an interface involved
        if iint<0: iint+=self.nint
        x0=self.interface[iint]-self.c[iint]*t #Characteristic traceback
        if iint==0: #At leftmost interface
            if self.bcLeft==2 and x0<self.xlower and \
              (mod(x0-self.xlower,self.xupper-self.xlower)+self.xlower)<self.interface[-1]:#Crosses periodic boundary
                tnew=t-(self.interface[iint]-self.interface[iint-1]+self.xupper-self.xlower)/self.c[iint]
                w2val=self.w1int(iint-1,tnew)*self.CR_rl[iint-1] + self.w2int(iint-1,tnew)*self.CT_lr[iint-1]
            else: #Not periodic BCs or doesn't reach rightmost interface
                w2val=pl.dot(self.l2[iint,:],self.q0(x0))*self.z[iint]
        elif x0>self.interface[iint-1]: #Doesn't reach next interface
            w2val=pl.dot(self.l2[iint,:],self.q0(x0))*self.z[iint]
        else: #Characteristic crosses interface to left (iint-1)
            tnew=t-(self.interface[iint]-self.interface[iint-1])/self.c[iint]
            w2val=self.w1int(iint-1,tnew)*self.CR_rl[iint-1] + self.w2int(iint-1,tnew)*self.CT_lr[iint-1]
        return w2val


    #==============================
    def getmat(self,x):
    #==============================
        """Determine index of material that x resides in"""
        if self.bcLeft==2: #Periodic BCs
            x=mod(x-self.xlower,self.xupper-self.xlower)+self.xlower
        try:
            thismat=max(pl.find(x-self.interface>0))+1
        except:
            #x lies in leftmost material
            thismat=0
        return thismat


    #==============================
    def whichint(self,x,t,dir):
    #==============================
        """Determine the first interface (if any) that characteristics from
           x cross"""
        mat=self.getmat(x)
        c=self.c[mat]*dir
        x0=x+c*t
        if self.bcLeft<>2 or not ((mat==0 and dir==-1) or (mat==self.nint-1 and dir==1)):
            if mat>0:
                if (x0-self.interface[mat-1])<0:
                    #Characteristic crosses interface to left
                    iint=mat-1
                    return iint
            if mat<self.nint:
                if (x0-self.interface[mat])>0:
                    #Characteristic crosses interface to right
                    iint=mat
                    return iint
        elif mat==0:  #Periodic BCs and leftmost material and dir=-1
            if x0<self.xlower:
                return -1
        else: #Periodic BCs and rightmost material and dir=1
            if x0>self.xupper:
                return -1

        return -100  #if we get to here, characteristic crosses no interfaces


    #==============================
    def q0(self,x):
    #==============================
        """Evaluate initial condition"""

        if self.bcLeft==2: #Periodic BCs
            x=pl.mod(x-self.xlower,self.xupper-self.xlower)+self.xlower

        # Eventually: a case statement for different values of self.ic

        #C5 Newton polynomial:
        x0=self.x0
        a=self.a
        q=pl.array([0.,0.])
        mat=self.getmat(x)
        z=self.z[mat]
        if abs(x-x0)<=a:
            q[0]=(x-x0-a)**6 * (x-x0+a)**6/a**12
            q[1]=q[0]/z #Purely right-going since Z=1

        return q

#==============================
#End of acMedium class
#==============================


