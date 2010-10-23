c
c
c =========================================================
       subroutine qinit(maxmx,meqn,mbc,mx,xlower,dx,q,maux,aux)
c =========================================================
c
c     # Set initial conditions for q.
c     # Pulse in pressure, zero velocity
c
c
      implicit double precision (a-h,o-z)
      dimension q(1-mbc:maxmx+mbc, meqn)
      dimension aux(1-mbc:maxmx+mbc, *)
      common /cqinit/ a,x0,ic
c
c
        write(*,*) a,ic
      beta=0.1d0
      do 150 i=1,mx
         xcell = xlower + (i-0.5d0)*dx

         go to (10,20,30,40,50,60,70,80,90,100) ic


   10    continue
c        # half ellipse:
         if (xcell.gt.-4d0 .and. xcell.lt.-2d0) then
             q(i,1) = dsqrt(1.d0 - (xcell+3.d0)**2)
            else
             q(i,1) = 0.d0
            endif
         q(i,2) = q(i,1)/aux(i,1)
         go to 150

   20    continue
c        # single discontinuity:
         if (xcell .lt. -2.d0) then
             q(i,1) = 1.d0
            else
             q(i,1) = 0.d0
            endif
         q(i,2) = q(i,1)/aux(i,1)
         go to 150

   30    continue
c        # Gaussian and square pulse:
c        q(i,1) = dexp(-beta*(xcell+2.0d0)**2)  
c        if (dabs(q(i,1)) .lt. 1d-30) q(i,1) = 0.d0
         q(i,1) = 0.d0
         if (xcell.gt.-5.d0 .and. xcell.lt.-3.d0) then
            q(i,1) = q(i,1) + 1.0d0
            endif
         q(i,2) = q(i,1)/aux(i,1)
         go to 150

   40    continue
c        # Gaussian only:
c         q(i,1) = dexp(-beta*(xcell+2.0d0)**2)  

c        #Initialization with exact cell averages
c         xl=xcell-0.5d0*dx
c         xr=xcell+0.5d0*dx
c         sqrtbeta=sqrt(beta)
c         q(i,1) = 0.5d0*sqrt(3.1415926d0/beta)/dx*
c     2   (erf(sqrt(beta)*(xr+2.0d0)) - erf(sqrt(beta)*(xl+2.0d0)))

c        #Initialization with gauss quadrature
         xl=xcell-dx/sqrt(12.d0);
         xr=xcell+dx/sqrt(12.d0);
         q(i,1) = 0.5d0*(dexp(-beta*(xl-x0)**2) 
     2                     + dexp(-beta*(xr-x0)**2))
         q(i,2) = q(i,1)/aux(i,1)
         go to 150

   50    continue
c        # Wave packet:
         q(i,1) = dexp(-beta*(xcell-x0)**2)*dcos(20.d0/beta*xcell)
         q(i,2) = q(i,1)/aux(i,1)
c        #Initialization with exact cell averages
c         q(i,1) = dexp(-beta*(xcell+2.0d0)**2)*dcos(20.d0/beta*xcell)
c         q(i,2) = q(i,1)
         go to 150

   60    continue
c        # C3 Newton polynomial
c         pi = 3.14159265358979323846d0
         x0 = -4.d0
         xrelr=xcell+0.5d0*dx-x0
         xrell=xcell-0.5d0*dx-x0
c        #Point values would be given by:
c        # q(i,1) = (xcell-x0-a)^4*(xcell-x0+a)^4/a**8.d0
c        #Initialization with exact cell averages:
         q(i,1) = (xrelr**9.d0/(9.d0*a**8.d0)-4.d0*xrelr**7.d0
     2             /(7.d0*a**6.d0)
     2             + 6.d0*xrelr**5.d0/(5.d0*a**4.d0)
     2             - 4.d0*xrelr**3.d0/(3.d0*a*a) + xrelr
     2           - (xrell**9.d0/(9.d0*a**8.d0)-4.d0*xrell**7.d0
     2             /(7.d0*a**6.d0)
     2             + 6.d0*xrell**5.d0/(5.d0*a**4.d0)
     2             - 4.d0*xrell**3.d0/(3.d0*a*a) + xrell))/dx
c        # Ought to make the zeroing more accurate
         if (dabs(xcell-x0) .gt. a) q(i,1) = 0.d0
         q(i,2) = q(i,1)/aux(i,1)
         go to 150

   70    continue
c        # C3 Newton polynomial modulated by a cosine
         x0 = -4.d0
         xr=xcell+0.5d0*dx-x0
         xl=xcell-0.5d0*dx-x0
c        #Point values would be given by:
c        # q(i,1) = cos((x-x0)/beta)*(xcell-x0-a)^4*(xcell-x0+a)^4/a**8.d0
c        #Initialization with exact cell averages
         yy=beta/a
         a1=-yy**2.d0*(8.d0+144.d0*yy**2.d0+2880.d0*yy**4.d0
     2        + 40320.d0*yy**6.d0)
         a3=(yy/a)**2.d0*(24.d0+480.d0*yy**2.d0+6720.d0*yy**4.d0)
         a5=-(yy/a**2.d0)**2.d0*(24.d0+336.d0*yy**2.d0)
         a7=8.d0*yy**2.d0/a**6.d0
         b0=beta*(1.d0+8.d0*yy**2.d0+144.d0*yy**4.d0
     2            + 2880.d0*yy**6.d0 + 40320.d0*yy**8.d0)
         b2=-(yy/a)*(4.d0+72.d0*yy**2.d0+1440.d0*yy**4.d0
     2               + 20160.d0*yy**6.d0)
         b4=(yy/a**3.d0)*(6.d0+120.d0*yy**2.d0+1680.d0*yy**4.d0)
         b6=-(yy/a**5.d0)*(4.d0+56.d0*yy**2.d0)
         b8=yy/a**7.d0
         q(i,1)=(((b0+b2*xr**2.d0+b4*xr**4.d0+b6*xr**6.d0+b8*xr**8.d0)
     2           *dsin(xr/beta)
     2          +(a1*xr+a3*xr**3.d0+a5*xr**5.d0+a7*xr**7.d0)
     2           *dcos(xr/beta)) -
     2          ((b0+b2*xl**2.d0+b4*xl**4.d0+b6*xl**6.d0+b8*xl**8.d0)
     2           *dsin(xl/beta)
     2          +(a1*xl+a3*xl**3.d0+a5*xl**5.d0+a7*xl**7.d0)
     2           *dcos(xl/beta)))/dx
c        # Ought to make the zeroing more accurate
         if (dabs(xcell-x0) .gt. a) q(i,1) = 0.d0
         q(i,2) = q(i,1)/aux(i,1)
         go to 150

   80    continue
c        # C3 Newton polynomial modulated by a cosine
c        # Point values for finite difference method
         x0 = -4.d0
         xnode=xcell-0.5d0*dx
         q(i,1) = dcos((xnode-x0)/beta)*(xnode-x0-a)**4.d0
     &                 *(xnode-x0+a)**4.d0/a**8.d0
c        # Ought to make the zeroing more accurate
         if (dabs(xnode-x0) .gt. a) q(i,1) = 0.d0
         q(i,2) = q(i,1)/aux(i,1)
         go to 150

  90     continue
c        C5 Newton polynomial:
c        Initialization with exact cell averages
!        x0 =-4.d0
         xrelr=xcell+0.5d0*dx-x0
         xrell=xcell-0.5d0*dx-x0
c        #Point values would be given by:
c        # q(i,1) = ((xcell-x0)-a)^4*((xcell-x0)+a)^4/a**8.d0
c        #Initialization with exact cell averages:
         q(i,1) = (xrelr**13.d0/(13.d0*a**12.d0)
     2             - 6.d0*xrelr**11.d0/(11.d0*a**10.d0)
     2             + 5.d0*xrelr**9.d0/(3.d0*a**8.d0)
     2             - 20.d0*xrelr**7.d0/(7.d0*a**6.d0)
     2             + 3.d0*xrelr**5.d0/(a**4.d0) 
     2             - 2.d0*xrelr**3.d0/(a*a) + xrelr
     2            - (xrell**13.d0/(13.d0*a**12.d0)
     2             - 6.d0*xrell**11.d0/(11.d0*a**10.d0)
     2             + 5.d0*xrell**9.d0/(3.d0*a**8.d0)
     2             - 20.d0*xrell**7.d0/(7.d0*a**6.d0)
     2             + 3.d0*xrell**5.d0/(a**4.d0) 
     2             - 2.d0*xrell**3.d0/(a*a) + xrell))/dx
c        # Ought to make the zeroing more accurate
         if (dabs(xcell-x0) .gt. a) q(i,1) = 0.d0
         q(i,2) = q(i,1)/aux(i,1)
         go to 150

  100    continue
c        # Zero:
         q(i,1) = 0.d0
         q(i,2) = 0.d0
         if(i.eq.mx/2) then
           q(i,1) = 1.e-10
           q(i,2) = q(i,1)/aux(i,1)
         endif
         go to 150



  150    continue
c
      return
      end
