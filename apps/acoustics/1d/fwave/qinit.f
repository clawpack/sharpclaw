c
c
c =========================================================
       subroutine qinit(maxmx,meqn,mbc,mx,xlower,dx,q,maux,aux)
c =========================================================
c
c     # Set initial conditions for q.
c
c
      implicit double precision (a-h,o-z)
      dimension q(1-mbc:mx+mbc, meqn)
      dimension aux(1-mbc:mx+mbc, 2)
      common /cqinit/ a,x0,ic
c
c
!       aux(i,1) is the density
!       aux(i,2) is the bulk modulus
!       eps = -p/K
!       m   = rho * u

      do 150 i=1,mx
	 xcell = xlower + (i-0.5d0)*dx

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
         q(i,1) = -q(i,1)/aux(i,2)  !Convert pressure to strain

         !For purely right-going pulse, m = eps*Z
         q(i,2) = -q(i,1) * dsqrt(aux(i,1)*aux(i,2)) 

  150    continue
c
      return
      end
