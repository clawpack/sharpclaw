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
      dimension q(1-mbc:maxmx+mbc, meqn)
      dimension aux(1-mbc:maxmx+mbc, *)
c
c
      do 150 i=1,mx
	 xcell = xlower + (i-0.5d0)*dx
!        q(i,1) = g0((xcell-100.5d0)/10.d0)
         q(i,1) = 0.d0
!        if (xcell .gt. 100.d0 .and. xcell .lt. 101.d0) q(i,1) = 1.d0
         q(i,2) = 0.d0
!        sig = 1.0*dexp(-((xcell-75.d0)/10.d0)**2.d0)
!        q(i,1) = dlog(sig+1)/aux(i,2)
  150    continue
c
      return
      end
