!     ============================================
      subroutine b4step(maxmx,mbc,mx,meqn,q,
     &            xlower,dx,t,dt,maux,aux)
!     ============================================
!
!     # called from claw1 before each call to step1.
!     # use to set time-dependent aux arrays or perform other tasks
!     # which must be done every time step.
!
!
!     
      implicit double precision (a-h,o-z)
      dimension q(1-mbc:maxmx+mbc, meqn)
      dimension aux(1-mbc:maxmx+mbc, *)
      common /comtimereverse/ trtime, itrdone

      if ((t+1.e-10).ge.trtime .and. itrdone.eq.0) then
        do i=1-mbc,mx+mbc
          q(i,2)=-q(i,2)
        enddo
        itrdone=1
      endif

!Output time traces at particular points

      i = 1
      sig = sigma(q(i,1),i,aux(i,2),aux(i,3))
      u = q(i,2)/aux(i,1)
      write(17,1005) t,q(i,1),q(i,2),sig,u
 1005 format(e16.8,4e16.8)

      i = mx
      sig = sigma(q(i,1),i,aux(i,2),aux(i,3))
      u = q(i,2)/aux(i,1)
      write(18,1005) t,q(i,1),q(i,2),sig,u

      i = 5
      sig = sigma(q(i,1),i,aux(i,2),aux(i,3))
      u = q(i,2)/aux(i,1)
      write(28,1005) t,q(i,1),q(i,2),sig,u

      return
      end
