!
!
! =========================================================
       subroutine qinit(maxmx,meqn,mbc,mx,xlower,dx,q,maux,aux)
! =========================================================
!
!     # Set initial conditions for q.
!
!
      implicit double precision (a-h,o-z)
      dimension q(1-mbc:maxmx+mbc, meqn)
      dimension aux(1-mbc:maxmx+mbc, *)

      integer :: mx(1)
      double precision :: xlower(1),dx(1)

!
!
      ql = 0.2d0
      qr = 0.1d0

      do i=1,mx(1)
        xcell = xlower(1) + (i-0.5d0)*dx(1)
        if (xcell .lt. 0.d0) then
            q(i,1) = ql
        else
            q(i,1) = qr
        endif

      enddo
!
      
      end subroutine qinit
