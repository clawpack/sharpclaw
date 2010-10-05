!     ============================================
      subroutine b4step(maxnx,mbc,nx,meqn,q, &
                 xlower,dx,time,dt,maux,aux)
!     ============================================
!
!     # called before each call to step
!     # use to set time-dependent aux arrays or perform other tasks.
!
!     # make velocity time dependent, reversing flow.
!
!     
      implicit double precision (a-h,o-z)
      integer :: maxnx, nx(2)
      dimension q(1-mbc:nx(1)+mbc,1-mbc:nx(2)+mbc, meqn)
      dimension aux(1-mbc:nx(1)+mbc,1-mbc:nx(2)+mbc, maux)
      double precision :: xlower(2), dx(2)
!      external psi
      common /comvt/ tperiod,pi2
!
      if (tperiod .eq. 0.d0) then
!        # special case --- indication that velocities specified in 
!        # setaux should be used for all time.
         return
         endif

       vt = dcos(pi2*(time+dt/2.d0)/tperiod)
!
       do i=1-mbc,nx(1)+mbc
          do j=1-mbc,nx(2)+mbc
!           # coordinates of lower left corner of grid cell:
            xll = xlower(1) + (i-1)*dx(1)
            yll = xlower(2) + (j-1)*dx(2)

!           # difference stream function psi to get normal velocities:
            aux(i,j,1) = -(psi(xll, yll+dx(2)) - psi(xll,yll)) / dx(2)
            aux(i,j,2) =  (psi(xll+dx(1), yll) - psi(xll,yll)) / dx(1)
!
!           # multiply by time-factor:
            aux(i,j,1) = vt * aux(i,j,1)
            aux(i,j,2) = vt * aux(i,j,2)
        end do
      end do
      return
      end
