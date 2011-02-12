
c
c
c     =====================================================
      subroutine tfluct(ixy,maxm,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,
     &                   s,adq)
c     =====================================================
! # Solve Riemann problems for the 1D advection equation q_t + (u*q)_x = 0.

! --------------------------------------------------------------------
! # In conservation form, with CELL-CENTERED velocities specified in
! # the auxiliary variable
! # aux(i,1)  =  u-velocity in cell i
! # BE CAREFUL, UNDERSTAND WELL THE PROBLEM AS DESCRIBED IN
! # Finite volume method for hyperbolic problems, R.J. LeVeque
! --------------------------------------------------------------------

c     # On input, ql contains the state vector at the left edge of each cell
c     #           qr contains the state vector at the right edge of each cell
c
c     # On output, fwave contains the waves as jumps in f,
c     #            s the speeds,
c     #
c     #            amdq = A^- Delta q, 
c     #            apdq = A^+ Delta q,
c     #                   the decomposition of the flux difference
c     #                       f(qr(i-1)) - f(ql(i))
c     #                   into leftgoing and rightgoing parts respectively.
c     #
c
c     # Note that the i'th Riemann problem has left state ql(i,:)
c     #                                    and right state qr(i,:)
c     # From the basic clawpack routines, this routine is called with ql = qr
c
c
      implicit double precision (a-h,o-z)
c
      integer :: mx(1)
      integer :: ixy, maxmx, meqn, mwaves, mbc
      double precision :: ql(1-mbc:mx(1)+mbc, meqn)
      double precision :: qr(1-mbc:mx(1)+mbc, meqn)
      double precision :: auxl(1-mbc:mx(1)+mbc, 1)
      double precision :: auxr(1-mbc:mx(1)+mbc, 1)
      double precision :: s(1-mbc:mx(1)+mbc, mwaves)
      double precision :: wave(1-mbc:mx(1)+mbc, meqn, mwaves)
      double precision :: adq(1-mbc:mx(1)+mbc, meqn)
  
      integer :: i
      double precision :: ui, uim, qi, qim, qstar

      do i=2-mbc,mx(1)+mbc
        ui = auxr(i,1)       ! Velocity at right edge of cell i
        uim = auxl(i,1)    ! Velocity at left edge of cell i
        qi = qr(i,1)         ! Density tracer at right edge of cell i
        qim = ql(i,1)      ! Density tracer at left edge of cell i
        if (ui .gt. 0.d0) then
            qstar = uim*qim/ui   ! The normal flux (f = uq) MUST BE CONTINUOUS,                                       ! i.e. uim * qim = ui * qi

            wave(i,1,1) = qi - qstar
            s(i,1) = ui
            adq(i,1) = ui*qi - uim*qim
        else
            qstar = ui*qi/uim
            wave(i,1,1) = qstar - qim
            s(i,1) = uim
            adq(i,1) = ui*qi - uim*qim
        endif
      enddo
 
      return
      end subroutine tfluct
