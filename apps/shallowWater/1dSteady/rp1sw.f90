! =============================================================================
  subroutine rp(ixy,maxmx,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
! =============================================================================
!
! # solve Riemann problems for the 1D shallow water equations
! #   (h)_t + (u h)_x = 0 
! #   (uh)_t + ( uuh + .5*gh^2 )_x = -gh(b)_x 
! # using Roe's approximate Riemann solver and f-wave approach.
! # Therefore the contribution of the source term in the discharge equation
! # is taken into account here.
!
! # On input, ql contains the state vector at the left edge of each cell
! #           qr contains the state vector at the right edge of each cell
! # On output, wave contains the f-waves, 
! #            s the speeds, 
! #            amdq the  left-going flux difference  A^- \Delta q
! #            apdq the right-going flux difference  A^+ \Delta q
!
! # Note that the i'th Riemann problem has left state qr(i-1,:)
! #                                    and right state ql(i,:)

!  implicit none
!  double precision :: ql(1-mbc:maxmx+mbc, meqn)
!  double precision :: qr(1-mbc:maxmx+mbc, meqn)
!  double precision :: s(1-mbc:maxmx+mbc, mwaves)
!  double precision :: wave(1-mbc:maxmx+mbc, meqn, mwaves)
!  double precision :: amdq(1-mbc:maxmx+mbc, meqn)
!  double precision :: apdq(1-mbc:maxmx+mbc, meqn)
!  double precision :: auxl(1-mbc:maxmx+mbc, *)
!  double precision :: auxr(1-mbc:maxmx+mbc, *)

!  double precision :: grav
!  common /param/ grav
  
 ! double precision :: delta(2)
 ! logical efix
 ! data efix /.true./    !# use entropy fix for transonic rarefactions
          

  return
  end subroutine rp



