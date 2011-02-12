! ============================================================================
  subroutine rp(ixy,maxmx,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,fwave,s,amdq,apdq)
! ============================================================================

! # Solve Riemann problems for the 1D advection equation q_t + (u*q)_x = 0.

! --------------------------------------------------------------------
! # In conservation form, with CELL-CENTERED velocities specified in
! # the auxiliary variable
! # aux(i,1)  =  u-velocity in cell i
! # BE CAREFUL, UNDERSTAND WELL THE PROBLEM AS DESCRIBED IN
! # Finite volume method for hyperbolic problems, R.J. LeVeque
! --------------------------------------------------------------------

! # On input, ql contains the state vector at the left edge of each cell
! #           qr contains the state vector at the right edge of each cell
! # On output, wave contains the waves,
! #            s the speeds,
! #            amdq the  left-going flux difference  A^- \Delta q
! #            apdq the right-going flux difference  A^+ \Delta q
!
!     # Note that the i'th Riemann problem has left state qr(i-1,:)
!     #                                    and right state ql(i,:)
!     # From the basic clawpack routine step1, rp is called with ql = qr = q.

  implicit none
 
  integer :: mx(1)
  integer :: ixy, maxmx, meqn, mwaves, mbc
  double precision :: ql(1-mbc:mx(1)+mbc, meqn)
  double precision :: qr(1-mbc:mx(1)+mbc, meqn)
  double precision :: auxl(1-mbc:mx(1)+mbc, 1)
  double precision :: auxr(1-mbc:mx(1)+mbc, 1)
  double precision :: s(1-mbc:mx(1)+mbc, mwaves)
  double precision :: fwave(1-mbc:mx(1)+mbc, meqn, mwaves)
  double precision :: amdq(1-mbc:mx(1)+mbc, meqn)
  double precision :: apdq(1-mbc:mx(1)+mbc, meqn)
  
  integer :: i
  double precision :: ui, uim, qi, qim, fi, fim
  

  do i=2-mbc,mx(1)+mbc
  	ui = auxl(i,1)       ! Velocity in cell i
	uim = auxr(i-1,1)    ! Velocity in cell i-1
	qi = ql(i,1)         ! Density tracer in cell i
	qim = qr(i-1,1)      ! Density tracer in cell i-1
	
	fim = uim*qim
	fi = ui*qi
	
	fwave(i,1,1) = fi - fim
	
    if (ui .gt. 0.d0) then
	    s(i,1) = ui
	    amdq(i,1) = 0.d0
	    apdq(i,1) = fwave(i,1,1)
	else
	    s(i,1) = uim
	    amdq(i,1) = fwave(i,1,1)
	    apdq(i,1) = 0.d0
    endif
    
  enddo
  
  
  
  return

  end subroutine rp