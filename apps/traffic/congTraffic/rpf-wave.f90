! =============================================================================
  subroutine rp(ixy,maxmx,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,fwave,s,amdq,apdq)
! =============================================================================
!
! # Solve Riemann problems for the non-linear traffic flow equation 
! # q_t + (u_max(x)*q*(1-q))_x = 0,
! # with variable speed limit umax, using the f-wave approach.
! # Since u_max=u_max(x), the problem has a spatially varying flux
! # function. 
!
!
! # On input, ql contains the state vector at the left edge of each cell
! #           qr contains the state vector at the right edge of each cell
!
! # On output, wave contains the waves, 
! #            s the speeds, 
! #            amdq the  left-going flux difference  A^- \Delta q
! #            apdq the right-going flux difference  A^+ \Delta q
!
! # Note that the i'th Riemann problem has left state qr(i-1,:)
! #                                    and right state ql(i,:)
! # From the basic clawpack routine step1, rp is called with ql = qr = q.

  implicit none
  
  integer :: ixy, maxmx, meqn, mwaves, mbc
  integer :: mx(1)
  double precision :: ql(1-mbc:maxmx+mbc, meqn)
  double precision :: qr(1-mbc:maxmx+mbc, meqn)
  double precision :: s(1-mbc:maxmx+mbc, mwaves)
  double precision :: fwave(1-mbc:maxmx+mbc, meqn, mwaves)
  double precision :: amdq(1-mbc:maxmx+mbc, meqn)
  double precision :: apdq(1-mbc:maxmx+mbc, meqn)
  double precision :: auxl(1-mbc:maxmx+mbc, 1)
  double precision :: auxr(1-mbc:maxmx+mbc, 1)
   
  integer :: i
  double precision :: fim1, fi, sim1, si, dq, f0
  
  
  
  do i=2-mbc,mx(1)+mbc

  	! # Flux in each cell and flux difference
    fim1 = auxr(i-1,1)*qr(i-1,1)*(1.d0 - qr(i-1,1))
    fi = auxl(i,1)*ql(i,1)*(1.d0 - ql(i,1))
    fwave(i,1,1) = fi - fim1

    ! # Characteristic speed in each cell:
	sim1 = auxr(i-1,1)*(1.d0 - 2.d0*qr(i-1,1))
	si = auxl(i,1)*(1.d0 - 2.d0*ql(i,1))

    if (sim1 .lt. 0.d0 .and. si .le. 0.d0) then
    	! # left-going
        s(i,1) = sim1
	    amdq(i,1) = fwave(i,1,1)
	    apdq(i,1) = 0.d0
    elseif (sim1 .ge. 0.d0 .and. si .gt. 0.d0) then
    	! # right-going
        s(i,1) = si
	    amdq(i,1) = 0.d0
	    apdq(i,1) = fwave(i,1,1)
	elseif (sim1 .lt. 0.d0 .and. si .gt. 0.d0) then
    	! # transonic rarefaction
        ! # split fwave between amdq and apdq:
        s(i,1) = 0.5d0*(sim1 + si)
        dq = ql(i,1) - qr(i-1,1)

        ! # Entropy fix:  (perhaps doesn't work for all cases!!!)
        ! # This assumes the flux in the transonic case should
        ! # correspond to q=0.5 on the side with the smaller umax value.
        f0 = dmin1(auxr(i-1,1),auxl(i,1))*0.25d0
        amdq(i,1) = f0 - fim1
        apdq(i,1) = fi - f0

    else
    	! # Transonic shock
        s(i,1) = 0.5d0*(sim1 + si)
        if (s(i,1) .lt. 0.d0) then 
        	amdq(i,1) = fwave(i,1,1)
            apdq(i,1) = 0.d0
        elseif (s(i,1) .gt. 0.d0) then 
            amdq(i,1) = 0.d0
            apdq(i,1) = fwave(i,1,1)
        else
	        amdq(i,1) = 0.5d0 * fwave(i,1,1) 
	        apdq(i,1) = 0.5d0 * fwave(i,1,1)
        endif
    endif
  enddo
  
  return
  end subroutine rp

