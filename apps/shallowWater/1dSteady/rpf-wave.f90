! ================================================================================
	subroutine rp(ixy,maxmx,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,fwave,s,amdq,apdq)
! ================================================================================
!
! # Solve Riemann problems for the 1D shallow water equations
! # (h)_t + (h*u)_x = 0
! # (hu)_t + (h*u^2 + 1/2*grav*h^2) = -grav*h*(b)_x
! # using f-wave algorithm and Roe's approximate Riemann solver.  
! #
! # With the f-wave approach the source term in the discharge equation 
! # is treated here. Therefore, the user do not have to provide 
! # a subroutine which compute the contribution of the source term to the residual.
! # Thus the default (empty) src1.f90 subroutine is called.
!
! # On input, ql contains the state vector at the left edge of each cell
! #           qr contains the state vector at the right edge of each cell
! # On output, wave contains the waves, 
! #            s the speeds, 
! #            amdq the  left-going flux difference  A^- \Delta q
! #            apdq the right-going flux difference  A^+ \Delta q
!
! # Note that the i'th Riemann problem has left state qr(i-1,:)
! #                                    and right state ql(i,:)
! # From the basic clawpack routine step1, rp is called with ql = qr = q.
!
! # Here meqn=mwaves=2 should be passed from the calling routine

	implicit none 
	
	integer :: ixy, maxmx, meqn, mwaves, mbc
	integer :: mx(1)
    double precision :: ql(1-mbc:maxmx+mbc, meqn)
    double precision :: qr(1-mbc:maxmx+mbc, meqn)
    double precision :: auxl(1-mbc:maxmx+mbc, *)
    double precision :: auxr(1-mbc:maxmx+mbc, *)
    double precision :: s(1-mbc:maxmx+mbc, mwaves)
    double precision :: fwave(1-mbc:maxmx+mbc, meqn, mwaves)
    double precision :: amdq(1-mbc:maxmx+mbc, meqn)
    double precision :: apdq(1-mbc:maxmx+mbc, meqn)
    
    
    double precision :: grav
	common /comrp/ grav
	
	double precision :: hl, ul, bl, hr, ur, br, hbar, uhat, chat
	double precision :: R(2,2)
	double precision :: L(2,2)
	double precision :: fluxDiff(2), beta(2)
	
	integer :: i, j, k
	
    
    do i=2-mbc,mx(1)+mbc
    	! # Left states
		hl = qr(i-1,1)
		ul = qr(i-1,2)/hl
		bl = auxr(i-1,1)
		
		! # Right states
		hr = ql(i,1)
		ur = ql(i,2)/hr
		br = auxl(i,1)
		
		! # Average states (they come from the Roe's linearization)
		hbar = 0.5*(hr+hl)
		uhat = (dsqrt(hl)*ul + dsqrt(hr)*ur)/(dsqrt(hl)+dsqrt(hr))
		chat = dsqrt(grav*hbar)
		
		! # Flux differences + discretized source term
		fluxDiff(1) = (hr*ur) - (hl*ul)
		fluxDiff(2) = (hr*ur*ur + 0.5*grav*hr*hr) - (hl*ul*ul + 0.5*grav*hl*hl)
		
		! # Wave speeds
		s(i,1) = uhat-chat
		s(i,2) = uhat+chat
		
		! # Right eigenvectors
		R(1,1) = 1.d0
		R(2,1) = uhat-chat
		
		R(1,2) = 1.d0
		R(2,2) = uhat+chat
	
	
	    L(1,1) = (chat+uhat)/(2.d0*chat)
		L(2,1) = (chat-uhat)/(2.d0*chat)
		
		L(1,2) = -1d0/(2.d0*chat)
		L(2,2) = 1.d0/(2.d0*chat)
		
		! # Coefficients beta
		do j=1,meqn
			beta(j) = 0.d0
			do k=1,meqn
				beta(j) = beta(j) + L(j,k)*fluxDiff(k)
			enddo
		enddo

		! # Flux waves
		do j=1,meqn
			do k=1,meqn
				fwave(i,j,k) = beta(k)*R(j,k)
			enddo
		enddo
		
		
		! # Fluctuations
		do j=1,meqn
			amdq(i,j) = 0.d0
			apdq(i,j) = 0.d0
			do k=1,meqn
				if (s(i,k) .lt. 0.d0) then
					amdq(i,j) = amdq(i,j) + fwave(i,j,k)
				elseif (s(i,k) .ge.0) then
					apdq(i,j) = apdq(i,j) + fwave(i,j,k)
				else
					amdq(i,j) = amdq(i,j) + 0.5*fwave(i,j,k)
					apdq(i,j) = apdq(i,j) + 0.5*fwave(i,j,k)
				endif
			enddo
		enddo
	enddo
	
	return
	end subroutine rp