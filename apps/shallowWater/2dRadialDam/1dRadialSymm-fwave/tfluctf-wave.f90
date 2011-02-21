! =======================================================================
	subroutine tfluct(ixy,maxmx,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,s,adq)
! =======================================================================
!
! # Compute the total fluctucations for the 2D shallow water equations
! # in radial symmetry
! #
! # (h)_t + (h*u)_r = -h*u/r
! # (hu)_t + (h*u^2 + 1/2*grav*h^2)_r = -h*u^2/r
! #
! # by solving the Riemann problem with the f-wave algorithm and the
! # Roe's approximate Riemann solver.
! # The initial condition for the Riemann problem are: left state  ql(i,:)
! #                                                    right state qr(i,:)
! #
! # With the f-wave approach the source term in the discharge equation 
! # is treated here. Therefore, the user do not have to provide 
! # a subroutine which compute the contribution of the source term to the residual.
! # Thus the default (empty) src1.f90 subroutine is called.
! #
! # On output, s contains the speeds, 
! #            adq the total fluctuation A \Delta q
!
!
! # From the basic clawpack routine step1, rp is called with ql = qr = q.

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
    double precision :: adq(1-mbc:maxmx+mbc, meqn)
    
    
    double precision :: grav
	common /comrp/ grav
	
	double precision :: hl, ul, hr, ur, hbar, uhat, chat
	double precision :: R(2,2)
	double precision :: L(2,2)
	double precision :: fluxDiff(2), beta(2)
	
	integer :: i, j, k
	
    
    do i=2-mbc,mx(1)+mbc
    	! # Left states
		hl = ql(i,1)
		ul = ql(i,2)/hl
		
		! # Right states
		hr = qr(i,1)
		ur = qr(i,2)/hr
		
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
		
		do j=1,meqn
			adq(i,j)=amdq(i,j) + apdq(i,j)
		enddo
			
		
		
	enddo
	
	return
	end subroutine tfluct
