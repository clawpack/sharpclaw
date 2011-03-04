! ================================================================================
	subroutine rp(ixy,maxmx,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,fwave,s,amdq,apdq)
! ================================================================================

! # Solve Riemann problems for the 2D shallow water equations
! #
! # (h)_t + (h*u)_x = 0
! # (hu)_t + (h*u^2 + 1/2*grav*h^2)_x + (h*u*v)_y = -grav*h*(b)_x
! # (hv)_t + (h*u*v)_x + (h*v^2 + 1/2*grav*h^2)_y = -grav*h*(b)_y
! #
! # using f-wave algorithm and Roe's approximate Riemann solver.  
! #
! # With the f-wave approach the source term in the discharge equation 
! # should be treated here. However, its contribution has been taken into account 
! # in the tfluctf-wave.f90 subroutine, which solves a Riemann problem at the interface.!
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



	implicit none 
	
	integer :: ixy, maxmx, meqn, mwaves, mbc
	integer :: mx
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
	
	double precision :: hl, ul, vl, bl, hr, ur, vr, br, hbar, uhat, vhat, chat, n_1, n_2
	double precision :: R(3,3)
	double precision :: L(3,3)
	double precision :: fluxDiff(3), beta(3)
	
	
	integer :: i, j, k
		
	! # Define normal and tangential directions of the current slice. This is needed
	! # to use the right components of the system that correspond to the momentum in
	! # the normal and tangential directions. Since this Riemann solver is used for
	! # both x- and y-slices we do not use names "n" for the normal component and "t"
	! # for the tangential component, but we simply define a vector n = (n_1,n_2)^T
	! # whose component depends on the current slice direction.
	! # 
	! # Look for instance: D. Ambrosi
	! #                    Approximation of shallow water equations by Roe's Riemann solver
	! #                    International Journal of numerical methods in fluids
	! #                    VOL. 20, 157-168 (1995)
	! #
    if (ixy .eq. 1) then
	  n_1 = 1.d0  ! # This means vertical interface,
	  n_2 = 0.d0  ! # e.g. normal vector pointing in the x direction
	else
	  n_1 = 0.d0 ! # This means horizontal interface,
	  n_2 = 1.d0 ! # e.g. normal vector pointing in the y direction
	endif

	!Now, it is possible to loop over the cells in this slice
    do i=2-mbc,mx+mbc
      	! # Left states
		hl = qr(i-1,1)
		ul = qr(i-1,2)/hl
		vl = qr(i-1,3)/hl
		bl = auxr(i-1,1)
		
		! # Right states
		hr = ql(i,1)
		ur = ql(i,2)/hr
		vr = ql(i,3)/hr
		br = auxl(i,1)
	
		! # Roe average states (Roe's linearization)
		hbar = 1.d0/2.d0*(hr + hl)
		uhat = (dsqrt(hr)*ur + dsqrt(hl)*ul)/(dsqrt(hr) + dsqrt(hl)) 
		vhat = (dsqrt(hr)*vr + dsqrt(hl)*vl)/(dsqrt(hr) + dsqrt(hl))
		chat = dsqrt(grav*hbar)
		
		! # Flux differences
		! ##################
		! # In order to compute these differences in an efficient way (avoid if statements)
		! # we compute the scalar product between the velocity vector at the interface and the
		! # normal vector "{n} = {n_1,n_2}^T" which has been defined above.
		fluxDiff(1) = hr*(ur*n_1 + vr*n_2) - hl*(ul*n_1 + vl*n_2) ! Easy to understand
		
		! # If the slice is in the y-direction we have:
		! # fluxDiff(2) = (hr*ur^2 + 1/2*grav*hr^2) - (hl*ul^2 + 1/2*grav*hl^2)
		! # fluxDiff(3) = (hr*ur*vr) - (hl*ul*vl)
		! #
		! # If the slice is in the x-direction we have:
		! # fluxDiff(2) = (hr*ur*vr) - (hl*ul*vl)
		! # fluxDiff(3) = (hr*vr^2 + 1/2*grav*hr^2) - (hl*vl^2 + 1/2*grav*hl^2)
		! #
		! # Using the vector component n_1 and n_2 defined above,
		! # this two possibilities can be achieved in the following way:
		fluxDiff(2) = (hr*ur*(ur*n_1 + vr*n_2)+0.5*grav*hr*hr*n_1)-(hl*ul*(ul*n_1 + vl*n_2)+0.5*grav*hl*hl*n_1)!+grav*hbar*(br-bl)*n_1
		
		fluxDiff(3) = (hr*vr*(ur*n_1 + vr*n_2)+0.5*grav*hr*hr*n_2)-(hl*vl*(ul*n_1 + vl*n_2)+0.5*grav*hl*hl*n_2)!+grav*hbar*(br-bl)*n_2
		
		! # Wave speeds
		s(i,1) = (uhat*n_1 + vhat*n_2) - chat
		s(i,2) = (uhat*n_1 + vhat*n_2)
		s(i,2) = (uhat*n_1 + vhat*n_2) + chat
		
		! # Right eigenvectors (columns)
		R(1,1) = 1.d0
		R(2,1) = uhat - chat*n_1
		R(3,1) = vhat - chat*n_2
		
		R(1,2) = 0.d0
		R(2,2) = -n_2
		R(3,2) = n_1
		
		R(1,3) = 1.d0
		R(2,3) = uhat + chat*n_1
		R(3,3) = vhat + chat*n_2
		
	
		! # Left eigenvectors (rows)
	    L(1,1) = ((uhat*n_1 + vhat*n_2) + chat)/(2.d0*chat)
	    L(1,2) = -n_1/(2.d0*chat)
	    L(1,3) = -n_2/(2.d0*chat)
	    
		L(2,1) = uhat*n_2 - vhat*n_1 
		L(2,2) = -n_2
		L(2,3) = n_1
		
		L(3,1) = (chat - (uhat*n_1 + vhat*n_2))/(2.d0*chat)
		L(3,2) = n_1/(2.d0*chat)
		L(3,3) = n_2/(2.d0*chat)
		
		
		! # Coefficients beta which multiply the right eigenvectors (see below)
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
				elseif (s(i,k) .ge. 0.d0) then
					apdq(i,j) = apdq(i,j) + fwave(i,j,k)
				else
					amdq(i,j) = amdq(i,j) + 1.d0/2.d0*fwave(i,j,k)
					apdq(i,j) = apdq(i,j) + 1.d0/2.d0*fwave(i,j,k)
		 		endif
			enddo
		enddo
		
	enddo
	
	return
	end subroutine rp
