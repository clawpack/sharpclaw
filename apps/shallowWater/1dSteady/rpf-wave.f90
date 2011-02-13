! =============================================================================
  subroutine rp(ixy,maxmx,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,fwave,s,amdq,apdq)
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
!
! # On output, wave contains the f-waves, 
! #            s the speeds, 
! #            amdq the  left-going flux difference  A^- \Delta q
! #            apdq the right-going flux difference  A^+ \Delta q
!
! # Note that the i'th Riemann problem has left state qr(i-1,:)
! #                                    and right state ql(i,:)

  implicit none
  
  integer :: ixy, maxmx, meqn, mwaves, mbc, mx(1)
  double precision :: ql(1-mbc:maxmx+mbc, meqn)
  double precision :: qr(1-mbc:maxmx+mbc, meqn)
  double precision :: s(1-mbc:maxmx+mbc, mwaves)
  double precision :: fwave(1-mbc:maxmx+mbc, meqn, mwaves)
  double precision :: amdq(1-mbc:maxmx+mbc, meqn)
  double precision :: apdq(1-mbc:maxmx+mbc, meqn)
  double precision :: auxl(1-mbc:maxmx+mbc, *)
  double precision :: auxr(1-mbc:maxmx+mbc, *)

  double precision :: grav
  common /param/ grav
  
  logical efix
  data efix /.true./    !# use entropy fix for transonic rarefactions
          
  double precision :: hl, ul, bl, hr, ur, br, hbar, chat, uhat
  double precision :: R(2,2)
  double precision :: L(2,2)
  double precision :: deltaf(2)
  double precision :: beta(2)
  
  
  integer :: i, j, k
  
  
  do i = 2-mbc, mx(1)+mbc
  	! Take left and right states at the interface
  	! Left state
  	hl = qr(i-1,1)
  	ul = qr(i-1,2)/hl
    bl = auxr(i-1,1)
  	
  	! Right state
  	hr = ql(i,1)
  	ur = ql(i,2)/hr
  	br = auxl(i,1)
  	
  	
  	! Contruct averaged quantities
  	hbar = (hl+hr)/2.d0
  	chat = dsqrt(grav*hbar) ! Roe average
  	uhat = (dsqrt(hl)*ul + dsqrt(hr)*ur)/(dsqrt(hl) + dsqrt(hr)) ! Roe average
  	
  	! Construct matrices of the left and theright vector
  	R(1,1) =  1.d0
  	R(2,1) = uhat - chat
  	R(1,2) = 1.d0
  	R(2,2) = uhat + chat
  	
  	L(1,1) = (chat + uhat)/(2.d0*chat)
  	L(1,2) = -1.d0/(2*chat)
  	
  	L(2,1) = (chat - uhat)/(2.d0*chat)
  	L(2,2) = 1.d0/(2*chat)	
  	
  	! Compute the flux differences
  	! This flux difference will be used to compute the f-waves
  	! Remember: flux for (h) : h*u
  	!           flux for (hu): h*u^2 + 1/2*g*h^2
  	deltaf(1) = (hr*ur) - (hl*ul)
  	deltaf(2) = (hr*ur*ur + 1/2*grav*hr*hr) - (hl*ul*ul + 1/2*grav*hl*hl) + 0.02*grav*hbar*(br-bl)
  	
  	
    ! Compute coefficients beta
  	do j = 1, meqn
  		beta(j) = 0.d0
  		do k = 1,meqn
  			beta(j) = beta(j) + L(j,k)*deltaf(k)
  	    enddo
    enddo
  	    
    ! Compute f-wave's and the waves speed
    do j = 1, meqn
  		do k = 1,meqn
  			fwave(i,j,k) = beta(k)*R(j,k)
  	    enddo
    enddo
    
    s(i,1) = uhat - chat
    s(i,2) = uhat + chat
   


 
    ! Finally compute the fluctuations
    if (efix)
    
    else
    	do  j = 1,meqn
    		amdq(i,j) = 0.d0
    		apdq(i,j) = 0.d0
    		do k = 1, meqn
    			if (s(i,k) .lt. 0.d0) then
    				amdq(i,j) = amdq(i,j) + fwave(i,j,k)
    			elseif (s(i,k) .gt. 0.d0) then
    				apdq(i,j) = apdq(i,j) + fwave(i,j,k)
    	    	else
            		amdq(i,j) = amdq(i,j) + 1/2*fwave(i,j,k)
					apdq(i,j) = apdq(i,j) + 1/2*fwave(i,j,k)
    			endif
    		enddo		
    	enddo
  endif
    	
    
  enddo
  	

  end subroutine rp



