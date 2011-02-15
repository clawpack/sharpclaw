! =========================================================
  	subroutine src(q,dq,aux,t,dt)
! =========================================================
! # Set problem parameters for the following 1D shallow water equations 
! #
! # (h)_t + (h*u)_x = 0
! # (hu)_t + (h*u^2 + 1/2*grav*h^2)_x = -grav*h*(b)_x
! #
! # where b=b(x) represents the bottom topography.
!
! # Note that in general src() should evaluate 
! # (delta t) * dq/dt = (delta t) * psi(t)
! # instantaneously (not integrate over a timestep)
! # This differs from the function src() in CLAWPACK


  	Use ClawParams
      
  	implicit none

  	double precision, intent(in) :: q(1-mbc:nx(1)+mbc, meqn)
	double precision, intent(out) :: dq(1-mbc:nx(1)+mbc, meqn)
  	double precision :: aux(1-mbc:nx(1)+mbc, *)
  	double precision :: dt, t

  	double precision :: grav
  	common /comrp/ grav

  	integer :: i
  
  	do i=1,nx(1)
    	dq(i,1) = 0.d0
    	dq(i,2) = -1.0d0*grav*q(i,1)*(aux(i,1)-aux(i-1,1))/dx(1)*dt !src = -g*h*(bottom)_x
  	enddo
  

  	return
  	end subroutine src

