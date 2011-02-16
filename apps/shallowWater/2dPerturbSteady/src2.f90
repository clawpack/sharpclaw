! =========================================================
  subroutine src(q,dq,aux,t,dt)
! =========================================================
! # Set problem parameters for the following 2D shallow water equations 
! #
! # (h)_t + (h*u)_x + (h*v)_y = 0
! # (hu)_t + (h*u^2 + 1/2*grav*h^2)_x + (h*u*v)_y = -grav*h*(b)_x
! # (hv)_t + (h*u*v)_x + (h*v^2 + 1/2*grav*h^2)_y = -grav*h*(b)_y
! #
! # where b=b(x,y) represents the bottom topography.
!
! # Note that in general src() should evaluate 
! # (delta t) * dq/dt = (delta t) * psi(t)
! # instantaneously (not integrate over a timestep)
! # This differs from the function src() in CLAWPACK


  	Use ClawParams
      
  	implicit none

  	double precision, intent(in) :: q(1-mbc:nx(1)+mbc, 1-mbc:nx(2)+mbc,  meqn)
	double precision, intent(out) :: dq(1-mbc:nx(1)+mbc, 1-mbc:nx(2)+mbc, meqn)
  	double precision :: aux(1-mbc:nx(1)+mbc, 1-mbc:nx(2)+mbc, maux)
  	double precision :: dt, t
  	
  	
  	double precision :: grav
  	common /comrp/ grav
  	
  	double precision :: b_x, b_y

  	integer :: i, j
  
  	do i=1,nx(1)
  		do j=1,nx(2)
    		dq(i,j,1) = 0.d0
    		!dq(i,j,2) = 0.d0
    		!dq(i,j,3) = 0.d0
    		b_x = (aux(i,j,1)-aux(i-1,j,1))/dx(1)
    		dq(i,j,2) = -1.0d0*grav*q(i,j,1)*b_x*dt !src = -g*h*(bottom)_x
    		b_y = (aux(i,j,1)-aux(i,j-1,1))/dx(2)
    		dq(i,j,3) = -1.0d0*grav*q(i,j,1)*b_y*dt !src = -g*h*(bottom)_y
    	enddo
  	enddo
  

  	return
  	end subroutine src