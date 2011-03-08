! ====================
	subroutine setprob
! ====================
!
! # Set auxiliary array for the following 2D shallow water equations 
! #
! # (h)_t + (h*u)_x + (h*v)_y = 0
! # (hu)_t + (h*u^2 + 1/2*grav*h^2)_x + (h*u*v)_y = -grav*h*(b)_x
! # (hv)_t + (h*u*v)_x + (h*v^2 + 1/2*grav*h^2)_y = -grav*h*(b)_y
! #
! # where b=b(x,y) represents the bottom topography.
!
! # The value of the gravity acceleration is stored in grav
! # The value of the perturbation for the initial solution is stored in eps

	implicit none
      
	integer :: iunit
	character*12 fname
  
	double precision :: grav, eps
  	common /comrp/ grav
  	common /perturb/ eps

  	iunit = 7
  	fname = 'setprob.data'

  	call opendatafile(iunit, fname)
                
	! # Gravitational constant:
  	read(7,*) grav

	! # Perturbation for the initial condition
  	read(7,*) eps

  	return
  	end subroutine setprob
