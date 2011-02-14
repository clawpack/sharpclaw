! ==================
  subroutine setprob
! ==================
!
! # Set problem parameters for the following 1D shallow water equations 
! #
! # (h)_t + (h*u)_x = 0
! # (hu)_t + (h*u^2 + 1/2*grav*h^2) = -grav*h*(b)_x
! #
! # where b=b(x) represents the bottom topography.
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
