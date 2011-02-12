! ==================
  subroutine setprob
! ==================

! # Set initial distribution of car density.
! # Set velocity of the cars.
  
  implicit none
  
  integer :: iunit
  character*12 fname
  
  double precision :: rho1,rho2,u1,u2
  common /param/ rho1,rho2,u1,u2

  iunit = 7
  fname = 'setprob.data'

  call opendatafile(iunit, fname)
                
! # Initial car density:
  read(7,*) rho1,rho2

! # Velocity
  read(7,*) u1,u2


  return
  end subroutine setprob
