! ====================
  subroutine setprob
! ====================

! # Set the velocity for scalar advection
! # This value is passed to the Riemann solver rp....f90 in a common block
  
  implicit none
  
  double precision :: u1,u2,rho1,rho2
  common /comtraf/ u1,u2,rho1,rho2
  
  character*12 fname
  integer :: iunit
  
  iunit = 7
  fname = 'setprob.data'

  call opendatafile(iunit, fname)

  read(7,*) u1,u2
  read(7,*) rho1,rho2

  end subroutine setprob