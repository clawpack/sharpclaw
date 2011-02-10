! ========================================
      subroutine setprob
! ========================================
      implicit double precision (a-h,o-z)
      character*12 fname
      common /param/ grav
      common /perturb/ eps

!    # Set some parameters for the 1D shallow water equations

      iunit = 7
      fname = 'setprob.data'

      call opendatafile(iunit, fname)
                
!     # Gravitational constant:
      read(7,*) grav

!     # Data for the initial condition. It is used in qinit.f90:
      read(7,*) eps


      return
      end subroutine setprob
