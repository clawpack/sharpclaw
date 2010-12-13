! ========================================
      subroutine setprob
! ========================================
      implicit double precision (a-h,o-z)
      character*12 fname
      common /comrp/ grav
      common /cdisc/ x0,y0,alf,beta,r0,idisc
      common /comic/ hin,hout

!    # Set the material parameters for the acoustic equations
!    # Passed to the Riemann solver rp1.f in a common block


      iunit = 7
      fname = 'setprob.data'

      call opendatafile(iunit, fname)
                
!     # gravitational constant:
      read(7,*) grav

!     # data for radial dam-break problem:
      idisc=2 ! Hard-coded (should it be an input parameter??)
      read(7,*) x0,y0
      read(7,*) r0
      read(7,*) hin
      read(7,*) hout


!
      return
      end subroutine setprob 
