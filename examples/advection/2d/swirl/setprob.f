      subroutine setprob
      implicit double precision (a-h,o-z)
      character*16 fname
      common /compsi/ pi
      common /comvt/ tperiod,pi2

      iunit = 7
      fname = 'setprob.data'
      
      call opendatafile(iunit, fname)
c
c     # compute pi, used in psi.f
      pi = 4.d0 * datan(1.d0)
c
c     # save 2*pi and tperiod in common block for use in b4step2:
c
      pi2 = 2.d0*pi
c
      read(7,*) tperiod

      return
      end
