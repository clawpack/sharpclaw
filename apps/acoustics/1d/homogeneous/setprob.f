      subroutine setprob
      implicit double precision (a-h,o-z)
      character*16 fname
      common /cparam/ rho,bulk,cc,zz   
      common /cqinit/ beta
c
c     # Set the material parameters for the acoustic equations
c     # Passed to the Riemann solver rp1.f in a common block
c
      iunit = 7
      fname = 'setprob.data'
      
      call opendatafile(iunit, fname)

c     # density:
      read(7,*) rho

c     # bulk modulus:
      read(7,*) bulk

c     # sound speed:
      cc = dsqrt(bulk/rho)

c     # impedance:
      zz = cc*rho

c     # beta for initial conditions:
      read(7,*) beta

      return
      end
