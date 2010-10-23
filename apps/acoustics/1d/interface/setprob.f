      subroutine setprob
      implicit double precision (a-h,o-z)
      character*16 fname
      common /cqinit/ a,x0,ic
      common /comaux/ rho(1000),c(1000),dinter(1000),ninter
c
c     # Set the material parameters for the acoustic equations
c
      iunit = 7
      fname = 'setprob.data'
      
      call opendatafile(iunit, fname)

c     # choice of initial data:
      read(7,*) ic
      read(7,*) a
      read(7,*) x0
c     # Number of interfaces:
      read(7,*) ninter
c     # interface location:
      read(7,*) (dinter(ni), ni=1,ninter)
c     # Piecewise constant medium with single interface at x=0
c     # Density and sound speed to left and right:
      read(7,*) (rho(ni), ni=1,ninter+1)
      read(7,*) (c(ni), ni=1,ninter+1)
c
      return
      end
