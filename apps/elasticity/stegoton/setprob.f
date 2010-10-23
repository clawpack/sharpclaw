      subroutine setprob
      implicit double precision (a-h,o-z)
      character*16 fname

      common /combc/ omega
      common /comwall/ pi,t1,a1,tw1
      common /combcu/ ubc(1:1000)
      common /combcn/ nt,ntmax
      common /comtimereverse/ trtime, trdone
      common /commat/ dKA, dKB, rhoA, rhoB

c
c     # Set the material parameters for the acoustic equations
c
      iunit = 7
      fname = 'setprob.data'
      
      call opendatafile(iunit, fname)

      pi = 4.d0*datan(1.d0)

c     # parameters for wall
      read(7,*) t1
      read(7,*) a1
      read(7,*) tw1
      read(7,*) dKA
      read(7,*) dKB
      read(7,*) rhoA
      read(7,*) rhoB
      read(7,*) trtime

c      open(unit=19,file='bc.data',status='old',form='formatted')
c      read(19,*) ntmax
c      do i=1,ntmax
c         read(19,*) ti,ubc(i)
c         enddo
c      nt = 0
c      write(6,*) 'Using bc.data...  ntmax = ', ntmax

      return
      end
