      subroutine setprob
      implicit double precision (a-h,o-z)
      character*16 fname

      common /comic/ beta,ic
      common /comrp/ u

c
      iunit = 7
      fname = 'setprob.data'
      
      call opendatafile(iunit, fname)


c     # parameters for wall
      read(7,*) u
      read(7,*) beta
      read(7,*) ic

      return
      end
