c
c
c =========================================================
       subroutine qinit(maxmx,meqn,mbc,mx,xlower,dx,q,maux,aux)
c =========================================================
c
c     # Set initial conditions for q.
c     # Read in wave from file wave.data
c
c
      implicit double precision (a-h,o-z)
      dimension q(1-mbc:maxmx+mbc, meqn)
      dimension aux(1-mbc:maxmx+mbc, *)

      open(unit=22,file='wave.data',status='old',form='formatted')
c
c
      do i=1,5
         read(22,*) xjunk
         enddo

      do 150 i=1,mx
         read(22,*) q(i,1),q(i,2)
c
c        # isolate one soliton:
         xcell = xlower + (i-0.5d0)*dx
c        if (xcell.lt.150.d0 .or. xcell.gt.165.d0) then
         if (xcell.lt.35.d0 .or. xcell.gt.50.d0) then
            q(i,1) = 0.d0
            q(i,2) = 0.d0
            endif
         write(23,*) i,q(i,1),q(i,2)
  150    continue
c
      return
      end
