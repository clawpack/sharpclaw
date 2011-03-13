!     ============================================
      subroutine setaux(maxmx,mbc,mx,xlower,dx,maux,aux)
!     ============================================
!
!     # set auxiliary arrays 
!     # variable coefficient acoustics
!     #  aux(i,1) = density Z in i'th cell
!     #  aux(i,2) = bulk modulus c in i'th cell
!
!     # Piecewise constant medium with single interface at x=0
!     # Density and sound speed to left and right are set in setprob.f
!
!     
      implicit double precision (a-h,o-z)
      dimension aux(1-mbc:mx+mbc, 2)
      common /comaux/ rho(1000),c(1000),dinter(1000),ninter

      open(unit=31,file='fort.aux',status='unknown',form='formatted')

       ni=1
       do i=1-mbc,mx+mbc
          xcellr = xlower + i*dx
          do while (xcellr.gt.dinter(ni).and.ni.lt.ninter+1)
              ni=ni+1
          enddo
          aux(i,1) = rho(ni)
          aux(i,2) = c(ni)**2.d0*rho(ni)
        enddo


        do i=1,mx
          write(31,701) aux(i,1), aux(i,2)
  701     format(2e16.6)
          enddo

       close(unit=31)

       return
       end
