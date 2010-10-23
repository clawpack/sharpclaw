!     ============================================
      subroutine setaux(maxmx,mbc,mx,xlower,dx,maux,aux)
!     ============================================
!
!     # set auxiliary arrays 
!     # variable coefficient acoustics
!     #  aux(i,1) = rho in ith cell
!     #  aux(i,2) = coefA (nonlinear bulk modulus)
!     #  aux(i,3) = coefB (linear bulk modulus)
!
!     # rho value for i-1 < x < i is determined by rho(i)
!     # input from setprob.rho in setprob.f
!
!     
!     TODO:
!       - Add capability for many layers
!       - Add smoothly varying media as run-time option

      implicit double precision (a-h,o-z)
      dimension aux(1-mbc:maxmx+mbc, 3)
      common /commat/ dKA, dKB, rhoA, rhoB

      open(unit=31,file='fort.aux',status='unknown',form='formatted')
!

       pi = datan(1.d0)*4.d0
       xupper = xlower + mx*dx
       do 100 i=1-mbc,mx+mbc
          xcell = xlower + (i-0.5d0)*dx
!
!         # discrete layers:
!         # layer A between j and j+frac1
!         # layer B between j+frac1 and j+1

          frac1 = 1.d0/2.d0
          frac2 = 3.d0/2.d0
          frac3 = 3.d0/2.d0

          ix = xcell
          xfrac = xcell - ix
          if (xfrac .lt. frac1) then
!     &          .and.xcell .lt. (xupper-xlower)*0.5d0) then
!             # layer A:
              aux(i,1) = rhoA
              aux(i,2) = dKA
              aux(i,3) = 0.d0
            elseif (xfrac .lt. frac2) then
!     &         .and. xcell .lt. (xupper-xlower)*0.5d0) then
!             # layer B:
              aux(i,1) = rhoB
              aux(i,2) = dKB
              aux(i,3) = 0.0d0
!           elseif (xfrac .lt. frac3) then
!             aux(i,1) = 2.0d0
!             aux(i,2) = 0.d0
!             aux(i,3) = 2.0d0
!            else
!             # linear material
!              aux(i,1) = 10.0d0
!              aux(i,2) = 0.d0
!              aux(i,3) = 10.d0
            endif

!         For sinotons:
!         aux(i,1)=2.5d0+1.5d0*(dsin(2.d0*pi*xcell))**1.d0
!         aux(i,2)=2.5d0+1.5d0*(dsin(2.d0*pi*(xcell+0.0d0)))**1.d0
!         aux(i,3)=0.d0
  100     continue

!      # make material uniform in ghost cells:
       do ibc=1,mbc
          aux(1-ibc,1) = aux(1,1)
          aux(mx+ibc,1) = aux(mx,1)
          aux(1-ibc,2) = aux(1,2)
          aux(mx+ibc,2) = aux(mx,2)
          aux(1-ibc,3) = aux(1,3)
          aux(mx+ibc,3) = aux(mx,3)
          enddo

	do i=1,mx
          write(31,701) aux(i,1),aux(i,2),aux(i,3)
  701     format(3e16.6)
          enddo

       close(unit=31)
!
       return
       end
