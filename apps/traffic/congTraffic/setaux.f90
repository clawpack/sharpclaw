!     ============================================
      subroutine setaux(maxmx,mbc,mx,xlower,dx,maux,aux)
!     ============================================
!
!     # Set auxiliary array for traffic flow with variable speed limit
!     
      implicit double precision (a-h,o-z)
      
      integer :: mx(1)
      double precision:: dx(1)
      double precision :: xlower(1)
      dimension aux(1-mbc:maxmx+mbc, 1)

      do i=1-mbc,mx(1)+mbc
        xcell = xlower(1) + (i-0.5d0)*dx(1)
        if (xcell .lt. 0.d0) then
             aux(i,1) = 2.0d0
         else
             aux(i,1) = 1.0d0
         endif
      enddo

       
       end subroutine setaux
