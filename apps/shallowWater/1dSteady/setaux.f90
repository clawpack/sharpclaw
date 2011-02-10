! ==================================================
  subroutine setaux(maxmx,mbc,mx,xlower,dx,maux,aux)
! ==================================================
!
! # Set auxiliary array for the shallow water equation with a bottom topography
! # B=B(x) 
! # aux array contains the bottom topography height at the cells faces
! # Why at the heigth at the interfaces? I find it more intuitive for the 
! # calculation of the bottom topography derivative. This derivative appears in 
! # the source term (RHS) of the shallow water for the discharge equation.
!     
  implicit none
      
  integer :: maxmx, mbc, mx(1), maux
  double precision:: dx(1)
  double precision :: xlower(1)
  double precision :: aux(1-mbc:maxmx+mbc, 1)

  double precision :: pi, xcell
  integer :: i
      
      
  pi = 4.d0*datan(1.d0)     

  do i=1-mbc,mx(1)+mbc
  	xcell = xlower(1) + (i-0.5d0)*dx(1)
    if (xcell .gt. 1.4 .and. xcell .lt. 1.6) then
        aux(i,1) = 0.25d0 * (dcos(10.d0*pi*(xcell-1.5d0)) + 1.d0)
    else
    	aux(i,1) = 0.d0
    endif
   write(*,*) aux(i,1)
  enddo
       
  end subroutine setaux
