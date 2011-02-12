! ==================================================
  subroutine setaux(maxmx,mbc,mx,xlower,dx,maux,aux)
! ==================================================
!
! # Set auxiliary array for traffic flow with variable speed limit umax.
! # The auxiliary array contains the distribution of the velocity
!
!     
  implicit none
      
  integer :: maxmx, mbc, maux
  integer :: mx(1)
  double precision:: xlower(1),dx(1)
  double precision :: aux(1-mbc:maxmx+mbc, 1)

  integer :: i
  double precision :: xCell
  
  do i=1-mbc,mx(1)+mbc
  	xCell = xlower(1) + (i-0.5d0)*dx(1)
  		if (xCell .le. 0.d0) then
  			aux(i,1) = 2.d0
        else
        	aux(i,1) = 1.d0
        endif
  enddo

  return     
  end subroutine setaux
