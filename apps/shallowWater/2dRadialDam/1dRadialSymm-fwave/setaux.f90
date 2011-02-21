! ==================================================
  subroutine setaux(maxmx,mbc,mx,xlower,dx,maux,aux)
! ==================================================
!
! # Set auxiliary array for the following 1D shallow water equations 
! #
! # (h)_t + (r*h*u)_r = 0
! # (hu)_t + (r*h*u^2 + r*1/2*grav*h^2)_r = -1/2*grav*h^2
! #
! # where r is the radial position.
!
! # The auxiliary array aux(i,1) contains the height of the bottom topography 
! # at the i-th cell center.
   
  implicit none
      
  integer :: maxmx, mbc, maux
  integer :: mx(1)
  double precision:: xlower(1),dx(1)
  double precision :: aux(1-mbc:mx(1)+mbc, *)
  
  
  integer :: i
  double precision :: xInterface

  do i=1-mbc,mx(1)+mbc
  	aux(i,1) = xlower(1) + (i-0.5d0)*dx(1)
  	aux(i,2) = dx(1)
  enddo
  
 

  return     
  end subroutine setaux