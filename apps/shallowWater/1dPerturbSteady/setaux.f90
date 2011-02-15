! ==================================================
  subroutine setaux(maxmx,mbc,mx,xlower,dx,maux,aux)
! ==================================================
!
! # Set auxiliary array for the following 1D shallow water equations 
! #
! # (h)_t + (h*u)_x = 0
! # (hu)_t + (h*u^2 + 1/2*grav*h^2)_x = -grav*h*(b)_x
! #
! # where b=b(x) represents the bottom topography.
!
! # The auxiliary array aux(i,1) contains the height of the bottom topography 
! # at the i-th cell center.
!
! # The settings used here are those reported in:
! # Y. Xing a, C.-W. Shu, High order well-balanced finite volume WENO schemes
! # and discontinuous Galerkin methods for a class of hyperbolic systems with source terms.
! # Journal of Computational Physics 214 (2006) 567–598. 
! # Section 6.1.3. A small perturbation of a steady state water.
   
  implicit none
      
  integer :: maxmx, mbc, maux
  integer :: mx(1)
  double precision:: xlower(1),dx(1)
  double precision :: aux(1-mbc:mx(1)+mbc, 1)
  
  
  integer :: i
  double precision :: xCell, pi
  
  
  pi = 4.d0*datan(1.d0)
  
  
  do i=1-mbc,mx(1)+mbc
  	xCell = xlower(1) + (i-0.5d0)*dx(1)
  	if (xCell .ge. 1.4d0 .and. xCell .le. 1.6d0) then
  		aux(i,1) = 0.25*(dcos(10*pi*(xCell-1.5d0)) + 1.d0)
    else
       	aux(i,1) = 0.d0
    endif
  enddo
  
 

  return     
  end subroutine setaux
