! =====================================================
	subroutine setaux(maxmx,mbc,mx,xlower,dx,maux,aux)
! =====================================================
!
! # Set auxiliary array for the following 2D shallow water equations 
! #
! # (h)_t + (h*u)_x + (h*v)_y = 0
! # (hu)_t + (h*u^2 + 1/2*grav*h^2)_x + (h*u*v)_y = -grav*h*(b)_x
! # (hv)_t + (h*u*v)_x + (h*v^2 + 1/2*grav*h^2)_y = -grav*h*(b)_y
! #
! # where b=b(x,y) represents the bottom topography.
!
! # The auxiliary array aux(i,1) contains the height of the bottom topography 
! # at the i-th cell center.
!
! # The settings used here are those reported in:
! # Y. Xing a, C.-W. Shu, High order well-balanced finite volume WENO schemes
! # and discontinuous Galerkin methods for a class of hyperbolic systems with source terms.
! # Journal of Computational Physics 214 (2006) 567–598. 
! # Section 6.2.3. A small perturbation of a two dimensional steady state water.
	
	implicit none
	
	integer :: maxmx, mbc, maux
	integer :: mx(2)
	double precision :: xlower(2), dx(2)
    double precision :: aux(1-mbc:mx(1)+mbc,1-mbc:mx(2)+mbc, *)
    
    integer :: i, j
	double precision :: xCell, yCell
	
	do i=1-mbc,mx(1)+mbc
		xCell = xlower(1) + (i-0.5d0)*dx(1)
    	do j=1-mbc,mx(2)+mbc
    		yCell = xlower(2) + (j-0.5d0)*dx(2)
    		aux(i,j,1) = 0.8d0*dexp(-5.d0*(xCell-0.9d0)**2-50.d0*(yCell-0.5d0)**2)
    	enddo
    enddo
    
    return
    end subroutine setaux
