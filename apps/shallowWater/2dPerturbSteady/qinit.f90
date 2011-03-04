! ========================================================
  subroutine qinit(maxmx,meqn,mbc,mx,xlower,dx,q,maux,aux)
! ========================================================
!
! # Set initial conditions for the following 2D shallow water equations:
! #
! # (h)_t + (h*u)_x + (h*v)_y = 0
! # (hu)_t + (h*u^2 + 1/2*grav*h^2)_x + (h*u*v)_y = -grav*h*(b)_x
! # (hv)_t + (h*u*v)_x + (h*v^2 + 1/2*grav*h^2)_y = -grav*h*(b)_y
! #
! # where b=b(x,y) represents the bottom topography.
!
! # The settings used here are those reported in:
! # Y. Xing a, C.-W. Shu, High order well-balanced finite volume WENO schemes
! # and discontinuous Galerkin methods for a class of hyperbolic systems with source terms.
! # Journal of Computational Physics 214 (2006) 567–598. 
! # Section 6.1.3. A small perturbation of a two dimensional steady state water.

  	implicit none
  
  	integer :: maxmx, meqn, mbc, maux
  	integer :: mx(2)
  	double precision :: xlower(2), dx(2)

  	double precision, intent(out) :: q(1-mbc:mx(1)+mbc, 1-mbc:mx(2)+mbc, meqn)
  	double precision :: aux(1-mbc:mx(1)+mbc,1-mbc:mx(2)+mbc, maux)
  
 	double precision :: eps
  	common /perturb/ eps
  
    integer :: i, j
	double precision :: xCell, yCell
	

	do i=1-mbc,mx(1)+mbc
		xCell = xlower(1) + (i-0.5d0)*dx(1)
    	do j=1-mbc,mx(2)+mbc	
    		yCell = xlower(2) + (j-0.5d0)*dx(2)
    		
    		if (xCell .ge. 0.05d0 .and. xCell .le. 0.15d0) then
				q(i,j,1) = 1.0d0 + eps!- aux(i,j,1) + eps
			else
				q(i,j,1) = 1.0d0 !- aux(i,j,1)
			endif

			!q(i,j,1) = 1.d0
			q(i,j,2) = 0.d0 !h*u
			q(i,j,3) = 0.d0 !h*v
			
		enddo
	enddo
	
	return
	end subroutine qinit
