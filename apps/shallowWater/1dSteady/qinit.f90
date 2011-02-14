! ========================================================
  subroutine qinit(maxmx,meqn,mbc,mx,xlower,dx,q,maux,aux)
! ========================================================
!
! # Set initial conditions for the following 1D shallow water equations:
! #
! # (h)_t + (h*u)_x = 0
! # (hu)_t + (h*u^2 + 1/2*grav*h^2) = -grav*h*(b)_x
! #
! # where b=b(x) represents the bottom topography.
!
! # The settings used here are those reported in:
! # Y. Xing a, C.-W. Shu, High order well-balanced finite volume WENO schemes
! # and discontinuous Galerkin methods for a class of hyperbolic systems with source terms.
! # Journal of Computational Physics 214 (2006) 567–598. 
! # Section 6.1.3. A small perturbation of a steady state water.


  implicit none
  
  integer :: maxmx, meqn, mbc, maux
  integer :: mx(1)
  double precision :: xlower(1), dx(1)

  double precision, intent(out) :: q(1-mbc:mx(1)+mbc, meqn)
  double precision :: aux(1-mbc:mx(1)+mbc, *)
  
  double precision :: eps
  common /perturb/ eps
  
  integer :: i
  double precision :: xCell
  
  do i=1,mx(1)
  	q(i,1) = 1.d0 - aux(i,1)
  	q(i,2) = 0.d0
  	
  	!xCell = xlower(1) + (i-0.5d0)*dx(1)
    !if (xCell .ge. 1.1d0 .and. xCell .le. 1.2d0) then
    !	q(i,1) = q(i,1) !+ eps
    !endif
  enddo

  return
  end subroutine qinit
