! =====================================================================
  subroutine qinit(maxmx,meqn,mbc,mx,xlower,dx,q,maux,aux)
! =====================================================================
!
! # Set initial solution for the height (h) and the discharge (hu).
!
! # The bottom topography is an hump specified in the aux array (see setaux.f90).
! # See: Yulong Xing,  and Chi-Wang Shu, "High order well-balanced finite 
!        volume WENO schemes and discontinuous Galerkin methods for a class
!        of hyperbolic systems with source terms". 
!        Journal of Computational Physics, Volume 214, Issue 2, 20 May 2006, 
!        Pages 567-598
!     
!
  implicit none
      
  integer :: maxmx, meqn, mbc, mx(1), maux
  double precision :: xlower(1),dx(1) 
  double precision :: q(1-mbc:maxmx+mbc, meqn)
  double precision :: aux(1-mbc:maxmx+mbc, 1)
  
  double precision :: eps
  common /perturb/ eps
 
  integer :: i
  double precision :: xcell

  write(*,*) eps
	
  do i=1,mx(1)
  	q(i,1) = 1.d0 - aux(i,1)
    q(i,2) = 0.d0
	xcell = xlower(1) + (i-0.5d0)*dx(1)
	if (xcell.gt.1.1d0 .and. xcell.lt.1.2d0) then
		q(i,1) = q(i,1) + eps
	endif 
 enddo
      
 return
 end subroutine qinit