! =========================================================
  subroutine qinit(maxmx,meqn,mbc,mx,xlower,dx,q,maux,aux)
! =========================================================
!
! # Set initial conditions for q, i.e. initial distribution of cars density

  implicit none 
  
  integer :: maxmx, meqn, mbc, maux
  integer :: mx(1)
  double precision :: xlower(1),dx(1)
  double precision :: q(1-mbc:maxmx+mbc, meqn)
  double precision :: aux(1-mbc:maxmx+mbc, *)
  
  double precision :: rho1, rho2
  common /param/ rho1,rho2

  integer :: i
  double precision :: xCell
  
  do i=1,mx(1)
  	xCell = xlower(1) + (i-0.5d0)*dx(1)
    if (xCell .lt. 0.d0) then
    	q(i,1) = rho1
    else
    	q(i,1) = rho2
    endif

  enddo
  
  return
  end subroutine qinit
