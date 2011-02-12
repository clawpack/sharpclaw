!========================================================
 subroutine qinit(maxmx,meqn,mbc,mx,xlower,dx,q,maux,aux)
!========================================================

! # Set initial conditions for q.


  implicit none
  
  integer :: maxmx, meqn, mbc, maux
  integer :: mx(1)
  double precision :: dx(1)
  double precision :: xlower(1)
  double precision :: q(1-mbc:mx(1)+mbc, meqn)
  double precision :: aux(1-mbc:mx(1)+mbc, *)
  
  double precision :: u1,u2,rho1,rho2
  common /comtraf/ u1,u2,rho1,rho2
  
  integer :: i
  double precision :: xcell
  
  write(*,*) u1,u2,rho1,rho2
      
  do i=1,mx(1)
  	xcell = xlower(1) + (i-0.5d0)*dx(1)
    if (xcell .lt. 0.0d0) then
    	q(i,1) = rho1
	else
	    q(i,1) = rho2
	endif
  enddo
  
  !#.and. xcell .gt. -20.0d0
     
     
  end subroutine qinit