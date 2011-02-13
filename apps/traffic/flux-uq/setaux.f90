!============================================
 subroutine setaux(maxmx,mbc,mx,xlower,dx,maux,aux)
!============================================

! # Set auxiliary arrays 
! # aux(i,1) = velocity u_i in i'th cell for advection
     
  implicit none
      
  integer :: mx(1)
  integer :: maxmx, mbc, maux
  double precision :: dx(1)
  double precision :: xlower(1)
  double precision :: aux(1-mbc:mx(1)+mbc, maux)
      
  double precision :: u1,u2,rho1,rho2
  common /comtraf/ u1,u2,rho1,rho2
  
  integer :: i
  double precision :: pi2,xedge
  
  pi2 = 8.d0*datan(1.d0)

  do i=1-mbc,mx(1)+mbc
  	xedge = xlower(1) + (i-0.5d0)*dx(1)

	if (xedge .lt. 0.d0) then
		aux(i,1) = u1
	else
	    aux(i,1) = u2
	endif
!	  	write(*,*) dx(1),i,xedge,aux(i,1)

  enddo

  end subroutine setaux
