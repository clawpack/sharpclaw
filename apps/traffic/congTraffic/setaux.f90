! ==================================================
  subroutine setaux(maxmx,mbc,mx,xlower,dx,maux,aux)
! ==================================================
!
! # Set auxiliary array for the the non-linear traffic flow equation 
! # q_t + (u_max(x)*q*(1-q))_x = 0.
!
! # The auxiliary array contains the distribution of the velocity, i.e. 
! # u_max=u_max(x)
!
!     
  implicit none
      
  integer :: maxmx, mbc, maux
  integer :: mx(1)
  double precision:: xlower(1),dx(1)
  double precision :: aux(1-mbc:mx(1)+mbc, 1)
  
  double precision :: rho1, rho2, u1, u2
  common /param/ rho1, rho2, u1, u2

  integer :: i
  double precision :: xCell
  
  !write(*,*) u1,u2
  
  do i=1-mbc,mx(1)+mbc
  	xCell = xlower(1) + (i-0.5d0)*dx(1)
  		if (xCell .lt. 0.d0) then
  			aux(i,1) = u1
        else
        	aux(i,1) = u2
        endif
        
         !write(*,*) i,aux(i,1)
  enddo
  
 

  return     
  end subroutine setaux
