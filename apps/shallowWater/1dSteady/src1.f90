! =========================================================
  subroutine src(q,dq,aux,t,dt)
! =========================================================
! src() evaluates (delta t) * dq/dt = (delta t) * psi(t)
! instantaneously (not integrate over a timestep)
! This differs from the function src() in CLAWPACK


  Use ClawParams
      
  implicit none

  double precision, intent(in) :: q(1-mbc:nx(1)+mbc, meqn)
  double precision :: aux(1-mbc:nx(1)+mbc, *)
  double precision :: dt, t
  double precision, intent(out) :: dq(1-mbc:nx(1)+mbc, meqn)

  double precision :: grav
  common /comrp/ grav

  integer :: i
  
  do i=1,nx(1)
    dq(i,1) = 0.d0
    dq(i,2) = -1.0d0*grav*q(i,1)*((aux(i+1,1)-aux(i-1,1))/(2*dx(1)))*dt !src = -g*h*(bottom)_x
  enddo
  

  return
  end subroutine src

