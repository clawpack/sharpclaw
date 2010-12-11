! =========================================================
      subroutine src(q,dq,aux,t,dt)
! =========================================================


      Use ClawParams
      
      implicit none

      double precision, intent(in) :: q(1-mbc:nx(1)+mbc, meqn)
      double precision :: aux,dt,t
      integer :: i

!     src() evaluates (delta t) * dq/dt = (delta t) * psi(t)
!     instantaneously (not integrate over a timestep)
!     This differs from the function src() in CLAWPACK

      double precision, intent(out) :: dq(1-mbc:nx(1)+mbc, meqn)
      double precision :: rad_pos


      !tmp = dx
      
      do i=1,nx(1)
        rad_pos = xlower(ndim) + (i-0.5d0)*dx(ndim)
        dq(i,1) = -1.d0*dt*q(i,1)*q(i,2)/rad_pos
        dq(i,2) = -1.0d0*dt*q(i,1)*q(i,2)**2/rad_pos
      enddo

      end subroutine src

