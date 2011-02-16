! =========================================================
      subroutine src(q,dq,aux,t,dt)
! =========================================================
!
!     Use this routine if there is no source term
!     Simply sets the dq vector to zero

!     Note that in general src() should evaluate 
!     (delta t) * dq/dt = (delta t) * psi(t)
!     instantaneously (not integrate over a timestep)
!     This differs from the function src() in CLAWPACK

      use ClawParams
      implicit none

      double precision, intent(in) :: q, aux, dt, t
      double precision :: dq(1-mbc:nx(1)+mbc, 1-mbc:nx(2)+mbc, meqn)
      integer :: i, j, m

      forall(i=1:nx(1), j=1:nx(2), m=1:meqn)
      	dq(i,j,m) = 0.d0
      end forall
            

      end subroutine src
