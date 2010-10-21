! =========================================================
      subroutine src(q,dq,aux,t,dt)
! =========================================================
      use ClawParams
      implicit none

      double precision, intent(in) :: q,aux,dt,t
      integer :: i,m
!
!     Use this routine if there is no source term
!     Simply zeros the dq vector

!     Note that in general src() should evaluate 
!     (delta t) * dq/dt = (delta t) * psi(t)
!     instantaneously (not integrate over a timestep)
!     This differs from the function src() in CLAWPACK

      double precision :: dq(1-mbc:nx(1)+mbc, meqn)


      forall (i=1:nx(1), m=1:meqn)
        dq(i,m) = 0.d0
      end forall

      end subroutine src
