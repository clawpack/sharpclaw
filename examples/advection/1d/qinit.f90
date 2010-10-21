! =========================================================
subroutine qinit(maxmx,meqn,mbc,mx,xlower,dx,q,maux,aux)
! =========================================================

    ! Set initial conditions for q.


    implicit double precision (a-h,o-z)
    integer, parameter :: ndim=1
    integer :: mx(ndim)
    double precision :: xlower(ndim), dx(ndim)
    double precision, intent(out) :: q(1-mbc:mx(1)+mbc, meqn)
    dimension aux(1-mbc:mx(1)+mbc, *)
    common /comic/ beta,ic


    pi2 = 8.d0*datan(1.d0)  !# = 2 * pi
    do i=1,mx(1)
        xcell = xlower(1) + (i-0.5d0)*dx(1)

        ! gaussian:
        q(i,1) = dexp(-beta * (xcell-0.3d0)**2)

    enddo

end subroutine qinit
