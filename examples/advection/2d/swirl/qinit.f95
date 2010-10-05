! =====================================================
subroutine qinit(maxmx,meqn,mbc,nx,xlower,dx,q,maux,aux)
! =====================================================

    ! Set initial conditions for q.
    ! Sample scalar equation with data that is piecewise constant with
    ! q = 1.0  if  0.1 < x < 0.6   and   0.1 < y < 0.6
    !     0.1  otherwise

    implicit none
    integer :: i,j
    double precision :: xi,yj
    integer, intent(in) :: maxmx,meqn,mbc,maux
    integer,intent(in) :: nx(2)
    double precision, intent(in) :: aux(1-mbc:nx(1)+mbc, 1-mbc:nx(2)+mbc, maux)
    double precision, intent(in) :: xlower(2),dx(2)
    double precision :: q(1-mbc:nx(1)+mbc, 1-mbc:nx(2)+mbc, meqn)

    do i=1,nx(1)
        xi = xlower(1) + (i-0.5d0)*dx(1)
        do j=1,nx(2)
            yj = xlower(2) + (j-0.5d0)*dx(2)
            if (xi.lt.0.5d0) then
            ! if (xi.lt.0.5d0 .and. xi.gt.0.1d0 .and. yj.gt.0.1d0 .and. &
            ! yj.lt.0.3d0) then
                q(i,j,1) = 1.d0
            else
                q(i,j,1) = 0.d0
            endif
        enddo
    enddo
end subroutine qinit
