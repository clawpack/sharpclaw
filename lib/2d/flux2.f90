! ==========================================================
subroutine flux2(q,g,dq,aux,dt,cfl,t,rp,tfluct)
! ==========================================================

    ! Evaluate (delta t) *dq/dt
    ! On entry, q should give the initial data for this step
    ! dimension-by-dimension version
    ! (No genuinely multi-d terms)

    ! g%dq1d is used to return increments to q from flux1
    ! See the flux1 documentation for more information.

    use ClawData
    use ClawParams
    implicit none

    type(griddat) g

    external rp,tfluct
    double precision, target, intent(in) :: q(1-mbc:nx(1)+mbc, 1-mbc:nx(2)+mbc, meqn)
    double precision, intent(inout) :: dq(1-mbc:nx(1)+mbc, 1-mbc:nx(2)+mbc, meqn)
    double precision, target, intent(in) :: aux(1-mbc:nx(1)+mbc, 1-mbc:nx(2)+mbc,maux)
    double precision, intent(in) :: dt,t
    double precision, intent(out) :: cfl
    integer :: i,j,m
    double precision :: cfl1d
    double precision, pointer :: auxp(:,:),q1dp(:,:)

    cfl = 0.d0
    dq=0.d0

    ! perform x-sweeps
    ! ==================

    do j = 0,nx(2)+1

        ! copy data along a slice into 1d arrays:
        q1dp => q(:,j,:)
        if (maux .gt. 0)  then
            auxp => aux(:,j,:)
        endif

        ! compute modification dq1d along this slice:
        call flux1(q1dp,g,g%dq1d,auxp,dt,cfl1d,t,rp,tfluct,1)
        cfl = dmax1(cfl,cfl1d)

        if (mcapa.eq.0) then
            ! no capa array.  Standard flux differencing:
            forall(i=1:nx(1), m=1:meqn)
                dq(i,j,m) = dq(i,j,m)+g%dq1d(i,m)
            end forall
        else
            ! with capa array.  Which is correct?
            forall(i=1:nx(1), m=1:meqn)
                dq(i,j,m) = dq(i,j,m)+g%dq1d(i,m)
            end forall
            ! dq(i,j,m) = dq(i,j,m)+g%dq1d(i,m)/aux(i,j,mcapa)
        endif
    enddo !end x sweeps


    ! perform y sweeps
    ! ==================

    do i = 0, nx(1)+1
        ! copy data along a slice into 1d arrays:
        q1dp => q(i,:,:)
        forall(j=1-mbc:nx(2)+mbc, m=1:meqn)
            g%q1d(j,m) = q(i,j,m)
        end forall

        if (maux .gt. 0)  then
            auxp => aux(i,:,:)
        endif

        call flux1(q1dp,g,g%dq1d,auxp,dt,cfl1d,t,rp,tfluct,2)
        cfl = dmax1(cfl,cfl1d)

        if (mcapa.eq.0) then
            ! no capa array.  Standard flux differencing:
            forall(j=1:nx(2),m=1:meqn)
                dq(i,j,m) = dq(i,j,m)+g%dq1d(j,m)
            end forall
        else
            ! with capa array.  Which is correct?
            ! dq(i,j,m) = dq(i,j,m)+g%dq1d(j,m)/aux(i,j,mcapa)
            dq(i,:,:) = dq(i,:,:)+g%dq1d
        endif
    enddo !end y sweeps

end subroutine flux2
