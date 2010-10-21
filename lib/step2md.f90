! ==========================================================
subroutine step2md(q,g,dq,aux,dt,cfl,t,rp,tfluct)
! ==========================================================

    ! Evaluate dq/dt * (Delta t)
    ! On entry, q should give the initial data for this step
    ! This semi-discrete implementation is genuinely multi-d
    ! and was adapted from code provided by Yulong Xing

    ! This value is passed into the Riemann solvers. The scaled fluctuations
    ! go into the array dq.

    ! dq(i,.) modifies Q of cell i

    ! WENO5 reconstruction with 3-pt. gauss quadrature
    ! Genuinely multi-d reconstruction

    use ClawData
    use ClawParams
    implicit none

    type(griddat) g

    external  rp,tfluct
    double precision, target, intent(in) :: q(1-mbc:nx(1)+mbc, 1-mbc:nx(2)+mbc, meqn)
    double precision, intent(inout) :: dq(1-mbc:nx(1)+mbc, 1-mbc:nx(2)+mbc, meqn)
    double precision :: qgauss(1-mbc:nx(1)+mbc, 1-mbc:nx(2)+mbc, meqn, 3)
    double precision :: q1dgauss(1-mbc:maxnx+mbc, meqn, 3)
    double precision :: qlgauss(1-mbc:maxnx+mbc, meqn, 3)
    double precision :: qrgauss(1-mbc:maxnx+mbc, meqn, 3)
    double precision :: aux(1-mbc:nx(1)+mbc, 1-mbc:nx(2)+mbc, *)
    
    double precision :: w(3,3,3), w5(5,3), wt(3)

    double precision, intent(in) :: dt,t
    double precision, intent(out) :: cfl
    integer i,j,m,ma,mw,ig,ixy
    double precision :: cfl1d, dtdx, dtdy

    ! Since this is genuinely multi-d, we can't just operate on
    ! 1D slices except at the lowest level, so flux2 has been
    ! moved up into here.

    ! More Gauss Quadrature weights
    wt(1) = 5.d0/18.d0
    wt(2) = 8.d0/18.d0
    wt(3) = 5.d0/18.d0

    cfl = 0.d0
    dtdx = dt/dx(1)
    dtdy = dt/dx(2)

    if (mcapa.eq.0) then
        ! no capa array:
        do i=1-mbc,nx(1)+mbc
            g%dtdx1d(i) = dtdx
        enddo
        do j=1-mbc,nx(2)+mbc
            g%dtdy1d(j) = dtdy
        enddo
    endif

    dq=0.d0

    ! =====================================================
    !  Evaluate A*q_x
    ! =====================================================

    ixy=2
    ! Loop over vertical slices
    do i = 1-mbc, nx(1)+mbc

        ! Store a vertical slice in q1d
        forall(j=1-mbc:nx(2)+mbc, m=1:meqn)
            g%q1d(j,m) = q(i,j,m)
        end forall

        if (mcapa.gt.0) then
            forall(j=1-mbc:nx(2)+mbc)
               g%dtdy1d(j) = dtdy / aux(i,j,mcapa)
            end forall
        endif

        if (maux .gt. 0)  then
            forall(ma=1:maux, j=1-mbc:nx(2)+mbc)
                g%aux2(j,ma) = aux(i,  j,ma)
            end forall
        endif

        ! Reconstruct horizontal cell averages along a vertical slice
        call weno_gauss(ixy,maxnx,meqn,mbc,nx(2),g%q1d,q1dgauss,g%aux2,g%aux2,w,w5)

        ! Store reconstructed vertical slice in qgauss
        forall(j=1-mbc:nx(2)+mbc, m=1:meqn, ig=1:3)
            qgauss(i,j,m,ig) = q1dgauss(j,m,ig)
        end forall
            
    enddo ! End of outer loop over x

    ! apply BCs to qgauss()
    ! call bc2(maxnx,maxmy,meqn,mbc,nx,my,xlower,ylower, &
    !           dx,dy,qgauss,maux,aux,t,dt,mthbc)

    ixy=1
    ! Loop over horizontal slices
    do j = 1-mbc, nx(2)+mbc

        ! Store a horizontal slice of qgauss in q1dgauss
        forall(i=1-mbc:nx(1)+mbc, m=1:meqn, ig=1:3)
            q1dgauss(i,m,ig) = qgauss(i,j,m,ig)
        end forall

        if (mcapa.gt.0)  then
            forall(i=1-mbc:nx(1)+mbc)
               g%dtdx1d(i) = dtdx / aux(i,j,mcapa)
            end forall
        endif

        if (maux .gt. 0)  then
            forall(ma=1:maux, i=1-mbc:nx(1)+mbc)
                  g%aux2(i,ma) = aux(i,j  ,ma)
            end forall
        endif

        ! Reconstruct left/right values at each gauss pt. on each
        ! interface along a horizontal slice
        call weno_lines(ixy,maxnx,meqn,mbc,nx,q1dgauss,&
                            qlgauss,qrgauss,g%aux2,g%aux2)

        ! solve Riemann problem at each interface gauss pt. and compute 
        ! fluctuations

        do ig=1,3

            forall (i=1-mbc:nx(1)+mbc, m=1:meqn)
                g%ql(i,m)=qlgauss(i,m,ig)
                g%qr(i,m)=qrgauss(i,m,ig)
            end forall
               
            call rpn2(ixy,maxnx,meqn,mwaves,mbc,nx(1),g%ql,g%qr,g%aux2,g%aux2,&
                        g%wave,g%s,g%amdq,g%apdq)

            ! solve internal Riemann problem in each cell and compute updates
            do i = 1-mbc,nx(1)+mbc
                do m=1,meqn
                    g%qr(i-1,m)=g%ql(i,m)
                    g%ql(i  ,m)=g%qr(i,m)
                enddo
                do m=1,maux
                    g%auxr(i-1,m)=g%aux2(i,m)
                    g%auxl(i  ,m)=g%aux2(i,m)
                enddo
            enddo

            call rpn2(ixy,maxnx,meqn,mwaves,mbc,nx(1),g%ql,g%qr,g%auxl,g%auxr,&
                        g%wave,g%s,g%amdq2,g%apdq2)

            if (mcapa.eq.0) then
                ! no capa array.  Standard flux differencing:
                forall (m=1:meqn, i=1:nx(1)+1)
                    dq(i,j,m) = dq(i,j,m) - wt(ig)*g%dtdx1d(i)*(g%apdq(i,m) + &
                                g%amdq(i+1,m) + g%apdq2(i,m) + g%amdq2(i,m))
                end forall
            else
                ! with capa array.
                forall (m=1:meqn, i=1:nx(1)+1)
                    dq(i,j,m) = dq(i,j,m) - wt(ig)*g%dtdx1d(i)*(g%apdq(i,m) + &
                        g%amdq(i+1,m) + g%apdq2(i,m) + g%amdq2(i,m))/aux(i,j,mcapa)
                end forall
            endif

        enddo ! Loop over gauss points

    enddo ! End of outer loop over y

    ! compute maximum wave speed for checking Courant number:
    ! This usually needs to be within loops above
    ! But this is more efficient and works for linear problems
    cfl1d = 0.d0
    do mw=1,mwaves
        do i=1,nx(1)+1
            ! if s>0 use dtdx1d(i) to compute CFL,
            ! if s<0 use dtdx1d(i-1) to compute CFL:
            cfl1d = dmax1(cfl1d, g%dtdx1d(i)*g%s(i,mw), &
                                -g%dtdx1d(i-1)*g%s(i,mw))
        enddo
    enddo
    cfl = dmax1(cfl,cfl1d)

! =====================================================
!  END OF A*q_x
! =====================================================

! =====================================================
!  START OF B*q_y
! =====================================================

    ixy=1
    ! Loop over horizontal slices
    do j = 1-mbc, nx(2)+mbc
        ! Store a horizontal slice in q1d
        forall(i=1-mbc:nx(1)+mbc, m=1:meqn)
            g%q1d(i,m) = q(i,j,m)
        end forall

        if (mcapa.gt.0)  then
            g%dtdx1d = dtdx / aux(:,j,mcapa)
        endif

        if (maux .gt. 0)  then
            forall(ma=1:maux, i=1-mbc:nx(1)+mbc)
                g%aux2(i,ma) = aux(i,j  ,ma)
            end forall
        endif

        ! Reconstruct vertical cell averages along a horizontal slice
        call weno_gauss(ixy,maxnx,meqn,mbc,nx(1),g%q1d,q1dgauss,&
                            g%aux2,g%aux2,w,w5)

        ! Store reconstructed horizontal slice in qgauss
        forall(i=1-mbc:nx(1)+mbc,m=1:meqn,ig=1:3)
            qgauss(i,j,m,ig) = q1dgauss(i,m,ig)
        end forall
            
    enddo !End of outer loop over y

    ! Need to apply BCs to qgauss() here
    ! call bc2(maxnx,maxmy,meqn,mbc,nx(1),my,xlower,ylower,
    !    &               dx,dy,qgauss,maux,aux,t,dt,mthbc)

    ixy=2
    ! Loop over vertical slices
    do i = 1-mbc, nx(1)+mbc
        ! Store a vertical slice of qgauss in q1dgauss
        forall(j=1-mbc:nx(2)+mbc, m=1:meqn, ig=1:3)
            q1dgauss(j,m,ig) = qgauss(i,j,m,ig)
        end forall

        if (mcapa.gt.0)  then
            g%dtdy1d = dtdy / aux(i,:,mcapa)
        endif

        if (maux .gt. 0)  then
            forall(ma=1:maux, j=1-mbc:nx(2)+mbc)
                g%aux2(j,ma) = aux(i,  j,ma)
            end forall
        endif

        ! Reconstruct left/right values at each gauss pt. on each
        ! interface along a vertical slice
        call weno_lines(ixy,maxnx,meqn,mbc,nx(2),q1dgauss,&
                            qlgauss,qrgauss,g%aux2,g%aux2)

        ! solve Riemann problem at each interface gauss pt. and compute 
        ! fluctuations

        do ig=1,3

            forall (j = 1-mbc:nx(2)+mbc, m=1:meqn)
                g%ql(j,m)=qlgauss(j,m,ig)
                g%qr(j,m)=qrgauss(j,m,ig)
            end forall
               
            call rpn2(ixy,maxnx,meqn,mwaves,mbc,nx(2),g%ql,g%qr,g%aux2,g%aux2, &
                      g%wave,g%s,g%amdq,g%apdq)

            ! solve internal Riemann problem in each cell and compute updates
            forall(j=1-mbc:nx(2)+mbc,m=1:meqn)
                g%qr(j-1,m)=g%ql(j,m)
                g%ql(j  ,m)=g%qr(j,m)
            end forall
            forall(j=1-mbc:nx(2)+mbc,m=1:maux)
                g%auxr(j-1,m)=g%aux2(j,m)
                g%auxl(j  ,m)=g%aux2(j,m)
            end forall

            call rpn2(ixy,maxnx,meqn,mwaves,mbc,nx(2),g%ql,g%qr,g%auxl,g%auxr, &
                         g%wave,g%s,g%amdq2,g%apdq2)

            if (mcapa.eq.0) then
                ! no capa array.  Standard flux differencing:
                forall(m=1:meqn, j=1:nx(2)+1)
                    dq(i,j,m) = dq(i,j,m) - wt(ig)*g%dtdy1d(j)*(g%apdq(j,m) + &
                                    g%amdq(j+1,m) + g%apdq2(j,m) + g%amdq2(j,m))
                end forall
            else
                ! with capa array.
                forall(m=1:meqn, j=1:nx(2)+1)
                    dq(i,j,m) = dq(i,j,m) - wt(ig)*g%dtdy1d(j)*(g%apdq(j,m) + &
                        g%amdq(j+1,m) + g%apdq2(j,m) + g%amdq2(j,m))/aux(i,j,mcapa)
                end forall
            endif

        enddo!Loop over gauss points

    enddo !End of outer loop over x

    ! compute maximum wave speed for checking Courant number:
    ! This usually needs to be within loops above
    ! But this is more efficient and works for linear problems
    cfl1d = 0.d0
    do 51 mw=1,mwaves
        do 51 j=1,nx(2)+1
            ! if s>0 use dtdy1d(i) to compute CFL,
            ! if s<0 use dtdy1d(i-1) to compute CFL:
            cfl1d = dmax1(cfl1d, g%dtdy1d(i)*g%s(i,mw), &
                                -g%dtdy1d(i-1)*g%s(i,mw))
   51    continue
    cfl = dmax1(cfl,cfl1d)

! =====================================================
!  END OF B*q_y
! =====================================================

end subroutine step2md
