! ===================================================================
subroutine flux1(q1d,g,dq1d,aux,dt,cfl,t,rp,tfluct,ixy)
! ===================================================================
!
!     # Evaluate (delta t) * dq(t)/dt
!
!     SharpClaw
!     Author: David Ketcheson
!
!     amdq, apdq, amdq2, apdq2, wave, and s are used locally:
!
!     amdq(1-mbc:mx+mbc, meqn) = left-going flux-differences
!     apdq(1-mbc:mx+mbc, meqn) = right-going flux-differences
!        e.g. amdq(i,m) = m'th component of A^- \Delta q from i'th Riemann
!                         problem (between cells i-1 and i).
!
!     wave(1-mbc:mx+mbc, meqn, mwaves) = waves from solution of
!                                           Riemann problems,
!            wave(i,m,mw) = mth component of jump in q across
!                           wave in family mw in Riemann problem between
!                           states i-1 and i.
!
!     s(1-mbc:mx+mbc, mwaves) = wave speeds,
!            s(i,mw) = speed of wave in family mw in Riemann problem between
!                      states i-1 and i.
!
!     t is the time at which we want to evaluate dq/dt, which may not
!      be the current simulation time
! ===================================================================

    USE ClawData
    use ClawParams
    USE reconstruct

    implicit double precision (a-h,o-z)

    type(griddat) g

    double precision :: q1d(1-mbc:nx(ixy)+mbc, meqn)
    dimension    dq1d(1-mbc:maxnx+mbc, meqn)
    double precision, target :: aux(1-mbc:nx(ixy)+mbc, maux)
    double precision, intent(out) :: cfl
    integer, intent(in) :: ixy
    integer mx
    double precision, pointer :: auxpl(:,:), auxpr(:,:)

    mx=nx(ixy)

    if (mcapa.gt.0) then
        g%dtdx = dt / (dx(ixy)*aux(:,mcapa))
    else
        g%dtdx = dt/dx(ixy)
    endif
    if (ndim.gt.1) dq1d=0.d0

    select case(lim_type)
        ! Non-limited reconstruction of components of q (simplest approach)
        case(0)
        select case(char_decomp)
            case(0)
                call q2qlqr_poly(q1d,g%ql,g%qr,mx)
            case(1)
                ! wave-based unlimited reconstruction
                call rp(ixy,maxnx,meqn,mwaves,mbc,mx,&
                        q1d,q1d,aux,aux,g%wave,g%s,g%amdq,g%apdq)
                call q2qlqr_poly_wave(q1d,g%ql,g%qr,g%wave,g%s,mx)
        end select
        case(1)
        select case(char_decomp)
            case(0)
                ! Fill in TVD reconstruction w/o char. decomp. here
                call tvd2(q1d,g%ql,g%qr,mthlim)
            case(1)
                ! wave-based second order reconstruction
                call rp(ixy,maxnx,meqn,mwaves,mbc,mx,&
                        q1d,q1d,aux,aux,g%wave,g%s,g%amdq,g%apdq)
                call tvd2_wave(q1d,g%ql,g%qr,g%wave,g%s,mthlim)
            case(2)
                ! characteristic-wise second order reconstruction
                call evec(mx,meqn,mbc,mx,q1d,aux,aux,g%evl,g%evr)
                call tvd2_char(q1d,g%ql,g%qr,mthlim,g%evl,g%evr)
        end select
        case(2)
        select case (char_decomp)
            case (0)
                ! no characteristic decomposition
                call weno5(q1d,g%ql,g%qr)
            case (1)
                ! wave-based reconstruction
                call rp(ixy,maxnx,meqn,mwaves,mbc,mx,&
                        q1d,q1d,aux,aux,g%wave,g%s,g%amdq,g%apdq)
                call weno5_wave(q1d,g%ql,g%qr,g%wave)
            case (2)
                ! characteristic-wise reconstruction
                call evec(mx,meqn,mbc,mx,q1d,aux,aux,g%evl,g%evr)
                call weno5_char(q1d,g%ql,g%qr,g%evl,g%evr)
            case (3)
                ! transmission-based reconstruction
                call evec(mx,meqn,mbc,mx,q1d,aux,aux,g%evl,g%evr)
                call weno5_trans(q1d,g%ql,g%qr,g%evl,g%evr)
            case default
                write(*,*) 'ERROR: Unrecognized characteristic decomposition &
                            option'
                write(*,*) 'You should set 0<=char_decomp<=3'
                stop
        end select
      end select


    ! solve Riemann problem at each interface 
    ! -----------------------------------------
    call rp(ixy,maxnx,meqn,mwaves,mbc,mx,g%ql,g%qr,aux,aux, &
              g%wave,g%s,g%amdq,g%apdq)

    ! compute maximum wave speed:
    cfl = 0.d0
    do mw=1,mwaves
        do i=1,mx+1
            ! if s>0 use dtdx(i) to compute CFL,
            ! if s<0 use dtdx(i-1) to compute CFL:
            cfl = dmax1(cfl, g%dtdx(i)*g%s(i,mw), -g%dtdx(i-1)*g%s(i,mw))
        enddo
    enddo

    ! Find total fluctuation within each cell
    if (tfluct_solver==1) then
        ! tfluct should be a special solver that uses the parameters aux(i)
        ! to solve a Riemann problem with left state ql(i)
        ! and right state qr(i), and returns a total fluctuation in amdq2
        ! NOTE that here amdq2 is really a total fluctuation (should be
        ! called adq); we do it this way just to avoid declaring more storage
        call tfluct(ixy,maxnx,meqn,mwaves,mbc,mx,g%ql,g%qr, &
                     aux,aux,g%s,g%amdq2)

        ! Modify q using fluctuations:
        ! Note this may not correspond to a conservative flux-differencing
        ! for equations not in conservation form.  It is conservative if
        ! adq = f(qr(i)) - f(ql(i)).

        forall (i=1:mx, m=1:meqn)
            dq1d(i,m) = dq1d(i,m) - g%dtdx(i)*(g%apdq(i,m) + &
                            g%amdq2(i,m) + g%amdq(i+1,m))
        end forall

    else
        ! Or we can just swap things around and use the usual Riemann solver
        ! This may be more convenient, but is less efficient.
        do i = 1-mbc+1,mx+mbc
            do m = 1, meqn
                g%qr(i-1,m) = g%ql(i,m)
                g%ql(i  ,m) = g%qr(i,m)
            enddo
        enddo
        auxpr => aux(1-mbc+1:,:)
        auxpl => aux
        call rp(ixy,maxnx,meqn,mwaves,mbc,mx,g%ql,g%qr, &
                 auxpr,auxpl,g%wave,g%s,g%amdq2,g%apdq2)

        forall(i=1:mx, m=1:meqn)
            dq1d(i,m) = dq1d(i,m)-g%dtdx(i)*(g%amdq(i+1,m)+ &
                        g%apdq(i,m)+g%amdq2(i,m)+g%apdq2(i,m))
        end forall
    endif

end subroutine flux1
