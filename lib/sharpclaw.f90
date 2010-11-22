! ======================================================================
subroutine sharpclaw(tstart,tend,cfl_maxused,cfl_last,dt_minused, &
                dt_maxused,dt_last,steps_taken,info,dt_initial,dt_max, &
                dt_variable,time_integrator,verbosity, &
                cfl_max,cfl_desired,max_steps,r1,r2,bc,rp,tfluct,src,b4step)
! ======================================================================

    USE ClawData
    USE ClawParams
    use Global
    implicit none

    external bc,rp,tfluct,src,b4step

    type(griddat) gdat
    type(rkreg) r1,r2
    
    double precision :: gp
    double precision :: tstart,tend,cflmax,cfl
    double precision :: dtmax,dtmin,dt,t,told
    integer :: info,maxn,n,stat
    
    double precision :: cfl_maxused, cfl_last
    double precision :: dt_minused, dt_maxused, dt_last
    double precision :: dt_initial,dt_max
    integer :: steps_taken
    double precision :: cfl_max, cfl_desired
    integer :: max_steps
    integer :: dt_variable,time_integrator,verbosity

    logical :: retake_step
    double precision :: dtcom,tcom
    common /comxt/ dtcom,tcom
    
    info = 0
    t = tstart
    maxn = max_steps
    dt = dt_initial   !# initial dt
    cflmax = 0.d0
    dtmin = dt
    dtmax = dt
    steps_taken = 0

    if (dt_variable == 0) then
        ! fixed size time steps.  Compute the number of steps:
        if (tend < tstart) then
            ! single step mode
            maxn = 1
        else
            maxn = (tend - tstart + 1d-10) / dt
            if (dabs(maxn*dt - (tend-tstart)) > &
                1d-5*(tend-tstart)) then
                ! dt doesn't divide time interval integer number of times
                info = 2
                return
            endif
        endif
    endif


    if ((dt_variable==1) .and. (cfl_desired>cfl_max)) then
        info = 3
        return
    endif

    maxnx=maxval(nx)

    allocate( gdat%s(1-mbc:maxnx+mbc, mwaves),STAT=stat )
    if (stat>0) write(6,*) '*** ERROR ***'
    allocate( gdat%wave(1-mbc:maxnx+mbc, meqn, mwaves),STAT=stat )
    if (stat>0) write(6,*) '*** ERROR ***'
    allocate( gdat%ql(1-mbc:maxnx+mbc, meqn),STAT=stat )
    if (stat>0) write(6,*) '*** ERROR ***'
    allocate( gdat%qr(1-mbc:maxnx+mbc, meqn),STAT=stat )
    if (stat>0) write(6,*) '*** ERROR ***'
    allocate( gdat%dtdx(1-mbc:maxnx+mbc),STAT=stat )
    if (stat>0) write(6,*) '*** ERROR ***'
    allocate( gdat%amdq(1-mbc:maxnx+mbc, meqn),STAT=stat )
    if (stat>0) write(6,*) '*** ERROR ***'
    allocate( gdat%apdq(1-mbc:maxnx+mbc, meqn),STAT=stat )
    if (stat>0) write(6,*) '*** ERROR ***'
    allocate( gdat%amdq2(1-mbc:maxnx+mbc, meqn),STAT=stat )
    if (stat>0) write(6,*) '*** ERROR ***'
    allocate( gdat%apdq2(1-mbc:maxnx+mbc, meqn),STAT=stat )
    if (stat>0) write(6,*) '*** ERROR ***'
    allocate( gdat%evl(1-mbc:maxnx+mbc,meqn,meqn),STAT=stat )
    if (stat>0) write(6,*) '*** ERROR ***'
    allocate( gdat%evr(1-mbc:maxnx+mbc,meqn,meqn),STAT=stat )
    if (stat>0) write(6,*) '*** ERROR ***'
    if (tfluct_solver==0) then
        allocate( gdat%auxl(1-mbc:maxnx+mbc,maux) )
        allocate( gdat%auxr(1-mbc:maxnx+mbc,maux) )
    endif

    ! Allocate extra arrays needed in 2D
    if (ndim>1) then
        allocate(gdat%q1d(1-mbc:maxnx+mbc,meqn))
        allocate(gdat%dq1d(1-mbc:maxnx+mbc,meqn))
        allocate(gdat%aux2(1-mbc:maxnx+mbc,maux))
        if (multid_recon==1) then
            allocate(gdat%dtdx1d(1-mbc:nx(1)+mbc))
            allocate(gdat%dtdy1d(1-mbc:nx(2)+mbc))

            allocate(gdat%qgauss(1-mbc:nx(1)+mbc, 1-mbc:nx(2)+mbc, meqn,3))
            allocate(gdat%q1dgauss(1-mbc:maxnx+mbc,meqn,3))
            allocate(gdat%qlgauss(1-mbc:maxnx+mbc,meqn,3))
            allocate(gdat%qrgauss(1-mbc:maxnx+mbc,meqn,3))


        endif
    endif

    retake_step = .FALSE.

    ! ======================================================================
    ! Main loop
    ! ======================================================================
    do n=1,maxn
        if (.not.retake_step) then
            told = t   !# time at beginning of time step.

            ! adjust dt to hit tend exactly if we're near end of computation
            ! (unless tend < tstart, which is a flag to take only a single step)
            if (told+dt>tend .and. tstart<tend) dt = tend - told

        !HACK This goes below bc really
            ! call user-supplied routine which might set aux arrays
            ! for this time step, for example.
        call b4step(maxnx,mbc,nx,meqn,q,xlower,dx,told,dt,maux,aux)
        !END HACK

            if (dt_variable==1) then
                ! save old q in case we need to retake step with smaller dt:
                qold=q
            endif
        endif
        retake_step = .FALSE.

        t = told + dt    !# time at end of step


        ! store dt and t in the common block comxt in case they are needed
        ! in the Riemann solvers (for variable coefficients)
         tcom = told
         dtcom = dt

        ! ------------------------------------------------------------------
        ! main steps in algorithm:
        ! ------------------------------------------------------------------
        ! extend data from grid to bordering boundary cells:
        select case(ndim)
            case(1)
                call bc(maxnx,meqn,mbc,nx,xlower,dx,q,maux,aux,told,dt,mthbc)
            case(2)
        end select


        ! take a step on the conservation law, including source terms
        select case(time_integrator)
            case (1) ! SSPRK22 time stepping (2-stage, 2nd order)
                call ts_ssp22(gdat,r1,dt,cfl,bc,src,rp,tfluct,told) 

            case (2) ! SSPRK42 time stepping (4-stage, 2nd order)
                call ts_ssp42(gdat,r1,dt,cfl,bc,src,rp,tfluct,told)

            case (3) ! SSPRK33 time stepping (3-stage, 3rd order)
                call ts_ssp33(gdat,r1,dt,cfl,bc,src,rp,tfluct,told)

            case (4) ! SSPRK10,4 stepping      
                call ts_ssp104(gdat,r1,r2,dt,cfl,bc,src,rp,tfluct,told)
            case default
                write(*,*) 'ERROR: Unknown time stepping option'
                write(*,*) 'time_integrator should be between 1 and 4'
                stop
        end select


        if (verbosity == 1) write(6,601) n,cfl,dt,t
  601   format('SClaw... Step',i4, &
                         '   Courant number =',f6.3,'  dt =',d12.4, &
                         '  t =',d12.4)

        !  choose new time step if variable time step
        if (dt_variable == 1) then
            if (cfl > 0.d0) then
                dt = dmin1(dt_max, dt * cfl_desired/cfl)
                dtmin = dmin1(dt,dtmin)
                dtmax = dmax1(dt,dtmax)
            else
                dt = dt_max
            endif
        endif

        ! check to see if the Courant number was too large:
        if (cfl .le. cfl_max) then
            ! accept this step
            cflmax = dmax1(cfl,cflmax)
        else
            ! reject this step
            t = told
            q=qold

            if (verbosity == 1) then
                write(6,602) 
  602           format('SClaw rejecting step... ', &
                            'Courant number too large')
            endif
            if (dt_variable==1) then
                ! if variable dt, go back and take a smaller step
                retake_step = .TRUE.
            else
                ! if fixed dt, give up and return
                cflmax = dmax1(cfl,cflmax)
                return
            endif
        endif

        !see if we are done:
        steps_taken = steps_taken + 1
        if (t >= tend) exit

    enddo

    ! too many timesteps?
    if (dt_variable==1 .and. t<tend .and. steps_taken == maxn) info = 11

    ! Courant number too large with fixed dt?
    if (dt_variable==0 .and. cflmax > cfl_max) info = 12

    tend = t
    cfl_maxused = cflmax
    cfl_last = cfl
    dt_minused = dtmin
    dt_maxused = dtmax
    dt_last = dt

    deallocate(gdat%s,gdat%wave,gdat%ql,gdat%qr,gdat%dtdx)
    deallocate(gdat%amdq,gdat%apdq,gdat%amdq2,gdat%apdq2,gdat%evl,gdat%evr)
    if (tfluct_solver==0) then
        deallocate(gdat%auxl,gdat%auxr)
    endif

end subroutine sharpclaw


! ==============================================================
subroutine ts_ssp22(g,r1,dt,cfl,bc,src,rp,tfluct,told)
! ==============================================================
    ! SSP22 time stepping (2nd order, 2 stages)
    USE ClawData
    use Global
    implicit double precision (a-h,o-z)
    type(griddat) g
    type(rkreg) r1
    external bc,rp,tfluct,src

    ! Constants for RK time stepping
    c2 =1.d0

    ! Store (dt)*L(q^n) = (dt)*d/dt(q^n(t^n) in r1
    call step(q,g,r1%qrk,aux,dt,cfl,told,rp,src,tfluct)

    ! q^(1) = q^n + dt*L(q^n) => q
    q= q+r1%qrk

    ! extend data from grid to bordering boundary cells:
    ! t^(1) = t^n + c2*dt
    that = told+c2*dt

    ! Store (dt)*L(q^1) = (dt)*d/dt(q^(1)(t^(1)) in r1
    call step(q,g,r1%qrk,aux,dt,cfl,that,rp,src,tfluct)

    ! q = 1/2*u^n + 1/2*u^(1) + 1/2*dt*L(q^1)
    q= 0.5d0*(qold + q + r1%qrk)

end subroutine ts_ssp22


! ==============================================================
subroutine ts_ssp42(g,r1,dt,cfl,bc,src,rp,tfluct,told)
! ==============================================================
    ! SSP42 time stepping (2nd order, 4 stages)
    USE ClawData                   
    use Global
    implicit double precision (a-h,o-z)
    type(griddat) g
    type(rkreg) r1
    external bc,rp,tfluct,src

    ! Constants for RK time stepping
    third=1.d0/3.d0
    a21=third
    c2 = third
    c3 = 2.d0*third
    c4 = 1.d0
    b4= 0.25d0
    
    ! Store L(q^n) = d/dt(q^n(t^n) in qrk1
    ! q^(1) = q^n + dt*L(q^n)

    call step(q,g,r1%qrk,aux,dt,cfl,told,rp,src,tfluct)

    ! q^(1) = q^n + 1/3*dt*L(q^n) => q
    q= q+a21*r1%qrk

    ! t^(1) = t^n + c2*dt
    that = told+c2*dt

    ! Store L(q^1) = d/dt(q^(1)(t^(1)) in qrk1
    call step(q,g,r1%qrk,aux,dt,cfl,that,rp,src,tfluct)

    ! q^(2) = q^(1) + 1/3*dt*L(q^(1)) => q
    q= q + third*r1%qrk

    ! t^(2) = t^n + c3*dt
    that = told+c3*dt

    ! Store L(q^2) = d/dt(q^(2)(t^(2)) in qrk1
    call step(q,g,r1%qrk,aux,dt,cfl,that,rp,src,tfluct)

    ! q^(3) = q^(2) + 1/3*dt*L(q^(2)) => q
    q= q + third*r1%qrk

    ! t^(3) = t^n + c4*dt
    that = told+c4*dt

    ! Store L(q^(3)) = d/dt(q^(3)(t^(3)) in qrk1
    call step(q,g,r1%qrk,aux,dt,cfl,that,rp,src,tfluct)

    ! q^n+1 = 1/4*q^n + 3/4*q^(3) + 1/4*dt*L(q^(3)) => q
    q=0.25d0*qold + 0.75d0*q+ b4*r1%qrk

end subroutine ts_ssp42


! ==============================================================
subroutine ts_ssp33(g,r1,dt,cfl,bc,src,rp,tfluct,told)
! ==============================================================
    ! SSP33 time stepping (3rd order, 3 stages)
    USE ClawData                   
    use Global
    implicit double precision (a-h,o-z)
    external bc,rp,tfluct,src
    type(griddat) g
    type(rkreg) r1

    ! Constants for RK time stepping
    c2 = 1.d0
    c3 = 0.5d0
    r1%qrk=0.d0

    ! Store L(q^n) = d/dt(q^n(t^n) in qrk1
    call step(q,g,r1%qrk,aux,dt,cfl,told,rp,src,tfluct)

    ! q^(1) = q^n + dt*L(q^n) => q
    q= q+r1%qrk

    ! t^(1) = t^n + c2*dt
    that = told+c2*dt

    ! Store L(q^1) = d/dt(q^(1)(t^(1)) in qrk1
    call step(q,g,r1%qrk,aux,dt,cfl,that,rp,src,tfluct)

    ! q^(2) = 3/4*q^n + 1/4*(q^(1) + dt*L(q^(1))) => q
    q= 0.75d0*qold +0.25d0*(q + r1%qrk)

    ! t^(2) = t^n + c3*dt
    that = told+c3*dt

    ! Store L(q^2) = d/dt(q^(2)(t^(2)) in qrk1
    call step(q,g,r1%qrk,aux,dt,cfl,that,rp,src,tfluct)

    ! q^n+1 = 1/3*q^n + 2/3*(q^(2) + dt*L(q^(2))) => q
    q=(qold + 2.d0*(q+ r1%qrk))/3.d0

end subroutine ts_ssp33


! ==============================================================
subroutine ts_ssp104(g,r1,r2,dt,cfl,bc,src,rp,tfluct,told)
! ==============================================================
    ! SSPRK10,4 stepping
    ! c=[0 1/6 1/3 1/2 2/3 1/3 1/2 2/3 5/6 1]
    ! c=[0 1/6 2/6 3/6 4/6 2/6 3/6 4/6 5/6 1]
    USE ClawData                   
    use Global
    implicit double precision (a-h,o-z)
    type(griddat) g
    type(rkreg) r1,r2
    external bc,rp,tfluct,src

    ! Stage 2
    ! Store (dt)*L(q^n) = (dt)*d/dt(q^n(t^n) in qrk1
    call step(q,g,r1%qrk,aux,dt,cfl,told,rp,src,tfluct)

    q= q+1.d0/6.d0*r1%qrk
    that=told+dt/6.d0

    ! Stage 3
    ! Store (dt)*L(q^1) = (dt)*d/dt(q^1(t^1) in qrk1
    call step(q,g,r1%qrk,aux,dt,cfl,that,rp,src,tfluct)

    q= q+1.d0/6.d0*r1%qrk
    that=told+dt/3.d0

    ! Stage 4
    ! Store (dt)*L(q^2) = (dt)*d/dt(q^2(t^2) in qrk1
    call step(q,g,r1%qrk,aux,dt,cfl,that,rp,src,tfluct)

    q= q+1.d0/6.d0*r1%qrk
    that=told+0.5d0*dt

    ! Stage 5
    ! Store (dt)*L(q^3) = (dt)*d/dt(q^3(t^3) in qrk1
    call step(q,g,r1%qrk,aux,dt,cfl,that,rp,src,tfluct)

    q= q+1.d0/6.d0*r1%qrk
    that=told+2.d0/3.d0*dt

    ! Stage 6
    ! Store (dt)*L(q^4) = (dt)*d/dt(q^4(t^4) in qrk1
    call step(q,g,r1%qrk,aux,dt,cfl,that,rp,src,tfluct)

    r2%qrk= q+1.d0/6.d0*r1%qrk
    r1%qrk= 0.04d0*qold+0.36d0*r2%qrk
    q= 0.6d0*qold+0.4d0*r2%qrk

    that=told+1.d0/3.d0*dt

    ! Stage 7
    ! Store (dt)*L(q^5) = (dt)*d/dt(q^5(t^5) in qrk2
    call step(q,g,r2%qrk,aux,dt,cfl,that,rp,src,tfluct)

    q= q+1.d0/6.d0*r2%qrk
    that=told+0.5d0*dt

    ! Stage 8
    ! Store (dt)*L(q^6) = (dt)*d/dt(q^6(t^6) in qrk2
    call step(q,g,r2%qrk,aux,dt,cfl,that,rp,src,tfluct)

    q= q+1.d0/6.d0*r2%qrk
    that=told+2.d0/3.d0*dt

    ! Stage 9
    ! Store (dt)*L(q^7) = (dt)*d/dt(q^7(t^7) in qrk2
    call step(q,g,r2%qrk,aux,dt,cfl,that,rp,src,tfluct)

    q= q+1.d0/6.d0*r2%qrk
    that=told+5.d0/6.d0*dt

    ! Stage 10
    ! Store (dt)*L(q^8) = (dt)*d/dt(q^8(t^8) in qrk2
    call step(q,g,r2%qrk,aux,dt,cfl,that,rp,src,tfluct)

    q= q+1.d0/6.d0*r2%qrk
    that=told+dt

    ! u^n+1
    ! Store (dt)*L(q^9) = (dt)*d/dt(q^8(t^8) in qrk2
    call step(q,g,r2%qrk,aux,dt,cfl,that,rp,src,tfluct)

    q= r1%qrk+0.6d0*q+0.1d0*r2%qrk

end subroutine ts_ssp104
