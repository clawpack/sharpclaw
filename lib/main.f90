! ===================================================================
!  Program: 	SharpClaw
!  File:   	main.f95
! ===================================================================

!     Documentation is available at
!                 http://www.amath.washington.edu/~claw/extensions/wenoclaw/
!
!     Author: David I. Ketcheson
!     Version of October, 2008 --  SharpClaw Version 0.3
!

program sharpclaw_main
    
    USE ClawData
    USE ClawParams
    use Global
    use reconstruct
    implicit none
    external bc,tfluct,rp,src,b4step
    
    integer :: mx,stat,i,n,nstepout,mw,nstop,iframe,info,md
    integer :: nout,outstyle
    double precision :: tfinal
    
    double precision :: t0,tstart,tend,dtout
    double precision, allocatable, dimension (:) :: tout
    
    character*16 fname

    double precision :: cfl_maxused, cfl_last
    double precision :: dt_minused, dt_maxused, dt_last
    double precision :: dt_initial,dt_max
    integer :: steps_taken
    double precision :: cfl_max, cfl_desired
    integer :: max_steps
    integer :: dt_variable,time_integrator,verbosity

    type(rkreg) r1,r2

    open(55,file='sharpclaw.data',status='old',form='formatted')

    ! Read the input in standard form from claw1ez.data:
    ! For a description of input parameters see the documentation at
    !    http://www.amath.washington.edu/~claw/extensions/wenoclaw

    fname = 'sharpclaw.data'
    call opendatafile(55,fname)

    read(55,*) ndim     ! Dimension of the problem (1 or 2)
    allocate(nx(ndim))
    do md=1,ndim
      read(55,*) nx(md)
    end do

    maxnx=maxval(nx)

    ! i/o variables
    read(55,*) nout
    read(55,*) outstyle
    select case(outstyle)
        case(1)
            read(55,*) tfinal
            nstepout = 1
        case(2)
            allocate(tout(nout))
            read(55,*) (tout(i), i=1,nout)
            nstepout = 1
        case(3)
          read(55,*) nstepout, nstop
          nout = nstop
    end select

    ! timestepping variables
    read(55,*) dt_initial
    read(55,*) dt_max
    read(55,*) cfl_max
    read(55,*) cfl_desired
    read(55,*) max_steps

    ! input parameters for clawpack routines
    read(55,*) dt_variable
    read(55,*) time_integrator
    read(55,*) verbosity
    read(55,*) src_term
    read(55,*) mcapa
    read(55,*) maux 
    read(55,*) tfluct_solver
    read(55,*) char_decomp

    read(55,*) meqn
    read(55,*) mwaves
    read(55,*) lim_type
    allocate(mthlim(mwaves))
    read(55,*) (mthlim(mw), mw=1,mwaves)

    ! physical domain:
    read(55,*) t0
    allocate(xlower(ndim))
    allocate(xupper(ndim))
    do md=1,ndim
      read(55,*) xlower(md)
      read(55,*) xupper(md)
    end do
!
!     # boundary conditions:
    read(55,*) mbc
    allocate(mthbc(2*ndim))
    do md=1,2*ndim
      read(55,*) mthbc(md)
    end do
    
    ! close the file
    close(55)

    ! Allocate q and aux
    call qalloc(mbc,nx,meqn,maux,r1,r2,time_integrator)

    ! Work array allocations
    ! all should eventually be here once the arrays are moved into modules
    call alloc_recon_workspace(maxnx,mbc,meqn,mwaves,lim_type,char_decomp)


    if ((mthbc(1).eq.2 .and. mthbc(2).ne.2) .or. &
        (mthbc(2).eq.2 .and. mthbc(1).ne.2)) then
        write(6,*) '*** ERROR ***  periodic boundary conditions'
        write(6,*) ' require mthbc(1) and mthbc(2) BOTH be set to 2'
        stop 
    endif

    ! grid spacing
    allocate(dx(ndim))
    do i=1,ndim
        dx(i) = (xupper(i) - xlower(i)) / float(nx(i))
    enddo

    ! time increments between outputting solution:
    if (outstyle .eq. 1) then
        dtout = (tfinal - t0)/float(nout)
    endif

    mx=nx(1)

    ! call users routine setprob to set any specific parameters
    ! or other initialization required.
    call setprob

    ! set aux array:
    if (maux .gt. 0)  then
        call setaux(mx,mbc,nx,xlower,dx,maux,aux)
    endif

    ! Set initial conditions:
    call qinit(mx,meqn,mbc,nx,xlower,dx,q,maux,aux)

    ! Output initial data
    call output(meqn,mbc,nx,xlower,dx,q,t0,0,aux,maux)
    write(6,601) 0, t0

        open(10,file='fort.info',status='unknown',form='formatted')
    ! =================================================================
    ! Main loop
    ! =================================================================
    tend = t0
    do n=1,nout
        tstart = tend
        select case(outstyle)
            case(1)
                tend = tstart + dtout
            case(2)
                tend = tout(n)
            case(3)
                tend = tstart - 1.d0  !# single-step mode
        end select

        call sharpclaw(tstart,tend,cfl_maxused,cfl_last,dt_minused, &
                dt_maxused,dt_last,steps_taken,info,dt_initial,dt_max, &
                dt_variable,time_integrator,verbosity, &
                cfl_max,cfl_desired,max_steps,r1,r2,bc,rp,tfluct,src,b4step)

        ! check to see if an error occured:
        if (info .ne. 0) then
            write(6,*) '*** ERROR in sharpclaw ***  info =',info
            if (info.eq.1) then
               write(6,*) '***   either mx > mx or mbc < 2'
               endif
            if (info.eq.2) then
               write(6,*) '***   dt does not divide (tend - tstart)'
               write(6,*) '***   and dt is fixed since dt_variable=0'
               endif
            if (info.eq.3) then
               write(6,*) '***   dt_variable=1 and cfl_desired > cfl_max'
               endif
            if (info.eq.4) then
               write(6,*) '***   mwork is too small'
               endif
            if (info.eq.11) then
               write(6,*) '***   Too many times steps, n > max_steps'
               endif
            if (info.eq.12) then
               write(6,*) &
                '***   The Courant number is greater than cfl_max'
               write(6,*) '***   and dt is fixed since dt_variable=0'
               endif

        endif

         dt_initial = dt_last

!        # output solution at this time
!        ------------------------------
!        # if outstyle=1 or 2, then nstepout=1 and we output every time
!        # we reach this point, since wclaw1 was called for the entire time
!        # increment between outputs.

!        # if outstyle=3 then we only output if we have taken nstepout
!        # time steps since the last output.

!        # iframe is the frame number used to form file names in out1
        iframe = n/nstepout
        if (iframe*nstepout == n) then
            call output(meqn,mbc,nx,xlower,dx,q,tend,iframe,aux,maux)
            write(6,601) iframe,tend
            write(10,1010) tend,info,dt_minused,dt_maxused, &
                 dt_last,cfl_maxused,cfl_last,steps_taken
        endif

!        # formats for writing out information about this call to claw:

  601    format('main: Frame ',i4, &
                 ' matlab plot files done at time t =', &
                 d12.4,/)

 1010    format('tend =',d15.4,/, &
             'info =',i5,/,'smallest dt =',d15.4,/,'largest dt =', &
             d15.4,/,'last dt =',d15.4,/,'largest cfl =', &
               d15.4,/,'last cfl =',d15.4,/,'steps taken =',i4,/)

    enddo

    close(10)
    deallocate(q,qold,aux,mthlim,nx,mthbc,xlower,xupper)
    if (outstyle == 2) then
        deallocate(tout)
    endif

    deallocate(r1%qrk)
    if (time_integrator==4) then
        deallocate(r2%qrk)
    endif

    call dealloc_recon_workspace(lim_type,char_decomp)


end program sharpclaw_main
