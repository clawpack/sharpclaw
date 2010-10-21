! =========================================================
subroutine output(meqn,mbc,mx,xlower,dx,q,t,iframe,aux,maux)
! =========================================================

    ! Output the results for a general system of conservation laws
    ! in 1 dimension

    ! Write the results to the file fort.q<iframe>
    ! Use format required by matlab script  plotclaw1.m
    ! The same format is used by the amrclaw package.
    ! Here it's adapted to output just the single grid.

    ! set outaux = .true. to also output the aux arrays to fort.a<iframe>

    implicit double precision (a-h,o-z)
    integer, parameter :: ndim=1
    double precision :: q(1-mbc:mx+mbc, meqn)
    double precision :: aux(1-mbc:mx+mbc, maux)
    double precision, intent(in) :: xlower, dx
    character*10 fname1, fname2, fname3
    logical outaux

    outaux = .false.

    ! Write the results to the file fort.q<iframe>
    ! Use format required by matlab script  plotclaw1.m
    ! The same format is used by the amrclaw package.  
    ! Here it's adapted to output just the single grid.

    ! first create the file name and open file

    fname1 = 'fort.qxxxx'
    fname2 = 'fort.txxxx'
    fname3 = 'fort.axxxx'
    nstp = iframe
    do ipos = 10, 7, -1
        idigit = mod(nstp,10)
        fname1(ipos:ipos) = char(ichar('0') + idigit)
        fname2(ipos:ipos) = char(ichar('0') + idigit)
        fname3(ipos:ipos) = char(ichar('0') + idigit)
        nstp = nstp / 10
    enddo

    open(unit=50,file=fname1,status='unknown',form='formatted')
    open(unit=60,file=fname2,status='unknown',form='formatted')


    ! the following parameters are used in amrclaw where there are
    ! multiple grids.  Here they are all set to 1:
    ngrids = 1
    mptr = 1
    level = 1

    write(50,1001) mptr,level,mx
 1001 format(i5,'                 grid_number',/, &
             i5,'                 AMR_level',/, &
             i5,'                 mx')

    write(50,1002) xlower,dx
 1002 format(e18.8,'    xlow', /, &
             e18.8,'    dx', /)

    do i=1,mx
        do m=1,meqn
            ! exponents with more than 2 digits cause problems reading
            ! into matlab... reset tiny values to zero:
            if (dabs(q(i,m)) .lt. 1d-99) q(i,m) = 0.d0
        enddo
        write(50,1005) (q(i,m), m=1,meqn)
 1005     format(4e16.8)

    enddo
    write(50,*) ' '
    write(50,*) ' '

    if (outaux) then 
        ! also output the aux arrays:
        open(unit=70,file=fname3,status='unknown',form='formatted')
        write(70,1001) mptr,level,mx
        write(70,1002) xlower,dx
        do i=1,mx
            do m=1,maux
                ! exponents with more than 2 digits cause problems reading
                ! into matlab... reset tiny values to zero:
                if (dabs(aux(i,m)) .lt. 1d-99) aux(i,m) = 0.d0
            enddo

            write(70,1005) (aux(i,m), m=1,maux)

        enddo
        write(70,*) ' '
        close(unit=70)
    endif

    write(60,1000) t,meqn,ngrids,maux,1

 1000 format(e26.16,'    time', /, &
             i5,'                 meqn'/, &
             i5,'                 ngrids'/, &
             i5,'                 maux'/, &
             i5,'                 ndim'/,/)

    close(unit=50)
    close(unit=60)

end subroutine output
