! ==============================================================
	subroutine output(meqn,mbc,nx,xlower,dx,q,t,iframe,aux,maux)
! ==============================================================
!
! # Output the results for a general system of conservation laws
! # in 2 dimensions
!
! # Here the quantity (h + bottom topo) is printed in the fort.qxxxx file
! # as a first variable, i.e. q[0].


	implicit double precision (a-h,o-z)
    integer, parameter :: ndim=2
    integer, intent(in) :: nx(ndim)
    double precision :: q(1-mbc:nx(1)+mbc, 1-mbc:nx(2)+mbc, meqn)
    double precision :: physVar(1-mbc:nx(1)+mbc, 1-mbc:nx(2)+mbc, meqn)
    double precision :: aux(1-mbc:nx(1)+mbc, 1-mbc:nx(2)+mbc, maux)
    character*10 fname1, fname2, fname3
    double precision, intent(in) :: xlower(ndim), dx(ndim)
    logical outaux
      
    outaux = .false.

!   # Write the results to the file fort.q<iframe>
!   # Use format required by matlab script  plotclaw2.m
!   # The same format is used by the amrclaw package.
!   # Here it's adapted to output just the single grid.
!   # first create the file name and open file
!
    fname1 = 'fort.qxxxx'
    fname2 = 'fort.txxxx'
    fname3 = 'fort.axxxx'
    nstp = iframe
    	do 55 ipos = 10, 7, -1
        	idigit = mod(nstp,10)
            fname1(ipos:ipos) = char(ichar('0') + idigit)
            fname2(ipos:ipos) = char(ichar('0') + idigit)
            fname3(ipos:ipos) = char(ichar('0') + idigit)
            nstp = nstp / 10
 55      continue

	open(unit=50,file=fname1,status='unknown',form='formatted')
    open(unit=60,file=fname2,status='unknown',form='formatted')

!
!   # the following parameters are used in amrclaw where there are
!   # multiple grids.  Here they are all set to 1:
    ngrids = 1
    mptr = 1
    level = 1

    write(50,1001) mptr,level,nx(1),nx(2)
 1001 format(i5,'                 grid_number',/, &
             i5,'                 AMR_level',/, &
             i5,'                 mx',/, &
             i5,'                 my')


      write(50,1002) xlower(1),xlower(2),dx(1),dx(2)
 1002 format(e26.16,'    xlow', /, &
             e26.16,'    ylow', /, &
             e26.16,'    dx', /, &
             e26.16,'    dy',/)
!

      do 20 j=1,nx(2)
        do 10 i=1,nx(1)
          do m=1,meqn
          	 physVar(i,j,1)=q(i,j,1)+aux(i,j,1)
!            # exponents with more than 2 digits cause problems reading
!            # into matlab... reset tiny values to zero:
             if (dabs(physVar(i,j,m)) .lt. 1d-99) physVar(i,j,m) = 0.d0
             enddo
!
          write(50,1005) (physVar(i,j,m), m=1,meqn)
 1005     format(4e26.16)
!
 10       continue
        write(50,*) ' '
 20     continue
      write(50,*) ' '
      
      
      if (outaux) then 
!     # also output the aux arrays:
      open(unit=70,file=fname3,status='unknown',form='formatted')
      write(70,1001) mptr,level,nx(1),nx(2)
      write(70,1002) xlower(1),xlower(2),dx(1),dx(2)
      do 120 j=1,nx(2)
         do 110 i=1,nx(1)
            do m=1,maux
!              # exponents with more than 2 digits cause problems reading
!              # into matlab... reset tiny values to zero:
               if (dabs(aux(i,j,m)) .lt. 1d-99) aux(i,j,m) = 0.d0
            enddo
!
            write(70,1005) (aux(i,j,m), m=1,maux)
!
  110       continue
         write(70,*) ' '
  120    continue
      write(70,*) ' '
      close(unit=70)
      endif


      write(60,1000) t,meqn,ngrids,maux,2
 1000 format(e26.16,'    time', /, &
             i5,'                 meqn'/, &
             i5,'                 ngrids'/, &
             i5,'                 maux'/, &
             i5,'                 ndim'/,/)
!

      close(unit=50)
      close(unit=60)

      return
      end
