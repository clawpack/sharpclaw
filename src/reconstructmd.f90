module reconstructmd
    integer, parameter :: ngauss=3
    double precision, allocatable, private :: uu(:,:),dq(:,:)
    double precision, allocatable, private :: uh(:,:,:),gg(:,:),hh(:,:),u(:,:,:)
    double precision, allocatable, private :: evl(:,:,:),evr(:,:,:)
    double precision, private :: w(ngauss,ngauss,ngauss), w5(5,ngauss)
    double precision, private  :: epweno = 1.e-36

    double precision, allocatable :: hgg(:,:,:), hff(:,:,:)
    double precision, allocatable :: du(:,:,:),uhh(:,:,:,:),ugg(:,:,:,:)


contains

    subroutine alloc_reconmd_workspace(maxnx,nx,mbc,meqn,mwaves,char_decomp)

        integer,intent(in) :: maxnx,mbc,meqn,mwaves,char_decomp
        integer,intent(in) :: nx(2)

        allocate(hgg(maxnx+2*mbc,ngauss,ngauss))
        allocate(hff(maxnx+2*mbc,meqn,ngauss))
        allocate(uhh(maxnx+2*mbc,meqn,2,ngauss))
        allocate(ugg(maxnx+2*mbc,meqn,2,ngauss))

        allocate(dq(maxnx+2*mbc,meqn))
        allocate(uu(maxnx+2*mbc,2))
        allocate(hh(maxnx+2*mbc,-2:2))

    end subroutine alloc_reconmd_workspace

    subroutine dealloc_reconmd_workspace(char_decomp)
        integer,intent(in) :: char_decomp

        deallocate(hgg)
        deallocate(hff)
        deallocate(uhh)
        deallocate(ugg)
        deallocate(dq)
        deallocate(uu)
        deallocate(hh)

    end subroutine dealloc_reconmd_workspace

    subroutine set_md_weno_weights()
        !Gauss quadrature weights
        !This should be generalized for different quadrature choices
        gp = dsqrt(0.6d0)

        !Correspondence to Yulong Xing's code:
        !w(1,*,*) = wa(*,*)
        !w(2,*,*) = wb(*,*)
        !w(3,*,*) = wc(*,*)

        w(1,1,1) = 1./30. - gp/4.
        w(1,1,2) = -1./15. + gp
        w(1,1,3) = 31./30. - 3.*gp/4.
        w(1,2,1) = 1./30. + gp/4.
        w(1,2,2) = 14./15.
        w(1,2,3) = 1./30. - gp/4. 
        w(1,3,1) = 31./30. + 3.*gp/4.
        w(1,3,2) = -1./15. - gp 
        w(1,3,3) = 1./30. + gp/4.
        
        w(2,1,1) = -1./24.
        w(2,1,2) = 1./12.
        w(2,1,3) = 23./24.
        w(2,2,1) = -1./24.
        w(2,2,2) = 13/12.
        w(2,2,3) = -1./24.
        w(2,3,1) = 23./24.
        w(2,3,2) = 1./12.
        w(2,3,3) = -1./24.
        
        w(3,1,1) = 1./30. + gp/4. 
        w(3,1,2) = -1./15. - gp 
        w(3,1,3) = 31./30. + 3.*gp/4.
        w(3,2,1) = 1./30. - gp/4. 
        w(3,2,2) = 14./15.
        w(3,2,3) = 1./30. + gp/4. 
        w(3,3,1) = 31./30. - 3.*gp/4.
        w(3,3,2) = -1./15. + gp
        w(3,3,3) = 1./30. - gp/4.
        
        !for point -sqrt(3/5)
        w5(1,1) = -11./240.*gp - 3./800.
        w5(2,1) = 41./120.*gp+29./600.
        w5(3,1) = 1093./1200.
        w5(4,1) = -41./120.*gp+29./600.
        w5(5,1) = 11./240.*gp - 3./800.
        !for point 0
        w5(1,2) = 3./640.
        w5(2,2) = -29./480.
        w5(3,2) = 1067./960.
        w5(4,2) = -29./480.
        w5(5,2) = 3./640.
        !for point sqrt(3/5)
        w5(1,3) = 11./240.*gp - 3./800.
        w5(2,3) = -41./120.*gp+29./600.
        w5(3,3) = 1093./1200.
        w5(4,3) = 41./120.*gp+29./600.
        w5(5,3) = -11./240.*gp - 3./800.

    end subroutine set_md_weno_weights


    ! ===================================================================
    subroutine weno_gauss(q1d,q_gauss,mx)
    ! ===================================================================

    ! fifth order WENO reconstruction
    ! 2D version for fully multi-d reconstruction
    ! using Gauss quadrature to calculate the waves                            
    ! works with weno5_lines2()
    ! this part just reconstructs along a slice and evaluates at
    ! Gauss quadrature "lines"
    ! Adapted from code of Yulong Xing

    ! Accepts q1d: cell averages of q along a slice
    ! Returns q_gauss: 1D cell averages of q along gauss lines normal to
    !  slice direction

    ! very inefficiently coded for now due to being adapted from the
    ! version with char. decomp.

    implicit double precision (a-h,o-z)
    double precision, intent(in) :: q1d(:,:)
    double precision, intent(out) :: q_gauss(:,:,:)
    integer, intent(in) :: mx

    integer, parameter :: mbc=3
    integer :: meqn, mx2

    mx2  = mx+2*mbc; meqn = size(q1d,2)


    ! loop over all equations (all components).  

    forall(m=1:meqn,i=2:mx2)
        ! compute and store the differences of the cell averages
        dq(i,m)=q1d(i+1,m)-q1d(i,m)
    end forall

    do m=1,meqn

        ! No characteristic projection
        do  i = mbc,mx2-mbc+1
            do m1 = -2,2
                hh(i,m1) = dq(i+m1,m)
            enddo !m1 loop

            do m1 = 1,ngauss
                !DK could just loop over gauss points here
                hgg(i,m1,1) = ( w(1,m1,1)*q1d(i-3+m1,m) &
                               +w(1,m1,2)*q1d(i-2+m1,m) &
                               +w(1,m1,3)*q1d(i-1+m1,m))
                hgg(i,m1,2) = ( w(2,m1,1)*q1d(i-3+m1,m) &
                               +w(2,m1,2)*q1d(i-2+m1,m) &
                               +w(2,m1,3)*q1d(i-1+m1,m))
                hgg(i,m1,3) = ( w(3,m1,1)*q1d(i-3+m1,m) &
                               +w(3,m1,2)*q1d(i-2+m1,m) &
                               +w(3,m1,3)*q1d(i-1+m1,m))
           enddo !m1 loop
        enddo !i loop

        ! m1=1: construct ql

        im=1
        ione=im
        inone=-im
        intwo=-2*im
  
        do i=3,mx2-2
  
            t1=im*(hh(i,intwo)-hh(i,inone))
            t2=im*(hh(i,inone)-hh(i,0    ))
            t3=im*(hh(i,0    )-hh(i,ione ))
  
            tt1=13.*t1**2+3.*(   hh(i,intwo)-3.*hh(i,inone))**2
            tt2=13.*t2**2+3.*(   hh(i,inone)+   hh(i,0    ))**2
            tt3=13.*t3**2+3.*(3.*hh(i,0    )-   hh(i,ione ))**2
       
            tt1=(epweno+tt1)**2
            tt2=(epweno+tt2)**2
            tt3=(epweno+tt3)**2

            d1 = w5(1,1) / w(1,1,1)
            d3 = w5(5,1) / w(1,3,3)
            d2 = 1 - d1 -d3

            s1 =d1*tt2*tt3
            s2 =d2*tt1*tt3
            s3 =d3*tt1*tt2
            t0 =1./(s1+s2+s3)
            s1 =s1*t0
            s3 =s3*t0
  
            hff(i,m,1) = s1*(hgg(i,1,1)-hgg(i,2,1)) + hgg(i,2,1) &
                         + s3*(hgg(i,3,1)-hgg(i,2,1))

            d1 = w5(1,2)/w(2,1,1)
            d3 = w5(5,2)/w(2,3,3)
            d2 = 1-d1-d3

            ! Negative weights here: special treatment
            d1p = 0.5d0*( d1 + 3.d0*dabs(d1) )
            d2p = 0.5d0*( d2 + 3.d0*dabs(d2) )
            d3p = 0.5d0*( d3 + 3.d0*dabs(d3) )
            d1n = d1p - d1
            d2n = d2p - d2
            d3n = d3p - d3
            sump = d1p+d2p+d3p
            d1pp = d1p / sump
            d2pp = d2p / sump
            d3pp = d3p / sump
            sumn = d1n+d2n+d3n
            d1nn = d1n / sumn
            d2nn = d2n / sumn
            d3nn = d3n / sumn
            ! positive part
            s1 = d1pp * tt2 * tt3
            s2 = d2pp * tt1 * tt3
            s3 = d3pp * tt1 * tt2
            t0 = 1. / ( s1 + s2 + s3 )
            s1 = s1 * t0
            s3 = s3 * t0

            hff(i,m,2) = sump* ( s1*(hgg(i,1,2)-hgg(i,2,2)) &
                      + hgg(i,2,2) + s3*(hgg(i,3,2)-hgg(i,2,2)) )
            ! negative part
            s1 = d1nn * tt2 * tt3
            s2 = d2nn * tt1 * tt3
            s3 = d3nn * tt1 * tt2
            t0 = 1. / ( s1 + s2 + s3 )
            s1 = s1 * t0
            s3 = s3 * t0

            hff(i,m,2) = hff(i,m,2) - sumn* (s1*(hgg(i,1,2)-hgg(i,2,2)) &
                    + hgg(i,2,2) + s3*(hgg(i,3,2)-hgg(i,2,2)) )

            d1 = w5(1,3) / w(3,1,1)
            d3 = w5(5,3) / w(3,3,3)
            d2= 1 - d1 - d3
            
            s1 = d1 * tt2 * tt3
            s2 = d2 * tt1 * tt3
            s3 = d3 * tt1 * tt2
            t0 = 1. / ( s1 + s2 + s3 )
            s1 = s1 * t0
            s3 = s3 * t0

            hff(i,m,3) = s1*(hgg(i,1,3)-hgg(i,2,3)) &
                     + hgg(i,2,3) + s3*(hgg(i,3,3)-hgg(i,2,3))


        enddo
    enddo
    !----------------- loop in "m"  ends  here  -------------------

    do ig =  1,3
        ! ig=1 ---> -sqrt(0.6), ig=2 ---> 0, ig=3 ---> sqrt(0.6)
        do m =  1, meqn
            do i = mbc,mx2-mbc+1
                q_gauss(i,m,ig) = hff(i,m,ig)
            enddo ! loop over cells (i)
        enddo ! loop over eqns (m)
    enddo !loop over gauss points (ig)

    return
    end subroutine weno_gauss

    ! ===================================================================
    subroutine weno_lines(mx,q_gauss,ql,qr)
    ! ===================================================================

    ! fifth order WENO reconstruction
    ! Adapted from code provided by Yulong Xing
    ! 2D version for fully multi-d reconstruction
    ! Using Gauss quadrature to calculate the waves                            
    ! works with weno5_gauss2()
    ! this part reconstructs along ngauss Gauss quadrature "lines"
    ! in a slice

    ! Takes q_gauss: 1D averages of q in each cell at gauss pts. in other
    !                 dimension
    ! Returns ql,qr: Reconstructed point values on either side of each
    !                 interface, at gauss quadrature points on that face

    !      I am assuming that mbc=3 and all the q values from -2 to mx+3
    !      are available

    ! very inefficiently coded for now due to adaptation from the
    ! version with char. decomp.

    implicit double precision (a-h,o-z)
    double precision, intent(in) :: q_gauss(:,:,:)
    double precision, intent(out) :: ql(:,:,:),qr(:,:,:)
    integer, intent(in) :: mx

    integer, parameter :: mbc=3
    integer :: meqn, mx2

    mx2  = mx+2*mbc; meqn = size(q_gauss,2)

    ! loop over all equations (all components).  

    ! compute and store the differences of the cell averages
    ! along each gauss line
    forall(ig=1:ngauss,m=1:meqn,i=2:mx2)
        hff(i,m,ig)=q_gauss(i,m,ig)-q_gauss(i-1,m,ig)
    endforall

    do m=1,meqn

        ! No characteristic projection

        do ig=1,ngauss

            forall(i=mbc:mx2-mbc+1,m1=-2:2)
                hh(i,m1) = hff(i+m1,m,ig)
            endforall

            ! the reconstruction
            ! m1=1: construct ql
            ! m1=2: construct qr

            do m1=1,2
        
                im=(-1)**(m1+1)
                ione=im; inone=-im; intwo=-2*im
  
                do i=1,mx+1
  
                    t1=im*(hh(i,intwo)-hh(i,inone))
                    t2=im*(hh(i,inone)-hh(i,0    ))
                    t3=im*(hh(i,0    )-hh(i,ione ))
  
                    tt1=13.*t1**2+3.*(   hh(i,intwo)-3.*hh(i,inone))**2
                    tt2=13.*t2**2+3.*(   hh(i,inone)+   hh(i,0    ))**2
                    tt3=13.*t3**2+3.*(3.*hh(i,0    )-   hh(i,ione ))**2
       
                    tt1=(epweno+tt1)**2
                    tt2=(epweno+tt2)**2
                    tt3=(epweno+tt3)**2
                    s1 =tt2*tt3
                    s2 =6.*tt1*tt3
                    s3 =3.*tt1*tt2
                    t0 =1./(s1+s2+s3)
                    s1 =s1*t0
                    s3 =s3*t0
  
                    ugg(i,m,m1,ig) = ( s1*(t2-t1) + (0.5*s3-0.25)*(t3-t2) ) /3.



                enddo !Loop over cells (i)
            enddo !Loop over interface side (m1)
        enddo !Loop over gauss lines
    enddo !Loop over equations (m)
    !----------------- loop in "m"  ends  here  -------------------

    forall(ig=1:ngauss,m1=1:2,m=1:meqn,i=2:mx2-mbc+1)
        uhh(i,m,m1,ig)=(-q_gauss(i-2,m,ig)+7.d0*(q_gauss(i-1,m,ig)+ &
                         q_gauss(i,m,ig))-q_gauss(i+1,m,ig) )/12.d0
        uhh(i,m,m1,ig) = uhh(i,m,m1,ig) + ugg(i,m,m1,ig)

    end forall

    forall(m=1:meqn,i=mbc+1:mx2-mbc,ig=1:ngauss)
        qr(i-1,m,ig)=uhh(i,m,1,ig)
        ql(i,m,ig)=uhh(i,m,2,ig)
    end forall
                
    return

    end subroutine weno_lines

end module reconstructmd
