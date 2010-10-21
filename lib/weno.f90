! ===================================================================
subroutine weno5(q,g,mx)
! ===================================================================

! This is a subroutine of fifth order WENO reconstruction, for CLAWpack. 
! Provided by Chi-Wang Shu

      use ClawData
      use ClawParams

      implicit double precision (a-h,o-z)

      type(griddat) g
      dimension  q(1-mbc:mx+mbc,meqn)

      epweno = 1.e-36

!     # loop over all equations (all components).  
!     # the reconstruction is performed component-wise, no attempt has been
!     # done to implement characteristic decompositions in the reconstruction

      do 100 m=1,meqn

!        # compute and store the differences of the cell averages

         do i=2-mbc,mx+mbc
            g%dq1m(i)=q(i,m)-q(i-1,m)
            enddo

!        # the reconstruction

         do 50 m1=1,2

!          # m1=1: construct ql
!          # m1=2: construct qr

           im=(-1)**(m1+1)
           ione=im
           inone=-im
           intwo=-2*im
  
           do 20 i=1,mx+1
  
             t1=im*(g%dq1m(i+intwo)-g%dq1m(i+inone))
             t2=im*(g%dq1m(i+inone)-g%dq1m(i      ))
             t3=im*(g%dq1m(i      )-g%dq1m(i+ione ))
  
             tt1=13.*t1**2+3.*(   g%dq1m(i+intwo)-3.*g%dq1m(i+inone))**2
             tt2=13.*t2**2+3.*(   g%dq1m(i+inone)+   g%dq1m(i      ))**2
             tt3=13.*t3**2+3.*(3.*g%dq1m(i      )-   g%dq1m(i+ione ))**2
       
             tt1=(epweno+tt1)**2
             tt2=(epweno+tt2)**2
             tt3=(epweno+tt3)**2
             s1 =tt2*tt3
             s2 =6.*tt1*tt3
             s3 =3.*tt1*tt2
             t0 =1./(s1+s2+s3)
             s1 =s1*t0
             s3 =s3*t0
  
             g%uu(i,m1) = (s1*(t2-t1)+(0.5*s3-0.25)*(t3-t2))/3. &
                    +(-q(i-2,m)+7.*(q(i-1,m)+q(i,m))-q(i+1,m))/12.

   20        continue
   50      continue

           g%qr(0:mx,m)=g%uu(1:mx+1,1)
           g%ql(1:mx+1,m)=g%uu(1:mx+1,2)

  100    continue

      return
      end

! ===================================================================
subroutine weno5_char(q,g,mx)
! ===================================================================

! This one uses characteristic decomposition
!  evl, evr are left and right eigenvectors at each interface

    use ClawData
    use ClawParams

    implicit double precision (a-h,o-z)

    type(griddat) g
    dimension q(1-mbc:maxnx+mbc,meqn)

    epweno = 1.e-36

    ! loop over all equations (all components).  
    ! the reconstruction is performed using characteristic decomposition

      do 99 m=1,meqn
!        # compute and store the differences of the cell averages
         do i=3-mbc,mx+mbc-2
            g%dq(i,m)=q(i,m)-q(i-1,m)
!           # Compute the part of the reconstruction that is
!              stencil-independent
            g%qr(i-1,m) = (-q(i-2,m)+7.*(q(i-1,m)+q(i,m))-q(i+1,m))/12.
            g%ql(i,m)   = g%qr(i-1,m)
            enddo
  99    continue

      do 100 ip=1,meqn

!       # Project the difference of the cell averages to the
!       # 'm'th characteristic field

        do m2 = -2,2
           do  i = 2-mbc,mx+mbc
              g%hh(i,m2) = 0.d0
              do m=1,meqn 
                g%hh(i,m2) = g%hh(i,m2)+ g%evl(i,ip,m)*g%dq(i+m2,m)
              enddo
           enddo
        enddo


!        # the reconstruction

         do 50 m1=1,2

!          # m1=1: construct ql
!          # m1=2: construct qr

           im=(-1)**(m1+1)
           ione=im
           inone=-im
           intwo=-2*im
  
           do 20 i=1,mx+1
  
             t1=im*(g%hh(i,intwo)-g%hh(i,inone))
             t2=im*(g%hh(i,inone)-g%hh(i,0    ))
             t3=im*(g%hh(i,0    )-g%hh(i,ione ))
  
             tt1=13.*t1**2+3.*(   g%hh(i,intwo)-3.*g%hh(i,inone))**2
             tt2=13.*t2**2+3.*(   g%hh(i,inone)+   g%hh(i,0    ))**2
             tt3=13.*t3**2+3.*(3.*g%hh(i,0    )-   g%hh(i,ione ))**2

             tt1=(epweno+tt1)**2
             tt2=(epweno+tt2)**2
             tt3=(epweno+tt3)**2
             s1 =tt2*tt3
             s2 =6.*tt1*tt3
             s3 =3.*tt1*tt2
             t0 =1./(s1+s2+s3)
             s1 =s1*t0
             s3 =s3*t0
  
             g%uu(i,m1) = ( s1*(t2-t1) + (0.5*s3-0.25)*(t3-t2) ) /3.

   20        continue !end loop over interfaces
   50      continue !end loop over which side of interface
! Project to the physical space:
           do m = 1,meqn
           do i = 2-mbc,mx+mbc
                 g%qr(i-1,m) = g%qr(i-1,m) + g%evr(i,m,ip)*g%uu(i,1)
                 g%ql(i  ,m) = g%ql(i  ,m) + g%evr(i,m,ip)*g%uu(i,2)
           enddo
           enddo
  100    continue !end loop over waves


      return
      end

! ===================================================================
subroutine weno5_trans(q,g,mx)
! ===================================================================

!      Transmission-based WENO reconstruction

    use ClawData
    use ClawParams

    implicit double precision (a-h,o-z)

    type(griddat) g

    external evec
    dimension q(1-mbc:maxnx+mbc,meqn)

    epweno = 1.e-36


!     # loop over all equations (all components).  
!     # the reconstruction is performed using characteristic decomposition

      do 98 m=1,meqn
!        # compute and store the differences of the cell averages
         do i=2-mbc,mx+mbc
            g%dq(i,m)=q(i,m)-q(i-1,m)
            enddo
  98    continue

!       #Compute left and right eigenvectors
      call evec(mx,meqn,mbc,mx,q,g%auxl,g%auxr,g%evl,g%evr)


!       # Find wave strengths at each interface
!       # 'm'th characteristic field
      do 99 mw=1,meqn
           do  i = 2-mbc,mx+mbc
              g%gg(i,mw) = 0.d0
              do m=1,meqn
                g%gg(i,mw) = g%gg(i,mw)+ g%evl(i,mw,m)*g%dq(i,m)
              enddo
           enddo
  99    continue

      do 100 mw=1,meqn


!       # Project the waves to the
!       # 'm'th characteristic field

        do m1 = -2,2
           do  i = 2-mbc,mx+mbc
              g%hh(i,m1) = 0.d0
              do m=1,meqn 
                g%hh(i,m1) = g%hh(i,m1)+g%evl(i,mw,m)* &
                            g%gg(i+m1,mw)*g%evr(i+m1,m,mw)
              enddo
           enddo
        enddo


!        # the reconstruction

         do 50 m1=1,2

!          # m1=1: construct ql
!          # m1=2: construct qr

           im=(-1)**(m1+1)
           ione=im
           inone=-im
           intwo=-2*im
  
           do 20 i=1,mx+1
  
             t1=im*(g%hh(i,intwo)-g%hh(i,inone))
             t2=im*(g%hh(i,inone)-g%hh(i,0    ))
             t3=im*(g%hh(i,0    )-g%hh(i,ione ))
  
             tt1=13.*t1**2+3.*(   g%hh(i,intwo)-3.*g%hh(i,inone))**2
             tt2=13.*t2**2+3.*(   g%hh(i,inone)+   g%hh(i,0    ))**2
             tt3=13.*t3**2+3.*(3.*g%hh(i,0    )-   g%hh(i,ione ))**2
       
             tt1=(epweno+tt1)**2
             tt2=(epweno+tt2)**2
             tt3=(epweno+tt3)**2
             s1 =tt2*tt3
             s2 =6.*tt1*tt3
             s3 =3.*tt1*tt2
             t0 =1./(s1+s2+s3)
             s1 =s1*t0
             s3 =s3*t0
  
             g%u(i,mw,m1) = ( s1*(t2-t1) + (0.5*s3-0.25)*(t3-t2) ) /3.

   20        continue
   50      continue
  100    continue

! Project to the physical space:

        do m1 =  1,2
           do m =  1, meqn
              do i = 2-mbc,mx+mbc
                 g%uh(i,m,m1) =( -q(i-2,m) + 7*( q(i-1,m)+q(i,m) ) &
                                         - q(i+1,m) )/12.
                 do mw=1,meqn 
                g%uh(i,m,m1) = g%uh(i,m,m1) + g%evr(i,m,mw)*g%u(i,mw,m1)
                 enddo
              enddo
           enddo
        enddo

      do m =  1, meqn
         do i=1,mx+1
           g%qr(i-1,m)=g%uh(i,m,1)
           g%ql(i,m)=g%uh(i,m,2)
         enddo
      enddo

      return
      end

! ===================================================================
subroutine weno5_wave(q,ql,qr,wave,s,mx)
! ===================================================================
!
!  Fifth order WENO reconstruction, based on waves
!  which are later interpreted as slopes.

    use ClawData
    use ClawParams

    implicit double precision (a-h,o-z)

    dimension q(1-mbc:maxnx+mbc,meqn)
    dimension ql(1-mbc:maxnx+mbc,meqn),qr(1-mbc:maxnx+mbc,meqn)
    dimension wave(1-mbc:maxnx+mbc,meqn,mwaves)
    dimension s(1-mbc:maxnx+mbc,mwaves)
    dimension u(2)

    epweno = 1.e-36

    ! loop over interfaces (i-1/2)
    do i=3-mbc,mx+mbc-2
        ! Compute the part of the reconstruction that is stencil-independent
        do m=1,meqn
            qr(i-1,m) = (-q(i-2,m)+7.*(q(i-1,m)+q(i,m))-q(i+1,m))/12.
            ql(i,m)   = qr(i-1,m)
        enddo
        ! the reconstruction is performed in terms of waves
        do mw=1,mwaves
            ! loop over which side of x_i-1/2 we're on
            do m1=1,2
                ! m1=1: construct q^-_{i-1/2}
                ! m1=2: construct q^+_{i-1/2}
                im=(-1)**(m1+1)
                ione=im
                inone=-im
                intwo=-2*im
  
                ! find thetas
                wnorm2 = wave(i      ,1,mw)*wave(i,1,mw)
                theta1 = wave(i+intwo,1,mw)*wave(i,1,mw)
                theta2 = wave(i+inone,1,mw)*wave(i,1,mw)
                theta3 = wave(i+ione ,1,mw)*wave(i,1,mw)
                do m=2,meqn
                    wnorm2 = wnorm2 + wave(i,      m,mw)*wave(i,m,mw)
                    theta1 = theta1 + wave(i+intwo,m,mw)*wave(i,m,mw)
                    theta2 = theta2 + wave(i+inone,m,mw)*wave(i,m,mw)
                    theta3 = theta3 + wave(i+ione ,m,mw)*wave(i,m,mw)
                enddo

                t1=im*(theta1-theta2)
                t2=im*(theta2-wnorm2)
                t3=im*(wnorm2-theta3)
  
                tt1=13.*t1**2+3.*(theta1   -3.*theta2)**2
                tt2=13.*t2**2+3.*(theta2   +   wnorm2)**2
                tt3=13.*t3**2+3.*(3.*wnorm2-   theta3)**2
       
                tt1=(epweno+tt1)**2
                tt2=(epweno+tt2)**2
                tt3=(epweno+tt3)**2
                s1 =tt2*tt3
                s2 =6.*tt1*tt3
                s3 =3.*tt1*tt2
                t0 =1./(s1+s2+s3)
                s1 =s1*t0
                s3 =s3*t0
  
                if(wnorm2.gt.1.e-14) then
                    u(m1) = ( s1*(t2-t1) + (0.5*s3-0.25)*(t3-t2) ) /3.
                    wnorm2=1.d0/wnorm2
                else
                    u(m1) = 0.d0
                    wnorm2=0.d0
                endif
            enddo !end loop over which side of interface
            do m=1,meqn
                qr(i-1,m) = qr(i-1,m) +  u(1)*wave(i,m,mw)*wnorm2
                ql(i  ,m) = ql(i  ,m) +  u(2)*wave(i,m,mw)*wnorm2
            enddo
        enddo !loop over waves
    enddo !loop over interfaces

end subroutine weno5_wave
