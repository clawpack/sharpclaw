! ===================================================================
       subroutine weno_gauss(ixy,q1d,q_gauss,auxl,auxr,w,w5)
! ===================================================================

! This is a subroutine of fifth order WENO reconstruction, for CLAWpack.
! Provided by Yulong Xing
! 2D version for fully multi-d reconstruction
! using Gauss quadrature to calculate the waves                            
! works with weno5_lines2()
! this part just reconstructs along a slice and evaluates at
! Gauss quadrature "lines"

! Accepts q1d: cell averages of q along a slice
! Returns q_gauss: 1D cell averages of q along gauss lines normal to
!  slice direction

! very inefficiently coded for now due to being adapted from the
! version with char. decomp.

      use ClawData
      use ClawParams

      implicit double precision (a-h,o-z)

      type(griddat) g

      parameter (ngauss = 3)
      dimension q1d(1-mbc:maxnx+mbc,meqn)
      dimension q_gauss(1-mbc:maxnx+mbc,meqn,ngauss)
      dimension w(ngauss,ngauss,ngauss), w5(5,ngauss)

      epweno = 1.e-36

!     # loop over all equations (all components).  

      do 99 m=1,meqn
!        # compute and store the differences of the cell averages
         do i=1-mbc,mx+mbc-1
            g%dq(i,m)=q1d(i+1,m)-q1d(i,m)
            enddo
  99    continue

      do 100 m=1,meqn

!       # No characteristic projection

        do  i = 2-mbc,mx+mbc
           do m1 = -2,2
              g%hh(i,m1) = 0.d0
                g%hh(i,m1) = g%dq(i+m1,m)
           enddo !m1 loop

           do m1 = 1,ngauss
              !DK could just loop over gauss points here
                 g%hgg(i,m1,1) = ( w(1,m1,1)*q1d(i-3+m1,m) &
                               +w(1,m1,2)*q1d(i-2+m1,m) &
                               +w(1,m1,3)*q1d(i-1+m1,m))
                 g%hgg(i,m1,2) = ( w(2,m1,1)*q1d(i-3+m1,m) &
                               +w(2,m1,2)*q1d(i-2+m1,m) &
                               +w(2,m1,3)*q1d(i-1+m1,m))
                 g%hgg(i,m1,3) = ( w(3,m1,1)*q1d(i-3+m1,m) &
                               +w(3,m1,2)*q1d(i-2+m1,m) &
                               +w(3,m1,3)*q1d(i-1+m1,m))

           enddo !m1 loop
        enddo !i loop


!        # the reconstruction


!          # m1=1: construct ql

           im=1
           ione=im
           inone=-im
           intwo=-2*im
  
           do 20 i=0,mx
  
             t1=im*(g%hh(i,intwo)-g%hh(i,inone))
             t2=im*(g%hh(i,inone)-g%hh(i,0    ))
             t3=im*(g%hh(i,0    )-g%hh(i,ione ))
  
             tt1=13.*t1**2+3.*(   g%hh(i,intwo)-3.*g%hh(i,inone))**2
             tt2=13.*t2**2+3.*(   g%hh(i,inone)+   g%hh(i,0    ))**2
             tt3=13.*t3**2+3.*(3.*g%hh(i,0    )-   g%hh(i,ione ))**2
       
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
  
             g%hff(i,m,1) = s1*(g%hgg(i,1,1)-g%hgg(i,2,1)) + g%hgg(i,2,1) &
                         + s3*(g%hgg(i,3,1)-g%hgg(i,2,1))

             d1 = w5(1,2)/w(2,1,1)
             d3 = w5(5,2)/w(2,3,3)
             d2 = 1-d1-d3

!  Negative weights here: special treatment
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
!         positive part
          s1 = d1pp * tt2 * tt3
          s2 = d2pp * tt1 * tt3
          s3 = d3pp * tt1 * tt2
          t0 = 1. / ( s1 + s2 + s3 )
          s1 = s1 * t0
          s3 = s3 * t0
!	    
          g%hff(i,m,2) = sump* ( s1*(g%hgg(i,1,2)-g%hgg(i,2,2)) &
                      + g%hgg(i,2,2) + s3*(g%hgg(i,3,2)-g%hgg(i,2,2)) )
!         negative part
          s1 = d1nn * tt2 * tt3
          s2 = d2nn * tt1 * tt3
          s3 = d3nn * tt1 * tt2
          t0 = 1. / ( s1 + s2 + s3 )
          s1 = s1 * t0
          s3 = s3 * t0
!	    
        g%hff(i,m,2) = g%hff(i,m,2) - sumn* (s1*(g%hgg(i,1,2)-g%hgg(i,2,2)) &
                    + g%hgg(i,2,2) + s3*(g%hgg(i,3,2)-g%hgg(i,2,2)) )
!
          d1 = w5(1,3) / w(3,1,1)
          d3 = w5(5,3) / w(3,3,3)
          d2= 1 - d1 - d3
 
          s1 = d1 * tt2 * tt3
          s2 = d2 * tt1 * tt3
          s3 = d3 * tt1 * tt2
          t0 = 1. / ( s1 + s2 + s3 )
          s1 = s1 * t0
          s3 = s3 * t0
!	    
          g%hff(i,m,3) = s1*(g%hgg(i,1,3)-g%hgg(i,2,3)) &
                     + g%hgg(i,2,3) + s3*(g%hgg(i,3,3)-g%hgg(i,2,3))


   20        continue
  100    continue
!----------------- loop in "m"  ends  here  -------------------

        do ig =  1,3
!       ig=1 ---> -sqrt(0.6), ig=2 ---> 0, ig=3 ---> sqrt(0.6)
           do m =  1, meqn
              do i = 2-mbc,mx+mbc
                q_gauss(i,m,ig) = g%hff(i,m,ig)
              enddo ! loop over cells (i)
           enddo ! loop over eqns (m)
        enddo !loop over gauss points (ig)

      return
      end
