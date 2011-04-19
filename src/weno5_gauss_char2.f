c
c ===================================================================
       subroutine weno5_gauss_char2(ixy,maxm,meqn,mbc,mx,q1d,q_gauss,
     2                          auxl,auxr,w,w5)
c ===================================================================

c This is a subroutine of fifth order WENO reconstruction, for CLAWpack.
c Provided by Yulong Xing
c This one uses characteristic decomposition 
c 2D version for fully multi-d reconstruction
c using Gauss quadrature to calculate the waves                            
c works with weno5_char_lines2()
c this part just reconstructs along a slice and evaluates at
c Gauss quadrature "lines"

c Accepts q1d: 2D cell averages of q along a slice
c Returns q_gauss: 1D cell averages of q along gauss lines normal to
c  slice direction

      implicit double precision (a-h,o-z)

      external  evecint2
      dimension q1d(1-mbc:maxm+mbc,meqn)
      dimension q_gauss(1-mbc:maxm+mbc,meqn,ngauss)
      dimension evl(1-mbc:maxm+mbc,meqn,meqn)
      dimension evr(1-mbc:maxm+mbc,meqn,meqn)
      dimension w(ngauss,ngauss,ngauss), w5(5,ngauss)

      parameter (mbcweno=3, mxweno=10000, meqnweno=6)
      parameter (ngauss = 3)
      dimension dq(1-mbcweno:mxweno+mbcweno,meqnweno)
      dimension hh(1-mbcweno:mxweno+mbcweno,-2:2)
      dimension hgg(1-mbcweno:mxweno+mbcweno,ngauss,ngauss)
      dimension hff(1-mbcweno:mxweno+mbcweno,meqnweno,ngauss)


      if (maxm .gt. mxweno) then
         write(6,*) 'Error -- increase mxweno'
         stop 
         endif

      if (meqn .gt. meqnweno) then
         write(6,*) 'Error -- increase meqnweno'
         stop 
         endif

      if (mbcweno .ne. mbc) then
         write(6,*) 'Error -- mbcweno in weno.f must agree with mbc'
         stop 
         endif

c     # epweno is a parameter used in WENO weights to prevent the denominator
c     # to become zero.  it does NOT need to be changed

      epweno = 1.e-29

c     # loop over all equations (all components).  
c     # the reconstruction is performed using characteristic decomposition


      do 99 m=1,meqn
c        # compute and store the differences of the cell averages
         do i=1-mbc,mx+mbc-1
            dq(i,m)=q1d(i+1,m)-q1d(i,m)
            enddo
  99    continue

c       #Compute left and right eigenvectors
      call evecint2(ixy,maxm,meqn,mbc,mx,q1d,auxl,auxr,evl,evr)

      do 100 m=1,meqn

c       # Project the difference of the cell averages to the
c       # 'm'th characteristic field

        do  i = 2-mbc,mx+mbc
           do m1 = -2,2
              hh(i,m1) = 0.d0
              do mm=1,meqn 
                hh(i,m1) = hh(i,m1)+ evl(i,m,mm)*dq(i+m1,mm)
              enddo !meqn loop
           enddo !m1 loop

           do m1 = 1,ngauss
              hgg(i,m1,1)=0.d0
              hgg(i,m1,2)=0.d0
              hgg(i,m1,3)=0.d0
              do mm=1,meqn
              !DK could just loop over gauss points here
                 hgg(i,m1,1) = hgg(i,m1,1) + evl(i,m,mm) *
     &                         ( w(1,m1,1)*q1d(i-3+m1,mm)
     &                          +w(1,m1,2)*q1d(i-2+m1,mm)
     &                          +w(1,m1,3)*q1d(i-1+m1,mm))
                 hgg(i,m1,2) = hgg(i,m1,2) + evl(i,m,mm) *
     &                         ( w(2,m1,1)*q1d(i-3+m1,mm)
     &                          +w(2,m1,2)*q1d(i-2+m1,mm)
     &                          +w(2,m1,3)*q1d(i-1+m1,mm))
                 hgg(i,m1,3) = hgg(i,m1,3) + evl(i,m,mm) *
     &                         ( w(3,m1,1)*q1d(i-3+m1,mm)
     &                          +w(3,m1,2)*q1d(i-2+m1,mm)
     &                          +w(3,m1,3)*q1d(i-1+m1,mm))

              enddo !meqn loop
           enddo !m1 loop

        enddo !i loop


c        # the reconstruction


c          # m1=1: construct ql

           im=1
           ione=im
           inone=-im
           intwo=-2*im
  
           do 20 i=0,mx
  
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
  
             hff(i,m,1) = s1*(hgg(i,1,1)-hgg(i,2,1)) + hgg(i,2,1)
     &                    + s3*(hgg(i,3,1)-hgg(i,2,1))

             d1 = w5(1,2)/w(2,1,1)
             d3 = w5(5,2)/w(2,3,3)
             d2 = 1-d1-d3

c  Negative weights here: special treatment
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
c         positive part
          s1 = d1pp * tt2 * tt3
          s2 = d2pp * tt1 * tt3
          s3 = d3pp * tt1 * tt2
          t0 = 1. / ( s1 + s2 + s3 )
          s1 = s1 * t0
          s3 = s3 * t0
c	    
          hff(i,m,2) = sump* ( s1*(hgg(i,1,2)-hgg(i,2,2)) 
     &                 + hgg(i,2,2) + s3*(hgg(i,3,2)-hgg(i,2,2)) )
c         negative part
          s1 = d1nn * tt2 * tt3
          s2 = d2nn * tt1 * tt3
          s3 = d3nn * tt1 * tt2
          t0 = 1. / ( s1 + s2 + s3 )
          s1 = s1 * t0
          s3 = s3 * t0
c	    
        hff(i,m,2) = hff(i,m,2) - sumn* (s1*(hgg(i,1,2)-hgg(i,2,2)) 
     &               + hgg(i,2,2) + s3*(hgg(i,3,2)-hgg(i,2,2)) )
c
          d1 = w5(1,3) / w(3,1,1)
          d3 = w5(5,3) / w(3,3,3)
          d2= 1 - d1 - d3
 
          s1 = d1 * tt2 * tt3
          s2 = d2 * tt1 * tt3
          s3 = d3 * tt1 * tt2
          t0 = 1. / ( s1 + s2 + s3 )
          s1 = s1 * t0
          s3 = s3 * t0
c	    
          hff(i,m,3) = s1*(hgg(i,1,3)-hgg(i,2,3)) 
     &                + hgg(i,2,3) + s3*(hgg(i,3,3)-hgg(i,2,3))


   20        continue
  100    continue
*----------------- loop in "m"  ends  here  -------------------

c Project to the physical space:

c       #Compute left and right eigenvectors
c      call evec(ixy,maxm,meqn,mwaves,mbc,mx,q,auxl,auxr,evl,evr)

        do ig =  1,3
c       ig=1 ---> -sqrt(0.6), ig=2 ---> 0, ig=3 ---> sqrt(0.6)
           do m =  1, meqn
              do i = 2-mbc,mx+mbc
                 q_gauss(i,m,ig)=0.d0
                 do mm=1,meqn 
                   q_gauss(i,m,ig) = q_gauss(i,m,ig) + 
     &                              evr(i,m,mm)*hff(i,mm,ig)

                 enddo ! loop over eqns (mm)
              enddo ! loop over cells (i)
           enddo ! loop over eqns (m)
        enddo !loop over gauss points (ig)

      return
      end
