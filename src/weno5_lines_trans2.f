c
c ===================================================================
      subroutine weno5_lines_trans2(ixy,maxm,meqn,mbc,mx,q_gauss,ql,qr,
     &                      auxl,auxr,evec)
c ===================================================================

c This is a subroutine of fifth order WENO reconstruction, for CLAWpack.
c Adapted from code provided by Yulong Xing
c Uses characteristic decomposition 
c 2D version for fully multi-d reconstruction
c Using Gauss quadrature to calculate the waves                            
c works with weno5_char_gauss2()
c this part reconstructs along ngauss Gauss quadrature "lines"
c in a slice

c Takes q_gauss: 1D averages of q in each cell at gauss pts. in other
c                 dimension
c Returns ql,qr: Reconstructed point values on either side of each
c                 interface, at gauss quadrature points on that face

c      I am assuming that mbc=3 and all the q values from -2 to mx+3
c      are available

      implicit double precision (a-h,o-z)

      external  evec
      parameter (ngauss = 3)

      dimension q_gauss(1-mbc:maxm+mbc,meqn,ngauss)
      dimension ql(1-mbc:maxm+mbc,meqn,ngauss)
      dimension qr(1-mbc:maxm+mbc,meqn,ngauss)
      dimension evl(1-mbc:maxm+mbc,meqn,meqn)
      dimension evr(1-mbc:maxm+mbc,meqn,meqn)

      parameter (mbcweno = 3, mxweno=10000,meqnweno=6)
      dimension qg(1-mbcweno:mxweno+mbcweno,meqnweno)
      dimension du(1-mbcweno:mxweno+mbcweno,meqnweno,ngauss)
      dimension gg(1-mbcweno:mxweno+mbcweno,meqnweno,ngauss)
      dimension uh(1-mbcweno:mxweno+mbcweno,meqnweno,2,ngauss)
      dimension  u(1-mbcweno:mxweno+mbcweno,meqnweno,2,ngauss)
      dimension hh(1-mbcweno:mxweno+mbcweno,-2:2)
      dimension gau_evl(1-mbcweno:mxweno+mbcweno,
     &                  meqnweno,meqnweno,ngauss)
      dimension gau_evr(1-mbcweno:mxweno+mbcweno,
     &                  meqnweno,meqnweno,ngauss)


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

c           # compute and store the differences of the cell averages
c           # along each gauss line
c           # Get the eigenvectors at the same time
      do ig = 1,ngauss
         do m=1,meqn
            do i=2-mbc,mx+mbc
               du(i,m,ig)=q_gauss(i,m,ig)-q_gauss(i-1,m,ig)
               qg(i,m)=q_gauss(i,m,ig)
            enddo !Loop over cells
         enddo !Loop over equations

         call evec(ixy,maxm,meqn,mwaves,mbc,mx,qg,auxl,auxr,evl,evr)

         do m=1,meqn
           do mm=1,meqn
             do i=2-mbc,mx+mbc
               gau_evl(i,m,mm,ig)=evl(i,m,mm)
               gau_evr(i,m,mm,ig)=evr(i,m,mm)
             enddo !Loop over cells
           enddo
         enddo !Loop over equations
 
      enddo !Loop over gauss points

c       # Find waves at each interface
c       # 'm'th characteristic field
      do 99 m=1,meqn
        do i=2-mbc,mx+mbc
          do ig=1,ngauss
          gg(i,m,ig)=0.d0
          do mm=1,meqn
            gg(i,m,ig) = gg(i,m,ig) + gau_evl(i,m,mm,ig)*du(i,m,ig)
          enddo
          enddo
        enddo
  99  continue

      do 100 m=1,meqn

c       # Project the waves to the
c       # 'm'th characteristic field

        do 70 ig=1,ngauss

           do  i = 2-mbc,mx+mbc
              do m1 = -2,2
                 hh(i,m1) = 0.d0
                 do mm=1,meqn !DK Maybe meqn should be mwave?
                   hh(i,m1) = hh(i,m1)+ gau_evl(i,m,mm,ig)*
     &                   gg(i+m1,m,ig)*gau_evr(i+m1,mm,m,ig)
                 enddo !meqn loop
              enddo !m1 loop
           enddo !i loop


c        # the reconstruction
c        # m1=1: construct ql
c        # m1=2: construct qr

         do 50 m1=1,2
        
           im=(-1)**(m1+1)
           ione=im
           inone=-im
           intwo=-2*im
  
           do 20 i=1,mx+1
  
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
  
             u(i,m,m1,ig) = ( s1*(t2-t1) + (0.5*s3-0.25)*(t3-t2) ) /3.



   20          continue !Loop over cells (i)
   50        continue !Loop over interface side (m1)
   70      continue !Loop over gauss lines
  100    continue !Loop over equations (m)
*----------------- loop in "m"  ends  here  -------------------

c Project to the physical space:

        do ig =  1,3
c       ig=1 ---> -sqrt(0.6), ig=2 ---> 0, ig=3 ---> sqrt(0.6)
           do m1 = 1,2
              do m =  1, meqn
                 do i = 2-mbc,mx+mbc
              uh(i,m,m1,ig)=(-q_gauss(i-2,m,ig)+7.d0*(q_gauss(i-1,m,ig)+
     &                        q_gauss(i,m,ig))-q_gauss(i+1,m,ig) )/12.d0

                    do mm=1,meqn !DK Maybe meqn should be mwave?
                      uh(i,m,m1,ig) = uh(i,m,m1,ig) + 
     &                                gau_evr(i,m,mm,ig)*u(i,mm,m1,ig)
                    enddo ! loop over eqns (mm)

                 enddo ! loop over cells (i)
              enddo ! loop over eqns (m)
           enddo ! loop over interface side (m1)
        enddo !loop over gauss points (ig)

      do m =  1, meqn
        do i=1,mx+1
          do ig =  1,3
            qr(i-1,m,ig)=uh(i,m,1,ig)
            ql(i,m,ig)=uh(i,m,2,ig)
          enddo
        enddo
      enddo
                
      return
      end
