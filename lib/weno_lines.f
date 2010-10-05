c
c ===================================================================
      subroutine weno_lines(ixy,q_gauss,ql,qr,auxl,auxr)
c ===================================================================

c This is a subroutine of fifth order WENO reconstruction, for CLAWpack.
c Adapted from code provided by Yulong Xing
c 2D version for fully multi-d reconstruction
c Using Gauss quadrature to calculate the waves                            
c works with weno5_gauss2()
c this part reconstructs along ngauss Gauss quadrature "lines"
c in a slice

c Takes q_gauss: 1D averages of q in each cell at gauss pts. in other
c                 dimension
c Returns ql,qr: Reconstructed point values on either side of each
c                 interface, at gauss quadrature points on that face

c      I am assuming that mbc=3 and all the q values from -2 to mx+3
c      are available

c very inefficiently coded for now due to adaptation from the
c version with char. decomp.

      implicit double precision (a-h,o-z)

      parameter (ngauss = 3)

      dimension q_gauss(1-mbc:maxnx+mbc,meqn,ngauss)
      dimension ql(1-mbc:maxnx+mbc,meqn,ngauss)
      dimension qr(1-mbc:maxnx+mbc,meqn,ngauss)

      epweno = 1.e-36

c     # loop over all equations (all components).  

c           # compute and store the differences of the cell averages
c           # along each gauss line
c           # Get the eigenvectors at the same time
      do ig = 1,ngauss
         do m=1,meqn
            do i=2-mbc,mx+mbc
               g%hff(i,m,ig)=q_gauss(i,m,ig)-q_gauss(i-1,m,ig)
            enddo !Loop over cells
         enddo !Loop over equations
      enddo !Loop over gauss points

      do 100 m=1,meqn

c       # No characteristic projection

        do 70 ig=1,ngauss

           do  i = 2-mbc,mx+mbc
              do m1 = -2,2
                g%hh(i,m1) = g%hff(i+m1,m,ig)
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
  
             ugg(i,m,m1,ig) = ( s1*(t2-t1) + (0.5*s3-0.25)*(t3-t2) ) /3.



   20          continue !Loop over cells (i)
   50        continue !Loop over interface side (m1)
   70      continue !Loop over gauss lines
  100    continue !Loop over equations (m)
*----------------- loop in "m"  ends  here  -------------------

        do ig =  1,3
c       ig=1 ---> -sqrt(0.6), ig=2 ---> 0, ig=3 ---> sqrt(0.6)
           do m1 = 1,2
              do m =  1, meqn
                 do i = 2-mbc,mx+mbc
              g%uhh(i,m,m1,ig)=(-q_gauss(i-2,m,ig)+7.d0*(q_gauss(i-1,m,ig)+
     &                        q_gauss(i,m,ig))-q_gauss(i+1,m,ig) )/12.d0

                      g%uhh(i,m,m1,ig) = g%uhh(i,m,m1,ig) + 
     &                                ugg(i,m,m1,ig)

                 enddo ! loop over cells (i)
              enddo ! loop over eqns (m)
           enddo ! loop over interface side (m1)
        enddo !loop over gauss points (ig)

      do m =  1, meqn
        do i=1,mx+1
          do ig =  1,3
            qr(i-1,m,ig)=g%uhh(i,m,1,ig)
            ql(i,m,ig)=g%uhh(i,m,2,ig)
          enddo
        enddo
      enddo
                
      return
      end
