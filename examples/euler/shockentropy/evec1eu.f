!     =====================================================
      subroutine evec(maxm,meqn,mbc,mx,q,auxl,auxr,evl,evr)
!     =====================================================
!
!     # Calculation of left and right eigenvectors
!       for the acoustics equations in 1d, with
!     # variable coefficients (heterogeneous media)
!
!     # auxl(i,1) should contain the impedance Z in cell i
!     # auxl(i,2) should contain the sound speed c in cell i
!
!     # On output, evl(i) and evr(i) contain left/right eigenvectors
!     # at interface i-1/2
!
!
      implicit double precision (a-h,o-z)
!
      dimension  auxl(1-mbc:maxm+mbc, 2)
      dimension  auxr(1-mbc:maxm+mbc, 2)
      dimension evl(1-mbc:maxm+mbc,meqn,meqn)
      dimension evr(1-mbc:maxm+mbc,meqn,meqn)
!
      parameter (max2 = 2002)  !# assumes at most 2000 grid points with mbc=2
      dimension u(-1:max2),enth(-1:max2),a(-1:max2)
!
!     # Compute Roe-averaged quantities:
!
      do 20 i=2-mbc,mx+mbc
         rhsqrtl = dsqrt(q(i-1,1))
         rhsqrtr = dsqrt(q(i,1))
         pl = gamma1*(q(i-1,3) - 0.5d0*(q(i-1,2)**2)/q(i-1,1))
         pr = gamma1*(q(i,3) - 0.5d0*(q(i,2)**2)/q(i,1))
         rhsq2 = rhsqrtl + rhsqrtr
         u(i) = (q(i-1,2)/rhsqrtl + q(i,2)/rhsqrtr) / rhsq2
         enth(i) = (((q(i-1,3)+pl)/rhsqrtl
     &             + (q(i,3)+pr)/rhsqrtr)) / rhsq2
         a2 = gamma1*(enth(i) - .5d0*u(i)**2)
         a(i) = dsqrt(a2)
   20    continue

!        1      1       1
!
! R =   u-c     u      u+c
!
!       H-uc   u^2/2   H+uc


      do 20 i = 2-mbc, mx+mbc
         evr(i,1,1) = 1.d0
         evr(i,1,2) = 1.d0
         evr(i,1,3) = 1.d0
         evr(i,2,1) = u(i)-a(i)
         evr(i,2,2) = u(i)
         evr(i,2,3) = u(i)+a(i)
         evr(i,3,1) = H(i)-u(i)*a(i)
         evr(i,3,2) = 0.5d0*u(i)**2
         evr(i,3,3) = H(i)+u(i)*a(i)

! Need to fill in evl!!!
   20    continue

      return
      end
