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
!      -Z_i-1  Z_i
! R =
!        1      1

      do 20 i = 2-mbc, mx+mbc
         evr(i,1,1) = -auxl(i-1,1)
         evr(i,1,2) = auxl(i,1)
         evr(i,2,1) = 1.d0
         evr(i,2,2) = 1.d0
         z=1.d0/(auxl(i-1,1)+auxl(i,1))
         evl(i,1,1) = -1.d0*z
         evl(i,1,2) = auxl(i,1)*z
         evl(i,2,1) = 1.d0*z
         evl(i,2,2) = auxl(i-1,1)*z
   20    continue

      return
      end
