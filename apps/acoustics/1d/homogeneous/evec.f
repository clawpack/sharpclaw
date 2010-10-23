
c
c
c     =====================================================
      subroutine evec(maxm,meqn,mbc,mx,q,auxl,auxr,evl,evr)
c     =====================================================
c
c     # Calculation of left and right eigenvectors
c       for the acoustics equations in 1d, with
c     # variable coefficients (heterogeneous media)
c
c     # auxl(i,1) should contain the impedance Z in cell i
c     # auxl(i,2) should contain the sound speed c in cell i
c
c     # On output, evl(i) and evr(i) contain left/right eigenvectors
c     # at interface i-1/2
c
c
      implicit double precision (a-h,o-z)
c
      dimension  auxl(1-mbc:maxm+mbc, 2)
      dimension  auxr(1-mbc:maxm+mbc, 2)
      dimension  q(1-mbc:maxm+mbc, meqn)
      dimension evl(1-mbc:maxm+mbc,meqn,meqn)
      dimension evr(1-mbc:maxm+mbc,meqn,meqn)
      common /cparam/ rho,bulk,cc,zz
c
c      -Z_i-1  Z_i
c R =
c        1      1

      do 20 i = 2-mbc, mx+mbc
         evr(i,1,1) = -zz
         evr(i,1,2) = zz
         evr(i,2,1) = 1.d0
         evr(i,2,2) = 1.d0
         z=1.d0/(zz+zz)
         evl(i,1,1) = -z
         evl(i,1,2) = zz*z
         evl(i,2,1) = 1.d0*z
         evl(i,2,2) = zz*z
   20    continue

      return
      end
