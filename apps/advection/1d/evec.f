
c
c
c     =====================================================
      subroutine evec(maxm,meqn,mbc,mx,q,auxl,auxr,evl,evr)
c     =====================================================
c
c     # Calculation of left and right eigenvectors
c       for the linear advection equation in 1d 
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
      common /comrp/ u
c

      do 20 i = 2-mbc, mx+mbc
         evr(i,1,1) = 1
         evl(i,1,1) = 1
   20    continue

      return
      end
