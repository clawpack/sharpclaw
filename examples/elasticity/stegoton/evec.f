c
c
c     =====================================================
      subroutine evec(maxm,meqn,mbc,mx,q,auxl,auxr,evl,evr)
c     =====================================================
c
c     # Calculation of left and right eigenvectors
c     # for the nonlinear elastic equations in 1d, with
c     #  variable coefficients
c
c     # aux(i,1) = rho(i)
c
c     # function sigma(eps,i) gives stress-strain relation in i'th cell
c     # function sigmap(eps,i) gives d/(d eps) of sigma
c     #    For linear:   sigma(eps,i) = K_i * eps, and sigmap(eps,i) = K_i

c     # On input, q contains the state vector in each cell
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
c
      do 20 i = 2-mbc, mx+mbc
         rhoi = auxl(i,1)
         rhoim = auxr(i-1,1)
         epsi = q(i,1)
         epsim = q(i-1,1)
         urhoi = q(i,2)
         urhoim = q(i-1,2)
                  
c #linearize on each side:

         bulki = sigmap(epsi,i)
         bulkim = sigmap(epsim,i-1)
         ci = dsqrt(bulki / rhoi)
         cim = dsqrt(bulkim / rhoim)
         zi = ci*rhoi
         zim = cim*rhoim

c  Right eigenvector matrix
c        1      1
c R =
c        Z_i-1  -Z_i

         evl(i,1,1) = 1.d0
         evl(i,1,2) = 1.d0
         evl(i,2,1) = zim
         evl(i,2,2) = -zi

c  Left eigenvector matrix
         zh=1.d0/(zim+zi);

         evr(i,1,1) = zi*zh
         evr(i,1,2) = zh
         evr(i,2,1) = zim*zh
         evr(i,2,2) = -zh
   20    continue

      return
      end
