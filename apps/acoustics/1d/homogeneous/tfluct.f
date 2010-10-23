
c
c
c     =====================================================
      subroutine tfluct(ixy,maxmx,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,
     &                  s,adq)
c     =====================================================
c
c     # Total fluctuation solver for the acoustics equations in 1d,
c
c     # On input, ql contains the state vector at the left edge of each cell
c     #           qr contains the state vector at the right edge of each cell
c
c     # On output, wave contains the waves,
c     #            s the speeds,
c     #
c     #            amdq = A^- Delta q, 
c     #            apdq = A^+ Delta q,
c     #                   the decomposition of the flux difference
c     #                       f(qr(i-1)) - f(ql(i))
c     #                   into leftgoing and rightgoing parts respectively.
c     #
c
c     # Note that the i'th Riemann problem has left state ql(i,:)
c     #                                    and right state qr(i,:)
c
      implicit double precision (a-h,o-z)
c
      dimension wave(1-mbc:maxmx+mbc, meqn, mwaves)
      dimension    s(1-mbc:maxmx+mbc, mwaves)
      dimension   ql(1-mbc:maxmx+mbc, meqn)
      dimension   qr(1-mbc:maxmx+mbc, meqn)
      dimension adq(1-mbc:maxmx+mbc, meqn)
c
c     local arrays
c     ------------
      dimension delta(2)
c
c     # density, bulk modulus, and sound speed, and impedence of medium:
c     # (should be set in setprob.f)
      common /cparam/ rho,bulk,cc,zz   
c
c
c     # split the jump in q at each interface into waves
c
c     # find a1 and a2, the coefficients of the 2 eigenvectors:
      do 20 i = 2-mbc, mx+mbc
         delta(1) = qr(i,1) - ql(i,1)
         delta(2) = qr(i,2) - ql(i,2)
         a1 = (-delta(1) + zz*delta(2)) / (2.d0*zz)
         a2 =  (delta(1) + zz*delta(2)) / (2.d0*zz)
c
c        # Compute the waves.
c
         wave(i,1,1) = -a1*zz
         wave(i,2,1) = a1
         s(i,1) = -cc
c
         wave(i,1,2) = a2*zz
         wave(i,2,2) = a2
         s(i,2) = cc
c
   20    continue
c
c
c     # compute the leftgoing and rightgoing flux differences:
c     # Note s(i,1) < 0   and   s(i,2) > 0.
c
      do 220 m=1,meqn
         do 220 i = 2-mbc, mx+mbc
            amdq = s(i,1)*wave(i,m,1)
            apdq = s(i,2)*wave(i,m,2)
            adq(i,m) = amdq+apdq
  220       continue
c
      return
      end
