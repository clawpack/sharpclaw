c
c
c =========================================================
      subroutine tfluct(ixy,maxmx,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,
     &           s,adq)
c =========================================================
c
c     # find total fluctucations for the 1D advection equation q_t + u*q_x = 0.
c     # For constant advection velocity u, passed in common block.
c
c     # The advection speed u is passed in the common block comrp
c     # On input, ql contains the state vector at the left edge of each cell
c     #           qr contains the state vector at the right edge of each cell
c     # On output, wave contains the waves,
c     #            s the speeds,
c     #            amdq the  left-going flux difference  A^- \Delta q
c     #            apdq the right-going flux difference  A^+ \Delta q
c
c     # Note that the i'th Riemann problem has left state ql(i,:)
c     #                                    and right state qr(i,:)
c
      implicit double precision (a-h,o-z)
      dimension   ql(1-mbc:maxmx+mbc, meqn)
      dimension   qr(1-mbc:maxmx+mbc, meqn)
      dimension    s(1-mbc:maxmx+mbc, mwaves)
      dimension wave(1-mbc:maxmx+mbc, meqn, mwaves)
      dimension adq(1-mbc:maxmx+mbc, meqn)
      common /comrp/ u
c
c
c
      do 30 i=2-mbc,mx+mbc
c
c        # Compute the wave and speed
c
         wave(i,1,1) = qr(i,1) - ql(i,1)
         s(i,1) = u
         amdq = dmin1(u, 0.d0) * wave(i,1,1)
         apdq = dmax1(u, 0.d0) * wave(i,1,1)
         adq(i,1) = amdq + apdq
   30    continue
c
      return
      end
