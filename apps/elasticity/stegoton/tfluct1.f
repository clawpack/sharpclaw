
c
c
c     =====================================================
      subroutine tfluct(ixy,maxm,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,
     &                   s,adq)
c     =====================================================
c
c     # Riemann solver for the nonlinear elastic equations in 1d,
c     #  variable coefficients
c     #   eps_t - (m/rho(x))_x = 0
c     #   m_t - sigma(eps,x)_x =0
c     # where eps=strain, m=rho*u=momentum
c
c     # aux(i,1) = rho(i)
c     # function sigma(eps,i) gives stress-strain relation in i'th cell
c     # function sigmap(eps,i) gives d/(d eps) of sigma
c     #    For linear:   sigma(eps,i) = K_i * eps, and sigmap(eps,i) = K_i
c
c     # On input, ql contains the state vector at the left edge of each cell
c     #           qr contains the state vector at the right edge of each cell
c
c     # On output, fwave contains the waves as jumps in f,
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
c     # From the basic clawpack routines, this routine is called with ql = qr
c
c
      implicit double precision (a-h,o-z)
c
      external sigma,sigmap
      dimension auxl(1-mbc:maxm+mbc, 3)
      dimension auxr(1-mbc:maxm+mbc, 3)
      dimension fwave(1-mbc:maxm+mbc, meqn, mwaves)
      dimension    s(1-mbc:maxm+mbc, mwaves)
      dimension   ql(1-mbc:maxm+mbc, meqn)
      dimension   qr(1-mbc:maxm+mbc, meqn)
      dimension adq(1-mbc:maxm+mbc, meqn)
c
c     local arrays
c     ------------
      dimension delta(2)
c
c
c     # split the jump in q at each interface into waves
c
      do 20 i = 2-mbc, mx+mbc
         rhoi = auxr(i,1)
         rhoim = auxl(i,1)
         epsi = qr(i,1)
         epsim = ql(i,1)
         urhoi = qr(i,2)
         urhoim = ql(i,2)

c        #linearize on each side:

         bulki = sigmap(epsi,i,auxr(i,2),auxr(i,3))
         bulkim = sigmap(epsim,i,auxl(i,2),auxl(i,3))
         ci = dsqrt(bulki / rhoi)
         cim = dsqrt(bulkim / rhoim)
       	 zi = ci*rhoi
         zim = cim*rhoim

         du = urhoi/rhoi - urhoim/rhoim
         dsig = sigma(epsi,i,auxr(i,2),auxr(i,3)) 
     &         - sigma(epsim,i,auxl(i,2),auxl(i,3))
         b1 = -(zi*du + dsig) / (zim + zi)
         b2 = -(zim*du - dsig) / (zim + zi)
         a1 = b1 / (-cim)
         a2 = b2 / ci

         estarim = epsim + a1
         estari = epsi - a2

         ui = urhoi / rhoi
         uim = urhoim / rhoim

         ustar = urhoim/rhoim + (estarim - epsim) * cim
c               = urhoi/rhoi - (estari - epsi) * ci
c               = uim - b1
        
c
c        # Compute the waves.
c
         fwave(i,1,1) = b1
         fwave(i,2,1) = b1 * zim
         s(i,1) = -cim
c
         fwave(i,1,2) = b2
         fwave(i,2,2) = b2*(-zi)
         s(i,2) = ci

   20    continue
c
c
c     # compute the leftgoing and rightgoing fluctuations:
c     # Note s(i,1) < 0   and   s(i,2) > 0.
c
      do 220 m=1,meqn
         do 220 i = 2-mbc, mx+mbc
            amdq = fwave(i,m,1)
            apdq = fwave(i,m,2)
            adq(i,m) = amdq+apdq
  220       continue
c
      return
      end
