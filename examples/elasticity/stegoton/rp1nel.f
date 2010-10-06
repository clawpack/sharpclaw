!     =====================================================
      subroutine rp(ixy,maxm,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,
     &                  fwave,s,amdq,apdq)
!     =====================================================
!
!     # Riemann solver for the nonlinear elastic equations in 1d,
!     #  variable coefficients
!     #   eps_t - (m/rho(x))_x = 0
!     #   m_t - sigma(eps,x)_x =0
!     # where eps=strain, m=rho*u=momentum
!
!     # aux(i,1) = rho(i)
!     # function sigma(eps,i) gives stress-strain relation in ith cell
!     # function sigmap(eps,i) gives d/(d eps) of sigma
!     #    For linear:   sigma(eps,i) = K_i * eps, and sigmap(eps,i) = K_i
!
!     # On input, ql contains the state vector at the left edge of each cell
!     #           qr contains the state vector at the right edge of each cell
!
!     # On output, fwave contains the waves as jumps in f,
!     #            s the speeds,
!     #
!     #            amdq = A^- Delta q, 
!     #            apdq = A^+ Delta q,
!     #                   the decomposition of the flux difference
!     #                       f(qr(i-1)) - f(ql(i))
!     #                   into leftgoing and rightgoing parts respectively.
!     #
!
!     # Note that the ith Riemann problem has left state qr(i-1,:)
!     #                                    and right state ql(i,:)
!     # From the basic clawpack routines, this routine is called with ql = qr
!
!
      implicit double precision (a-h,o-z)
!
      dimension auxl(1-mbc:maxm+mbc, 3)
      dimension auxr(1-mbc:maxm+mbc, 3)
      dimension fwave(1-mbc:maxm+mbc, meqn, mwaves)
      dimension    s(1-mbc:maxm+mbc, mwaves)
      dimension   ql(1-mbc:maxm+mbc, meqn)
      dimension   qr(1-mbc:maxm+mbc, meqn)
      dimension apdq(1-mbc:maxm+mbc, meqn)
      dimension amdq(1-mbc:maxm+mbc, meqn)
!
!     local arrays
!     ------------
      dimension delta(2)
!
!
!     # split the jump in q at each interface into waves
!
      do 20 i = 2-mbc, mx+mbc
         rhoi = auxl(i,1)
         rhoim = auxr(i-1,1)
         epsi = ql(i,1)
         epsim = qr(i-1,1)
         urhoi = ql(i,2)
         urhoim = qr(i-1,2)

!        #linearize on each side:

         bulki = sigmap(epsi,i,auxl(i,2),auxl(i,3))
         bulkim = sigmap(epsim,i-1,auxr(i-1,2),auxr(i-1,3))
         ci = dsqrt(bulki / rhoi)
         cim = dsqrt(bulkim / rhoim)
         zi = ci*rhoi
         zim = cim*rhoim

         du = urhoi/rhoi - urhoim/rhoim
         dsig = sigma(epsi,i,auxl(i,2),auxl(i,3)) 
     &          - sigma(epsim,i-1,auxr(i-1,2),auxr(i-1,3))
         b1 = -(zi*du + dsig) / (zim + zi)
         b2 = -(zim*du - dsig) / (zim + zi)
!         a1 = b1 / (-cim)
!         a2 = b2 / ci

!         estarim = epsim + a1
!         estari = epsi - a2

!         ui = urhoi / rhoi
!         uim = urhoim / rhoim

!         ustar = urhoim/rhoim + (estarim - epsim) * cim
!               = urhoi/rhoi - (estari - epsi) * ci
!               = uim - b1
        
!
!        # Compute the waves.
!
         fwave(i,1,1) = b1
         fwave(i,2,1) = b1 * zim
         s(i,1) = -cim
!
         fwave(i,1,2) = b2
         fwave(i,2,2) = b2*(-zi)
         s(i,2) = ci

   20    continue
!
!
!     # compute the leftgoing and rightgoing fluctuations:
!     # Note s(i,1) < 0   and   s(i,2) > 0.
!
      do 220 m=1,meqn
         do 220 i = 2-mbc, mx+mbc
            amdq(i,m) = fwave(i,m,1)
            apdq(i,m) = fwave(i,m,2)
  220       continue
!
      return
      end


!     --------------------------------------------
      double precision function sigma(eps,i,coefA,coefB)
!     --------------------------------------------
      implicit double precision (a-h,o-z)

!     # stress-strain relation in ith cell


!     # nonlinear in both layers:
!     sigma = (coefB + coefA)*(dexp(eps) - 1.d0)

!     # nonlinear layer 0:
!     sigma = coefB*eps + coefA*(dexp(eps) - 1.d0)

!     # nonlinear layer 1:
!     sigma = coefB*(dexp(eps) - 1.d0) + coefA*eps

      if (coefA.gt.0.d0) then
!         # layer A:
!         sigma = coefA*eps 
          sigma = dexp(coefA*eps) - 1.d0
!          sigma = coefA*eps + 0.5d0*eps**2
        else
!         # layer B:
          sigma = coefB*eps 
!         sigma = dexp(coefB*eps) - 1.d0
        endif

      return
      end


!     --------------------------------------------
      double precision function sigmap(eps,i,coefA,coefB)
!     --------------------------------------------
      implicit double precision (a-h,o-z)

!     # derivative of stress-strain relation in ith cell

!     # nonlinear in both layers:
!     sigmap = (coefB + coefA)*dexp(eps)

!     # nonlinear layer 0:
!     sigmap = coefB + coefA*dexp(eps)

!     # nonlinear layer 1:
!     sigmap = coefB*dexp(eps) + coefA

      if (coefA.gt.0.d0) then
!         # layer A:
!         sigmap = coefA
          sigmap = coefA*dexp(coefA*eps)
!         sigmap = coefA + eps
        else
!         # layer B:
          sigmap = coefB
!         sigmap = coefB*dexp(coefB*eps)
        endif


      return
      end
