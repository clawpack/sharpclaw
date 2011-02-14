c =========================================================
      subroutine rp(ixy,maxmx,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,
     &		 wave,s,amdq,apdq)
c =========================================================
c
c # Solve Riemann problems for the 1D shallow water equations
c # (h)_t + (h*u)_x = 0
c # (hu)_t + (h*u^2 + 1/2*grav*h^2) = -grav*h*(b)_x
c # using q-wave algorithm and Roe's approximate Riemann solver 
c # with entropy fix for transonic rarefractions.  
c #
c # Note that with the q-wave approach the source term in the discharge equation 
c # is treated with a fractional step method. Therefore, the user should provide 
c # a subroutine which compute the contribution of the source term to the residual.
c # In this case, the aforementioned subroutine is src1.f90
c
c # On input, ql contains the state vector at the left edge of each cell
c #           qr contains the state vector at the right edge of each cell
c # On output, wave contains the waves, 
c #            s the speeds, 
c #            amdq the  left-going flux difference  A^- \Delta q
c #            apdq the right-going flux difference  A^+ \Delta q
c
c # Note that the i'th Riemann problem has left state qr(i-1,:)
c #                                    and right state ql(i,:)
c # From the basic clawpack routine step1, rp is called with ql = qr = q.
c
c # Here meqn=mwaves=2 should be passed from the calling routine

      implicit double precision (a-h,o-z)
      dimension   ql(1-mbc:maxmx+mbc, meqn)
      dimension   qr(1-mbc:maxmx+mbc, meqn)
      dimension    s(1-mbc:maxmx+mbc, mwaves)
      dimension wave(1-mbc:maxmx+mbc, meqn, mwaves)
      dimension amdq(1-mbc:maxmx+mbc, meqn)
      dimension apdq(1-mbc:maxmx+mbc, meqn)
      dimension auxl(1-mbc:maxmx+mbc, *)
      dimension auxr(1-mbc:maxmx+mbc, *)
c
c     # local storage
c     ---------------
      dimension delta(2)
      logical efix
      common /comrp/ grav
      data efix /.false./    !# use entropy fix for transonic rarefactions
          
c
      do 30 i=2-mbc,mx+mbc
c
c         
c     # compute  Roe-averaged quantities: 
         ubar = (qr(i-1,2)/dsqrt(qr(i-1,1)) + ql(i,2)/dsqrt(ql(i,1)))/
     .       ( dsqrt(qr(i-1,1)) + dsqrt(ql(i,1)) )
         cbar=dsqrt(0.5d0*grav*(qr(i-1,1) + ql(i,1)))
         
c     # delta(1)=h(i)-h(i-1) and  delta(2)=hu(i)-hu(i-1)
      delta(1) = ql(i,1) - qr(i-1,1)
      delta(2) = ql(i,2) - qr(i-1,2)

c # compute coeffs in the evector expansion of delta(1),delta(2)
      a1 = 0.5d0*(-delta(2) + (ubar + cbar) * delta(1))/cbar
      a2 = 0.5d0*( delta(2) - (ubar - cbar) * delta(1))/cbar

c     # finally, compute the waves.
         wave(i,1,1) = a1
         wave(i,2,1) = a1*(ubar - cbar)
         s(i,1) = ubar - cbar
         
         wave(i,1,2) = a2
         wave(i,2,2) = a2*(ubar + cbar)
         s(i,2) = ubar + cbar
         
   30 continue

c     # No entropy fix
c     ----------------------------------------------
c     # amdq = SUM s*wave   over left-going waves
c     # apdq = SUM s*wave   over right-going waves

      do 100 m=1,2
         do 100 i=2-mbc, mx+mbc
            amdq(i,m) = 0.d0
            apdq(i,m) = 0.d0
            do 90 mw=1,mwaves
               if (s(i,mw) .lt. 0.d0) then
                   amdq(i,m) = amdq(i,m) + s(i,mw)*wave(i,m,mw)
                 else
                   apdq(i,m) = apdq(i,m) + s(i,mw)*wave(i,m,mw)
                 endif
   90          continue
  100       continue
      go to 900


  900 continue
      return
      end
