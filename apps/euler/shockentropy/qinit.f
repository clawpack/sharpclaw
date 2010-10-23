c
c
c =========================================================
       subroutine qinit(maxmx,meqn,mbc,mx,xlower,dx,q,maux,aux)
c =========================================================
c
c     # Set initial conditions for q.
c     # Smooth entropy wave hitting a shock
c
c
      implicit double precision (a-h,o-z)
      dimension q(1-mbc:maxmx+mbc, meqn)
      dimension aux(1-mbc:maxmx+mbc, *)
      common /param/ gamma, gamma1
c
c
       sloc = -4.d0
c
c      # data in left state:
       rhol = 3.857143d0
       ul = 2.629369d0
       pl = 10.333333d0
       rhoul = rhol*ul
       el = pl/gamma1 + 0.5d0*rhol*ul**2

c      # data in right state:
       ur = 0.d0
       pr = 1.d0
c

      do 150 i=1,mx
             xcell = xlower + (i-0.5d0)*dx
             if (xcell .lt. sloc) then
                 q(i,1) = rhol
                 q(i,2) = rhoul
                 q(i,3) = el
                else
                 q(i,1) = 1.d0 + 0.2d0*sin(5.d0*xcell)
                 q(i,2) = q(i,1)*ur
                 q(i,3) = pr/gamma1 + 0.5d0*q(i,1)*ur**2
                endif

  150    continue
c
      return
      end
