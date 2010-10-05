c     ============================================
      subroutine setaux(maxmx,mbc,nx,xlower,dx,maux,aux)
c     ============================================
c
c     # set auxiliary arrays 

c     #   aux(i,j,1) is edge velocity at "left" boundary of grid point (i,j)
c     #   aux(i,j,2) is edge velocity at "bottom" boundary of grid point (i,j)

c
c     
      implicit double precision (a-h,o-z)
      dimension nx(2),xlower(2),dx(2)
      dimension aux(1-mbc:nx(1)+mbc,1-mbc:nx(2)+mbc, maux)
c
c     # constant velocities which are used if tperiod=0 is specified
c     # in setprob.data

      do 20 i=1-mbc,nx(1)+mbc
         do 20 j=1-mbc,nx(2)+mbc

c           # coordinates of lower left corner of grid cell:
            xll = xlower(1) + (i-1)*dx(1)
            yll = xlower(2) + (j-1)*dx(2)

c           # difference stream function psi to get normal velocities:
            aux(i,j,1) = -(psi(xll, yll+dx(2)) - psi(xll,yll)) / dx(2)
            aux(i,j,2) =  (psi(xll+dx(1), yll) - psi(xll,yll)) / dx(1)
   20       continue

c
       return
       end
