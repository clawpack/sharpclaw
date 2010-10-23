! =================================================================
subroutine bc(maxmx,meqn,mbc,mx,xlower,dx,q,maux,aux,t,dt,mthbc)
! =================================================================
!
!     # Standard boundary condition choices for wclaw1
!
!     # At each boundary  k = 1 (left),  2 (right):
!     #   mthbc(k) =  0  for user-supplied BC's (must be inserted!)
!     #            =  1  for zero-order extrapolation
!     #            =  2  for periodic boundary conditions
!     #            =  3  for solid walls, assuming this can be implemented
!     #                  by reflecting the data about the boundary and then
!     #                  negating the 2'nd component of q.
!     ------------------------------------------------
!
!     # Extend the data from the computational region
!     #      i = 1, 2, ..., mx2
!     # to the virtual cells outside the region, with
!     #      i = 1-ibc  and   i = mx+ibc   for ibc=1,...,mbc
!
      implicit double precision (a-h,o-z)
      dimension q(1-mbc:mx+mbc, meqn)
      dimension aux(1-mbc:mx+mbc, *)

      dimension mthbc(2)

!
!
!-------------------------------------------------------
!     # left boundary:
!-------------------------------------------------------
      go to (100,110,120,130) mthbc(1)+1
!
  100 continue
!     # user-specified boundary conditions go here in place of error output
      a1=1.d0
      t1=1.d0
      tw1=1.d0
      do 105 m=1,meqn
         do 105 ibc=1,mbc
               q(1-ibc,m) = q(ibc,m)
  105       continue
!     # wall velocity:  (make negative for expansion)
      vwall = -a1*g0((t-t1)/tw1)
!     # adjust the normal momentum:
      do 106 ibc=1,mbc
            q(1-ibc,2) = vwall
            q(1-ibc,1) = 0.0d0*vwall*aux(ibc,1) - q(ibc,2)
  106    continue
      go to 199
!
  110 continue
!     # zero-order extrapolation:
      do 115 m=1,meqn
         do 115 ibc=1,mbc
               q(1-ibc,m) = q(1,m)
  115       continue
      go to 199

  120 continue
!     # periodic:  
      do 125 m=1,meqn
         do 125 ibc=1,mbc
               q(1-ibc,m) = q(mx+1-ibc,m)
  125       continue
      go to 199

  130 continue
!     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do 135 m=1,meqn
         do 135 ibc=1,mbc
               q(1-ibc,m) = q(ibc,m)
  135       continue
!     # negate the normal velocity:
      do 136 ibc=1,mbc
            q(1-ibc,2) = -q(ibc,2)
  136    continue
      go to 199

  199 continue

!
!-------------------------------------------------------
!     # right boundary:
!-------------------------------------------------------
      go to (200,210,220,230) mthbc(2)+1
!
  200 continue
!     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(2)=0 and no BCs specified in bc2'
      stop
      go to 299

  210 continue
!     # zero-order extrapolation:
      do 215 m=1,meqn
         do 215 ibc=1,mbc
               q(mx+ibc,m) = q(mx,m)
  215       continue
      go to 299

  220 continue
!     # periodic:  
      do 225 m=1,meqn
         do 225 ibc=1,mbc
               q(mx+ibc,m) = q(ibc,m)
  225       continue
      go to 299

  230 continue
!     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do 235 m=1,meqn
         do 235 ibc=1,mbc
               q(mx+ibc,m) = q(mx+1-ibc,m)
  235       continue
      do 236 ibc=1,mbc
            q(mx+ibc,2) = -q(mx+1-ibc,2)
  236    continue
      go to 299

  299 continue
!
      return
      end

!     ===============================
      double precision function g0(t)
!     ===============================

      implicit double precision (a-h,o-z)

      if (dabs(t) .lt. 1.d0)  then
        g0 = dcos(pi*t)
      else
        g0 = 0.d0
      endif

      return
      end

