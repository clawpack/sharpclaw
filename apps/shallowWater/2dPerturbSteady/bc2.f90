! ==========================================================
subroutine bc(meqn,mbc,nx,xlower,dx,q,maux,aux,t,dt,mthbc)
! ==========================================================
    ! Standard boundary condition choices for claw2

    ! At each boundary  k = 1 (left),  2 (right),  3 (top), 4 (bottom):
    !   mthbc(k) =  0  for user-supplied BC's (must be inserted!)
    !            =  1  for zero-order extrapolation
    !            =  2  for periodic boundary coniditions
    !            =  3  for solid walls, assuming this can be implemented
    !                  by reflecting the data about the boundary and then
    !                  negating the 2'nd (for k=1,2) or 3'rd (for k=3,4)
    !                  component of q.

    ! Extend the data from the interior cells (1:nx(1), 1:nx(2))
    ! to the ghost cells outside the region:
    !   (i, 1-jbc)   for jbc = 1,mbc,  i = 1-mbc, nx(1)+mbc
    !   (i, nx(2)+jbc)  for jbc = 1,mbc,  i = 1-mbc, nx(1)+mbc
    !   (1-ibc, j)   for ibc = 1,mbc,  j = 1-mbc, nx(2)+mbc
    !   (nx(1)+ibc, j)  for ibc = 1,mbc,  j = 1-mbc, nx(2)+mbc

    implicit double precision (a-h,o-z)
    integer :: nx(2)
    dimension q(1-mbc:nx(1)+mbc, 1-mbc:nx(2)+mbc, meqn)
    dimension aux(1-mbc:nx(1)+mbc, 1-mbc:nx(2)+mbc, *)
    dimension mthbc(4),xlower(2),dx(2)

!-------------------------------------------------------
!     left boundary:
!-------------------------------------------------------
      go to (100,110,120,130) mthbc(1)+1

  100 continue
        ! user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(1)=0 and no BCs specified in bc2'
      stop
      go to 199
!
  110 continue
        ! zero-order extrapolation:
      do 115 m=2,meqn
         do 115 ibc=1,mbc
            do 115 j = 1-mbc, nx(2)+mbc
               q(1-ibc,j,1) = q(1,j,1) + aux(1,j,1) - aux(1-ibc,j,1)
               q(1-ibc,j,m) = q(1,j,m)
  115       continue
      go to 199

  120 continue
        ! periodic:  
      do 125 m=1,meqn
         do 125 ibc=1,mbc
            do 125 j = 1-mbc, nx(2)+mbc
               q(1-ibc,j,m) = q(nx(1)+1-ibc,j,m)
  125       continue
      go to 199

  130 continue
        ! solid wall (assumes 2'nd component is velocity or momentum in x):
      do 135 m=1,meqn
         do 135 ibc=1,mbc
            do 135 j = 1-mbc, nx(2)+mbc
               q(1-ibc,j,m) = q(ibc,j,m)
  135       continue
        ! negate the normal velocity:
      do 136 ibc=1,mbc
         do 136 j = 1-mbc, nx(2)+mbc
            q(1-ibc,j,2) = -q(ibc,j,2)
  136    continue
      go to 199

  199 continue

!-------------------------------------------------------
!     # right boundary:
!-------------------------------------------------------
      go to (200,210,220,230) mthbc(2)+1

  200 continue
        ! user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(2)=0 and no BCs specified in bc2'
      stop
      go to 299

  210 continue
!     # zero-order extrapolation:
      do 215 m=2,meqn
         do 215 ibc=1,mbc
            do 215 j = 1-mbc, nx(2)+mbc
               q(nx(1)+ibc,j,1) = q(nx(1),j,1) + aux(nx(1),j,1) - aux(nx(1)+ibc,j,1)  
               q(nx(1)+ibc,j,m) = q(nx(1),j,m)
  215       continue
      go to 299
      

  220 continue
!     # periodic:  
      do 225 m=1,meqn
         do 225 ibc=1,mbc
            do 225 j = 1-mbc, nx(2)+mbc
               q(nx(1)+ibc,j,m) = q(ibc,j,m)
  225       continue
      go to 299

  230 continue
!     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do 235 m=1,meqn
         do 235 ibc=1,mbc
            do 235 j = 1-mbc, nx(2)+mbc
               q(nx(1)+ibc,j,m) = q(nx(1)+1-ibc,j,m)
  235       continue
!     # negate the normal velocity:
      do 236 ibc=1,mbc
         do 236 j = 1-mbc, nx(2)+mbc
            q(nx(1)+ibc,j,2) = -q(nx(1)+1-ibc,j,2)
  236    continue
      go to 299

  299 continue
!
!-------------------------------------------------------
!     # bottom boundary:
!-------------------------------------------------------
      go to (300,310,320,330) mthbc(3)+1
!
  300 continue
!     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(3)=0 and no BCs specified in bc2'
      stop
      go to 399
!
  310 continue
!     # zero-order extrapolation:
      do 315 m=2,meqn
         do 315 jbc=1,mbc
            do 315 i = 1-mbc, nx(1)+mbc
               q(i,1-jbc,1) = q(i,1,1) + aux(i,1,1) - aux(i,1-jbc,1)
               q(i,1-jbc,m) = q(i,1,m)
  315       continue
      go to 399

  320 continue
!     # periodic:  
      do 325 m=1,meqn
         do 325 jbc=1,mbc
            do 325 i = 1-mbc, nx(1)+mbc
               q(i,1-jbc,m) = q(i,nx(2)+1-jbc,m)
  325       continue
      go to 399

  330 continue
!     # solid wall (assumes 3'rd component is velocity or momentum in y):
      do 335 m=1,meqn
         do 335 jbc=1,mbc
            do 335 i = 1-mbc, nx(1)+mbc
               q(i,1-jbc,m) = q(i,jbc,m)
  335       continue
!     # negate the normal velocity:
      do 336 jbc=1,mbc
         do 336 i = 1-mbc, nx(1)+mbc
            q(i,1-jbc,3) = -q(i,jbc,3)
  336    continue
      go to 399

  399 continue
!
!-------------------------------------------------------
!     # top boundary:
!-------------------------------------------------------
      go to (400,410,420,430) mthbc(4)+1
!
  400 continue
!     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(4)=0 and no BCs specified in bc2'
      stop
      go to 499

  410 continue
!     # zero-order extrapolation:
      do 415 m=1,meqn
         do 415 jbc=1,mbc
            do 415 i = 1-mbc, nx(1)+mbc
               q(i,nx(2)+jbc,1) = q(i,nx(2),1) + aux(i,nx(2),1) - aux(i,nx(2)+jbc,1)
               q(i,nx(2)+jbc,m) = q(i,nx(2),m)
  415       continue
      go to 499

  420 continue
!     # periodic:  
      do 425 m=1,meqn
         do 425 jbc=1,mbc
            do 425 i = 1-mbc, nx(1)+mbc
               q(i,nx(2)+jbc,m) = q(i,jbc,m)
  425       continue
      go to 499

  430 continue
!     # solid wall (assumes 3'rd component is velocity or momentum in y):
      do 435 m=1,meqn
         do 435 jbc=1,mbc
            do 435 i = 1-mbc, nx(1)+mbc
               q(i,nx(2)+jbc,m) = q(i,nx(2)+1-jbc,m)
  435       continue
!     # negate the normal velocity:
      do 436 jbc=1,mbc
         do 436 i = 1-mbc, nx(1)+mbc
            q(i,nx(2)+jbc,3) = -q(i,nx(2)+1-jbc,3)
  436    continue
      go to 499

  499 continue

      return
      end
