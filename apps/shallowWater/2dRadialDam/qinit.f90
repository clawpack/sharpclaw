!     ========================================================
       subroutine qinit(ndim,meqn,mbc,nx,xlower,dx,q,maux,aux)
!     ========================================================
!
!     # Set initial conditions for q.
!     # Shallow water with radial dam break problem, h = hin inside
!     # circle specified in fdisc.f
!
       implicit double precision (a-h,o-z)
       integer :: ndim
       integer, intent(in) :: nx(ndim)
       double precision :: xlower(ndim), dx(ndim)
       dimension q(1-mbc:nx(1)+mbc, 1-mbc:nx(2)+mbc, meqn)
       common /comic/ hin,hout
       
       do 20 i=1,nx(1)
          xlow = xlower(1) + (i-1.d0)*dx(1)
          do 20 j=1,nx(2)
             ylow = xlower(2) + (j-1.d0)*dx(2)
	     call cellave(xlow,ylow,dx(1),dx(2),win)
	     q(i,j,1) = hin*win + hout*(1.d0-win)
	     q(i,j,2) = 0.d0
	     q(i,j,3) = 0.d0
  20         continue
       return
       end
