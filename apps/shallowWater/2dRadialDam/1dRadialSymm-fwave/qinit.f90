! =============================================================
       subroutine qinit(ndim,meqn,mbc,nx,xlower,dx,q,maux,aux)
! =============================================================

     ! Set initial conditions for q.


      implicit double precision (a-h,o-z)
      integer :: ndim
      integer :: nx(ndim)
      double precision :: xlower(ndim), dx(ndim)

      double precision, intent(out) :: q(1-mbc:nx(1)+mbc, meqn)
      dimension aux(1-mbc:nx(1)+mbc, *)
      common/cdisc/ x0,y0,alf,beta,r0,idisc
      common /comic/ hin,hout

      pi = 4.d0*datan(1.d0)

      do i=1,nx(1)
        xcell = xlower(1) + (i-0.5d0)*dx(1)
        circBound = xlower(1)+r0
        if (xcell .lt. circBound) then
             h = hin
        else
             h = hout
        endif
        q(i,1) = h
        q(i,2) = 0.d0
      enddo

end subroutine qinit
