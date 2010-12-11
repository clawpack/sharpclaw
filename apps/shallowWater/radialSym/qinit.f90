! =========================================================
       subroutine qinit(maxmx,meqn,mbc,mx,xlower,dx,q,maux,aux)
! =========================================================

     ! Set initial conditions for q.


      implicit double precision (a-h,o-z)
      integer, parameter :: ndim=1
      integer :: mx(ndim)
      double precision :: xlower(ndim), dx(ndim)

      double precision, intent(out) :: q(1-mbc:mx(1)+mbc, meqn)
      dimension aux(1-mbc:mx(1)+mbc, *)
      common/cdisc/ x0,y0,alf,beta,r0,idisc
      common /comic/ hin,hout

      pi = 4.d0*datan(1.d0)
      width = 0.2d0

     !write(*,*) r0,hin,hout
      do i=1,mx(1)
        xcell = xlower(1) + (i-0.5d0)*dx(1)
        if (xcell .lt. r0) then
             h = hin
        else
             h = hout
        endif
        q(i,1) = h
        q(i,2) = 0.d0
      enddo

end subroutine qinit
