! =============================================================================
  subroutine tfluct(ixy,maxmx,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,s,adq)
! =============================================================================
!
! # This subroutine computes the total fluctuation (wave approach)
! # for the the non-linear traffic flow equation 
! # q_t + (u_max(x)*q*(1-q))_x = 0,
! # for the f-wave algorithm.
!
! # Note that u_max=u_max(x). Thus, the problem has a spatially varying flux
! # function. 
! 
!
! # On input, ql contains the state vector at the left edge of each cell
! #           qr contains the state vector at the right edge of each cell
!
! # On output, s contains the f-wave speeds, 
! #            adq the total fluctuation, i.e. the flux difference 
! #            expressed in terms of flux-waves.
!
! # Note that the i'th Riemann problem has left state ql(i,:)
! #                                    and right state qr(i,:)

  implicit none
  
  integer :: ixy, maxmx, meqn, mwaves, mbc
  integer :: mx(1)
  double precision :: ql(1-mbc:maxmx+mbc, meqn)
  double precision :: qr(1-mbc:maxmx+mbc, meqn)
  double precision :: s(1-mbc:maxmx+mbc, mwaves)
  double precision :: fwave(1-mbc:maxmx+mbc, meqn, mwaves)
  double precision :: adq(1-mbc:maxmx+mbc, meqn)
  double precision :: auxl(1-mbc:maxmx+mbc, 1)
  double precision :: auxr(1-mbc:maxmx+mbc, 1)
   
  integer :: i
  double precision :: fim1, fi, sim1, si, dq, f0
  
  do i=2-mbc,mx(1)+mbc

  	! # Flux in each cell and flux difference
    fim1 = auxl(i,1)*ql(i,1)*(1.d0 - ql(i,1))
    fi = auxr(i,1)*qr(i,1)*(1.d0 - qr(i,1))
    fwave(i,1,1) = fi - fim1

    ! # The total fluctuation for this case is simply the flux difference! 
    adq(i,1) = fwave(i,1,1)
    
  enddo
  
  return
  end subroutine tfluct

