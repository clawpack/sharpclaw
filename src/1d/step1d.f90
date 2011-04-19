! ===================================================================
subroutine step(q,g,dq,aux,dt,cfl,t,rp,src,tfluct)
! ===================================================================
    ! Program: SharpClaw
    ! File: 1d/step1d.f95    
    ! Author: David Ketcheson

    ! This is the one-dimensional step routine:
    ! Apply boundary conditions to q
    ! Apply source terms, if anyway
    ! Then compute the fluxes

    use ClawData
    use ClawParams
    implicit none
    type(griddat) g
    double precision, intent(out) :: dq,cfl
    double precision, intent(in) :: q,aux,dt,t
    external rp,src,tfluct

    call bc(maxnx,meqn,mbc,nx,xlower,dx,q,maux,aux,t,dt,mthbc)

    ! Evaluate source terms: Evaluate (delta t) * psi(q,t) and store in dq
    ! if there is no source term, the default routine simply sets dq=0
    ! BUT IT MUST STILL BE CALLED!
    call src(q,dq,aux,t,dt)

    !Now add the hyperbolic part to dq
    call flux1(q,g,dq,aux,dt,cfl,t,rp,tfluct,1)

end subroutine step
