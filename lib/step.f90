! ===================================================================
subroutine step(q,g,dq,aux,dt,cfl,t,rp,src,tfluct)
! ===================================================================
    ! Program: SharpClaw
    ! File: step.f95    
    ! Author: David Ketcheson

    ! Apply boundary conditions to q
    ! Then call appropriate subroutine to compute (delta t) * dq/dt
    use ClawData
    use ClawParams
    implicit none
    type(griddat) g
    double precision, intent(out) :: dq,cfl
    double precision, intent(in) :: q,aux,dt,t
    external rp,src,tfluct

    select case(ndim)
        case(1)
            call bc(maxnx,meqn,mbc,nx,xlower,dx,q,maux,aux,t,dt,mthbc)
        case(2)
            call bc(meqn,mbc,nx,xlower,dx,q,maux,aux,t,dt,mthbc)
    end select

    ! Evaluate source terms: Evaluate (delta t) * psi(q,t) and store in dq
    ! if there is no source term, the default routine simply sets dq=0
    ! BUT IT MUST STILL BE CALLED!
    call src(q,dq,aux,t,dt)

    !Now add the hyperbolic part to dq
    select case(ndim)
        case(1)
            call step1(q,g,dq,aux,dt,cfl,t,rp,tfluct,1)
        case(2)
            select case(multid_recon)
                case(0)
                    call step2(q,g,dq,aux,dt,cfl,t,rp,tfluct)
                case(1)
                    ! call step2md(q,g,dq,aux,dt,cfl,t,rp,tfluct)
            end select
    end select

end subroutine step
