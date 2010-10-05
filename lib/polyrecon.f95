! Non-limited polynomial reconstruction routines
! For linear wave propagation problems

! ===================================================================
      subroutine q2qlqr_poly(q,ql,qr,mx)
! ===================================================================

! Reconstruction based on upwind-centered differencing

      use ClawData
      use ClawParams

      implicit none

      double precision ::  q(1-mbc:mx+mbc,meqn), &
                          ql(1-mbc:mx+mbc,meqn), &
                          qr(1-mbc:mx+mbc,meqn)
      double precision :: du(1-mbc:mx+mbc,meqn)
      double precision :: theta(-3:3), c(-3:2)
      integer i,j,k,mw,m,mx
      double precision :: ur,ul,invnorm

      k=(mthlim(1)-1)/2

      select case (mthlim(1))
        case (3)
          c(-1)=1.d0/6.d0
          c( 0)=2.d0/6.d0
        case (7)
          c(-3)=  3.d0/420.d0
          c(-2)=-22.d0/420.d0
          c(-1)= 79.d0/420.d0
          c( 0)=180.d0/420.d0
          c( 1)=-34.d0/420.d0
          c( 2)=  4.d0/420.d0
        case default
          write(*,*) 'ERROR: Invalid reconstruction order'
          stop
      end select

      forall (i=2-mbc:mx+mbc, m=1:meqn)
        du(i,m)=q(i,m)-q(i-1,m)
      end forall

      do 20 i=k+1,mx-k-1

        do m=1,meqn
          ! Stencil-independent part of reconstruction
          qr(i-1,m) = q(i-1,m)
          ql(i  ,m) = q(i  ,m)
        enddo

        do m=1,meqn ! loop over equations
  
          do j=-k,k
            theta(j)=du(i+j,m)
          enddo
  
          ul=0.d0
          ur=0.d0

          do j=-k,k-1
            ur = ur + c(j)*theta( j)
            ul = ul + c(j)*theta(-j)
            !write(*,*) i,j,theta(j),du(i+j,m),ur
          enddo
 
          ! reconstruct q^-_\imh
          qr(i-1,m) = qr(i-1,m) + ur
          ! reconstruct q^+_\imh
          ql(i  ,m) = ql(i  ,m) - ul

        enddo !loop over equations

   20     continue  !loop over interfaces


      return
      end


! ===================================================================
      subroutine q2qlqr_poly_fwave(q,ql,qr,wave,s,mx)
! ===================================================================

! Reconstruction based on upwind-centered differencing
! With fwave-slope reconstruction

      use ClawData
      use ClawParams

      implicit none

      double precision ::  q(1-mbc:mx+mbc,meqn), &
                          ql(1-mbc:mx+mbc,meqn), &
                          qr(1-mbc:mx+mbc,meqn)
      double precision :: wave(1-mbc:mx+mbc,meqn,mwaves)
      double precision :: s(   1-mbc:mx+mbc,       mwaves)
      double precision :: theta(-3:3), c(-3:2)
      integer i,j,k,mw,m,mx
      double precision :: ur,ul,invwavenorm

      rord=7  !Order of finite differencing for reconstruction
                !This should be a run-time parameter
      k=(rord-1)/2

      select case (rord)
        case (3)
          c(-1)=1.d0/6.d0
          c( 0)=2.d0/6.d0
        case (7)
          c(-3)=  3.d0/420.d0
          c(-2)=-22.d0/420.d0
          c(-1)= 79.d0/420.d0
          c( 0)=180.d0/420.d0
          c( 1)=-34.d0/420.d0
          c( 2)=  4.d0/420.d0
        case default
          write(*,*) 'ERROR: Invalid reconstruction order'
          stop
      end select

!     convert fwaves to waves by dividing by the sound speed
      do i=1-mbc+1,mx+mbc
      do mw=1,mwaves
      do m=1,meqn
        wave(i,m,mw)=wave(i,m,mw)/s(i,mw)
      enddo
      enddo
      enddo

      do 20 i=1,mx+1

        do m=1,meqn
          ! Stencil-independent part of reconstruction
          qr(i-1,m) = q(i-1,m)
          ql(i  ,m) = q(i  ,m)
        enddo

        do mw=1,mwaves ! loop over waves
  
          do j=-k,k
            theta(j)=wave(i+j,1,mw)*wave(i,1,mw)
            do m=2,meqn
              theta(j)=theta(j) + wave(i+j,m,mw)*wave(i,m,mw)
            enddo
          enddo
  
          if (theta(0).gt.1.e-14) then
            invwavenorm = 1.d0/theta(0)
            ur=0.d0
            ul=0.d0

            do j=-k,k-1
              ur = ur + c(j)*theta( j)
              ul = ul + c(j)*theta(-j)
            enddo
 
            do m=1,meqn ! Add in contribution from this wave family
              ! reconstruct q^-_\imh
             qr(i-1,m) = qr(i-1,m) + ur*invwavenorm*wave(i,m,mw)
              ! reconstruct q^+_\imh
             ql(i  ,m) = ql(i  ,m) - ul*invwavenorm*wave(i,m,mw)
            enddo
          endif

        enddo !loop over waves

   20     continue  !loop over interfaces


      return
      end

! ===================================================================
      subroutine q2qlqr_poly_wave(q,ql,qr,wave,s,mx)
! ===================================================================

! Reconstruction based on upwind-centered differencing
! With wave-slope reconstruction

      use ClawData
      use ClawParams

      implicit none

      double precision :: q(1-mbc:maxnx+mbc,meqn), &
                ql(1-mbc:maxnx+mbc,meqn), &
                qr(1-mbc:maxnx+mbc,meqn)
      double precision :: wave(1-mbc:maxnx+mbc,meqn,mwaves)
      double precision :: s(   1-mbc:maxnx+mbc,       mwaves)
      double precision :: theta(-3:3)
      double precision :: c(-3:2)
      integer i,j,k,mw,m,mx
      double precision :: ur,ul,wavenormed

      k=(mthlim(1)-1)/2

      select case (mthlim(1))
        case (3)
          c(-1)=1.d0/6.d0
          c( 0)=2.d0/6.d0
        case (7)
          c(-3)=  3.d0/420.d0
          c(-2)=-22.d0/420.d0
          c(-1)= 79.d0/420.d0
          c( 0)=180.d0/420.d0
          c( 1)=-34.d0/420.d0
          c( 2)=  4.d0/420.d0
        case default
          write(*,*) 'ERROR: Invalid reconstruction order'
          stop
      end select

      do 20 i=1-mbc+k,mx+mbc-k+1

        do m=1,meqn
          ! Stencil-independent part of reconstruction
          qr(i-1,m) = q(i-1,m)
          ql(i  ,m) = q(i  ,m)
        enddo

        do mw=1,mwaves ! loop over waves
  
          do j=-k,k
            theta(j)=wave(i+j,1,mw)*wave(i,1,mw)
            do m=2,meqn
              theta(j)=theta(j) + wave(i+j,m,mw)*wave(i,m,mw)
            enddo
          enddo
  
          if (theta(0).gt.1.e-14) then
            wavenormed = 1.d0/theta(0)
            ur=0.d0
            ul=0.d0

            do j=-k,k-1
              ur = ur + c(j)*theta( j)
              ul = ul + c(j)*theta(-j)
            enddo
 
            do m=1,meqn ! Add in contribution from this wave family
              ! reconstruct q^-_\imh
             qr(i-1,m) = qr(i-1,m) + ur*wavenormed*wave(i,m,mw)
              ! reconstruct q^+_\imh
             ql(i  ,m) = ql(i  ,m) - ul*wavenormed*wave(i,m,mw)
            enddo
          endif

        enddo !loop over waves

   20     continue  !loop over interfaces


      return
      end
