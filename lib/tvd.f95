! ===================================================================
subroutine tvd2(q,g,mx)
! ===================================================================
! Second order TVD reconstruction for WENOCLAW

      use ClawData
      use ClawParams

      implicit double precision (a-h,o-z)

      type(griddat) g
      dimension  q(1-mbc:mx+mbc,meqn)

!     # loop over all equations (all components).  
!     # the reconstruction is performed component-wise

      do m=1,meqn

!        # compute and store the differences of the cell averages

         do i=2-mbc,mx+mbc
            dqm=dqp
            dqp=q(i+1,m)-q(i,m)
            r=dqp/dqm

            select case(mthlim(m))
                case(1)
                    ! minmod
                    qlimitr = dmax1(0.d0, dmin1(1.d0, r))
                case(2)
                    ! superbee
                    qlimitr = dmax1(0.d0, dmin1(1.d0, 2.d0*r), dmin1(2.d0, r))
                case(3)
                    ! van Leer
                    qlimitr = (r + dabs(r)) / (1.d0 + dabs(r))
                case(4)
                    ! monotonized centered
                    c = (1.d0 + r)/2.d0
                    qlimitr = dmax1(0.d0, dmin1(c, 2.d0, 2.d0*r))
                case(5)
                    ! Cada & Torrilhon simple
                    beta=2.d0
                    xgamma=2.d0
                    alpha=1.d0/3.d0
                    pp=(2.d0+r)/3.d0
                    amax = dmax1(-alpha*r,0.d0,dmin1(beta*r,pp,xgamma))
                    qlimitr = dmax1(0.d0, dmin1(pp,amax))
            end select

           g%qr(i,m) = q(i,m) + 0.5d0*qlimitr*dqm
           g%ql(i,m) = q(i,m) - 0.5d0*qlimitr*dqm

         enddo
      enddo

      return
      end

! ===================================================================
subroutine tvd2_char(q,g,mx)
! ===================================================================

! Second order TVD reconstruction for WENOCLAW
! This one uses characteristic decomposition

!  g%evl, g%evr are left and right eigenvectors at each interface

    use ClawData
    use ClawParams

    implicit double precision (a-h,o-z)

    type(griddat) g

    dimension q(1-mbc:maxnx+mbc,meqn)

    dimension dq(1-mbc:maxnx+mbc,meqn)
    dimension  u(1-mbc:maxnx+mbc,meqn,2)
    dimension hh(1-mbc:maxnx+mbc,-1:1)


    ! loop over all equations (all components).  
    ! the reconstruction is performed using characteristic decomposition

    ! compute and store the differences of the cell averages
    do m=1,meqn
        do i=2-mbc,mx+mbc
            dq(i,m)=q(i,m)-q(i-1,m)
        enddo
    enddo

    do m=1,meqn

        ! Project the difference of the cell averages to the
        ! 'm'th characteristic field
        do m1 = -1,1
            do  i = 2-mbc,mx+mbc
                hh(i,m1) = 0.d0
                do mm=1,meqn
                    hh(i,m1) = hh(i,m1)+ g%evl(i,m,mm)*dq(i+m1,mm)
                enddo
            enddo
        enddo


        ! the reconstruction
        do m1=1,2
            im=(-1)**(m1+1)
            ! m1=1: construct qr_i-1
            ! m1=2: construct ql_i

            do i=1,mx+1
                ! dqp=hh(i,m1-1)
                ! dqm=hh(i,m1-2)
                if (dabs(hh(i,m1-2)).gt.1.e-14) then
                    r=hh(i,m1-1)/hh(i,m1-2)
                else
                    r=0.d0
                endif
           
                select case(mthlim(m))
                    case(1)
                        ! minmod
                        slimitr = dmax1(0.d0, dmin1(1.d0, r))
                    case(2)
                        ! superbee
                        slimitr = dmax1(0.d0, dmin1(1.d0, 2.d0*r), dmin1(2.d0, r))
                    case(3)
                        ! van Leer
                        slimitr = (r + dabs(r)) / (1.d0 + dabs(r))
                    case(4)
                        ! monotonized centered
                        c = (1.d0 + r)/2.d0
                        slimitr = dmax1(0.d0, dmin1(c, 2.d0, 2.d0*r))
                end select
    
                 u(i,m,m1) = im*0.5d0*slimitr*hh(i,m1-2)

            enddo
        enddo
    enddo

    ! Project to the physical space:
    do m =  1, meqn
        do i = 2-mbc,mx+mbc
            g%qr(i-1,m)=q(i-1,m)
            g%ql(i  ,m)=q(i  ,m)
            do mm=1,meqn 
                g%qr(i-1,m) = g%qr(i-1,m) + g%evr(i,m,mm)*u(i,mm,1)
                g%ql(i  ,m) = g%ql(i  ,m) + g%evr(i,m,mm)*u(i,mm,2)
            enddo
        enddo
    enddo
end subroutine tvd2_char

! ===================================================================
subroutine tvd2_wave(q,g,mx)
! ===================================================================
! Second order TVD reconstruction for WENOCLAW
! This one uses projected waves

    use ClawData
    use ClawParams

    implicit double precision (a-h,o-z)

    type(griddat) g

    dimension q(1-mbc:maxnx+mbc,meqn)


    do i=1,mx+1
        do m =  1, meqn
            g%qr(i-1,m) = q(i-1,m)
            g%ql(i  ,m) = q(i  ,m)
        enddo 
    enddo

    ! loop over all equations (all components).  
    ! the reconstruction is performed using characteristic decomposition

    do 200 mw=1,mwaves
        if (mthlim(mw).eq.0) go to 200
        dotr = 0.d0
        do 190 i=0,mx+1
          wnorm2=0.d0
          dotl=dotr
          dotr=0.d0
          do 5 m=1,meqn
            wnorm2 = wnorm2 + g%wave(i,m,mw)**2
            dotr = dotr + g%wave(i,m,mw)*g%wave(i+1,m,mw)
   5        continue
          if (i.eq.0) go to 190
          if (wnorm2.eq.0.d0) go to 190
            
          if (g%s(i,mw).gt.0.d0) then
            r = dotl / wnorm2
          else
            r = dotr / wnorm2
          endif

                select case(mthlim(mw))
                    case(1)
                        ! minmod
                        wlimitr = dmax1(0.d0, dmin1(1.d0, r))
                    case(2)
                        ! superbee
                        wlimitr = dmax1(0.d0, dmin1(1.d0, 2.d0*r), dmin1(2.d0, r))
                    case(3)
                        ! van Leer
                        wlimitr = (r + dabs(r)) / (1.d0 + dabs(r))
                    case(4)
                        ! monotonized centered
                        c = (1.d0 + r)/2.d0
                        wlimitr = dmax1(0.d0, dmin1(c, 2.d0, 2.d0*r))
                    case(5)
                        ! Cada & Torrilhon simple
                        beta=2.d0
                        xgamma=2.d0
                        alpha=1.d0/3.d0
                        pp=(2.d0+r)/3.d0
                        amax = dmax1(-alpha*r,0.d0,dmin1(beta*r,pp,xgamma))
                        wlimitr = dmax1(0.d0, dmin1(pp,amax))
                end select


             g%uu(i,mw) = 0.5d0*wlimitr

             do m =  1, meqn
               g%qr(i-1,m) = g%qr(i-1,m) + g%wave(i,m,mw)*g%uu(i,mw)
               g%ql(i  ,m) = g%ql(i  ,m) - g%wave(i,m,mw)*g%uu(i,mw)
             enddo ! end loop over equations

  190      continue !end loop over interfaces
  200    continue !end loop over waves

      return
      end
