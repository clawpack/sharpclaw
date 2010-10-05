      program driver
c
c  Generic driver routine for wclaw2
c
c  WENOCLAW Version 0.2
c
c
      implicit double precision (a-h,o-z)

c     # set parameters for maximum array sizes used in declarations
c     # these must be increased for larger problems.
c
c
      parameter (maxmx =   100)
      parameter (maxmy =   100)
      parameter (mwork =  47700)
      parameter (mrkwork =  11236)

      parameter (mbc = 3)
      parameter (meqn = 1)
      parameter (mwaves = 1)
      parameter (maux = 2)

      dimension    q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      dimension  aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)
      dimension mthlim(mwaves)
      dimension work(mwork),rkwork(mrkwork)
c
      call wclaw2ez(maxmx,maxmy,meqn,mwaves,mbc,maux,mwork,mthlim,
     &           q,work,mrkwork,rkwork,aux)

      stop 
      end
