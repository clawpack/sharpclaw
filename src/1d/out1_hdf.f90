!
!
! =========================================================
      subroutine out1(maxmx,meqn,mbc,mx,xlower,dx,q,t,iframe,aux,maux)
! =========================================================
!
!
!     # Output the results for a general system of conservation laws
!     # in 1 dimension to a hdf file as a scientific data set.
!     # See http://hdf.ncsa.uiuc.edu/ for more info on HDF.
!     # The results in this output file can be plotted in MATLAB
!     # using the "plotclaw1" script.
!
!     # Revised 2003 by Peter Blossey
!     # Adapted code written by Sorin Mitran to standard F77 CLAWPACK style
!
      implicit double precision (a-h,o-z)
!
      parameter   (nDim = 1)
      dimension    q(1-mbc:maxmx+mbc, meqn)
      dimension aux(1-mbc:maxmx+mbc, maux)
      character*14 fname
      character*13 qname2
      character*9  qname
!
!     # HDF: Declare variables that describe datasets and HDF files.
!
      integer    sd_id, sds_id, sds_rank
      integer    sds_dims, sds_start, sds_edges, sds_stride
      dimension  sds_dims(nDim), sds_start(nDim)
      dimension  sds_edges(nDim), sds_stride(nDim) 
      dimension  qbuf(21), qout(mx)
!
!     # HDF: Declare external HDF functions
!
      integer  sfstart, sfcreate, sfwdata, sfscompress, sfendacc, sfend
      external sfstart, sfcreate, sfwdata, sfscompress, sfendacc, sfend
!
!     # HDF: Set up HDF constants
!
      integer 	DFACC_READ, DFACC_WRITE, DFACC_CREATE 
      parameter(DFACC_READ = 1, DFACC_WRITE = 2, DFACC_CREATE = 4)

      integer   DFNT_FLOAT64, DFNT_INT32
      parameter(DFNT_FLOAT64 = 6, DFNT_INT32 = 24)

      integer   SUCCEED, FAIL
      parameter(SUCCEED = 0, FAIL = -1)
!
!     # HDF: Set up compression constants for HDF file.
!
      integer   COMP_CODE_DEFLATE, DEFLATE_LEVEL
      parameter (COMP_CODE_DEFLATE = 4, DEFLATE_LEVEL = 6)
!
!     # Write the results to the file fort.q<iframe>
!     # Use format required by matlab script  plotclaw2.m
!     # The same format is used by the amrclaw package.
!     # Here it's adapted to output just the single grid.
!
!     # first create the file name and open file
!
      fname = 'fort.q' &
           // char(ichar('0') + mod(iframe/1000,10))  &
           // char(ichar('0') + mod(iframe/100,10))  &
           // char(ichar('0') + mod(iframe/10,10))  &
           // char(ichar('0') + mod(iframe,10)) &
           // '.hdf'
!
!     # Specify grid number and create a string which will describe this
!     # grid in HDF file.  This could be an input for simulations with
!     # multiple grids, as in AMRCLAW.  
!
      ngrids_out=1
      qname = 'grid_' &
           // char(ichar('0') + mod(ngrids_out/1000,10))  &
           // char(ichar('0') + mod(ngrids_out/100,10))  &
           // char(ichar('0') + mod(ngrids_out/10,10))  &
           // char(ichar('0') + mod(ngrids_out,10))
!
!     # HDF: create hdf file.
!
      sd_id = sfstart(fname, DFACC_CREATE)
      if (sd_id.eq.FAIL) THEN
         WRITE(*,*) 'Failed to create HDF file (call to sfstart)'
         STOP
      end if
!
!     # HDF: create a data set for parameters describing q in HDF file.
!
      sds_rank = 1
      sds_dims(1) = 21
      
      sds_id = sfcreate(sd_id,qname,DFNT_FLOAT64,sds_rank,sds_dims)
      if (sds_id.eq.FAIL) THEN
         WRITE(*,*) 'Failed to create scientific data set in HDF file'
         STOP
      end if
!
!     # HDF: set up parameters describing data set.
!
      sds_start(1)  = 0
      sds_edges(1)  = sds_dims(1)
      sds_stride(1) = 1        
      qbuf(1) = ngrids_out
      qbuf(2) = nDim
      qbuf(3) = t
      qbuf(4) = meqn
      qbuf(5) = 1.
      qbuf(6) = mx
      qbuf(7) = 0.
      qbuf(8) = 0.
      qbuf(9) = 0.
      qbuf(10) = xlower
      qbuf(11) = 0.
      qbuf(12) = 0.
      qbuf(13) = 0.
      qbuf(14) = xlower+mx*dx
      qbuf(15) = 0.
      qbuf(16) = 0.
      qbuf(17) = 0.
      qbuf(18) = dx
      qbuf(19) = 0.
      qbuf(20) = 0.
      qbuf(21) = 0.
      istat = sfwdata(sds_id,sds_start,sds_stride,sds_edges,qbuf)  
      istat = sfendacc(sds_id)  
!
!     # Loop over fields in q
!
      do m = 1,meqn
!
!     # HDF: create a data set for parameters describing q in HDF file.
!
         qname2 = qname // '_' &
              // char(ichar('0') + mod(m/100,10))  &
              // char(ichar('0') + mod(m/10,10))  &
              // char(ichar('0') + mod(m,10))
!
         sds_rank = nDim
         sds_dims(1) = mx
!      
         sds_id = sfcreate(sd_id,qname,DFNT_FLOAT64,sds_rank,sds_dims)
         if (sds_id.eq.FAIL) THEN
            WRITE(*,*) 'Failed to create data set in HDF file'
            STOP
         end if
!
!     # HDF: set up parameters describing data set.
!
         sds_start(1)  = 0
         sds_edges(1)  = sds_dims(1)
         sds_stride(1) = 1        

!     # Copy current field of q into qout.
!
         do i = 1,mx
            qout(i) = q(i,m)
         end do
!
!     # HDF: set compression mode and write data to hdf file.
!
         istat=sfscompress(sds_id,COMP_CODE_DEFLATE,DEFLATE_LEVEL)
         istat = sfwdata(sds_id,sds_start,sds_stride,sds_edges,qout)
         if (istat.eq.FAIL) then
            WRITE(*,*) 'Failed to write SDS (call to sfwdata)'
            STOP
         end if
!
!     # HDF: Close the data set
!
         istat = sfendacc(sds_id)  
         if (istat.eq.FAIL) then
            WRITE(*,*) 'Failed to close SDS (call to sfendacc)'
            STOP
         end if
      end do
!
!     # HDF: Close HDF file.
!
      istat = sfend(sd_id)
      if (istat.eq.FAIL) then
         WRITE(*,*) 'Failed to close ', fname, ' (call to sfend)'
         STOP
      end if

      return
      end
