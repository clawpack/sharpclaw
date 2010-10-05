c
c  =============================================================================
c  SharpClaw -- Semi-discrete High order Accurate Package for Conservation LAWs
c  =============================================================================
c
c  Written by David Ketcheson based on the CLAWPACK software of R.J. LeVeque.
c
c  Description:
c  ------------
c     SharpClaw is a library of Fortran routines designed to solve
c     hyperbolic systems of PDEs with arbitrarily high order accuracy. 
c     The numerical method used to solve these equations is 
c     described in a paper at
c     http://www.amath.washington.edu/~ketch/home/CV_files/ketch_hyp2006.pdf
c
c  Directory structure:
c  --------------------
c     In this directory the following sub-directories are available:
c     $SharpClaw/1d/lib       --  contains the SharpClaw library files
c     $SharpClaw/1d/examples  --  example problems
c     $SharpClaw/1d/utils     --  Perl scripts for modifying CLAWPACK setups
c                                to run in SharpClaw
c     $SharpClaw/2d/lib       --  contains the 2D SharpClaw library files
c     $SharpClaw/2d/examples  --  2D example problems
c
c =========================================================================

c
c  Copyright by the authors
c
c  This software is made available for research and instructional use only.
c  You may copy and use this software without charge for these non-commercial
c  purposes, provided that the copyright notice and associated text is
c  reproduced on all copies.  For all other uses (including distribution of
c  modified versions), please contact the author at the address given below.
c
c  *** This software is made available "as is" without any assurance that it
c  *** will work for your purposes.  The software may in fact have defects, so
c  *** use the software at your own risk.
c
c  --------------------------------------
c    CLAWPACK Version 4.2,  
c    Webpage: http://www.amath.washington.edu/~claw
c  --------------------------------------
c    Author:  Randall J. LeVeque
c             Applied Mathematics
c             Box 352420
c             University of Washington,
c             Seattle, WA 98195-2420
c             rjl@amath.washington.edu
c =========================================================================
