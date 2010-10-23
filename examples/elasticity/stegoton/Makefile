
# Makefile for Clawpack code in this directory.
# This version only sets the local files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = $(CLAW)/util/Makefile.common

# See the above file for details and a list of make options, or type
#   make .help
# at the unix prompt.

# Adjust these variables if desired:
# ----------------------------------

CLAW_PKG = SharpClaw                  # Clawpack package to use
CLAW_EXE = xsclaw                    # Executable to create
CLAW_setrun_file = setrun.py        # File containing function to make data
CLAW_OUTDIR = _output               # Directory for output
CLAW_setplot_file = setplot.py      # File containing function to set plots
CLAW_PLOTDIR = _plots               # Directory for plots

# Environment variable FC should be set to fortran compiler, e.g. gfortran
FC ?= gfortran   # default if not set as environment variable
# Add any desired compiler flags such as -g here:
FFLAGS = -O3


# ---------------------------------
# List of sources for this program:
# ---------------------------------

    
CLAW_SOURCES = \
  qinit.f \
  rp1nel.f \
  tfluct1.f \
  setaux.f \
  bc1.f \
  out1nel.f \
  b4step1.f \
  setprob.f 


# Clawpack library to be used:
CLAW_LIB = $(SCLAW)/lib
    
CLAW_LIBSOURCES = \
  $(CLAW_LIB)/1d/Global1.f90 \
  $(CLAW_LIB)/ClawData.f90 \
  $(CLAW_LIB)/ClawParams.f90 \
  $(CLAW_LIB)/reconstruct.f90 \
  $(CLAW_LIB)/main.f90 \
  $(CLAW_LIB)/sharpclaw.f90 \
  $(CLAW_LIB)/1d/qalloc1.f90 \
  $(CLAW_LIB)/step.f90 \
  $(CLAW_LIB)/step1.f90 \
  $(CLAW_LIB)/step2.f90 \
  $(CLAW_LIB)/polyrecon.f90 \
  $(CLAW_LIB)/1d/src1.f90 \
  $(CLAW_LIB)/evec.f90 \
  $(CLAW_LIB)/opendatafile.f

#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)


### DO NOT remove this line - make depends on it ###
    
