
# Makefile for Clawpack code in this directory.
# This version only sets the local files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = $(CLAW)/util/Makefile.common

# See the above file for details and a list of make options, or type
#   make .help
# at the unix prompt.

# Adjust these variables if desired:
# ----------------------------------

CLAW_PKG = SharpClaw                # Clawpack package to use
CLAW_EXE = xsclaw                   # Executable to create
CLAW_setrun_file = setrun.py        # File containing function to make data
CLAW_OUTDIR = _output               # Directory for output
CLAW_setplot_file = setplot.py      # File containing function to set plots
CLAW_PLOTDIR = _plots               # Directory for plots

# Environment variable FC should be set to fortran compiler, e.g. gfortran
FC ?= gfortran   # default if not set as environment variable
# Add any desired compiler flags such as -g here:
FFLAGS = 


# ---------------------------------
# List of sources for this program:
# ---------------------------------

    
CLAW_SOURCES = \
  qinit.f90 \
  rpf-wave.f \
  setaux.f90 \



# Clawpack library to be used:
CLAW_LIB = $(SCLAW)/lib
    
CLAW_LIBSOURCES = \
  $(CLAW_LIB)/1d/Global1.f90 \
  $(CLAW_LIB)/ClawData.f90 \
  $(CLAW_LIB)/ClawParams.f90 \
  $(CLAW_LIB)/reconstruct.f90 \
  $(CLAW_LIB)/reconstructmd.f90 \
  $(CLAW_LIB)/main.f90 \
  $(CLAW_LIB)/sharpclaw.f90 \
  $(CLAW_LIB)/1d/qalloc1.f90 \
  $(CLAW_LIB)/1d/step1d.f90 \
  $(CLAW_LIB)/flux1.f90 \
  $(CLAW_LIB)/polyrecon.f90 \
  $(CLAW_LIB)/1d/bc1.f90 \
  $(CLAW_LIB)/b4step.f90 \
  $(CLAW_LIB)/1d/src1.f90 \
  $(CLAW_LIB)/1d/out1.f90 \
  $(CLAW_LIB)/opendatafile.f \
  $(CLAW_LIB)/tfluct.f \
  $(CLAW_LIB)/evec.f90 \
  $(CLAW_LIB)/setprob.f90


#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)


### DO NOT remove this line - make depends on it ###
    
