#
#		Makefile for the wenoclaw code:
#
#	To compile the wclaw1 library files type:  make lib
#       (from the lib sub-directory)
#
#	To make an executable, type:              make xwclaw
#       (from the application sub-directory)
#
#       To make an executable that generates
#       output in HDF (version 4) format, type:   make xwclawhdf
#       (from the application sub-directory)
#
#	To compile a single file.f90 type:	  make file.o
#
#	To combine all source files type:	  make program
#       (from the application sub-directory)
#       (sometimes easier to debug all in one)
#
#	To combine all library source:    make backup
#       (date is appended)
#
.SUFFIXES:
.SUFFIXES: .f .f90 .o .mod .F .f90
GFORT = gfortran
FFLAGS = -O 
LFLAGS = 
F77 = gfortran -c
F95 = $(GFORT) -c
LINK	= $(GFORT)

.f90.o: ; $(F95) $(FFLAGS) $*.f90 -o $*.o 
.f90.mod: ; $(F95) $(FFLAGS) $*.f90 -o $*.o 

LIB2DOBJECTS = \
  Global2.o\
  ../ClawParams.o\
  ../ClawData.o\
  ../reconstruct.o\
  ../reconstructmd.o\
  ../main.o\
  ../sharpclaw.o \
  step2d.o \
  ../flux1.o \
  flux2.o \
  flux2md.o \
  ../polyrecon.o\
  ../setprob.o \
  ../setaux.o \
  ../b4step.o \
  ../evec.o \
  ../tfluct.o \
  qalloc2.o \
  bc2.o \
  src2.o \
  out2.o 

SOURCES = \
  Global2.f90\
  ../ClawParams.f90\
  ../ClawData.f90\
  ../reconstruct.f90\
  ../reconstructmd.f90\
  ../main.f90\
  ../sharpclaw.f90 \
  step2d.f90 \
  ../flux1.f90 \
  flux2.f90 \
  flux2md.f90 \
  ../polyrecon.f90\
  ../setprob.f90 \
  ../setaux.f90 \
  ../b4step.f90 \
  ../evec.f90 \
  ../tfluct.f90 \
  qalloc2.f90 \
  bc2.f90 \
  src2.f90 \
  out2.f90 

lib2d: $(LIB2DOBJECTS)

backup : $(SOURCES)
	cat $(SOURCES) $(IFS) > clawpack2lib.`date +%y%m%d` 

clean:
	rm -f *.o ../*.o *.mod ../*.mod
### DO NOT remove this line - make depends on it ###
