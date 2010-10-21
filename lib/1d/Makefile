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

LIB1DOBJECTS = \
  Global1.o\
  ../ClawParams.o\
  ../ClawData.o\
  ../main.o\
  ../sharpclaw.o \
  ../step.o \
  ../step1.o \
  ../polyrecon.o\
  ../weno.o\
  ../tvd.o\
  ../setprob.o \
  ../setaux.o \
  ../b4step.o \
  ../evec.o \
  ../tfluct.o \
  ../step2.o \
  qalloc1.o \
  bc1.o \
  src1.o \
  out1.o \
  out1_hdf.o

SOURCES = \
  Global1.f90\
  ../ClawParams.f90\
  ../ClawData.f90\
  ../main.f90\
  ../sharpclaw.f90 \
  ../step.f90 \
  ../step1.f90 \
  ../polyrecon.f90\
  ../weno.f90\
  ../tvd.f90\
  ../setprob.f90 \
  ../setaux.f90 \
  ../b4step.f90 \
  ../evec.f90 \
  ../tfluct.f90 \
  ../step2.f90 \
  qalloc1.f90 \
  bc1.f90 \
  src1.f90 \
  out1.f90 \
  out1_hdf.f90 

lib1d: $(LIB1DOBJECTS)

backup : $(SOURCES)
	cat $(SOURCES) $(IFS) > clawpack1lib.`date +%y%m%d` 

clean:
	rm -f *.o ../*.o *.mod ../*.mod
### DO NOT remove this line - make depends on it ###
