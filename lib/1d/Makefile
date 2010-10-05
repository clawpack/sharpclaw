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
#	To compile a single file.f95 type:	  make file.o
#
#	To combine all source files type:	  make program
#       (from the application sub-directory)
#       (sometimes easier to debug all in one)
#
#	To combine all library source:    make backup
#       (date is appended)
#
.SUFFIXES:
.SUFFIXES: .f .f95 .o .mod .F .f95
GFORT = gfortran
FFLAGS = -O 
LFLAGS = 
F77 = gfortran -c
F95 = $(GFORT) -c
LINK	= $(GFORT)

.f95.o: ; $(F95) $(FFLAGS) $*.f95 -o $*.o 
.f95.mod: ; $(F95) $(FFLAGS) $*.f95 -o $*.o 

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
  Global1.f95\
  ../ClawParams.f95\
  ../ClawData.f95\
  ../main.f95\
  ../sharpclaw.f95 \
  ../step.f95 \
  ../step1.f95 \
  ../polyrecon.f95\
  ../weno.f95\
  ../tvd.f95\
  ../setprob.f95 \
  ../setaux.f95 \
  ../b4step.f95 \
  ../evec.f95 \
  ../tfluct.f95 \
  ../step2.f95 \
  qalloc1.f95 \
  bc1.f95 \
  src1.f95 \
  out1.f95 \
  out1_hdf.f95 

lib1d: $(LIB1DOBJECTS)

backup : $(SOURCES)
	cat $(SOURCES) $(IFS) > clawpack1lib.`date +%y%m%d` 

clean:
	rm -f *.o ../*.o *.mod ../*.mod
### DO NOT remove this line - make depends on it ###
