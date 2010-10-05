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

LIB2DOBJECTS = \
  Global2.o\
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
  ../step2md.o \
  ../weno_gauss.o \
  ../weno_lines.o \
  qalloc2.o \
  bc2.o \
  src2.o \
  out2.o 

SOURCES = \
  Global2.f95\
  ../ClawParams.f95\
  ../ClawData.f95\
  ../main.f95\
  ../sharpclaw.f95 \
  ../step.f95 \
  ../step2.f95 \
  ../step2md.f95 \
  ../weno_gauss.f95 \
  ../weno_lines.f95 \
  ../polyrecon.f95\
  ../weno.f95\
  ../tvd.f95\
  ../setprob.f95 \
  ../setaux.f95 \
  ../b4step.f95 \
  ../evec.f95 \
  ../tfluct.f95 \
  ../step2.f95 \
  qalloc2.f95 \
  bc2.f95 \
  src2.f95 \
  out2.f95 

lib2d: $(LIB2DOBJECTS)

backup : $(SOURCES)
	cat $(SOURCES) $(IFS) > clawpack2lib.`date +%y%m%d` 

clean:
	rm -f *.o ../*.o *.mod ../*.mod
### DO NOT remove this line - make depends on it ###
