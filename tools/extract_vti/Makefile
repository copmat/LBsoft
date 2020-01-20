#
# makefile
#

FC=gfortran
#FFLAGS= -O0 -g -fcheck=all -fbacktrace -c
#LDFLAGS= -O0 -g -fcheck=all -fbacktrace -o
FFLAGS= -O1 -c
LDFLAGS= -O1 -o
EX=main.x
SHELL=/bin/sh

#qshep2d.o:qshep2d.f90
#	$(FC) $(FFLAGS) qshep2d.f90 $(GSLFLAGS)

all: main

main:parse_mod.o extract_vti.o
	$(FC) $(LDFLAGS) $(EX) parse_mod.o extract_vti.o

parse_mod.o:parse_mod.f90
	$(FC) $(FFLAGS) parse_mod.f90

extract_vti.o:extract_vti.f90
	$(FC) $(FFLAGS) extract_vti.f90

clean:
	rm -rf *.mod *.o main.x

