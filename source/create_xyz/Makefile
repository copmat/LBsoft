#
# makefile
#

FC=gfortran
FFLAGS= -O2 -cpp -c
LDFLAGS= -O2 -cpp -o
EX=create_xyz.x
SHELL=/bin/sh

all: main

main:create_xyz.o
	$(FC) $(LDFLAGS) $(EX) create_xyz.o
	
create_xyz.o:create_xyz.f90
	$(FC) $(FFLAGS) create_xyz.f90
	
clean:
	rm -rf *.mod *.o $(EX)

