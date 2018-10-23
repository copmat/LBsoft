#
# makefile
#

FC=gfortran
FFLAGS= -O1 -cpp -c
LDFLAGS= -O1 -cpp -o
EX=create_xyz.x
SHELL=/bin/sh

all: main

main:create_xyz.o
	$(FC) $(LDFLAGS) $(EX) create_xyz.o
	
create_xyz.o:create_xyz.f90
	$(FC) $(FFLAGS) create_xyz.f90
	
clean:
	rm -rf *.mod *.o $(EX)

