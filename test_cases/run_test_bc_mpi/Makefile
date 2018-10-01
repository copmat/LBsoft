#
# makefile
#

FC=gfortran
FFLAGS= -O1 -cpp -c
LDFLAGS= -O1 -cpp -o
EX=compare.x
SHELL=/bin/sh

all: main

main:compare.o
	$(FC) $(LDFLAGS) $(EX) compare.o
	
compare.o:compare.f90
	$(FC) $(FFLAGS) compare.f90
	
clean:
	rm -rf *.mod *.o $(EX)

