#
# makefile
#

FC=gfortran
FFLAGS= -O0 -cpp -c
LDFLAGS= -O0 -cpp -o
EX=test.x
SHELL=/bin/sh

all: main

main:run_test.o
	$(FC) $(LDFLAGS) $(EX) run_test.o
	
run_test.o:run_test.f90
	$(FC) $(FFLAGS) run_test.f90
	
clean:
	rm -rf *.mod *.o $(EX)

