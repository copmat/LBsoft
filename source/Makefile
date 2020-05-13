#
# makefile
#

BINROOT = ../execute
FC=undefined
CC=undefined
FFLAGS=undefined
LDFLAGS=undefined
TYPE=undefined
EX=main.x
EXP=main_mpi.x
EXE = $(BINROOT)/$(EX)
EXEP = $(BINROOT)/$(EXP)
SHELL=/bin/sh

def:	all

all:
	@echo "Error - you should specify a target machine!"
	@echo "Possible choices are:              "
	@echo "                                   "
	@echo "gfortran                           "
	@echo "gfortran-mpi                       "
	@echo "intel                              "
	@echo "intel-mpi                          "
	@echo "intel-mpi-skylake                  "
	@echo "intel-mpi-knl                      "
	@echo "                                   "
	@echo "Possible choices for debugging are:"
	@echo "                                   "
	@echo "gfortran-debugger                  "
	@echo "gfortran-mpidebugger               "
	@echo "intel-debugger                     "
	@echo "intel-mpidebugger                  "
	@echo "                                   "
	@echo "Please examine Makefile for further details "

help: all
	
gfortran:
	$(MAKE) CC=gcc \
	FC=gfortran \
	CFLAGS="-O2 -cpp -I $(CURDIR) -c" \
	FFLAGS="-O2 -cpp -I $(CURDIR) -c" \
	LDFLAGS="-O2 -cpp -I $(CURDIR) -o" \
	TYPE="seq" \
	EX=$(EX) BINROOT=$(BINROOT) seq
	
gfortran-mpi:
	$(MAKE) CC=mpicc \
	FC=mpif90 \
	CFLAGS="-O2 -cpp -I $(CURDIR) -c" \
	FFLAGS="-O2 -cpp -I $(CURDIR) -c" \
	LDFLAGS="-O2 -cpp -I $(CURDIR) -o" \
	TYPE=mpi \
	EX=$(EXP) BINROOT=$(BINROOT) mpi
	
intel:
	$(MAKE) CC=icc \
	FC=ifort \
	CFLAGS="-O2 -I $(CURDIR) -c" \
	FFLAGS="-O2 -cpp -DINTEL -I $(CURDIR) -c" \
	LDFLAGS="-O2 -cpp -DINTEL -I $(CURDIR) -o" \
	TYPE="seq" \
	EX=$(EX) BINROOT=$(BINROOT) seq
	
intel-mpi:
	$(MAKE) CC=mpiicc \
	FC=mpiifort \
	CFLAGS="-O2 -I $(CURDIR) -c" \
	FFLAGS="-O2 -cpp -DINTEL -I $(CURDIR) -c" \
	LDFLAGS="-O2 -cpp -DINTEL -I $(CURDIR) -o" \
	TYPE=mpi \
	EX=$(EXP) BINROOT=$(BINROOT) mpi
	
intel-mpi-fast:
	$(MAKE) CC=mpiicc \
	FC=mpiifort \
	CFLAGS="-O2 -I $(CURDIR) -c" \
	FFLAGS="-fast -cpp -DINTEL -I $(CURDIR) -c" \
	LDFLAGS="-cpp -DINTEL -I $(CURDIR) -o" \
	TYPE=mpi \
	EX=$(EXP) BINROOT=$(BINROOT) mpi
	
intel-mpi-skylake:
	$(MAKE) CC=mpiicc \
	FC=mpiifort \
	CFLAGS="-O2 -I $(CURDIR) -xSKYLAKE-AVX512 -c" \
	FFLAGS="-O2 -cpp -DINTEL -I $(CURDIR) -xSKYLAKE-AVX512 -c" \
	LDFLAGS="-O2 -cpp -DINTEL -I $(CURDIR) -xSKYLAKE-AVX512 -o" \
	TYPE=mpi \
	EX=$(EXP) BINROOT=$(BINROOT) mpi

intel-mpi-knl:
	@mpiifort -v
	@echo "----------------------------"
	$(MAKE) CC=mpiicc \
       FC=mpiifort \
       CFLAGS="-O2 -I $(CURDIR) -xMIC-AVX512 -c" \
       FFLAGS="-Ofast -cpp -DINTEL -I $(CURDIR) -xMIC-AVX512 -c" \
       LDFLAGS="-Ofast -cpp -DINTEL -I $(CURDIR) -xMIC-AVX512 -o" \
       TYPE=mpi \
       EX=$(EXP) BINROOT=$(BINROOT) mpi


gfortran-debugger:
	$(MAKE) CC=gcc \
	FC=gfortran \
	CFLAGS="-O0 -g -cpp -I $(CURDIR) -c" \
	FFLAGS="-O0 -g -fcheck=all -fbacktrace \
	-ffpe-trap=invalid,zero,overflow,underflow,denormal \
	-cpp -I $(CURDIR) -c" \
	LDFLAGS="-O0 -g -fcheck=all -fbacktrace \
	-ffpe-trap=invalid,zero,overflow,underflow,denormal \
	-cpp -I $(CURDIR) -o" \
	TYPE=seq \
	EX=$(EX) BINROOT=$(BINROOT) seq

gfortran-mpidebugger:
	$(MAKE) CC=mpicc \
	FC=mpif90 \
	CFLAGS="-O2 -g -cpp -I $(CURDIR) -c" \
	FFLAGS="-O2 -g -fcheck=all -fbacktrace \
	-ffpe-trap=invalid,zero,overflow,underflow,denormal \
	-cpp -I $(CURDIR) -c" \
	LDFLAGS="-O2 -g -fcheck=all -fbacktrace \
	-ffpe-trap=invalid,zero,overflow,underflow,denormal \
	-cpp -I $(CURDIR) -o" \
	TYPE=mpi \
	EX=$(EXP) BINROOT=$(BINROOT) mpi
	
intel-debugger:
	$(MAKE) CC=icc \
	FC=ifort \
	CFLAGS="-O0 -g -I $(CURDIR) -c" \
	FFLAGS="-O0 -g -check all -traceback -check bounds -cpp -DINTEL -I $(CURDIR) -c" \
	LDFLAGS="-O0 -g -check all -traceback -check bounds -cpp -DINTEL -I $(CURDIR) -o" \
	TYPE=seq \
	EX=$(EX) BINROOT=$(BINROOT) seq
	
intel-mpidebugger:
	$(MAKE) CC=mpiicc \
	FC=mpiifort \
	CFLAGS="-O0 -g -I $(CURDIR) -c" \
	FFLAGS="-O0 -g -check all -traceback -check bounds -cpp -DINTEL -I $(CURDIR) -c" \
	LDFLAGS="-O0 -g -check all -traceback -check bounds -cpp -DINTEL -I $(CURDIR) -o" \
	TYPE=mpi \
	EX=$(EXP) BINROOT=$(BINROOT) mpi
	
seq:get_mem.o fluids_lattices_mod.o serial_version_mod.o \
	profiling_mod.o parse_mod.o error_mod.o utility_mod.o \
	fluids_bc_mod.o fluids_mod.o particles_mod.o write_output_mod.o \
	integrator_mod.o statistic_mod.o io_mod.o main.o
	$(FC) $(LDFLAGS) $(EX) get_mem.o fluids_lattices_mod.o \
	serial_version_mod.o profiling_mod.o parse_mod.o error_mod.o \
	utility_mod.o fluids_bc_mod.o fluids_mod.o particles_mod.o \
	write_output_mod.o integrator_mod.o statistic_mod.o io_mod.o main.o
	mv $(EX) $(EXE)
	
mpi:get_mem.o fluids_lattices_mod.o parallel_version_mod.o \
	profiling_mod.o parse_mod.o error_mod.o utility_mod.o \
	fluids_bc_mod.o fluids_mod.o particles_mod.o write_output_mod.o \
	integrator_mod.o statistic_mod.o io_mod.o main.o
	$(FC) $(LDFLAGS) $(EX) get_mem.o fluids_lattices_mod.o \
	parallel_version_mod.o profiling_mod.o parse_mod.o error_mod.o \
	utility_mod.o fluids_bc_mod.o fluids_mod.o particles_mod.o \
	write_output_mod.o integrator_mod.o statistic_mod.o io_mod.o main.o
	mv $(EXP) $(EXEP)
	
get_mem.o:get_mem.c
	$(CC) $(CFLAGS) get_mem.c
	
fluids_lattices_mod.o:fluids_lattices_mod.f90
	$(FC) $(FFLAGS) fluids_lattices_mod.f90
	
serial_version_mod.o:serial_version_mod.f90
	$(FC) $(FFLAGS) serial_version_mod.f90
	
parallel_version_mod.o:parallel_version_mod.f90
	$(FC) $(FFLAGS) parallel_version_mod.f90
	
profiling_mod.o:profiling_mod.f90
	$(FC) $(FFLAGS) profiling_mod.f90
	
parse_mod.o:parse_mod.f90
	$(FC) $(FFLAGS) parse_mod.f90
	
error_mod.o:error_mod.f90
	$(FC) $(FFLAGS) error_mod.f90
	
utility_mod.o:utility_mod.f90
	$(FC) $(FFLAGS) utility_mod.f90
	
fluids_bc_mod.o:fluids_bc_mod.f90
	$(FC) $(FFLAGS) fluids_bc_mod.f90
	
fluids_mod.o:fluids_mod.f90
	$(FC) $(FFLAGS) fluids_mod.f90
	
particles_mod.o:particles_mod.f90
	$(FC) $(FFLAGS) particles_mod.f90
	
write_output_mod.o:write_output_mod.f90
	$(FC) $(FFLAGS) write_output_mod.f90
	
integrator_mod.o:integrator_mod.f90
	$(FC) $(FFLAGS) integrator_mod.f90
	
statistic_mod.o:statistic_mod.f90
	$(FC) $(FFLAGS) statistic_mod.f90
	
io_mod.o:io_mod.f90
	$(FC) $(FFLAGS) io_mod.f90
	
qsort_mod.o:	qsort_mod.f90
	$(FC) $(FFLAGS) qsort_mod.f90
	
main.o:main.f90
	$(FC) $(FFLAGS) main.f90
	
clean:
	rm -rf *.mod *.o

