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

def:	gfortran-mpidebugger

all:
	@echo "Error - you should specify a target machine!"
	@echo "Possible choices are:              "
	@echo "                                   "
	@echo "gfortran                           "
	@echo "gfortran-mpi                       "
	@echo "intel                              "
	@echo "intel-mpi                          "
	@echo "intel-mpi-skylake                  "
	@echo "intel-mpidebugger		  "
	@echo "                                   "
	@echo "Possible choices for debugging are:"
	@echo "                                   "
	@echo "gfortran-debugger                  "
	@echo "gfortran-mpidebugger               "
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
	FFLAGS="-O2 -cpp -DMPI -I $(CURDIR) -c" \
	LDFLAGS="-O2 -cpp -DMPI -I $(CURDIR) -o" \
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
	FFLAGS="-O2 -cpp -DINTEL -DMPI -I $(CURDIR) -c" \
	LDFLAGS="-O2 -cpp -DINTEL -DMPI -I $(CURDIR) -o" \
	TYPE=mpi \
	EX=$(EXP) BINROOT=$(BINROOT) mpi
	
intel-mpi-fast:
	$(MAKE) CC=mpiicc \
	FC=mpiifort \
	CFLAGS="-O2 -I $(CURDIR) -c" \
	FFLAGS="-fast -cpp -DINTEL -DMPI -I $(CURDIR) -c" \
	LDFLAGS="-cpp -DINTEL -DMPI -I $(CURDIR) -o" \
	TYPE=mpi \
	EX=$(EXP) BINROOT=$(BINROOT) mpi
	
intel-mpi-skylake:
	$(MAKE) CC=mpiicc \
	FC=mpiifort \
	CFLAGS="-O2 -I $(CURDIR) -xSKYLAKE-AVX512 -c" \
	FFLAGS="-O2 -cpp -DINTEL -DMPI -I $(CURDIR) -xSKYLAKE-AVX512 -c" \
	LDFLAGS="-O2 -cpp -DINTEL -DMPI -I $(CURDIR) -xSKYLAKE-AVX512 -o" \
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
	-cpp -DMPI -I $(CURDIR) -c" \
	LDFLAGS="-O2 -g -fcheck=all -fbacktrace \
	-ffpe-trap=invalid,zero,overflow,underflow,denormal \
	-cpp -DMPI -I $(CURDIR) -o" \
	TYPE=mpi \
	EX=$(EXP) BINROOT=$(BINROOT) mpi
	
intel-mpidebugger:
	$(MAKE) CC=mpiicc \
	FC=mpiifort \
	CFLAGS="-O0 -g -I $(CURDIR) -c" \
	FFLAGS="-O0 -g -check all -cpp -DINTEL -DMPI -I $(CURDIR) -c" \
	LDFLAGS="-O0 -g -check all -cpp -DINTEL -DMPI -I $(CURDIR) -o" \
	TYPE=mpi \
	EX=$(EXP) BINROOT=$(BINROOT) mpi
	
seq:IR_Precision.o Lib_Pack_Data.o Lib_Base64.o \
	Lib_VTK_IO.o get_mem.o serial_version_mod.o \
	profiling_mod.o parse_mod.o error_mod.o aop_mod.o utility_mod.o \
	comsetup_mod.o fluids_mod.o particles_mod.o write_output_mod.o \
	integrator_mod.o statistic_mod.o io_mod.o main.o
	$(FC) $(LDFLAGS) $(EX) IR_Precision.o \
	Lib_Pack_Data.o Lib_Base64.o Lib_VTK_IO.o get_mem.o \
	serial_version_mod.o profiling_mod.o parse_mod.o error_mod.o \
	aop_mod.o utility_mod.o comsetup_mod.o fluids_mod.o particles_mod.o \
	write_output_mod.o integrator_mod.o statistic_mod.o io_mod.o main.o
	mv $(EX) $(EXE)
	
mpi:IR_Precision.o Lib_Pack_Data.o Lib_Base64.o \
	Lib_VTK_IO.o get_mem.o parallel_version_mod.o \
	profiling_mod.o parse_mod.o error_mod.o aop_mod.o utility_mod.o \
	comsetup_mod.o mpi_comm.o fluids_mod.o particles_mod.o write_output_mod.o \
	integrator_mod.o statistic_mod.o io_mod.o main.o
	$(FC) $(LDFLAGS) $(EX) IR_Precision.o \
	Lib_Pack_Data.o Lib_Base64.o Lib_VTK_IO.o get_mem.o \
	parallel_version_mod.o profiling_mod.o parse_mod.o error_mod.o \
	aop_mod.o utility_mod.o comsetup_mod.o mpi_comm.o fluids_mod.o particles_mod.o \
	write_output_mod.o integrator_mod.o statistic_mod.o io_mod.o main.o
	mv $(EXP) $(EXEP)
	
IR_Precision.o:IR_Precision.f90
	$(FC) $(FFLAGS) IR_Precision.f90
	
Lib_Pack_Data.o:Lib_Pack_Data.f90
	$(FC) $(FFLAGS) Lib_Pack_Data.f90
	
Lib_Base64.o:Lib_Base64.f90
	$(FC) $(FFLAGS) Lib_Base64.f90
	
Lib_VTK_IO.o:Lib_VTK_IO.f90
	$(FC) $(FFLAGS) Lib_VTK_IO.f90
	
get_mem.o:get_mem.c
	$(CC) $(CFLAGS) get_mem.c
	
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
	
aop_mod.o:aop_mod.f90
	$(FC) $(FFLAGS) aop_mod.f90
	
utility_mod.o:utility_mod.f90
	$(FC) $(FFLAGS) utility_mod.f90
	
comsetup_mod.o:comsetup_mod.f90
	$(FC) $(FFLAGS) comsetup_mod.f90
	
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
	
mpi_comm.o:	mpi_comm.f90 aop_mod.o comsetup_mod.o
	$(FC) $(FFLAGS) mpi_comm.f90
	
qsort_mod.o:	qsort_mod.f90
	$(FC) $(FFLAGS) qsort_mod.f90
	
main.o:main.f90
	$(FC) $(FFLAGS) main.f90
	
clean:
	rm -rf *.mod *.o

