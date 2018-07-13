#
# makefile
#

CC=gcc
FC=gfortran
CFLAGS= -O0 -g -c
FFLAGS= -O0 -g -fbacktrace -ffpe-trap=zero -fcheck=all -cpp -I $(CURDIR) -c
LDFLAGS= -O0 -g -fbacktrace -ffpe-trap=zero -fcheck=all -cpp -I $(CURDIR) -o
FFLAGSVTK=-O0 -g -fcheck=all -cpp -DPRC=4 -DIPRC=4 -c
#CFLAGS= -O2 -c
#FFLAGS= -O2 -cpp -I $(CURDIR) -c
#LDFLAGS= -O2 -cpp -I $(CURDIR) -o
#FFLAGSVTK= -O2 -cpp -DPRC=4 -DIPRC=4 -c
EX=../execute/main.x
SHELL=/bin/sh

all: main

main:IR_Precision.o Lib_Pack_Data.o Lib_Base64.o \
	Lib_VTK_IO.o get_mem.o serial_version_mod.o \
	profiling_mod.o parse_mod.o error_mod.o utility_mod.o \
	fluids_mod.o write_output_mod.o integrator_mod.o io_mod.o main.o
	$(FC) $(LDFLAGS) $(EX) IR_Precision.o \
	Lib_Pack_Data.o Lib_Base64.o Lib_VTK_IO.o get_mem.o \
	serial_version_mod.o profiling_mod.o parse_mod.o error_mod.o \
	utility_mod.o fluids_mod.o write_output_mod.o integrator_mod.o io_mod.o main.o
	
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
	
utility_mod.o:utility_mod.f90
	$(FC) $(FFLAGS) utility_mod.f90
	
fluids_mod.o:fluids_mod.f90
	$(FC) $(FFLAGS) fluids_mod.f90
	
write_output_mod.o:write_output_mod.f90
	$(FC) $(FFLAGS) write_output_mod.f90
	
integrator_mod.o:integrator_mod.f90
	$(FC) $(FFLAGS) integrator_mod.f90
	
io_mod.o:io_mod.f90
	$(FC) $(FFLAGS) io_mod.f90
	
main.o:main.f90
	$(FC) $(FFLAGS) main.f90
	
clean:
	rm -rf *.mod *.o $(EX)

