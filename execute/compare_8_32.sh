#!/bin/bash


./compareSer_mpi_pops_bounds.sh $1 0 1  16 1  16 1  16
./compareSer_mpi_pops_bounds.sh $1 1 17 32 1  16 1  16
./compareSer_mpi_pops_bounds.sh $1 2 1  16 17 32 1  16
./compareSer_mpi_pops_bounds.sh $1 3 17 32 17 32 1  16

./compareSer_mpi_pops_bounds.sh $1 4 1  16 1  16 17 32
./compareSer_mpi_pops_bounds.sh $1 5 17 32 1  16 17 32
./compareSer_mpi_pops_bounds.sh $1 6 1  16 17 32 17 32
./compareSer_mpi_pops_bounds.sh $1 7 17 32 17 32 17 32
