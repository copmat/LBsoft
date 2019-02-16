#!/bin/bash


./compareSer_mpi_pops_bounds.sh $1 0 1  8  1  32  1 32
./compareSer_mpi_pops_bounds.sh $1 1 9 16  1  32  1 32
./compareSer_mpi_pops_bounds.sh $1 2 17 24 17 32 1 32
./compareSer_mpi_pops_bounds.sh $1 3 25 32 17 32 1 32
