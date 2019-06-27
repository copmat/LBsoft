#!/bin/bash
#SBATCH --nodes=2 --tasks-per-node=32 --cpus-per-task=1
#SBATCH --time 1:00:00
#SBATCH -A Pra17_4374
#SBATCH -p knl_usr_dbg
#SBATCH --job-name=lb256
#SBATCH --mail-type=ALL
#SBATCH --mail-user=fabio.bonaccorso@roma2.infn.it


# Source global definitions
module load env-knl
module purge; module load autoload intelmpi/2018--binary

ulimit -s unlimited

mpiifort -v
echo "------------------------"
mpirun ./main_mpi.x
