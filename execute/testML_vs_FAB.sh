#!/bin/bash

if [ "$#" -lt 1 ]; then
        echo "Usage: $0 Num_Procs"
        exit
fi


exe_ml=./main_mpi.x.ml
exefab=./main_mpi.x.fab

for test in 2D_Poiseuille_xy 2D_Poiseuille_xz 2D_Poiseuille_yz 3D_Particle_pbc 3D_Rotating_Particle 3D_Spinodal run_test_dec_xyz
do
	echo "TEST: $test"

	cp ../test_cases/$test/input.* .

	rm -fr fab ml; mkdir fab ml

	mpirun -np $1 $exe_ml |& sed -ne "/nstep/,/MEMORY MONITOR*/p" ; mv output* ml 
	mpirun -np $1 $exefab |& sed -ne "/nstep/,/MEMORY MONITOR*/p" ; mv output* fab
	diff -rq ml/ fab/
		
	if [ $? != 0 ]; then
		echo "Test: $test NOK"
		exit
	fi

done
