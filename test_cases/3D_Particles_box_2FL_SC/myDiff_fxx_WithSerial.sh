#!/bin/bash

i0=$(printf "%06d\n" $1)
i1=$(printf "%06d\n" $2)

for num in $(seq -w $i0 $i1)
do
	for name in fxx_particle_delete_fluids fxx_particle_create_fluids fxx_compute_inter_force fxx_compute_psi_sc_particles fxx_mergeforce fxb_apply_particle_bounce_back
	do
		sort -n -k 2 $name.iter"$num".*.dat > $name.iter"$num".000000.dat.mpi
		file="$name.iter"$num".000000.dat"
		maxDiff.sh $file".mpi" ser/$file 4 5 6 | awk '{if ($4 > 5.0e-15) {print; err=1}} END {exit (1-err)}' && echo $file
		#tkdiff $file master/$file
	done
	echo "-------------------------------"
done
