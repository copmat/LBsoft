#!/bin/bash

args=$*

for it in $args
do
	for i in $(seq 1 8)
	do
		cat forceInt.atom00000"$i".step000001.iter"$it".*.dat | ./trisort.sh 3 4 5 > atom"$i".step1.iter"$it".mpi
		cat forceInt.atom00000"$i".step000002.iter"$it".*.dat | ./trisort.sh 4 5 6 > atom"$i".step2.iter"$it".mpi
		cat forceInt.atom00000"$i".step000003.iter"$it".*.dat | ./trisort.sh 3 4 5 > atom"$i".step3.iter"$it".mpi

                cat part_delfluid.atom00000"$i".iter"$it".*.dat | ./trisort.sh 3 4 5 > part_delfluid.atom00000"$i".iter"$it".mpi
                cat part_createfluid.atom00000"$i".iter"$it".*.dat | ./trisort.sh 3 4 5 > part_createfluid.atom00000"$i".iter"$it".mpi

	done

	cat mergeforce.iter"$it".*.dat | ./trisort.sh 2 4 1  | uniq > mergeforce.iter"$it".000000.dat.mpi
	cat force_particle_bb.iter"$it".*.dat | ./trisort.sh 2 4 1 > force_particle_bb.iter"$it".000000.dat.mpi
	cat nve_lf.iter"$it".*.dat | ./trisort.sh 2 4 1 > nve_lf.iter"$it".000000.dat.mpi
done
