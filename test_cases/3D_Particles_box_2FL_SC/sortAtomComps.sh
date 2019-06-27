#!/bin/bash

args=$*

natoms=8

for numit in $args
do
	it=$(printf "%06d\n" $numit)

	for i in $(seq 1 $natoms)
	do
		cat forceInt.atom00000"$i".step000001.iter"$it".*.dat | ./trisort.sh 3 4 5 > atom"$i".step1.iter"$it".mpi
		cat forceInt.atom00000"$i".step000002.iter"$it".*.dat | ./trisort.sh 3 4 5 > atom"$i".step2.iter"$it".mpi
		cat forceInt.atom00000"$i".step000003.iter"$it".*.dat | ./trisort.sh 3 4 5 > atom"$i".step3.iter"$it".mpi

                cat forceInt_SC.atom00000"$i".iter"$it".*.dat | ./quadsort.sh 2 3 4 5 > atom"$i".SC.iter"$it".mpi
                cat forceInt_SC_grad.atom00000"$i".iter"$it".*.dat | ./quadsort.sh 1 2 3 4 > atom"$i".SC_grad.iter"$it".mpi

                cat part_delfluid.atom00000"$i".iter"$it".*.dat | ./trisort.sh 2 3 4 > part_delfluid.atom00000"$i".iter"$it".mpi
                cat part_createfluid.atom00000"$i".iter"$it".*.dat | ./trisort.sh 2 3 4 > part_createfluid.atom00000"$i".iter"$it".mpi
	done

	cat mergeforce.iter"$it".*.dat | ./trisort.sh 2 4 1  | uniq > mergeforce.iter"$it".000000.dat.mpi
	cat force_particle_bb.iter"$it".*.dat | ./trisort.sh 2 4 1 > force_particle_bb.iter"$it".000000.dat.mpi
	cat nve_lf.iter"$it".*.dat | ./trisort.sh 2 4 1 > nve_lf.iter"$it".000000.dat.mpi
done
