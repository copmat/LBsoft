#!/bin/bash

args=$*

for it in $args
do
	for i in $(seq 1 8)
	do
		cat forceInt.atom00000"$i".step000001.iter"$it".*.dat | ./trisort.sh 3 4 5 > atom"$i".step1.iter"$it".mpi
		cat forceInt.atom00000"$i".step000002.iter"$it".*.dat | ./trisort.sh 4 5 6 > atom"$i".step2.iter"$it".mpi
		cat forceInt.atom00000"$i".step000003.iter"$it".*.dat | ./trisort.sh 3 4 5 > atom"$i".step3.iter"$it".mpi
	done

	cat mergeforce.iter"$it".*.dat | ./trisort.sh 2 4 1  | uniq > mergeforce.iter"$it".000000.dat.sort
	cat md.timeadv1.iter"$it".*.dat | ./trisort.sh 2 4 1 > md.timeadv1.iter"$it".000000.dat.sort
	cat nve_lf.iter"$it".*.dat | ./trisort.sh 2 4 1 > nve_lf.iter"$it".000000.dat.sort
done
