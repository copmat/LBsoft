#!/bin/bash

for it in $(seq 1 $1)
do
	for i in $(seq 1 8)
	do
		cat forceInt.atom00000"$i".step000001.iter00000"$it".*.dat | ./trisort.sh 3 4 5 > atom"$i".step1.iter"$it".mpi
		cat forceInt.atom00000"$i".step000002.iter00000"$it".*.dat | ./trisort.sh 4 5 6 > atom"$i".step2.iter"$it".mpi
		cat forceInt.atom00000"$i".step000003.iter00000"$it".*.dat | ./trisort.sh 3 4 5 > atom"$i".step3.iter"$it".mpi
	done

	cat nve_lf.iter00000"$it".*.dat | ./trisort.sh 2 4 1 > nve_lf.iter00000"$it".000000.dat.sort
done
