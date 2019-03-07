#!/bin/bash

it=$1
echo "Diffing iter:$it"

echo "Diffing step1"
for i in $(seq 1 8); do diff -q atom"$i".step1.iter"$it".sort ../atom"$i".step1.iter"$it".mpi ; done

echo "Diffing step2"
for i in $(seq 1 8); do diff -q atom"$i".step2.iter"$it".sort ../atom"$i".step2.iter"$it".mpi ; done

echo "Diffing step3"
for i in $(seq 1 8); do diff -q atom"$i".step3.iter"$it".sort ../atom"$i".step3.iter"$it".mpi ; done

if [[ "$2" != "-n" ]]; then
echo "Diffing presort"
for i in $(seq 1 8); do diff -q presort.atom00000"$i".iter"$it".000000.dat ../presort.atom00000"$i".iter"$it".000000.dat ; done

echo "Diffing presort_pcf"
for i in $(seq 1 8); do diff -q presort_pcf.atom00000"$i".iter"$it".000000.dat ../presort_pcf.atom00000"$i".iter"$it".000000.dat ; done
fi

echo "Diffing part_createfluid"
for i in $(seq 1 8); do diff -q part_createfluid.atom00000"$i".iter"$it".sort ../part_createfluid.atom00000"$i".iter"$it".mpi; done

echo "Diffing part_delfluid"
for i in $(seq 1 8); do diff -q part_delfluid.atom00000"$i".iter"$it".sort ../part_delfluid.atom00000"$i".iter"$it".mpi; done

echo "Diffing mergeforce"
diff -q mergeforce.iter"$it".000000.dat.sort ../mergeforce.iter"$it".000000.dat.mpi

echo "Diffing force_particle_bb"
diff -q force_particle_bb.iter"$it".000000.dat.sort ../force_particle_bb.iter"$it".000000.dat.mpi

echo "Diffing nve_lf"
diff -q nve_lf.iter"$it".000000.dat.sort ../nve_lf.iter"$it".000000.dat.mpi
