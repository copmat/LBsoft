#!/bin/bash

it=$1
echo "Diffing iter:$it"

echo "Diffing step1"
for i in $(seq 1 8); do diff -q atom"$i".step1.iter"$it".sort ../atom"$i".step1.iter"$it".mpi ; done

echo "Diffing step2"
for i in $(seq 1 8); do diff -q atom"$i".step2.iter"$it".sort ../atom"$i".step2.iter"$it".mpi ; done

echo "Diffing step3"
for i in $(seq 1 8); do diff -q atom"$i".step3.iter"$it".sort ../atom"$i".step3.iter"$it".mpi ; done

echo "Diffing presort"
for i in $(seq 1 8); do diff -q presort.atom00000"$i".iter"$it".000000.dat ../presort.atom00000"$i".iter"$it".000000.dat ; done

echo "Diffing mergeforce"
diff -q mergeforce.iter"$it".000000.dat.sort ../mergeforce.iter"$it".000000.dat.sort

echo "Diffing force_particle_bb"
diff -q force_particle_bb.iter"$it".000000.dat.sort ../force_particle_bb.iter"$it".000000.dat.sort

echo "Diffing nve_lf"
diff -q nve_lf.iter"$it".000000.dat.sort ../nve_lf.iter"$it".000000.dat.sort
