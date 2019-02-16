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
for i in $(seq 1 8); do diff -q presort.atom00000"$i".iter00000"$it".rank000000.dat ../presort.atom00000"$i".iter00000"$it".rank000000.dat ; done

echo "Diffing mergeforce"
diff -q mergeforce.iter00000"$it".000000.dat ../mergeforce.iter00000"$it".000000.dat

echo "Diffing nve_lf"
diff -q nve_lf.iter00000"$it".000000.dat.sort ../nve_lf.iter00000"$it".000000.dat.sort
