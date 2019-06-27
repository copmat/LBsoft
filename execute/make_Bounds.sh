#!/bin/bash

for i in $(seq -w 000000 $1)
do
	bound=$(head -n2 output$i.map | tail -n1 | awk '{print $1,$2,$3,$4,$5,$6}')
	y=$(( $i + 0 ))
	echo './compareSer_mpi_pops_bounds.sh $1 ' "$y $bound"
done
