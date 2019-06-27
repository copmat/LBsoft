#!/bin/bash

i0=$(printf "%06d\n" $1)
i1=$(printf "%06d\n" $2)

for num in $(seq -w $i0 $i1)
do
	for file in $(ls -rt fx*_*[^l].iter"$num".000000.dat)
	do
		maxDiff.sh $file master/$file 3 4 5 | awk '{if ($4 > 1.0e-15) {print; err=1}} END {exit (1-err)}' && echo $file
		#tkdiff $file master/$file
	done
done
