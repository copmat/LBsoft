#!/bin/bash

numRows=$(wc -l $1 | cut -d' ' -f1)

field=$2

cat $1 $3 | awk -v f=$field -v rows=$numRows '
function myabs(x) {
	if (x>=0) return x
	return -x
}

BEGIN { max = -1.0E+300 }
{
	if (NR<=rows) {
		buf[NR] = $f
		next
	}

	val = myabs($f - buf[NR-rows])
	if (max < val) {
		max = val
		if (verbose) print NR, $f, buf[NR-rows], "max =", max
	}
}
END { print max }
'
