#!/bin/bash

numRows=$(wc -l $1 | cut -d' ' -f1)

field=$2
verbose=$4

cat $1 $3 | awk -v f=$field -v rows=$numRows -v verbose=$verbose '
function myabs(x) {
	if (x>=0) return x
	return -x
}

BEGIN { max = 0 }
{
	if (NR<=rows) {
		buf[NR] = $f
		next
	}
	
	// Skip bc-points
	if ($1<1 || $1>32) next
	if ($2<1 || $2>32) next
	if ($3<1 || $3>32) next

	val = myabs($f - buf[NR-rows])
	if (verbose && val>1.0E-4) print $0, "<->", buf[NR-rows]
	if (max < val) {
		max = val
		maxline = NR-rows
		maxVals = $1 " " $2 " " $3 " "  $f "<-> " buf[NR-rows]
	}
}
END {
	if (max < 1.0E-14) maxVals=""
	print "col=",f, "maxErr=", max, "@line:", maxline, maxVals }
'
