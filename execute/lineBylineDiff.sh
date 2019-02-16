#!/bin/bash

numRows=$(wc -l $1 | cut -d' ' -f1)

field=$2
verbose=$4

cat $1 $3 | awk -v f=$field -v rows=$numRows -v verbose=$verbose '
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
	if (verbose && val>0) print $1,$2,$3, $f, buf[NR-rows]
	if (max < val) {
		max = val
		maxline = NR-rows
	}
}
END { print "col=",f, "maxErr=", max, "@line:", maxline }
'
