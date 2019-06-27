#!/bin/bash

x=$1
y=$2
z=$3

awk -v x=$x -v y=$y -v z=$z '
function myabs(x) {
	if (x>=0) return x
	return -x
}

BEGIN { max = 0 }
{
	if (myabs($1-x)>1) next
	if (myabs($2-y)>1) next
	if (myabs($3-z)>1) next

	print
}
'
