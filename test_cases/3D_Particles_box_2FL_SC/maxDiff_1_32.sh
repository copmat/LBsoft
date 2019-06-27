#!/bin/bash

curDir=$(dirname $0)

file1=$1
file2=$2

shift; shift

if [[ $1 == "-v" ]];  then
	verbose=1
	shift
fi

args=$*
for i in $args
do
  $curDir/lineBylineDiff_1_32.sh $file1 $i $file2 $verbose
done
