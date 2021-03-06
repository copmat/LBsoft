#!/bin/bash

file1=$1
file2=$2

shift; shift

if [ "$1" == "-v" ]; then
   verb=1
   shift
fi

args=$*

curDir=$(dirname $0)

for i in $args
do
  $curDir/lineBylineDiff.sh $file1 $i $file2 $verb
done
