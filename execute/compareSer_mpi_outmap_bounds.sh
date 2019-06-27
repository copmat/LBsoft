#!/bin/sh

tmpfile1=/dev/shm/tmp.file1
tmpfile2=/dev/shm/tmp.file2

name=$1
num=$(printf "%06d\n" $2)

shift; shift

xm=$1
xp=$2
ym=$3
yp=$4
zm=$5
zp=$6

cond="$xm<=\$1 && \$1<=$xp && $ym<=\$2 && \$2<=$yp && $zm<=\$3 && \$3<=$zp"

echo "Diffing output.map, ser vs mpi.$num"

awk "$cond" output000000.map > $tmpfile1
awk "$cond" ../output"$num".map > $tmpfile2

maxDiff.sh $tmpfile1 $tmpfile2 $(seq 4 9)
