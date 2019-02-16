#!/bin/sh

tool=tkdiff
tool=./myless.sh
tool="diff -q"

tmpfile1=/dev/shm/tmp.file1
tmpfile2=/dev/shm/tmp.file2
tmpfile3=/dev/shm/tmp.file3
tmpfile4=/dev/shm/tmp.file4

name=$1
num=$2

shift; shift

xm=$1
xp=$2
ym=$3
yp=$4
zm=$5
zp=$6

cond="$xm<=\$1 && \$1<=$xp && $ym<=\$2 && \$2<=$yp && $zm<=\$3 && \$3<=$zp"

echo "Diffing $name, ser vs mpi.$num"

echo "Diffing pops"
awk "$cond" "$name".000000.dat > $tmpfile1
awk "$cond" ../"$name".00000"$num".dat > $tmpfile2
$tool $tmpfile1 $tmpfile2

echo "Diffing rho,vel,isfluid"
awk "$cond" "$name".000000.dat1 > $tmpfile3
awk "$cond" ../"$name".00000"$num".dat1 > $tmpfile4
$tool $tmpfile3 $tmpfile4

#rm -f $tmpfile1 $tmpfile2
