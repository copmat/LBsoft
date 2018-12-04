#!/bin/sh

curdir=$(pwd)

for i in 2D_* 3D_*
do
	cd $curdir

	echo $i
	cd $i
	ls

	../../../execute/main_mpi.x 
	mv global.map global.map.orig
	mv output000000.map output000000.map.orig 
    if [ -f restart.xyz ]; then
       mv restart.xyz restart.xyz.orig
    fi
done
