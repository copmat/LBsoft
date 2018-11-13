#!/bin/sh

curdir=$(pwd)

for i in 2D_* 3D_Spinodal
do
	cd $curdir

	echo $i
	cd $i
	ls

	../../../execute/main.x 
	mv global.map global.map.orig
	mv output000000.map output000000.map.orig 

done
