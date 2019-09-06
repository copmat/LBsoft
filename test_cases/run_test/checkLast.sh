#!/bin/bash

curDir=$(pwd)

tests="2D_Poiseuille_gradP_xy 2D_Poiseuille_gradP_xz 2D_Poiseuille_gradP_yz 2D_Poiseuille_xy 2D_Poiseuille_xz 2D_Poiseuille_yz 2D_Shear_xy 2D_Shear_xz 2D_Shear_yz 3D_Spinodal 3D_Particle_pbc 3D_Rotating_Particle 3D_Shear_Particle 3D_Particles_box 3D_Particle_SC" 


for i in $tests
do
	cd $curDir; cd ../$i
	echo $i

	if [ -e output000000.map ]; then
		grep "components" input.dat
		grep "MD" input.dat

		echo "Diffs in output"
		$curDir/maxDiff.sh output000000.map output000000.map.orig 4 5 6 7 | awk 'BEGIN {ORS=" "} {print $4}END {ORS="";print "\n"}'

		if [ -e input.xyz ]; then
			NF=$(head -n4 input.xyz | awk '{if ($1 =="C") print NF}' | head -n1)
			echo "Diffs in restart.xyz"
			$curDir/maxDiff.sh restart.xyz restart.xyz.orig $(seq 2 $NF) | awk 'BEGIN {ORS=" "} {print $4}END {ORS="";print "\n"}'
		fi

	else
		echo "Dead stuff in $i !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	fi

	echo ""
done
