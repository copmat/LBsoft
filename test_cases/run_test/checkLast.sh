#!/bin/bash

exe=/home/nodal2/IIT/LBsoft/execute/main_mpi.x

curDir=$(pwd)
declare -A ris

tests="2D_Poiseuille_gradP_xy 2D_Poiseuille_gradP_xz 2D_Poiseuille_gradP_yz 2D_Poiseuille_xy 2D_Poiseuille_xz 2D_Poiseuille_yz 2D_Shear_xy 2D_Shear_xz 2D_Shear_yz 3D_Spinodal 3D_Particle_pbc 3D_Rotating_Particle 3D_Shear_Particle 3D_Particles_box 3D_Particle_SC" 


for i in $tests
do
	if [ -e output000000.map ]; then
		echo $i
		cd $curDir; cd ../$i
		grep "components" input.dat
		grep "MD" input.dat

		echo "Diffs in output"
		maxDiff.sh output000000.map output000000.map.orig 4 5 6 7 | awk '{ok[NR] = $4} END{PROCINFO["sorted_in"] = "@ind_num_asc"; OFMT="%14.3f"; for (i in ok) line = line ok[i] " "; print line}'

		if [ -e input.xyz ]; then
			NF=$(head -n4 input.xyz | awk '{if ($1 =="C") print NF}' | head -n1)
			echo "Diffs in restart.xyz"
			maxDiff.sh restart.xyz restart.xyz.orig $(seq 2 $NF) | awk '{ok[NR] = $4} END{PROCINFO["sorted_in"] = "@ind_num_asc"; OFMT="%14.3f"; for (i in ok) line = line ok[i] " "; print line}'
		fi

		echo ""
	fi
done
