#!/bin/bash


curDir=$(pwd)
exe=$curDir/../../execute/main_mpi.x

ls -l $exe

tests="2D_Poiseuille_gradP_xy 2D_Poiseuille_gradP_xz 2D_Poiseuille_gradP_yz 2D_Poiseuille_xy 2D_Poiseuille_xz 2D_Poiseuille_yz 2D_Shear_xy 2D_Shear_xz 2D_Shear_yz 3D_Spinodal 3D_Particle_pbc 3D_Rotating_Particle 3D_Shear_Particle 3D_Particles_box 3D_Particle_SC" 


declare -A comp
declare -A MD
declare -A ris
declare -A ris2

for i in $tests
do
	echo $i
	cd $curDir; cd ../$i

	$exe |& grep -A 100 nstep

	comp[$i]=$(grep "components" input.dat)
	MD[$i]=$(grep "MD" input.dat)

	ris[$i]=$(maxDiff.sh output000000.map output000000.map.orig 4 5 6 7 | awk '{ok[NR] = $4} END{PROCINFO["sorted_in"] = "@ind_num_asc"; for (i in ok) line = line ok[i] " "; print line}')

	if [ -e input.xyz ]; then
                NF=$(head -n4 input.xyz | awk '{if ($1 =="C") print NF}' | head -n1)
		fields=$(seq 2 $NF)
		ris2[$i]=$(maxDiff.sh restart.xyz restart.xyz.orig $fields | awk '{ok[NR] = $4} END{PROCINFO["sorted_in"] = "@ind_num_asc";  for (i in ok) line = line ok[i] " "; print line}')
        fi
done


echo ""
echo ""

for i in $tests
do
	echo "Test:$i"
	echo "${comp[$i]}"
	if [ "${MD[$i]}" != "" ]; then
		echo "${MD[$i]}"
	fi

        echo "Diffs in output"
	echo "${ris[$i]}"
	if [ "${MD[$i]}" != "" ]; then
		echo "Diffs in restart.xyz"
		echo "${ris2[$i]}"
	fi

	echo ""
done
