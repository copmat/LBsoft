./diffPops.sh compute_omega $1
#./diffPops.sh aft_inter_part_and_grid $1
#./diffPops.sh aft_compute_densities_wall $1
./diffPops.sh compute_psi_sc $1
#./diffPops.sh aft_compute_sc_particles $1
./diffPops.sh compute_force_sc $1

./diffPops.sh collision_fluids $1
./diffPops.sh aft_apply_bounceback_pop $1
