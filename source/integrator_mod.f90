
#include "default_macro.h"
 module integrator_mod
 
!***********************************************************************
!     
!     LBsoft module containing integrators data routines
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification july 2018
!     
!***********************************************************************
 
 use fluids_lattices_mod
 use version_mod,     only : idrank,finalize_world,i4back,ownern
 use error_mod
 use profiling_mod,   only : start_timing2,end_timing2, &
                       ldiagnostic
 use fluids_bc_mod,   only: ixpbc,iypbc,izpbc,isfluid, &
                        compute_densities_wall,driver_bc_isfluid, &
                        driver_bc_new_isfluid,update_isfluid,pimage
 use fluids_mod,      only : initialize_fluid_force,compute_fluid_force_sc, &
                        driver_bc_densities,compute_psi_sc,&
                        driver_collision_fluids,compute_omega, &
                        moments_fluids,rhoR,rhoB, &
                        lpair_SC,lsingle_fluid, &
                        driver_streaming_fluids,test_fake_pops,&
                        probe_pops_in_node, &
                        driver_bc_velocities,lbc_halfway, &
                        driver_apply_bounceback_halfway_pop, &
                        probe_red_moments_in_node,print_all_pops, &
                        print_all_pops_center,driver_bc_pop_selfcomm, &
                        print_all_pops_area_shpere, &
                        nx,ny,nz, &
                        driver_bc_all_pops, print_all_pops2, &
                        rescale_fluid_mass,lmass_rescale,lColourG

 use particles_mod,    only : parlst,lparticles, &
                        vertest,initialize_particle_force, &
                        driver_inter_force,nve_lf, &
                        integrate_particles_vv,merge_particle_energies,&
                        initialize_particle_energy,lvv, &
                        apply_particle_bounce_back, &
                        store_old_pos_vel_part,xxx,yyy,zzz, &
                        inter_part_and_grid,force_particle_bounce_back,&
                        merge_particle_force,build_new_isfluid, &
                        compute_mean_particle_force,rdimx,rdimy,rdimz, &
                        nsphere,compute_psi_sc_particles, &
                        spherelist,restore_particles

 use write_output_mod, only : write_vtk_frame,write_xyz
 
 implicit none

 private
 
 integer, public, protected, save :: integrator
 integer, public, protected, save :: nstep=0
 integer, public, protected, save :: nstepmax=0
 logical, public, protected, save :: lintegrator=.false.
 
 ! is a restore run?
 logical, public, protected, save :: l_restore=.false.
 
 real(kind=PRC), public, save :: tstep = ONE
 real(kind=PRC), public, save :: initime = ZERO
 real(kind=PRC), public, save :: endtime = ZERO
 logical, public, save :: lendtime
 
 public :: set_nstep
 public :: update_nstep
 public :: driver_integrator
 public :: set_nstepmax
 public :: set_restore
 public :: get_restore
 
 contains
 
 subroutine set_nstepmax(itemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the nstepmax for the 
!     integration loop
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification july 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: itemp
  
  nstepmax = itemp
  
  return
  
 end subroutine set_nstepmax
 
 subroutine set_nstep(itemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the nstep counter
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification july 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: itemp
  
  nstep = itemp
  
  return
  
 end subroutine set_nstep
 
 subroutine update_nstep
 
!***********************************************************************
!     
!     LBsoft subroutine for updating the nstep counter
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification july 2018
!     
!***********************************************************************
 
  implicit none
  
  nstep = nstep + 1
  
  return
  
 end subroutine update_nstep
 
 subroutine set_restore(l_temp)
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the l_restore flag for the
!     restarting procedure
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification october 2019
!     
!***********************************************************************
  
  implicit none
  logical, intent(in) :: l_temp

  l_restore = l_temp
  
  return
  
 end subroutine set_restore

 subroutine get_restore(l_temp)
 
!***********************************************************************
!     
!     LBsoft subroutine for inquiring the l_restore flag for the
!     restarting procedure
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification october 2019
!     
!***********************************************************************
  
  implicit none
  logical, intent(out)  :: l_temp
  
  l_temp = l_restore
  
  return
  
 end subroutine get_restore
 
 subroutine driver_integrator(mytime)
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the time integration of system
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification july 2018
!     
!***********************************************************************
 
  implicit none
  real(kind=PRC), intent(inout) :: mytime
  real(kind=PRC) :: new_time
  real(kind=PRC)::myrho,myu,myv,myw
  logical :: ltest
  integer :: itemp,jtemp,ktemp,i,j,k,l,m,ii,jj,kk,i2,j2,k2,i3,j3,k3
  integer(kind=IPRC) :: i4
  logical :: newlst
  logical :: debug = .false., debug1 = .false.


  !debug  = nstep >= 147
  !debug1  = nstep >= 265
  
  new_time = real(nstep,kind=PRC)*tstep
  
  if(lmass_rescale)then
    if(ldiagnostic)call start_timing2("LB","rescale_fluid_mass")
    call rescale_fluid_mass(nstep)
    if(ldiagnostic)call end_timing2("LB","rescale_fluid_mass")
  endif
  
  if(lparticles)then
    if(ldiagnostic)call start_timing2("LB","build_new_isfluid")
    call initialize_particle_force
    call build_new_isfluid(nstep, debug)
    if(ldiagnostic)call end_timing2("LB","build_new_isfluid")
  endif
  
  if(ldiagnostic)call start_timing2("LB","initialize_force")
  call initialize_fluid_force
  if(ldiagnostic)call end_timing2("LB","initialize_force")
  
  if(ldiagnostic)call start_timing2("LB","moments_fluids")
  call moments_fluids(nstep)
  if(ldiagnostic)call end_timing2("LB","moments_fluids")

  if(ldiagnostic)call start_timing2("LB","driver_bc_densities")
  call driver_bc_densities
  if(ldiagnostic)call end_timing2("LB","driver_bc_densities")
  
  if(ldiagnostic)call start_timing2("LB","driver_bc_velocities")
  call driver_bc_velocities
  if(ldiagnostic)call end_timing2("LB","driver_bc_velocities")
  if (debug1) call print_all_pops2(131, "moments_", nstep)
  
  if(ldiagnostic)call start_timing2("LB","compute_omega")
  call compute_omega
  if(ldiagnostic)call end_timing2("LB","compute_omega")
  
  ! if (debug1) call print_all_pops2(131, "aft_compute_omega", nstep)
  
  if(lparticles)then
    if(ldiagnostic)call start_timing2("MD","inter_part_and_grid")
    call inter_part_and_grid(nstep, lparticles, debug)
    if(ldiagnostic)call end_timing2("MD","inter_part_and_grid")
    if (debug1) call print_all_pops2(131, "rmdel_", nstep)
    
    if(ldiagnostic)call start_timing2("MD","driver_bc_densities")
    call driver_bc_densities
    if(ldiagnostic)call end_timing2("MD","driver_bc_densities")
    
    if(ldiagnostic)call start_timing2("MD","driver_bc_velocities")
    call driver_bc_velocities
    if(ldiagnostic)call end_timing2("MD","driver_bc_velocities")
    
    if(ldiagnostic)call start_timing2("MD","update_isfluid")
    call update_isfluid
    if(ldiagnostic)call end_timing2("MD","update_isfluid")
    
    if(ldiagnostic)call start_timing2("MD","driver_bc_isfluid")
    call driver_bc_isfluid
    if(ldiagnostic)call end_timing2("MD","driver_bc_isfluid")
    ! if (debug1) call print_all_pops2(131, "aft_inter_part_and_grid", nstep)
    
    newlst=.false.
    if(ldiagnostic)call start_timing2("MD","vertest")
    call vertest(nstep, newlst)
    if(ldiagnostic)call end_timing2("MD","vertest")
    
    if(newlst)then
      if(ldiagnostic)call start_timing2("MD","driver_nlist")
      call parlst(nstep,newlst)
      if(ldiagnostic)call end_timing2("MD","driver_nlist")
    endif
    
    if(ldiagnostic)call start_timing2("MD","driver_inter_f")
    call initialize_particle_energy
    call driver_inter_force(nstep, .true.)
    if(ldiagnostic)call end_timing2("MD","driver_inter_f")
  endif

  if(lpair_SC .or. lColourG .or. lparticles)then
    if(ldiagnostic)call start_timing2("LB","compute_densities_wall")
    call compute_densities_wall(nstep,lsingle_fluid,rhoR,rhoB)
    if(ldiagnostic)call end_timing2("LB","compute_densities_wall")
    
    if(ldiagnostic)call start_timing2("LB","driver_bc_densities")
    call driver_bc_densities
    if(ldiagnostic)call end_timing2("LB","driver_bc_densities")
    if (debug1) call print_all_pops2(131, "wall_", nstep)
  endif
    
  if(lpair_SC)then
    if(ldiagnostic)call start_timing2("LB","compute_psi_sc")
    call compute_psi_sc
    if(ldiagnostic)call end_timing2("LB","compute_psi_sc")
    ! if (debug1) call print_all_pops2(131, "aft_compute_psi_sc", nstep)

    if(lparticles)then
      if(ldiagnostic)call start_timing2("MD","compute_sc_particles")
      call compute_psi_sc_particles(nstep, debug)
      if(ldiagnostic)call end_timing2("MD","compute_sc_particles")

      ! if (debug1) call print_all_pops2(131, "aft_compute_sc_particles", nstep)
    endif

    if(ldiagnostic)call start_timing2("LB","compute_force_sc")
    call compute_fluid_force_sc
    if(ldiagnostic)call end_timing2("LB","compute_force_sc")
    
    ! if (debug1) call print_all_pops2(131, "aft_compute_force_sc", nstep)
  endif

  if(ldiagnostic)call start_timing2("IO","write_vtk_frame")
   call write_vtk_frame(nstep)
  if(ldiagnostic)call end_timing2("IO","write_vtk_frame")
  if (debug1) call print_all_pops2(131, "aft_vtk", nstep)

  if(ldiagnostic)call start_timing2("LB","collision_fluids")
  call driver_collision_fluids(nstep)
  if(ldiagnostic)call end_timing2("LB","collision_fluids")
  if (debug1) call print_all_pops2(131, "cg_", nstep)
  
  if(lparticles)then
    if(ldiagnostic)call start_timing2("MD","driver_bc_all_pops")
    call driver_bc_all_pops
    if(ldiagnostic)call end_timing2("MD","driver_bc_all_pops")

    if(ldiagnostic)call start_timing2("IO","write_xyz")
    call write_xyz(nstep)
    if(ldiagnostic)call end_timing2("IO","write_xyz")
    
    if(ldiagnostic)call start_timing2("MD","apply_part_bback")
    call apply_particle_bounce_back(nstep, debug, debug1)
    if(ldiagnostic)call end_timing2("MD","apply_part_bback")
    if (debug1) call print_all_pops2(131, "partBB_", nstep)
    
    if(ldiagnostic)call start_timing2("MD","merge_part_force")
    call merge_particle_force(nstep, debug)
    if(ldiagnostic)call end_timing2("MD","merge_part_force")
    
    if(ldiagnostic)call start_timing2("MD","force_particle_bb")
    call force_particle_bounce_back(nstep, debug)
    if(ldiagnostic)call end_timing2("MD","force_particle_bb")
    
    if(ldiagnostic)call start_timing2("MD","mean_partf")
    call compute_mean_particle_force
    if(ldiagnostic)call end_timing2("MD","mean_partf")
    
    if(ldiagnostic)call start_timing2("MD","store_old_part")
    call store_old_pos_vel_part
    if(ldiagnostic)call end_timing2("MD","store_old_part")

    if(ldiagnostic)call start_timing2("MD","integrate_lf")
    call nve_lf(nstep, debug)
    call merge_particle_energies
    if(ldiagnostic)call end_timing2("MD","integrate_lf")
    
    if(ldiagnostic)call start_timing2("MD","restore_partic")
    call restore_particles(nstep)
    if(ldiagnostic)call end_timing2("MD","restore_partic")
    
    if(ldiagnostic)call start_timing2("MD","driver_bc_all_pops")
    call driver_bc_all_pops
    if(ldiagnostic)call end_timing2("MD","driver_bc_all_pops")
  endif
  
  if(lbc_halfway)then
    if(ldiagnostic)call start_timing2("LB","apply_bback_pop_hf")
    call driver_apply_bounceback_halfway_pop(nstep)
    if(ldiagnostic)call end_timing2("LB","apply_bback_pop_hf")
  endif
  if (debug1) call print_all_pops2(131, "bb_", nstep)

  if(ldiagnostic)call start_timing2("LB","streaming_fluids")
  call driver_streaming_fluids(lparticles)
  if(ldiagnostic)call end_timing2("LB","streaming_fluids")
  
  
  if(debug1)call print_all_pops2(131, "step_", nstep)
  
  mytime = new_time
  
  return
  
 end subroutine driver_integrator

 end module integrator_mod

