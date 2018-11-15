
#include <default_macro.h>
 module integrator_mod
 
!***********************************************************************
!     
!     LBsoft module containing integrators data routines
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification july 2018
!     
!***********************************************************************

 use version_mod,     only : idrank,finalize_world
 use error_mod
 use profiling_mod,   only : start_timing2,end_timing2, &
                       ldiagnostic
 use lbempi_mod,      only : i4back,ownern
 use fluids_mod,      only : initialize_fluid_force,compute_fluid_force_sc, &
                        driver_bc_densities,&
                        driver_collision_fluids,compute_omega, &
                        moments_fluids,driver_copy_densities_wall, &
                        lpair_SC,driver_apply_bounceback_pop, &
                        driver_streaming_fluids,aoptpR,test_fake_pops,&
                        probe_pops_in_node,ex,ey,ez,isfluid, &
                        driver_bc_velocities,lbc_halfway, &
                        driver_apply_bounceback_halfway_pop
 use particles_mod,    only : driver_neighborhood_list,lparticles, &
                        vertest,initialize_particle_force, &
                        driver_inter_force,integrate_particles_lf, &
                        integrate_particles_vv,merge_particle_energies,&
                        initialize_particle_energy,lvv, &
                        apply_particle_bounce_back, &
                        store_old_pos_vel_part, &
                        inter_part_and_grid,force_particle_bounce_back,&
                        merge_particle_force
 use write_output_mod, only : write_vtk_frame,write_xyz
 
 implicit none

 private
 
 integer, public, protected, save :: integrator
 integer, public, protected, save :: nstep=0
 integer, public, protected, save :: nstepmax=0
 logical, public, protected, save :: lintegrator=.false.
 
 real(kind=PRC), public, save :: tstep = ONE
 real(kind=PRC), public, save :: initime = ZERO
 real(kind=PRC), public, save :: endtime = ZERO
 logical, public, save :: lendtime
 
 public :: set_nstep
 public :: update_nstep
 public :: driver_integrator
 public :: set_nstepmax
 
 contains
 
 subroutine set_nstepmax(itemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the nstepmax for the 
!     integration loop
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
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
!     licensed under Open Software License v. 3.0 (OSL-3.0)
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
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification july 2018
!     
!***********************************************************************
 
  implicit none
  
  nstep = nstep + 1
  
  return
  
 end subroutine update_nstep
 
 subroutine driver_integrator(mytime)
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the time integration of system
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification july 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout) :: mytime
  
  real(kind=PRC) :: new_time
  
  real(kind=PRC)::myrho,myu,myv,myw
  logical :: ltest
  integer :: itemp,jtemp,ktemp,l
  integer(kind=IPRC) :: i4
  logical :: newlst
  
  
  
  new_time = real(nstep,kind=PRC)*tstep
  
  if(ldiagnostic)call start_timing2("LB","initialize_force")
  call initialize_fluid_force
  if(ldiagnostic)call end_timing2("LB","initialize_force")
  
  if(ldiagnostic)call start_timing2("LB","moments_fluids")
  call moments_fluids
  if(ldiagnostic)call end_timing2("LB","moments_fluids")
  
  if(lparticles)then
    
    newlst=.false.
    if(ldiagnostic)call start_timing2("MD","vertest")
    call vertest(newlst,tstep)
    if(ldiagnostic)call end_timing2("MD","vertest")
      
    if(ldiagnostic)call start_timing2("MD","driver_nlist")
    call driver_neighborhood_list(newlst,nstep)
    if(ldiagnostic)call end_timing2("MD","driver_nlist")
    
    if(ldiagnostic)call start_timing2("MD","driver_inter_f")
    call initialize_particle_energy
    call initialize_particle_force
    call driver_inter_force
    if(ldiagnostic)call end_timing2("MD","driver_inter_f")
    
    if(ldiagnostic)call start_timing2("MD","inter_part_and_grid")
    call inter_part_and_grid
    if(ldiagnostic)call end_timing2("MD","inter_part_and_grid")
    
  endif
  
  if(ldiagnostic)call start_timing2("LB","driver_bc_densities")
  call driver_bc_densities
  if(ldiagnostic)call end_timing2("LB","driver_bc_densities")
  
  if(ldiagnostic)call start_timing2("LB","driver_bc_velocities")
  call driver_bc_velocities
  if(ldiagnostic)call end_timing2("LB","driver_bc_velocities")
  
  if(lpair_SC)then
    if(ldiagnostic)call start_timing2("LB","driver_densities_wall")
    call driver_copy_densities_wall
    if(ldiagnostic)call end_timing2("LB","driver_densities_wall")
    
    if(ldiagnostic)call start_timing2("LB","compute_force_sc")
    call compute_fluid_force_sc
    if(ldiagnostic)call end_timing2("LB","compute_force_sc")
  endif
  
  if(ldiagnostic)call start_timing2("LB","compute_omega")
  call compute_omega
  if(ldiagnostic)call end_timing2("LB","compute_omega")
  
  if(ldiagnostic)call start_timing2("IO","write_vtk_frame")
  call write_vtk_frame(nstep)
  if(ldiagnostic)call end_timing2("IO","write_vtk_frame")
  
  if(ldiagnostic)call start_timing2("LB","collision_fluids")
  call driver_collision_fluids
  if(ldiagnostic)call end_timing2("LB","collision_fluids")
  
  if(lbc_halfway)then
    if(ldiagnostic)call start_timing2("LB","apply_bback_pop_hf")
    call driver_apply_bounceback_halfway_pop
    if(ldiagnostic)call end_timing2("LB","apply_bback_pop_hf")
  endif
  
  if(lparticles)then
    if(ldiagnostic)call start_timing2("IO","write_xyz")
    call write_xyz(nstep)
    if(ldiagnostic)call end_timing2("IO","write_xyz")
    
    if(ldiagnostic)call start_timing2("LB","apply_part_bback")
    call apply_particle_bounce_back
    if(ldiagnostic)call end_timing2("LB","apply_part_bback")
    
    call merge_particle_force
    
    call force_particle_bounce_back
    call initialize_particle_force
    call store_old_pos_vel_part
    if(lvv)then
      if(ldiagnostic)call start_timing2("MD","integrate_vv")
      call integrate_particles_vv(1,nstep)
      if(ldiagnostic)call end_timing2("MD","integrate_vv")
      
      if(ldiagnostic)call start_timing2("MD","integrate_vv")
      call integrate_particles_vv(2,nstep)
      call merge_particle_energies
      if(ldiagnostic)call end_timing2("MD","integrate_vv")
    else
      if(ldiagnostic)call start_timing2("MD","integrate_lf")
      call integrate_particles_lf(nstep)
      call merge_particle_energies
      if(ldiagnostic)call end_timing2("MD","integrate_lf")
    endif
  endif
  
  if(ldiagnostic)call start_timing2("LB","streaming_fluids")
  call driver_streaming_fluids(lparticles)
  if(ldiagnostic)call end_timing2("LB","streaming_fluids")
  
  if(.not.lbc_halfway)then
    if(ldiagnostic)call start_timing2("LB","apply_bback_pop")
    call driver_apply_bounceback_pop
    if(ldiagnostic)call end_timing2("LB","apply_bback_pop")
  endif
  
  mytime = new_time
  
  return
  
 end subroutine driver_integrator

 end module integrator_mod

