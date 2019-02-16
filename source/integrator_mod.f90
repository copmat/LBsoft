
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
                        driver_bc_densities,compute_psi_sc,&
                        driver_collision_fluids,compute_omega, &
                        moments_fluids,compute_densities_wall, &
                        lpair_SC,driver_apply_bounceback_pop, &
                        driver_streaming_fluids,aoptpR,test_fake_pops,&
                        probe_pops_in_node,ex,ey,ez,isfluid, &
                        driver_bc_velocities,lbc_halfway, &
                        driver_apply_bounceback_halfway_pop, &
                        probe_red_moments_in_node,print_all_pops, &
                        print_all_pops_center,driver_bc_pop_selfcomm, &
                        print_all_pops_area_shpere,ex,ey,ez,pimage, &
                        ixpbc,iypbc,izpbc,nx,ny,nz,opp, &
                        aoptpR, driver_bc_pops, print_all_pops2, &
                        driver_bc_pops_NOK

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
  integer :: itemp,jtemp,ktemp,i,j,k,l,m,ii,jj,kk,i2,j2,k2,i3,j3,k3
  integer(kind=IPRC) :: i4
  logical :: newlst
  logical :: debug
  

  debug = nstep <= 10
  debug = modulo(nstep, 500) == 0


  new_time = real(nstep,kind=PRC)*tstep
  
  if(lparticles)then
    if(ldiagnostic)call start_timing2("LB","build_new_isfluid")
    call initialize_particle_force
    call build_new_isfluid(nstep)
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
  
  if(ldiagnostic)call start_timing2("LB","compute_omega")
  call compute_omega
  if(ldiagnostic)call end_timing2("LB","compute_omega")
  
  if (debug) call print_all_pops2(131, "aft_compute_omega", nstep, aoptpR)

  if(lparticles)then
    if(ldiagnostic)call start_timing2("MD","inter_part_and_grid")
    call inter_part_and_grid(nstep, lparticles)
    if(ldiagnostic)call end_timing2("MD","inter_part_and_grid")
    
    if (debug) call print_all_pops2(131, "aft_inter_part_and_grid", nstep, aoptpR)

    newlst=.false.
    if(ldiagnostic)call start_timing2("MD","vertest")
    call vertest(newlst)
    if(ldiagnostic)call end_timing2("MD","vertest")
    
    if(ldiagnostic)call start_timing2("MD","driver_nlist")
    if(newlst) call parlst
    if(ldiagnostic)call end_timing2("MD","driver_nlist")
    
    if(ldiagnostic)call start_timing2("MD","driver_inter_f")
    call initialize_particle_energy
    call driver_inter_force(nstep)
    if(ldiagnostic)call end_timing2("MD","driver_inter_f")
  endif

  if (debug) call print_all_pops2(131, "aft_inter_force", nstep, aoptpR)

  if(lpair_SC .or. lparticles)then
    ! Solved in inter_part_and_grid) call driver_bc_pops

    if(ldiagnostic)call start_timing2("LB","compute_densities_wall")
    call compute_densities_wall
    if(ldiagnostic)call end_timing2("LB","compute_densities_wall")
  endif
!  if (debug) call print_all_pops2(131, "aft_densities_wall", nstep, aoptpR)
    
  if(lpair_SC)then
    if(ldiagnostic)call start_timing2("LB","compute_psi_sc")
    call compute_psi_sc
    if(ldiagnostic)call end_timing2("LB","compute_psi_sc")
    
    if(lparticles)then
      if(ldiagnostic)call start_timing2("MD","compute_sc_particles")
      call compute_psi_sc_particles(nstep)
      if(ldiagnostic)call end_timing2("MD","compute_sc_particles")
    endif
    
    if(ldiagnostic)call start_timing2("LB","compute_force_sc")
    call compute_fluid_force_sc
    if(ldiagnostic)call end_timing2("LB","compute_force_sc")
  endif
  
  if(ldiagnostic)call start_timing2("LB","collision_fluids")
  call driver_collision_fluids
  if(ldiagnostic)call end_timing2("LB","collision_fluids")
  
  if(lbc_halfway)then
    if(ldiagnostic)call start_timing2("LB","apply_bback_pop_hf")
    call driver_apply_bounceback_halfway_pop
    if(ldiagnostic)call end_timing2("LB","apply_bback_pop_hf")
  endif
  
  if (debug) call print_all_pops2(131, "bef_driver_bc_pops", nstep, aoptpR)

  if(lparticles)then
    ! call driver_bc_pops(lparticles)
    call driver_bc_pops_NOK

    if (debug) call print_all_pops2(131, "aft_driver_bc_pops", nstep, aoptpR)

    if(ldiagnostic)call start_timing2("IO","write_xyz")
    call write_xyz(nstep)
    if(ldiagnostic)call end_timing2("IO","write_xyz")
    
    if(ldiagnostic)call start_timing2("LB","apply_part_bback")
    call apply_particle_bounce_back(nstep, debug)
    if(ldiagnostic)call end_timing2("LB","apply_part_bback")
    
    if (debug) call print_all_pops2(131, "aft_part_bb", nstep, aoptpR)

    call merge_particle_force(nstep, debug)
    
    call force_particle_bounce_back(nstep, debug)
    
    call compute_mean_particle_force
    
    call store_old_pos_vel_part

    if(ldiagnostic)call start_timing2("MD","integrate_lf")
    call nve_lf(nstep, debug)
    call merge_particle_energies
    if(ldiagnostic)call end_timing2("MD","integrate_lf")

    call restore_particles

    ! call driver_bc_pops(lparticles)
    call driver_bc_pops_NOK
  endif

  if(ldiagnostic)call start_timing2("IO","write_vtk_frame")
   call write_vtk_frame(nstep)
  if(ldiagnostic)call end_timing2("IO","write_vtk_frame")

  if(ldiagnostic)call start_timing2("LB","streaming_fluids")
  call driver_streaming_fluids(lparticles)
  if(ldiagnostic)call end_timing2("LB","streaming_fluids")
  
  if (debug) call print_all_pops2(131, "bef_apply_bounceback_pop", nstep, aoptpR)
  if(.not.lbc_halfway)then
    if(ldiagnostic)call start_timing2("LB","apply_bback_pop")
    call driver_apply_bounceback_pop
    if(ldiagnostic)call end_timing2("LB","apply_bback_pop")
  endif
  
  if (debug) call print_all_pops2(131, "aft_apply_bounceback_pop", nstep, aoptpR)

  mytime = new_time
  
  return
  
 end subroutine driver_integrator

 end module integrator_mod

