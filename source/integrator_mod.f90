
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
                        aoptpR, driver_bc_pops, print_all_pops2

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

 include 'mpif.h'

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

  real(kind(1.d0)),save :: time1, time2, elap=0,elap1=0
  

  !debug = nstep <= 10
  !debug = modulo(nstep, 500) == 0
  debug = .false.


  new_time = real(nstep,kind=PRC)*tstep
  
  if(ldiagnostic)call start_timing2("LB","collision_fluids")
  call driver_collision_fluids
  if(ldiagnostic)call end_timing2("LB","collision_fluids")

  elap = elap + time2-time1

  if(ldiagnostic)call start_timing2("LB","streaming_fluids")
  call driver_streaming_fluids(lparticles)
  if(ldiagnostic)call end_timing2("LB","streaming_fluids")

  mytime = new_time
  
  return
  
 end subroutine driver_integrator

 end module integrator_mod

