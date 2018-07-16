
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

 use version_mod,     only : idrank
 use error_mod
 use profiling_mod,  only : start_timing2,end_timing2, &
                       ldiagnostic
 use fluids_mod,      only : initialize_fluid_force,compute_fluid_force_sc, &
                        collision_fluids,driver_bc_hvars,driver_bc_pops,&
                        collision_fluids_unique_omega,lunique_omega, &
                        compute_omega_bimix,streaming_fluids, &
                        moments_fluids
 use write_output_mod, only : write_vtk_frame
 
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
 
  new_time = real(nstep,kind=PRC)*tstep
  
  if(ldiagnostic)call start_timing2("LB","moments_fluids")
  call moments_fluids
  if(ldiagnostic)call end_timing2("LB","moments_fluids")
  
  if(ldiagnostic)call start_timing2("LB","driver_bc_hvars")
  call driver_bc_hvars
  if(ldiagnostic)call end_timing2("LB","driver_bc_hvars")
  
  if(ldiagnostic)call start_timing2("LB","initialize_force")
  call initialize_fluid_force
  if(ldiagnostic)call end_timing2("LB","initialize_force")
  
  if(ldiagnostic)call start_timing2("LB","compute_force_sc")
  call compute_fluid_force_sc
  if(ldiagnostic)call end_timing2("LB","compute_force_sc")
  
  if(ldiagnostic)call start_timing2("IO","write_vtk_frame")
  call write_vtk_frame(nstep)
  if(ldiagnostic)call end_timing2("IO","write_vtk_frame")
  
  if(ldiagnostic)call start_timing2("LB","collision_fluids")
  if(lunique_omega)then
    call collision_fluids_unique_omega
  else
    call compute_omega_bimix
    call collision_fluids
  endif
  if(ldiagnostic)call end_timing2("LB","collision_fluids")
  
  if(ldiagnostic)call start_timing2("LB","driver_bc_pops")
  call driver_bc_pops
  if(ldiagnostic)call end_timing2("LB","driver_bc_pops")
  
  if(ldiagnostic)call start_timing2("LB","streaming_fluids")
  call streaming_fluids
  if(ldiagnostic)call end_timing2("LB","streaming_fluids")
  
  mytime = new_time
  
  return
  
 end subroutine driver_integrator

 end module integrator_mod

