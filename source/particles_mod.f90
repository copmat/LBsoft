
#include <default_macro.h>
 module particles_mod
 
!***********************************************************************
!     
!     LBsoft module containing variable and subroutines of the
!     particle managing
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
 use version_mod,    only : idrank,mxrank,or_world_larr,finalize_world,get_sync_world
 use error_mod
 use aop_mod
 use utility_mod, only : Pi,modulvec,cross,dot,gauss,ibuffservice, &
                   allocate_array_ibuffservice,buffservice, &
                   allocate_array_buffservice,lbuffservice, &
                   allocate_array_lbuffservice, &
                   nbuffservice3d,buffservice3d, &
                   rand_noseeded,linit_seed,gauss_noseeded,write_fmtnumb

 use lbempi_mod,  only : commspop, commrpop, i4find, i4back, &
                   ownern,deallocate_ownern,commexch_dens, &
                   commwait_dens,comm_init_isfluid,commexch_vel_component, &
                   commwait_vel_component

 
 implicit none
 
 private
 
 public :: initialize_particles
 
 contains
 
 subroutine initialize_particles
 
  implicit none
  
  
  
  return
 
 end subroutine initialize_particles
 
 end module particles_mod
