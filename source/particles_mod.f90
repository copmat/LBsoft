
#include <default_macro.h>
 module particles_mod
 
!***********************************************************************
!     
!     LBsoft module containing variable and subroutines of the
!     particle managing
!
!     few variable names and subroutines are inspired from
!     DL POLY CLASSICS distributed under BSD licence from
!     the daresbury laboratory (in primis we like to acknowledge 
!     prof. w. smith)
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
 use version_mod, only : idrank,mxrank,or_world_larr,finalize_world, &
                   get_sync_world
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
 
 !key for periodic boundary conditions
 integer, public, protected, save :: imcon=0
 
 !cell vectors
 double precision, dimension(9), public, protected, save :: cell
 
 !key for activate the body rotation
 logical, public, protected, save :: lrotate=.true.
 
 !position of particles
 double precision, allocatable, public, protected, save :: xxx(:)
 double precision, allocatable, public, protected, save :: yyy(:)
 double precision, allocatable, public, protected, save :: zzz(:)
 !components of angular velocity
 double precision, allocatable, public, protected, save :: vxx(:)
 double precision, allocatable, public, protected, save :: vyy(:)
 double precision, allocatable, public, protected, save :: vzz(:)
 !components of force
 double precision, allocatable, public, protected, save :: fxx(:)
 double precision, allocatable, public, protected, save :: fyy(:)
 double precision, allocatable, public, protected, save :: fzz(:)
 !radius of particles
 double precision, allocatable, public, protected, save :: rdim(:)
 !particle mass
 double precision, allocatable, public, protected, save :: weight(:)
 !rotational inertia in body fixed frame
 double precision, allocatable, public, protected, save :: rotin(:)
 !components of angular velocity
 double precision, allocatable, public, protected, save :: oxx(:)
 double precision, allocatable, public, protected, save :: oyy(:)
 double precision, allocatable, public, protected, save :: ozz(:)
 !components of quaternion vector
 double precision, allocatable, public, protected, save :: q0(:)
 double precision, allocatable, public, protected, save :: q1(:)
 double precision, allocatable, public, protected, save :: q2(:)
 double precision, allocatable, public, protected, save :: q3(:)
 !components of torque on rigid body
 double precision, allocatable, public, protected, save :: tqx(:)
 double precision, allocatable, public, protected, save :: tqy(:)
 double precision, allocatable, public, protected, save :: tqz(:)
 
 public :: initialize_particles
 
 contains
 
 subroutine initialize_particles
 
  implicit none
  
  
  
  return
 
 end subroutine initialize_particles
 
 end module particles_mod
