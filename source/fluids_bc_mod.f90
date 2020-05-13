
#include <default_macro.h>

 module fluids_bc_mod
 
!***********************************************************************
!     
!     LBsoft module containing variable and subroutines of the
!     fluid boundary conditions
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification May 2020
!     
!***********************************************************************
 
 use fluids_lattices_mod
 use version_mod,    only : idrank,mxrank,or_world_larr,finalize_world, &
                   get_sync_world,sum_world_farr,commspop,commrpop, &
                   i4find,i4back,ownern,deallocate_ownern,commexch_dens, &
                   commwait_dens,comm_init_isfluid,commexch_vel_component, &
                   commwait_vel_component,commexch_isfluid, &
                   commwait_isfluid,commexch_single_halo, &
                   commwait_single_halo,collective_readFile_pops, &
                   collective_writeFile_pops,collective_readFile, &
                   collective_readFile_int,collective_writeFile, &
                   collective_writeFile_int
 use error_mod
 use utility_mod, only : Pi,modulvec,cross,dot,gauss,ibuffservice, &
                   allocate_array_ibuffservice,buffservice,conv_rad, &
                   allocate_array_buffservice,lbuffservice, &
                   allocate_array_lbuffservice,xcross,ycross,zcross, &
                   nbuffservice3d,buffservice3d,int_cube_sphere,fcut, &
                   rand_noseeded,linit_seed,gauss_noseeded,write_fmtnumb, &
                   openLogFile

 
 implicit none
 
 private
 
 integer, save, protected, public :: ibctype=0
 integer, save, public :: ixpbc, iypbc, izpbc
 
 public :: set_boundary_conditions_type
 
 contains
 
 subroutine set_boundary_conditions_type(itemp1,itemp2,itemp3)
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the ibctype variable taking into
!     account the periodic bc
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification May 2020
!     
!***********************************************************************
   
  implicit none
  
  integer, intent(in) :: itemp1,itemp2,itemp3
  integer ::i,j,k
  integer, dimension(0:1,0:1,0:1), parameter :: iselbct = &
   reshape((/0,1,2,3,4,5,6,7/),(/2,2,2/))
 ! 0=F    1=T  periodic
 !     itemp1      itemp2      itemp3     ibctype
 !          0           0           0           0
 !          1           0           0           1
 !          0           1           0           2
 !          1           1           0           3
 !          0           0           1           4
 !          1           0           1           5
 !          0           1           1           6
 !          1           1           1           7
  
  
  if( itemp1 .ne. 0 .and. itemp1 .ne. 1)then
    call error(9)
  endif
  
  if( itemp2 .ne. 0 .and. itemp2 .ne. 1)then
    call error(9)
  endif
  
  if( itemp3 .ne. 0 .and. itemp3 .ne. 1)then
    call error(9)
  endif
  
  ibctype=iselbct(itemp1,itemp2,itemp3)
  
  return
  
 end subroutine set_boundary_conditions_type
 
 end module fluids_bc_mod
