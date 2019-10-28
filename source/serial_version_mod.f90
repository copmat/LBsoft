
#include <default_macro.h>
 module version_mod
 
!***********************************************************************
!     
!     LBsoft module containing subroutines which replace all the 
!     communication tasks providing a serial version of the code
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 use aop_mod 
 
 implicit none
 
 private
 
 integer, public, save :: idrank=0
 integer, public, save :: mxrank=1
 
#if LATTICE==319
 integer, parameter :: links=18
  !lattice vectors
 integer, dimension(0:links), parameter :: &
  ex = (/ 0, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1, 0, 0, 0, 0/)
   !      0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18
 integer, dimension(0:links), parameter :: &
  ey = (/ 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1/)
 integer, dimension(0:links), parameter :: &
  ez = (/ 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1/)
 integer, dimension(0:links), parameter, public :: &
  opp =(/ 0, 2, 1, 4, 3, 6, 5, 8, 7,10, 9,12,11,14,13,16,15,18,17/)
#endif
 
 integer, protected, public :: nprocz, nprocy, nprocx
 integer, save :: nxy2,nx2,nbuff_comm,nx_comm,ny_comm,nz_comm
 
 integer, protected, public :: domdec
 
 integer, save, dimension(3) :: ipbc,globalDims
 integer(kind=2), public, protected, allocatable :: ownern(:)
 integer, dimension(:), allocatable, public, protected :: gminx,gmaxx, &
   gminy,gmaxy,gminz,gmaxz
   
 logical, save :: ldo_second=.false.
 
 public :: print_version
 public :: get_rank_world
 public :: get_size_world
 public :: get_sync_world
 public :: init_world
 public :: finalize_world
 public :: abort_world
 public :: time_world
 public :: bcast_world_i
 public :: bcast_world_l
 public :: bcast_world_f
 public :: bcast_world_iarr
 public :: bcast_world_larr
 public :: bcast_world_farr
 public :: bcast_world_carr
 public :: sum_world_iarr
 public :: sum_world_farr
 public :: min_world_iarr
 public :: min_world_farr
 public :: max_world_iarr
 public :: max_world_farr
 public :: and_world_larr
 public :: or_world_larr
 public :: or1_world_larr
 public :: irecv_world_i
 public :: irecv_world_f
 public :: irecv_world_l
 public :: isend_world_i
 public :: isend_world_f
 public :: isend_world_l
 public :: irecv_world_iarr
 public :: irecv_world_farr
 public :: irecv_world_larr
 public :: isend_world_iarr
 public :: isend_world_farr
 public :: isend_world_larr
 public :: wait_world
 public :: waitall_world
 public :: sum_world_qarr
 public :: open_file_vtk_par
 public :: close_file_vtk_par
 public :: print_header_vtk_par
 public :: print_footer_vtk_par
 public :: setupcom
 public :: set_domdec
 public :: set_domain
 public :: i4back
 public :: i4find
 public :: deallocate_ownern
 public :: commexch_dens
 public :: commwait_dens
 public :: commexch_vel_component
 public :: commwait_vel_component
 public :: commexch_single_halo
 public :: commwait_single_halo
 public :: commexch_isfluid
 public :: commwait_isfluid
 public :: comm_init_isfluid
 public :: create_findneigh_list_single_halo
 public :: create_findneigh_list_hvar_isfluid
 public :: create_findneigh_list_pops
 public :: ownernfind
 public :: ownernfind_arr
 public :: driving_print_binary_1d_vtk
 public :: driving_print_binary_3d_vtk
 public :: collective_readFile_pops
 public :: collective_writeFile_pops
 
 public :: collective_readFile
 public :: collective_readFile_int
 public :: collective_writeFile
 public :: collective_writeFile_int
 
 public :: commspop
 public :: commrpop
 
 contains
 
 subroutine print_version(iu)
 
!***********************************************************************
!     
!     LBsoft subroutine for printing the code version which is
!     currently in use
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  character(len=50), intent(out) :: iu
  
  write(iu,'(a50)')'The code is running in serial mode                '
  
  return
  
 end subroutine print_version
 
 subroutine get_rank_world()
 
!***********************************************************************
!     
!     LBsoft subroutine to determine identity of processing node 
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  idrank=0
  
  return
  
 end subroutine get_rank_world
 
 subroutine get_size_world()
 
!***********************************************************************
!     
!     LBsoft subroutine to determine the number of processing nodes
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  mxrank=1
  
  return
  
 end subroutine get_size_world
 
 subroutine get_sync_world()
 
!***********************************************************************
!     
!     LBsoft subroutine to synchronize the processing nodes
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  return
  
 end subroutine get_sync_world
 
 subroutine init_world()
 
!***********************************************************************
!     
!     LBsoft subroutine to initialize the MPI work
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************

  implicit none
  
  return
 
 end subroutine init_world
 
 subroutine finalize_world()
 
!***********************************************************************
!     
!     LBsoft subroutine to finalize the MPI work
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  return
  
 end subroutine finalize_world
 
 subroutine abort_world()
 
!***********************************************************************
!     
!     LBsoft subroutine to abort the MPI work
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  integer ier
  
  return
  
 end subroutine abort_world
 
 subroutine time_world(timecpu)
 
!***********************************************************************
!     
!     LBsoft subroutine to use the CPU clock
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC), intent(out) :: timecpu
  
  call cpu_time(timecpu)
  
  return
  
 end subroutine time_world
 
 subroutine bcast_world_i(buffer)
 
!***********************************************************************
!     
!     LBsoft subroutine to broadcast an integer number to all other 
!     nodes
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(inout) :: buffer
  
  return
  
 end subroutine bcast_world_i
 
 subroutine bcast_world_l(buffer)
 
!***********************************************************************
!     
!     LBsoft subroutine to broadcast a logical variable to all other 
!     nodes
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  logical, intent(inout) :: buffer
  
  return
  
 end subroutine bcast_world_l
 
 subroutine bcast_world_f(buffer)
 
!***********************************************************************
!     
!     LBsoft subroutine to broadcast a float number to all 
!     other nodes
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC), intent(inout) :: buffer
  
  return
  
 end subroutine bcast_world_f
 
 subroutine bcast_world_iarr(buffer,narr)
 
!***********************************************************************
!     
!     LBsoft subroutine to broadcast an integer array to all other 
!     nodes
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(inout), dimension(narr) :: buffer
  integer, intent(in) :: narr
  
  return
  
 end subroutine bcast_world_iarr
 
 subroutine bcast_world_larr(buffer,narr)
 
!***********************************************************************
!     
!     LBsoft subroutine to broadcast a logical array to all other 
!     nodes
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  logical, intent(inout), dimension(narr) :: buffer
  integer, intent(in) :: narr
  
  return
  
 end subroutine bcast_world_larr
 
 subroutine bcast_world_farr(buffer,narr)
 
!***********************************************************************
!     
!     LBsoft subroutine to broadcast a float array to all 
!     other nodes
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC), intent(inout), dimension(narr) :: buffer
  integer, intent(in) :: narr
  
  return
  
 end subroutine bcast_world_farr
 
 subroutine bcast_world_carr(mxlen,argument,narr)
 
!***********************************************************************
!     
!     LBsoft subroutine to broadcast a character array to all 
!     other nodes
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: mxlen
  character(len=mxlen), intent(in), dimension(narr) :: argument
  integer, intent(in) :: narr
  
  return
  
 end subroutine bcast_world_carr
 
 subroutine sum_world_iarr(argument,narr,buffersub)
 
!***********************************************************************
!     
!     LBsoft global summation subroutine for a integer array
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  integer, intent(inout), dimension(narr), optional :: buffersub
  
  integer ier
  
  if(present(buffersub))buffersub(1:narr)=argument(1:narr)
  
  return
  
 end subroutine sum_world_iarr
 
 subroutine sum_world_farr(argument,narr,buffersub)
 
!***********************************************************************
!     
!     LBsoft global summation subroutine for a float array
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC), intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  real(kind=PRC), intent(inout), dimension(narr), optional :: buffersub
  
  integer ier
  
  if(present(buffersub))buffersub(1:narr)=argument(1:narr)
  
  return
  
 end subroutine sum_world_farr
 
 subroutine min_world_iarr(argument,narr,buffersub)
 
!***********************************************************************
!     
!     LBsoft global minimum subroutine for an integer array
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  integer, intent(inout), dimension(narr), optional :: buffersub
  
  integer ier
  
  if(present(buffersub))buffersub(1:narr)=argument(1:narr)
  
  return
  
 end subroutine min_world_iarr
 
 subroutine min_world_farr(argument,narr,buffersub)
 
!***********************************************************************
!     
!     LBsoft global minimum subroutine for a float array
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC), intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  real(kind=PRC), intent(inout), dimension(narr), optional :: buffersub
  
  integer ier
  
  if(present(buffersub))buffersub(1:narr)=argument(1:narr)
  
  return
  
 end subroutine min_world_farr
 
 subroutine max_world_iarr(argument,narr,buffersub)
 
!***********************************************************************
!     
!     LBsoft global maximum subroutine for an integer array
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  integer, intent(inout), dimension(narr), optional :: buffersub
  
  integer ier
  
  if(present(buffersub))buffersub(1:narr)=argument(1:narr)
  
  return
  
 end subroutine max_world_iarr
 
 subroutine max_world_farr(argument,narr,buffersub)
 
!***********************************************************************
!     
!     LBsoft global maximum subroutine for a float array
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC), intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  real(kind=PRC), intent(inout), dimension(narr), optional :: buffersub
  
  integer ier
  
  if(present(buffersub))buffersub(1:narr)=argument(1:narr)
  
  return
  
 end subroutine max_world_farr
 
 subroutine and_world_larr(argument,narr,buffersub)
 
!***********************************************************************
!     
!     LBsoft global 'logical and' subroutine for a logical array
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  logical, intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  logical, intent(inout), dimension(narr), optional :: buffersub
  
  integer ier
  
  if(present(buffersub))buffersub(1:narr)=argument(1:narr)
  
  return
  
 end subroutine and_world_larr
 
 subroutine or_world_larr(argument,narr,buffersub)
 
!***********************************************************************
!     
!     LBsoft global 'logical or' subroutine for a logical array
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  logical, intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  logical, intent(inout), dimension(narr), optional :: buffersub
  
  integer ier
  
  if(present(buffersub))buffersub(1:narr)=argument(1:narr)
  
  return
  
 end subroutine or_world_larr
 
 subroutine or1_world_larr(argument,narr)
 
!***********************************************************************
!     
!     LBsoft global 'logical or' subroutine for a logical array
!     of kind 1
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  logical(kind=1), intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  
  integer ier
  
  return
  
 end subroutine or1_world_larr
 
 subroutine irecv_world_i(argument,idsource,tag,irequest)
 
!***********************************************************************
!     
!     LBsoft global subroutine for receiving an integer value
!     by a nonblocking receive operation
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(out) :: argument
  integer, intent(in) :: idsource,tag
  integer, intent(out) :: irequest
  integer :: ierr
  
  return
  
 end subroutine irecv_world_i
 
 subroutine irecv_world_f(argument,idsource,tag,irequest)
 
!***********************************************************************
!     
!     LBsoft global subroutine for receiving a float value
!     by a nonblocking receive operation
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC), intent(out) :: argument
  integer, intent(in) :: idsource,tag
  integer, intent(out) :: irequest
  integer :: ierr
  
  return
  
 end subroutine irecv_world_f
 
 subroutine irecv_world_l(argument,idsource,tag,irequest)
 
!***********************************************************************
!     
!     LBsoft global subroutine for receiving a logical value
!     by a nonblocking receive operation
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
  
  implicit none
  
  logical, intent(out) :: argument
  integer, intent(in) :: idsource,tag
  integer, intent(out) :: irequest
  integer :: ierr
  
  return
  
 end subroutine irecv_world_l
 
 subroutine isend_world_i(argument,iddest,tag,irequest)
 
!***********************************************************************
!     
!     LBsoft global subroutine for sending an integer value
!     by a nonblocking send operation
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: argument
  integer, intent(in) :: iddest,tag
  integer, intent(out) :: irequest
  integer :: ierr
  
  return
  
 end subroutine isend_world_i
 
 subroutine isend_world_f(argument,iddest,tag,irequest)
 
!***********************************************************************
!     
!     LBsoft global subroutine for sending a float value
!     by a nonblocking send operation
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC), intent(in) :: argument
  integer, intent(in) :: iddest,tag
  integer, intent(out) :: irequest
  integer :: ierr
  
  return
  
 end subroutine isend_world_f
 
 subroutine isend_world_l(argument,iddest,tag,irequest)
 
!***********************************************************************
!     
!     LBsoft global subroutine for sending a logical value
!     by a nonblocking send operation
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
  
  implicit none
  
  logical, intent(in) :: argument
  integer, intent(in) :: iddest,tag
  integer, intent(out) :: irequest
  integer :: ierr
  
  return
  
 end subroutine isend_world_l
 
 subroutine irecv_world_iarr(argument,narr,idsource,tag,irequest)
 
!***********************************************************************
!     
!     LBsoft global subroutine for receiving an integer array
!     by a nonblocking receive operation
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: narr
  integer, dimension(narr), intent(out) :: argument
  integer, intent(in) :: idsource,tag
  integer, intent(out) :: irequest
  integer :: ierr
  
  return
  
 end subroutine irecv_world_iarr
 
 subroutine irecv_world_farr(argument,narr,idsource,tag,irequest)
 
!***********************************************************************
!     
!     LBsoft global subroutine for receiving a float array
!     by a nonblocking receive operation
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: narr
  real(kind=PRC), dimension(narr), intent(out) :: argument
  integer, intent(in) :: idsource,tag
  integer, intent(out) :: irequest
  integer :: ierr
  
  return
  
 end subroutine irecv_world_farr
 
 subroutine irecv_world_larr(argument,narr,idsource,tag,irequest)
 
!***********************************************************************
!     
!     LBsoft global subroutine for receiving a logical array
!     by a nonblocking receive operation
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: narr
  logical, dimension(narr), intent(out) :: argument
  integer, intent(in) :: idsource,tag
  integer, intent(out) :: irequest
  integer :: ierr
  
  return
  
 end subroutine irecv_world_larr
 
 subroutine isend_world_iarr(argument,narr,iddest,tag,irequest)
 
!***********************************************************************
!     
!     LBsoft global subroutine for sending an integer array
!     by a nonblocking send operation
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: narr
  integer, dimension(narr), intent(in) :: argument
  integer, intent(in) :: iddest,tag
  integer, intent(out) :: irequest
  integer :: ierr
  
  return
  
 end subroutine isend_world_iarr
 
 subroutine isend_world_farr(argument,narr,iddest,tag,irequest)
 
!***********************************************************************
!     
!     LBsoft global subroutine for sending a float array
!     by a nonblocking send operation
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: narr
  real(kind=PRC), dimension(narr), intent(in) :: argument
  integer, intent(in) :: iddest,tag
  integer, intent(out) :: irequest
  integer :: ierr
  
  return
  
 end subroutine isend_world_farr
 
 subroutine isend_world_larr(argument,narr,iddest,tag,irequest)
 
!***********************************************************************
!     
!     LBsoft global subroutine for sending a logical array
!     by a nonblocking send operation
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: narr
  logical, dimension(narr), intent(in) :: argument
  integer, intent(in) :: iddest,tag
  integer, intent(out) :: irequest
  integer :: ierr
  
  return
  
 end subroutine isend_world_larr
 
 subroutine wait_world(argument)
 
!***********************************************************************
!     
!     LBsoft global subroutine for waiting a nonblocking operation 
!     to complete.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(inout) :: argument
  integer :: ierr
  
  return
  
 end subroutine wait_world
 
 subroutine waitall_world(argument,narr)
 
!***********************************************************************
!     
!     LBsoft global subroutine for waiting a collection of nonblocking  
!     operations to complete
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: narr
  integer, dimension(narr), intent(inout) :: argument
  integer :: ierr
  
  return
  
 end subroutine waitall_world
 
 subroutine sum_world_qarr(argument,narr,buffersub)
 
!***********************************************************************
!     
!     LBsoft global summation subroutine for a float array
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC*2), intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  real(kind=PRC*2), intent(inout), dimension(narr), optional :: buffersub
  
  integer ier
  
  if(present(buffersub))buffersub(1:narr)=argument(1:narr)
  
  return
  
 end subroutine sum_world_qarr
 
 subroutine open_file_vtk_par(iotest,nn,myname,E_IO)
 
!***********************************************************************
!     
!     LBsoft subroutine for opening the vtk legacy file
!     in serial IO
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: iotest,nn
  character(len=nn) :: myname
  integer, intent(out) :: e_io
  
  return
  
 endsubroutine open_file_vtk_par
 
 subroutine close_file_vtk_par(iotest,e_io)
 
!***********************************************************************
!     
!     LBsoft subroutine for closing the vtk legacy file
!     in serial IO
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: iotest
  integer, intent(out) :: e_io
  
  return
  
 endsubroutine close_file_vtk_par
 
 subroutine print_header_vtk_par(iotest,offsetsub,nn,header,E_IO)
 
!***********************************************************************
!     
!     LBsoft subroutine for writing the header part of
!     in VTK legacy file in serial IO
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: iotest,offsetsub,nn
  character(len=500) :: header
  integer, intent(out) :: E_IO
  
  return
  
 endsubroutine print_header_vtk_par
 
 subroutine print_footer_vtk_par(iotest,offsetsub,footer,E_IO)
 
!***********************************************************************
!     
!     LBsoft subroutine for writing the footer part of
!     in VTK legacy file in serial IO
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: iotest,offsetsub
  character(len=30) :: footer
  integer, intent(out) :: E_IO
  
  return
  
 endsubroutine print_footer_vtk_par
 
 subroutine setupcom(nx,ny,nz,nbuff,ibctype,ixpbc,iypbc,izpbc,minx,maxx, &
   miny,maxy,minz,maxz,lsingle_fluid)
  
!***********************************************************************
!     
!     LBsoft subroutine for the initial setup of the domain 
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Bernaschi
!     last modification October 2019
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: nx,ny,nz,nbuff,ibctype
  integer, intent(inout) :: ixpbc,iypbc,izpbc,minx,maxx,miny,maxy, &
   minz,maxz
  
  LOGICAL, intent(IN) :: lsingle_fluid
  
  integer :: i
  integer :: ownernlb,ownernub

  
  nx_comm=nx
  ny_comm=ny
  nz_comm=nz
  
  globalDims(1)=nx
  globalDims(2)=ny
  globalDims(3)=nz

  nxy2=(nx+2*nbuff)*(ny+2*nbuff)
  nx2=nx+(2*nbuff)
  
  nbuff_comm=nbuff
  
  minx=1
  maxx=nx
  miny=1
  maxy=ny
  minz=1
  maxz=nz
  
  ldo_second=(.not. lsingle_fluid)
  
  
    !     itemp1      itemp2      itemp3     ibctype
    !          0           0           0           0
    !          1           0           0           1
    !          0           1           0           2
    !          1           1           0           3
    !          0           0           1           4
    !          1           0           1           5
    !          0           1           1           6
    !          1           1           1           7

    

    ! ibctype 0 corresponds to a lattice with no PBC
    ixpbc=0
    iypbc=0
    izpbc=0
    IF(ibctype.EQ.1) THEN
       ixpbc=1
    ENDIF
    IF(ibctype.EQ.2) THEN
       iypbc=1
    ENDIF
    IF(ibctype.EQ.3) THEN
       ixpbc=1
       iypbc=1
    ENDIF
    IF(ibctype.EQ.4) THEN
       izpbc=1
    ENDIF
    IF(ibctype.EQ.5) THEN
       ixpbc=1
       izpbc=1
    ENDIF
    IF(ibctype.EQ.6) THEN
       iypbc=1
       izpbc=1
    ENDIF
    IF(ibctype.EQ.7) THEN
       ixpbc=1
       iypbc=1
       izpbc=1
    ENDIF
    
  ipbc(1)=ixpbc
  ipbc(2)=iypbc
  ipbc(3)=izpbc
  
  !ownernlb=(1-nbuff)+((1-nbuff)*nx2)+((1-nbuff)*nxy2)
  ownernlb=1
  ownernub=(nx+2*nbuff)*(ny+2*nbuff)*(nz+2*nbuff)
  
  ALLOCATE(ownern(ownernlb:ownernub))
  do i=ownernlb,ownernub
     ownern(i)=0
  enddo
  
  allocate(gminx(-1:mxrank-1))
  allocate(gmaxx(-1:mxrank-1))
  allocate(gminy(-1:mxrank-1))
  allocate(gmaxy(-1:mxrank-1))
  allocate(gminz(-1:mxrank-1))
  allocate(gmaxz(-1:mxrank-1))
  
  
  gminx(-1:mxrank-1)=0
  gmaxx(-1:mxrank-1)=0
  gminy(-1:mxrank-1)=0
  gmaxy(-1:mxrank-1)=0
  gminz(-1:mxrank-1)=0
  gmaxz(-1:mxrank-1)=0
     
  do i=0,mxrank-1
    if(i==idrank)then
      gminx(i)=minx
      gmaxx(i)=maxx
      gminy(i)=miny
      gmaxy(i)=maxy
      gminz(i)=minz
      gmaxz(i)=maxz
    endif
  enddo
  call sum_world_iarr(gminx,mxrank+1)
  call sum_world_iarr(gmaxx,mxrank+1)
  call sum_world_iarr(gminy,mxrank+1)
  call sum_world_iarr(gmaxy,mxrank+1)
  call sum_world_iarr(gminz,mxrank+1)
  call sum_world_iarr(gmaxz,mxrank+1)
  
  !here store the pbc condition
  gminx(-1)=ixpbc
  gminy(-1)=iypbc
  gminz(-1)=izpbc
  !here store the total dimension
  gmaxx(-1)=nx
  gmaxy(-1)=ny
  gmaxz(-1)=nz
  
  return
  
 end subroutine setupcom
 
 subroutine set_domdec(itemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the protected value of domdec
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: itemp
  
  domdec=itemp
  
  return
  
 end subroutine set_domdec
 
  subroutine set_domain(itemp1,itemp2,itemp3)
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the domain decomposition in input
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: itemp1,itemp2,itemp3
  
  nprocx=itemp1
  nprocy=itemp2
  nprocz=itemp3
  
  return
  
 end subroutine set_domain
 
 subroutine create_findneigh_list_single_halo(nx,ny,nz,nbuff,ibctype, &
   ixpbc,iypbc,izpbc,minx,maxx,miny,maxy,minz,maxz)

!***********************************************************************
!     
!     LBsoft subroutine for driving the creation of the list of 
!     the HVAR exchange by MPI sending and receiving
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************

  implicit none
 
  INTEGER, intent(in) :: nx,ny,nz,nbuff,ibctype
  INTEGER, intent(inout) :: ixpbc,iypbc,izpbc,minx,maxx,miny,maxy, &
   minz,maxz
   
  return
  
 end subroutine create_findneigh_list_single_halo
 
 subroutine create_findneigh_list_hvar_isfluid(nx,ny,nz,nbuff,ibctype, &
   ixpbc,iypbc,izpbc,minx,maxx,miny,maxy,minz,maxz)

!***********************************************************************
!     
!     LBsoft subroutine for driving the creation of the list of 
!     the HVAR exchange by MPI sending and receiving
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************

  implicit none
 
  INTEGER, intent(in) :: nx,ny,nz,nbuff,ibctype
  INTEGER, intent(inout) :: ixpbc,iypbc,izpbc,minx,maxx,miny,maxy, &
   minz,maxz
   
  return
  
 end subroutine create_findneigh_list_hvar_isfluid

 subroutine create_findneigh_list_pops(nx,ny,nz,nbuff,ibctype,isfluid, &
   ixpbc,iypbc,izpbc,minx,maxx,miny,maxy,minz,maxz)

!***********************************************************************
!     
!     LBsoft subroutine for driving the creation of the list of 
!     the populations exchange by MPI sending and receiving
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
  
  implicit none
 
  INTEGER, intent(in) :: nx,ny,nz,nbuff,ibctype
  INTEGER(kind=1), intent(in), allocatable, dimension(:,:,:) :: isfluid
  INTEGER, intent(inout) :: ixpbc,iypbc,izpbc,minx,maxx,miny,maxy, &
   minz,maxz
   
  return
  
 end subroutine
 
 subroutine comm_init_isfluid(temp)

!***********************************************************************
!     
!     LBsoft subroutine to exchange ISFLUID
!     amond the MPI ranks
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************

  IMPLICIT NONE

  INTEGER(kind=1), dimension(:,:,:), allocatable :: temp
  
  return
  
 end subroutine comm_init_isfluid
 
 subroutine commexch_isfluid(temp)
 
!***********************************************************************
!     
!     LBsoft subroutine for exchanging ISFLUID variables
!     among MPI ranks
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
   
   IMPLICIT NONE
   
   INTEGER(kind=1), dimension(:,:,:), allocatable :: temp
   
   return
   
   END subroutine commexch_isfluid
   
 subroutine commwait_isfluid(temp)
    
!***********************************************************************
!     
!     LBsoft subroutine for waiting the exchange ISFLUID 
!     variables among MPI ranks
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
   
   IMPLICIT NONE
   
   INTEGER(kind=1), dimension(:,:,:), allocatable :: temp
   
   return
   
 END subroutine commwait_isfluid
 
 subroutine commexch_dens(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for exchanging hydro-dynamics variables
!     among MPI ranks
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
   
   IMPLICIT NONE
   
   REAL(KIND=PRC), dimension(:,:,:), allocatable :: dtemp
   
   return
   
  END subroutine commexch_dens
  
 subroutine commwait_dens(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for waiting the exchange hydro-dynamics 
!     variables among MPI ranks
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
   
   IMPLICIT NONE
   
   REAL(KIND=PRC), dimension(:,:,:), allocatable :: dtemp
   
   return
   
 END subroutine commwait_dens
 
 subroutine commexch_vel_component(dtemp)

!***********************************************************************
!     
!     LBsoft subroutine for exchanging velocities variables
!     among MPI ranks
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
   
   IMPLICIT NONE
   
   REAL(KIND=PRC), dimension(:,:,:), allocatable :: dtemp
   
   return
   
 END SUBROUTINE commexch_vel_component
   
 SUBROUTINE commwait_vel_component(dtemp)

!***********************************************************************
!     
!     LBsoft subroutine for waiting the exchange velocities 
!     variables among MPI ranks
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
   
   IMPLICIT NONE
   
   REAL(KIND=PRC), dimension(:,:,:), allocatable :: dtemp
   
   return
   
 END SUBROUTINE commwait_vel_component
 
 subroutine commexch_single_halo(dtemp)

!***********************************************************************
!     
!     LBsoft subroutine for exchanging velocities variables
!     among MPI ranks
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
   
   IMPLICIT NONE
   
   REAL(KIND=PRC), dimension(:,:,:), allocatable :: dtemp
   
   return
   
 END SUBROUTINE commexch_single_halo
   
 SUBROUTINE commwait_single_halo(dtemp)

!***********************************************************************
!     
!     LBsoft subroutine for waiting the exchange velocities 
!     variables among MPI ranks
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
   
   IMPLICIT NONE
   
   REAL(KIND=PRC), dimension(:,:,:), allocatable :: dtemp
   
   return
   
 END SUBROUTINE commwait_single_halo

 pure function ownernfind_arr(i,j,k,nxs,nys,nzs,nbuffs,ownerns)
  
!***********************************************************************
!     
!     LBsoft subroutine for finding the owner of a node provided 
!     in i,j,k triplet-form taken from ownerns array
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
  
 implicit none
   
  integer, intent(in) :: i,j,k,nxs,nys,nzs,nbuffs
  integer(kind=2), allocatable, intent(in) :: ownerns(:)
  integer(kind=IPRC) :: i4
   
  integer :: ownernfind_arr
    
  ownernfind_arr=0
   
  return
   
 end function ownernfind_arr
 
 pure function ownernfind(i,j,k,mxranksub,gminxs,gmaxxs, &
     gminys,gmaxys,gminzs,gmaxzs) result(i4find)
 
!***********************************************************************
!     
!     LBsoft subroutine for finding the owner of a node provided 
!     in i,j,k triplet-form computed on-fly
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none
  integer,intent(in) :: i !< i-location on mesh
  integer,intent(in) :: j !< j-location on mesh
  integer,intent(in) :: k !< k-location on mesh
  integer,intent(in) :: mxranksub
  !this is necessary to make it pure
  integer,intent(in), dimension(-1:mxranksub-1) :: gminxs,gmaxxs, &
   gminys,gmaxys,gminzs,gmaxzs
     
  integer :: i4find
  integer :: l,itemp,jtemp,ktemp
  logical :: ltest(1),lt(3)
    
  itemp=i
  jtemp=j
  ktemp=k
    
    
  i4find=0
    
  return
  
 end function ownernfind
 
 subroutine commspop(pop_ptr)
 
!***********************************************************************
!     
!     LBsoft subroutine for sending the population over the MPI
!     network
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************

  implicit none

  type(REALPTR), dimension(0:links):: pop_ptr
  
  if(mxrank==1)return
  
  return

 end subroutine commspop
 
 subroutine commrpop(pop_ptr,lparticles,isfluid)
 
!***********************************************************************
!     
!     LBsoft subroutine for receiving the population over the MPI
!     network
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************

  implicit none

  type(REALPTR), dimension(0:links):: pop_ptr
  logical, intent(in) :: lparticles
  integer(kind=1), allocatable, dimension(:,:,:), intent(in) :: isfluid
  
  if(mxrank==1)return
 

  return
 
 end subroutine commrpop
 
 function i4back(i,j,k) result(myout)
  
!***********************************************************************
!     
!     LBsoft subroutine for finding a unique integer number i4-form
!     for each i,j,k triplet
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
  
  implicit none
  
  integer,INTENT(in) :: i !< i-location on mesh
  integer,INTENT(in) :: j !< j-location on mesh
  integer,INTENT(in) :: k !< k-location on mesh
  integer(kind=IPRC) :: myout

  myout = int(k+nbuff_comm-1,kind=IPRC)*nxy2 + &
   int(j+nbuff_comm-1,kind=IPRC)*nx2 + int(i+nbuff_comm,kind=IPRC)
  
  return
  
 end function i4back
 
 function i4find(i4sub)

!***********************************************************************
!     
!     LBsoft subroutine for returning the i,j,k triplet-form of 
!     a i4-form
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
  
  implicit none
  integer(kind=IPRC), intent(in) :: i4sub
  integer, dimension(3) :: i4find
    
  i4find(3) = ((i4sub-1)/nxy2)-nbuff_comm+1
  i4find(2) = ((i4sub-1 - (i4find(3)+nbuff_comm-1)*nxy2)/nx2)-nbuff_comm+1
  i4find(1) = (i4sub-nbuff_comm+1 - (i4find(3)+nbuff_comm-1)*nxy2 - (i4find(2)+nbuff_comm-1)*nx2)-1
    
  return
  
 end function i4find
 
 SUBROUTINE deallocate_ownern
  
!***********************************************************************
!     
!     LBsoft subroutine to deallocate ownern array
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification September 2018
!     
!***********************************************************************
  
   implicit none
  
   DEALLOCATE(ownern)
  
   return

 end subroutine deallocate_ownern
 
 subroutine collective_writeFile(it, mynamefile, myvar, nbuff)
  
!***********************************************************************
!     
!     LBsoft subroutine for writing a scalar field with PRC
!     precision in binary format in serial IO
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
  
  implicit none
  REAL(KIND=PRC), dimension(:,:,:), allocatable :: myvar
  integer, intent(in) :: it, nbuff
  character(len=120),intent(in) :: mynamefile
  integer :: ierr
  
  integer, parameter :: fh=914
  
  open(unit=fh,file=trim(mynamefile),action='write', &
   status='replace',form='unformatted',access='stream',iostat=ierr)
  
  call write_binary_1d_bynary(fh,globalDims,myvar)
  
  close(unit=fh,iostat=ierr)
  
  return
  
 end subroutine collective_writeFile

 subroutine collective_writeFile_int(it, mynamefile, myvar, nbuff)
 
!***********************************************************************
!     
!     LBsoft subroutine for writing a scalar field with INTEGER
!     of 1 byte in binary format in serial IO
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none
  integer(kind=1), allocatable, dimension(:,:,:) :: myvar
  integer, intent(in) :: it, nbuff
  character(len=120),intent(in) :: mynamefile
  integer :: ierr

  integer, parameter :: fh=915
  
  open(unit=fh,file=trim(mynamefile),action='write', &
   status='replace',form='unformatted',access='stream',iostat=ierr)
  
  call write_binary_1d_bynary_int(fh,globalDims,myvar)
  
  close(unit=fh,iostat=ierr)
  
  return
  
 end subroutine collective_writeFile_int
 
 subroutine collective_readFile(it, mynamefile, myvar, nbuff)

!***********************************************************************
!     
!     LBsoft subroutine for reading a scalar field with PRC
!     precision in binary format in serial IO
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************

  implicit none
  REAL(KIND=PRC), dimension(:,:,:), allocatable :: myvar
  integer, intent(in) :: it, nbuff
  character(len=120),intent(in) :: mynamefile
  integer :: ierr
  
  integer, parameter :: fh=912
  
  open(unit=fh,file=trim(mynamefile),action='read', &
   status='old',form='unformatted',access='stream',iostat=ierr)
  
  call read_binary_1d_bynary(fh,globalDims,myvar)
  
  close(unit=fh,iostat=ierr)
  
  return
 
 end subroutine collective_readFile

 subroutine collective_readFile_int(it, mynamefile, myvar, nbuff)
    
!***********************************************************************
!     
!     LBsoft subroutine for reading a scalar field with INTEGER
!     of 1 byte in binary format in serial IO
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
  implicit none
  integer(kind=1), allocatable, dimension(:,:,:) :: myvar
  integer, intent(in) :: it, nbuff
  character(len=120),intent(in) :: mynamefile
  integer :: ierr
  
  integer, parameter :: fh=913
  
  open(unit=fh,file=trim(mynamefile),action='read', &
   status='old',form='unformatted',access='stream',iostat=ierr)
  
  call read_binary_1d_bynary_int(fh,globalDims,myvar)
  
  close(unit=fh,iostat=ierr)
  
  return
  
 end subroutine collective_readFile_int
 
 subroutine collective_readFile_pops(it, mynamefile, f00,f01, &
  f02,f03,f04,f05,f06,f07,f08,f09,f10,f11,f12,f13,f14,f15,f16,f17,f18, &
  nbuff)

!***********************************************************************
!     
!     LBsoft subroutine for reading the pops with PRC
!     precision in binary format in serial IO
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
  
  implicit none
  REAL(kind=PRC), allocatable, dimension(:,:,:)  :: f00,f01, &
   f02,f03,f04,f05,f06,f07,f08,f09,f10,f11,f12,f13,f14,f15,f16,f17,f18
  integer, intent(in) :: it, nbuff
  character(len=120),intent(in) :: mynamefile
  integer :: ierr
  
  integer, parameter :: fh=911
  
  open(unit=fh,file=trim(mynamefile),action='read', &
   status='old',form='unformatted',access='stream',iostat=ierr)
  
  
  call read_binary_1d_bynary(fh,globalDims,f00)
  call read_binary_1d_bynary(fh,globalDims,f01)
  call read_binary_1d_bynary(fh,globalDims,f02)
  call read_binary_1d_bynary(fh,globalDims,f03)
  call read_binary_1d_bynary(fh,globalDims,f04)
  call read_binary_1d_bynary(fh,globalDims,f05)
  call read_binary_1d_bynary(fh,globalDims,f06)
  call read_binary_1d_bynary(fh,globalDims,f07)
  call read_binary_1d_bynary(fh,globalDims,f08)
  call read_binary_1d_bynary(fh,globalDims,f09)
  call read_binary_1d_bynary(fh,globalDims,f10)
  call read_binary_1d_bynary(fh,globalDims,f11)
  call read_binary_1d_bynary(fh,globalDims,f12)
  call read_binary_1d_bynary(fh,globalDims,f13)
  call read_binary_1d_bynary(fh,globalDims,f14)
  call read_binary_1d_bynary(fh,globalDims,f15)
  call read_binary_1d_bynary(fh,globalDims,f16)
  call read_binary_1d_bynary(fh,globalDims,f17)
  call read_binary_1d_bynary(fh,globalDims,f18)
  

  close(unit=fh,iostat=ierr)
  
  return
  
 end subroutine collective_readFile_pops
 
 subroutine read_binary_1d_bynary(myio,mydim,myvar)
 
!***********************************************************************
!     
!     LBsoft subroutine for reading a scalar field with PRC
!     precision in binary format in serial IO
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: myio
  integer, intent(in), dimension(3) :: mydim
  real(kind=PRC), allocatable, dimension(:,:,:)  :: myvar
  
  integer :: i,j,k
  
  do k=1,mydim(3)
    do j=1,mydim(2)
      do i=1,mydim(1)
        read(myio)myvar(i,j,k)
      enddo
    enddo
  enddo
  
  return
 
 end subroutine read_binary_1d_bynary
 
 subroutine read_binary_1d_bynary_int(myio,mydim,myvar)
 
!***********************************************************************
!     
!     LBsoft subroutine for reading a scalar field of INTEGER
!     with 1 byte in binary format in serial IO
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: myio
  integer, intent(in), dimension(3) :: mydim
  integer(kind=1), allocatable, dimension(:,:,:)  :: myvar
  
  integer :: i,j,k
  
  do k=1,mydim(3)
    do j=1,mydim(2)
      do i=1,mydim(1)
        read(myio)myvar(i,j,k)
      enddo
    enddo
  enddo
  
  return
 
 end subroutine read_binary_1d_bynary_int
 
 subroutine collective_writeFile_pops(it, mynamefile, f00,f01, &
  f02,f03,f04,f05,f06,f07,f08,f09,f10,f11,f12,f13,f14,f15,f16,f17,f18, &
  nbuff)

!***********************************************************************
!     
!     LBsoft subroutine for writing the pops with PRC
!     precision in binary format in serial IO
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
  
  implicit none
  REAL(kind=PRC), allocatable, dimension(:,:,:)  :: f00,f01, &
   f02,f03,f04,f05,f06,f07,f08,f09,f10,f11,f12,f13,f14,f15,f16,f17,f18
  integer, intent(in) :: it, nbuff
  character(len=120),intent(in) :: mynamefile
  integer :: ierr
  
  integer, parameter :: fh=911
  
  open(unit=fh,file=trim(mynamefile),action='write', &
   status='replace',form='unformatted',access='stream',iostat=ierr)
  
  
  call write_binary_1d_bynary(fh,globalDims,f00)
  call write_binary_1d_bynary(fh,globalDims,f01)
  call write_binary_1d_bynary(fh,globalDims,f02)
  call write_binary_1d_bynary(fh,globalDims,f03)
  call write_binary_1d_bynary(fh,globalDims,f04)
  call write_binary_1d_bynary(fh,globalDims,f05)
  call write_binary_1d_bynary(fh,globalDims,f06)
  call write_binary_1d_bynary(fh,globalDims,f07)
  call write_binary_1d_bynary(fh,globalDims,f08)
  call write_binary_1d_bynary(fh,globalDims,f09)
  call write_binary_1d_bynary(fh,globalDims,f10)
  call write_binary_1d_bynary(fh,globalDims,f11)
  call write_binary_1d_bynary(fh,globalDims,f12)
  call write_binary_1d_bynary(fh,globalDims,f13)
  call write_binary_1d_bynary(fh,globalDims,f14)
  call write_binary_1d_bynary(fh,globalDims,f15)
  call write_binary_1d_bynary(fh,globalDims,f16)
  call write_binary_1d_bynary(fh,globalDims,f17)
  call write_binary_1d_bynary(fh,globalDims,f18)
  

  close(unit=fh,iostat=ierr)
  
  return
  
 end subroutine collective_writeFile_pops
 
 subroutine write_binary_1d_bynary(myio,mydim,myvar)
 
!***********************************************************************
!     
!     LBsoft subroutine for writing a scalar field with PRC
!     precision in binary format in serial IO
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: myio
  integer, intent(in), dimension(3) :: mydim
  real(kind=PRC), allocatable, dimension(:,:,:)  :: myvar
  
  integer :: i,j,k
  
  do k=1,mydim(3)
    do j=1,mydim(2)
      do i=1,mydim(1)
        write(myio)myvar(i,j,k)
      enddo
    enddo
  enddo
  
  return
 
 end subroutine write_binary_1d_bynary
 
 subroutine write_binary_1d_bynary_int(myio,mydim,myvar)
 
!***********************************************************************
!     
!     LBsoft subroutine for writing a scalar field of INTEGER
!     with 1 byte in binary format in serial IO
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: myio
  integer, intent(in), dimension(3) :: mydim
  integer(kind=1), allocatable, dimension(:,:,:)  :: myvar
  
  integer :: i,j,k
  
  do k=1,mydim(3)
    do j=1,mydim(2)
      do i=1,mydim(1)
        write(myio)myvar(i,j,k)
      enddo
    enddo
  enddo
  
  return
 
 end subroutine write_binary_1d_bynary_int
  
 subroutine driving_print_binary_1d_vtk(iotest,headoff,nbyte,nn,myvar,e_io)
 
!***********************************************************************
!     
!     LBsoft subroutine for writing a scalar field with single
!     precision in VTK legacy binary format in serial IO
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: iotest,headoff,nbyte,nn
  real(4), intent(in), dimension(:,:,:), allocatable :: myvar
  integer, intent(inout) :: e_io
  
  call print_binary_1d_vtk_serial(iotest,headoff,nbyte,nn, &
    reshape(myvar,(/nn/)),e_io)
  
  return
  
 end subroutine driving_print_binary_1d_vtk
 
 subroutine print_binary_1d_vtk_serial(iotest,headoff,nbyte,nn,myvar1d,E_IO)
 
!***********************************************************************
!     
!     LBsoft subroutine for writing a scalar field with single
!     precision in VTK legacy binary format in serial IO
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: iotest,headoff,nbyte,nn
  real(4), intent(in), dimension(nn) :: myvar1d
  integer, intent(out) :: E_IO
  
  integer :: ii
  
  write(iotest,iostat=E_IO)int(nbyte,kind=4),(myvar1d(ii),ii=1,nn)
  
  return
  
 end subroutine print_binary_1d_vtk_serial
 
 subroutine driving_print_binary_3d_vtk(iotest,headoff,nbyte,nn,myvar1, &
  myvar2,myvar3,e_io)
  
!***********************************************************************
!     
!     LBsoft subroutine for writing a vector field with single
!     precision in VTK legacy binary format in serial IO
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: iotest,headoff,nbyte,nn
  real(4), intent(in), dimension(:,:,:), allocatable :: myvar1,myvar2,myvar3
  integer, intent(inout) :: e_io
  
  call print_binary_3d_vtk_serial(iotest,headoff,nbyte,nn, &
   reshape(myvar1,(/nn/)),reshape(myvar2,(/nn/)), &
   reshape(myvar3,(/nn/)),e_io)
  
  return
  
 end subroutine driving_print_binary_3d_vtk
 
 subroutine print_binary_3d_vtk_serial(iotest,headoff,nbyte,nn,myvarx, &
  myvary,myvarz,E_IO)
  
!***********************************************************************
!     
!     LBsoft subroutine for writing a vector field with single
!     precision in VTK legacy binary format in serial IO
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: iotest,headoff,nbyte,nn
  real(4), intent(in), dimension(nn) :: myvarx,myvary,myvarz
  integer, intent(out) :: E_IO
  
  integer :: ii
  
  write(iotest,iostat=E_IO)int(nbyte,kind=4),(myvarx(ii),myvary(ii), &
   myvarz(ii),ii=1,nn)
  
  return
  
 end subroutine print_binary_3d_vtk_serial
 
 end module version_mod
