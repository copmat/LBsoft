
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
 
 implicit none
 
 private
 
 integer, public, save :: idrank=0
 integer, public, save :: mxrank=1
 
 integer, save :: myid,numprocs,nxy2,nx2
 
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
 
 end module version_mod
