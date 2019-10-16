
#include <default_macro.h>
#define CHECKTYPE 0
#define ADDSYNC 0
 module version_mod
 
!***********************************************************************
!     
!     LBsoft module containing subroutines which manage all the 
!     communication tasks for the code parallel version and assign
!     a subset of beads to a specific node
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 use, intrinsic ::  iso_c_binding
 
 implicit none
 
 include 'mpif.h'
 
 private
 
 integer, public, save :: idrank=0
 integer, public, save :: mxrank=1
 
#if PRC==4
 integer, parameter :: MY_FLOAT=MPI_REAL4
#elif PRC==8
 integer, parameter :: MY_FLOAT=MPI_REAL8
#endif

#if LPRC==4
 integer, parameter :: MY_LOGICAL=MPI_LOGICAL
#endif

#if IPRC==4
 integer, parameter :: MY_INTEGER=MPI_INTEGER4
#endif
 
 integer, parameter :: MY_CHAR=MPI_CHARACTER
 
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
 public :: or_world_larr, or1_world_larr
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
 
 integer, allocatable, dimension(:), save :: ibuffer
 real(kind=PRC), allocatable, dimension(:), save :: fbuffer
 logical, allocatable, dimension(:), save :: lbuffer
 real(kind=c_long_double), allocatable, dimension(:), save :: qbuffer

 integer, save :: nibuffer,nfbuffer,ndbuffer,nlbuffer,nqbuffer
 logical, save :: aibuffer=.false.
 logical, save :: afbuffer=.false.
 logical, save :: adbuffer=.false.
 logical, save :: albuffer=.false.
 logical, save :: aqbuffer=.false.

 logical, allocatable, dimension(:), save :: lbuffer1
 logical, save :: albuffer1=.false.
 integer, save :: nlbuffer1

 integer, parameter :: nincrement=100
 
 contains
 
 subroutine allocate_ibuffer(narr)
 
!***********************************************************************
!     
!     LBsoft subroutine for allocating the buffer array of integer
!     type
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: narr
  
  if(aibuffer)then
    if(narr>nibuffer)then
      deallocate(ibuffer)
      nibuffer=narr+nincrement
      allocate(ibuffer(nibuffer))
    endif
  else
    nibuffer=narr+nincrement
    allocate(ibuffer(nibuffer))
    aibuffer=.true.
  endif
  
  return
  
 end subroutine allocate_ibuffer
 
 subroutine allocate_fbuffer(narr)
 
!***********************************************************************
!     
!     LBsoft subroutine for allocating the buffer array of float
!     type
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: narr
  
  if(afbuffer)then
    if(narr>nfbuffer)then
      deallocate(fbuffer)
      nfbuffer=narr+nincrement
      allocate(fbuffer(nfbuffer))
    endif
  else
    nfbuffer=narr+nincrement
    allocate(fbuffer(nfbuffer))
    afbuffer=.true.
  endif
  
  return
  
 end subroutine allocate_fbuffer
 
 subroutine allocate_lbuffer(narr)
 
!***********************************************************************
!     
!     LBsoft subroutine for allocating the buffer array of logical
!     type
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: narr
  
  if(albuffer)then
    if(narr>nlbuffer)then
      deallocate(lbuffer)
      nlbuffer=narr+nincrement
      allocate(lbuffer(nlbuffer))
    endif
  else
    nlbuffer=narr+nincrement
    allocate(lbuffer(nlbuffer))
    albuffer=.true.
  endif
  
  return
  
 end subroutine allocate_lbuffer
 
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
  

  if(mxrank==1)then
    write(iu,'(a40,i4,a6)')'The code is running in parallel mode on ',mxrank, &
       ' CPUs '
  else
    write(iu,'(a40,i4,a6)')'The code is running in parallel mode on ',mxrank, &
       ' CPUs '
  endif
  
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
  
  integer ier

  call MPI_COMM_RANK(MPI_COMM_WORLD,idrank,ier)
  
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
  
  
  integer ier
  
  call MPI_COMM_SIZE(MPI_COMM_WORLD,mxrank,ier)
  
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
  
  integer ier
  
  call MPI_BARRIER(MPI_COMM_WORLD,ier)
  
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
  
  integer ier
  
  call MPI_INIT(ier)
  
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
  
  integer ier
  
  call MPI_FINALIZE(ier)
  
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
 
  call MPI_ABORT(MPI_COMM_WORLD,30,ier)
  
  return
  
 end subroutine abort_world
 
 subroutine time_world(timecpu)
 
!***********************************************************************
!     
!     LBsoft subroutine to use the MPI CPU clock
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification January 2017
!     
!***********************************************************************
 
  implicit none
  
  
  
  real(kind=PRC), intent(out) :: timecpu
  
  if(idrank==0)call cpu_time(timecpu)
  call bcast_world_f(timecpu)
  
  return
  
 end subroutine time_world
 
 subroutine bcast_world_i(argument)
 
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
  
  integer, intent(inout) :: argument
  
  integer ier
  
#if CHECKTYPE==1
  if(kind(argument).ne. IPRC )then
    write(6,'(2i8,a,i8)')idrank,kind(argument), &
     ' in bcast_world_i expected',IPRC
    call flush(6)
  endif
#endif
  
  call MPI_BCAST(argument,1,MY_INTEGER,0,MPI_COMM_WORLD,ier)
  
#if ADDSYNC==1
  call get_sync_world
#endif
  
  return
  
 end subroutine bcast_world_i
 
 subroutine bcast_world_l(argument)
 
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
  
  logical, intent(inout) :: argument
    
  integer ier
  
#if CHECKTYPE==1
  if(kind(argument).ne. LPRC )then
    write(6,'(2i8,a,i8)')idrank,kind(argument), &
     ' in bcast_world_l expected',LPRC
    call flush(6)
  endif
#endif
  
  call MPI_BCAST(argument,1,MY_LOGICAL,0,MPI_COMM_WORLD,ier)
  
#if ADDSYNC==1
  call get_sync_world
#endif
  
  return
  
 end subroutine bcast_world_l
 
 subroutine bcast_world_f(argument)
 
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
  
  real(kind=PRC), intent(inout) :: argument
  
  integer ier
  
#if CHECKTYPE==1
  if(kind(argument).ne. PRC )then
    write(6,'(2i8,a,i8)')idrank,kind(argument), &
     ' in bcast_world_f expected',PRC
    call flush(6)
  endif
#endif
  
  call MPI_BCAST(argument,1,MY_FLOAT,0,MPI_COMM_WORLD,ier)
  
#if ADDSYNC==1
  call get_sync_world
#endif
  
  return
  
 end subroutine bcast_world_f
 
 subroutine bcast_world_iarr(argument,narr)
 
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
  
  integer, intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  
  integer ier
  
  call MPI_BCAST(argument,narr,MY_INTEGER,0,MPI_COMM_WORLD,ier)
  
#if ADDSYNC==1
  call get_sync_world
#endif
  
  return
  
 end subroutine bcast_world_iarr
 
 subroutine bcast_world_larr(argument,narr)
 
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
  
  
  
  logical, intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
    
  integer ier
  
  call MPI_BCAST(argument,narr,MY_LOGICAL,0,MPI_COMM_WORLD,ier)
  
#if ADDSYNC==1
  call get_sync_world
#endif
  
  return
  
 end subroutine bcast_world_larr
 
 subroutine bcast_world_farr(argument,narr)
  
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
  
  real(kind=PRC), intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  
  integer ier
  
  call MPI_BCAST(argument,narr,MY_FLOAT,0,MPI_COMM_WORLD,ier)
  
#if ADDSYNC==1
  call get_sync_world
#endif
  
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
  character(len=mxlen), intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  
  integer ier
  
  call MPI_BCAST(argument,mxlen*narr,MY_CHAR,0,MPI_COMM_WORLD,ier)
  
#if ADDSYNC==1
  call get_sync_world
#endif
  
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
  
  if(present(buffersub))then
    
    call MPI_ALLREDUCE(argument,buffersub,narr,MY_INTEGER, &
      MPI_SUM,MPI_COMM_WORLD,ier)
    
  else
    
    call allocate_ibuffer(narr)
    
    call MPI_ALLREDUCE(argument,ibuffer,narr,MY_INTEGER, &
      MPI_SUM,MPI_COMM_WORLD,ier)
    
    argument(1:narr)=ibuffer(1:narr)
    
  endif
  
#if ADDSYNC==1
  call get_sync_world
#endif
  
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
  
  if(present(buffersub))then
    
    call MPI_ALLREDUCE(argument,buffersub,narr,MY_FLOAT, &
      MPI_SUM,MPI_COMM_WORLD,ier)
    
  else
    
    call allocate_fbuffer(narr)
    
    call MPI_ALLREDUCE(argument,fbuffer,narr,MY_FLOAT, &
     MPI_SUM,MPI_COMM_WORLD,ier)
    
    argument(1:narr)=fbuffer(1:narr)
    
  endif
  
#if ADDSYNC==1
  call get_sync_world
#endif
  
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
  
  if(present(buffersub))then
    
    call MPI_ALLREDUCE(argument,buffersub,narr,MY_INTEGER, &
      MPI_MIN,MPI_COMM_WORLD,ier)
    
  else
    
    call allocate_ibuffer(narr)
    
    call MPI_ALLREDUCE(argument,ibuffer,narr,MY_INTEGER, &
      MPI_MIN,MPI_COMM_WORLD,ier)
      
    argument(1:narr)=ibuffer(1:narr)
    
  endif
  
#if ADDSYNC==1
  call get_sync_world
#endif
  
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
  
  if(present(buffersub))then
    
    call MPI_ALLREDUCE(argument,buffersub,narr,MY_FLOAT, &
      MPI_MIN,MPI_COMM_WORLD,ier)
    
  else
    
    call allocate_fbuffer(narr)
    
    call MPI_ALLREDUCE(argument,fbuffer,narr,MY_FLOAT, &
      MPI_MIN,MPI_COMM_WORLD,ier)
    
    argument(1:narr)=fbuffer(1:narr)
    
  endif
  
#if ADDSYNC==1
  call get_sync_world
#endif
  
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
  
  if(present(buffersub))then
    
    call MPI_ALLREDUCE(argument,buffersub,narr,MY_INTEGER, &
      MPI_MAX,MPI_COMM_WORLD,ier)
    
  else
    
    call allocate_ibuffer(narr)
    
    call MPI_ALLREDUCE(argument,ibuffer,narr,MY_INTEGER, &
      MPI_MAX,MPI_COMM_WORLD,ier)
      
    argument(1:narr)=ibuffer(1:narr)
    
  endif
  
#if ADDSYNC==1
  call get_sync_world
#endif
  
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
  
  if(present(buffersub))then
    
    call MPI_ALLREDUCE(argument,buffersub,narr,MY_FLOAT, &
      MPI_MAX,MPI_COMM_WORLD,ier)
    
  else
    
    call allocate_fbuffer(narr)
    
    call MPI_ALLREDUCE(argument,fbuffer,narr,MY_FLOAT, &
      MPI_MAX,MPI_COMM_WORLD,ier)
    
    argument(1:narr)=fbuffer(1:narr)
    
  endif
  
#if ADDSYNC==1
  call get_sync_world
#endif
  
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
  
  if(present(buffersub))then
    
    call MPI_ALLREDUCE(argument,buffersub,narr,MY_LOGICAL, &
      MPI_LAND,MPI_COMM_WORLD,ier)
    
  else
    
    call allocate_lbuffer(narr)
    
    call MPI_ALLREDUCE(argument,lbuffer,narr,MY_LOGICAL, &
      MPI_LAND,MPI_COMM_WORLD,ier)
      
    argument(1:narr)=lbuffer(1:narr)
    
  endif
  
#if ADDSYNC==1
  call get_sync_world
#endif
  
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
  
  if(present(buffersub))then
    
    call MPI_ALLREDUCE(argument,buffersub,narr,MY_LOGICAL, &
      MPI_LOR,MPI_COMM_WORLD,ier)
    
  else
    
    call allocate_lbuffer(narr)
    
    call MPI_ALLREDUCE(argument,lbuffer,narr,MY_LOGICAL, &
      MPI_LOR,MPI_COMM_WORLD,ier)
      
    argument(1:narr)=lbuffer(1:narr)
    
  endif
  
#if ADDSYNC==1
  call get_sync_world
#endif
  
  return
  
 end subroutine or_world_larr
 
 subroutine allocate_lbuffer1(narr)
  implicit none
  integer, intent(in) :: narr


  if(albuffer1)then
    if(narr>nlbuffer1)then
      deallocate(lbuffer1)
      nlbuffer1=narr+nincrement
      allocate(lbuffer1(nlbuffer1))
    endif
  else
    nlbuffer1=narr+nincrement
    allocate(lbuffer1(nlbuffer1))
    albuffer1=.true.
  endif
 end subroutine allocate_lbuffer1

 subroutine or1_world_larr(argument,narr,buffersub)
  implicit none
  logical(kind=1), intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  logical, intent(inout), dimension(narr), optional :: buffersub
  logical, allocatable, dimension(:) :: temp

  integer ier

  if(present(buffersub))then
    call MPI_ALLREDUCE(argument,buffersub,narr,MY_LOGICAL, &
      MPI_LOR,MPI_COMM_WORLD,ier)
      argument(1:narr) = buffersub(1:narr)
  else
    call allocate_lbuffer(narr)
    allocate(temp(narr))
    temp(1:narr) = argument(1:narr)

    call MPI_ALLREDUCE(temp,lbuffer,narr,MY_LOGICAL, &
      MPI_LOR,MPI_COMM_WORLD,ier)

    argument(1:narr) = lbuffer(1:narr)
    deallocate(temp)
  endif

#if ADDSYNC==1
  call get_sync_world
#endif
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
  
  if(mxrank==1)return
  
  CALL MPI_IRECV(argument,1,MY_INTEGER, &
   idsource,tag,MPI_COMM_WORLD,irequest,ierr)
  
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
  
  if(mxrank==1)return
  
  CALL MPI_IRECV(argument,1,MY_FLOAT, &
   idsource,tag,MPI_COMM_WORLD,irequest,ierr)
  
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
  
  if(mxrank==1)return
  
  CALL MPI_IRECV(argument,1,MY_LOGICAL, &
   idsource,tag,MPI_COMM_WORLD,irequest,ierr)
  
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
  
  if(mxrank==1)return
  
  CALL MPI_ISEND(argument,1,MY_INTEGER, &
   iddest,tag,MPI_COMM_WORLD,irequest,ierr)
  
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
  
  if(mxrank==1)return
  
  CALL MPI_ISEND(argument,1,MY_FLOAT, &
   iddest,tag,MPI_COMM_WORLD,irequest,ierr)
  
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
  
  if(mxrank==1)return
  
  CALL MPI_ISEND(argument,1,MY_LOGICAL, &
   iddest,tag,MPI_COMM_WORLD,irequest,ierr)
  
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
  
  if(mxrank==1)return
  
  CALL MPI_IRECV(argument,narr,MY_INTEGER, &
   idsource,tag,MPI_COMM_WORLD,irequest,ierr)
  
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
  
  if(mxrank==1)return
  
  CALL MPI_IRECV(argument,narr,MY_FLOAT, &
   idsource,tag,MPI_COMM_WORLD,irequest,ierr)
  
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
  
  if(mxrank==1)return
  
  CALL MPI_IRECV(argument,narr,MY_LOGICAL, &
   idsource,tag,MPI_COMM_WORLD,irequest,ierr)
  
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
  
  if(mxrank==1)return
  
  CALL MPI_ISEND(argument,narr,MY_INTEGER, &
   iddest,tag,MPI_COMM_WORLD,irequest,ierr)
  
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
  
  if(mxrank==1)return
  
  CALL MPI_ISEND(argument,narr,MY_FLOAT, &
   iddest,tag,MPI_COMM_WORLD,irequest,ierr)
  
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
  
  if(mxrank==1)return
  
  CALL MPI_ISEND(argument,narr,MY_LOGICAL, &
   iddest,tag,MPI_COMM_WORLD,irequest,ierr)
  
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
  integer, dimension(MPI_STATUS_SIZE) :: istatus
  integer :: ierr
  
  if(mxrank==1)return
   
  CALL MPI_WAIT(argument,istatus,ierr)
  
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
  integer, dimension(MPI_STATUS_SIZE,narr) :: istatus
  integer :: ierr
  
  if(mxrank==1)return
   
  CALL MPI_WAITALL(narr,argument,istatus,ierr)
  
  return
  
 end subroutine waitall_world

 subroutine allocate_qbuffer(narr)
  implicit none
  integer, intent(in) :: narr

  if(aqbuffer)then
    if(narr>nqbuffer)then
      deallocate(qbuffer)
      nqbuffer=narr+nincrement
      allocate(qbuffer(nqbuffer))
    endif
  else
    nqbuffer=narr+nincrement
    allocate(qbuffer(nqbuffer))
    aqbuffer=.true.
  endif
 end subroutine allocate_qbuffer

 subroutine sum_world_qarr(argument,narr,buffersub)
  
!***********************************************************************
!     
!     LBsoft global summation subroutine for a float array
!     in quadruple-precision
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC*2), intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  real(kind=PRC*2), intent(inout), dimension(narr), optional :: buffersub
  integer ier
  real(kind=c_long_double), dimension(narr) :: converted


  converted(1:narr) = argument(1:narr)

  if(present(buffersub))then
    call MPI_ALLREDUCE(argument,buffersub,narr,MPI_REAL16, &
       MPI_SUM,MPI_COMM_WORLD,ier)
   else
    call allocate_qbuffer(narr)
    call MPI_ALLREDUCE(converted,qbuffer,narr,MPI_REAL16, &
     MPI_SUM,MPI_COMM_WORLD,ier)
    argument(1:narr)=qbuffer(1:narr)
  endif
  
  return

 end subroutine sum_world_qarr
 
 subroutine open_file_vtk_par(iotest,nn,myname,e_io)
 
!***********************************************************************
!     
!     LBsoft subroutine for opening the vtk legacy file
!     in parallel IO
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
  
  call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(myname), &
                       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                       MPI_INFO_NULL,iotest,e_io)
  return
  
 endsubroutine open_file_vtk_par
 
 subroutine close_file_vtk_par(iotest,e_io)
 
!***********************************************************************
!     
!     LBsoft subroutine for closing the vtk legacy file
!     in parallel IO
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: iotest
  integer, intent(out) :: e_io
  
  call MPI_FILE_CLOSE(iotest, e_io)
  
  return
  
 endsubroutine close_file_vtk_par
 
 subroutine print_header_vtk_par(iotest,offsetsub,nn,header,e_io)
 
!***********************************************************************
!     
!     LBsoft subroutine for writing the header part of
!     in VTK legacy file in parallel IO
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: iotest,offsetsub,nn
  character(len=500) :: header
  integer, intent(out) :: e_io
  
  integer(kind=MPI_OFFSET_KIND) :: myoffset
  
  myoffset=int(offsetsub,kind=MPI_OFFSET_KIND)
  
  if(idrank==0)call MPI_File_write_at(iotest,myoffset,header(1:nn),nn, &
   MPI_CHARACTER,MPI_STATUS_IGNORE,e_io)
  
  return
  
 endsubroutine print_header_vtk_par
 
 subroutine print_footer_vtk_par(iotest,offsetsub,footer,e_io)
 
!***********************************************************************
!     
!     LBsoft subroutine for writing the footer part of
!     in VTK legacy file in parallel IO
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: iotest,offsetsub
  character(len=30) :: footer
  integer, intent(out) :: e_io
  
  integer(kind=MPI_OFFSET_KIND) :: myoffset
  
  myoffset=int(offsetsub,kind=MPI_OFFSET_KIND)
  
  if(idrank==0)call MPI_File_write_at(iotest,myoffset,footer(1:30),30, &
   MPI_CHARACTER,MPI_STATUS_IGNORE,e_io)
  
  return
  
 endsubroutine print_footer_vtk_par
  
 end module version_mod
