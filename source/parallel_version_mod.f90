
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
!     last modification January 2017
!     
!***********************************************************************
 
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
 public :: sum_world_iarr
 public :: sum_world_farr
 public :: min_world_iarr
 public :: min_world_farr
 public :: max_world_iarr
 public :: max_world_farr
 public :: and_world_larr
 public :: or_world_larr
 
 integer, allocatable, dimension(:), save :: ibuffer
 real(kind=PRC), allocatable, dimension(:), save :: fbuffer
 logical, allocatable, dimension(:), save :: lbuffer
 integer, save :: nibuffer,nfbuffer,ndbuffer,nlbuffer
 logical, save :: aibuffer=.false.
 logical, save :: afbuffer=.false.
 logical, save :: adbuffer=.false.
 logical, save :: albuffer=.false.
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
 
 subroutine allocate_dbuffer(narr)
 
!***********************************************************************
!     
!     LBsoft subroutine for allocating the buffer array of double
!     precision type
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: narr
  
  if(adbuffer)then
    if(narr>ndbuffer)then
      deallocate(dbuffer)
      ndbuffer=narr+nincrement
      allocate(dbuffer(ndbuffer))
    endif
  else
    ndbuffer=narr+nincrement
    allocate(dbuffer(ndbuffer))
    adbuffer=.true.
  endif
  
  return
  
 end subroutine allocate_dbuffer
 
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
  
 end module version_mod
