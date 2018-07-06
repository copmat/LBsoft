
#include <default_macro.h>
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
 
 private
 
 integer, public, save :: idrank=0
 integer, public, save :: mxrank=1
 
 interface bcast_world
   module procedure bcast_world_i
   module procedure bcast_world_l
   module procedure bcast_world_f
   module procedure bcast_world_d
   module procedure bcast_world_iarr
   module procedure bcast_world_larr
   module procedure bcast_world_farr
   module procedure bcast_world_darr
 end interface bcast_world
 
 interface sum_world
   module procedure sum_world_iarr
   module procedure sum_world_farr
   module procedure sum_world_darr
 end interface sum_world
 
 interface min_world
   module procedure min_world_iarr
   module procedure min_world_farr
   module procedure min_world_darr
 end interface min_world
 
 interface max_world
   module procedure max_world_iarr
   module procedure max_world_farr
   module procedure max_world_darr
 end interface max_world
 
 public :: print_version
 public :: get_rank_world
 public :: get_size_world
 public :: get_sync_world
 public :: init_world
 public :: finalize_world
 public :: abort_world
 public :: time_world
 public :: bcast_world
 public :: sum_world
 public :: min_world
 public :: max_world
 public :: and_world_larr
 public :: or_world_larr
 
 integer, allocatable, dimension(:), save :: ibuffer
 real(4), allocatable, dimension(:), save :: fbuffer
 real(8), allocatable, dimension(:), save :: dbuffer
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
!     LBsoft subroutine for allocating the buffer array of single
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
  
  include 'mpif.h'
  
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
  
  include 'mpif.h'
  
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
  
  include 'mpif.h'
  
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
  
  include 'mpif.h'
  
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
  
  include 'mpif.h'
  
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
  
  include 'mpif.h'
  
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
  
  include 'mpif.h'
  
  real(kind=PRC), intent(out) :: timecpu
  
  if(idrank==0)call cpu_time(timecpu)
  call bcast_world(timecpu)
  
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
  
  include 'mpif.h'
  
  integer, intent(inout) :: argument
  
  integer ier
  
  call MPI_BCAST(argument,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  
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
  
  include 'mpif.h'
  
  logical, intent(inout) :: argument
  
  integer ier
  
  call MPI_BCAST(argument,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)
  
  return
  
 end subroutine bcast_world_l
 
 subroutine bcast_world_f(argument)
 
!***********************************************************************
!     
!     LBsoft subroutine to broadcast an single precision number to all 
!     other nodes
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  include 'mpif.h'
  
  real(4), intent(inout) :: argument
  
  integer ier
  
  call MPI_BCAST(argument,1,MPI_REAL4,0,MPI_COMM_WORLD,ier)
  
  return
  
 end subroutine bcast_world_f
 
 subroutine bcast_world_d(argument)
 
!***********************************************************************
!     
!     LBsoft subroutine to broadcast an double precision number to all 
!     other nodes
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  include 'mpif.h'
  
  real(8), intent(inout) :: argument
  
  integer ier
  
  call MPI_BCAST(argument,1,MPI_REAL8,0,MPI_COMM_WORLD,ier)
  
  return
  
 end subroutine bcast_world_d
 
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
  
  include 'mpif.h'
  
  integer, intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  
  integer ier
  
  call MPI_BCAST(argument,narr,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  
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
  
  include 'mpif.h'
  
  logical, intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  
  integer ier
  
  call MPI_BCAST(argument,narr,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)
  
  return
  
 end subroutine bcast_world_larr
 
 subroutine bcast_world_farr(argument,narr)
  
!***********************************************************************
!     
!     LBsoft subroutine to broadcast an single precision array to all 
!     other nodes
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  include 'mpif.h'
  
  real(4), intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  
  integer ier
  
  call MPI_BCAST(argument,narr,MPI_REAL4,0,MPI_COMM_WORLD,ier)
  
  return
  
 end subroutine bcast_world_farr
 
 subroutine bcast_world_darr(argument,narr)
  
!***********************************************************************
!     
!     LBsoft subroutine to broadcast an double precision array to all 
!     other nodes
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  include 'mpif.h'
  
  real(8), intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  
  integer ier
  
  call MPI_BCAST(argument,narr,MPI_REAL8,0,MPI_COMM_WORLD,ier)
  
  return
  
 end subroutine bcast_world_darr
 
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
  
  include 'mpif.h'
  
  integer, intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  integer, intent(inout), dimension(narr), optional :: buffersub
  
  integer ier
  
  if(present(buffersub))then
    
    call MPI_ALLREDUCE(argument,buffersub,narr,MPI_INTEGER, &
      MPI_SUM,MPI_COMM_WORLD,ier)
    
  else
    
    call allocate_ibuffer(narr)
    
    call MPI_ALLREDUCE(argument,ibuffer,narr,MPI_INTEGER, &
      MPI_SUM,MPI_COMM_WORLD,ier)
    
    argument(1:narr)=ibuffer(1:narr)
    
  endif
  
  return
  
 end subroutine sum_world_iarr
 
 subroutine sum_world_farr(argument,narr,buffersub)
 
!***********************************************************************
!     
!     LBsoft global summation subroutine for a double precision array
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  include 'mpif.h'
  
  real(4), intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  real(4), intent(inout), dimension(narr), optional :: buffersub
  
  integer ier
  
  if(present(buffersub))then
    
    call MPI_ALLREDUCE(argument,buffersub,narr,MPI_REAL4, &
      MPI_SUM,MPI_COMM_WORLD,ier)
    
  else
    
    call allocate_fbuffer(narr)
    
    call MPI_ALLREDUCE(argument,fbuffer,narr,MPI_REAL4, &
      MPI_SUM,MPI_COMM_WORLD,ier)
    
    argument(1:narr)=fbuffer(1:narr)
    
  endif
  
  return
  
 end subroutine sum_world_farr
 
 subroutine sum_world_darr(argument,narr,buffersub)
 
!***********************************************************************
!     
!     LBsoft global summation subroutine for a double precision array
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  include 'mpif.h'
  
  real(8), intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  real(8), intent(inout), dimension(narr), optional :: buffersub
  
  integer ier
  
  if(present(buffersub))then
    
    call MPI_ALLREDUCE(argument,buffersub,narr,MPI_REAL8, &
      MPI_SUM,MPI_COMM_WORLD,ier)
    
  else
    
    call allocate_dbuffer(narr)
    
    call MPI_ALLREDUCE(argument,dbuffer,narr,MPI_REAL8, &
      MPI_SUM,MPI_COMM_WORLD,ier)
    
    argument(1:narr)=dbuffer(1:narr)
    
  endif
  
  return
  
 end subroutine sum_world_darr
 
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
  
  include 'mpif.h'
  
  integer, intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  integer, intent(inout), dimension(narr), optional :: buffersub
  
  integer ier
  
  if(present(buffersub))then
    
    call MPI_ALLREDUCE(argument,buffersub,narr,MPI_INTEGER, &
      MPI_MIN,MPI_COMM_WORLD,ier)
    
  else
    
    call allocate_ibuffer(narr)
    
    call MPI_ALLREDUCE(argument,ibuffer,narr,MPI_INTEGER, &
      MPI_MIN,MPI_COMM_WORLD,ier)
      
    argument(1:narr)=ibuffer(1:narr)
    
  endif
  
  return
  
 end subroutine min_world_iarr
 
 subroutine min_world_farr(argument,narr,buffersub)
 
!***********************************************************************
!     
!     LBsoft global minimum subroutine for a double precision array
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  include 'mpif.h'
  
  real(4), intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  real(4), intent(inout), dimension(narr), optional :: buffersub
  
  integer ier
  
  if(present(buffersub))then
    
    call MPI_ALLREDUCE(argument,buffersub,narr,MPI_REAL4, &
      MPI_MIN,MPI_COMM_WORLD,ier)
    
  else
    
    call allocate_fbuffer(narr)
    
    call MPI_ALLREDUCE(argument,fbuffer,narr,MPI_REAL4, &
      MPI_MIN,MPI_COMM_WORLD,ier)
    
    argument(1:narr)=fbuffer(1:narr)
    
  endif
  
  return
  
 end subroutine min_world_farr
 
 subroutine min_world_darr(argument,narr,buffersub)
 
!***********************************************************************
!     
!     LBsoft global minimum subroutine for a double precision array
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  include 'mpif.h'
  
  real(8), intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  real(8), intent(inout), dimension(narr), optional :: buffersub
  
  integer ier
  
  if(present(buffersub))then
  
    call MPI_ALLREDUCE(argument,buffersub,narr,MPI_REAL8, &
      MPI_MIN,MPI_COMM_WORLD,ier)
  
  else
  
    call allocate_dbuffer(narr)
  
    call MPI_ALLREDUCE(argument,dbuffer,narr,MPI_REAL8, &
      MPI_MIN,MPI_COMM_WORLD,ier)
    
    argument(1:narr)=dbuffer(1:narr)
  
  endif
  
  return
  
 end subroutine min_world_darr
 
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
  
  include 'mpif.h'
  
  integer, intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  integer, intent(inout), dimension(narr), optional :: buffersub
  
  integer ier
  
  if(present(buffersub))then
    
    call MPI_ALLREDUCE(argument,buffersub,narr,MPI_INTEGER, &
      MPI_MAX,MPI_COMM_WORLD,ier)
    
  else
    
    call allocate_ibuffer(narr)
    
    call MPI_ALLREDUCE(argument,ibuffer,narr,MPI_INTEGER, &
      MPI_MAX,MPI_COMM_WORLD,ier)
      
    argument(1:narr)=ibuffer(1:narr)
    
  endif
  
  return
  
 end subroutine max_world_iarr
 
 subroutine max_world_farr(argument,narr,buffersub)
 
!***********************************************************************
!     
!     LBsoft global maximum subroutine for a double precision array
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  include 'mpif.h'
  
  real(4), intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  real(4), intent(inout), dimension(narr), optional :: buffersub
  
  integer ier
  
  if(present(buffersub))then
    
    call MPI_ALLREDUCE(argument,buffersub,narr,MPI_REAL4, &
      MPI_MAX,MPI_COMM_WORLD,ier)
    
  else
    
    call allocate_fbuffer(narr)
    
    call MPI_ALLREDUCE(argument,fbuffer,narr,MPI_REAL4, &
      MPI_MAX,MPI_COMM_WORLD,ier)
    
    argument(1:narr)=fbuffer(1:narr)
    
  endif
  
  return
  
 end subroutine max_world_farr
 
 subroutine max_world_darr(argument,narr,buffersub)
 
!***********************************************************************
!     
!     LBsoft global maximum subroutine for a double precision array
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  include 'mpif.h'
  
  real(8), intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  real(8), intent(inout), dimension(narr), optional :: buffersub
  
  integer ier
  
  if(present(buffersub))then
  
    call MPI_ALLREDUCE(argument,buffersub,narr,MPI_REAL8, &
      MPI_MAX,MPI_COMM_WORLD,ier)
  
  else
  
    call allocate_dbuffer(narr)
  
    call MPI_ALLREDUCE(argument,dbuffer,narr,MPI_REAL8, &
      MPI_MAX,MPI_COMM_WORLD,ier)
    
    argument(1:narr)=dbuffer(1:narr)
  
  endif
  
  return
  
 end subroutine max_world_darr
 
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
  
  include 'mpif.h'
  
  logical, intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  logical, intent(inout), dimension(narr), optional :: buffersub
  
  integer ier
  
  if(present(buffersub))then
    
    call MPI_ALLREDUCE(argument,buffersub,narr,MPI_LOGICAL, &
      MPI_LAND,MPI_COMM_WORLD,ier)
    
  else
    
    call allocate_lbuffer(narr)
    
    call MPI_ALLREDUCE(argument,lbuffer,narr,MPI_LOGICAL, &
      MPI_LAND,MPI_COMM_WORLD,ier)
      
    argument(1:narr)=lbuffer(1:narr)
    
  endif
  
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
  
  include 'mpif.h'
  
  logical, intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  logical, intent(inout), dimension(narr), optional :: buffersub
  
  integer ier
  
  if(present(buffersub))then
    
    call MPI_ALLREDUCE(argument,buffersub,narr,MPI_LOGICAL, &
      MPI_LOR,MPI_COMM_WORLD,ier)
    
  else
    
    call allocate_lbuffer(narr)
    
    call MPI_ALLREDUCE(argument,lbuffer,narr,MPI_LOGICAL, &
      MPI_LOR,MPI_COMM_WORLD,ier)
      
    argument(1:narr)=lbuffer(1:narr)
    
  endif
  
  return
  
 end subroutine or_world_larr
  
 end module version_mod
