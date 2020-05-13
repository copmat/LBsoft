
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
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 use, intrinsic ::  iso_c_binding
 use fluids_lattices_mod
 
 implicit none
 
 include 'mpif.h'
 
 private
 
 integer, public, save :: idrank=0
 integer, public, save :: mxrank=1
 
#if PRC==4
 integer, parameter :: MYFLOAT=MPI_REAL4
#elif PRC==8
 integer, parameter :: MYFLOAT=MPI_REAL8
#endif

#if LPRC==4
 integer, parameter :: MYLOGICAL=MPI_LOGICAL
#endif

#if IPRC==4
 integer, parameter :: MYINTEGER=MPI_INTEGER4
#endif
 
 integer, parameter :: MYCHAR=MPI_CHARACTER
 
 integer, parameter :: MYINT=MPI_INTEGER
 integer, parameter :: MYINT1=MPI_INTEGER1
 integer, parameter :: MYINT2=MPI_INTEGER2
 integer, parameter :: MYCOMM=MPI_COMM_WORLD
 
 integer, parameter :: maxneigh=32
 integer, parameter :: npop=19
 integer, parameter :: movetag=393000, &
  findntag=524000,findptag=655000,findwtag=786000
  
 integer, parameter :: findntag_hvar=findntag*2
 integer, parameter :: tag_hvar=movetag*2
 integer, parameter :: tag_hvar2=movetag*3
 integer, parameter :: tag_isfluid=movetag*4
 integer, parameter :: findntag_sing=findntag*2
 integer, parameter :: tag_sing=movetag*2
 
 !define LARGEINT 1073741824
 integer, parameter :: LARGEINT=1073741824
 
 
 integer :: rc, ierr
 integer, public, protected :: lnofsite, itpppsize
 
 integer, protected, public :: domdec
 integer, parameter :: syncsend=1
 integer, save :: requestm(0:maxneigh-1), requestp(0:maxneigh-1), &
  requestw(0:maxneigh-1), requesti(0:maxneigh-1), requesto(0:maxneigh-1), &
  requesth(0:maxneigh-1), requesth_wall(0:maxneigh-1), requests(0:maxneigh-1), &
  request_hvar(0:maxneigh-1),requestm_hvar(0:maxneigh-1), &
  request_hvar2(0:maxneigh-1),requestm_hvar2(0:maxneigh-1), &
  request_isfluid(0:maxneigh-1),requestm_isfluid(0:maxneigh-1), &
  request_sing(0:maxneigh-1),requestm_sing(0:maxneigh-1)
  
 logical :: firstmove=.true.
  
 logical :: allpbc
 
 integer :: n_pe2recv_fluid, n_pe2send_fluid
 integer :: i_pe2send_fluid(0:maxneigh-1), n_pop2send_fluid(0:maxneigh-1)
 integer :: i_pe2recv_fluid(0:maxneigh-1), n_pop2recv_fluid(0:maxneigh-1)
  
 integer :: n_pe2recv_fluid_hvar, n_pe2send_fluid_hvar
 integer :: i_pe2send_fluid_hvar(0:maxneigh-1), n_var2send_fluid(0:maxneigh-1)
 integer :: i_pe2recv_fluid_hvar(0:maxneigh-1), n_var2recv_fluid(0:maxneigh-1)
 
 integer :: n_pe2recv_fluid_sing, n_pe2send_fluid_sing
 integer :: i_pe2send_fluid_sing(0:maxneigh-1), n_sing2send_fluid(0:maxneigh-1)
 integer :: i_pe2recv_fluid_sing(0:maxneigh-1), n_sing2recv_fluid(0:maxneigh-1)

 integer, protected, public :: nprocz, nprocy, nprocx
 integer, save :: nxy2, nx2,nx_comm,ny_comm,nz_comm
 integer :: nbuff_comm
 integer, save, dimension(3) :: ipbc

 integer, pointer :: countnpp(:)
 integer(kind=2), public, protected, allocatable :: ownern(:)

 integer(kind=IPRC), allocatable :: i_pop2send_fluid(:,:,:), i_pop2recv_fluid(:,:,:)
 integer(kind=IPRC), allocatable :: i_var2send_fluid(:,:,:), i_var2recv_fluid(:,:,:)
 integer(kind=IPRC), allocatable :: i_sing2send_fluid(:,:,:), i_sing2recv_fluid(:,:,:)
 integer(kind=IPRC), allocatable :: i_pop2send_wall(:,:,:), i_pop2recv_wall(:,:,:)
 integer(kind=IPRC), allocatable :: i_pop2send_inlet(:,:,:), i_pop2recv_inlet(:,:,:)
 integer(kind=IPRC), allocatable :: i_pop2send_outlet(:,:,:), i_pop2recv_outlet(:,:,:)
  
 real(kind=PRC), pointer :: buffpops(:,:), buffpopr(:,:)
 real(kind=PRC), pointer :: bufftags(:,:), bufftagr(:,:)
  
 real(kind=PRC), pointer :: buffs_hvar(:,:), buffr_hvar(:,:)
 real(kind=PRC), pointer :: buffs_hvar2(:,:), buffr_hvar2(:,:)
 
 real(kind=PRC), pointer :: buffs_sing(:,:), buffr_sing(:,:)
  
 integer(KIND=1), pointer :: buffs_isfluid(:,:), buffr_isfluid(:,:)
  
 logical, save :: ldo_second=.false.
  
 integer, dimension(:), allocatable, public, protected :: gminx,gmaxx, &
   gminy,gmaxy,gminz,gmaxz
  
  
 real(kind=PRC),allocatable,dimension(:,:),public :: buff_transf_r,buff_transf_s
  
  
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
 
 integer :: globalDims(3)
 integer :: ldims(3)
 integer :: mymin(3)
 integer :: mymax(3)
 
 public :: print_version
 public :: get_rank_world
 public :: get_size_world
 public :: get_sync_world
 public :: init_world
 public :: finalize_world
 public :: abort_world
 public :: time_world
 public :: wtime
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
 public :: driving_print_binary_1d_vtk
 public :: driving_print_binary_3d_vtk
 public :: collective_readFile_pops
 public :: collective_writeFile_pops
 
 public :: collective_readFile
 public :: collective_readFile_int
 public :: collective_writeFile
 public :: collective_writeFile_int
 
 public :: ownernfind
 public :: ownernfind_arr
 public :: commspop
 public :: commrpop
 
 contains
 
 subroutine allocate_ibuffer(narr)
 
!***********************************************************************
!     
!     LBsoft subroutine for allocating the buffer array of integer
!     type
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
 
 function wtime()

!***********************************************************************
!     
!     LBsoft subroutine for computing the wall-clock time
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     modified by: M. Lauricella
!     last modification January 2020
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC) wtime
  
  wtime =  MPI_Wtime()
  
  return
  
 end function wtime
 
 subroutine bcast_world_i(argument)
 
!***********************************************************************
!     
!     LBsoft subroutine to broadcast an integer number to all other 
!     nodes
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
  
  call MPI_BCAST(argument,1,MYINTEGER,0,MPI_COMM_WORLD,ier)
  
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
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
  
  call MPI_BCAST(argument,1,MYLOGICAL,0,MPI_COMM_WORLD,ier)
  
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
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
  
  call MPI_BCAST(argument,1,MYFLOAT,0,MPI_COMM_WORLD,ier)
  
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
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  
  integer ier
  
  call MPI_BCAST(argument,narr,MYINTEGER,0,MPI_COMM_WORLD,ier)
  
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
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  
  
  logical, intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
    
  integer ier
  
  call MPI_BCAST(argument,narr,MYLOGICAL,0,MPI_COMM_WORLD,ier)
  
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
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC), intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  
  integer ier
  
  call MPI_BCAST(argument,narr,MYFLOAT,0,MPI_COMM_WORLD,ier)
  
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
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: mxlen
  character(len=mxlen), intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  
  integer ier
  
  call MPI_BCAST(argument,mxlen*narr,MYCHAR,0,MPI_COMM_WORLD,ier)
  
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
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
    
    call MPI_ALLREDUCE(argument,buffersub,narr,MYINTEGER, &
      MPI_SUM,MPI_COMM_WORLD,ier)
    
  else
    
    call allocate_ibuffer(narr)
    
    call MPI_ALLREDUCE(argument,ibuffer,narr,MYINTEGER, &
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
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
    
    call MPI_ALLREDUCE(argument,buffersub,narr,MYFLOAT, &
      MPI_SUM,MPI_COMM_WORLD,ier)
    
  else
    
    call allocate_fbuffer(narr)
    
    call MPI_ALLREDUCE(argument,fbuffer,narr,MYFLOAT, &
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
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
    
    call MPI_ALLREDUCE(argument,buffersub,narr,MYINTEGER, &
      MPI_MIN,MPI_COMM_WORLD,ier)
    
  else
    
    call allocate_ibuffer(narr)
    
    call MPI_ALLREDUCE(argument,ibuffer,narr,MYINTEGER, &
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
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
    
    call MPI_ALLREDUCE(argument,buffersub,narr,MYFLOAT, &
      MPI_MIN,MPI_COMM_WORLD,ier)
    
  else
    
    call allocate_fbuffer(narr)
    
    call MPI_ALLREDUCE(argument,fbuffer,narr,MYFLOAT, &
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
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
    
    call MPI_ALLREDUCE(argument,buffersub,narr,MYINTEGER, &
      MPI_MAX,MPI_COMM_WORLD,ier)
    
  else
    
    call allocate_ibuffer(narr)
    
    call MPI_ALLREDUCE(argument,ibuffer,narr,MYINTEGER, &
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
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
    
    call MPI_ALLREDUCE(argument,buffersub,narr,MYFLOAT, &
      MPI_MAX,MPI_COMM_WORLD,ier)
    
  else
    
    call allocate_fbuffer(narr)
    
    call MPI_ALLREDUCE(argument,fbuffer,narr,MYFLOAT, &
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
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
    
    call MPI_ALLREDUCE(argument,buffersub,narr,MYLOGICAL, &
      MPI_LAND,MPI_COMM_WORLD,ier)
    
  else
    
    call allocate_lbuffer(narr)
    
    call MPI_ALLREDUCE(argument,lbuffer,narr,MYLOGICAL, &
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
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
    
    call MPI_ALLREDUCE(argument,buffersub,narr,MYLOGICAL, &
      MPI_LOR,MPI_COMM_WORLD,ier)
    
  else
    
    call allocate_lbuffer(narr)
    
    call MPI_ALLREDUCE(argument,lbuffer,narr,MYLOGICAL, &
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

 subroutine or1_world_larr(argument,narr)
 
!***********************************************************************
!     
!     LBsoft global 'logical or' subroutine for a logical array
!     of kind 1
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  logical(kind=1), intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  logical, allocatable, dimension(:) :: temp

  integer ier

  call allocate_lbuffer(narr)
  allocate(temp(narr))
  temp(1:narr) = argument(1:narr)

  call MPI_ALLREDUCE(temp,lbuffer,narr,MYLOGICAL, &
    MPI_LOR,MPI_COMM_WORLD,ier)

  argument(1:narr) = lbuffer(1:narr)
  deallocate(temp)

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
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
  
  CALL MPI_IRECV(argument,1,MYINTEGER, &
   idsource,tag,MPI_COMM_WORLD,irequest,ierr)
  
  return
  
 end subroutine irecv_world_i
 
 subroutine irecv_world_f(argument,idsource,tag,irequest)
 
!***********************************************************************
!     
!     LBsoft global subroutine for receiving a float value
!     by a nonblocking receive operation
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
  
  CALL MPI_IRECV(argument,1,MYFLOAT, &
   idsource,tag,MPI_COMM_WORLD,irequest,ierr)
  
  return
  
 end subroutine irecv_world_f
 
 subroutine irecv_world_l(argument,idsource,tag,irequest)
 
!***********************************************************************
!     
!     LBsoft global subroutine for receiving a logical value
!     by a nonblocking receive operation
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
  
  CALL MPI_IRECV(argument,1,MYLOGICAL, &
   idsource,tag,MPI_COMM_WORLD,irequest,ierr)
  
  return
  
 end subroutine irecv_world_l
 
 subroutine isend_world_i(argument,iddest,tag,irequest)
 
!***********************************************************************
!     
!     LBsoft global subroutine for sending an integer value
!     by a nonblocking send operation
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
  
  CALL MPI_ISEND(argument,1,MYINTEGER, &
   iddest,tag,MPI_COMM_WORLD,irequest,ierr)
  
  return
  
 end subroutine isend_world_i
 
 subroutine isend_world_f(argument,iddest,tag,irequest)
 
!***********************************************************************
!     
!     LBsoft global subroutine for sending a float value
!     by a nonblocking send operation
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
  
  CALL MPI_ISEND(argument,1,MYFLOAT, &
   iddest,tag,MPI_COMM_WORLD,irequest,ierr)
  
  return
  
 end subroutine isend_world_f
 
 subroutine isend_world_l(argument,iddest,tag,irequest)
 
!***********************************************************************
!     
!     LBsoft global subroutine for sending a logical value
!     by a nonblocking send operation
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
  
  CALL MPI_ISEND(argument,1,MYLOGICAL, &
   iddest,tag,MPI_COMM_WORLD,irequest,ierr)
  
  return
  
 end subroutine isend_world_l
 
 subroutine irecv_world_iarr(argument,narr,idsource,tag,irequest)
 
!***********************************************************************
!     
!     LBsoft global subroutine for receiving an integer array
!     by a nonblocking receive operation
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
  
  CALL MPI_IRECV(argument,narr,MYINTEGER, &
   idsource,tag,MPI_COMM_WORLD,irequest,ierr)
  
  return
  
 end subroutine irecv_world_iarr
 
 subroutine irecv_world_farr(argument,narr,idsource,tag,irequest)
 
!***********************************************************************
!     
!     LBsoft global subroutine for receiving a float array
!     by a nonblocking receive operation
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
  
  CALL MPI_IRECV(argument,narr,MYFLOAT, &
   idsource,tag,MPI_COMM_WORLD,irequest,ierr)
  
  return
  
 end subroutine irecv_world_farr
 
 subroutine irecv_world_larr(argument,narr,idsource,tag,irequest)
 
!***********************************************************************
!     
!     LBsoft global subroutine for receiving a logical array
!     by a nonblocking receive operation
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
  
  CALL MPI_IRECV(argument,narr,MYLOGICAL, &
   idsource,tag,MPI_COMM_WORLD,irequest,ierr)
  
  return
  
 end subroutine irecv_world_larr
 
 subroutine isend_world_iarr(argument,narr,iddest,tag,irequest)
 
!***********************************************************************
!     
!     LBsoft global subroutine for sending an integer array
!     by a nonblocking send operation
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
  
  CALL MPI_ISEND(argument,narr,MYINTEGER, &
   iddest,tag,MPI_COMM_WORLD,irequest,ierr)
  
  return
  
 end subroutine isend_world_iarr
 
 subroutine isend_world_farr(argument,narr,iddest,tag,irequest)
 
!***********************************************************************
!     
!     LBsoft global subroutine for sending a float array
!     by a nonblocking send operation
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
  
  CALL MPI_ISEND(argument,narr,MYFLOAT, &
   iddest,tag,MPI_COMM_WORLD,irequest,ierr)
  
  return
  
 end subroutine isend_world_farr
 
 subroutine isend_world_larr(argument,narr,iddest,tag,irequest)
 
!***********************************************************************
!     
!     LBsoft global subroutine for sending a logical array
!     by a nonblocking send operation
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
  
  CALL MPI_ISEND(argument,narr,MYLOGICAL, &
   iddest,tag,MPI_COMM_WORLD,irequest,ierr)
  
  return
  
 end subroutine isend_world_larr
 
 subroutine wait_world(argument)
 
!***********************************************************************
!     
!     LBsoft global subroutine for waiting a nonblocking operation 
!     to complete.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
  
  call MPI_File_write_at(iotest,myoffset,header(1:nn),nn, &
   MPI_CHARACTER,MPI_STATUS_IGNORE,e_io)
  
  return
  
 endsubroutine print_header_vtk_par
 
 subroutine print_footer_vtk_par(iotest,offsetsub,footer,e_io)
 
!***********************************************************************
!     
!     LBsoft subroutine for writing the footer part of
!     in VTK legacy file in parallel IO
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
  
  call MPI_File_write_at(iotest,myoffset,footer(1:30),30, &
   MPI_CHARACTER,MPI_STATUS_IGNORE,e_io)
  
  return
  
 end subroutine print_footer_vtk_par
 
 subroutine setupcom(nx,ny,nz,nbuff,ibctype,ixpbc,iypbc,izpbc,minx,maxx, &
   miny,maxy,minz,maxz,lsingle_fluid)
  
!***********************************************************************
!     
!     LBsoft subroutine for the initial setup of the domain 
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Bernaschi
!     last modification October 2019
!     
!***********************************************************************
  
  INTEGER, intent(in) :: nx,ny,nz,nbuff,ibctype
  INTEGER, intent(inout) :: ixpbc,iypbc,izpbc,minx,maxx,miny,maxy, &
   minz,maxz
  LOGICAL, intent(IN) :: lsingle_fluid
   
  
  
  INTEGER :: i,j,k,itemp
  INTEGER :: ownernlb,ownernub
  logical :: ltest(1)
  INTEGER(kind=IPRC) :: i4
  
  nx_comm=nx
  ny_comm=ny
  nz_comm=nz
  
  globalDims(1)=nx
  globalDims(2)=ny
  globalDims(3)=nz
  
  nxy2=(nx+2*nbuff)*(ny+2*nbuff)
  nx2=nx+(2*nbuff)
  nbuff_comm=nbuff
  ALLOCATE(countnpp(0:mxrank-1))
  DO i=0,mxrank-1
     countnpp(i)=0
  ENDDO
  ownernlb=1 !(1-nbuff)+((1-nbuff)*nx2)+((1-nbuff)*nxy2)
  ownernub=(nx+2*nbuff)*(ny+2*nbuff)*(nz+2*nbuff)
!max  write(0,*),'NX=',nx,'NY=',ny,'NZ=',nz,'ownern lb=',ownernlb,'ownern up=',ownernub
  ALLOCATE(ownern(ownernlb:ownernub))
  do i=ownernlb,ownernub
     ownern(i)=-1
  enddo
  
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
    
    call cartdeco(nx,ny,nz,nbuff,.true.,ownern,ownernlb,minx,maxx,miny,maxy, &
     minz,maxz)
    
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
    

!   check if node ownerships are assigned in natural order
    ltest(1)=.false.
    do k=minz-nbuff,maxz+nbuff
      do j=miny-nbuff,maxy+nbuff
        do i=minx-nbuff,maxx+nbuff
          i4=i4back(i,j,k)
          itemp=ownernfind(i,j,k,mxrank,gminx,gmaxx, &
           gminy,gmaxy,gminz,gmaxz)
          if(itemp.ne.int(ownern(i4)))then
            write(6,*)'wrong assignment of node ownership'
            write(6,*)'expected',idrank,i,j,k,itemp
            write(6,*)'obtained',idrank,i,j,k,int(ownern(i4))
            call flush(6)
            ltest(1)=.true.
            goto 121
          endif
        enddo
      enddo
    enddo
 121 continue
    call or_world_larr(ltest,1)
    if(ltest(1))then
      call finalize_world
      stop
    endif

    ldims(1)=maxx-minx+1
    ldims(2)=maxy-miny+1
    ldims(3)=maxz-minz+1
    
    mymin(1)=minx
    mymin(2)=miny
    mymin(3)=minz
    mymax(1)=maxx
    mymax(2)=maxy
    mymax(3)=maxz
    
    return

 end subroutine setupcom
 
 SUBROUTINE cartdeco(nx,ny,nz,nbuff,first, ownern, ownernlb,minx,maxx,miny,maxy, &
   minz,maxz)

!***********************************************************************
!     
!     LBsoft subroutine for defining the domain decomposition for
!     MPI ranks 
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Bernaschi
!     last modification October 2019
!     
!***********************************************************************

    IMPLICIT NONE
    
    integer, intent(in) :: nx,ny,nz,nbuff
    LOGICAL,INTENT(in) :: first
    ! INTEGER(kind=IPRC) :: tempown(0:*)
    ! INTEGER :: order(0:*)
    INTEGER,INTENT(in):: ownernlb
    INTEGER(kind=2) :: ownern(ownernlb:*)
    INTEGER, intent(inout) :: minx,maxx,miny,maxy,minz,maxz

    INTEGER(kind=IPRC) :: i4,i4target
    INTEGER :: i,j,k,itemp,jtemp,ktemp
    INTEGER :: ix,iy,iz
    INTEGER :: currproc, currpoff=0, curroffz, curroffy, curroffx
    INTEGER :: slice=0, slicez=0, slicey=0, slicex=0
    INTEGER :: m, iandex
    LOGICAL :: ltest(1)

   
    IF(first) THEN
       
       minx=LARGEINT
       miny=LARGEINT
       minz=LARGEINT
       maxx=-LARGEINT
       maxy=-LARGEINT
       maxz=-LARGEINT
       
       IF(mxrank.EQ.1) THEN
          domdec=1
          nprocx=1
          nprocy=1
          nprocz=1
       ENDIF

       IF(domdec.LT.1.OR.domdec.GT.7) THEN

          IF(mxrank.GT.1) CALL MPI_FINALIZE(ierr)

          IF(idrank.EQ.0) THEN
             PRINT*,'decomposition:', domdec
          ENDIF
          STOP
       ENDIF

       IF(domdec.GT.3.AND.MOD(mxrank,2).NE.0) THEN

          IF(mxrank.GT.1) CALL MPI_FINALIZE(ierr)

          IF(idrank.EQ.0) THEN
             PRINT*,'Invalid number of processors',mxrank,&
                  'for this decomposition'
          ENDIF
          STOP
       ENDIF
       IF(domdec.GT.3.AND.mxrank.LT.4) THEN

          IF(mxrank.GT.1) CALL MPI_FINALIZE(ierr)

          IF(idrank.EQ.0) THEN
             PRINT*,&
                  'At least 4 processors are required for this decomposition'
          ENDIF
          STOP
       ENDIF

       IF(domdec.EQ.7) THEN
          IF(mxrank.LT.8) THEN

             IF(mxrank.GT.1) CALL MPI_FINALIZE(ierr)

             IF(idrank.EQ.0) THEN
                PRINT*,&
                     'At least 8 processors are required for this decomposition'
             ENDIF
             STOP
          ENDIF
          IF(mxrank.EQ.8) THEN
             nprocz=2
             nprocy=2
             nprocx=2
          ELSE IF(mxrank.EQ.16) THEN
             nprocz=2
             nprocy=2
             nprocx=4
          ELSE IF(mxrank.EQ.32) THEN
             nprocz=2
             nprocy=4
             nprocx=4
          ELSE IF(mxrank.EQ.64) THEN
             nprocz=4
             nprocy=4
             nprocx=4
          ELSE IF(mxrank.EQ.128) THEN
             nprocz=4
             nprocy=4
             nprocx=8
          ELSE IF(mxrank.EQ.256) THEN
             nprocz=4
             nprocy=8
             nprocx=8
          ELSE IF(mxrank.EQ.512) THEN
             nprocz=8
             nprocy=8
             nprocx=8
          ELSE IF(mxrank.EQ.1024) THEN
             nprocz=8
             nprocy=8
             nprocx=16
          ELSE IF(mxrank.EQ.2048) THEN
             nprocz=8
             nprocy=16
             nprocx=16
          ELSE IF(mxrank.EQ.4096) THEN
             nprocz=16
             nprocy=16
             nprocx=16
          ELSE IF(mxrank.EQ.8192) THEN
             nprocz=16
             nprocy=16
             nprocx=32
          ELSE IF(mxrank.EQ.16384) THEN
             nprocz=16
             nprocy=32
             nprocx=32
          ELSE IF(mxrank.EQ.32768) THEN
             nprocz=32
             nprocy=32
             nprocx=32
          ELSE IF(idrank.EQ.0) THEN
             PRINT*,"Invalid number of processors",mxrank

             IF(mxrank.GT.1) CALL MPI_FINALIZE(ierr)

             STOP
          ENDIF
       ENDIF

       m=0
    ENDIF
    
    currproc=0
    curroffz=0
    curroffy=0
    slice=0; slicez=0; slicey=0; slicex=0
    IF(domdec.EQ.1) THEN      !along z
       slice=(nz)/mxrank
       IF(MOD((nz),mxrank).NE.0) THEN
          slice=slice+1
       ENDIF
    ENDIF
    IF(domdec.EQ.7) THEN      !along zyx
       slicez=(nz)/nprocz
       IF(MOD((nz),nprocz).NE.0) THEN
          slicez=slicez+1
       ENDIF
    ENDIF
    DO iz=1, nz
       IF(domdec.EQ.1) THEN   !along z
          IF(slice.EQ.0.AND.iz.LT.(nz)) THEN
             currproc=currproc+1
             slice=((nz)-(iz))/(mxrank-currproc)-1
             IF(MOD((nz)-(iz),(mxrank-currproc)).NE.0) THEN
                slice=slice+1
             ENDIF
          ELSE
             slice=slice-1
          ENDIF
       ELSE IF(domdec.EQ.4.OR.domdec.EQ.5) THEN !along zy and zx
          IF(iz.GT.(nz)/2) THEN
             currpoff=mxrank/2
          ELSE
             currpoff=0
          ENDIF
          IF(domdec.EQ.4) THEN !along zy
             currproc=currpoff
             slice=(ny)/(mxrank/2)
             IF(MOD((ny),(mxrank/2)).NE.0) THEN
                slice=slice+1
             ENDIF
          ENDIF
       ELSE IF(domdec.EQ.7) THEN !along zyx
          IF(slicez.EQ.0.AND.iz.LT.(nz)) THEN
             curroffz=curroffz+1
             slicez=((nz)-(iz))/(nprocz-curroffz)
             IF(MOD((nz)-(iz),(nprocz-curroffz)).NE.0) THEN
                slicez=slicez+1
             ENDIF
          ENDIF
          slicez=slicez-1
          curroffy=0
          slicey=(ny)/nprocy
          IF(MOD((ny),nprocy).NE.0) THEN
             slicey=slicey+1
          ENDIF
       ELSE IF(domdec.EQ.2) THEN !along y
          currproc=0
          slice=(ny)/mxrank
          IF(MOD((ny),mxrank).NE.0) THEN
             slice=slice+1
          ENDIF
       ENDIF
       DO iy=1, ny
          IF(domdec.EQ.2) THEN
             IF(slice.EQ.0.AND.iy.LT.(ny)) THEN
                currproc=currproc+1
                slice=((ny)-(iy))/(mxrank-currproc)-1
                IF(MOD((ny)-(iy),(mxrank-currproc)).NE.0) THEN
                   slice=slice+1
                ENDIF
             ELSE
                slice=slice-1
             ENDIF
          ELSE IF(domdec.EQ.4) THEN !along zy
             IF(slice.EQ.0.AND.iy.LT.(ny)) THEN
                currproc=currproc+1
                slice=((ny)-(iy))/ &
                     ((mxrank/2)-(currproc-currpoff))-1
                IF(MOD((ny)-(iy),&
                     (mxrank/2)-(currproc-currpoff)).NE.0) THEN
                   slice=slice+1
                ENDIF
             ELSE
                slice=slice-1
             ENDIF
          ELSE IF(domdec.EQ.6) THEN !along yx
             IF(iy.GT.(ny)/2) THEN
                currpoff=mxrank/2
             ELSE
                currpoff=0
             ENDIF
             currproc=currpoff
             slice=(nx)/(mxrank/2)
             IF(MOD((nx),(mxrank/2)).NE.0) THEN
                slice=slice+1
             ENDIF
          ELSE IF(domdec.EQ.7) THEN !along zyx
             IF(slicey.EQ.0.AND.iy.LT.(ny)) THEN
                curroffy=curroffy+1
                slicey=((ny)-(iy))/(nprocy-curroffy)
                IF(MOD((ny)-(iy),(nprocy-curroffy)).NE.0) THEN
                   slicey=slicey+1
                ENDIF
             ENDIF
             slicey=slicey-1
             curroffx=0
             slicex=(nx)/nprocx
             IF(MOD((nx),nprocx).NE.0) THEN
                slicex=slicex+1
             ENDIF
          ELSE IF(domdec.EQ.3) THEN !along x
             currproc=0
             slice=(nx)/mxrank
             IF(MOD((nx),mxrank).NE.0) THEN
                slice=slice+1
             ENDIF
          ELSE IF(domdec.EQ.5) THEN !along zx
             currproc=currpoff
             slice=(nx)/(mxrank/2)
             IF(MOD((nx),(mxrank/2)).NE.0) THEN
                slice=slice+1
             ENDIF
          ENDIF
          DO ix=1, nx
             IF(domdec.EQ.7) THEN
                currproc=curroffx+curroffy*nprocx+curroffz*nprocx*nprocy
             ENDIF
             IF(currproc.LT.0.OR.currproc.GE.mxrank) THEN
                PRINT*,'Invalid currproc',currproc,nprocx,nprocy,nprocz,&
                     curroffx, curroffy,curroffz
                STOP
             ENDIF
             ! i4=iz*nxy2 + iy*nx2 + ix
             i4=i4back(ix,iy,iz)
             IF(domdec.EQ.3) THEN !along x
                IF(slice.EQ.0.AND.ix.LT.(nx)) THEN
                   currproc=currproc+1
                   slice=((nx)-(ix))/(mxrank-currproc)-1
                   IF(MOD((nx)-(ix),(mxrank-currproc)).NE.0) THEN
                      slice=slice+1
                   ENDIF
                ELSE
                   slice=slice-1
                ENDIF
             ELSE IF(domdec.EQ.5.OR.domdec.EQ.6) THEN !along zx or along yx
                IF(slice.EQ.0.AND.ix.LT.(nx)) THEN
                   currproc=currproc+1
                   slice=((nx)-(ix))/ &
                        ((mxrank/2)-(currproc-currpoff))-1
                   IF(MOD((nx)-(ix),&
                        (mxrank/2)-(currproc-currpoff)).NE.0) THEN
                      slice=slice+1
                   ENDIF
                ELSE
                   slice=slice-1
                ENDIF
             ELSE IF(domdec.EQ.7) THEN !along zyx
                IF(slicex.EQ.1.AND.ix.LT.(nx)) THEN
                   curroffx=curroffx+1
                   slicex=((nx)-(ix))/ &
                        (nprocx-curroffx)
                   IF(MOD((nx)-(ix),(nprocx-curroffx)).NE.0) THEN
                      slicex=slicex+1
                   ENDIF
                ELSE
                   slicex=slicex-1
                ENDIF
             ENDIF
             IF(first) THEN
                countnpp(currproc)=countnpp(currproc)+1
                ownern(i4)=currproc
                IF(currproc.eq.idrank) THEN ! check the local boundaries
                   IF(ix.LT.minx) THEN
                      minx=ix
                   ENDIF
                   IF(iy.LT.miny) THEN
                      miny=iy
                   ENDIF
                   IF(iz.LT.minz) THEN
                      minz=iz
                   ENDIF
                   IF(ix.GT.maxx) THEN
                      maxx=ix
                   ENDIF
                   IF(iy.GT.maxy) THEN
                      maxy=iy
                   ENDIF
                   IF(iz.GT.maxz) THEN
                      maxz=iz
                   ENDIF
                ENDIF
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    
    if(minx.le.1) then
        minx=1
     endif
     if(miny.le.1) then
        miny=1
     endif
     if(minz.le.1) then
        minz=1
     endif
     if(maxx.ge.nx) then
        maxx=nx
     endif
     if(maxy.ge.ny) then
        maxy=ny
     endif
     if(maxz.ge.nz) then
        maxz=nz
     endif
    
#if 0
    do i=0,mxrank-1
      if(i==idrank)then
        write(6,*)idrank,minx,maxx,miny,maxy,minz,maxz
        call flush(6)
      endif
      call get_sync_world
    enddo
#endif
    
    ltest(1)=.false.
    
    do k=1-nbuff,nz+nbuff
       do j=1-nbuff,ny+nbuff
         do i=1-nbuff,nx+nbuff
           if(i<1.or.j<1.or.k<1.or.i>nx.or.j>ny.or.k>nz)then
             itemp=i
             jtemp=j
             ktemp=k
             if(ipbc(1)==1)then
               if(i<1)itemp=itemp+nx
               if(i>nx)itemp=itemp-nx
             else
               if(i<1)itemp=1
               if(i>nx)itemp=nx
             endif
             if(ipbc(2)==1)then
               if(j<1)jtemp=jtemp+ny
               if(j>ny)jtemp=jtemp-ny
             else
               if(j<1)jtemp=1
               if(j>ny)jtemp=ny
             endif
             if(ipbc(3)==1)then
               if(k<1)ktemp=ktemp+nz
               if(k>nz)ktemp=ktemp-nz
             else
               if(k<1)ktemp=1
               if(k>nz)ktemp=nz
             endif
             if(itemp<1.or.jtemp<1.or.ktemp<1.or.itemp>nx.or.jtemp>ny.or.ktemp>nz)then
               write(6,*)'ERROR in building the ownern frame'
               ltest(1)=.true.
               goto 121
             endif
             i4target=i4back(i,j,k)
             i4=i4back(itemp,jtemp,ktemp)
             if(i4target<0 .or. i4<0)then
               write(6,*)'ERROR in cartdeco id ',idrank,' ijk ', &
                i,j,k,' i4target ',i4target,' ijktemp ',itemp,jtemp,ktemp,i4
               ltest(1)=.true.
               goto 121
             else
               ownern(i4target)=ownern(i4)
               currproc=ownern(i4)
               countnpp(currproc)=countnpp(currproc)+1
             endif
           endif
         enddo
       enddo
     enddo
121  continue
     call or_world_larr(ltest,1)
     if(ltest(1))then
       call finalize_world
       stop
     endif
    

!   check the natural order of the subdoumains
    ltest(1)=.false.
    do k=1,nprocz
      do j=1,nprocy
        do i=1,nprocx
          ix=idnodeback(i,j,k)
          if(ix==idrank)then
            i4=i4back(minx,miny,minz)
            if(ownern(i4).ne.idrank)then
              write(6,*)'wrong order of subdomains'
              write(6,*)idrank,minx,maxx,miny,maxy,minz,maxz
              write(6,*)'expected',idrank,idnodefind(idrank),idrank
              write(6,*)'obtained',idrank,idnodefind(int(ownern(i4))),ownern(i4)
              call flush(6)
              ltest(1)=.true.
            endif
          endif
        enddo
      enddo
    enddo
    call or_world_larr(ltest,1)
    if(ltest(1))then
      call finalize_world
      stop
    endif

    
    
  end subroutine cartdeco
 
 subroutine set_domdec(itemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the protected value of domdec
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
 
 function idnodefind(i4sub) result(i4find)

!***********************************************************************
!     
!     LBsoft subroutine for finding the id of MPI rank provided 
!     the liquid node in i4-form
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Bernaschi
!     last modification October 2019
!     
!***********************************************************************
  
  implicit none
    
  integer, intent(in) :: i4sub
  integer, dimension(3) :: i4find
    
  i4find(3) = (i4sub/nprocx*nprocy)
  i4find(2) = ((i4sub - i4find(3)*nprocx*nprocy)/nprocx)
  i4find(1) = (i4sub - i4find(3)*nprocx*nprocy - i4find(2)*nprocx)
    
  return
  
 end function idnodefind
 
 function idnodeback(i,j,k) result(myout)
 
!***********************************************************************
!     
!     LBsoft subroutine for finding the id of MPI rank provided 
!     the node of the cartesian domain decomposition
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Bernaschi
!     last modification October 2019
!     
!***********************************************************************
    
    implicit none
    integer,intent(in) :: i !< i-location on mesh
    integer,intent(in) :: j !< j-location on mesh
    integer,intent(in) :: k !< k-location on mesh
    integer :: myout
     
    myout = i+j*nprocx+k*nprocx*nprocy-1

 end function idnodeback
  
 pure function ownernfind_arr(i,j,k,nxs,nys,nzs,nbuffs,ownerns)
  
!***********************************************************************
!     
!     LBsoft subroutine for finding the owner of a node provided 
!     in i,j,k triplet-form taken from ownerns array
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
  
 implicit none
   
  integer, intent(in) :: i,j,k,nxs,nys,nzs,nbuffs
  integer(kind=2), allocatable, intent(in) :: ownerns(:)
  integer(kind=IPRC) :: i4
   
  integer :: ownernfind_arr
   
  i4=int(k+nbuffs-1,kind=IPRC)*(nxs+2*nbuffs)*(nys+2*nbuffs)+ &
   int(j+nbuffs-1,kind=IPRC)*(nxs+2*nbuffs)+int(i+nbuffs,kind=IPRC)
    
  ownernfind_arr=int(ownerns(i4))
   
  return
   
 end function ownernfind_arr
 
 pure function ownernfind(i,j,k,mxranksub,gminxs,gmaxxs, &
     gminys,gmaxys,gminzs,gmaxzs) result(i4find)
 
!***********************************************************************
!     
!     LBsoft subroutine for finding the owner of a node provided 
!     in i,j,k triplet-form computed on-fly
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
    
  if(gminxs(-1)==1)then
    if(i<1)itemp=itemp+gmaxxs(-1)
    if(i>gmaxxs(-1))itemp=itemp-gmaxxs(-1)
  else
    if(i<1)itemp=1
    if(i>gmaxxs(-1))itemp=gmaxxs(-1)
  endif
  if(gminys(-1)==1)then
     if(j<1)jtemp=jtemp+gmaxys(-1)
     if(j>gmaxys(-1))jtemp=jtemp-gmaxys(-1)
  else
  if(j<1)jtemp=1
    if(j>gmaxys(-1))jtemp=gmaxys(-1)
  endif
  if(gminzs(-1)==1)then
    if(k<1)ktemp=ktemp+gmaxzs(-1)
    if(k>gmaxzs(-1))ktemp=ktemp-gmaxzs(-1)
  else
    if(k<1)ktemp=1
    if(k>gmaxzs(-1))ktemp=gmaxzs(-1)
  endif
    
  do l=0,mxranksub-1
    lt(1:3)=.false.
    if(itemp.ge.gminxs(l) .and. itemp.le.gmaxxs(l))lt(1)=.true.
    if(jtemp.ge.gminys(l) .and. jtemp.le.gmaxys(l))lt(2)=.true.
    if(ktemp.ge.gminzs(l) .and. ktemp.le.gmaxzs(l))lt(3)=.true.
    if(all(lt))then
      i4find=l
      return
    endif
  enddo
    
  i4find=-1
    
  return
  
 end function ownernfind
 
 subroutine create_findneigh_list_single_halo(nx,ny,nz,nbuff,ibctype, &
   ixpbc,iypbc,izpbc,minx,maxx,miny,maxy,minz,maxz)
   
!***********************************************************************
!     
!     LBsoft subroutine for driving the creation of the list of 
!     the HVAR exchange by MPI sending and receiving
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Bernaschi
!     last modification October 2019
!     
!***********************************************************************

  implicit none
 
  INTEGER, intent(in) :: nx,ny,nz,nbuff,ibctype
  INTEGER, intent(inout) :: ixpbc,iypbc,izpbc,minx,maxx,miny,maxy, &
   minz,maxz
  

  call findneigh_single_halo(nx,ny,nz,ibctype,ixpbc,iypbc,izpbc, &
   minx,maxx,miny,maxy,minz,maxz)
   
  return
  
 end subroutine create_findneigh_list_single_halo
 
 subroutine create_findneigh_list_hvar_isfluid(nx,ny,nz,nbuff,ibctype, &
   ixpbc,iypbc,izpbc,minx,maxx,miny,maxy,minz,maxz)
   
!***********************************************************************
!     
!     LBsoft subroutine for driving the creation of the list of 
!     the HVAR exchange by MPI sending and receiving
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Bernaschi
!     last modification October 2019
!     
!***********************************************************************

  implicit none
 
  INTEGER, intent(in) :: nx,ny,nz,nbuff,ibctype
  INTEGER, intent(inout) :: ixpbc,iypbc,izpbc,minx,maxx,miny,maxy, &
   minz,maxz
  

  call findneigh_hvar_isfluid(nx,ny,nz,nbuff,ibctype,ixpbc,iypbc,izpbc, &
   minx,maxx,miny,maxy,minz,maxz)
   
  return
  
 end subroutine create_findneigh_list_hvar_isfluid

 subroutine create_findneigh_list_pops(nx,ny,nz,nbuff,ibctype,isfluid, &
   ixpbc,iypbc,izpbc,minx,maxx,miny,maxy,minz,maxz)
!***********************************************************************
!     
!     LBsoft subroutine for driving the creation of the list of 
!     the populations exchange by MPI sending and receiving
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Bernaschi
!     last modification October 2019
!     
!***********************************************************************

  implicit none
  
  INTEGER, intent(in) :: nx,ny,nz,nbuff,ibctype
  INTEGER(kind=1), intent(in), allocatable, dimension(:,:,:) :: isfluid
  INTEGER, intent(inout) :: ixpbc,iypbc,izpbc,minx,maxx,miny,maxy, &
   minz,maxz
  
  
  call findneigh(nx,ny,nz,nbuff,ibctype,isfluid,ixpbc,iypbc,izpbc,minx, &
   maxx,miny,maxy,minz,maxz)
   
  return
  
 end subroutine create_findneigh_list_pops
 
 SUBROUTINE findneigh_single_halo(nx,ny,nz,ibctype,ixpbc,iypbc,izpbc, &
   minx,maxx,miny,maxy,minz,maxz)

!***********************************************************************
!     
!     LBsoft subroutine for creating the list of a single node of
!     halo to exchange data by MPI sending and receiving
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Bernaschi
!     last modification October 2019
!     
!***********************************************************************

    IMPLICIT NONE

    INTEGER :: request(0:maxneigh-1)
    
    INTEGER :: status(MPI_STATUS_SIZE,0:maxneigh-1)
    
    
    INTEGER, intent(in) :: nx,ny,nz,ibctype
    INTEGER, intent(inout) :: ixpbc,iypbc,izpbc,minx,maxx,miny,maxy, &
   minz,maxz
    
    INTEGER :: lminx, lmaxx, lminy, lmaxy, lminz, lmaxz
    INTEGER(kind=IPRC) :: i4,i4temp,i4shift
    INTEGER :: i,j,jy,k,ip,maxsender,maxreceiver
    INTEGER :: ishift,jshift,kshift
    INTEGER(kind=IPRC) :: nghb
    INTEGER :: n, l, m,ks,mym(3)
    INTEGER :: ifold,jfold,kfold
    INTEGER :: itemp,jtemp,ktemp
    INTEGER :: itemp2,jtemp2,ktemp2
    INTEGER :: firstminpe,firstmaxpe
    INTEGER :: i4_send,i4_recv

    INTEGER :: allneigh(0:maxneigh*mxrank-1)

    INTEGER :: sender(0:mxrank-1)
    INTEGER :: receiver(0:mxrank-1)
    
    logical, parameter :: lverbose=.false.
    logical :: ltemp(1)
    
    integer, parameter :: nsing=1
    
    if(mxrank==1)return
    
    DO i=0,mxrank-1
       !sender(i)=0
       receiver(i)=0
    ENDDO
    DO i=0,maxneigh-1
       i_pe2send_fluid_sing(i)=-1
       n_sing2send_fluid(i)=-1
    ENDDO
    DO i=0,maxneigh-1
       i_pe2recv_fluid_sing(i)=-1
       n_sing2recv_fluid(i)=-1
    ENDDO
    n_pe2recv_fluid_sing=0
    n_pe2send_fluid_sing=0
    
    if(minx.eq.1) then
       if(ixpbc.eq.1) then
          lminx=1
       else
          lminx=0
       endif
    else
       lminx=minx
    endif
    if(miny.eq.1) then
       if(iypbc.eq.1) then
          lminy=1
       else
          lminy=0
       endif
    else
       lminy=miny
    endif
    if(minz.eq.1) then
       if(izpbc.eq.1) then
          lminz=1
       else
          lminz=0
       endif
    else
       lminz=minz
    endif
    if(maxx.eq.nx) then
       if(ixpbc.eq.1) then
          lmaxx=nx
       else
          lmaxx=nx+1
       endif
    else
       lmaxx=maxx
    endif
    if(maxy.eq.ny) then
       if(iypbc.eq.1) then
          lmaxy=ny
       else
          lmaxy=ny+1
       endif
    else
       lmaxy=maxy
    endif
    if(maxz.eq.nz) then
       if(izpbc.eq.1) then
          lmaxz=nz
       else
          lmaxz=nz+1
       endif
    else
       lmaxz=maxz
    endif
    
    
    if(lverbose)then
      call flush(6)
      call get_sync_world
      do i=0,mxrank-1
        if(i==idrank)then
          write(6,*)'id=',idrank,'lminx=',lminx,'lmaxx=',lmaxx,'lminy=',lminy, &
              'lmaxy=',lmaxy,'lminz=',lminz,'lmaxz=',lmaxz  
        endif
        call flush(6)
        call get_sync_world
      enddo
    endif
    
    ! receiver(i): count no. of neighboring fluid nodes that talk with i-th PE
    !determine what I have to receive from each node/PE around me
     ltemp(1)=.false.
     do k=minz-nsing,maxz+nsing
       do j=miny-nsing,maxy+nsing
         do i=minx-nsing,maxx+nsing
           if(i<minx.or.j<miny.or.k<minz.or.i>maxx.or.j>maxy.or.k>maxz)then
                itemp=i
                jtemp=j
                ktemp=k
                itemp2=i
                jtemp2=j
                ktemp2=k
                !apply periodic conditions if necessary
                if(ixpbc.eq.1) then
                   if(itemp<1) then
                      itemp=itemp+nx
                   endif
                   if(itemp>nx) then
                      itemp=itemp-nx
                   endif
                endif
                if(iypbc.eq.1) then
                   if(jtemp<1) then
                      jtemp=jtemp+ny
                   endif
                   if(jtemp>ny) then
                      jtemp=jtemp-ny
                   endif
                endif
                if(izpbc.eq.1) then
                   if(ktemp<1) then
                      ktemp=ktemp+nz
                   endif
                   if(ktemp>nz) then
                      ktemp=ktemp-nz
                   endif
                endif
                i4=i4back(itemp,jtemp,ktemp)
                if(i4<1)then
                  write(6,*)idrank,' something bad ',itemp2,jtemp2,ktemp2, &
                   itemp,jtemp,ktemp,minx,maxx
                  ltemp(1)=.true.
                endif
                if(ownern(i4)<0 .or.ownern(i4)>(mxrank-1))then
                  write(6,*)idrank,' something bad ',itemp2,jtemp2,ktemp2, &
                   itemp,jtemp,ktemp,i4,ownern(i4)
                  ltemp(1)=.true.
                endif
                IF(ownern(i4).NE.idrank) THEN
                   receiver(ownern(i4))=receiver(ownern(i4))+1
                ENDIF
             endif
          enddo
       enddo
    enddo
    call or_world_larr(ltemp,1)
    if(ltemp(1))then
      call finalize_world
      stop
    endif
    
    
    if(lverbose)then
      call flush(6)
      call get_sync_world
      do i=0,mxrank-1
        if(i==idrank)then
          write(6,*)'id=',idrank,'receiver=',receiver(0:mxrank-1)
        endif
        call flush(6)
        call get_sync_world
      enddo
    endif
    
    !determine the max number of fluid nodes for node/PE that I have to receive 
    j=0
    maxreceiver=0
    DO i=0,mxrank-1

       IF(receiver(i).GT.0) THEN ! select procs that talks to this one

          IF(j.GE.maxneigh) THEN
             PRINT *,'Task',idrank,'Invalid index for myreceiver',j
             STOP
          ENDIF

          IF(receiver(i).GT.maxreceiver) THEN
             maxreceiver=receiver(i)  ! get the max no. of nodes to talk to
          ENDIF

          i_pe2recv_fluid_sing(j)=i    ! index=process # rank of the list of fluid nodes to talk to (that should be sent)
          n_sing2recv_fluid(j)=0        ! no. of fluid nodes to talk to (initialization)
          receiver(i)=j                ! id of PE to talk in local ordering
          j=j+1                        ! local ordering index
       ELSE
          receiver(i)=-1
       ENDIF

    ENDDO
    n_pe2recv_fluid_sing=j
    
    !allocate buffer
    IF(maxreceiver.GT.0.AND.n_pe2recv_fluid_sing.GT.0) THEN
       ALLOCATE(i_sing2recv_fluid(0:1,0:maxreceiver-1,0:n_pe2recv_fluid_sing-1))
       ALLOCATE(buffr_sing(0:maxreceiver-1,0:n_pe2recv_fluid_sing-1))
    ENDIF
    
    !initialization
    DO i=0,maxneigh-1
       i_pe2send_fluid_sing(i)=-1
    ENDDO

    CALL MPI_ALLGATHER(i_pe2recv_fluid_sing,maxneigh,MPI_INTEGER,allneigh, &
     maxneigh,MPI_INTEGER, MPI_COMM_WORLD, ierr)

    !allneigh = i_pe2recv_fluid_sing

    
    !determine to who I have to send
    j=0
    DO i=0,mxrank*maxneigh-1
       IF(allneigh(i)==idrank) THEN
          IF(j.GE.maxneigh) THEN
             PRINT *,'Task',idrank,'Invalid index for i_pe2send_fluid_sing',j
             STOP
          ENDIF
          i_pe2send_fluid_sing(j)=i/maxneigh
          j=j+1
       ENDIF
    ENDDO
    n_pe2send_fluid_sing=j
    
    if(lverbose)then
      call flush(6)
      call get_sync_world
      do i=0,mxrank-1
        if(i==idrank)then
          write(6,*)'id=',idrank,'n_pe2send_fluid_sing=', &
           n_pe2send_fluid_sing,i_pe2send_fluid_sing(0:n_pe2send_fluid_sing-1)
          write(6,*)'id=',idrank,'n_pe2recv_fluid_sing=', &
           n_pe2recv_fluid_sing,i_pe2recv_fluid_sing(0:n_pe2recv_fluid_sing-1)
          write(6,*)'id=',idrank,'maxreceiver=',maxreceiver
        endif
        call flush(6)
        call get_sync_world
      enddo
    endif
    
    
    !synopsis:
    !n_pe2send_fluid = towards how many I have to send
    !n_pe2recv_fluid = from how many I have to receive
    !i_pe2send_fluid = list of the nodes/PEs from whom I have to receive
    !i_pe2recv_fluid = list of the nodes/PEs from whom I have to send
    
    ! construct the following lists:
    !
    ! i_pe2send_fluid : list of PEs to send fluid info
    ! i_pe2recv_fluid : list of PEs to receive fluid info
    ! i_pop2send      : list of populations/sites to send fluid info
    ! i_pop2recv      : list of populations/sites to receive fluid info
    ltemp(1)=.false.

    do k=minz-nsing,maxz+nsing
      do jy=miny-nsing,maxy+nsing
        do i=minx-nsing,maxx+nsing
           if(i<minx.or.jy<miny.or.k<minz.or.i>maxx.or.jy>maxy.or.k>maxz)then
                itemp=i
                jtemp=jy
                ktemp=k
                itemp2=i
                jtemp2=jy
                ktemp2=k
                !apply periodic conditions if necessary
                if(ixpbc.eq.1) then
                   if(itemp<1) then
                      itemp=itemp+nx
                   endif
                   if(itemp>nx) then
                      itemp=itemp-nx
                   endif
                endif
                if(iypbc.eq.1) then
                   if(jtemp<1) then
                      jtemp=jtemp+ny
                   endif
                   if(jtemp>ny) then
                      jtemp=jtemp-ny
                   endif
                endif
                if(izpbc.eq.1) then
                   if(ktemp<1) then
                      ktemp=ktemp+nz
                   endif
                   if(ktemp>nz) then
                      ktemp=ktemp-nz
                   endif
                endif
                i4=i4back(itemp,jtemp,ktemp)

                IF(ownern(i4).NE.idrank) THEN   ! select nodes with ownership /= idrank
                   j=receiver(ownern(i4))  !find the unique ID of the cast operation
                   if(i4back(itemp,jtemp,ktemp)<0)then
                     write(6,*)'id',idrank,' something bad ',itemp,jtemp,ktemp
                     ltemp(1)=.true.
                   endif
                   if(i4back(itemp2,jtemp2,ktemp2)<0)then
                     write(6,*)'id',idrank,' csomething bad ',itemp2,jtemp2,ktemp2
                     ltemp(1)=.true.
                   endif
                   i_sing2recv_fluid(0,n_sing2recv_fluid(j),j)= &
                    i4back(itemp,jtemp,ktemp) !unique ID of the sending node
                   i_sing2recv_fluid(1,n_sing2recv_fluid(j),j)= &
                    i4back(itemp2,jtemp2,ktemp2) !unique ID of the receiving node
                   n_sing2recv_fluid(j)=n_sing2recv_fluid(j)+1 !total number of pops which have to be sent for the ID cast OP
                ENDIF
             endif
          enddo
       enddo
    enddo
    call or_world_larr(ltemp,1)
    if(ltemp(1))then
      call finalize_world
      stop
    endif
    !send the number of pops which should be received
    !determine how much I have to received

    DO i=0, n_pe2send_fluid_sing-1
       CALL MPI_IRECV(n_sing2send_fluid(i),1,MPI_INTEGER, &
            & i_pe2send_fluid_sing(i),i_pe2send_fluid_sing(i)+findntag_sing, &
            & MPI_COMM_WORLD,request(i),ierr)
    ENDDO
    DO i=0, n_pe2recv_fluid_sing-1
       CALL MPI_SEND(n_sing2recv_fluid(i),1,MPI_INTEGER, &
            & i_pe2recv_fluid_sing(i),idrank+findntag_sing, &
            & MPI_COMM_WORLD,ierr)
    ENDDO
    IF(n_pe2send_fluid_sing.GT.0) THEN
       CALL MPI_WAITALL(n_pe2send_fluid_sing,request,status,ierr)
    ENDIF

    
    if(lverbose)then
      call flush(6)
      call get_sync_world
      do i=0,mxrank-1
        if(i==idrank)then
          write(6,*)'id=',idrank,'i_pe2send_fluid_sing=', &
           i_pe2send_fluid_sing(0:n_pe2send_fluid_sing-1)
          write(6,*)'id=',idrank,'n_sing2send_fluid=', &
           n_sing2send_fluid(0:n_pe2send_fluid_sing-1)
        endif
        call flush(6)
        call get_sync_world
      enddo
      call flush(6)
      call get_sync_world
      do i=0,mxrank-1
        if(i==idrank)then
          write(6,*)'id=',idrank,'i_pe2recv_fluid_sing=', &
           i_pe2recv_fluid_sing(0:n_pe2recv_fluid_sing-1)
          write(6,*)'id=',idrank,'n_sing2recv_fluid=', &
           n_sing2recv_fluid(0:n_pe2send_fluid_sing-1)
        endif
        call flush(6)
        call get_sync_world
      enddo
    endif
    
   
    
    !determine the maximum numper of pops which should be received
    maxsender=0
    DO i=0,n_pe2send_fluid_sing-1
       IF(n_sing2send_fluid(i).GT.maxsender) THEN
          maxsender=n_sing2send_fluid(i)
       ENDIF
    ENDDO
    !allocate to receive
    IF(maxsender.GT.0.AND.n_pe2send_fluid_sing.GT.0) THEN
       ALLOCATE(i_sing2send_fluid(0:1,0:maxsender-1,0:n_pe2send_fluid_sing-1))

       ALLOCATE(buffs_sing(0:maxsender-1,0:n_pe2send_fluid_sing-1))

    ENDIF
   
    DO i=0, n_pe2send_fluid_sing-1
       CALL MPI_IRECV(i_sing2send_fluid(0,0,i),2*n_sing2send_fluid(i),MYINT, &
            & i_pe2send_fluid_sing(i),i_pe2send_fluid_sing(i)+findntag_sing, &
            & MPI_COMM_WORLD,request(i),ierr)
    ENDDO
    DO i=0, n_pe2recv_fluid_sing-1
       CALL MPI_SEND(i_sing2recv_fluid(0,0,i),2*n_sing2recv_fluid(i),MYINT, &
            & i_pe2recv_fluid_sing(i),idrank+findntag_sing, &
            & MPI_COMM_WORLD,ierr)
    ENDDO
    IF(n_pe2recv_fluid_sing.GT.0) THEN
       CALL MPI_WAITALL(n_pe2send_fluid_sing,request,status,ierr)
    ENDIF
    
    
    call get_sync_world
    
    ltemp(1)=.false.
    DO i=0, n_pe2send_fluid_sing-1
       DO j=0,n_sing2send_fluid(i)-1
          i4_send=i_sing2send_fluid(0,j,i)
          i4_recv=i_sing2send_fluid(1,j,i)
          if(i4_send<1)then
            write(6,*)idrank,i,j,n_pe2send_fluid_sing,n_sing2send_fluid(i), &
             i_pe2send_fluid_sing(i)
          endif
          if(ownern(i4_send).NE.idrank)then
            write(6,*)'idrank = ',idrank,&
             'Horror in task hsing. I have to send ',i4_send, &
              ' belonging to ',ownern(i4_send)
            ltemp(1)=.true.
          endif
       ENDDO
    ENDDO
    call or_world_larr(ltemp,1)
    if(ltemp(1))then
      call finalize_world
      stop
    endif
    
    ltemp(1)=.false.
    DO i=0, n_pe2recv_fluid_sing-1
       DO j=0,n_sing2recv_fluid(i)-1
          i4_send=i_sing2recv_fluid(0,j,i)
          i4_recv=i_sing2recv_fluid(1,j,i)
          if(ownern(i4_send).EQ.idrank)then
            write(6,*)'idrank = ',idrank,&
             'Horror in task sing. I have to receiving ',i4_send, &
              ' belonging to me ',ownern(i4_send)
            ltemp(1)=.true.
          endif
        ENDDO
    ENDDO
    call or_world_larr(ltemp,1)
    if(ltemp(1))then
      call finalize_world
      stop
    endif
    
    

 END SUBROUTINE findneigh_single_halo
 
 SUBROUTINE findneigh_hvar_isfluid(nx,ny,nz,nbuff,ibctype,ixpbc,iypbc,izpbc, &
   minx,maxx,miny,maxy,minz,maxz)

!***********************************************************************
!     
!     LBsoft subroutine for creating the list of 
!     the HVAR exchanged by MPI sending and receiving
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Bernaschi
!     last modification October 2019
!     
!***********************************************************************

    IMPLICIT NONE

    INTEGER :: request(0:maxneigh-1)
    
    INTEGER :: status(MPI_STATUS_SIZE,0:maxneigh-1)
    
    
    INTEGER, intent(in) :: nx,ny,nz,nbuff,ibctype
    INTEGER, intent(inout) :: ixpbc,iypbc,izpbc,minx,maxx,miny,maxy, &
   minz,maxz
    
    INTEGER :: lminx, lmaxx, lminy, lmaxy, lminz, lmaxz
    INTEGER(kind=IPRC) :: i4,i4temp,i4shift
    INTEGER :: i,j,jy,k,ip,maxsender,maxreceiver
    INTEGER :: ishift,jshift,kshift
    INTEGER(kind=IPRC) :: nghb
    INTEGER :: n, l, m,ks,mym(3)
    INTEGER :: ifold,jfold,kfold
    INTEGER :: itemp,jtemp,ktemp
    INTEGER :: itemp2,jtemp2,ktemp2
    INTEGER :: firstminpe,firstmaxpe
    INTEGER :: i4_send,i4_recv

    INTEGER :: allneigh(0:maxneigh*mxrank-1)

    INTEGER :: sender(0:mxrank-1)
    INTEGER :: receiver(0:mxrank-1)
    
    logical, parameter :: lverbose=.false.
    logical :: ltemp(1)
    
    if(mxrank==1)return
    
    DO i=0,mxrank-1
       !sender(i)=0
       receiver(i)=0
    ENDDO
    DO i=0,maxneigh-1
       i_pe2send_fluid_hvar(i)=-1
       n_var2send_fluid(i)=-1
    ENDDO
    DO i=0,maxneigh-1
       i_pe2recv_fluid_hvar(i)=-1
       n_var2recv_fluid(i)=-1
    ENDDO
    n_pe2recv_fluid_hvar=0
    n_pe2send_fluid_hvar=0
    
    if(minx.eq.1) then
       if(ixpbc.eq.1) then
          lminx=1
       else
          lminx=0
       endif
    else
       lminx=minx
    endif
    if(miny.eq.1) then
       if(iypbc.eq.1) then
          lminy=1
       else
          lminy=0
       endif
    else
       lminy=miny
    endif
    if(minz.eq.1) then
       if(izpbc.eq.1) then
          lminz=1
       else
          lminz=0
       endif
    else
       lminz=minz
    endif
    if(maxx.eq.nx) then
       if(ixpbc.eq.1) then
          lmaxx=nx
       else
          lmaxx=nx+1
       endif
    else
       lmaxx=maxx
    endif
    if(maxy.eq.ny) then
       if(iypbc.eq.1) then
          lmaxy=ny
       else
          lmaxy=ny+1
       endif
    else
       lmaxy=maxy
    endif
    if(maxz.eq.nz) then
       if(izpbc.eq.1) then
          lmaxz=nz
       else
          lmaxz=nz+1
       endif
    else
       lmaxz=maxz
    endif
    
    
    if(lverbose)then
      call flush(6)
      call get_sync_world
      do i=0,mxrank-1
        if(i==idrank)then
          write(6,*)'id=',idrank,'lminx=',lminx,'lmaxx=',lmaxx,'lminy=',lminy, &
              'lmaxy=',lmaxy,'lminz=',lminz,'lmaxz=',lmaxz  
        endif
        call flush(6)
        call get_sync_world
      enddo
    endif
    
    ! receiver(i): count no. of neighboring fluid nodes that talk with i-th PE
    !determine what I have to receive from each node/PE around me
     ltemp(1)=.false.
     do k=minz-nbuff,maxz+nbuff
       do j=miny-nbuff,maxy+nbuff
         do i=minx-nbuff,maxx+nbuff
           if(i<minx.or.j<miny.or.k<minz.or.i>maxx.or.j>maxy.or.k>maxz)then
                itemp=i
                jtemp=j
                ktemp=k
                itemp2=i
                jtemp2=j
                ktemp2=k
                !apply periodic conditions if necessary
                if(ixpbc.eq.1) then
                   if(itemp<1) then
                      itemp=itemp+nx
                   endif
                   if(itemp>nx) then
                      itemp=itemp-nx
                   endif
                endif
                if(iypbc.eq.1) then
                   if(jtemp<1) then
                      jtemp=jtemp+ny
                   endif
                   if(jtemp>ny) then
                      jtemp=jtemp-ny
                   endif
                endif
                if(izpbc.eq.1) then
                   if(ktemp<1) then
                      ktemp=ktemp+nz
                   endif
                   if(ktemp>nz) then
                      ktemp=ktemp-nz
                   endif
                endif
                i4=i4back(itemp,jtemp,ktemp)
                if(i4<1)then
                  write(6,*)idrank,' something bad ',itemp2,jtemp2,ktemp2, &
                   itemp,jtemp,ktemp,minx,maxx
                  ltemp(1)=.true.
                endif
                if(ownern(i4)<0 .or.ownern(i4)>(mxrank-1))then
                  write(6,*)idrank,' something bad ',itemp2,jtemp2,ktemp2, &
                   itemp,jtemp,ktemp,i4,ownern(i4)
                  ltemp(1)=.true.
                endif
                IF(ownern(i4).NE.idrank) THEN
                   receiver(ownern(i4))=receiver(ownern(i4))+1
                ENDIF
             endif
          enddo
       enddo
    enddo
    call or_world_larr(ltemp,1)
    if(ltemp(1))then
      call finalize_world
      stop
    endif
    
    
    if(lverbose)then
      call flush(6)
      call get_sync_world
      do i=0,mxrank-1
        if(i==idrank)then
          write(6,*)'id=',idrank,'receiver=',receiver(0:mxrank-1)
        endif
        call flush(6)
        call get_sync_world
      enddo
    endif
    
    !determine the max number of fluid nodes for node/PE that I have to receive 
    j=0
    maxreceiver=0
    DO i=0,mxrank-1

       IF(receiver(i).GT.0) THEN ! select procs that talks to this one

          IF(j.GE.maxneigh) THEN
             PRINT *,'Task',idrank,'Invalid index for myreceiver',j
             STOP
          ENDIF

          IF(receiver(i).GT.maxreceiver) THEN
             maxreceiver=receiver(i)  ! get the max no. of nodes to talk to
          ENDIF

          i_pe2recv_fluid_hvar(j)=i    ! index=process # rank of the list of fluid nodes to talk to (that should be sent)
          n_var2recv_fluid(j)=0        ! no. of fluid nodes to talk to (initialization)
          receiver(i)=j                ! id of PE to talk in local ordering
          j=j+1                        ! local ordering index
       ELSE
          receiver(i)=-1
       ENDIF

    ENDDO
    n_pe2recv_fluid_hvar=j
    
    !allocate buffer
    IF(maxreceiver.GT.0.AND.n_pe2recv_fluid_hvar.GT.0) THEN
       ALLOCATE(i_var2recv_fluid(0:1,0:maxreceiver-1,0:n_pe2recv_fluid_hvar-1))

       ALLOCATE(buffr_hvar(0:maxreceiver-1,0:n_pe2recv_fluid_hvar-1))
       if(ldo_second)ALLOCATE(buffr_hvar2(0:maxreceiver-1,0:n_pe2recv_fluid_hvar-1))
       ALLOCATE(buffr_isfluid(0:maxreceiver-1,0:n_pe2recv_fluid_hvar-1))
    ENDIF
    
    !initialization
    DO i=0,maxneigh-1
       i_pe2send_fluid_hvar(i)=-1
    ENDDO

    CALL MPI_ALLGATHER(i_pe2recv_fluid_hvar,maxneigh,MPI_INTEGER,allneigh, &
     maxneigh,MPI_INTEGER, MPI_COMM_WORLD, ierr)

    !allneigh = i_pe2recv_fluid_hvar

    
    !determine to who I have to send
    j=0
    DO i=0,mxrank*maxneigh-1
       IF(allneigh(i)==idrank) THEN
          IF(j.GE.maxneigh) THEN
             PRINT *,'Task',idrank,'Invalid index for i_pe2send_fluid_hvar',j
             STOP
          ENDIF
          i_pe2send_fluid_hvar(j)=i/maxneigh
          j=j+1
       ENDIF
    ENDDO
    n_pe2send_fluid_hvar=j
    
    if(lverbose)then
      call flush(6)
      call get_sync_world
      do i=0,mxrank-1
        if(i==idrank)then
          write(6,*)'id=',idrank,'n_pe2send_fluid_hvar=', &
           n_pe2send_fluid_hvar,i_pe2send_fluid_hvar(0:n_pe2send_fluid_hvar-1)
          write(6,*)'id=',idrank,'n_pe2recv_fluid_hvar=', &
           n_pe2recv_fluid_hvar,i_pe2recv_fluid_hvar(0:n_pe2recv_fluid_hvar-1)
          write(6,*)'id=',idrank,'maxreceiver=',maxreceiver
        endif
        call flush(6)
        call get_sync_world
      enddo
    endif
    
    
    !synopsis:
    !n_pe2send_fluid = towards how many I have to send
    !n_pe2recv_fluid = from how many I have to receive
    !i_pe2send_fluid = list of the nodes/PEs from whom I have to receive
    !i_pe2recv_fluid = list of the nodes/PEs from whom I have to send
    
    ! construct the following lists:
    !
    ! i_pe2send_fluid : list of PEs to send fluid info
    ! i_pe2recv_fluid : list of PEs to receive fluid info
    ! i_pop2send      : list of populations/sites to send fluid info
    ! i_pop2recv      : list of populations/sites to receive fluid info
    ltemp(1)=.false.

    do k=minz-nbuff,maxz+nbuff
      do jy=miny-nbuff,maxy+nbuff
        do i=minx-nbuff,maxx+nbuff
           if(i<minx.or.jy<miny.or.k<minz.or.i>maxx.or.jy>maxy.or.k>maxz)then
                itemp=i
                jtemp=jy
                ktemp=k
                itemp2=i
                jtemp2=jy
                ktemp2=k
                !apply periodic conditions if necessary
                if(ixpbc.eq.1) then
                   if(itemp<1) then
                      itemp=itemp+nx
                   endif
                   if(itemp>nx) then
                      itemp=itemp-nx
                   endif
                endif
                if(iypbc.eq.1) then
                   if(jtemp<1) then
                      jtemp=jtemp+ny
                   endif
                   if(jtemp>ny) then
                      jtemp=jtemp-ny
                   endif
                endif
                if(izpbc.eq.1) then
                   if(ktemp<1) then
                      ktemp=ktemp+nz
                   endif
                   if(ktemp>nz) then
                      ktemp=ktemp-nz
                   endif
                endif
                i4=i4back(itemp,jtemp,ktemp)

                IF(ownern(i4).NE.idrank) THEN   ! select nodes with ownership /= idrank
                   j=receiver(ownern(i4))  !find the unique ID of the cast operation
                   if(i4back(itemp,jtemp,ktemp)<0)then
                     write(6,*)'id',idrank,' something bad ',itemp,jtemp,ktemp
                     ltemp(1)=.true.
                   endif
                   if(i4back(itemp2,jtemp2,ktemp2)<0)then
                     write(6,*)'id',idrank,' csomething bad ',itemp2,jtemp2,ktemp2
                     ltemp(1)=.true.
                   endif
                   i_var2recv_fluid(0,n_var2recv_fluid(j),j)= &
                    i4back(itemp,jtemp,ktemp) !unique ID of the sending node
                   i_var2recv_fluid(1,n_var2recv_fluid(j),j)= &
                    i4back(itemp2,jtemp2,ktemp2) !unique ID of the receiving node
                   n_var2recv_fluid(j)=n_var2recv_fluid(j)+1 !total number of pops which have to be sent for the ID cast OP
                ENDIF
             endif
          enddo
       enddo
    enddo
    call or_world_larr(ltemp,1)
    if(ltemp(1))then
      call finalize_world
      stop
    endif
    !send the number of pops which should be received
    !determine how much I have to received

    DO i=0, n_pe2send_fluid_hvar-1
       CALL MPI_IRECV(n_var2send_fluid(i),1,MPI_INTEGER, &
            & i_pe2send_fluid_hvar(i),i_pe2send_fluid_hvar(i)+findntag_hvar, &
            & MPI_COMM_WORLD,request(i),ierr)
    ENDDO
    DO i=0, n_pe2recv_fluid_hvar-1
       CALL MPI_SEND(n_var2recv_fluid(i),1,MPI_INTEGER, &
            & i_pe2recv_fluid_hvar(i),idrank+findntag_hvar, &
            & MPI_COMM_WORLD,ierr)
    ENDDO
    IF(n_pe2send_fluid_hvar.GT.0) THEN
       CALL MPI_WAITALL(n_pe2send_fluid_hvar,request,status,ierr)
    ENDIF

    
    if(lverbose)then
      call flush(6)
      call get_sync_world
      do i=0,mxrank-1
        if(i==idrank)then
          write(6,*)'id=',idrank,'i_pe2send_fluid_hvar=', &
           i_pe2send_fluid_hvar(0:n_pe2send_fluid_hvar-1)
          write(6,*)'id=',idrank,'n_var2send_fluid=', &
           n_var2send_fluid(0:n_pe2send_fluid_hvar-1)
        endif
        call flush(6)
        call get_sync_world
      enddo
      call flush(6)
      call get_sync_world
      do i=0,mxrank-1
        if(i==idrank)then
          write(6,*)'id=',idrank,'i_pe2recv_fluid_hvar=', &
           i_pe2recv_fluid_hvar(0:n_pe2recv_fluid_hvar-1)
          write(6,*)'id=',idrank,'n_var2recv_fluid=', &
           n_var2recv_fluid(0:n_pe2send_fluid_hvar-1)
        endif
        call flush(6)
        call get_sync_world
      enddo
    endif
    
   
    
    !determine the maximum numper of pops which should be received
    maxsender=0
    DO i=0,n_pe2send_fluid_hvar-1
       IF(n_var2send_fluid(i).GT.maxsender) THEN
          maxsender=n_var2send_fluid(i)
       ENDIF
    ENDDO
    !allocate to receive
    IF(maxsender.GT.0.AND.n_pe2send_fluid_hvar.GT.0) THEN
       ALLOCATE(i_var2send_fluid(0:1,0:maxsender-1,0:n_pe2send_fluid_hvar-1))

       ALLOCATE(buffs_hvar(0:maxsender-1,0:n_pe2send_fluid_hvar-1))
       if(ldo_second)ALLOCATE(buffs_hvar2(0:maxsender-1,0:n_pe2send_fluid_hvar-1))
       ALLOCATE(buffs_isfluid(0:maxsender-1,0:n_pe2send_fluid_hvar-1))

    ENDIF
    
    DO i=0, n_pe2send_fluid_hvar-1
       CALL MPI_IRECV(i_var2send_fluid(0,0,i),2*n_var2send_fluid(i),MYINT, &
            & i_pe2send_fluid_hvar(i),i_pe2send_fluid_hvar(i)+findntag_hvar, &
            & MPI_COMM_WORLD,request(i),ierr)
    ENDDO
    DO i=0, n_pe2recv_fluid_hvar-1
       CALL MPI_SEND(i_var2recv_fluid(0,0,i),2*n_var2recv_fluid(i),MYINT, &
            & i_pe2recv_fluid_hvar(i),idrank+findntag_hvar, &
            & MPI_COMM_WORLD,ierr)
    ENDDO
    IF(n_pe2recv_fluid_hvar.GT.0) THEN
       CALL MPI_WAITALL(n_pe2send_fluid_hvar,request,status,ierr)
    ENDIF
    
    
    call get_sync_world
    
    ltemp(1)=.false.
    DO i=0, n_pe2send_fluid_hvar-1
       DO j=0,n_var2send_fluid(i)-1
          i4_send=i_var2send_fluid(0,j,i)
          i4_recv=i_var2send_fluid(1,j,i)
          if(i4_send<1)then
            write(6,*)idrank,i,j,n_pe2send_fluid_hvar,n_var2send_fluid(i),i_pe2send_fluid_hvar(i)
          endif
          if(ownern(i4_send).NE.idrank)then
            write(6,*)'idrank = ',idrank,&
             'Horror in task hvar. I have to send ',i4_send, &
              ' belonging to ',ownern(i4_send)
            ltemp(1)=.true.
          endif
       ENDDO
    ENDDO
    call or_world_larr(ltemp,1)
    if(ltemp(1))then
      call finalize_world
      stop
    endif
    
    ltemp(1)=.false.
    DO i=0, n_pe2recv_fluid_hvar-1
       DO j=0,n_var2recv_fluid(i)-1
          i4_send=i_var2recv_fluid(0,j,i)
          i4_recv=i_var2recv_fluid(1,j,i)
          if(ownern(i4_send).EQ.idrank)then
            write(6,*)'idrank = ',idrank,&
             'Horror in task hvar. I have to receiving ',i4_send, &
              ' belonging to me ',ownern(i4_send)
            ltemp(1)=.true.
          endif
        ENDDO
    ENDDO
    call or_world_larr(ltemp,1)
    if(ltemp(1))then
      call finalize_world
      stop
    endif
    
    

 END SUBROUTINE findneigh_hvar_isfluid
 
 SUBROUTINE findneigh(nx,ny,nz,nbuff,ibctype,isfluid,ixpbc,iypbc,izpbc,minx,maxx, &
       miny,maxy,minz,maxz)

!***********************************************************************
!     
!     LBsoft subroutine for creating the list of 
!     the populations exchanged by MPI sending and receiving
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Bernaschi
!     last modification October 2019
!     
!***********************************************************************
 
    IMPLICIT NONE

    INTEGER :: request(0:maxneigh-1)
    
    INTEGER :: status(MPI_STATUS_SIZE,0:maxneigh-1)
    

    INTEGER, intent(in) :: nx,ny,nz,nbuff,ibctype
    INTEGER(kind=1), intent(in), allocatable, dimension(:,:,:) :: isfluid
    INTEGER, intent(inout) :: ixpbc,iypbc,izpbc,minx,maxx,miny,maxy, &
         minz,maxz

    INTEGER :: lminx, lmaxx, lminy, lmaxy, lminz, lmaxz
    INTEGER(kind=IPRC) :: i4,i4temp,i4shift
    INTEGER :: i,j,jy,k,ip,maxsender,maxreceiver
    INTEGER :: ishift,jshift,kshift
    INTEGER(kind=IPRC) :: nghb
    INTEGER :: n, l, m,ks,mym(3)
    INTEGER :: ifold,jfold,kfold
    INTEGER :: itemp,jtemp,ktemp
    INTEGER :: itemp2,jtemp2,ktemp2
    INTEGER :: firstminpe,firstmaxpe

    INTEGER :: allneigh(0:maxneigh*mxrank-1)

    INTEGER :: sender(0:mxrank-1)
    
    logical, parameter :: lverbose=.false.
    logical :: ltemp(1)
    
    if(mxrank==1)return

    DO i=0,mxrank-1
       sender(i)=0
    ENDDO
    DO i=0,maxneigh-1
       i_pe2send_fluid(i)=-1
       n_pop2send_fluid(i)=-1
    ENDDO
    n_pe2recv_fluid=0
    n_pe2send_fluid=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(minx.lt.1) then
       if(ixpbc.eq.1) then
          lminx=1
       else
          lminx=0
       endif
    else
       lminx=minx
    endif
    if(miny.lt.1) then
       if(iypbc.eq.1) then
          lminy=1
       else
          lminy=0
       endif
    else
       lminy=miny
    endif
    if(minz.lt.1) then
       if(izpbc.eq.1) then
          lminz=1
       else
          lminz=0
       endif
    else
       lminz=minz
    endif
    if(maxx.gt.nx) then
       if(ixpbc.eq.1) then
          lmaxx=nx
       else
          lmaxx=nx+1
       endif
    else
       lmaxx=maxx
    endif
    if(maxy.gt.ny) then
       if(iypbc.eq.1) then
          lmaxy=ny
       else
          lmaxy=ny+1
       endif
    else
       lmaxy=maxy
    endif
    if(maxz.gt.nz) then
       if(izpbc.eq.1) then
          lmaxz=nz
       else
          lmaxz=nz+1
       endif
    else
       lmaxz=maxz
    endif
    
    
    
    if(lverbose)write(6,*)'id=',idrank,'lminx=',lminx,'lmaxx=',lmaxx, &
     'lminy=',lminy,'lmaxy=',lmaxy,'lminz=',lminz,'lmaxz=',lmaxz  
    
    ! sender(i): count no. of neighboring fluid nodes that talk with i-th PE
    !determine what I have to send to each node/PE around me
    do l=1,npop-1
       ishift=ex(l)
       jshift=ey(l)
       kshift=ez(l)
       do k=lminz,lmaxz
          do j=lminy,lmaxy
             do i=lminx,lmaxx
                
                itemp=i+ishift;
                jtemp=j+jshift;
                ktemp=k+kshift;
                itemp2=i+ishift;
                jtemp2=j+jshift;
                ktemp2=k+kshift;
                !apply periodic conditions if necessary
                if(ixpbc.eq.1) then
                   if(itemp.eq.0) then
                      itemp=nx
                   endif
                   if(itemp.eq.(nx+1)) then
                      itemp=1
                   endif
                endif
                if(iypbc.eq.1) then
                   if(jtemp.eq.0) then
                      jtemp=ny
                   endif
                   if(jtemp.eq.(ny+1)) then
                      jtemp=1
                   endif
                endif
                if(izpbc.eq.1) then
                   if(ktemp.eq.0) then
                      ktemp=nz
                   endif
                   if(ktemp.eq.(nz+1)) then
                      ktemp=1
                   endif
                endif
                i4=i4back(itemp,jtemp,ktemp)
                !max                if(idrank.eq.0) then
                !max                   write(50,*)'i',itemp,'j',jtemp,'k',ktemp,'owner',ownern(i4)
                !max                endif
                !max                if(idrank.eq.1) then
                !max                   write(51,*)'i',itemp,'j',jtemp,'k',ktemp,'owner',ownern(i4)
                !max                endif   
                IF(ownern(i4).NE.idrank .AND. isfluid(i,j,k)/=3 .AND. &
                 isfluid(i,j,k)/=0) THEN
                   sender(ownern(i4))=sender(ownern(i4))+1
                   !max                   write(0,*)'ownern ',ownern(i4),'for node ',i4,'ix=',i,'iy=',j,'iz=',k,'l=',l
                   !max                   write(0,*)'itemp=',itemp,'jtemp=',jtemp,'ktemp=',ktemp
                ENDIF
             enddo
          enddo
       enddo
    enddo
    if(lverbose)write(6,*)'id=',idrank,'sender=',sender(0:1)
    !determine the max number of fluid nodes for node/PE to whom I have to communicate 
    ltemp(1)=.false.
    j=0
    maxsender=0
    DO i=0,mxrank-1

       IF(sender(i).GT.0) THEN ! select procs that talks to this one

          IF(j.GE.maxneigh) THEN
             PRINT *,'Task',idrank,'Invalid index for myreceiver',j
             ltemp(1)=.true.
          ENDIF

          IF(sender(i).GT.maxsender) THEN
             maxsender=sender(i)  ! get the max no. of nodes to talk to
          ENDIF

          i_pe2send_fluid(j)=i    ! index=process # rank of the list of fluid nodes to talk to (that should be sent)
          n_pop2send_fluid(j)=0   ! no. of fluid nodes to talk to  (initialization)
          sender(i)=j             ! id of PE to talk in local ordering
          j=j+1                   ! local ordering index
       ELSE
          sender(i)=-1
       ENDIF

    ENDDO
    n_pe2send_fluid=j
    call or_world_larr(ltemp,1)
    if(ltemp(1))then
      call finalize_world
      stop
    endif
    
    !allocate buffer
    IF(maxsender.GT.0.AND.n_pe2send_fluid.GT.0) THEN
       ALLOCATE(i_pop2send_fluid(0:1,0:maxsender-1,0:n_pe2send_fluid-1))

       ALLOCATE(buffpops(0:maxsender-1,0:n_pe2send_fluid-1))

    ENDIF
    !max    write(0,*)'id=',idrank,'maxsender=',maxsender,'n_pe2send_fluid=',n_pe2send_fluid
    DO i=0,maxneigh-1
       i_pe2recv_fluid(i)=-1
    ENDDO
    
    CALL MPI_ALLGATHER(i_pe2send_fluid,maxneigh,MPI_INTEGER,allneigh, &
     maxneigh,MPI_INTEGER, MPI_COMM_WORLD, ierr)

    !allneigh = i_pe2send_fluid

    !determine from who I have to receive
    ltemp(1)=.false.
    j=0
    DO i=0,mxrank*maxneigh-1
       IF(allneigh(i)==idrank) THEN
          IF(j.GE.maxneigh) THEN
             PRINT *,'Task',idrank,'Invalid index for i_pe2recv_fluid',j
             ltemp(1)=.true.
          ENDIF
          i_pe2recv_fluid(j)=i/maxneigh
          j=j+1
       ENDIF
    ENDDO
    n_pe2recv_fluid=j
    call or_world_larr(ltemp,1)
    if(ltemp(1))then
      call finalize_world
      stop
    endif
    !synopsis:
    !n_pe2send_fluid = towards how many I have to send
    !n_pe2recv_fluid = from how many I have to receive
    !i_pe2send_fluid = list of the nodes/PEs from whom I have to receive
    !i_pe2recv_fluid = list of the nodes/PEs from whom I have to send
    
    if(lverbose)then
      write(6,*)'id=',idrank,'n_pe2send_fluid=',n_pe2send_fluid, &
       i_pe2send_fluid(0:n_pe2send_fluid-1)
      write(6,*)'id=',idrank,'n_pe2recv_fluid=',n_pe2recv_fluid, &
       i_pe2recv_fluid(0:n_pe2recv_fluid-1)
    endif
    

    ! construct the following lists:
    !
    ! i_pe2send_fluid : list of PEs to send fluid info
    ! i_pe2recv_fluid : list of PEs to receive fluid info
    ! i_pop2send      : list of populations/sites to send fluid info
    ! i_pop2recv      : list of populations/sites to receive fluid info

    do l=1,npop-1
       ishift=ex(l)
       jshift=ey(l)
       kshift=ez(l)
       do k=minz,maxz
          do jy=miny,maxy
             do i=minx,maxx
                itemp=i+ishift;
                jtemp=jy+jshift;
                ktemp=k+kshift;
                itemp2=i+ishift;
                jtemp2=jy+jshift;
                ktemp2=k+kshift;
                !apply periodic conditions if necessary
                if(ixpbc.eq.1) then
                   if(itemp.eq.0) then
                      itemp=nx
                   endif
                   if(itemp.eq.(nx+1)) then
                      itemp=1
                   endif
                endif
                if(iypbc.eq.1) then
                   if(jtemp.eq.0) then
                      jtemp=ny
                   endif
                   if(jtemp.eq.(ny+1)) then
                      jtemp=1
                   endif
                endif
                if(izpbc.eq.1) then
                   if(ktemp.eq.0) then
                      ktemp=nz
                   endif
                   if(ktemp.eq.(nz+1)) then
                      ktemp=1
                   endif
                endif
                i4=i4back(itemp,jtemp,ktemp)

                IF(ownern(i4).NE.idrank .AND. isfluid(i,jy,k)/=3 .AND. &
                 isfluid(i,jy,k)/=0) THEN   ! select nodes with ownership /= idrank
                   j=sender(ownern(i4))  !find the unique ID of the cast operation
                   i_pop2send_fluid(0,n_pop2send_fluid(j),j)=l !direction of the population
                   i_pop2send_fluid(1,n_pop2send_fluid(j),j)=i4back(i,jy,k) !unique ID of the node
                   n_pop2send_fluid(j)=n_pop2send_fluid(j)+1 !total number of pops which have to be sent for the ID cast OP
                ENDIF
             enddo
          enddo
       enddo
    enddo
    !send the number of pops which should be received
    !determine how much I have to received
    
    DO i=0, n_pe2recv_fluid-1
       CALL MPI_IRECV(n_pop2recv_fluid(i),1,MPI_INTEGER, &
            & i_pe2recv_fluid(i),i_pe2recv_fluid(i)+findntag, &
            & MPI_COMM_WORLD,request(i),ierr)
    ENDDO
    DO i=0, n_pe2send_fluid-1
       CALL MPI_SEND(n_pop2send_fluid(i),1,MPI_INTEGER, &
            & i_pe2send_fluid(i),idrank+findntag, &
            & MPI_COMM_WORLD,ierr)
    ENDDO
    IF(n_pe2recv_fluid.GT.0) THEN
       CALL MPI_WAITALL(n_pe2recv_fluid,request,status,ierr)
    ENDIF
    
    !determine the maximum numper of pops which should be received
    maxreceiver=0
    DO i=0,n_pe2recv_fluid-1
       IF(n_pop2recv_fluid(i).GT.maxreceiver) THEN
          maxreceiver=n_pop2recv_fluid(i)
       ENDIF
    ENDDO
    !allocate to receive
    IF(maxreceiver.GT.0.AND.n_pe2recv_fluid.GT.0) THEN
       ALLOCATE(i_pop2recv_fluid(0:1,0:maxreceiver-1,0:n_pe2recv_fluid-1))

       ALLOCATE(buffpopr(0:maxreceiver-1,0:n_pe2recv_fluid-1))

    ENDIF
    
    DO i=0, n_pe2recv_fluid-1
       CALL MPI_IRECV(i_pop2recv_fluid(0,0,i),2*n_pop2recv_fluid(i),MYINT, &
            & i_pe2recv_fluid(i),i_pe2recv_fluid(i)+findntag, &
            & MPI_COMM_WORLD,request(i),ierr)
    ENDDO
    DO i=0, n_pe2send_fluid-1
       CALL MPI_SEND(i_pop2send_fluid(0,0,i),2*n_pop2send_fluid(i),MYINT, &
            & i_pe2send_fluid(i),idrank+findntag, &
            & MPI_COMM_WORLD,ierr)
    ENDDO
    IF(n_pe2recv_fluid.GT.0) THEN
       CALL MPI_WAITALL(n_pe2recv_fluid,request,status,ierr)
    ENDIF
    
    ltemp(1)=.false.
    DO i=0, n_pe2recv_fluid-1
       DO j=0,n_pop2recv_fluid(i)-1
          ishift=ex(i_pop2recv_fluid(0,j,i))
          jshift=ey(i_pop2recv_fluid(0,j,i))
          kshift=ez(i_pop2recv_fluid(0,j,i))
          i4=i_pop2recv_fluid(1,j,i)
          !ktemp = (i4/INT(nxy2,KIND=IPRC))
          !jtemp = ((i4 - ktemp*INT(nxy2,KIND=IPRC))/INT(nx2,KIND=IPRC))
          !itemp = (i4 - ktemp*INT(nxy2,KIND=IPRC) - jtemp*INT(nx2,KIND=IPRC))
          mym=i4find(i4)
          ktemp = mym(3)
          jtemp = mym(2)
          itemp = mym(1)
          ktemp = ktemp+kshift
          jtemp = jtemp+jshift
          itemp = itemp+ishift
          !apply periodic conditions if necessary
          if(ixpbc.eq.1) then
             if(itemp.eq.0) then
                itemp=nx
             endif
             if(itemp.eq.(nx+1)) then
                itemp=1
             endif
          endif
          if(iypbc.eq.1) then
             if(jtemp.eq.0) then
                jtemp=ny
             endif
             if(jtemp.eq.(ny+1)) then
                jtemp=1
             endif
          endif
          if(izpbc.eq.1) then
             if(ktemp.eq.0) then
                ktemp=nz
             endif
             if(ktemp.eq.(nz+1)) then
                ktemp=1
             endif
          endif
          i4temp=i4back(itemp,jtemp,ktemp)
          i_pop2recv_fluid(1,j,i)=i4temp
          if(ownern(i4temp).ne.idrank) then
             write(50+idrank,*)'Horror in task',idrank,i4,i4temp,itemp,jtemp,&
                  &                        ktemp,ishift,jshift,kshift
             ltemp(1)=.true.
          endif
       ENDDO
    ENDDO
    call or_world_larr(ltemp,1)
    if(ltemp(1))then
      call finalize_world
      stop
    endif


  END SUBROUTINE findneigh
 
 subroutine comm_init_isfluid(temp)

!***********************************************************************
!     
!     LBsoft subroutine to exchange ISFLUID
!     amond the MPI ranks
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
  
  IMPLICIT NONE

  INTEGER(kind=1), dimension(:,:,:), allocatable :: temp
  
  INTEGER :: request_sub(0:maxneigh-1)
  INTEGER :: requestm_sub(0:maxneigh-1)
  INTEGER :: status(MPI_STATUS_SIZE,0:maxneigh-1)
  INTEGER :: i,j,mym(3)
  INTEGER(kind=IPRC) :: i4, itemp, jtemp, ktemp

  
  if(mxrank==1)return
  
  do i=0, n_pe2recv_fluid_hvar-1
    CALL MPI_IRECV(buffr_isfluid(0,i),n_var2recv_fluid(i),MYINT1, &
      &  i_pe2recv_fluid_hvar(i),i_pe2recv_fluid_hvar(i)+tag_isfluid, &
        MPI_COMM_WORLD,requestm_sub(i),ierr)
  enddo
  
  
   do i=0, n_pe2send_fluid_hvar-1
      do j=0,n_var2send_fluid(i)-1
         i4=i_var2send_fluid(0,j,i)
         mym=i4find(i4)
         ktemp = mym(3)
         jtemp = mym(2)
         itemp = mym(1)
         buffs_isfluid(j,i)=temp(itemp,jtemp,ktemp)
      enddo

      CALL MPI_ISEND(buffs_isfluid(0,i),n_var2send_fluid(i),MYINT1, &
       & i_pe2send_fluid_hvar(i),idrank+tag_isfluid,MPI_COMM_WORLD,request_sub(i),ierr)
     

   enddo

   if(n_pe2recv_fluid_hvar.gt.0) then
      CALL MPI_WAITALL(n_pe2recv_fluid_hvar,requestm_sub,status,ierr)
   endif
   do i=0, n_pe2recv_fluid_hvar-1
      do j=0,n_var2recv_fluid(i)-1
         i4=i_var2recv_fluid(1,j,i)
         mym=i4find(i4)
         ktemp = mym(3)
         jtemp = mym(2)
         itemp = mym(1)
         temp(itemp,jtemp,ktemp)=buffr_isfluid(j,i)
      enddo
   enddo
   
   if(n_pe2send_fluid_hvar.gt.0) then
      call MPI_WAITALL(n_pe2send_fluid_hvar,request_sub,status,ierr)
   endif
   
 end subroutine comm_init_isfluid
 
 subroutine commexch_isfluid(temp)

!***********************************************************************
!     
!     LBsoft subroutine for exchanging ISFLUID variables
!     among MPI ranks
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  INTEGER(kind=1), dimension(:,:,:), allocatable :: temp
  
  INTEGER :: status(MPI_STATUS_SIZE,0:maxneigh-1)
  INTEGER :: i,j,mym(3)
  INTEGER(kind=IPRC) :: i4, itemp, jtemp, ktemp

  
  if(mxrank==1)return
  
  do i=0, n_pe2recv_fluid_hvar-1
    CALL MPI_IRECV(buffr_isfluid(0,i),n_var2recv_fluid(i),MYINT1, &
      &  i_pe2recv_fluid_hvar(i),i_pe2recv_fluid_hvar(i)+tag_isfluid, &
        MPI_COMM_WORLD,requestm_isfluid(i),ierr)
  enddo
  
  
   do i=0, n_pe2send_fluid_hvar-1
      do j=0,n_var2send_fluid(i)-1
         i4=i_var2send_fluid(0,j,i)
         mym=i4find(i4)
         ktemp = mym(3)
         jtemp = mym(2)
         itemp = mym(1)
         buffs_isfluid(j,i)=temp(itemp,jtemp,ktemp)
      enddo

      CALL MPI_ISEND(buffs_isfluid(0,i),n_var2send_fluid(i),MYINT1, &
       & i_pe2send_fluid_hvar(i),idrank+tag_isfluid,MPI_COMM_WORLD,request_isfluid(i),ierr)
     

   enddo

 
   
 end subroutine commexch_isfluid
   
 subroutine commwait_isfluid(temp)

!***********************************************************************
!     
!     LBsoft subroutine for waiting the exchange ISFLUID 
!     variables among MPI ranks
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
   
  INTEGER(kind=1), dimension(:,:,:), allocatable :: temp
  
  INTEGER :: status(MPI_STATUS_SIZE,0:maxneigh-1)
  INTEGER :: i,j,mym(3)
  INTEGER(kind=IPRC) :: i4, itemp, jtemp, ktemp
  
  if(mxrank==1)return

   if(n_pe2recv_fluid_hvar.gt.0) then
      CALL MPI_WAITALL(n_pe2recv_fluid_hvar,requestm_isfluid,status,ierr)
   endif
   do i=0, n_pe2recv_fluid_hvar-1
      do j=0,n_var2recv_fluid(i)-1
         i4=i_var2recv_fluid(1,j,i)
         mym=i4find(i4)
         ktemp = mym(3)
         jtemp = mym(2)
         itemp = mym(1)
         temp(itemp,jtemp,ktemp)=buffr_isfluid(j,i)
      enddo
   enddo
   
   if(n_pe2send_fluid_hvar.gt.0) then
      CALL MPI_WAITALL(n_pe2send_fluid_hvar,request_isfluid,status,ierr)
   endif
   

 end subroutine commwait_isfluid
 
 subroutine commexch_dens(dtemp,dtemp2)

!***********************************************************************
!     
!     LBsoft subroutine for exchanging density variables
!     among MPI ranks
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************

  IMPLICIT NONE

  REAL(KIND=PRC), dimension(:,:,:), allocatable :: dtemp,dtemp2
  INTEGER :: status(MPI_STATUS_SIZE,0:maxneigh-1)
  INTEGER :: i,j,mym(3)
  INTEGER(kind=IPRC) :: i4, itemp, jtemp, ktemp

  if(mxrank==1)return
  
  do i=0, n_pe2recv_fluid_hvar-1
    CALL MPI_IRECV(buffr_hvar(0,i),n_var2recv_fluid(i),MYFLOAT, &
      &  i_pe2recv_fluid_hvar(i),i_pe2recv_fluid_hvar(i)+tag_hvar, &
        MPI_COMM_WORLD,requestm_hvar(i),ierr)
  enddo
  
  


   do i=0, n_pe2send_fluid_hvar-1
      do j=0,n_var2send_fluid(i)-1
         i4=i_var2send_fluid(0,j,i)
         mym=i4find(i4)
         ktemp = mym(3)
         jtemp = mym(2)
         itemp = mym(1)
         buffs_hvar(j,i)=dtemp(itemp,jtemp,ktemp)
      enddo

      CALL MPI_ISEND(buffs_hvar(0,i),n_var2send_fluid(i),MYFLOAT, &
       & i_pe2send_fluid_hvar(i),idrank+tag_hvar,MPI_COMM_WORLD,request_hvar(i),ierr)
     

   enddo
   
   if(ldo_second)then
    do i=0, n_pe2recv_fluid_hvar-1
      CALL MPI_IRECV(buffr_hvar2(0,i),n_var2recv_fluid(i),MYFLOAT, &
        &  i_pe2recv_fluid_hvar(i),i_pe2recv_fluid_hvar(i)+tag_hvar2, &
          MPI_COMM_WORLD,requestm_hvar2(i),ierr)
    enddo
    do i=0, n_pe2send_fluid_hvar-1
      do j=0,n_var2send_fluid(i)-1
         i4=i_var2send_fluid(0,j,i)
         mym=i4find(i4)
         ktemp = mym(3)
         jtemp = mym(2)
         itemp = mym(1)
         buffs_hvar2(j,i)=dtemp2(itemp,jtemp,ktemp)
         
      enddo

      CALL MPI_ISEND(buffs_hvar2(0,i),n_var2send_fluid(i),MYFLOAT, &
       & i_pe2send_fluid_hvar(i),idrank+tag_hvar2,MPI_COMM_WORLD,request_hvar2(i),ierr)
   enddo
  endif
 
  return
 
 end subroutine commexch_dens
 
 subroutine commwait_dens(dtemp,dtemp2)
 
!***********************************************************************
!     
!     LBsoft subroutine for waiting the exchange densities 
!     variables among MPI ranks
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************

  IMPLICIT NONE

  REAL(KIND=PRC), dimension(:,:,:), allocatable :: dtemp,dtemp2
  INTEGER :: request(0:maxneigh-1)
  INTEGER :: status(MPI_STATUS_SIZE,0:maxneigh-1)
  INTEGER :: i,j,mym(3)
  INTEGER(kind=IPRC) :: i4, itemp, jtemp, ktemp
  
  if(mxrank==1)return

   if(n_pe2recv_fluid_hvar.gt.0) then
      CALL MPI_WAITALL(n_pe2recv_fluid_hvar,requestm_hvar,status,ierr)
   endif
   do i=0, n_pe2recv_fluid_hvar-1
      do j=0,n_var2recv_fluid(i)-1
         i4=i_var2recv_fluid(1,j,i)
         mym=i4find(i4)
         ktemp = mym(3)
         jtemp = mym(2)
         itemp = mym(1)
         dtemp(itemp,jtemp,ktemp)=buffr_hvar(j,i)
      enddo
   enddo
   
   if(n_pe2send_fluid_hvar.gt.0) then
      CALL MPI_WAITALL(n_pe2send_fluid_hvar,request_hvar,status,ierr)
   endif
   
   if(ldo_second)then
     if(n_pe2recv_fluid_hvar.gt.0) then
       CALL MPI_WAITALL(n_pe2recv_fluid_hvar,requestm_hvar2,status,ierr)
     endif
     do i=0, n_pe2recv_fluid_hvar-1
       do j=0,n_var2recv_fluid(i)-1
         i4=i_var2recv_fluid(1,j,i)
         mym=i4find(i4)
         ktemp = mym(3)
         jtemp = mym(2)
         itemp = mym(1)
         dtemp2(itemp,jtemp,ktemp)=buffr_hvar2(j,i)
       enddo
     enddo
   
     if(n_pe2send_fluid_hvar.gt.0) then
       CALL MPI_WAITALL(n_pe2send_fluid_hvar,request_hvar2,status,ierr)
     endif
   endif

 end subroutine commwait_dens
 
 subroutine commexch_vel_component(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for exchanging velocities variables
!     among MPI ranks
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************

  IMPLICIT NONE

  REAL(KIND=PRC), dimension(:,:,:), allocatable :: dtemp
  INTEGER :: status(MPI_STATUS_SIZE,0:maxneigh-1)
  INTEGER :: i,j,mym(3)
  INTEGER(kind=IPRC) :: i4, itemp, jtemp, ktemp

  if(mxrank==1)return
  
  do i=0, n_pe2recv_fluid_hvar-1
    CALL MPI_IRECV(buffr_hvar(0,i),n_var2recv_fluid(i),MYFLOAT, &
      &  i_pe2recv_fluid_hvar(i),i_pe2recv_fluid_hvar(i)+tag_hvar, &
        MPI_COMM_WORLD,requestm_hvar(i),ierr)
  enddo
  
  
   do i=0, n_pe2send_fluid_hvar-1
      do j=0,n_var2send_fluid(i)-1
         i4=i_var2send_fluid(0,j,i)
         mym=i4find(i4)
         ktemp = mym(3)
         jtemp = mym(2)
         itemp = mym(1)
         buffs_hvar(j,i)=dtemp(itemp,jtemp,ktemp)
      enddo

      CALL MPI_ISEND(buffs_hvar(0,i),n_var2send_fluid(i),MYFLOAT, &
       & i_pe2send_fluid_hvar(i),idrank+tag_hvar,MPI_COMM_WORLD,request_hvar(i),ierr)
     

   enddo
   

 end subroutine commexch_vel_component
   
 subroutine commwait_vel_component(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for waiting the exchange velocities 
!     variables among MPI ranks
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************

  IMPLICIT NONE

  REAL(KIND=PRC), dimension(:,:,:), allocatable :: dtemp
  INTEGER :: request(0:maxneigh-1)
  INTEGER :: status(MPI_STATUS_SIZE,0:maxneigh-1)
  INTEGER :: i,j,mym(3)
  INTEGER(kind=IPRC) :: i4, itemp, jtemp, ktemp
  
  if(mxrank==1)return

   if(n_pe2recv_fluid_hvar.gt.0) then
      CALL MPI_WAITALL(n_pe2recv_fluid_hvar,requestm_hvar,status,ierr)
   endif
   do i=0, n_pe2recv_fluid_hvar-1
      do j=0,n_var2recv_fluid(i)-1
         i4=i_var2recv_fluid(1,j,i)
         mym=i4find(i4)
         ktemp = mym(3)
         jtemp = mym(2)
         itemp = mym(1)
         dtemp(itemp,jtemp,ktemp)=buffr_hvar(j,i)
      enddo
   enddo
   
   if(n_pe2send_fluid_hvar.gt.0) then
      CALL MPI_WAITALL(n_pe2send_fluid_hvar,request_hvar,status,ierr)
   endif
   

 end subroutine commwait_vel_component
 
 subroutine commexch_single_halo(dtemp,nbuff, &
  ix,fx,iy,fy,iz,fz)
 
!***********************************************************************
!     
!     LBsoft subroutine for exchanging velocities variables
!     among MPI ranks
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************

  IMPLICIT NONE
  
  integer, intent(in) :: nbuff,ix,fx,iy,fy,iz,fz
  real(kind=PRC), &
   dimension(ix-nbuff:fx+nbuff,iy-nbuff:fy+nbuff,iz-nbuff:fz+nbuff) :: &
   dtemp
  INTEGER :: status(MPI_STATUS_SIZE,0:maxneigh-1)
  INTEGER :: i,j,mym(3)
  INTEGER(kind=IPRC) :: i4, itemp, jtemp, ktemp

  if(mxrank==1)return
  
  do i=0, n_pe2recv_fluid_sing-1
    CALL MPI_IRECV(buffr_sing(0,i),n_sing2recv_fluid(i),MYFLOAT, &
      &  i_pe2recv_fluid_sing(i),i_pe2recv_fluid_sing(i)+tag_sing, &
        MPI_COMM_WORLD,requestm_sing(i),ierr)
  enddo
  
  
   do i=0, n_pe2send_fluid_sing-1
      do j=0,n_sing2send_fluid(i)-1
         i4=i_sing2send_fluid(0,j,i)
         mym=i4find(i4)
         ktemp = mym(3)
         jtemp = mym(2)
         itemp = mym(1)
         buffs_sing(j,i)=dtemp(itemp,jtemp,ktemp)
      enddo

      CALL MPI_ISEND(buffs_sing(0,i),n_sing2send_fluid(i),MYFLOAT, &
       & i_pe2send_fluid_sing(i),idrank+tag_sing,MPI_COMM_WORLD,request_sing(i),ierr)
     

   enddo
   

 end subroutine commexch_single_halo
   
 subroutine commwait_single_halo(dtemp,nbuff,ix,fx,iy,fy,iz,fz)
 
!***********************************************************************
!     
!     LBsoft subroutine for waiting the exchange velocities 
!     variables among MPI ranks
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************

  IMPLICIT NONE

  integer, intent(in) :: nbuff,ix,fx,iy,fy,iz,fz
  real(kind=PRC), &
   dimension(ix-nbuff:fx+nbuff,iy-nbuff:fy+nbuff,iz-nbuff:fz+nbuff) :: &
   dtemp
  INTEGER :: request(0:maxneigh-1)
  INTEGER :: status(MPI_STATUS_SIZE,0:maxneigh-1)
  INTEGER :: i,j,mym(3)
  INTEGER(kind=IPRC) :: i4, itemp, jtemp, ktemp
  
  if(mxrank==1)return

   if(n_pe2recv_fluid_sing.gt.0) then
      CALL MPI_WAITALL(n_pe2recv_fluid_sing,requestm_sing,status,ierr)
   endif
   do i=0, n_pe2recv_fluid_sing-1
      do j=0,n_sing2recv_fluid(i)-1
         i4=i_sing2recv_fluid(1,j,i)
         mym=i4find(i4)
         ktemp = mym(3)
         jtemp = mym(2)
         itemp = mym(1)
         dtemp(itemp,jtemp,ktemp)=buffr_sing(j,i)
      enddo
   enddo
   
   if(n_pe2send_fluid_sing.gt.0) then
      CALL MPI_WAITALL(n_pe2send_fluid_sing,request_sing,status,ierr)
   endif
   

 end subroutine commwait_single_halo
 
 subroutine commspop(fluidsub)
 
!***********************************************************************
!     
!     LBsoft subroutine for sending the population over the MPI
!     network
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************

  implicit none

  real(kind=PRC), allocatable, dimension(:,:,:,:)  :: fluidsub
  integer :: request(0:maxneigh-1)
  integer :: istatus(MPI_STATUS_SIZE,0:maxneigh-1)
  integer :: i,j,mym(3)
  integer(kind=IPRC) :: i4, itemp, jtemp, ktemp
  
  if(mxrank==1)return
  
  if(firstmove) then
    do i=0, n_pe2recv_fluid-1
      call MPI_IRECV(buffpopr(0,i),n_pop2recv_fluid(i),MYFLOAT, &
       i_pe2recv_fluid(i),i_pe2recv_fluid(i)+movetag, &
       MPI_COMM_WORLD,requestm(i),ierr)
    enddo
    firstmove=.false.
  endif
  do i=0, n_pe2send_fluid-1
    do j=0,n_pop2send_fluid(i)-1
      i4=i_pop2send_fluid(1,j,i)
      mym=i4find(i4)
      ktemp = mym(3)
      jtemp = mym(2)
      itemp = mym(1)
      buffpops(j,i)=fluidsub(itemp,jtemp,ktemp,i_pop2send_fluid(0,j,i))
    enddo
    if(syncsend.eq.1) then
      call MPI_SEND(buffpops(0,i),n_pop2send_fluid(i),MYFLOAT, &
       i_pe2send_fluid(i),idrank+movetag,MPI_COMM_WORLD,ierr)
    else
      call MPI_ISEND(buffpops(0,i),n_pop2send_fluid(i),MYFLOAT, &
       i_pe2send_fluid(i),idrank+movetag,MPI_COMM_WORLD,request(i),ierr)
    endif
  enddo
   
  return

 end subroutine commspop
 
 subroutine commrpop(fluidsub,lparticles,isfluid)
 
!***********************************************************************
!     
!     LBsoft subroutine for receiving the population over the MPI
!     network
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************

  implicit none

  real(kind=PRC), allocatable, dimension(:,:,:,:)  :: fluidsub
  logical, intent(in) :: lparticles
  integer(kind=1), allocatable, dimension(:,:,:), intent(in) :: isfluid
  
  integer :: request(0:maxneigh-1)
  integer :: istatus(MPI_STATUS_SIZE,0:maxneigh-1)
  integer :: i,j,mym(3),idir,iopp
  integer(kind=IPRC) :: i4, itemp, jtemp, ktemp
  
  
  if(mxrank==1)return
  if(n_pe2recv_fluid.gt.0) then
    call MPI_WAITALL(n_pe2recv_fluid,requestm,istatus,ierr)
   endif
   if(lparticles)then
     do i=0, n_pe2recv_fluid-1
       do j=0,n_pop2recv_fluid(i)-1
         i4=i_pop2recv_fluid(1,j,i)
         mym=i4find(i4)
         ktemp = mym(3)
         jtemp = mym(2)
         itemp = mym(1)
         idir=i_pop2recv_fluid(0,j,i)
         iopp=opp(idir)
         if(isfluid(itemp+ex(iopp),jtemp+ey(iopp),ktemp+ez(iopp))<2 .or. &
          isfluid(itemp+ex(iopp),jtemp+ey(iopp),ktemp+ez(iopp))>5)then
           fluidsub(itemp,jtemp,ktemp,idir)=buffpopr(j,i)
         endif
       enddo
     enddo
   else
     do i=0, n_pe2recv_fluid-1
       do j=0,n_pop2recv_fluid(i)-1
         i4=i_pop2recv_fluid(1,j,i)
         mym=i4find(i4)
         ktemp = mym(3)
         jtemp = mym(2)
         itemp = mym(1)
         fluidsub(itemp,jtemp,ktemp,i_pop2recv_fluid(0,j,i))=buffpopr(j,i)
       enddo
     enddo
   endif
   
   do i=0, n_pe2recv_fluid-1
     call MPI_IRECV(buffpopr(0,i),n_pop2recv_fluid(i),MYFLOAT, &
      i_pe2recv_fluid(i),i_pe2recv_fluid(i)+movetag, &
      MPI_COMM_WORLD,requestm(i),ierr)
   enddo

   if(syncsend.eq.0 .and. n_pe2send_fluid.gt.0) then
     call MPI_WAITALL(n_pe2send_fluid,request,istatus,ierr)
   endif

  return
 
 end subroutine commrpop
 
 function i4back(i,j,k) result(myout)
  
!***********************************************************************
!     
!     LBsoft subroutine for finding a unique integer number i4-form
!     for each i,j,k triplet
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification September 2018
!     
!***********************************************************************
  
   implicit none
  
   DEALLOCATE(ownern)
  
   return

 end subroutine deallocate_ownern
 
 subroutine driving_print_binary_1d_vtk(iotest,headoff,nbyte,nn,myvar,e_io)
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the writing of a scalar field 
!     with single precision in VTK legacy binary format in parallel IO
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: iotest,headoff,nbyte,nn
  real(4), intent(in), dimension(:,:,:), allocatable  :: myvar
  integer, intent(inout) :: e_io
  
  if(mxrank==1)then
    call print_binary_1d_vtk_serial(iotest,headoff,nbyte,nn, &
      reshape(myvar,(/nn/)),e_io)
  else
    call print_binary_1d_vtk_col(iotest,headoff,nbyte,nn, &
      myvar,e_io)
  endif
  
  return
  
 end subroutine driving_print_binary_1d_vtk
 
 subroutine print_binary_1d_vtk_serial(iotest,headoff,nbyte,nn,myvar1d,E_IO)
 
!***********************************************************************
!     
!     LBsoft subroutine for writing a scalar field with single
!     precision in VTK legacy binary format in serial IO
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
 
 subroutine collective_writeFile(it,mynamefile,myvar,nbuff, &
  ix,fx,iy,fy,iz,fz)
  
!***********************************************************************
!     
!     LBsoft subroutine for writing a scalar field with PRC
!     precision in binary format in parallel IO
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: it,nbuff,ix,fx,iy,fy,iz,fz
  real(kind=PRC), &
   dimension(ix-nbuff:fx+nbuff,iy-nbuff:fy+nbuff,iz-nbuff:fz+nbuff), &
   intent(in)  :: myvar
  character(len=120),intent(in) :: mynamefile
  integer :: ierr,fh, filetype,mystarts(3),i
  integer :: memDims(3),memOffs(3), imemtype
  integer (kind=MPI_OFFSET_KIND) offset
  INTEGER STATUS(MPI_STATUS_SIZE)


  call MPI_FILE_OPEN(MPI_COMM_WORLD, mynamefile, &
   MPI_MODE_CREATE+MPI_MODE_WRONLY, &
   MPI_INFO_NULL, fh, ierr)
  !if (idrank== 0) print *, "MPI_FILE_OPEN)fh=",fh,"ierr=",ierr

  ! We need global sizes: globalDims
  ! We need local sizes: lDims
  ! We need start indices: mymin-1
  mystarts = mymin - 1

  call MPI_Type_create_subarray(3, globalDims, lDims, mystarts, &
   MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, filetype, ierr);
  !if (idrank== 0) print *, "MPI_Type_create_subarray)ierr=",ierr, globalDims, lDims, mystarts

  call MPI_Type_commit(filetype, ierr);

  call MPI_File_Set_View(fh, 0_MPI_OFFSET_KIND, MPI_DOUBLE, filetype, &
   "native", MPI_INFO_NULL, ierr)
  !if (idrank== 0) print *, "MPI_File_Set_View)ierr=",ierr

  ! We need full local sizes: memDims
  memDims = lDims + 2*nbuff
  memOffs = [ nbuff, nbuff, nbuff ]

  call MPI_TYPE_CREATE_SUBARRAY(3, memDims, lDims, memOffs, &
   MPI_ORDER_FORTRAN, MPI_DOUBLE, imemtype, ierr)
  !if (idrank== 0) print *, "MPI_Type_create_subarray)ierr=",ierr, memDims, lDims, memOffs

  call MPI_TYPE_COMMIT(imemtype, ierr)

  call MPI_FILE_WRITE_ALL(fh, myvar, 1, imemtype, status, ierr)
  !if (idrank== 0) print *, "MPI_FILE_WRITE_ALL)ierr=",ierr

  call MPI_FILE_CLOSE(fh, ierr)
  
  return
  
 end subroutine collective_writeFile

 subroutine collective_writeFile_int(it,mynamefile,myvar,nbuff, &
  ix,fx,iy,fy,iz,fz)
 
!***********************************************************************
!     
!     LBsoft subroutine for writing a scalar field with INTEGER
!     of 1 byte in binary format in parallel IO
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: it,nbuff,ix,fx,iy,fy,iz,fz
  integer(kind=1), &
   dimension(ix-nbuff:fx+nbuff,iy-nbuff:fy+nbuff,iz-nbuff:fz+nbuff), &
   intent(in)  :: myvar
  character(len=120),intent(in) :: mynamefile
  integer :: ierr,fh, filetype,mystarts(3),i
  integer :: memDims(3),memOffs(3), imemtype
  integer (kind=MPI_OFFSET_KIND) offset
  INTEGER STATUS(MPI_STATUS_SIZE)


  call MPI_FILE_OPEN(MPI_COMM_WORLD, mynamefile, &
   MPI_MODE_CREATE+MPI_MODE_WRONLY, &
   MPI_INFO_NULL, fh, ierr)
  !if (id_rank== 0) print *, "MPI_FILE_OPEN)fh=",fh,"ierr=",ierr

  ! We need global sizes: globalDims
  ! We need local sizes: lDims
  ! We need start indices: mymin-1
  mystarts = mymin - 1

  call MPI_Type_create_subarray(3, globalDims, lDims, mystarts, &
   MPI_ORDER_FORTRAN, MPI_CHAR, filetype, ierr);
  !if (id_rank== 0) print *, "MPI_Type_create_subarray)ierr=",ierr, globalDims, lDims, mystarts

  call MPI_Type_commit(filetype, ierr);

  call MPI_File_Set_View(fh, 0_MPI_OFFSET_KIND, MPI_CHAR, filetype, &
   "native", MPI_INFO_NULL, ierr)
  !if (id_rank== 0) print *, "MPI_File_Set_View)ierr=",ierr

  ! We need full local sizes: memDims
  memDims = lDims + 2*nbuff
  memOffs = [ nbuff, nbuff, nbuff ]

  call MPI_TYPE_CREATE_SUBARRAY(3, memDims, lDims, memOffs, &
   MPI_ORDER_FORTRAN, MPI_CHAR, imemtype, ierr)
  !if (id_rank== 0) print *, "MPI_Type_create_subarray)ierr=",ierr, memDims, lDims, memOffs

  call MPI_TYPE_COMMIT(imemtype, ierr)

  call MPI_FILE_WRITE_ALL(fh, myvar, 1, imemtype, status, ierr)
  !if (id_rank== 0) print *, "MPI_FILE_WRITE_ALL)ierr=",ierr

  call MPI_FILE_CLOSE(fh, ierr)
  
  return
  
 end subroutine collective_writeFile_int
 
 subroutine collective_readFile(it,mynamefile,myvar,nbuff, &
  ix,fx,iy,fy,iz,fz)

!***********************************************************************
!     
!     LBsoft subroutine for reading a scalar field with PRC
!     precision in binary format in parallel IO
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************

  implicit none
  
  integer, intent(in) :: it,nbuff,ix,fx,iy,fy,iz,fz
  real(kind=PRC), &
   dimension(ix-nbuff:fx+nbuff,iy-nbuff:fy+nbuff,iz-nbuff:fz+nbuff) :: &
   myvar
  character(len=120),intent(in) :: mynamefile
  integer :: ierr,fh, filetype,mystarts(3),i
  integer :: memDims(3),memOffs(3), imemtype
  integer (kind=MPI_OFFSET_KIND) offset
  INTEGER STATUS(MPI_STATUS_SIZE)


  call MPI_FILE_OPEN(MPI_COMM_WORLD, mynamefile, &
    MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)
  !if (id_rank== 0) print *, "MPI_FILE_OPEN)fh=",fh,"ierr=",ierr

  ! We need global sizes: globalDims
  ! We need local sizes: lDims
  ! We need start indices: mymin-1
  mystarts = mymin - 1

  call MPI_Type_create_subarray(3, globalDims, lDims, mystarts, &
    MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, filetype, ierr);
  !if (id_rank== 0) print *, "MPI_Type_create_subarray)ierr=",ierr, globalDims, lDims, mystarts

  call MPI_Type_commit(filetype, ierr);

  call MPI_File_Set_View(fh, 0_MPI_OFFSET_KIND, MPI_DOUBLE, filetype, &
    "native", MPI_INFO_NULL, ierr)
  !if (id_rank== 0) print *, "MPI_File_Set_View)ierr=",ierr

  ! We need full local sizes: memDims
  memDims = lDims + 2*nbuff
  memOffs = [ nbuff, nbuff, nbuff ]

  call MPI_TYPE_CREATE_SUBARRAY(3, memDims, lDims, memOffs, &
   MPI_ORDER_FORTRAN, MPI_DOUBLE, imemtype, ierr)
  !if (id_rank== 0) print *, "MPI_Type_create_subarray)ierr=",ierr, memDims, lDims, memOffs

  call MPI_TYPE_COMMIT(imemtype, ierr)

  call MPI_FILE_READ_ALL(fh, myvar, 1, imemtype, status, ierr)
  !if (id_rank== 0) print *, "MPI_FILE_WRITE_ALL)ierr=",ierr

  call MPI_FILE_CLOSE(fh, ierr)
  
  return
 
 end subroutine collective_readFile

 subroutine collective_readFile_int(it,mynamefile,myvar,nbuff, &
  ix,fx,iy,fy,iz,fz)
    
!***********************************************************************
!     
!     LBsoft subroutine for reading a scalar field with INTEGER
!     of 1 byte in binary format in parallel IO
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
  implicit none
  
  integer, intent(in) :: it,nbuff,ix,fx,iy,fy,iz,fz
  integer(kind=1), &
   dimension(ix-nbuff:fx+nbuff,iy-nbuff:fy+nbuff,iz-nbuff:fz+nbuff) :: &
   myvar
  character(len=120),intent(in) :: mynamefile
  integer :: ierr,fh, filetype,mystarts(3),i
  integer :: memDims(3),memOffs(3), imemtype
  integer (kind=MPI_OFFSET_KIND) offset
  INTEGER STATUS(MPI_STATUS_SIZE)


  call MPI_FILE_OPEN(MPI_COMM_WORLD, mynamefile, &
    MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)
  !if (id_rank== 0) print *, "MPI_FILE_OPEN)fh=",fh,"ierr=",ierr

  ! We need global sizes: globalDims
  ! We need local sizes: lDims
  ! We need start indices: mymin-1
  mystarts = mymin - 1

  call MPI_Type_create_subarray(3, globalDims, lDims, mystarts, &
    MPI_ORDER_FORTRAN, MPI_CHAR, filetype, ierr);

  call MPI_Type_commit(filetype, ierr);

  call MPI_File_Set_View(fh, 0_MPI_OFFSET_KIND, MPI_CHAR, filetype, &
    "native", MPI_INFO_NULL, ierr)

  ! We need full local sizes: memDims
  memDims = lDims + 2*nbuff
  memOffs = [ nbuff, nbuff, nbuff ]

  call MPI_TYPE_CREATE_SUBARRAY(3, memDims, lDims, memOffs, &
   MPI_ORDER_FORTRAN, MPI_CHAR, imemtype, ierr)

  call MPI_TYPE_COMMIT(imemtype, ierr)

  call MPI_FILE_READ_ALL(fh, myvar, 1, imemtype, status, ierr)

  call MPI_FILE_CLOSE(fh, ierr)
  
  return
  
 end subroutine collective_readFile_int
 
 subroutine collective_readFile_pops(it,mynamefile,f00,f01, &
  f02,f03,f04,f05,f06,f07,f08,f09,f10,f11,f12,f13,f14,f15,f16,f17,f18, &
  nbuff,ix,fx,iy,fy,iz,fz)

!***********************************************************************
!     
!     LBsoft subroutine for reading a scalar field with PRC
!     precision in binary format in parallel IO
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
  
  implicit none
  integer, intent(in) :: it,nbuff,ix,fx,iy,fy,iz,fz
  real(kind=PRC), &
   dimension(ix-nbuff:fx+nbuff,iy-nbuff:fy+nbuff,iz-nbuff:fz+nbuff) :: &
   f00,f01, &
   f02,f03,f04,f05,f06,f07,f08,f09,f10,f11,f12,f13,f14,f15,f16,f17,f18
  character(len=120),intent(in) :: mynamefile
  integer :: ierr,fh, filetype,mystarts(3),i
  integer :: memDims(3),memOffs(3), imemtype
  integer (kind=MPI_OFFSET_KIND) offset
  INTEGER STATUS(MPI_STATUS_SIZE)


  call MPI_FILE_OPEN(MPI_COMM_WORLD, mynamefile, &
   MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)
  !if (id_rank== 0) print *, "MPI_FILE_OPEN)fh=",fh,"ierr=",ierr

  ! We need global sizes: globalDims
  ! We need local sizes: lDims
  ! We need start indices: mymin-1
  mystarts = mymin - 1

  call MPI_Type_create_subarray(3, globalDims, lDims, mystarts, &
   MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, filetype, ierr);

  call MPI_Type_commit(filetype, ierr);

  call MPI_File_Set_View(fh, 0_MPI_OFFSET_KIND, MPI_DOUBLE, filetype, &
    "native", MPI_INFO_NULL, ierr)

  ! We need full local sizes: memDims
  memDims = lDims + 2*nbuff
  memOffs = [ nbuff, nbuff, nbuff ]

  call MPI_TYPE_CREATE_SUBARRAY(3, memDims, lDims, memOffs, &
   MPI_ORDER_FORTRAN, MPI_DOUBLE, imemtype, ierr)

  call MPI_TYPE_COMMIT(imemtype, ierr)

  call MPI_FILE_READ_ALL(fh, f00, 1, imemtype, status, ierr)
  call MPI_FILE_READ_ALL(fh, f01, 1, imemtype, status, ierr)
  call MPI_FILE_READ_ALL(fh, f02, 1, imemtype, status, ierr)
  call MPI_FILE_READ_ALL(fh, f03, 1, imemtype, status, ierr)
  call MPI_FILE_READ_ALL(fh, f04, 1, imemtype, status, ierr)
  call MPI_FILE_READ_ALL(fh, f05, 1, imemtype, status, ierr)
  call MPI_FILE_READ_ALL(fh, f06, 1, imemtype, status, ierr)
  call MPI_FILE_READ_ALL(fh, f07, 1, imemtype, status, ierr)
  call MPI_FILE_READ_ALL(fh, f08, 1, imemtype, status, ierr)
  call MPI_FILE_READ_ALL(fh, f09, 1, imemtype, status, ierr)
  call MPI_FILE_READ_ALL(fh, f10, 1, imemtype, status, ierr)
  call MPI_FILE_READ_ALL(fh, f11, 1, imemtype, status, ierr)
  call MPI_FILE_READ_ALL(fh, f12, 1, imemtype, status, ierr)
  call MPI_FILE_READ_ALL(fh, f13, 1, imemtype, status, ierr)
  call MPI_FILE_READ_ALL(fh, f14, 1, imemtype, status, ierr)
  call MPI_FILE_READ_ALL(fh, f15, 1, imemtype, status, ierr)
  call MPI_FILE_READ_ALL(fh, f16, 1, imemtype, status, ierr)
  call MPI_FILE_READ_ALL(fh, f17, 1, imemtype, status, ierr)
  call MPI_FILE_READ_ALL(fh, f18, 1, imemtype, status, ierr)

  call MPI_FILE_CLOSE(fh, ierr)
  
  return
  
 end subroutine collective_readFile_pops
 
 subroutine collective_writeFile_pops(it, mynamefile, f00,f01, &
  f02,f03,f04,f05,f06,f07,f08,f09,f10,f11,f12,f13,f14,f15,f16,f17,f18, &
  nbuff,ix,fx,iy,fy,iz,fz)
  
!***********************************************************************
!     
!     LBsoft subroutine for writing a scalar field with PRC
!     precision in binary format in parallel IO
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: it,nbuff,ix,fx,iy,fy,iz,fz
  real(kind=PRC), &
   dimension(ix-nbuff:fx+nbuff,iy-nbuff:fy+nbuff,iz-nbuff:fz+nbuff) :: &
   f00,f01, &
   f02,f03,f04,f05,f06,f07,f08,f09,f10,f11,f12,f13,f14,f15,f16,f17,f18
  character(len=120),intent(in) :: mynamefile
  integer :: ierr,fh, filetype,mystarts(3),i
  integer :: memDims(3),memOffs(3), imemtype
  integer (kind=MPI_OFFSET_KIND) offset
  integer :: istatus(MPI_STATUS_SIZE)


  call MPI_FILE_OPEN(MPI_COMM_WORLD, mynamefile, &
   MPI_MODE_CREATE+MPI_MODE_WRONLY, &
   MPI_INFO_NULL, fh, ierr)
  !if (id_rank== 0) print *, "MPI_FILE_OPEN)fh=",fh,"ierr=",ierr

  ! We need global sizes: globalDims
  ! We need local sizes: lDims
  ! We need start indices: mymin-1
  mystarts = mymin - 1

  call MPI_Type_create_subarray(3, globalDims, lDims, mystarts, &
   MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, filetype, ierr);
  !if (idrank== 0) print *, "MPI_Type_create_subarray)ierr=",ierr, globalDims, lDims, mystarts

  call MPI_Type_commit(filetype, ierr);

  call MPI_File_Set_View(fh, 0_MPI_OFFSET_KIND, MPI_DOUBLE, filetype, &
   "native", MPI_INFO_NULL, ierr)
  !if (idrank== 0) print *, "MPI_File_Set_View)ierr=",ierr
  
  ! We need full local sizes: memDims
  memDims = lDims + 2*nbuff
  memOffs = [ nbuff, nbuff, nbuff ]
  
  call MPI_TYPE_CREATE_SUBARRAY(3, memDims, lDims, memOffs, MPI_ORDER_FORTRAN, MPI_DOUBLE, imemtype, ierr)
  !if (idrank== 0) print *, "MPI_Type_create_subarray)ierr=",ierr, memDims, lDims, memOffs
  
  call MPI_TYPE_COMMIT(imemtype, ierr)
  
  call MPI_FILE_WRITE_ALL(fh, f00, 1, imemtype, istatus, ierr)
  call MPI_FILE_WRITE_ALL(fh, f01, 1, imemtype, istatus, ierr)
  call MPI_FILE_WRITE_ALL(fh, f02, 1, imemtype, istatus, ierr)
  call MPI_FILE_WRITE_ALL(fh, f03, 1, imemtype, istatus, ierr)
  call MPI_FILE_WRITE_ALL(fh, f04, 1, imemtype, istatus, ierr)
  call MPI_FILE_WRITE_ALL(fh, f05, 1, imemtype, istatus, ierr)
  call MPI_FILE_WRITE_ALL(fh, f06, 1, imemtype, istatus, ierr)
  call MPI_FILE_WRITE_ALL(fh, f07, 1, imemtype, istatus, ierr)
  call MPI_FILE_WRITE_ALL(fh, f08, 1, imemtype, istatus, ierr)
  call MPI_FILE_WRITE_ALL(fh, f09, 1, imemtype, istatus, ierr)
  call MPI_FILE_WRITE_ALL(fh, f10, 1, imemtype, istatus, ierr)
  call MPI_FILE_WRITE_ALL(fh, f11, 1, imemtype, istatus, ierr)
  call MPI_FILE_WRITE_ALL(fh, f12, 1, imemtype, istatus, ierr)
  call MPI_FILE_WRITE_ALL(fh, f13, 1, imemtype, istatus, ierr)
  call MPI_FILE_WRITE_ALL(fh, f14, 1, imemtype, istatus, ierr)
  call MPI_FILE_WRITE_ALL(fh, f15, 1, imemtype, istatus, ierr)
  call MPI_FILE_WRITE_ALL(fh, f16, 1, imemtype, istatus, ierr)
  call MPI_FILE_WRITE_ALL(fh, f17, 1, imemtype, istatus, ierr)
  call MPI_FILE_WRITE_ALL(fh, f18, 1, imemtype, istatus, ierr)
  !if (idrank== 0) print *, "MPI_FILE_WRITE_ALL)ierr=",ierr

  call MPI_FILE_CLOSE(fh, ierr)
  
  return
      
 end subroutine collective_writeFile_pops
  
 subroutine print_binary_1d_vtk_col(iotest,headoff,nbyte,nn,myvar,e_io)
 
!***********************************************************************
!     
!     LBsoft subroutine for writing a scalar field with single
!     precision in VTK legacy binary format in parallel IO
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: iotest,headoff,nbyte,nn
  real(4), intent(in), dimension(:,:,:) :: myvar
  integer, intent(out) :: e_io
  
  integer(kind=MPI_OFFSET_KIND) :: myoffset
  
  integer, parameter :: byteint=4
  
  integer :: filetypesub
  integer, parameter :: nbuffsub=0
  
  integer :: memDims(3),memOffs(3),mystarts(3),imemtype
  
  myoffset=int(headoff,kind=MPI_OFFSET_KIND)
  
  mystarts = mymin - 1
  
  if(idrank==0)call MPI_File_write_at(iotest,myoffset,int(nbyte,kind=4),1, &
     MPI_INTEGER,MPI_STATUS_IGNORE,e_io)
  myoffset = myoffset + 1 * byteint
  
  call MPI_Type_create_subarray(3,globalDims,lDims,mystarts, &
    MPI_ORDER_FORTRAN,MPI_REAL4,filetypesub,e_io)
        
  call MPI_Type_commit(filetypesub, e_io)
   
  call MPI_File_Set_View(iotest,myoffset,MPI_REAL4,filetypesub, &
    "native",MPI_INFO_NULL,e_io)
  ! We need full local sizes: memDims
  memDims = lDims + 2*nbuffsub
  memOffs = [ nbuffsub, nbuffsub, nbuffsub ]

  call MPI_TYPE_CREATE_SUBARRAY(3,memDims,lDims,memOffs, &
   MPI_ORDER_FORTRAN,MPI_REAL4,imemtype,e_io)

  call MPI_TYPE_COMMIT(imemtype,e_io)

  call MPI_FILE_WRITE_ALL(iotest,myvar,1,imemtype,MPI_STATUS_IGNORE,e_io)
  
  return
  
 end subroutine print_binary_1d_vtk_col
 
 subroutine driving_print_binary_3d_vtk(iotest,headoff,nbyte,nn,myvar1, &
  myvar2,myvar3,e_io)
  
!***********************************************************************
!     
!     LBsoft subroutine for driving the writing of a vector field 
!     with single precision in VTK legacy binary format in parallel IO
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: iotest,headoff,nbyte,nn
  real(4), intent(in), dimension(:,:,:), allocatable :: myvar1,myvar2,myvar3
  integer, intent(inout) :: e_io
  
  integer :: i,j,k,l
  integer :: nx1,nx2,ny1,ny2,nz1,nz2
  real(kind=4), allocatable, dimension(:,:,:,:) :: service4
  
  if(mxrank==1)then
    call print_binary_3d_vtk_serial(iotest,headoff,nbyte,nn, &
     reshape(myvar1,(/nn/)),reshape(myvar2,(/nn/)), &
     reshape(myvar3,(/nn/)),e_io)
  
    return
  endif
  
  nx1=mymin(1)
  nx2=mymax(1)
  ny1=mymin(2)
  ny2=mymax(2)
  nz1=mymin(3)
  nz2=mymax(3)
  
  allocate(service4(1:3,nx1:nx2,ny1:ny2,nz1:nz2))
  service4(1:3,nx1:nx2,ny1:ny2,nz1:nz2)=ZERO
  do k=nz1,nz2
    do j=ny1,ny2
      do i=nx1,nx2
        service4(1,i,j,k)=myvar1(i,j,k)
        service4(2,i,j,k)=myvar2(i,j,k)
        service4(3,i,j,k)=myvar3(i,j,k)
      enddo
    enddo
  enddo
  
  call print_binary_3d_vtk_col(iotest,headoff,nbyte,nn,service4,e_io)
     
  deallocate(service4)
  
  return
  
 end subroutine driving_print_binary_3d_vtk
 
 subroutine print_binary_3d_vtk_serial(iotest,headoff,nbyte,nn,myvarx, &
  myvary,myvarz,E_IO)
  
!***********************************************************************
!     
!     LBsoft subroutine for writing a vector field with single
!     precision in VTK legacy binary format in serial IO
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
 
 subroutine print_binary_3d_vtk_col(iotest,headoff,nbyte,nn,myvar,e_io)
  
!***********************************************************************
!     
!     LBsoft subroutine for writing a vector field with single
!     precision in VTK legacy binary format in parallel IO
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: iotest,headoff,nbyte,nn
  real(4), intent(in), dimension(:,:,:,:) :: myvar
  integer, intent(out) :: e_io
  
  integer(kind=MPI_OFFSET_KIND) :: myoffset
  
  integer, parameter :: byteint=4
  
  integer :: filetypesub
  integer, parameter :: nbuffsub=0
  
  integer :: imemtype
  
  integer, dimension(4) :: velglobalDims,velldims,velmystarts, &
   velmemDims,velmemOffs
  
  myoffset=int(headoff,kind=MPI_OFFSET_KIND)
  
  velglobalDims(1)=3
  velglobalDims(2:4)=globalDims(1:3)
  velldims(1)=3
  velldims(2:4)=lDims(1:3)
  velmystarts(1) = 0
  velmystarts(2:4) = mymin(1:3)-1
  
  if(idrank==0)call MPI_File_write_at(iotest,myoffset,int(nbyte,kind=4),1, &
     MPI_INTEGER,MPI_STATUS_IGNORE,e_io)
  myoffset = myoffset + 1 * byteint
  
  call MPI_Type_create_subarray(4,velglobalDims,velldims,velmystarts, &
    MPI_ORDER_FORTRAN,MPI_REAL4,filetypesub,e_io)
        
  call MPI_Type_commit(filetypesub, e_io)
   
  call MPI_File_Set_View(iotest,myoffset,MPI_REAL4,filetypesub, &
    "native",MPI_INFO_NULL,e_io)
  ! We need full local sizes: memDims
  velmemDims(1) = vellDims(1)
  velmemDims(2:4) = vellDims(2:4) + 2*nbuffsub
  velmemOffs = [ 0, nbuffsub, nbuffsub, nbuffsub ]

  call MPI_TYPE_CREATE_SUBARRAY(4,velmemDims,velldims,velmemOffs, &
   MPI_ORDER_FORTRAN,MPI_REAL4,imemtype,e_io)

  call MPI_TYPE_COMMIT(imemtype,e_io)

  call MPI_FILE_WRITE_ALL(iotest,myvar,1,imemtype,MPI_STATUS_IGNORE,e_io)
  
  return
  
 end subroutine print_binary_3d_vtk_col
  
 end module version_mod
