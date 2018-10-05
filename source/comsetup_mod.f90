!TBD: Massimo August 12 2018
!(1)check the limits of each subdomain
!(2)check for periodic boundary conditions (done, to be tested)
!(3)remove local nbuff and use the variable defined in fluids module
!(4)check for bounceback conditions
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include <default_macro.h>

MODULE lbempi_mod
  use version_mod, only : idrank, mxrank,finalize_world,or_world_larr, &
   get_sync_world,sum_world_iarr
  use aop_mod 
  IMPLICIT NONE

  INTEGER :: myid=0, numprocs=1
  
  
#if LATTICE==319
 integer, parameter, private :: links=18
  !lattice vectors
 integer, dimension(0:links), parameter, private :: &
  ex = (/ 0, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1, 0, 0, 0, 0/)
   !      0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18
 integer, dimension(0:links), parameter, private :: &
  ey = (/ 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1/)
 integer, dimension(0:links), parameter, private :: &
  ez = (/ 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1/)
#endif


#if defined(MPI)

  ! USE mpi 
  include 'mpif.h'
  INTEGER, SAVE :: MYFLOAT=MPI_DOUBLE_PRECISION
  INTEGER, SAVE :: MYINT=MPI_INTEGER
  INTEGER, SAVE :: MYINT1=MPI_INTEGER1
  INTEGER, SAVE :: MYINT2=MPI_INTEGER2
  INTEGER, SAVE :: MYCOMM=MPI_COMM_WORLD

#endif
!!! temporary parameter
  
!!!

  INTEGER, PARAMETER :: maxneigh=32
  INTEGER, PARAMETER :: npop=19
  INTEGER, PARAMETER :: movetag=393000, &
           & findntag=524000,findptag=655000,findwtag=786000
           
           
  INTEGER, PARAMETER :: findntag_hvar=findntag*2
  INTEGER, PARAMETER :: tag_hvar=movetag*2
  INTEGER, PARAMETER :: tag_hvar2=movetag*3
  INTEGER, PARAMETER :: tag_isfluid=movetag*4
  
  integer rc, ierr
  integer,PUBLIC :: lnofsite, itpppsize

  integer, protected, public :: domdec
  integer, parameter :: syncsend=1
  INTEGER,SAVE :: requestm(0:maxneigh-1), requestp(0:maxneigh-1), &
       requestw(0:maxneigh-1), requesti(0:maxneigh-1), requesto(0:maxneigh-1), &
       requesth(0:maxneigh-1), requesth_wall(0:maxneigh-1), requests(0:maxneigh-1), &
       request_hvar(0:maxneigh-1),requestm_hvar(0:maxneigh-1), &
       request_hvar2(0:maxneigh-1),requestm_hvar2(0:maxneigh-1)
  LOGICAL :: firstmove=.true.
  
  LOGICAL :: allpbc

  INTEGER :: n_pe2recv_fluid, n_pe2send_fluid
  INTEGER :: i_pe2send_fluid(0:maxneigh-1), n_pop2send_fluid(0:maxneigh-1)
  INTEGER :: i_pe2recv_fluid(0:maxneigh-1), n_pop2recv_fluid(0:maxneigh-1)
  
  INTEGER :: n_pe2recv_fluid_hvar, n_pe2send_fluid_hvar
  INTEGER :: i_pe2send_fluid_hvar(0:maxneigh-1), n_var2send_fluid(0:maxneigh-1)
  INTEGER :: i_pe2recv_fluid_hvar(0:maxneigh-1), n_var2recv_fluid(0:maxneigh-1)

  integer, protected, public :: nprocz, nprocy, nprocx
  INTEGER :: nxy2, nx2
  integer :: nbuff_comm

  integer, POINTER :: countnpp(:)
  INTEGER(KIND=2), ALLOCATABLE :: ownern(:)

  INTEGER(kind=IPRC), ALLOCATABLE :: i_pop2send_fluid(:,:,:), i_pop2recv_fluid(:,:,:)
  INTEGER(kind=IPRC), ALLOCATABLE :: i_var2send_fluid(:,:,:), i_var2recv_fluid(:,:,:)
  INTEGER(kind=IPRC), ALLOCATABLE :: i_pop2send_wall(:,:,:), i_pop2recv_wall(:,:,:)
  INTEGER(kind=IPRC), ALLOCATABLE :: i_pop2send_inlet(:,:,:), i_pop2recv_inlet(:,:,:)
  INTEGER(kind=IPRC), ALLOCATABLE :: i_pop2send_outlet(:,:,:), i_pop2recv_outlet(:,:,:)
  
  real(kind=PRC), POINTER :: buffpops(:,:), buffpopr(:,:)
  real(kind=PRC), POINTER :: bufftags(:,:), bufftagr(:,:)
  
  real(kind=PRC), POINTER :: buffs_hvar(:,:), buffr_hvar(:,:)
  real(kind=PRC), POINTER :: buffs_hvar2(:,:), buffr_hvar2(:,:)
  
  INTEGER(KIND=1), POINTER :: buffs_isfluid(:,:), buffr_isfluid(:,:)
  
  logical, save :: ldo_second=.false.
  
  integer, dimension(:), allocatable, public, protected :: gminx,gmaxx, &
   gminy,gmaxy,gminz,gmaxz
  
  public :: commspop, commrpop
  PUBLIC :: buff_transf_r,buff_transf_s
  REAL(kind=PRC),ALLOCATABLE,DIMENSION(:,:) :: buff_transf_r,buff_transf_s
#if defined(FNMD)
    INTEGER, PARAMETER :: mddirnpop=27
    INTEGER :: mddir(0:mddirnpop-npop-1)
#endif
  
  public :: setupcom
  public :: set_domdec
  public :: set_domain
  public :: i4back
  public :: i4find
  public :: comm_hvar
  public :: deallocate_ownern
  public :: commexch_dens
  public :: commwait_dens
  public :: comm_init_isfluid
  public :: create_findneigh_list_hvar
  public :: create_findneigh_list_pops
  public :: commexch_vel_component
  public :: commwait_vel_component
  
CONTAINS
#define LARGEINT 1073741824

 subroutine set_domdec(itemp)
 
  implicit none
  
  integer, intent(in) :: itemp
  
  domdec=itemp
  
  return
  
 end subroutine set_domdec
 
 subroutine set_domain(itemp1,itemp2,itemp3)
 
  implicit none
  
  integer, intent(in) :: itemp1,itemp2,itemp3
  
  nprocx=itemp1
  nprocy=itemp2
  nprocz=itemp3
  
  return
  
 end subroutine set_domain

#ifdef MPI
 SUBROUTINE setupcom(nx,ny,nz,nbuff,ibctype,ixpbc,iypbc,izpbc,minx,maxx, &
   miny,maxy,minz,maxz,lsingle_fluid)
  
 
  
  INTEGER, intent(in) :: nx,ny,nz,nbuff,ibctype
  INTEGER, intent(inout) :: ixpbc,iypbc,izpbc,minx,maxx,miny,maxy, &
   minz,maxz
  LOGICAL, intent(IN) :: lsingle_fluid
   
  
  
  INTEGER :: i
  INTEGER :: ownernlb,ownernub

  
  myid=idrank
  numprocs=mxrank
  

  nxy2=(nx+2*nbuff)*(ny+2*nbuff)
  nx2=nx+(2*nbuff)
  nbuff_comm=nbuff
  ALLOCATE(countnpp(0:numprocs-1))
  DO i=0,numprocs-1
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
  
  call cartdeco(nx,ny,nz,nbuff,.true.,ownern,ownernlb,minx,maxx,miny,maxy, &
   minz,maxz)
  
  
  
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
    
    allocate(gminx(0:mxrank-1))
    allocate(gmaxx(0:mxrank-1))
    allocate(gminy(0:mxrank-1))
    allocate(gmaxy(0:mxrank-1))
    allocate(gminz(0:mxrank-1))
    allocate(gmaxz(0:mxrank-1))
    
    
    gminx(0:mxrank-1)=0
    gmaxx(0:mxrank-1)=0
    gminy(0:mxrank-1)=0
    gmaxy(0:mxrank-1)=0
    gminz(0:mxrank-1)=0
    gmaxz(0:mxrank-1)=0
       
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
    call sum_world_iarr(gminx,mxrank)
    call sum_world_iarr(gmaxx,mxrank)
    call sum_world_iarr(gminy,mxrank)
    call sum_world_iarr(gmaxy,mxrank)
    call sum_world_iarr(gminz,mxrank)
    call sum_world_iarr(gmaxz,mxrank)
    
    return

END SUBROUTINE setupcom


subroutine create_findneigh_list_hvar(nx,ny,nz,nbuff,ibctype,ixpbc,iypbc,izpbc,minx,maxx, &
   miny,maxy,minz,maxz)

  implicit none
 
  INTEGER, intent(in) :: nx,ny,nz,nbuff,ibctype
  INTEGER, intent(inout) :: ixpbc,iypbc,izpbc,minx,maxx,miny,maxy, &
   minz,maxz
  

  call findneigh_hvar(nx,ny,nz,nbuff,ibctype,ixpbc,iypbc,izpbc,minx,maxx, &
   miny,maxy,minz,maxz)
   
  return
  
end subroutine

subroutine create_findneigh_list_pops(nx,ny,nz,nbuff,ibctype,isfluid,ixpbc,iypbc,izpbc,minx,maxx, &
   miny,maxy,minz,maxz)

  implicit none
  
  INTEGER, intent(in) :: nx,ny,nz,nbuff,ibctype
  INTEGER(kind=1), intent(in), allocatable, dimension(:,:,:) :: isfluid
  INTEGER, intent(inout) :: ixpbc,iypbc,izpbc,minx,maxx,miny,maxy, &
   minz,maxz
  
  
 !max  write(0,*)'countnpp(',myid,')=',countnpp(myid)
  call findneigh(nx,ny,nz,nbuff,ibctype,isfluid,ixpbc,iypbc,izpbc,minx,maxx,miny,maxy, &
   minz,maxz)
   
  return
  
end subroutine
  
#else
  subroutine setupcom(nx,ny,nz,nbuff,ibctype,ixpbc,iypbc,izpbc,minx,maxx, &
   miny,maxy,minz,maxz,lsingle_fluid)
  
!***********************************************************************
!     
!     LBsoft subroutine for the initial setup of the domain 
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification September 2018
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: nx,ny,nz,nbuff,ibctype
  integer, intent(inout) :: ixpbc,iypbc,izpbc,minx,maxx,miny,maxy, &
   minz,maxz
  
  LOGICAL, intent(IN) :: lsingle_fluid
  
  integer :: i
  integer :: ownernlb,ownernub

  
  myid=idrank
  numprocs=mxrank

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
  
  !ownernlb=(1-nbuff)+((1-nbuff)*nx2)+((1-nbuff)*nxy2)
  ownernlb=1
  ownernub=(nx+2*nbuff)*(ny+2*nbuff)*(nz+2*nbuff)
  
  ALLOCATE(ownern(ownernlb:ownernub))
  do i=ownernlb,ownernub
     ownern(i)=0
  enddo
  
  allocate(gminx(0:mxrank-1))
  allocate(gmaxx(0:mxrank-1))
  allocate(gminy(0:mxrank-1))
  allocate(gmaxy(0:mxrank-1))
  allocate(gminz(0:mxrank-1))
  allocate(gmaxz(0:mxrank-1))
  
  
  gminx(0:mxrank-1)=0
  gmaxx(0:mxrank-1)=0
  gminy(0:mxrank-1)=0
  gmaxy(0:mxrank-1)=0
  gminz(0:mxrank-1)=0
  gmaxz(0:mxrank-1)=0
     
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
  call sum_world_iarr(gminx,mxrank)
  call sum_world_iarr(gmaxx,mxrank)
  call sum_world_iarr(gminy,mxrank)
  call sum_world_iarr(gmaxy,mxrank)
  call sum_world_iarr(gminz,mxrank)
  call sum_world_iarr(gmaxz,mxrank)
  
  return
  
 end subroutine setupcom
 
 subroutine create_findneigh_list_hvar(nx,ny,nz,nbuff,ibctype,ixpbc,iypbc,izpbc,minx,maxx, &
   miny,maxy,minz,maxz)

  implicit none
 
  INTEGER, intent(in) :: nx,ny,nz,nbuff,ibctype
  INTEGER, intent(inout) :: ixpbc,iypbc,izpbc,minx,maxx,miny,maxy, &
   minz,maxz
  

   
  return
  
end subroutine

 subroutine create_findneigh_list_pops(nx,ny,nz,nbuff,ibctype,isfluid,ixpbc,iypbc,izpbc,minx,maxx, &
   miny,maxy,minz,maxz)

  implicit none
 
  INTEGER, intent(in) :: nx,ny,nz,nbuff,ibctype
  INTEGER(kind=1), intent(in), allocatable, dimension(:,:,:) :: isfluid
  INTEGER, intent(inout) :: ixpbc,iypbc,izpbc,minx,maxx,miny,maxy, &
   minz,maxz
  

   
  return
  
end subroutine

#endif
 
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
 
  SUBROUTINE cartdeco(nx,ny,nz,nbuff,first, ownern, ownernlb,minx,maxx,miny,maxy, &
   minz,maxz)

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
       
#if defined(FNMD)
       IF(numprocs.GT.1) THEN
          mddir(0)=+(nx+2*nbuff)*(ny+2*nbuff)+1+(nx+2*nbuff) !z+1, x+1, y+1
          mddir(1)=-(nx+2*nbuff)*(ny+2*nbuff)+1+(nx+2*nbuff) !z-1, x+1, y+1
          mddir(2)=+(nx+2*nbuff)*(ny+2*nbuff)-1+(nx+2*nbuff) !z+1, x-1, y+1
          mddir(3)=-(nx+2*nbuff)*(ny+2*nbuff)-1+(nx+2*nbuff) !z-1, x-1, y+1
          mddir(4)=+(nx+2*nbuff)*(ny+2*nbuff)+1-(nx+2*nbuff) !z+1, x+1, y-1
          mddir(5)=-(nx+2*nbuff)*(ny+2*nbuff)+1-(nx+2*nbuff) !z-1, x+1, y-1
          mddir(6)=+(nx+2*nbuff)*(ny+2*nbuff)-1-(nx+2*nbuff) !z+1, x-1. y-1
          mddir(7)=-(nx+2*nbuff)*(ny+2*nbuff)-1-(nx+2*nbuff) !z-1, x-1, y-1
       ENDIF
#endif
       minx=LARGEINT
       miny=LARGEINT
       minz=LARGEINT
       maxx=-LARGEINT
       maxy=-LARGEINT
       maxz=-LARGEINT
       
       IF(numprocs.EQ.1) THEN
          domdec=1
          nprocx=1
          nprocy=1
          nprocz=1
       ENDIF

       IF(domdec.LT.1.OR.domdec.GT.7) THEN
#if defined(MPI)
          IF(numprocs.GT.1) CALL MPI_FINALIZE(ierr)
#endif
          IF(myid.EQ.0) THEN
             PRINT*,'decomposition:', domdec
          ENDIF
          STOP
       ENDIF

       IF(domdec.GT.3.AND.MOD(numprocs,2).NE.0) THEN
#if defined(MPI)
          IF(numprocs.GT.1) CALL MPI_FINALIZE(ierr)
#endif
          IF(myid.EQ.0) THEN
             PRINT*,'Invalid number of processors',numprocs,&
                  'for this decomposition'
          ENDIF
          STOP
       ENDIF
       IF(domdec.GT.3.AND.numprocs.LT.4) THEN
#if defined(MPI)
          IF(numprocs.GT.1) CALL MPI_FINALIZE(ierr)
#endif
          IF(myid.EQ.0) THEN
             PRINT*,&
                  'At least 4 processors are required for this decomposition'
          ENDIF
          STOP
       ENDIF

       IF(domdec.EQ.7) THEN
          IF(numprocs.LT.8) THEN
#if defined(MPI)
             IF(numprocs.GT.1) CALL MPI_FINALIZE(ierr)
#endif
             IF(myid.EQ.0) THEN
                PRINT*,&
                     'At least 8 processors are required for this decomposition'
             ENDIF
             STOP
          ENDIF
          IF(numprocs.EQ.8) THEN
             nprocz=2
             nprocy=2
             nprocx=2
          ELSE IF(numprocs.EQ.16) THEN
             nprocz=2
             nprocy=2
             nprocx=4
          ELSE IF(numprocs.EQ.32) THEN
             nprocz=2
             nprocy=4
             nprocx=4
          ELSE IF(numprocs.EQ.64) THEN
             nprocz=4
             nprocy=4
             nprocx=4
          ELSE IF(numprocs.EQ.128) THEN
             nprocz=4
             nprocy=4
             nprocx=8
          ELSE IF(numprocs.EQ.256) THEN
             nprocz=4
             nprocy=8
             nprocx=8
          ELSE IF(numprocs.EQ.512) THEN
             nprocz=8
             nprocy=8
             nprocx=8
          ELSE IF(numprocs.EQ.1024) THEN
             nprocz=8
             nprocy=8
             nprocx=16
          ELSE IF(numprocs.EQ.2048) THEN
             nprocz=8
             nprocy=16
             nprocx=16
          ELSE IF(numprocs.EQ.4096) THEN
             nprocz=16
             nprocy=16
             nprocx=16
          ELSE IF(numprocs.EQ.8192) THEN
             nprocz=16
             nprocy=16
             nprocx=32
          ELSE IF(numprocs.EQ.16384) THEN
             nprocz=16
             nprocy=32
             nprocx=32
          ELSE IF(numprocs.EQ.32768) THEN
             nprocz=32
             nprocy=32
             nprocx=32
          ELSE IF(myid.EQ.0) THEN
             PRINT*,"Invalid number of processors",numprocs
#if defined(MPI)
             IF(numprocs.GT.1) CALL MPI_FINALIZE(ierr)
#endif
             STOP
          ENDIF
       ENDIF

       m=0
    ENDIF
    
    currproc=0
    curroffz=0
    curroffy=0
    IF(domdec.EQ.1) THEN      !along z
       slice=(nz)/numprocs
       IF(MOD((nz),numprocs).NE.0) THEN
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
          IF(slice.EQ.1.AND.iz.LT.(nz)) THEN
             currproc=currproc+1
             slice=((nz)-(iz))/(numprocs-currproc)
             IF(MOD((nz)-(iz),(numprocs-currproc)).NE.0) THEN
                slice=slice+1
             ENDIF
          ELSE
             slice=slice-1
          ENDIF
       ELSE IF(domdec.EQ.4.OR.domdec.EQ.5) THEN !along zy and zx
          IF(iz.GT.(nz)/2) THEN
             currpoff=numprocs/2
          ELSE
             currpoff=0
          ENDIF
          IF(domdec.EQ.4) THEN !along zy
             currproc=currpoff
             slice=(ny)/(numprocs/2)
             IF(MOD((ny),(numprocs/2)).NE.0) THEN
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
          slice=(ny)/numprocs
          IF(MOD((ny),numprocs).NE.0) THEN
             slice=slice+1
          ENDIF
       ENDIF
       DO iy=1, ny
          IF(domdec.EQ.2) THEN
             IF(slice.EQ.1.AND.iy.LT.(ny)) THEN
                currproc=currproc+1
                slice=((ny)-(iy))/(numprocs-currproc)
                IF(MOD((ny)-(iy),(numprocs-currproc)).NE.0) THEN
                   slice=slice+1
                ENDIF
             ELSE
                slice=slice-1
             ENDIF
          ELSE IF(domdec.EQ.4) THEN !along zy
             IF(slice.EQ.1.AND.iy.LT.(ny)) THEN
                currproc=currproc+1
                slice=((ny)-(iy))/ &
                     ((numprocs/2)-(currproc-currpoff))
                IF(MOD((ny)-(iy),&
                     (numprocs/2)-(currproc-currpoff)).NE.0) THEN
                   slice=slice+1
                ENDIF
             ELSE
                slice=slice-1
             ENDIF
          ELSE IF(domdec.EQ.6) THEN !along yx
             IF(iy.GT.(ny)/2) THEN
                currpoff=numprocs/2
             ELSE
                currpoff=0
             ENDIF
             currproc=currpoff
             slice=(nx)/(numprocs/2)
             IF(MOD((nx),(numprocs/2)).NE.0) THEN
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
             slice=(nx)/numprocs
             IF(MOD((nx),numprocs).NE.0) THEN
                slice=slice+1
             ENDIF
          ELSE IF(domdec.EQ.5) THEN !along zx
             currproc=currpoff
             slice=(nx)/(numprocs/2)
             IF(MOD((nx),(numprocs/2)).NE.0) THEN
                slice=slice+1
             ENDIF
          ENDIF
          DO ix=1, nx
             IF(domdec.EQ.7) THEN
                currproc=curroffx+curroffy*nprocx+curroffz*nprocx*nprocy
             ENDIF
             IF(currproc.LT.0.OR.currproc.GE.numprocs) THEN
                PRINT*,'Invalid currproc',currproc,nprocx,nprocy,nprocz,&
                     curroffx, curroffy,curroffz
                STOP
             ENDIF
             ! i4=iz*nxy2 + iy*nx2 + ix
             i4=i4back(ix,iy,iz)
             IF(domdec.EQ.3) THEN !along x
                IF(slice.EQ.1.AND.ix.LT.(nx)) THEN
                   currproc=currproc+1
                   slice=((nx)-(ix))/(numprocs-currproc)
                   IF(MOD((nx)-(ix),(numprocs-currproc)).NE.0) THEN
                      slice=slice+1
                   ENDIF
                ELSE
                   slice=slice-1
                ENDIF
             ELSE IF(domdec.EQ.5.OR.domdec.EQ.6) THEN !along zx or along yx
                IF(slice.EQ.1.AND.ix.LT.(nx)) THEN
                   currproc=currproc+1
                   slice=((nx)-(ix))/ &
                        ((numprocs/2)-(currproc-currpoff))
                   IF(MOD((nx)-(ix),&
                        (numprocs/2)-(currproc-currpoff)).NE.0) THEN
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
!max            write(0,*)'ownern(',i4,')=',currproc,'ix=',ix,'iy=',iy,'iz=',iz
                IF(currproc.eq.myid) THEN ! check the local boundaries
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
    
    ltest(1)=.false.
    
    do k=1-nbuff,nz+nbuff
       do j=1-nbuff,ny+nbuff
         do i=1-nbuff,nx+nbuff
           if(i<minx.or.j<miny.or.k<minz.or.i>maxx.or.j>maxy.or.k>maxz)then
             itemp=i
             jtemp=j
             ktemp=k
             if(i<1)itemp=1
             if(i>nx)itemp=nx
             if(j<1)jtemp=1
             if(j>ny)jtemp=ny
             if(k<1)ktemp=1
             if(k>nz)ktemp=nz
             i4target=i4back(i,j,k)
             i4=i4back(itemp,jtemp,ktemp)
             if(i4target<0 .or. i4<0)then
               write(6,*)'ERROR in cartdeco id ',idrank,' ijk ', &
                i,j,k,' i4target ',i4target,' ijktemp ',itemp,jtemp,ktemp,i4
               ltest(1)=.true.
             else
               ownern(i4target)=ownern(i4)
               currproc=ownern(i4)
               countnpp(currproc)=countnpp(currproc)+1
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
     
    
    
  END SUBROUTINE cartdeco

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !     max Possible domain decomposition:
  !     max 1D= along z=1, y=2, x=3
  !     max 2D= along zy=4, zx=5, yx=6
  !     max 3D= along zyx=7
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !> Returns the i4-form of a i,j,k triplet
  !! Assumes (nx+2*nbuff) (ny+2*nbuff) (nz+2*nbuff) strides
  FUNCTION i4back(i,j,k) RESULT(out)
     INTEGER,INTENT(in) :: i !< i-location on mesh
     INTEGER,INTENT(in) :: j !< j-location on mesh
     INTEGER,INTENT(in) :: k !< k-location on mesh
     INTEGER(kind=IPRC) :: out
#if 0
     out = INT(k,kind=IPRC)*nxy2 + INT(j,kind=IPRC)*nx2 + INT(i,kind=IPRC)
#else
     out = INT(k+nbuff_comm-1,kind=IPRC)*nxy2 + INT(j+nbuff_comm-1,kind=IPRC)*nx2 + INT(i+nbuff_comm,kind=IPRC)
     !myindex=(k+nbuff-1)*(ny+2*nbuff)*(nx+2*nbuff) + (j+nbuff-1)*(nx+2*nbuff) + i+nbuff
#endif
  END FUNCTION i4back
  
 !> Returns the i,j,k triplet-form of a i4-form
  FUNCTION i4find(i4sub)
    
    implicit none
    integer(kind=IPRC), intent(in) :: i4sub
    integer, dimension(3) :: i4find
    
#if 0
    i4find(3) = (i4sub/nxy2)
    i4find(2) = ((i4sub - i4find(3)*nxy2)/nx2)
    i4find(1) = (i4sub - i4find(3)*nxy2 - i4find(2)*nx2)
#else
    i4find(3) = ((i4sub-1)/nxy2)-nbuff_comm+1
    i4find(2) = ((i4sub-1 - (i4find(3)+nbuff_comm-1)*nxy2)/nx2)-nbuff_comm+1
    i4find(1) = (i4sub-nbuff_comm+1 - (i4find(3)+nbuff_comm-1)*nxy2 - (i4find(2)+nbuff_comm-1)*nx2)-1
#endif
    
    return
  END FUNCTION i4find
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE findneigh(nx,ny,nz,nbuff,ibctype,isfluid,ixpbc,iypbc,izpbc,minx,maxx, &
       miny,maxy,minz,maxz)

    IMPLICIT NONE

    INTEGER :: request(0:maxneigh-1)
#if defined(MPI)
    INTEGER :: status(MPI_STATUS_SIZE,0:maxneigh-1)
#endif

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

    INTEGER :: allneigh(0:maxneigh*numprocs-1)

    INTEGER :: sender(0:numprocs-1)
    
    logical, parameter :: lverbose=.false.
    logical :: ltemp(1)
    
    if(numprocs==1)return

    DO i=0,numprocs-1
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
    
    
    
    if(lverbose)write(6,*)'id=',idrank,'lminx=',lminx,'lmaxx=',lmaxx,'lminy=',lminy, &
              'lmaxy=',lmaxy,'lminz=',lminz,'lmaxz=',lmaxz  
    
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
                IF(ownern(i4).NE.myid.AND.isfluid(i,j,k)/=3) THEN
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
    DO i=0,numprocs-1

       IF(sender(i).GT.0) THEN ! select procs that talks to this one

          IF(j.GE.maxneigh) THEN
             PRINT *,'Task',myid,'Invalid index for myreceiver',j
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
    !max    write(0,*)'id=',myid,'maxsender=',maxsender,'n_pe2send_fluid=',n_pe2send_fluid
    DO i=0,maxneigh-1
       i_pe2recv_fluid(i)=-1
    ENDDO
#if defined(MPI)
    CALL MPI_ALLGATHER(i_pe2send_fluid,maxneigh,MPI_INTEGER,allneigh,maxneigh,MPI_INTEGER, MPI_COMM_WORLD, ierr)
#else
    allneigh = i_pe2send_fluid
#endif
    !determine from who I have to receive
    ltemp(1)=.false.
    j=0
    DO i=0,numprocs*maxneigh-1
       IF(allneigh(i)==myid) THEN
          IF(j.GE.maxneigh) THEN
             PRINT *,'Task',myid,'Invalid index for i_pe2recv_fluid',j
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
      write(6,*)'id=',idrank,'n_pe2send_fluid=',n_pe2send_fluid,i_pe2send_fluid(0:n_pe2send_fluid-1)
      write(6,*)'id=',idrank,'n_pe2recv_fluid=',n_pe2recv_fluid,i_pe2recv_fluid(0:n_pe2recv_fluid-1)
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

                IF(ownern(i4).NE.myid.AND.isfluid(i,jy,k)/=3) THEN   ! select nodes with ownership /= myid
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
#if defined(MPI)
    DO i=0, n_pe2recv_fluid-1
       CALL MPI_IRECV(n_pop2recv_fluid(i),1,MPI_INTEGER, &
            & i_pe2recv_fluid(i),i_pe2recv_fluid(i)+findntag, &
            & MPI_COMM_WORLD,request(i),ierr)
    ENDDO
    DO i=0, n_pe2send_fluid-1
       CALL MPI_SEND(n_pop2send_fluid(i),1,MPI_INTEGER, &
            & i_pe2send_fluid(i),myid+findntag, &
            & MPI_COMM_WORLD,ierr)
    ENDDO
    IF(n_pe2recv_fluid.GT.0) THEN
       CALL MPI_WAITALL(n_pe2recv_fluid,request,status,ierr)
    ENDIF
#endif
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
#if defined(MPI)
    DO i=0, n_pe2recv_fluid-1
       CALL MPI_IRECV(i_pop2recv_fluid(0,0,i),2*n_pop2recv_fluid(i),MYINT, &
            & i_pe2recv_fluid(i),i_pe2recv_fluid(i)+findntag, &
            & MPI_COMM_WORLD,request(i),ierr)
    ENDDO
    DO i=0, n_pe2send_fluid-1
       CALL MPI_SEND(i_pop2send_fluid(0,0,i),2*n_pop2send_fluid(i),MYINT, &
            & i_pe2send_fluid(i),myid+findntag, &
            & MPI_COMM_WORLD,ierr)
    ENDDO
    IF(n_pe2recv_fluid.GT.0) THEN
       CALL MPI_WAITALL(n_pe2recv_fluid,request,status,ierr)
    ENDIF
#endif
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
          if(ownern(i4temp).ne.myid) then
             write(50+myid,*)'Horror in task',myid,i4,i4temp,itemp,jtemp,&
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
  

  
  SUBROUTINE findneigh_hvar(nx,ny,nz,nbuff,ibctype,ixpbc,iypbc,izpbc,minx,maxx, &
   miny,maxy,minz,maxz)

    IMPLICIT NONE

    INTEGER :: request(0:maxneigh-1)
#if defined(MPI)
    INTEGER :: status(MPI_STATUS_SIZE,0:maxneigh-1)
#endif
    
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

    INTEGER :: allneigh(0:maxneigh*numprocs-1)

    INTEGER :: sender(0:numprocs-1)
    INTEGER :: receiver(0:numprocs-1)
    
    logical, parameter :: lverbose=.false.
    logical :: ltemp(1)
    
    if(numprocs==1)return
    
    DO i=0,numprocs-1
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
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
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
                  write(6,*)idrank,' something bad ',itemp2,jtemp2,ktemp2,itemp,jtemp,ktemp,minx,maxx
                  ltemp(1)=.true.
                endif
                !max                if(idrank.eq.0) then
                !max                   write(50,*)'i',itemp,'j',jtemp,'k',ktemp,'owner',ownern(i4)
                !max                endif
                !max                if(idrank.eq.1) then
                !max                   write(51,*)'i',itemp,'j',jtemp,'k',ktemp,'owner',ownern(i4)
                !max                endif  
                if(ownern(i4)<0 .or.ownern(i4)>(mxrank-1))then
                  write(6,*)idrank,' something bad ',itemp2,jtemp2,ktemp2,itemp,jtemp,ktemp,i4,ownern(i4)
                  ltemp(1)=.true.
                endif
                IF(ownern(i4).NE.myid) THEN
                   !sender(ownern(i4))=sender(ownern(i4))+1
                   receiver(ownern(i4))=receiver(ownern(i4))+1
                   !max                   write(0,*)'ownern ',ownern(i4),'for node ',i4,'ix=',i,'iy=',j,'iz=',k,'l=',l
                   !max                   write(0,*)'itemp=',itemp,'jtemp=',jtemp,'ktemp=',ktemp
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
          write(6,*)'id=',idrank,'receiver=',receiver(0:numprocs-1)
        endif
        call flush(6)
        call get_sync_world
      enddo
    endif
    
    !determine the max number of fluid nodes for node/PE that I have to receive 
    j=0
    maxreceiver=0
    DO i=0,numprocs-1

       IF(receiver(i).GT.0) THEN ! select procs that talks to this one

          IF(j.GE.maxneigh) THEN
             PRINT *,'Task',myid,'Invalid index for myreceiver',j
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
    
!max    write(0,*)'id=',myid,'maxreceiver=',maxreceiver,'n_pe2recv_fluid_hvar=',n_pe2recv_fluid_hvar
    !initialization
    DO i=0,maxneigh-1
       i_pe2send_fluid_hvar(i)=-1
    ENDDO
#if defined(MPI)
    CALL MPI_ALLGATHER(i_pe2recv_fluid_hvar,maxneigh,MPI_INTEGER,allneigh,maxneigh,MPI_INTEGER, MPI_COMM_WORLD, ierr)
#else
    allneigh = i_pe2recv_fluid_hvar
#endif
    
    !determine to who I have to send
    j=0
    DO i=0,numprocs*maxneigh-1
       IF(allneigh(i)==myid) THEN
          IF(j.GE.maxneigh) THEN
             PRINT *,'Task',myid,'Invalid index for i_pe2send_fluid_hvar',j
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
          write(6,*)'id=',idrank,'n_pe2send_fluid_hvar=',n_pe2send_fluid_hvar,i_pe2send_fluid_hvar(0:n_pe2send_fluid_hvar-1)
          write(6,*)'id=',idrank,'n_pe2recv_fluid_hvar=',n_pe2recv_fluid_hvar,i_pe2recv_fluid_hvar(0:n_pe2recv_fluid_hvar-1)
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

                IF(ownern(i4).NE.myid) THEN   ! select nodes with ownership /= myid
                   j=receiver(ownern(i4))  !find the unique ID of the cast operation
                   if(i4back(itemp,jtemp,ktemp)<0)then
                     write(6,*)'id',myid,' something bad ',itemp,jtemp,ktemp
                     ltemp(1)=.true.
                   endif
                   if(i4back(itemp2,jtemp2,ktemp2)<0)then
                     write(6,*)'id',myid,' csomething bad ',itemp2,jtemp2,ktemp2
                     ltemp(1)=.true.
                   endif
                   i_var2recv_fluid(0,n_var2recv_fluid(j),j)=i4back(itemp,jtemp,ktemp) !unique ID of the sending node
                   i_var2recv_fluid(1,n_var2recv_fluid(j),j)=i4back(itemp2,jtemp2,ktemp2) !unique ID of the receiving node
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
#if defined(MPI)
    DO i=0, n_pe2send_fluid_hvar-1
       CALL MPI_IRECV(n_var2send_fluid(i),1,MPI_INTEGER, &
            & i_pe2send_fluid_hvar(i),i_pe2send_fluid_hvar(i)+findntag_hvar, &
            & MPI_COMM_WORLD,request(i),ierr)
    ENDDO
    DO i=0, n_pe2recv_fluid_hvar-1
       CALL MPI_SEND(n_var2recv_fluid(i),1,MPI_INTEGER, &
            & i_pe2recv_fluid_hvar(i),myid+findntag_hvar, &
            & MPI_COMM_WORLD,ierr)
    ENDDO
    IF(n_pe2send_fluid_hvar.GT.0) THEN
       CALL MPI_WAITALL(n_pe2send_fluid_hvar,request,status,ierr)
    ENDIF
#endif
    
    if(lverbose)then
      call flush(6)
      call get_sync_world
      do i=0,mxrank-1
        if(i==idrank)then
          write(6,*)'id=',idrank,'i_pe2send_fluid_hvar=',i_pe2send_fluid_hvar(0:n_pe2send_fluid_hvar-1)
          write(6,*)'id=',idrank,'n_var2send_fluid=',n_var2send_fluid(0:n_pe2send_fluid_hvar-1)
        endif
        call flush(6)
        call get_sync_world
      enddo
      call flush(6)
      call get_sync_world
      do i=0,mxrank-1
        if(i==idrank)then
          write(6,*)'id=',idrank,'i_pe2recv_fluid_hvar=',i_pe2recv_fluid_hvar(0:n_pe2recv_fluid_hvar-1)
          write(6,*)'id=',idrank,'n_var2recv_fluid=',n_var2recv_fluid(0:n_pe2send_fluid_hvar-1)
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
#if defined(MPI)
    DO i=0, n_pe2send_fluid_hvar-1
       CALL MPI_IRECV(i_var2send_fluid(0,0,i),2*n_var2send_fluid(i),MYINT, &
            & i_pe2send_fluid_hvar(i),i_pe2send_fluid_hvar(i)+findntag_hvar, &
            & MPI_COMM_WORLD,request(i),ierr)
    ENDDO
    DO i=0, n_pe2recv_fluid_hvar-1
       CALL MPI_SEND(i_var2recv_fluid(0,0,i),2*n_var2recv_fluid(i),MYINT, &
            & i_pe2recv_fluid_hvar(i),myid+findntag_hvar, &
            & MPI_COMM_WORLD,ierr)
    ENDDO
    IF(n_pe2recv_fluid.GT.0) THEN
       CALL MPI_WAITALL(n_pe2send_fluid_hvar,request,status,ierr)
    ENDIF
#endif
    
    call get_sync_world
    
    ltemp(1)=.false.
    DO i=0, n_pe2send_fluid_hvar-1
       DO j=0,n_var2send_fluid(i)-1
          i4_send=i_var2send_fluid(0,j,i)
          i4_recv=i_var2send_fluid(1,j,i)
          if(i4_send<1)then
            write(6,*)idrank,i,j,n_pe2send_fluid_hvar,n_var2send_fluid(i),i_pe2send_fluid_hvar(i)
          endif
          if(ownern(i4_send).NE.myid)then
            write(6,*)'idrank = ',idrank,&
             'Horror in task hvar. I have to send ',i4_send,' belonging to ',ownern(i4_send)
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
          if(ownern(i4_send).EQ.myid)then
            write(6,*)'idrank = ',idrank,&
             'Horror in task hvar. I have to receiving ',i4_send,' belonging to me ',ownern(i4_send)
            ltemp(1)=.true.
          endif
        ENDDO
    ENDDO
    call or_world_larr(ltemp,1)
    if(ltemp(1))then
      call finalize_world
      stop
    endif
    
    

  END SUBROUTINE findneigh_hvar
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if defined(MPI)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE commspop(pop_ptr)

  IMPLICIT NONE

  ! REAL(KIND=PRC) :: pop_ptr(0:,0:)
  ! REAL(KIND=PRC) :: pop_ptr(0:,0:)
  !max REAL(KIND=PRC) :: pop_ptr(0:,0:)
  TYPE(REALPTR), dimension(0:links):: pop_ptr
!max  REAL(KIND=PRC) :: pop_ptr(0:npop-1,0:*)
  INTEGER :: request(0:maxneigh-1)
  INTEGER :: status(MPI_STATUS_SIZE,0:maxneigh-1)
  INTEGER :: i,j,mym(3)
  INTEGER(kind=IPRC) :: i4, itemp, jtemp, ktemp
  
  if(numprocs==1)return
  
   if(firstmove) then
      do i=0, n_pe2recv_fluid-1
         CALL MPI_IRECV(buffpopr(0,i),n_pop2recv_fluid(i),MYFLOAT, &
                     &  i_pe2recv_fluid(i),i_pe2recv_fluid(i)+movetag, &
                        MPI_COMM_WORLD,requestm(i),ierr)
      enddo
      firstmove=.false.
   endif
   do i=0, n_pe2send_fluid-1
      do j=0,n_pop2send_fluid(i)-1
         !max         buffpops(j,i)=pop_ptr(i_pop2send_fluid(0,j,i),i_pop2send_fluid(1,j,i))
         i4=i_pop2send_fluid(1,j,i)
         !ktemp = i4/INT(nxy2,KIND=IPRC)
         !jtemp = (i4 - ktemp*INT(nxy2,KIND=IPRC))/INT(nx2,KIND=IPRC)
         !itemp = i4 - ktemp*INT(nxy2,KIND=IPRC) - jtemp*INT(nx2,KIND=IPRC)
         mym=i4find(i4)
         ktemp = mym(3)
         jtemp = mym(2)
         itemp = mym(1)
!max         write(50+myid,*)'id=',myid,'l=',i_pop2send_fluid(0,j,i),'ix=',itemp,'iy=',jtemp,'iz=',ktemp
         buffpops(j,i)=pop_ptr(i_pop2send_fluid(0,j,i))%p(itemp,jtemp,ktemp)
      enddo
      if(syncsend.eq.1) then
         CALL MPI_SEND(buffpops(0,i),n_pop2send_fluid(i),MYFLOAT, &
                     & i_pe2send_fluid(i),myid+movetag,MPI_COMM_WORLD,ierr)
      else
         CALL MPI_ISEND(buffpops(0,i),n_pop2send_fluid(i),MYFLOAT, &
                     & i_pe2send_fluid(i),myid+movetag,MPI_COMM_WORLD,request(i),ierr)
      endif

   enddo

   END SUBROUTINE commspop
 
   SUBROUTINE commrpop(pop_ptr)

  IMPLICIT NONE

  ! REAL(KIND=PRC) :: pop_ptr(0:,0:)
  ! REAL(KIND=PRC) :: pop_ptr(0:,0:)
  !max REAL(KIND=PRC) :: pop_ptr(0:,0:)
  TYPE(REALPTR), dimension(0:links):: pop_ptr
!max  REAL(KIND=PRC) :: pop_ptr(0:npop-1,0:*)
  INTEGER :: request(0:maxneigh-1)
  INTEGER :: status(MPI_STATUS_SIZE,0:maxneigh-1)
  INTEGER :: i,j,mym(3)
  INTEGER(kind=IPRC) :: i4, itemp, jtemp, ktemp
  
  
  if(numprocs==1)return
  
   if(n_pe2recv_fluid.gt.0) then
      CALL MPI_WAITALL(n_pe2recv_fluid,requestm,status,ierr)
   endif
   do i=0, n_pe2recv_fluid-1
      do j=0,n_pop2recv_fluid(i)-1
         i4=i_pop2recv_fluid(1,j,i)
         !ktemp = i4/INT(nxy2,KIND=IPRC)
         !jtemp = (i4 - ktemp*INT(nxy2,KIND=IPRC))/INT(nx2,KIND=IPRC)
         !itemp = i4 - ktemp*INT(nxy2,KIND=IPRC) - jtemp*INT(nx2,KIND=IPRC)
         mym=i4find(i4)
         ktemp = mym(3)
         jtemp = mym(2)
         itemp = mym(1)
         pop_ptr(i_pop2recv_fluid(0,j,i))%p(itemp,jtemp,ktemp)=buffpopr(j,i)
      enddo
   enddo

   do i=0, n_pe2recv_fluid-1
      CALL MPI_IRECV(buffpopr(0,i),n_pop2recv_fluid(i),MYFLOAT, &
                   & i_pe2recv_fluid(i),i_pe2recv_fluid(i)+movetag, &
                     MPI_COMM_WORLD,requestm(i),ierr)
   enddo

   if(syncsend.eq.0.and.n_pe2send_fluid.gt.0) then
      CALL MPI_WAITALL(n_pe2send_fluid,request,status,ierr)
   endif

   END SUBROUTINE commrpop
   
   SUBROUTINE commexch_dens(dtemp,dtemp2)

  IMPLICIT NONE

  ! REAL(KIND=PRC) :: pop_ptr(0:,0:)
  ! REAL(KIND=PRC) :: pop_ptr(0:,0:)
  !max REAL(KIND=PRC) :: pop_ptr(0:,0:)
  REAL(KIND=PRC), dimension(:,:,:), allocatable :: dtemp,dtemp2
!max  REAL(KIND=PRC) :: pop_ptr(0:npop-1,0:*)
  INTEGER :: request(0:maxneigh-1)
  INTEGER :: status(MPI_STATUS_SIZE,0:maxneigh-1)
  INTEGER :: i,j,mym(3)
  INTEGER(kind=IPRC) :: i4, itemp, jtemp, ktemp

  if(numprocs==1)return
  
  do i=0, n_pe2recv_fluid_hvar-1
    CALL MPI_IRECV(buffr_hvar(0,i),n_var2recv_fluid(i),MYFLOAT, &
      &  i_pe2recv_fluid_hvar(i),i_pe2recv_fluid_hvar(i)+tag_hvar, &
        MPI_COMM_WORLD,requestm_hvar(i),ierr)
  enddo
  
  


   do i=0, n_pe2send_fluid_hvar-1
      do j=0,n_var2send_fluid(i)-1
         !max         buffpops(j,i)=pop_ptr(i_pop2send_fluid(0,j,i),i_pop2send_fluid(1,j,i))
         i4=i_var2send_fluid(0,j,i)
         !ktemp = i4/INT(nxy2,KIND=IPRC)
         !jtemp = (i4 - ktemp*INT(nxy2,KIND=IPRC))/INT(nx2,KIND=IPRC)
         !itemp = i4 - ktemp*INT(nxy2,KIND=IPRC) - jtemp*INT(nx2,KIND=IPRC)
         mym=i4find(i4)
         ktemp = mym(3)
         jtemp = mym(2)
         itemp = mym(1)
!max         write(50+myid,*)'id=',myid,'l=',i_pop2send_fluid(0,j,i),'ix=',itemp,'iy=',jtemp,'iz=',ktemp
         buffs_hvar(j,i)=dtemp(itemp,jtemp,ktemp)
         !buffpops(j,i)=pop_ptr(i_pop2send_fluid(0,j,i))%p(itemp,jtemp,ktemp)
      enddo

      CALL MPI_ISEND(buffs_hvar(0,i),n_var2send_fluid(i),MYFLOAT, &
       & i_pe2send_fluid_hvar(i),myid+tag_hvar,MPI_COMM_WORLD,request_hvar(i),ierr)
     

   enddo
   
   if(ldo_second)then
    do i=0, n_pe2recv_fluid_hvar-1
      CALL MPI_IRECV(buffr_hvar2(0,i),n_var2recv_fluid(i),MYFLOAT, &
        &  i_pe2recv_fluid_hvar(i),i_pe2recv_fluid_hvar(i)+tag_hvar2, &
          MPI_COMM_WORLD,requestm_hvar2(i),ierr)
    enddo
    do i=0, n_pe2send_fluid_hvar-1
      do j=0,n_var2send_fluid(i)-1
         !max         buffpops(j,i)=pop_ptr(i_pop2send_fluid(0,j,i),i_pop2send_fluid(1,j,i))
         i4=i_var2send_fluid(0,j,i)
         !ktemp = i4/INT(nxy2,KIND=IPRC)
         !jtemp = (i4 - ktemp*INT(nxy2,KIND=IPRC))/INT(nx2,KIND=IPRC)
         !itemp = i4 - ktemp*INT(nxy2,KIND=IPRC) - jtemp*INT(nx2,KIND=IPRC)
         mym=i4find(i4)
         ktemp = mym(3)
         jtemp = mym(2)
         itemp = mym(1)
!max         write(50+myid,*)'id=',myid,'l=',i_pop2send_fluid(0,j,i),'ix=',itemp,'iy=',jtemp,'iz=',ktemp
         buffs_hvar2(j,i)=dtemp2(itemp,jtemp,ktemp)
         !buffpops(j,i)=pop_ptr(i_pop2send_fluid(0,j,i))%p(itemp,jtemp,ktemp)
      enddo

      CALL MPI_ISEND(buffs_hvar2(0,i),n_var2send_fluid(i),MYFLOAT, &
       & i_pe2send_fluid_hvar(i),myid+tag_hvar2,MPI_COMM_WORLD,request_hvar2(i),ierr)
   enddo
  endif

   END SUBROUTINE commexch_dens
   
  SUBROUTINE commexch_vel_component(dtemp)

  IMPLICIT NONE

  ! REAL(KIND=PRC) :: pop_ptr(0:,0:)
  ! REAL(KIND=PRC) :: pop_ptr(0:,0:)
  !max REAL(KIND=PRC) :: pop_ptr(0:,0:)
  REAL(KIND=PRC), dimension(:,:,:), allocatable :: dtemp
!max  REAL(KIND=PRC) :: pop_ptr(0:npop-1,0:*)
  INTEGER :: request(0:maxneigh-1)
  INTEGER :: status(MPI_STATUS_SIZE,0:maxneigh-1)
  INTEGER :: i,j,mym(3)
  INTEGER(kind=IPRC) :: i4, itemp, jtemp, ktemp

  if(numprocs==1)return
  
  do i=0, n_pe2recv_fluid_hvar-1
    CALL MPI_IRECV(buffr_hvar(0,i),n_var2recv_fluid(i),MYFLOAT, &
      &  i_pe2recv_fluid_hvar(i),i_pe2recv_fluid_hvar(i)+tag_hvar, &
        MPI_COMM_WORLD,requestm_hvar(i),ierr)
  enddo
  
  


   do i=0, n_pe2send_fluid_hvar-1
      do j=0,n_var2send_fluid(i)-1
         !max         buffpops(j,i)=pop_ptr(i_pop2send_fluid(0,j,i),i_pop2send_fluid(1,j,i))
         i4=i_var2send_fluid(0,j,i)
         !ktemp = i4/INT(nxy2,KIND=IPRC)
         !jtemp = (i4 - ktemp*INT(nxy2,KIND=IPRC))/INT(nx2,KIND=IPRC)
         !itemp = i4 - ktemp*INT(nxy2,KIND=IPRC) - jtemp*INT(nx2,KIND=IPRC)
         mym=i4find(i4)
         ktemp = mym(3)
         jtemp = mym(2)
         itemp = mym(1)
!max         write(50+myid,*)'id=',myid,'l=',i_pop2send_fluid(0,j,i),'ix=',itemp,'iy=',jtemp,'iz=',ktemp
         buffs_hvar(j,i)=dtemp(itemp,jtemp,ktemp)
         !buffpops(j,i)=pop_ptr(i_pop2send_fluid(0,j,i))%p(itemp,jtemp,ktemp)
      enddo

      CALL MPI_ISEND(buffs_hvar(0,i),n_var2send_fluid(i),MYFLOAT, &
       & i_pe2send_fluid_hvar(i),myid+tag_hvar,MPI_COMM_WORLD,request_hvar(i),ierr)
     

   enddo
   

   END SUBROUTINE commexch_vel_component
 
   SUBROUTINE commwait_dens(dtemp,dtemp2)

  IMPLICIT NONE

  ! REAL(KIND=PRC) :: pop_ptr(0:,0:)
  ! REAL(KIND=PRC) :: pop_ptr(0:,0:)
  !max REAL(KIND=PRC) :: pop_ptr(0:,0:)
  !TYPE(REALPTR), dimension(0:links):: pop_ptr
  REAL(KIND=PRC), dimension(:,:,:), allocatable :: dtemp,dtemp2
!max  REAL(KIND=PRC) :: pop_ptr(0:npop-1,0:*)
  INTEGER :: request(0:maxneigh-1)
  INTEGER :: status(MPI_STATUS_SIZE,0:maxneigh-1)
  INTEGER :: i,j,mym(3)
  INTEGER(kind=IPRC) :: i4, itemp, jtemp, ktemp
  
  if(numprocs==1)return

   if(n_pe2recv_fluid_hvar.gt.0) then
      CALL MPI_WAITALL(n_pe2recv_fluid_hvar,requestm_hvar,status,ierr)
   endif
   do i=0, n_pe2recv_fluid_hvar-1
      do j=0,n_var2recv_fluid(i)-1
         i4=i_var2recv_fluid(1,j,i)
         !ktemp = i4/INT(nxy2,KIND=IPRC)
         !jtemp = (i4 - ktemp*INT(nxy2,KIND=IPRC))/INT(nx2,KIND=IPRC)
         !itemp = i4 - ktemp*INT(nxy2,KIND=IPRC) - jtemp*INT(nx2,KIND=IPRC)
         mym=i4find(i4)
         ktemp = mym(3)
         jtemp = mym(2)
         itemp = mym(1)
         !pop_ptr(i_pop2recv_fluid(0,j,i))%p(itemp,jtemp,ktemp)=buffpopr(j,i)
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
         !ktemp = i4/INT(nxy2,KIND=IPRC)
         !jtemp = (i4 - ktemp*INT(nxy2,KIND=IPRC))/INT(nx2,KIND=IPRC)
         !itemp = i4 - ktemp*INT(nxy2,KIND=IPRC) - jtemp*INT(nx2,KIND=IPRC)
         mym=i4find(i4)
         ktemp = mym(3)
         jtemp = mym(2)
         itemp = mym(1)
         !pop_ptr(i_pop2recv_fluid(0,j,i))%p(itemp,jtemp,ktemp)=buffpopr(j,i)
         dtemp2(itemp,jtemp,ktemp)=buffr_hvar2(j,i)
       enddo
     enddo
   
     if(n_pe2send_fluid_hvar.gt.0) then
       CALL MPI_WAITALL(n_pe2send_fluid_hvar,request_hvar2,status,ierr)
     endif
   endif

   END SUBROUTINE commwait_dens
   
   SUBROUTINE commwait_vel_component(dtemp)

  IMPLICIT NONE

  ! REAL(KIND=PRC) :: pop_ptr(0:,0:)
  ! REAL(KIND=PRC) :: pop_ptr(0:,0:)
  !max REAL(KIND=PRC) :: pop_ptr(0:,0:)
  !TYPE(REALPTR), dimension(0:links):: pop_ptr
  REAL(KIND=PRC), dimension(:,:,:), allocatable :: dtemp
!max  REAL(KIND=PRC) :: pop_ptr(0:npop-1,0:*)
  INTEGER :: request(0:maxneigh-1)
  INTEGER :: status(MPI_STATUS_SIZE,0:maxneigh-1)
  INTEGER :: i,j,mym(3)
  INTEGER(kind=IPRC) :: i4, itemp, jtemp, ktemp
  
  if(numprocs==1)return

   if(n_pe2recv_fluid_hvar.gt.0) then
      CALL MPI_WAITALL(n_pe2recv_fluid_hvar,requestm_hvar,status,ierr)
   endif
   do i=0, n_pe2recv_fluid_hvar-1
      do j=0,n_var2recv_fluid(i)-1
         i4=i_var2recv_fluid(1,j,i)
         !ktemp = i4/INT(nxy2,KIND=IPRC)
         !jtemp = (i4 - ktemp*INT(nxy2,KIND=IPRC))/INT(nx2,KIND=IPRC)
         !itemp = i4 - ktemp*INT(nxy2,KIND=IPRC) - jtemp*INT(nx2,KIND=IPRC)
         mym=i4find(i4)
         ktemp = mym(3)
         jtemp = mym(2)
         itemp = mym(1)
         !pop_ptr(i_pop2recv_fluid(0,j,i))%p(itemp,jtemp,ktemp)=buffpopr(j,i)
         dtemp(itemp,jtemp,ktemp)=buffr_hvar(j,i)
      enddo
   enddo
   
   if(n_pe2send_fluid_hvar.gt.0) then
      CALL MPI_WAITALL(n_pe2send_fluid_hvar,request_hvar,status,ierr)
   endif
   

   END SUBROUTINE commwait_vel_component
   
   SUBROUTINE comm_hvar(dtemp,dtemp2)
   
   IMPLICIT NONE
   
   REAL(KIND=PRC), dimension(:,:,:), allocatable :: dtemp,dtemp2
   
   if(numprocs==1)return
   
   call commexch_dens(dtemp,dtemp2)
   call commwait_dens(dtemp,dtemp2)
   
   
   return
   
   END SUBROUTINE comm_hvar
   
   
   
   
  SUBROUTINE comm_init_isfluid(temp)

  IMPLICIT NONE

  INTEGER(kind=1), dimension(:,:,:), allocatable :: temp
  
  INTEGER :: request_sub(0:maxneigh-1)
  INTEGER :: requestm_sub(0:maxneigh-1)
  INTEGER :: status(MPI_STATUS_SIZE,0:maxneigh-1)
  INTEGER :: i,j,mym(3)
  INTEGER(kind=IPRC) :: i4, itemp, jtemp, ktemp

  
  if(numprocs==1)return
  
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
       & i_pe2send_fluid_hvar(i),myid+tag_isfluid,MPI_COMM_WORLD,request_sub(i),ierr)
     

   enddo

   if(n_pe2recv_fluid_hvar.gt.0) then
      CALL MPI_WAITALL(n_pe2recv_fluid_hvar,requestm_sub,status,ierr)
   endif
   do i=0, n_pe2recv_fluid_hvar-1
      do j=0,n_var2recv_fluid(i)-1
         i4=i_var2recv_fluid(1,j,i)
         !ktemp = i4/INT(nxy2,KIND=IPRC)
         !jtemp = (i4 - ktemp*INT(nxy2,KIND=IPRC))/INT(nx2,KIND=IPRC)
         !itemp = i4 - ktemp*INT(nxy2,KIND=IPRC) - jtemp*INT(nx2,KIND=IPRC)
         mym=i4find(i4)
         ktemp = mym(3)
         jtemp = mym(2)
         itemp = mym(1)
         !pop_ptr(i_pop2recv_fluid(0,j,i))%p(itemp,jtemp,ktemp)=buffpopr(j,i)
         temp(itemp,jtemp,ktemp)=buffr_isfluid(j,i)
      enddo
   enddo
   
   if(n_pe2send_fluid_hvar.gt.0) then
      CALL MPI_WAITALL(n_pe2send_fluid_hvar,request_sub,status,ierr)
   endif
   
   END SUBROUTINE comm_init_isfluid
   
#else

  SUBROUTINE commspop(pop_ptr)

  IMPLICIT NONE

  ! REAL(KIND=PRC) :: pop_ptr(0:,0:)
  ! REAL(KIND=PRC) :: pop_ptr(0:,0:)
  !max REAL(KIND=PRC) :: pop_ptr(0:,0:)
  TYPE(REALPTR), dimension(0:links):: pop_ptr
  
  return
  
  end subroutine commspop 
  
  
  SUBROUTINE commrpop(pop_ptr)

  IMPLICIT NONE

  TYPE(REALPTR), dimension(0:links):: pop_ptr
  
  return
  
  end subroutine commrpop
  
  SUBROUTINE comm_hvar(dtemp)
   
   IMPLICIT NONE
   
   REAL(KIND=PRC), dimension(:,:,:), allocatable :: dtemp
   
   return
   
   END SUBROUTINE comm_hvar
   
   SUBROUTINE commexch_dens(dtemp)
   
   IMPLICIT NONE
   
   REAL(KIND=PRC), dimension(:,:,:), allocatable :: dtemp
   
   return
   
   END SUBROUTINE commexch_dens
   
    SUBROUTINE commwait_dens(dtemp)
   
   IMPLICIT NONE
   
   REAL(KIND=PRC), dimension(:,:,:), allocatable :: dtemp
   
   return
   
   END SUBROUTINE commwait_dens
   
   SUBROUTINE commexch_vel_component(dtemp)
   
   IMPLICIT NONE
   
   REAL(KIND=PRC), dimension(:,:,:), allocatable :: dtemp
   
   return
   
   END SUBROUTINE commexch_vel_component
   
    SUBROUTINE commwait_vel_component(dtemp)
   
   IMPLICIT NONE
   
   REAL(KIND=PRC), dimension(:,:,:), allocatable :: dtemp
   
   return
   
   END SUBROUTINE commwait_vel_component
   
  SUBROUTINE comm_init_isfluid(temp)

  IMPLICIT NONE

  INTEGER(kind=1), dimension(:,:,:), allocatable :: temp
  
  return
  
  END SUBROUTINE comm_init_isfluid
  
#endif


   END MODULE lbempi_mod
