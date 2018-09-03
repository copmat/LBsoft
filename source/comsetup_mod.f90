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
#if defined(CUDA) && defined(USEPINME)
! Interface to cudaMallocHost and cudaFree
MODULE cuda_alloc

  USE iso_c_binding

  IMPLICIT NONE
  ! Define the C pointer as type (C_PTR)
  TYPE(C_PTR) :: cptr_locu, cptr_locv, cptr_locw, cptr_locrho
  TYPE(C_PTR) :: cptr_locju, cptr_locjv, cptr_locjw
  TYPE(C_PTR) :: cptr_locfu, cptr_locfv, cptr_locfw

  INTERFACE

     ! cudaMallocHost
     INTEGER (C_INT) FUNCTION cudaMallocHost(buffer, size)  bind(C,name="cudaMallocHost")
       USE iso_c_binding
       IMPLICIT NONE
       TYPE (C_PTR)  :: buffer
       INTEGER (C_INT), value :: size
     END FUNCTION cudaMallocHost

     ! cudaFreeHost
     INTEGER (C_INT) FUNCTION cudaFreeHost(buffer)  bind(C,name="cudaFreeHost")
       USE iso_c_binding
       IMPLICIT NONE
       TYPE (C_PTR), value :: buffer
     END FUNCTION cudaFreeHost
  END INTERFACE

END MODULE cuda_alloc

#endif

MODULE lbempi_mod
  use version_mod, only : idrank, mxrank
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
  INTEGER, SAVE :: MYCOMM=MPI_COMM_WORLD

#endif
!!! temporary parameter
  
!!!

  INTEGER, PARAMETER :: maxneigh=32
  INTEGER, PARAMETER :: npop=19
  INTEGER, PARAMETER :: movetag=393000, &
           & findntag=524000,findptag=655000,findwtag=786000
  
  integer rc, ierr
  integer,PUBLIC :: lnofsite, itpppsize

  INTEGER :: domdec, syncsend
  INTEGER,SAVE :: requestm(0:maxneigh-1), requestp(0:maxneigh-1), &
       requestw(0:maxneigh-1), requesti(0:maxneigh-1), requesto(0:maxneigh-1), &
       requesth(0:maxneigh-1), requesth_wall(0:maxneigh-1), requests(0:maxneigh-1)
  LOGICAL :: firstmove=.true.
  LOGICAL :: allpbc

  INTEGER :: n_pe2recv_fluid, n_pe2send_fluid
  INTEGER :: i_pe2send_fluid(0:maxneigh-1), n_pop2send_fluid(0:maxneigh-1)
  INTEGER :: i_pe2recv_fluid(0:maxneigh-1), n_pop2recv_fluid(0:maxneigh-1)

  INTEGER :: nprocz, nprocy, nprocx
  INTEGER :: nxy2, nx2

  integer, POINTER :: countnpp(:)
  INTEGER(KIND=2), ALLOCATABLE :: ownern(:)

  INTEGER(kind=IPRC), ALLOCATABLE :: i_pop2send_fluid(:,:,:), i_pop2recv_fluid(:,:,:)
  INTEGER(kind=IPRC), ALLOCATABLE :: i_pop2send_wall(:,:,:), i_pop2recv_wall(:,:,:)
  INTEGER(kind=IPRC), ALLOCATABLE :: i_pop2send_inlet(:,:,:), i_pop2recv_inlet(:,:,:)
  INTEGER(kind=IPRC), ALLOCATABLE :: i_pop2send_outlet(:,:,:), i_pop2recv_outlet(:,:,:)

  real(kind=PRC), POINTER :: buffpops(:,:), buffpopr(:,:)
  real(kind=PRC), POINTER :: bufftags(:,:), bufftagr(:,:)
  
  public :: commspop, commrpop
  PUBLIC :: buff_transf_r,buff_transf_s
  REAL(kind=PRC),ALLOCATABLE,DIMENSION(:,:) :: buff_transf_r,buff_transf_s
#if defined(FNMD)
    INTEGER, PARAMETER :: mddirnpop=27
    INTEGER :: mddir(0:mddirnpop-npop-1)
#endif

CONTAINS
#define LARGEINT 1073741824
SUBROUTINE setupcom(nx,ny,nz,ibctype,ixpbc,iypbc,izpbc,minx,maxx, &
   miny,maxy,minz,maxz)
  
 
  
  INTEGER, intent(in) :: nx,ny,nz,ibctype
  INTEGER, intent(inout) :: ixpbc,iypbc,izpbc,minx,maxx,miny,maxy, &
   minz,maxz
  
  INTEGER :: nbuff=2
  INTEGER :: i
  INTEGER :: ownernlb,ownernub

  domdec=7
  myid=idrank
  numprocs=mxrank
  syncsend=1

  nxy2=(nx+2*nbuff)*(ny+2*nbuff)
  nx2=nx+(2*nbuff)
  ALLOCATE(countnpp(0:numprocs-1))
  DO i=0,numprocs-1
     countnpp(i)=0
  ENDDO
  ownernlb=(1-nbuff)+((1-nbuff)*nx2)+((1-nbuff)*nxy2)
  ownernub=(nx+2*nbuff)*(ny+2*nbuff)*(nz+2*nbuff)
!max  write(0,*),'NX=',nx,'NY=',ny,'NZ=',nz,'ownern lb=',ownernlb,'ownern up=',ownernub
  ALLOCATE(ownern(ownernlb:ownernub))
  do i=ownernlb,ownernub
     ownern(i)=-1
  enddo
  call cartdeco(nx,ny,nz,.true.,ownern,ownernlb,minx,maxx,miny,maxy, &
   minz,maxz)
!max  write(0,*)'countnpp(',myid,')=',countnpp(myid)
  call findneigh(nx,ny,nz,ibctype,ixpbc,iypbc,izpbc,minx,maxx,miny,maxy, &
   minz,maxz)

  DEALLOCATE(ownern)
END SUBROUTINE setupcom

  SUBROUTINE cartdeco(nx,ny,nz,first, ownern, ownernlb,minx,maxx,miny,maxy, &
   minz,maxz)

    IMPLICIT NONE
    
    integer, intent(in) :: nx,ny,nz
    LOGICAL,INTENT(in) :: first
    ! INTEGER(kind=IPRC) :: tempown(0:*)
    ! INTEGER :: order(0:*)
    INTEGER,INTENT(in):: ownernlb
    INTEGER(kind=2) :: ownern(ownernlb:*)
    INTEGER, intent(inout) :: minx,maxx,miny,maxy,minz,maxz

    INTEGER(kind=IPRC) :: i4
    INTEGER :: i,j,k
    INTEGER :: ix,iy,iz
    INTEGER :: currproc, currpoff=0, curroffz, curroffy, curroffx
    INTEGER :: slice=0, slicez=0, slicey=0, slicex=0
    INTEGER :: m, iandex

    INTEGER :: nbuff=2
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
       slice=(nz+2*nbuff)/numprocs
       IF(MOD((nz+2*nbuff),numprocs).NE.0) THEN
          slice=slice+1
       ENDIF
    ENDIF
    IF(domdec.EQ.7) THEN      !along zyx
       slicez=(nz+2*nbuff)/nprocz
       IF(MOD((nz+2*nbuff),nprocz).NE.0) THEN
          slicez=slicez+1
       ENDIF
    ENDIF
    DO iz=1-nbuff, nz+nbuff
       IF(domdec.EQ.1) THEN   !along z
          IF(slice.EQ.1.AND.iz.LT.(nz+nbuff)) THEN
             currproc=currproc+1
             slice=((nz+2*nbuff)-(iz+nbuff))/(numprocs-currproc)
             IF(MOD((nz+2*nbuff)-(iz+nbuff),(numprocs-currproc)).NE.0) THEN
                slice=slice+1
             ENDIF
          ELSE
             slice=slice-1
          ENDIF
       ELSE IF(domdec.EQ.4.OR.domdec.EQ.5) THEN !along zy and zx
          IF(iz.GT.(nz+nbuff)/2) THEN
             currpoff=numprocs/2
          ELSE
             currpoff=0
          ENDIF
          IF(domdec.EQ.4) THEN !along zy
             currproc=currpoff
             slice=(ny+2*nbuff)/(numprocs/2)
             IF(MOD((ny+2*nbuff),(numprocs/2)).NE.0) THEN
                slice=slice+1
             ENDIF
          ENDIF
       ELSE IF(domdec.EQ.7) THEN !along zyx
          IF(slicez.EQ.0.AND.iz.LT.(nz+nbuff)) THEN
             curroffz=curroffz+1
             slicez=((nz+2*nbuff)-(iz))/(nprocz-curroffz)
             IF(MOD((nz+2*nbuff)-(iz),(nprocz-curroffz)).NE.0) THEN
                slicez=slicez+1
             ENDIF
          ENDIF
          slicez=slicez-1
          curroffy=0
          slicey=(ny+2*nbuff)/nprocy
          IF(MOD((ny+2*nbuff),nprocy).NE.0) THEN
             slicey=slicey+1
          ENDIF
       ELSE IF(domdec.EQ.2) THEN !along y
          currproc=0
          slice=(ny+2*nbuff)/numprocs
          IF(MOD((ny+2*nbuff),numprocs).NE.0) THEN
             slice=slice+1
          ENDIF
       ENDIF
       DO iy=1-nbuff, ny+nbuff
          IF(domdec.EQ.2) THEN
             IF(slice.EQ.1.AND.iy.LT.(ny+nbuff)) THEN
                currproc=currproc+1
                slice=((ny+2*nbuff)-(iy+nbuff))/(numprocs-currproc)
                IF(MOD((ny+2*nbuff)-(iy+nbuff),(numprocs-currproc)).NE.0) THEN
                   slice=slice+1
                ENDIF
             ELSE
                slice=slice-1
             ENDIF
          ELSE IF(domdec.EQ.4) THEN !along zy
             IF(slice.EQ.1.AND.iy.LT.(ny+nbuff)) THEN
                currproc=currproc+1
                slice=((ny+2*nbuff)-(iy+nbuff))/ &
                     ((numprocs/2)-(currproc-currpoff))
                IF(MOD((ny+2*nbuff)-(iy+nbuff),&
                     (numprocs/2)-(currproc-currpoff)).NE.0) THEN
                   slice=slice+1
                ENDIF
             ELSE
                slice=slice-1
             ENDIF
          ELSE IF(domdec.EQ.6) THEN !along yx
             IF(iy.GT.(ny+nbuff)/2) THEN
                currpoff=numprocs/2
             ELSE
                currpoff=0
             ENDIF
             currproc=currpoff
             slice=(nx+2*nbuff)/(numprocs/2)
             IF(MOD((nx+2*nbuff),(numprocs/2)).NE.0) THEN
                slice=slice+1
             ENDIF
          ELSE IF(domdec.EQ.7) THEN !along zyx
             IF(slicey.EQ.0.AND.iy.LT.(ny+nbuff)) THEN
                curroffy=curroffy+1
                slicey=((ny+2*nbuff)-(iy))/(nprocy-curroffy)
                IF(MOD((ny+2*nbuff)-(iy),(nprocy-curroffy)).NE.0) THEN
                   slicey=slicey+1
                ENDIF
             ENDIF
             slicey=slicey-1
             curroffx=0
             slicex=(nx+2*nbuff)/nprocx
             IF(MOD((nx+2*nbuff),nprocx).NE.0) THEN
                slicex=slicex+1
             ENDIF
          ELSE IF(domdec.EQ.3) THEN !along x
             currproc=0
             slice=(nx+2*nbuff)/numprocs
             IF(MOD((nx+2*nbuff),numprocs).NE.0) THEN
                slice=slice+1
             ENDIF
          ELSE IF(domdec.EQ.5) THEN !along zx
             currproc=currpoff
             slice=(nx+2*nbuff)/(numprocs/2)
             IF(MOD((nx+2*nbuff),(numprocs/2)).NE.0) THEN
                slice=slice+1
             ENDIF
          ENDIF
          DO ix=1-nbuff, nx+nbuff
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
                IF(slice.EQ.1.AND.ix.LT.(nx+nbuff)) THEN
                   currproc=currproc+1
                   slice=((nx+2*nbuff)-(ix+nbuff))/(numprocs-currproc)
                   IF(MOD((nx+2*nbuff)-(ix+nbuff),(numprocs-currproc)).NE.0) THEN
                      slice=slice+1
                   ENDIF
                ELSE
                   slice=slice-1
                ENDIF
             ELSE IF(domdec.EQ.5.OR.domdec.EQ.6) THEN !along zx or along yx
                IF(slice.EQ.1.AND.ix.LT.(nx+nbuff)) THEN
                   currproc=currproc+1
                   slice=((nx+2*nbuff)-(ix+nbuff))/ &
                        ((numprocs/2)-(currproc-currpoff))
                   IF(MOD((nx+2*nbuff)-(ix+nbuff),&
                        (numprocs/2)-(currproc-currpoff)).NE.0) THEN
                      slice=slice+1
                   ENDIF
                ELSE
                   slice=slice-1
                ENDIF
             ELSE IF(domdec.EQ.7) THEN !along zyx
                IF(slicex.EQ.1.AND.ix.LT.(nx+nbuff)) THEN
                   curroffx=curroffx+1
                   slicex=((nx+2*nbuff)-(ix+nbuff))/ &
                        (nprocx-curroffx)
                   IF(MOD((nx+2*nbuff)-(ix+nbuff),(nprocx-curroffx)).NE.0) THEN
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
     out = INT(k,kind=IPRC)*nxy2 + INT(j,kind=IPRC)*nx2 + INT(i,kind=IPRC)
  END FUNCTION i4back
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE findneigh(nx,ny,nz,ibctype,ixpbc,iypbc,izpbc,minx,maxx, &
   miny,maxy,minz,maxz)

#if defined(CUDA) && defined(USEPINME)
    ! USE iso_c_binding
    ! USE cuda_alloc
    IMPLICIT NONE
    
    
    
    ! Define the floating point kind to be  single_precision
    INTEGER, PARAMETER :: fp_kind = KIND(0.0)
    INTEGER :: res
    ! The allocation is performed by C function call.
#else
    IMPLICIT NONE
#endif
    INTEGER :: request(0:maxneigh-1)
#if defined(MPI)
    INTEGER :: status(MPI_STATUS_SIZE,0:maxneigh-1)
#endif
    
    INTEGER, intent(in) :: nx,ny,nz,ibctype
    INTEGER, intent(inout) :: ixpbc,iypbc,izpbc,minx,maxx,miny,maxy, &
   minz,maxz
    
    INTEGER :: lminx, lmaxx, lminy, lmaxy, lminz, lmaxz
    INTEGER(kind=IPRC) :: i4,i4temp,i4shift
    INTEGER :: i,j,jy,k,ip,maxsender,maxreceiver
    INTEGER :: ishift,jshift,kshift
    INTEGER(kind=IPRC) :: nghb
    INTEGER :: n, l, m,ks
    INTEGER :: ifold,jfold,kfold
    INTEGER :: itemp,jtemp,ktemp
    INTEGER :: firstminpe,firstmaxpe

    INTEGER :: allneigh(0:maxneigh*numprocs-1)

    INTEGER :: sender(0:numprocs-1)

    !     itemp1      itemp2      itemp3     ibctype
    !          0           0           0           0
    !          1           0           0           1
    !          0           1           0           2
    !          1           1           0           3
    !          0           0           1           4
    !          1           0           1           5
    !          0           1           1           6
    !          1           1           1           7

    INTEGER :: nbuff=2

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
       lminx=1
    else
       lminx=minx
    endif
    if(miny.lt.1) then
       lminy=1
    else
       lminy=miny
    endif
    if(minz.lt.1) then
       lminz=1
    else
       lminz=minz
    endif
    if(maxx.gt.nx) then
       lmaxx=nx
    else
       lmaxx=maxx
    endif
    if(maxy.gt.ny) then
       lmaxy=ny
    else
       lmaxy=maxy
    endif
    if(maxz.gt.nz) then
       lmaxz=nz
    else
       lmaxz=maxz
    endif
    write(0,*)'id=',idrank,'lminx=',lminx,'lmaxx=',lmaxx,'lminy=',lminy, &
              'lmaxy=',lmaxy,'lminz=',lminz,'lmaxz=',lmaxz  
    ! sender(i): count no. of neighboring fluid nodes that talk with i-th PE
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
                IF(ownern(i4).NE.myid) THEN
                   sender(ownern(i4))=sender(ownern(i4))+1
                   !max                   write(0,*)'ownern ',ownern(i4),'for node ',i4,'ix=',i,'iy=',j,'iz=',k,'l=',l
                   !max                   write(0,*)'itemp=',itemp,'jtemp=',jtemp,'ktemp=',ktemp
                ENDIF
             enddo
          enddo
       enddo
    enddo

    j=0
    maxsender=0
    DO i=0,numprocs-1

       IF(sender(i).GT.0) THEN ! select procs that talks to this one

          IF(j.GE.maxneigh) THEN
             PRINT *,'Task',myid,'Invalid index for myreceiver',j
             STOP
          ENDIF

          IF(sender(i).GT.maxsender) THEN
             maxsender=sender(i)  ! get the max no. of nodes to talk to
          ENDIF

          i_pe2send_fluid(j)=i    ! list of fluid nodes to talk to
          n_pop2send_fluid(j)=0   ! no. of fluid nodes to talk to
          sender(i)=j             ! id of PE to talk in local ordering
          j=j+1                   ! local ordering index
       ELSE
          sender(i)=-1
       ENDIF

    ENDDO
    n_pe2send_fluid=j
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
    j=0
    DO i=0,numprocs*maxneigh-1
       IF(allneigh(i)==myid) THEN
          IF(j.GE.maxneigh) THEN
             PRINT *,'Task',myid,'Invalid index for i_pe2recv_fluid',j
             STOP
          ENDIF
          i_pe2recv_fluid(j)=i/maxneigh
          j=j+1
       ENDIF
    ENDDO
    n_pe2recv_fluid=j

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
       do k=lminz,lmaxz
          do jy=lminy,lmaxy
             do i=lminx,lmaxx
                itemp=i+ishift;
                jtemp=jy+jshift;
                ktemp=k+kshift;
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

                IF(ownern(i4).NE.myid) THEN   ! select nodes with ownership /= myid
                   j=sender(ownern(i4))
                   i_pop2send_fluid(0,n_pop2send_fluid(j),j)=l
                   i_pop2send_fluid(1,n_pop2send_fluid(j),j)=i4back(i,jy,k)
                   n_pop2send_fluid(j)=n_pop2send_fluid(j)+1
                ENDIF
             enddo
          enddo
       enddo
    enddo

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
    maxreceiver=0
    DO i=0,n_pe2recv_fluid-1
       IF(n_pop2recv_fluid(i).GT.maxreceiver) THEN
          maxreceiver=n_pop2recv_fluid(i)
       ENDIF
    ENDDO
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
#if 0
    ! this part is not necessary because the domain is regular and all nodes
    ! are fluid
    ! convert the i4-like send list (in i_pop2send_fluid) into ifl-like form

    DO i=0, n_pe2send_fluid-1
       DO j=0,n_pop2send_fluid(i)-1

          i_pop2send_fluid(1,j,i)=searchind(itppp,firstfree-1,i_pop2send_fluid(1,j,i))

          IF(i_pop2send_fluid(1,j,i).LT.0) THEN

             i_pop2send_fluid(1,j,i)=searchind(itppp,firstfree-1,nv+1)

             IF(i_pop2send_fluid(1,j,i).LT.0) THEN
                PRINT*,"(",myid,") could not find dummy node firstfree is",firstfree
             ENDIF
          ENDIF
       ENDDO
    ENDDO

    ! match the received nodes (in i4 form) with the pool of my nodes
    ! convert the i4-like receive list (in i_pop2recv_fluid) into ifl-like form

    DO i=0, n_pe2recv_fluid-1
       DO j=0,n_pop2recv_fluid(i)-1

          i_pop2recv_fluid(1,j,i)=searchind(itppp,firstfree-1,i_pop2recv_fluid(1,j,i))

          IF(i_pop2recv_fluid(1,j,i).LT.0) THEN

             i_pop2recv_fluid(1,j,i)=searchind(itppp,firstfree-1,nv+1)

             IF(i_pop2recv_fluid(1,j,i).LT.0) THEN
                PRINT*,"(",myid,") could not find dummy node firstfree is",firstfree
             ENDIF
          ENDIF
       ENDDO
    ENDDO
#else
    DO i=0, n_pe2recv_fluid-1
       DO j=0,n_pop2recv_fluid(i)-1
          ishift=ex(i_pop2recv_fluid(0,j,i))
          jshift=ey(i_pop2recv_fluid(0,j,i))
          kshift=ez(i_pop2recv_fluid(0,j,i))
          i4=i_pop2recv_fluid(1,j,i)
          ktemp = (i4/INT(nxy2,KIND=IPRC))
          jtemp = ((i4 - ktemp*INT(nxy2,KIND=IPRC))/INT(nx2,KIND=IPRC))
          itemp = (i4 - ktemp*INT(nxy2,KIND=IPRC) - jtemp*INT(nx2,KIND=IPRC))
          ktemp = ktemp+kshift
          jtemp = jtemp+jshift
          itemp = itemp+ishift
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
          endif
       ENDDO
    ENDDO
#endif

  END SUBROUTINE findneigh
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
  INTEGER :: i,j
  INTEGER(kind=IPRC) :: i4, itemp, jtemp, ktemp

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
         ktemp = i4/INT(nxy2,KIND=IPRC)
         jtemp = (i4 - ktemp*INT(nxy2,KIND=IPRC))/INT(nx2,KIND=IPRC)
         itemp = i4 - ktemp*INT(nxy2,KIND=IPRC) - jtemp*INT(nx2,KIND=IPRC)
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
  INTEGER :: i,j
  INTEGER(kind=IPRC) :: i4, itemp, jtemp, ktemp

   if(n_pe2recv_fluid.gt.0) then
      CALL MPI_WAITALL(n_pe2recv_fluid,requestm,status,ierr)
   endif
   do i=0, n_pe2recv_fluid-1
      do j=0,n_pop2recv_fluid(i)-1
         i4=i_pop2recv_fluid(1,j,i)
         ktemp = i4/INT(nxy2,KIND=IPRC)
         jtemp = (i4 - ktemp*INT(nxy2,KIND=IPRC))/INT(nx2,KIND=IPRC)
         itemp = i4 - ktemp*INT(nxy2,KIND=IPRC) - jtemp*INT(nx2,KIND=IPRC)
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
#else

  SUBROUTINE commspop(pop_ptr)

  IMPLICIT NONE

  ! REAL(KIND=PRC) :: pop_ptr(0:,0:)
  ! REAL(KIND=PRC) :: pop_ptr(0:,0:)
  !max REAL(KIND=PRC) :: pop_ptr(0:,0:)
  TYPE(REALPTR), dimension(0:links):: pop_ptr
  
  return
  
  end subroutine   
  
  
  SUBROUTINE commrpop(pop_ptr)

  IMPLICIT NONE

  TYPE(REALPTR), dimension(0:links):: pop_ptr
  
  return
  
  end subroutine

#endif


   END MODULE lbempi_mod
