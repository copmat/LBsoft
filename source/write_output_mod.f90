
#include <default_macro.h>
 module write_output_mod
 
!***********************************************************************
!     
!     LBsoft module for output VTK routines
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
 use version_mod,    only : idrank,mxrank,get_sync_world
 use lbempi_mod,     only : gminx,gmaxx,gminy,gmaxy,gminz,gmaxz,ownern,&
  i4back
 use error_mod
 use utility_mod,    only : write_fmtnumb,ltest_mode
 use fluids_mod,     only : nx,ny,nz,rhoR,rhoB,u,v,w,lsingle_fluid, &
  minx, maxx, miny, maxy, minz, maxz
  
  private
  
  character(len=*), parameter :: filenamevtk='output'
  integer, save, public, protected :: ivtkevery=50
  logical, save, public, protected :: lvtkfile=.false.
  logical, save, public, protected :: lvtkownern=.true.
  
  public :: write_vtk_frame
  public :: set_value_ivtkevery
  public :: write_test_map
  
 contains
 
 subroutine set_value_ivtkevery(ltemp,itemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for settimg time interval value
!     used to print the vtk file
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification july 2018
!     
!***********************************************************************
 
  implicit none
  
  logical, intent(in) :: ltemp
  integer, intent(in) :: itemp
  
  lvtkfile=ltemp
  ivtkevery=itemp
  
  return
  
 end subroutine set_value_ivtkevery
 
 subroutine write_vtk_frame(nstepsub)
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the serial/parallel writing
!     in structured VTK legacy binary format
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: nstepsub
  
  if(mxrank==1)then
    call write_vtk_frame_serial(nstepsub)
  else
    call write_vtk_frame_parallel(nstepsub)
    !if(mxrank/=1)call error(11)
  endif
  
  return
  
 end subroutine write_vtk_frame
  
 subroutine write_vtk_frame_serial(nstepsub)
 
!***********************************************************************
!     
!     LBsoft subroutine for writing the hydrodynamic variables
!     in structured VTK legacy binary format in serial IO
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification july 2018
!     
!***********************************************************************
     
  use IR_Precision
  use Lib_VTK_IO, only : VTK_INI_XML,VTK_GEO_XML,VTK_DAT_XML,VTK_VAR_XML,VTK_END_XML,&
   VTK_FLD_XML,PVTK_INI_XML,PVTK_GEO_XML,PVTK_DAT_XML,PVTK_VAR_XML,PVTK_END_XML
 
 
  implicit none
  
  integer, intent(in) :: nstepsub
  
  real(kind=4), allocatable, dimension(:,:,:), save ::x,y,z
  integer :: e_io,nnt
  real(kind=4), allocatable, dimension(:,:,:), save :: service1, &
   service2,service3
  
  character(len=120) :: sevt
  integer :: l,ii,jj,kk,i,j,k,myindex
 
 
  integer :: nx1,nx2,ny1,ny2,nz1,nz2,nn,nxmin,nymin,nnn
  
  
  logical, save :: lfirst=.true.
  
  if(.not. lvtkfile)return
  
  if(mxrank/=1)call error(11)
  
  if(idrank/=0)return
  
  if(mod(nstepsub,ivtkevery)/=0)return
  
  nx1=1
  nx2=nx
  ny1=1
  ny2=ny
  nz1=1
  nz2=nz
  
  nnn=nx*ny*nz
  nn=(nx2-nx1+1)*(ny2-ny1+1)*(nz2-nz1+1)
  
  if(lfirst)then
    lfirst=.false.
    allocate(x(nx1:nx2,ny1:ny2,nz1:nz2))
    allocate(y(nx1:nx2,ny1:ny2,nz1:nz2))
    allocate(z(nx1:nx2,ny1:ny2,nz1:nz2))
    allocate(service1(nx1:nx2,ny1:ny2,nz1:nz2))
    allocate(service2(nx1:nx2,ny1:ny2,nz1:nz2))
    allocate(service3(nx1:nx2,ny1:ny2,nz1:nz2))
    
    forall(i=nx1:nx2,j=ny1:ny2,k=nz1:nz2)
      x(i,j,k) = real(i,kind=R4P)
      y(i,j,k) = real(j,kind=R4P)
      z(i,j,k) = real(k,kind=R4P)
    end forall
  endif
  
  sevt=repeat(' ',120)
  sevt=trim(filenamevtk)// &
   trim(write_fmtnumb(nstepsub))//'.vts'
 
  E_IO = VTK_INI_XML(output_format='binary', &
   filename=trim(sevt),mesh_topology='StructuredGrid', &
   nx1=nx1,nx2=nx2,ny1=ny1, &
   ny2=ny2,nz1=nz1,nz2=nz2)
   
  E_IO = VTK_GEO_XML(nx1=nx1,nx2=nx2,ny1=ny1, &
   ny2=ny2,nz1=nz1,nz2=nz2,NN=nn, &
   X=reshape(x,(/nn/) ),Y=reshape(y,(/nn/)),Z=reshape(z,(/nn/)))
   
  E_IO = VTK_DAT_XML(var_location = 'node', &
   var_block_action = 'open')
 
 
 
 
  service1(nx1:nx2,ny1:ny2,nz1:nz2)=rhoR(1:nx,1:ny,1:nz)
  E_IO = VTK_VAR_XML(NC_NN=nn,varname ='density1', &
   var=reshape(service1,(/nn/)))
  
  if(.not. lsingle_fluid)then
    service2(nx1:nx2,ny1:ny2,nz1:nz2)=rhoB(1:nx,1:ny,1:nz)
    E_IO = VTK_VAR_XML(NC_NN=nn,varname ='density2', &
     var=reshape(service2,(/nn/)))
  endif
  
 
  service1(nx1:nx2,ny1:ny2,nz1:nz2)=u(1:nx,1:ny,1:nz)
  service2(nx1:nx2,ny1:ny2,nz1:nz2)=v(1:nx,1:ny,1:nz)
  service3(nx1:nx2,ny1:ny2,nz1:nz2)=w(1:nx,1:ny,1:nz)
  E_IO = VTK_VAR_XML(NC_NN=nn,varname='velocity', &
   varX=reshape(service1,(/nn/)), &
   varY=reshape(service2,(/nn/)), &
   varZ=reshape(service3,(/nn/)))
   
  E_IO = VTK_DAT_XML(var_location = 'node', &
   var_block_action = 'close')
  
  E_IO = VTK_GEO_XML()
  E_IO = VTK_END_XML()
 
 
  return
     
 end subroutine write_vtk_frame_serial
 
 subroutine initialize_vtk
 
!***********************************************************************
!     
!     LBsoft subroutine for creating the folders containing the files
!     in structured VTK legacy binary format in parallel IO
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification july 2018
!     
!***********************************************************************
  
  implicit none
     
     integer :: i,l
     character(len=255) :: path,makedirectory
     logical :: lexist
     character :: delimiter
     
     
     
     path=repeat(' ',255)
     call getcwd(path)
     !call get_environment_variable('DELIMITER',delimiter)
     path=trim(path)
     delimiter=path(1:1)
     if(delimiter==' ')delimiter='/'
     
     
     makedirectory=repeat(' ',255)
     makedirectory = 'rank'//trim(write_fmtnumb(idrank))//delimiter
#ifdef INTEL
     inquire(directory=trim(makedirectory),exist=lexist)
#else
     inquire(file=trim(makedirectory),exist=lexist)
#endif
     if(.not. lexist)then
       makedirectory=repeat(' ',255)
       makedirectory = 'mkdir rank'//trim(write_fmtnumb(idrank))
       call system(makedirectory)
     endif
     call get_sync_world()
     
     return
 end subroutine
 
 subroutine write_vtk_frame_parallel(nstepsub)
 
!***********************************************************************
!     
!     LBsoft subroutine for writing the hydrodynamic variables
!     in structured VTK legacy binary format in parallel IO
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification july 2018
!     
!***********************************************************************
     
  use IR_Precision
  use Lib_VTK_IO, only : VTK_INI_XML,VTK_GEO_XML,VTK_DAT_XML,VTK_VAR_XML,VTK_END_XML,&
   VTK_FLD_XML,PVTK_INI_XML,PVTK_GEO_XML,PVTK_DAT_XML,PVTK_VAR_XML,PVTK_END_XML
 
 
  implicit none
  
  integer, intent(in) :: nstepsub
  
  real(kind=4), allocatable, dimension(:,:,:), save ::x,y,z
  integer :: e_io,nnt
  real(kind=4), allocatable, dimension(:,:,:), save :: service1, &
   service2,service3
  
  character(len=120) :: sevt
  integer :: l,ii,jj,kk,i,j,k,myindex
 
 
  integer :: nx1,nx2,ny1,ny2,nz1,nz2,nn,nxmin,nymin,nnn
  integer :: mf(0:mxrank-1)
  
  integer(4),dimension(0:mxrank-1) :: nx1_p,nx2_p, &
   ny1_p,ny2_p,nz1_p,nz2_p,nn_p
  
  logical, save :: lfirst=.true.
  
  logical, save :: lfirstdel=.true.
  character*1, save :: delimiter
  INTEGER(kind=IPRC) :: i4
     
  if(lfirstdel)then
    call getcwd(sevt)
    sevt=trim(sevt)
    delimiter=sevt(1:1)
    if(delimiter==' ')delimiter='/'
    lfirstdel=.false.
    call initialize_vtk
   endif
  
  if(.not. lvtkfile)return
  
  if(mxrank==1)call error(11)
  
  if(mod(nstepsub,ivtkevery)/=0)return
  
  nx1=1
  nx2=nx
  ny1=1
  ny2=ny
  nz1=1
  nz2=nz
  
  do i=0,mxrank-1
    mf(i)=i
  enddo
  
  do i=0,mxrank-1
    nx1_p(i)=gminx(i)
    nx2_p(i)=gmaxx(i)
    ny1_p(i)=gminy(i)
    ny2_p(i)=gmaxy(i)
    nz1_p(i)=gminz(i)
    nz2_p(i)=gmaxz(i)
    if(nx2_p(i)<nx)then
      nx2_p(i)=nx2_p(i)+1
    endif
    if(ny2_p(i)<ny)then
      ny2_p(i)=ny2_p(i)+1
    endif
    if(nz2_p(i)<nz)then
      nz2_p(i)=nz2_p(i)+1
    endif
    nn_p(i)=(nx2_p(i)-nx1_p(i)+1)*(ny2_p(i)-ny1_p(i)+1)*(nz2_p(i)-nz1_p(i)+1)
  enddo
  
  if(idrank==0)then
    sevt=repeat(' ',120)
    sevt=trim(filenamevtk)//trim(write_fmtnumb(nstepsub))//'.pvts'
    ! pvts
    E_IO = PVTK_INI_XML(filename = trim(sevt), &
     mesh_topology = 'PStructuredGrid', &
     nx1=nx1,nx2=nx2,ny1=ny1,ny2=ny2,nz1=nz1,nz2=nz2, tp='Float32')
    do i=0,mxrank-1
      sevt=repeat(' ',120)
      sevt='rank'//trim(write_fmtnumb(i))//delimiter//trim(filenamevtk)// &
       trim(write_fmtnumb(nstepsub))//'_part'//trim(write_fmtnumb(i))//'.vts'
      E_IO = PVTK_GEO_XML(nx1=nx1_p(i),nx2=nx2_p(i),ny1=ny1_p(i), &
        ny2=ny2_p(i),nz1=nz1_p(i),nz2=nz2_p(i),source=trim(sevt))
    enddo
    E_IO = PVTK_DAT_XML(var_location = 'node', &
     var_block_action = 'open')
    E_IO = PVTK_VAR_XML(varname='density1',tp='Float32')
    if(.not. lsingle_fluid)then
      E_IO = PVTK_VAR_XML(varname='density2',tp='Float32')
    endif
    if(lvtkownern)then
      E_IO = PVTK_VAR_XML(varname='ownern',tp='Float32')
    endif
    E_IO = PVTK_VAR_XML(Nc=3_I4P,varname='velocity',tp='Float32')
    E_IO = PVTK_DAT_XML(var_location = 'node', &
     var_block_action = 'close')
    E_IO = PVTK_END_XML()
  endif
  
  nnn=nx*ny*nz
  nn=(nx2-nx1+1)*(ny2-ny1+1)*(nz2-nz1+1)
  
  if(lfirst)then
    lfirst=.false.
    allocate(x(nx1_p(idrank):nx2_p(idrank),ny1_p(idrank):ny2_p(idrank), &
     nz1_p(idrank):nz2_p(idrank)))
    allocate(y(nx1_p(idrank):nx2_p(idrank),ny1_p(idrank):ny2_p(idrank), &
     nz1_p(idrank):nz2_p(idrank)))
    allocate(z(nx1_p(idrank):nx2_p(idrank),ny1_p(idrank):ny2_p(idrank), &
     nz1_p(idrank):nz2_p(idrank)))
    
    forall(i=nx1_p(idrank):nx2_p(idrank),j=ny1_p(idrank):ny2_p(idrank), &
     k=nz1_p(idrank):nz2_p(idrank))
      x(i,j,k) = real(i,kind=R4P)
      y(i,j,k) = real(j,kind=R4P)
      z(i,j,k) = real(k,kind=R4P)
    end forall
       
    allocate(service1(nx1_p(idrank):nx2_p(idrank),ny1_p(idrank):ny2_p(idrank), &
     nz1_p(idrank):nz2_p(idrank)))
    allocate(service2(nx1_p(idrank):nx2_p(idrank),ny1_p(idrank):ny2_p(idrank), &
     nz1_p(idrank):nz2_p(idrank)))
    allocate(service3(nx1_p(idrank):nx2_p(idrank),ny1_p(idrank):ny2_p(idrank), &
     nz1_p(idrank):nz2_p(idrank)))
  endif
  
  sevt=repeat(' ',120)
  sevt='rank'//trim(write_fmtnumb(idrank))//delimiter//trim(filenamevtk)// &
   trim(write_fmtnumb(nstepsub))//'_part'//trim(write_fmtnumb(idrank))//'.vts'
     
  E_IO = VTK_INI_XML(cf=mf(idrank),output_format='binary', &
   filename=trim(sevt),mesh_topology='StructuredGrid', &
   nx1=nx1_p(idrank),nx2=nx2_p(idrank),ny1=ny1_p(idrank), &
   ny2=ny2_p(idrank),nz1=nz1_p(idrank),nz2=nz2_p(idrank))
      
  E_IO = VTK_GEO_XML(cf=mf(idrank),nx1=nx1_p(idrank),nx2=nx2_p(idrank),ny1=ny1_p(idrank), &
   ny2=ny2_p(idrank),nz1=nz1_p(idrank),nz2=nz2_p(idrank),NN=nn_p(idrank), &
   X=reshape(x,(/nn_p(idrank)/) ),Y=reshape(y,(/nn_p(idrank)/)),Z=reshape(z,(/nn_p(idrank)/)))
       
  E_IO = VTK_DAT_XML(cf=mf(idrank),var_location = 'node', &
   var_block_action = 'open')
     
  service1(nx1_p(idrank):nx2_p(idrank),ny1_p(idrank):ny2_p(idrank), &
   nz1_p(idrank):nz2_p(idrank))= &
   rhoR(nx1_p(idrank):nx2_p(idrank),ny1_p(idrank):ny2_p(idrank), &
   nz1_p(idrank):nz2_p(idrank))
  E_IO = VTK_VAR_XML(cf=mf(idrank),NC_NN=nn_p(idrank),varname ='density1', &
   var=reshape(service1,(/nn_p(idrank)/)))
     
  if(.not. lsingle_fluid)then
    service1(nx1_p(idrank):nx2_p(idrank),ny1_p(idrank):ny2_p(idrank), &
     nz1_p(idrank):nz2_p(idrank))= &
     rhoB(nx1_p(idrank):nx2_p(idrank),ny1_p(idrank):ny2_p(idrank), &
     nz1_p(idrank):nz2_p(idrank))
    E_IO = VTK_VAR_XML(cf=mf(idrank),NC_NN=nn_p(idrank),varname ='density2', &
     var=reshape(service1,(/nn_p(idrank)/)))
  endif
  if(lvtkownern)then
    do k=nz1_p(idrank),nz2_p(idrank)
      do j=ny1_p(idrank),ny2_p(idrank)
        do i=nx1_p(idrank),nx2_p(idrank)
          i4=i4back(i,j,k)
          service1(i,j,k)= ownern(i4)
        enddo
      enddo
    enddo
    E_IO = VTK_VAR_XML(cf=mf(idrank),NC_NN=nn_p(idrank),varname ='ownern', &
     var=reshape(service1,(/nn_p(idrank)/)))
  endif
     
  service1(nx1_p(idrank):nx2_p(idrank),ny1_p(idrank):ny2_p(idrank), &
   nz1_p(idrank):nz2_p(idrank))= &
   u(nx1_p(idrank):nx2_p(idrank),ny1_p(idrank):ny2_p(idrank), &
   nz1_p(idrank):nz2_p(idrank))
  service2(nx1_p(idrank):nx2_p(idrank),ny1_p(idrank):ny2_p(idrank), &
   nz1_p(idrank):nz2_p(idrank))= &
   v(nx1_p(idrank):nx2_p(idrank),ny1_p(idrank):ny2_p(idrank), &
   nz1_p(idrank):nz2_p(idrank))
  service3(nx1_p(idrank):nx2_p(idrank),ny1_p(idrank):ny2_p(idrank), &
   nz1_p(idrank):nz2_p(idrank))= &
   w(nx1_p(idrank):nx2_p(idrank),ny1_p(idrank):ny2_p(idrank), &
   nz1_p(idrank):nz2_p(idrank))
     
  E_IO = VTK_VAR_XML(cf=mf(idrank),NC_NN=nn_p(idrank),varname='velocity', &
   varX=reshape(service1,(/nn_p(idrank)/)), &
   varY=reshape(service2,(/nn_p(idrank)/)), &
   varZ=reshape(service3,(/nn_p(idrank)/)))
      
     
  E_IO = VTK_DAT_XML(cf=mf(idrank),var_location = 'node', &
   var_block_action = 'close')
      
  E_IO = VTK_GEO_XML(cf=mf(idrank))
  E_IO = VTK_END_XML(cf=mf(idrank))
  
  return
     
 end subroutine write_vtk_frame_parallel
 
 subroutine write_test_map()
  
!***********************************************************************
!     
!     LBsoft subroutine for writing the hydrodynamic variables
!     in ASCII format for diagnostic purposes
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification july 2018
!     
!***********************************************************************

  implicit none
  
  integer :: i,j,k,l
  
  integer, parameter :: iomap=250
  
  character(len=256) :: file_loc_proc

  if(.not. ltest_mode)return
  
  if(idrank==0)then 
    open(iomap-1,file='global.map',status='replace',action='write')
    write(iomap-1,'(4i10)')mxrank,nx,ny,nz
    close(iomap-1)
  endif

  l=idrank
  file_loc_proc = 'output'//trim(write_fmtnumb(l))//'.map'
  open(iomap+l,file=trim(file_loc_proc),status='replace',action='write')
  write(iomap+l,'(i10)')mxrank
  if(lsingle_fluid)then
    write(iomap+l,'(8i10)')minx,maxx,miny,maxy,minz,maxz, PRC, 4
    do k=minz,maxz
      do j=miny,maxy
        do i=minx,maxx
          write(iomap+l,'(3i8,4f24.16)')i,j,k,rhoR(i,j,k), &
           u(i,j,k),v(i,j,k),w(i,j,k)
        enddo
      enddo
    enddo
  else
    write(iomap+l,'(8i10)')minx,maxx,miny,maxy,minz,maxz, PRC, 5
    do k=minz,maxz
      do j=miny,maxy
        do i=minx,maxx
          write(iomap+l,'(3i8,5f24.16)')i,j,k,rhoR(i,j,k),rhoB(i,j,k), &
           u(i,j,k),v(i,j,k),w(i,j,k)
        enddo
      enddo
    enddo
  endif
  
  close(iomap+l)
  
 end subroutine write_test_map
 
 end module write_output_mod
