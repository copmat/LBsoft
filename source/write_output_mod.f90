
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
  
 use version_mod,    only : idrank,mxrank
 use error_mod
 use utility_mod,    only : write_fmtnumb,ltest_mode
 use fluids_mod,     only : nx,ny,nz,rhoR,rhoB,u,v,w,lsingle_fluid, &
  minx, maxx, miny, maxy, minz, maxz
  
  private
  
  character(len=*), parameter :: filenamevtk='output'
  integer, save, public, protected :: ivtkevery=50
  logical, save, public, protected :: lvtkfile=.false.
  
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
!     LBsoft subroutine for writing the hydrodynamic variables
!     in structured VTK legacy binary format in serial
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
     
 end subroutine write_vtk_frame
 
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
    write(iomap-1,'(i10)')mxrank,nx,ny,nz
    close(iomap-1)
  endif

  l=idrank
  file_loc_proc = 'output'//trim(write_fmtnumb(l))//'.map'
  open(iomap+l,file=trim(file_loc_proc),status='replace',action='write')
  write(iomap+l,'(i10)')mxrank
  if(lsingle_fluid)then
    write(iomap+l,'(6i10)')nx,ny,nz, PRC, 4
    do k=minz,maxz
      do j=miny,maxy
        do i=minx,maxx
          write(iomap+l,'(3i8,4f20.10)')i,j,k,rhoR(i,j,k), &
           u(i,j,k),v(i,j,k),w(i,j,k)
        enddo
      enddo
    enddo
  else
    write(iomap+l,'(6i10)')nx,ny,nz, PRC, 5
    do k=minz,maxz
      do j=miny,maxy
        do i=minx,maxx
          write(iomap+l,'(3i8,5f20.10)')i,j,k,rhoR(i,j,k),rhoB(i,j,k), &
           u(i,j,k),v(i,j,k),w(i,j,k)
        enddo
      enddo
    enddo
  endif
  
  close(iomap+l)
  
 end subroutine write_test_map
 
 end module write_output_mod
