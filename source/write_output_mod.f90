
#include <default_macro.h>
 module write_output_mod
 
!***********************************************************************
!     
!     LBsoft module for output VTK routines
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
 use version_mod,    only : idrank,mxrank,get_sync_world,finalize_world, &
  bcast_world_i,open_file_vtk_par,close_file_vtk_par, &
  print_header_vtk_par,print_footer_vtk_par,gminx,gmaxx,gminy,gmaxy, &
  gminz,gmaxz,ownern,i4back,driving_print_binary_1d_vtk, &
  driving_print_binary_3d_vtk
 use error_mod
 use utility_mod,    only : write_fmtnumb,ltest_mode, &
  test_little_endian,space_fmtnumb12,space_fmtnumb
  


 use fluids_mod,     only : nx,ny,nz,rhoR,rhoB,u,v,w,lsingle_fluid, &
  minx, maxx, miny, maxy, minz, maxz,isfluid, dumpHvar,bin_oneFile,dump_oneFile

 use particles_mod,  only : natms,xxx,yyy,zzz,lparticles,cell, &
  ishape,lrotate,natms,natms_tot,q0,q1,q2,q3,vxx,vyy,vzz, &
  take_rotversorx,take_rotversory,take_rotversorz, &
  clean_fluid_inside_particle,q2eul, atmbook, natms_ext, dumppart_oneFile
 
  private
  
  logical, save :: lframe_isfluid=.true.
  integer, parameter :: iotest=180
  integer, parameter :: ioxyz=190
  integer, parameter :: mxln=120
  character(len=*), parameter :: filenamevtk='out'
  integer, save, public, protected :: ivtkevery=100
  logical, save, public, protected :: lvtkfile=.false.
  logical, save, public, protected :: lvtkownern=.false.
  integer, save, public, protected :: ixyzevery=100
  logical, save, public, protected :: lxyzfile=.false.
  integer, save, public, protected :: ibinevery=100
  integer, save, public, protected :: idumpevery=10000
  logical, save, public, protected :: lbinevery=.false.
  logical, save, public, protected :: lelittle=.true.
  character(len=mxln), save :: dir_out
  character(len=mxln), save :: dir_out_rank
  character(len=mxln), allocatable, save :: dir_out_grank(:)
  character(len=mxln), save :: extentvtk
  character(len=mxln), save :: extentvtk_isfluid
  character(len=8), save, allocatable, dimension(:) :: namevarvtk
  character(len=500), save, allocatable, dimension(:) :: headervtk
  character(len=30), save, allocatable, dimension(:) :: footervtk
  integer, save, allocatable, dimension(:) :: ndimvtk
  integer, save, allocatable, dimension(:) :: vtkoffset
  integer, save, allocatable, dimension(:) :: ndatavtk
  integer, save, allocatable, dimension(:) :: nheadervtk
  integer, save :: nfilevtk
  integer, save, allocatable, dimension(:) :: varlistvtk
  
  character, save :: delimiter
  
  public :: write_vtk_frame
  public :: init_output
  public :: set_value_ivtkevery
  public :: write_test_map
  public :: write_xyz
  public :: write_xyz_close
  public :: set_value_ixyzevery
  public :: write_particle_xyz
  public :: set_value_ibinevery
  public :: dumpForOutput
  public :: set_value_dumpevery
  public :: read_restart_file
  public :: write_vtk_isfluid
  
 contains
 
 subroutine init_output(nprintlistvtk,printlistvtk)
 
!***********************************************************************
!     
!     LBsoft subroutine for creating the folders containing the files
!     in image VTK legacy binary format in parallel IO
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: nprintlistvtk
  integer, intent(in), allocatable, dimension(:) :: printlistvtk
  character(len=255) :: path,makedirectory
  logical :: lexist
  
  integer :: i,nn,indent,myoffset,new_myoffset,iend,nn_isfluid
  integer, parameter :: byter4=4
  integer, parameter :: byteint=4
  
  if((.not. lvtkfile).and.(.not. lxyzfile) .and. (.not. lbinevery))return
  
  call test_little_endian(lelittle)
  
  path = repeat(' ',255)
  call getcwd(path)
  
  !call get_environment_variable('DELIMITER',delimiter)
  path = trim(path)
  delimiter = path(1:1)
  if (delimiter==' ') delimiter='/'


  if(lbinevery) then
     makedirectory=repeat(' ',255)
     makedirectory = 'dumpBin'
     dir_out=trim(makedirectory)
#ifdef INTEL
     inquire(directory=trim(makedirectory),exist=lexist)
#else
     inquire(file=trim(makedirectory),exist=lexist)
#endif
   
     if(.not. lexist)then
       if(idrank==0) then
         makedirectory=repeat(' ',255)
         makedirectory = 'mkdir dumpBin'
         call system(makedirectory)
       endif
     endif
  endif
  
  makedirectory=repeat(' ',255)
  makedirectory = 'output'//delimiter
  dir_out=trim(makedirectory)
#ifdef INTEL
  inquire(directory=trim(makedirectory),exist=lexist)
#else
  inquire(file=trim(makedirectory),exist=lexist)
#endif
  
  if(.not. lexist)then
    if(idrank==0) then
      makedirectory=repeat(' ',255)
      makedirectory = 'mkdir output'
      call system(makedirectory)
    endif
  endif

  call get_sync_world()
  
  call write_xyz_open(trim(dir_out)//'output.xyz')
  
  if(.not. lvtkfile)return
  
  if(lframe_isfluid)then
    if(mxrank>1)then
      call warning(72)
      lframe_isfluid=.false.
    endif
  endif
  
  makedirectory=repeat(' ',255)
  makedirectory=trim(path)//delimiter//'output'//delimiter
  
  extentvtk =  space_fmtnumb(1) // ' ' // space_fmtnumb(nx) // ' ' &
        // space_fmtnumb(1) // ' ' // space_fmtnumb(ny) // ' ' &
        // space_fmtnumb(1) // ' ' // space_fmtnumb(nz)
  
  if(lframe_isfluid)then
    extentvtk_isfluid =  space_fmtnumb(0) // ' ' // space_fmtnumb(nx+1) // ' ' &
        // space_fmtnumb(0) // ' ' // space_fmtnumb(ny+1) // ' ' &
        // space_fmtnumb(0) // ' ' // space_fmtnumb(nz+1)
  else
    extentvtk_isfluid =  space_fmtnumb(1) // ' ' // space_fmtnumb(nx) // ' ' &
        // space_fmtnumb(1) // ' ' // space_fmtnumb(ny) // ' ' &
        // space_fmtnumb(1) // ' ' // space_fmtnumb(nz)
  endif
  nfilevtk=nprintlistvtk
  
  allocate(varlistvtk(0:nfilevtk))
  allocate(namevarvtk(0:nfilevtk))
  allocate(ndimvtk(0:nfilevtk))
  allocate(headervtk(0:nfilevtk))
  allocate(footervtk(0:nfilevtk))
  allocate(nheadervtk(0:nfilevtk))
  allocate(vtkoffset(0:nfilevtk))
  allocate(ndatavtk(0:nfilevtk))
  varlistvtk(0)=-1
  varlistvtk(1:nfilevtk)=printlistvtk(1:nfilevtk)
  do i=0,nfilevtk
    select case(varlistvtk(i))
    case(-1)
      namevarvtk(i)='isfluid '
      ndimvtk(i)=1
    case(1)
      namevarvtk(i)='rho1    '
      ndimvtk(i)=1
    case(2)
      namevarvtk(i)='rho2    '
      ndimvtk(i)=1
    case(3)
      namevarvtk(i)='phase   '
      ndimvtk(i)=1
    case(4)
      namevarvtk(i)='vel     '
      ndimvtk(i)=3
    case(5)
      namevarvtk(i)='part    '
      ndimvtk(i)=1
    case default
      call error(38)
    end select
  enddo
  
  nn=nx*ny*nz
  if(lframe_isfluid)then
    nn_isfluid=(nx+2)*(ny+2)*(nz+2)
  else
    nn_isfluid=nn
  endif
  
  do i=0,nfilevtk
    myoffset=0
    indent=0
    if(i==0)then
      call header_vtk(headervtk(i),namevarvtk(i),extentvtk_isfluid, &
       ndimvtk(i),0,iend,myoffset, new_myoffset,indent)
      vtkoffset(i)=new_myoffset
      myoffset=new_myoffset+byteint+ndimvtk(i)*nn_isfluid*byter4
      ndatavtk(i)=ndimvtk(i)*nn_isfluid*byter4
      nheadervtk(i)=iend
      call footer_vtk(footervtk(i),0,iend,myoffset, &
       new_myoffset,indent)
    else
      call header_vtk(headervtk(i),namevarvtk(i),extentvtk,ndimvtk(i),0, &
       iend,myoffset,new_myoffset,indent)
      vtkoffset(i)=new_myoffset
      myoffset=new_myoffset+byteint+ndimvtk(i)*nn*byter4
      ndatavtk(i)=ndimvtk(i)*nn*byter4
     nheadervtk(i)=iend
      call footer_vtk(footervtk(i),0,iend,myoffset, &
       new_myoffset,indent)
    endif
  enddo

#if 0
  call chdir(makedirectory)
  
  if(idrank==0)then
    allocate(dir_out_grank(0:mxrank-1))
    do i=0,mxrank-1
      makedirectory=repeat(' ',255)
      makedirectory = 'rank'//trim(write_fmtnumb(i))//delimiter
      dir_out_grank(i)=trim(dir_out)//trim(makedirectory)
    enddo
  endif
  
  makedirectory=repeat(' ',255)
  makedirectory = 'rank'//trim(write_fmtnumb(idrank))//delimiter
  dir_out_rank=trim(dir_out)//trim(makedirectory)
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
  
  call chdir(path)
  
#endif
  
  call get_sync_world
  
  return

 end subroutine init_output
 
 subroutine write_vtp_file(iunitsub,fnamesub,istepsub)
 
!***********************************************************************
!     
!     LBsoft subroutine for writing the POLYDATA VTK file
!     for the particles
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: iunitsub
  character(len=*), intent(in) :: fnamesub
  integer, intent(in) :: istepsub
  
  integer :: iatm,imio,jatm
  
  
  character(len=120) :: sevt
  integer, allocatable, dimension(:,:) :: list
  
  real(kind=PRC) :: mydist,miomax,vmod
  real(kind=PRC), dimension(3) :: uxx,uyy,uzz
  
  integer :: i,j,k,ii,jj,kk
  
  103 format (A)
  104 format (ES16.5)
  105 format (A,1X,I12,I12,I12)
  106 format (A,I12,A)
  107 format (A,I12)
  108 format (ES16.5,ES16.5,ES16.5)
  109 format (I12,I12,I12)
  110 format (I12,I12)
  112 format (A,I12,I12)
  
  if(idrank/=0)return
  
  sevt=repeat(' ',120)
  sevt=trim(dir_out)//trim(fnamesub)// &
   trim(write_fmtnumb(istepsub))//'.vtk'
   
  
  open(unit=iotest,file=trim(sevt),status='replace',action='write')
  
  write(iotest,103) '# vtk DataFile Version 2.0'
  write(iotest,103) 'Field Emission Device - Charge Density Plot'
  write(iotest,103) 'ASCII'
  write(iotest,103) 'DATASET POLYDATA'
  write(iotest,106) 'POINTS ',natms_tot,' float'
  do iatm=1,natms_tot
    write(iotest,108)xxx(iatm),yyy(iatm),zzz(iatm)
  enddo
  
  write(iotest,*)
  write(iotest,107) 'POINT_DATA ', natms_tot
  write(iotest,103) 'SCALARS Type float 1'
  write(iotest,103) 'LOOKUP_TABLE default'
  
  do iatm=1,natms_tot
    write(iotest,*)1 !ltype(iatm)
  enddo
  
  if(lrotate)then
    
    write(iotest,103) 'VECTORS Vectx float'
    do iatm=1,natms_tot
      write(iotest,*)take_rotversorx(q0(iatm),q1(iatm),q2(iatm),q3(iatm))
    enddo
    
    write(iotest,103) 'VECTORS Vecty float'
    do iatm=1,natms_tot
      write(iotest,*)take_rotversory(q0(iatm),q1(iatm),q2(iatm),q3(iatm))
    enddo
    
    write(iotest,103) 'VECTORS Vectz float'
    do iatm=1,natms_tot
      write(iotest,*)take_rotversorz(q0(iatm),q1(iatm),q2(iatm),q3(iatm))
    enddo
    
  endif
  
  close(iotest)
   
  return
  
 end subroutine write_vtp_file
 
 subroutine write_xyz_open(filename)
 
!***********************************************************************
!     
!     LBsoft subroutine for opening the file in xyz format
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
  
  implicit none
  
  character(len=*), intent(in) :: filename
  
  if(idrank/=0)return
  if(.not. lxyzfile)return
  
  open(ioxyz,file=trim(filename),status='replace',action='write')
  
  return
    
 end subroutine write_xyz_open
 
 subroutine write_xyz(istepsub)
 
!***********************************************************************
!     
!     LBsoft subroutine for writing a frame in the xyz file
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: istepsub
  
  integer :: j
  character(len=8), parameter :: mystring8='C       '
  character(len=13), parameter :: mystring13='time step    '
  
  if(idrank/=0)return
  if(.not. lxyzfile)return
  if(mod(istepsub,ixyzevery)/=0)return

  write(ioxyz,'(i8)') natms_tot
  write(ioxyz,'(a13,i12,3f16.8)')mystring13,istepsub,cell(1),cell(5),cell(9)
  do j=1,natms_tot
    write(ioxyz,"(a8,3f16.6)")mystring8,xxx(j),yyy(j),zzz(j)
  end do

  return

    
 end subroutine write_xyz

 subroutine write_xyz_close()
  
!***********************************************************************
!     
!     LBsoft subroutine for closing the file in xyz format
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
  
  implicit none
  
  if(idrank/=0)return
  if(.not. lxyzfile)return
  close(ioxyz)
  
  return
    
 end subroutine write_xyz_close

 subroutine set_value_ivtkevery(ltemp,itemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for settimg time interval value
!     used to print the vtk file
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
 
 subroutine set_value_ixyzevery(ltemp,itemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for settimg time interval value
!     used to print the xyz file
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
  
  logical, intent(in) :: ltemp
  integer, intent(in) :: itemp
  
  lxyzfile=ltemp
  ixyzevery=itemp
  
  return
  
 end subroutine set_value_ixyzevery
 
 subroutine write_vtk_frame(nstepsub,wantRestore)
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the serial/parallel writing
!     in structured VTK legacy binary format
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: nstepsub
  logical, intent(in), optional :: wantRestore
  
  integer, parameter :: vtkoutput=55
  integer :: i
  
  if((.not. lvtkfile).and.(.not. lxyzfile))return
  
  if(mod(nstepsub,ivtkevery)/=0)return
  if(present(wantRestore))then
    if(.not.wantRestore)call clean_fluid_inside_particle
  else
    call clean_fluid_inside_particle
  endif
  
  if(mxrank==1)then
    do i=1,nfilevtk
      select case(varlistvtk(i))
      case(1)
        call write_vtk_frame_serial(vtkoutput+i,nstepsub,ndatavtk(i), &
         nheadervtk(i),headervtk(i),vtkoffset(i),footervtk(i), &
         namevarvtk(i),ndimvtk(i),varlistvtk(i),rhoR)
      case(2)
        call write_vtk_frame_serial(vtkoutput+i,nstepsub,ndatavtk(i), &
         nheadervtk(i),headervtk(i),vtkoffset(i),footervtk(i), &
         namevarvtk(i),ndimvtk(i),varlistvtk(i),rhoB)
      case(3)
        call write_vtk_frame_serial(vtkoutput+i,nstepsub,ndatavtk(i), &
         nheadervtk(i),headervtk(i),vtkoffset(i),footervtk(i), &
         namevarvtk(i),ndimvtk(i),varlistvtk(i),rhoR,rhoB)
      case(4)
        call write_vtk_frame_serial(vtkoutput+i,nstepsub,ndatavtk(i), &
         nheadervtk(i),headervtk(i),vtkoffset(i),footervtk(i), &
         namevarvtk(i),ndimvtk(i),varlistvtk(i),u,v,w)
      case(5)
        if(lparticles)then
          call write_vtp_file(181,'outatm',nstepsub)
        endif
      case default
        call error(40)
      end select
    enddo
  else
    do i=1,nfilevtk
      select case(varlistvtk(i))
      case(1)
        call write_vtk_frame_parallel(vtkoutput+i,nstepsub,ndatavtk(i), &
         nheadervtk(i),headervtk(i),vtkoffset(i),footervtk(i), &
         namevarvtk(i),ndimvtk(i),varlistvtk(i),rhoR)
      case(2)
        call write_vtk_frame_parallel(vtkoutput+i,nstepsub,ndatavtk(i), &
         nheadervtk(i),headervtk(i),vtkoffset(i),footervtk(i), &
         namevarvtk(i),ndimvtk(i),varlistvtk(i),rhoB)
      case(3)
        call write_vtk_frame_parallel(vtkoutput+i,nstepsub,ndatavtk(i), &
         nheadervtk(i),headervtk(i),vtkoffset(i),footervtk(i), &
         namevarvtk(i),ndimvtk(i),varlistvtk(i),rhoR,rhoB)
      case(4)
        call write_vtk_frame_parallel(vtkoutput+i,nstepsub,ndatavtk(i), &
         nheadervtk(i),headervtk(i),vtkoffset(i),footervtk(i), &
         namevarvtk(i),ndimvtk(i),varlistvtk(i),u,v,w)
      case(5)
        if(lparticles)then
          call write_vtp_file(181,'outatm',nstepsub)
        endif
      case default
        call error(40)
      end select
    enddo
  endif
  
  return
  
 end subroutine write_vtk_frame
 
 subroutine write_vtk_isfluid(nstepsub)
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the serial/parallel writing
!     in structured VTK legacy binary format
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: nstepsub
  
  integer, parameter :: vtkoutput=55
  integer :: i
  
  if(.not. lvtkfile)return
  
  
  if(mxrank==1)then
    i=0
    call write_vtk_isfluid_serial(vtkoutput+i,nstepsub,ndatavtk(i), &
     nheadervtk(i),headervtk(i),vtkoffset(i),footervtk(i), &
     namevarvtk(i),ndimvtk(i),varlistvtk(i))
  else
    i=0
    call write_vtk_isfluid_parallel(vtkoutput+i,nstepsub,ndatavtk(i), &
     nheadervtk(i),headervtk(i),vtkoffset(i),footervtk(i), &
     namevarvtk(i),ndimvtk(i),varlistvtk(i))
  endif
  
  return
  
 end subroutine write_vtk_isfluid
 
 subroutine write_vtk_isfluid_serial(ioprint,nstepsub,ndata,nheader,header, &
  headoff,footer,filenamevar,ndim,ntype)
 
!***********************************************************************
!     
!     LBsoft subroutine for writing the isfluid variable
!     in structured VTK legacy binary format in serial IO
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification November 2019
!     
!***********************************************************************
 
 
  implicit none
  
  integer, intent(in) :: ioprint,nstepsub,ndata,nheader,headoff,ndim, &
   ntype
  character(len=8), intent(in) :: filenamevar
  character(len=500), intent(in) :: header
  character(len=30), intent(in) :: footer
  integer :: e_io,nnt
  real(kind=4), allocatable, dimension(:,:,:) :: service1
  
  character(len=120) :: sevt
  integer :: l,ii,jj,kk,i,j,k,myindex,iost,endoff
  
  integer :: nx1,nx2,ny1,ny2,nz1,nz2,nn,nxmin,nymin,nnn
  
  character(len=120) :: s_buffer
  
  character(len=500) :: mystring500
  integer, parameter :: byter4=4
  integer, parameter :: byteint=4
  real(kind=PRC) :: miomax,vmod
  real(kind=4) :: r4temp,b4temp
  
  if(.not. lvtkfile)return
  
  if(mxrank/=1)call error(11)
  
  if(idrank/=0)return
  
  if(lframe_isfluid)then
    nx1=0
    nx2=nx+1
    ny1=0
    ny2=ny+1
    nz1=0
    nz2=nz+1
  else
    nx1=1
    nx2=nx
    ny1=1
    ny2=ny
    nz1=1
    nz2=nz
  endif
  
  nn=(nx2-nx1+1)*(ny2-ny1+1)*(nz2-nz1+1)
  nnn=nn
  
  allocate(service1(nx1:nx2,ny1:ny2,nz1:nz2))
  service1(nx1:nx2,ny1:ny2,nz1:nz2)=ZERO
  
  if(ndim==1)then
  !requested scalar field in output
    if(ntype==-1)then
       !requested phase field in output
       !compute phase field on fly
       do k=nz1,nz2
        do j=ny1,ny2
          do i=nx1,nx2
            service1(i,j,k)=real(isfluid(i,j,k),kind=4)
          enddo
        enddo
      enddo
    else
      call error(42)
    endif
  elseif(ndim==3)then
    !requested vector field in output
    call error(42)
  endif
  
  sevt=repeat(' ',120)
  
  sevt = trim(dir_out) // trim(filenamevtk)//'_'//trim(filenamevar)// &
   '_'//trim(write_fmtnumb(nstepsub)) // '.vti'
  
  call open_file_vtk(ioprint,120,sevt,e_io)
  
  call print_header_vtk(ioprint,0,nheader,header,E_IO)
  
  call driving_print_binary_1d_vtk(ioprint,headoff,ndata,nn,service1, &
   iost)
  
  endoff=headoff+ndata+byteint
  call print_footer_vtk(ioprint,endoff,footer,E_IO)
  
  call close_file_vtk(ioprint,e_io)
  
  if(allocated(service1))deallocate(service1)
 
  return
     
 end subroutine write_vtk_isfluid_serial
 
 subroutine write_vtk_isfluid_parallel(ioprint,nstepsub,ndata,nheader,header, &
  headoff,footer,filenamevar,ndim,ntype)
 
!***********************************************************************
!     
!     LBsoft subroutine for writing the isfluid variable
!     in structured VTK legacy binary format in parallel IO
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification November 2019
!     
!***********************************************************************
 
 
  implicit none
  
  integer, intent(in) :: ioprint,nstepsub,ndata,nheader,headoff,ndim, &
   ntype
  character(len=8), intent(in) :: filenamevar
  character(len=500), intent(in) :: header
  character(len=30), intent(in) :: footer
  integer :: e_io,nnt
  real(kind=4), allocatable, dimension(:,:,:) :: service1
  
  character(len=120) :: sevt
  integer :: l,ii,jj,kk,i,j,k,myindex,iost,endoff
  
  integer :: nx1,nx2,ny1,ny2,nz1,nz2,nn,nxmin,nymin,nnn
  
  character(len=120) :: s_buffer
  
  character(len=500) :: mystring500
  integer, parameter :: byter4=4
  integer, parameter :: byteint=4
  real(kind=PRC) :: miomax,vmod
  real(kind=4) :: r4temp,b4temp
  
  if(.not. lvtkfile)return
  if(lframe_isfluid)call error(44)
  
  nx1=minx
  nx2=maxx
  ny1=miny
  ny2=maxy
  nz1=minz
  nz2=maxz
  
  nnn=nx*ny*nz
  nn=(nx2-nx1+1)*(ny2-ny1+1)*(nz2-nz1+1)
  
  if(ndim==1)then
    allocate(service1(nx1:nx2,ny1:ny2,nz1:nz2))
    service1(nx1:nx2,ny1:ny2,nz1:nz2)=ZERO
  !requested scalar field in output
    if(ntype==-1)then
       !requested phase field in output
       !compute phase field on fly
       do k=minz,maxz
        do j=miny,maxy
          do i=minx,maxx
            service1(i,j,k)=real(isfluid(i,j,k),kind=4)
          enddo
        enddo
      enddo
    else
      call error(43)
    endif
  elseif(ndim==3)then
    !requested vector field in output
    call error(43)
  endif
  
  sevt=repeat(' ',120)
  
  sevt = trim(dir_out) // trim(filenamevtk)//'_'//trim(filenamevar)// &
   '_'//trim(write_fmtnumb(nstepsub)) // '.vti'
  
  call open_file_vtk_par(ioprint,120,sevt,e_io)
  
  if(idrank==0)call print_header_vtk_par(ioprint,0,nheader,header,E_IO)
  
  endoff=headoff+ndata+byteint
  
  if(idrank==0)call print_footer_vtk_par(ioprint,endoff,footer,E_IO)
  
  call driving_print_binary_1d_vtk(ioprint,headoff,ndata,nn,service1, &
   iost)
  
  call close_file_vtk_par(ioprint,e_io)
  
  if(allocated(service1))deallocate(service1)
  
  return
     
 end subroutine write_vtk_isfluid_parallel
  
 subroutine write_vtk_frame_serial(ioprint,nstepsub,ndata,nheader,header, &
  headoff,footer,filenamevar,ndim,ntype,myvar,myvary,myvarz)
 
!***********************************************************************
!     
!     LBsoft subroutine for writing the hydrodynamic variables
!     in structured VTK legacy binary format in serial IO
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
 
  implicit none
  
  integer, intent(in) :: ioprint,nstepsub,ndata,nheader,headoff,ndim, &
   ntype
  character(len=8), intent(in) :: filenamevar
  character(len=500), intent(in) :: header
  character(len=30), intent(in) :: footer
  real(kind=PRC), allocatable, dimension(:,:,:) :: myvar
  real(kind=PRC), allocatable, dimension(:,:,:), optional :: myvary,myvarz
  integer :: e_io,nnt
  real(kind=4), allocatable, dimension(:,:,:) :: service1 , &
   service2,service3
  
  character(len=120) :: sevt
  integer :: l,ii,jj,kk,i,j,k,myindex,iost,endoff
  
  integer :: nx1,nx2,ny1,ny2,nz1,nz2,nn,nxmin,nymin,nnn
  
  character(len=120) :: s_buffer
  
  character(len=500) :: mystring500
  integer, parameter :: byter4=4
  integer, parameter :: byteint=4
  real(kind=PRC) :: miomax,vmod
  real(kind=4) :: r4temp,b4temp
  
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
  
  allocate(service1(nx1:nx2,ny1:ny2,nz1:nz2))
  service1(nx1:nx2,ny1:ny2,nz1:nz2)=ZERO
  
  if(ndim==1)then
  !requested scalar field in output
    if(ntype==3)then
       !requested phase field in output
       !compute phase field on fly
       do k=minz,maxz
        do j=miny,maxy
          do i=minx,maxx
            if(isfluid(i,j,k)==3 .or. isfluid(i,j,k)==0)cycle
            r4temp=real(myvar(i,j,k),kind=4)
            b4temp=real(myvary(i,j,k),kind=4)
            service1(i,j,k)=(r4temp-b4temp)/(r4temp+b4temp)
          enddo
        enddo
      enddo
    else
      do k=minz,maxz
        do j=miny,maxy
          do i=minx,maxx
            if(isfluid(i,j,k)==3 .or. isfluid(i,j,k)==0)cycle
            service1(i,j,k)=real(myvar(i,j,k),kind=4)
          enddo
        enddo
      enddo
    endif
  elseif(ndim==3)then
    !requested vector field in output
    allocate(service2(nx1:nx2,ny1:ny2,nz1:nz2))
    allocate(service3(nx1:nx2,ny1:ny2,nz1:nz2))
    service2(nx1:nx2,ny1:ny2,nz1:nz2)=ZERO
    service3(nx1:nx2,ny1:ny2,nz1:nz2)=ZERO
    do k=minz,maxz
      do j=miny,maxy
        do i=minx,maxx
          if(isfluid(i,j,k)==3 .or. isfluid(i,j,k)==0)cycle
          service1(i,j,k)=real(myvar(i,j,k),kind=4)
          service2(i,j,k)=real(myvary(i,j,k),kind=4)
          service3(i,j,k)=real(myvarz(i,j,k),kind=4)
        enddo
      enddo
    enddo
  endif
  
  sevt=repeat(' ',120)
  
  sevt = trim(dir_out) // trim(filenamevtk)//'_'//trim(filenamevar)// &
   '_'//trim(write_fmtnumb(nstepsub)) // '.vti'
  
  call open_file_vtk(ioprint,120,sevt,e_io)
  
  call print_header_vtk(ioprint,0,nheader,header,E_IO)
  
  if(ndim==1)then
    call driving_print_binary_1d_vtk(ioprint,headoff,ndata,nn,service1, &
     iost)
  elseif(ndim==3)then
    call driving_print_binary_3d_vtk(ioprint,headoff,ndata,nn,service1, &
     service2,service3,iost)
  endif
  
  endoff=headoff+ndata+byteint
  call print_footer_vtk(ioprint,endoff,footer,E_IO)
  
  call close_file_vtk(ioprint,e_io)
  
  if(allocated(service1))deallocate(service1)
  if(allocated(service2))deallocate(service2)
  if(allocated(service3))deallocate(service3)
 
  return
     
 end subroutine write_vtk_frame_serial
 
 subroutine write_vtk_frame_parallel(ioprint,nstepsub,ndata,nheader,header, &
  headoff,footer,filenamevar,ndim,ntype,myvar,myvary,myvarz)
 
!***********************************************************************
!     
!     LBsoft subroutine for writing the hydrodynamic variables
!     in structured VTK legacy binary format in parallel IO
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
 
  implicit none
  
  integer, intent(in) :: ioprint,nstepsub,ndata,nheader,headoff,ndim, &
   ntype
  character(len=8), intent(in) :: filenamevar
  character(len=500), intent(in) :: header
  character(len=30), intent(in) :: footer
  real(kind=PRC), allocatable, dimension(:,:,:) :: myvar
  real(kind=PRC), allocatable, dimension(:,:,:), optional :: myvary,myvarz
  integer :: e_io,nnt
  real(kind=4), allocatable, dimension(:,:,:) :: service1,service2,service3
  
  character(len=120) :: sevt
  integer :: l,ii,jj,kk,i,j,k,myindex,iost,endoff
  
  integer :: nx1,nx2,ny1,ny2,nz1,nz2,nn,nxmin,nymin,nnn
  
  character(len=120) :: s_buffer
  
  character(len=500) :: mystring500
  integer, parameter :: byter4=4
  integer, parameter :: byteint=4
  real(kind=PRC) :: miomax,vmod
  real(kind=4) :: r4temp,b4temp
  
  if(.not. lvtkfile)return
  
  if(mod(nstepsub,ivtkevery)/=0)return
  
  nx1=minx
  nx2=maxx
  ny1=miny
  ny2=maxy
  nz1=minz
  nz2=maxz
  
  nnn=nx*ny*nz
  nn=(nx2-nx1+1)*(ny2-ny1+1)*(nz2-nz1+1)
  
  
  
  if(ndim==1)then
    allocate(service1(nx1:nx2,ny1:ny2,nz1:nz2))
    service1(nx1:nx2,ny1:ny2,nz1:nz2)=ZERO
  !requested scalar field in output
    if(ntype==3)then
       !requested phase field in output
       !compute phase field on fly
       do k=minz,maxz
        do j=miny,maxy
          do i=minx,maxx
            if(isfluid(i,j,k)==3 .or. isfluid(i,j,k)==0)cycle
            r4temp=real(myvar(i,j,k),kind=4)
            b4temp=real(myvary(i,j,k),kind=4)
            service1(i,j,k)=(r4temp-b4temp)/(r4temp+b4temp)
          enddo
        enddo
      enddo
    else
      do k=minz,maxz
        do j=miny,maxy
          do i=minx,maxx
            if(isfluid(i,j,k)==3 .or. isfluid(i,j,k)==0)cycle
            service1(i,j,k)=real(myvar(i,j,k),kind=4)
          enddo
        enddo
      enddo
    endif
  elseif(ndim==3)then
    !requested vector field in output
    allocate(service1(nx1:nx2,ny1:ny2,nz1:nz2))
    allocate(service2(nx1:nx2,ny1:ny2,nz1:nz2))
    allocate(service3(nx1:nx2,ny1:ny2,nz1:nz2))
    service1(nx1:nx2,ny1:ny2,nz1:nz2)=ZERO
    service2(nx1:nx2,ny1:ny2,nz1:nz2)=ZERO
    service3(nx1:nx2,ny1:ny2,nz1:nz2)=ZERO
    do k=minz,maxz
      do j=miny,maxy
        do i=minx,maxx
          if(isfluid(i,j,k)==3 .or. isfluid(i,j,k)==0)cycle
          service1(i,j,k)=real(myvar(i,j,k),kind=4)
          service2(i,j,k)=real(myvary(i,j,k),kind=4)
          service3(i,j,k)=real(myvarz(i,j,k),kind=4)
        enddo
      enddo
    enddo
  endif
  
  sevt=repeat(' ',120)
  
  sevt = trim(dir_out) // trim(filenamevtk)//'_'//trim(filenamevar)// &
   '_'//trim(write_fmtnumb(nstepsub)) // '.vti'
  
  call open_file_vtk_par(ioprint,120,sevt,e_io)
  
  if(idrank==0)call print_header_vtk_par(ioprint,0,nheader,header,E_IO)
  
  endoff=headoff+ndata+byteint
  
  if(idrank==0)call print_footer_vtk_par(ioprint,endoff,footer,E_IO)
  
  if(ndim==1)then
    call driving_print_binary_1d_vtk(ioprint,headoff,ndata,nn,service1, &
     iost)
  elseif(ndim==3)then
    call driving_print_binary_3d_vtk(ioprint,headoff,ndata,nn,service1, &
     service2,service3,iost)
  endif
  
  call close_file_vtk_par(ioprint,e_io)
  
  if(allocated(service1))deallocate(service1)
  if(allocated(service2))deallocate(service2)
  if(allocated(service3))deallocate(service3)
  
  return
     
 end subroutine write_vtk_frame_parallel
 
 subroutine open_file_vtk(iotest,nn,myname,E_IO)
 
!***********************************************************************
!     
!     LBsoft subroutine for opening the vtk legacy file
!     in serial IO
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
  
  open(unit=iotest,file=trim(myname),&
   form='unformatted',access='stream',action='write',status='replace', &
   iostat=e_io)
  
  return
  
 endsubroutine open_file_vtk
 
 subroutine close_file_vtk(iotest,e_io)
 
!***********************************************************************
!     
!     LBsoft subroutine for closing the vtk legacy file
!     in serial IO
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: iotest
  integer, intent(out) :: e_io
   
  close(unit=iotest,iostat=e_io)
  
  return
  
 endsubroutine close_file_vtk
 
 subroutine print_header_vtk(iotest,offsetsub,nn,header,E_IO)
 
!***********************************************************************
!     
!     LBsoft subroutine for writing the header part of
!     in VTK legacy file in serial IO
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: iotest,offsetsub,nn
  character(len=500) :: header
  integer, intent(out) :: E_IO
  
  write(iotest, iostat=E_IO)header(1:nn)
  
  return
  
 endsubroutine print_header_vtk
 
 subroutine print_footer_vtk(iotest,offsetsub,footer,E_IO)
 
!***********************************************************************
!     
!     LBsoft subroutine for writing the footer part of
!     in VTK legacy file in serial IO
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: iotest,offsetsub
  character(len=30) :: footer
  integer, intent(out) :: E_IO
  
  write(iotest, iostat=E_IO)footer(1:30)
  
  return
  
 endsubroutine print_footer_vtk
 
 subroutine initialize_vtk
 
!***********************************************************************
!     
!     LBsoft subroutine for creating the folders containing the files
!     in structured VTK legacy binary format in parallel IO
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
 
 subroutine write_test_map()
  
!***********************************************************************
!     
!     LBsoft subroutine for writing the hydrodynamic variables
!     in ASCII format for diagnostic purposes
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification july 2018
!     
!***********************************************************************

  implicit none
  
  integer :: i,j,k,l
  
  integer, parameter :: iomap=250
  
  character(len=256) :: file_loc_proc

  if(.not. ltest_mode)return
  
  call clean_fluid_inside_particle
  
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
          write(iomap+l,'(3i8,4f24.16,I2)')i,j,k,rhoR(i,j,k), &
           u(i,j,k),v(i,j,k),w(i,j,k),isfluid(i,j,k)
        enddo
      enddo
    enddo
  else
    write(iomap+l,'(8i10)')minx,maxx,miny,maxy,minz,maxz, PRC, 5
    do k=minz,maxz
      do j=miny,maxy
        do i=minx,maxx
          write(iomap+l,'(3i8,5f24.16,I2)')i,j,k,rhoR(i,j,k),rhoB(i,j,k), &
           u(i,j,k),v(i,j,k),w(i,j,k),isfluid(i,j,k)
        enddo
      enddo
    enddo
  endif
  
  close(iomap+l)
  
 end subroutine write_test_map
 
 subroutine write_particle_xyz()
  
!***********************************************************************
!     
!     LBsoft subroutine for writing the particle variables
!     in xyz format for diagnostic purposes
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************

  implicit none
  integer :: i,j,myi
  integer, parameter :: ioxyz=280
  character(len=8), parameter :: mystring8='C       '
  real(kind=PRC) :: dtemp(3), qtemp(0:3)
  
  
  if(.not.lparticles)return
  
  if(idrank==0)then 
    open(ioxyz,file='restart.xyz',status='replace',action='write')
    write(ioxyz,'(i10)')natms_tot
    if(lrotate)then
      write(ioxyz,'(a)')'read list vxx vyy vzz phi theta psi '
    else
      write(ioxyz,'(a)')'read list vxx vyy vzz '
    endif
    close(ioxyz)
  endif
  
  do j=0,mxrank-1
    if(j==idrank)then
      open(ioxyz,file='restart.xyz',status='old',action='write', &
       position='append')
      if(lrotate)then
        do myi=1,natms
          i = atmbook(myi)
          qtemp = [ q0(i),q1(i),q2(i),q3(i) ]
          call q2eul(qtemp,dtemp(1),dtemp(3),dtemp(2))
          write(ioxyz,'(a8,9g16.8,i10)')mystring8,xxx(i),yyy(i),zzz(i), &
           vxx(i),vyy(i),vzz(i),dtemp(1:3),i
        enddo
      else
        do myi=1,natms
          i = atmbook(myi)
          write(ioxyz,'(a8,6g16.8,i10)')mystring8,xxx(i),yyy(i),zzz(i), &
           vxx(i),vyy(i),vzz(i),i
        enddo
      endif
      close(ioxyz)
    endif
    call get_sync_world
  enddo
  
  return
  
 end subroutine write_particle_xyz


 subroutine set_value_ibinevery(ltemp,itemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for settimg time interval value
!     used to print the raw data in binary files
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none
  logical, intent(in) :: ltemp
  integer, intent(in) :: itemp
  
  ibinevery=itemp
  lbinevery = .true.
  
  return
  
 end subroutine set_value_ibinevery

 subroutine set_value_dumpevery(itemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for settimg time interval value
!     used to print the restart files
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none
  integer, intent(in) :: itemp
  
  idumpevery=itemp
  
  return
  
 end subroutine set_value_dumpevery

 subroutine dumpForOutput(nstep,mytime,lprint)
 
!***********************************************************************
!     
!     LBsoft subroutine for writing raw output data for the restarting
!     procedure and binary files
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!*********************************************************************** 
 
  implicit none
  integer, intent(in) :: nstep
  real(kind=PRC), intent(in) :: mytime
  logical, intent(in) :: lprint
  character(len=120) :: mynamefile
  
  if(mod(nstep,idumpevery)==0 .or. lprint) then
    call dump_oneFile(nstep)
    if(lparticles) call dumppart_oneFile(nstep)
    call write_restart_file(idumpevery,135,'dumpGlobal.save.dat', &
     nstep,mytime)
  endif
  
  if(lbinevery)then
    if(mod(nstep,ibinevery)/=0) return
    
    call bin_oneFile(nstep)
  
    if(.not.lparticles)return

    if(idrank==0) then
       mynamefile=repeat(' ',120)
       mynamefile='dumpBin'//delimiter//'dumpParticles.' // write_fmtnumb(nstep) // '.dat'
       open(unit=133,file=trim(mynamefile),form='unformatted',status='replace')
       write(133) xxx(1:natms_tot)
       write(133) yyy(1:natms_tot)
       write(133) zzz(1:natms_tot)
       write(133) vxx(1:natms_tot)
       write(133) vyy(1:natms_tot)
       write(133) vzz(1:natms_tot)
       
       if(lrotate)then
        write(133) q0(1:natms_tot)
        write(133) q1(1:natms_tot)
        write(133) q2(1:natms_tot)
        write(133) q3(1:natms_tot)
       endif
       close(133)
     endif
  endif
  
  return
  
 end subroutine dumpForOutput
 
 subroutine write_restart_file(idumpevery,myunit,myname,nstep,mytime)
  
!***********************************************************************
!     
!     LBsoft subroutine for writing data for the restarting
!     procedure
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: idumpevery,myunit,nstep
  
  character(len=*), intent(in) :: myname
  
  real(kind=PRC), intent(in) :: mytime
 
  
  if(idrank==0)then
    open(unit=myunit,file=trim(myname),status='replace', &
     action='write',form='unformatted')
    
    write(myunit)nstep
    
    close(myunit)
    
  endif
  
  return
  
 end subroutine write_restart_file
 
 subroutine read_restart_file(myunit,myname,nstep_temp)
  
!***********************************************************************
!     
!     LBsoft subroutine for writing data for the restarting
!     procedure
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: myunit
  
  character(len=*), intent(in) :: myname
  
  integer, intent(inout) :: nstep_temp
  
  logical :: lexist
  
  inquire(file=trim(myname),exist=lexist)
  if(.not. lexist)then
    call error(35)
  endif
  
  if(idrank==0)then
    open(unit=myunit,file=trim(myname),status='old', &
     action='read',form='unformatted')
    read(myunit)nstep_temp
    close(myunit)
  endif
  
  call bcast_world_i(nstep_temp)
  
  call warning(53,real(nstep_temp,kind=PRC))
  
  return
  
 end subroutine read_restart_file
 
 subroutine header_vtk(mystring500,namevar,extent,ncomp,iinisub,iend,myoffset, &
  new_myoffset,indent)
 
  implicit none
  
  character(len=8),intent(in) :: namevar
  character(len=120),intent(in) :: extent
  integer, intent(in) :: ncomp,iinisub,myoffset
  integer, intent(out) :: iend,new_myoffset
  integer, intent(inout) :: indent
  
  !namevar='density1'
  
  character(len=500), intent(out) :: mystring500
  ! End-character for binary-record finalize.
  character(1), parameter:: end_rec = char(10) 
  character(1) :: string1
  character(len=*),parameter :: topology='ImageData' 
  integer :: ioffset,nele,bytechar,byteint,byter4,byter8,iini
  
  iini=iinisub
  bytechar=kind(end_rec)
  byteint=kind(iini)
  byter4  = 4
  byter8  = 8
  
  mystring500=repeat(' ',500)
  
  iend=iini
  
  iini=iend+1
  nele=22
  iend=iend+nele
  mystring500(iini:iend)='<?xml version="1.0"?>'//end_rec
  
  new_myoffset=myoffset
  new_myoffset = new_myoffset + nele * bytechar
 
  
  iini=iend+1
  nele=67
  iend=iend+nele
  if(lelittle)then  
    mystring500(iini:iend) = '<VTKFile type="'//trim(topology)// &
     '" version="0.1" byte_order="LittleEndian">'//end_rec
  else
    mystring500(iini:iend) = '<VTKFile type="'//trim(topology)// &
     '" version="0.1" byte_order="BigEndian">   '//end_rec
  endif
  
  new_myoffset = new_myoffset + 67 * bytechar
 
  
  indent = indent + 2
  iini=iend+1
  nele=70
  iend=iend+nele
  mystring500(iini:iend) = repeat(' ',indent)//'<'//trim(topology)//' WholeExtent="'//&
                 trim(extent)//'">'//end_rec
  

  new_myoffset = new_myoffset + 70 * bytechar
 
  
  indent = indent + 2
  iini=iend+1
  nele=63
  iend=iend+nele
  mystring500(iini:iend) = repeat(' ',indent)//'<Piece Extent="'//trim(extent)//'">'//end_rec
  
  new_myoffset = new_myoffset + 63 * bytechar
 
  
! initializing offset pointer
  ioffset = 0 
  
  indent = indent + 2
  iini=iend+1
  nele=18
  iend=iend+nele
  mystring500(iini:iend)=repeat(' ',indent)//'<PointData>'//end_rec
  
  new_myoffset = new_myoffset + 18 * bytechar
  
  indent = indent + 2
  iini=iend+1
  nele=115
  iend=iend+nele
  
  if(ncomp/=1 .and. ncomp/=3)call error(37)
  write(string1,'(i1)')ncomp
  mystring500(iini:iend)=repeat(' ',indent)//'<DataArray type="Float32" Name="'// &
   namevar//'" NumberOfComponents="'//string1// '" '//&
   'format="appended" offset="'//space_fmtnumb12(ioffset)//'"/>'//end_rec
  
  new_myoffset = new_myoffset + 115 * bytechar
 
  
  indent = indent - 2
  iini=iend+1
  nele=19
  iend=iend+nele
  mystring500(iini:iend)=repeat(' ',indent)//'</PointData>'//end_rec
  
  new_myoffset = new_myoffset + 19 * bytechar
  
  
  indent = indent - 2
  iini=iend+1
  nele=13
  iend=iend+nele
  mystring500(iini:iend)=repeat(' ',indent)//'</Piece>'//end_rec
  
  
  new_myoffset = new_myoffset + 13 * bytechar
 
  
  indent = indent - 2
  iini=iend+1
  nele=15
  iend=iend+nele
  mystring500(iini:iend)=repeat(' ',indent)//'</'//trim(topology)//'>'//end_rec

  new_myoffset = new_myoffset + 15 * bytechar
 

  iini=iend+1
  nele=32
  iend=iend+nele
  mystring500(iini:iend)=repeat(' ',indent)//'<AppendedData encoding="raw">'//end_rec
  
  new_myoffset = new_myoffset + 32 * bytechar
  
  iini=iend+1
  nele=1
  iend=iend+nele
  mystring500(iini:iend)='_'
  
  new_myoffset = new_myoffset + 1 * bytechar
  
  return
  
 end subroutine header_vtk
 
 subroutine footer_vtk(mystring30,iinisub,iend,myoffset, &
  new_myoffset,indent)
 
  implicit none
  
  integer, intent(in) :: iinisub,myoffset
  integer, intent(out) :: iend,new_myoffset
  integer, intent(inout) :: indent
  
  
  character(len=30), intent(out) :: mystring30
  ! End-character for binary-record finalize.
  character(1), parameter:: end_rec = char(10) 
  character(1) :: string1
  character(len=*),parameter :: topology='ImageData' 
  integer :: ioffset,nele,bytechar,byteint,byter4,byter8,iini
  
  iini=iinisub
  bytechar=kind(end_rec)
  byteint=kind(iini)
  byter4  = 4
  byter8  = 8
  
  mystring30=repeat(' ',30)
  
  iend=iini
  
  iini=iend+1
  nele=1
  iend=iend+nele
  mystring30(iini:iend)=end_rec
  
  new_myoffset = myoffset
  new_myoffset = new_myoffset + 1 * bytechar
 
  
  
  iini=iend+1
  nele=18
  iend=iend+nele
  mystring30(iini:iend)=repeat(' ',indent)//'</AppendedData>'//end_rec
  
  new_myoffset = new_myoffset + 18 * bytechar
  
  iini=iend+1
  nele=11
  iend=iend+nele
  mystring30(iini:iend)='</VTKFile>'//end_rec
  
  if(iend/=30)then
    call error(39)
  endif
  
  return
  
 end subroutine footer_vtk
 
 end module write_output_mod
