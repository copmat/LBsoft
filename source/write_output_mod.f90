
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
  
 use version_mod,    only : idrank,mxrank,get_sync_world,finalize_world
 use lbempi_mod,     only : gminx,gmaxx,gminy,gmaxy,gminz,gmaxz,ownern,&
  i4back
 use error_mod
 use utility_mod,    only : write_fmtnumb,ltest_mode
 use fluids_mod,     only : nx,ny,nz,rhoR,rhoB,u,v,w,lsingle_fluid, &
  minx, maxx, miny, maxy, minz, maxz,isfluid, dumpHvar
 use particles_mod,  only : natms,xxx,yyy,zzz,lparticles,cell, &
  ishape,lrotate,natms,natms_tot,q0,q1,q2,q3,vxx,vyy,vzz, &
  take_rotversorx,take_rotversory,take_rotversorz, &
  clean_fluid_inside_particle,q2eul, atmbook, natms_ext
 
  private
  
  integer, parameter :: iotest=180
  integer, parameter :: ioxyz=190
  integer, parameter :: mxln=120
  character(len=*), parameter :: filenamevtk='out'
  integer, save, public, protected :: ivtkevery=50
  logical, save, public, protected :: lvtkfile=.false.
  logical, save, public, protected :: lvtkstruct=.false.
  logical, save, public, protected :: lvtkownern=.false.
  integer, save, public, protected :: ixyzevery=50
  logical, save, public, protected :: lxyzfile=.false.
  integer, save, public, protected :: istatevery=5000000
  logical, save, public, protected :: lstatevery=.false.
  character(len=mxln), save :: dir_out
  character(len=mxln), save :: dir_out_rank
  character(len=mxln), allocatable, save :: dir_out_grank(:)
  
  character, save :: delimiter
  
  public :: write_vtk_frame
  public :: init_output
  public :: set_value_ivtkevery
  public :: write_test_map
  public :: write_xyz
  public :: write_xyz_close
  public :: set_value_ixyzevery
  public :: write_particle_xyz
  public :: set_value_istatevery, dumpForStats
  public :: writePVD
  
 contains
 
 subroutine init_output()
 
!***********************************************************************
!     
!     LBsoft subroutine for creating the folders containing the files
!     in image VTK legacy binary format in parallel IO
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
  
  
  character(len=255) :: path,makedirectory
  logical :: lexist
  
  integer :: i
  
  if((.not. lvtkfile).and.(.not. lxyzfile) .and. (.not. lstatevery))return
  
  path = repeat(' ',255)
  call getcwd(path)
  
  !call get_environment_variable('DELIMITER',delimiter)
  path = trim(path)
  delimiter = path(1:1)
  if (delimiter==' ') delimiter='/'
  
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

  makedirectory=repeat(' ',255)
  makedirectory = 'dumpStat'
  dir_out=trim(makedirectory)
#ifdef INTEL
  inquire(directory=trim(makedirectory),exist=lexist)
#else
  inquire(file=trim(makedirectory),exist=lexist)
#endif

  if(.not. lexist)then
    if(idrank==0) then
      makedirectory=repeat(' ',255)
      makedirectory = 'mkdir dumpStat'
      call system(makedirectory)
    endif
  endif
  
  call get_sync_world()
  
  call write_xyz_open(trim(dir_out)//'output.xyz')
  
  if(.not. lvtkfile)return
  
  makedirectory=repeat(' ',255)
  makedirectory=trim(path)//delimiter//'output'//delimiter
  

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
  
  call get_sync_world
  
  return

 end subroutine init_output

 subroutine writeImageDataVTI(fname)
 
!***********************************************************************
!     
!     LBsoft subroutine for writing the hydrodynamic variables
!     in image VTK legacy binary format
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: F. Bonaccorso
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
  character(len=120),intent(in) :: fname
  character(len=120) :: fnameFull,extent
  integer i,j,k
  
1000 format ("     ", f15.8)
1001 format (f20.8)
1003 format (3f20.8)
  
  fnameFull = trim(dir_out_rank)// trim(fname) // '_' // trim(write_fmtnumb(idrank)) //'.vti'
  open(unit=iotest,file=trim(fnameFull),status='replace',action='write')

  extent =  trim(write_fmtnumb(minx)) // ' ' // trim(write_fmtnumb(maxx+1)) // ' ' &
        // trim(write_fmtnumb(miny)) // ' ' // trim(write_fmtnumb(maxy+1)) // ' ' &
        // trim(write_fmtnumb(minz)) // ' ' // trim(write_fmtnumb(maxz+1))

  write(iotest,'(a)') '<VTKFile type="ImageData" version="1.0">'
  write(iotest,'(a)') ' <ImageData WholeExtent="' // trim(extent) // '" >'
  write(iotest,'(a)') ' <Piece Extent="' // trim(extent) // '">'
  write(iotest,'(a)') '   <CellData>'
  write(iotest,'(a)') '    <DataArray type="Float32" Name="density1" format="ascii" >'

  do k=minz,maxz
     do j=miny,maxy
       do i=minx,maxx
         write(iotest,fmt=1001) rhoR(i,j,k)
       enddo
     enddo
  enddo
  write(iotest,'(a)') '    </DataArray>'
  if(.not. lsingle_fluid)then
  write(iotest,'(a)') '    <DataArray type="Float32" Name="density2" format="ascii" >'

  do k=minz,maxz
     do j=miny,maxy
       do i=minx,maxx
         write(iotest,fmt=1001) rhoB(i,j,k)
       enddo
     enddo
  enddo
  write(iotest,'(a)') '    </DataArray>'
  endif
  write(iotest,'(2a)') '    <DataArray type="Float32" Name="velocity" ', &
   'NumberOfComponents="3" format="ascii" >'

  do k=minz,maxz
     do j=miny,maxy
       do i=minx,maxx
         write(iotest,fmt=1003) u(i,j,k),v(i,j,k),w(i,j,k)
       enddo
     enddo
  enddo
  write(iotest,'(a)') '    </DataArray>'
  
  write(iotest,'(a)') '   </CellData>'
  write(iotest,'(a)') ' </Piece>'
  write(iotest,'(a)') ' </ImageData>'
  write(iotest,'(a)') '</VTKFile >'
  close(iotest)
  
  return
  
 end subroutine writeImageDataVTI

 subroutine writeImageDataPVTI(fname)
 
!***********************************************************************
!     
!     LBsoft subroutine for writing the parallel driving file
!     for image VTK legacy output
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: F. Bonaccorso
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
  
  character(len=120),intent(in) :: fname
  character(len=120) :: fnameFull,extent
  character(len=255) :: makedirectory
  integer i
  
  if(idrank/=0)return
  
  fnameFull = trim(dir_out) // trim(fname) // '.pvti'
  open(unit=iotest,file=trim(fnameFull),status='replace',action='write')

  extent =  '1 ' // trim(write_fmtnumb(nx+1)) // ' ' &
        // '1 ' // trim(write_fmtnumb(ny+1)) // ' ' &
        // '1 ' // trim(write_fmtnumb(nz+1))

  write(iotest,'(a)') '<VTKFile type="PImageData" version="1.0">'
  write(iotest,'(a)') '  <PImageData WholeExtent="' // trim(extent) // '">'
  write(iotest,'(a)') '   <PCellData>'
  write(iotest,'(a)') '    <PDataArray type="Float32" Name="density1" />'
  if(.not. lsingle_fluid)then
  write(iotest,'(a)') '    <PDataArray type="Float32" Name="density2" />'
  endif
  write(iotest,'(2a)')'    <PDataArray type="Float32" Name="velocity" ', &
   'NumberOfComponents="3"/>'
  write(iotest,'(a)') '   </PCellData>'

  do i=0,mxrank-1
   extent =  trim(write_fmtnumb(gminx(i))) // ' ' // trim(write_fmtnumb(gmaxx(i)+1)) // ' ' &
        // trim(write_fmtnumb(gminy(i))) // ' ' // trim(write_fmtnumb(gmaxy(i)+1)) // ' ' &
        // trim(write_fmtnumb(gminz(i))) // ' ' // trim(write_fmtnumb(gmaxz(i)+1))
    makedirectory=repeat(' ',255)
    makedirectory = 'rank'//trim(write_fmtnumb(i))//delimiter
    write(iotest,'(a)') '    <Piece Extent="' // trim(extent) // '" Source="' // &
    trim(makedirectory)//trim(fname) // '_' // trim(write_fmtnumb(i)) //'.vti" />'
  enddo

  write(iotest,'(a)') '  </PImageData>'
  write(iotest,'(a)') '</VTKFile>'
  close(iotest)
 
  return
 
 end subroutine writeImageDataPVTI

 subroutine writeImageDataPVD(fname)
 
!***********************************************************************
!     
!     LBsoft subroutine for writing the collection driving PVD file
!     for VTK legacy output
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: F. Bonaccorso
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
  
  character(len=120),intent(in) :: fname

  character(len=120) :: fnameFull

  fnameFull = trim(fname) // '.pvd'

  open(unit=iotest,file=trim(fnameFull),status='replace',action='write')

  write(iotest,'(a)') '<VTKFile type="Collection" version="1.0">'
  write(iotest,'(a)') '  <Collection>'
  write(iotest,'(a)') '   <DataSet part="0" file="output/' // trim(fname) // '.pvti"/>'
  write(iotest,'(a)') '  </Collection>'
  write(iotest,'(a)') '</VTKFile>'
  close(iotest)
  
  return
  
 end subroutine writeImageDataPVD
 
 subroutine write_vtp_file(iunitsub,fnamesub,istepsub)
 
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
  
  sevt=repeat(' ',120)
  sevt=trim(dir_out)//trim(fnamesub)// &
   trim(write_fmtnumb(istepsub))//'.vtk'
   
  
  open(unit=iotest,file=trim(sevt),status='replace',action='write')
  
  write(iotest,103) '# vtk DataFile Version 2.0'
  write(iotest,103) 'Field Emission Device - Charge Density Plot'
  write(iotest,103) 'ASCII'
  write(iotest,103) 'DATASET POLYDATA'
  write(iotest,106) 'POINTS ',natms,' float'
  do iatm=1,natms
    write(iotest,108)xxx(iatm),yyy(iatm),zzz(iatm)
  enddo
  
  write(iotest,*)
  write(iotest,107) 'POINT_DATA ', natms
  write(iotest,103) 'SCALARS Type float 1'
  write(iotest,103) 'LOOKUP_TABLE default'
  
  do iatm=1,natms
    write(iotest,*)1 !ltype(iatm)
  enddo
  
  if(lrotate)then
    
    write(iotest,103) 'VECTORS Vectx float'
    do iatm=1,natms
      write(iotest,*)take_rotversorx(q0(iatm),q1(iatm),q2(iatm),q3(iatm))
    enddo
    
    write(iotest,103) 'VECTORS Vecty float'
    do iatm=1,natms
      write(iotest,*)take_rotversory(q0(iatm),q1(iatm),q2(iatm),q3(iatm))
    enddo
    
    write(iotest,103) 'VECTORS Vectz float'
    do iatm=1,natms
      write(iotest,*)take_rotversorz(q0(iatm),q1(iatm),q2(iatm),q3(iatm))
    enddo
    
  endif
  
  close(iotest)
   
  return
  
 end subroutine

 subroutine writePVD(nstepsub)
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the write of the PVD file
!     for VTK legacy output
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: F. Bonaccorso
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: nstepsub
  character(len=120) :: fname

  if(mod(nstepsub,ivtkevery)/=0)return
  
  fname=repeat(' ',120)
  fname = trim(filenamevtk) // trim(write_fmtnumb(nstepsub))
  
  if(idrank==0)then
    !! call writeImageDataPVD(fname)
    call writeImageDataPVTI(fname)
  endif
  
  call writeImageDataVTI(fname)
  
  return
  
 end subroutine writePVD
 
 subroutine write_xyz_open(filename)
 
!***********************************************************************
!     
!     LBsoft subroutine for opening the file in xyz format
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
  
  implicit none
  
  character(len=*), intent(in) :: filename
  
  if(.not. lxyzfile)return
  
  open(ioxyz,file=trim(filename),status='replace',action='write')
  
  return
    
 end subroutine write_xyz_open
 
 subroutine write_xyz(istepsub)
 
!***********************************************************************
!     
!     LBsoft subroutine for writing a frame in the xyz file
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: istepsub
  
  integer :: j
  character(len=8), parameter :: mystring8='C       '
  character(len=13), parameter :: mystring13='time step    '
  
  if(.not. lxyzfile)return
  if(mod(istepsub,ixyzevery)/=0)return

  write(ioxyz,'(i8)') natms
  write(ioxyz,'(a13,i12,3f16.8)')mystring13,istepsub,cell(1),cell(5),cell(9)
  do j=1,natms
    write(ioxyz,"(a8,3f16.6)")mystring8,xxx(j),yyy(j),zzz(j)
  end do

  return

    
 end subroutine write_xyz

 subroutine write_xyz_close()
  
!***********************************************************************
!     
!     LBsoft subroutine for closing the file in xyz format
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
  
  implicit none
  
  if(.not. lxyzfile)return
  close(ioxyz)
  
  return
    
 end subroutine write_xyz_close

 subroutine set_value_ivtkevery(ltemp,ltemp2,itemp)
 
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
  
  logical, intent(in) :: ltemp,ltemp2
  integer, intent(in) :: itemp
  
  lvtkfile=ltemp
  lvtkstruct=ltemp2
  ivtkevery=itemp
  
  return
  
 end subroutine set_value_ivtkevery
 
 subroutine set_value_ixyzevery(ltemp,itemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for settimg time interval value
!     used to print the xyz file
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
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
  
  if((.not. lvtkfile).and.(.not. lxyzfile))return
  
  if(mod(nstepsub,ivtkevery)/=0)return
  
  call clean_fluid_inside_particle
  
  if(lparticles)then
    if(mxrank==1)then
      call write_vtp_file(181,'outatm',nstepsub)
!    else
!      call error(11)
!       write (6,*) "Skipping in Parallel"
    endif
  endif
  
  if(lvtkstruct)then
    if(mxrank==1)then
      call write_vtk_frame_serial(nstepsub)
    else
      call write_vtk_frame_parallel(nstepsub)
    endif
  else
    call writePVD(nstepsub)
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
  
  real(kind=PRC) :: miomax,vmod
  
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
  !sevt=trim(filenamevtk)// &
  ! trim(write_fmtnumb(nstepsub))//'.vts'
  sevt = trim(dir_out) // trim(filenamevtk)// &
   trim(write_fmtnumb(nstepsub)) // '.vts'
 
  E_IO = VTK_INI_XML(output_format='raw', &
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
   
#if 0
  service1(nx1:nx2,ny1:ny2,nz1:nz2)=real(isfluid(1:nx,1:ny,1:nz),kind=R4P)
  E_IO = VTK_VAR_XML(NC_NN=nn,varname ='isfluid', &
   var=reshape(service1,(/nn/))) 
#endif
   
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
    !call initialize_vtk
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
    !sevt=trim(filenamevtk)//trim(write_fmtnumb(nstepsub))//'.pvts'
    sevt=trim(dir_out)//trim(filenamevtk)//trim(write_fmtnumb(nstepsub)) &
     // '.pvts'
    ! pvts
    E_IO = PVTK_INI_XML(filename = trim(sevt), &
     mesh_topology = 'PStructuredGrid', &
     nx1=nx1,nx2=nx2,ny1=ny1,ny2=ny2,nz1=nz1,nz2=nz2, tp='Float32')
    do i=0,mxrank-1
      sevt=repeat(' ',120)
      !sevt='rank'//trim(write_fmtnumb(i))//delimiter//trim(filenamevtk)// &
      ! trim(write_fmtnumb(nstepsub))//'_part'//trim(write_fmtnumb(i))//'.vts'
      sevt='rank'//trim(write_fmtnumb(i))//delimiter//trim(filenamevtk)// &
       trim(write_fmtnumb(nstepsub))//'_'//trim(write_fmtnumb(i))//'.vts'
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
  !sevt='rank'//trim(write_fmtnumb(idrank))//delimiter//trim(filenamevtk)// &
  ! trim(write_fmtnumb(nstepsub))//'_part'//trim(write_fmtnumb(idrank))//'.vts'
  sevt=trim(dir_out_rank)//trim(filenamevtk)//trim(write_fmtnumb(nstepsub)) &
   // '_' // trim(write_fmtnumb(idrank)) //'.vts'
     
  E_IO = VTK_INI_XML(cf=mf(idrank),output_format='raw', &
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
!     licensed under Open Software License v. 3.0 (OSL-3.0)
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


 subroutine set_value_istatevery(itemp)
  implicit none
  integer, intent(in) :: itemp
  
  istatevery=itemp
  lstatevery = .true.
 end subroutine set_value_istatevery

 subroutine dumpForStats(nstep)
  implicit none
  integer, intent(in) :: nstep
  character(len=120) :: mynamefile
  
  if(mod(nstep,istatevery)/=0)return
  
  call dumpHvar(nstep)

  if(.not.lparticles)return

  if(idrank==0) then
     mynamefile=repeat(' ',120)
     mynamefile='dumpStat'//delimiter//'dumpParticles.' // write_fmtnumb(nstep) // '.dat'
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
 end subroutine dumpForStats
 
 end module write_output_mod
