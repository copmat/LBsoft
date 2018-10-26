
#include <default_macro.h>
 module particles_mod
 
!***********************************************************************
!     
!     LBsoft module containing variable and subroutines of the
!     particle managing
!
!     few variable names and subroutines are inspired to
!     DL POLY CLASSICS distributed under BSD licence from
!     the daresbury laboratory (in primis we like to acknowledge 
!     prof. w. smith)
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
 use version_mod, only : idrank,mxrank,or_world_larr,finalize_world, &
                   get_sync_world,isend_world_farr,bcast_world_iarr, &
                   waitall_world,irecv_world_farr,isend_world_iarr, &
                   irecv_world_iarr,max_world_iarr,sum_world_iarr, &
                   sum_world_farr
 use error_mod
 use aop_mod
 use utility_mod, only : Pi,modulvec,cross,dot,gauss,ibuffservice, &
                   allocate_array_ibuffservice,nbuffservice,buffservice, &
                   allocate_array_buffservice,lbuffservice, &
                   allocate_array_lbuffservice, &
                   nbuffservice3d,buffservice3d,allocate_array_bdf, &
                   xdf,ydf,zdf,rand_noseeded,linit_seed,gauss_noseeded,&
                   write_fmtnumb,dcell,invert

 use lbempi_mod,  only : commspop, commrpop, i4find, i4back, &
                   ownern,deallocate_ownern,commexch_dens, &
                   commwait_dens,comm_init_isfluid,commexch_vel_component, &
                   commwait_vel_component,ownernfind,ownernfind_arr, &
                   ownern,gminx,gmaxx,gminy,gmaxy,gminz,gmaxz
 
 use fluids_mod,  only : nx,ny,nz,nbuff,minx,maxx,miny,maxy,minz,maxz

 
 implicit none
 
 private
 
 !key for periodic boundary conditions
 integer, public, protected, save :: imcon=0
 
 !number of particles in the subdomain (process)
 integer, public, protected, save :: natms=0
 !number of particles in the subdomain (process) including the halo
 integer, public, protected, save :: natms_ext=0
 !total number of particles in all the domain
 integer, public, protected, save :: natms_tot=0
 
 !expected total number of particles in each the subdomain
 integer, public, protected, save :: mxatms=0
 
 !expected number of particles in each neighborhood list
 integer, public, protected, save :: mxlist=0
 
 !expected local size of the neighborhood list
 integer, public, protected, save :: msatms=0
 
 !expected number of particles in each neighborhood listcell
 integer, public, protected, save :: mslistcell=0
 
 !number of subdomains in parallel decomposition
 integer, public, protected, save :: nbig_cells=0
 
 !max number of parameters for pair force fields
 integer, public, parameter :: mxpvdw=2
 
 !max number of pair force fields
 integer, public, parameter :: mxvdw=1
 
 !number of pair force fields defined in input
 integer, public, protected, save :: ntpvdw=0
 
 !activate particle part
 logical, public, protected, save :: lparticles=.false.
 
 !set unique mass
 logical, public, protected, save :: lumass=.false.
 
 !time step of the lagrangian integrator
 real(kind=PRC), public, protected, save :: tstepatm=ONE
 
 !value of unique mass
 real(kind=PRC), public, protected, save :: umass=ONE
 
 !cell vectors
 real(kind=PRC), dimension(9), public, protected, save :: cell
 
 !expected variance of the particle density given as ratio on natms_tot
 real(kind=PRC), public, protected, save :: densvar=ELEVEN/TEN
 real(kind=PRC), public, protected, save :: cx=ZERO
 real(kind=PRC), public, protected, save :: cy=ZERO
 real(kind=PRC), public, protected, save :: cz=ZERO
 
 !kinetic energy
 real(kind=PRC), public, protected, save :: engke
 
 !rotational kinetic energy
 real(kind=PRC), public, protected, save :: engrot
 
 !configurational energy
 real(kind=PRC), public, protected, save :: engcfg
 
 !Potential cut-off
 real(kind=PRC), public, protected, save :: rcut=ZERO
 !Verlet neighbour list shell width
 real(kind=PRC), public, protected, save :: delr=ZERO
 
 !key for activate the body rotation
 logical, public, protected, save :: lrotate=.false.
 
 !key for activate the body rotation
 logical, public, protected, save :: lspherical=.true.
 
 !there are enough cells for the link approach within the same process?
 logical, public, protected, save :: lnolink=.false.
 
 !book of global particle ID
 integer, allocatable, public, protected, save :: atmbook(:)
 
 !number of particle ID in the link cells within the same process
 integer, allocatable, public, protected, save :: nlinklist(:)
 !list of particle ID in the link cells within the same process
 integer, allocatable, public, protected, save :: linklist(:,:)
 
 !number of link cell within the same process
 integer, save :: ncells=0
 
 !global minimum number of link cell within the same process
 integer, save :: ncellsmin=0
 
 !global maximum number of link cell within the same process
 integer, save :: ncellsmax=0
 
 !verlet list
 integer, allocatable, save :: lentry(:)
 integer, allocatable, save :: list(:,:)
 
 !store the type of pair force field
 integer, allocatable, public, protected, save :: ltpvdw(:)
 
 !position of particles
 real(kind=PRC), allocatable, public, protected, save :: xxx(:)
 real(kind=PRC), allocatable, public, protected, save :: yyy(:)
 real(kind=PRC), allocatable, public, protected, save :: zzz(:)
 !components of angular velocity
 real(kind=PRC), allocatable, public, protected, save :: vxx(:)
 real(kind=PRC), allocatable, public, protected, save :: vyy(:)
 real(kind=PRC), allocatable, public, protected, save :: vzz(:)
 !components of force
 real(kind=PRC), allocatable, public, protected, save :: fxx(:)
 real(kind=PRC), allocatable, public, protected, save :: fyy(:)
 real(kind=PRC), allocatable, public, protected, save :: fzz(:)
 !radius of particles
 real(kind=PRC), allocatable, public, protected, save :: rdim(:)
 !radius of particles
 real(kind=PRC), allocatable, public, protected, save :: rdimx(:)
 real(kind=PRC), allocatable, public, protected, save :: rdimy(:)
 real(kind=PRC), allocatable, public, protected, save :: rdimz(:)
 !particle mass
 real(kind=PRC), allocatable, public, protected, save :: weight(:)
 !rotational inertia in body fixed frame
 real(kind=PRC), allocatable, public, protected, save :: rotin(:)
 !components of angular velocity
 real(kind=PRC), allocatable, public, protected, save :: oxx(:)
 real(kind=PRC), allocatable, public, protected, save :: oyy(:)
 real(kind=PRC), allocatable, public, protected, save :: ozz(:)
 !components of quaternion vector
 real(kind=PRC), allocatable, public, protected, save :: q0(:)
 real(kind=PRC), allocatable, public, protected, save :: q1(:)
 real(kind=PRC), allocatable, public, protected, save :: q2(:)
 real(kind=PRC), allocatable, public, protected, save :: q3(:)
 !components of torque on rigid body
 real(kind=PRC), allocatable, public, protected, save :: tqx(:)
 real(kind=PRC), allocatable, public, protected, save :: tqy(:)
 real(kind=PRC), allocatable, public, protected, save :: tqz(:)
 
 !parameters of all the pair force fields
 real(kind=PRC), allocatable, public, protected, save :: prmvdw(:,:)
 
 public :: allocate_particles
 public :: set_natms_tot
 public :: set_lparticles
 public :: set_densvar
 public :: set_rcut
 public :: set_delr
 public :: set_lrotate
 public :: initialize_map_particles
 public :: rotmat_2_quat
 public :: driver_neighborhood_list
 public :: set_ntpvdw
 public :: allocate_field_array
 public :: set_field_array
 public :: vertest
 public :: set_umass
 public :: initialize_particle_force
 public :: driver_inter_force
 public :: initialize_lf
 public :: nve_lf
 
 contains
 
 subroutine set_natms_tot(itemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the global number of particles 
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: itemp
  
  natms_tot=itemp
  
  return
  
 end subroutine set_natms_tot
 
 subroutine set_lparticles(ltemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the lparticles protected variable
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
  
  logical, intent(in) :: ltemp
  
  lparticles=ltemp
  
  return
  
 end subroutine set_lparticles
 
 subroutine set_densvar(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the densvar protected variable
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: dtemp
  
  densvar=dtemp
  
  return
  
 end subroutine set_densvar
 
 subroutine set_rcut(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the rcut protected variable
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: dtemp
  
  rcut=dtemp
  
  return
  
 end subroutine set_rcut
 
 subroutine set_delr(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the delr protected variable
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: dtemp
  
  delr=dtemp
  
  return
  
 end subroutine set_delr
 
 subroutine set_lrotate(ltemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the lrotate protected variable
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
  
  logical, intent(in) :: ltemp
  
  lrotate=ltemp
  
  return
  
 end subroutine set_lrotate
 
 subroutine set_umass(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the ntpvdw protected variable
!     ntpvdw denotes the total number of pair force fields 
!     defined in input
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: dtemp
  
  lumass=.true.
  umass=dtemp
  
  return
  
 end subroutine set_umass
 
 subroutine set_ntpvdw(itemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the ntpvdw protected variable
!     ntpvdw denotes the total number of pair force fields 
!     defined in input
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: itemp
  
  ntpvdw=itemp
  
  return
  
 end subroutine set_ntpvdw
 
 subroutine allocate_field_array
 
!***********************************************************************
!     
!     LBsoft subroutine for allocate the prmvdw protected variable
!     prmvdw stores all the parameters of the pair force fields 
!     defined in input
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
  
  integer :: fail(2)
  
  allocate(prmvdw(mxpvdw,ntpvdw),stat=fail(1))
  allocate(ltpvdw(ntpvdw),stat=fail(2))
  
  call sum_world_iarr(fail,2)
  if(any(fail.ne.0))call error(26)
  
  return
  
 end subroutine allocate_field_array
 
 subroutine set_field_array(iarr1,iarr2,itemp1,dtemp2)
 
!***********************************************************************
!     
!     LBsoft subroutine for allocate the prmvdw protected variable
!     prmvdw stores all the parameters of the pair force fields 
!     defined in input
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: iarr1,iarr2
  integer, dimension(iarr1), intent(in) :: itemp1
  real(kind=PRC), dimension(iarr2,iarr1), intent(in) :: dtemp2
  
  integer :: i
  
  do i=1,iarr1
    ltpvdw(i)=itemp1(i)
    prmvdw(1:iarr2,i)=dtemp2(1:iarr2,i)
  enddo
  
  return
  
 end subroutine set_field_array
 
 subroutine allocate_particles(ibctype,tstep)
 
!***********************************************************************
!     
!     LBsoft subroutine for allocate the particle variables
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: ibctype
  real(kind=PRC), intent(in) :: tstep
  
  integer :: i
  integer, parameter :: nistatmax=100
  integer, dimension(nistatmax) :: istat
  real(kind=PRC) :: cellx,celly,cellz,volm,dens,ratio
  logical :: ltest(1)
  integer :: ilx,ily,ilz
  
  if(.not. lparticles)return
  
  tstepatm=tstep
  
  imcon=ibctype
  
  mxatms=ceiling(real(natms_tot,kind=PRC)*densvar)
  volm=real(nx+1,kind=PRC)*real(nx+1,kind=PRC)*real(nx+1,kind=PRC)
  dens=dble(mxatms)/volm
  ratio=(ONE+HALF)*dens*(FOUR*pi/THREE)*(rcut+delr)**THREE
  mxlist=min(nint(ratio),(mxatms+1)/2)
  mxlist=max(mxlist,4)
  
  msatms=ceiling(real(natms_tot/mxrank,kind=PRC)*densvar)
  
  cellx=huge(ONE)
  celly=huge(ONE)
  cellz=huge(ONE)
  do i=0,mxrank-1
    cellx=min(real(gmaxx(i)-gminx(i),kind=PRC)+ONE,cellx)
    celly=min(real(gmaxy(i)-gminy(i),kind=PRC)+ONE,celly)
    cellz=min(real(gmaxz(i)-gminz(i),kind=PRC)+ONE,cellz)
  enddo
  
  ilx=int(cellx/(rcut+delr))
  ily=int(celly/(rcut+delr))
  ilz=int(cellz/(rcut+delr))
  !check there are enough link cells
  lnolink=.false.
  if(ilx.lt.3)lnolink=.true.
  if(ily.lt.3)lnolink=.true.
  if(ilz.lt.3)lnolink=.true.
  ncellsmin=ilx*ily*ilz
  
  if(lnolink)then
    call warning(21)
  else
    !cordination number of spheres in bcc and fcc is 12, thus....
    mslistcell=max(((mxatms/mxrank)/ncellsmin+1),12)
  endif
  
#if 0
  if(.not. lnolink)then
    lnolink=.true.
    if(idrank==0)then
      write(6,'(a)')'ATTENTION: link cell within the sub domain is under developing!'
      write(6,'(a)')'ATTENTION: link cell will be applied only over the sub domains!'
    endif
  endif
#endif
  
  cellx=ZERO
  celly=ZERO
  cellz=ZERO
  do i=0,mxrank-1
    cellx=max(real(gmaxx(i)-gminx(i),kind=PRC)+ONE,cellx)
    celly=max(real(gmaxy(i)-gminy(i),kind=PRC)+ONE,celly)
    cellz=max(real(gmaxz(i)-gminz(i),kind=PRC)+ONE,cellz)
  enddo
  
  ilx=int(cellx/(rcut+delr))
  ily=int(celly/(rcut+delr))
  ilz=int(cellz/(rcut+delr))
  ncellsmax=ilx*ily*ilz
  
  nbig_cells=mxrank
  
  cx=real(nx+1,kind=PRC)*HALF
  cy=real(ny+1,kind=PRC)*HALF
  cz=real(nz+1,kind=PRC)*HALF
  
  cell(1)=real(nx,kind=PRC)
  cell(2)=ZERO
  cell(3)=ZERO
  cell(4)=ZERO
  cell(5)=real(ny,kind=PRC)
  cell(6)=ZERO
  cell(7)=ZERO
  cell(8)=ZERO
  cell(9)=real(nz,kind=PRC)
  
  istat(1:nistatmax)=0
  
  allocate(atmbook(mxatms),stat=istat(100))
  
  allocate(xxx(mxatms),stat=istat(1))
  allocate(yyy(mxatms),stat=istat(2))
  allocate(zzz(mxatms),stat=istat(3))
  
  allocate(vxx(mxatms),stat=istat(4))
  allocate(vyy(mxatms),stat=istat(5))
  allocate(vzz(mxatms),stat=istat(6))
  
  allocate(fxx(mxatms),stat=istat(7))
  allocate(fyy(mxatms),stat=istat(8))
  allocate(fzz(mxatms),stat=istat(9))
  
  if(lspherical)then
    allocate(rdim(mxatms),stat=istat(10))
  else
    allocate(rdimx(mxatms),stat=istat(23))
    allocate(rdimy(mxatms),stat=istat(24))
    allocate(rdimz(mxatms),stat=istat(25))
  endif
  allocate(weight(mxatms),stat=istat(11))
  allocate(rotin(mxatms),stat=istat(12))
  
  allocate(oxx(mxatms),stat=istat(13))
  allocate(oyy(mxatms),stat=istat(14))
  allocate(ozz(mxatms),stat=istat(15))
  
  allocate(q0(mxatms),stat=istat(16))
  allocate(q1(mxatms),stat=istat(17))
  allocate(q2(mxatms),stat=istat(18))
  allocate(q3(mxatms),stat=istat(19))
  
  allocate(tqx(mxatms),stat=istat(20))
  allocate(tqy(mxatms),stat=istat(21))
  allocate(tqz(mxatms),stat=istat(22))
  
#if 1
  if(mxrank>1)then
    if(idrank==0)write(6,'(a)')'ATTENTION: MD in parallel is under developing!'
    call error(-1)
  endif
  
  if(lrotate)then
    if(idrank==0)write(6,'(a)')'ATTENTION: MD rotation integrator is under developing!'
    call error(-1)
  endif
  
  !at the moment only verlet list
  if(.not. lnolink)then
    if(idrank==0)write(6,'(a)')'ATTENTION: at the moment link cell is not implemente!'
    if(idrank==0)write(6,'(a)')'ATTENTION: automatic switch to Verlet list mode!'
    lnolink=.true.
  endif
#endif
  
  if(lnolink)then
    allocate(lentry(msatms),stat=istat(23))
    allocate(list(mxlist,msatms),stat=istat(24))
  else
    allocate(nlinklist(ncellsmax),stat=istat(23))
    allocate(linklist(mslistcell,ncellsmax),stat=istat(24))
  endif
  
  ltest=.false.
  if(any(istat.ne.0))then
    do i=1,nistatmax
      if(istat(i).ne.0)exit
    enddo
    call warning(2,dble(i))
    ltest=.true.
  endif

  call or_world_larr(ltest,1)
  if(ltest(1))call error(20)
  
  atmbook(1:mxatms)=0
  
  xxx(1:mxatms)=ZERO
  yyy(1:mxatms)=ZERO
  zzz(1:mxatms)=ZERO
  
  vxx(1:mxatms)=ZERO
  vyy(1:mxatms)=ZERO
  vzz(1:mxatms)=ZERO
  
  fxx(1:mxatms)=ZERO
  fyy(1:mxatms)=ZERO
  fzz(1:mxatms)=ZERO
  
  if(lspherical)then
    rdim(1:mxatms)=ZERO
  else
    rdimx(1:mxatms)=ZERO
    rdimy(1:mxatms)=ZERO
    rdimz(1:mxatms)=ZERO
  endif
  weight(1:mxatms)=ZERO
  rotin(1:mxatms)=ZERO
  
  oxx(1:mxatms)=ZERO
  oyy(1:mxatms)=ZERO
  ozz(1:mxatms)=ZERO
  
  q0(1:mxatms)=ZERO
  q1(1:mxatms)=ZERO
  q2(1:mxatms)=ZERO
  q3(1:mxatms)=ZERO
  
  tqx(1:mxatms)=ZERO
  tqy(1:mxatms)=ZERO
  tqz(1:mxatms)=ZERO
  
  if(lnolink)then
    lentry(1:msatms)=0
    list(1:mxlist,1:msatms)=0
  else
    nlinklist(1:ncellsmax)=0
    linklist(1:mslistcell,1:ncellsmax)=0
  endif
  
  return
 
 end subroutine allocate_particles
 
 subroutine initialize_map_particles(nxyzlist_sub,xyzlist_sub, &
  xxs,yys,zzs,ots)
  
!***********************************************************************
!     
!     LBsoft subroutine for mapping all the particles and send
!     their values to the appropiate process following the LB
!     domain decomposition
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: nxyzlist_sub
  integer, intent(in), allocatable, dimension(:)  :: xyzlist_sub
  real(kind=PRC), allocatable, dimension(:) :: xxs,yys,zzs
  real(kind=PRC), allocatable, dimension(:,:) :: ots
  
  integer :: i,j,k,ids,sub_i,idrank_sub
  logical :: ltest(1),lqinput,ltinput
  integer, dimension(0:mxrank-1) :: isend_nparticle
  integer :: isend_nvar
  integer, dimension(4+nxyzlist_sub,1:mxrank-1) :: irequest_send
  integer, dimension(4+nxyzlist_sub) :: irequest_recv
  real(kind=PRC), dimension(9) :: rot,newrot
  real(kind=PRC) :: matrixmio(3,3),dmio(3)
  real(kind=PRC), dimension(3) :: xsubm,ysubm,zsubm
#ifdef CHECKQUAT
  real(kind=PRC), parameter :: toll=real(1.d-4,kind=PRC)
  logical :: ltestrot(1)=.false.
#endif
  
  if(.not. lparticles)return
  
  ltest=.false.
  
  isend_nvar=4+nxyzlist_sub
  
  if(idrank==0)then
    
    call pbc_images_centered(imcon,natms_tot,cell,cx,cy,cz,xxs,yys,zzs)
    
    do idrank_sub=0,mxrank-1
      sub_i=0
      do i=1,natms_tot
#ifndef DEOWERN
        ids=ownernfind_arr(nint(xxs(i)),nint(yys(i)),nint(zzs(i)), &
         nx,ny,nz,nbuff,ownern)
#else
        ids=ownernfind(nint(xxs(i)),nint(yys(i)),nint(zzs(i)), &
         mxrank,gminx,gmaxx,gminy,gmaxy,gminz,gmaxz)
#endif
        if(ids==idrank_sub)then
          sub_i=sub_i+1
          if(sub_i>mxatms)then
            ltest(1)=.true.
            goto 460
          endif
        endif
      enddo
      isend_nparticle(idrank_sub)=sub_i
    enddo
  endif
  
 460 continue
 
  call or_world_larr(ltest,1)
  if(ltest(1))then
    call warning(15,densvar)
    call error(21)
  endif
  
  call bcast_world_iarr(isend_nparticle,mxrank)
  
  natms=isend_nparticle(idrank)
  
  if(idrank==0)then
    
    do idrank_sub=1,mxrank-1
      sub_i=0
      do i=1,natms_tot
#ifndef DEOWERN
        ids=ownernfind_arr(nint(xxs(i)),nint(yys(i)),nint(zzs(i)), &
         nx,ny,nz,nbuff,ownern)
#else
        ids=ownernfind(nint(xxs(i)),nint(yys(i)),nint(zzs(i)), &
         mxrank,gminx,gmaxx,gminy,gmaxy,gminz,gmaxz)
#endif
        if(ids==idrank_sub)then
          sub_i=sub_i+1
          atmbook(sub_i)=i
          xxx(sub_i)=xxs(i)
          yyy(sub_i)=yys(i)
          zzz(sub_i)=zzs(i)
          if(nxyzlist_sub>0)then
            do j=1,nxyzlist_sub
              k=xyzlist_sub(j)
              select case(k)
              case(2)
                weight(sub_i)=ots(j,i)
              case(3)
                rdimx(sub_i)=ots(j,i)
              case(4)
                rdimy(sub_i)=ots(j,i)
              case(5)
                rdimz(sub_i)=ots(j,i)
              case(6)
                rdim(sub_i)=ots(j,i)
              case(7)
                vxx(sub_i)=ots(j,i)
              case(8)
                vyy(sub_i)=ots(j,i)
              case(9)
                vzz(sub_i)=ots(j,i)
              case(10)
                oxx(sub_i)=ots(j,i)
              case(11)
                oyy(sub_i)=ots(j,i)
              case(12)
                ozz(sub_i)=ots(j,i)
              case(13)
                fxx(sub_i)=ots(j,i)
              case(14)
                fyy(sub_i)=ots(j,i)
              case(15)
                fzz(sub_i)=ots(j,i)
              case(16)
                tqx(sub_i)=ots(j,i)
              case(17)
                tqy(sub_i)=ots(j,i)
              case(18)
                tqz(sub_i)=ots(j,i)
              case(19)
                q1(sub_i)=ots(j,i)
              case(20)
                q2(sub_i)=ots(j,i)
              case(21)
                q3(sub_i)=ots(j,i)
              case(22)
                q0(sub_i)=ots(j,i)
              case(23)
                q1(sub_i)=ots(j,i)
              case(24)
                q2(sub_i)=ots(j,i)
              case(25)
                q3(sub_i)=ots(j,i)
              end select
            enddo
          endif
        endif
      enddo
      call isend_world_iarr(atmbook,isend_nparticle(idrank_sub),idrank_sub, &
       120+idrank_sub+1,irequest_send(1,idrank_sub))
      call isend_world_farr(xxx,isend_nparticle(idrank_sub),idrank_sub, &
       120+idrank_sub+2,irequest_send(2,idrank_sub))
      call isend_world_farr(yyy,isend_nparticle(idrank_sub),idrank_sub, &
       120+idrank_sub+3,irequest_send(3,idrank_sub))
      call isend_world_farr(zzz,isend_nparticle(idrank_sub),idrank_sub, &
       120+idrank_sub+4,irequest_send(4,idrank_sub))
      do j=1,nxyzlist_sub
        k=xyzlist_sub(j)
        select case(k)
        case(2)
          call isend_world_farr(weight,isend_nparticle(idrank_sub),idrank_sub, &
           120+idrank_sub+4+j,irequest_send(4+j,idrank_sub))
        case(3)
          call isend_world_farr(rdimx,isend_nparticle(idrank_sub),idrank_sub, &
           120+idrank_sub+4+j,irequest_send(4+j,idrank_sub))
        case(4)
          call isend_world_farr(rdimy,isend_nparticle(idrank_sub),idrank_sub, &
           120+idrank_sub+4+j,irequest_send(4+j,idrank_sub))
        case(5)
          call isend_world_farr(rdimz,isend_nparticle(idrank_sub),idrank_sub, &
           120+idrank_sub+4+j,irequest_send(4+j,idrank_sub))
        case(6)
          call isend_world_farr(rdim,isend_nparticle(idrank_sub),idrank_sub, &
           120+idrank_sub+4+j,irequest_send(4+j,idrank_sub))
        case(7)
          call isend_world_farr(vxx,isend_nparticle(idrank_sub),idrank_sub, &
           120+idrank_sub+4+j,irequest_send(4+j,idrank_sub))
        case(8)
          call isend_world_farr(vyy,isend_nparticle(idrank_sub),idrank_sub, &
           120+idrank_sub+4+j,irequest_send(4+j,idrank_sub))
        case(9)
          call isend_world_farr(vzz,isend_nparticle(idrank_sub),idrank_sub, &
           120+idrank_sub+4+j,irequest_send(4+j,idrank_sub))
        case(10)
          call isend_world_farr(oxx,isend_nparticle(idrank_sub),idrank_sub, &
           120+idrank_sub+4+j,irequest_send(4+j,idrank_sub))
        case(11)
          call isend_world_farr(oyy,isend_nparticle(idrank_sub),idrank_sub, &
           120+idrank_sub+4+j,irequest_send(4+j,idrank_sub))
        case(12)
          call isend_world_farr(ozz,isend_nparticle(idrank_sub),idrank_sub, &
           120+idrank_sub+4+j,irequest_send(4+j,idrank_sub))
        case(13)
          call isend_world_farr(fxx,isend_nparticle(idrank_sub),idrank_sub, &
           120+idrank_sub+4+j,irequest_send(4+j,idrank_sub))
        case(14)
          call isend_world_farr(fyy,isend_nparticle(idrank_sub),idrank_sub, &
           120+idrank_sub+4+j,irequest_send(4+j,idrank_sub))
        case(15)
          call isend_world_farr(fzz,isend_nparticle(idrank_sub),idrank_sub, &
           120+idrank_sub+4+j,irequest_send(4+j,idrank_sub))
        case(16)
          call isend_world_farr(tqx,isend_nparticle(idrank_sub),idrank_sub, &
           120+idrank_sub+4+j,irequest_send(4+j,idrank_sub))
        case(17)
          call isend_world_farr(tqy,isend_nparticle(idrank_sub),idrank_sub, &
           120+idrank_sub+4+j,irequest_send(4+j,idrank_sub))
        case(18)
          call isend_world_farr(tqz,isend_nparticle(idrank_sub),idrank_sub, &
           120+idrank_sub+4+j,irequest_send(4+j,idrank_sub))
        case(19)
          call isend_world_farr(q1,isend_nparticle(idrank_sub),idrank_sub, &
           120+idrank_sub+4+j,irequest_send(4+j,idrank_sub))
        case(20)
          call isend_world_farr(q2,isend_nparticle(idrank_sub),idrank_sub, &
           120+idrank_sub+4+j,irequest_send(4+j,idrank_sub))
        case(21)
          call isend_world_farr(q3,isend_nparticle(idrank_sub),idrank_sub, &
           120+idrank_sub+4+j,irequest_send(4+j,idrank_sub))
        case(22)
          call isend_world_farr(q0,isend_nparticle(idrank_sub),idrank_sub, &
           120+idrank_sub+4+j,irequest_send(4+j,idrank_sub))
        case(23)
          call isend_world_farr(q1,isend_nparticle(idrank_sub),idrank_sub, &
           120+idrank_sub+4+j,irequest_send(4+j,idrank_sub))
        case(24)
          call isend_world_farr(q2,isend_nparticle(idrank_sub),idrank_sub, &
           120+idrank_sub+4+j,irequest_send(4+j,idrank_sub))
        case(25)
          call isend_world_farr(q3,isend_nparticle(idrank_sub),idrank_sub, &
           120+idrank_sub+4+j,irequest_send(4+j,idrank_sub))
        end select
      enddo
      call waitall_world(irequest_send(1:isend_nvar,idrank_sub),isend_nvar)
    enddo
    
    sub_i=0
    do i=1,natms_tot
#ifndef DEOWERN
      ids=ownernfind_arr(nint(xxs(i)),nint(yys(i)),nint(zzs(i)), &
       nx,ny,nz,nbuff,ownern)
#else
      ids=ownernfind(nint(xxs(i)),nint(yys(i)),nint(zzs(i)), &
       mxrank,gminx,gmaxx,gminy,gmaxy,gminz,gmaxz)
#endif
      if(ids==0)then
        sub_i=sub_i+1
        atmbook(sub_i)=i
        xxx(sub_i)=xxs(i)
        yyy(sub_i)=yys(i)
        zzz(sub_i)=zzs(i)
        if(nxyzlist_sub>0)then
          do j=1,nxyzlist_sub
            k=xyzlist_sub(j)
            select case(k)
            case(2)
              weight(sub_i)=ots(j,i)
            case(3)
              rdimx(sub_i)=ots(j,i)
            case(4)
              rdimy(sub_i)=ots(j,i)
            case(5)
              rdimz(sub_i)=ots(j,i)
            case(6)
              rdim(sub_i)=ots(j,i)
            case(7)
              vxx(sub_i)=ots(j,i)
            case(8)
              vyy(sub_i)=ots(j,i)
            case(9)
              vzz(sub_i)=ots(j,i)
            case(10)
              oxx(sub_i)=ots(j,i)
            case(11)
              oyy(sub_i)=ots(j,i)
            case(12)
              ozz(sub_i)=ots(j,i)
            case(13)
              fxx(sub_i)=ots(j,i)
            case(14)
              fyy(sub_i)=ots(j,i)
            case(15)
              fzz(sub_i)=ots(j,i)
            case(16)
              tqx(sub_i)=ots(j,i)
            case(17)
              tqy(sub_i)=ots(j,i)
            case(18)
              tqz(sub_i)=ots(j,i)
            case(19)
              q1(sub_i)=ots(j,i)
            case(20)
              q2(sub_i)=ots(j,i)
            case(21)
              q3(sub_i)=ots(j,i)
            case(22)
              q0(sub_i)=ots(j,i)
            case(23)
              q1(sub_i)=ots(j,i)
            case(24)
              q2(sub_i)=ots(j,i)
            case(25)
              q3(sub_i)=ots(j,i)
            end select
          enddo
        endif
      endif
    enddo
  else !if not idrank=0
    call irecv_world_iarr(atmbook,isend_nparticle(idrank),0, &
     120+idrank+1,irequest_recv(1))
    call irecv_world_farr(xxx,isend_nparticle(idrank),0, &
     120+idrank+2,irequest_recv(2))
    call irecv_world_farr(yyy,isend_nparticle(idrank),0, &
     120+idrank+3,irequest_recv(3))
    call irecv_world_farr(zzz,isend_nparticle(idrank),0, &
     120+idrank+4,irequest_recv(4))
    do j=1,nxyzlist_sub
      k=xyzlist_sub(j)
      select case(k)
      case(2)
        call irecv_world_farr(weight,isend_nparticle(idrank),0, &
         120+idrank+4+j,irequest_recv(4+j))
      case(3)
        call irecv_world_farr(rdimx,isend_nparticle(idrank),0, &
         120+idrank+4+j,irequest_recv(4+j))
      case(4)
        call irecv_world_farr(rdimy,isend_nparticle(idrank),0, &
         120+idrank+4+j,irequest_recv(4+j))
      case(5)
        call irecv_world_farr(rdimz,isend_nparticle(idrank),0, &
         120+idrank+4+j,irequest_recv(4+j))
      case(6)
        call irecv_world_farr(rdim,isend_nparticle(idrank),0, &
         120+idrank+4+j,irequest_recv(4+j))
      case(7)
        call irecv_world_farr(vxx,isend_nparticle(idrank),0, &
         120+idrank+4+j,irequest_recv(4+j))
      case(8)
        call irecv_world_farr(vyy,isend_nparticle(idrank),0, &
         120+idrank+4+j,irequest_recv(4+j))
      case(9)
        call irecv_world_farr(vzz,isend_nparticle(idrank),0, &
         120+idrank+4+j,irequest_recv(4+j))
      case(10)
        call irecv_world_farr(oxx,isend_nparticle(idrank),0, &
         120+idrank+4+j,irequest_recv(4+j))
      case(11)
        call irecv_world_farr(oyy,isend_nparticle(idrank),0, &
         120+idrank+4+j,irequest_recv(4+j))
      case(12)
        call irecv_world_farr(ozz,isend_nparticle(idrank),0, &
         120+idrank+4+j,irequest_recv(4+j))
      case(13)
        call irecv_world_farr(fxx,isend_nparticle(idrank),0, &
         120+idrank+4+j,irequest_recv(4+j))
      case(14)
        call irecv_world_farr(fyy,isend_nparticle(idrank),0, &
         120+idrank+4+j,irequest_recv(4+j))
      case(15)
        call irecv_world_farr(fzz,isend_nparticle(idrank),0, &
         120+idrank+4+j,irequest_recv(4+j))
      case(16)
        call irecv_world_farr(tqx,isend_nparticle(idrank),0, &
         120+idrank+4+j,irequest_recv(4+j))
      case(17)
        call irecv_world_farr(tqy,isend_nparticle(idrank),0, &
         120+idrank+4+j,irequest_recv(4+j))
      case(18)
        call irecv_world_farr(tqz,isend_nparticle(idrank),0, &
         120+idrank+4+j,irequest_recv(4+j))
      case(19)
        call irecv_world_farr(q1,isend_nparticle(idrank),0, &
         120+idrank+4+j,irequest_recv(4+j))
      case(20)
        call irecv_world_farr(q2,isend_nparticle(idrank),0, &
         120+idrank+4+j,irequest_recv(4+j))
      case(21)
        call irecv_world_farr(q3,isend_nparticle(idrank),0, &
         120+idrank+4+j,irequest_recv(4+j))
      case(22)
        call irecv_world_farr(q0,isend_nparticle(idrank),0, &
         120+idrank+4+j,irequest_recv(4+j))
      case(23)
        call irecv_world_farr(q1,isend_nparticle(idrank),0, &
         120+idrank+4+j,irequest_recv(4+j))
      case(24)
        call irecv_world_farr(q2,isend_nparticle(idrank),0, &
         120+idrank+4+j,irequest_recv(4+j))
      case(25)
        call irecv_world_farr(q3,isend_nparticle(idrank),0, &
         120+idrank+4+j,irequest_recv(4+j))
      end select
    enddo
    call waitall_world(irequest_recv(1:isend_nvar),isend_nvar)
  endif
  
  call get_sync_world
  
! since the data have been sent, we deallocate
  if(allocated(xxs))deallocate(xxs)
  if(allocated(yys))deallocate(yys)
  if(allocated(zzs))deallocate(zzs)
  if(allocated(ots))deallocate(ots)
  
#if 0
  call print_all_particles(100,'mioprima',1)
#endif
  
  lqinput=.false.
  ltinput=.false.
  if(nxyzlist_sub>0)then
    do j=1,nxyzlist_sub
      k=xyzlist_sub(j)
      if(k>=13 .and. k<=21)ltinput=.true.
      if(k>=22)lqinput=.true.
    enddo
  endif

  if((.not. ltinput) .and. (.not. lqinput))then
    call warning(19)
    !transform the rotational matrix in quaternions using an uniform
    !distribution
    do i=1,natms
      xsubm(1:3)=ZERO
      ysubm(1:3)=ZERO
      zsubm(1:3)=ZERO
      xsubm(1)=ONE
      ysubm(2)=ONE
      zsubm(3)=ONE
      call random_number(dmio(1))
      call random_number(dmio(2))
      call random_number(dmio(3))
      dmio(1:3)=dmio(1:3)*TWO*Pi
      do j=1,3
        call matrix_rotaxis(j,dmio(j),matrixmio)
        call rotate_vect(1,3,xsubm,ysubm,zsubm,(/ZERO,ZERO,ZERO/), &
         matrixmio)
      enddo
      rot(1)=xsubm(1)
      rot(2)=xsubm(2)
      rot(3)=xsubm(3)
      rot(4)=ysubm(1)
      rot(5)=ysubm(2)
      rot(6)=ysubm(3)
      rot(7)=zsubm(1)
      rot(8)=zsubm(2)
      rot(9)=zsubm(3)
      call rotmat_2_quat(rot,q0(i),q1(i),q2(i),q3(i))
#ifdef CHECKQUAT
      call quat_2_rotmat(q0(i),q1(i),q2(i),q3(i),newrot)
      do j=1,9
        if(abs(newrot(j)-rot(j))>toll)then
          ltestrot=.true.
        endif
      enddo
#endif
    enddo
  else
    if(ltinput)then
      !transform the rotational matrix in quaternions
      do i=1,natms
        rot(1)=fxx(i)
        rot(2)=fyy(i)
        rot(3)=fzz(i)
        rot(4)=tqx(i)
        rot(5)=tqy(i)
        rot(6)=tqz(i)
        rot(7)=q1(i)
        rot(8)=q2(i)
        rot(9)=q3(i)
        call rotmat_2_quat(rot,q0(i),q1(i),q2(i),q3(i))
        fxx(i)=ZERO
        fyy(i)=ZERO
        fzz(i)=ZERO
        tqx(i)=ZERO
        tqy(i)=ZERO
        tqz(i)=ZERO
      enddo
    endif
  endif
  
#ifdef CHECKQUAT
  call or_world_larr(ltestrot,1)
  if(ltestrot(1))then
    call error(22)
  endif
#endif

  if(lumass)then
    forall(i=1:natms)
      weight(i)=umass
    end forall
  endif
  
  return
  
 end subroutine initialize_map_particles
 
 subroutine driver_neighborhood_list(newlst,nstepsub)
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the neighborhood list construction
!
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
  
  logical, intent(in) :: newlst
  integer, intent(in) :: nstepsub
  
  if(lnolink)then
    call parlst(newlst)
  else
    call map_linkcell
  endif
  
 end subroutine driver_neighborhood_list
 
 subroutine vertest(newlst,tstep)
      
!***********************************************************************
!     
!     LBsoft subroutine to test for updating of Verlet list
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
  
  implicit none
  
  logical, intent(out) :: newlst
  
  logical, save :: newjob=.true.
  integer :: i,j,k,moved(1),ibuff
  integer :: fail=0
  real(kind=PRC) :: rmax,dr,tstep
      
  real(kind=PRC), allocatable, save :: xold(:),yold(:),zold(:)
  
  if((natms+mxrank-1)/mxrank.gt.msatms)call error(24)
      
  if(newjob)then
!   set up initial arrays 
    allocate (xold(msatms),yold(msatms),zold(msatms),stat=fail)
    if(fail.ne.0)call error(25)
    forall(i=1:natms)
      xold(i)=ZERO
      yold(i)=ZERO
      zold(i)=ZERO
    end forall
    newjob=.false.
    newlst=.true.
  else
    
!   integrate velocities 
    forall(i=1:natms)
      xold(i)=xold(i)+vxx(i)
      yold(i)=yold(i)+vyy(i)
      zold(i)=zold(i)+vzz(i)
    end forall
    
!   maximum displacement 
    rmax=(delr/TWO)**TWO
    
!   test atomic displacements
    moved(1)=0
    do i=1,natms
      dr=tstep**2*(xold(i)**2+yold(i)**2+zold(i)**2)
      if(dr.gt.rmax)moved(1)=moved(1)+1
    enddo
        
!   global sum of moved atoms
    call sum_world_iarr(moved,1)
        
!   test for new verlet list
    newlst=(moved(1).ge.2)
        
!     update stored positions
    if(newlst)then
      forall(i=1:natms)
        xold(i)=ZERO
        yold(i)=ZERO
        zold(i)=ZERO
      end forall
    endif
  endif
      
  return
  
 end subroutine vertest
 
 subroutine parlst(newlst)
  
!***********************************************************************
!     
!     LBsoft subroutine to compute Verlet Nlist
!
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
  
  implicit none

! Arguments
      
  logical,intent(in) :: newlst
  
! separation vectors and powers thereof
  
  real(kind=PRC) :: rsq,xm,ym,zm,rsqcut
  
! Loop counters
  integer :: iatm,ii,i
  integer :: k,jatm
  integer :: itype,ilentry

  logical :: lchk(1)
  integer :: ibig(1),idum

  if(newlst)then

    lchk=.false.
    ibig=0
    
    rsqcut = (rcut+delr)**TWO
    
    call allocate_array_bdf(natms)
    
    lentry(1:msatms)=0
    
    ii = 0
    do iatm = 1,natms
      ii=iatm
      k=0
      do jatm = 1,natms
        if(jatm<=iatm)cycle
        k=k+1
        xdf(k)=xxx(jatm)-xxx(iatm)
        ydf(k)=yyy(jatm)-yyy(iatm)
        zdf(k)=zzz(jatm)-zzz(iatm)
      enddo
      
      call pbc_images(imcon,k,cell,xdf,ydf,zdf)
      
      k=0
      ilentry=0
      do jatm = 1,natms
        if(jatm<=iatm)cycle
        k=k+1
        rsq=xdf(k)**TWO+ydf(k)**TWO+zdf(k)**TWO
        if(rsq<=rsqcut) then
          ilentry=ilentry+1
          if(ilentry.gt.mxlist)then
            lchk=.true.
            ibig(1)=max(ilentry,ibig(1))
            write(6,*)iatm,jatm,ilentry,rsq,rsqcut
            exit
          else
            list(ilentry,ii)=jatm
          endif
        endif
        lentry(ii)=ilentry
      enddo
      
    enddo
  
    
!   terminate job if neighbour list array exceeded
    call or_world_larr(lchk,1)
    if(lchk(1))then
      call max_world_iarr(ibig,1)
      call warning(23,real(ibig(1),kind=PRC))
      call error(21)
    endif
  endif

  return

 end subroutine parlst
 
 subroutine initialize_particle_force
 
!***********************************************************************
!     
!     LBsoft subroutine to initialize the force array
!
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************

  implicit none
  
  integer :: i
  
  forall(i=1:natms)
    fxx(i)=ZERO
    fyy(i)=ZERO
    fzz(i)=ZERO
  end forall
  
  return
  
 end subroutine initialize_particle_force
 
 subroutine driver_inter_force
 
!***********************************************************************
!     
!     LBsoft subroutine to drive the inter-particle force computation
!
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************

  implicit none
  real(kind=PRC), parameter :: tol=real(1.d-4,kind=PRC)
  integer :: i
  
  engcfg=ZERO
  
  call compute_inter_force(lentry,list)
  
  return
  
 end subroutine driver_inter_force
 
 subroutine compute_inter_force(lentrysub,listsub)
  
!***********************************************************************
!     
!     LBsoft subroutine to compute the inter-particle force
!
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
  
  implicit none

  integer, allocatable, dimension(:), intent(in) :: lentrysub
  integer, allocatable, dimension(:,:), intent(in) :: listsub
  
! separation vectors and powers thereof
  
  real(kind=PRC) :: rsq,xm,ym,zm,rsqcut,rrr,eps,sig,vvv,ggg,rmin,vmin, &
   kappa,rlimit,gmin
  
! Loop counters
  integer :: iatm,ii,i,ivdw
  integer :: k,jatm
  integer :: itype,ilentry
  
  real(kind=PRC), parameter :: s2rmin=TWO**(ONE/SIX)
  
  call allocate_array_bdf(natms)
  
  select case(ltpvdw(1))
  case (0)
    return
  case (1)
  
  ivdw=1
  eps=prmvdw(1,ivdw)
  sig=prmvdw(2,ivdw)
  rmin=s2rmin*sig
  rsqcut = (rmin)**TWO
  vmin=FOUR*eps*(sig/rmin)**SIX*((sig/rmin)**SIX-ONE)
  
  ii = 0
  do iatm = 1,natms
    ii=iatm
    do k = 1,lentrysub(ii)
      jatm=listsub(k,ii)
      xdf(k)=xxx(jatm)-xxx(iatm)
      ydf(k)=yyy(jatm)-yyy(iatm)
      zdf(k)=zzz(jatm)-zzz(iatm)
    enddo
      
    call pbc_images(imcon,lentrysub(ii),cell,xdf,ydf,zdf)
    
    do k = 1,lentrysub(ii)
      jatm=listsub(k,ii)
      rsq=xdf(k)**TWO+ydf(k)**TWO+zdf(k)**TWO
      if(rsq<=rsqcut) then
        rrr=sqrt(rsq)
        vvv=FOUR*eps*(sig/rrr)**SIX*((sig/rrr)**SIX-ONE)+vmin
        ggg=TWENTYFOUR*eps/rrr*(sig/rrr)**SIX*(TWO*(sig/rrr)**SIX-ONE)
        engcfg=engcfg+vvv
        fxx(iatm)=fxx(iatm)+ggg*xdf(k)/rrr
        fxx(jatm)=fxx(jatm)-ggg*xdf(k)/rrr
        fyy(iatm)=fyy(iatm)+ggg*ydf(k)/rrr
        fyy(jatm)=fyy(jatm)-ggg*ydf(k)/rrr
        fzz(iatm)=fzz(iatm)+ggg*zdf(k)/rrr
        fzz(jatm)=fzz(jatm)-ggg*zdf(k)/rrr
      endif
    enddo
  enddo
  
  case (2)
  
    ivdw=1
    eps=prmvdw(1,ivdw)
    sig=prmvdw(2,ivdw)
    
    rsqcut = (rcut)**TWO
    
    ii = 0
    do iatm = 1,natms
      ii=iatm
      do k = 1,lentrysub(ii)
        jatm=listsub(k,ii)
        xdf(k)=xxx(jatm)-xxx(iatm)
        ydf(k)=yyy(jatm)-yyy(iatm)
        zdf(k)=zzz(jatm)-zzz(iatm)
      enddo
        
      call pbc_images(imcon,lentrysub(ii),cell,xdf,ydf,zdf)
      
      do k = 1,lentrysub(ii)
        jatm=listsub(k,ii)
        rsq=xdf(k)**TWO+ydf(k)**TWO+zdf(k)**TWO
        if(rsq<=rsqcut) then
          rrr=sqrt(rsq)
          vvv=FOUR*eps*(sig/rrr)**SIX*((sig/rrr)**SIX-ONE)
          ggg=TWENTYFOUR*eps/rrr*(sig/rrr)**SIX*(TWO*(sig/rrr)**SIX-ONE)
          engcfg=engcfg+vvv
          fxx(iatm)=fxx(iatm)+ggg*xdf(k)/rrr
          fxx(jatm)=fxx(jatm)-ggg*xdf(k)/rrr
          fyy(iatm)=fyy(iatm)+ggg*ydf(k)/rrr
          fyy(jatm)=fyy(jatm)-ggg*ydf(k)/rrr
          fzz(iatm)=fzz(iatm)+ggg*zdf(k)/rrr
          fzz(jatm)=fzz(jatm)-ggg*zdf(k)/rrr
        endif
      enddo
    enddo
    
  case (3)
  
    ivdw=1
    kappa=prmvdw(1,ivdw)
    rmin=prmvdw(2,ivdw)
    rlimit=rmin-NINE/TEN
    rsqcut = (rmin)**TWO
    vmin=kappa*(rmin-rlimit)**(FIVE*HALF)
    gmin=FIVE*HALF*kappa*(rmin-rlimit)**(THREE*HALF)
    
  ii = 0
  do iatm = 1,natms
    ii=iatm
    do k = 1,lentrysub(ii)
      jatm=listsub(k,ii)
      xdf(k)=xxx(jatm)-xxx(iatm)
      ydf(k)=yyy(jatm)-yyy(iatm)
      zdf(k)=zzz(jatm)-zzz(iatm)
    enddo
      
    call pbc_images(imcon,lentrysub(ii),cell,xdf,ydf,zdf)
    
    do k = 1,lentrysub(ii)
      jatm=listsub(k,ii)
      rsq=xdf(k)**TWO+ydf(k)**TWO+zdf(k)**TWO
      if(rsq<=rsqcut) then
        rrr=sqrt(rsq)
        if(rrr<=rlimit)then
          vvv=vmin
          ggg=gmin
        else
          vvv=kappa*(rmin-rrr)**(FIVE*HALF)
          ggg=FIVE*HALF*kappa*(rmin-rrr)**(THREE*HALF)
        endif
        engcfg=engcfg+vvv
        fxx(iatm)=fxx(iatm)+ggg*xdf(k)/rrr
        fxx(jatm)=fxx(jatm)-ggg*xdf(k)/rrr
        fyy(iatm)=fyy(iatm)+ggg*ydf(k)/rrr
        fyy(jatm)=fyy(jatm)-ggg*ydf(k)/rrr
        fzz(iatm)=fzz(iatm)+ggg*zdf(k)/rrr
        fzz(jatm)=fzz(jatm)-ggg*zdf(k)/rrr
      endif
    enddo
  enddo
    
  case default
    call error(27)
  end select
  
  return

 end subroutine compute_inter_force
 
 subroutine initialize_lf()

!***********************************************************************
!     
!     LBsoft subroutine to initialize the velocity at half timestep
!     back in order to apply the verlet leapfrog integrator
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************

  implicit none
  
  integer :: i
  
! report the atoms velocity to half timestep back 
  forall(i=1:natms)      
    vxx(i)=vxx(i)-HALF*tstepatm/weight(i)*fxx(i)
    vyy(i)=vyy(i)-HALF*tstepatm/weight(i)*fyy(i)
    vzz(i)=vzz(i)-HALF*tstepatm/weight(i)*fzz(i)
  end forall
  
  return
      
 end subroutine initialize_lf
 
 subroutine nve_lf()

!***********************************************************************
!     
!     LBsoft subroutine for integrating newtonian EOM by Verlet leapfrog
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************

  implicit none

  integer, parameter :: nfailmax=10
  integer, dimension(nfailmax) :: fail
  integer :: i,k
  logical :: ltest(1)
  
  real(kind=PRC), allocatable :: uxx(:),uyy(:),uzz(:)    
  real(kind=PRC), allocatable :: xxo(:),yyo(:),zzo(:)
  real(kind=PRC), allocatable :: vxo(:),vyo(:),vzo(:)
      
! allocate working arrays
  fail(1:nfailmax)=0
  
  allocate(uxx(mxatms),stat=fail(1))
  allocate(uyy(mxatms),stat=fail(2))
  allocate(uzz(mxatms),stat=fail(3))
  allocate(xxo(mxatms),stat=fail(4))
  allocate(yyo(mxatms),stat=fail(5))
  allocate(zzo(mxatms),stat=fail(6))
  allocate(vxo(mxatms),stat=fail(7))
  allocate(vyo(mxatms),stat=fail(8))
  allocate(vzo(mxatms),stat=fail(9))
  
  ltest=.false.
  if(any(fail.ne.0))then
    do i=1,nfailmax
      if(fail(i).ne.0)exit
    enddo
    call warning(2,dble(i))
    ltest=.true.
  endif

  call or_world_larr(ltest,1)
  if(ltest(1))call error(28)

! store initial values of position and velocity    
  forall(i=1:natms)
    xxo(i)=xxx(i)
    yyo(i)=yyy(i)
    zzo(i)=zzz(i)
    vxo(i)=vxx(i)
    vyo(i)=vyy(i)
    vzo(i)=vzz(i)
  end forall
      
! move atoms by leapfrog algorithm    
  forall(i=1:natms)
!   update velocities       
    uxx(i)=vxx(i)+tstepatm/weight(i)*fxx(i)
    uyy(i)=vyy(i)+tstepatm/weight(i)*fyy(i)
    uzz(i)=vzz(i)+tstepatm/weight(i)*fzz(i)
!   update positions
    xxx(i)=xxo(i)+tstepatm*uxx(i)
    yyy(i)=yyo(i)+tstepatm*uyy(i)
    zzz(i)=zzo(i)+tstepatm*uzz(i)    
  end forall
  
#if 1
  do i=1,natms
    if(abs(uxx(i))>ONE .or. abs(uyy(i))>ONE .or. abs(uzz(i))>ONE )then
      write(6,*)'ERROR - numerical instability'
      call error(-1)
    endif
  enddo
#endif
  
! calculate full timestep velocity
  forall(i=1:natms)
    vxx(i)=HALF*(vxx(i)+uxx(i))
    vyy(i)=HALF*(vyy(i)+uyy(i))
    vzz(i)=HALF*(vzz(i)+uzz(i))
  end forall
      
! calculate kinetic energy    
  call getkin(engke)
    
! periodic boundary condition
  call pbc_images_centered(imcon,natms,cell,cx,cy,cz,xxx,yyy,zzz)
    
! updated velocity     
  forall(i=1:natms)
    vxx(i)=uxx(i)
    vyy(i)=uyy(i)
    vzz(i)=uzz(i)
  end forall
  

! deallocate work arrays
  deallocate (uxx,uyy,uzz,stat=fail(2))
  deallocate (xxo,yyo,zzo,stat=fail(3))
  deallocate (vxo,vyo,vzo,stat=fail(4))
      
  return
      
 end subroutine nve_lf
 
 subroutine getkin(engkes)

!***********************************************************************
!
!     LBsoft subroutine to calculate system kinetic energy
!
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC), intent(out) :: engkes
  
  integer :: i
  real(kind=PRC) :: engke(1)
  
  engke(1)=ZERO
  
  do i=1,natms
    engke(1)=engke(1)+weight(i)*(vxx(i)**TWO+vyy(i)**TWO+vzz(i)**TWO)
  enddo

  call sum_world_farr(engke,1)
  
  engkes=HALF*engke(1)

  return
  
 end subroutine getkin
 
 subroutine rotmat_2_quat(rot,q0s,q1s,q2s,q3s)
 
!***********************************************************************
!     
!     LBsoft subroutine for determining the quaternion associated
!     to a rotational matrix
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in), dimension(9) :: rot
  real(kind=PRC), intent(out) :: q0s,q1s,q2s,q3s
  
  real(kind=PRC) :: aq,bq,cq,dq,eq,fq,gq,hq,rnorm
  
  real(kind=PRC), parameter :: toll=real(1.d-4,kind=PRC)
 
 !determine quaternions from rotational matrix
  aq=rot(1)+rot(5)
  bq=rot(2)-rot(4)
  cq=rot(6)-rot(8)
  dq=rot(2)+rot(4)
  eq=rot(3)+rot(7)
  fq=rot(6)+rot(8)
  gq=rot(3)-rot(7)
  hq=rot(1)-rot(5)
          
  q0s=HALF*sqrt(aq+sqrt(aq*aq+bq*bq))
          
  if(q0s.gt.toll)then
    q1s=-FOURTH*cq/q0s
    q2s=FOURTH*gq/q0s
    q3s=-FOURTH*bq/q0s
  else
    q1s=HALF*sqrt(hq+sqrt(hq*hq+dq*dq))
    if(q1s.gt.toll)then
      q2s=FOURTH*dq/q1s
      q3s=FOURTH*eq/q1s
    else
      q2s=HALF*sqrt(-hq+sqrt(hq*hq+dq*dq))
      if(q2s.gt.toll)then
        q3s=FOURTH*fq/q2s
      else
        q3s=ONE
      endif
    endif
  endif
  
 !normalize quaternions
  rnorm=ONE/sqrt(q0s**TWO+q1s**TWO+q2s**TWO+q3s**TWO)
  q0s=rnorm*q0s
  q1s=rnorm*q1s
  q2s=rnorm*q2s
  q3s=rnorm*q3s
  
  return
  
 end subroutine rotmat_2_quat
 
 subroutine quat_2_rotmat(q0,q1,q2,q3,rot)
 
!***********************************************************************
!     
!     LBsoft subroutine for determining the rotational matrix from
!     quaternion using x convention for euler angles
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
      
  real(kind=PRC), intent(in) :: q0,q1,q2,q3
  real(kind=PRC), intent(out), dimension(9) :: rot
      
  rot(1)=q0**TWO+q1**TWO-q2**TWO-q3**TWO
  rot(2)=TWO*(q1*q2-q0*q3)
  rot(3)=TWO*(q1*q3+q0*q2)
  rot(4)=TWO*(q1*q2+q0*q3)
  rot(5)=q0**TWO-q1**TWO+q2**TWO-q3**TWO
  rot(6)=TWO*(q2*q3-q0*q1)
  rot(7)=TWO*(q1*q3-q0*q2)
  rot(8)=TWO*(q2*q3+q0*q1)
  rot(9)=q0**TWO-q1**TWO-q2**TWO+q3**TWO
  
  return
  
 end subroutine quat_2_rotmat
 
 subroutine rotate_vect(isub,lsub,xsubmy,ysubmy,zsubmy,center,matrix)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying a rotation of an array
!     depending on a rotation matrix
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: isub,lsub
  real(kind=PRC), dimension(*), intent(inout) :: xsubmy,ysubmy,zsubmy
  real(kind=PRC), dimension(3), intent(in) :: center
  real(kind=PRC), dimension(3,3), intent(in) :: matrix
  
  integer :: i
  real(kind=PRC), dimension(3) :: dvec1,dvec2
  
  do i=isub,lsub
    dvec1(1)=xsubmy(i)-center(1)
    dvec1(2)=ysubmy(i)-center(2)
    dvec1(3)=zsubmy(i)-center(3)
    call rotation(dvec1,matrix,dvec2)
    xsubmy(i)=dvec2(1)+center(1)
    ysubmy(i)=dvec2(2)+center(2)
    zsubmy(i)=dvec2(3)+center(3)
  enddo
  
  return
  
 end subroutine rotate_vect
 
 subroutine rotation(vec1,matrix,vec2)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying a rotation of a single vector
!     depending on a rotation matrix
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: vec1(3),matrix(3,3)
  real(kind=PRC), dimension(3), intent(out) :: vec2
  
  integer :: i,j
  
  vec2(:)=ZERO
  do i=1,3
    do j=1,3
      vec2(i)=matrix(i,j)*vec1(j)+vec2(i)
    enddo
  enddo
  
  return
  
 end subroutine rotation
 
 subroutine matrix_rotaxis(iaxis,theta,matrix)
 
!***********************************************************************
!     
!     LBsoft subroutine for building the rotation matrix
!     of a given angle and axis
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: iaxis
  real(kind=PRC), intent(in) :: theta
  real(kind=PRC), dimension(3,3), intent(out) :: matrix
  
  real(kind=PRC) :: ct,st
  
  ct=cos(theta)
  st=sin(theta)
  
  matrix(:,:)= ZERO
  
  select case(iaxis)
    case(1)
      matrix(1,1)= ONE
      matrix(2,2)=ct
      matrix(3,2)=-st
      matrix(2,3)=st
      matrix(3,3)=ct
    case(2)
      matrix(1,1)=ct
      matrix(3,1)=st
      matrix(2,2)= ONE
      matrix(1,3)=-st
      matrix(3,3)=ct
    case(3)
      matrix(1,1)=ct
      matrix(2,1)=-st
      matrix(1,2)=st
      matrix(2,2)=ct
      matrix(3,3)= ONE
    case default
      matrix(1:3,1:3)= ZERO
  end select
  
  return
  
 end subroutine matrix_rotaxis 
 
 subroutine map_linkcell
 
!***********************************************************************
!     
!     LBsoft subroutine for building the the linklist of the linkcells
!     contained in the subdomain of the process
!     note: if rcut is too big (or cellsub too small), the linkcell 
!     cannot be applied
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC) :: rcell(9),cellsub(9),det,rcsq,xdc,ydc,zdc
  logical :: loutbound(1)
  
  integer :: i,ilx,ily,ilz,icell,ix,iy,iz
  
  if(lnolink)return
  
  call allocate_array_bdf(natms)
        
!  real space cut off 
   rcsq=(rcut+delr)**TWO
  
   cellsub(1:9)=ZERO
   cellsub(1)=real(maxx-minx,kind=PRC)+ONE
   cellsub(5)=real(maxy-miny,kind=PRC)+ONE
   cellsub(9)=real(maxz-minz,kind=PRC)+ONE
  
   call invert(cellsub,rcell,det)
   
   ilx=int(cellsub(1)/(rcut+delr))
   ily=int(cellsub(5)/(rcut+delr))
   ilz=int(cellsub(9)/(rcut+delr))
!     check there are enough link cells
   lnolink=.false.
   if(ilx.lt.3)lnolink=.true.
   if(ily.lt.3)lnolink=.true.
   if(ilz.lt.3)lnolink=.true.
   
   ncells=ilx*ily*ilz
        
!  link-cell cutoff for reduced space
   xdc=real(ilx,kind=PRC)
   ydc=real(ily,kind=PRC)
   zdc=real(ilz,kind=PRC)
        
!  reduced space coordinates
   forall(i=1:natms)
     xdf(i)=rcell(1)*(xxx(i)-minx-HALF)
     ydf(i)=rcell(5)*(yyy(i)-miny-HALF)
     zdf(i)=rcell(9)*(zzz(i)-minz-HALF)
   end forall
   
!  fill up linklist
   nlinklist(1:ncellsmax)=0
   loutbound=.false.
   do i=1,natms
     ix=min(int(xdc*xdf(i)),ilx-1)
     iy=min(int(ydc*ydf(i)),ily-1)
     iz=min(int(zdc*zdf(i)),ilz-1)
     icell=1+ix+ilx*(iy+ily*iz)
     nlinklist(icell)=nlinklist(icell)+1
     if(nlinklist(icell)>mslistcell)then
       loutbound=.true.
       exit
     else
       linklist(icell,nlinklist(icell))=i
     endif
   enddo
   
   call or_world_larr(loutbound,1)
   if(loutbound(1))then
     call warning(22,real(mslistcell,kind=PRC))
     call warning(15,densvar)
     call error(23)
   endif
   
   return
   
  end subroutine map_linkcell
 
 subroutine pbc_images_centered(imcons,natmssub,cells,cx,cy,cz,xxs,yys,zzs)
      
!***********************************************************************
!     
!     LBsoft subroutine for calculating the minimum image
!     of particles within a specified cell
!     
!     dl_poly subroutine for calculating the minimum image
!     of atom pairs within a specified MD cell
!     
!     pbc conditions
!     0=F    1=T  periodic
!     itemp1      itemp2      itemp3       imcon
!          0           0           0           0
!          1           0           0           1
!          0           1           0           2
!          1           1           0           3
!          0           0           1           4
!          1           0           1           5
!          0           1           1           6
!          1           1           1           7
!     
!     note: in all cases the centre of the cell is specified by cx cy cz 
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: imcons,natmssub
  real(kind=PRC), intent(in) :: cells(9),cx,cy,cz
  real(kind=PRC), allocatable, dimension(:) :: xxs,yys,zzs

  integer :: i
  real(kind=PRC) aaa,bbb,ccc
  
  select case(imcons)
  case(0)
    return
  case(1)
    aaa=ONE/cells(1)
    forall(i=1:natmssub)
      xxs(i)=xxs(i)-cx
      xxs(i)=xxs(i)-cell(1)*nint(aaa*xxs(i))+cx
    end forall
  case(2)
    bbb=ONE/cells(5)
    forall(i=1:natmssub)
      yys(i)=yys(i)-cy
      yys(i)=yys(i)-cell(5)*nint(bbb*yys(i))+cy
    end forall
  case(3)
    aaa=ONE/cells(1)
    bbb=ONE/cells(5)
    forall(i=1:natmssub)
      xxs(i)=xxs(i)-cx
      yys(i)=yys(i)-cy
      xxs(i)=xxs(i)-cell(1)*nint(aaa*xxs(i))+cx
      yys(i)=yys(i)-cell(5)*nint(bbb*yys(i))+cy
    end forall
  case(4)
    ccc=ONE/cells(9)
    forall(i=1:natmssub)
      zzs(i)=zzs(i)-cz
      zzs(i)=zzs(i)-cell(9)*nint(ccc*zzs(i))+cz
    end forall
  case(5)
    aaa=ONE/cells(1)
    ccc=ONE/cells(9)
    forall(i=1:natmssub)
      xxs(i)=xxs(i)-cx
      zzs(i)=zzs(i)-cz
      xxs(i)=xxs(i)-cell(1)*nint(aaa*xxs(i))+cx
      zzs(i)=zzs(i)-cell(9)*nint(ccc*zzs(i))+cz
    end forall
  case(6)
    bbb=ONE/cells(5)
    ccc=ONE/cells(9)
    forall(i=1:natmssub)
      yys(i)=yys(i)-cy
      zzs(i)=zzs(i)-cz
      yys(i)=yys(i)-cell(5)*nint(bbb*yys(i))+cy
      zzs(i)=zzs(i)-cell(9)*nint(ccc*zzs(i))+cz
    end forall
  case(7)
    aaa=ONE/cells(1)
    bbb=ONE/cells(5)
    ccc=ONE/cells(9)
    forall(i=1:natmssub)
      xxs(i)=xxs(i)-cx
      yys(i)=yys(i)-cy
      zzs(i)=zzs(i)-cz
      xxs(i)=xxs(i)-cell(1)*nint(aaa*xxs(i))+cx
      yys(i)=yys(i)-cell(5)*nint(bbb*yys(i))+cy
      zzs(i)=zzs(i)-cell(9)*nint(ccc*zzs(i))+cz
    end forall
  end select
  
  return
      
 end subroutine pbc_images_centered
 
 subroutine pbc_images(imcons,natmssub,cells,xxs,yys,zzs)
      
!***********************************************************************
!     
!     LBsoft subroutine for calculating the minimum image
!     of particles within a specified cell
!     
!     dl_poly subroutine for calculating the minimum image
!     of atom pairs within a specified MD cell
!     
!     pbc conditions
!     0=F    1=T  periodic
!     itemp1      itemp2      itemp3       imcon
!          0           0           0           0
!          1           0           0           1
!          0           1           0           2
!          1           1           0           3
!          0           0           1           4
!          1           0           1           5
!          0           1           1           6
!          1           1           1           7
!     
!     note: in all cases the centre of the cell is at (0,0,0)
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: imcons,natmssub
  real(kind=PRC), intent(in) :: cells(9)
  real(kind=PRC), allocatable, dimension(:) :: xxs,yys,zzs

  integer :: i
  real(kind=PRC) aaa,bbb,ccc
  
  select case(imcons)
  case(0)
    return
  case(1)
    aaa=ONE/cells(1)
    forall(i=1:natmssub)
      xxs(i)=xxs(i)-cell(1)*nint(aaa*xxs(i))
    end forall
  case(2)
    bbb=ONE/cells(5)
    forall(i=1:natmssub)
      yys(i)=yys(i)-cell(5)*nint(bbb*yys(i))
    end forall
  case(3)
    aaa=ONE/cells(1)
    bbb=ONE/cells(5)
    forall(i=1:natmssub)
      xxs(i)=xxs(i)-cell(1)*nint(aaa*xxs(i))
      yys(i)=yys(i)-cell(5)*nint(bbb*yys(i))
    end forall
  case(4)
    ccc=ONE/cells(9)
    forall(i=1:natmssub)
      zzs(i)=zzs(i)-cell(9)*nint(ccc*zzs(i))
    end forall
  case(5)
    aaa=ONE/cells(1)
    ccc=ONE/cells(9)
    forall(i=1:natmssub)
      xxs(i)=xxs(i)-cell(1)*nint(aaa*xxs(i))
      zzs(i)=zzs(i)-cell(9)*nint(ccc*zzs(i))
    end forall
  case(6)
    bbb=ONE/cells(5)
    ccc=ONE/cells(9)
    forall(i=1:natmssub)
      yys(i)=yys(i)-cell(5)*nint(bbb*yys(i))
      zzs(i)=zzs(i)-cell(9)*nint(ccc*zzs(i))
    end forall
  case(7)
    aaa=ONE/cells(1)
    bbb=ONE/cells(5)
    ccc=ONE/cells(9)
    forall(i=1:natmssub)
      xxs(i)=xxs(i)-cell(1)*nint(aaa*xxs(i))
      yys(i)=yys(i)-cell(5)*nint(bbb*yys(i))
      zzs(i)=zzs(i)-cell(9)*nint(ccc*zzs(i))
    end forall
  end select
  
  return
      
 end subroutine pbc_images
 
 subroutine print_all_particles(iosub,filenam,itersub)
 
!***********************************************************************
!     
!     LBsoft subroutine for writing all the particles in ASCII format 
!     for diagnostic purposes always ordered also in parallel
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: iosub,itersub
  character(len=*), intent(in) :: filenam
  
  character(len=120) :: mynamefile
  integer :: i,j
  
  mynamefile=repeat(' ',120)
  mynamefile=trim(filenam)//write_fmtnumb(itersub)//'.xyz'
  
  if(idrank==0) then
    open(unit=iosub*idrank+23,file=trim(mynamefile),status='replace')
    write(iosub*idrank+23,'(i8)')natms_tot
    write(iosub*idrank+23,'(i8)')itersub
    close(iosub*idrank+23)
  endif
  
  do i=1,natms_tot
    do j=1,natms
      if(atmbook(j)==i)then
        open(unit=iosub*idrank+23,file=trim(mynamefile),status='old', &
         position='append')
        write(iosub*idrank+23,'(a8,6f16.8)')'C       ', &
         xxx(j),yyy(j),zzz(j),vxx(j),vyy(j),vzz(j)
        close(iosub*idrank+23)
      endif
    enddo
    call get_sync_world
  enddo
  
  return
  
 end subroutine print_all_particles
 
 end module particles_mod
