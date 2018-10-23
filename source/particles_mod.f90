
#include <default_macro.h>
 module particles_mod
 
!***********************************************************************
!     
!     LBsoft module containing variable and subroutines of the
!     particle managing
!
!     few variable names and subroutines are inspired from
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
                   irecv_world_iarr
 use error_mod
 use aop_mod
 use utility_mod, only : Pi,modulvec,cross,dot,gauss,ibuffservice, &
                   allocate_array_ibuffservice,nbuffservice,buffservice, &
                   allocate_array_buffservice,lbuffservice, &
                   allocate_array_lbuffservice, &
                   nbuffservice3d,buffservice3d, &
                   rand_noseeded,linit_seed,gauss_noseeded,write_fmtnumb

 use lbempi_mod,  only : commspop, commrpop, i4find, i4back, &
                   ownern,deallocate_ownern,commexch_dens, &
                   commwait_dens,comm_init_isfluid,commexch_vel_component, &
                   commwait_vel_component,ownernfind,ownernfind_arr, &
                   ownern,gminx,gmaxx,gminy,gmaxy,gminz,gmaxz
 
 use fluids_mod,  only : nx,ny,nz,nbuff

 
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
 integer, public, protected, save :: mslist=0
 
 !number of subdomains in parallel decomposition
 integer, public, protected, save :: nbig_cells=0
 
 !activate particle part
 logical, public, protected, save :: lparticles=.false.
 
 !cell vectors
 real(kind=PRC), dimension(9), public, protected, save :: cell
 
 !expected variance of the particle density given as ratio on natms_tot
 real(kind=PRC), public, protected, save :: densvar=ELEVEN/TEN
 real(kind=PRC), public, protected, save :: cx=ZERO
 real(kind=PRC), public, protected, save :: cy=ZERO
 real(kind=PRC), public, protected, save :: cz=ZERO
 
 !Potential cut-off
 real(kind=PRC), public, protected, save :: rcut=ZERO
 !Verlet neighbour list shell width
 real(kind=PRC), public, protected, save :: delr=ZERO
 
 !key for activate the body rotation
 logical, public, protected, save :: lrotate=.true.
 
 !key for activate the body rotation
 logical, public, protected, save :: lspherical=.true.
 
 !book of global particle ID
 integer, allocatable, public, protected, save :: atmbook(:)
 
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
 
 public :: allocate_particles
 public :: set_natms_tot
 public :: set_lparticles
 public :: set_densvar
 public :: set_rcut
 public :: set_delr
 public :: initialize_map_particles
 public :: rotmat_2_quat
 
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
 
 subroutine allocate_particles(ibctype)
 
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
  
  integer :: i
  integer, parameter :: nistatmax=100
  integer, dimension(nistatmax) :: istat
  
  logical :: ltest(1)
  
  if(.not. lparticles)return
  
  imcon=ibctype
  
  mxatms=ceiling(real(natms_tot,kind=PRC)*densvar)
  mslist=(mxatms/mxrank+1)*((mxatms-1)/mxrank+1)
  
  nbig_cells=mxrank
  
  cx=real(nx+1,kind=PRC)*HALF
  cy=real(ny+1,kind=PRC)*HALF
  cz=real(nz+1,kind=PRC)*HALF
  
  cell(1)=real(nx+1,kind=PRC)
  cell(2)=ZERO
  cell(3)=ZERO
  cell(4)=ZERO
  cell(5)=real(ny+1,kind=PRC)
  cell(6)=ZERO
  cell(7)=ZERO
  cell(8)=ZERO
  cell(9)=real(nz+1,kind=PRC)
  
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
  logical :: ltest(1)
  integer, dimension(0:mxrank-1) :: isend_nparticle
  integer :: isend_nvar
  integer, dimension(4+nxyzlist_sub,1:mxrank-1) :: irequest_send
  integer, dimension(4+nxyzlist_sub) :: irequest_recv
  real(kind=PRC), dimension(9) :: rot
  real(kind=PRC) :: matrixmio(3,3),dmio(3)
  real(kind=PRC), dimension(3) :: xsubm,ysubm,zsubm
  
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
  
  ltest=.false.
  if(nxyzlist_sub>0)then
    do j=1,nxyzlist_sub
      k=xyzlist_sub(j)
      if(k>=13 .and. k<=21)ltest=.true.
    enddo
  endif

  if(.not. ltest(1))then
    call warning(19)
    !transform the rotational matrix in quaternions using an uniform
    !distribution
    do i=1,natms
      xsubm(1:3)=ZERO
      ysubm(1:3)=ZERO
      zsubm(1:3)=ZERO
      xsubm(1)=ONE
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
    enddo
  else
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
  
  if(idrank==0)write(6,'(a)')'ATTENTION: MD is under developing!'
  call error(-1)
  
  return
  
 end subroutine initialize_map_particles
 
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
 
! subroutine parlink &
!       (newlst,natms,nneut,idnode,mxnode,imcon,rcut,delr)
      
!!***********************************************************************
!!     
!!     dl_poly subroutine for constructing the verlet neighbour
!!     list based on link-cell method with neutral groups
!!     frozen atoms taken into account
!!     
!!     to be used with the link version of exclude :excludeneu_link
!!     parallel replicated data version
!!     
!!     copyright - daresbury laboratory 1996
!!     author    - t. forester january 1996.
!!     
!!***********************************************************************
      
!      implicit none
      
!      logical lchk,newlst,linc,lfrzi,ldo,swop
!      integer natms,nneut,idnode,mxnode,imcon,idum,ibig(1)
!      integer i,nsbcll,ilx,ily,ilz,ncells,ix,iy,iz
!      integer icell,j,ic,ii,kc,jx,jy,jz,jc,ineu,ik,jneu,ineua,jneua
!      integer ika,jneua1,i1,j1
!      real(8) rcut,delr,rcsq,xm,ym,zm,det,xdc,ydc,zdc,tx,ty,tz
!      real(8) cx,cy,cz,sxd,syd,szd,rsq,xd,yd,zd
      
!  real(kind=PRC), allocatable, dimension(:) :: uxx,uyy,uzz
      
      
!  integer, parameter, dimension(27) :: &
!   nix=(/ 0, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1, 1,-1/)
!  integer, parameter, dimension(27) :: &
!   niy=(/ 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1, 1,-1,-1, 1, 1,-1/)
!  integer, parameter, dimension(27) :: &
!   niz=(/ 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1,-1, 1/)
  
!  integer, save :: fail=0
!  logical, save :: newjob=.true.
      
!!     allocate work arrays
      
!  allocate (uxx(mxatms),uyy(mxatms),uzz(mxatms),stat=fail)
!  if(fail.ne.0)call error(1900)
      
!  lchk=.true.
!  ibig=0
!  if(newlst)then
   
!!   zero link arrays
!    forall(i=1:natms)
!      link(i)=0
!    end forall
        
!!     construct pair force neighbour list
        
!        do i=1,mslist
          
!          lentry(i)=0
          
!        enddo
        
!!     real space cut off 
        
!        rcsq=(rcut+delr)**2
!!     
!!     create mock cell vector for non-periodic system
        
!        if(imcon.eq.0.or.imcon.eq.6)then
          
!!     find maximum x,y,z postions
          
!          xm=0.d0
!          ym=0.d0
!          zm=0.d0
          
!          do i=1,natms
            
!            xm=max(xm,abs(xxx(i)))
!            ym=max(ym,abs(yyy(i)))
!            zm=max(zm,abs(zzz(i)))
            
!          enddo
          
          
!        endif
        
!        call dcell(cell,celprp)
!        call invert(cell,rcell,det)
        
!!     number of subcells
        
!        nsbcll=14
        
!        ilx=int(celprp(7)/(rcut+delr))
!        ily=int(celprp(8)/(rcut+delr))
!        ilz=int(celprp(9)/(rcut+delr))
!!     
!!     check there are enough link cells
        
!        linc=.false.
!        if(ilx.lt.3)linc=.true.
!        if(ily.lt.3)linc=.true.
!        if(ilz.lt.3)linc=.true.
!        if(linc) call error(305)
        
!        ncells=ilx*ily*ilz
!        if(ncells.gt.mxcell)then
          
!          if(idnode.eq.0) write(nrite,*) 'mxcell must be >= ',ncells
!          call  error(392)
          
!        endif
        
!!     calculate link cell indices
        
!        do i=1,ncells
          
!          lct(i)=0
          
!        enddo
        
!!     link-cell cutoff for reduced space
        
!        xdc=dble(ilx)
!        ydc=dble(ily)
!        zdc=dble(ilz)
        
!!     reduced space coordinates
        
!        if(newjob)then
          
!          newjob=.false.
!          call images(imcon,idnode,mxnode,natms,cell,xxx,yyy,zzz)
          
!          if(mxnode.gt.1)  call merge
!     x      (idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)
          
!        endif
        
!        do i=1,natms
          
!          tx=xxx(i)
!          ty=yyy(i)
!          tz=zzz(i)
          
!          uxx(i)=(rcell(1)*tx+rcell(4)*ty+rcell(7)*tz)+0.5d0
!          uyy(i)=(rcell(2)*tx+rcell(5)*ty+rcell(8)*tz)+0.5d0
!          uzz(i)=(rcell(3)*tx+rcell(6)*ty+rcell(9)*tz)+0.5d0
          
!        enddo
        
!!     link neighbours 
        
!        do i=1,natms
          
!          ix=min(int(xdc*uxx(i)),ilx-1)
!          iy=min(int(ydc*uyy(i)),ily-1)
!          iz=min(int(zdc*uzz(i)),ilz-1)
          
!          icell=1+ix+ilx*(iy+ily*iz)
          
!          j=lct(icell)
!          lct(icell)=i
!          link(i)=j
          
!        enddo
        
!!     set control variables for loop over subcells
        
!        ix=1
!        iy=1
!        iz=1
        
!!     primary loop over subcells
        
!        do ic=1,ncells
          
!          ii=lct(ic)
!          if(ii.gt.0)then
            
!!     secondary loop over subcells
            
!            do kc=1,nsbcll
              
!              i=ii
              
!              cx=0.d0
!              cy=0.d0
!              cz=0.d0
!              jx=ix+nix(kc)
!              jy=iy+niy(kc)
!              jz=iz+niz(kc)
              
!!     minimum image convention
              
!              if(jx.gt.ilx)then
                
!                jx=jx-ilx
!                cx=1.d0
                
!              elseif(jx.lt.1)then
                
!                jx=jx+ilx
!                cx=-1.d0
                
!              endif
              
!              if(jy.gt.ily)then
                
!                jy=jy-ily
!                cy=1.d0
                
!              elseif(jy.lt.1)then
                
!                jy=jy+ily
!                cy=-1.d0
                
!              endif
              
!              if(jz.gt.ilz)then
                
!                jz=jz-ilz
!                cz=1.d0
                
!              elseif(jz.lt.1)then
                
!                jz=jz+ilz
!                cz=-1.d0
                
!              endif
              
!!     index of neighbouring cell
              
!              jc=jx+ilx*((jy-1)+ily*(jz-1))
!              j=lct(jc)
              
!!     ignore if empty
              
!              if(j.gt.0)then
                
!                do while(i.ne.0)
                  
!!     test if site is of interest to this node
                  
!                  ineu=lstneu(i)
!                  ik=0
                  
!!     i's  group index for this processor
                  
!                  if(mod(ineu-1,mxnode).eq.idnode)ik=((ineu-1)/mxnode)+1
                  
!!     test if i is a frozen atom
                  
!                  lfrzi=(lstfrz(i).ne.0)
                  
!                  if(ic.eq.jc) j=link(i)
!                  if(j.gt.0)then
                    
!                    do while(j.ne.0)
                      
!                      jneu=lstneu(j)
                      
!!     swop tests for switching of group indices,
!!     ldo for 'doing' interaction
                      
!                      swop=.false.
!                      ldo=(ik.gt.0)
!                      jneua=jneu
!                      ineua=ineu
!                      ika=ik
                      
!!     keep only Brode-Ahlrichs pairs
                      
!                      if(jneua.ge.ineua)then
                        
!                        if(jneua-ineua.gt.nneut/2)then 
                          
!                          swop=(mod(jneu-1,mxnode).eq.idnode)
!                          if(swop)then 
!                            ldo=((nneut+ineua-jneua).le.(nneut-1)/2)
!                          else
!                            ldo=.false.
!                          endif
                          
!                        endif
                        
!                      elseif(nneut+jneua-ineua.gt.(nneut-1)/2)then
                        
!                        swop=(mod(jneu-1,mxnode).eq.idnode)
!                        if(swop)then
!                          ldo=((ineua-jneua).le.nneut/2)
!                        else
!                          ldo=.false.
!                        endif
                        
!                      endif
                      
!                      if(swop.and.ldo)then
!                        jneua=ineu
!                        ineua=jneu
!                        ika=((jneu-1)/mxnode)+1
!                      endif
                      
!!     test of frozen atom pairs
                      
!                      if(lfrzi.and.ldo)ldo=(lstfrz(j).eq.0)
                      
!!     check we haven't already included this group in the list ...
                      
!                      jneua1=0
!                      do while(ldo.and.jneua1.lt.min(lentry(ika),mxlist))
                        
!                        jneua1=jneua1+1
!                        if(list(ika,jneua1).eq.jneua) ldo=.false.
                        
!                      enddo
                      
!                      if(ldo)then
                        
!!     distance in real space : minimum image applied
                        
!                        sxd=uxx(j)-uxx(i)+cx
!                        syd=uyy(j)-uyy(i)+cy
!                        szd=uzz(j)-uzz(i)+cz
                        
!                        xd=cell(1)*sxd+cell(4)*syd+cell(7)*szd
!                        yd=cell(2)*sxd+cell(5)*syd+cell(8)*szd
!                        zd=cell(3)*sxd+cell(6)*syd+cell(9)*szd
                        
!                        rsq=xd*xd+yd*yd+zd*zd
                        
!!     test of distance
                        
!                        if(rsq.lt.rcsq)then
                          
!                          lentry(ika)=lentry(ika)+1
!                          if(lentry(ika).gt.mxlist)then
                            
!                            ibig=max(ibig,lentry(ika))
!                            lchk=.false.
                            
!                          else
                            
!                            list(ika,lentry(ika))=jneua
                            
!                          endif
                          
!                        endif
                        
!                      endif
                      
!                      j=link(j)
                      
!                    enddo
                    
!                  endif
                  
!                  j=lct(jc)
!                  i=link(i)
                  
!                enddo
                
!              endif
              
!            enddo
            
!          endif
          
!          ix=ix+1
!          if(ix.gt.ilx)then
            
!            ix=1
!            iy=iy+1
            
!            if(iy.gt.ily)then
              
!              iy=1
!              iz=iz+1
              
!            endif
            
!          endif
          
!        enddo
        
!!     terminate job if neighbour list array exceeded
        
!        if(mxnode.gt.1) call gstate(lchk)
        
!        if(.not.lchk)then
          
!          call max_world_iarr(ibig,1)
!          if(idnode.eq.0)then
!            write(nrite,*)'mxlist must be at least ',ibig
!            write(nrite,*)'mxlist is currently ',mxlist
!          endif
!          call error(107)
          
!        endif
        
!!     sort list into order ..
!!     use link as a work array
        
!        ik=0
!        do i=1+idnode,nneut,mxnode
          
!          ik=ik+1
!          do j=1,lentry(ik)
            
!            link(j)=list(ik,j)
            
!          enddo
!          call shellsort(lentry(ik),link)
          
!!     ensure Brode-Ahlrichs ordering
          
!          i1=lentry(ik)+1
!          j1=0
!          do j=1,lentry(ik)
            
!            if(link(j).ge.i)then
              
!              j1=j1+1
!              list(ik,j1)=link(j)
!              i1=min(i1,j)
              
!            endif
            
!          enddo
          
!          do j=1,i1-1
            
!            j1=j1+1
!            list(ik,j1)=link(j)
            
!          enddo
          
!        enddo
        
!      endif
      
!!     deallocate work arrays 
      
!      deallocate (uxx,uyy,uzz,stat=fail)
      
!  return
  
! end subroutine parlinkneu
 
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
