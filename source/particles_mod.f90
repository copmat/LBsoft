
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
                   sum_world_farr, bcast_world_i, bcast_world_farr, &
                   or1_world_larr

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
 
 use fluids_mod,  only : nx,ny,nz,nbuff,minx,maxx,miny,maxy,minz,maxz, &
                   set_lbc_halfway,lbc_halfway,cssq,links,ex,ey,ez, &
                   init_particle_2_isfluid,push_comm_isfluid, &
                   ixpbc,iypbc,izpbc,isfluid,particle_bounce_back, &
                   initialize_new_isfluid,update_isfluid, &
                   driver_bc_isfluid,mapping_new_isfluid,&
                   particle_delete_fluids,particle_create_fluids, &
                   erase_fluids_in_particles,lunique_omega,omega, &
                   omega_to_viscosity,viscR,pimage,opp, &
                   compute_sc_particle_interact, driver_bc_pops, driver_bc_densities, &
                   driver_bc_pops_NOK

 
 implicit none
 
 private
 
 !key for periodic boundary conditions
 integer, public, protected, save :: imcon=0
 
 !key for integrator type
 integer, public, protected, save :: keyint=1
 
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
 
 !number of subdomains in parallel decomposition
 integer, public, protected, save :: nbig_cells=0
 
 !number of particle types with different features
 integer, public, protected, save :: ntype=0
 
 !max number of parameters for pair force fields
 integer, public, parameter :: mxpvdw=3
 
 !max number of pair force fields
 integer, public, parameter :: mxvdw=1
 
 !max number of different particle type
 integer, public, parameter :: mxntype=4
 
 !number of pair force fields defined in input
 integer, public, protected, save :: ntpvdw=0
 
 !number of max iteration in the quaternion update algorithm
 integer, public, protected, save :: mxquat=20
 
 !activate particle part
 logical, public, protected, save :: lparticles=.false.
 
 !activate reflecting side wall
 logical, public, protected, save :: lsidewall_md=.false.
 
 !activate lubrication force
 logical, public, protected, save :: llubrication=.false.
 
 !kbT factor
 real(kind=PRC), public, protected, save :: tempboltz=cssq
 
 !time step of the lagrangian integrator
 real(kind=PRC), public, protected, save :: tstepatm=ONE
 
 !total degree of freedom of the particles
 real(kind=PRC), public, protected, save :: degfre=ZERO
 
 !rotational degree of freedom of a single particle
 real(kind=PRC), public, protected, save :: degrot=ZERO
 
 !tollerance in the quaternion update algorithm
 real(kind=PRC), public, protected, save :: quattol=real(1.d-12,kind=PRC)
 
 !square tollerance in the quaternion update algorithm 
 real(kind=PRC), public, protected, save :: sqquattol
 
 !initiali particle temperature in KbT unit
 real(kind=PRC), public, protected, save :: init_temp=ONE
 
 !cell vectors
 real(kind=PRC), dimension(9), public, protected, save :: cell
 
 !expected variance of the particle density given as ratio on natms_tot
 real(kind=PRC), public, protected, save :: densvar=ELEVEN/TEN
 real(kind=PRC), public, protected, save :: cx=ZERO
 real(kind=PRC), public, protected, save :: cy=ZERO
 real(kind=PRC), public, protected, save :: cz=ZERO
 
 !kinetic energy
 real(kind=PRC), public, protected, save :: engke=ZERO
 
 !rotational kinetic energy
 real(kind=PRC), public, protected, save :: engrot=ZERO
 
 !configurational energy
 real(kind=PRC), public, protected, save :: engcfg=ZERO
 
 !total energy
 real(kind=PRC), public, protected, save :: engtot=ZERO
 
 !Potential cut-off
 real(kind=PRC), public, protected, save :: rcut=ZERO
 
 !Verlet neighbour list shell width
 real(kind=PRC), public, protected, save :: delr=ZERO
 
 !force constant of reflecting side wall (hertzian force field)
 real(kind=PRC), public, protected, save :: sidewall_md_k=ZERO
 
 !force distance from reflecting side wall (hertzian force field)
 real(kind=PRC), public, protected, save :: sidewall_md_rdist=ONE
 
 !mean particle force along directions x y z
 real(kind=PRC), public, protected, save :: meanpfx=ZERO
 real(kind=PRC), public, protected, save :: meanpfy=ZERO
 real(kind=PRC), public, protected, save :: meanpfz=ZERO
 
 !external particle force along directions x y z
 real(kind=PRC), save, protected, public :: ext_fxx=ZERO
 real(kind=PRC), save, protected, public :: ext_fyy=ZERO
 real(kind=PRC), save, protected, public :: ext_fzz=ZERO
 
 !external particle torque along directions x y z
 real(kind=PRC), save, protected, public :: ext_tqx=ZERO
 real(kind=PRC), save, protected, public :: ext_tqy=ZERO
 real(kind=PRC), save, protected, public :: ext_tqz=ZERO
 
 !key for set initial particle temperature
 logical, public, protected, save :: linit_temp=.false.
 
 !key for activate the velocity Verlet
 logical, public, protected, save :: lvv=.false.
 
 !key for activate the body rotation
 logical, public, protected, save :: lrotate=.false.
 
 !there are enough cells for the link approach within the same process?
 logical, public, protected, save :: lnolink=.false.
 
 !book of global particle ID
 integer, allocatable, public, protected, save :: atmbook(:)
 

 !number of link cell within the same process
 integer, save :: ncells=0
 
 !global minimum number of link cell within the same process
 integer, save :: ncellsmin=0
 
 !number of reflecting nodes in the spherical particle
 integer, save, protected, public :: nsphere
 
 !coordinate list of reflecting nodes in the spherical particle
 integer, allocatable, dimension(:,:), save, protected, public :: spherelist
 
 !radial distance list of reflecting nodes in the spherical particle
 real(kind=PRC), allocatable, dimension(:), save :: spheredist
 
 !unit vector list of reflecting nodes in the spherical particle
 real(kind=PRC), allocatable, dimension(:,:), save :: sphereuvec
 
 !number of dead nodes in the spherical particle
 integer, save :: nspheredead
 
 !coordinate list of dead nodes in the spherical particle
 integer, allocatable, dimension(:,:), save :: spherelistdead
 
 !verlet list
 integer, allocatable, save :: lentry(:)
 integer, allocatable, save :: list(:,:)
 
 !store the type of pair force field
 integer, allocatable, public, protected, save :: ltpvdw(:)
 
 !mask to associated the type of pair force field to two particle types
 integer, allocatable, public, protected, save :: mskvdw(:,:)
 
 !flag for particle moving grid overlap
 logical(kind=1), allocatable, public, protected, save :: lmove(:)
 
 !flag for particle leaving the subdomain (for parallel version)
 logical(kind=1), allocatable, public, protected, save :: lmove_dom(:)
 
 !string denoting the particle type
 character(len=8), allocatable, public, protected, save :: atmnamtype(:)
 
 !key for particle shape
 integer, allocatable, public, protected, save :: ishape(:)
 
 !set unique mass
 logical, allocatable, public, protected, save :: lumass(:)
 
 !set unique radius
 logical, allocatable, public, protected, save :: lurdim(:)
 
 !type of particles
 integer, allocatable, public, protected, save :: ltype(:)
 
 !position of particles
 real(kind=PRC), allocatable, public, protected, save :: xxx(:)
 real(kind=PRC), allocatable, public, protected, save :: yyy(:)
 real(kind=PRC), allocatable, public, protected, save :: zzz(:)
 
 !old position of particles
 real(kind=PRC), allocatable, public, protected, save :: xxo(:)
 real(kind=PRC), allocatable, public, protected, save :: yyo(:)
 real(kind=PRC), allocatable, public, protected, save :: zzo(:)
 
 !components of linear velocity
 real(kind=PRC), allocatable, public, protected, save :: vxx(:)
 real(kind=PRC), allocatable, public, protected, save :: vyy(:)
 real(kind=PRC), allocatable, public, protected, save :: vzz(:)
 
 !old components of linear velocity
 real(kind=PRC), allocatable, public, protected, save :: vxo(:)
 real(kind=PRC), allocatable, public, protected, save :: vyo(:)
 real(kind=PRC), allocatable, public, protected, save :: vzo(:)
 
 !components of force
 real(kind=PRC), allocatable, public, protected, save :: fxx(:)
 real(kind=PRC), allocatable, public, protected, save :: fyy(:)
 real(kind=PRC), allocatable, public, protected, save :: fzz(:)
 
 !components of linear force due to bounce back
 real(kind=PRC), allocatable, public, protected, save :: fxb(:)
 real(kind=PRC), allocatable, public, protected, save :: fyb(:)
 real(kind=PRC), allocatable, public, protected, save :: fzb(:)
 
 !old components of linear force due to bounce back
 real(kind=PRC), allocatable, public, protected, save :: fxbo(:)
 real(kind=PRC), allocatable, public, protected, save :: fybo(:)
 real(kind=PRC), allocatable, public, protected, save :: fzbo(:)
 
 !components of torque due to bounce back
 real(kind=PRC), allocatable, public, protected, save :: txb(:)
 real(kind=PRC), allocatable, public, protected, save :: tyb(:)
 real(kind=PRC), allocatable, public, protected, save :: tzb(:)
 
 !old components of torque due to bounce back
 real(kind=PRC), allocatable, public, protected, save :: txbo(:)
 real(kind=PRC), allocatable, public, protected, save :: tybo(:)
 real(kind=PRC), allocatable, public, protected, save :: tzbo(:)
 
 !radius of particles
 real(kind=PRC), allocatable, public, protected, save :: rdimx(:)
 real(kind=PRC), allocatable, public, protected, save :: rdimy(:)
 real(kind=PRC), allocatable, public, protected, save :: rdimz(:)
 !particle mass
 real(kind=PRC), allocatable, public, protected, save :: weight(:)
 !inverse particle mass
 real(kind=PRC), allocatable, public, protected, save :: rmass(:)
 !rotational inertia in body fixed frame
 real(kind=PRC), allocatable, public, protected, save :: rotin(:)
 !components of angular velocity (NOTE: in body referement system)
 real(kind=PRC), allocatable, public, protected, save :: oxx(:)
 real(kind=PRC), allocatable, public, protected, save :: oyy(:)
 real(kind=PRC), allocatable, public, protected, save :: ozz(:)
 !components of quaternion vector
 real(kind=PRC), allocatable, public, protected, save :: q0(:)
 real(kind=PRC), allocatable, public, protected, save :: q1(:)
 real(kind=PRC), allocatable, public, protected, save :: q2(:)
 real(kind=PRC), allocatable, public, protected, save :: q3(:)
 !components of torque on rigid body (NOTE: in world referement system)
 real(kind=PRC), allocatable, public, protected, save :: tqx(:)
 real(kind=PRC), allocatable, public, protected, save :: tqy(:)
 real(kind=PRC), allocatable, public, protected, save :: tqz(:)
 
 !rotational inertia about x, y, z
 real(kind=PRC), allocatable, public, protected, save :: rotinx(:)
 real(kind=PRC), allocatable, public, protected, save :: rotiny(:)
 real(kind=PRC), allocatable, public, protected, save :: rotinz(:)
 
 !parameters of all the pair force fields
 real(kind=PRC), allocatable, public, protected, save :: prmvdw(:,:)
 

 public :: allocate_particles
 public :: set_natms_tot
 public :: set_lparticles
 public :: set_densvar
 public :: set_rcut
 public :: set_delr
 public :: set_ishape
 public :: set_lrotate
 public :: set_lvv
 public :: set_lsidewall_md
 public :: set_sidewall_md_k
 public :: set_sidewall_md_rdist
 public :: initialize_map_particles
 public :: rotmat_2_quat
 public :: parlst
 public :: set_ntpvdw
 public :: allocate_field_array
 public :: set_field_array
 public :: vertest
 public :: set_umass
 public :: set_rdim
 public :: initialize_particle_force
 public :: initialize_particle_energy
 public :: driver_inter_force
 public :: initialize_integrator_lf
 public :: nve_lf
 public :: integrate_particles_vv
 public :: merge_particle_energies
 public :: set_init_temp
 public :: init_particles_fluid_interaction
 public :: apply_particle_bounce_back
 public :: store_old_pos_vel_part
 public :: inter_part_and_grid
 public :: force_particle_bounce_back
 public :: merge_particle_force
 public :: build_new_isfluid
 public :: set_atmnamtype
 public :: set_ntype
 public :: find_type
 public :: allocate_particle_features
 public :: take_rotversor
 public :: take_rotversorx
 public :: take_rotversory
 public :: take_rotversorz
 public :: clean_fluid_inside_particle
 public :: compute_mean_particle_force
 public :: set_llubrication
 public :: set_value_ext_force_particles
 public :: set_value_ext_torque_particles
 public :: q2eul
 public :: compute_psi_sc_particles, restore_particles
 
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
 
 subroutine set_ntype(itemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the ntype of particles 
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: itemp
  
  ntype=itemp
  
  return
  
 end subroutine set_ntype
 
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
 
 subroutine set_ishape(itemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the ishape protected variable
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in), dimension(1:mxntype) :: itemp
  
  if(.not. lparticles)return
  
  ishape(1:mxntype)=itemp(1:mxntype)
  
  where(ishape(1:mxntype)==0)
    lurdim(1:mxntype)=.true.
  end where
  
  return
  
 end subroutine set_ishape
 
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
 
 subroutine set_lvv(ltemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the lvv protected variable
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
  
  logical, intent(in) :: ltemp
  
  lvv=ltemp
  
  return
  
 end subroutine set_lvv
 
 subroutine set_umass(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the umass protected variable
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in), dimension(1:mxntype) :: dtemp
  
  if(.not. lparticles)return
  
  lumass(1:ntype)=.true.
  weight(1:mxntype)=dtemp(1:mxntype)
  
  return
  
 end subroutine set_umass
 
 subroutine set_rdim(dtempx,dtempy,dtempz)
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the urdim protected variable
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in), dimension(1:mxntype) :: dtempx,dtempy,dtempz
  
  if(.not. lparticles)return
  
  rdimx(1:mxntype)=dtempx(1:mxntype)
  rdimy(1:mxntype)=dtempy(1:mxntype)
  rdimz(1:mxntype)=dtempz(1:mxntype)
  
  return
  
 end subroutine set_rdim
 
 subroutine set_init_temp(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the ntpvdw protected variable
!     init_temp denotes the initiali particle temperature in KbT unit
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: dtemp
  
  linit_temp=.true.
  init_temp=dtemp
  
  return
  
 end subroutine set_init_temp
 
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
 
 subroutine set_atmnamtype(ctemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for allocate the atmnamtype protected variable
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
 
  implicit none
  
  character(len=8), intent(in), dimension(mxntype) :: ctemp
  
  integer :: i
  
  if(.not. lparticles)return
  
  do i=1,mxntype
    atmnamtype(i)=ctemp(i)
  enddo
  
  return
  
 end subroutine set_atmnamtype
 
 subroutine set_lsidewall_md(ltemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the lsidewall_md protected variable
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
 
  implicit none
  
  logical, intent(in) :: ltemp
  
  lsidewall_md=ltemp
  
  return
  
 end subroutine set_lsidewall_md
 
 subroutine set_sidewall_md_k(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the sidewall_md_k protected variable
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification December 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: dtemp
  
  sidewall_md_k=dtemp
  
  return
  
 end subroutine set_sidewall_md_k
 
 subroutine set_sidewall_md_rdist(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the sidewall_md_rdist 
!     protected variable
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification December 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: dtemp
  
  sidewall_md_rdist=dtemp
  
  return
  
 end subroutine set_sidewall_md_rdist
 
 subroutine set_llubrication(ltemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the llubrication protected variable
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification December 2018
!     
!***********************************************************************
 
  implicit none
  
  logical, intent(in) :: ltemp
  
  llubrication=ltemp
  
  return
  
 end subroutine set_llubrication
 
 subroutine set_value_ext_force_particles(dtemp1,dtemp2,dtemp3)
 
!***********************************************************************
!     
!     LBsoft subroutine for set the external force values acting
!     on the particles 
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification December 2018
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC), intent(in) :: dtemp1,dtemp2,dtemp3
  
  ext_fxx = dtemp1
  ext_fyy = dtemp2
  ext_fzz = dtemp3
  
  return
  
 end subroutine set_value_ext_force_particles
 
 subroutine set_value_ext_torque_particles(dtemp1,dtemp2,dtemp3)
 
!***********************************************************************
!     
!     LBsoft subroutine for set the external torque values acting
!     on the particles 
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification December 2018
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC), intent(in) :: dtemp1,dtemp2,dtemp3
  
  ext_tqx = dtemp1
  ext_tqy = dtemp2
  ext_tqz = dtemp3
  
  return
  
 end subroutine set_value_ext_torque_particles
 
 subroutine compute_mean_particle_force()
 
!***********************************************************************
!     
!     LBsoft subroutine for computing the mean particle force
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
  
  implicit none
  
  integer :: i, myi
  real(kind=PRC) :: dtemp(3)
  
  dtemp(1:3)=ZERO
  
  do myi=1,natms
    i = atmbook(myi)
    dtemp(1)=dtemp(1)+fxx(i)
    dtemp(2)=dtemp(2)+fyy(i)
    dtemp(3)=dtemp(3)+fzz(i)
  enddo
  
  if(mxrank>1)then
    call sum_world_farr(dtemp,3)
  endif
  
  meanpfx=dtemp(1)/real(natms_tot,kind=PRC)
  meanpfy=dtemp(2)/real(natms_tot,kind=PRC)
  meanpfz=dtemp(3)/real(natms_tot,kind=PRC)
 end subroutine compute_mean_particle_force
 

 subroutine allocate_particle_features()
 
!***********************************************************************
!     
!     LBsoft subroutine for allocate the atmnamtype protected variable
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
 
  implicit none
  
  integer :: i
  integer, parameter :: nistatmax=100
  integer, dimension(nistatmax) :: istat
  logical :: ltest(1)
  
  istat(1:nistatmax)=0
  
  
  allocate(rdimx(mxntype),stat=istat(1))
  allocate(rdimy(mxntype),stat=istat(2))
  allocate(rdimz(mxntype),stat=istat(3))
  
  allocate(weight(mxntype),stat=istat(4))
  
  allocate(rotin(mxntype),stat=istat(5))
  
  allocate(rmass(mxntype),stat=istat(6))
  allocate(atmnamtype(mxntype),stat=istat(7))
  
  allocate(ishape(mxntype),stat=istat(8))
  
  allocate(lumass(mxntype),stat=istat(9))
  allocate(lurdim(mxntype),stat=istat(10))
  
  if(lrotate)then
    allocate(rotinx(mxntype),stat=istat(11))
    allocate(rotiny(mxntype),stat=istat(12))
    allocate(rotinz(mxntype),stat=istat(13))
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
  if(ltest(1))call error(32)
 end subroutine allocate_particle_features

 
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
  
  integer :: fail(3)
  
  allocate(prmvdw(mxpvdw,ntpvdw),stat=fail(1))
  allocate(ltpvdw(ntpvdw),stat=fail(2))
  allocate(mskvdw(mxntype,mxntype),stat=fail(3))
  
  call sum_world_iarr(fail,3)
  if(any(fail.ne.0))call error(26)
  
  prmvdw(1:mxpvdw,1:ntpvdw)=ZERO
  ltpvdw(1:ntpvdw)=0
  mskvdw(1:mxntype,1:mxntype)=0
  
  return
  
 end subroutine allocate_field_array
 
 subroutine set_field_array(iarr1,iarr2,itemp1,dtemp2,itemp2)
 
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
  integer, dimension(mxntype,mxntype), intent(in) :: itemp2
  
  integer :: i
  
  if(.not. lparticles)return
  
  do i=1,iarr1
    ltpvdw(i)=itemp1(i)
    prmvdw(1:iarr2,i)=dtemp2(1:iarr2,i)
  enddo
  
  mskvdw(1:mxntype,1:mxntype)=itemp2(1:mxntype,1:mxntype)
  
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
  
  
  tstepatm=tstep
  
  imcon=ibctype
  
  mxatms=ceiling(real(natms_tot,kind=PRC)*densvar)
  mxatms = natms_tot    !! In this version, all processes have all atoms from begin to end!


  volm=real(nx+1,kind=PRC)*real(nx+1,kind=PRC)*real(nx+1,kind=PRC)
  dens=dble(mxatms)/volm
  ratio=(ONE+HALF)*dens*(FOUR*pi/THREE)*(rcut+delr)**THREE
  mxlist=min(nint(ratio),(mxatms+1)/2)
  mxlist=max(mxlist,1000)
  
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
  
  allocate(lmove_dom(mxatms),stat=istat(98))
  allocate(lmove(mxatms),stat=istat(99))
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
  
  if(lrotate)then
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
  
  endif
  
  allocate(ltype(mxatms),stat=istat(23))
  
  allocate(xxo(mxatms),stat=istat(26))
  allocate(yyo(mxatms),stat=istat(27))
  allocate(zzo(mxatms),stat=istat(28))
  
  allocate(vxo(mxatms),stat=istat(29))
  allocate(vyo(mxatms),stat=istat(30))
  allocate(vzo(mxatms),stat=istat(31))
  
  allocate(fxb(mxatms),stat=istat(32))
  allocate(fyb(mxatms),stat=istat(33))
  allocate(fzb(mxatms),stat=istat(34))
  
  allocate(fxbo(mxatms),stat=istat(35))
  allocate(fybo(mxatms),stat=istat(36))
  allocate(fzbo(mxatms),stat=istat(37))
  
  if(lrotate)then
    allocate(txb(mxatms),stat=istat(38))
    allocate(tyb(mxatms),stat=istat(39))
    allocate(tzb(mxatms),stat=istat(40))
    
    allocate(txbo(mxatms),stat=istat(41))
    allocate(tybo(mxatms),stat=istat(42))
    allocate(tzbo(mxatms),stat=istat(43))
  endif
  
  
!  if(mxrank>1)then
!    if(idrank==0)write(6,'(a)')'ATTENTION: MD in parallel is under developing!'
!    call error(-1)
!  endif
  
#if 1
  if(any(ishape(1:ntype)/=0))then
    if(idrank==0)write(6,'(a)')'ATTENTION: non spherical particle part is under developing!'
    call error(-1)
  endif
  
  !at the moment only verlet list
  if(.not. lnolink)then
    if(idrank==0)write(6,'(a)')'ATTENTION: at the moment link cell is not implemented!'
    if(idrank==0)write(6,'(a)')'ATTENTION: automatic switch to Verlet list mode!'
    lnolink=.true.
  endif
  
#endif
  
  allocate(lentry(mxatms),stat=istat(23))
  allocate(list(mxlist,mxatms),stat=istat(24))

  
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
  
  lmove(1:mxatms)=.false.
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
  
  if(lrotate)then
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
  endif
  
  ltype(1:mxatms)=0
  
  xxo(1:mxatms)=ZERO
  yyo(1:mxatms)=ZERO
  zzo(1:mxatms)=ZERO
  
  vxo(1:mxatms)=ZERO
  vyo(1:mxatms)=ZERO
  vzo(1:mxatms)=ZERO
  
  fxb(1:mxatms)=ZERO
  fyb(1:mxatms)=ZERO
  fzb(1:mxatms)=ZERO
  
  fxbo(1:mxatms)=ZERO
  fybo(1:mxatms)=ZERO
  fzbo(1:mxatms)=ZERO
  
  if(lrotate)then
    txb(1:mxatms)=ZERO
    tyb(1:mxatms)=ZERO
    tzb(1:mxatms)=ZERO
    
    txbo(1:mxatms)=ZERO
    tybo(1:mxatms)=ZERO
    tzbo(1:mxatms)=ZERO
  endif
  
  lentry(1:mxatms)=0
  list(1:mxlist,1:mxatms)=0
 end subroutine allocate_particles
 

 subroutine find_type(carr,nch,itype)
  
!***********************************************************************
!     
!     LBsoft subroutine for finding the atom type from its label
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: nch
  character(len=nch), intent(in) :: carr
  integer, intent(out) :: itype
  integer :: i
  
  character(len=8) :: myservice
  
  myservice(1:8)=trim(carr(1:nch))//repeat(' ',8-nch)
  
  itype=0
  do i=1,mxntype
    if(trim(atmnamtype(i))==trim(myservice))then
      itype=i
      return
    endif
  enddo
  
  return
  
 end subroutine find_type
 
 subroutine initialize_map_particles(nxyzlist_sub,xyzlist_sub, &
  ltypes,xxs,yys,zzs,ots,lvelocity)
  
!***********************************************************************
!     
!     LBsoft subroutine for mapping all the particles and send
!     their values to the appropiate process following the LB
!     domain decomposition
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
 
  implicit none
  integer, intent(in) :: nxyzlist_sub
  integer, intent(in), allocatable, dimension(:)  :: xyzlist_sub
  integer, allocatable, dimension(:) :: ltypes
  real(kind=PRC), allocatable, dimension(:) :: xxs,yys,zzs
  real(kind=PRC), allocatable, dimension(:,:) :: ots
  logical, intent(in) :: lvelocity
  
  integer :: i,j,k,l,iatm,ids,idrank_sub,itype, myi, num_ext
  logical :: ltest(1),lqinput
  integer, dimension(0:mxrank-1) :: isend_nparticle
  integer :: isend_nvar
  integer, dimension(4+nxyzlist_sub,1:mxrank-1) :: irequest_send
  integer, dimension(4+nxyzlist_sub) :: irequest_recv
  real(kind=PRC), dimension(9) :: rot,newrot
  real(kind=PRC) :: matrixmio(3,3),dmio(3),dmio2(3),rnorm
  real(kind=PRC), dimension(3) :: xsubm,ysubm,zsubm
  real(kind=PRC), dimension(0:3) :: qs
#ifdef CHECKQUAT
  real(kind=PRC), parameter :: toll=real(1.d-4,kind=PRC)
  logical :: ltestrot(1)=.false.
#endif
  
  
  isend_nvar=4+nxyzlist_sub
    
  if(nxyzlist_sub>0)then
    do j=1,nxyzlist_sub
      k=xyzlist_sub(j)
      select case(k)
      case(2)
        call warning(33,ZERO,'mass')
        lumass=.true.
      case(3)
        call warning(33,ZERO,'radx')
      case(4)
        call warning(33,ZERO,'rady')
      case(5)
        call warning(33,ZERO,'radz')
      case(6)
        call warning(33,ZERO,'rad')
      end select
    enddo
  endif
  
  
  if(idrank==0)then
      do i=1,natms_tot
          ltype(i)=ltypes(i)
          xxx(i)=xxs(i)
          yyy(i)=yys(i)
          zzz(i)=zzs(i)

          if(nxyzlist_sub>0)then
            do j=1,nxyzlist_sub
              k=xyzlist_sub(j)
              select case(k)
              case(7)
                vxx(i)=ots(j,i)
              case(8)
                vyy(i)=ots(j,i)
              case(9)
                vzz(i)=ots(j,i)
              case(10)
                oxx(i)=ots(j,i)
              case(11)
                oyy(i)=ots(j,i)
              case(12)
                ozz(i)=ots(j,i)
              case(13)
                q1(i)=ots(j,i)
              case(14)
                q2(i)=ots(j,i)
              case(15)
                q3(i)=ots(j,i)
              end select
            enddo
          endif
      enddo
  endif

  call bcast_world_i(natms_tot)
  call get_sync_world

  call bcast_world_iarr(ltype,natms_tot)
  call bcast_world_farr(xxx,natms_tot)
  call bcast_world_farr(yyy,natms_tot)
  call bcast_world_farr(zzz,natms_tot)

! compute inverse mass
  do itype=1,ntype
    rmass(itype)=ONE/weight(itype)
  enddo

! Initialize linear momentum if not given in xyz file
  if (idrank==0) then
    natms = natms_tot
    if(linit_temp)call init_velocity
  endif

  do j=1,nxyzlist_sub
    k=xyzlist_sub(j)
    select case(k)
    case(7)
      call bcast_world_farr(vxx,natms_tot)
    case(8)
      call bcast_world_farr(vyy,natms_tot)
    case(9)
      call bcast_world_farr(vzz,natms_tot)
    case(10)
      call bcast_world_farr(oxx,natms_tot)
    case(11)
      call bcast_world_farr(oyy,natms_tot)
    case(12)
      call bcast_world_farr(ozz,natms_tot)
    case(13)
      call bcast_world_farr(q1,natms_tot)
    case(14)
      call bcast_world_farr(q2,natms_tot)
    case(15)
      call bcast_world_farr(q3,natms_tot)
    end select
  enddo

  call get_sync_world
  
! since the data have been sent, we deallocate
  if(allocated(ltypes))deallocate(ltypes)
  if(allocated(xxs))deallocate(xxs)
  if(allocated(yys))deallocate(yys)
  if(allocated(zzs))deallocate(zzs)
  if(allocated(ots))deallocate(ots)
  

  ! Count my atoms, put them in atmbook list
  ! Also count halo atoms, at end of atmbook
  natms = 0; num_ext = 0

  do i=1,natms_tot
#ifndef DEOWERN
    ids=ownernfind_arr(nint(xxx(i)),nint(yyy(i)),nint(zzz(i)), &
        nx,ny,nz,nbuff,ownern)
#else
    ids=ownernfind(nint(xxx(i)),nint(yyy(i)),nint(zzz(i)), &
        mxrank,gminx,gmaxx,gminy,gmaxy,gminz,gmaxz)
#endif
    if (idrank==ids) then
        natms = natms + 1
        atmbook(natms) = i
    else
        call checkIfExtAtom(i, num_ext)
    endif
  enddo

  ! ATTENTION DOMAIN DECOMPOSITION SHOULD BE ADDED HERE AND natms_ext PROPERLY SET
  ! 1) First parallelization with all atoms on all processes FB
  do iatm=1, num_ext
    atmbook(natms+iatm) = atmbook(mxatms - num_ext + iatm)
  enddo
  natms_ext= natms + num_ext

#if 1
  call print_all_particles(100,'atomSetup',1)
#endif
  
! total degree of freedom of particles
  if(lrotate)then
    degfre=SIX*real(natms,kind=PRC)
  else
    degfre=THREE*real(natms,kind=PRC)
  endif

  
! start the initialization of the rotation part
  
  if(lrotate)then
    !initialize moments of inertia
    where(ishape(1:ntype)==1)
      !oblate
      rotinx(1:ntype)=ONE/FIVE*weight(1:ntype)*(rdimy(1:ntype)+rdimz(1:ntype))**TWO
      rotiny(1:ntype)=ONE/FIVE*weight(1:ntype)*(rdimx(1:ntype)+rdimz(1:ntype))**TWO
      rotinz(1:ntype)=ONE/FIVE*weight(1:ntype)*(rdimx(1:ntype)+rdimy(1:ntype))**TWO
    elsewhere
      !sphere
      rotinx(1:ntype)=TWO/FIVE*weight(1:ntype)*(rdimx(1:ntype))**TWO
      rotiny(1:ntype)=TWO/FIVE*weight(1:ntype)*(rdimx(1:ntype))**TWO
      rotinz(1:ntype)=TWO/FIVE*weight(1:ntype)*(rdimx(1:ntype))**TWO
    end where
  endif
  
  lqinput=.false.
  if(nxyzlist_sub>0)then
    do j=1,nxyzlist_sub
      k=xyzlist_sub(j)
      if(k>=13 .and. k<=15)lqinput=.true.
    enddo
  endif
  
  if(lrotate)then
    if(.not. lqinput)then
      call warning(19)
      !transform the rotational matrix in quaternions using an uniform
      !distribution
      do myi=1,natms_ext
        i = atmbook(myi)
        xsubm(1:3)=ZERO
        ysubm(1:3)=ZERO
        zsubm(1:3)=ZERO
        xsubm(1)=ONE
        ysubm(2)=ONE
        zsubm(3)=ONE
        call random_number(dmio(1))
        call random_number(dmio(2))
        call random_number(dmio(3))
        dmio(1)=(dmio(1)-HALF)*TWO*Pi
        dmio(2)=(dmio(2)-HALF)*Pi
        dmio(3)=(dmio(3)-HALF)*TWO*Pi
        call eul2q(dmio(1),dmio(3),dmio(2),qs)
        q0(i)=qs(0)
        q1(i)=qs(1)
        q2(i)=qs(2)
        q3(i)=qs(3)
#ifdef CHECKQUAT
        call q2eul(qs,dmio2(1),dmio2(3),dmio2(2))
        do j=1,3
          if(abs(dmio(j)-dmio2(j))>toll)then
            ltestrot=.true.
          endif
        enddo
#endif
      enddo
    else
      !transform the Euler angles to quaternions
      do myi=1,natms_ext
        i = atmbook(myi)
        dmio(1)=q1(i)
        dmio(2)=q2(i)
        dmio(3)=q3(i)
        call eul2q(dmio(1),dmio(3),dmio(2),qs)
        q0(i)=qs(0)
        q1(i)=qs(1)
        q2(i)=qs(2)
        q3(i)=qs(3)
#ifdef CHECKQUAT
        call q2eul(qs,dmio2(1),dmio2(3),dmio2(2))
        do j=1,3
          if(abs(dmio(j)-dmio2(j))>toll)then
            ltestrot=.true.
          endif
        enddo
#endif
      enddo
    endif
  
#ifdef CHECKQUAT
    call or_world_larr(ltestrot,1)
    if(ltestrot(1))then
      call error(22)
    endif
#endif
    
    do myi=1,natms_ext
      i = atmbook(myi)
!      write (6,*) __FILE__,__LINE__, "i=", i, myi<=natms, xxx(i),yyy(i),zzz(i), vxx(i),vyy(i),vzz(i)
      rnorm=ONE/sqrt(q0(i)**TWO+q1(i)**TWO+q2(i)**TWO+q3(i)**TWO)
      q0(i)=q0(i)*rnorm
      q1(i)=q1(i)*rnorm
      q2(i)=q2(i)*rnorm
      q3(i)=q3(i)*rnorm
    enddo
    
  endif
  
  if(lparticles)then
    if(.not. lbc_halfway)then
      call set_lbc_halfway(.true.)
      call warning(26)
    endif
    
  endif
 end subroutine initialize_map_particles
 

 subroutine checkIfExtAtom(i, num_ext)
  implicit none
  integer, intent(in) :: i
  integer, intent(inout) :: num_ext
  real(kind=PRC) :: radius
  integer :: extrema(3,8), xm,ym,zm, xp,yp,zp
  integer :: ids, c, onEnterNum



  onEnterNum = num_ext

  radius = floor(rdimx(1))

  xm = pimage(ixpbc, nint(xxx(i)-radius), nx)
  xp = pimage(ixpbc, nint(xxx(i)+radius), nx)
  ym = pimage(iypbc, nint(yyy(i)-radius), ny)
  yp = pimage(iypbc, nint(yyy(i)+radius), ny)
  zm = pimage(izpbc, nint(zzz(i)-radius), nz)
  zp = pimage(izpbc, nint(zzz(i)+radius), nz)

  extrema(:,1) = [ xm, ym, zm ]
  extrema(:,2) = [ xm, ym, zp ]
  extrema(:,3) = [ xm, yp, zm ]
  extrema(:,4) = [ xm, yp, zp ]
  extrema(:,5) = [ xp, ym, zm ]
  extrema(:,6) = [ xp, ym, zp ]
  extrema(:,7) = [ xp, yp, zm ]
  extrema(:,8) = [ xp, yp, zp ]

  do c = 1,8
!    if (extrema(1,c) < 0) extrema(1,c)=0
!    if (extrema(2,c) < 0) extrema(2,c)=0
!    if (extrema(3,c) < 0) extrema(3,c)=0
!
!    if (extrema(1,c) > nx+nbuff) extrema(1,c)=nx+nbuff
!    if (extrema(2,c) > ny+nbuff) extrema(2,c)=ny+nbuff
!    if (extrema(3,c) > nz+nbuff) extrema(3,c)=nz+nbuff

    CYCLE_OUT_INTERVAL(extrema(1,c), minx-1, maxx+1)
    CYCLE_OUT_INTERVAL(extrema(2,c), miny-1, maxy+1)
    CYCLE_OUT_INTERVAL(extrema(3,c), minz-1, maxz+1)

    atmbook(mxatms - num_ext) = i
    num_ext = num_ext + 1
    exit
  enddo

 end subroutine checkIfExtAtom


 subroutine init_particles_fluid_interaction
  
!***********************************************************************
!     
!     LBsoft subroutine to initialize particle templates
!     and particle fluid interactions
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
 
  implicit none
  integer :: myi, iatm,i,j,k,l
  
  if(.not. lparticles)return
  
  call spherical_template(rdimx(1),nsphere,spherelist,spheredist, &
   nspheredead,spherelistdead)
 

 ! initialize isfluid according to the particle presence
  do myi=1,natms_ext
    iatm = atmbook(myi)
    i=nint(xxx(iatm))
    j=nint(yyy(iatm))
    k=nint(zzz(iatm))
    call init_particle_2_isfluid(i,j,k,nsphere,spherelist,spheredist, &
     nspheredead,spherelistdead, myi <= natms)
  enddo
  
 !push the isfluid communication if necessary
  call push_comm_isfluid
  
#if 0
  k=floor(rdimx(1))
  write(6,*)nsphere
  do i=1,nsphere
    write(6,*)i,spherelist(1:3,i),spheredist(i)
    if(spherelist(1,i)>k .or. spherelist(2,i)>k .or. spherelist(3,i)>k)then
      write(6,*)'error'
      exit
    endif
    if(spherelist(1,i)<-k .or. spherelist(2,i)<-k .or. spherelist(3,i)<-k)then
      write(6,*)'error'
      exit
    endif
  enddo
  
  write(6,*)nspheredead
  do i=1,nspheredead
    write(6,*)i,spherelistdead(1:3,i)
    if(spherelistdead(1,i)>k .or. spherelistdead(2,i)>k .or. spherelistdead(3,i)>k)then
      write(6,*)'error'
      exit
    endif
    if(spherelistdead(1,i)<-k .or. spherelistdead(2,i)<-k .or. spherelistdead(3,i)<-k)then
      write(6,*)'error'
      exit
    endif
  enddo
  
  call finalize_world
  stop
#endif
  
 end subroutine init_particles_fluid_interaction
 

 subroutine apply_particle_bounce_back(nstep)
 
!***********************************************************************
!     
!     LBsoft subroutine to apply the particle bounce back on
!     fluids
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
 
  implicit none
  integer, intent(in) :: nstep
  integer :: myi, iatm,i,j,k,itype
  real(kind=PRC) :: myrot(9),oat(0:3),qtemp(0:3),qversor(0:3)


  do myi=1,natms_ext
    iatm = atmbook(myi)
    fxb(iatm)=ZERO
    fyb(iatm)=ZERO
    fzb(iatm)=ZERO
  enddo
  
  if(lrotate)then
  
  do myi=1,natms_ext
    iatm = atmbook(myi)
!    write (6,*) __FILE__,__LINE__, "iatm=", iatm
    txb(iatm)=ZERO
    tyb(iatm)=ZERO
    tzb(iatm)=ZERO

    i=nint(xxx(iatm))
    j=nint(yyy(iatm))
    k=nint(zzz(iatm))
    itype=ltype(iatm)
    !transform oxx oyy ozz from body ref to world ref
    qtemp(0)=q0(iatm)
    qtemp(1)=q1(iatm)
    qtemp(2)=q2(iatm)
    qtemp(3)=q3(iatm)
    qversor(0)=ZERO
    qversor(1)=oxx(iatm)
    qversor(2)=oyy(iatm)
    qversor(3)=ozz(iatm)
    oat=qtrimult(qtemp,qversor,qconj(qtemp))
    
    call particle_bounce_back(.true., nstep,iatm,myi<=natms, lrotate,i,j,k,nsphere, &
     spherelist,spheredist,rdimx(itype),rdimy(itype),rdimz(itype), &
     xxx(iatm),yyy(iatm),zzz(iatm), &
     vxx(iatm),vyy(iatm),vzz(iatm), &
     fxb(iatm),fyb(iatm),fzb(iatm),oat(1),oat(2),oat(3), &
     txb(iatm),tyb(iatm),tzb(iatm))

     write (6,*) __FILE__,__LINE__, myi <= natms, "iatm=", iatm, &
        "f=", fxb(iatm),fyb(iatm),fzb(iatm), &
        "t=", txb(iatm),tyb(iatm),tzb(iatm)
  enddo
  
  else
  
  do myi=1,natms_ext
    iatm = atmbook(myi)

    i=nint(xxx(iatm))
    j=nint(yyy(iatm))
    k=nint(zzz(iatm))
    itype=ltype(iatm)
    call particle_bounce_back(.false., nstep,iatm,myi<=natms, lrotate,i,j,k,nsphere, &
     spherelist,spheredist,rdimx(itype),rdimy(itype),rdimz(itype), &
     xxx(iatm),yyy(iatm),zzz(iatm), &
     vxx(iatm),vyy(iatm),vzz(iatm), &
     fxb(iatm),fyb(iatm),fzb(iatm))

!     write (6,*) __FILE__,__LINE__, myi <= natms, "iatm=", iatm, &
!        "f=", fxb(iatm),fyb(iatm),fzb(iatm)
  enddo

  endif

 end subroutine apply_particle_bounce_back
 

 subroutine force_particle_bounce_back()
 
!***********************************************************************
!     
!     LBsoft subroutine to compute the force of bounceback acting
!     on the particles at time t (bounce force is computed at t+1/2)
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
 
  implicit none
  integer :: iatm, myi
  

  !f(t) = (f(t+1/2)+f(t-1/2))/2
  do myi=1,natms
    iatm = atmbook(myi)
    fxx(iatm)=fxx(iatm)+(fxb(iatm)+fxbo(iatm))*HALF
    fyy(iatm)=fyy(iatm)+(fyb(iatm)+fybo(iatm))*HALF
    fzz(iatm)=fzz(iatm)+(fzb(iatm)+fzbo(iatm))*HALF
  enddo
  
  do myi=1,natms
    iatm = atmbook(myi)
    fxbo(iatm)=fxb(iatm)
    fybo(iatm)=fyb(iatm)
    fzbo(iatm)=fzb(iatm)
  enddo
  
  if(.not. lrotate)return
  
!  forall(iatm=1:natms,abs(txb(iatm))>0.1d0)txb(iatm)=txb(iatm)/abs(txb(iatm))*0.1d0
!  forall(iatm=1:natms,abs(tyb(iatm))>0.1d0)tyb(iatm)=tyb(iatm)/abs(tyb(iatm))*0.1d0
!  forall(iatm=1:natms,abs(tzb(iatm))>0.1d0)tzb(iatm)=tzb(iatm)/abs(tzb(iatm))*0.1d0
  
  !f(t) = (f(t+1/2)+f(t-1/2))/2
  do myi=1,natms
    iatm = atmbook(myi)
    tqx(iatm)=tqx(iatm)+(txb(iatm)+txbo(iatm))*HALF
    tqy(iatm)=tqy(iatm)+(tyb(iatm)+tybo(iatm))*HALF
    tqz(iatm)=tqz(iatm)+(tzb(iatm)+tzbo(iatm))*HALF
  enddo
  
  do myi=1,natms
    iatm = atmbook(myi)
    txbo(iatm)=txb(iatm)
    tybo(iatm)=tyb(iatm)
    tzbo(iatm)=tzb(iatm)
  enddo
 end subroutine force_particle_bounce_back

 
 subroutine merge_particle_force()
 
!***********************************************************************
!     
!     LBsoft subroutine to merge all the particle forces
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
  
  implicit none
  integer :: iatm, myi


  call sum_world_farr(fxb, natms_tot)
  call sum_world_farr(fyb, natms_tot)
  call sum_world_farr(fzb, natms_tot)
  call sum_world_farr(txb, natms_tot)
  call sum_world_farr(tyb, natms_tot)
  call sum_world_farr(tzb, natms_tot)
  
  do iatm=1,natms_tot
    write (6,*) __FILE__,__LINE__, "iatm=", iatm, &
        "f=", fxb(iatm),fyb(iatm),fzb(iatm), &
        "t=", txb(iatm),tyb(iatm),tzb(iatm)
  enddo

  do myi=natms+1,natms_ext
    iatm = atmbook(myi)
    fxb(iatm) = 0
    fyb(iatm) = 0
    fzb(iatm) = 0
    txb(iatm) = 0
    tyb(iatm) = 0
    tzb(iatm) = 0
  enddo

!  do myi=1,natms_ext
!    iatm = atmbook(myi)
!    write (6,*) __FILE__,__LINE__, myi <= natms, "iatm=", iatm, &
!        "f=", fxb(iatm),fyb(iatm),fzb(iatm), &
!        "t=", txb(iatm),tyb(iatm),tzb(iatm)
!  enddo

 end subroutine merge_particle_force

  
 subroutine check_moving_particles
  
!***********************************************************************
!     
!     LBsoft subroutine to initialize particle templates
!     and particle fluid interactions
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
 
  implicit none
  integer :: iatm,myi, i
  

  lmove(1:natms_tot) = .false.

  do myi=1,natms
    iatm=atmbook(myi)
    if (nint(xxx(iatm)) /= nint(xxo(iatm)) .or. &
        nint(yyy(iatm)) /= nint(yyo(iatm)) .or. &
        nint(zzz(iatm)) /= nint(zzo(iatm))) then
      lmove(iatm) = .true.
    endif
  enddo
  
  call or1_world_larr(lmove, natms_tot)

!  SERIAL VERSION
! ! check all the particle centers if they are moved
! where(nint(xxx(1:natms_ext))/=nint(xxo(1:natms_ext)).or. &
!  nint(yyy(1:natms_ext))/=nint(yyo(1:natms_ext)) .or. &
!  nint(zzz(1:natms_ext))/=nint(zzo(1:natms_ext)))
!   lmove(1:natms_ext)=.true.
! elsewhere
!   lmove(1:natms_ext)=.false.
! end where
! if(mxrank>1)then
!   !check if the particle is leaving the sub domain
!   forall(i=1:natms,lmove(i))
!     lmove_dom(i)=(nint(xxx(i))<minx .or. nint(xxx(i))>maxx .or. &
!      nint(yyy(i))<miny .or. nint(yyy(i))>maxy .or. &
!      nint(zzz(i))<minz .or. nint(zzz(i))>maxz)
!   end forall
!   !check if the particle is entering in the sub domain
!   !NOTE: particle in the halo are from natms+1 up to natms_ext
!   if(natms_ext>=natms+1)then
!     forall(i=natms+1:natms_ext,lmove(i))
!       lmove_dom(i)=(nint(xxx(i))>=minx .and. nint(xxx(i))<=maxx .and. &
!        nint(yyy(i))>=miny .and. nint(yyy(i))<=maxy .and. &
!        nint(zzz(i))>=minz .and. nint(zzz(i))<=maxz)
!     end forall
!   endif
! endif
!endif
  
!  write (6,*) __FILE__,__LINE__, "lmove=", lmove(1:natms_ext)
 end subroutine check_moving_particles
 

 subroutine build_new_isfluid(nstep)
  
!***********************************************************************
!     
!     LBsoft subroutine to build the new_isfluid array if necessary
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
  
  implicit none
  integer, intent(in) :: nstep
 

  call initialize_new_isfluid
    
  call check_moving_particles
  
  call mapping_new_isfluid(natms,natms_ext,atmbook,nsphere, &
     spherelist,spheredist,nspheredead,spherelistdead,lmove, &
     xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz,xxo,yyo,zzo)
     
  call particle_delete_fluids(nstep,natms_ext,atmbook,nsphere, &
     spherelist,spheredist,nspheredead,spherelistdead,lmove,lrotate, &
     ltype,xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz,tqx,tqy,tqz,xxo,yyo,zzo, &
     rdimx,rdimy,rdimz)
  
 end subroutine build_new_isfluid
 
 subroutine inter_part_and_grid(nstep, lparticles)
  
!***********************************************************************
!     
!     LBsoft subroutine to manage moving particle effects
!     to fluid nodes and create-delete fluid nodes if necessary
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
 
  implicit none
  integer, intent(in) :: nstep
  logical, intent(in) :: lparticles
  
  call particle_create_fluids(nstep,natms_ext,atmbook,nsphere, &
     spherelist,spheredist,nspheredead,spherelistdead,lmove,lrotate, &
     ltype,xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz,tqx,tqy,tqz,xxo,yyo,zzo, &
     rdimx,rdimy,rdimz)

  call driver_bc_densities
  ! call driver_bc_pops(.false.)
  call driver_bc_pops_NOK

  call update_isfluid
  call driver_bc_isfluid
 end subroutine inter_part_and_grid
 

 subroutine clean_fluid_inside_particle
   
!***********************************************************************
!     
!     LBsoft subroutine to remove fluid inside the particles
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
 
  implicit none
  
  if(.not. lparticles)return
 
  call erase_fluids_in_particles(natms,atmbook, nsphere,spherelist, &
   spheredist,nspheredead,spherelistdead,ltype,xxx,yyy,zzz, &
   vxx,vyy,vzz,rdimx,rdimy,rdimz)
  
 end subroutine clean_fluid_inside_particle
 
 subroutine compute_psi_sc_particles(nstep)
  
!***********************************************************************
!     
!     LBsoft subroutine to compute the Shan Chen pseudo potential
!     in presence of particles
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification January 2019
!     
!***********************************************************************
  
  implicit none
  integer, intent(in) :: nstep
  integer :: iatm,i,j,k,itype
  real(kind=PRC) :: myrot(9),oat(0:3),qtemp(0:3),qversor(0:3)
  

  if(lrotate)then
  
  do iatm=1,natms_ext
    i=nint(xxx(iatm))
    j=nint(yyy(iatm))
    k=nint(zzz(iatm))
    itype=ltype(iatm)
    !transform xversor from body ref to world ref
    qtemp(0)=q0(iatm)
    qtemp(1)=q1(iatm)
    qtemp(2)=q2(iatm)
    qtemp(3)=q3(iatm)
    qversor(0:3)=ZERO
    qversor(1)=ONE
    oat=qtrimult(qtemp,qversor,qconj(qtemp))
    call compute_sc_particle_interact(nstep,iatm, iatm<=natms, lrotate,i,j,k,nsphere, &
     spherelist,spheredist,rdimx(itype),rdimy(itype),rdimz(itype), &
     xxx(iatm),yyy(iatm),zzz(iatm), &
     vxx(iatm),vyy(iatm),vzz(iatm), &
     fxx(iatm),fyy(iatm),fzz(iatm),oat(1),oat(2),oat(3), &
     tqx(iatm),tqy(iatm),tqz(iatm))
  enddo
  
  else
  
  do iatm=1,natms
    i=nint(xxx(iatm))
    j=nint(yyy(iatm))
    k=nint(zzz(iatm))
    itype=ltype(iatm)
    call compute_sc_particle_interact(nstep,iatm,iatm<=natms, lrotate,i,j,k,nsphere, &
     spherelist,spheredist,rdimx(itype),rdimy(itype),rdimz(itype), &
     xxx(iatm),yyy(iatm),zzz(iatm), &
     vxx(iatm),vyy(iatm),vzz(iatm), &
     fxx(iatm),fyy(iatm),fzz(iatm))
  enddo
  
  endif
 end subroutine compute_psi_sc_particles
 
 
 subroutine vertest(newlst)
      
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
  integer :: myi, i,moved(1)
  integer :: fail=0
  real(kind=PRC) :: rmax,dr, checkSpace
  real(kind=PRC), allocatable, save :: xold(:),yold(:),zold(:)
  

  checkSpace = (natms+mxrank-1)/mxrank
  if(checkSpace > mxatms) call error(24)
      
  if(newjob)then
!   set up initial arrays 
    allocate (xold(mxatms),yold(mxatms),zold(mxatms),stat=fail)
    if(fail.ne.0)call error(25)

    do myi=1,natms
      i = atmbook(myi)
      xold(i)=ZERO
      yold(i)=ZERO
      zold(i)=ZERO
    enddo

    newjob=.false.
    newlst=.true.

    return
  endif
    
!   integrate velocities 
    do myi=1,natms
      i = atmbook(myi)
      xold(i)=xold(i)+vxx(i)
      yold(i)=yold(i)+vyy(i)
      zold(i)=zold(i)+vzz(i)
    enddo
    
!   maximum displacement 
    rmax=(delr/TWO)**TWO
    
!   test atomic displacements
    moved(1)=0
    do myi=1,natms
      i = atmbook(myi)
      dr=tstepatm**2*(xold(i)**2+yold(i)**2+zold(i)**2)
      if(dr.gt.rmax)moved(1)=moved(1)+1
    enddo
        
!   global sum of moved atoms
    call sum_world_iarr(moved,1)
        
!   test for new verlet list
    newlst=(moved(1).ge.2)
        
!   update stored positions
    if(newlst)then
      do myi=1,natms
        i = atmbook(myi)
        xold(i)=ZERO
        yold(i)=ZERO
        zold(i)=ZERO
      enddo
    endif
 end subroutine vertest

 
 subroutine parlst
  
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
! separation vectors and powers thereof
  real(kind=PRC) :: rsq,xm,ym,zm,rsqcut
  
! Loop counters
  integer :: iatm,ii,i, myi, myj
  integer :: k,jatm
  integer :: itype,ilentry
  logical :: lchk(1)
  integer :: ibig(1),idum



    lchk=.false.
    ibig=0
    
    rsqcut = (rcut+delr)**TWO
    
    call allocate_array_bdf(natms)
    
    lentry(1:mxatms)=0
    
    ii = 0
    do myi = 1,natms
      iatm = atmbook(myi)
      ii=myi
      k=0
      do myj = myi+1,natms
        jatm = atmbook(myj)
        k=k+1
        xdf(k)=xxx(jatm)-xxx(iatm)
        ydf(k)=yyy(jatm)-yyy(iatm)
        zdf(k)=zzz(jatm)-zzz(iatm)
      enddo
      
      call pbc_images(imcon,k,cell,xdf,ydf,zdf)
      
      k=0
      ilentry=0
      do myj = myi+1,natms
        jatm = atmbook(myj)
        k=k+1
        rsq=xdf(k)**TWO + ydf(k)**TWO + zdf(k)**TWO
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
      
!      write(6,*) myi, 'list for',iatm, ':',list(1:ilentry,ii)
    enddo
  
    
!   terminate job if neighbour list array exceeded
    call or_world_larr(lchk,1)
    if(lchk(1))then
      call max_world_iarr(ibig,1)
      call warning(23,real(ibig(1),kind=PRC))
      call error(21)
    endif
 end subroutine parlst
 

 subroutine initialize_particle_force
 
!***********************************************************************
!     
!     LBsoft subroutine to initialize the force arrays
!
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************

  implicit none
  integer :: i, myi
  
  do myi = 1,natms
    i = atmbook(myi)
    fxx(i)=ext_fxx
    fyy(i)=ext_fyy
    fzz(i)=ext_fzz
  enddo
  
  if(.not. lrotate)return

  !initialize also the torque of forces
  do myi = 1,natms
    i = atmbook(myi)
	tqx(i)=ext_tqx
    tqy(i)=ext_tqy
    tqz(i)=ext_tqz
  enddo
 end subroutine initialize_particle_force
 

 subroutine initialize_particle_energy
 
!***********************************************************************
!     
!     LBsoft subroutine to initialize the particle energy
!
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************

  implicit none
  
  
  engcfg=ZERO
  engke=ZERO
  engtot=ZERO
 end subroutine initialize_particle_energy
 

 subroutine driver_inter_force(nstepsub)
 
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
  integer, intent(in) :: nstepsub
  real(kind=PRC), parameter :: tol=real(1.d-4,kind=PRC)
  integer :: i
  
  call compute_inter_force(lentry,list)
  call compute_sidewall_force(nstepsub)
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
   kappa,rlimit,gmin,rpar,lubfactor,rparcut,rparcap,rsrparcut,mxrsqcut,&
   ux,uy,uz,visc
  
! Loop counters
  integer :: iatm,ii,i,j,ivdw,itype,jtype,iimax,xcm,ycm,zcm
  integer :: k,jatm, myi
  integer :: ilentry
  
  real(kind=PRC), parameter :: s2rmin=TWO**(ONE/SIX)
  
  call allocate_array_bdf(natms)
  
  if(lunique_omega)then
    visc=viscR
  endif
  
  do ivdw=1,ntpvdw
  
    select case(ltpvdw(ivdw))
    case (0)
      return

    case (1)
    
    !Weeks-Chandler-Andersen
    
    eps=prmvdw(1,ivdw)
    sig=prmvdw(2,ivdw)
    rmin=s2rmin*sig
    rsqcut = (rmin)**TWO
    vmin=FOUR*eps*(sig/rmin)**SIX*((sig/rmin)**SIX-ONE)
    do j=1,mxntype
      do k=1,mxntype
        if(mskvdw(j,k)==ivdw)then
          itype=j
          jtype=k
        endif
      enddo
    enddo
    rpar=rdimx(itype)+rdimx(jtype)
    lubfactor=SIX*Pi*rpar**FOUR/(FOUR*rpar**TWO)
    rparcut=rpar+TWO/THREE
    rsrparcut=rparcut**TWO
    mxrsqcut=max(rsrparcut,rsqcut)
    rparcap=rpar+ONE/TEN
    
    do iatm = 1,natms
      itype=ltype(iatm)
      if(all(mskvdw(1:ntpvdw,itype)/=ivdw))cycle
      ii = 0
      do k = 1,lentrysub(iatm)
        jatm=listsub(k,iatm)
        jtype=ltype(jatm)
        if(mskvdw(itype,jtype)/=ivdw)cycle
        ii=ii+1
        xdf(ii)=xxx(jatm)-xxx(iatm)
        ydf(ii)=yyy(jatm)-yyy(iatm)
        zdf(ii)=zzz(jatm)-zzz(iatm)
      enddo
      iimax=ii
      
      call pbc_images(imcon,iimax,cell,xdf,ydf,zdf)
      
      ii = 0
      do k = 1,lentrysub(iatm)
        jatm=listsub(k,iatm)
        jtype=ltype(jatm)
        if(mskvdw(itype,jtype)/=ivdw)cycle
        ii=ii+1
        rsq=xdf(ii)**TWO+ydf(ii)**TWO+zdf(ii)**TWO
        if(rsq<=mxrsqcut)then
          rrr=sqrt(rsq)
          if(rrr<=rmin)then
            vvv=FOUR*eps*(sig/rrr)**SIX*((sig/rrr)**SIX-ONE)+vmin
            ggg=TWENTYFOUR*eps/rrr*(sig/rrr)**SIX*(TWO*(sig/rrr)**SIX-ONE)
            engcfg=engcfg+vvv
            fxx(iatm)=fxx(iatm)-ggg*xdf(ii)/rrr
            fxx(jatm)=fxx(jatm)+ggg*xdf(ii)/rrr
            fyy(iatm)=fyy(iatm)-ggg*ydf(ii)/rrr
            fyy(jatm)=fyy(jatm)+ggg*ydf(ii)/rrr
            fzz(iatm)=fzz(iatm)-ggg*zdf(ii)/rrr
            fzz(jatm)=fzz(jatm)+ggg*zdf(ii)/rrr
          endif
          if(llubrication)then
            if(rrr<=rparcut)then
              ux=xdf(ii)/rrr
              uy=ydf(ii)/rrr
              uz=zdf(ii)/rrr
              if(rrr<=rparcap)rrr=rparcap
              if(.not. lunique_omega)then
                xcm=nint(xxx(iatm)+xdf(ii)*HALF)
                ycm=nint(yyy(iatm)+ydf(ii)*HALF)
                zcm=nint(zzz(iatm)+zdf(ii)*HALF)
                xcm=pimage(ixpbc,xcm,nx)
                ycm=pimage(iypbc,ycm,ny)
                zcm=pimage(izpbc,zcm,nz)
                visc=omega_to_viscosity(omega(xcm,ycm,zcm))
              endif
              fxx(iatm)=fxx(iatm)-lubfactor*(rdimx(itype)**FOUR)*ux* &
               (ux*(vxx(iatm)-vxx(jatm)))*(ONE/(rrr-rpar)-ONE)
              fxx(jatm)=fxx(jatm)+lubfactor*(rdimx(itype)**FOUR)*ux* &
               (ux*(vxx(iatm)-vxx(jatm)))*(ONE/(rrr-rpar)-ONE)
              fyy(iatm)=fyy(iatm)-lubfactor*(rdimx(itype)**FOUR)*uy* &
               (uy*(vyy(iatm)-vyy(jatm)))*(ONE/(rrr-rpar)-ONE)
              fyy(jatm)=fyy(jatm)+lubfactor*(rdimx(jtype)**FOUR)*uy* &
               (uy*(vyy(iatm)-vyy(jatm)))*(ONE/(rrr-rpar)-ONE)
              fzz(iatm)=fzz(iatm)-lubfactor*(rdimx(jtype)**FOUR)*uz* &
               (uz*(vzz(iatm)-vzz(jatm)))*(ONE/(rrr-rpar)-ONE)
              fzz(jatm)=fzz(jatm)+lubfactor*(rdimx(jtype)**FOUR)*uz* &
               (uz*(vzz(iatm)-vzz(jatm)))*(ONE/(rrr-rpar)-ONE)
            endif 
          endif
        endif
      enddo
    enddo
  
    case (2)
    
    !Lennard-Jones
    
    eps=prmvdw(1,ivdw)
    sig=prmvdw(2,ivdw)
    rsqcut = (rcut)**TWO
    
    do j=1,mxntype
      do k=1,mxntype
        if(mskvdw(j,k)==ivdw)then
          itype=j
          jtype=k
        endif
      enddo
    enddo
    rpar=rdimx(itype)+rdimx(jtype)
    lubfactor=SIX*Pi*rpar**FOUR/(FOUR*rpar**TWO)
    rparcut=rpar+TWO/THREE
    rsrparcut=rparcut**TWO
    mxrsqcut=max(rsrparcut,rsqcut)
    rparcap=rpar+ONE/TEN
    
    
    
    
    do iatm = 1,natms
      itype=ltype(iatm)
      if(all(mskvdw(1:ntpvdw,itype)/=ivdw))cycle
      ii = 0
      do k = 1,lentrysub(iatm)
        jatm=listsub(k,iatm)
        jtype=ltype(jatm)
        if(mskvdw(itype,jtype)/=ivdw)cycle
        ii=ii+1
        xdf(ii)=xxx(jatm)-xxx(iatm)
        ydf(ii)=yyy(jatm)-yyy(iatm)
        zdf(ii)=zzz(jatm)-zzz(iatm)
      enddo
      iimax=ii
      
      call pbc_images(imcon,iimax,cell,xdf,ydf,zdf)
      
      ii = 0
      do k = 1,lentrysub(iatm)
        jatm=listsub(k,iatm)
        jtype=ltype(jatm)
        if(mskvdw(itype,jtype)/=ivdw)cycle
        ii=ii+1
        rsq=xdf(ii)**TWO+ydf(ii)**TWO+zdf(ii)**TWO
        if(rsq<=mxrsqcut) then
          rrr=sqrt(rsq)
          vvv=FOUR*eps*(sig/rrr)**SIX*((sig/rrr)**SIX-ONE)
          ggg=TWENTYFOUR*eps/rrr*(sig/rrr)**SIX*(TWO*(sig/rrr)**SIX-ONE)
          engcfg=engcfg+vvv
          fxx(iatm)=fxx(iatm)-ggg*xdf(ii)/rrr
          fxx(jatm)=fxx(jatm)+ggg*xdf(ii)/rrr
          fyy(iatm)=fyy(iatm)-ggg*ydf(ii)/rrr
          fyy(jatm)=fyy(jatm)+ggg*ydf(ii)/rrr
          fzz(iatm)=fzz(iatm)-ggg*zdf(ii)/rrr
          fzz(jatm)=fzz(jatm)+ggg*zdf(ii)/rrr
          if(llubrication)then
            if(rrr<=rparcut)then
              ux=xdf(ii)/rrr
              uy=ydf(ii)/rrr
              uz=zdf(ii)/rrr
              if(rrr<=rparcap)rrr=rparcap
              if(.not. lunique_omega)then
                xcm=nint(xxx(iatm)+xdf(ii)*HALF)
                ycm=nint(yyy(iatm)+ydf(ii)*HALF)
                zcm=nint(zzz(iatm)+zdf(ii)*HALF)
                xcm=pimage(ixpbc,xcm,nx)
                ycm=pimage(iypbc,ycm,ny)
                zcm=pimage(izpbc,zcm,nz)
                visc=omega_to_viscosity(omega(xcm,ycm,zcm))
              endif
              fxx(iatm)=fxx(iatm)-lubfactor*(rdimx(itype)**FOUR)*ux* &
               (ux*(vxx(iatm)-vxx(jatm)))*(ONE/(rrr-rpar)-ONE)
              fxx(jatm)=fxx(jatm)+lubfactor*(rdimx(itype)**FOUR)*ux* &
               (ux*(vxx(iatm)-vxx(jatm)))*(ONE/(rrr-rpar)-ONE)
              fyy(iatm)=fyy(iatm)-lubfactor*(rdimx(itype)**FOUR)*uy* &
               (uy*(vyy(iatm)-vyy(jatm)))*(ONE/(rrr-rpar)-ONE)
              fyy(jatm)=fyy(jatm)+lubfactor*(rdimx(jtype)**FOUR)*uy* &
               (uy*(vyy(iatm)-vyy(jatm)))*(ONE/(rrr-rpar)-ONE)
              fzz(iatm)=fzz(iatm)-lubfactor*(rdimx(jtype)**FOUR)*uz* &
               (uz*(vzz(iatm)-vzz(jatm)))*(ONE/(rrr-rpar)-ONE)
              fzz(jatm)=fzz(jatm)+lubfactor*(rdimx(jtype)**FOUR)*uz* &
               (uz*(vzz(iatm)-vzz(jatm)))*(ONE/(rrr-rpar)-ONE)
            endif 
          endif
        endif
      enddo
    enddo
    
    case (3)
  
    kappa=prmvdw(1,ivdw)
    rmin=prmvdw(2,ivdw)
    rlimit=prmvdw(3,ivdw)
    rsqcut = (rmin)**TWO
    vmin=kappa*(rmin-rlimit)**(FIVE*HALF)
    gmin=FIVE*HALF*kappa*(rmin-rlimit)**(THREE*HALF)
    do j=1,mxntype
      do k=1,mxntype
        if(mskvdw(j,k)==ivdw)then
          itype=j
          jtype=k
        endif
      enddo
    enddo
    rpar=rdimx(itype)+rdimx(jtype)
    lubfactor=SIX*Pi/(rpar**TWO)
    rparcut=rpar+TWO/THREE
    rsrparcut=rparcut**TWO
    mxrsqcut=max(rsrparcut,rsqcut)
    rparcap=rpar+ONE/TEN
    
  
    do myi = 1,natms
      iatm = atmbook(myi)
      itype=ltype(iatm)
      if(all(mskvdw(1:ntpvdw,itype)/=ivdw))cycle
      ii = 0
      do k = 1,lentrysub(myi)
        jatm=listsub(k,myi)
        jtype=ltype(jatm)
        if(mskvdw(itype,jtype)/=ivdw)cycle
        ii=ii+1
        xdf(ii)=xxx(jatm)-xxx(iatm)
        ydf(ii)=yyy(jatm)-yyy(iatm)
        zdf(ii)=zzz(jatm)-zzz(iatm)
      enddo
      iimax=ii
    
      call pbc_images(imcon,iimax,cell,xdf,ydf,zdf)
      
      ii = 0
      do k = 1,lentrysub(myi)
        jatm=listsub(k,myi)
        jtype=ltype(jatm)
        if(mskvdw(itype,jtype)/=ivdw)cycle
        ii=ii+1
        rsq=xdf(ii)**TWO+ydf(ii)**TWO+zdf(ii)**TWO

        ! FABIO: MPI HACK
        if (rsq <= ZERO) then
            write (6,*) __FILE__,__LINE__, "idrank=",idrank, "rsq=", rsq
            rsq = 0.5
        endif

        if(rsq<=mxrsqcut) then
          rrr=sqrt(rsq)
          if(rrr<=rmin)then
            if(rrr<=rlimit)then
              vvv=vmin+(rlimit-rrr)*gmin
              ggg=gmin
            else
              vvv=kappa*(rmin-rrr)**(FIVE*HALF)
              ggg=FIVE*HALF*kappa*(rmin-rrr)**(THREE*HALF)
            endif
            engcfg=engcfg+vvv
            fxx(iatm)=fxx(iatm)-ggg*xdf(ii)/rrr
            fxx(jatm)=fxx(jatm)+ggg*xdf(ii)/rrr
            fyy(iatm)=fyy(iatm)-ggg*ydf(ii)/rrr
            fyy(jatm)=fyy(jatm)+ggg*ydf(ii)/rrr
            fzz(iatm)=fzz(iatm)-ggg*zdf(ii)/rrr
            fzz(jatm)=fzz(jatm)+ggg*zdf(ii)/rrr
          endif

          if(llubrication)then
            if(rrr<=rparcut)then
              ux=xdf(ii)/rrr
              uy=ydf(ii)/rrr
              uz=zdf(ii)/rrr
              if(rrr<=rparcap)rrr=rparcap
              if(.not. lunique_omega)then
                xcm=nint(xxx(iatm)+xdf(ii)*HALF)
                ycm=nint(yyy(iatm)+ydf(ii)*HALF)
                zcm=nint(zzz(iatm)+zdf(ii)*HALF)
                xcm=pimage(ixpbc,xcm,nx)
                ycm=pimage(iypbc,ycm,ny)
                zcm=pimage(izpbc,zcm,nz)
                visc=omega_to_viscosity(omega(xcm,ycm,zcm))
              endif
              fxx(iatm)=fxx(iatm)-lubfactor*(rdimx(itype)**FOUR)*ux* &
               (ux*(vxx(iatm)-vxx(jatm)))*(ONE/(rrr-rpar)-ONE)
              fxx(jatm)=fxx(jatm)+lubfactor*(rdimx(itype)**FOUR)*ux* &
               (ux*(vxx(iatm)-vxx(jatm)))*(ONE/(rrr-rpar)-ONE)
              fyy(iatm)=fyy(iatm)-lubfactor*(rdimx(itype)**FOUR)*uy* &
               (uy*(vyy(iatm)-vyy(jatm)))*(ONE/(rrr-rpar)-ONE)
              fyy(jatm)=fyy(jatm)+lubfactor*(rdimx(jtype)**FOUR)*uy* &
               (uy*(vyy(iatm)-vyy(jatm)))*(ONE/(rrr-rpar)-ONE)
              fzz(iatm)=fzz(iatm)-lubfactor*(rdimx(jtype)**FOUR)*uz* &
               (uz*(vzz(iatm)-vzz(jatm)))*(ONE/(rrr-rpar)-ONE)
              fzz(jatm)=fzz(jatm)+lubfactor*(rdimx(jtype)**FOUR)*uz* &
               (uz*(vzz(iatm)-vzz(jatm)))*(ONE/(rrr-rpar)-ONE)
            endif 

          endif

        endif
      enddo
    enddo
    
    case default
      call error(27)
    end select
  
  enddo
  
  return

 end subroutine compute_inter_force
 
 
 
 subroutine merge_particle_energies()
 
  implicit none
  
  real(kind=PRC) :: fsum(2)
  
! all reduce sum
  if(mxrank>1)then
    fsum(1)=engke
    fsum(2)=engcfg
    call sum_world_farr(fsum,2)
    engke=fsum(1)
    engcfg=fsum(2)
  endif
  
! calculate total energy  
  engtot=engke+engcfg+engrot
  
  return
  
 end subroutine merge_particle_energies
 
 subroutine compute_sidewall_force(nstepsub)
  
!***********************************************************************
!     
!     LBsoft subroutine to compute the interaction force with
!     the side wall if necessary
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: nstepsub
  
  integer :: iatm,itype, myi
  real(kind=PRC) :: inflimit,suplimit(3),vmin,gmin,vvv,ggg,rrr,rlimit, &
   suprcup(3)
  logical, dimension(1) :: ltest
  
  if(.not. lsidewall_md)return
  
  ltest(1)=.false.
  inflimit=sidewall_md_rdist+HALF
  rlimit=ONE+HALF
  vmin=sidewall_md_k*(ONE**(FIVE*HALF))
  gmin=FIVE*HALF*sidewall_md_k*(ONE**(THREE*HALF))
  suplimit(1)=cell(1)-sidewall_md_rdist-HALF
  suplimit(2)=cell(5)-sidewall_md_rdist-HALF
  suplimit(3)=cell(9)-sidewall_md_rdist-HALF
  suprcup(1)=cell(1)-ONE-HALF
  suprcup(2)=cell(5)-ONE-HALF
  suprcup(3)=cell(9)-ONE-HALF
  
  do myi=1,natms
    iatm = atmbook(myi)
    itype=ltype(iatm)

    if(ixpbc==0) then
      if(xxx(iatm)-rdimx(itype)<=inflimit)then
        rrr=xxx(iatm)-rdimx(itype)
        if(rrr<=rlimit)then
          vvv=vmin+(rlimit-rrr)*gmin
          ggg=gmin
          ltest(1)=.true.
        else
          vvv=sidewall_md_k*(inflimit-rrr)**(FIVE*HALF)
          ggg=FIVE*HALF*sidewall_md_k*(inflimit-rrr)**(THREE*HALF)
        endif
        engcfg=engcfg+vvv
        fxx(iatm)=fxx(iatm)+ggg
      endif

      if(xxx(iatm)+rdimx(itype)>=suplimit(1))then
        rrr=xxx(iatm)+rdimx(itype)
        if(rrr>=suprcup(1))then
          vvv=vmin+(rrr-suprcup(1))*gmin
          ggg=gmin
          ltest(1)=.true.
        else
          vvv=sidewall_md_k*(rrr-suplimit(1))**(FIVE*HALF)
          ggg=FIVE*HALF*sidewall_md_k*(rrr-suplimit(1))**(THREE*HALF)
        endif
        engcfg=engcfg+vvv
        fxx(iatm)=fxx(iatm)-ggg
      endif
    endif
    
    if(iypbc==0) then
      if(yyy(iatm)-rdimy(itype)<=inflimit)then
        rrr=yyy(iatm)-rdimy(itype)
        if(rrr<=rlimit)then
          vvv=vmin+(rlimit-rrr)*gmin
          ggg=gmin
          ltest(1)=.true.
        else
          vvv=sidewall_md_k*(inflimit-rrr)**(FIVE*HALF)
          ggg=FIVE*HALF*sidewall_md_k*(inflimit-rrr)**(THREE*HALF)
        endif
        engcfg=engcfg+vvv
        fyy(iatm)=fyy(iatm)+ggg
      endif

      if(yyy(iatm)+rdimy(itype)>=suplimit(2))then
        rrr=yyy(iatm)+rdimy(itype)
        if(rrr>=suprcup(2))then
          vvv=vmin+(rrr-suprcup(2))*gmin
          ggg=gmin
          ltest(1)=.true.
          write(6,*)nstepsub,yyy(iatm),yyy(iatm)+rdimy(itype),-ggg
        else
          vvv=sidewall_md_k*(rrr-suplimit(2))**(FIVE*HALF)
          ggg=FIVE*HALF*sidewall_md_k*(rrr-suplimit(2))**(THREE*HALF)
        endif
        engcfg=engcfg+vvv
        fyy(iatm)=fyy(iatm)-ggg
      endif
    endif
    
    if(izpbc==0) then
      if(zzz(iatm)-rdimz(itype)<=inflimit)then
        rrr=zzz(iatm)-rdimz(itype)
        if(rrr<=rlimit)then
          vvv=vmin+(rlimit-rrr)*gmin
          ggg=gmin
          ltest(1)=.true.
        else
          vvv=sidewall_md_k*(inflimit-rrr)**(FIVE*HALF)
          ggg=FIVE*HALF*sidewall_md_k*(inflimit-rrr)**(THREE*HALF)
        endif
        engcfg=engcfg+vvv
        fzz(iatm)=fzz(iatm)+ggg
      endif

      if(zzz(iatm)+rdimz(itype)>=suplimit(3))then
        rrr=zzz(iatm)+rdimz(itype)
        if(rrr>=suprcup(3))then
          vvv=vmin+(rrr-suprcup(3))*gmin
          ggg=gmin
          ltest(1)=.true.
        else
          vvv=sidewall_md_k*(rrr-suplimit(3))**(FIVE*HALF)
          ggg=FIVE*HALF*sidewall_md_k*(rrr-suplimit(3))**(THREE*HALF)
        endif
        engcfg=engcfg+vvv
        fzz(iatm)=fzz(iatm)-ggg
      endif
    endif
    
  enddo
  
  call or_world_larr(ltest,1)
  if(ltest(1))call warning(49)
  
 end subroutine compute_sidewall_force
 

 subroutine init_velocity()
      
!***********************************************************************
!     
!     dl_poly subroutine for resetting the system velocities
!     
!     copyright - daresbury laboratory
!     author    - w. smith    may 2007
!     
!***********************************************************************

  implicit none
  integer :: i, myi
  real(kind=PRC) :: temp,tolnce,sigma
  real(kind=PRC), parameter :: tmin=real(1.d-2,kind=PRC)
      
! set atomic velocities from gaussian distribution
      
  ! MPI: This is called only from rank==0, so skip atmbook
  do i=1,natms
    sigma=sqrt(tempboltz*init_temp*rmass(ltype(i)))
    vxx(i)=sigma*gauss()
    vyy(i)=sigma*gauss()
    vzz(i)=sigma*gauss()
  enddo
  
  if(init_temp>tmin)call vscaleg()
 end subroutine init_velocity
  
  
 subroutine store_old_pos_vel_part
 
!*********************************************************************
!     
!     LBsoft subroutine for storing the old position and 
!     velocity arrays of particles
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
  
  integer :: i, myi
  
  ! store initial values of position and velocity    
  do myi=1,natms
    i = atmbook(myi)
    xxo(i)=xxx(i)
    yyo(i)=yyy(i)
    zzo(i)=zzz(i)
    vxo(i)=vxx(i)
    vyo(i)=vyy(i)
    vzo(i)=vzz(i)
  enddo
 end subroutine store_old_pos_vel_part
  

 subroutine vscaleg()

!*********************************************************************
!     
!     LBsoft subroutine for scaling the velocity arrays to the
!     desired temperature
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
      
  implicit none
  integer :: i, myi
  real(kind=PRC) :: buffer(1),sumke,ratio_temp,mytemp,ratio_sigma
  

! compute kinetic energy
  call getkin(sumke)
  
  if(mxrank>1)then
    buffer(1)=sumke
    call sum_world_farr(buffer(1),1)
    sumke=buffer(1)
  endif
  
  mytemp=TWO*(sumke+engrot)/(tempboltz*degfre)
! calculate temperature ratio
  ratio_temp=init_temp/mytemp  
  ratio_sigma=sqrt(init_temp/mytemp)
  
  do myi=1,natms
     i = atmbook(myi)
     vxx(i)=ratio_sigma*vxx(i)
     vyy(i)=ratio_sigma*vyy(i)
     vzz(i)=ratio_sigma*vzz(i)
  enddo
 end subroutine vscaleg

 
 subroutine initialize_integrator_lf()

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
  integer :: i, myi,ltypeloc
  
  
  if(lvv)return
  
  sqquattol=quattol**TWO
  
  select case(keyint)
  case (1) 
! report the atoms velocity to half timestep back for leap frog
    do myi=1,natms
      i = atmbook(myi)
      ltypeloc = ltype(i)
      vxx(i)=vxx(i)-HALF*tstepatm*rmass(ltypeloc)*fxx(i)
      vyy(i)=vyy(i)-HALF*tstepatm*rmass(ltypeloc)*fyy(i)
      vzz(i)=vzz(i)-HALF*tstepatm*rmass(ltypeloc)*fzz(i)
    enddo

  case default
    return
  end select
 end subroutine initialize_integrator_lf
 

 subroutine nve_lf(nstepsub)

!***********************************************************************
!     
!     LBsoft subroutine for integrating newtonian EOM by 
!     Verlet leapfrog
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************

  implicit none
  integer, intent(in) :: nstepsub
  integer, parameter :: nfailmax=10
  integer, dimension(nfailmax) :: fail
  integer :: i,j,k,itype,itq, myi
  logical :: ltest(1)
  real(kind=PRC) :: trx,try,trz,delx,dely,delz,engrotbuff(1)
  
! versor of particle frame 
  real(kind=PRC), dimension(3) :: uxx,uyy,yzz
  
! working variables
  real(kind=PRC) :: oqx,oqy,oqz
  real(kind=PRC) :: opx,opy,opz
  
! working arrays
  real(kind=PRC), allocatable :: bxx(:),byy(:),bzz(:)
  

! allocate working arrays
  fail(1:nfailmax)=0
  
  allocate(bxx(mxatms),stat=fail(1))
  allocate(byy(mxatms),stat=fail(2))
  allocate(bzz(mxatms),stat=fail(3))
  
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
      
! move atoms by leapfrog algorithm    
  do myi=1,natms
    i = atmbook(myi)
!   update velocities       
    bxx(i)=vxx(i)+tstepatm*rmass(ltype(i))*fxx(i)
    byy(i)=vyy(i)+tstepatm*rmass(ltype(i))*fyy(i)
    bzz(i)=vzz(i)+tstepatm*rmass(ltype(i))*fzz(i)
!   update positions
    xxx(i)=xxo(i)+tstepatm*bxx(i)
    yyy(i)=yyo(i)+tstepatm*byy(i)
    zzz(i)=zzo(i)+tstepatm*bzz(i)
  enddo
  
#if 1
  do myi=1,natms
    i = atmbook(myi)
    if(abs(bxx(i))>ONE .or. abs(byy(i))>ONE .or. abs(bzz(i))>ONE )then
      write(6,*)'WARNING - velocity cap applied at step ',nstepsub
      if(bxx(i)>cssq)bxx(i)=cssq
      if(bxx(i)<cssq)bxx(i)=-cssq
      if(byy(i)>cssq)byy(i)=cssq
      if(byy(i)<cssq)byy(i)=-cssq
      if(bzz(i)>cssq)bzz(i)=cssq
      if(bzz(i)<cssq)bzz(i)=-cssq
      !call error(-1)
    endif
  enddo
#endif
  
! calculate full timestep velocity
  do myi=1,natms
    i = atmbook(myi)
    vxx(i)=HALF*(vxx(i)+bxx(i))
    vyy(i)=HALF*(vyy(i)+byy(i))
    vzz(i)=HALF*(vzz(i)+bzz(i))
  enddo
  
! calculate kinetic energy at time step n
  call getkin(engke)
  
! periodic boundary condition
  call pbc_images_centered(imcon,natms,cell,cx,cy,cz,xxx,yyy,zzz)
    
! restore free atom half step velocity   
  do myi=1,natms
    i = atmbook(myi)
    vxx(i)=bxx(i)
    vyy(i)=byy(i)
    vzz(i)=bzz(i)
  enddo
  
  
! rigid body motion
  if(lrotate)then
    
!   initialiaze the accumulator of rotational kinetic energy
    engrotbuff(1)=ZERO
    
#ifdef SVANBERG
    call rot_midstep_impl_lf(engrotbuff(1))
#else
    call rot_fullstep_impl_lf(engrotbuff(1))
#endif
     
!   compute the total rotational kinetic energy
    call sum_world_farr(engrotbuff,1)
    engrot=HALF*engrotbuff(1)
  endif
  
! deallocate work arrays
  deallocate (bxx,byy,bzz,stat=fail(2))
 end subroutine nve_lf

 
 subroutine rot_midstep_impl_lf(engrotsum)
  implicit none
  real(kind=PRC), intent(inout) :: engrotsum
  integer :: i,itq,itype, myi
  real(kind=PRC), parameter :: onefive=real(1.5d0,kind=PRC)
  
! working variables
  real(kind=PRC) :: delx,dely,delz
  real(kind=PRC) :: oqx,oqy,oqz
  real(kind=PRC) :: opx,opy,opz
  real(kind=PRC) :: trx,try,trz
  real(kind=PRC) :: eps,rnorm,myrot(9)
  real(kind=PRC) :: qn0,qn1,qn2,qn3
  real(kind=PRC) :: qn0a,qn1a,qn2a,qn3a
  real(kind=PRC) :: dq0,dq1,dq2,dq3
  

  do myi=1,natms
      i = atmbook(myi)
      itype=ltype(i)
      
      !current rotational matrix 
      call quat_2_rotmat(q0(i),q1(i),q2(i),q3(i),myrot)
      
!     store current angular velocity which are stile at n-1/2
      opx=oxx(i)
      opy=oyy(i)
      opz=ozz(i)
      
!     compute angular velocity increment of one step
      trx=(tqx(i)*myrot(1)+tqy(i)*myrot(4)+tqz(i)*myrot(7))/rotinx(itype)+ &
       (rotiny(itype)-rotinz(itype))*opy*opz/rotinx(itype)
      try=(tqx(i)*myrot(2)+tqy(i)*myrot(5)+tqz(i)*myrot(8))/rotiny(itype)+ &
       (rotinz(itype)-rotinx(itype))*opz*opx/rotiny(itype)
      trz=(tqx(i)*myrot(3)+tqy(i)*myrot(6)+tqz(i)*myrot(9))/rotinz(itype)+ &
       (rotinx(itype)-rotiny(itype))*opx*opy/rotinz(itype)
      
      delx=tstepatm*trx
      dely=tstepatm*try
      delz=tstepatm*trz
      
!     angular velocity at time step n
      opx=oxx(i)+delx*HALF
      opy=oyy(i)+dely*HALF
      opz=ozz(i)+delz*HALF
      
      
!     add its contribution to the total rotational kinetic energy at time step n
      engrotsum=engrotsum+(rotinx(itype)*opx**TWO+ &
       rotiny(itype)*opy**TWO+ &
       rotinz(itype)*opz**TWO)
       
!     compute derivative at n
      dq0=(-q1(i)*opx-q2(i)*opy-q3(i)*opz)*HALF
      dq1=( q0(i)*opx-q3(i)*opy+q2(i)*opz)*HALF
      dq2=( q3(i)*opx+q0(i)*opy-q1(i)*opz)*HALF
      dq3=(-q2(i)*opx+q1(i)*opy+q0(i)*opz)*HALF
      
!     update at t+1/2
      qn0=q0(i)+dq0*tstepatm*HALF
      qn1=q1(i)+dq1*tstepatm*HALF
      qn2=q2(i)+dq2*tstepatm*HALF
      qn3=q3(i)+dq3*tstepatm*HALF
      
      rnorm=ONE/sqrt(qn0**TWO+qn1**TWO+qn2**TWO+qn3**TWO)
      qn0=qn0*rnorm
      qn1=qn1*rnorm
      qn2=qn2*rnorm
      qn3=qn3*rnorm
      
!     angular velocity at timestep n+1/2
      oxx(i)=oxx(i)+delx
      oyy(i)=oyy(i)+dely
      ozz(i)=ozz(i)+delz
      
      
!     store angular velocity at timestep n+1/2
      opx=oxx(i)
      opy=oyy(i)
      opz=ozz(i)
      
!     assign new quaternions
      
!     iteration of new quaternions (lab fixed)
  
      itq=0
      eps=real(1.0d9,kind=PRC)
  
      do while((itq.lt.mxquat).and.(eps.gt.sqquattol))
        itq=itq+1
        
!       compute derivative at n+1/2
        dq0=(-qn1*opx-qn2*opy-qn3*opz)*HALF
        dq1=( qn0*opx-qn3*opy+qn2*opz)*HALF
        dq2=( qn3*opx+qn0*opy-qn1*opz)*HALF
        dq3=(-qn2*opx+qn1*opy+qn0*opz)*HALF
        
!       update at t+1/2
        qn0a=qn0+dq0*tstepatm*HALF
        qn1a=qn1+dq1*tstepatm*HALF
        qn2a=qn2+dq2*tstepatm*HALF
        qn3a=qn3+dq3*tstepatm*HALF
        
        rnorm=ONE/sqrt(qn0a**TWO+qn1a**TWO+qn2a**TWO+qn3a**TWO)
        qn0a=qn0a*rnorm
        qn1a=qn1a*rnorm
        qn2a=qn2a*rnorm
        qn3a=qn3a*rnorm
        
!      convergence test 
        eps=((qn0a-qn0)**TWO+(qn1a-qn1)**TWO+(qn2a-qn2)**TWO+ &
         (qn3a-qn3)**TWO)*tstepatm**TWO
              
        qn0=qn0a
        qn1=qn1a
        qn2=qn2a
        qn3=qn3a
      enddo
      
! store new quaternions
      q0(i)=q0(i)+dq0*tstepatm
      q1(i)=q1(i)+dq1*tstepatm
      q2(i)=q2(i)+dq2*tstepatm
      q3(i)=q3(i)+dq3*tstepatm
      
      rnorm=ONE/sqrt(q0(i)**TWO+q1(i)**TWO+q2(i)**TWO+q3(i)**TWO)
      q0(i)=q0(i)*rnorm
      q1(i)=q1(i)*rnorm
      q2(i)=q2(i)*rnorm
      q3(i)=q3(i)*rnorm
    enddo
 end subroutine rot_midstep_impl_lf

 
 subroutine rot_fullstep_impl_lf(engrotsum)
 
  implicit none
  
  real(kind=PRC), intent(inout) :: engrotsum
  
  integer :: i,itq,itype
  
  real(kind=PRC), parameter :: onefive=real(1.5d0,kind=PRC)
  
! working variables
  real(kind=PRC) :: delx,dely,delz
  real(kind=PRC) :: oqx,oqy,oqz
  real(kind=PRC) :: opx,opy,opz
  real(kind=PRC) :: trx,try,trz
  real(kind=PRC) :: bxx,byy,bzz
  
  real(kind=PRC) :: eps,rnorm,myrot(9)
  
  real(kind=PRC) :: qn0,qn1,qn2,qn3
  real(kind=PRC) :: qn0a,qn1a,qn2a,qn3a
  real(kind=PRC) :: qn0b,qn1b,qn2b,qn3b
  
  do i=1,natms
    
    itype=ltype(i)
    
    !current rotational matrix 
    call quat_2_rotmat(q0(i),q1(i),q2(i),q3(i),myrot)
      
!   store current angular velocity which are stile at n-1/2
    opx=oxx(i)
    opy=oyy(i)
    opz=ozz(i)
    
!   compute angular velocity increment of one step
    
    trx=(tqx(i)*myrot(1)+tqy(i)*myrot(4)+tqz(i)*myrot(7))/rotinx(itype)+ &
     (rotiny(itype)-rotinz(itype))*opy*opz/rotinx(itype)
    try=(tqx(i)*myrot(2)+tqy(i)*myrot(5)+tqz(i)*myrot(8))/rotiny(itype)+ &
     (rotinz(itype)-rotinx(itype))*opz*opx/rotiny(itype)
    trz=(tqx(i)*myrot(3)+tqy(i)*myrot(6)+tqz(i)*myrot(9))/rotinz(itype)+ &
     (rotinx(itype)-rotiny(itype))*opx*opy/rotinz(itype)
    
    delx=tstepatm*trx
    dely=tstepatm*try
    delz=tstepatm*trz
    
!   angular velocity at time step n
    
    opx=oxx(i)+delx*HALF
    opy=oyy(i)+dely*HALF
    opz=ozz(i)+delz*HALF
    
!   angular velocity at time step n+1/2
    
    bxx=oxx(i)+delx
    byy=oyy(i)+dely
    bzz=ozz(i)+delz
    
!   angular velocity at time step n+1  (needed for quat algorithm)
     
    oqx=oxx(i)+delx*onefive
    oqy=oyy(i)+dely*onefive
    oqz=ozz(i)+delz*onefive
    
!   angular velocity at timestep n
    
    oxx(i)=oxx(i)+HALF*delx
    oyy(i)=oyy(i)+HALF*dely
    ozz(i)=ozz(i)+HALF*delz
    
!   add its contribution to the total rotational kinetic energy at time step n
    
    engrotsum=engrotsum+(rotinx(itype)*oxx(i)**TWO+ &
     rotiny(itype)*oyy(i)**TWO+ &
     rotinz(itype)*ozz(i)**TWO)
     
!   restore half step angular velocity
    
    opx=oxx(i)
    opy=oyy(i)
    opz=ozz(i)
    oxx(i)=bxx
    oyy(i)=byy
    ozz(i)=bzz
    
!   assign new quaternions
    
    qn0=q0(i)+(-q1(i)*opx-q2(i)*opy-q3(i)*opz)*tstepatm*HALF
    qn1=q1(i)+( q0(i)*opx-q3(i)*opy+q2(i)*opz)*tstepatm*HALF
    qn2=q2(i)+( q3(i)*opx+q0(i)*opy-q1(i)*opz)*tstepatm*HALF
    qn3=q3(i)+(-q2(i)*opx+q1(i)*opy+q0(i)*opz)*tstepatm*HALF
  
    qn0b=ZERO
    qn1b=ZERO
    qn2b=ZERO
    qn3b=ZERO
      
    itq=0
    eps=real(1.0d9,kind=PRC)
  
    do while((itq.lt.mxquat).and.(eps.gt.sqquattol))
      
      itq=itq+1
      
      qn0a=HALF*(-q1(i)*opx-q2(i)*opy-q3(i)*opz)+ &
       HALF*(-qn1*oqx-qn2*oqy-qn3*oqz)
      qn1a=HALF*(  q0(i)*opx-q3(i)*opy+q2(i)*opz)+ &
       HALF*( qn0*oqx-qn3*oqy+qn2*oqz)
      qn2a=HALF*(  q3(i)*opx+q0(i)*opy-q1(i)*opz)+ &
       HALF*( qn3*oqx+qn0*oqy-qn1*oqz)
      qn3a=HALF*(-q2(i)*opx+q1(i)*opy+q0(i)*opz)+ &
       HALF*(-qn2*oqx+qn1*oqy+qn0*oqz)
      
      qn0=q0(i)+HALF*qn0a*tstepatm
      qn1=q1(i)+HALF*qn1a*tstepatm
      qn2=q2(i)+HALF*qn2a*tstepatm
      qn3=q3(i)+HALF*qn3a*tstepatm
      
      rnorm=ONE/sqrt(qn0**TWO+qn1**TWO+qn2**TWO+qn3**TWO)
      qn0=qn0*rnorm
      qn1=qn1*rnorm
      qn2=qn2*rnorm
      qn3=qn3*rnorm
        
!  convergence test 
       
      eps=((qn0a-qn0b)**TWO+(qn1a-qn1b)**TWO+(qn2a-qn2b)**TWO+ &
       (qn3a-qn3b)**TWO)*tstepatm**TWO
      
      qn0b=qn0a
      qn1b=qn1a
      qn2b=qn2a
      qn3b=qn3a
          
    enddo
      
    q0(i)=qn0
    q1(i)=qn1
    q2(i)=qn2
    q3(i)=qn3
      
    
  enddo
  
  return
  
 end subroutine rot_fullstep_impl_lf
 
 function take_rotversor(idir,qa,qb,qc,qd)
 
!***********************************************************************
!     
!     LBsoft subroutine for determining the versor along the idir
!     direction from the quaternion using x convention for euler angles
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: idir
  
  real(kind=PRC), intent(in) :: qa,qb,qc,qd
  
  real(kind=PRC) :: q0,q1,q2,q3
  
  real(kind=PRC), dimension(3) :: take_rotversor
  
#ifdef SVANBERG
  
  real(kind=PRC), dimension(0:3) :: qversor,qtemp,qtemp2
  
  qtemp(0)=qa
  qtemp(1)=qb
  qtemp(2)=qc
  qtemp(3)=qd
  
  select case(idir)
  case(1)
    
    qversor = (/ ZERO , ONE , ZERO , ZERO /)
    
  case(2)
    
    qversor = (/ ZERO , ZERO , ONE , ZERO /)
    
  case(3)
    
    qversor = (/ ZERO , ZERO , ZERO , ONE /)
    
  case default
    
    !so you want to make me angry!
    take_rotversor=(/ZERO,ZERO,ZERO/)
    return
    
  end select
  
  qtemp2=qtrimult(qtemp,qversor,qconj(qtemp))
  take_rotversor=qtemp2(1:3)

#else
  
  q0=qa
  q1=-qb
  q2=-qc
  q3=-qd
  
  select case(idir)
  case(1)
  
    take_rotversor(1)=q0**TWO+q1**TWO-q2**TWO-q3**TWO
    take_rotversor(2)=TWO*(q1*q2-q0*q3)
    take_rotversor(3)=TWO*(q1*q3+q0*q2)
  
  case(2)
  
    take_rotversor(1)=TWO*(q1*q2+q0*q3)
    take_rotversor(2)=q0**TWO-q1**TWO+q2**TWO-q3**TWO
    take_rotversor(3)=TWO*(q2*q3-q0*q1)
  
  case(3)
  
    take_rotversor(1)=TWO*(q1*q3-q0*q2)
    take_rotversor(2)=TWO*(q2*q3+q0*q1)
    take_rotversor(3)=q0**TWO-q1**TWO-q2**TWO+q3**TWO
  
  case default
  
    !so you want to make me angry!
    take_rotversor=(/ZERO,ZERO,ZERO/)
    
  end select
  
#endif
  
  return
  
 end function take_rotversor
 
 function take_rotversorx(qa,qb,qc,qd)
 
!***********************************************************************
!     
!     LBsoft subroutine for determining the versor along x from
!     the quaternion using x convention for euler angles
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
 
  implicit none
      
  real(kind=PRC), intent(in) :: qa,qb,qc,qd
  
  real(kind=PRC) :: q0,q1,q2,q3
  
  real(kind=PRC), dimension(3) :: take_rotversorx
  
#ifdef SVANBERG
  
  real(kind=PRC), dimension(0:3), parameter :: qversor= &
   (/ ZERO , ONE , ZERO , ZERO /)
  
  real(kind=PRC), dimension(0:3) :: qtemp,qtemp2
  
  qtemp(0)=qa
  qtemp(1)=qb
  qtemp(2)=qc
  qtemp(3)=qd
  
  qtemp2=qtrimult(qtemp,qversor,qconj(qtemp))
  take_rotversorx=qtemp2(1:3)
  
#else
  
  q0=qa
  q1=-qb
  q2=-qc
  q3=-qd
  
  take_rotversorx(1)=q0**TWO+q1**TWO-q2**TWO-q3**TWO
  take_rotversorx(2)=TWO*(q1*q2-q0*q3)
  take_rotversorx(3)=TWO*(q1*q3+q0*q2)
  
#endif
  
  return
  
 end function take_rotversorx
 
 function take_rotversory(qa,qb,qc,qd)
 
!***********************************************************************
!     
!     LBsoft subroutine for determining the versor along y from
!     the quaternion using x convention for euler angles
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
 
  implicit none
      
  real(kind=PRC), intent(in) :: qa,qb,qc,qd
  
  real(kind=PRC) :: q0,q1,q2,q3
  
  real(kind=PRC), dimension(3) :: take_rotversory
  
#ifdef SVANBERG

  real(kind=PRC), dimension(0:3), parameter :: qversor= &
   (/ ZERO , ZERO , ONE , ZERO /)
  
  real(kind=PRC), dimension(0:3) :: qtemp,qtemp2
  
  qtemp(0)=qa
  qtemp(1)=qb
  qtemp(2)=qc
  qtemp(3)=qd
  
  qtemp2=qtrimult(qtemp,qversor,qconj(qtemp))
  take_rotversory=qtemp2(1:3)

#else
  
  q0=qa
  q1=-qb
  q2=-qc
  q3=-qd
  
  take_rotversory(1)=TWO*(q1*q2+q0*q3)
  take_rotversory(2)=q0**TWO-q1**TWO+q2**TWO-q3**TWO
  take_rotversory(3)=TWO*(q2*q3-q0*q1)
  
#endif
  
  return
  
 end function take_rotversory
 
 function take_rotversorz(qa,qb,qc,qd)
 
!***********************************************************************
!     
!     LBsoft subroutine for determining the versor along z from
!     the quaternion using x convention for euler angles
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
 
  implicit none
      
  real(kind=PRC), intent(in) :: qa,qb,qc,qd
  
  real(kind=PRC) :: q0,q1,q2,q3
  
  real(kind=PRC), dimension(3) :: take_rotversorz
  
#ifdef SVANBERG

  real(kind=PRC), dimension(0:3), parameter :: qversor= &
   (/ ZERO , ZERO , ZERO , ONE /)
  
  real(kind=PRC), dimension(0:3) :: qtemp,qtemp2
  
  qtemp(0)=qa
  qtemp(1)=qb
  qtemp(2)=qc
  qtemp(3)=qd
  
  qtemp2=qtrimult(qtemp,qversor,qconj(qtemp))
  take_rotversorz=qtemp2(1:3)

#else
  
  q0=qa
  q1=-qb
  q2=-qc
  q3=-qd
  
  take_rotversorz(1)=TWO*(q1*q3-q0*q2)
  take_rotversorz(2)=TWO*(q2*q3+q0*q1)
  take_rotversorz(3)=q0**TWO-q1**TWO-q2**TWO+q3**TWO
  
#endif
  
  return
  
 end function take_rotversorz
 
 pure function qsum(qsa,qsb)
 
!***********************************************************************
!     
!     LBsoft subroutine to compute the sum of two quaternions
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification December 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in), dimension(0:3) :: qsa,qsb
  
  real(kind=PRC), dimension(0:3) :: qsum
  
  qsum(0:3)=qsa(0:3)+qsb(0:3)
  
  return
  
 end function qsum
 
 pure function qdiff(qsa,qsb)
 
!***********************************************************************
!     
!     LBsoft subroutine to compute the difference of two quaternions
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification December 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in), dimension(0:3) :: qsa,qsb
  
  real(kind=PRC), dimension(0:3) :: qdiff
  
  qdiff(0:3)=qsa(0:3)-qsb(0:3)
  
  return
  
 end function qdiff
 
 pure function qmult(qsa,qsb)
 
!***********************************************************************
!     
!     LBsoft subroutine to compute the product of two quaternions.
!     NOTE: it is not commutative!
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification December 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in), dimension(0:3) :: qsa,qsb
  
  real(kind=PRC), dimension(0:3) :: qmult
  
  qmult(0)=qsa(0)*qsb(0)-qsa(1)*qsb(1)-qsa(2)*qsb(2)-qsa(3)*qsb(3)
  qmult(1:3)=qsa(0)*qsb(1:3)+qsb(0)*qsa(1:3)+cross(qsa(1:3),qsb(1:3))
  
  return
  
 end function qmult
 
 pure function qtrimult(qsa,qsb,qsc)
 
!***********************************************************************
!     
!     LBsoft subroutine to compute the triple product of two
!     quaternions.
!     NOTE: it is not commutative!
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification December 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in), dimension(0:3) :: qsa,qsb,qsc
  
  real(kind=PRC), dimension(0:3) :: qtrimult,qtemps
  
  qtemps=qmult(qsb,qsc)
  qtrimult=qmult(qsa,qtemps)
  
  return
  
 end function qtrimult
 
 pure function qzeroreal(qsa)
 
!***********************************************************************
!     
!     LBsoft subroutine to set to zero the real part of a
!     quaternions.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification December 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in), dimension(0:3) :: qsa
  
  real(kind=PRC), dimension(0:3) :: qzeroreal
  
  qzeroreal(0)=ZERO
  qzeroreal(1:3)=qsa(1:3)
  
  return
  
 end function qzeroreal
 
 pure function qzeroimg(qsa)
 
!***********************************************************************
!     
!     LBsoft subroutine to set to zero the imaginary part of a
!     quaternions.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification December 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in), dimension(0:3) :: qsa
  
  real(kind=PRC), dimension(0:3) :: qzeroimg
  
  qzeroimg(0)=qsa(0)
  qzeroimg(1:3)=ZERO
  
  return
  
 end function qzeroimg
 
 pure function qdot(qsa,qsb)
 
!***********************************************************************
!     
!     LBsoft subroutine to compute the dot product of two quaternions
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification December 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in), dimension(0:3) :: qsa,qsb
  
  real(kind=PRC) :: qdot
  
  qdot=qsa(0)*qsb(0)+qsa(1)*qsb(1)+qsa(2)*qsb(2)+qsa(3)*qsb(3)
  
  return
  
 end function qdot
 
 pure function qmodule(qs)
 
!***********************************************************************
!     
!     LBsoft subroutine to compute the module of a quaternion
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification December 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in), dimension(0:3) :: qs
  
  real(kind=PRC) :: qmodule
  
  qmodule = sqrt(qs(0)**TWO+qs(1)**TWO+qs(2)**TWO+qs(3)**TWO)
  
  return
  
 end function qmodule
 
 pure function qang(ang,vect)
 
!***********************************************************************
!     
!     LBsoft subroutine to compute the quaternion associated to
!     a rotation of angle [ang] around the unit vector [vect]
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification December 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: ang
  real(kind=PRC), intent(in), dimension(1:3) :: vect
  
  real(kind=PRC), dimension(0:3) :: qang
  
  real(kind=PRC) :: lmod,uvect(1:3)
  
  lmod = sqrt(vect(1)**TWO+vect(2)**TWO+vect(3)**TWO)
  uvect(1:3)=vect(1:3)/lmod
  
  qang(0)=cos(ang*HALF)
  qang(1:3)=sin(ang*HALF)*uvect(1:3)
  
  lmod = sqrt(qang(0)**TWO+qang(1)**TWO+qang(2)**TWO+qang(3)**TWO)
  
  qang(0:3)=qang(0:3)/lmod
  
  return
  
 end function qang
 
 subroutine qnorm(qs)
 
!***********************************************************************
!     
!     LBsoft subroutine to normalize a quaternion
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification December 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), dimension(0:3) :: qs
  
  real(kind=PRC) :: lmod
  
  lmod = sqrt(qs(0)**TWO+qs(1)**TWO+qs(2)**TWO+qs(3)**TWO)
  
  qs(0:3)=qs(0:3)/lmod
  
  return
  
 end subroutine qnorm
 
 pure function qconj(qs)
 
!***********************************************************************
!     
!     LBsoft subroutine to compute the conjugate of a quaternion
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification December 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in), dimension(0:3) :: qs
  
  real(kind=PRC), dimension(0:3) :: qconj
  
  qconj(0)=qs(0)
  qconj(1:3)=-qs(1:3)
  
  return
  
 end function qconj
 
 pure function  qinv(qs)
 
!***********************************************************************
!     
!     LBsoft subroutine to compute the inverse of a quaternion
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification December 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in), dimension(0:3) :: qs
  
  real(kind=PRC), dimension(0:3) :: qinv
  
  real(kind=PRC) :: lmod
  
  lmod = sqrt(qs(0)**TWO+qs(1)**TWO+qs(2)**TWO+qs(3)**TWO)
  
  qinv(0)=qs(0)/lmod
  qinv(1:3)=-qs(1:3)/lmod
  
  return
  
 end function qinv
 
 subroutine eul2q(phis,psis,thetas,qs)
 
!***********************************************************************
!     
!     LBsoft subroutine to compute the quaternion associated to
!     the Euler angles
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification December 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: phis,psis,thetas
  real(kind=PRC), intent(out), dimension(0:3) :: qs
  
  real(kind=PRC) :: cy,sy,cp,sp,cr,sr
  
  cy = cos(phis * HALF)
  sy = sin(phis * HALF)
  cp = cos(psis * HALF)
  sp = sin(psis * HALF)
  cr = cos(thetas * HALF)
  sr = sin(thetas * HALF)
  
  qs(0) = cy * cp * cr + sy * sp * sr
  qs(1) = cy * cp * sr - sy * sp * cr
  qs(2) = sy * cp * sr + cy * sp * cr
  qs(3) = sy * cp * cr - cy * sp * sr
  
  return
  
 end subroutine eul2q
 
 subroutine q2eul(qs,phis,psis,thetas)
 
!***********************************************************************
!     
!     LBsoft subroutine to compute the three Euler angles associated to
!     a quaternion
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification December 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in), dimension(0:3) :: qs
  real(kind=PRC), intent(out) :: phis,psis,thetas
  
  real(kind=PRC) :: sinr_cosp,cosr_cosp,siny_cosp,cosy_cosp,sinp
  
  sinr_cosp = TWO * (qs(0) * qs(1) + qs(2) * qs(3))
  cosr_cosp = ONE- TWO * (qs(1)**TWO + qs(2)**TWO)
  thetas = atan2(sinr_cosp, cosr_cosp)

! psis (y-axis rotation)
  sinp = TWO * (qs(0) * qs(2) - qs(3) * qs(1))
  if(abs(sinp)>=ONE)then
    psis = sign(Pi*HALF, sinp) 
  else
    psis = asin(sinp)
  endif

! phis (z-axis rotation)
  siny_cosp = TWO * (qs(0) * qs(3) + qs(1) * qs(2))
  cosy_cosp = ONE - TWO * (qs(2)**TWO + qs(3)**TWO)
  phis = atan2(siny_cosp, cosy_cosp)
  
  return
  
 end subroutine q2eul
 
 subroutine q2mat(qs,mat)
 
!***********************************************************************
!     
!     LBsoft subroutine to compute the rotation matrix associated to
!     a quaternion
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification December 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in), dimension(0:3) :: qs
  real(kind=PRC), intent(out), dimension(3,3) :: mat
  
  mat(1,1) = qs(0)**TWO+qs(1)**TWO-qs(2)**TWO-qs(3)**TWO
  mat(2,1) = TWO*(qs(1)*qs(2)+qs(0)*qs(3))
  mat(3,1) = TWO*(qs(1)*qs(3)-qs(0)*qs(2))
  mat(1,2) = TWO*(qs(1)*qs(2)-qs(0)*qs(3))
  mat(2,2) = qs(0)**TWO-qs(1)**TWO+qs(2)**TWO-qs(3)**TWO
  mat(3,2) = TWO*(qs(2)*qs(3)+qs(0)*qs(1))
  mat(1,3) = TWO*(qs(1)*qs(3)+qs(0)*qs(2))
  mat(2,3) = TWO*(qs(2)*qs(3)-qs(0)*qs(1))
  mat(3,3) = qs(0)**TWO-qs(1)**TWO-qs(2)**TWO+qs(3)**TWO
  
  return
  
 end subroutine q2mat
 
 subroutine integrate_particles_vv(isw,nstepsub)
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the particle integration by
!     velocity Verlet
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: isw,nstepsub
  
  select case(keyint)
  case (1) 
    call nve_vv(isw)
  case default
    call nve_vv(isw)
  end select
  
  return
  
 end subroutine integrate_particles_vv
 
 subroutine nve_vv(isw)

!***********************************************************************
!     
!     LBsoft subroutine for integrating newtonian EOM by 
!     velocity Verlet 
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************

  implicit none
  
  integer, intent(in) :: isw
  
  integer :: i
  
! update velocities for first and second stages
  forall(i=1:natms)
    vxx(i)=vxx(i)+HALF*tstepatm*fxx(i)*rmass(ltype(i))
    vyy(i)=vyy(i)+HALF*tstepatm*fyy(i)*rmass(ltype(i))
    vzz(i)=vzz(i)+HALF*tstepatm*fzz(i)*rmass(ltype(i))
  end forall
  
  select case(isw)
  case(1)
  
!   update positions
    forall(i=1:natms)
      xxx(i)=xxx(i)+tstepatm*vxx(i)
      yyy(i)=yyy(i)+tstepatm*vyy(i)
      zzz(i)=zzz(i)+tstepatm*vzz(i)
    end forall
    
!   periodic boundary condition
    call pbc_images_centered(imcon,natms,cell,cx,cy,cz,xxx,yyy,zzz)
  
  case(2)
  
!   calculate kinetic energy    
    call getkin(engke)
    
  end select
      
  return
      
 end subroutine nve_vv
 
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
  integer :: i, myi

  
  engkes=ZERO
  
  do myi=1,natms
    i = atmbook(myi)
    engkes=engkes+weight(ltype(i))*(vxx(i)**TWO+vyy(i)**TWO+vzz(i)**TWO)
  enddo
  
  engkes=HALF*engkes
 end subroutine getkin

 
 subroutine rotmat_2_quat(rot,q0s,q1s,q2s,q3s)
 
!***********************************************************************
!     
!     LBsoft subroutine for determining the quaternion associated
!     to a rotational matrix within a given tollerance value
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
    integer :: i,myi

    mynamefile=repeat(' ',120)
    mynamefile=trim(filenam)//write_fmtnumb(itersub)//'_'//write_fmtnumb(idrank)//'.xyz'
  
    open(unit=23,file=trim(mynamefile),status='replace')
    write(23,'(i8)')natms
    write(23,'(i8)')itersub

    do myi=1,natms
        i = atmbook(myi)
        write(23,'(a8,6f16.8)')'C       ', xxx(i),yyy(i),zzz(i),vxx(i),vyy(i),vzz(i)
    enddo
    close(23)

 end subroutine print_all_particles
 

 subroutine spherical_template(rdimsub,outnum,outlist,outdist, &
  outnumdead,outlistdead)
 
!***********************************************************************
!     
!     LBsoft subroutine for creating the spherical template (grid approx of a sphere)
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification Jan 2019
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: rdimsub
  integer, intent(out) :: outnum,outnumdead
  integer, allocatable, dimension(:,:) :: outlist,outlistdead
  real(kind=PRC), allocatable, dimension(:) :: outdist
  
  integer :: i,j,k,l,m,rmax,rmin,ishift,jshift,kshift
  integer :: ioppshift,joppshift,koppshift
  
  real(kind=PRC) :: vec(3),rdist,sqrcut
  integer, allocatable, dimension(:,:,:) :: issub
  
  integer, allocatable, dimension(:) :: icount
  
  sqrcut=rdimsub**TWO
  
  rmax=floor(rdimsub)+1
  rmin=floor(rdimsub)
  allocate(issub(-rmax:rmax,-rmax:rmax,-rmax:rmax))
  issub(-rmax:rmax,-rmax:rmax,-rmax:rmax)=1
  
  do k=-rmax,rmax
    vec(3)=real(k,kind=PRC)
    do j=-rmax,rmax
      vec(2)=real(j,kind=PRC)
      do i=-rmax,rmax
        vec(1)=real(i,kind=PRC)
        rdist=vec(1)**TWO+vec(2)**TWO+vec(3)**TWO
        if(rdist<=sqrcut)then
          issub(i,j,k)=3
          if(i>rmin .or. i<-rmin .or. &
           j>rmin .or. j<-rmin .or. &
           k>rmin .or. k<-rmin)then
            call error(29)
          endif
        endif
      enddo
    enddo
  enddo
  
  do k=-rmin,rmin
    do j=-rmin,rmin
      do i=-rmin,rmin
        if(issub(i,j,k)==3)then
          do l=1,links
            if(issub(i+ex(l),j+ey(l),k+ez(l))==1)then
              issub(i,j,k)=2
              exit
            endif
          enddo
        endif
      enddo
    enddo
  enddo
  
  l=0
  do k=-rmax,rmax
    do j=-rmax,rmax
      do i=-rmax,rmax
        if(issub(i,j,k)==2)then
          l=l+1
        endif
      enddo
    enddo
  enddo
  
  outnum=l
  if(allocated(outlist))deallocate(outlist)
  allocate(outlist(3,outnum))
  if(allocated(outdist))deallocate(outdist)
  allocate(outdist(outnum))
  l=0
  do k=-rmax,rmax
    vec(3)=real(k,kind=PRC)
    do j=-rmax,rmax
      vec(2)=real(j,kind=PRC)
      do i=-rmax,rmax
        vec(1)=real(i,kind=PRC)
        if(issub(i,j,k)==2)then
          rdist=vec(1)**TWO+vec(2)**TWO+vec(3)**TWO
          l=l+1
          outlist(1,l)=i
          outlist(2,l)=j
          outlist(3,l)=k
          outdist(l)=sqrt(rdist)
        endif
      enddo
    enddo
  enddo


  l=0
  do k=-rmax,rmax
    do j=-rmax,rmax
      do i=-rmax,rmax
        if(issub(i,j,k)==3)then
          l=l+1
        endif
      enddo
    enddo
  enddo
  
  outnumdead=l
  if(allocated(outlistdead))deallocate(outlistdead)
  allocate(outlistdead(3,outnumdead))
  
  l=0
  do k=-rmax,rmax
    do j=-rmax,rmax
      do i=-rmax,rmax
        if(issub(i,j,k)==3)then
          l=l+1
          outlistdead(1,l)=i
          outlistdead(2,l)=j
          outlistdead(3,l)=k
        endif
      enddo
    enddo
  enddo
  
  allocate(icount(outnum))
  icount(1:outnum)=0
  do m=1,outnum
    i=outlist(1,m)
    j=outlist(2,m)
    k=outlist(3,m)
    do l=1,links
      ishift=i+ex(l)
      jshift=j+ey(l)
      kshift=k+ez(l)
      ioppshift=i+ex(opp(l))
      joppshift=j+ey(opp(l))
      koppshift=k+ez(opp(l))
      if(issub(ishift,jshift,kshift)/=1)cycle
      if(issub(ioppshift,joppshift,koppshift)==1)then
        icount(m)=icount(m)+1
      endif
    enddo
  enddo

 end subroutine spherical_template
 
 subroutine restore_particles
  implicit none
  integer :: i, myi, iatm, ids, num_ext
  logical(kind=1), dimension(natms_tot) :: mine


!  if (mxrank==1) return

  mine(:) = .false.
  do myi=1,natms
      i = atmbook(myi)
      mine(i) = .true.
  enddo

  forall(i=1:natms_tot, .not. mine(i))
    ! Erase atoms by others
    xxx(i) = 0
    yyy(i) = 0
    zzz(i) = 0
    vxx(i) = 0
    vyy(i) = 0
    vzz(i) = 0

    q0(i) = 0
    q1(i) = 0
    q2(i) = 0
    q3(i) = 0

    xxo(i) = 0
    yyo(i) = 0
    zzo(i) = 0
    vxo(i) = 0
    vyo(i) = 0
    vzo(i) = 0
  end forall

  call sum_world_farr(xxx,natms_tot)
  call sum_world_farr(yyy,natms_tot)
  call sum_world_farr(zzz,natms_tot)
  call sum_world_farr(vxx,natms_tot)
  call sum_world_farr(vyy,natms_tot)
  call sum_world_farr(vzz,natms_tot)

  call sum_world_farr(xxo,natms_tot)
  call sum_world_farr(yyo,natms_tot)
  call sum_world_farr(zzo,natms_tot)
  call sum_world_farr(vxo,natms_tot)
  call sum_world_farr(vyo,natms_tot)
  call sum_world_farr(vzo,natms_tot)

  call sum_world_farr(q0,natms_tot)
  call sum_world_farr(q1,natms_tot)
  call sum_world_farr(q2,natms_tot)
  call sum_world_farr(q3,natms_tot)


  ! Count my atoms, put them in atmbook list
  ! Also count halo atoms, at end of atmbook
  natms = 0; num_ext = 0

  do i=1,natms_tot
#ifndef DEOWERN
    ids=ownernfind_arr(nint(xxx(i)),nint(yyy(i)),nint(zzz(i)), &
        nx,ny,nz,nbuff,ownern)
#else
    ids=ownernfind(nint(xxx(i)),nint(yyy(i)),nint(zzz(i)), &
        mxrank,gminx,gmaxx,gminy,gmaxy,gminz,gmaxz)
#endif
    if (idrank==ids) then
        natms = natms + 1
        atmbook(natms) = i
    else
        call checkIfExtAtom(i, num_ext)
    endif
  enddo

  ! 1) First parallelization with all atoms on all processes FB
  do iatm=1, num_ext
    atmbook(natms+iatm) = atmbook(mxatms - num_ext + iatm)
  enddo
  natms_ext= natms + num_ext

!  do myi=1,natms_ext
!      i = atmbook(myi)
!      write (6,*) __FILE__,__LINE__, "i=", i, myi<=natms, xxx(i),yyy(i),zzz(i)
!  enddo

 end subroutine restore_particles

 end module particles_mod
