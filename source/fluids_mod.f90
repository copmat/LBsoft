
#include <default_macro.h>
 module fluids_mod
 
!***********************************************************************
!     
!     LBsoft module containing variable and subroutines of the
!     fluid components
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
 use version_mod,    only : idrank,or_world_larr
 use error_mod
 use utility_mod, only : Pi,modulvec,cross,dot,gauss,ibuffservice, &
                   allocate_array_ibuffservice,buffservice, &
                   allocate_array_buffservice,lbuffservice, &
                   allocate_array_lbuffservice, &
                   buffservice3d,allocate_array_buffservice3d
 
 implicit none
 
 private
 
 integer, save, protected, public :: nx=0
 integer, save, protected, public :: ny=0
 integer, save, protected, public :: nz=0
 
 integer, save, protected, public :: idistselect=0
 
 integer, save, protected, public :: ibctype=0
 
 logical, save, protected, public :: lalloc_hvars=.false.
 logical, save, protected, public :: lalloc_pops=.false.
 
 logical, save, protected, public :: lforce_add=.false.
 
 real(kind=PRC), save, protected, public :: t_LB = ONE
 
 real(kind=PRC), save, protected, public :: beta   = ZERO
 real(kind=PRC), save, protected, public :: akl    = ZERO
 real(kind=PRC), save, protected, public :: awall  = ZERO
 real(kind=PRC), save, protected, public :: tauR   = ZERO
 real(kind=PRC), save, protected, public :: tauB   = ZERO
 real(kind=PRC), save, protected, public :: omegaR = ZERO
 real(kind=PRC), save, protected, public :: omegaB = ZERO
 real(kind=PRC), save, protected, public :: viscR  = ZERO
 real(kind=PRC), save, protected, public :: viscB  = ZERO
 real(kind=PRC), save, protected, public :: meanR  = ZERO
 real(kind=PRC), save, protected, public :: meanB  = ZERO
 real(kind=PRC), save, protected, public :: stdevR  = ZERO
 real(kind=PRC), save, protected, public :: stdevB  = ZERO
 real(kind=PRC), save, protected, public :: initial_u = ZERO
 real(kind=PRC), save, protected, public :: initial_v = ZERO
 real(kind=PRC), save, protected, public :: initial_w = ZERO
 real(kind=PRC), save, protected, public :: ext_fu = ZERO
 real(kind=PRC), save, protected, public :: ext_fv = ZERO
 real(kind=PRC), save, protected, public :: ext_fw = ZERO
 
 real(kind=PRC), save, protected, public :: pair_SC = ZERO
 
 
 integer, save, protected, public, allocatable, dimension(:,:,:) :: isfluid
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: rhoR
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: rhoB
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: u
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: v
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: w
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: omega
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: fuR
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: fvR
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: fwR
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: fuB
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: fvB
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: fwB
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: gradpsixR
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: gradpsiyR
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: gradpsizR
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: gradpsixB
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: gradpsiyB
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: gradpsizB
 
#if LATTICE==319
 
 integer, parameter, public :: latt_dim=3
 integer, parameter, public :: links=18
 
 integer, parameter, public :: nbuff=2
 
 real(kind=PRC), parameter, public :: cssq = ( ONE / THREE )
 
 real(kind=PRC), parameter :: p0 = ( TWELVE / THIRTYSIX )
 real(kind=PRC), parameter :: p1 = ( TWO/ THIRTYSIX )
 real(kind=PRC), parameter :: p2 = ( ONE / THIRTYSIX )
 real(kind=PRC), dimension(0:links), parameter, public :: &
  p = (/p0,p1,p1,p1,p1,p1,p1,p2,p2,p2,p2,p2,p2,p2,p2,p2,p2,p2,p2/)
 
 !lattice vectors
 integer, dimension(0:links), parameter, public :: &
  ex = (/0,1,-1,0,0,0,0,1,-1,-1, 1, 1,-1,-1, 1, 0, 0, 0, 0/)
   !     0,1, 2,3,4,5,6,7, 8, 9,10,11,12,13,14,15,16,17,18
 integer, dimension(0:links), parameter, public :: &
  ey = (/0,0,0,1,-1,0,0,1,-1,1,-1, 0, 0, 0, 0, 1,-1,-1, 1/)
 integer, dimension(0:links), parameter, public :: &
  ez = (/0,0,0,0,0,1,-1,0,0,0, 0, 1,-1, 1,-1, 1,-1, 1,-1/)
 integer, dimension(0:links), parameter, public :: &
  opp =(/0,2,1,4,3,6,5, 8,7,10,9,12,11,14,13,16,15,18,17/)
  
 real(kind=PRC), dimension(0:links), parameter, public :: &
  dex = real(ex,kind=PRC)
 real(kind=PRC), dimension(0:links), parameter, public :: &
  dey = real(ey,kind=PRC)
 real(kind=PRC), dimension(0:links), parameter, public :: &
  dez = real(ez,kind=PRC)
 
#endif 



#if LATTICE==319

 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: & 
  f00R,f01R,f02R,f03R,f04R,f05R,f06R,f07R,f08R
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: &
  f00B,f01B,f02B,f03B,f04B,f05B,f06B,f07B,f08B
 
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: & 
  f09R,f10R,f11R,f12R,f13R,f14R,f15R,f16R,f17R,f18R
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: &
  f09B,f10B,f11B,f12B,f13B,f14B,f15B,f16B,f17B,f18B

#endif
 
 public :: set_random_dens_fluids
 public :: set_initial_dens_fluids
 public :: set_initial_vel_fluids
 public :: set_mean_value_dens_fluids
 public :: set_stdev_value_dens_fluids
 public :: set_initial_dist_type
 public :: set_initial_dim_box
 public :: set_mean_value_vel_fluids
 public :: allocate_fluids
 public :: initialize_fluids
 public :: set_boundary_conditions_type
 public :: driver_bc_pops
 public :: initialize_fluid_force
 public :: set_value_ext_force_fluids
 public :: compute_fluid_force_sc
 public :: collision_fluids
 public :: driver_bc_hvars
 
 contains
 
 subroutine allocate_fluids
 
!***********************************************************************
!     
!     LBsoft subroutine for allocating arrays which 
!     describe the fluids
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  integer, parameter :: nistatmax=100
  integer, dimension(nistatmax) :: istat
  integer :: myzero,mynx,myny,mynz
  
  integer :: i
  logical, dimension(1) :: ltest=.false.
  
  myzero=1-nbuff
  mynx=nx+nbuff
  myny=ny+nbuff
  mynz=nz+nbuff
  
  istat=0
  
  allocate(rhoR(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(1))
  allocate(rhoB(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(2))
  
  allocate(u(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(3))
  allocate(v(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(4))
  allocate(w(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(5))
  
  allocate(fuR(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(6))
  allocate(fvR(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(7))
  allocate(fwR(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(8))
  
  allocate(fuB(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(9))
  allocate(fvB(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(10))
  allocate(fwB(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(11))
  
  allocate(gradpsixR(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(12))
  allocate(gradpsiyR(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(13))
  allocate(gradpsizR(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(14))
  
  allocate(gradpsixB(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(15))
  allocate(gradpsiyB(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(16))
  allocate(gradpsizB(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(17))
  
  allocate(omega(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(18))
  
  allocate(f00R(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(30))
  allocate(f01R(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(31))
  allocate(f02R(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(32))
  allocate(f03R(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(33))
  allocate(f04R(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(34))
  allocate(f05R(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(35))
  allocate(f06R(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(36))
  allocate(f07R(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(37))
  allocate(f08R(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(38))
  allocate(f09R(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(39))
  allocate(f10R(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(40))
  allocate(f11R(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(41))
  allocate(f12R(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(42))
  allocate(f13R(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(43))
  allocate(f14R(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(44))
  allocate(f15R(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(45))
  allocate(f16R(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(46))
  allocate(f17R(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(47))
  allocate(f18R(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(48))
  
  allocate(f00B(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(60))
  allocate(f01B(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(61))
  allocate(f02B(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(62))
  allocate(f03B(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(63))
  allocate(f04B(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(64))
  allocate(f05B(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(65))
  allocate(f06B(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(66))
  allocate(f07B(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(67))
  allocate(f08B(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(68))
  allocate(f09B(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(69))
  allocate(f10B(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(60))
  allocate(f11B(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(61))
  allocate(f12B(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(62))
  allocate(f13B(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(63))
  allocate(f14B(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(64))
  allocate(f15B(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(65))
  allocate(f16B(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(66))
  allocate(f17B(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(67))
  allocate(f18B(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(68))
  
  ltest=.false.
  if(any(istat.ne.0))then
    do i=1,nistatmax
      if(istat(i).ne.0)exit
    enddo
    call warning(2,dble(i))
    ltest=.true.
  endif
  
  call or_world_larr(ltest,1)
  if(ltest(1))call error(4)
  
  call allocate_array_buffservice3d(myzero,mynx,myzero,myny,myzero,mynz)
  
  return
  
 end subroutine allocate_fluids
 
 subroutine set_boundary_conditions_type(itemp1,itemp2,itemp3)
 
  implicit none
  
  integer, intent(in) :: itemp1,itemp2,itemp3
  
  integer, dimension(0:1,0:1,0:1) :: iselbct = &
   reshape((/0,1,2,3,4,5,6,7/),(/2,2,2/))
  
  if( itemp1 .ne. 0 .and. itemp1 .ne. 1)then
    call error(9)
  endif
  
  if( itemp2 .ne. 0 .and. itemp2 .ne. 1)then
    call error(9)
  endif
  
  if( itemp3 .ne. 0 .and. itemp3 .ne. 1)then
    call error(9)
  endif
  
  ibctype=iselbct(itemp1,itemp2,itemp3)
  
  return
  
 end subroutine set_boundary_conditions_type
 
 subroutine initialize_fluids
 
!***********************************************************************
!     
!     LBsoft subroutine for initialize the fluids and read the
!     restart file if requested
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer :: i,j,k
  
  select case(idistselect)
  case(1)
    call set_random_dens_fluids
  case default
    call set_initial_dens_fluids
  end select
  
  call set_initial_vel_fluids
  
  call driver_bc_hvars
  
  call set_initial_pop_fluids
  
  call driver_bc_pops
  
  !restart to be added
  
  return
  
 end subroutine initialize_fluids
 
 subroutine initialize_fluid_force
 
!***********************************************************************
!     
!     LBsoft subroutine for initialize the force fluids
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer :: i,j,k
  
  fuR(1:nx,1:ny,1:nz)=ext_fu
  fvR(1:nx,1:ny,1:nz)=ext_fv
  fwR(1:nx,1:ny,1:nz)=ext_fw
  
  fuB(1:nx,1:ny,1:nz)=ext_fu
  fvB(1:nx,1:ny,1:nz)=ext_fv
  fwB(1:nx,1:ny,1:nz)=ext_fw
  
  return
 
 end subroutine initialize_fluid_force
 
 subroutine set_random_dens_fluids
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the initial density of fluids 
!     following a random gaussian distribution
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  use utility_mod, only: gauss
  
  implicit none
  
  rhoR(1:nx,1:ny,1:nz)=meanR+stdevR*gauss()
  rhoB(1:nx,1:ny,1:nz)=meanB+stdevB*gauss()
  
  return
  
 end subroutine set_random_dens_fluids
 
 subroutine set_initial_dens_fluids
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the initial density of fluids 
!     equal to a given value
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  rhoR(1:nx,1:ny,1:nz)=meanR
  rhoB(1:nx,1:ny,1:nz)=meanB
  
  return
  
 end subroutine set_initial_dens_fluids
 
 subroutine set_initial_vel_fluids
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the initial velocity of fluids 
!     equal to given values
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  u(1:nx,1:ny,1:nz)=initial_u
  v(1:nx,1:ny,1:nz)=initial_v
  w(1:nx,1:ny,1:nz)=initial_w
  
  return
  
 end subroutine set_initial_vel_fluids
 
 subroutine set_initial_pop_fluids
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the initial populations of fluids 
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  integer :: i,j,k
  
  !red fluid
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f00R(i,j,k)=equil_pop00(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f01R(i,j,k)=equil_pop01(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f02R(i,j,k)=equil_pop02(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f03R(i,j,k)=equil_pop03(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f04R(i,j,k)=equil_pop04(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f05R(i,j,k)=equil_pop05(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f06R(i,j,k)=equil_pop06(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f07R(i,j,k)=equil_pop07(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f08R(i,j,k)=equil_pop08(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f09R(i,j,k)=equil_pop09(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f10R(i,j,k)=equil_pop10(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f11R(i,j,k)=equil_pop11(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f12R(i,j,k)=equil_pop12(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f13R(i,j,k)=equil_pop13(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f14R(i,j,k)=equil_pop14(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f15R(i,j,k)=equil_pop15(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f16R(i,j,k)=equil_pop16(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f17R(i,j,k)=equil_pop17(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f18R(i,j,k)=equil_pop18(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))
  end forall
  
  !blue fluid
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f00B(i,j,k)=equil_pop00(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f01B(i,j,k)=equil_pop01(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f02B(i,j,k)=equil_pop02(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f03B(i,j,k)=equil_pop03(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f04B(i,j,k)=equil_pop04(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f05B(i,j,k)=equil_pop05(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f06B(i,j,k)=equil_pop06(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f07B(i,j,k)=equil_pop07(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f08B(i,j,k)=equil_pop08(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f09B(i,j,k)=equil_pop09(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f10B(i,j,k)=equil_pop10(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f11B(i,j,k)=equil_pop11(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f12B(i,j,k)=equil_pop12(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f13B(i,j,k)=equil_pop13(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f14B(i,j,k)=equil_pop14(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f15B(i,j,k)=equil_pop15(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f16B(i,j,k)=equil_pop16(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f17B(i,j,k)=equil_pop17(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f18B(i,j,k)=equil_pop18(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))
  end forall
  
  return
  
 end subroutine set_initial_pop_fluids
 
 subroutine set_mean_value_dens_fluids(dtemp1,dtemp2)
 
!***********************************************************************
!     
!     LBsoft subroutine for change the mean density of fluids 
!     inside this module
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC), intent(in) :: dtemp1,dtemp2
  
  meanR = dtemp1
  meanB = dtemp2
  
  return
  
 end subroutine set_mean_value_dens_fluids
 
 subroutine set_stdev_value_dens_fluids(dtemp1,dtemp2)
 
!***********************************************************************
!     
!     LBsoft subroutine for change the standard devviation 
!     in the density of fluids inside this module
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC), intent(in) :: dtemp1,dtemp2
  
  stdevR = dtemp1
  stdevB = dtemp2
  
  return
  
 end subroutine set_stdev_value_dens_fluids
 
 subroutine set_initial_value_vel_fluids(dtemp1,dtemp2,dtemp3)
 
!***********************************************************************
!     
!     LBsoft subroutine for change the initial velocity values of fluids 
!     inside this module
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC), intent(in) :: dtemp1,dtemp2,dtemp3
  
  initial_u = dtemp1
  initial_v = dtemp2
  initial_w = dtemp3
  
  return
  
 end subroutine set_initial_value_vel_fluids
 
 subroutine set_value_ext_force_fluids(dtemp1,dtemp2,dtemp3)
 
!***********************************************************************
!     
!     LBsoft subroutine for set the external force values acting
!     on the fluids 
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC), intent(in) :: dtemp1,dtemp2,dtemp3
  
  ext_fu = dtemp1
  ext_fv = dtemp2
  ext_fw = dtemp3
  
  return
  
 end subroutine set_value_ext_force_fluids
 
 subroutine set_initial_dist_type(itemp1)
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the initial distribution type for
!     the density fluid initialization
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: itemp1
  
  idistselect = itemp1
  
  return
  
 end subroutine set_initial_dist_type
 
 subroutine set_initial_dim_box(itemp1,itemp2,itemp3)
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the initial dimensions of the
!     simulation box
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: itemp1,itemp2,itemp3
  
  nx = itemp1
  ny = itemp2
  nz = itemp3
  
  return
  
 end subroutine set_initial_dim_box
 
 subroutine set_mean_value_vel_fluids(dtemp1,dtemp2,dtemp3)
 
!***********************************************************************
!     
!     LBsoft subroutine for change the mean velocity of fluids 
!     inside this module
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC), intent(in) :: dtemp1,dtemp2,dtemp3
  
  initial_u = dtemp1
  initial_v = dtemp2
  initial_w = dtemp3
  
  return
  
 end subroutine set_mean_value_vel_fluids
 
 subroutine driver_bc_hvars
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the boundary conditions
!     to hydrodynamic variables
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  select case(ibctype)
  case(7)
    call apply_pbc_hvars
  case default
    call apply_pbc_hvars
  end select
  
  return
  
 end subroutine driver_bc_hvars
 
 subroutine apply_pbc_hvars
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the full periodic boundary 
!     conditions to hydrodynamic variables
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  !sides 6
  
  !apply pbc at north 1
  call apply_pbc_north(rhoR)
  call apply_pbc_north(rhoB)
  
  call apply_pbc_north(u)
  call apply_pbc_north(v)
  call apply_pbc_north(w)
  
  !apply pbc at south 2
  call apply_pbc_south(rhoR)
  call apply_pbc_south(rhoB)
  
  call apply_pbc_south(u)
  call apply_pbc_south(v)
  call apply_pbc_south(w)
  
  !apply pbc at east 3
  call apply_pbc_east(rhoR)
  call apply_pbc_east(rhoB)
  
  call apply_pbc_east(u)
  call apply_pbc_east(v)
  call apply_pbc_east(w)
  
  !apply pbc at west 4
  call apply_pbc_west(rhoR)
  call apply_pbc_west(rhoB)
  
  call apply_pbc_west(u)
  call apply_pbc_west(v)
  call apply_pbc_west(w)
  
  !apply pbc at front 5
  call apply_pbc_front(rhoR)
  call apply_pbc_front(rhoB)
  
  call apply_pbc_front(u)
  call apply_pbc_front(v)
  call apply_pbc_front(w)
  
  !apply pbc at rear 6
  call apply_pbc_rear(rhoR)
  call apply_pbc_rear(rhoB)
  
  call apply_pbc_rear(u)
  call apply_pbc_rear(v)
  call apply_pbc_rear(w)
  
  !edges 12
  
  !apply pbc at front east 1
  call apply_pbc_edge_front_east(rhoR)
  call apply_pbc_edge_front_east(rhoB)
  
  call apply_pbc_edge_front_east(u)
  call apply_pbc_edge_front_east(v)
  call apply_pbc_edge_front_east(w)
  
  !apply pbc at front west 2
  call apply_pbc_edge_front_west(rhoR)
  call apply_pbc_edge_front_west(rhoB)
  
  call apply_pbc_edge_front_west(u)
  call apply_pbc_edge_front_west(v)
  call apply_pbc_edge_front_west(w)
  
  !apply pbc at north east 3
  call apply_pbc_edge_north_east(rhoR)
  call apply_pbc_edge_north_east(rhoB)
  
  call apply_pbc_edge_north_east(u)
  call apply_pbc_edge_north_east(v)
  call apply_pbc_edge_north_east(w)
  
  !apply pbc at north front 4
  call apply_pbc_edge_north_front(rhoR)
  call apply_pbc_edge_north_front(rhoB)
  
  call apply_pbc_edge_north_front(u)
  call apply_pbc_edge_north_front(v)
  call apply_pbc_edge_north_front(w)
  
  !apply pbc at north rear 5
  call apply_pbc_edge_north_rear(rhoR)
  call apply_pbc_edge_north_rear(rhoB)
  
  call apply_pbc_edge_north_rear(u)
  call apply_pbc_edge_north_rear(v)
  call apply_pbc_edge_north_rear(w)
  
  !apply pbc at north west 6
  call apply_pbc_edge_north_west(rhoR)
  call apply_pbc_edge_north_west(rhoB)
  
  call apply_pbc_edge_north_west(u)
  call apply_pbc_edge_north_west(v)
  call apply_pbc_edge_north_west(w)
  
  !apply pbc at rear east 7
  call apply_pbc_edge_rear_east(rhoR)
  call apply_pbc_edge_rear_east(rhoB)
  
  call apply_pbc_edge_rear_east(u)
  call apply_pbc_edge_rear_east(v)
  call apply_pbc_edge_rear_east(w)
  
  !apply pbc at rear west 8
  call apply_pbc_edge_rear_west(rhoR)
  call apply_pbc_edge_rear_west(rhoB)
  
  call apply_pbc_edge_rear_west(u)
  call apply_pbc_edge_rear_west(v)
  call apply_pbc_edge_rear_west(w)
  
  !apply pbc at south east 9
  call apply_pbc_edge_south_east(rhoR)
  call apply_pbc_edge_south_east(rhoB)
  
  call apply_pbc_edge_south_east(u)
  call apply_pbc_edge_south_east(v)
  call apply_pbc_edge_south_east(w)
  
  !apply pbc at south front 10
  call apply_pbc_edge_south_front(rhoR)
  call apply_pbc_edge_south_front(rhoB)
  
  call apply_pbc_edge_south_front(u)
  call apply_pbc_edge_south_front(v)
  call apply_pbc_edge_south_front(w)
  
  !apply pbc at south rear 11
  call apply_pbc_edge_south_rear(rhoR)
  call apply_pbc_edge_south_rear(rhoB)
  
  call apply_pbc_edge_south_rear(u)
  call apply_pbc_edge_south_rear(v)
  call apply_pbc_edge_south_rear(w)
  
  
  !apply pbc at south west 12
  call apply_pbc_edge_south_west(rhoR)
  call apply_pbc_edge_south_west(rhoB)
  
  call apply_pbc_edge_south_west(u)
  call apply_pbc_edge_south_west(v)
  call apply_pbc_edge_south_west(w)
  
  !corner 8
  
  !apply pbc at north east front 1
  call apply_pbc_corner_north_east_front(rhoR)
  call apply_pbc_corner_north_east_front(rhoB)
  
  call apply_pbc_corner_north_east_front(u)
  call apply_pbc_corner_north_east_front(v)
  call apply_pbc_corner_north_east_front(w)
  
  !apply pbc at north east rear 2
  call apply_pbc_corner_north_east_rear(rhoR)
  call apply_pbc_corner_north_east_rear(rhoB)
  
  call apply_pbc_corner_north_east_rear(u)
  call apply_pbc_corner_north_east_rear(v)
  call apply_pbc_corner_north_east_rear(w)
  
  !apply pbc at north west rear 3
  call apply_pbc_corner_north_west_rear(rhoR)
  call apply_pbc_corner_north_west_rear(rhoB)
  
  call apply_pbc_corner_north_west_rear(u)
  call apply_pbc_corner_north_west_rear(v)
  call apply_pbc_corner_north_west_rear(w)
  
  !apply pbc at north west front 4
  call apply_pbc_corner_north_west_front(rhoR)
  call apply_pbc_corner_north_west_front(rhoB)
  
  call apply_pbc_corner_north_west_front(u)
  call apply_pbc_corner_north_west_front(v)
  call apply_pbc_corner_north_west_front(w)
  
  !apply pbc at south east front 5
  call apply_pbc_corner_south_east_front(rhoR)
  call apply_pbc_corner_south_east_front(rhoB)
  
  call apply_pbc_corner_south_east_front(u)
  call apply_pbc_corner_south_east_front(v)
  call apply_pbc_corner_south_east_front(w)
  
  !apply pbc at south west front 6
  call apply_pbc_corner_south_west_front(rhoR)
  call apply_pbc_corner_south_west_front(rhoB)
  
  call apply_pbc_corner_south_west_front(u)
  call apply_pbc_corner_south_west_front(v)
  call apply_pbc_corner_south_west_front(w)
  
  !apply pbc at south west rear 7
  call apply_pbc_corner_south_west_rear(rhoR)
  call apply_pbc_corner_south_west_rear(rhoB)
  
  call apply_pbc_corner_south_west_rear(u)
  call apply_pbc_corner_south_west_rear(v)
  call apply_pbc_corner_south_west_rear(w)
  
  !apply pbc at south east rear 8
  call apply_pbc_corner_south_east_rear(rhoR)
  call apply_pbc_corner_south_east_rear(rhoB)
  
  call apply_pbc_corner_south_east_rear(u)
  call apply_pbc_corner_south_east_rear(v)
  call apply_pbc_corner_south_east_rear(w)
  
  return
  
 end subroutine apply_pbc_hvars
 
 subroutine driver_bc_pops
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the boundary conditions
!     to fluid populations
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  select case(ibctype)
  case(7)
    call apply_pbc_pops
  case default
    call apply_pbc_pops
  end select
  
  return
  
 end subroutine driver_bc_pops
 
 subroutine apply_pbc_pops 
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the full periodic boundary 
!     conditions to fluid populations
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  !sides 6
  
  !apply pbc at north 1
  !red fluid
  call apply_pbc_north(f00R)
  call apply_pbc_north(f01R)
  call apply_pbc_north(f02R)
  call apply_pbc_north(f03R)
  call apply_pbc_north(f04R)
  call apply_pbc_north(f05R)
  call apply_pbc_north(f06R)
  call apply_pbc_north(f07R)
  call apply_pbc_north(f08R)
  call apply_pbc_north(f09R)
  call apply_pbc_north(f10R)
  call apply_pbc_north(f11R)
  call apply_pbc_north(f12R)
  call apply_pbc_north(f13R)
  call apply_pbc_north(f14R)
  call apply_pbc_north(f15R)
  call apply_pbc_north(f16R)
  call apply_pbc_north(f17R)
  call apply_pbc_north(f18R)
  !blue fluid
  call apply_pbc_north(f00B)
  call apply_pbc_north(f01B)
  call apply_pbc_north(f02B)
  call apply_pbc_north(f03B)
  call apply_pbc_north(f04B)
  call apply_pbc_north(f05B)
  call apply_pbc_north(f06B)
  call apply_pbc_north(f07B)
  call apply_pbc_north(f08B)
  call apply_pbc_north(f09B)
  call apply_pbc_north(f10B)
  call apply_pbc_north(f11B)
  call apply_pbc_north(f12B)
  call apply_pbc_north(f13B)
  call apply_pbc_north(f14B)
  call apply_pbc_north(f15B)
  call apply_pbc_north(f16B)
  call apply_pbc_north(f17B)
  call apply_pbc_north(f18B)
  
  !apply pbc at south 2
  !red fluid
  call apply_pbc_south(f00R)
  call apply_pbc_south(f01R)
  call apply_pbc_south(f02R)
  call apply_pbc_south(f03R)
  call apply_pbc_south(f04R)
  call apply_pbc_south(f05R)
  call apply_pbc_south(f06R)
  call apply_pbc_south(f07R)
  call apply_pbc_south(f08R)
  call apply_pbc_south(f09R)
  call apply_pbc_south(f10R)
  call apply_pbc_south(f11R)
  call apply_pbc_south(f12R)
  call apply_pbc_south(f13R)
  call apply_pbc_south(f14R)
  call apply_pbc_south(f15R)
  call apply_pbc_south(f16R)
  call apply_pbc_south(f17R)
  call apply_pbc_south(f18R)
  !blue fluid
  call apply_pbc_south(f00B)
  call apply_pbc_south(f01B)
  call apply_pbc_south(f02B)
  call apply_pbc_south(f03B)
  call apply_pbc_south(f04B)
  call apply_pbc_south(f05B)
  call apply_pbc_south(f06B)
  call apply_pbc_south(f07B)
  call apply_pbc_south(f08B)
  call apply_pbc_south(f09B)
  call apply_pbc_south(f10B)
  call apply_pbc_south(f11B)
  call apply_pbc_south(f12B)
  call apply_pbc_south(f13B)
  call apply_pbc_south(f14B)
  call apply_pbc_south(f15B)
  call apply_pbc_south(f16B)
  call apply_pbc_south(f17B)
  call apply_pbc_south(f18B)
  
  !apply pbc at east 3
  !red fluid
  call apply_pbc_east(f00R)
  call apply_pbc_east(f01R)
  call apply_pbc_east(f02R)
  call apply_pbc_east(f03R)
  call apply_pbc_east(f04R)
  call apply_pbc_east(f05R)
  call apply_pbc_east(f06R)
  call apply_pbc_east(f07R)
  call apply_pbc_east(f08R)
  call apply_pbc_east(f09R)
  call apply_pbc_east(f10R)
  call apply_pbc_east(f11R)
  call apply_pbc_east(f12R)
  call apply_pbc_east(f13R)
  call apply_pbc_east(f14R)
  call apply_pbc_east(f15R)
  call apply_pbc_east(f16R)
  call apply_pbc_east(f17R)
  call apply_pbc_east(f18R)
  !blue fluid
  call apply_pbc_east(f00B)
  call apply_pbc_east(f01B)
  call apply_pbc_east(f02B)
  call apply_pbc_east(f03B)
  call apply_pbc_east(f04B)
  call apply_pbc_east(f05B)
  call apply_pbc_east(f06B)
  call apply_pbc_east(f07B)
  call apply_pbc_east(f08B)
  call apply_pbc_east(f09B)
  call apply_pbc_east(f10B)
  call apply_pbc_east(f11B)
  call apply_pbc_east(f12B)
  call apply_pbc_east(f13B)
  call apply_pbc_east(f14B)
  call apply_pbc_east(f15B)
  call apply_pbc_east(f16B)
  call apply_pbc_east(f17B)
  call apply_pbc_east(f18B)
  
  !apply pbc at west 4
  !red fluid
  call apply_pbc_west(f00R)
  call apply_pbc_west(f01R)
  call apply_pbc_west(f02R)
  call apply_pbc_west(f03R)
  call apply_pbc_west(f04R)
  call apply_pbc_west(f05R)
  call apply_pbc_west(f06R)
  call apply_pbc_west(f07R)
  call apply_pbc_west(f08R)
  call apply_pbc_west(f09R)
  call apply_pbc_west(f10R)
  call apply_pbc_west(f11R)
  call apply_pbc_west(f12R)
  call apply_pbc_west(f13R)
  call apply_pbc_west(f14R)
  call apply_pbc_west(f15R)
  call apply_pbc_west(f16R)
  call apply_pbc_west(f17R)
  call apply_pbc_west(f18R)
  !blue fluid
  call apply_pbc_west(f00B)
  call apply_pbc_west(f01B)
  call apply_pbc_west(f02B)
  call apply_pbc_west(f03B)
  call apply_pbc_west(f04B)
  call apply_pbc_west(f05B)
  call apply_pbc_west(f06B)
  call apply_pbc_west(f07B)
  call apply_pbc_west(f08B)
  call apply_pbc_west(f09B)
  call apply_pbc_west(f10B)
  call apply_pbc_west(f11B)
  call apply_pbc_west(f12B)
  call apply_pbc_west(f13B)
  call apply_pbc_west(f14B)
  call apply_pbc_west(f15B)
  call apply_pbc_west(f16B)
  call apply_pbc_west(f17B)
  call apply_pbc_west(f18B)
  
  !apply pbc at front 5
  !red fluid
  call apply_pbc_front(f00R)
  call apply_pbc_front(f01R)
  call apply_pbc_front(f02R)
  call apply_pbc_front(f03R)
  call apply_pbc_front(f04R)
  call apply_pbc_front(f05R)
  call apply_pbc_front(f06R)
  call apply_pbc_front(f07R)
  call apply_pbc_front(f08R)
  call apply_pbc_front(f09R)
  call apply_pbc_front(f10R)
  call apply_pbc_front(f11R)
  call apply_pbc_front(f12R)
  call apply_pbc_front(f13R)
  call apply_pbc_front(f14R)
  call apply_pbc_front(f15R)
  call apply_pbc_front(f16R)
  call apply_pbc_front(f17R)
  call apply_pbc_front(f18R)
  !blue fluid
  call apply_pbc_front(f00B)
  call apply_pbc_front(f01B)
  call apply_pbc_front(f02B)
  call apply_pbc_front(f03B)
  call apply_pbc_front(f04B)
  call apply_pbc_front(f05B)
  call apply_pbc_front(f06B)
  call apply_pbc_front(f07B)
  call apply_pbc_front(f08B)
  call apply_pbc_front(f09B)
  call apply_pbc_front(f10B)
  call apply_pbc_front(f11B)
  call apply_pbc_front(f12B)
  call apply_pbc_front(f13B)
  call apply_pbc_front(f14B)
  call apply_pbc_front(f15B)
  call apply_pbc_front(f16B)
  call apply_pbc_front(f17B)
  call apply_pbc_front(f18B)
  
  !apply pbc at rear 6
  !red fluid
  call apply_pbc_rear(f00R)
  call apply_pbc_rear(f01R)
  call apply_pbc_rear(f02R)
  call apply_pbc_rear(f03R)
  call apply_pbc_rear(f04R)
  call apply_pbc_rear(f05R)
  call apply_pbc_rear(f06R)
  call apply_pbc_rear(f07R)
  call apply_pbc_rear(f08R)
  call apply_pbc_rear(f09R)
  call apply_pbc_rear(f10R)
  call apply_pbc_rear(f11R)
  call apply_pbc_rear(f12R)
  call apply_pbc_rear(f13R)
  call apply_pbc_rear(f14R)
  call apply_pbc_rear(f15R)
  call apply_pbc_rear(f16R)
  call apply_pbc_rear(f17R)
  call apply_pbc_rear(f18R)
  !blue fluid
  call apply_pbc_rear(f00B)
  call apply_pbc_rear(f01B)
  call apply_pbc_rear(f02B)
  call apply_pbc_rear(f03B)
  call apply_pbc_rear(f04B)
  call apply_pbc_rear(f05B)
  call apply_pbc_rear(f06B)
  call apply_pbc_rear(f07B)
  call apply_pbc_rear(f08B)
  call apply_pbc_rear(f09B)
  call apply_pbc_rear(f10B)
  call apply_pbc_rear(f11B)
  call apply_pbc_rear(f12B)
  call apply_pbc_rear(f13B)
  call apply_pbc_rear(f14B)
  call apply_pbc_rear(f15B)
  call apply_pbc_rear(f16B)
  call apply_pbc_rear(f17B)
  call apply_pbc_rear(f18B)
  
  !edges 12
  
  !apply pbc at front east 1
  !red fluid
  call apply_pbc_edge_front_east(f00R)
  call apply_pbc_edge_front_east(f01R)
  call apply_pbc_edge_front_east(f02R)
  call apply_pbc_edge_front_east(f03R)
  call apply_pbc_edge_front_east(f04R)
  call apply_pbc_edge_front_east(f05R)
  call apply_pbc_edge_front_east(f06R)
  call apply_pbc_edge_front_east(f07R)
  call apply_pbc_edge_front_east(f08R)
  call apply_pbc_edge_front_east(f09R)
  call apply_pbc_edge_front_east(f10R)
  call apply_pbc_edge_front_east(f11R)
  call apply_pbc_edge_front_east(f12R)
  call apply_pbc_edge_front_east(f13R)
  call apply_pbc_edge_front_east(f14R)
  call apply_pbc_edge_front_east(f15R)
  call apply_pbc_edge_front_east(f16R)
  call apply_pbc_edge_front_east(f17R)
  call apply_pbc_edge_front_east(f18R)
  !blue fluid
  call apply_pbc_edge_front_east(f00B)
  call apply_pbc_edge_front_east(f01B)
  call apply_pbc_edge_front_east(f02B)
  call apply_pbc_edge_front_east(f03B)
  call apply_pbc_edge_front_east(f04B)
  call apply_pbc_edge_front_east(f05B)
  call apply_pbc_edge_front_east(f06B)
  call apply_pbc_edge_front_east(f07B)
  call apply_pbc_edge_front_east(f08B)
  call apply_pbc_edge_front_east(f09B)
  call apply_pbc_edge_front_east(f10B)
  call apply_pbc_edge_front_east(f11B)
  call apply_pbc_edge_front_east(f12B)
  call apply_pbc_edge_front_east(f13B)
  call apply_pbc_edge_front_east(f14B)
  call apply_pbc_edge_front_east(f15B)
  call apply_pbc_edge_front_east(f16B)
  call apply_pbc_edge_front_east(f17B)
  call apply_pbc_edge_front_east(f18B)
  
  !apply pbc at front west 2
  !red fluid
  call apply_pbc_edge_front_west(f00R)
  call apply_pbc_edge_front_west(f01R)
  call apply_pbc_edge_front_west(f02R)
  call apply_pbc_edge_front_west(f03R)
  call apply_pbc_edge_front_west(f04R)
  call apply_pbc_edge_front_west(f05R)
  call apply_pbc_edge_front_west(f06R)
  call apply_pbc_edge_front_west(f07R)
  call apply_pbc_edge_front_west(f08R)
  call apply_pbc_edge_front_west(f09R)
  call apply_pbc_edge_front_west(f10R)
  call apply_pbc_edge_front_west(f11R)
  call apply_pbc_edge_front_west(f12R)
  call apply_pbc_edge_front_west(f13R)
  call apply_pbc_edge_front_west(f14R)
  call apply_pbc_edge_front_west(f15R)
  call apply_pbc_edge_front_west(f16R)
  call apply_pbc_edge_front_west(f17R)
  call apply_pbc_edge_front_west(f18R)
  !blue fluid
  call apply_pbc_edge_front_west(f00B)
  call apply_pbc_edge_front_west(f01B)
  call apply_pbc_edge_front_west(f02B)
  call apply_pbc_edge_front_west(f03B)
  call apply_pbc_edge_front_west(f04B)
  call apply_pbc_edge_front_west(f05B)
  call apply_pbc_edge_front_west(f06B)
  call apply_pbc_edge_front_west(f07B)
  call apply_pbc_edge_front_west(f08B)
  call apply_pbc_edge_front_west(f09B)
  call apply_pbc_edge_front_west(f10B)
  call apply_pbc_edge_front_west(f11B)
  call apply_pbc_edge_front_west(f12B)
  call apply_pbc_edge_front_west(f13B)
  call apply_pbc_edge_front_west(f14B)
  call apply_pbc_edge_front_west(f15B)
  call apply_pbc_edge_front_west(f16B)
  call apply_pbc_edge_front_west(f17B)
  call apply_pbc_edge_front_west(f18B)
  
  !apply pbc at north east 3
  !red fluid
  call apply_pbc_edge_north_east(f00R)
  call apply_pbc_edge_north_east(f01R)
  call apply_pbc_edge_north_east(f02R)
  call apply_pbc_edge_north_east(f03R)
  call apply_pbc_edge_north_east(f04R)
  call apply_pbc_edge_north_east(f05R)
  call apply_pbc_edge_north_east(f06R)
  call apply_pbc_edge_north_east(f07R)
  call apply_pbc_edge_north_east(f08R)
  call apply_pbc_edge_north_east(f09R)
  call apply_pbc_edge_north_east(f10R)
  call apply_pbc_edge_north_east(f11R)
  call apply_pbc_edge_north_east(f12R)
  call apply_pbc_edge_north_east(f13R)
  call apply_pbc_edge_north_east(f14R)
  call apply_pbc_edge_north_east(f15R)
  call apply_pbc_edge_north_east(f16R)
  call apply_pbc_edge_north_east(f17R)
  call apply_pbc_edge_north_east(f18R)
  !blue fluid
  call apply_pbc_edge_north_east(f00B)
  call apply_pbc_edge_north_east(f01B)
  call apply_pbc_edge_north_east(f02B)
  call apply_pbc_edge_north_east(f03B)
  call apply_pbc_edge_north_east(f04B)
  call apply_pbc_edge_north_east(f05B)
  call apply_pbc_edge_north_east(f06B)
  call apply_pbc_edge_north_east(f07B)
  call apply_pbc_edge_north_east(f08B)
  call apply_pbc_edge_north_east(f09B)
  call apply_pbc_edge_north_east(f10B)
  call apply_pbc_edge_north_east(f11B)
  call apply_pbc_edge_north_east(f12B)
  call apply_pbc_edge_north_east(f13B)
  call apply_pbc_edge_north_east(f14B)
  call apply_pbc_edge_north_east(f15B)
  call apply_pbc_edge_north_east(f16B)
  call apply_pbc_edge_north_east(f17B)
  call apply_pbc_edge_north_east(f18B)
  
  !apply pbc at north front 4
  call apply_pbc_edge_north_front(f00R)
  call apply_pbc_edge_north_front(f01R)
  call apply_pbc_edge_north_front(f02R)
  call apply_pbc_edge_north_front(f03R)
  call apply_pbc_edge_north_front(f04R)
  call apply_pbc_edge_north_front(f05R)
  call apply_pbc_edge_north_front(f06R)
  call apply_pbc_edge_north_front(f07R)
  call apply_pbc_edge_north_front(f08R)
  call apply_pbc_edge_north_front(f09R)
  call apply_pbc_edge_north_front(f10R)
  call apply_pbc_edge_north_front(f11R)
  call apply_pbc_edge_north_front(f12R)
  call apply_pbc_edge_north_front(f13R)
  call apply_pbc_edge_north_front(f14R)
  call apply_pbc_edge_north_front(f15R)
  call apply_pbc_edge_north_front(f16R)
  call apply_pbc_edge_north_front(f17R)
  call apply_pbc_edge_north_front(f18R)
  !blue fluid
  call apply_pbc_edge_north_front(f00B)
  call apply_pbc_edge_north_front(f01B)
  call apply_pbc_edge_north_front(f02B)
  call apply_pbc_edge_north_front(f03B)
  call apply_pbc_edge_north_front(f04B)
  call apply_pbc_edge_north_front(f05B)
  call apply_pbc_edge_north_front(f06B)
  call apply_pbc_edge_north_front(f07B)
  call apply_pbc_edge_north_front(f08B)
  call apply_pbc_edge_north_front(f09B)
  call apply_pbc_edge_north_front(f10B)
  call apply_pbc_edge_north_front(f11B)
  call apply_pbc_edge_north_front(f12B)
  call apply_pbc_edge_north_front(f13B)
  call apply_pbc_edge_north_front(f14B)
  call apply_pbc_edge_north_front(f15B)
  call apply_pbc_edge_north_front(f16B)
  call apply_pbc_edge_north_front(f17B)
  call apply_pbc_edge_north_front(f18B)
  
  !apply pbc at north rear 5
  !red fluid
  call apply_pbc_edge_north_rear(f00R)
  call apply_pbc_edge_north_rear(f01R)
  call apply_pbc_edge_north_rear(f02R)
  call apply_pbc_edge_north_rear(f03R)
  call apply_pbc_edge_north_rear(f04R)
  call apply_pbc_edge_north_rear(f05R)
  call apply_pbc_edge_north_rear(f06R)
  call apply_pbc_edge_north_rear(f07R)
  call apply_pbc_edge_north_rear(f08R)
  call apply_pbc_edge_north_rear(f09R)
  call apply_pbc_edge_north_rear(f10R)
  call apply_pbc_edge_north_rear(f11R)
  call apply_pbc_edge_north_rear(f12R)
  call apply_pbc_edge_north_rear(f13R)
  call apply_pbc_edge_north_rear(f14R)
  call apply_pbc_edge_north_rear(f15R)
  call apply_pbc_edge_north_rear(f16R)
  call apply_pbc_edge_north_rear(f17R)
  call apply_pbc_edge_north_rear(f18R)
  !blue fluid
  call apply_pbc_edge_north_rear(f00B)
  call apply_pbc_edge_north_rear(f01B)
  call apply_pbc_edge_north_rear(f02B)
  call apply_pbc_edge_north_rear(f03B)
  call apply_pbc_edge_north_rear(f04B)
  call apply_pbc_edge_north_rear(f05B)
  call apply_pbc_edge_north_rear(f06B)
  call apply_pbc_edge_north_rear(f07B)
  call apply_pbc_edge_north_rear(f08B)
  call apply_pbc_edge_north_rear(f09B)
  call apply_pbc_edge_north_rear(f10B)
  call apply_pbc_edge_north_rear(f11B)
  call apply_pbc_edge_north_rear(f12B)
  call apply_pbc_edge_north_rear(f13B)
  call apply_pbc_edge_north_rear(f14B)
  call apply_pbc_edge_north_rear(f15B)
  call apply_pbc_edge_north_rear(f16B)
  call apply_pbc_edge_north_rear(f17B)
  call apply_pbc_edge_north_rear(f18B)
  
  !apply pbc at north west 6
  !red fluid
  call apply_pbc_edge_north_west(f00R)
  call apply_pbc_edge_north_west(f01R)
  call apply_pbc_edge_north_west(f02R)
  call apply_pbc_edge_north_west(f03R)
  call apply_pbc_edge_north_west(f04R)
  call apply_pbc_edge_north_west(f05R)
  call apply_pbc_edge_north_west(f06R)
  call apply_pbc_edge_north_west(f07R)
  call apply_pbc_edge_north_west(f08R)
  call apply_pbc_edge_north_west(f09R)
  call apply_pbc_edge_north_west(f10R)
  call apply_pbc_edge_north_west(f11R)
  call apply_pbc_edge_north_west(f12R)
  call apply_pbc_edge_north_west(f13R)
  call apply_pbc_edge_north_west(f14R)
  call apply_pbc_edge_north_west(f15R)
  call apply_pbc_edge_north_west(f16R)
  call apply_pbc_edge_north_west(f17R)
  call apply_pbc_edge_north_west(f18R)
  !blue fluid
  call apply_pbc_edge_north_west(f00B)
  call apply_pbc_edge_north_west(f01B)
  call apply_pbc_edge_north_west(f02B)
  call apply_pbc_edge_north_west(f03B)
  call apply_pbc_edge_north_west(f04B)
  call apply_pbc_edge_north_west(f05B)
  call apply_pbc_edge_north_west(f06B)
  call apply_pbc_edge_north_west(f07B)
  call apply_pbc_edge_north_west(f08B)
  call apply_pbc_edge_north_west(f09B)
  call apply_pbc_edge_north_west(f10B)
  call apply_pbc_edge_north_west(f11B)
  call apply_pbc_edge_north_west(f12B)
  call apply_pbc_edge_north_west(f13B)
  call apply_pbc_edge_north_west(f14B)
  call apply_pbc_edge_north_west(f15B)
  call apply_pbc_edge_north_west(f16B)
  call apply_pbc_edge_north_west(f17B)
  call apply_pbc_edge_north_west(f18B)
  
  !apply pbc at rear east 7
  !red fluid
  call apply_pbc_edge_rear_east(f00R)
  call apply_pbc_edge_rear_east(f01R)
  call apply_pbc_edge_rear_east(f02R)
  call apply_pbc_edge_rear_east(f03R)
  call apply_pbc_edge_rear_east(f04R)
  call apply_pbc_edge_rear_east(f05R)
  call apply_pbc_edge_rear_east(f06R)
  call apply_pbc_edge_rear_east(f07R)
  call apply_pbc_edge_rear_east(f08R)
  call apply_pbc_edge_rear_east(f09R)
  call apply_pbc_edge_rear_east(f10R)
  call apply_pbc_edge_rear_east(f11R)
  call apply_pbc_edge_rear_east(f12R)
  call apply_pbc_edge_rear_east(f13R)
  call apply_pbc_edge_rear_east(f14R)
  call apply_pbc_edge_rear_east(f15R)
  call apply_pbc_edge_rear_east(f16R)
  call apply_pbc_edge_rear_east(f17R)
  call apply_pbc_edge_rear_east(f18R)
  !blue fluid
  call apply_pbc_edge_rear_east(f00B)
  call apply_pbc_edge_rear_east(f01B)
  call apply_pbc_edge_rear_east(f02B)
  call apply_pbc_edge_rear_east(f03B)
  call apply_pbc_edge_rear_east(f04B)
  call apply_pbc_edge_rear_east(f05B)
  call apply_pbc_edge_rear_east(f06B)
  call apply_pbc_edge_rear_east(f07B)
  call apply_pbc_edge_rear_east(f08B)
  call apply_pbc_edge_rear_east(f09B)
  call apply_pbc_edge_rear_east(f10B)
  call apply_pbc_edge_rear_east(f11B)
  call apply_pbc_edge_rear_east(f12B)
  call apply_pbc_edge_rear_east(f13B)
  call apply_pbc_edge_rear_east(f14B)
  call apply_pbc_edge_rear_east(f15B)
  call apply_pbc_edge_rear_east(f16B)
  call apply_pbc_edge_rear_east(f17B)
  call apply_pbc_edge_rear_east(f18B)
  
  !apply pbc at rear west 8
  !red fluid
  call apply_pbc_edge_rear_west(f00R)
  call apply_pbc_edge_rear_west(f01R)
  call apply_pbc_edge_rear_west(f02R)
  call apply_pbc_edge_rear_west(f03R)
  call apply_pbc_edge_rear_west(f04R)
  call apply_pbc_edge_rear_west(f05R)
  call apply_pbc_edge_rear_west(f06R)
  call apply_pbc_edge_rear_west(f07R)
  call apply_pbc_edge_rear_west(f08R)
  call apply_pbc_edge_rear_west(f09R)
  call apply_pbc_edge_rear_west(f10R)
  call apply_pbc_edge_rear_west(f11R)
  call apply_pbc_edge_rear_west(f12R)
  call apply_pbc_edge_rear_west(f13R)
  call apply_pbc_edge_rear_west(f14R)
  call apply_pbc_edge_rear_west(f15R)
  call apply_pbc_edge_rear_west(f16R)
  call apply_pbc_edge_rear_west(f17R)
  call apply_pbc_edge_rear_west(f18R)
  !blue fluid
  call apply_pbc_edge_rear_west(f00B)
  call apply_pbc_edge_rear_west(f01B)
  call apply_pbc_edge_rear_west(f02B)
  call apply_pbc_edge_rear_west(f03B)
  call apply_pbc_edge_rear_west(f04B)
  call apply_pbc_edge_rear_west(f05B)
  call apply_pbc_edge_rear_west(f06B)
  call apply_pbc_edge_rear_west(f07B)
  call apply_pbc_edge_rear_west(f08B)
  call apply_pbc_edge_rear_west(f09B)
  call apply_pbc_edge_rear_west(f10B)
  call apply_pbc_edge_rear_west(f11B)
  call apply_pbc_edge_rear_west(f12B)
  call apply_pbc_edge_rear_west(f13B)
  call apply_pbc_edge_rear_west(f14B)
  call apply_pbc_edge_rear_west(f15B)
  call apply_pbc_edge_rear_west(f16B)
  call apply_pbc_edge_rear_west(f17B)
  call apply_pbc_edge_rear_west(f18B)
  
  !apply pbc at south east 9
  !red fluid
  call apply_pbc_edge_south_east(f00R)
  call apply_pbc_edge_south_east(f01R)
  call apply_pbc_edge_south_east(f02R)
  call apply_pbc_edge_south_east(f03R)
  call apply_pbc_edge_south_east(f04R)
  call apply_pbc_edge_south_east(f05R)
  call apply_pbc_edge_south_east(f06R)
  call apply_pbc_edge_south_east(f07R)
  call apply_pbc_edge_south_east(f08R)
  call apply_pbc_edge_south_east(f09R)
  call apply_pbc_edge_south_east(f10R)
  call apply_pbc_edge_south_east(f11R)
  call apply_pbc_edge_south_east(f12R)
  call apply_pbc_edge_south_east(f13R)
  call apply_pbc_edge_south_east(f14R)
  call apply_pbc_edge_south_east(f15R)
  call apply_pbc_edge_south_east(f16R)
  call apply_pbc_edge_south_east(f17R)
  call apply_pbc_edge_south_east(f18R)
  !blue fluid
  call apply_pbc_edge_south_east(f00B)
  call apply_pbc_edge_south_east(f01B)
  call apply_pbc_edge_south_east(f02B)
  call apply_pbc_edge_south_east(f03B)
  call apply_pbc_edge_south_east(f04B)
  call apply_pbc_edge_south_east(f05B)
  call apply_pbc_edge_south_east(f06B)
  call apply_pbc_edge_south_east(f07B)
  call apply_pbc_edge_south_east(f08B)
  call apply_pbc_edge_south_east(f09B)
  call apply_pbc_edge_south_east(f10B)
  call apply_pbc_edge_south_east(f11B)
  call apply_pbc_edge_south_east(f12B)
  call apply_pbc_edge_south_east(f13B)
  call apply_pbc_edge_south_east(f14B)
  call apply_pbc_edge_south_east(f15B)
  call apply_pbc_edge_south_east(f16B)
  call apply_pbc_edge_south_east(f17B)
  call apply_pbc_edge_south_east(f18B)
  
  !apply pbc at south front 10
  !red fluid
  call apply_pbc_edge_south_front(f00R)
  call apply_pbc_edge_south_front(f01R)
  call apply_pbc_edge_south_front(f02R)
  call apply_pbc_edge_south_front(f03R)
  call apply_pbc_edge_south_front(f04R)
  call apply_pbc_edge_south_front(f05R)
  call apply_pbc_edge_south_front(f06R)
  call apply_pbc_edge_south_front(f07R)
  call apply_pbc_edge_south_front(f08R)
  call apply_pbc_edge_south_front(f09R)
  call apply_pbc_edge_south_front(f10R)
  call apply_pbc_edge_south_front(f11R)
  call apply_pbc_edge_south_front(f12R)
  call apply_pbc_edge_south_front(f13R)
  call apply_pbc_edge_south_front(f14R)
  call apply_pbc_edge_south_front(f15R)
  call apply_pbc_edge_south_front(f16R)
  call apply_pbc_edge_south_front(f17R)
  call apply_pbc_edge_south_front(f18R)
  !blue fluid
  call apply_pbc_edge_south_front(f00B)
  call apply_pbc_edge_south_front(f01B)
  call apply_pbc_edge_south_front(f02B)
  call apply_pbc_edge_south_front(f03B)
  call apply_pbc_edge_south_front(f04B)
  call apply_pbc_edge_south_front(f05B)
  call apply_pbc_edge_south_front(f06B)
  call apply_pbc_edge_south_front(f07B)
  call apply_pbc_edge_south_front(f08B)
  call apply_pbc_edge_south_front(f09B)
  call apply_pbc_edge_south_front(f10B)
  call apply_pbc_edge_south_front(f11B)
  call apply_pbc_edge_south_front(f12B)
  call apply_pbc_edge_south_front(f13B)
  call apply_pbc_edge_south_front(f14B)
  call apply_pbc_edge_south_front(f15B)
  call apply_pbc_edge_south_front(f16B)
  call apply_pbc_edge_south_front(f17B)
  call apply_pbc_edge_south_front(f18B)
  
  !apply pbc at south rear 11
  !red fluid
  call apply_pbc_edge_south_rear(f00R)
  call apply_pbc_edge_south_rear(f01R)
  call apply_pbc_edge_south_rear(f02R)
  call apply_pbc_edge_south_rear(f03R)
  call apply_pbc_edge_south_rear(f04R)
  call apply_pbc_edge_south_rear(f05R)
  call apply_pbc_edge_south_rear(f06R)
  call apply_pbc_edge_south_rear(f07R)
  call apply_pbc_edge_south_rear(f08R)
  call apply_pbc_edge_south_rear(f09R)
  call apply_pbc_edge_south_rear(f10R)
  call apply_pbc_edge_south_rear(f11R)
  call apply_pbc_edge_south_rear(f12R)
  call apply_pbc_edge_south_rear(f13R)
  call apply_pbc_edge_south_rear(f14R)
  call apply_pbc_edge_south_rear(f15R)
  call apply_pbc_edge_south_rear(f16R)
  call apply_pbc_edge_south_rear(f17R)
  call apply_pbc_edge_south_rear(f18R)
  !blue fluid
  call apply_pbc_edge_south_rear(f00B)
  call apply_pbc_edge_south_rear(f01B)
  call apply_pbc_edge_south_rear(f02B)
  call apply_pbc_edge_south_rear(f03B)
  call apply_pbc_edge_south_rear(f04B)
  call apply_pbc_edge_south_rear(f05B)
  call apply_pbc_edge_south_rear(f06B)
  call apply_pbc_edge_south_rear(f07B)
  call apply_pbc_edge_south_rear(f08B)
  call apply_pbc_edge_south_rear(f09B)
  call apply_pbc_edge_south_rear(f10B)
  call apply_pbc_edge_south_rear(f11B)
  call apply_pbc_edge_south_rear(f12B)
  call apply_pbc_edge_south_rear(f13B)
  call apply_pbc_edge_south_rear(f14B)
  call apply_pbc_edge_south_rear(f15B)
  call apply_pbc_edge_south_rear(f16B)
  call apply_pbc_edge_south_rear(f17B)
  call apply_pbc_edge_south_rear(f18B)
  
  !apply pbc at south west 12
  !red fluid
  call apply_pbc_edge_south_west(f00R)
  call apply_pbc_edge_south_west(f01R)
  call apply_pbc_edge_south_west(f02R)
  call apply_pbc_edge_south_west(f03R)
  call apply_pbc_edge_south_west(f04R)
  call apply_pbc_edge_south_west(f05R)
  call apply_pbc_edge_south_west(f06R)
  call apply_pbc_edge_south_west(f07R)
  call apply_pbc_edge_south_west(f08R)
  call apply_pbc_edge_south_west(f09R)
  call apply_pbc_edge_south_west(f10R)
  call apply_pbc_edge_south_west(f11R)
  call apply_pbc_edge_south_west(f12R)
  call apply_pbc_edge_south_west(f13R)
  call apply_pbc_edge_south_west(f14R)
  call apply_pbc_edge_south_west(f15R)
  call apply_pbc_edge_south_west(f16R)
  call apply_pbc_edge_south_west(f17R)
  call apply_pbc_edge_south_west(f18R)
  !blue fluid
  call apply_pbc_edge_south_west(f00B)
  call apply_pbc_edge_south_west(f01B)
  call apply_pbc_edge_south_west(f02B)
  call apply_pbc_edge_south_west(f03B)
  call apply_pbc_edge_south_west(f04B)
  call apply_pbc_edge_south_west(f05B)
  call apply_pbc_edge_south_west(f06B)
  call apply_pbc_edge_south_west(f07B)
  call apply_pbc_edge_south_west(f08B)
  call apply_pbc_edge_south_west(f09B)
  call apply_pbc_edge_south_west(f10B)
  call apply_pbc_edge_south_west(f11B)
  call apply_pbc_edge_south_west(f12B)
  call apply_pbc_edge_south_west(f13B)
  call apply_pbc_edge_south_west(f14B)
  call apply_pbc_edge_south_west(f15B)
  call apply_pbc_edge_south_west(f16B)
  call apply_pbc_edge_south_west(f17B)
  call apply_pbc_edge_south_west(f18B)
  
  !corner 8
  
  !apply pbc at north east front 1
  !red fluid
  call apply_pbc_corner_north_east_front(f00R)
  call apply_pbc_corner_north_east_front(f01R)
  call apply_pbc_corner_north_east_front(f02R)
  call apply_pbc_corner_north_east_front(f03R)
  call apply_pbc_corner_north_east_front(f04R)
  call apply_pbc_corner_north_east_front(f05R)
  call apply_pbc_corner_north_east_front(f06R)
  call apply_pbc_corner_north_east_front(f07R)
  call apply_pbc_corner_north_east_front(f08R)
  call apply_pbc_corner_north_east_front(f09R)
  call apply_pbc_corner_north_east_front(f10R)
  call apply_pbc_corner_north_east_front(f11R)
  call apply_pbc_corner_north_east_front(f12R)
  call apply_pbc_corner_north_east_front(f13R)
  call apply_pbc_corner_north_east_front(f14R)
  call apply_pbc_corner_north_east_front(f15R)
  call apply_pbc_corner_north_east_front(f16R)
  call apply_pbc_corner_north_east_front(f17R)
  call apply_pbc_corner_north_east_front(f18R)
  !blue fluid
  call apply_pbc_corner_north_east_front(f00B)
  call apply_pbc_corner_north_east_front(f01B)
  call apply_pbc_corner_north_east_front(f02B)
  call apply_pbc_corner_north_east_front(f03B)
  call apply_pbc_corner_north_east_front(f04B)
  call apply_pbc_corner_north_east_front(f05B)
  call apply_pbc_corner_north_east_front(f06B)
  call apply_pbc_corner_north_east_front(f07B)
  call apply_pbc_corner_north_east_front(f08B)
  call apply_pbc_corner_north_east_front(f09B)
  call apply_pbc_corner_north_east_front(f10B)
  call apply_pbc_corner_north_east_front(f11B)
  call apply_pbc_corner_north_east_front(f12B)
  call apply_pbc_corner_north_east_front(f13B)
  call apply_pbc_corner_north_east_front(f14B)
  call apply_pbc_corner_north_east_front(f15B)
  call apply_pbc_corner_north_east_front(f16B)
  call apply_pbc_corner_north_east_front(f17B)
  call apply_pbc_corner_north_east_front(f18B)
  
  !apply pbc at north east rear 2
  !red fluid
  call apply_pbc_corner_north_east_rear(f00R)
  call apply_pbc_corner_north_east_rear(f01R)
  call apply_pbc_corner_north_east_rear(f02R)
  call apply_pbc_corner_north_east_rear(f03R)
  call apply_pbc_corner_north_east_rear(f04R)
  call apply_pbc_corner_north_east_rear(f05R)
  call apply_pbc_corner_north_east_rear(f06R)
  call apply_pbc_corner_north_east_rear(f07R)
  call apply_pbc_corner_north_east_rear(f08R)
  call apply_pbc_corner_north_east_rear(f09R)
  call apply_pbc_corner_north_east_rear(f10R)
  call apply_pbc_corner_north_east_rear(f11R)
  call apply_pbc_corner_north_east_rear(f12R)
  call apply_pbc_corner_north_east_rear(f13R)
  call apply_pbc_corner_north_east_rear(f14R)
  call apply_pbc_corner_north_east_rear(f15R)
  call apply_pbc_corner_north_east_rear(f16R)
  call apply_pbc_corner_north_east_rear(f17R)
  call apply_pbc_corner_north_east_rear(f18R)
  !blue fluid
  call apply_pbc_corner_north_east_rear(f00B)
  call apply_pbc_corner_north_east_rear(f01B)
  call apply_pbc_corner_north_east_rear(f02B)
  call apply_pbc_corner_north_east_rear(f03B)
  call apply_pbc_corner_north_east_rear(f04B)
  call apply_pbc_corner_north_east_rear(f05B)
  call apply_pbc_corner_north_east_rear(f06B)
  call apply_pbc_corner_north_east_rear(f07B)
  call apply_pbc_corner_north_east_rear(f08B)
  call apply_pbc_corner_north_east_rear(f09B)
  call apply_pbc_corner_north_east_rear(f10B)
  call apply_pbc_corner_north_east_rear(f11B)
  call apply_pbc_corner_north_east_rear(f12B)
  call apply_pbc_corner_north_east_rear(f13B)
  call apply_pbc_corner_north_east_rear(f14B)
  call apply_pbc_corner_north_east_rear(f15B)
  call apply_pbc_corner_north_east_rear(f16B)
  call apply_pbc_corner_north_east_rear(f17B)
  call apply_pbc_corner_north_east_rear(f18B)
  
  !apply pbc at north west rear 3
  !red fluid
  call apply_pbc_corner_north_west_rear(f00R)
  call apply_pbc_corner_north_west_rear(f01R)
  call apply_pbc_corner_north_west_rear(f02R)
  call apply_pbc_corner_north_west_rear(f03R)
  call apply_pbc_corner_north_west_rear(f04R)
  call apply_pbc_corner_north_west_rear(f05R)
  call apply_pbc_corner_north_west_rear(f06R)
  call apply_pbc_corner_north_west_rear(f07R)
  call apply_pbc_corner_north_west_rear(f08R)
  call apply_pbc_corner_north_west_rear(f09R)
  call apply_pbc_corner_north_west_rear(f10R)
  call apply_pbc_corner_north_west_rear(f11R)
  call apply_pbc_corner_north_west_rear(f12R)
  call apply_pbc_corner_north_west_rear(f13R)
  call apply_pbc_corner_north_west_rear(f14R)
  call apply_pbc_corner_north_west_rear(f15R)
  call apply_pbc_corner_north_west_rear(f16R)
  call apply_pbc_corner_north_west_rear(f17R)
  call apply_pbc_corner_north_west_rear(f18R)
  !blue fluid
  call apply_pbc_corner_north_west_rear(f00B)
  call apply_pbc_corner_north_west_rear(f01B)
  call apply_pbc_corner_north_west_rear(f02B)
  call apply_pbc_corner_north_west_rear(f03B)
  call apply_pbc_corner_north_west_rear(f04B)
  call apply_pbc_corner_north_west_rear(f05B)
  call apply_pbc_corner_north_west_rear(f06B)
  call apply_pbc_corner_north_west_rear(f07B)
  call apply_pbc_corner_north_west_rear(f08B)
  call apply_pbc_corner_north_west_rear(f09B)
  call apply_pbc_corner_north_west_rear(f10B)
  call apply_pbc_corner_north_west_rear(f11B)
  call apply_pbc_corner_north_west_rear(f12B)
  call apply_pbc_corner_north_west_rear(f13B)
  call apply_pbc_corner_north_west_rear(f14B)
  call apply_pbc_corner_north_west_rear(f15B)
  call apply_pbc_corner_north_west_rear(f16B)
  call apply_pbc_corner_north_west_rear(f17B)
  call apply_pbc_corner_north_west_rear(f18B)
  
  !apply pbc at north west front 4
  !red fluid
  call apply_pbc_corner_north_west_front(f00R)
  call apply_pbc_corner_north_west_front(f01R)
  call apply_pbc_corner_north_west_front(f02R)
  call apply_pbc_corner_north_west_front(f03R)
  call apply_pbc_corner_north_west_front(f04R)
  call apply_pbc_corner_north_west_front(f05R)
  call apply_pbc_corner_north_west_front(f06R)
  call apply_pbc_corner_north_west_front(f07R)
  call apply_pbc_corner_north_west_front(f08R)
  call apply_pbc_corner_north_west_front(f09R)
  call apply_pbc_corner_north_west_front(f10R)
  call apply_pbc_corner_north_west_front(f11R)
  call apply_pbc_corner_north_west_front(f12R)
  call apply_pbc_corner_north_west_front(f13R)
  call apply_pbc_corner_north_west_front(f14R)
  call apply_pbc_corner_north_west_front(f15R)
  call apply_pbc_corner_north_west_front(f16R)
  call apply_pbc_corner_north_west_front(f17R)
  call apply_pbc_corner_north_west_front(f18R)
  !blue fluid
  call apply_pbc_corner_north_west_front(f00B)
  call apply_pbc_corner_north_west_front(f01B)
  call apply_pbc_corner_north_west_front(f02B)
  call apply_pbc_corner_north_west_front(f03B)
  call apply_pbc_corner_north_west_front(f04B)
  call apply_pbc_corner_north_west_front(f05B)
  call apply_pbc_corner_north_west_front(f06B)
  call apply_pbc_corner_north_west_front(f07B)
  call apply_pbc_corner_north_west_front(f08B)
  call apply_pbc_corner_north_west_front(f09B)
  call apply_pbc_corner_north_west_front(f10B)
  call apply_pbc_corner_north_west_front(f11B)
  call apply_pbc_corner_north_west_front(f12B)
  call apply_pbc_corner_north_west_front(f13B)
  call apply_pbc_corner_north_west_front(f14B)
  call apply_pbc_corner_north_west_front(f15B)
  call apply_pbc_corner_north_west_front(f16B)
  call apply_pbc_corner_north_west_front(f17B)
  call apply_pbc_corner_north_west_front(f18B)
  
  !apply pbc at south east front 5
  !red fluid
  call apply_pbc_corner_south_east_front(f00R)
  call apply_pbc_corner_south_east_front(f01R)
  call apply_pbc_corner_south_east_front(f02R)
  call apply_pbc_corner_south_east_front(f03R)
  call apply_pbc_corner_south_east_front(f04R)
  call apply_pbc_corner_south_east_front(f05R)
  call apply_pbc_corner_south_east_front(f06R)
  call apply_pbc_corner_south_east_front(f07R)
  call apply_pbc_corner_south_east_front(f08R)
  call apply_pbc_corner_south_east_front(f09R)
  call apply_pbc_corner_south_east_front(f10R)
  call apply_pbc_corner_south_east_front(f11R)
  call apply_pbc_corner_south_east_front(f12R)
  call apply_pbc_corner_south_east_front(f13R)
  call apply_pbc_corner_south_east_front(f14R)
  call apply_pbc_corner_south_east_front(f15R)
  call apply_pbc_corner_south_east_front(f16R)
  call apply_pbc_corner_south_east_front(f17R)
  call apply_pbc_corner_south_east_front(f18R)
  !blue fluid
  call apply_pbc_corner_south_east_front(f00B)
  call apply_pbc_corner_south_east_front(f01B)
  call apply_pbc_corner_south_east_front(f02B)
  call apply_pbc_corner_south_east_front(f03B)
  call apply_pbc_corner_south_east_front(f04B)
  call apply_pbc_corner_south_east_front(f05B)
  call apply_pbc_corner_south_east_front(f06B)
  call apply_pbc_corner_south_east_front(f07B)
  call apply_pbc_corner_south_east_front(f08B)
  call apply_pbc_corner_south_east_front(f09B)
  call apply_pbc_corner_south_east_front(f10B)
  call apply_pbc_corner_south_east_front(f11B)
  call apply_pbc_corner_south_east_front(f12B)
  call apply_pbc_corner_south_east_front(f13B)
  call apply_pbc_corner_south_east_front(f14B)
  call apply_pbc_corner_south_east_front(f15B)
  call apply_pbc_corner_south_east_front(f16B)
  call apply_pbc_corner_south_east_front(f17B)
  call apply_pbc_corner_south_east_front(f18B)
  
  !apply pbc at south west front 6
  !red fluid
  call apply_pbc_corner_south_west_front(f00R)
  call apply_pbc_corner_south_west_front(f01R)
  call apply_pbc_corner_south_west_front(f02R)
  call apply_pbc_corner_south_west_front(f03R)
  call apply_pbc_corner_south_west_front(f04R)
  call apply_pbc_corner_south_west_front(f05R)
  call apply_pbc_corner_south_west_front(f06R)
  call apply_pbc_corner_south_west_front(f07R)
  call apply_pbc_corner_south_west_front(f08R)
  call apply_pbc_corner_south_west_front(f09R)
  call apply_pbc_corner_south_west_front(f10R)
  call apply_pbc_corner_south_west_front(f11R)
  call apply_pbc_corner_south_west_front(f12R)
  call apply_pbc_corner_south_west_front(f13R)
  call apply_pbc_corner_south_west_front(f14R)
  call apply_pbc_corner_south_west_front(f15R)
  call apply_pbc_corner_south_west_front(f16R)
  call apply_pbc_corner_south_west_front(f17R)
  call apply_pbc_corner_south_west_front(f18R)
  !blue fluid
  call apply_pbc_corner_south_west_front(f00B)
  call apply_pbc_corner_south_west_front(f01B)
  call apply_pbc_corner_south_west_front(f02B)
  call apply_pbc_corner_south_west_front(f03B)
  call apply_pbc_corner_south_west_front(f04B)
  call apply_pbc_corner_south_west_front(f05B)
  call apply_pbc_corner_south_west_front(f06B)
  call apply_pbc_corner_south_west_front(f07B)
  call apply_pbc_corner_south_west_front(f08B)
  call apply_pbc_corner_south_west_front(f09B)
  call apply_pbc_corner_south_west_front(f10B)
  call apply_pbc_corner_south_west_front(f11B)
  call apply_pbc_corner_south_west_front(f12B)
  call apply_pbc_corner_south_west_front(f13B)
  call apply_pbc_corner_south_west_front(f14B)
  call apply_pbc_corner_south_west_front(f15B)
  call apply_pbc_corner_south_west_front(f16B)
  call apply_pbc_corner_south_west_front(f17B)
  call apply_pbc_corner_south_west_front(f18B)
  
  !apply pbc at south west rear 7
  !red fluid
  call apply_pbc_corner_south_west_rear(f00R)
  call apply_pbc_corner_south_west_rear(f01R)
  call apply_pbc_corner_south_west_rear(f02R)
  call apply_pbc_corner_south_west_rear(f03R)
  call apply_pbc_corner_south_west_rear(f04R)
  call apply_pbc_corner_south_west_rear(f05R)
  call apply_pbc_corner_south_west_rear(f06R)
  call apply_pbc_corner_south_west_rear(f07R)
  call apply_pbc_corner_south_west_rear(f08R)
  call apply_pbc_corner_south_west_rear(f09R)
  call apply_pbc_corner_south_west_rear(f10R)
  call apply_pbc_corner_south_west_rear(f11R)
  call apply_pbc_corner_south_west_rear(f12R)
  call apply_pbc_corner_south_west_rear(f13R)
  call apply_pbc_corner_south_west_rear(f14R)
  call apply_pbc_corner_south_west_rear(f15R)
  call apply_pbc_corner_south_west_rear(f16R)
  call apply_pbc_corner_south_west_rear(f17R)
  call apply_pbc_corner_south_west_rear(f18R)
  !blue fluid
  call apply_pbc_corner_south_west_rear(f00B)
  call apply_pbc_corner_south_west_rear(f01B)
  call apply_pbc_corner_south_west_rear(f02B)
  call apply_pbc_corner_south_west_rear(f03B)
  call apply_pbc_corner_south_west_rear(f04B)
  call apply_pbc_corner_south_west_rear(f05B)
  call apply_pbc_corner_south_west_rear(f06B)
  call apply_pbc_corner_south_west_rear(f07B)
  call apply_pbc_corner_south_west_rear(f08B)
  call apply_pbc_corner_south_west_rear(f09B)
  call apply_pbc_corner_south_west_rear(f10B)
  call apply_pbc_corner_south_west_rear(f11B)
  call apply_pbc_corner_south_west_rear(f12B)
  call apply_pbc_corner_south_west_rear(f13B)
  call apply_pbc_corner_south_west_rear(f14B)
  call apply_pbc_corner_south_west_rear(f15B)
  call apply_pbc_corner_south_west_rear(f16B)
  call apply_pbc_corner_south_west_rear(f17B)
  call apply_pbc_corner_south_west_rear(f18B)
  
  !apply pbc at south east rear 8
  !red fluid
  call apply_pbc_corner_south_east_rear(f00R)
  call apply_pbc_corner_south_east_rear(f01R)
  call apply_pbc_corner_south_east_rear(f02R)
  call apply_pbc_corner_south_east_rear(f03R)
  call apply_pbc_corner_south_east_rear(f04R)
  call apply_pbc_corner_south_east_rear(f05R)
  call apply_pbc_corner_south_east_rear(f06R)
  call apply_pbc_corner_south_east_rear(f07R)
  call apply_pbc_corner_south_east_rear(f08R)
  call apply_pbc_corner_south_east_rear(f09R)
  call apply_pbc_corner_south_east_rear(f10R)
  call apply_pbc_corner_south_east_rear(f11R)
  call apply_pbc_corner_south_east_rear(f12R)
  call apply_pbc_corner_south_east_rear(f13R)
  call apply_pbc_corner_south_east_rear(f14R)
  call apply_pbc_corner_south_east_rear(f15R)
  call apply_pbc_corner_south_east_rear(f16R)
  call apply_pbc_corner_south_east_rear(f17R)
  call apply_pbc_corner_south_east_rear(f18R)
  !blue fluid
  call apply_pbc_corner_south_east_rear(f00B)
  call apply_pbc_corner_south_east_rear(f01B)
  call apply_pbc_corner_south_east_rear(f02B)
  call apply_pbc_corner_south_east_rear(f03B)
  call apply_pbc_corner_south_east_rear(f04B)
  call apply_pbc_corner_south_east_rear(f05B)
  call apply_pbc_corner_south_east_rear(f06B)
  call apply_pbc_corner_south_east_rear(f07B)
  call apply_pbc_corner_south_east_rear(f08B)
  call apply_pbc_corner_south_east_rear(f09B)
  call apply_pbc_corner_south_east_rear(f10B)
  call apply_pbc_corner_south_east_rear(f11B)
  call apply_pbc_corner_south_east_rear(f12B)
  call apply_pbc_corner_south_east_rear(f13B)
  call apply_pbc_corner_south_east_rear(f14B)
  call apply_pbc_corner_south_east_rear(f15B)
  call apply_pbc_corner_south_east_rear(f16B)
  call apply_pbc_corner_south_east_rear(f17B)
  call apply_pbc_corner_south_east_rear(f18B)
  
  return
  
 end subroutine apply_pbc_pops
 
 subroutine apply_pbc_east(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the periodic boundary 
!     condition at the east side
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk
  
  forall(kk=1:nbuff,j=1:ny,k=1:nz)
    dtemp(nx+kk,j,k)=dtemp(kk,j,k)
  end forall
  
  return
  
 end subroutine apply_pbc_east
 
 subroutine apply_pbc_west(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the periodic boundary 
!     condition at the west side
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk

  forall(kk=1:nbuff,j=1:ny,k=1:nz)
    dtemp(1-kk,j,k)=dtemp(nx+1-kk,j,k)
  end forall
  
  return
  
 end subroutine apply_pbc_west
 
 subroutine apply_pbc_north(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the periodic boundary 
!     condition at the north side
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk
  
  forall(i=1:nx,j=1:ny,kk=1:nbuff)
    dtemp(i,j,nz+kk)= dtemp(i,j,kk)
  end forall
  
  return
  
 end subroutine apply_pbc_north
 
 subroutine apply_pbc_south(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the periodic boundary 
!     condition at the south side
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk
  
  forall(i=1:nx,j=1:ny,kk=1:nbuff)
    dtemp(i,j,1-kk)=dtemp(i,j,nz+1-kk)
  end forall
  
  return
  
 end subroutine apply_pbc_south
 
 subroutine apply_pbc_front(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the periodic boundary 
!     condition at the front side
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk
  
  forall(i=1:nx,kk=1:nbuff,k=1:nz)
    dtemp(i,1-kk,k)= dtemp(i,ny+1-kk,k)
  end forall

  return
  
 end subroutine apply_pbc_front
 
 subroutine apply_pbc_rear(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the periodic boundary 
!     condition at the rear side
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk
  
  forall(i=1:nx,kk=1:nbuff,k=1:nz)
    dtemp(i,ny+kk,k)= dtemp(i,kk,k)
  end forall

  return
  
 end subroutine apply_pbc_rear
 
 subroutine apply_pbc_edge_front_east(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the periodic boundary 
!     condition at the front east edge
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk,kkk
  
  forall(kkk=1:nbuff,kk=1:nbuff,k=1:nz)
    dtemp(nx+kkk,1-kk,k)=dtemp(kkk,ny+1-kk,k)
  end forall

  return
  
 end subroutine apply_pbc_edge_front_east
 
 subroutine apply_pbc_edge_front_west(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the periodic boundary 
!     condition at the front west edge
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk,kkk
  
  forall(kkk=1:nbuff,kk=1:nbuff,k=1:nz)
    dtemp(1-kkk,1-kk,k)=dtemp(nx+1-kkk,ny+1-kk,k)
  end forall

  return
  
 end subroutine apply_pbc_edge_front_west
 
 subroutine apply_pbc_edge_north_east(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the periodic boundary 
!     condition at the north east edge
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk,kkk
  
  forall(kkk=1:nbuff,j=1:ny,kk=1:nbuff)
    dtemp(nx+kkk,j,nz+kk)=dtemp(kkk,j,kk)
  end forall

  return
  
 end subroutine apply_pbc_edge_north_east
 
 subroutine apply_pbc_edge_north_front(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the periodic boundary 
!     condition at the north front edge
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk,kkk
  
  forall(i=1:nx,kkk=1:nbuff,kk=1:nbuff)
    dtemp(i,1-kkk,nz+kk)=dtemp(i,ny+1-kkk,kk)
  end forall

  return
  
 end subroutine apply_pbc_edge_north_front
 
 subroutine apply_pbc_edge_north_rear(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the periodic boundary 
!     condition at the north rear edge
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk,kkk
  
  forall(i=1:nx,kkk=1:nbuff,kk=1:nbuff)
    dtemp(i,ny+kkk,nz+kk)= dtemp(i,kkk,kk)
  end forall

  return
  
 end subroutine apply_pbc_edge_north_rear
 
 subroutine apply_pbc_edge_north_west(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the periodic boundary 
!     condition at the north west edge
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk,kkk
  
  forall(kkk=1:nbuff,j=1:ny,kk=1:nbuff)
    dtemp(1-kkk,j,nz+kk)= dtemp(nx+1-kkk,j,kk)
  end forall
  
  return
  
 end subroutine apply_pbc_edge_north_west
 
 subroutine apply_pbc_edge_rear_east(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the periodic boundary 
!     condition at the rear east edge
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk,kkk
  
  forall(kkk=1:nbuff,kk=1:nbuff,k=1:nz)
    dtemp(nx+kkk,ny+kk,k)= dtemp(kkk,kk,k)
  end forall
  
  return
  
 end subroutine apply_pbc_edge_rear_east
 
 subroutine apply_pbc_edge_rear_west(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the periodic boundary 
!     condition at the rear west edge
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk,kkk
  
  forall(kkk=1:nbuff,kk=1:nbuff,k=1:nz)
    dtemp(1-kkk,ny+kk,k)=dtemp(nx+1-kkk,kk,k)
  end forall
  
  return
  
 end subroutine apply_pbc_edge_rear_west
 
 subroutine apply_pbc_edge_south_east(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the periodic boundary 
!     condition at the south east edge
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk,kkk
  
  forall(kkk=1:nbuff,j=1:ny,kk=1:nbuff)
    dtemp(nx+kkk,j,1-kk)=dtemp(kkk,j,nz+1-kk)
  end forall
  
  return
  
 end subroutine apply_pbc_edge_south_east
 
 subroutine apply_pbc_edge_south_front(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the periodic boundary 
!     condition at the south east edge
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk,kkk
  
  forall(i=1:nx,kkk=1:nbuff,kk=1:nbuff)
    dtemp(i,1-kkk,1-kk)=dtemp(i,ny+1-kkk,nz+1-kk)
  end forall
  
  return
  
 end subroutine apply_pbc_edge_south_front
 
 subroutine apply_pbc_edge_south_rear(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the periodic boundary 
!     condition at the south rear edge
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk,kkk
  
  forall(i=1:nx,kkk=1:nbuff,kk=1:nbuff)
    dtemp(i,ny+kkk,1-kk)=dtemp(i,kkk,nz+1-kk)
  end forall
  
  return
  
 end subroutine apply_pbc_edge_south_rear
 
 subroutine apply_pbc_edge_south_west(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the periodic boundary 
!     condition at the south rear edge
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk,kkk
  
  forall(kkk=1:nbuff,j=1:ny,kk=1:nbuff)
    dtemp(1-kkk,j,1-kk)=dtemp(nx+1-kkk,j,nz+1-kk)
  end forall
  
  return
  
 end subroutine apply_pbc_edge_south_west
 
 subroutine apply_pbc_corner_north_east_front(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the periodic boundary 
!     condition at the north east front corner
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk,kkk,kkkk
  
  forall(kk=1:nbuff,kkk=1:nbuff,kkkk=1:nbuff)
    dtemp(nx+kk,1-kkk,nz+kkkk)=dtemp(kk,ny+1-kkk,kkkk)
  end forall
  
  return
  
 end subroutine apply_pbc_corner_north_east_front
 
 subroutine apply_pbc_corner_north_east_rear(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the periodic boundary 
!     condition at the north east rear corner
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk,kkk,kkkk
  
  forall(kk=1:nbuff,kkk=1:nbuff,kkkk=1:nbuff)
    dtemp(nx+kk,ny+kkk,nz+kkkk)=dtemp(kk,kkk,kkkk)
  end forall
  
  return
  
 end subroutine apply_pbc_corner_north_east_rear
 
 subroutine apply_pbc_corner_north_west_rear(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the periodic boundary 
!     condition at the north west rear corner
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk,kkk,kkkk
  
  forall(kk=1:nbuff,kkk=1:nbuff,kkkk=1:nbuff)
    dtemp(1-kk,ny+kkk,nz+kkkk)=dtemp(nx+1-kk,kkk,kkkk)
  end forall
  
  return
  
 end subroutine apply_pbc_corner_north_west_rear
 
 subroutine apply_pbc_corner_north_west_front(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the periodic boundary 
!     condition at the north west front corner
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk,kkk,kkkk
  
  forall(kk=1:nbuff,kkk=1:nbuff,kkkk=1:nbuff)
    dtemp(1-kk,1-kkk,nz+kkkk)=dtemp(nx+1-kk,ny+1-kkk,kkkk)
  end forall
  
  return
  
 end subroutine apply_pbc_corner_north_west_front
 
 subroutine apply_pbc_corner_south_east_front(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the periodic boundary 
!     condition at the south east front corner
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk,kkk,kkkk
  
  forall(kk=1:nbuff,kkk=1:nbuff,kkkk=1:nbuff)
    dtemp(nx+kk,1-kkk,1-kkkk)=dtemp(kk,ny+1-kkk,nz+1-kkkk)
  end forall
  
  return
  
 end subroutine apply_pbc_corner_south_east_front
 
 subroutine apply_pbc_corner_south_west_front(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the periodic boundary 
!     condition at the south west front corner
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk,kkk,kkkk
  
  forall(kk=1:nbuff,kkk=1:nbuff,kkkk=1:nbuff)
    dtemp(1-kk,1-kkk,1-kkkk)=dtemp(nx+1-kk,ny+1-kkk,nz+1-kkkk)
  end forall
  
  return
  
 end subroutine apply_pbc_corner_south_west_front
 
 subroutine apply_pbc_corner_south_west_rear(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the periodic boundary 
!     condition at the south west rear corner
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk,kkk,kkkk
  
  forall(kk=1:nbuff,kkk=1:nbuff,kkkk=1:nbuff)
    dtemp(1-kk,ny+kkk,1-kkkk)=dtemp(nx+1-kk,kkk,nz+1-kkkk)
  end forall
  
  return
  
 end subroutine apply_pbc_corner_south_west_rear
 
 subroutine apply_pbc_corner_south_east_rear(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the periodic boundary 
!     condition at the south east rear corner
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk,kkk,kkkk
  
  forall(kk=1:nbuff,kkk=1:nbuff,kkkk=1:nbuff)
    dtemp(nx+kk,ny+kkk,1-kkkk)=dtemp(kk,kkk,nz+1-kkkk)
  end forall
  
  return
  
 end subroutine apply_pbc_corner_south_east_rear
 
 pure function equil_pop00(myrho,myu,myv,myw)
 
!***********************************************************************
!     
!     LBsoft function for Boltzmann equilibrium distribution 
!     for the population zero
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: myrho,myu,myv,myw
  
  real(kind=PRC) :: equil_pop00
  
  real(kind=PRC) :: uv
  
  integer, parameter :: myl = 0
  
  real(kind=PRC), parameter :: mycssq = cssq
  real(kind=PRC), parameter :: myp = p(myl)
  
  real(kind=PRC), parameter :: mydex = dex(myl)
  real(kind=PRC), parameter :: mydey = dey(myl)
  real(kind=PRC), parameter :: mydez = dez(myl)
  
  equil_pop00=myrho*myp*(ONE-(HALF/mycssq)*(myu**TWO+myv**TWO+myw**TWO))
 
  return
  
 end function equil_pop00
 
 pure function equil_pop01(myrho,myu,myv,myw)
 
!***********************************************************************
!     
!     LBsoft function for Boltzmann equilibrium distribution 
!     for the population one
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: myrho,myu,myv,myw
  
  real(kind=PRC) :: equil_pop01
  
  real(kind=PRC) :: uv
  
  integer, parameter :: myl = 1
  
  real(kind=PRC), parameter :: mycssq = cssq
  real(kind=PRC), parameter :: myp = p(myl)
  
  real(kind=PRC), parameter :: mydex = dex(myl)
  real(kind=PRC), parameter :: mydey = dey(myl)
  real(kind=PRC), parameter :: mydez = dez(myl)
  
  !occhio gli zeri su dex dey dez andrebbero tolti con pazienza su tutti gli equil_pop#
  
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop01=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
 
  return
  
 end function equil_pop01
 
 pure function equil_pop02(myrho,myu,myv,myw)
 
!***********************************************************************
!     
!     LBsoft function for Boltzmann equilibrium distribution 
!     for the population two
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: myrho,myu,myv,myw
  
  real(kind=PRC) :: equil_pop02
  
  real(kind=PRC) :: uv
  
  integer, parameter :: myl = 2
  
  real(kind=PRC), parameter :: mycssq = cssq
  real(kind=PRC), parameter :: myp = p(myl)
  
  real(kind=PRC), parameter :: mydex = dex(myl)
  real(kind=PRC), parameter :: mydey = dey(myl)
  real(kind=PRC), parameter :: mydez = dez(myl)
  
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop02=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
 
  return
  
 end function equil_pop02
 
 pure function equil_pop03(myrho,myu,myv,myw)
 
!***********************************************************************
!     
!     LBsoft function for Boltzmann equilibrium distribution 
!     for the population three
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: myrho,myu,myv,myw
  
  real(kind=PRC) :: equil_pop03
  
  real(kind=PRC) :: uv
  
  integer, parameter :: myl = 3
  
  real(kind=PRC), parameter :: mycssq = cssq
  real(kind=PRC), parameter :: myp = p(myl)
  
  real(kind=PRC), parameter :: mydex = dex(myl)
  real(kind=PRC), parameter :: mydey = dey(myl)
  real(kind=PRC), parameter :: mydez = dez(myl)
  
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop03=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
 
  return
  
 end function equil_pop03
 
 pure function equil_pop04(myrho,myu,myv,myw)
 
!***********************************************************************
!     
!     LBsoft function for Boltzmann equilibrium distribution 
!     for the population four
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: myrho,myu,myv,myw
  
  real(kind=PRC) :: equil_pop04
  
  real(kind=PRC) :: uv
  
  integer, parameter :: myl = 4
  
  real(kind=PRC), parameter :: mycssq = cssq
  real(kind=PRC), parameter :: myp = p(myl)
  
  real(kind=PRC), parameter :: mydex = dex(myl)
  real(kind=PRC), parameter :: mydey = dey(myl)
  real(kind=PRC), parameter :: mydez = dez(myl)
  
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop04=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
 
  return
  
 end function equil_pop04
 
 pure function equil_pop05(myrho,myu,myv,myw)
 
!***********************************************************************
!     
!     LBsoft function for Boltzmann equilibrium distribution 
!     for the population five
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: myrho,myu,myv,myw
  
  real(kind=PRC) :: equil_pop05
  
  real(kind=PRC) :: uv
  
  integer, parameter :: myl = 5
  
  real(kind=PRC), parameter :: mycssq = cssq
  real(kind=PRC), parameter :: myp = p(myl)
  
  real(kind=PRC), parameter :: mydex = dex(myl)
  real(kind=PRC), parameter :: mydey = dey(myl)
  real(kind=PRC), parameter :: mydez = dez(myl)
  
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop05=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
 
  return
  
 end function equil_pop05
 
 pure function equil_pop06(myrho,myu,myv,myw)
 
!***********************************************************************
!     
!     LBsoft function for Boltzmann equilibrium distribution 
!     for the population six
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: myrho,myu,myv,myw
  
  real(kind=PRC) :: equil_pop06
  
  real(kind=PRC) :: uv
  
  integer, parameter :: myl = 6
  
  real(kind=PRC), parameter :: mycssq = cssq
  real(kind=PRC), parameter :: myp = p(myl)
  
  real(kind=PRC), parameter :: mydex = dex(myl)
  real(kind=PRC), parameter :: mydey = dey(myl)
  real(kind=PRC), parameter :: mydez = dez(myl)
  
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop06=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
 
  return
  
 end function equil_pop06
 
 pure function equil_pop07(myrho,myu,myv,myw)
 
!***********************************************************************
!     
!     LBsoft function for Boltzmann equilibrium distribution 
!     for the population seven
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: myrho,myu,myv,myw
  
  real(kind=PRC) :: equil_pop07
  
  real(kind=PRC) :: uv
  
  integer, parameter :: myl = 7
  
  real(kind=PRC), parameter :: mycssq = cssq
  real(kind=PRC), parameter :: myp = p(myl)
  
  real(kind=PRC), parameter :: mydex = dex(myl)
  real(kind=PRC), parameter :: mydey = dey(myl)
  real(kind=PRC), parameter :: mydez = dez(myl)
  
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop07=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
 
  return
  
 end function equil_pop07
 
 pure function equil_pop08(myrho,myu,myv,myw)
 
!***********************************************************************
!     
!     LBsoft function for Boltzmann equilibrium distribution 
!     for the population eight
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: myrho,myu,myv,myw
  
  real(kind=PRC) :: equil_pop08
  
  real(kind=PRC) :: uv
  
  integer, parameter :: myl = 8
  
  real(kind=PRC), parameter :: mycssq = cssq
  real(kind=PRC), parameter :: myp = p(myl)
  
  real(kind=PRC), parameter :: mydex = dex(myl)
  real(kind=PRC), parameter :: mydey = dey(myl)
  real(kind=PRC), parameter :: mydez = dez(myl)
  
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop08=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
 
  return
  
 end function equil_pop08
 
 pure function equil_pop09(myrho,myu,myv,myw)
 
!***********************************************************************
!     
!     LBsoft function for Boltzmann equilibrium distribution 
!     for the population nine
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: myrho,myu,myv,myw
  
  real(kind=PRC) :: equil_pop09
  
  real(kind=PRC) :: uv
  
  integer, parameter :: myl = 9
  
  real(kind=PRC), parameter :: mycssq = cssq
  real(kind=PRC), parameter :: myp = p(myl)
  
  real(kind=PRC), parameter :: mydex = dex(myl)
  real(kind=PRC), parameter :: mydey = dey(myl)
  real(kind=PRC), parameter :: mydez = dez(myl)
  
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop09=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
 
  return
  
 end function equil_pop09
 
 pure function equil_pop10(myrho,myu,myv,myw)
 
!***********************************************************************
!     
!     LBsoft function for Boltzmann equilibrium distribution 
!     for the population ten
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: myrho,myu,myv,myw
  
  real(kind=PRC) :: equil_pop10
  
  real(kind=PRC) :: uv
  
  integer, parameter :: myl = 10
  
  real(kind=PRC), parameter :: mycssq = cssq
  real(kind=PRC), parameter :: myp = p(myl)
  
  real(kind=PRC), parameter :: mydex = dex(myl)
  real(kind=PRC), parameter :: mydey = dey(myl)
  real(kind=PRC), parameter :: mydez = dez(myl)
  
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop10=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
 
  return
  
 end function equil_pop10
 
 pure function equil_pop11(myrho,myu,myv,myw)
 
!***********************************************************************
!     
!     LBsoft function for Boltzmann equilibrium distribution 
!     for the population eleven
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: myrho,myu,myv,myw
  
  real(kind=PRC) :: equil_pop11
  
  real(kind=PRC) :: uv
  
  integer, parameter :: myl = 11
  
  real(kind=PRC), parameter :: mycssq = cssq
  real(kind=PRC), parameter :: myp = p(myl)
  
  real(kind=PRC), parameter :: mydex = dex(myl)
  real(kind=PRC), parameter :: mydey = dey(myl)
  real(kind=PRC), parameter :: mydez = dez(myl)
  
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop11=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
 
  return
  
 end function equil_pop11
 
 pure function equil_pop12(myrho,myu,myv,myw)
 
!***********************************************************************
!     
!     LBsoft function for Boltzmann equilibrium distribution 
!     for the population twelve
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: myrho,myu,myv,myw
  
  real(kind=PRC) :: equil_pop12
  
  real(kind=PRC) :: uv
  
  integer, parameter :: myl = 12
  
  real(kind=PRC), parameter :: mycssq = cssq
  real(kind=PRC), parameter :: myp = p(myl)
  
  real(kind=PRC), parameter :: mydex = dex(myl)
  real(kind=PRC), parameter :: mydey = dey(myl)
  real(kind=PRC), parameter :: mydez = dez(myl)
  
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop12=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
 
  return
  
 end function equil_pop12
 
 pure function equil_pop13(myrho,myu,myv,myw)
 
!***********************************************************************
!     
!     LBsoft function for Boltzmann equilibrium distribution 
!     for the population thirteen
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: myrho,myu,myv,myw
  
  real(kind=PRC) :: equil_pop13
  
  real(kind=PRC) :: uv
  
  integer, parameter :: myl = 13
  
  real(kind=PRC), parameter :: mycssq = cssq
  real(kind=PRC), parameter :: myp = p(myl)
  
  real(kind=PRC), parameter :: mydex = dex(myl)
  real(kind=PRC), parameter :: mydey = dey(myl)
  real(kind=PRC), parameter :: mydez = dez(myl)
  
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop13=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
 
  return
  
 end function equil_pop13
 
 pure function equil_pop14(myrho,myu,myv,myw)
 
!***********************************************************************
!     
!     LBsoft function for Boltzmann equilibrium distribution 
!     for the population fourteen
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: myrho,myu,myv,myw
  
  real(kind=PRC) :: equil_pop14
  
  real(kind=PRC) :: uv
  
  integer, parameter :: myl = 14
  
  real(kind=PRC), parameter :: mycssq = cssq
  real(kind=PRC), parameter :: myp = p(myl)
  
  real(kind=PRC), parameter :: mydex = dex(myl)
  real(kind=PRC), parameter :: mydey = dey(myl)
  real(kind=PRC), parameter :: mydez = dez(myl)
  
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop14=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
 
  return
  
 end function equil_pop14
 
 pure function equil_pop15(myrho,myu,myv,myw)
 
!***********************************************************************
!     
!     LBsoft function for Boltzmann equilibrium distribution 
!     for the population fiveteen
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: myrho,myu,myv,myw
  
  real(kind=PRC) :: equil_pop15
  
  real(kind=PRC) :: uv
  
  integer, parameter :: myl = 15
  
  real(kind=PRC), parameter :: mycssq = cssq
  real(kind=PRC), parameter :: myp = p(myl)
  
  real(kind=PRC), parameter :: mydex = dex(myl)
  real(kind=PRC), parameter :: mydey = dey(myl)
  real(kind=PRC), parameter :: mydez = dez(myl)
  
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop15=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
 
  return
  
 end function equil_pop15
 
 pure function equil_pop16(myrho,myu,myv,myw)
 
!***********************************************************************
!     
!     LBsoft function for Boltzmann equilibrium distribution 
!     for the population sixteen
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: myrho,myu,myv,myw
  
  real(kind=PRC) :: equil_pop16
  
  real(kind=PRC) :: uv
  
  integer, parameter :: myl = 16
  
  real(kind=PRC), parameter :: mycssq = cssq
  real(kind=PRC), parameter :: myp = p(myl)
  
  real(kind=PRC), parameter :: mydex = dex(myl)
  real(kind=PRC), parameter :: mydey = dey(myl)
  real(kind=PRC), parameter :: mydez = dez(myl)
  
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop16=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
 
  return
  
 end function equil_pop16
 
 pure function equil_pop17(myrho,myu,myv,myw)
 
!***********************************************************************
!     
!     LBsoft function for Boltzmann equilibrium distribution 
!     for the population seventeen
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: myrho,myu,myv,myw
  
  real(kind=PRC) :: equil_pop17
  
  real(kind=PRC) :: uv
  
  integer, parameter :: myl = 17
  
  real(kind=PRC), parameter :: mycssq = cssq
  real(kind=PRC), parameter :: myp = p(myl)
  
  real(kind=PRC), parameter :: mydex = dex(myl)
  real(kind=PRC), parameter :: mydey = dey(myl)
  real(kind=PRC), parameter :: mydez = dez(myl)
  
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop17=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
 
  return
  
 end function equil_pop17
 
 pure function equil_pop18(myrho,myu,myv,myw)
 
!***********************************************************************
!     
!     LBsoft function for Boltzmann equilibrium distribution 
!     for the population eighteen
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: myrho,myu,myv,myw
  
  real(kind=PRC) :: equil_pop18
  
  real(kind=PRC) :: uv
  
  integer, parameter :: myl = 18
  
  real(kind=PRC), parameter :: mycssq = cssq
  real(kind=PRC), parameter :: myp = p(myl)
  
  real(kind=PRC), parameter :: mydex = dex(myl)
  real(kind=PRC), parameter :: mydey = dey(myl)
  real(kind=PRC), parameter :: mydez = dez(myl)
  
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop18=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
 
  return
  
 end function equil_pop18
 
 subroutine compute_fluid_force_sc
 
!***********************************************************************
!     
!     LBsoft subroutine for computing the Shan Chen pair interaction
!     forces
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer :: i,j,k
  
  !red fluid
  call compute_grad_on_lattice(rhoR,gradpsixR,gradpsiyR,gradpsizR)
  !blue fluid
  call compute_grad_on_lattice(rhoB,gradpsixB,gradpsiyB,gradpsizB)
  
  !red fluid
  forall(i=1:nx,j=1:ny,k=1:nz)
    fuR(i,j,k) = fuR(i,j,k) - pair_SC*rhoR(i,j,k)*gradpsixB(i,j,k)
  end forall
  forall(i=1:nx,j=1:ny,k=1:nz)
    fvR(i,j,k) = fvR(i,j,k) - pair_SC*rhoR(i,j,k)*gradpsiyB(i,j,k)
  end forall
  forall(i=1:nx,j=1:ny,k=1:nz)
    fwR(i,j,k) = fwR(i,j,k) - pair_SC*rhoR(i,j,k)*gradpsizB(i,j,k)
  end forall
  !blue fluid
  forall(i=1:nx,j=1:ny,k=1:nz)
    fuB(i,j,k) = fuB(i,j,k) - pair_SC*rhoB(i,j,k)*gradpsixR(i,j,k)
  end forall
  forall(i=1:nx,j=1:ny,k=1:nz)
    fvB(i,j,k) = fvB(i,j,k) - pair_SC*rhoB(i,j,k)*gradpsiyR(i,j,k)
  end forall
  forall(i=1:nx,j=1:ny,k=1:nz)
    fwB(i,j,k) = fwB(i,j,k) - pair_SC*rhoB(i,j,k)*gradpsizB(i,j,k)
  end forall
  
  return
  
 end subroutine compute_fluid_force_sc
 
 subroutine convert_fluid_force_to_velshifted
 
!***********************************************************************
!     
!     LBsoft subroutine for converting the forces to the correspondent
!     velocity shifted for the forces
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer :: i,j,k
  
  !red fluid
  forall(i=1:nx,j=1:ny,k=1:nz)
    fuR(i,j,k) = fuR(i,j,k)*t_LB / rhoR(i,j,k) + u(i,j,k)
  end forall
  forall(i=1:nx,j=1:ny,k=1:nz)
    fvR(i,j,k) = fvR(i,j,k)*t_LB / rhoR(i,j,k) + v(i,j,k)
  end forall
  forall(i=1:nx,j=1:ny,k=1:nz)
    fwR(i,j,k) = fwR(i,j,k)*t_LB / rhoR(i,j,k) + w(i,j,k)
  end forall
  !blue fluid
  forall(i=1:nx,j=1:ny,k=1:nz)
    fuB(i,j,k) = fuB(i,j,k)*t_LB / rhoB(i,j,k) + u(i,j,k)
  end forall
  forall(i=1:nx,j=1:ny,k=1:nz)
    fvB(i,j,k) = fvB(i,j,k)*t_LB / rhoB(i,j,k) + v(i,j,k)
  end forall
  forall(i=1:nx,j=1:ny,k=1:nz)
    fwB(i,j,k) = fwB(i,j,k)*t_LB / rhoB(i,j,k) + w(i,j,k)
  end forall
  
  return
  
 end subroutine convert_fluid_force_to_velshifted
 
 subroutine compute_grad_on_lattice(myarr,mygradx,mygrady,mygradz)
 
!***********************************************************************
!     
!     LBsoft subroutine for computing the gradient of a scalar
!     over the lattice
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC), allocatable, dimension(:,:,:), intent(in) :: myarr
  real(kind=PRC), allocatable, dimension(:,:,:), intent(inout) :: &
   mygradx,mygrady,mygradz
   
  integer :: i,j,k,l
  
#if LATTICE==319
  !occhio gli zeri su dex dey dez andrebbero tolti con pazienza
  forall(i=1:nx,j=1:ny,k=1:nz)
    mygradx(i,j,k)= &
     myarr(i+ex(0),j+ey(0),k+ez(0))*p(0)*dex(0)+ & !00
     myarr(i+ex(1),j+ey(1),k+ez(1))*p(1)*dex(1)+ & !01
     myarr(i+ex(2),j+ey(2),k+ez(2))*p(2)*dex(2)+ & !02
     myarr(i+ex(3),j+ey(3),k+ez(3))*p(3)*dex(3)+ & !03
     myarr(i+ex(4),j+ey(4),k+ez(4))*p(4)*dex(4)+ & !04
     myarr(i+ex(5),j+ey(5),k+ez(5))*p(5)*dex(5)+ & !05
     myarr(i+ex(6),j+ey(6),k+ez(6))*p(6)*dex(6)+ & !06
     myarr(i+ex(7),j+ey(7),k+ez(7))*p(7)*dex(7)+ & !07
     myarr(i+ex(8),j+ey(8),k+ez(8))*p(8)*dex(8)+ & !08
     myarr(i+ex(9),j+ey(9),k+ez(9))*p(9)*dex(9)+ & !09
     myarr(i+ex(10),j+ey(10),k+ez(10))*p(10)*dex(10)+ & !10
     myarr(i+ex(11),j+ey(11),k+ez(11))*p(11)*dex(11)+ & !11
     myarr(i+ex(12),j+ey(12),k+ez(12))*p(12)*dex(12)+ & !12
     myarr(i+ex(13),j+ey(13),k+ez(13))*p(13)*dex(13)+ & !13
     myarr(i+ex(14),j+ey(14),k+ez(14))*p(14)*dex(14)+ & !14
     myarr(i+ex(15),j+ey(15),k+ez(15))*p(15)*dex(15)+ & !15
     myarr(i+ex(16),j+ey(16),k+ez(16))*p(16)*dex(16)+ & !16
     myarr(i+ex(17),j+ey(17),k+ez(17))*p(17)*dex(17)+ & !17
     myarr(i+ex(18),j+ey(18),k+ez(18))*p(18)*dex(18)    !18
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    mygrady(i,j,k)= &
     myarr(i+ex(0),j+ey(0),k+ez(0))*p(0)*dey(0)+ & !00
     myarr(i+ex(1),j+ey(1),k+ez(1))*p(1)*dey(1)+ & !01
     myarr(i+ex(2),j+ey(2),k+ez(2))*p(2)*dey(2)+ & !02
     myarr(i+ex(3),j+ey(3),k+ez(3))*p(3)*dey(3)+ & !03
     myarr(i+ex(4),j+ey(4),k+ez(4))*p(4)*dey(4)+ & !04
     myarr(i+ex(5),j+ey(5),k+ez(5))*p(5)*dey(5)+ & !05
     myarr(i+ex(6),j+ey(6),k+ez(6))*p(6)*dey(6)+ & !06
     myarr(i+ex(7),j+ey(7),k+ez(7))*p(7)*dey(7)+ & !07
     myarr(i+ex(8),j+ey(8),k+ez(8))*p(8)*dey(8)+ & !08
     myarr(i+ex(9),j+ey(9),k+ez(9))*p(9)*dey(9)+ & !09
     myarr(i+ex(10),j+ey(10),k+ez(10))*p(10)*dey(10)+ & !10
     myarr(i+ex(11),j+ey(11),k+ez(11))*p(11)*dey(11)+ & !11
     myarr(i+ex(12),j+ey(12),k+ez(12))*p(12)*dey(12)+ & !12
     myarr(i+ex(13),j+ey(13),k+ez(13))*p(13)*dey(13)+ & !13
     myarr(i+ex(14),j+ey(14),k+ez(14))*p(14)*dey(14)+ & !14
     myarr(i+ex(15),j+ey(15),k+ez(15))*p(15)*dey(15)+ & !15
     myarr(i+ex(16),j+ey(16),k+ez(16))*p(16)*dey(16)+ & !16
     myarr(i+ex(17),j+ey(17),k+ez(17))*p(17)*dey(17)+ & !17
     myarr(i+ex(18),j+ey(18),k+ez(18))*p(18)*dey(18)    !18
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    mygradz(i,j,k)= &
     myarr(i+ex(0),j+ey(0),k+ez(0))*p(0)*dez(0)+ & !00
     myarr(i+ex(1),j+ey(1),k+ez(1))*p(1)*dez(1)+ & !01
     myarr(i+ex(2),j+ey(2),k+ez(2))*p(2)*dez(2)+ & !02
     myarr(i+ex(3),j+ey(3),k+ez(3))*p(3)*dez(3)+ & !03
     myarr(i+ex(4),j+ey(4),k+ez(4))*p(4)*dez(4)+ & !04
     myarr(i+ex(5),j+ey(5),k+ez(5))*p(5)*dez(5)+ & !05
     myarr(i+ex(6),j+ey(6),k+ez(6))*p(6)*dez(6)+ & !06
     myarr(i+ex(7),j+ey(7),k+ez(7))*p(7)*dez(7)+ & !07
     myarr(i+ex(8),j+ey(8),k+ez(8))*p(8)*dez(8)+ & !08
     myarr(i+ex(9),j+ey(9),k+ez(9))*p(9)*dez(9)+ & !09
     myarr(i+ex(10),j+ey(10),k+ez(10))*p(10)*dez(10)+ & !10
     myarr(i+ex(11),j+ey(11),k+ez(11))*p(11)*dez(11)+ & !11
     myarr(i+ex(12),j+ey(12),k+ez(12))*p(12)*dez(12)+ & !12
     myarr(i+ex(13),j+ey(13),k+ez(13))*p(13)*dez(13)+ & !13
     myarr(i+ex(14),j+ey(14),k+ez(14))*p(14)*dez(14)+ & !14
     myarr(i+ex(15),j+ey(15),k+ez(15))*p(15)*dez(15)+ & !15
     myarr(i+ex(16),j+ey(16),k+ez(16))*p(16)*dez(16)+ & !16
     myarr(i+ex(17),j+ey(17),k+ez(17))*p(17)*dez(17)+ & !17
     myarr(i+ex(18),j+ey(18),k+ez(18))*p(18)*dez(18)    !18
  end forall
  
#endif
  
  return
  
 end subroutine compute_grad_on_lattice
 
 subroutine collision_fluids
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the collisional step
!     on the Boltzmann populations 
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  integer :: i,j,k
  
  if(lforce_add)then
    
    call convert_fluid_force_to_velshifted
    
    !red fluid
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f00R(i,j,k)=(ONE-omega(i,j,k))*f00R(i,j,k)+ (omega(i,j,k)-ONE)* &
       equil_pop00(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))+ &
       equil_pop00(rhoR(i,j,k),fuR(i,j,k),fvR(i,j,k),fwR(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f01R(i,j,k)=(ONE-omega(i,j,k))*f01R(i,j,k)+ (omega(i,j,k)-ONE)* &
       equil_pop01(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))+ &
       equil_pop01(rhoR(i,j,k),fuR(i,j,k),fvR(i,j,k),fwR(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f02R(i,j,k)=(ONE-omega(i,j,k))*f02R(i,j,k)+ (omega(i,j,k)-ONE)* &
       equil_pop02(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))+ &
       equil_pop02(rhoR(i,j,k),fuR(i,j,k),fvR(i,j,k),fwR(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f03R(i,j,k)=(ONE-omega(i,j,k))*f03R(i,j,k)+ (omega(i,j,k)-ONE)* &
       equil_pop03(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))+ &
       equil_pop03(rhoR(i,j,k),fuR(i,j,k),fvR(i,j,k),fwR(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f04R(i,j,k)=(ONE-omega(i,j,k))*f04R(i,j,k)+ (omega(i,j,k)-ONE)* &
       equil_pop04(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))+ &
       equil_pop04(rhoR(i,j,k),fuR(i,j,k),fvR(i,j,k),fwR(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f05R(i,j,k)=(ONE-omega(i,j,k))*f05R(i,j,k)+ (omega(i,j,k)-ONE)* &
       equil_pop05(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))+ &
       equil_pop05(rhoR(i,j,k),fuR(i,j,k),fvR(i,j,k),fwR(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f06R(i,j,k)=(ONE-omega(i,j,k))*f06R(i,j,k)+ (omega(i,j,k)-ONE)* &
       equil_pop06(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))+ &
       equil_pop06(rhoR(i,j,k),fuR(i,j,k),fvR(i,j,k),fwR(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f07R(i,j,k)=(ONE-omega(i,j,k))*f07R(i,j,k)+ (omega(i,j,k)-ONE)* &
       equil_pop07(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))+ &
       equil_pop07(rhoR(i,j,k),fuR(i,j,k),fvR(i,j,k),fwR(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f08R(i,j,k)=(ONE-omega(i,j,k))*f08R(i,j,k)+ (omega(i,j,k)-ONE)* &
       equil_pop08(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))+ &
       equil_pop08(rhoR(i,j,k),fuR(i,j,k),fvR(i,j,k),fwR(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f09R(i,j,k)=(ONE-omega(i,j,k))*f09R(i,j,k)+ (omega(i,j,k)-ONE)* &
       equil_pop09(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))+ &
       equil_pop09(rhoR(i,j,k),fuR(i,j,k),fvR(i,j,k),fwR(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f10R(i,j,k)=(ONE-omega(i,j,k))*f10R(i,j,k)+ (omega(i,j,k)-ONE)* &
       equil_pop10(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))+ &
       equil_pop10(rhoR(i,j,k),fuR(i,j,k),fvR(i,j,k),fwR(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f11R(i,j,k)=(ONE-omega(i,j,k))*f11R(i,j,k)+ (omega(i,j,k)-ONE)* &
       equil_pop11(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))+ &
       equil_pop11(rhoR(i,j,k),fuR(i,j,k),fvR(i,j,k),fwR(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f12R(i,j,k)=(ONE-omega(i,j,k))*f12R(i,j,k)+ (omega(i,j,k)-ONE)* &
       equil_pop12(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))+ &
       equil_pop12(rhoR(i,j,k),fuR(i,j,k),fvR(i,j,k),fwR(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f13R(i,j,k)=(ONE-omega(i,j,k))*f13R(i,j,k)+ (omega(i,j,k)-ONE)* &
       equil_pop13(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))+ &
       equil_pop13(rhoR(i,j,k),fuR(i,j,k),fvR(i,j,k),fwR(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f14R(i,j,k)=(ONE-omega(i,j,k))*f14R(i,j,k)+ (omega(i,j,k)-ONE)* &
       equil_pop14(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))+ &
       equil_pop14(rhoR(i,j,k),fuR(i,j,k),fvR(i,j,k),fwR(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f15R(i,j,k)=(ONE-omega(i,j,k))*f15R(i,j,k)+ (omega(i,j,k)-ONE)* &
       equil_pop15(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))+ &
       equil_pop15(rhoR(i,j,k),fuR(i,j,k),fvR(i,j,k),fwR(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f16R(i,j,k)=(ONE-omega(i,j,k))*f16R(i,j,k)+ (omega(i,j,k)-ONE)* &
       equil_pop16(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))+ &
       equil_pop16(rhoR(i,j,k),fuR(i,j,k),fvR(i,j,k),fwR(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f17R(i,j,k)=(ONE-omega(i,j,k))*f17R(i,j,k)+ (omega(i,j,k)-ONE)* &
       equil_pop17(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))+ &
       equil_pop17(rhoR(i,j,k),fuR(i,j,k),fvR(i,j,k),fwR(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f18R(i,j,k)=(ONE-omega(i,j,k))*f18R(i,j,k)+ (omega(i,j,k)-ONE)* &
       equil_pop18(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))+ &
       equil_pop18(rhoR(i,j,k),fuR(i,j,k),fvR(i,j,k),fwR(i,j,k))
    end forall
    
    !blue fluid
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f00B(i,j,k)=(ONE-omega(i,j,k))*f00B(i,j,k)+ (omega(i,j,k)-ONE)* &
       equil_pop00(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))+ &
       equil_pop00(rhoB(i,j,k),fuB(i,j,k),fvB(i,j,k),fwB(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f01B(i,j,k)=(ONE-omega(i,j,k))*f01B(i,j,k)+ (omega(i,j,k)-ONE)* &
       equil_pop01(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))+ &
       equil_pop01(rhoB(i,j,k),fuB(i,j,k),fvB(i,j,k),fwB(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f02B(i,j,k)=(ONE-omega(i,j,k))*f02B(i,j,k)+ (omega(i,j,k)-ONE)* &
       equil_pop02(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))+ &
       equil_pop02(rhoB(i,j,k),fuB(i,j,k),fvB(i,j,k),fwB(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f03B(i,j,k)=(ONE-omega(i,j,k))*f03B(i,j,k)+ (omega(i,j,k)-ONE)* &
       equil_pop03(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))+ &
       equil_pop03(rhoB(i,j,k),fuB(i,j,k),fvB(i,j,k),fwB(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f04B(i,j,k)=(ONE-omega(i,j,k))*f04B(i,j,k)+ (omega(i,j,k)-ONE)* &
       equil_pop04(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))+ &
       equil_pop04(rhoB(i,j,k),fuB(i,j,k),fvB(i,j,k),fwB(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f05B(i,j,k)=(ONE-omega(i,j,k))*f05B(i,j,k)+ (omega(i,j,k)-ONE)* &
       equil_pop05(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))+ &
       equil_pop05(rhoB(i,j,k),fuB(i,j,k),fvB(i,j,k),fwB(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f06B(i,j,k)=(ONE-omega(i,j,k))*f06B(i,j,k)+ (omega(i,j,k)-ONE)* &
       equil_pop06(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))+ &
       equil_pop06(rhoB(i,j,k),fuB(i,j,k),fvB(i,j,k),fwB(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f07B(i,j,k)=(ONE-omega(i,j,k))*f07B(i,j,k)+ (omega(i,j,k)-ONE)* &
       equil_pop07(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))+ &
       equil_pop07(rhoB(i,j,k),fuB(i,j,k),fvB(i,j,k),fwB(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f08B(i,j,k)=(ONE-omega(i,j,k))*f08B(i,j,k)+ (omega(i,j,k)-ONE)* &
       equil_pop08(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))+ &
       equil_pop08(rhoB(i,j,k),fuB(i,j,k),fvB(i,j,k),fwB(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f09B(i,j,k)=(ONE-omega(i,j,k))*f09B(i,j,k)+ (omega(i,j,k)-ONE)* &
       equil_pop09(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))+ &
       equil_pop09(rhoB(i,j,k),fuB(i,j,k),fvB(i,j,k),fwB(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f10B(i,j,k)=(ONE-omega(i,j,k))*f10B(i,j,k)+ (omega(i,j,k)-ONE)* &
       equil_pop10(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))+ &
       equil_pop10(rhoB(i,j,k),fuB(i,j,k),fvB(i,j,k),fwB(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f11B(i,j,k)=(ONE-omega(i,j,k))*f11B(i,j,k)+ (omega(i,j,k)-ONE)* &
       equil_pop11(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))+ &
       equil_pop11(rhoB(i,j,k),fuB(i,j,k),fvB(i,j,k),fwB(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f12B(i,j,k)=(ONE-omega(i,j,k))*f12B(i,j,k)+ (omega(i,j,k)-ONE)* &
       equil_pop12(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))+ &
       equil_pop12(rhoB(i,j,k),fuB(i,j,k),fvB(i,j,k),fwB(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f13B(i,j,k)=(ONE-omega(i,j,k))*f13B(i,j,k)+ (omega(i,j,k)-ONE)* &
       equil_pop13(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))+ &
       equil_pop13(rhoB(i,j,k),fuB(i,j,k),fvB(i,j,k),fwB(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f14B(i,j,k)=(ONE-omega(i,j,k))*f14B(i,j,k)+ (omega(i,j,k)-ONE)* &
       equil_pop14(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))+ &
       equil_pop14(rhoB(i,j,k),fuB(i,j,k),fvB(i,j,k),fwB(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f15B(i,j,k)=(ONE-omega(i,j,k))*f15B(i,j,k)+ (omega(i,j,k)-ONE)* &
       equil_pop15(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))+ &
       equil_pop15(rhoB(i,j,k),fuB(i,j,k),fvB(i,j,k),fwB(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f16B(i,j,k)=(ONE-omega(i,j,k))*f16B(i,j,k)+ (omega(i,j,k)-ONE)* &
       equil_pop16(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))+ &
       equil_pop16(rhoB(i,j,k),fuB(i,j,k),fvB(i,j,k),fwB(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f17B(i,j,k)=(ONE-omega(i,j,k))*f17B(i,j,k)+ (omega(i,j,k)-ONE)* &
       equil_pop17(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))+ &
       equil_pop17(rhoB(i,j,k),fuB(i,j,k),fvB(i,j,k),fwB(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f18B(i,j,k)=(ONE-omega(i,j,k))*f18B(i,j,k)+ (omega(i,j,k)-ONE)* &
       equil_pop18(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))+ &
       equil_pop18(rhoB(i,j,k),fuB(i,j,k),fvB(i,j,k),fwB(i,j,k))
    end forall
     
  else

    forall(i=1:nx,j=1:ny,k=1:nz)
      f00R(i,j,k)=f00R(i,j,k)+omega(i,j,k)* &
       (equil_pop00(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))-f00R(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f01R(i,j,k)=f01R(i,j,k)+omega(i,j,k)* &
       (equil_pop01(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))-f01R(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f02R(i,j,k)=f02R(i,j,k)+omega(i,j,k)* &
       (equil_pop02(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))-f02R(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f03R(i,j,k)=f03R(i,j,k)+omega(i,j,k)* &
       (equil_pop03(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))-f03R(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f04R(i,j,k)=f04R(i,j,k)+omega(i,j,k)* &
       (equil_pop04(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))-f04R(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f05R(i,j,k)=f05R(i,j,k)+omega(i,j,k)* &
       (equil_pop05(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))-f05R(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f06R(i,j,k)=f06R(i,j,k)+omega(i,j,k)* &
       (equil_pop06(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))-f06R(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f07R(i,j,k)=f07R(i,j,k)+omega(i,j,k)* &
       (equil_pop07(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))-f07R(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f08R(i,j,k)=f08R(i,j,k)+omega(i,j,k)* &
       (equil_pop08(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))-f08R(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f09R(i,j,k)=f09R(i,j,k)+omega(i,j,k)* &
       (equil_pop09(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))-f09R(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f10R(i,j,k)=f10R(i,j,k)+omega(i,j,k)* &
       (equil_pop10(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))-f10R(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f11R(i,j,k)=f11R(i,j,k)+omega(i,j,k)* &
       (equil_pop11(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))-f11R(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f12R(i,j,k)=f12R(i,j,k)+omega(i,j,k)* &
       (equil_pop12(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))-f12R(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f13R(i,j,k)=f13R(i,j,k)+omega(i,j,k)* &
       (equil_pop13(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))-f13R(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f14R(i,j,k)=f14R(i,j,k)+omega(i,j,k)* &
       (equil_pop14(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))-f14R(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f15R(i,j,k)=f15R(i,j,k)+omega(i,j,k)* &
       (equil_pop15(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))-f15R(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f16R(i,j,k)=f16R(i,j,k)+omega(i,j,k)* &
       (equil_pop16(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))-f16R(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f17R(i,j,k)=f17R(i,j,k)+omega(i,j,k)* &
       (equil_pop17(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))-f17R(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f18R(i,j,k)=f18R(i,j,k)+omega(i,j,k)* &
       (equil_pop18(rhoR(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))-f18R(i,j,k))
    end forall
    
    !blue fluid
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f00B(i,j,k)=f00B(i,j,k)+omega(i,j,k)* &
       (equil_pop00(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))-f00B(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f01B(i,j,k)=f01B(i,j,k)+omega(i,j,k)* &
       (equil_pop01(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))-f01B(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f02B(i,j,k)=f02B(i,j,k)+omega(i,j,k)* &
       (equil_pop02(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))-f02B(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f03B(i,j,k)=f03B(i,j,k)+omega(i,j,k)* &
       (equil_pop03(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))-f03B(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f04B(i,j,k)=f04B(i,j,k)+omega(i,j,k)* &
       (equil_pop04(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))-f04B(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f05B(i,j,k)=f05B(i,j,k)+omega(i,j,k)* &
       (equil_pop05(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))-f05B(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f06B(i,j,k)=f06B(i,j,k)+omega(i,j,k)* &
       (equil_pop06(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))-f06B(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f07B(i,j,k)=f07B(i,j,k)+omega(i,j,k)* &
       (equil_pop07(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))-f07B(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f08B(i,j,k)=f08B(i,j,k)+omega(i,j,k)* &
       (equil_pop08(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))-f08B(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f09B(i,j,k)=f09B(i,j,k)+omega(i,j,k)* &
       (equil_pop09(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))-f09B(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f10B(i,j,k)=f10B(i,j,k)+omega(i,j,k)* &
       (equil_pop10(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))-f10B(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f11B(i,j,k)=f11B(i,j,k)+omega(i,j,k)* &
       (equil_pop11(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))-f11B(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f12B(i,j,k)=f12B(i,j,k)+omega(i,j,k)* &
       (equil_pop12(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))-f12B(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f13B(i,j,k)=f13B(i,j,k)+omega(i,j,k)* &
       (equil_pop13(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))-f13B(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f14B(i,j,k)=f14B(i,j,k)+omega(i,j,k)* &
       (equil_pop14(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))-f14B(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f15B(i,j,k)=f15B(i,j,k)+omega(i,j,k)* &
       (equil_pop15(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))-f15B(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f16B(i,j,k)=f16B(i,j,k)+omega(i,j,k)* &
       (equil_pop16(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))-f16B(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f17B(i,j,k)=f17B(i,j,k)+omega(i,j,k)* &
       (equil_pop17(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))-f17B(i,j,k))
    end forall
    
    forall(i=1:nx,j=1:ny,k=1:nz)
      f18B(i,j,k)=f18B(i,j,k)+omega(i,j,k)* &
       (equil_pop18(rhoB(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k))-f18B(i,j,k))
    end forall
  
  endif
  
  return
  
 end subroutine collision_fluids
 
 end module fluids_mod
