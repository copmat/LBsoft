
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
                   buffservice3d,allocate_array_buffservice3d, &
                   rand_noseeded,linit_seed,gauss_noseeded
 
 implicit none
 
 private
 
 integer, save, protected, public :: nx=0
 integer, save, protected, public :: ny=0
 integer, save, protected, public :: nz=0
 
 integer, save, protected, public :: LBintegrator=0
 
 integer, save, protected, public :: idistselect=0
 
 integer, save, protected, public :: ibctype=0
 
 integer, save, protected, public :: bc_type_east=0
 integer, save, protected, public :: bc_type_west=0
 integer, save, protected, public :: bc_type_front=0
 integer, save, protected, public :: bc_type_rear=0
 integer, save, protected, public :: bc_type_north=0
 integer, save, protected, public :: bc_type_south=0
 
 logical, save, protected, public :: lalloc_hvars=.false.
 logical, save, protected, public :: lalloc_pops=.false.
 
 logical, save, protected, public :: lforce_add=.false.
 logical, save, protected, public :: lShanChen=.false.
 logical, save, protected, public :: lsing_SC=.false.
 logical, save, protected, public :: lpair_SC=.false.
 
 logical, save, protected, public :: lunique_omega=.false.
 
 logical, save, protected, public :: lsingle_fluid=.false.
 
 logical, save, protected, public :: lvorticity=.false.
 
 logical, save, protected, public :: lwall_SC=.false.
 
 real(kind=PRC), save, protected, public :: t_LB = ONE
 
 real(kind=PRC), save, protected, public :: beta         = ZERO
 real(kind=PRC), save, protected, public :: akl          = ZERO
 real(kind=PRC), save, protected, public :: awall        = ZERO
 real(kind=PRC), save, protected, public :: tauR         = ZERO
 real(kind=PRC), save, protected, public :: tauB         = ZERO
 real(kind=PRC), save, protected, public :: unique_omega = ZERO
 real(kind=PRC), save, protected, public :: viscR        = ZERO
 real(kind=PRC), save, protected, public :: viscB        = ZERO
 real(kind=PRC), save, protected, public :: meanR        = ZERO
 real(kind=PRC), save, protected, public :: meanB        = ZERO
 real(kind=PRC), save, protected, public :: stdevR       = ZERO
 real(kind=PRC), save, protected, public :: stdevB       = ZERO
 real(kind=PRC), save, protected, public :: initial_u    = ZERO
 real(kind=PRC), save, protected, public :: initial_v    = ZERO
 real(kind=PRC), save, protected, public :: initial_w    = ZERO
 real(kind=PRC), save, protected, public :: ext_fu       = ZERO
 real(kind=PRC), save, protected, public :: ext_fv       = ZERO
 real(kind=PRC), save, protected, public :: ext_fw       = ZERO
 
 real(kind=PRC), save, protected, public :: bc_rhoR_east = ZERO
 real(kind=PRC), save, protected, public :: bc_rhoR_west = ZERO
 real(kind=PRC), save, protected, public :: bc_rhoR_front= ZERO
 real(kind=PRC), save, protected, public :: bc_rhoR_rear = ZERO
 real(kind=PRC), save, protected, public :: bc_rhoR_north= ZERO
 real(kind=PRC), save, protected, public :: bc_rhoR_south= ZERO
 
 real(kind=PRC), save, protected, public :: bc_rhoB_east = ZERO
 real(kind=PRC), save, protected, public :: bc_rhoB_west = ZERO
 real(kind=PRC), save, protected, public :: bc_rhoB_front= ZERO
 real(kind=PRC), save, protected, public :: bc_rhoB_rear = ZERO
 real(kind=PRC), save, protected, public :: bc_rhoB_north= ZERO
 real(kind=PRC), save, protected, public :: bc_rhoB_south= ZERO
 
 real(kind=PRC), save, protected, public :: bc_u_east = ZERO
 real(kind=PRC), save, protected, public :: bc_u_west = ZERO
 real(kind=PRC), save, protected, public :: bc_u_front= ZERO
 real(kind=PRC), save, protected, public :: bc_u_rear = ZERO
 real(kind=PRC), save, protected, public :: bc_u_north= ZERO
 real(kind=PRC), save, protected, public :: bc_u_south= ZERO
 
 real(kind=PRC), save, protected, public :: bc_v_east = ZERO
 real(kind=PRC), save, protected, public :: bc_v_west = ZERO
 real(kind=PRC), save, protected, public :: bc_v_front= ZERO
 real(kind=PRC), save, protected, public :: bc_v_rear = ZERO
 real(kind=PRC), save, protected, public :: bc_v_north= ZERO
 real(kind=PRC), save, protected, public :: bc_v_south= ZERO
 
 real(kind=PRC), save, protected, public :: bc_w_east = ZERO
 real(kind=PRC), save, protected, public :: bc_w_west = ZERO
 real(kind=PRC), save, protected, public :: bc_w_front= ZERO
 real(kind=PRC), save, protected, public :: bc_w_rear = ZERO
 real(kind=PRC), save, protected, public :: bc_w_north= ZERO
 real(kind=PRC), save, protected, public :: bc_w_south= ZERO
 
 real(kind=PRC), save, protected, public :: pair_SC = ZERO
 
 real(kind=PRC), save, protected, public :: wallR_SC = ONE
 real(kind=PRC), save, protected, public :: wallB_SC = ONE
 
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
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: psiR
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: psiB
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: gradpsixR
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: gradpsiyR
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: gradpsizR
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: gradpsixB
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: gradpsiyB
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: gradpsizB
 
#if LATTICE==319
 
 character(len=6), parameter, public :: latt_name="d3q19 "
 
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
  ex = (/ 0, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1, 0, 0, 0, 0/)
   !      0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18
 integer, dimension(0:links), parameter, public :: &
  ey = (/ 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1/)
 integer, dimension(0:links), parameter, public :: &
  ez = (/ 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1/)
 integer, dimension(0:links), parameter, public :: &
  opp =(/ 0, 2, 1, 4, 3, 6, 5, 8, 7,10, 9,12,11,14,13,16,15,18,17/)
  
 real(kind=PRC), dimension(0:links), parameter, public :: &
  dex = real(ex,kind=PRC)
 real(kind=PRC), dimension(0:links), parameter, public :: &
  dey = real(ey,kind=PRC)
 real(kind=PRC), dimension(0:links), parameter, public :: &
  dez = real(ez,kind=PRC)
 
#endif 

 real(kind=PRC), save, public, allocatable, dimension(:,:,:) :: & 
  f00R,f01R,f02R,f03R,f04R,f05R,f06R,f07R,f08R
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: &
  f00B,f01B,f02B,f03B,f04B,f05B,f06B,f07B,f08B
 
 real(kind=PRC), save, public, allocatable, dimension(:,:,:) :: & 
  f09R,f10R,f11R,f12R,f13R,f14R,f15R,f16R,f17R,f18R
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: &
  f09B,f10B,f11B,f12B,f13B,f14B,f15B,f16B,f17B,f18B
 
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
 public :: driver_collision_fluids
 public :: driver_bc_densities
 public :: driver_bc_velocities
 public :: set_fluid_force_sc
 public :: set_value_viscosity
 public :: set_value_tau
 public :: compute_omega
 public :: driver_streaming_fluids
 public :: moments_fluids
 public :: driver_reflect_densities
 public :: set_lsingle_fluid
 public :: driver_apply_bounceback_pop
 public :: probe_red_moments_in_node
 public :: probe_blue_moments_in_node
 public :: set_value_bc_east
 public :: set_value_bc_west
 public :: set_value_bc_front
 public :: set_value_bc_rear
 public :: set_value_bc_north
 public :: set_value_bc_south
 public :: set_fluid_wall_sc
 public :: set_LBintegrator_type
 
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
  
  integer :: i,j,k
  logical, dimension(1) :: ltest=.false.
  
  myzero=1-nbuff
  mynx=nx+nbuff
  if(mod(nx,2)==0)mynx=mynx+1
  myny=ny+nbuff
  mynz=nz+nbuff
  
  istat=0
  
  allocate(rhoR(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(1))
  
  allocate(u(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(3))
  allocate(v(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(4))
  allocate(w(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(5))
  
  allocate(fuR(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(6))
  allocate(fvR(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(7))
  allocate(fwR(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(8))
  
  if(lShanChen)then
    allocate(gradpsixR(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(12))
    allocate(gradpsiyR(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(13))
    allocate(gradpsizR(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(14))
    allocate(psiR(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(19))
  endif
  
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
  
  if(.not. lsingle_fluid)then
  
    allocate(rhoB(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(2))
  
    allocate(fuB(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(9))
    allocate(fvB(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(10))
    allocate(fwB(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(11))
    
    if(lShanChen)then
      allocate(gradpsixB(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(15))
      allocate(gradpsiyB(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(16))
      allocate(gradpsizB(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(17))
      allocate(psiB(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(20))
    endif
    
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
  if(ltest(1))call error(4)
  
  call allocate_array_buffservice3d(myzero,mynx,myzero,myny,myzero,mynz)
  
  return
  
 end subroutine allocate_fluids
 
 subroutine set_boundary_conditions_type(itemp1,itemp2,itemp3)
 
  implicit none
  
  integer, intent(in) :: itemp1,itemp2,itemp3
  integer ::i,j,k
  integer, dimension(0:1,0:1,0:1), parameter :: iselbct = &
   reshape((/0,1,2,3,4,5,6,7/),(/2,2,2/))
 ! 0=F    1=T  periodic
 !     itemp1      itemp2      itemp3     ibctype
 !          0           0           0           0
 !          1           0           0           1
 !          0           1           0           2
 !          1           1           0           3
 !          0           0           1           4
 !          1           0           1           5
 !          0           1           1           6
 !          1           1           1           7
  
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
  case(2)
    call set_uniform_dens_fluids
  case default
    call set_initial_dens_fluids
  end select
  
  call set_initial_vel_fluids
  
  call driver_bc_densities
  
  if(lvorticity)call driver_bc_velocities
  
  call driver_set_initial_pop_fluids
  
  call driver_bc_pops
  
  call driver_reflect_pops
  
  if(lunique_omega)call set_unique_omega
  
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
  
  if(.not. lforce_add)return
  
  forall(i=1:nx,j=1:ny,k=1:nz)fuR(i,j,k)=ext_fu
  forall(i=1:nx,j=1:ny,k=1:nz)fvR(i,j,k)=ext_fv
  forall(i=1:nx,j=1:ny,k=1:nz)fwR(i,j,k)=ext_fw
  
  if(lsingle_fluid)return
  
  forall(i=1:nx,j=1:ny,k=1:nz)fuB(i,j,k)=ext_fu
  forall(i=1:nx,j=1:ny,k=1:nz)fvB(i,j,k)=ext_fv
  forall(i=1:nx,j=1:ny,k=1:nz)fwB(i,j,k)=ext_fw
  
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
  
  implicit none
  
  integer :: i,j,k
  
  if(linit_seed)then
  
    do k=1,nz
      do j=1,ny
        do i=1,nx
          rhoR(i,j,k)=meanR+stdevR*gauss()
        enddo
      enddo
    enddo
    
    if(lsingle_fluid)return
    
    do k=1,nz
      do j=1,ny
        do i=1,nx
          rhoB(i,j,k)=meanB+stdevB*gauss()
        enddo
      enddo
    enddo
  
  else
  
    forall (i=1:nx,j=1:ny,k=1:nz)
      rhoR(i,j,k)=meanR+stdevR*gauss_noseeded(i,j,k,1)
    end forall
    
    if(lsingle_fluid)return
    
    forall (i=1:nx,j=1:ny,k=1:nz)
      rhoB(i,j,k)=meanB+stdevB*gauss_noseeded(i,j,k,100)
    end forall
  
  endif
  
  return
  
 end subroutine set_random_dens_fluids
 
 subroutine set_uniform_dens_fluids
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the initial density of fluids 
!     following a random uniform distribution
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  integer :: i,j,k
  
  real(kind=PRC) :: dtemp1
  
  if(linit_seed)then
  
    do k=1,nz
      do j=1,ny
        do i=1,nx
          call random_number(dtemp1)
          rhoR(i,j,k)=meanR+stdevR*dtemp1
        enddo
      enddo
    enddo
    
    if(lsingle_fluid)return
    
    do k=1,nz
      do j=1,ny
        do i=1,nx
          call random_number(dtemp1)
          rhoB(i,j,k)=meanB+stdevB*dtemp1
        enddo
      enddo
    enddo
  
  else
  
    forall (i=1:nx,j=1:ny,k=1:nz)
      rhoR(i,j,k)=meanR+stdevR*rand_noseeded(i,j,k,1)
    end forall
    
    if(lsingle_fluid)return
    
    forall (i=1:nx,j=1:ny,k=1:nz)
      rhoB(i,j,k)=meanB+stdevB*rand_noseeded(i,j,k,100)
    end forall
  
  endif
  
  return
  
 end subroutine set_uniform_dens_fluids
 
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
  
  integer :: i,j,k
  
  forall(i=1:nx,j=1:ny,k=1:nz)rhoR(i,j,k)=meanR
  if(lsingle_fluid)return
  forall(i=1:nx,j=1:ny,k=1:nz)rhoB(i,j,k)=meanB
  
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
  
  integer :: i,j,k
  
  forall(i=1:nx,j=1:ny,k=1:nz)u(i,j,k)=initial_u
  forall(i=1:nx,j=1:ny,k=1:nz)v(i,j,k)=initial_v
  forall(i=1:nx,j=1:ny,k=1:nz)w(i,j,k)=initial_w
  
  return
  
 end subroutine set_initial_vel_fluids
 
 subroutine driver_set_initial_pop_fluids
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the setting of the initial 
!     fluid populations
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  call set_initial_pop_fluids(rhoR,u,v,w,f00R,f01R,&
  f02R,f03R,f04R,f05R,f06R,f07R,f08R,f09R,f10R, &
  f11R,f12R,f13R,f14R,f15R,f16R,f17R,f18R)
  
  if(lsingle_fluid)return
  
  call set_initial_pop_fluids(rhoB,u,v,w,f00B,f01B,&
  f02B,f03B,f04B,f05B,f06B,f07B,f08B,f09B,f10B, &
  f11B,f12B,f13B,f14B,f15B,f16B,f17B,f18B)
  
  return
  
 end subroutine driver_set_initial_pop_fluids
 
 subroutine set_initial_pop_fluids(rhosub,usub,vsub,wsub,f00sub,f01sub,&
  f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub,f09sub,f10sub, &
  f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub)
 
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
  
  real(kind=PRC), allocatable, dimension(:,:,:)  :: rhosub,usub,vsub, &
   wsub,f00sub,f01sub,f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub,&
   f09sub,f10sub,f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub
  
  integer :: i,j,k
  
  !red fluid
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f00sub(i,j,k)=equil_pop00(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k), &
     wsub(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f01sub(i,j,k)=equil_pop01(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k), &
     wsub(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f02sub(i,j,k)=equil_pop02(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k), &
     wsub(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f03sub(i,j,k)=equil_pop03(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k), &
     wsub(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f04sub(i,j,k)=equil_pop04(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k), &
     wsub(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f05sub(i,j,k)=equil_pop05(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k), &
     wsub(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f06sub(i,j,k)=equil_pop06(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k), &
     wsub(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f07sub(i,j,k)=equil_pop07(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k), &
     wsub(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f08sub(i,j,k)=equil_pop08(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k), &
     wsub(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f09sub(i,j,k)=equil_pop09(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k), &
     wsub(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f10sub(i,j,k)=equil_pop10(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k), &
     wsub(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f11sub(i,j,k)=equil_pop11(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k), &
     wsub(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f12sub(i,j,k)=equil_pop12(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k), &
     wsub(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f13sub(i,j,k)=equil_pop13(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k), &
     wsub(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f14sub(i,j,k)=equil_pop14(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k), &
     wsub(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f15sub(i,j,k)=equil_pop15(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k), &
     wsub(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f16sub(i,j,k)=equil_pop16(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k), &
     wsub(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f17sub(i,j,k)=equil_pop17(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k), &
     wsub(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f18sub(i,j,k)=equil_pop18(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k), &
     wsub(i,j,k))
  end forall
  
  return
  
 end subroutine set_initial_pop_fluids
 
 subroutine set_value_bc_east(itemp1,dtemp1,dtemp2,dtemp3,dtemp4,dtemp5)
 
!***********************************************************************
!     
!     LBsoft subroutine for set the conditions of fluids 
!     at the east open boundary
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: itemp1
  real(kind=PRC), intent(in) :: dtemp1,dtemp2,dtemp3,dtemp4,dtemp5
  
  bc_type_east = itemp1
  bc_rhoR_east = dtemp1
  bc_rhoB_east = dtemp2
  bc_u_east = dtemp3
  bc_v_east = dtemp4
  bc_w_east = dtemp5
  
  return
  
 end subroutine set_value_bc_east
 
 subroutine set_value_bc_west(itemp1,dtemp1,dtemp2,dtemp3,dtemp4,dtemp5)
 
!***********************************************************************
!     
!     LBsoft subroutine for set the conditions of fluids 
!     at the west open boundary
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: itemp1
  real(kind=PRC), intent(in) :: dtemp1,dtemp2,dtemp3,dtemp4,dtemp5
  
  bc_type_west = itemp1
  bc_rhoR_west = dtemp1
  bc_rhoB_west = dtemp2
  bc_u_west = dtemp3
  bc_v_west = dtemp4
  bc_w_west = dtemp5
  
  return
  
 end subroutine set_value_bc_west
 
 subroutine set_value_bc_front(itemp1,dtemp1,dtemp2,dtemp3,dtemp4,dtemp5)
 
!***********************************************************************
!     
!     LBsoft subroutine for set the conditions of fluids 
!     at the front open boundary
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: itemp1
  real(kind=PRC), intent(in) :: dtemp1,dtemp2,dtemp3,dtemp4,dtemp5
  
  bc_type_front = itemp1
  bc_rhoR_front = dtemp1
  bc_rhoB_front = dtemp2
  bc_u_front = dtemp3
  bc_v_front = dtemp4
  bc_w_front = dtemp5
  
  return
  
 end subroutine set_value_bc_front
 
 subroutine set_value_bc_rear(itemp1,dtemp1,dtemp2,dtemp3,dtemp4,dtemp5)
 
!***********************************************************************
!     
!     LBsoft subroutine for set the conditions of fluids 
!     at the rear open boundary
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: itemp1
  real(kind=PRC), intent(in) :: dtemp1,dtemp2,dtemp3,dtemp4,dtemp5
  
  bc_type_rear = itemp1
  bc_rhoR_rear = dtemp1
  bc_rhoB_rear = dtemp2
  bc_u_rear = dtemp3
  bc_v_rear = dtemp4
  bc_w_rear = dtemp5
  
  return
  
 end subroutine set_value_bc_rear
 
 subroutine set_value_bc_north(itemp1,dtemp1,dtemp2,dtemp3,dtemp4,dtemp5)
 
!***********************************************************************
!     
!     LBsoft subroutine for set the conditions of fluids 
!     at the north open boundary
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: itemp1
  real(kind=PRC), intent(in) :: dtemp1,dtemp2,dtemp3,dtemp4,dtemp5
  
  bc_type_north = itemp1
  bc_rhoR_north = dtemp1
  bc_rhoB_north = dtemp2
  bc_u_north = dtemp3
  bc_v_north = dtemp4
  bc_w_north = dtemp5
  
  return
  
 end subroutine set_value_bc_north
 
 subroutine set_value_bc_south(itemp1,dtemp1,dtemp2,dtemp3,dtemp4,dtemp5)
 
!***********************************************************************
!     
!     LBsoft subroutine for set the conditions of fluids 
!     at the south open boundary
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: itemp1
  real(kind=PRC), intent(in) :: dtemp1,dtemp2,dtemp3,dtemp4,dtemp5
  
  bc_type_south = itemp1
  bc_rhoR_south = dtemp1
  bc_rhoB_south = dtemp2
  bc_u_south = dtemp3
  bc_v_south = dtemp4
  bc_w_south = dtemp5
  
  return
  
 end subroutine set_value_bc_south
 
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
  
  lforce_add=.true.
  ext_fu = dtemp1
  ext_fv = dtemp2
  ext_fw = dtemp3
  
  return
  
 end subroutine set_value_ext_force_fluids
 
 subroutine set_LBintegrator_type(itemp1)
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the LB integrator type for
!     the collisional step
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: itemp1
  
  LBintegrator = itemp1
  
  return
  
 end subroutine set_LBintegrator_type
 
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
 
 subroutine set_lsingle_fluid(ltemp1)
 
!***********************************************************************
!     
!     LBsoft subroutine for set the value of lsingle_fluid for select 
!     the single fluid mode
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  logical, intent(in) :: ltemp1
  
  lsingle_fluid=ltemp1
  
  return
  
 end subroutine set_lsingle_fluid
 
 subroutine set_fluid_force_sc(ltemp1,dtemp1)
 
!***********************************************************************
!     
!     LBsoft subroutine for set the value of the pair ShanChen 
!     force constant
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  logical, intent(in) :: ltemp1
  real(kind=PRC), intent(in) :: dtemp1
  
  lforce_add=.true.
  lpair_SC = ltemp1
  lShanChen = (lShanChen .or. lpair_SC)
  pair_SC = dtemp1
  
  return
  
 end subroutine set_fluid_force_sc 
 
 subroutine set_fluid_wall_sc(dtemp1,dtemp2)
 
!***********************************************************************
!     
!     LBsoft subroutine for set the value of the pair ShanChen 
!     wall coupling constant
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC), intent(in) :: dtemp1,dtemp2
  
  wallR_SC = dtemp1
  wallB_SC = dtemp2
  lwall_SC = (wallR_SC/=ONE .or. wallB_SC/=ONE )
  
  return
  
 end subroutine set_fluid_wall_sc
 
 subroutine set_value_viscosity(dtemp1,dtemp2,temp_lsingle_fluid)
 
!***********************************************************************
!     
!     LBsoft subroutine for set the viscosity value and the relaxation 
!     time tau
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC), intent(in) :: dtemp1,dtemp2
  logical, intent(in) :: temp_lsingle_fluid
  
  if(dtemp1==dtemp2 .or. temp_lsingle_fluid)then
    lunique_omega=.true.
    viscR=dtemp1
    viscB=dtemp2
    tauR=viscR/cssq+HALF
    tauB=viscB/cssq+HALF
    unique_omega = viscosity_to_omega(viscR)
  else
    lunique_omega=.false.
    viscR=dtemp1
    viscB=dtemp2
    tauR=viscR/cssq+HALF
    tauB=viscB/cssq+HALF
  endif
  
  return
  
 end subroutine set_value_viscosity
 
 subroutine set_value_tau(dtemp1,dtemp2,temp_lsingle_fluid)
 
!***********************************************************************
!     
!     LBsoft subroutine for set the relaxation time tau value and 
!     the viscosity
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC), intent(in) :: dtemp1,dtemp2
  logical, intent(in) :: temp_lsingle_fluid
  
  if(dtemp1==dtemp2 .or. temp_lsingle_fluid)then
    lunique_omega=.true.
    tauR=dtemp1
    tauB=dtemp2
    viscR=cssq*(tauR-HALF)
    viscB=cssq*(tauB-HALF)
    unique_omega = viscosity_to_omega(viscR)
  else
    lunique_omega=.false.
    tauR=dtemp1
    tauB=dtemp2
    viscR=cssq*(tauR-HALF)
    viscB=cssq*(tauB-HALF)
  endif
  
  return
  
 end subroutine set_value_tau
 
 subroutine set_unique_omega
 
!***********************************************************************
!     
!     LBsoft subroutine for set a unique omega value
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  integer :: i,j,k
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    omega(i,j,k)=unique_omega
  end forall
  
  return
  
 end subroutine set_unique_omega
 
 pure function viscosity_to_omega(dtemp1)
 
!***********************************************************************
!     
!     LBsoft function to convert the viscosity to omega
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: dtemp1
  
  real(kind=PRC), parameter :: myfactor = TWO / cssq
  
  real(kind=PRC) :: viscosity_to_omega
  
  viscosity_to_omega = TWO / ( myfactor * dtemp1 + ONE )
  
  return
 
 end function viscosity_to_omega
 
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
    fwB(i,j,k) = fwB(i,j,k) - pair_SC*rhoB(i,j,k)*gradpsizR(i,j,k)
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
  
  if(lsingle_fluid)return
  
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
  !tolti gli zeri su dex dey dez con pazienza
  forall(i=1:nx,j=1:ny,k=1:nz)
    mygradx(i,j,k)= &
     myarr(i+ex(1),j+ey(1),k+ez(1))*p(1)*dex(1)+ & !01
     myarr(i+ex(2),j+ey(2),k+ez(2))*p(2)*dex(2)+ & !02
     myarr(i+ex(7),j+ey(7),k+ez(7))*p(7)*dex(7)+ & !07
     myarr(i+ex(8),j+ey(8),k+ez(8))*p(8)*dex(8)+ & !08
     myarr(i+ex(9),j+ey(9),k+ez(9))*p(9)*dex(9)+ & !09
     myarr(i+ex(10),j+ey(10),k+ez(10))*p(10)*dex(10)+ & !10
     myarr(i+ex(11),j+ey(11),k+ez(11))*p(11)*dex(11)+ & !11
     myarr(i+ex(12),j+ey(12),k+ez(12))*p(12)*dex(12)+ & !12
     myarr(i+ex(13),j+ey(13),k+ez(13))*p(13)*dex(13)+ & !13
     myarr(i+ex(14),j+ey(14),k+ez(14))*p(14)*dex(14)    !14
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    mygrady(i,j,k)= &
     myarr(i+ex(3),j+ey(3),k+ez(3))*p(3)*dey(3)+ & !03
     myarr(i+ex(4),j+ey(4),k+ez(4))*p(4)*dey(4)+ & !04
     myarr(i+ex(7),j+ey(7),k+ez(7))*p(7)*dey(7)+ & !07
     myarr(i+ex(8),j+ey(8),k+ez(8))*p(8)*dey(8)+ & !08
     myarr(i+ex(9),j+ey(9),k+ez(9))*p(9)*dey(9)+ & !09
     myarr(i+ex(10),j+ey(10),k+ez(10))*p(10)*dey(10)+ & !10
     myarr(i+ex(15),j+ey(15),k+ez(15))*p(15)*dey(15)+ & !15
     myarr(i+ex(16),j+ey(16),k+ez(16))*p(16)*dey(16)+ & !16
     myarr(i+ex(17),j+ey(17),k+ez(17))*p(17)*dey(17)+ & !17
     myarr(i+ex(18),j+ey(18),k+ez(18))*p(18)*dey(18)    !18
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    mygradz(i,j,k)= &
     myarr(i+ex(5),j+ey(5),k+ez(5))*p(5)*dez(5)+ & !05
     myarr(i+ex(6),j+ey(6),k+ez(6))*p(6)*dez(6)+ & !06
     myarr(i+ex(11),j+ey(11),k+ez(11))*p(11)*dez(11)+ & !11
     myarr(i+ex(12),j+ey(12),k+ez(12))*p(12)*dez(12)+ & !12
     myarr(i+ex(13),j+ey(13),k+ez(13))*p(13)*dez(13)+ & !13
     myarr(i+ex(14),j+ey(14),k+ez(14))*p(14)*dez(14)+ & !14
     myarr(i+ex(15),j+ey(15),k+ez(15))*p(15)*dez(15)+ & !15
     myarr(i+ex(16),j+ey(16),k+ez(16))*p(16)*dez(16)+ & !16
     myarr(i+ex(17),j+ey(17),k+ez(17))*p(17)*dez(17)+ & !17
     myarr(i+ex(18),j+ey(18),k+ez(18))*p(18)*dez(18)    !18
  end forall
#else
  !occhio gli zeri su dex dey dez andrebbero tolti con pazienza
  forall(i=1:nx,j=1:ny,k=1:nz)
    mygradx(i,j,k)= &
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
 
 subroutine driver_collision_fluids
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the collisional step
!     on the Boltzmann populations 
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  select case(LBintegrator)
  case (0)
    call collision_fluids_BGK(rhoR,u,v,w,omega,f00R,f01R,&
     f02R,f03R,f04R,f05R,f06R,f07R,f08R,f09R,f10R, &
     f11R,f12R,f13R,f14R,f15R,f16R,f17R,f18R)
    if(lsingle_fluid)return
    call collision_fluids_BGK(rhoB,u,v,w,omega,f00B,f01B,&
     f02B,f03B,f04B,f05B,f06B,f07B,f08B,f09B,f10B, &
     f11B,f12B,f13B,f14B,f15B,f16B,f17B,f18B)
  case (1)
    call convert_fluid_force_to_velshifted
    call collision_fluids_EDM(rhoR,u,v,w,fuR,fvR,fwR,omega,f00R,f01R,&
     f02R,f03R,f04R,f05R,f06R,f07R,f08R,f09R,f10R, &
     f11R,f12R,f13R,f14R,f15R,f16R,f17R,f18R)
    if(lsingle_fluid)return
    call collision_fluids_EDM(rhoB,u,v,w,fuB,fvB,fwB,omega,f00B,f01B,&
     f02B,f03B,f04B,f05B,f06B,f07B,f08B,f09B,f10B, &
     f11B,f12B,f13B,f14B,f15B,f16B,f17B,f18B)
  case default
    call error(14)
  end select
  
  return
  
 end subroutine driver_collision_fluids
 
 subroutine collision_fluids_BGK(rhosub,usub,vsub,wsub,omegas,f00sub,f01sub,&
  f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub,f09sub,f10sub, &
  f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub)
 
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
  
  real(kind=PRC), allocatable, dimension(:,:,:)  :: rhosub,usub,vsub, &
   wsub,omegas,f00sub,f01sub,f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub,&
   f09sub,f10sub,f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub
  
  integer :: i,j,k
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f00sub(i,j,k)=f00sub(i,j,k)+omegas(i,j,k)* &
     (equil_pop00(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k),wsub(i,j,k))- &
     f00sub(i,j,k))
  end forall
    
  forall(i=1:nx,j=1:ny,k=1:nz)
    f01sub(i,j,k)=f01sub(i,j,k)+omegas(i,j,k)* &
     (equil_pop01(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k),wsub(i,j,k))- &
     f01sub(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f02sub(i,j,k)=f02sub(i,j,k)+omegas(i,j,k)* &
     (equil_pop02(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k),wsub(i,j,k))- &
     f02sub(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f03sub(i,j,k)=f03sub(i,j,k)+omegas(i,j,k)* &
     (equil_pop03(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k),wsub(i,j,k))- &
     f03sub(i,j,k))
   end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f04sub(i,j,k)=f04sub(i,j,k)+omegas(i,j,k)* &
     (equil_pop04(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k),wsub(i,j,k))- &
     f04sub(i,j,k))
  end forall
    
  forall(i=1:nx,j=1:ny,k=1:nz)
    f05sub(i,j,k)=f05sub(i,j,k)+omegas(i,j,k)* &
     (equil_pop05(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k),wsub(i,j,k))- &
     f05sub(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f06sub(i,j,k)=f06sub(i,j,k)+omegas(i,j,k)* &
     (equil_pop06(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k),wsub(i,j,k))- &
     f06sub(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f07sub(i,j,k)=f07sub(i,j,k)+omegas(i,j,k)* &
     (equil_pop07(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k),wsub(i,j,k))- &
     f07sub(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f08sub(i,j,k)=f08sub(i,j,k)+omegas(i,j,k)* &
     (equil_pop08(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k),wsub(i,j,k))- &
     f08sub(i,j,k))
  end forall
    
  forall(i=1:nx,j=1:ny,k=1:nz)
    f09sub(i,j,k)=f09sub(i,j,k)+omegas(i,j,k)* &
     (equil_pop09(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k),wsub(i,j,k))- &
     f09sub(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f10sub(i,j,k)=f10sub(i,j,k)+omegas(i,j,k)* &
     (equil_pop10(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k),wsub(i,j,k))- &
     f10sub(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f11sub(i,j,k)=f11sub(i,j,k)+omegas(i,j,k)* &
     (equil_pop11(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k),wsub(i,j,k))- &
     f11sub(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f12sub(i,j,k)=f12sub(i,j,k)+omegas(i,j,k)* &
     (equil_pop12(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k),wsub(i,j,k))- &
     f12sub(i,j,k))
  end forall
    
  forall(i=1:nx,j=1:ny,k=1:nz)
    f13sub(i,j,k)=f13sub(i,j,k)+omegas(i,j,k)* &
     (equil_pop13(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k),wsub(i,j,k))- &
     f13sub(i,j,k))
  end forall
    
  forall(i=1:nx,j=1:ny,k=1:nz)
    f14sub(i,j,k)=f14sub(i,j,k)+omegas(i,j,k)* &
     (equil_pop14(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k),wsub(i,j,k))- &
     f14sub(i,j,k))
  end forall
    
  forall(i=1:nx,j=1:ny,k=1:nz)
    f15sub(i,j,k)=f15sub(i,j,k)+omegas(i,j,k)* &
     (equil_pop15(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k),wsub(i,j,k))- &
     f15sub(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f16sub(i,j,k)=f16sub(i,j,k)+omegas(i,j,k)* &
     (equil_pop16(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k),wsub(i,j,k))- &
     f16sub(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f17sub(i,j,k)=f17sub(i,j,k)+omegas(i,j,k)* &
     (equil_pop17(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k),wsub(i,j,k))- &
     f17sub(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f18sub(i,j,k)=f18sub(i,j,k)+omegas(i,j,k)* &
     (equil_pop18(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k),wsub(i,j,k))- &
     f18sub(i,j,k))
  end forall
  
  return
  
 end subroutine collision_fluids_BGK
 
 subroutine collision_fluids_EDM(rhosub,usub,vsub,wsub,fusub,fvsub, &
  fwsub,omegas,f00sub,f01sub,f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,&
  f08sub,f09sub,f10sub,f11sub,f12sub,f13sub,f14sub,f15sub,f16sub, &
  f17sub,f18sub)
 
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
  
  real(kind=PRC), allocatable, dimension(:,:,:)  :: rhosub,usub,vsub, &
   wsub,fusub,fvsub,fwsub,omegas,f00sub,f01sub,f02sub,f03sub,f04sub, &
   f05sub,f06sub,f07sub,f08sub,f09sub,f10sub,f11sub,f12sub,f13sub, &
   f14sub,f15sub,f16sub,f17sub,f18sub
  
  integer :: i,j,k
    
  forall(i=1:nx,j=1:ny,k=1:nz)
    f00sub(i,j,k)=(ONE-omegas(i,j,k))*f00sub(i,j,k)+(omegas(i,j,k)-ONE)* &
     equil_pop00(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k),wsub(i,j,k))+ &
     equil_pop00(rhosub(i,j,k),fusub(i,j,k),fvsub(i,j,k),fwsub(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f01sub(i,j,k)=(ONE-omegas(i,j,k))*f01sub(i,j,k)+(omegas(i,j,k)-ONE)* &
     equil_pop01(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k),wsub(i,j,k))+ &
     equil_pop01(rhosub(i,j,k),fusub(i,j,k),fvsub(i,j,k),fwsub(i,j,k))
  end forall
    
  forall(i=1:nx,j=1:ny,k=1:nz)
    f02sub(i,j,k)=(ONE-omegas(i,j,k))*f02sub(i,j,k)+(omegas(i,j,k)-ONE)* &
     equil_pop02(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k),wsub(i,j,k))+ &
     equil_pop02(rhosub(i,j,k),fusub(i,j,k),fvsub(i,j,k),fwsub(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f03sub(i,j,k)=(ONE-omegas(i,j,k))*f03sub(i,j,k)+(omegas(i,j,k)-ONE)* &
     equil_pop03(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k),wsub(i,j,k))+ &
     equil_pop03(rhosub(i,j,k),fusub(i,j,k),fvsub(i,j,k),fwsub(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f04sub(i,j,k)=(ONE-omegas(i,j,k))*f04sub(i,j,k)+(omegas(i,j,k)-ONE)* &
     equil_pop04(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k),wsub(i,j,k))+ &
     equil_pop04(rhosub(i,j,k),fusub(i,j,k),fvsub(i,j,k),fwsub(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f05sub(i,j,k)=(ONE-omegas(i,j,k))*f05sub(i,j,k)+(omegas(i,j,k)-ONE)* &
     equil_pop05(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k),wsub(i,j,k))+ &
     equil_pop05(rhosub(i,j,k),fusub(i,j,k),fvsub(i,j,k),fwsub(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f06sub(i,j,k)=(ONE-omegas(i,j,k))*f06sub(i,j,k)+(omegas(i,j,k)-ONE)* &
     equil_pop06(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k),wsub(i,j,k))+ &
     equil_pop06(rhosub(i,j,k),fusub(i,j,k),fvsub(i,j,k),fwsub(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f07sub(i,j,k)=(ONE-omegas(i,j,k))*f07sub(i,j,k)+(omegas(i,j,k)-ONE)* &
     equil_pop07(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k),wsub(i,j,k))+ &
     equil_pop07(rhosub(i,j,k),fusub(i,j,k),fvsub(i,j,k),fwsub(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f08sub(i,j,k)=(ONE-omegas(i,j,k))*f08sub(i,j,k)+(omegas(i,j,k)-ONE)* &
     equil_pop08(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k),wsub(i,j,k))+ &
     equil_pop08(rhosub(i,j,k),fusub(i,j,k),fvsub(i,j,k),fwsub(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f09sub(i,j,k)=(ONE-omegas(i,j,k))*f09sub(i,j,k)+(omegas(i,j,k)-ONE)* &
     equil_pop09(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k),wsub(i,j,k))+ &
     equil_pop09(rhosub(i,j,k),fusub(i,j,k),fvsub(i,j,k),fwsub(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f10sub(i,j,k)=(ONE-omegas(i,j,k))*f10sub(i,j,k)+(omegas(i,j,k)-ONE)* &
     equil_pop10(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k),wsub(i,j,k))+ &
     equil_pop10(rhosub(i,j,k),fusub(i,j,k),fvsub(i,j,k),fwsub(i,j,k))
  end forall
    
  forall(i=1:nx,j=1:ny,k=1:nz)
    f11sub(i,j,k)=(ONE-omegas(i,j,k))*f11sub(i,j,k)+(omegas(i,j,k)-ONE)* &
     equil_pop11(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k),wsub(i,j,k))+ &
     equil_pop11(rhosub(i,j,k),fusub(i,j,k),fvsub(i,j,k),fwsub(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f12sub(i,j,k)=(ONE-omegas(i,j,k))*f12sub(i,j,k)+(omegas(i,j,k)-ONE)* &
     equil_pop12(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k),wsub(i,j,k))+ &
     equil_pop12(rhosub(i,j,k),fusub(i,j,k),fvsub(i,j,k),fwsub(i,j,k))
  end forall
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    f13sub(i,j,k)=(ONE-omegas(i,j,k))*f13sub(i,j,k)+(omegas(i,j,k)-ONE)* &
     equil_pop13(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k),wsub(i,j,k))+ &
     equil_pop13(rhosub(i,j,k),fusub(i,j,k),fvsub(i,j,k),fwsub(i,j,k))
  end forall
    
  forall(i=1:nx,j=1:ny,k=1:nz)
    f14sub(i,j,k)=(ONE-omegas(i,j,k))*f14sub(i,j,k)+(omegas(i,j,k)-ONE)* &
     equil_pop14(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k),wsub(i,j,k))+ &
     equil_pop14(rhosub(i,j,k),fusub(i,j,k),fvsub(i,j,k),fwsub(i,j,k))
  end forall
    
  forall(i=1:nx,j=1:ny,k=1:nz)
    f15sub(i,j,k)=(ONE-omegas(i,j,k))*f15sub(i,j,k)+(omegas(i,j,k)-ONE)* &
     equil_pop15(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k),wsub(i,j,k))+ &
     equil_pop15(rhosub(i,j,k),fusub(i,j,k),fvsub(i,j,k),fwsub(i,j,k))
  end forall
    
  forall(i=1:nx,j=1:ny,k=1:nz)
    f16sub(i,j,k)=(ONE-omegas(i,j,k))*f16sub(i,j,k)+(omegas(i,j,k)-ONE)* &
     equil_pop16(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k),wsub(i,j,k))+ &
     equil_pop16(rhosub(i,j,k),fusub(i,j,k),fvsub(i,j,k),fwsub(i,j,k))
  end forall
    
  forall(i=1:nx,j=1:ny,k=1:nz)
    f17sub(i,j,k)=(ONE-omegas(i,j,k))*f17sub(i,j,k)+(omegas(i,j,k)-ONE)* &
     equil_pop17(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k),wsub(i,j,k))+ &
     equil_pop17(rhosub(i,j,k),fusub(i,j,k),fvsub(i,j,k),fwsub(i,j,k))
  end forall
    
  forall(i=1:nx,j=1:ny,k=1:nz)
    f18sub(i,j,k)=(ONE-omegas(i,j,k))*f18sub(i,j,k)+(omegas(i,j,k)-ONE)* &
     equil_pop18(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k),wsub(i,j,k))+ &
     equil_pop18(rhosub(i,j,k),fusub(i,j,k),fvsub(i,j,k),fwsub(i,j,k))
  end forall
  
  return
  
 end subroutine collision_fluids_EDM
 
 subroutine compute_omega
 
!***********************************************************************
!     
!     LBsoft subroutine for computing the relative omega of the
!     fluid system
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  integer :: i,j,k
  
  if(lunique_omega)return
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    omega(i,j,k)=viscosity_to_omega( &
     ONE/((rhoR(i,j,k)/(rhoB(i,j,k)+rhoR(i,j,k)))*(ONE/viscR) + &
     (rhoB(i,j,k)/(rhoB(i,j,k)+rhoR(i,j,k)))*(ONE/viscB)) )
  end forall
  
  return
  
 end subroutine compute_omega
 
 subroutine driver_streaming_fluids
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the streaming step
!     on the Boltzmann populations
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  
  call streaming_fluids(f00R,f01R,f02R,f03R,f04R, &
   f05R,f06R,f07R,f08R,f09R,f10R,f11R,f12R,f13R, &
   f14R,f15R,f16R,f17R,f18R)
   
  if(lsingle_fluid)return
  
  call streaming_fluids(f00B,f01B,f02B,f03B,f04B, &
   f05B,f06B,f07B,f08B,f09B,f10B,f11B,f12B,f13B, &
   f14B,f15B,f16B,f17B,f18B)
  
  
  return
  
 end subroutine driver_streaming_fluids
 
 subroutine streaming_fluids(f00sub,f01sub,f02sub,f03sub,f04sub, &
   f05sub,f06sub,f07sub,f08sub,f09sub,f10sub,f11sub,f12sub,f13sub, &
   f14sub,f15sub,f16sub,f17sub,f18sub)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the streaming step
!     on the Boltzmann populations
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC), allocatable, dimension(:,:,:)  :: f00sub,f01sub, &
   f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub,f09sub,f10sub, &
   f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub
  
  integer :: i,j,k,l,ishift,jshift,kshift
  

  l=1
  ishift=ex(l)
  jshift=ey(l)
  kshift=ez(l)
  forall(i=0:nx+1,j=0:ny+1,k=0:nz+1)
    buffservice3d(i+ishift,j+jshift,k+kshift) = f01sub(i,j,k)
  end forall
  forall(i=0:nx+1,j=0:ny+1,k=0:nz+1)f01sub(i,j,k) = buffservice3d(i,j,k)
  
  l=2
  ishift=ex(l)
  jshift=ey(l)
  kshift=ez(l)
  forall(i=0:nx+1,j=0:ny+1,k=0:nz+1)
    buffservice3d(i+ishift,j+jshift,k+kshift) = f02sub(i,j,k)
  end forall
  forall(i=0:nx+1,j=0:ny+1,k=0:nz+1)f02sub(i,j,k) = buffservice3d(i,j,k)
  
  l=3
  ishift=ex(l)
  jshift=ey(l)
  kshift=ez(l)
  forall(i=0:nx+1,j=0:ny+1,k=0:nz+1)
    buffservice3d(i+ishift,j+jshift,k+kshift) = f03sub(i,j,k)
  end forall
  forall(i=0:nx+1,j=0:ny+1,k=0:nz+1)f03sub(i,j,k) = buffservice3d(i,j,k)
  
  l=4
  ishift=ex(l)
  jshift=ey(l)
  kshift=ez(l)
  forall(i=0:nx+1,j=0:ny+1,k=0:nz+1)
    buffservice3d(i+ishift,j+jshift,k+kshift) = f04sub(i,j,k)
  end forall
  forall(i=0:nx+1,j=0:ny+1,k=0:nz+1)f04sub(i,j,k) = buffservice3d(i,j,k)
  
  l=5
  ishift=ex(l)
  jshift=ey(l)
  kshift=ez(l)
  forall(i=0:nx+1,j=0:ny+1,k=0:nz+1)
    buffservice3d(i+ishift,j+jshift,k+kshift) = f05sub(i,j,k)
  end forall
  forall(i=0:nx+1,j=0:ny+1,k=0:nz+1)f05sub(i,j,k) = buffservice3d(i,j,k)
  
  l=6
  ishift=ex(l)
  jshift=ey(l)
  kshift=ez(l)
  forall(i=0:nx+1,j=0:ny+1,k=0:nz+1)
    buffservice3d(i+ishift,j+jshift,k+kshift) = f06sub(i,j,k)
  end forall
  forall(i=0:nx+1,j=0:ny+1,k=0:nz+1)f06sub(i,j,k) = buffservice3d(i,j,k)
  
  l=7
  ishift=ex(l)
  jshift=ey(l)
  kshift=ez(l)
  forall(i=0:nx+1,j=0:ny+1,k=0:nz+1)
    buffservice3d(i+ishift,j+jshift,k+kshift) = f07sub(i,j,k)
  end forall
  forall(i=0:nx+1,j=0:ny+1,k=0:nz+1)f07sub(i,j,k) = buffservice3d(i,j,k)
  
  l=8
  ishift=ex(l)
  jshift=ey(l)
  kshift=ez(l)
  forall(i=0:nx+1,j=0:ny+1,k=0:nz+1)
    buffservice3d(i+ishift,j+jshift,k+kshift) = f08sub(i,j,k)
  end forall
  forall(i=0:nx+1,j=0:ny+1,k=0:nz+1)f08sub(i,j,k) = buffservice3d(i,j,k)
  
  l=9
  ishift=ex(l)
  jshift=ey(l)
  kshift=ez(l)
  forall(i=0:nx+1,j=0:ny+1,k=0:nz+1)
    buffservice3d(i+ishift,j+jshift,k+kshift) = f09sub(i,j,k)
  end forall
  forall(i=0:nx+1,j=0:ny+1,k=0:nz+1)f09sub(i,j,k) = buffservice3d(i,j,k)
  
  l=10
  ishift=ex(l)
  jshift=ey(l)
  kshift=ez(l)
  forall(i=0:nx+1,j=0:ny+1,k=0:nz+1)
    buffservice3d(i+ishift,j+jshift,k+kshift) = f10sub(i,j,k)
  end forall
  forall(i=0:nx+1,j=0:ny+1,k=0:nz+1)f10sub(i,j,k) = buffservice3d(i,j,k)
  
  l=11
  ishift=ex(l)
  jshift=ey(l)
  kshift=ez(l)
  forall(i=0:nx+1,j=0:ny+1,k=0:nz+1)
    buffservice3d(i+ishift,j+jshift,k+kshift) = f11sub(i,j,k)
  end forall
  forall(i=0:nx+1,j=0:ny+1,k=0:nz+1)f11sub(i,j,k) = buffservice3d(i,j,k)
  
  l=12
  ishift=ex(l)
  jshift=ey(l)
  kshift=ez(l)
  forall(i=0:nx+1,j=0:ny+1,k=0:nz+1)
    buffservice3d(i+ishift,j+jshift,k+kshift) = f12sub(i,j,k)
  end forall
  forall(i=0:nx+1,j=0:ny+1,k=0:nz+1)f12sub(i,j,k) = buffservice3d(i,j,k)
  
  l=13
  ishift=ex(l)
  jshift=ey(l)
  kshift=ez(l)
  forall(i=0:nx+1,j=0:ny+1,k=0:nz+1)
    buffservice3d(i+ishift,j+jshift,k+kshift) = f13sub(i,j,k)
  end forall
  forall(i=0:nx+1,j=0:ny+1,k=0:nz+1)f13sub(i,j,k) = buffservice3d(i,j,k)
  
  l=14
  ishift=ex(l)
  jshift=ey(l)
  kshift=ez(l)
  forall(i=0:nx+1,j=0:ny+1,k=0:nz+1)
    buffservice3d(i+ishift,j+jshift,k+kshift) = f14sub(i,j,k)
  end forall
  forall(i=0:nx+1,j=0:ny+1,k=0:nz+1)f14sub(i,j,k) = buffservice3d(i,j,k)
  
  l=15
  ishift=ex(l)
  jshift=ey(l)
  kshift=ez(l)
  forall(i=0:nx+1,j=0:ny+1,k=0:nz+1)
    buffservice3d(i+ishift,j+jshift,k+kshift) = f15sub(i,j,k)
  end forall
  forall(i=0:nx+1,j=0:ny+1,k=0:nz+1)f15sub(i,j,k) = buffservice3d(i,j,k)
  
  l=16
  ishift=ex(l)
  jshift=ey(l)
  kshift=ez(l)
  forall(i=0:nx+1,j=0:ny+1,k=0:nz+1)
    buffservice3d(i+ishift,j+jshift,k+kshift) = f16sub(i,j,k)
  end forall
  forall(i=0:nx+1,j=0:ny+1,k=0:nz+1)f16sub(i,j,k) = buffservice3d(i,j,k)
  
  l=17
  ishift=ex(l)
  jshift=ey(l)
  kshift=ez(l)
  forall(i=0:nx+1,j=0:ny+1,k=0:nz+1)
    buffservice3d(i+ishift,j+jshift,k+kshift) = f17sub(i,j,k)
  end forall
  forall(i=0:nx+1,j=0:ny+1,k=0:nz+1)f17sub(i,j,k) = buffservice3d(i,j,k)
  
  l=18
  ishift=ex(l)
  jshift=ey(l)
  kshift=ez(l)
  forall(i=0:nx+1,j=0:ny+1,k=0:nz+1)
    buffservice3d(i+ishift,j+jshift,k+kshift) = f18sub(i,j,k)
  end forall
  forall(i=0:nx+1,j=0:ny+1,k=0:nz+1)f18sub(i,j,k) = buffservice3d(i,j,k)
  
  return
  
 end subroutine streaming_fluids
 
 subroutine moments_fluids
 
!***********************************************************************
!     
!     LBsoft subroutine for computing the moments on the 
!     Boltzmann populations and estimate the hydrodynamic variables
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  integer :: i,j,k,l
  
  real(kind=PRC) :: ddx,ddy,ddz
  
  !compute density and accumulate mass flux
  
  forall(i=1:nx,j=1:ny,k=1:nz)
    u(i,j,k)    = ZERO
    v(i,j,k)    = ZERO
    w(i,j,k)    = ZERO
  end forall
  
  !red fluid
  
  if(lsingle_fluid)then
    
    l=0
    ddx=dex(l)
    ddy=dey(l)
    ddz=dez(l)
    forall(i=1:nx,j=1:ny,k=1:nz)
      rhoR(i,j,k) = f00R(i,j,k)
    end forall
  
    l=1
    ddx=dex(l)
    ddy=dey(l)
    ddz=dez(l)
    forall(i=1:nx,j=1:ny,k=1:nz)
      rhoR(i,j,k) = f01R(i,j,k) + rhoR(i,j,k)
      u(i,j,k)    = f01R(i,j,k)*ddx + u(i,j,k)
    end forall
    
    l=2
    ddx=dex(l)
    ddy=dey(l)
    ddz=dez(l)
    forall(i=1:nx,j=1:ny,k=1:nz)
      rhoR(i,j,k) = f02R(i,j,k) + rhoR(i,j,k)
      u(i,j,k)    = f02R(i,j,k)*ddx + u(i,j,k)
    end forall
  
    l=3
    ddx=dex(l)
    ddy=dey(l)
    ddz=dez(l)
    forall(i=1:nx,j=1:ny,k=1:nz)
      rhoR(i,j,k) = f03R(i,j,k) + rhoR(i,j,k)
      v(i,j,k)    = f03R(i,j,k)*ddy + v(i,j,k)
    end forall
    
    l=4
    ddx=dex(l)
    ddy=dey(l)
    ddz=dez(l)
    forall(i=1:nx,j=1:ny,k=1:nz)
      rhoR(i,j,k) = f04R(i,j,k) + rhoR(i,j,k)
      v(i,j,k)    = f04R(i,j,k)*ddy + v(i,j,k)
    end forall
  
    l=5
    ddx=dex(l)
    ddy=dey(l)
    ddz=dez(l)
    forall(i=1:nx,j=1:ny,k=1:nz)
      rhoR(i,j,k) = f05R(i,j,k) + rhoR(i,j,k)
      w(i,j,k)    = f05R(i,j,k)*ddz + w(i,j,k)
    end forall
    
    l=6
    ddx=dex(l)
    ddy=dey(l)
    ddz=dez(l)
    forall(i=1:nx,j=1:ny,k=1:nz)
      rhoR(i,j,k) = f06R(i,j,k) + rhoR(i,j,k)
      w(i,j,k)    = f06R(i,j,k)*ddz + w(i,j,k)
    end forall
  
    l=7
    ddx=dex(l)
    ddy=dey(l)
    ddz=dez(l)
    forall(i=1:nx,j=1:ny,k=1:nz)
      rhoR(i,j,k) = f07R(i,j,k) + rhoR(i,j,k)
      u(i,j,k)    = f07R(i,j,k)*ddx + u(i,j,k)
      v(i,j,k)    = f07R(i,j,k)*ddy + v(i,j,k)
    end forall
    
    l=8
    ddx=dex(l)
    ddy=dey(l)
    ddz=dez(l)
    forall(i=1:nx,j=1:ny,k=1:nz)
      rhoR(i,j,k) = f08R(i,j,k) + rhoR(i,j,k)
      u(i,j,k)    = f08R(i,j,k)*ddx + u(i,j,k)
      v(i,j,k)    = f08R(i,j,k)*ddy + v(i,j,k)
    end forall
    
    l=9
    ddx=dex(l)
    ddy=dey(l)
    ddz=dez(l)
    forall(i=1:nx,j=1:ny,k=1:nz)
      rhoR(i,j,k) = f09R(i,j,k) + rhoR(i,j,k)
      u(i,j,k)    = f09R(i,j,k)*ddx + u(i,j,k)
      v(i,j,k)    = f09R(i,j,k)*ddy + v(i,j,k)
    end forall
    
    l=10
    ddx=dex(l)
    ddy=dey(l)
    ddz=dez(l)
    forall(i=1:nx,j=1:ny,k=1:nz)
      rhoR(i,j,k) = f10R(i,j,k) + rhoR(i,j,k)
      u(i,j,k)    = f10R(i,j,k)*ddx + u(i,j,k)
      v(i,j,k)    = f10R(i,j,k)*ddy + v(i,j,k)
    end forall
    
    l=11
    ddx=dex(l)
    ddy=dey(l)
    ddz=dez(l)
    forall(i=1:nx,j=1:ny,k=1:nz)
      rhoR(i,j,k) = f11R(i,j,k) + rhoR(i,j,k)
      u(i,j,k)    = f11R(i,j,k)*ddx + u(i,j,k)
      w(i,j,k)    = f11R(i,j,k)*ddz + w(i,j,k)
    end forall
    
    l=12
    ddx=dex(l)
    ddy=dey(l)
    ddz=dez(l)
    forall(i=1:nx,j=1:ny,k=1:nz)
      rhoR(i,j,k) = f12R(i,j,k) + rhoR(i,j,k)
      u(i,j,k)    = f12R(i,j,k)*ddx + u(i,j,k)
      w(i,j,k)    = f12R(i,j,k)*ddz + w(i,j,k)
    end forall
    
    l=13
    ddx=dex(l)
    ddy=dey(l)
    ddz=dez(l)
    forall(i=1:nx,j=1:ny,k=1:nz)
      rhoR(i,j,k) = f13R(i,j,k) + rhoR(i,j,k)
      u(i,j,k)    = f13R(i,j,k)*ddx + u(i,j,k)
      w(i,j,k)    = f13R(i,j,k)*ddz + w(i,j,k)
    end forall
    
    l=14
    ddx=dex(l)
    ddy=dey(l)
    ddz=dez(l)
    forall(i=1:nx,j=1:ny,k=1:nz)
      rhoR(i,j,k) = f14R(i,j,k) + rhoR(i,j,k)
      u(i,j,k)    = f14R(i,j,k)*ddx + u(i,j,k)
      w(i,j,k)    = f14R(i,j,k)*ddz + w(i,j,k)
    end forall
    
    l=15
    ddx=dex(l)
    ddy=dey(l)
    ddz=dez(l)
    forall(i=1:nx,j=1:ny,k=1:nz)
      rhoR(i,j,k) = f15R(i,j,k) + rhoR(i,j,k)
      v(i,j,k)    = f15R(i,j,k)*ddy + v(i,j,k)
      w(i,j,k)    = f15R(i,j,k)*ddz + w(i,j,k)
    end forall
  
    l=16
    ddx=dex(l)
    ddy=dey(l)
    ddz=dez(l)
    forall(i=1:nx,j=1:ny,k=1:nz)
      rhoR(i,j,k) = f16R(i,j,k) + rhoR(i,j,k)
      v(i,j,k)    = f16R(i,j,k)*ddy + v(i,j,k)
      w(i,j,k)    = f16R(i,j,k)*ddz + w(i,j,k)
    end forall
  
    l=17
    ddx=dex(l)
    ddy=dey(l)
    ddz=dez(l)
    forall(i=1:nx,j=1:ny,k=1:nz)
      rhoR(i,j,k) = f17R(i,j,k) + rhoR(i,j,k)
      v(i,j,k)    = f17R(i,j,k)*ddy + v(i,j,k)
      w(i,j,k)    = f17R(i,j,k)*ddz + w(i,j,k)
    end forall
  
    l=18
    ddx=dex(l)
    ddy=dey(l)
    ddz=dez(l)
    forall(i=1:nx,j=1:ny,k=1:nz)
      rhoR(i,j,k) = f18R(i,j,k) + rhoR(i,j,k)
      v(i,j,k)    = f18R(i,j,k)*ddy + v(i,j,k)
      w(i,j,k)    = f18R(i,j,k)*ddz + w(i,j,k)
    end forall
    
    !compute speed from mass flux
    forall(i=1:nx,j=1:ny,k=1:nz)
      u(i,j,k) = u(i,j,k)/rhoR(i,j,k)
      v(i,j,k) = v(i,j,k)/rhoR(i,j,k)
      w(i,j,k) = w(i,j,k)/rhoR(i,j,k)
    end forall
    return
  endif
  
  l=0
  ddx=dex(l)/tauR
  ddy=dey(l)/tauR
  ddz=dez(l)/tauR
  forall(i=1:nx,j=1:ny,k=1:nz)
    rhoR(i,j,k) = f00R(i,j,k)
  end forall
  
  l=1
  ddx=dex(l)/tauR
  ddy=dey(l)/tauR
  ddz=dez(l)/tauR
  forall(i=1:nx,j=1:ny,k=1:nz)
    rhoR(i,j,k) = f01R(i,j,k) + rhoR(i,j,k)
    u(i,j,k)    = f01R(i,j,k)*ddx + u(i,j,k)
  end forall
  
  l=2
  ddx=dex(l)/tauR
  ddy=dey(l)/tauR
  ddz=dez(l)/tauR
  forall(i=1:nx,j=1:ny,k=1:nz)
    rhoR(i,j,k) = f02R(i,j,k) + rhoR(i,j,k)
    u(i,j,k)    = f02R(i,j,k)*ddx + u(i,j,k)
  end forall
  
  l=3
  ddx=dex(l)/tauR
  ddy=dey(l)/tauR
  ddz=dez(l)/tauR
  forall(i=1:nx,j=1:ny,k=1:nz)
    rhoR(i,j,k) = f03R(i,j,k) + rhoR(i,j,k)
    v(i,j,k)    = f03R(i,j,k)*ddy + v(i,j,k)
  end forall
  
  l=4
  ddx=dex(l)/tauR
  ddy=dey(l)/tauR
  ddz=dez(l)/tauR
  forall(i=1:nx,j=1:ny,k=1:nz)
    rhoR(i,j,k) = f04R(i,j,k) + rhoR(i,j,k)
    v(i,j,k)    = f04R(i,j,k)*ddy + v(i,j,k)
  end forall
  
  l=5
  ddx=dex(l)/tauR
  ddy=dey(l)/tauR
  ddz=dez(l)/tauR
  forall(i=1:nx,j=1:ny,k=1:nz)
    rhoR(i,j,k) = f05R(i,j,k) + rhoR(i,j,k)
    w(i,j,k)    = f05R(i,j,k)*ddz + w(i,j,k)
  end forall
  
  l=6
  ddx=dex(l)/tauR
  ddy=dey(l)/tauR
  ddz=dez(l)/tauR
  forall(i=1:nx,j=1:ny,k=1:nz)
    rhoR(i,j,k) = f06R(i,j,k) + rhoR(i,j,k)
    w(i,j,k)    = f06R(i,j,k)*ddz + w(i,j,k)
  end forall
  
  l=7
  ddx=dex(l)/tauR
  ddy=dey(l)/tauR
  ddz=dez(l)/tauR
  forall(i=1:nx,j=1:ny,k=1:nz)
    rhoR(i,j,k) = f07R(i,j,k) + rhoR(i,j,k)
    u(i,j,k)    = f07R(i,j,k)*ddx + u(i,j,k)
    v(i,j,k)    = f07R(i,j,k)*ddy + v(i,j,k)
  end forall
  
  l=8
  ddx=dex(l)/tauR
  ddy=dey(l)/tauR
  ddz=dez(l)/tauR
  forall(i=1:nx,j=1:ny,k=1:nz)
    rhoR(i,j,k) = f08R(i,j,k) + rhoR(i,j,k)
    u(i,j,k)    = f08R(i,j,k)*ddx + u(i,j,k)
    v(i,j,k)    = f08R(i,j,k)*ddy + v(i,j,k)
  end forall
  
  l=9
  ddx=dex(l)/tauR
  ddy=dey(l)/tauR
  ddz=dez(l)/tauR
  forall(i=1:nx,j=1:ny,k=1:nz)
    rhoR(i,j,k) = f09R(i,j,k) + rhoR(i,j,k)
    u(i,j,k)    = f09R(i,j,k)*ddx + u(i,j,k)
    v(i,j,k)    = f09R(i,j,k)*ddy + v(i,j,k)
  end forall
  
  l=10
  ddx=dex(l)/tauR
  ddy=dey(l)/tauR
  ddz=dez(l)/tauR
  forall(i=1:nx,j=1:ny,k=1:nz)
    rhoR(i,j,k) = f10R(i,j,k) + rhoR(i,j,k)
    u(i,j,k)    = f10R(i,j,k)*ddx + u(i,j,k)
    v(i,j,k)    = f10R(i,j,k)*ddy + v(i,j,k)
  end forall
  
  l=11
  ddx=dex(l)/tauR
  ddy=dey(l)/tauR
  ddz=dez(l)/tauR
  forall(i=1:nx,j=1:ny,k=1:nz)
    rhoR(i,j,k) = f11R(i,j,k) + rhoR(i,j,k)
    u(i,j,k)    = f11R(i,j,k)*ddx + u(i,j,k)
    w(i,j,k)    = f11R(i,j,k)*ddz + w(i,j,k)
  end forall
  
  l=12
  ddx=dex(l)/tauR
  ddy=dey(l)/tauR
  ddz=dez(l)/tauR
  forall(i=1:nx,j=1:ny,k=1:nz)
    rhoR(i,j,k) = f12R(i,j,k) + rhoR(i,j,k)
    u(i,j,k)    = f12R(i,j,k)*ddx + u(i,j,k)
    w(i,j,k)    = f12R(i,j,k)*ddz + w(i,j,k)
  end forall
  
  l=13
  ddx=dex(l)/tauR
  ddy=dey(l)/tauR
  ddz=dez(l)/tauR
  forall(i=1:nx,j=1:ny,k=1:nz)
    rhoR(i,j,k) = f13R(i,j,k) + rhoR(i,j,k)
    u(i,j,k)    = f13R(i,j,k)*ddx + u(i,j,k)
    w(i,j,k)    = f13R(i,j,k)*ddz + w(i,j,k)
  end forall
  
  l=14
  ddx=dex(l)/tauR
  ddy=dey(l)/tauR
  ddz=dez(l)/tauR
  forall(i=1:nx,j=1:ny,k=1:nz)
    rhoR(i,j,k) = f14R(i,j,k) + rhoR(i,j,k)
    u(i,j,k)    = f14R(i,j,k)*ddx + u(i,j,k)
    w(i,j,k)    = f14R(i,j,k)*ddz + w(i,j,k)
  end forall
  
  l=15
  ddx=dex(l)/tauR
  ddy=dey(l)/tauR
  ddz=dez(l)/tauR
  forall(i=1:nx,j=1:ny,k=1:nz)
    rhoR(i,j,k) = f15R(i,j,k) + rhoR(i,j,k)
    v(i,j,k)    = f15R(i,j,k)*ddy + v(i,j,k)
    w(i,j,k)    = f15R(i,j,k)*ddz + w(i,j,k)
  end forall
  
  l=16
  ddx=dex(l)/tauR
  ddy=dey(l)/tauR
  ddz=dez(l)/tauR
  forall(i=1:nx,j=1:ny,k=1:nz)
    rhoR(i,j,k) = f16R(i,j,k) + rhoR(i,j,k)
    v(i,j,k)    = f16R(i,j,k)*ddy + v(i,j,k)
    w(i,j,k)    = f16R(i,j,k)*ddz + w(i,j,k)
  end forall
  
  l=17
  ddx=dex(l)/tauR
  ddy=dey(l)/tauR
  ddz=dez(l)/tauR
  forall(i=1:nx,j=1:ny,k=1:nz)
    rhoR(i,j,k) = f17R(i,j,k) + rhoR(i,j,k)
    v(i,j,k)    = f17R(i,j,k)*ddy + v(i,j,k)
    w(i,j,k)    = f17R(i,j,k)*ddz + w(i,j,k)
  end forall
  
  l=18
  ddx=dex(l)/tauR
  ddy=dey(l)/tauR
  ddz=dez(l)/tauR
  forall(i=1:nx,j=1:ny,k=1:nz)
    rhoR(i,j,k) = f18R(i,j,k) + rhoR(i,j,k)
    v(i,j,k)    = f18R(i,j,k)*ddy + v(i,j,k)
    w(i,j,k)    = f18R(i,j,k)*ddz + w(i,j,k)
  end forall
  
  !blue fluid
  
  l=0
  ddx=dex(l)/tauB
  ddy=dey(l)/tauB
  ddz=dez(l)/tauB
  forall(i=1:nx,j=1:ny,k=1:nz)
    rhoB(i,j,k) = f00B(i,j,k)
  end forall
  
  l=1
  ddx=dex(l)/tauB
  ddy=dey(l)/tauB
  ddz=dez(l)/tauB
  forall(i=1:nx,j=1:ny,k=1:nz)
    rhoB(i,j,k) = f01B(i,j,k) + rhoB(i,j,k)
    u(i,j,k)    = f01B(i,j,k)*ddx + u(i,j,k)
  end forall
  
  l=2
  ddx=dex(l)/tauB
  ddy=dey(l)/tauB
  ddz=dez(l)/tauB
  forall(i=1:nx,j=1:ny,k=1:nz)
    rhoB(i,j,k) = f02B(i,j,k) + rhoB(i,j,k)
    u(i,j,k)    = f02B(i,j,k)*ddx + u(i,j,k)
  end forall
  
  l=3
  ddx=dex(l)/tauB
  ddy=dey(l)/tauB
  ddz=dez(l)/tauB
  forall(i=1:nx,j=1:ny,k=1:nz)
    rhoB(i,j,k) = f03B(i,j,k) + rhoB(i,j,k)
    v(i,j,k)    = f03B(i,j,k)*ddy + v(i,j,k)
  end forall
  
  l=4
  ddx=dex(l)/tauB
  ddy=dey(l)/tauB
  ddz=dez(l)/tauB
  forall(i=1:nx,j=1:ny,k=1:nz)
    rhoB(i,j,k) = f04B(i,j,k) + rhoB(i,j,k)
    v(i,j,k)    = f04B(i,j,k)*ddy + v(i,j,k)
  end forall
  
  l=5
  ddx=dex(l)/tauB
  ddy=dey(l)/tauB
  ddz=dez(l)/tauB
  forall(i=1:nx,j=1:ny,k=1:nz)
    rhoB(i,j,k) = f05B(i,j,k) + rhoB(i,j,k)
    w(i,j,k)    = f05B(i,j,k)*ddz + w(i,j,k)
  end forall
  
  l=6
  ddx=dex(l)/tauB
  ddy=dey(l)/tauB
  ddz=dez(l)/tauB
  forall(i=1:nx,j=1:ny,k=1:nz)
    rhoB(i,j,k) = f06B(i,j,k) + rhoB(i,j,k)
    w(i,j,k)    = f06B(i,j,k)*ddz + w(i,j,k)
  end forall
  
  l=7
  ddx=dex(l)/tauB
  ddy=dey(l)/tauB
  ddz=dez(l)/tauB
  forall(i=1:nx,j=1:ny,k=1:nz)
    rhoB(i,j,k) = f07B(i,j,k) + rhoB(i,j,k)
    u(i,j,k)    = f07B(i,j,k)*ddx + u(i,j,k)
    v(i,j,k)    = f07B(i,j,k)*ddy + v(i,j,k)
  end forall
  
  l=8
  ddx=dex(l)/tauB
  ddy=dey(l)/tauB
  ddz=dez(l)/tauB
  forall(i=1:nx,j=1:ny,k=1:nz)
    rhoB(i,j,k) = f08B(i,j,k) + rhoB(i,j,k)
    u(i,j,k)    = f08B(i,j,k)*ddx + u(i,j,k)
    v(i,j,k)    = f08B(i,j,k)*ddy + v(i,j,k)
  end forall
  
  l=9
  ddx=dex(l)/tauB
  ddy=dey(l)/tauB
  ddz=dez(l)/tauB
  forall(i=1:nx,j=1:ny,k=1:nz)
    rhoB(i,j,k) = f09B(i,j,k) + rhoB(i,j,k)
    u(i,j,k)    = f09B(i,j,k)*ddx + u(i,j,k)
    v(i,j,k)    = f09B(i,j,k)*ddy + v(i,j,k)
  end forall
  
  l=10
  ddx=dex(l)/tauB
  ddy=dey(l)/tauB
  ddz=dez(l)/tauB
  forall(i=1:nx,j=1:ny,k=1:nz)
    rhoB(i,j,k) = f10B(i,j,k) + rhoB(i,j,k)
    u(i,j,k)    = f10B(i,j,k)*ddx + u(i,j,k)
    v(i,j,k)    = f10B(i,j,k)*ddy + v(i,j,k)
  end forall
  
  l=11
  ddx=dex(l)/tauB
  ddy=dey(l)/tauB
  ddz=dez(l)/tauB
  forall(i=1:nx,j=1:ny,k=1:nz)
    rhoB(i,j,k) = f11B(i,j,k) + rhoB(i,j,k)
    u(i,j,k)    = f11B(i,j,k)*ddx + u(i,j,k)
    w(i,j,k)    = f11B(i,j,k)*ddz + w(i,j,k)
  end forall
  
  l=12
  ddx=dex(l)/tauB
  ddy=dey(l)/tauB
  ddz=dez(l)/tauB
  forall(i=1:nx,j=1:ny,k=1:nz)
    rhoB(i,j,k) = f12B(i,j,k) + rhoB(i,j,k)
    u(i,j,k)    = f12B(i,j,k)*ddx + u(i,j,k)
    w(i,j,k)    = f12B(i,j,k)*ddz + w(i,j,k)
  end forall
  
  l=13
  ddx=dex(l)/tauB
  ddy=dey(l)/tauB
  ddz=dez(l)/tauB
  forall(i=1:nx,j=1:ny,k=1:nz)
    rhoB(i,j,k) = f13B(i,j,k) + rhoB(i,j,k)
    u(i,j,k)    = f13B(i,j,k)*ddx + u(i,j,k)
    w(i,j,k)    = f13B(i,j,k)*ddz + w(i,j,k)
  end forall
  
  l=14
  ddx=dex(l)/tauB
  ddy=dey(l)/tauB
  ddz=dez(l)/tauB
  forall(i=1:nx,j=1:ny,k=1:nz)
    rhoB(i,j,k) = f14B(i,j,k) + rhoB(i,j,k)
    u(i,j,k)    = f14B(i,j,k)*ddx + u(i,j,k)
    w(i,j,k)    = f14B(i,j,k)*ddz + w(i,j,k)
  end forall
  
  l=15
  ddx=dex(l)/tauB
  ddy=dey(l)/tauB
  ddz=dez(l)/tauB
  forall(i=1:nx,j=1:ny,k=1:nz)
    rhoB(i,j,k) = f15B(i,j,k) + rhoB(i,j,k)
    v(i,j,k)    = f15B(i,j,k)*ddy + v(i,j,k)
    w(i,j,k)    = f15B(i,j,k)*ddz + w(i,j,k)
  end forall
  
  l=16
  ddx=dex(l)/tauB
  ddy=dey(l)/tauB
  ddz=dez(l)/tauB
  forall(i=1:nx,j=1:ny,k=1:nz)
    rhoB(i,j,k) = f16B(i,j,k) + rhoB(i,j,k)
    v(i,j,k)    = f16B(i,j,k)*ddy + v(i,j,k)
    w(i,j,k)    = f16B(i,j,k)*ddz + w(i,j,k)
  end forall
  
  l=17
  ddx=dex(l)/tauB
  ddy=dey(l)/tauB
  ddz=dez(l)/tauB
  forall(i=1:nx,j=1:ny,k=1:nz)
    rhoB(i,j,k) = f17B(i,j,k) + rhoB(i,j,k)
    v(i,j,k)    = f17B(i,j,k)*ddy + v(i,j,k)
    w(i,j,k)    = f17B(i,j,k)*ddz + w(i,j,k)
  end forall
  
  l=18
  ddx=dex(l)/tauB
  ddy=dey(l)/tauB
  ddz=dez(l)/tauB
  forall(i=1:nx,j=1:ny,k=1:nz)
    rhoB(i,j,k) = f18B(i,j,k) + rhoB(i,j,k)
    v(i,j,k)    = f18B(i,j,k)*ddy + v(i,j,k)
    w(i,j,k)    = f18B(i,j,k)*ddz + w(i,j,k)
  end forall
  
  !compute speed from mass flux
  forall(i=1:nx,j=1:ny,k=1:nz)
    u(i,j,k) = u(i,j,k)/(rhoR(i,j,k)/tauR + rhoB(i,j,k)/tauB)
    v(i,j,k) = v(i,j,k)/(rhoR(i,j,k)/tauR + rhoB(i,j,k)/tauB)
    w(i,j,k) = w(i,j,k)/(rhoR(i,j,k)/tauR + rhoB(i,j,k)/tauB)
  end forall
    
  return
 
 end subroutine moments_fluids
 
 subroutine probe_red_moments_in_node(i,j,k,dtemp1,dtemp2,dtemp3,dtemp4)
 
!***********************************************************************
!     
!     LBsoft subroutine for probing the Red hydrodynamic moments
!     at the point i j k
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: i,j,k
  real(kind=PRC), intent(out) :: dtemp1,dtemp2,dtemp3,dtemp4
  
  integer :: l
  
  dtemp1=ZERO
  dtemp2=ZERO
  dtemp3=ZERO
  dtemp4=ZERO
  
  l=0
  dtemp1 = f00R(i,j,k) + dtemp1 
  dtemp2    = f00R(i,j,k)*dex(l) + dtemp2
  dtemp3    = f00R(i,j,k)*dey(l) + dtemp3
  dtemp4    = f00R(i,j,k)*dez(l) + dtemp4
  
  l=1
  dtemp1 = f01R(i,j,k) + dtemp1
  dtemp2    = f01R(i,j,k)*dex(l) + dtemp2
  dtemp3    = f01R(i,j,k)*dey(l) + dtemp3
  dtemp4    = f01R(i,j,k)*dez(l) + dtemp4
  
  l=2
  dtemp1 = f02R(i,j,k) + dtemp1
  dtemp2    = f02R(i,j,k)*dex(l) + dtemp2
  dtemp3    = f02R(i,j,k)*dey(l) + dtemp3
  dtemp4    = f02R(i,j,k)*dez(l) + dtemp4
  
  l=3
  dtemp1 = f03R(i,j,k) + dtemp1
  dtemp2    = f03R(i,j,k)*dex(l) + dtemp2
  dtemp3    = f03R(i,j,k)*dey(l) + dtemp3
  dtemp4    = f03R(i,j,k)*dez(l) + dtemp4
  
  l=4
  dtemp1 = f04R(i,j,k) + dtemp1
  dtemp2    = f04R(i,j,k)*dex(l) + dtemp2
  dtemp3    = f04R(i,j,k)*dey(l) + dtemp3
  dtemp4    = f04R(i,j,k)*dez(l) + dtemp4
  
  l=5
  dtemp1 = f05R(i,j,k) + dtemp1
  dtemp2    = f05R(i,j,k)*dex(l) + dtemp2
  dtemp3    = f05R(i,j,k)*dey(l) + dtemp3
  dtemp4    = f05R(i,j,k)*dez(l) + dtemp4
  
  l=6
  dtemp1 = f06R(i,j,k) + dtemp1
  dtemp2    = f06R(i,j,k)*dex(l) + dtemp2
  dtemp3    = f06R(i,j,k)*dey(l) + dtemp3
  dtemp4    = f06R(i,j,k)*dez(l) + dtemp4
  
  l=7
  dtemp1 = f07R(i,j,k) + dtemp1
  dtemp2    = f07R(i,j,k)*dex(l) + dtemp2
  dtemp3    = f07R(i,j,k)*dey(l) + dtemp3
  dtemp4    = f07R(i,j,k)*dez(l) + dtemp4
  
  l=8
  dtemp1 = f08R(i,j,k) + dtemp1
  dtemp2    = f08R(i,j,k)*dex(l) + dtemp2
  dtemp3    = f08R(i,j,k)*dey(l) + dtemp3
  dtemp4    = f08R(i,j,k)*dez(l) + dtemp4
  
  l=9
  dtemp1 = f09R(i,j,k) + dtemp1
  dtemp2    = f09R(i,j,k)*dex(l) + dtemp2
  dtemp3    = f09R(i,j,k)*dey(l) + dtemp3
  dtemp4    = f09R(i,j,k)*dez(l) + dtemp4
  
  l=10
  dtemp1 = f10R(i,j,k) + dtemp1
  dtemp2    = f10R(i,j,k)*dex(l) + dtemp2
  dtemp3    = f10R(i,j,k)*dey(l) + dtemp3
  dtemp4    = f10R(i,j,k)*dez(l) + dtemp4
  
  l=11
  dtemp1 = f11R(i,j,k) + dtemp1
  dtemp2    = f11R(i,j,k)*dex(l) + dtemp2
  dtemp3    = f11R(i,j,k)*dey(l) + dtemp3
  dtemp4    = f11R(i,j,k)*dez(l) + dtemp4
  
  l=12
  dtemp1 = f12R(i,j,k) + dtemp1
  dtemp2    = f12R(i,j,k)*dex(l) + dtemp2
  dtemp3    = f12R(i,j,k)*dey(l) + dtemp3
  dtemp4    = f12R(i,j,k)*dez(l) + dtemp4
  
  l=13
  dtemp1 = f13R(i,j,k) + dtemp1
  dtemp2    = f13R(i,j,k)*dex(l) + dtemp2
  dtemp3    = f13R(i,j,k)*dey(l) + dtemp3
  dtemp4    = f13R(i,j,k)*dez(l) + dtemp4
  
  l=14
  dtemp1 = f14R(i,j,k) + dtemp1
  dtemp2    = f14R(i,j,k)*dex(l) + dtemp2
  dtemp3    = f14R(i,j,k)*dey(l) + dtemp3
  dtemp4    = f14R(i,j,k)*dez(l) + dtemp4
  
  l=15
  dtemp1 = f15R(i,j,k) + dtemp1
  dtemp2    = f15R(i,j,k)*dex(l) + dtemp2
  dtemp3    = f15R(i,j,k)*dey(l) + dtemp3
  dtemp4    = f15R(i,j,k)*dez(l) + dtemp4
  
  l=16
  dtemp1 = f16R(i,j,k) + dtemp1
  dtemp2    = f16R(i,j,k)*dex(l) + dtemp2
  dtemp3    = f16R(i,j,k)*dey(l) + dtemp3
  dtemp4    = f16R(i,j,k)*dez(l) + dtemp4
  
  l=17
  dtemp1 = f17R(i,j,k) + dtemp1
  dtemp2    = f17R(i,j,k)*dex(l) + dtemp2
  dtemp3    = f17R(i,j,k)*dey(l) + dtemp3
  dtemp4    = f17R(i,j,k)*dez(l) + dtemp4
  
  l=18
  dtemp1 = f18R(i,j,k) + dtemp1
  dtemp2    = f18R(i,j,k)*dex(l) + dtemp2
  dtemp3    = f18R(i,j,k)*dey(l) + dtemp3
  dtemp4    = f18R(i,j,k)*dez(l) + dtemp4
  
  return
  
 end subroutine probe_red_moments_in_node
 
 subroutine probe_blue_moments_in_node(i,j,k,dtemp1,dtemp2,dtemp3,dtemp4)
 
!***********************************************************************
!     
!     LBsoft subroutine for probing the Blue hydrodynamic moments
!     at the point i j k
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: i,j,k
  real(kind=PRC), intent(out) :: dtemp1,dtemp2,dtemp3,dtemp4
  
  integer :: l
  
  dtemp1=ZERO
  dtemp2=ZERO
  dtemp3=ZERO
  dtemp4=ZERO
  
  l=0
  dtemp1 = f00B(i,j,k) + dtemp1 
  dtemp2    = f00B(i,j,k)*dex(l) + dtemp2
  dtemp3    = f00B(i,j,k)*dey(l) + dtemp3
  dtemp4    = f00B(i,j,k)*dez(l) + dtemp4
  
  l=1
  dtemp1 = f01B(i,j,k) + dtemp1
  dtemp2    = f01B(i,j,k)*dex(l) + dtemp2
  dtemp3    = f01B(i,j,k)*dey(l) + dtemp3
  dtemp4    = f01B(i,j,k)*dez(l) + dtemp4
  
  l=2
  dtemp1 = f02B(i,j,k) + dtemp1
  dtemp2    = f02B(i,j,k)*dex(l) + dtemp2
  dtemp3    = f02B(i,j,k)*dey(l) + dtemp3
  dtemp4    = f02B(i,j,k)*dez(l) + dtemp4
  
  l=3
  dtemp1 = f03B(i,j,k) + dtemp1
  dtemp2    = f03B(i,j,k)*dex(l) + dtemp2
  dtemp3    = f03B(i,j,k)*dey(l) + dtemp3
  dtemp4    = f03B(i,j,k)*dez(l) + dtemp4
  
  l=4
  dtemp1 = f04B(i,j,k) + dtemp1
  dtemp2    = f04B(i,j,k)*dex(l) + dtemp2
  dtemp3    = f04B(i,j,k)*dey(l) + dtemp3
  dtemp4    = f04B(i,j,k)*dez(l) + dtemp4
  
  l=5
  dtemp1 = f05B(i,j,k) + dtemp1
  dtemp2    = f05B(i,j,k)*dex(l) + dtemp2
  dtemp3    = f05B(i,j,k)*dey(l) + dtemp3
  dtemp4    = f05B(i,j,k)*dez(l) + dtemp4
  
  l=6
  dtemp1 = f06B(i,j,k) + dtemp1
  dtemp2    = f06B(i,j,k)*dex(l) + dtemp2
  dtemp3    = f06B(i,j,k)*dey(l) + dtemp3
  dtemp4    = f06B(i,j,k)*dez(l) + dtemp4
  
  l=7
  dtemp1 = f07B(i,j,k) + dtemp1
  dtemp2    = f07B(i,j,k)*dex(l) + dtemp2
  dtemp3    = f07B(i,j,k)*dey(l) + dtemp3
  dtemp4    = f07B(i,j,k)*dez(l) + dtemp4
  
  l=8
  dtemp1 = f08B(i,j,k) + dtemp1
  dtemp2    = f08B(i,j,k)*dex(l) + dtemp2
  dtemp3    = f08B(i,j,k)*dey(l) + dtemp3
  dtemp4    = f08B(i,j,k)*dez(l) + dtemp4
  
  l=9
  dtemp1 = f09B(i,j,k) + dtemp1
  dtemp2    = f09B(i,j,k)*dex(l) + dtemp2
  dtemp3    = f09B(i,j,k)*dey(l) + dtemp3
  dtemp4    = f09B(i,j,k)*dez(l) + dtemp4
  
  l=10
  dtemp1 = f10B(i,j,k) + dtemp1
  dtemp2    = f10B(i,j,k)*dex(l) + dtemp2
  dtemp3    = f10B(i,j,k)*dey(l) + dtemp3
  dtemp4    = f10B(i,j,k)*dez(l) + dtemp4
  
  l=11
  dtemp1 = f11B(i,j,k) + dtemp1
  dtemp2    = f11B(i,j,k)*dex(l) + dtemp2
  dtemp3    = f11B(i,j,k)*dey(l) + dtemp3
  dtemp4    = f11B(i,j,k)*dez(l) + dtemp4
  
  l=12
  dtemp1 = f12B(i,j,k) + dtemp1
  dtemp2    = f12B(i,j,k)*dex(l) + dtemp2
  dtemp3    = f12B(i,j,k)*dey(l) + dtemp3
  dtemp4    = f12B(i,j,k)*dez(l) + dtemp4
  
  l=13
  dtemp1 = f13B(i,j,k) + dtemp1
  dtemp2    = f13B(i,j,k)*dex(l) + dtemp2
  dtemp3    = f13B(i,j,k)*dey(l) + dtemp3
  dtemp4    = f13B(i,j,k)*dez(l) + dtemp4
  
  l=14
  dtemp1 = f14B(i,j,k) + dtemp1
  dtemp2    = f14B(i,j,k)*dex(l) + dtemp2
  dtemp3    = f14B(i,j,k)*dey(l) + dtemp3
  dtemp4    = f14B(i,j,k)*dez(l) + dtemp4
  
  l=15
  dtemp1 = f15B(i,j,k) + dtemp1
  dtemp2    = f15B(i,j,k)*dex(l) + dtemp2
  dtemp3    = f15B(i,j,k)*dey(l) + dtemp3
  dtemp4    = f15B(i,j,k)*dez(l) + dtemp4
  
  l=16
  dtemp1 = f16B(i,j,k) + dtemp1
  dtemp2    = f16B(i,j,k)*dex(l) + dtemp2
  dtemp3    = f16B(i,j,k)*dey(l) + dtemp3
  dtemp4    = f16B(i,j,k)*dez(l) + dtemp4
  
  l=17
  dtemp1 = f17B(i,j,k) + dtemp1
  dtemp2    = f17B(i,j,k)*dex(l) + dtemp2
  dtemp3    = f17B(i,j,k)*dey(l) + dtemp3
  dtemp4    = f17B(i,j,k)*dez(l) + dtemp4
  
  l=18
  dtemp1 = f18B(i,j,k) + dtemp1
  dtemp2    = f18B(i,j,k)*dex(l) + dtemp2
  dtemp3    = f18B(i,j,k)*dey(l) + dtemp3
  dtemp4    = f18B(i,j,k)*dez(l) + dtemp4
  
  return
  
 end subroutine probe_blue_moments_in_node
 
!*************START PART TO DEFINE THE EQUILIBRIUM FUNCTS***************
 
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
  
#if LATTICE==319
  uv=(ONE/mycssq)*(myu*mydex)
  equil_pop01=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#else
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop01=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#endif
  
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
  
#if LATTICE==319
  uv=(ONE/mycssq)*(myu*mydex)
  equil_pop02=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#else
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop02=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#endif
  
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
  
#if LATTICE==319
  uv=(ONE/mycssq)*(myv*mydey)
  equil_pop03=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#else
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop03=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#endif
  
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
  
#if LATTICE==319
  uv=(ONE/mycssq)*(myv*mydey)
  equil_pop04=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#else
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop04=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#endif
  
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
  
#if LATTICE==319
  uv=(ONE/mycssq)*(myw*mydez)
  equil_pop05=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#else
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop05=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#endif
  
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
  
#if LATTICE==319
  uv=(ONE/mycssq)*(myw*mydez)
  equil_pop06=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#else
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop06=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#endif
  
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
  
#if LATTICE==319
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey)
  equil_pop07=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#else
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop07=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#endif
  
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
  
#if LATTICE==319
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey)
  equil_pop08=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#else
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop08=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#endif
  
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
  
#if LATTICE==319
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey)
  equil_pop09=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#else
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop09=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#endif
  
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
  
#if LATTICE==319
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey)
  equil_pop10=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#else
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop10=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#endif
  
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
  
#if LATTICE==319
  uv=(ONE/mycssq)*(myu*mydex + myw*mydez)
  equil_pop11=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#else
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop11=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#endif
  
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
  
#if LATTICE==319
  uv=(ONE/mycssq)*(myu*mydex  + myw*mydez)
  equil_pop12=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#else
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop12=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#endif
  
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
  
#if LATTICE==319
  uv=(ONE/mycssq)*(myu*mydex + myw*mydez)
  equil_pop13=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#else
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop13=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#endif
  
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
  
#if LATTICE==319
  uv=(ONE/mycssq)*(myu*mydex + myw*mydez)
  equil_pop14=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#else
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop14=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#endif
  
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
  
#if LATTICE==319
  uv=(ONE/mycssq)*(myv*mydey + myw*mydez)
  equil_pop15=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#else
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop15=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#endif
  
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
  
#if LATTICE==319
  uv=(ONE/mycssq)*(myv*mydey + myw*mydez)
  equil_pop16=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#else
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop16=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#endif
  
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
  
#if LATTICE==319
  uv=(ONE/mycssq)*(myv*mydey + myw*mydez)
  equil_pop17=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#else
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop17=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#endif
  
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
  
#if LATTICE==319
  uv=(ONE/mycssq)*(myv*mydey + myw*mydez)
  equil_pop18=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#else
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop18=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#endif
  
  return
  
 end function equil_pop18

!**************END PART TO DEFINE THE EQUILIBRIUM FUNCTS****************

!*****************START PART TO MANAGE THE PERIODIC BC******************

subroutine driver_bc_densities
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the boundary conditions
!     to fluid densities
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  select case(ibctype)
  case(0) ! 0 0 0
    return
  case(1) ! 1 0 0
    call apply_pbc_densities_along_x
  case(2) ! 0 1 0
    call apply_pbc_densities_along_y
  case(3) ! 1 1 0 
    call apply_pbc_densities_along_xy
  case(4) ! 0 0 1 
    call apply_pbc_densities_along_z
  case(5) ! 1 0 1 
    call apply_pbc_densities_along_xz
  case(6) ! 0 1 1 
    call apply_pbc_densities_along_yz
  case(7) ! 1 1 1
    call apply_pbc_densities
  case default
    call error(12)
  end select
  
  return
  
 end subroutine driver_bc_densities
 
 subroutine driver_bc_velocities
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the boundary conditions
!     to fluid densities
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  select case(ibctype)
  case(0) ! 0 0 0
    return
  case(1) ! 1 0 0
    call apply_pbc_velocities_along_x
  case(2) ! 0 1 0
    call apply_pbc_velocities_along_y
  case(3) ! 1 1 0 
    call apply_pbc_velocities_along_xy
  case(4) ! 0 0 1 
    call apply_pbc_velocities_along_z
  case(5) ! 1 0 1 
    call apply_pbc_velocities_along_xz
  case(6) ! 0 1 1 
    call apply_pbc_velocities_along_yz
  case(7) ! 1 1 1
    call apply_pbc_velocities
  case default
    call error(12)
  end select
  
  return
  
 end subroutine driver_bc_velocities
 
 subroutine apply_pbc_densities
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the full periodic boundary 
!     conditions to fluid densities
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  !sides 6
  
  !apply pbc at north 1   !z
  call apply_pbc_north(rhoR)
  
  !apply pbc at south 2   !z
  call apply_pbc_south(rhoR)
  
  !apply pbc at east 3    !x
  call apply_pbc_east(rhoR)
  
  !apply pbc at west 4    !x
  call apply_pbc_west(rhoR)
  
  !apply pbc at front 5   !y
  call apply_pbc_front(rhoR)
  
  !apply pbc at rear 6   !y
  call apply_pbc_rear(rhoR)
  
  !edges 12
  
  !apply pbc at front east 1   !xy
  call apply_pbc_edge_front_east(rhoR)
  
  !apply pbc at front west 2   !xy
  call apply_pbc_edge_front_west(rhoR)
  
  !apply pbc at north east 3   !xz
  call apply_pbc_edge_north_east(rhoR)
  
  !apply pbc at north front 4   !yz
  call apply_pbc_edge_north_front(rhoR)
  
  !apply pbc at north rear 5   !yz
  call apply_pbc_edge_north_rear(rhoR)
  
  !apply pbc at north west 6   !xz
  call apply_pbc_edge_north_west(rhoR)
  
  !apply pbc at rear east 7   !xy
  call apply_pbc_edge_rear_east(rhoR)
  
  !apply pbc at rear west 8   !xy
  call apply_pbc_edge_rear_west(rhoR)
  
  !apply pbc at south east 9  !xz
  call apply_pbc_edge_south_east(rhoR)
  
  !apply pbc at south front 10  !yz
  call apply_pbc_edge_south_front(rhoR)
  
  !apply pbc at south rear 11  !yz
  call apply_pbc_edge_south_rear(rhoR)
  
  !apply pbc at south west 12  !xz
  call apply_pbc_edge_south_west(rhoR)
  
  !corner 8
  
  !apply pbc at north east front 1
  call apply_pbc_corner_north_east_front(rhoR)
  
  !apply pbc at north east rear 2
  call apply_pbc_corner_north_east_rear(rhoR)
  
  !apply pbc at north west rear 3
  call apply_pbc_corner_north_west_rear(rhoR)
  
  !apply pbc at north west front 4
  call apply_pbc_corner_north_west_front(rhoR)
  
  !apply pbc at south east front 5
  call apply_pbc_corner_south_east_front(rhoR)
  
  !apply pbc at south west front 6
  call apply_pbc_corner_south_west_front(rhoR)
  
  !apply pbc at south west rear 7
  call apply_pbc_corner_south_west_rear(rhoR)
  
  !apply pbc at south east rear 8
  call apply_pbc_corner_south_east_rear(rhoR)
  
  if(lsingle_fluid)return
  
  !sides 6
  
  !apply pbc at north 1   !z
  call apply_pbc_north(rhoB)
  
  !apply pbc at south 2   !z
  call apply_pbc_south(rhoB)
  
  !apply pbc at east 3    !x
  call apply_pbc_east(rhoB)
  
  !apply pbc at west 4    !x
  call apply_pbc_west(rhoB)
  
  !apply pbc at front 5   !y
  call apply_pbc_front(rhoB)
  
  !apply pbc at rear 6   !y
  call apply_pbc_rear(rhoB)
  
  !edges 12
  
  !apply pbc at front east 1   !xy
  call apply_pbc_edge_front_east(rhoB)
  
  !apply pbc at front west 2   !xy
  call apply_pbc_edge_front_west(rhoB)
  
  !apply pbc at north east 3   !xz
  call apply_pbc_edge_north_east(rhoB)
  
  !apply pbc at north front 4   !yz
  call apply_pbc_edge_north_front(rhoB)
  
  !apply pbc at north rear 5   !yz
  call apply_pbc_edge_north_rear(rhoB)
  
  !apply pbc at north west 6   !xz
  call apply_pbc_edge_north_west(rhoB)
  
  !apply pbc at rear east 7   !xy
  call apply_pbc_edge_rear_east(rhoB)
  
  !apply pbc at rear west 8   !xy
  call apply_pbc_edge_rear_west(rhoB)
  
  !apply pbc at south east 9  !xz
  call apply_pbc_edge_south_east(rhoB)
  
  !apply pbc at south front 10  !yz
  call apply_pbc_edge_south_front(rhoB)
  
  !apply pbc at south rear 11  !yz
  call apply_pbc_edge_south_rear(rhoB)
  
  !apply pbc at south west 12  !xz
  call apply_pbc_edge_south_west(rhoB)
  
  !corner 8
  
  !apply pbc at north east front 1
  call apply_pbc_corner_north_east_front(rhoB)
  
  !apply pbc at north east rear 2
  call apply_pbc_corner_north_east_rear(rhoB)
  
  !apply pbc at north west rear 3
  call apply_pbc_corner_north_west_rear(rhoB)
  
  !apply pbc at north west front 4
  call apply_pbc_corner_north_west_front(rhoB)
  
  !apply pbc at south east front 5
  call apply_pbc_corner_south_east_front(rhoB)
  
  !apply pbc at south west front 6
  call apply_pbc_corner_south_west_front(rhoB)
  
  !apply pbc at south west rear 7
  call apply_pbc_corner_south_west_rear(rhoB)
  
  !apply pbc at south east rear 8
  call apply_pbc_corner_south_east_rear(rhoB)
  
  
  
  return
  
 end subroutine apply_pbc_densities
 
 subroutine apply_pbc_velocities
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the full periodic boundary 
!     conditions to fluid velocities
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  !sides 6
  
  !apply pbc at north 1   !z
  call apply_pbc_north(u)
  call apply_pbc_north(v)
  call apply_pbc_north(w)
  
  !apply pbc at south 2   !z
  call apply_pbc_south(u)
  call apply_pbc_south(v)
  call apply_pbc_south(w)
  
  !apply pbc at east 3    !x
  call apply_pbc_east(u)
  call apply_pbc_east(v)
  call apply_pbc_east(w)
  
  !apply pbc at west 4    !x
  call apply_pbc_west(u)
  call apply_pbc_west(v)
  call apply_pbc_west(w)
  
  !apply pbc at front 5   !y
  call apply_pbc_front(u)
  call apply_pbc_front(v)
  call apply_pbc_front(w)
  
  !apply pbc at rear 6   !y
  call apply_pbc_rear(u)
  call apply_pbc_rear(v)
  call apply_pbc_rear(w)
  
  !edges 12
  
  !apply pbc at front east 1   !xy
  call apply_pbc_edge_front_east(u)
  call apply_pbc_edge_front_east(v)
  call apply_pbc_edge_front_east(w)
  
  !apply pbc at front west 2   !xy
  call apply_pbc_edge_front_west(u)
  call apply_pbc_edge_front_west(v)
  call apply_pbc_edge_front_west(w)
  
  !apply pbc at north east 3   !xz
  call apply_pbc_edge_north_east(u)
  call apply_pbc_edge_north_east(v)
  call apply_pbc_edge_north_east(w)
  
  !apply pbc at north front 4   !yz
  call apply_pbc_edge_north_front(u)
  call apply_pbc_edge_north_front(v)
  call apply_pbc_edge_north_front(w)
  
  !apply pbc at north rear 5   !yz
  call apply_pbc_edge_north_rear(u)
  call apply_pbc_edge_north_rear(v)
  call apply_pbc_edge_north_rear(w)
  
  !apply pbc at north west 6   !xz
  call apply_pbc_edge_north_west(u)
  call apply_pbc_edge_north_west(v)
  call apply_pbc_edge_north_west(w)
  
  !apply pbc at rear east 7   !xy
  call apply_pbc_edge_rear_east(u)
  call apply_pbc_edge_rear_east(v)
  call apply_pbc_edge_rear_east(w)
  
  !apply pbc at rear west 8   !xy
  call apply_pbc_edge_rear_west(u)
  call apply_pbc_edge_rear_west(v)
  call apply_pbc_edge_rear_west(w)
  
  !apply pbc at south east 9  !xz
  call apply_pbc_edge_south_east(u)
  call apply_pbc_edge_south_east(v)
  call apply_pbc_edge_south_east(w)
  
  !apply pbc at south front 10  !yz
  call apply_pbc_edge_south_front(u)
  call apply_pbc_edge_south_front(v)
  call apply_pbc_edge_south_front(w)
  
  !apply pbc at south rear 11  !yz
  call apply_pbc_edge_south_rear(u)
  call apply_pbc_edge_south_rear(v)
  call apply_pbc_edge_south_rear(w)
  
  !apply pbc at south west 12  !xz
  call apply_pbc_edge_south_west(u)
  call apply_pbc_edge_south_west(v)
  call apply_pbc_edge_south_west(w)
  
  !corner 8
  
  !apply pbc at north east front 1
  call apply_pbc_corner_north_east_front(u)
  call apply_pbc_corner_north_east_front(v)
  call apply_pbc_corner_north_east_front(w)
  
  !apply pbc at north east rear 2
  call apply_pbc_corner_north_east_rear(u)
  call apply_pbc_corner_north_east_rear(v)
  call apply_pbc_corner_north_east_rear(w)
  
  !apply pbc at north west rear 3
  call apply_pbc_corner_north_west_rear(u)
  call apply_pbc_corner_north_west_rear(v)
  call apply_pbc_corner_north_west_rear(w)
  
  !apply pbc at north west front 4
  call apply_pbc_corner_north_west_front(u)
  call apply_pbc_corner_north_west_front(v)
  call apply_pbc_corner_north_west_front(w)
  
  !apply pbc at south east front 5
  call apply_pbc_corner_south_east_front(u)
  call apply_pbc_corner_south_east_front(v)
  call apply_pbc_corner_south_east_front(w)
  
  !apply pbc at south west front 6
  call apply_pbc_corner_south_west_front(u)
  call apply_pbc_corner_south_west_front(v)
  call apply_pbc_corner_south_west_front(w)
  
  !apply pbc at south west rear 7
  call apply_pbc_corner_south_west_rear(u)
  call apply_pbc_corner_south_west_rear(v)
  call apply_pbc_corner_south_west_rear(w)
  
  !apply pbc at south east rear 8
  call apply_pbc_corner_south_east_rear(u)
  call apply_pbc_corner_south_east_rear(v)
  call apply_pbc_corner_south_east_rear(w)
  
  return
  
 end subroutine apply_pbc_velocities
 
 subroutine apply_pbc_densities_along_x
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the periodic boundary 
!     conditions to fluid densities along x
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  !sides 2 (6)
  
  !apply pbc at east 3    !x
  call apply_pbc_east(rhoR)
  
  !apply pbc at west 4    !x
  call apply_pbc_west(rhoR)
  
  if(lsingle_fluid)return
  
  !sides 2 (6)
  
  !apply pbc at east 3    !x
  call apply_pbc_east(rhoB)
  
  !apply pbc at west 4    !x
  call apply_pbc_west(rhoB)
    
  return
  
 end subroutine apply_pbc_densities_along_x
 
 subroutine apply_pbc_velocities_along_x
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the periodic boundary 
!     conditions to fluid velocities along x
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  !sides 2 (6)
  
  !apply pbc at east 3    !x
  call apply_pbc_east(u)
  call apply_pbc_east(v)
  call apply_pbc_east(w)
  
  !apply pbc at west 4    !x
  call apply_pbc_west(u)
  call apply_pbc_west(v)
  call apply_pbc_west(w)
  
  return
  
 end subroutine apply_pbc_velocities_along_x
 
 subroutine apply_pbc_densities_along_y
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the periodic boundary 
!     conditions to fluid densities along y
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  !sides 2 (6)
  
  !apply pbc at front 5   !y
  call apply_pbc_front(rhoR)
  
  !apply pbc at rear 6    !y
  call apply_pbc_rear(rhoR)
  
  if(lsingle_fluid)return
  
  !sides 2 (6)
  
  !apply pbc at front 5   !y
  call apply_pbc_front(rhoB)
  
  !apply pbc at rear 6    !y
  call apply_pbc_rear(rhoB)
  
  return
  
 end subroutine apply_pbc_densities_along_y
 
 subroutine apply_pbc_velocities_along_y
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the periodic boundary 
!     conditions to fluid velocities along y
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  !sides 2 (6)
  
  !apply pbc at front 5   !y
  call apply_pbc_front(u)
  call apply_pbc_front(v)
  call apply_pbc_front(w)
  
  !apply pbc at rear 6    !y
  call apply_pbc_rear(u)
  call apply_pbc_rear(v)
  call apply_pbc_rear(w)
  
  
  return
  
 end subroutine apply_pbc_velocities_along_y
 
 subroutine apply_pbc_densities_along_z
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the full periodic boundary 
!     conditions to fluid densities along z
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  !sides 6 (4)
  
  !apply pbc at north 1   !z
  call apply_pbc_north(rhoR)
  
  !apply pbc at south 2   !z
  call apply_pbc_south(rhoR)
  
  if(lsingle_fluid)return
  
  !sides 6 (4)
  
  !apply pbc at north 1   !z
  call apply_pbc_north(rhoB)
  
  !apply pbc at south 2   !z
  call apply_pbc_south(rhoB)
  
  return
  
 end subroutine apply_pbc_densities_along_z
 
 subroutine apply_pbc_velocities_along_z
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the full periodic boundary 
!     conditions to fluid velocities along z
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  !sides 2 (6)
  
  !apply pbc at north 1   !z
  call apply_pbc_north(u)
  call apply_pbc_north(v)
  call apply_pbc_north(w)
  
  !apply pbc at south 2   !z
  call apply_pbc_south(u)
  call apply_pbc_south(v)
  call apply_pbc_south(w)
  
  return
  
 end subroutine apply_pbc_velocities_along_z
 
 subroutine apply_pbc_densities_along_xy
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the periodic boundary 
!     conditions to fluid densities along x and y
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  !sides 4 (6)
  
  !apply pbc at east 3    !x
  call apply_pbc_east(rhoR)
  
  !apply pbc at west 4    !x
  call apply_pbc_west(rhoR)
  
  !apply pbc at front 5   !y
  call apply_pbc_front(rhoR)
  
  !apply pbc at rear 6    !y
  call apply_pbc_rear(rhoR)
  
  !edges 4 (12)
  
  !apply pbc at front east 1   !xy
  call apply_pbc_edge_front_east(rhoR)
  
  !apply pbc at front west 2   !xy
  call apply_pbc_edge_front_west(rhoR)
  
  !apply pbc at rear east 7   !xy
  call apply_pbc_edge_rear_east(rhoR)
  
  !apply pbc at rear west 8   !xy
  call apply_pbc_edge_rear_west(rhoR)
  
  if(lsingle_fluid)return
  
  !sides 4 (6)
  
  !apply pbc at east 3    !x
  call apply_pbc_east(rhoB)
  
  !apply pbc at west 4    !x
  call apply_pbc_west(rhoB)
  
  !apply pbc at front 5   !y
  call apply_pbc_front(rhoB)
  
  !apply pbc at rear 6    !y
  call apply_pbc_rear(rhoB)
  
  !edges 4 (12)
  
  !apply pbc at front east 1   !xy
  call apply_pbc_edge_front_east(rhoB)
  
  !apply pbc at front west 2   !xy
  call apply_pbc_edge_front_west(rhoB)
  
  !apply pbc at rear east 7   !xy
  call apply_pbc_edge_rear_east(rhoB)
  
  !apply pbc at rear west 8   !xy
  call apply_pbc_edge_rear_west(rhoB)
  
  
  return
  
 end subroutine apply_pbc_densities_along_xy
 
 subroutine apply_pbc_velocities_along_xy
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the periodic boundary 
!     conditions to fluid velocities along x and y
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  !sides 4 (6)
  
  !apply pbc at east 3    !x
  call apply_pbc_east(u)
  call apply_pbc_east(v)
  call apply_pbc_east(w)
  
  !apply pbc at west 4    !x
  call apply_pbc_west(u)
  call apply_pbc_west(v)
  call apply_pbc_west(w)
  
  !apply pbc at front 5   !y
  call apply_pbc_front(u)
  call apply_pbc_front(v)
  call apply_pbc_front(w)
  
  !apply pbc at rear 6    !y
  call apply_pbc_rear(u)
  call apply_pbc_rear(v)
  call apply_pbc_rear(w)
  
  !edges 4 (12)
  
  !apply pbc at front east 1   !xy
  call apply_pbc_edge_front_east(u)
  call apply_pbc_edge_front_east(v)
  call apply_pbc_edge_front_east(w)
  
  !apply pbc at front west 2   !xy
  call apply_pbc_edge_front_west(u)
  call apply_pbc_edge_front_west(v)
  call apply_pbc_edge_front_west(w)
  
  !apply pbc at rear east 7   !xy
  call apply_pbc_edge_rear_east(u)
  call apply_pbc_edge_rear_east(v)
  call apply_pbc_edge_rear_east(w)
  
  !apply pbc at rear west 8   !xy
  call apply_pbc_edge_rear_west(u)
  call apply_pbc_edge_rear_west(v)
  call apply_pbc_edge_rear_west(w)
  
  
  return
  
 end subroutine apply_pbc_velocities_along_xy
 
 subroutine apply_pbc_densities_along_xz
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the full periodic boundary 
!     conditions to fluid densities along x and z
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  !sides 6 (4)
  
  !apply pbc at north 1   !z
  call apply_pbc_north(rhoR)
  
  !apply pbc at south 2   !z
  call apply_pbc_south(rhoR)
  
  !apply pbc at east 3    !x
  call apply_pbc_east(rhoR)
  
  !apply pbc at west 4    !x
  call apply_pbc_west(rhoR)
  
  !edges 4 (12)
  
  !apply pbc at north east 3   !xz
  call apply_pbc_edge_north_east(rhoR)
  
  !apply pbc at north west 6   !xz
  call apply_pbc_edge_north_west(rhoR)
  
  !apply pbc at south east 9  !xz
  call apply_pbc_edge_south_east(rhoR)
  
  !apply pbc at south west 12  !xz
  call apply_pbc_edge_south_west(rhoR)
  
  if(lsingle_fluid)return
  
  !sides 6 (4)
  
  !apply pbc at north 1   !z
  call apply_pbc_north(rhoB)
  
  !apply pbc at south 2   !z
  call apply_pbc_south(rhoB)
  
  !apply pbc at east 3    !x
  call apply_pbc_east(rhoB)
  
  !apply pbc at west 4    !x
  call apply_pbc_west(rhoB)
  
  !edges 4 (12)
  
  !apply pbc at north east 3   !xz
  call apply_pbc_edge_north_east(rhoB)
  
  !apply pbc at north west 6   !xz
  call apply_pbc_edge_north_west(rhoB)
  
  !apply pbc at south east 9  !xz
  call apply_pbc_edge_south_east(rhoB)
  
  !apply pbc at south west 12  !xz
  call apply_pbc_edge_south_west(rhoB)
  
  return
  
 end subroutine apply_pbc_densities_along_xz
 
 subroutine apply_pbc_velocities_along_xz
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the full periodic boundary 
!     conditions to fluid velocities along x and z
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  !sides 6 (4)
  
  !apply pbc at north 1   !z
  call apply_pbc_north(u)
  call apply_pbc_north(v)
  call apply_pbc_north(w)
  
  !apply pbc at south 2   !z
  call apply_pbc_south(u)
  call apply_pbc_south(v)
  call apply_pbc_south(w)
  
  !apply pbc at east 3    !x
  call apply_pbc_east(u)
  call apply_pbc_east(v)
  call apply_pbc_east(w)
  
  !apply pbc at west 4    !x
  call apply_pbc_west(u)
  call apply_pbc_west(v)
  call apply_pbc_west(w)
  
  !edges 4 (12)
  
  !apply pbc at north east 3   !xz
  call apply_pbc_edge_north_east(u)
  call apply_pbc_edge_north_east(v)
  call apply_pbc_edge_north_east(w)
  
  !apply pbc at north west 6   !xz
  call apply_pbc_edge_north_west(u)
  call apply_pbc_edge_north_west(v)
  call apply_pbc_edge_north_west(w)
  
  !apply pbc at south east 9  !xz
  call apply_pbc_edge_south_east(u)
  call apply_pbc_edge_south_east(v)
  call apply_pbc_edge_south_east(w)
  
  !apply pbc at south west 12  !xz
  call apply_pbc_edge_south_west(u)
  call apply_pbc_edge_south_west(v)
  call apply_pbc_edge_south_west(w)
  
  return
  
 end subroutine apply_pbc_velocities_along_xz
 
 subroutine apply_pbc_densities_along_yz
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the full periodic boundary 
!     conditions to fluid densities along y and z
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  !sides 4 (6)
  
  !apply pbc at north 1   !z
  call apply_pbc_north(rhoR)
  
  !apply pbc at south 2   !z
  call apply_pbc_south(rhoR)
  
  !apply pbc at front 5   !y
  call apply_pbc_front(rhoR)
  
  !apply pbc at rear 6    !y
  call apply_pbc_rear(rhoR)
  
  !edges 4 (12)
  
  !apply pbc at north front 4   !yz
  call apply_pbc_edge_north_front(rhoR)
  
  !apply pbc at north rear 5   !yz
  call apply_pbc_edge_north_rear(rhoR)
  
  !apply pbc at south front 10  !yz
  call apply_pbc_edge_south_front(rhoR)
  
  !apply pbc at south rear 11  !yz
  call apply_pbc_edge_south_rear(rhoR)
  
  if(lsingle_fluid)return
  
  !sides 4 (6)
  
  !apply pbc at north 1   !z
  call apply_pbc_north(rhoB)
  
  !apply pbc at south 2   !z
  call apply_pbc_south(rhoB)
  
  !apply pbc at front 5   !y
  call apply_pbc_front(rhoB)
  
  !apply pbc at rear 6    !y
  call apply_pbc_rear(rhoB)
  
  !edges 4 (12)
  
  !apply pbc at north front 4   !yz
  call apply_pbc_edge_north_front(rhoB)
  
  !apply pbc at north rear 5   !yz
  call apply_pbc_edge_north_rear(rhoB)
  
  !apply pbc at south front 10  !yz
  call apply_pbc_edge_south_front(rhoB)
  
  !apply pbc at south rear 11  !yz
  call apply_pbc_edge_south_rear(rhoB)
  
  return
  
 end subroutine apply_pbc_densities_along_yz
 
 subroutine apply_pbc_velocities_along_yz
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the full periodic boundary 
!     conditions to fluid velocities along y and z
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  !sides 4 (6)
  
  !apply pbc at north 1   !z
  call apply_pbc_north(u)
  call apply_pbc_north(v)
  call apply_pbc_north(w)
  
  !apply pbc at south 2   !z
  call apply_pbc_south(u)
  call apply_pbc_south(v)
  call apply_pbc_south(w)
  
  !apply pbc at front 5   !y
  call apply_pbc_front(u)
  call apply_pbc_front(v)
  call apply_pbc_front(w)
  
  !apply pbc at rear 6    !y
  call apply_pbc_rear(u)
  call apply_pbc_rear(v)
  call apply_pbc_rear(w)
  
  !edges 4 (12)
  
  !apply pbc at north front 4   !yz
  call apply_pbc_edge_north_front(u)
  call apply_pbc_edge_north_front(v)
  call apply_pbc_edge_north_front(w)
  
  !apply pbc at north rear 5   !yz
  call apply_pbc_edge_north_rear(u)
  call apply_pbc_edge_north_rear(v)
  call apply_pbc_edge_north_rear(w)
  
  !apply pbc at south front 10  !yz
  call apply_pbc_edge_south_front(u)
  call apply_pbc_edge_south_front(v)
  call apply_pbc_edge_south_front(w)
  
  !apply pbc at south rear 11  !yz
  call apply_pbc_edge_south_rear(u)
  call apply_pbc_edge_south_rear(v)
  call apply_pbc_edge_south_rear(w)
  
  return
  
 end subroutine apply_pbc_velocities_along_yz
 
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
  case(0) ! 0 0 0
    return
  case(1) ! 1 0 0 
    call driver_pbc_pops_along_x
  case(2) ! 0 1 0 
    call driver_pbc_pops_along_y
  case(3) ! 1 1 0 
    call driver_pbc_pops_along_xy
  case(4) ! 0 0 2 
    call driver_pbc_pops_along_z
  case(5) ! 1 0 1 
    call driver_pbc_pops_along_xz
  case(6) ! 0 1 1 
    call driver_pbc_pops_along_yz
  case(7) ! 1 1 1
    call driver_pbc_pops
  case default
    call error(12)
  end select
  
  return
  
 end subroutine driver_bc_pops
 
 subroutine driver_pbc_pops 
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the full periodic boundary 
!     conditions to fluid populations
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************

  implicit none
 
  call apply_pbc_pops(f00R,f01R,f02R,f03R,f04R,f05R,f06R,f07R,f08R, &
   f09R,f10R,f11R,f12R,f13R,f14R,f15R,f16R,f17R,f18R) 
  
  if(lsingle_fluid)return
 
  call apply_pbc_pops(f00B,f01B,f02B,f03B,f04B,f05B,f06B,f07B,f08B, &
   f09B,f10B,f11B,f12B,f13B,f14B,f15B,f16B,f17B,f18B)
   
  return
  
 end subroutine driver_pbc_pops
 
 subroutine apply_pbc_pops(f00sub,f01sub,f02sub,f03sub,f04sub,f05sub, &
  f06sub,f07sub,f08sub,f09sub,f10sub,f11sub,f12sub,f13sub,f14sub, &
  f15sub,f16sub,f17sub,f18sub) 
 
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
  
  real(kind=PRC), allocatable, dimension(:,:,:)  :: f00sub,f01sub, &
  f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub,f09sub,f10sub, &
  f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub
  
  !sides 6
  
  !apply pbc at north 1 !z
  !red fluid
  call apply_pbc_north(f00sub)
  call apply_pbc_north(f01sub)
  call apply_pbc_north(f02sub)
  call apply_pbc_north(f03sub)
  call apply_pbc_north(f04sub)
  call apply_pbc_north(f05sub)
  call apply_pbc_north(f06sub)
  call apply_pbc_north(f07sub)
  call apply_pbc_north(f08sub)
  call apply_pbc_north(f09sub)
  call apply_pbc_north(f10sub)
  call apply_pbc_north(f11sub)
  call apply_pbc_north(f12sub)
  call apply_pbc_north(f13sub)
  call apply_pbc_north(f14sub)
  call apply_pbc_north(f15sub)
  call apply_pbc_north(f16sub)
  call apply_pbc_north(f17sub)
  call apply_pbc_north(f18sub)
  
  !apply pbc at south 2 !z
  !red fluid
  call apply_pbc_south(f00sub)
  call apply_pbc_south(f01sub)
  call apply_pbc_south(f02sub)
  call apply_pbc_south(f03sub)
  call apply_pbc_south(f04sub)
  call apply_pbc_south(f05sub)
  call apply_pbc_south(f06sub)
  call apply_pbc_south(f07sub)
  call apply_pbc_south(f08sub)
  call apply_pbc_south(f09sub)
  call apply_pbc_south(f10sub)
  call apply_pbc_south(f11sub)
  call apply_pbc_south(f12sub)
  call apply_pbc_south(f13sub)
  call apply_pbc_south(f14sub)
  call apply_pbc_south(f15sub)
  call apply_pbc_south(f16sub)
  call apply_pbc_south(f17sub)
  call apply_pbc_south(f18sub)
  
  !apply pbc at east 3 !x
  !red fluid
  call apply_pbc_east(f00sub)
  call apply_pbc_east(f01sub)
  call apply_pbc_east(f02sub)
  call apply_pbc_east(f03sub)
  call apply_pbc_east(f04sub)
  call apply_pbc_east(f05sub)
  call apply_pbc_east(f06sub)
  call apply_pbc_east(f07sub)
  call apply_pbc_east(f08sub)
  call apply_pbc_east(f09sub)
  call apply_pbc_east(f10sub)
  call apply_pbc_east(f11sub)
  call apply_pbc_east(f12sub)
  call apply_pbc_east(f13sub)
  call apply_pbc_east(f14sub)
  call apply_pbc_east(f15sub)
  call apply_pbc_east(f16sub)
  call apply_pbc_east(f17sub)
  call apply_pbc_east(f18sub)
  
  !apply pbc at west 4 !x
  !red fluid
  call apply_pbc_west(f00sub)
  call apply_pbc_west(f01sub)
  call apply_pbc_west(f02sub)
  call apply_pbc_west(f03sub)
  call apply_pbc_west(f04sub)
  call apply_pbc_west(f05sub)
  call apply_pbc_west(f06sub)
  call apply_pbc_west(f07sub)
  call apply_pbc_west(f08sub)
  call apply_pbc_west(f09sub)
  call apply_pbc_west(f10sub)
  call apply_pbc_west(f11sub)
  call apply_pbc_west(f12sub)
  call apply_pbc_west(f13sub)
  call apply_pbc_west(f14sub)
  call apply_pbc_west(f15sub)
  call apply_pbc_west(f16sub)
  call apply_pbc_west(f17sub)
  call apply_pbc_west(f18sub)
  
  !apply pbc at front 5 !y
  !red fluid
  call apply_pbc_front(f00sub)
  call apply_pbc_front(f01sub)
  call apply_pbc_front(f02sub)
  call apply_pbc_front(f03sub)
  call apply_pbc_front(f04sub)
  call apply_pbc_front(f05sub)
  call apply_pbc_front(f06sub)
  call apply_pbc_front(f07sub)
  call apply_pbc_front(f08sub)
  call apply_pbc_front(f09sub)
  call apply_pbc_front(f10sub)
  call apply_pbc_front(f11sub)
  call apply_pbc_front(f12sub)
  call apply_pbc_front(f13sub)
  call apply_pbc_front(f14sub)
  call apply_pbc_front(f15sub)
  call apply_pbc_front(f16sub)
  call apply_pbc_front(f17sub)
  call apply_pbc_front(f18sub)
  
  !apply pbc at rear 6 !y
  !red fluid
  call apply_pbc_rear(f00sub)
  call apply_pbc_rear(f01sub)
  call apply_pbc_rear(f02sub)
  call apply_pbc_rear(f03sub)
  call apply_pbc_rear(f04sub)
  call apply_pbc_rear(f05sub)
  call apply_pbc_rear(f06sub)
  call apply_pbc_rear(f07sub)
  call apply_pbc_rear(f08sub)
  call apply_pbc_rear(f09sub)
  call apply_pbc_rear(f10sub)
  call apply_pbc_rear(f11sub)
  call apply_pbc_rear(f12sub)
  call apply_pbc_rear(f13sub)
  call apply_pbc_rear(f14sub)
  call apply_pbc_rear(f15sub)
  call apply_pbc_rear(f16sub)
  call apply_pbc_rear(f17sub)
  call apply_pbc_rear(f18sub)
  
  !edges 12
  
  !apply pbc at front east 1 !xy
  !red fluid
  call apply_pbc_edge_front_east(f00sub)
  call apply_pbc_edge_front_east(f01sub)
  call apply_pbc_edge_front_east(f02sub)
  call apply_pbc_edge_front_east(f03sub)
  call apply_pbc_edge_front_east(f04sub)
  call apply_pbc_edge_front_east(f05sub)
  call apply_pbc_edge_front_east(f06sub)
  call apply_pbc_edge_front_east(f07sub)
  call apply_pbc_edge_front_east(f08sub)
  call apply_pbc_edge_front_east(f09sub)
  call apply_pbc_edge_front_east(f10sub)
  call apply_pbc_edge_front_east(f11sub)
  call apply_pbc_edge_front_east(f12sub)
  call apply_pbc_edge_front_east(f13sub)
  call apply_pbc_edge_front_east(f14sub)
  call apply_pbc_edge_front_east(f15sub)
  call apply_pbc_edge_front_east(f16sub)
  call apply_pbc_edge_front_east(f17sub)
  call apply_pbc_edge_front_east(f18sub)
  
  !apply pbc at front west 2 !xy
  !red fluid
  call apply_pbc_edge_front_west(f00sub)
  call apply_pbc_edge_front_west(f01sub)
  call apply_pbc_edge_front_west(f02sub)
  call apply_pbc_edge_front_west(f03sub)
  call apply_pbc_edge_front_west(f04sub)
  call apply_pbc_edge_front_west(f05sub)
  call apply_pbc_edge_front_west(f06sub)
  call apply_pbc_edge_front_west(f07sub)
  call apply_pbc_edge_front_west(f08sub)
  call apply_pbc_edge_front_west(f09sub)
  call apply_pbc_edge_front_west(f10sub)
  call apply_pbc_edge_front_west(f11sub)
  call apply_pbc_edge_front_west(f12sub)
  call apply_pbc_edge_front_west(f13sub)
  call apply_pbc_edge_front_west(f14sub)
  call apply_pbc_edge_front_west(f15sub)
  call apply_pbc_edge_front_west(f16sub)
  call apply_pbc_edge_front_west(f17sub)
  call apply_pbc_edge_front_west(f18sub)
  
  !apply pbc at north east 3 !xz
  !red fluid
  call apply_pbc_edge_north_east(f00sub)
  call apply_pbc_edge_north_east(f01sub)
  call apply_pbc_edge_north_east(f02sub)
  call apply_pbc_edge_north_east(f03sub)
  call apply_pbc_edge_north_east(f04sub)
  call apply_pbc_edge_north_east(f05sub)
  call apply_pbc_edge_north_east(f06sub)
  call apply_pbc_edge_north_east(f07sub)
  call apply_pbc_edge_north_east(f08sub)
  call apply_pbc_edge_north_east(f09sub)
  call apply_pbc_edge_north_east(f10sub)
  call apply_pbc_edge_north_east(f11sub)
  call apply_pbc_edge_north_east(f12sub)
  call apply_pbc_edge_north_east(f13sub)
  call apply_pbc_edge_north_east(f14sub)
  call apply_pbc_edge_north_east(f15sub)
  call apply_pbc_edge_north_east(f16sub)
  call apply_pbc_edge_north_east(f17sub)
  call apply_pbc_edge_north_east(f18sub)
  
  !apply pbc at north front 4 !yz
  !red fluid
  call apply_pbc_edge_north_front(f00sub)
  call apply_pbc_edge_north_front(f01sub)
  call apply_pbc_edge_north_front(f02sub)
  call apply_pbc_edge_north_front(f03sub)
  call apply_pbc_edge_north_front(f04sub)
  call apply_pbc_edge_north_front(f05sub)
  call apply_pbc_edge_north_front(f06sub)
  call apply_pbc_edge_north_front(f07sub)
  call apply_pbc_edge_north_front(f08sub)
  call apply_pbc_edge_north_front(f09sub)
  call apply_pbc_edge_north_front(f10sub)
  call apply_pbc_edge_north_front(f11sub)
  call apply_pbc_edge_north_front(f12sub)
  call apply_pbc_edge_north_front(f13sub)
  call apply_pbc_edge_north_front(f14sub)
  call apply_pbc_edge_north_front(f15sub)
  call apply_pbc_edge_north_front(f16sub)
  call apply_pbc_edge_north_front(f17sub)
  call apply_pbc_edge_north_front(f18sub)
  
  !apply pbc at north rear 5 !yz
  !red fluid
  call apply_pbc_edge_north_rear(f00sub)
  call apply_pbc_edge_north_rear(f01sub)
  call apply_pbc_edge_north_rear(f02sub)
  call apply_pbc_edge_north_rear(f03sub)
  call apply_pbc_edge_north_rear(f04sub)
  call apply_pbc_edge_north_rear(f05sub)
  call apply_pbc_edge_north_rear(f06sub)
  call apply_pbc_edge_north_rear(f07sub)
  call apply_pbc_edge_north_rear(f08sub)
  call apply_pbc_edge_north_rear(f09sub)
  call apply_pbc_edge_north_rear(f10sub)
  call apply_pbc_edge_north_rear(f11sub)
  call apply_pbc_edge_north_rear(f12sub)
  call apply_pbc_edge_north_rear(f13sub)
  call apply_pbc_edge_north_rear(f14sub)
  call apply_pbc_edge_north_rear(f15sub)
  call apply_pbc_edge_north_rear(f16sub)
  call apply_pbc_edge_north_rear(f17sub)
  call apply_pbc_edge_north_rear(f18sub)
  
  !apply pbc at north west 6 !xz
  !red fluid
  call apply_pbc_edge_north_west(f00sub)
  call apply_pbc_edge_north_west(f01sub)
  call apply_pbc_edge_north_west(f02sub)
  call apply_pbc_edge_north_west(f03sub)
  call apply_pbc_edge_north_west(f04sub)
  call apply_pbc_edge_north_west(f05sub)
  call apply_pbc_edge_north_west(f06sub)
  call apply_pbc_edge_north_west(f07sub)
  call apply_pbc_edge_north_west(f08sub)
  call apply_pbc_edge_north_west(f09sub)
  call apply_pbc_edge_north_west(f10sub)
  call apply_pbc_edge_north_west(f11sub)
  call apply_pbc_edge_north_west(f12sub)
  call apply_pbc_edge_north_west(f13sub)
  call apply_pbc_edge_north_west(f14sub)
  call apply_pbc_edge_north_west(f15sub)
  call apply_pbc_edge_north_west(f16sub)
  call apply_pbc_edge_north_west(f17sub)
  call apply_pbc_edge_north_west(f18sub)
  
  !apply pbc at rear east 7 !xy
  !red fluid
  call apply_pbc_edge_rear_east(f00sub)
  call apply_pbc_edge_rear_east(f01sub)
  call apply_pbc_edge_rear_east(f02sub)
  call apply_pbc_edge_rear_east(f03sub)
  call apply_pbc_edge_rear_east(f04sub)
  call apply_pbc_edge_rear_east(f05sub)
  call apply_pbc_edge_rear_east(f06sub)
  call apply_pbc_edge_rear_east(f07sub)
  call apply_pbc_edge_rear_east(f08sub)
  call apply_pbc_edge_rear_east(f09sub)
  call apply_pbc_edge_rear_east(f10sub)
  call apply_pbc_edge_rear_east(f11sub)
  call apply_pbc_edge_rear_east(f12sub)
  call apply_pbc_edge_rear_east(f13sub)
  call apply_pbc_edge_rear_east(f14sub)
  call apply_pbc_edge_rear_east(f15sub)
  call apply_pbc_edge_rear_east(f16sub)
  call apply_pbc_edge_rear_east(f17sub)
  call apply_pbc_edge_rear_east(f18sub)
  
  !apply pbc at rear west 8 !xy
  !red fluid
  call apply_pbc_edge_rear_west(f00sub)
  call apply_pbc_edge_rear_west(f01sub)
  call apply_pbc_edge_rear_west(f02sub)
  call apply_pbc_edge_rear_west(f03sub)
  call apply_pbc_edge_rear_west(f04sub)
  call apply_pbc_edge_rear_west(f05sub)
  call apply_pbc_edge_rear_west(f06sub)
  call apply_pbc_edge_rear_west(f07sub)
  call apply_pbc_edge_rear_west(f08sub)
  call apply_pbc_edge_rear_west(f09sub)
  call apply_pbc_edge_rear_west(f10sub)
  call apply_pbc_edge_rear_west(f11sub)
  call apply_pbc_edge_rear_west(f12sub)
  call apply_pbc_edge_rear_west(f13sub)
  call apply_pbc_edge_rear_west(f14sub)
  call apply_pbc_edge_rear_west(f15sub)
  call apply_pbc_edge_rear_west(f16sub)
  call apply_pbc_edge_rear_west(f17sub)
  call apply_pbc_edge_rear_west(f18sub)
  
  !apply pbc at south east 9 !xz
  !red fluid
  call apply_pbc_edge_south_east(f00sub)
  call apply_pbc_edge_south_east(f01sub)
  call apply_pbc_edge_south_east(f02sub)
  call apply_pbc_edge_south_east(f03sub)
  call apply_pbc_edge_south_east(f04sub)
  call apply_pbc_edge_south_east(f05sub)
  call apply_pbc_edge_south_east(f06sub)
  call apply_pbc_edge_south_east(f07sub)
  call apply_pbc_edge_south_east(f08sub)
  call apply_pbc_edge_south_east(f09sub)
  call apply_pbc_edge_south_east(f10sub)
  call apply_pbc_edge_south_east(f11sub)
  call apply_pbc_edge_south_east(f12sub)
  call apply_pbc_edge_south_east(f13sub)
  call apply_pbc_edge_south_east(f14sub)
  call apply_pbc_edge_south_east(f15sub)
  call apply_pbc_edge_south_east(f16sub)
  call apply_pbc_edge_south_east(f17sub)
  call apply_pbc_edge_south_east(f18sub)
  
  !apply pbc at south front 10 !yz
  !red fluid
  call apply_pbc_edge_south_front(f00sub)
  call apply_pbc_edge_south_front(f01sub)
  call apply_pbc_edge_south_front(f02sub)
  call apply_pbc_edge_south_front(f03sub)
  call apply_pbc_edge_south_front(f04sub)
  call apply_pbc_edge_south_front(f05sub)
  call apply_pbc_edge_south_front(f06sub)
  call apply_pbc_edge_south_front(f07sub)
  call apply_pbc_edge_south_front(f08sub)
  call apply_pbc_edge_south_front(f09sub)
  call apply_pbc_edge_south_front(f10sub)
  call apply_pbc_edge_south_front(f11sub)
  call apply_pbc_edge_south_front(f12sub)
  call apply_pbc_edge_south_front(f13sub)
  call apply_pbc_edge_south_front(f14sub)
  call apply_pbc_edge_south_front(f15sub)
  call apply_pbc_edge_south_front(f16sub)
  call apply_pbc_edge_south_front(f17sub)
  call apply_pbc_edge_south_front(f18sub)
  
  !apply pbc at south rear 11 !yz
  !red fluid
  call apply_pbc_edge_south_rear(f00sub)
  call apply_pbc_edge_south_rear(f01sub)
  call apply_pbc_edge_south_rear(f02sub)
  call apply_pbc_edge_south_rear(f03sub)
  call apply_pbc_edge_south_rear(f04sub)
  call apply_pbc_edge_south_rear(f05sub)
  call apply_pbc_edge_south_rear(f06sub)
  call apply_pbc_edge_south_rear(f07sub)
  call apply_pbc_edge_south_rear(f08sub)
  call apply_pbc_edge_south_rear(f09sub)
  call apply_pbc_edge_south_rear(f10sub)
  call apply_pbc_edge_south_rear(f11sub)
  call apply_pbc_edge_south_rear(f12sub)
  call apply_pbc_edge_south_rear(f13sub)
  call apply_pbc_edge_south_rear(f14sub)
  call apply_pbc_edge_south_rear(f15sub)
  call apply_pbc_edge_south_rear(f16sub)
  call apply_pbc_edge_south_rear(f17sub)
  call apply_pbc_edge_south_rear(f18sub)
  
  !apply pbc at south west 12 !xz
  !red fluid
  call apply_pbc_edge_south_west(f00sub)
  call apply_pbc_edge_south_west(f01sub)
  call apply_pbc_edge_south_west(f02sub)
  call apply_pbc_edge_south_west(f03sub)
  call apply_pbc_edge_south_west(f04sub)
  call apply_pbc_edge_south_west(f05sub)
  call apply_pbc_edge_south_west(f06sub)
  call apply_pbc_edge_south_west(f07sub)
  call apply_pbc_edge_south_west(f08sub)
  call apply_pbc_edge_south_west(f09sub)
  call apply_pbc_edge_south_west(f10sub)
  call apply_pbc_edge_south_west(f11sub)
  call apply_pbc_edge_south_west(f12sub)
  call apply_pbc_edge_south_west(f13sub)
  call apply_pbc_edge_south_west(f14sub)
  call apply_pbc_edge_south_west(f15sub)
  call apply_pbc_edge_south_west(f16sub)
  call apply_pbc_edge_south_west(f17sub)
  call apply_pbc_edge_south_west(f18sub)
  
  !corner 8
  
  !apply pbc at north east front 1
  !red fluid
  call apply_pbc_corner_north_east_front(f00sub)
  call apply_pbc_corner_north_east_front(f01sub)
  call apply_pbc_corner_north_east_front(f02sub)
  call apply_pbc_corner_north_east_front(f03sub)
  call apply_pbc_corner_north_east_front(f04sub)
  call apply_pbc_corner_north_east_front(f05sub)
  call apply_pbc_corner_north_east_front(f06sub)
  call apply_pbc_corner_north_east_front(f07sub)
  call apply_pbc_corner_north_east_front(f08sub)
  call apply_pbc_corner_north_east_front(f09sub)
  call apply_pbc_corner_north_east_front(f10sub)
  call apply_pbc_corner_north_east_front(f11sub)
  call apply_pbc_corner_north_east_front(f12sub)
  call apply_pbc_corner_north_east_front(f13sub)
  call apply_pbc_corner_north_east_front(f14sub)
  call apply_pbc_corner_north_east_front(f15sub)
  call apply_pbc_corner_north_east_front(f16sub)
  call apply_pbc_corner_north_east_front(f17sub)
  call apply_pbc_corner_north_east_front(f18sub)
  
  !apply pbc at north east rear 2
  !red fluid
  call apply_pbc_corner_north_east_rear(f00sub)
  call apply_pbc_corner_north_east_rear(f01sub)
  call apply_pbc_corner_north_east_rear(f02sub)
  call apply_pbc_corner_north_east_rear(f03sub)
  call apply_pbc_corner_north_east_rear(f04sub)
  call apply_pbc_corner_north_east_rear(f05sub)
  call apply_pbc_corner_north_east_rear(f06sub)
  call apply_pbc_corner_north_east_rear(f07sub)
  call apply_pbc_corner_north_east_rear(f08sub)
  call apply_pbc_corner_north_east_rear(f09sub)
  call apply_pbc_corner_north_east_rear(f10sub)
  call apply_pbc_corner_north_east_rear(f11sub)
  call apply_pbc_corner_north_east_rear(f12sub)
  call apply_pbc_corner_north_east_rear(f13sub)
  call apply_pbc_corner_north_east_rear(f14sub)
  call apply_pbc_corner_north_east_rear(f15sub)
  call apply_pbc_corner_north_east_rear(f16sub)
  call apply_pbc_corner_north_east_rear(f17sub)
  call apply_pbc_corner_north_east_rear(f18sub)
  
  !apply pbc at north west rear 3
  !red fluid
  call apply_pbc_corner_north_west_rear(f00sub)
  call apply_pbc_corner_north_west_rear(f01sub)
  call apply_pbc_corner_north_west_rear(f02sub)
  call apply_pbc_corner_north_west_rear(f03sub)
  call apply_pbc_corner_north_west_rear(f04sub)
  call apply_pbc_corner_north_west_rear(f05sub)
  call apply_pbc_corner_north_west_rear(f06sub)
  call apply_pbc_corner_north_west_rear(f07sub)
  call apply_pbc_corner_north_west_rear(f08sub)
  call apply_pbc_corner_north_west_rear(f09sub)
  call apply_pbc_corner_north_west_rear(f10sub)
  call apply_pbc_corner_north_west_rear(f11sub)
  call apply_pbc_corner_north_west_rear(f12sub)
  call apply_pbc_corner_north_west_rear(f13sub)
  call apply_pbc_corner_north_west_rear(f14sub)
  call apply_pbc_corner_north_west_rear(f15sub)
  call apply_pbc_corner_north_west_rear(f16sub)
  call apply_pbc_corner_north_west_rear(f17sub)
  call apply_pbc_corner_north_west_rear(f18sub)
  
  !apply pbc at north west front 4
  !red fluid
  call apply_pbc_corner_north_west_front(f00sub)
  call apply_pbc_corner_north_west_front(f01sub)
  call apply_pbc_corner_north_west_front(f02sub)
  call apply_pbc_corner_north_west_front(f03sub)
  call apply_pbc_corner_north_west_front(f04sub)
  call apply_pbc_corner_north_west_front(f05sub)
  call apply_pbc_corner_north_west_front(f06sub)
  call apply_pbc_corner_north_west_front(f07sub)
  call apply_pbc_corner_north_west_front(f08sub)
  call apply_pbc_corner_north_west_front(f09sub)
  call apply_pbc_corner_north_west_front(f10sub)
  call apply_pbc_corner_north_west_front(f11sub)
  call apply_pbc_corner_north_west_front(f12sub)
  call apply_pbc_corner_north_west_front(f13sub)
  call apply_pbc_corner_north_west_front(f14sub)
  call apply_pbc_corner_north_west_front(f15sub)
  call apply_pbc_corner_north_west_front(f16sub)
  call apply_pbc_corner_north_west_front(f17sub)
  call apply_pbc_corner_north_west_front(f18sub)
  
  !apply pbc at south east front 5
  !red fluid
  call apply_pbc_corner_south_east_front(f00sub)
  call apply_pbc_corner_south_east_front(f01sub)
  call apply_pbc_corner_south_east_front(f02sub)
  call apply_pbc_corner_south_east_front(f03sub)
  call apply_pbc_corner_south_east_front(f04sub)
  call apply_pbc_corner_south_east_front(f05sub)
  call apply_pbc_corner_south_east_front(f06sub)
  call apply_pbc_corner_south_east_front(f07sub)
  call apply_pbc_corner_south_east_front(f08sub)
  call apply_pbc_corner_south_east_front(f09sub)
  call apply_pbc_corner_south_east_front(f10sub)
  call apply_pbc_corner_south_east_front(f11sub)
  call apply_pbc_corner_south_east_front(f12sub)
  call apply_pbc_corner_south_east_front(f13sub)
  call apply_pbc_corner_south_east_front(f14sub)
  call apply_pbc_corner_south_east_front(f15sub)
  call apply_pbc_corner_south_east_front(f16sub)
  call apply_pbc_corner_south_east_front(f17sub)
  call apply_pbc_corner_south_east_front(f18sub)
  
  !apply pbc at south west front 6
  !red fluid
  call apply_pbc_corner_south_west_front(f00sub)
  call apply_pbc_corner_south_west_front(f01sub)
  call apply_pbc_corner_south_west_front(f02sub)
  call apply_pbc_corner_south_west_front(f03sub)
  call apply_pbc_corner_south_west_front(f04sub)
  call apply_pbc_corner_south_west_front(f05sub)
  call apply_pbc_corner_south_west_front(f06sub)
  call apply_pbc_corner_south_west_front(f07sub)
  call apply_pbc_corner_south_west_front(f08sub)
  call apply_pbc_corner_south_west_front(f09sub)
  call apply_pbc_corner_south_west_front(f10sub)
  call apply_pbc_corner_south_west_front(f11sub)
  call apply_pbc_corner_south_west_front(f12sub)
  call apply_pbc_corner_south_west_front(f13sub)
  call apply_pbc_corner_south_west_front(f14sub)
  call apply_pbc_corner_south_west_front(f15sub)
  call apply_pbc_corner_south_west_front(f16sub)
  call apply_pbc_corner_south_west_front(f17sub)
  call apply_pbc_corner_south_west_front(f18sub)
  
  !apply pbc at south west rear 7
  !red fluid
  call apply_pbc_corner_south_west_rear(f00sub)
  call apply_pbc_corner_south_west_rear(f01sub)
  call apply_pbc_corner_south_west_rear(f02sub)
  call apply_pbc_corner_south_west_rear(f03sub)
  call apply_pbc_corner_south_west_rear(f04sub)
  call apply_pbc_corner_south_west_rear(f05sub)
  call apply_pbc_corner_south_west_rear(f06sub)
  call apply_pbc_corner_south_west_rear(f07sub)
  call apply_pbc_corner_south_west_rear(f08sub)
  call apply_pbc_corner_south_west_rear(f09sub)
  call apply_pbc_corner_south_west_rear(f10sub)
  call apply_pbc_corner_south_west_rear(f11sub)
  call apply_pbc_corner_south_west_rear(f12sub)
  call apply_pbc_corner_south_west_rear(f13sub)
  call apply_pbc_corner_south_west_rear(f14sub)
  call apply_pbc_corner_south_west_rear(f15sub)
  call apply_pbc_corner_south_west_rear(f16sub)
  call apply_pbc_corner_south_west_rear(f17sub)
  call apply_pbc_corner_south_west_rear(f18sub)
  
  !apply pbc at south east rear 8
  !red fluid
  call apply_pbc_corner_south_east_rear(f00sub)
  call apply_pbc_corner_south_east_rear(f01sub)
  call apply_pbc_corner_south_east_rear(f02sub)
  call apply_pbc_corner_south_east_rear(f03sub)
  call apply_pbc_corner_south_east_rear(f04sub)
  call apply_pbc_corner_south_east_rear(f05sub)
  call apply_pbc_corner_south_east_rear(f06sub)
  call apply_pbc_corner_south_east_rear(f07sub)
  call apply_pbc_corner_south_east_rear(f08sub)
  call apply_pbc_corner_south_east_rear(f09sub)
  call apply_pbc_corner_south_east_rear(f10sub)
  call apply_pbc_corner_south_east_rear(f11sub)
  call apply_pbc_corner_south_east_rear(f12sub)
  call apply_pbc_corner_south_east_rear(f13sub)
  call apply_pbc_corner_south_east_rear(f14sub)
  call apply_pbc_corner_south_east_rear(f15sub)
  call apply_pbc_corner_south_east_rear(f16sub)
  call apply_pbc_corner_south_east_rear(f17sub)
  call apply_pbc_corner_south_east_rear(f18sub)
  
  return
  
 end subroutine apply_pbc_pops
 
 subroutine driver_pbc_pops_along_x
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the periodic boundary 
!     conditions to fluid populations along x
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
 
  call apply_pbc_pops_along_x(f00R,f01R,f02R,f03R,f04R,f05R,f06R,f07R, &
   f08R,f09R,f10R,f11R,f12R,f13R,f14R,f15R,f16R,f17R,f18R) 
  
  if(lsingle_fluid)return
 
  call apply_pbc_pops_along_x(f00B,f01B,f02B,f03B,f04B,f05B,f06B,f07B, &
   f08B,f09B,f10B,f11B,f12B,f13B,f14B,f15B,f16B,f17B,f18B)
   
  return
  
 end subroutine driver_pbc_pops_along_x
 
 subroutine apply_pbc_pops_along_x(f00sub,f01sub,f02sub,f03sub,f04sub, &
  f05sub,f06sub,f07sub,f08sub,f09sub,f10sub,f11sub,f12sub,f13sub, &
  f14sub,f15sub,f16sub,f17sub,f18sub)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the periodic boundary 
!     conditions to fluid populations along x
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), allocatable, dimension(:,:,:)  :: f00sub,f01sub, &
  f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub,f09sub,f10sub, &
  f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub
  
  !sides 2 (6)
  
  !apply pbc at east 3 !x
  !red fluid
  call apply_pbc_east(f00sub)
  call apply_pbc_east(f01sub)
  call apply_pbc_east(f02sub)
  call apply_pbc_east(f03sub)
  call apply_pbc_east(f04sub)
  call apply_pbc_east(f05sub)
  call apply_pbc_east(f06sub)
  call apply_pbc_east(f07sub)
  call apply_pbc_east(f08sub)
  call apply_pbc_east(f09sub)
  call apply_pbc_east(f10sub)
  call apply_pbc_east(f11sub)
  call apply_pbc_east(f12sub)
  call apply_pbc_east(f13sub)
  call apply_pbc_east(f14sub)
  call apply_pbc_east(f15sub)
  call apply_pbc_east(f16sub)
  call apply_pbc_east(f17sub)
  call apply_pbc_east(f18sub)
  
  !apply pbc at west 4 !x
  !red fluid
  call apply_pbc_west(f00sub)
  call apply_pbc_west(f01sub)
  call apply_pbc_west(f02sub)
  call apply_pbc_west(f03sub)
  call apply_pbc_west(f04sub)
  call apply_pbc_west(f05sub)
  call apply_pbc_west(f06sub)
  call apply_pbc_west(f07sub)
  call apply_pbc_west(f08sub)
  call apply_pbc_west(f09sub)
  call apply_pbc_west(f10sub)
  call apply_pbc_west(f11sub)
  call apply_pbc_west(f12sub)
  call apply_pbc_west(f13sub)
  call apply_pbc_west(f14sub)
  call apply_pbc_west(f15sub)
  call apply_pbc_west(f16sub)
  call apply_pbc_west(f17sub)
  call apply_pbc_west(f18sub)
  
  return
  
 end subroutine apply_pbc_pops_along_x
 
 subroutine driver_pbc_pops_along_y
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the periodic boundary 
!     conditions to fluid populations along y
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
 
  call apply_pbc_pops_along_y(f00R,f01R,f02R,f03R,f04R,f05R,f06R,f07R, &
   f08R,f09R,f10R,f11R,f12R,f13R,f14R,f15R,f16R,f17R,f18R) 
  
  if(lsingle_fluid)return
 
  call apply_pbc_pops_along_y(f00B,f01B,f02B,f03B,f04B,f05B,f06B,f07B, &
   f08B,f09B,f10B,f11B,f12B,f13B,f14B,f15B,f16B,f17B,f18B)
   
  return
  
 end subroutine driver_pbc_pops_along_y
 
 subroutine apply_pbc_pops_along_y(f00sub,f01sub,f02sub,f03sub,f04sub, &
  f05sub,f06sub,f07sub,f08sub,f09sub,f10sub,f11sub,f12sub,f13sub, &
  f14sub,f15sub,f16sub,f17sub,f18sub)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the periodic boundary 
!     conditions to fluid populations along y
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), allocatable, dimension(:,:,:)  :: f00sub,f01sub, &
  f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub,f09sub,f10sub, &
  f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub
  
  !sides 2 (6)
  
  !apply pbc at front 5 !y
  !red fluid
  call apply_pbc_front(f00sub)
  call apply_pbc_front(f01sub)
  call apply_pbc_front(f02sub)
  call apply_pbc_front(f03sub)
  call apply_pbc_front(f04sub)
  call apply_pbc_front(f05sub)
  call apply_pbc_front(f06sub)
  call apply_pbc_front(f07sub)
  call apply_pbc_front(f08sub)
  call apply_pbc_front(f09sub)
  call apply_pbc_front(f10sub)
  call apply_pbc_front(f11sub)
  call apply_pbc_front(f12sub)
  call apply_pbc_front(f13sub)
  call apply_pbc_front(f14sub)
  call apply_pbc_front(f15sub)
  call apply_pbc_front(f16sub)
  call apply_pbc_front(f17sub)
  call apply_pbc_front(f18sub)
  
  !apply pbc at rear 6 !y
  !red fluid
  call apply_pbc_rear(f00sub)
  call apply_pbc_rear(f01sub)
  call apply_pbc_rear(f02sub)
  call apply_pbc_rear(f03sub)
  call apply_pbc_rear(f04sub)
  call apply_pbc_rear(f05sub)
  call apply_pbc_rear(f06sub)
  call apply_pbc_rear(f07sub)
  call apply_pbc_rear(f08sub)
  call apply_pbc_rear(f09sub)
  call apply_pbc_rear(f10sub)
  call apply_pbc_rear(f11sub)
  call apply_pbc_rear(f12sub)
  call apply_pbc_rear(f13sub)
  call apply_pbc_rear(f14sub)
  call apply_pbc_rear(f15sub)
  call apply_pbc_rear(f16sub)
  call apply_pbc_rear(f17sub)
  call apply_pbc_rear(f18sub)
  
  return
  
 end subroutine apply_pbc_pops_along_y
 
 subroutine driver_pbc_pops_along_z
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the periodic boundary 
!     conditions to fluid populations along z
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
 
  call apply_pbc_pops_along_z(f00R,f01R,f02R,f03R,f04R,f05R,f06R,f07R, &
   f08R,f09R,f10R,f11R,f12R,f13R,f14R,f15R,f16R,f17R,f18R) 
  
  if(lsingle_fluid)return
 
  call apply_pbc_pops_along_z(f00B,f01B,f02B,f03B,f04B,f05B,f06B,f07B, &
   f08B,f09B,f10B,f11B,f12B,f13B,f14B,f15B,f16B,f17B,f18B)
   
  return
  
 end subroutine driver_pbc_pops_along_z
 
 subroutine apply_pbc_pops_along_z(f00sub,f01sub,f02sub,f03sub,f04sub, &
  f05sub,f06sub,f07sub,f08sub,f09sub,f10sub,f11sub,f12sub,f13sub, &
  f14sub,f15sub,f16sub,f17sub,f18sub)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the periodic boundary 
!     conditions to fluid populations along z
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), allocatable, dimension(:,:,:)  :: f00sub,f01sub, &
  f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub,f09sub,f10sub, &
  f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub
  
  !sides 2 (6)
  
  !apply pbc at north 1 !z
  !red fluid
  call apply_pbc_north(f00sub)
  call apply_pbc_north(f01sub)
  call apply_pbc_north(f02sub)
  call apply_pbc_north(f03sub)
  call apply_pbc_north(f04sub)
  call apply_pbc_north(f05sub)
  call apply_pbc_north(f06sub)
  call apply_pbc_north(f07sub)
  call apply_pbc_north(f08sub)
  call apply_pbc_north(f09sub)
  call apply_pbc_north(f10sub)
  call apply_pbc_north(f11sub)
  call apply_pbc_north(f12sub)
  call apply_pbc_north(f13sub)
  call apply_pbc_north(f14sub)
  call apply_pbc_north(f15sub)
  call apply_pbc_north(f16sub)
  call apply_pbc_north(f17sub)
  call apply_pbc_north(f18sub)
  
  !apply pbc at south 2 !z
  !red fluid
  call apply_pbc_south(f00sub)
  call apply_pbc_south(f01sub)
  call apply_pbc_south(f02sub)
  call apply_pbc_south(f03sub)
  call apply_pbc_south(f04sub)
  call apply_pbc_south(f05sub)
  call apply_pbc_south(f06sub)
  call apply_pbc_south(f07sub)
  call apply_pbc_south(f08sub)
  call apply_pbc_south(f09sub)
  call apply_pbc_south(f10sub)
  call apply_pbc_south(f11sub)
  call apply_pbc_south(f12sub)
  call apply_pbc_south(f13sub)
  call apply_pbc_south(f14sub)
  call apply_pbc_south(f15sub)
  call apply_pbc_south(f16sub)
  call apply_pbc_south(f17sub)
  call apply_pbc_south(f18sub)
  
  return
  
 end subroutine apply_pbc_pops_along_z
 
 subroutine driver_pbc_pops_along_xy
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the periodic boundary 
!     conditions to fluid populations along x and y
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
 
  call apply_pbc_pops_along_xy(f00R,f01R,f02R,f03R,f04R,f05R,f06R,f07R,&
   f08R,f09R,f10R,f11R,f12R,f13R,f14R,f15R,f16R,f17R,f18R) 
  
  if(lsingle_fluid)return
 
  call apply_pbc_pops_along_xy(f00B,f01B,f02B,f03B,f04B,f05B,f06B,f07B,&
   f08B,f09B,f10B,f11B,f12B,f13B,f14B,f15B,f16B,f17B,f18B)
   
  return
  
 end subroutine driver_pbc_pops_along_xy
 
 subroutine apply_pbc_pops_along_xy(f00sub,f01sub,f02sub,f03sub,f04sub,&
  f05sub,f06sub,f07sub,f08sub,f09sub,f10sub,f11sub,f12sub,f13sub, &
  f14sub,f15sub,f16sub,f17sub,f18sub)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the periodic boundary 
!     conditions to fluid populations along x and y
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), allocatable, dimension(:,:,:)  :: f00sub,f01sub, &
  f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub,f09sub,f10sub, &
  f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub
  
  !sides 4 (6)
  
  !apply pbc at east 3 !x
  !red fluid
  call apply_pbc_east(f00sub)
  call apply_pbc_east(f01sub)
  call apply_pbc_east(f02sub)
  call apply_pbc_east(f03sub)
  call apply_pbc_east(f04sub)
  call apply_pbc_east(f05sub)
  call apply_pbc_east(f06sub)
  call apply_pbc_east(f07sub)
  call apply_pbc_east(f08sub)
  call apply_pbc_east(f09sub)
  call apply_pbc_east(f10sub)
  call apply_pbc_east(f11sub)
  call apply_pbc_east(f12sub)
  call apply_pbc_east(f13sub)
  call apply_pbc_east(f14sub)
  call apply_pbc_east(f15sub)
  call apply_pbc_east(f16sub)
  call apply_pbc_east(f17sub)
  call apply_pbc_east(f18sub)
  
  !apply pbc at west 4 !x
  !red fluid
  call apply_pbc_west(f00sub)
  call apply_pbc_west(f01sub)
  call apply_pbc_west(f02sub)
  call apply_pbc_west(f03sub)
  call apply_pbc_west(f04sub)
  call apply_pbc_west(f05sub)
  call apply_pbc_west(f06sub)
  call apply_pbc_west(f07sub)
  call apply_pbc_west(f08sub)
  call apply_pbc_west(f09sub)
  call apply_pbc_west(f10sub)
  call apply_pbc_west(f11sub)
  call apply_pbc_west(f12sub)
  call apply_pbc_west(f13sub)
  call apply_pbc_west(f14sub)
  call apply_pbc_west(f15sub)
  call apply_pbc_west(f16sub)
  call apply_pbc_west(f17sub)
  call apply_pbc_west(f18sub)
  
  !apply pbc at front 5 !y
  !red fluid
  call apply_pbc_front(f00sub)
  call apply_pbc_front(f01sub)
  call apply_pbc_front(f02sub)
  call apply_pbc_front(f03sub)
  call apply_pbc_front(f04sub)
  call apply_pbc_front(f05sub)
  call apply_pbc_front(f06sub)
  call apply_pbc_front(f07sub)
  call apply_pbc_front(f08sub)
  call apply_pbc_front(f09sub)
  call apply_pbc_front(f10sub)
  call apply_pbc_front(f11sub)
  call apply_pbc_front(f12sub)
  call apply_pbc_front(f13sub)
  call apply_pbc_front(f14sub)
  call apply_pbc_front(f15sub)
  call apply_pbc_front(f16sub)
  call apply_pbc_front(f17sub)
  call apply_pbc_front(f18sub)
  
  !apply pbc at rear 6 !y
  !red fluid
  call apply_pbc_rear(f00sub)
  call apply_pbc_rear(f01sub)
  call apply_pbc_rear(f02sub)
  call apply_pbc_rear(f03sub)
  call apply_pbc_rear(f04sub)
  call apply_pbc_rear(f05sub)
  call apply_pbc_rear(f06sub)
  call apply_pbc_rear(f07sub)
  call apply_pbc_rear(f08sub)
  call apply_pbc_rear(f09sub)
  call apply_pbc_rear(f10sub)
  call apply_pbc_rear(f11sub)
  call apply_pbc_rear(f12sub)
  call apply_pbc_rear(f13sub)
  call apply_pbc_rear(f14sub)
  call apply_pbc_rear(f15sub)
  call apply_pbc_rear(f16sub)
  call apply_pbc_rear(f17sub)
  call apply_pbc_rear(f18sub)
  
  !edges 4 (12)
  
  !apply pbc at front east 1 !xy
  !red fluid
  call apply_pbc_edge_front_east(f00sub)
  call apply_pbc_edge_front_east(f01sub)
  call apply_pbc_edge_front_east(f02sub)
  call apply_pbc_edge_front_east(f03sub)
  call apply_pbc_edge_front_east(f04sub)
  call apply_pbc_edge_front_east(f05sub)
  call apply_pbc_edge_front_east(f06sub)
  call apply_pbc_edge_front_east(f07sub)
  call apply_pbc_edge_front_east(f08sub)
  call apply_pbc_edge_front_east(f09sub)
  call apply_pbc_edge_front_east(f10sub)
  call apply_pbc_edge_front_east(f11sub)
  call apply_pbc_edge_front_east(f12sub)
  call apply_pbc_edge_front_east(f13sub)
  call apply_pbc_edge_front_east(f14sub)
  call apply_pbc_edge_front_east(f15sub)
  call apply_pbc_edge_front_east(f16sub)
  call apply_pbc_edge_front_east(f17sub)
  call apply_pbc_edge_front_east(f18sub)
  
  !apply pbc at front west 2 !xy
  !red fluid
  call apply_pbc_edge_front_west(f00sub)
  call apply_pbc_edge_front_west(f01sub)
  call apply_pbc_edge_front_west(f02sub)
  call apply_pbc_edge_front_west(f03sub)
  call apply_pbc_edge_front_west(f04sub)
  call apply_pbc_edge_front_west(f05sub)
  call apply_pbc_edge_front_west(f06sub)
  call apply_pbc_edge_front_west(f07sub)
  call apply_pbc_edge_front_west(f08sub)
  call apply_pbc_edge_front_west(f09sub)
  call apply_pbc_edge_front_west(f10sub)
  call apply_pbc_edge_front_west(f11sub)
  call apply_pbc_edge_front_west(f12sub)
  call apply_pbc_edge_front_west(f13sub)
  call apply_pbc_edge_front_west(f14sub)
  call apply_pbc_edge_front_west(f15sub)
  call apply_pbc_edge_front_west(f16sub)
  call apply_pbc_edge_front_west(f17sub)
  call apply_pbc_edge_front_west(f18sub)
  
  !apply pbc at rear east 7 !xy
  !red fluid
  call apply_pbc_edge_rear_east(f00sub)
  call apply_pbc_edge_rear_east(f01sub)
  call apply_pbc_edge_rear_east(f02sub)
  call apply_pbc_edge_rear_east(f03sub)
  call apply_pbc_edge_rear_east(f04sub)
  call apply_pbc_edge_rear_east(f05sub)
  call apply_pbc_edge_rear_east(f06sub)
  call apply_pbc_edge_rear_east(f07sub)
  call apply_pbc_edge_rear_east(f08sub)
  call apply_pbc_edge_rear_east(f09sub)
  call apply_pbc_edge_rear_east(f10sub)
  call apply_pbc_edge_rear_east(f11sub)
  call apply_pbc_edge_rear_east(f12sub)
  call apply_pbc_edge_rear_east(f13sub)
  call apply_pbc_edge_rear_east(f14sub)
  call apply_pbc_edge_rear_east(f15sub)
  call apply_pbc_edge_rear_east(f16sub)
  call apply_pbc_edge_rear_east(f17sub)
  call apply_pbc_edge_rear_east(f18sub)
  
  !apply pbc at rear west 8 !xy
  !red fluid
  call apply_pbc_edge_rear_west(f00sub)
  call apply_pbc_edge_rear_west(f01sub)
  call apply_pbc_edge_rear_west(f02sub)
  call apply_pbc_edge_rear_west(f03sub)
  call apply_pbc_edge_rear_west(f04sub)
  call apply_pbc_edge_rear_west(f05sub)
  call apply_pbc_edge_rear_west(f06sub)
  call apply_pbc_edge_rear_west(f07sub)
  call apply_pbc_edge_rear_west(f08sub)
  call apply_pbc_edge_rear_west(f09sub)
  call apply_pbc_edge_rear_west(f10sub)
  call apply_pbc_edge_rear_west(f11sub)
  call apply_pbc_edge_rear_west(f12sub)
  call apply_pbc_edge_rear_west(f13sub)
  call apply_pbc_edge_rear_west(f14sub)
  call apply_pbc_edge_rear_west(f15sub)
  call apply_pbc_edge_rear_west(f16sub)
  call apply_pbc_edge_rear_west(f17sub)
  call apply_pbc_edge_rear_west(f18sub)
  
  return
  
 end subroutine apply_pbc_pops_along_xy
 
 subroutine driver_pbc_pops_along_xz
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the periodic boundary 
!     conditions to fluid populations along x and z
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
 
  call apply_pbc_pops_along_xz(f00R,f01R,f02R,f03R,f04R,f05R,f06R,f07R,&
   f08R,f09R,f10R,f11R,f12R,f13R,f14R,f15R,f16R,f17R,f18R) 
  
  if(lsingle_fluid)return
 
  call apply_pbc_pops_along_xz(f00B,f01B,f02B,f03B,f04B,f05B,f06B,f07B,&
   f08B,f09B,f10B,f11B,f12B,f13B,f14B,f15B,f16B,f17B,f18B)
   
  return
  
 end subroutine driver_pbc_pops_along_xz
 
 subroutine apply_pbc_pops_along_xz(f00sub,f01sub,f02sub,f03sub,f04sub,f05sub, &
  f06sub,f07sub,f08sub,f09sub,f10sub,f11sub,f12sub,f13sub,f14sub, &
  f15sub,f16sub,f17sub,f18sub)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the periodic boundary 
!     conditions to fluid populations along x and z
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), allocatable, dimension(:,:,:)  :: f00sub,f01sub, &
  f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub,f09sub,f10sub, &
  f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub
  
  !sides 4 (6)
  
  !apply pbc at north 1 !z
  !red fluid
  call apply_pbc_north(f00sub)
  call apply_pbc_north(f01sub)
  call apply_pbc_north(f02sub)
  call apply_pbc_north(f03sub)
  call apply_pbc_north(f04sub)
  call apply_pbc_north(f05sub)
  call apply_pbc_north(f06sub)
  call apply_pbc_north(f07sub)
  call apply_pbc_north(f08sub)
  call apply_pbc_north(f09sub)
  call apply_pbc_north(f10sub)
  call apply_pbc_north(f11sub)
  call apply_pbc_north(f12sub)
  call apply_pbc_north(f13sub)
  call apply_pbc_north(f14sub)
  call apply_pbc_north(f15sub)
  call apply_pbc_north(f16sub)
  call apply_pbc_north(f17sub)
  call apply_pbc_north(f18sub)
  
  !apply pbc at south 2 !z
  !red fluid
  call apply_pbc_south(f00sub)
  call apply_pbc_south(f01sub)
  call apply_pbc_south(f02sub)
  call apply_pbc_south(f03sub)
  call apply_pbc_south(f04sub)
  call apply_pbc_south(f05sub)
  call apply_pbc_south(f06sub)
  call apply_pbc_south(f07sub)
  call apply_pbc_south(f08sub)
  call apply_pbc_south(f09sub)
  call apply_pbc_south(f10sub)
  call apply_pbc_south(f11sub)
  call apply_pbc_south(f12sub)
  call apply_pbc_south(f13sub)
  call apply_pbc_south(f14sub)
  call apply_pbc_south(f15sub)
  call apply_pbc_south(f16sub)
  call apply_pbc_south(f17sub)
  call apply_pbc_south(f18sub)
  
  !apply pbc at east 3 !x
  !red fluid
  call apply_pbc_east(f00sub)
  call apply_pbc_east(f01sub)
  call apply_pbc_east(f02sub)
  call apply_pbc_east(f03sub)
  call apply_pbc_east(f04sub)
  call apply_pbc_east(f05sub)
  call apply_pbc_east(f06sub)
  call apply_pbc_east(f07sub)
  call apply_pbc_east(f08sub)
  call apply_pbc_east(f09sub)
  call apply_pbc_east(f10sub)
  call apply_pbc_east(f11sub)
  call apply_pbc_east(f12sub)
  call apply_pbc_east(f13sub)
  call apply_pbc_east(f14sub)
  call apply_pbc_east(f15sub)
  call apply_pbc_east(f16sub)
  call apply_pbc_east(f17sub)
  call apply_pbc_east(f18sub)
  
  !apply pbc at west 4 !x
  !red fluid
  call apply_pbc_west(f00sub)
  call apply_pbc_west(f01sub)
  call apply_pbc_west(f02sub)
  call apply_pbc_west(f03sub)
  call apply_pbc_west(f04sub)
  call apply_pbc_west(f05sub)
  call apply_pbc_west(f06sub)
  call apply_pbc_west(f07sub)
  call apply_pbc_west(f08sub)
  call apply_pbc_west(f09sub)
  call apply_pbc_west(f10sub)
  call apply_pbc_west(f11sub)
  call apply_pbc_west(f12sub)
  call apply_pbc_west(f13sub)
  call apply_pbc_west(f14sub)
  call apply_pbc_west(f15sub)
  call apply_pbc_west(f16sub)
  call apply_pbc_west(f17sub)
  call apply_pbc_west(f18sub)
  
  !edges 4 (12)
  
  !apply pbc at north east 3 !xz
  !red fluid
  call apply_pbc_edge_north_east(f00sub)
  call apply_pbc_edge_north_east(f01sub)
  call apply_pbc_edge_north_east(f02sub)
  call apply_pbc_edge_north_east(f03sub)
  call apply_pbc_edge_north_east(f04sub)
  call apply_pbc_edge_north_east(f05sub)
  call apply_pbc_edge_north_east(f06sub)
  call apply_pbc_edge_north_east(f07sub)
  call apply_pbc_edge_north_east(f08sub)
  call apply_pbc_edge_north_east(f09sub)
  call apply_pbc_edge_north_east(f10sub)
  call apply_pbc_edge_north_east(f11sub)
  call apply_pbc_edge_north_east(f12sub)
  call apply_pbc_edge_north_east(f13sub)
  call apply_pbc_edge_north_east(f14sub)
  call apply_pbc_edge_north_east(f15sub)
  call apply_pbc_edge_north_east(f16sub)
  call apply_pbc_edge_north_east(f17sub)
  call apply_pbc_edge_north_east(f18sub)
  
  !apply pbc at north west 6 !xz
  !red fluid
  call apply_pbc_edge_north_west(f00sub)
  call apply_pbc_edge_north_west(f01sub)
  call apply_pbc_edge_north_west(f02sub)
  call apply_pbc_edge_north_west(f03sub)
  call apply_pbc_edge_north_west(f04sub)
  call apply_pbc_edge_north_west(f05sub)
  call apply_pbc_edge_north_west(f06sub)
  call apply_pbc_edge_north_west(f07sub)
  call apply_pbc_edge_north_west(f08sub)
  call apply_pbc_edge_north_west(f09sub)
  call apply_pbc_edge_north_west(f10sub)
  call apply_pbc_edge_north_west(f11sub)
  call apply_pbc_edge_north_west(f12sub)
  call apply_pbc_edge_north_west(f13sub)
  call apply_pbc_edge_north_west(f14sub)
  call apply_pbc_edge_north_west(f15sub)
  call apply_pbc_edge_north_west(f16sub)
  call apply_pbc_edge_north_west(f17sub)
  call apply_pbc_edge_north_west(f18sub)
  
  !apply pbc at south east 9 !xz
  !red fluid
  call apply_pbc_edge_south_east(f00sub)
  call apply_pbc_edge_south_east(f01sub)
  call apply_pbc_edge_south_east(f02sub)
  call apply_pbc_edge_south_east(f03sub)
  call apply_pbc_edge_south_east(f04sub)
  call apply_pbc_edge_south_east(f05sub)
  call apply_pbc_edge_south_east(f06sub)
  call apply_pbc_edge_south_east(f07sub)
  call apply_pbc_edge_south_east(f08sub)
  call apply_pbc_edge_south_east(f09sub)
  call apply_pbc_edge_south_east(f10sub)
  call apply_pbc_edge_south_east(f11sub)
  call apply_pbc_edge_south_east(f12sub)
  call apply_pbc_edge_south_east(f13sub)
  call apply_pbc_edge_south_east(f14sub)
  call apply_pbc_edge_south_east(f15sub)
  call apply_pbc_edge_south_east(f16sub)
  call apply_pbc_edge_south_east(f17sub)
  call apply_pbc_edge_south_east(f18sub)
  
  !apply pbc at south west 12 !xz
  !red fluid
  call apply_pbc_edge_south_west(f00sub)
  call apply_pbc_edge_south_west(f01sub)
  call apply_pbc_edge_south_west(f02sub)
  call apply_pbc_edge_south_west(f03sub)
  call apply_pbc_edge_south_west(f04sub)
  call apply_pbc_edge_south_west(f05sub)
  call apply_pbc_edge_south_west(f06sub)
  call apply_pbc_edge_south_west(f07sub)
  call apply_pbc_edge_south_west(f08sub)
  call apply_pbc_edge_south_west(f09sub)
  call apply_pbc_edge_south_west(f10sub)
  call apply_pbc_edge_south_west(f11sub)
  call apply_pbc_edge_south_west(f12sub)
  call apply_pbc_edge_south_west(f13sub)
  call apply_pbc_edge_south_west(f14sub)
  call apply_pbc_edge_south_west(f15sub)
  call apply_pbc_edge_south_west(f16sub)
  call apply_pbc_edge_south_west(f17sub)
  call apply_pbc_edge_south_west(f18sub)
  
  return
  
 end subroutine apply_pbc_pops_along_xz
 
 subroutine driver_pbc_pops_along_yz
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the periodic boundary 
!     conditions to fluid populations along y and z
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
 
  call apply_pbc_pops_along_yz(f00R,f01R,f02R,f03R,f04R,f05R,f06R,f07R,&
   f08R,f09R,f10R,f11R,f12R,f13R,f14R,f15R,f16R,f17R,f18R) 
  
  if(lsingle_fluid)return
 
  call apply_pbc_pops_along_yz(f00B,f01B,f02B,f03B,f04B,f05B,f06B,f07B,&
   f08B,f09B,f10B,f11B,f12B,f13B,f14B,f15B,f16B,f17B,f18B)
   
  return
  
 end subroutine driver_pbc_pops_along_yz
 
 subroutine apply_pbc_pops_along_yz(f00sub,f01sub,f02sub,f03sub,f04sub, &
  f05sub,f06sub,f07sub,f08sub,f09sub,f10sub,f11sub,f12sub,f13sub, &
  f14sub,f15sub,f16sub,f17sub,f18sub)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the periodic boundary 
!     conditions to fluid populations along y and z
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), allocatable, dimension(:,:,:)  :: f00sub,f01sub, &
  f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub,f09sub,f10sub, &
  f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub
  
  !sides 4 (6)
  
  !apply pbc at north 1 !z
  !red fluid
  call apply_pbc_north(f00sub)
  call apply_pbc_north(f01sub)
  call apply_pbc_north(f02sub)
  call apply_pbc_north(f03sub)
  call apply_pbc_north(f04sub)
  call apply_pbc_north(f05sub)
  call apply_pbc_north(f06sub)
  call apply_pbc_north(f07sub)
  call apply_pbc_north(f08sub)
  call apply_pbc_north(f09sub)
  call apply_pbc_north(f10sub)
  call apply_pbc_north(f11sub)
  call apply_pbc_north(f12sub)
  call apply_pbc_north(f13sub)
  call apply_pbc_north(f14sub)
  call apply_pbc_north(f15sub)
  call apply_pbc_north(f16sub)
  call apply_pbc_north(f17sub)
  call apply_pbc_north(f18sub)
  
  !apply pbc at south 2 !z
  !red fluid
  call apply_pbc_south(f00sub)
  call apply_pbc_south(f01sub)
  call apply_pbc_south(f02sub)
  call apply_pbc_south(f03sub)
  call apply_pbc_south(f04sub)
  call apply_pbc_south(f05sub)
  call apply_pbc_south(f06sub)
  call apply_pbc_south(f07sub)
  call apply_pbc_south(f08sub)
  call apply_pbc_south(f09sub)
  call apply_pbc_south(f10sub)
  call apply_pbc_south(f11sub)
  call apply_pbc_south(f12sub)
  call apply_pbc_south(f13sub)
  call apply_pbc_south(f14sub)
  call apply_pbc_south(f15sub)
  call apply_pbc_south(f16sub)
  call apply_pbc_south(f17sub)
  call apply_pbc_south(f18sub)
  
  !apply pbc at front 5 !y
  !red fluid
  call apply_pbc_front(f00sub)
  call apply_pbc_front(f01sub)
  call apply_pbc_front(f02sub)
  call apply_pbc_front(f03sub)
  call apply_pbc_front(f04sub)
  call apply_pbc_front(f05sub)
  call apply_pbc_front(f06sub)
  call apply_pbc_front(f07sub)
  call apply_pbc_front(f08sub)
  call apply_pbc_front(f09sub)
  call apply_pbc_front(f10sub)
  call apply_pbc_front(f11sub)
  call apply_pbc_front(f12sub)
  call apply_pbc_front(f13sub)
  call apply_pbc_front(f14sub)
  call apply_pbc_front(f15sub)
  call apply_pbc_front(f16sub)
  call apply_pbc_front(f17sub)
  call apply_pbc_front(f18sub)
  
  !apply pbc at rear 6 !y
  !red fluid
  call apply_pbc_rear(f00sub)
  call apply_pbc_rear(f01sub)
  call apply_pbc_rear(f02sub)
  call apply_pbc_rear(f03sub)
  call apply_pbc_rear(f04sub)
  call apply_pbc_rear(f05sub)
  call apply_pbc_rear(f06sub)
  call apply_pbc_rear(f07sub)
  call apply_pbc_rear(f08sub)
  call apply_pbc_rear(f09sub)
  call apply_pbc_rear(f10sub)
  call apply_pbc_rear(f11sub)
  call apply_pbc_rear(f12sub)
  call apply_pbc_rear(f13sub)
  call apply_pbc_rear(f14sub)
  call apply_pbc_rear(f15sub)
  call apply_pbc_rear(f16sub)
  call apply_pbc_rear(f17sub)
  call apply_pbc_rear(f18sub)
  
  !edges 4 (12)
  
  !apply pbc at north front 4 !yz
  !red fluid
  call apply_pbc_edge_north_front(f00sub)
  call apply_pbc_edge_north_front(f01sub)
  call apply_pbc_edge_north_front(f02sub)
  call apply_pbc_edge_north_front(f03sub)
  call apply_pbc_edge_north_front(f04sub)
  call apply_pbc_edge_north_front(f05sub)
  call apply_pbc_edge_north_front(f06sub)
  call apply_pbc_edge_north_front(f07sub)
  call apply_pbc_edge_north_front(f08sub)
  call apply_pbc_edge_north_front(f09sub)
  call apply_pbc_edge_north_front(f10sub)
  call apply_pbc_edge_north_front(f11sub)
  call apply_pbc_edge_north_front(f12sub)
  call apply_pbc_edge_north_front(f13sub)
  call apply_pbc_edge_north_front(f14sub)
  call apply_pbc_edge_north_front(f15sub)
  call apply_pbc_edge_north_front(f16sub)
  call apply_pbc_edge_north_front(f17sub)
  call apply_pbc_edge_north_front(f18sub)
  
  !apply pbc at north rear 5 !yz
  !red fluid
  call apply_pbc_edge_north_rear(f00sub)
  call apply_pbc_edge_north_rear(f01sub)
  call apply_pbc_edge_north_rear(f02sub)
  call apply_pbc_edge_north_rear(f03sub)
  call apply_pbc_edge_north_rear(f04sub)
  call apply_pbc_edge_north_rear(f05sub)
  call apply_pbc_edge_north_rear(f06sub)
  call apply_pbc_edge_north_rear(f07sub)
  call apply_pbc_edge_north_rear(f08sub)
  call apply_pbc_edge_north_rear(f09sub)
  call apply_pbc_edge_north_rear(f10sub)
  call apply_pbc_edge_north_rear(f11sub)
  call apply_pbc_edge_north_rear(f12sub)
  call apply_pbc_edge_north_rear(f13sub)
  call apply_pbc_edge_north_rear(f14sub)
  call apply_pbc_edge_north_rear(f15sub)
  call apply_pbc_edge_north_rear(f16sub)
  call apply_pbc_edge_north_rear(f17sub)
  call apply_pbc_edge_north_rear(f18sub)
  
  !apply pbc at south front 10 !yz
  !red fluid
  call apply_pbc_edge_south_front(f00sub)
  call apply_pbc_edge_south_front(f01sub)
  call apply_pbc_edge_south_front(f02sub)
  call apply_pbc_edge_south_front(f03sub)
  call apply_pbc_edge_south_front(f04sub)
  call apply_pbc_edge_south_front(f05sub)
  call apply_pbc_edge_south_front(f06sub)
  call apply_pbc_edge_south_front(f07sub)
  call apply_pbc_edge_south_front(f08sub)
  call apply_pbc_edge_south_front(f09sub)
  call apply_pbc_edge_south_front(f10sub)
  call apply_pbc_edge_south_front(f11sub)
  call apply_pbc_edge_south_front(f12sub)
  call apply_pbc_edge_south_front(f13sub)
  call apply_pbc_edge_south_front(f14sub)
  call apply_pbc_edge_south_front(f15sub)
  call apply_pbc_edge_south_front(f16sub)
  call apply_pbc_edge_south_front(f17sub)
  call apply_pbc_edge_south_front(f18sub)
  
  !apply pbc at south rear 11 !yz
  !red fluid
  call apply_pbc_edge_south_rear(f00sub)
  call apply_pbc_edge_south_rear(f01sub)
  call apply_pbc_edge_south_rear(f02sub)
  call apply_pbc_edge_south_rear(f03sub)
  call apply_pbc_edge_south_rear(f04sub)
  call apply_pbc_edge_south_rear(f05sub)
  call apply_pbc_edge_south_rear(f06sub)
  call apply_pbc_edge_south_rear(f07sub)
  call apply_pbc_edge_south_rear(f08sub)
  call apply_pbc_edge_south_rear(f09sub)
  call apply_pbc_edge_south_rear(f10sub)
  call apply_pbc_edge_south_rear(f11sub)
  call apply_pbc_edge_south_rear(f12sub)
  call apply_pbc_edge_south_rear(f13sub)
  call apply_pbc_edge_south_rear(f14sub)
  call apply_pbc_edge_south_rear(f15sub)
  call apply_pbc_edge_south_rear(f16sub)
  call apply_pbc_edge_south_rear(f17sub)
  call apply_pbc_edge_south_rear(f18sub)
  
  return
  
 end subroutine apply_pbc_pops_along_yz
 
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
 
!******************END PART TO MANAGE THE PERIODIC BC*******************
 
!*****************START PART TO MANAGE THE BOUNCEBACK*******************
 
 subroutine driver_apply_bounceback_pop
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the bounce back of the fluid
!     populations if requested from the boundary conditions
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  select case(ibctype)
  case(0) ! 0 0 0 
    call apply_bounceback_pops_all
  case(1) ! 1 0 0 
    call apply_bounceback_pops_along_yz
  case(2) ! 0 1 0 
    call apply_bounceback_pops_along_xz
  case(3) ! 1 1 0 
    call apply_bounceback_pops_along_z
  case(4) ! 0 0 1 
    call apply_bounceback_pops_along_xy
  case(5) ! 1 0 1 
    call apply_bounceback_pops_along_y
  case(6) ! 0 1 1 
    call apply_bounceback_pops_along_x
  case(7) ! 1 1 1
    return 
  case default
    call error(12)
  end select
  
  return
  
 end subroutine driver_apply_bounceback_pop
 
 subroutine apply_bounceback_pops_all 
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the bounce back to all 
!     the side of box
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
#if LATTICE==319
  
  !sides 6
  
  !apply bounceback at north 1   !z
  call driver_bounceback_north(1,nx,1,ny)
  
  !apply bounceback at south 2   !z
  call driver_bounceback_south(1,nx,1,ny)
  
  !apply bounceback at east 3    !x
  call driver_bounceback_east(1,ny,1,nz)
  
  !apply bounceback at west 4    !x
  call driver_bounceback_west(1,ny,1,nz)
  
  !apply bounceback at front 5   !y
  call driver_bounceback_front(1,nx,1,nz)
  
  !apply bounceback at rear 6   !y
  call driver_bounceback_rear(1,nx,1,nz)
  
  !edges 12
  
  !apply bounceback at front east 1   !xy
  call driver_bounceback_edge_front_east(1,nz)
  
  !apply bounceback at front west 2   !xy
  call driver_bounceback_edge_front_west(1,nz)
  
  !apply bounceback at north east 3   !xz
  call driver_bounceback_edge_north_east(1,ny)
  
  !apply bounceback at north front 4   !yz
  call driver_bounceback_edge_north_front(1,nx)
  
  !apply bounceback at north rear 5   !yz
  call driver_bounceback_edge_north_rear(1,nx)
  
  !apply bounceback at north west 6   !xz
  call driver_bounceback_edge_north_west(1,ny)
  
  !apply bounceback at rear east 7   !xy
  call driver_bounceback_edge_rear_east(1,nz)
  
  !apply bounceback at rear west 8   !xy
  call driver_bounceback_edge_rear_west(1,nz)
  
  !apply bounceback at south east 9  !xz
  call driver_bounceback_edge_south_east(1,ny)
  
  !apply bounceback at south front 10  !yz
  call driver_bounceback_edge_south_front(1,nx)
  
  !apply bounceback at south rear 11  !yz
  call driver_bounceback_edge_south_rear(1,nx)
  
  !apply bounceback at south west 12  !xz
  call driver_bounceback_edge_south_west(1,ny)
  
  !corner 8
  
  !apply bounceback at north east front 1
  call driver_bounceback_corner_north_east_front
  
  !apply bounceback at north east rear 2
  call driver_bounceback_corner_north_east_rear
  
  !apply bounceback at north west rear 3
  call driver_bounceback_corner_north_west_rear
  
  !apply bounceback at north west front 4
  call driver_bounceback_corner_north_west_front
  
  !apply bounceback at south east front 5
  call driver_bounceback_corner_south_east_front
  
  !apply bounceback at south west front 6
  call driver_bounceback_corner_south_west_front
  
  !apply bounceback at south west rear 7
  call driver_bounceback_corner_south_west_rear
  
  !apply bounceback at south east rear 8
  call driver_bounceback_corner_south_east_rear
  
#else
  
  call warning(9,ZERO,latt_name)
  call error(12)
  
#endif
  
  return
  
 end subroutine apply_bounceback_pops_all
 
 subroutine apply_bounceback_pops_along_yz
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the bounce back along y and z
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
#if LATTICE==319
  
  !sides 4 (6)
  
  !apply bounceback at north 1   !z
  call driver_bounceback_north(1-nbuff,nx+nbuff,1,ny) !frame_x
  
  !apply bounceback at south 2   !z
  call driver_bounceback_south(1-nbuff,nx+nbuff,1,ny) !frame_x
  
  !apply bounceback at front 5   !y
  call driver_bounceback_front(1-nbuff,nx+nbuff,1,nz) !frame_x
  
  !apply bounceback at rear 6   !y
  call driver_bounceback_rear(1-nbuff,nx+nbuff,1,nz) !frame_x
  
  !edges 4 (12)
  
  !apply bounceback at north front 4   !yz
  call driver_bounceback_edge_north_front(1-nbuff,nx+nbuff) !frame
  
  !apply bounceback at north rear 5   !yz
  call driver_bounceback_edge_north_rear(1-nbuff,nx+nbuff) !frame
  
  !apply bounceback at south front 10  !yz
  call driver_bounceback_edge_south_front(1-nbuff,nx+nbuff) !frame
  
  !apply bounceback at south rear 11  !yz
  call driver_bounceback_edge_south_rear(1-nbuff,nx+nbuff) !frame
  
#else
  
  call warning(9,ZERO,latt_name)
  call error(12)
  
#endif
  
  return
  
 end subroutine apply_bounceback_pops_along_yz
 
 subroutine apply_bounceback_pops_along_xz
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the bounce back along x and z
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
#if LATTICE==319
  
  !sides 4 (6)
  
  !apply bounceback at north 1   !z
  call driver_bounceback_north(1,nx,1-nbuff,ny+nbuff) !frame_y
  
  !apply bounceback at south 2   !z
  call driver_bounceback_south(1,nx,1-nbuff,ny+nbuff) !frame_y
  
  !apply bounceback at east 3    !x
  call driver_bounceback_east(1-nbuff,ny+nbuff,1,nz) !frame_y
  
  !apply bounceback at west 4    !x
  call driver_bounceback_west(1-nbuff,ny+nbuff,1,nz) !frame_y
  
  !edges 4 (12)
  
  !apply bounceback at north east 3   !xz
  call driver_bounceback_edge_north_east(1-nbuff,ny+nbuff) !frame
  
  !apply bounceback at north west 6   !xz
  call driver_bounceback_edge_north_west(1-nbuff,ny+nbuff) !frame
  
  !apply bounceback at south east 9  !xz
  call driver_bounceback_edge_south_east(1-nbuff,ny+nbuff) !frame
  
  !apply bounceback at south west 12  !xz
  call driver_bounceback_edge_south_west(1-nbuff,ny+nbuff) !frame
  
#else
  
  call warning(9,ZERO,latt_name)
  call error(12)
  
#endif
  
  return
  
 end subroutine apply_bounceback_pops_along_xz
 
 subroutine apply_bounceback_pops_along_xy
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the bounce back along x and z
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
#if LATTICE==319
  
  !sides 4 (6)
  
  !apply bounceback at east 3    !x
  call driver_bounceback_east(1,ny,1-nbuff,nz+nbuff) !frame_z
  
  !apply bounceback at west 4    !x
  call driver_bounceback_west(1,ny,1-nbuff,nz+nbuff) !frame_z
  
  !apply bounceback at front 5   !y
  call driver_bounceback_front(1,nx,1-nbuff,nz+nbuff) !frame_z
  
  !apply bounceback at rear 6   !y
  call driver_bounceback_rear(1,nx,1-nbuff,nz+nbuff) !frame_z
  
  !edges 4 (12)
  
  !apply bounceback at front east 1   !xy
  call driver_bounceback_edge_front_east(1-nbuff,nz+nbuff) !frame
  
  !apply bounceback at front west 2   !xy
  call driver_bounceback_edge_front_west(1-nbuff,nz+nbuff) !frame
  
  !apply bounceback at rear east 7   !xy
  call driver_bounceback_edge_rear_east(1-nbuff,nz+nbuff) !frame
  
  !apply bounceback at rear west 8   !xy
  call driver_bounceback_edge_rear_west(1-nbuff,nz+nbuff) !frame
  
#else
  
  call warning(9,ZERO,latt_name)
  call error(12)
  
#endif
  
  return
  
 end subroutine apply_bounceback_pops_along_xy
 
 subroutine apply_bounceback_pops_along_z 
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the bounce back 
!     conditions to density variables along z
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  
#if LATTICE==319
  
  !sides 6
  
  !apply pbc at north 1 !z
  call driver_bounceback_north(1-nbuff,nx+nbuff,1-nbuff,ny+nbuff) !frame
  
  !apply pbc at south 2 !z
  call driver_bounceback_south(1-nbuff,nx+nbuff,1-nbuff,ny+nbuff) !frame
  
#else
  
  call warning(9,ZERO,latt_name)
  call error(12)
  
#endif
  
  return
  
 end subroutine apply_bounceback_pops_along_z
 
 subroutine apply_bounceback_pops_along_y 
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the bounce back 
!     conditions to density variables along y
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
#if LATTICE==319
  
  !sides 6
  
  !apply pbc at front 5 !y
  call driver_bounceback_front(1-nbuff,nx+nbuff,1-nbuff,nz+nbuff) !frame
  
  !apply pbc at rear 6 !y
  call driver_bounceback_rear(1-nbuff,nx+nbuff,1-nbuff,nz+nbuff) !frame
  
#else
  
  call warning(9,ZERO,latt_name)
  call error(12)
  
#endif
  
  return
  
 end subroutine apply_bounceback_pops_along_y
 
 subroutine apply_bounceback_pops_along_x 
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the bounce back
!     conditions to density variables along x
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
#if LATTICE==319
  
  !sides 6
  
  !apply pbc at east 3 !x
  call driver_bounceback_east(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff) !frame
  
  !apply pbc at west 4 !x
  call driver_bounceback_west(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff) !frame
  
#else
  
  call warning(9,ZERO,latt_name)
  call error(12)
  
#endif
  
  return
  
 end subroutine apply_bounceback_pops_along_x
 
 subroutine driver_bounceback_east(sy,ey,sz,ez)
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the bounceback
!     at the east side along the space defined in input sy,ey,sz,ez
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: sy,ey,sz,ez
  
  call apply_bounceback_east(sy,ey,sz,ez,bc_rhoR_east,bc_u_east, &
   bc_v_east,bc_w_east,rhoR,u,v,w,f00R,f01R,f02R,f03R,f04R,f05R,f06R, &
   f07R,f08R,f09R,f10R,f11R,f12R,f13R,f14R,f15R,f16R,f17R,f18R)
  
  if(lsingle_fluid)return
 
  call apply_bounceback_east(sy,ey,sz,ez,bc_rhoB_east,bc_u_east, &
   bc_v_east,bc_w_east,rhoB,u,v,w,f00B,f01B,f02B,f03B,f04B,f05B,f06B, &
   f07B,f08B,f09B,f10B,f11B,f12B,f13B,f14B,f15B,f16B,f17B,f18B)
   
  return
  
 end subroutine driver_bounceback_east
 
 subroutine apply_bounceback_east(sy,ey,sz,ez,rho_east,u_east, &
  v_east,w_east,rho_s,u_s,v_s,w_s,f00sub,f01sub,f02sub,f03sub,f04sub, &
  f05sub,f06sub,f07sub,f08sub,f09sub,f10sub,f11sub,f12sub,f13sub, &
  f14sub,f15sub,f16sub,f17sub,f18sub)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the bounceback
!     at the east side along the space defined in input sy,ey,sz,ez
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: sy,ey,sz,ez
  real(kind=PRC), intent(in) :: rho_east,u_east,v_east,w_east
  real(kind=PRC), allocatable, dimension(:,:,:)  :: rho_s,u_s,v_s,w_s, &
   f00sub,f01sub,f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub, &
   f09sub,f10sub,f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub
  
  integer :: i,j,k
  integer, parameter :: kk = 1
  
  real(kind=PRC), parameter :: pref_bouzidi=TWO/cssq
  real(kind=PRC), parameter :: cssq2 = ( HALF / cssq )
  real(kind=PRC), parameter :: cssq4 = ( HALF / (cssq*cssq) )
  
#if LATTICE==319
  
  select case (bc_type_east)  
  case (0)
  
    ! swap pop1 and pop2 along x
    forall(j=sy:ey,k=sz:ez)
      buffservice3d(nx+kk,j,k)=f01sub(nx+kk,j,k)
      f01sub(nx+kk,j,k)=f02sub(nx+kk,j,k)
      f02sub(nx+kk,j,k)=buffservice3d(nx+kk,j,k)
    end forall
    
    ! swap pop7 and pop8 along x
    forall(j=sy:ey,k=sz:ez)
      buffservice3d(nx+kk,j,k)=f07sub(nx+kk,j,k)
      f07sub(nx+kk,j,k)=f08sub(nx+kk,j,k)
      f08sub(nx+kk,j,k)=buffservice3d(nx+kk,j,k)
    end forall
    
    ! swap pop9 and pop10 along x
    forall(j=sy:ey,k=sz:ez)
      buffservice3d(nx+kk,j,k)=f09sub(nx+kk,j,k)
      f09sub(nx+kk,j,k)=f10sub(nx+kk,j,k)
      f10sub(nx+kk,j,k)=buffservice3d(nx+kk,j,k)
    end forall
    
    ! swap pop11 and pop12 along x
    forall(j=sy:ey,k=sz:ez)
      buffservice3d(nx+kk,j,k)=f11sub(nx+kk,j,k)
      f11sub(nx+kk,j,k)=f12sub(nx+kk,j,k)
      f12sub(nx+kk,j,k)=buffservice3d(nx+kk,j,k)
    end forall
    
    ! swap pop13 and pop14 along x
    forall(j=sy:ey,k=sz:ez)
      buffservice3d(nx+kk,j,k)=f13sub(nx+kk,j,k)
      f13sub(nx+kk,j,k)=f14sub(nx+kk,j,k)
      f14sub(nx+kk,j,k)=buffservice3d(nx+kk,j,k)
    end forall
    
  case (1)
  
    ! swap pop1 and pop2 along x
    forall(j=sy:ey,k=sz:ez)
      buffservice3d(nx+kk,j,k)=-f01sub(nx+kk,j,k)
      f01sub(nx+kk,j,k)=-f02sub(nx+kk,j,k) + p(1)*TWO*rho_east* &
       (ONE+(dex(1)*u_s(nx-kk+1,j,k))**TWO/cssq4 - u_s(nx-kk+1,j,k)**TWO/cssq2)
      f02sub(nx+kk,j,k)=buffservice3d(nx+kk,j,k) + p(2)*TWO*rho_east* &
       (ONE+(dex(2)*u_s(nx-kk+1,j,k))**TWO/cssq4 - u_s(nx-kk+1,j,k)**TWO/cssq2)
    end forall
    
    ! swap pop7 and pop8 along x
    forall(j=sy:ey,k=sz:ez)
      buffservice3d(nx+kk,j,k)=-f07sub(nx+kk,j,k)
      f07sub(nx+kk,j,k)=-f08sub(nx+kk,j,k) + p(7)*TWO*rho_east* &
       (ONE+(dex(7)*u_s(nx-kk+1,j,k))**TWO/cssq4 - u_s(nx-kk+1,j,k)**TWO/cssq2)
      f08sub(nx+kk,j,k)=buffservice3d(nx+kk,j,k) + p(8)*TWO*rho_east* &
       (ONE+(dex(8)*u_s(nx-kk+1,j,k))**TWO/cssq4 - u_s(nx-kk+1,j,k)**TWO/cssq2)
    end forall
    
    ! swap pop9 and pop10 along x
    forall(j=sy:ey,k=sz:ez)
      buffservice3d(nx+kk,j,k)=-f09sub(nx+kk,j,k)
      f09sub(nx+kk,j,k)=-f10sub(nx+kk,j,k) + p(9)*TWO*rho_east* &
       (ONE+(dex(9)*u_s(nx-kk+1,j,k))**TWO/cssq4 - u_s(nx-kk+1,j,k)**TWO/cssq2)
      f10sub(nx+kk,j,k)=buffservice3d(nx+kk,j,k) + p(10)*TWO*rho_east* &
       (ONE+(dex(10)*u_s(nx-kk+1,j,k))**TWO/cssq4 - u_s(nx-kk+1,j,k)**TWO/cssq2)
    end forall
    
    ! swap pop11 and pop12 along x
    forall(j=sy:ey,k=sz:ez)
      buffservice3d(nx+kk,j,k)=-f11sub(nx+kk,j,k)
      f11sub(nx+kk,j,k)=-f12sub(nx+kk,j,k) + p(11)*TWO*rho_east* &
       (ONE+(dex(11)*u_s(nx-kk+1,j,k))**TWO/cssq4 - u_s(nx-kk+1,j,k)**TWO/cssq2)
      f12sub(nx+kk,j,k)=buffservice3d(nx+kk,j,k) + p(12)*TWO*rho_east* &
       (ONE+(dex(12)*u_s(nx-kk+1,j,k))**TWO/cssq4 - u_s(nx-kk+1,j,k)**TWO/cssq2)
    end forall
    
    ! swap pop13 and pop14 along x
    forall(j=sy:ey,k=sz:ez)
      buffservice3d(nx+kk,j,k)=-f13sub(nx+kk,j,k)
      f13sub(nx+kk,j,k)=-f14sub(nx+kk,j,k) + p(13)*TWO*rho_east* &
       (ONE+(dex(13)*u_s(nx-kk+1,j,k))**TWO/cssq4 - u_s(nx-kk+1,j,k)**TWO/cssq2)
      f14sub(nx+kk,j,k)=buffservice3d(nx+kk,j,k) + p(14)*TWO*rho_east* &
       (ONE+(dex(14)*u_s(nx-kk+1,j,k))**TWO/cssq4 - u_s(nx-kk+1,j,k)**TWO/cssq2)
    end forall
    
  case (2)
  
    ! swap pop1 and pop2 along x
    forall(j=sy:ey,k=sz:ez)
      buffservice3d(nx+kk,j,k)=f01sub(nx+kk,j,k)
      f01sub(nx+kk,j,k)=f02sub(nx+kk,j,k) + p(1)*pref_bouzidi*rho_s(nx-kk+1,j,k)* &
       (dex(1)*u_east+dey(1)*v_east+dez(1)*w_east)
      f02sub(nx+kk,j,k)=buffservice3d(nx+kk,j,k) + p(2)*pref_bouzidi*rho_s(nx-kk+1,j,k)* &
       (dex(2)*u_east+dey(2)*v_east+dez(2)*w_east)
    end forall
    
    ! swap pop7 and pop8 along x
    forall(j=sy:ey,k=sz:ez)
      buffservice3d(nx+kk,j,k)=f07sub(nx+kk,j,k)
      f07sub(nx+kk,j,k)=f08sub(nx+kk,j,k) + p(7)*pref_bouzidi*rho_s(nx-kk+1,j,k)* &
       (dex(7)*u_east+dey(7)*v_east+dez(7)*w_east)
      f08sub(nx+kk,j,k)=buffservice3d(nx+kk,j,k) + p(8)*pref_bouzidi*rho_s(nx-kk+1,j,k)* &
       (dex(8)*u_east+dey(8)*v_east+dez(8)*w_east)
    end forall
    
    ! swap pop9 and pop10 along x
    forall(j=sy:ey,k=sz:ez)
      buffservice3d(nx+kk,j,k)=f09sub(nx+kk,j,k)
      f09sub(nx+kk,j,k)=f10sub(nx+kk,j,k) + p(9)*pref_bouzidi*rho_s(nx-kk+1,j,k)* &
       (dex(9)*u_east+dey(9)*v_east+dez(9)*w_east)
      f10sub(nx+kk,j,k)=buffservice3d(nx+kk,j,k) + p(10)*pref_bouzidi*rho_s(nx-kk+1,j,k)* &
       (dex(10)*u_east+dey(10)*v_east+dez(10)*w_east)
    end forall
    
    ! swap pop11 and pop12 along x
    forall(j=sy:ey,k=sz:ez)
      buffservice3d(nx+kk,j,k)=f11sub(nx+kk,j,k)
      f11sub(nx+kk,j,k)=f12sub(nx+kk,j,k) + p(11)*pref_bouzidi*rho_s(nx-kk+1,j,k)* &
       (dex(11)*u_east+dey(11)*v_east+dez(11)*w_east)
      f12sub(nx+kk,j,k)=buffservice3d(nx+kk,j,k) + p(12)*pref_bouzidi*rho_s(nx-kk+1,j,k)* &
       (dex(12)*u_east+dey(12)*v_east+dez(12)*w_east)
    end forall
    
    ! swap pop13 and pop14 along x
    forall(j=sy:ey,k=sz:ez)
      buffservice3d(nx+kk,j,k)=f13sub(nx+kk,j,k)
      f13sub(nx+kk,j,k)=f14sub(nx+kk,j,k) + p(13)*pref_bouzidi*rho_s(nx-kk+1,j,k)* &
       (dex(13)*u_east+dey(13)*v_east+dez(13)*w_east)
      f14sub(nx+kk,j,k)=buffservice3d(nx+kk,j,k) + p(14)*pref_bouzidi*rho_s(nx-kk+1,j,k)* &
       (dex(14)*u_east+dey(14)*v_east+dez(14)*w_east)
    end forall
    
  
  case (3)
  
    !red fluid
    forall(j=sy:ey,k=sz:ez)
      f00sub(nx+kk,j,k)=equil_pop00(rho_east,u_east,v_east,w_east)
    end forall
    
    forall(j=sy:ey,k=sz:ez)
      f01sub(nx+kk,j,k)=equil_pop01(rho_east,u_east,v_east,w_east)
    end forall
    
    forall(j=sy:ey,k=sz:ez)
      f02sub(nx+kk,j,k)=equil_pop02(rho_east,u_east,v_east,w_east)
    end forall
    
    forall(j=sy:ey,k=sz:ez)
      f03sub(nx+kk,j,k)=equil_pop03(rho_east,u_east,v_east,w_east)
    end forall
    
    forall(j=sy:ey,k=sz:ez)
      f04sub(nx+kk,j,k)=equil_pop04(rho_east,u_east,v_east,w_east)
    end forall
    
    forall(j=sy:ey,k=sz:ez)
      f05sub(nx+kk,j,k)=equil_pop05(rho_east,u_east,v_east,w_east)
    end forall
    
    forall(j=sy:ey,k=sz:ez)
      f06sub(nx+kk,j,k)=equil_pop06(rho_east,u_east,v_east,w_east)
    end forall
    
    forall(j=sy:ey,k=sz:ez)
      f07sub(nx+kk,j,k)=equil_pop07(rho_east,u_east,v_east,w_east)
    end forall
    
    forall(j=sy:ey,k=sz:ez)
      f08sub(nx+kk,j,k)=equil_pop08(rho_east,u_east,v_east,w_east)
    end forall
    
    forall(j=sy:ey,k=sz:ez)
      f09sub(nx+kk,j,k)=equil_pop09(rho_east,u_east,v_east,w_east)
    end forall
    
    forall(j=sy:ey,k=sz:ez)
      f10sub(nx+kk,j,k)=equil_pop10(rho_east,u_east,v_east,w_east)
    end forall
    
    forall(j=sy:ey,k=sz:ez)
      f11sub(nx+kk,j,k)=equil_pop11(rho_east,u_east,v_east,w_east)
    end forall
    
    forall(j=sy:ey,k=sz:ez)
      f12sub(nx+kk,j,k)=equil_pop12(rho_east,u_east,v_east,w_east)
    end forall
    
    forall(j=sy:ey,k=sz:ez)
      f13sub(nx+kk,j,k)=equil_pop13(rho_east,u_east,v_east,w_east)
    end forall
    
    forall(j=sy:ey,k=sz:ez)
      f14sub(nx+kk,j,k)=equil_pop14(rho_east,u_east,v_east,w_east)
    end forall
    
    forall(j=sy:ey,k=sz:ez)
      f15sub(nx+kk,j,k)=equil_pop15(rho_east,u_east,v_east,w_east)
    end forall
    
    forall(j=sy:ey,k=sz:ez)
      f16sub(nx+kk,j,k)=equil_pop16(rho_east,u_east,v_east,w_east)
    end forall
    
    forall(j=sy:ey,k=sz:ez)
      f17sub(nx+kk,j,k)=equil_pop17(rho_east,u_east,v_east,w_east)
    end forall
    
    forall(j=sy:ey,k=sz:ez)
      f18sub(nx+kk,j,k)=equil_pop18(rho_east,u_east,v_east,w_east)
    end forall
  
  case default
    call error(13)
  end select
#else
  
  call warning(9,ZERO,latt_name)
  call error(12)
  
#endif
  
  return
  
 end subroutine apply_bounceback_east
 
 subroutine driver_bounceback_west(sy,ey,sz,ez)
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the bounceback
!     at the west side along the space defined in input sy,ey,sz,ez
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: sy,ey,sz,ez
  
  call apply_bounceback_west(sy,ey,sz,ez,bc_rhoR_west,bc_u_west, &
   bc_v_west,bc_w_west,rhoR,u,v,w,f00R,f01R,f02R,f03R,f04R,f05R,f06R, &
   f07R,f08R,f09R,f10R,f11R,f12R,f13R,f14R,f15R,f16R,f17R,f18R)
  
  if(lsingle_fluid)return
 
  call apply_bounceback_west(sy,ey,sz,ez,bc_rhoB_west,bc_u_west, &
   bc_v_west,bc_w_west,rhoB,u,v,w,f00B,f01B,f02B,f03B,f04B,f05B,f06B, &
   f07B,f08B,f09B,f10B,f11B,f12B,f13B,f14B,f15B,f16B,f17B,f18B)
   
  return
  
 end subroutine driver_bounceback_west
 
 subroutine apply_bounceback_west(sy,ey,sz,ez,rho_west,u_west, &
  v_west,w_west,rho_s,u_s,v_s,w_s,f00sub,f01sub,f02sub,f03sub,f04sub, &
  f05sub,f06sub,f07sub,f08sub,f09sub,f10sub,f11sub,f12sub,f13sub, &
  f14sub,f15sub,f16sub,f17sub,f18sub)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the bounceback
!     at the west side along the space defined in input sy,ey,sz,ez
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: sy,ey,sz,ez
  real(kind=PRC), intent(in) :: rho_west,u_west,v_west,w_west
  real(kind=PRC), allocatable, dimension(:,:,:)  :: rho_s,u_s,v_s,w_s, &
   f00sub,f01sub,f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub, &
   f09sub,f10sub,f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub
  
  integer :: i,j,k
  integer, parameter :: kk = 1
  
  real(kind=PRC), parameter :: pref_bouzidi=TWO/cssq
  
  real(kind=PRC), parameter :: cssq2 = ( HALF / cssq )
  real(kind=PRC), parameter :: cssq4 = ( HALF / (cssq*cssq) )
  
#if LATTICE==319
  
  select case (bc_type_west)  
  case (0)
  
    ! swap pop1 and pop2 along x
    forall(j=sy:ey,k=sz:ez)
      buffservice3d(1-kk,j,k)=f01sub(1-kk,j,k)
      f01sub(1-kk,j,k)=f02sub(1-kk,j,k)
      f02sub(1-kk,j,k)=buffservice3d(1-kk,j,k)
    end forall
    
    ! swap pop7 and pop8 along x
    forall(j=sy:ey,k=sz:ez)
      buffservice3d(1-kk,j,k)=f07sub(1-kk,j,k)
      f07sub(1-kk,j,k)=f08sub(1-kk,j,k)
      f08sub(1-kk,j,k)=buffservice3d(1-kk,j,k)
    end forall
    
    ! swap pop9 and pop10 along x
    forall(j=sy:ey,k=sz:ez)
      buffservice3d(1-kk,j,k)=f09sub(1-kk,j,k)
      f09sub(1-kk,j,k)=f10sub(1-kk,j,k)
      f10sub(1-kk,j,k)=buffservice3d(1-kk,j,k)
    end forall
   
    ! swap pop11 and pop12 along x
    forall(j=sy:ey,k=sz:ez)
      buffservice3d(1-kk,j,k)=f11sub(1-kk,j,k)
      f11sub(1-kk,j,k)=f12sub(1-kk,j,k)
      f12sub(1-kk,j,k)=buffservice3d(1-kk,j,k)
    end forall
    
    ! swap pop13 and pop14 along x
    forall(j=sy:ey,k=sz:ez)
      buffservice3d(1-kk,j,k)=f13sub(1-kk,j,k)
      f13sub(1-kk,j,k)=f14sub(1-kk,j,k)
      f14sub(1-kk,j,k)=buffservice3d(1-kk,j,k)
    end forall
    
  case (1)
  
    ! swap pop1 and pop2 along x
    forall(j=sy:ey,k=sz:ez)
      buffservice3d(1-kk,j,k)=-f01sub(1-kk,j,k)
      f01sub(1-kk,j,k)=-f02sub(1-kk,j,k) + p(1)*TWO*rho_west* &
       (ONE+(dex(1)*u_s(kk,j,k))**TWO/cssq4 - u_s(kk,j,k)**TWO/cssq2)
      f02sub(1-kk,j,k)=buffservice3d(1-kk,j,k) + p(2)*TWO*rho_west* &
       (ONE+(dex(2)*u_s(kk,j,k))**TWO/cssq4 - u_s(kk,j,k)**TWO/cssq2)
    end forall
    
    ! swap pop7 and pop8 along x
    forall(j=sy:ey,k=sz:ez)
      buffservice3d(1-kk,j,k)=-f07sub(1-kk,j,k)
      f07sub(1-kk,j,k)=-f08sub(1-kk,j,k) + p(7)*TWO*rho_west* &
       (ONE+(dex(7)*u_s(kk,j,k))**TWO/cssq4 - u_s(kk,j,k)**TWO/cssq2)
      f08sub(1-kk,j,k)=buffservice3d(1-kk,j,k) + p(8)*TWO*rho_west* &
       (ONE+(dex(8)*u_s(kk,j,k))**TWO/cssq4 - u_s(kk,j,k)**TWO/cssq2)
    end forall
    
    ! swap pop9 and pop10 along x
    forall(j=sy:ey,k=sz:ez)
      buffservice3d(1-kk,j,k)=-f09sub(1-kk,j,k)
      f09sub(1-kk,j,k)=-f10sub(1-kk,j,k) + p(9)*TWO*rho_west* &
       (ONE+(dex(9)*u_s(kk,j,k))**TWO/cssq4 - u_s(kk,j,k)**TWO/cssq2)
      f10sub(1-kk,j,k)=buffservice3d(1-kk,j,k) + p(10)*TWO*rho_west* &
       (ONE+(dex(10)*u_s(kk,j,k))**TWO/cssq4 - u_s(kk,j,k)**TWO/cssq2)
    end forall
   
    ! swap pop11 and pop12 along x
    forall(j=sy:ey,k=sz:ez)
      buffservice3d(1-kk,j,k)=-f11sub(1-kk,j,k)
      f11sub(1-kk,j,k)=-f12sub(1-kk,j,k) + p(11)*TWO*rho_west* &
       (ONE+(dex(11)*u_s(kk,j,k))**TWO/cssq4 - u_s(kk,j,k)**TWO/cssq2)
      f12sub(1-kk,j,k)=buffservice3d(1-kk,j,k) + p(12)*TWO*rho_west* &
       (ONE+(dex(12)*u_s(kk,j,k))**TWO/cssq4 - u_s(kk,j,k)**TWO/cssq2)
    end forall
    
    ! swap pop13 and pop14 along x
    forall(j=sy:ey,k=sz:ez)
      buffservice3d(1-kk,j,k)=-f13sub(1-kk,j,k)
      f13sub(1-kk,j,k)=-f14sub(1-kk,j,k) + p(13)*TWO*rho_west* &
       (ONE+(dex(13)*u_s(kk,j,k))**TWO/cssq4 - u_s(kk,j,k)**TWO/cssq2)
      f14sub(1-kk,j,k)=buffservice3d(1-kk,j,k) + p(14)*TWO*rho_west* &
       (ONE+(dex(14)*u_s(kk,j,k))**TWO/cssq4 - u_s(kk,j,k)**TWO/cssq2)
    end forall
  
  case (2)
  
    ! swap pop1 and pop2 along x
    forall(j=sy:ey,k=sz:ez)
      buffservice3d(1-kk,j,k)=f01sub(1-kk,j,k)
      f01sub(1-kk,j,k)=f02sub(1-kk,j,k) + p(1)*pref_bouzidi*rho_s(kk,j,k)* &
       (dex(1)*u_west+dey(1)*v_west+dez(1)*w_west)
      f02sub(1-kk,j,k)=buffservice3d(1-kk,j,k) + p(2)*pref_bouzidi*rho_s(kk,j,k)* &
       (dex(2)*u_west+dey(2)*v_west+dez(2)*w_west)
    end forall
    
    ! swap pop7 and pop8 along x
    forall(j=sy:ey,k=sz:ez)
      buffservice3d(1-kk,j,k)=f07sub(1-kk,j,k)
      f07sub(1-kk,j,k)=f08sub(1-kk,j,k) + p(7)*pref_bouzidi*rho_s(kk,j,k)* &
       (dex(7)*u_west+dey(7)*v_west+dez(7)*w_west)
      f08sub(1-kk,j,k)=buffservice3d(1-kk,j,k) + p(8)*pref_bouzidi*rho_s(kk,j,k)* &
       (dex(8)*u_west+dey(8)*v_west+dez(8)*w_west)
    end forall
    
    ! swap pop9 and pop10 along x
    forall(j=sy:ey,k=sz:ez)
      buffservice3d(1-kk,j,k)=f09sub(1-kk,j,k)
      f09sub(1-kk,j,k)=f10sub(1-kk,j,k) + p(9)*pref_bouzidi*rho_s(kk,j,k)* &
       (dex(9)*u_west+dey(9)*v_west+dez(9)*w_west)
      f10sub(1-kk,j,k)=buffservice3d(1-kk,j,k) + p(10)*pref_bouzidi*rho_s(kk,j,k)* &
       (dex(10)*u_west+dey(10)*v_west+dez(10)*w_west)
    end forall
   
    ! swap pop11 and pop12 along x
    forall(j=sy:ey,k=sz:ez)
      buffservice3d(1-kk,j,k)=f11sub(1-kk,j,k)
      f11sub(1-kk,j,k)=f12sub(1-kk,j,k) + p(11)*pref_bouzidi*rho_s(kk,j,k)* &
       (dex(11)*u_west+dey(11)*v_west+dez(11)*w_west)
      f12sub(1-kk,j,k)=buffservice3d(1-kk,j,k) + p(12)*pref_bouzidi*rho_s(kk,j,k)* &
       (dex(12)*u_west+dey(12)*v_west+dez(12)*w_west)
    end forall
    
    ! swap pop13 and pop14 along x
    forall(j=sy:ey,k=sz:ez)
      buffservice3d(1-kk,j,k)=f13sub(1-kk,j,k)
      f13sub(1-kk,j,k)=f14sub(1-kk,j,k) + p(13)*pref_bouzidi*rho_s(kk,j,k)* &
       (dex(13)*u_west+dey(13)*v_west+dez(13)*w_west)
      f14sub(1-kk,j,k)=buffservice3d(1-kk,j,k) + p(14)*pref_bouzidi*rho_s(kk,j,k)* &
       (dex(14)*u_west+dey(14)*v_west+dez(14)*w_west)
    end forall

  case (3)
  
    !red fluid
    forall(j=sy:ey,k=sz:ez)
      f00sub(1-kk,j,k)=equil_pop00(rho_west,u_west,v_west,w_west)
    end forall
    
    forall(j=sy:ey,k=sz:ez)
      f01sub(1-kk,j,k)=equil_pop01(rho_west,u_west,v_west,w_west)
    end forall
    
    forall(j=sy:ey,k=sz:ez)
      f02sub(1-kk,j,k)=equil_pop02(rho_west,u_west,v_west,w_west)
    end forall
    
    forall(j=sy:ey,k=sz:ez)
      f03sub(1-kk,j,k)=equil_pop03(rho_west,u_west,v_west,w_west)
    end forall
    
    forall(j=sy:ey,k=sz:ez)
      f04sub(1-kk,j,k)=equil_pop04(rho_west,u_west,v_west,w_west)
    end forall
    
    forall(j=sy:ey,k=sz:ez)
      f05sub(1-kk,j,k)=equil_pop05(rho_west,u_west,v_west,w_west)
    end forall
    
    forall(j=sy:ey,k=sz:ez)
      f06sub(1-kk,j,k)=equil_pop06(rho_west,u_west,v_west,w_west)
    end forall
    
    forall(j=sy:ey,k=sz:ez)
      f07sub(1-kk,j,k)=equil_pop07(rho_west,u_west,v_west,w_west)
    end forall
    
    forall(j=sy:ey,k=sz:ez)
      f08sub(1-kk,j,k)=equil_pop08(rho_west,u_west,v_west,w_west)
    end forall
    
    forall(j=sy:ey,k=sz:ez)
      f09sub(1-kk,j,k)=equil_pop09(rho_west,u_west,v_west,w_west)
    end forall
    
    forall(j=sy:ey,k=sz:ez)
      f10sub(1-kk,j,k)=equil_pop10(rho_west,u_west,v_west,w_west)
    end forall
    
    forall(j=sy:ey,k=sz:ez)
      f11sub(1-kk,j,k)=equil_pop11(rho_west,u_west,v_west,w_west)
    end forall
    
    forall(j=sy:ey,k=sz:ez)
      f12sub(1-kk,j,k)=equil_pop12(rho_west,u_west,v_west,w_west)
    end forall
    
    forall(j=sy:ey,k=sz:ez)
      f13sub(1-kk,j,k)=equil_pop13(rho_west,u_west,v_west,w_west)
    end forall
    
    forall(j=sy:ey,k=sz:ez)
      f14sub(1-kk,j,k)=equil_pop14(rho_west,u_west,v_west,w_west)
    end forall
    
    forall(j=sy:ey,k=sz:ez)
      f15sub(1-kk,j,k)=equil_pop15(rho_west,u_west,v_west,w_west)
    end forall
    
    forall(j=sy:ey,k=sz:ez)
      f16sub(1-kk,j,k)=equil_pop16(rho_west,u_west,v_west,w_west)
    end forall
    
    forall(j=sy:ey,k=sz:ez)
      f17sub(1-kk,j,k)=equil_pop17(rho_west,u_west,v_west,w_west)
    end forall
    
    forall(j=sy:ey,k=sz:ez)
      f18sub(1-kk,j,k)=equil_pop18(rho_west,u_west,v_west,w_west)
    end forall
    
  case default
    call error(13)
  end select
  
#else
  
  call warning(9,ZERO,latt_name)
  call error(12)
  
#endif
  
  return
  
 end subroutine apply_bounceback_west
 
 subroutine driver_bounceback_north(sx,ex,sy,ey)
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the bounceback
!     at the north side along the space defined in input sx,ex,sy,ey
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: sx,ex,sy,ey
  
  call apply_bounceback_north(sx,ex,sy,ey,bc_rhoR_north,bc_u_north, &
   bc_v_north,bc_w_north,rhoR,u,v,w,f00R,f01R,f02R,f03R,f04R,f05R,f06R, &
   f07R,f08R,f09R,f10R,f11R,f12R,f13R,f14R,f15R,f16R,f17R,f18R)
  
  if(lsingle_fluid)return
 
  call apply_bounceback_north(sx,ex,sy,ey,bc_rhoB_north,bc_u_north, &
   bc_v_north,bc_w_north,rhoB,u,v,w,f00B,f01B,f02B,f03B,f04B,f05B,f06B, &
   f07B,f08B,f09B,f10B,f11B,f12B,f13B,f14B,f15B,f16B,f17B,f18B)
   
  return
  
 end subroutine driver_bounceback_north
 
 subroutine apply_bounceback_north(sx,ex,sy,ey,rho_north,u_north, &
  v_north,w_north,rho_s,u_s,v_s,w_s,f00sub,f01sub,f02sub,f03sub,f04sub, &
  f05sub,f06sub,f07sub,f08sub,f09sub,f10sub,f11sub,f12sub,f13sub, &
  f14sub,f15sub,f16sub,f17sub,f18sub)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the bounceback
!     at the north side along the space defined in input sx,ex,sy,ey
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: sx,ex,sy,ey
  real(kind=PRC), intent(in) :: rho_north,u_north,v_north,w_north
  real(kind=PRC), allocatable, dimension(:,:,:)  :: rho_s,u_s,v_s,w_s, &
   f00sub,f01sub,f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub, &
   f09sub,f10sub,f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub
  
  integer :: i,j,k
  integer, parameter :: kk = 1
  
  real(kind=PRC), parameter :: pref_bouzidi=TWO/cssq
  
  real(kind=PRC), parameter :: cssq2 = ( HALF / cssq )
  real(kind=PRC), parameter :: cssq4 = ( HALF / (cssq*cssq) )
  
#if LATTICE==319
  
  select case (bc_type_north)  
  case (0)
  
    ! swap pop5 and pop6 along z
    forall(i=sx:ex,j=sy:ey)
      buffservice3d(i,j,nz+kk)=f05sub(i,j,nz+kk)
      f05sub(i,j,nz+kk)=f06sub(i,j,nz+kk)
      f06sub(i,j,nz+kk)=buffservice3d(i,j,nz+kk)
    end forall
  
    ! swap pop11 and pop12 along z
    forall(i=sx:ex,j=sy:ey)
      buffservice3d(i,j,nz+kk)=f11sub(i,j,nz+kk)
      f11sub(i,j,nz+kk)=f12sub(i,j,nz+kk)
      f12sub(i,j,nz+kk)=buffservice3d(i,j,nz+kk)
    end forall
  
    ! swap pop13 and pop14 along z
    forall(i=sx:ex,j=sy:ey)
      buffservice3d(i,j,nz+kk)=f13sub(i,j,nz+kk)
      f13sub(i,j,nz+kk)=f14sub(i,j,nz+kk)
      f14sub(i,j,nz+kk)=buffservice3d(i,j,nz+kk)
    end forall
  
    ! swap pop15 and pop16 along z
    forall(i=sx:ex,j=sy:ey)
      buffservice3d(i,j,nz+kk)=f15sub(i,j,nz+kk)
      f15sub(i,j,nz+kk)=f16sub(i,j,nz+kk)
      f16sub(i,j,nz+kk)=buffservice3d(i,j,nz+kk)
    end forall
  
    ! swap pop17 and pop18 along z
    forall(i=sx:ex,j=sy:ey)
      buffservice3d(i,j,nz+kk)=f17sub(i,j,nz+kk)
      f17sub(i,j,nz+kk)=f18sub(i,j,nz+kk)
      f18sub(i,j,nz+kk)=buffservice3d(i,j,nz+kk)
    end forall
    
  case (1)
  
    ! swap pop5 and pop6 along z
    forall(i=sx:ex,j=sy:ey)
      buffservice3d(i,j,nz+kk)=-f05sub(i,j,nz+kk)
      f05sub(i,j,nz+kk)=-f06sub(i,j,nz+kk) + TWO*p(5)*rho_north* &
       (ONE+(dez(5)*w_s(i,j,nz-kk+1))**TWO/cssq4 - w_s(i,j,nz-kk+1)**TWO/cssq2)
      f06sub(i,j,nz+kk)=buffservice3d(i,j,nz+kk) + TWO*p(6)*rho_north* &
       (ONE+(dez(6)*w_s(i,j,nz-kk+1))**TWO/cssq4 - w_s(i,j,nz-kk+1)**TWO/cssq2)
    end forall
  
    ! swap pop11 and pop12 along z
    forall(i=sx:ex,j=sy:ey)
      buffservice3d(i,j,nz+kk)=-f11sub(i,j,nz+kk)
      f11sub(i,j,nz+kk)=-f12sub(i,j,nz+kk) + TWO*p(11)*rho_north* &
       (ONE+(dez(11)*w_s(i,j,nz-kk+1))**TWO/cssq4 - w_s(i,j,nz-kk+1)**TWO/cssq2)
      f12sub(i,j,nz+kk)=buffservice3d(i,j,nz+kk) + TWO*p(12)*rho_north* &
       (ONE+(dez(12)*w_s(i,j,nz-kk+1))**TWO/cssq4 - w_s(i,j,nz-kk+1)**TWO/cssq2)
    end forall
  
    ! swap pop13 and pop14 along z
    forall(i=sx:ex,j=sy:ey)
      buffservice3d(i,j,nz+kk)=-f13sub(i,j,nz+kk)
      f13sub(i,j,nz+kk)=-f14sub(i,j,nz+kk) + TWO*p(13)*rho_north* &
       (ONE+(dez(13)*w_s(i,j,nz-kk+1))**TWO/cssq4 - w_s(i,j,nz-kk+1)**TWO/cssq2)
      f14sub(i,j,nz+kk)=buffservice3d(i,j,nz+kk) + TWO*p(14)*rho_north* &
       (ONE+(dez(14)*w_s(i,j,nz-kk+1))**TWO/cssq4 - w_s(i,j,nz-kk+1)**TWO/cssq2)
    end forall
  
    ! swap pop15 and pop16 along z
    forall(i=sx:ex,j=sy:ey)
      buffservice3d(i,j,nz+kk)=-f15sub(i,j,nz+kk)
      f15sub(i,j,nz+kk)=-f16sub(i,j,nz+kk) + TWO*p(15)*rho_north* &
       (ONE+(dez(15)*w_s(i,j,nz-kk+1))**TWO/cssq4 - w_s(i,j,nz-kk+1)**TWO/cssq2)
      f16sub(i,j,nz+kk)=buffservice3d(i,j,nz+kk) + TWO*p(16)*rho_north* &
       (ONE+(dez(16)*w_s(i,j,nz-kk+1))**TWO/cssq4 - w_s(i,j,nz-kk+1)**TWO/cssq2)
    end forall
  
    ! swap pop17 and pop18 along z
    forall(i=sx:ex,j=sy:ey)
      buffservice3d(i,j,nz+kk)=-f17sub(i,j,nz+kk)
      f17sub(i,j,nz+kk)=-f18sub(i,j,nz+kk) + TWO*p(17)*rho_north* &
       (ONE+(dez(17)*w_s(i,j,nz-kk+1))**TWO/cssq4 - w_s(i,j,nz-kk+1)**TWO/cssq2)
      f18sub(i,j,nz+kk)=buffservice3d(i,j,nz+kk) + TWO*p(18)*rho_north* &
       (ONE+(dez(18)*w_s(i,j,nz-kk+1))**TWO/cssq4 - w_s(i,j,nz-kk+1)**TWO/cssq2)
    end forall
  
  case (2)
  
    ! swap pop5 and pop6 along z
    forall(i=sx:ex,j=sy:ey)
      buffservice3d(i,j,nz+kk)=f05sub(i,j,nz+kk)
      f05sub(i,j,nz+kk)=f06sub(i,j,nz+kk) + p(5)*pref_bouzidi*rho_s(i,j,nz-kk+1)* &
       (dex(5)*u_north+dey(5)*v_north+dez(5)*w_north)
      f06sub(i,j,nz+kk)=buffservice3d(i,j,nz+kk) + p(6)*pref_bouzidi*rho_s(i,j,nz-kk+1)* &
       (dex(6)*u_north+dey(6)*v_north+dez(6)*w_north)
    end forall
  
    ! swap pop11 and pop12 along z
    forall(i=sx:ex,j=sy:ey)
      buffservice3d(i,j,nz+kk)=f11sub(i,j,nz+kk)
      f11sub(i,j,nz+kk)=f12sub(i,j,nz+kk) + p(11)*pref_bouzidi*rho_s(i,j,nz-kk+1)* &
       (dex(11)*u_north+dey(11)*v_north+dez(11)*w_north)
      f12sub(i,j,nz+kk)=buffservice3d(i,j,nz+kk) + p(12)*pref_bouzidi*rho_s(i,j,nz-kk+1)* &
       (dex(12)*u_north+dey(12)*v_north+dez(12)*w_north)
    end forall
  
    ! swap pop13 and pop14 along z
    forall(i=sx:ex,j=sy:ey)
      buffservice3d(i,j,nz+kk)=f13sub(i,j,nz+kk)
      f13sub(i,j,nz+kk)=f14sub(i,j,nz+kk) + p(13)*pref_bouzidi*rho_s(i,j,nz-kk+1)* &
       (dex(13)*u_north+dey(13)*v_north+dez(13)*w_north)
      f14sub(i,j,nz+kk)=buffservice3d(i,j,nz+kk) + p(14)*pref_bouzidi*rho_s(i,j,nz-kk+1)* &
       (dex(14)*u_north+dey(14)*v_north+dez(14)*w_north)
    end forall
  
    ! swap pop15 and pop16 along z
    forall(i=sx:ex,j=sy:ey)
      buffservice3d(i,j,nz+kk)=f15sub(i,j,nz+kk)
      f15sub(i,j,nz+kk)=f16sub(i,j,nz+kk) + p(15)*pref_bouzidi*rho_s(i,j,nz-kk+1)* &
       (dex(15)*u_north+dey(15)*v_north+dez(15)*w_north)
      f16sub(i,j,nz+kk)=buffservice3d(i,j,nz+kk) + p(16)*pref_bouzidi*rho_s(i,j,nz-kk+1)* &
       (dex(16)*u_north+dey(16)*v_north+dez(16)*w_north)
    end forall
  
    ! swap pop17 and pop18 along z
    forall(i=sx:ex,j=sy:ey)
      buffservice3d(i,j,nz+kk)=f17sub(i,j,nz+kk)
      f17sub(i,j,nz+kk)=f18sub(i,j,nz+kk) + p(17)*pref_bouzidi*rho_s(i,j,nz-kk+1)* &
       (dex(17)*u_north+dey(17)*v_north+dez(17)*w_north)
      f18sub(i,j,nz+kk)=buffservice3d(i,j,nz+kk) + p(18)*pref_bouzidi*rho_s(i,j,nz-kk+1)* &
       (dex(18)*u_north+dey(18)*v_north+dez(18)*w_north)
    end forall
  
  case (3)
  
    !red fluid
    
    forall(i=sx:ex,j=sy:ey)
      f01sub(i,j,nz+kk)=equil_pop01(rho_north,u_north,v_north,w_north)
    end forall
    
    forall(i=sx:ex,j=sy:ey)
      f02sub(i,j,nz+kk)=equil_pop02(rho_north,u_north,v_north,w_north)
    end forall
    
    forall(i=sx:ex,j=sy:ey)
      f03sub(i,j,nz+kk)=equil_pop03(rho_north,u_north,v_north,w_north)
    end forall
    
    forall(i=sx:ex,j=sy:ey)
      f04sub(i,j,nz+kk)=equil_pop04(rho_north,u_north,v_north,w_north)
    end forall
    
    forall(i=sx:ex,j=sy:ey)
      f05sub(i,j,nz+kk)=equil_pop05(rho_north,u_north,v_north,w_north)
    end forall
    
    forall(i=sx:ex,j=sy:ey)
      f06sub(i,j,nz+kk)=equil_pop06(rho_north,u_north,v_north,w_north)
    end forall
    
    forall(i=sx:ex,j=sy:ey)
      f07sub(i,j,nz+kk)=equil_pop07(rho_north,u_north,v_north,w_north)
    end forall
    
    forall(i=sx:ex,j=sy:ey)
      f08sub(i,j,nz+kk)=equil_pop08(rho_north,u_north,v_north,w_north)
    end forall
    
    forall(i=sx:ex,j=sy:ey)
      f09sub(i,j,nz+kk)=equil_pop09(rho_north,u_north,v_north,w_north)
    end forall
    
    forall(i=sx:ex,j=sy:ey)
      f10sub(i,j,nz+kk)=equil_pop10(rho_north,u_north,v_north,w_north)
    end forall
    
    forall(i=sx:ex,j=sy:ey)
      f11sub(i,j,nz+kk)=equil_pop11(rho_north,u_north,v_north,w_north)
    end forall
    
    forall(i=sx:ex,j=sy:ey)
      f12sub(i,j,nz+kk)=equil_pop12(rho_north,u_north,v_north,w_north)
    end forall
    
    forall(i=sx:ex,j=sy:ey)
      f13sub(i,j,nz+kk)=equil_pop13(rho_north,u_north,v_north,w_north)
    end forall
    
    forall(i=sx:ex,j=sy:ey)
      f14sub(i,j,nz+kk)=equil_pop14(rho_north,u_north,v_north,w_north)
    end forall
    
    forall(i=sx:ex,j=sy:ey)
      f15sub(i,j,nz+kk)=equil_pop15(rho_north,u_north,v_north,w_north)
    end forall
    
    forall(i=sx:ex,j=sy:ey)
      f16sub(i,j,nz+kk)=equil_pop16(rho_north,u_north,v_north,w_north)
    end forall
    
    forall(i=sx:ex,j=sy:ey)
      f17sub(i,j,nz+kk)=equil_pop17(rho_north,u_north,v_north,w_north)
    end forall
    
    forall(i=sx:ex,j=sy:ey)
      f18sub(i,j,nz+kk)=equil_pop18(rho_north,u_north,v_north,w_north)
    end forall
    
  case default
    call error(13)
  end select
  
#else
  
  call warning(9,ZERO,latt_name)
  call error(12)
  
#endif
  
  return
  
 end subroutine apply_bounceback_north
 
 subroutine driver_bounceback_south(sx,ex,sy,ey)
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the bounceback
!     at the south side along the space defined in input sx,ex,sy,ey
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: sx,ex,sy,ey
  
  call apply_bounceback_south(sx,ex,sy,ey,bc_rhoR_south,bc_u_south, &
   bc_v_south,bc_w_south,rhoR,u,v,w,f00R,f01R,f02R,f03R,f04R,f05R,f06R, &
   f07R,f08R,f09R,f10R,f11R,f12R,f13R,f14R,f15R,f16R,f17R,f18R)
  
  if(lsingle_fluid)return
 
  call apply_bounceback_south(sx,ex,sy,ey,bc_rhoB_south,bc_u_south, &
   bc_v_south,bc_w_south,rhoB,u,v,w,f00B,f01B,f02B,f03B,f04B,f05B,f06B, &
   f07B,f08B,f09B,f10B,f11B,f12B,f13B,f14B,f15B,f16B,f17B,f18B)
   
  return
  
 end subroutine driver_bounceback_south
 
 subroutine apply_bounceback_south(sx,ex,sy,ey,rho_south,u_south, &
  v_south,w_south,rho_s,u_s,v_s,w_s,f00sub,f01sub,f02sub,f03sub,f04sub, &
  f05sub,f06sub,f07sub,f08sub,f09sub,f10sub,f11sub,f12sub,f13sub, &
  f14sub,f15sub,f16sub,f17sub,f18sub)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the bounceback 
!     at the south side along the space defined in input sx,ex,sy,ey
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: sx,ex,sy,ey
  real(kind=PRC), intent(in) :: rho_south,u_south,v_south,w_south
  real(kind=PRC), allocatable, dimension(:,:,:)  :: rho_s,u_s,v_s,w_s, &
   f00sub,f01sub,f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub, &
   f09sub,f10sub,f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub
  
  integer :: i,j,k
  integer, parameter :: kk = 1
  
  real(kind=PRC), parameter :: pref_bouzidi=TWO/cssq
  
  real(kind=PRC), parameter :: cssq2 = ( HALF / cssq )
  real(kind=PRC), parameter :: cssq4 = ( HALF / (cssq*cssq) )
  
#if LATTICE==319
  
  select case (bc_type_south)  
  case (0)
  
    ! swap pop5 and pop6 along z
    forall(i=sx:ex,j=sy:ey)
      buffservice3d(i,j,1-kk)=f05sub(i,j,1-kk)
      f05sub(i,j,1-kk)=f06sub(i,j,1-kk)
      f06sub(i,j,1-kk)=buffservice3d(i,j,1-kk)
    end forall
    
    ! swap pop11 and pop12 along z
    forall(i=sx:ex,j=sy:ey)
      buffservice3d(i,j,1-kk)=f11sub(i,j,1-kk)
      f11sub(i,j,1-kk)=f12sub(i,j,1-kk)
      f12sub(i,j,1-kk)=buffservice3d(i,j,1-kk)
    end forall
    
    ! swap pop13 and pop14 along z
    forall(i=sx:ex,j=sy:ey)
      buffservice3d(i,j,1-kk)=f13sub(i,j,1-kk)
      f13sub(i,j,1-kk)=f14sub(i,j,1-kk)
      f14sub(i,j,1-kk)=buffservice3d(i,j,1-kk)
    end forall
    
    ! swap pop15 and pop16 along z
    forall(i=sx:ex,j=sy:ey)
      buffservice3d(i,j,1-kk)=f15sub(i,j,1-kk)
      f15sub(i,j,1-kk)=f16sub(i,j,1-kk)
      f16sub(i,j,1-kk)=buffservice3d(i,j,1-kk)
    end forall
    
    ! swap pop17 and pop18 along z
    forall(i=sx:ex,j=sy:ey)
      buffservice3d(i,j,1-kk)=f17sub(i,j,1-kk)
      f17sub(i,j,1-kk)=f18sub(i,j,1-kk)
      f18sub(i,j,1-kk)=buffservice3d(i,j,1-kk)
    end forall
    
  case (1)
  
    ! swap pop5 and pop6 along z
    forall(i=sx:ex,j=sy:ey)
      buffservice3d(i,j,1-kk)=-f05sub(i,j,1-kk)
      f05sub(i,j,1-kk)=-f06sub(i,j,1-kk) + TWO*p(5)*rho_south* &
       (ONE+(dez(5)*w_s(i,j,kk))**TWO/cssq4 - w_s(i,j,kk)**TWO/cssq2)
      f06sub(i,j,1-kk)=buffservice3d(i,j,1-kk) + TWO*p(6)*rho_south* &
       (ONE+(dez(6)*w_s(i,j,kk))**TWO/cssq4 - w_s(i,j,kk)**TWO/cssq2)
    end forall
    
    ! swap pop11 and pop12 along z
    forall(i=sx:ex,j=sy:ey)
      buffservice3d(i,j,1-kk)=-f11sub(i,j,1-kk)
      f11sub(i,j,1-kk)=-f12sub(i,j,1-kk) + TWO*p(11)*rho_south* &
       (ONE+(dez(11)*w_s(i,j,kk))**TWO/cssq4 - w_s(i,j,kk)**TWO/cssq2)
      f12sub(i,j,1-kk)=buffservice3d(i,j,1-kk) + TWO*p(12)*rho_south* &
       (ONE+(dez(12)*w_s(i,j,kk))**TWO/cssq4 - w_s(i,j,kk)**TWO/cssq2)
    end forall
    
    ! swap pop13 and pop14 along z
    forall(i=sx:ex,j=sy:ey)
      buffservice3d(i,j,1-kk)=-f13sub(i,j,1-kk)
      f13sub(i,j,1-kk)=-f14sub(i,j,1-kk) + TWO*p(13)*rho_south* &
       (ONE+(dez(13)*w_s(i,j,kk))**TWO/cssq4 - w_s(i,j,kk)**TWO/cssq2)
      f14sub(i,j,1-kk)=buffservice3d(i,j,1-kk) + TWO*p(14)*rho_south* &
       (ONE+(dez(14)*w_s(i,j,kk))**TWO/cssq4 - w_s(i,j,kk)**TWO/cssq2)
    end forall
    
    ! swap pop15 and pop16 along z
    forall(i=sx:ex,j=sy:ey)
      buffservice3d(i,j,1-kk)=-f15sub(i,j,1-kk)
      f15sub(i,j,1-kk)=-f16sub(i,j,1-kk) + TWO*p(15)*rho_south* &
       (ONE+(dez(15)*w_s(i,j,kk))**TWO/cssq4 - w_s(i,j,kk)**TWO/cssq2)
      f16sub(i,j,1-kk)=buffservice3d(i,j,1-kk) + TWO*p(16)*rho_south* &
       (ONE+(dez(16)*w_s(i,j,kk))**TWO/cssq4 - w_s(i,j,kk)**TWO/cssq2)
    end forall
    
    ! swap pop17 and pop18 along z
    forall(i=sx:ex,j=sy:ey)
      buffservice3d(i,j,1-kk)=-f17sub(i,j,1-kk)
      f17sub(i,j,1-kk)=-f18sub(i,j,1-kk) + TWO*p(17)*rho_south* &
       (ONE+(dez(17)*w_s(i,j,kk))**TWO/cssq4 - w_s(i,j,kk)**TWO/cssq2)
      f18sub(i,j,1-kk)=buffservice3d(i,j,1-kk) + TWO*p(18)*rho_south* &
       (ONE+(dez(18)*w_s(i,j,kk))**TWO/cssq4 - w_s(i,j,kk)**TWO/cssq2)
    end forall
    
  case (2)
  
    ! swap pop5 and pop6 along z
    forall(i=sx:ex,j=sy:ey)
      buffservice3d(i,j,1-kk)=f05sub(i,j,1-kk)
      f05sub(i,j,1-kk)=f06sub(i,j,1-kk) + p(5)*pref_bouzidi*rho_s(i,j,kk)* &
       (dex(5)*u_south+dey(5)*v_south+dez(5)*w_south)
      f06sub(i,j,1-kk)=buffservice3d(i,j,1-kk) + p(6)*pref_bouzidi*rho_s(i,j,kk)* &
       (dex(6)*u_south+dey(6)*v_south+dez(6)*w_south)
    end forall
    
    ! swap pop11 and pop12 along z
    forall(i=sx:ex,j=sy:ey)
      buffservice3d(i,j,1-kk)=f11sub(i,j,1-kk)
      f11sub(i,j,1-kk)=f12sub(i,j,1-kk) + p(11)*pref_bouzidi*rho_s(i,j,kk)* &
       (dex(11)*u_south+dey(11)*v_south+dez(11)*w_south)
      f12sub(i,j,1-kk)=buffservice3d(i,j,1-kk) + p(12)*pref_bouzidi*rho_s(i,j,kk)* &
       (dex(12)*u_south+dey(12)*v_south+dez(12)*w_south)
    end forall
    
    ! swap pop13 and pop14 along z
    forall(i=sx:ex,j=sy:ey)
      buffservice3d(i,j,1-kk)=f13sub(i,j,1-kk)
      f13sub(i,j,1-kk)=f14sub(i,j,1-kk) + p(13)*pref_bouzidi*rho_s(i,j,kk)* &
       (dex(13)*u_south+dey(13)*v_south+dez(13)*w_south)
      f14sub(i,j,1-kk)=buffservice3d(i,j,1-kk) + p(14)*pref_bouzidi*rho_s(i,j,kk)* &
       (dex(14)*u_south+dey(14)*v_south+dez(14)*w_south)
    end forall
    
    ! swap pop15 and pop16 along z
    forall(i=sx:ex,j=sy:ey)
      buffservice3d(i,j,1-kk)=f15sub(i,j,1-kk)
      f15sub(i,j,1-kk)=f16sub(i,j,1-kk) + p(15)*pref_bouzidi*rho_s(i,j,kk)* &
       (dex(15)*u_south+dey(15)*v_south+dez(15)*w_south)
      f16sub(i,j,1-kk)=buffservice3d(i,j,1-kk) + p(16)*pref_bouzidi*rho_s(i,j,kk)* &
       (dex(16)*u_south+dey(16)*v_south+dez(16)*w_south)
    end forall
    
    ! swap pop17 and pop18 along z
    forall(i=sx:ex,j=sy:ey)
      buffservice3d(i,j,1-kk)=f17sub(i,j,1-kk)
      f17sub(i,j,1-kk)=f18sub(i,j,1-kk) + p(17)*pref_bouzidi*rho_s(i,j,kk)* &
       (dex(17)*u_south+dey(17)*v_south+dez(17)*w_south)
      f18sub(i,j,1-kk)=buffservice3d(i,j,1-kk) + p(18)*pref_bouzidi*rho_s(i,j,kk)* &
       (dex(18)*u_south+dey(18)*v_south+dez(18)*w_south)
    end forall
    
  case (3)
    
    !red fluid
    forall(i=sx:ex,j=sy:ey)
      f00sub(i,j,1-kk)=equil_pop00(rho_south,u_south,v_south,w_south)
    end forall
    
    forall(i=sx:ex,j=sy:ey)
      f01sub(i,j,1-kk)=equil_pop01(rho_south,u_south,v_south,w_south)
    end forall
    
    forall(i=sx:ex,j=sy:ey)
      f02sub(i,j,1-kk)=equil_pop02(rho_south,u_south,v_south,w_south)
    end forall
    
    forall(i=sx:ex,j=sy:ey)
      f03sub(i,j,1-kk)=equil_pop03(rho_south,u_south,v_south,w_south)
    end forall
    
    forall(i=sx:ex,j=sy:ey)
      f04sub(i,j,1-kk)=equil_pop04(rho_south,u_south,v_south,w_south)
    end forall
    
    forall(i=sx:ex,j=sy:ey)
      f05sub(i,j,1-kk)=equil_pop05(rho_south,u_south,v_south,w_south)
    end forall
    
    forall(i=sx:ex,j=sy:ey)
      f06sub(i,j,1-kk)=equil_pop06(rho_south,u_south,v_south,w_south)
    end forall
    
    forall(i=sx:ex,j=sy:ey)
      f07sub(i,j,1-kk)=equil_pop07(rho_south,u_south,v_south,w_south)
    end forall
    
    forall(i=sx:ex,j=sy:ey)
      f08sub(i,j,1-kk)=equil_pop08(rho_south,u_south,v_south,w_south)
    end forall
    
    forall(i=sx:ex,j=sy:ey)
      f09sub(i,j,1-kk)=equil_pop09(rho_south,u_south,v_south,w_south)
    end forall
    
    forall(i=sx:ex,j=sy:ey)
      f10sub(i,j,1-kk)=equil_pop10(rho_south,u_south,v_south,w_south)
    end forall
    
    forall(i=sx:ex,j=sy:ey)
      f11sub(i,j,1-kk)=equil_pop11(rho_south,u_south,v_south,w_south)
    end forall
    
    forall(i=sx:ex,j=sy:ey)
      f12sub(i,j,1-kk)=equil_pop12(rho_south,u_south,v_south,w_south)
    end forall
    
    forall(i=sx:ex,j=sy:ey)
      f13sub(i,j,1-kk)=equil_pop13(rho_south,u_south,v_south,w_south)
    end forall
    
    forall(i=sx:ex,j=sy:ey)
      f14sub(i,j,1-kk)=equil_pop14(rho_south,u_south,v_south,w_south)
    end forall
    
    forall(i=sx:ex,j=sy:ey)
      f15sub(i,j,1-kk)=equil_pop15(rho_south,u_south,v_south,w_south)
    end forall
    
    forall(i=sx:ex,j=sy:ey)
      f16sub(i,j,1-kk)=equil_pop16(rho_south,u_south,v_south,w_south)
    end forall
    
    forall(i=sx:ex,j=sy:ey)
      f17sub(i,j,1-kk)=equil_pop17(rho_south,u_south,v_south,w_south)
    end forall
    
    forall(i=sx:ex,j=sy:ey)
      f18sub(i,j,1-kk)=equil_pop18(rho_south,u_south,v_south,w_south)
    end forall
  
  case default
    call error(13)
  end select
  
#else
  
  call warning(9,ZERO,latt_name)
  call error(12)
  
#endif
  
  return
  
 end subroutine apply_bounceback_south
 
 subroutine driver_bounceback_front(sx,ex,sz,ez)
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the bounceback
!     at the front side along the space defined in input sx,ex,sz,ez
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: sx,ex,sz,ez
  
  call apply_bounceback_front(sx,ex,sz,ez,bc_rhoR_front,bc_u_front, &
   bc_v_front,bc_w_front,rhoR,u,v,w,f00R,f01R,f02R,f03R,f04R,f05R,f06R, &
   f07R,f08R,f09R,f10R,f11R,f12R,f13R,f14R,f15R,f16R,f17R,f18R)
  
  if(lsingle_fluid)return
 
  call apply_bounceback_front(sx,ex,sz,ez,bc_rhoB_front,bc_u_front, &
   bc_v_front,bc_w_front,rhoB,u,v,w,f00B,f01B,f02B,f03B,f04B,f05B,f06B, &
   f07B,f08B,f09B,f10B,f11B,f12B,f13B,f14B,f15B,f16B,f17B,f18B)
   
  return
  
 end subroutine driver_bounceback_front
 
 subroutine apply_bounceback_front(sx,ex,sz,ez,rho_front,u_front, &
  v_front,w_front,rho_s,u_s,v_s,w_s,f00sub,f01sub,f02sub,f03sub,f04sub, &
  f05sub,f06sub,f07sub,f08sub,f09sub,f10sub,f11sub,f12sub,f13sub, &
  f14sub,f15sub,f16sub,f17sub,f18sub)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the bounceback
!     at the front side along the space defined in input sx,ex,sz,ez
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: sx,ex,sz,ez
  real(kind=PRC), intent(in) :: rho_front,u_front,v_front,w_front
  real(kind=PRC), allocatable, dimension(:,:,:)  :: rho_s,u_s,v_s,w_s, &
   f00sub,f01sub,f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub, &
   f09sub,f10sub,f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub
  
  integer :: i,j,k
  integer, parameter :: kk = 1
  
  real(kind=PRC), parameter :: pref_bouzidi=TWO/cssq
  
  real(kind=PRC), parameter :: cssq2 = ( HALF / cssq )
  real(kind=PRC), parameter :: cssq4 = ( HALF / (cssq*cssq) )
  
#if LATTICE==319
  
  select case (bc_type_front)  
  case (0)
  
    ! swap pop3 and pop4 along y
    forall(i=sx:ex,k=sz:ez)
      buffservice3d(i,1-kk,k)=f03sub(i,1-kk,k)
      f03sub(i,1-kk,k)=f04sub(i,1-kk,k)
      f04sub(i,1-kk,k)=buffservice3d(i,1-kk,k)
    end forall
    
    ! swap pop7 and pop8 along y
    forall(i=sx:ex,k=sz:ez)
      buffservice3d(i,1-kk,k)=f07sub(i,1-kk,k)
      f07sub(i,1-kk,k)=f08sub(i,1-kk,k)
      f08sub(i,1-kk,k)=buffservice3d(i,1-kk,k)
    end forall
    
    ! swap pop9 and pop10 along y
    forall(i=sx:ex,k=sz:ez)
      buffservice3d(i,1-kk,k)=f09sub(i,1-kk,k)
      f09sub(i,1-kk,k)=f10sub(i,1-kk,k)
      f10sub(i,1-kk,k)=buffservice3d(i,1-kk,k)
    end forall
    
    ! swap pop15 and pop16 along y
    forall(i=sx:ex,k=sz:ez)
      buffservice3d(i,1-kk,k)=f15sub(i,1-kk,k)
      f15sub(i,1-kk,k)=f16sub(i,1-kk,k)
      f16sub(i,1-kk,k)=buffservice3d(i,1-kk,k)
    end forall
    
    ! swap pop17 and pop18 along y
    forall(i=sx:ex,k=sz:ez)
      buffservice3d(i,1-kk,k)=f17sub(i,1-kk,k)
      f17sub(i,1-kk,k)=f18sub(i,1-kk,k)
      f18sub(i,1-kk,k)=buffservice3d(i,1-kk,k)
    end forall
    
  case (1)
  
    ! swap pop3 and pop4 along y
    forall(i=sx:ex,k=sz:ez)
      buffservice3d(i,1-kk,k)=-f03sub(i,1-kk,k)
      f03sub(i,1-kk,k)=-f04sub(i,1-kk,k) + TWO*p(3)*rho_front* &
       (ONE+(dey(3)*v_s(i,kk,k))**TWO/cssq4 - v_s(i,kk,k)**TWO/cssq2)
      f04sub(i,1-kk,k)=buffservice3d(i,1-kk,k) + TWO*p(4)*rho_front* &
       (ONE+(dey(4)*v_s(i,kk,k))**TWO/cssq4 - v_s(i,kk,k)**TWO/cssq2)
    end forall
    
    ! swap pop7 and pop8 along y
    forall(i=sx:ex,k=sz:ez)
      buffservice3d(i,1-kk,k)=-f07sub(i,1-kk,k)
      f07sub(i,1-kk,k)=-f08sub(i,1-kk,k) + TWO*p(7)*rho_front* &
       (ONE+(dey(7)*v_s(i,kk,k))**TWO/cssq4 - v_s(i,kk,k)**TWO/cssq2)
      f08sub(i,1-kk,k)=buffservice3d(i,1-kk,k) + TWO*p(8)*rho_front* &
       (ONE+(dey(8)*v_s(i,kk,k))**TWO/cssq4 - v_s(i,kk,k)**TWO/cssq2)
    end forall
    
    ! swap pop9 and pop10 along y
    forall(i=sx:ex,k=sz:ez)
      buffservice3d(i,1-kk,k)=-f09sub(i,1-kk,k)
      f09sub(i,1-kk,k)=-f10sub(i,1-kk,k) + TWO*p(9)*rho_front* &
       (ONE+(dey(9)*v_s(i,kk,k))**TWO/cssq4 - v_s(i,kk,k)**TWO/cssq2)
      f10sub(i,1-kk,k)=buffservice3d(i,1-kk,k) + TWO*p(10)*rho_front* &
       (ONE+(dey(10)*v_s(i,kk,k))**TWO/cssq4 - v_s(i,kk,k)**TWO/cssq2)
    end forall
    
    ! swap pop15 and pop16 along y
    forall(i=sx:ex,k=sz:ez)
      buffservice3d(i,1-kk,k)=-f15sub(i,1-kk,k)
      f15sub(i,1-kk,k)=-f16sub(i,1-kk,k) + TWO*p(15)*rho_front* &
       (ONE+(dey(15)*v_s(i,kk,k))**TWO/cssq4 - v_s(i,kk,k)**TWO/cssq2)
      f16sub(i,1-kk,k)=buffservice3d(i,1-kk,k) + TWO*p(16)*rho_front* &
       (ONE+(dey(16)*v_s(i,kk,k))**TWO/cssq4 - v_s(i,kk,k)**TWO/cssq2)
    end forall
    
    ! swap pop17 and pop18 along y
    forall(i=sx:ex,k=sz:ez)
      buffservice3d(i,1-kk,k)=-f17sub(i,1-kk,k)
      f17sub(i,1-kk,k)=-f18sub(i,1-kk,k) + TWO*p(17)*rho_front* &
       (ONE+(dey(17)*v_s(i,kk,k))**TWO/cssq4 - v_s(i,kk,k)**TWO/cssq2)
      f18sub(i,1-kk,k)=buffservice3d(i,1-kk,k) + TWO*p(18)*rho_front* &
       (ONE+(dey(18)*v_s(i,kk,k))**TWO/cssq4 - v_s(i,kk,k)**TWO/cssq2)
    end forall
  
  case (2)
  
    ! swap pop3 and pop4 along y
    forall(i=sx:ex,k=sz:ez)
      buffservice3d(i,1-kk,k)=f03sub(i,1-kk,k)
      f03sub(i,1-kk,k)=f04sub(i,1-kk,k) + p(3)*pref_bouzidi*rho_s(i,kk,k)* &
       (dex(3)*u_front+dey(3)*v_front+dez(3)*w_front)
      f04sub(i,1-kk,k)=buffservice3d(i,1-kk,k) + p(4)*pref_bouzidi*rho_s(i,kk,k)* &
       (dex(4)*u_front+dey(4)*v_front+dez(4)*w_front)
    end forall
    
    ! swap pop7 and pop8 along y
    forall(i=sx:ex,k=sz:ez)
      buffservice3d(i,1-kk,k)=f07sub(i,1-kk,k)
      f07sub(i,1-kk,k)=f08sub(i,1-kk,k) + p(7)*pref_bouzidi*rho_s(i,kk,k)* &
       (dex(7)*u_front+dey(7)*v_front+dez(7)*w_front)
      f08sub(i,1-kk,k)=buffservice3d(i,1-kk,k) + p(8)*pref_bouzidi*rho_s(i,kk,k)* &
       (dex(8)*u_front+dey(8)*v_front+dez(8)*w_front)
    end forall
    
    ! swap pop9 and pop10 along y
    forall(i=sx:ex,k=sz:ez)
      buffservice3d(i,1-kk,k)=f09sub(i,1-kk,k)
      f09sub(i,1-kk,k)=f10sub(i,1-kk,k) + p(9)*pref_bouzidi*rho_s(i,kk,k)* &
       (dex(9)*u_front+dey(9)*v_front+dez(9)*w_front)
      f10sub(i,1-kk,k)=buffservice3d(i,1-kk,k) + p(10)*pref_bouzidi*rho_s(i,kk,k)* &
       (dex(10)*u_front+dey(10)*v_front+dez(10)*w_front)
    end forall
    
    ! swap pop15 and pop16 along y
    forall(i=sx:ex,k=sz:ez)
      buffservice3d(i,1-kk,k)=f15sub(i,1-kk,k)
      f15sub(i,1-kk,k)=f16sub(i,1-kk,k) + p(15)*pref_bouzidi*rho_s(i,kk,k)* &
       (dex(15)*u_front+dey(15)*v_front+dez(15)*w_front)
      f16sub(i,1-kk,k)=buffservice3d(i,1-kk,k) + p(16)*pref_bouzidi*rho_s(i,kk,k)* &
       (dex(16)*u_front+dey(16)*v_front+dez(16)*w_front)
    end forall
    
    ! swap pop17 and pop18 along y
    forall(i=sx:ex,k=sz:ez)
      buffservice3d(i,1-kk,k)=f17sub(i,1-kk,k)
      f17sub(i,1-kk,k)=f18sub(i,1-kk,k) + p(17)*pref_bouzidi*rho_s(i,kk,k)* &
       (dex(17)*u_front+dey(17)*v_front+dez(17)*w_front)
      f18sub(i,1-kk,k)=buffservice3d(i,1-kk,k) + p(18)*pref_bouzidi*rho_s(i,kk,k)* &
       (dex(18)*u_front+dey(18)*v_front+dez(18)*w_front)
    end forall
  
  case (3)
    
    !red fluid
    forall(i=sx:ex,k=sz:ez)
      f00sub(i,1-kk,k)=equil_pop00(rho_front,u_front,v_front,w_front)
    end forall
    
    forall(i=sx:ex,k=sz:ez)
      f01sub(i,1-kk,k)=equil_pop01(rho_front,u_front,v_front,w_front)
    end forall
    
    forall(i=sx:ex,k=sz:ez)
      f02sub(i,1-kk,k)=equil_pop02(rho_front,u_front,v_front,w_front)
    end forall
    
    forall(i=sx:ex,k=sz:ez)
      f03sub(i,1-kk,k)=equil_pop03(rho_front,u_front,v_front,w_front)
    end forall
    
    forall(i=sx:ex,k=sz:ez)
      f04sub(i,1-kk,k)=equil_pop04(rho_front,u_front,v_front,w_front)
    end forall
    
    forall(i=sx:ex,k=sz:ez)
      f05sub(i,1-kk,k)=equil_pop05(rho_front,u_front,v_front,w_front)
    end forall
    
    forall(i=sx:ex,k=sz:ez)
      f06sub(i,1-kk,k)=equil_pop06(rho_front,u_front,v_front,w_front)
    end forall
    
    forall(i=sx:ex,k=sz:ez)
      f07sub(i,1-kk,k)=equil_pop07(rho_front,u_front,v_front,w_front)
    end forall
    
    forall(i=sx:ex,k=sz:ez)
      f08sub(i,1-kk,k)=equil_pop08(rho_front,u_front,v_front,w_front)
    end forall
    
    forall(i=sx:ex,k=sz:ez)
      f09sub(i,1-kk,k)=equil_pop09(rho_front,u_front,v_front,w_front)
    end forall
    
    forall(i=sx:ex,k=sz:ez)
      f10sub(i,1-kk,k)=equil_pop10(rho_front,u_front,v_front,w_front)
    end forall
    
    forall(i=sx:ex,k=sz:ez)
      f11sub(i,1-kk,k)=equil_pop11(rho_front,u_front,v_front,w_front)
    end forall
    
    forall(i=sx:ex,k=sz:ez)
      f12sub(i,1-kk,k)=equil_pop12(rho_front,u_front,v_front,w_front)
    end forall
    
    forall(i=sx:ex,k=sz:ez)
      f13sub(i,1-kk,k)=equil_pop13(rho_front,u_front,v_front,w_front)
    end forall
    
    forall(i=sx:ex,k=sz:ez)
      f14sub(i,1-kk,k)=equil_pop14(rho_front,u_front,v_front,w_front)
    end forall
    
    forall(i=sx:ex,k=sz:ez)
      f15sub(i,1-kk,k)=equil_pop15(rho_front,u_front,v_front,w_front)
    end forall
    
    forall(i=sx:ex,k=sz:ez)
      f16sub(i,1-kk,k)=equil_pop16(rho_front,u_front,v_front,w_front)
    end forall
    
    forall(i=sx:ex,k=sz:ez)
      f17sub(i,1-kk,k)=equil_pop17(rho_front,u_front,v_front,w_front)
    end forall
    
    forall(i=sx:ex,k=sz:ez)
      f18sub(i,1-kk,k)=equil_pop18(rho_front,u_front,v_front,w_front)
    end forall
    
  case default
    call error(13)
  end select
  
#else
  
  call warning(9,ZERO,latt_name)
  call error(12)
  
#endif
  
  return
  
 end subroutine apply_bounceback_front
 
 subroutine driver_bounceback_rear(sx,ex,sz,ez)
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the bounceback
!     at the rear side along the space defined in input sx,ex,sz,ez
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: sx,ex,sz,ez
  
  call apply_bounceback_rear(sx,ex,sz,ez,bc_rhoR_rear,bc_u_rear, &
   bc_v_rear,bc_w_rear,rhoR,u,v,w,f00R,f01R,f02R,f03R,f04R,f05R,f06R, &
   f07R,f08R,f09R,f10R,f11R,f12R,f13R,f14R,f15R,f16R,f17R,f18R)
  
  if(lsingle_fluid)return
 
  call apply_bounceback_rear(sx,ex,sz,ez,bc_rhoB_rear,bc_u_rear, &
   bc_v_rear,bc_w_rear,rhoB,u,v,w,f00B,f01B,f02B,f03B,f04B,f05B,f06B, &
   f07B,f08B,f09B,f10B,f11B,f12B,f13B,f14B,f15B,f16B,f17B,f18B)
   
  return
  
 end subroutine driver_bounceback_rear
 
 subroutine apply_bounceback_rear(sx,ex,sz,ez,rho_rear,u_rear, &
  v_rear,w_rear,rho_s,u_s,v_s,w_s,f00sub,f01sub,f02sub,f03sub,f04sub, &
  f05sub,f06sub,f07sub,f08sub,f09sub,f10sub,f11sub,f12sub,f13sub, &
  f14sub,f15sub,f16sub,f17sub,f18sub)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the bounceback 
!     at the rear side along the space defined in input sx,ex,sz,ez
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: sx,ex,sz,ez
  real(kind=PRC), intent(in) :: rho_rear,u_rear,v_rear,w_rear
  real(kind=PRC), allocatable, dimension(:,:,:)  :: rho_s,u_s,v_s,w_s, &
   f00sub,f01sub,f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub, &
   f09sub,f10sub,f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub
  
  integer :: i,j,k
  integer, parameter :: kk = 1
  
  real(kind=PRC), parameter :: pref_bouzidi=TWO/cssq
  
  real(kind=PRC), parameter :: cssq2 = ( HALF / cssq )
  real(kind=PRC), parameter :: cssq4 = ( HALF / (cssq*cssq) )
  
#if LATTICE==319
  
  select case (bc_type_rear)  
  case (0)
  
    ! swap pop3 and pop4 along y
    forall(i=sx:ex,k=sz:ez)
      buffservice3d(i,ny+kk,k)=f03sub(i,ny+kk,k)
      f03sub(i,ny+kk,k)=f04sub(i,ny+kk,k)
      f04sub(i,ny+kk,k)=buffservice3d(i,ny+kk,k)
    end forall
    
    ! swap pop7 and pop8 along y
    forall(i=sx:ex,k=sz:ez)
      buffservice3d(i,ny+kk,k)=f07sub(i,ny+kk,k)
      f07sub(i,ny+kk,k)=f08sub(i,ny+kk,k)
      f08sub(i,ny+kk,k)=buffservice3d(i,ny+kk,k)
    end forall
    
    ! swap pop9 and pop10 along y
    forall(i=sx:ex,k=sz:ez)
      buffservice3d(i,ny+kk,k)=f09sub(i,ny+kk,k)
      f09sub(i,ny+kk,k)=f10sub(i,ny+kk,k)
      f10sub(i,ny+kk,k)=buffservice3d(i,ny+kk,k)
    end forall
    
    ! swap pop15 and pop16 along y
    forall(i=sx:ex,k=sz:ez)
      buffservice3d(i,ny+kk,k)=f15sub(i,ny+kk,k)
      f15sub(i,ny+kk,k)=f16sub(i,ny+kk,k)
      f16sub(i,ny+kk,k)=buffservice3d(i,ny+kk,k)
    end forall
    
    ! swap pop17 and pop18 along y
    forall(i=sx:ex,k=sz:ez)
      buffservice3d(i,ny+kk,k)=f17sub(i,ny+kk,k)
      f17sub(i,ny+kk,k)=f18sub(i,ny+kk,k)
      f18sub(i,ny+kk,k)=buffservice3d(i,ny+kk,k)
    end forall
    
  case (1)
  
    ! swap pop3 and pop4 along y
    forall(i=sx:ex,k=sz:ez)
      buffservice3d(i,ny+kk,k)=-f03sub(i,ny+kk,k)
      f03sub(i,ny+kk,k)=-f04sub(i,ny+kk,k) + TWO*p(3)*rho_rear* &
       (ONE+(dey(3)*v_s(i,ny-kk+1,k))**TWO/cssq4 - v_s(i,ny-kk+1,k)**TWO/cssq2)
      f04sub(i,ny+kk,k)=buffservice3d(i,ny+kk,k) + TWO*p(4)*rho_rear* &
       (ONE+(dey(4)*v_s(i,ny-kk+1,k))**TWO/cssq4 - v_s(i,ny-kk+1,k)**TWO/cssq2)
    end forall
    
    ! swap pop7 and pop8 along y
    forall(i=sx:ex,k=sz:ez)
      buffservice3d(i,ny+kk,k)=-f07sub(i,ny+kk,k)
      f07sub(i,ny+kk,k)=-f08sub(i,ny+kk,k) + TWO*p(7)*rho_rear* &
       (ONE+(dey(7)*v_s(i,ny-kk+1,k))**TWO/cssq4 - v_s(i,ny-kk+1,k)**TWO/cssq2)
      f08sub(i,ny+kk,k)=buffservice3d(i,ny+kk,k) + TWO*p(8)*rho_rear* &
       (ONE+(dey(8)*v_s(i,ny-kk+1,k))**TWO/cssq4 - v_s(i,ny-kk+1,k)**TWO/cssq2)
    end forall
    
    ! swap pop9 and pop10 along y
    forall(i=sx:ex,k=sz:ez)
      buffservice3d(i,ny+kk,k)=-f09sub(i,ny+kk,k)
      f09sub(i,ny+kk,k)=-f10sub(i,ny+kk,k) + TWO*p(9)*rho_rear* &
       (ONE+(dey(9)*v_s(i,ny-kk+1,k))**TWO/cssq4 - v_s(i,ny-kk+1,k)**TWO/cssq2)
      f10sub(i,ny+kk,k)=buffservice3d(i,ny+kk,k) + TWO*p(10)*rho_rear* &
       (ONE+(dey(10)*v_s(i,ny-kk+1,k))**TWO/cssq4 - v_s(i,ny-kk+1,k)**TWO/cssq2)
    end forall
    
    ! swap pop15 and pop16 along y
    forall(i=sx:ex,k=sz:ez)
      buffservice3d(i,ny+kk,k)=-f15sub(i,ny+kk,k)
      f15sub(i,ny+kk,k)=-f16sub(i,ny+kk,k) + TWO*p(15)*rho_rear* &
       (ONE+(dey(15)*v_s(i,ny-kk+1,k))**TWO/cssq4 - v_s(i,ny-kk+1,k)**TWO/cssq2)
      f16sub(i,ny+kk,k)=buffservice3d(i,ny+kk,k) + TWO*p(16)*rho_rear* &
       (ONE+(dey(16)*v_s(i,ny-kk+1,k))**TWO/cssq4 - v_s(i,ny-kk+1,k)**TWO/cssq2)
    end forall
    
    ! swap pop17 and pop18 along y
    forall(i=sx:ex,k=sz:ez)
      buffservice3d(i,ny+kk,k)=-f17sub(i,ny+kk,k)
      f17sub(i,ny+kk,k)=-f18sub(i,ny+kk,k) + TWO*p(17)*rho_rear* &
       (ONE+(dey(17)*v_s(i,ny-kk+1,k))**TWO/cssq4 - v_s(i,ny-kk+1,k)**TWO/cssq2)
      f18sub(i,ny+kk,k)=buffservice3d(i,ny+kk,k) + TWO*p(18)*rho_rear* &
       (ONE+(dey(18)*v_s(i,ny-kk+1,k))**TWO/cssq4 - v_s(i,ny-kk+1,k)**TWO/cssq2)
    end forall
    
  case (2)
  
    ! swap pop3 and pop4 along y
    forall(i=sx:ex,k=sz:ez)
      buffservice3d(i,ny+kk,k)=f03sub(i,ny+kk,k)
      f03sub(i,ny+kk,k)=f04sub(i,ny+kk,k) + p(3)*pref_bouzidi*rho_s(i,ny-kk+1,k)* &
       (dex(3)*u_rear+dey(3)*v_rear+dez(3)*w_rear)
      f04sub(i,ny+kk,k)=buffservice3d(i,ny+kk,k) + p(4)*pref_bouzidi*rho_s(i,ny-kk+1,k)* &
       (dex(4)*u_rear+dey(4)*v_rear+dez(4)*w_rear)
    end forall
    
    ! swap pop7 and pop8 along y
    forall(i=sx:ex,k=sz:ez)
      buffservice3d(i,ny+kk,k)=f07sub(i,ny+kk,k)
      f07sub(i,ny+kk,k)=f08sub(i,ny+kk,k) + p(7)*pref_bouzidi*rho_s(i,ny-kk+1,k)* &
       (dex(7)*u_rear+dey(7)*v_rear+dez(7)*w_rear)
      f08sub(i,ny+kk,k)=buffservice3d(i,ny+kk,k) + p(8)*pref_bouzidi*rho_s(i,ny-kk+1,k)* &
       (dex(8)*u_rear+dey(8)*v_rear+dez(8)*w_rear)
    end forall
    
    ! swap pop9 and pop10 along y
    forall(i=sx:ex,k=sz:ez)
      buffservice3d(i,ny+kk,k)=f09sub(i,ny+kk,k)
      f09sub(i,ny+kk,k)=f10sub(i,ny+kk,k) + p(9)*pref_bouzidi*rho_s(i,ny-kk+1,k)* &
       (dex(9)*u_rear+dey(9)*v_rear+dez(9)*w_rear)
      f10sub(i,ny+kk,k)=buffservice3d(i,ny+kk,k) + p(10)*pref_bouzidi*rho_s(i,ny-kk+1,k)* &
       (dex(10)*u_rear+dey(10)*v_rear+dez(10)*w_rear)
    end forall
    
    ! swap pop15 and pop16 along y
    forall(i=sx:ex,k=sz:ez)
      buffservice3d(i,ny+kk,k)=f15sub(i,ny+kk,k)
      f15sub(i,ny+kk,k)=f16sub(i,ny+kk,k) + p(15)*pref_bouzidi*rho_s(i,ny-kk+1,k)* &
       (dex(15)*u_rear+dey(15)*v_rear+dez(15)*w_rear)
      f16sub(i,ny+kk,k)=buffservice3d(i,ny+kk,k) + p(16)*pref_bouzidi*rho_s(i,ny-kk+1,k)* &
       (dex(16)*u_rear+dey(16)*v_rear+dez(16)*w_rear)
    end forall
    
    ! swap pop17 and pop18 along y
    forall(i=sx:ex,k=sz:ez)
      buffservice3d(i,ny+kk,k)=f17sub(i,ny+kk,k)
      f17sub(i,ny+kk,k)=f18sub(i,ny+kk,k) + p(17)*pref_bouzidi*rho_s(i,ny-kk+1,k)* &
       (dex(17)*u_rear+dey(17)*v_rear+dez(17)*w_rear)
      f18sub(i,ny+kk,k)=buffservice3d(i,ny+kk,k) + p(18)*pref_bouzidi*rho_s(i,ny-kk+1,k)* &
       (dex(18)*u_rear+dey(18)*v_rear+dez(18)*w_rear)
    end forall
    
  case (3)
    
    !red fluid
    
    forall(i=sx:ex,k=sz:ez)
      f01sub(i,ny+kk,k)=equil_pop01(rho_rear,u_rear,v_rear,w_rear)
    end forall
    
    forall(i=sx:ex,k=sz:ez)
      f02sub(i,ny+kk,k)=equil_pop02(rho_rear,u_rear,v_rear,w_rear)
    end forall
    
    forall(i=sx:ex,k=sz:ez)
      f03sub(i,ny+kk,k)=equil_pop03(rho_rear,u_rear,v_rear,w_rear)
    end forall
    
    forall(i=sx:ex,k=sz:ez)
      f04sub(i,ny+kk,k)=equil_pop04(rho_rear,u_rear,v_rear,w_rear)
    end forall
    
    forall(i=sx:ex,k=sz:ez)
      f05sub(i,ny+kk,k)=equil_pop05(rho_rear,u_rear,v_rear,w_rear)
    end forall
    
    forall(i=sx:ex,k=sz:ez)
      f06sub(i,ny+kk,k)=equil_pop06(rho_rear,u_rear,v_rear,w_rear)
    end forall
    
    forall(i=sx:ex,k=sz:ez)
      f07sub(i,ny+kk,k)=equil_pop07(rho_rear,u_rear,v_rear,w_rear)
    end forall
    
    forall(i=sx:ex,k=sz:ez)
      f08sub(i,ny+kk,k)=equil_pop08(rho_rear,u_rear,v_rear,w_rear)
    end forall
    
    forall(i=sx:ex,k=sz:ez)
      f09sub(i,ny+kk,k)=equil_pop09(rho_rear,u_rear,v_rear,w_rear)
    end forall
    
    forall(i=sx:ex,k=sz:ez)
      f10sub(i,ny+kk,k)=equil_pop10(rho_rear,u_rear,v_rear,w_rear)
    end forall
    
    forall(i=sx:ex,k=sz:ez)
      f11sub(i,ny+kk,k)=equil_pop11(rho_rear,u_rear,v_rear,w_rear)
    end forall
    
    forall(i=sx:ex,k=sz:ez)
      f12sub(i,ny+kk,k)=equil_pop12(rho_rear,u_rear,v_rear,w_rear)
    end forall
    
    forall(i=sx:ex,k=sz:ez)
      f13sub(i,ny+kk,k)=equil_pop13(rho_rear,u_rear,v_rear,w_rear)
    end forall
    
    forall(i=sx:ex,k=sz:ez)
      f14sub(i,ny+kk,k)=equil_pop14(rho_rear,u_rear,v_rear,w_rear)
    end forall
    
    forall(i=sx:ex,k=sz:ez)
      f15sub(i,ny+kk,k)=equil_pop15(rho_rear,u_rear,v_rear,w_rear)
    end forall
    
    forall(i=sx:ex,k=sz:ez)
      f16sub(i,ny+kk,k)=equil_pop16(rho_rear,u_rear,v_rear,w_rear)
    end forall
    
    forall(i=sx:ex,k=sz:ez)
      f17sub(i,ny+kk,k)=equil_pop17(rho_rear,u_rear,v_rear,w_rear)
    end forall
    
    forall(i=sx:ex,k=sz:ez)
      f18sub(i,ny+kk,k)=equil_pop18(rho_rear,u_rear,v_rear,w_rear)
    end forall
    
  case default
    call error(13)
  end select
  
#else
  
  call warning(9,ZERO,latt_name)
  call error(12)
  
#endif

  return
  
 end subroutine apply_bounceback_rear
 
 subroutine driver_bounceback_edge_front_east(sz,ez)
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the bounceback boundary 
!     condition at the front east edge along the space 
!     defined in input sz,ez
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: sz,ez
  
  call apply_bounceback_edge_front_east(sz,ez,f00R,f01R,f02R,f03R,f04R,&
   f05R,f06R,f07R,f08R,f09R,f10R,f11R,f12R,f13R,f14R,f15R,f16R,f17R, &
   f18R)
  
  if(lsingle_fluid)return
 
  call apply_bounceback_edge_front_east(sz,ez,f00B,f01B,f02B,f03B,f04B,&
   f05B,f06B,f07B,f08B,f09B,f10B,f11B,f12B,f13B,f14B,f15B,f16B,f17B, &
   f18B)
   
  return
  
 end subroutine driver_bounceback_edge_front_east
 
 subroutine apply_bounceback_edge_front_east(sz,ez,f00sub,f01sub, &
  f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub,f09sub,f10sub, &
  f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the bounceback boundary 
!     condition at the front east edge along the space 
!     defined in input sz,ez
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: sz,ez
  
  real(kind=PRC), allocatable, dimension(:,:,:) :: &
   f00sub,f01sub,f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub, &
   f09sub,f10sub,f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub
  
  integer :: i,j,k
  integer, parameter :: kk = 1
  integer, parameter :: kkk = 1
  
#if LATTICE==319

  ! swap pop1 and pop2 x
  forall(k=sz:ez)
    buffservice3d(nx+kkk,1-kk,k)=f01sub(nx+kkk,1-kk,k)
    f01sub(nx+kkk,1-kk,k)=f02sub(nx+kkk,1-kk,k)
    f02sub(nx+kkk,1-kk,k)=buffservice3d(nx+kkk,1-kk,k)
  end forall
  
  ! swap pop3 and pop4 y
  forall(k=sz:ez)
    buffservice3d(nx+kkk,1-kk,k)=f03sub(nx+kkk,1-kk,k)
    f03sub(nx+kkk,1-kk,k)=f04sub(nx+kkk,1-kk,k)
    f04sub(nx+kkk,1-kk,k)=buffservice3d(nx+kkk,1-kk,k)
  end forall
  
  ! swap pop5 and pop6 z
  forall(k=sz:ez)
    buffservice3d(nx+kkk,1-kk,k)=f05sub(nx+kkk,1-kk,k)
    f05sub(nx+kkk,1-kk,k)=f06sub(nx+kkk,1-kk,k)
    f06sub(nx+kkk,1-kk,k)=buffservice3d(nx+kkk,1-kk,k)
  end forall
  
  ! swap pop7 and pop8 y
  forall(k=sz:ez)
    buffservice3d(nx+kkk,1-kk,k)=f07sub(nx+kkk,1-kk,k)
    f07sub(nx+kkk,1-kk,k)=f08sub(nx+kkk,1-kk,k)
    f08sub(nx+kkk,1-kk,k)=buffservice3d(nx+kkk,1-kk,k)
  end forall
  
  ! swap pop9 and pop10
  forall(k=sz:ez)
    buffservice3d(nx+kkk,1-kk,k)=f09sub(nx+kkk,1-kk,k)
    f09sub(nx+kkk,1-kk,k)=f10sub(nx+kkk,1-kk,k)
    f10sub(nx+kkk,1-kk,k)=buffservice3d(nx+kkk,1-kk,k)
  end forall
  
  ! swap pop11 and pop12
  forall(k=sz:ez)
    buffservice3d(nx+kkk,1-kk,k)=f11sub(nx+kkk,1-kk,k)
    f11sub(nx+kkk,1-kk,k)=f12sub(nx+kkk,1-kk,k)
    f12sub(nx+kkk,1-kk,k)=buffservice3d(nx+kkk,1-kk,k)
  end forall
  
  ! swap pop13 and pop14
  forall(k=sz:ez)
    buffservice3d(nx+kkk,1-kk,k)=f13sub(nx+kkk,1-kk,k)
    f13sub(nx+kkk,1-kk,k)=f14sub(nx+kkk,1-kk,k)
    f14sub(nx+kkk,1-kk,k)=buffservice3d(nx+kkk,1-kk,k)
  end forall
  
  ! swap pop15 and pop16
  forall(k=sz:ez)
    buffservice3d(nx+kkk,1-kk,k)=f15sub(nx+kkk,1-kk,k)
    f15sub(nx+kkk,1-kk,k)=f16sub(nx+kkk,1-kk,k)
    f16sub(nx+kkk,1-kk,k)=buffservice3d(nx+kkk,1-kk,k)
  end forall
  
  ! swap pop17 and pop18
  forall(k=sz:ez)
    buffservice3d(nx+kkk,1-kk,k)=f17sub(nx+kkk,1-kk,k)
    f17sub(nx+kkk,1-kk,k)=f18sub(nx+kkk,1-kk,k)
    f18sub(nx+kkk,1-kk,k)=buffservice3d(nx+kkk,1-kk,k)
  end forall
  
#else
  
  call warning(9,ZERO,latt_name)
  call error(12)
  
#endif
  
  return
  
 end subroutine apply_bounceback_edge_front_east
 
 subroutine driver_bounceback_edge_front_west(sz,ez)
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the bounceback boundary 
!     condition at the front west edge along the space 
!     defined in input sz,ez
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: sz,ez
  
  call apply_bounceback_edge_front_west(sz,ez,f00R,f01R,f02R,f03R,f04R,&
   f05R,f06R,f07R,f08R,f09R,f10R,f11R,f12R,f13R,f14R,f15R,f16R,f17R, &
   f18R)
  
  if(lsingle_fluid)return
 
  call apply_bounceback_edge_front_west(sz,ez,f00B,f01B,f02B,f03B,f04B,&
   f05B,f06B,f07B,f08B,f09B,f10B,f11B,f12B,f13B,f14B,f15B,f16B,f17B, &
   f18B)
   
  return
  
 end subroutine driver_bounceback_edge_front_west
 
 subroutine apply_bounceback_edge_front_west(sz,ez,f00sub,f01sub, &
  f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub,f09sub,f10sub, &
  f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the bounceback boundary 
!     condition at the front west edge along the space 
!     defined in input sz,ez
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: sz,ez
  
  real(kind=PRC), allocatable, dimension(:,:,:) :: &
   f00sub,f01sub,f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub, &
   f09sub,f10sub,f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub
  
  integer :: i,j,k
  integer, parameter :: kk = 1
  integer, parameter :: kkk = 1
  
  
#if LATTICE==319
  
  ! swap pop1 and pop2 x
  forall(k=sz:ez)
    buffservice3d(1-kkk,1-kk,k)=f01sub(1-kkk,1-kk,k)
    f01sub(1-kkk,1-kk,k)=f02sub(1-kkk,1-kk,k)
    f02sub(1-kkk,1-kk,k)=buffservice3d(1-kkk,1-kk,k)
  end forall
  
  ! swap pop3 and pop4 y
  forall(k=sz:ez)
    buffservice3d(1-kkk,1-kk,k)=f03sub(1-kkk,1-kk,k)
    f03sub(1-kkk,1-kk,k)=f04sub(1-kkk,1-kk,k)
    f04sub(1-kkk,1-kk,k)=buffservice3d(1-kkk,1-kk,k)
  end forall
  
  ! swap pop5 and pop6 z
  forall(k=sz:ez)
    buffservice3d(1-kkk,1-kk,k)=f05sub(1-kkk,1-kk,k)
    f05sub(1-kkk,1-kk,k)=f06sub(1-kkk,1-kk,k)
    f06sub(1-kkk,1-kk,k)=buffservice3d(1-kkk,1-kk,k)
  end forall
  
  ! swap pop7 and pop8 y
  forall(k=sz:ez)
    buffservice3d(1-kkk,1-kk,k)=f07sub(1-kkk,1-kk,k)
    f07sub(1-kkk,1-kk,k)=f08sub(1-kkk,1-kk,k)
    f08sub(1-kkk,1-kk,k)=buffservice3d(1-kkk,1-kk,k)
  end forall
  
  ! swap pop9 and pop10
  forall(k=sz:ez)
    buffservice3d(1-kkk,1-kk,k)=f09sub(1-kkk,1-kk,k)
    f09sub(1-kkk,1-kk,k)=f10sub(1-kkk,1-kk,k)
    f10sub(1-kkk,1-kk,k)=buffservice3d(1-kkk,1-kk,k)
  end forall
  
  ! swap pop11 and pop12
  forall(k=sz:ez)
    buffservice3d(1-kkk,1-kk,k)=f11sub(1-kkk,1-kk,k)
    f11sub(1-kkk,1-kk,k)=f12sub(1-kkk,1-kk,k)
    f12sub(1-kkk,1-kk,k)=buffservice3d(1-kkk,1-kk,k)
  end forall
  
  ! swap pop13 and pop14
  forall(k=sz:ez)
    buffservice3d(1-kkk,1-kk,k)=f13sub(1-kkk,1-kk,k)
    f13sub(1-kkk,1-kk,k)=f14sub(1-kkk,1-kk,k)
    f14sub(1-kkk,1-kk,k)=buffservice3d(1-kkk,1-kk,k)
  end forall
  
  ! swap pop15 and pop16
  forall(k=sz:ez)
    buffservice3d(1-kkk,1-kk,k)=f15sub(1-kkk,1-kk,k)
    f15sub(1-kkk,1-kk,k)=f16sub(1-kkk,1-kk,k)
    f16sub(1-kkk,1-kk,k)=buffservice3d(1-kkk,1-kk,k)
  end forall
  
  ! swap pop17 and pop18
  forall(k=sz:ez)
    buffservice3d(1-kkk,1-kk,k)=f17sub(1-kkk,1-kk,k)
    f17sub(1-kkk,1-kk,k)=f18sub(1-kkk,1-kk,k)
    f18sub(1-kkk,1-kk,k)=buffservice3d(1-kkk,1-kk,k)
  end forall
  
#else
  
  call warning(9,ZERO,latt_name)
  call error(12)
  
#endif
  return
  
 end subroutine apply_bounceback_edge_front_west
 
 subroutine driver_bounceback_edge_north_east(sy,ey)
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the bounceback boundary 
!     condition at the north east edge along the space 
!     defined in input sy,ey
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: sy,ey
  
  call apply_bounceback_edge_north_east(sy,ey,f00R,f01R,f02R,f03R,f04R,&
   f05R,f06R,f07R,f08R,f09R,f10R,f11R,f12R,f13R,f14R,f15R,f16R,f17R, &
   f18R)
  
  if(lsingle_fluid)return
 
  call apply_bounceback_edge_north_east(sy,ey,f00B,f01B,f02B,f03B,f04B,&
   f05B,f06B,f07B,f08B,f09B,f10B,f11B,f12B,f13B,f14B,f15B,f16B,f17B, &
   f18B)
   
  return
  
 end subroutine driver_bounceback_edge_north_east
 
 subroutine apply_bounceback_edge_north_east(sy,ey,f00sub,f01sub, &
  f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub,f09sub,f10sub, &
  f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the bounceback boundary 
!     condition at the north east edge along the space 
!     defined in input sy,ey
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: sy,ey
  real(kind=PRC), allocatable, dimension(:,:,:) :: &
   f00sub,f01sub,f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub, &
   f09sub,f10sub,f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub
  
  integer :: i,j,k
  integer, parameter :: kk = 1
  integer, parameter :: kkk = 1
  
#if LATTICE==319

  ! swap pop1 and pop2 x
  forall(j=sy:ey)
    buffservice3d(nx+kkk,j,nz+kk)=f01sub(nx+kkk,j,nz+kk)
    f01sub(nx+kkk,j,nz+kk)=f02sub(nx+kkk,j,nz+kk)
    f02sub(nx+kkk,j,nz+kk)=buffservice3d(nx+kkk,j,nz+kk)
  end forall
  
  ! swap pop3 and pop4 y
  forall(j=sy:ey)
    buffservice3d(nx+kkk,j,nz+kk)=f03sub(nx+kkk,j,nz+kk)
    f03sub(nx+kkk,j,nz+kk)=f04sub(nx+kkk,j,nz+kk)
    f04sub(nx+kkk,j,nz+kk)=buffservice3d(nx+kkk,j,nz+kk)
  end forall
  
  ! swap pop5 and pop6 z
  forall(j=sy:ey)
    buffservice3d(nx+kkk,j,nz+kk)=f05sub(nx+kkk,j,nz+kk)
    f05sub(nx+kkk,j,nz+kk)=f06sub(nx+kkk,j,nz+kk)
    f06sub(nx+kkk,j,nz+kk)=buffservice3d(nx+kkk,j,nz+kk)
  end forall
  
  ! swap pop7 and pop8
  forall(j=sy:ey)
    buffservice3d(nx+kkk,j,nz+kk)=f07sub(nx+kkk,j,nz+kk)
    f07sub(nx+kkk,j,nz+kk)=f08sub(nx+kkk,j,nz+kk)
    f08sub(nx+kkk,j,nz+kk)=buffservice3d(nx+kkk,j,nz+kk)
  end forall
  
  ! swap pop9 and pop10
  forall(j=sy:ey)
    buffservice3d(nx+kkk,j,nz+kk)=f09sub(nx+kkk,j,nz+kk)
    f09sub(nx+kkk,j,nz+kk)=f10sub(nx+kkk,j,nz+kk)
    f10sub(nx+kkk,j,nz+kk)=buffservice3d(nx+kkk,j,nz+kk)
  end forall
  
  ! swap pop11 and pop12
  forall(j=sy:ey)
    buffservice3d(nx+kkk,j,nz+kk)=f11sub(nx+kkk,j,nz+kk)
    f11sub(nx+kkk,j,nz+kk)=f12sub(nx+kkk,j,nz+kk)
    f12sub(nx+kkk,j,nz+kk)=buffservice3d(nx+kkk,j,nz+kk)
  end forall
  
  ! swap pop13 and pop14
  forall(j=sy:ey)
    buffservice3d(nx+kkk,j,nz+kk)=f13sub(nx+kkk,j,nz+kk)
    f13sub(nx+kkk,j,nz+kk)=f14sub(nx+kkk,j,nz+kk)
    f14sub(nx+kkk,j,nz+kk)=buffservice3d(nx+kkk,j,nz+kk)
  end forall
  
  ! swap pop15 and pop16
  forall(j=sy:ey)
    buffservice3d(nx+kkk,j,nz+kk)=f15sub(nx+kkk,j,nz+kk)
    f15sub(nx+kkk,j,nz+kk)=f16sub(nx+kkk,j,nz+kk)
    f16sub(nx+kkk,j,nz+kk)=buffservice3d(nx+kkk,j,nz+kk)
  end forall
  
  ! swap pop17 and pop18
  forall(j=sy:ey)
    buffservice3d(nx+kkk,j,nz+kk)=f17sub(nx+kkk,j,nz+kk)
    f17sub(nx+kkk,j,nz+kk)=f18sub(nx+kkk,j,nz+kk)
    f18sub(nx+kkk,j,nz+kk)=buffservice3d(nx+kkk,j,nz+kk)
  end forall
  
#else
  
  call warning(9,ZERO,latt_name)
  call error(12)
  
#endif
  return
  
 end subroutine apply_bounceback_edge_north_east
 
 subroutine driver_bounceback_edge_north_front(sx,ex)
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the bounceback boundary 
!     condition at the north front edge along the space 
!     defined in input sx,ex
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: sx,ex
  
  call apply_bounceback_edge_north_front(sx,ex,f00R,f01R,f02R,f03R,f04R,&
   f05R,f06R,f07R,f08R,f09R,f10R,f11R,f12R,f13R,f14R,f15R,f16R,f17R, &
   f18R)
  
  if(lsingle_fluid)return
 
  call apply_bounceback_edge_north_front(sx,ex,f00B,f01B,f02B,f03B,f04B,&
   f05B,f06B,f07B,f08B,f09B,f10B,f11B,f12B,f13B,f14B,f15B,f16B,f17B, &
   f18B)
   
  return
  
 end subroutine driver_bounceback_edge_north_front
 
 subroutine apply_bounceback_edge_north_front(sx,ex,f00sub,f01sub, &
  f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub,f09sub,f10sub, &
  f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the bounceback boundary 
!     condition at the north front edge along the space 
!     defined in input sx,ex
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: sx,ex
  real(kind=PRC), allocatable, dimension(:,:,:) :: &
   f00sub,f01sub,f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub, &
   f09sub,f10sub,f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub
  
  integer :: i,j,k
  integer, parameter :: kk = 1
  integer, parameter :: kkk = 1
  
#if LATTICE==319
  
  ! swap pop1 and pop2 x
  forall(i=sx:ex)
    buffservice3d(i,1-kkk,nz+kk)=f01sub(i,1-kkk,nz+kk)
    f01sub(i,1-kkk,nz+kk)=f02sub(i,1-kkk,nz+kk)
    f02sub(i,1-kkk,nz+kk)=buffservice3d(i,1-kkk,nz+kk)
  end forall
  
  ! swap pop3 and pop4 y
  forall(i=sx:ex)
    buffservice3d(i,1-kkk,nz+kk)=f03sub(i,1-kkk,nz+kk)
    f03sub(i,1-kkk,nz+kk)=f04sub(i,1-kkk,nz+kk)
    f04sub(i,1-kkk,nz+kk)=buffservice3d(i,1-kkk,nz+kk)
  end forall
  
  ! swap pop5 and pop6 z
  forall(i=sx:ex)
    buffservice3d(i,1-kkk,nz+kk)=f05sub(i,1-kkk,nz+kk)
    f05sub(i,1-kkk,nz+kk)=f06sub(i,1-kkk,nz+kk)
    f06sub(i,1-kkk,nz+kk)=buffservice3d(i,1-kkk,nz+kk)
  end forall
  
  ! swap pop7 and pop8
  forall(i=sx:ex)
    buffservice3d(i,1-kkk,nz+kk)=f07sub(i,1-kkk,nz+kk)
    f07sub(i,1-kkk,nz+kk)=f08sub(i,1-kkk,nz+kk)
    f08sub(i,1-kkk,nz+kk)=buffservice3d(i,1-kkk,nz+kk)
  end forall
  
  ! swap pop9 and pop10
  forall(i=sx:ex)
    buffservice3d(i,1-kkk,nz+kk)=f09sub(i,1-kkk,nz+kk)
    f09sub(i,1-kkk,nz+kk)=f10sub(i,1-kkk,nz+kk)
    f10sub(i,1-kkk,nz+kk)=buffservice3d(i,1-kkk,nz+kk)
  end forall
  
  ! swap pop11 and pop12
  forall(i=sx:ex)
    buffservice3d(i,1-kkk,nz+kk)=f11sub(i,1-kkk,nz+kk)
    f11sub(i,1-kkk,nz+kk)=f12sub(i,1-kkk,nz+kk)
    f12sub(i,1-kkk,nz+kk)=buffservice3d(i,1-kkk,nz+kk)
  end forall
  
  ! swap pop13 and pop14
  forall(i=sx:ex)
    buffservice3d(i,1-kkk,nz+kk)=f13sub(i,1-kkk,nz+kk)
    f13sub(i,1-kkk,nz+kk)=f14sub(i,1-kkk,nz+kk)
    f14sub(i,1-kkk,nz+kk)=buffservice3d(i,1-kkk,nz+kk)
  end forall
  
  ! swap pop15 and pop16
  forall(i=sx:ex)
    buffservice3d(i,1-kkk,nz+kk)=f15sub(i,1-kkk,nz+kk)
    f15sub(i,1-kkk,nz+kk)=f16sub(i,1-kkk,nz+kk)
    f16sub(i,1-kkk,nz+kk)=buffservice3d(i,1-kkk,nz+kk)
  end forall
  
  ! swap pop17 and pop18
  forall(i=sx:ex)
    buffservice3d(i,1-kkk,nz+kk)=f17sub(i,1-kkk,nz+kk)
    f17sub(i,1-kkk,nz+kk)=f18sub(i,1-kkk,nz+kk)
    f18sub(i,1-kkk,nz+kk)=buffservice3d(i,1-kkk,nz+kk)
  end forall
  
#else
  
  call warning(9,ZERO,latt_name)
  call error(12)
  
#endif
  return
  
 end subroutine apply_bounceback_edge_north_front
 
 subroutine driver_bounceback_edge_north_rear(sx,ex)
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the bounceback boundary 
!     condition at the north rear edge along the space 
!     defined in input sx,ex
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: sx,ex
  
  call apply_bounceback_edge_north_rear(sx,ex,f00R,f01R,f02R,f03R,f04R,&
   f05R,f06R,f07R,f08R,f09R,f10R,f11R,f12R,f13R,f14R,f15R,f16R,f17R, &
   f18R)
  
  if(lsingle_fluid)return
 
  call apply_bounceback_edge_north_rear(sx,ex,f00B,f01B,f02B,f03B,f04B,&
   f05B,f06B,f07B,f08B,f09B,f10B,f11B,f12B,f13B,f14B,f15B,f16B,f17B, &
   f18B)
   
  return
  
 end subroutine driver_bounceback_edge_north_rear
 
 subroutine apply_bounceback_edge_north_rear(sx,ex,f00sub,f01sub, &
  f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub,f09sub,f10sub, &
  f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the bounceback boundary 
!     condition at the north rear edge
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: sx,ex
  real(kind=PRC), allocatable, dimension(:,:,:) :: &
   f00sub,f01sub,f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub, &
   f09sub,f10sub,f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub
  
  integer :: i,j,k
  integer, parameter :: kk = 1
  integer, parameter :: kkk = 1
  
#if LATTICE==319
  
  ! swap pop1 and pop2 y
  forall(i=sx:ex)
    buffservice3d(i,ny+kkk,nz+kk)=f01sub(i,ny+kkk,nz+kk)
    f01sub(i,ny+kkk,nz+kk)=f02sub(i,ny+kkk,nz+kk)
    f02sub(i,ny+kkk,nz+kk)=buffservice3d(i,ny+kkk,nz+kk)
  end forall
  
 ! swap pop3 and pop4 y
  forall(i=sx:ex)
    buffservice3d(i,ny+kkk,nz+kk)=f03sub(i,ny+kkk,nz+kk)
    f03sub(i,ny+kkk,nz+kk)=f04sub(i,ny+kkk,nz+kk)
    f04sub(i,ny+kkk,nz+kk)=buffservice3d(i,ny+kkk,nz+kk)
  end forall
  
  ! swap pop5 and pop6 z
  forall(i=sx:ex)
    buffservice3d(i,ny+kkk,nz+kk)=f05sub(i,ny+kkk,nz+kk)
    f05sub(i,ny+kkk,nz+kk)=f06sub(i,ny+kkk,nz+kk)
    f06sub(i,ny+kkk,nz+kk)=buffservice3d(i,ny+kkk,nz+kk)
  end forall
  
  ! swap pop7 and pop8
  forall(i=sx:ex)
    buffservice3d(i,ny+kkk,nz+kk)=f07sub(i,ny+kkk,nz+kk)
    f07sub(i,ny+kkk,nz+kk)=f08sub(i,ny+kkk,nz+kk)
    f08sub(i,ny+kkk,nz+kk)=buffservice3d(i,ny+kkk,nz+kk)
  end forall
  
  ! swap pop9 and pop10
  forall(i=sx:ex)
    buffservice3d(i,ny+kkk,nz+kk)=f09sub(i,ny+kkk,nz+kk)
    f09sub(i,ny+kkk,nz+kk)=f10sub(i,ny+kkk,nz+kk)
    f10sub(i,ny+kkk,nz+kk)=buffservice3d(i,ny+kkk,nz+kk)
  end forall
  
  ! swap pop11 and pop12
  forall(i=sx:ex)
    buffservice3d(i,ny+kkk,nz+kk)=f11sub(i,ny+kkk,nz+kk)
    f11sub(i,ny+kkk,nz+kk)=f12sub(i,ny+kkk,nz+kk)
    f12sub(i,ny+kkk,nz+kk)=buffservice3d(i,ny+kkk,nz+kk)
  end forall
  
  ! swap pop13 and pop14
  forall(i=sx:ex)
    buffservice3d(i,ny+kkk,nz+kk)=f13sub(i,ny+kkk,nz+kk)
    f13sub(i,ny+kkk,nz+kk)=f14sub(i,ny+kkk,nz+kk)
    f14sub(i,ny+kkk,nz+kk)=buffservice3d(i,ny+kkk,nz+kk)
  end forall
  
  ! swap pop15 and pop16
  forall(i=sx:ex)
    buffservice3d(i,ny+kkk,nz+kk)=f15sub(i,ny+kkk,nz+kk)
    f15sub(i,ny+kkk,nz+kk)=f16sub(i,ny+kkk,nz+kk)
    f16sub(i,ny+kkk,nz+kk)=buffservice3d(i,ny+kkk,nz+kk)
  end forall
  
  ! swap pop17 and pop18
  forall(i=sx:ex)
    buffservice3d(i,ny+kkk,nz+kk)=f17sub(i,ny+kkk,nz+kk)
    f17sub(i,ny+kkk,nz+kk)=f18sub(i,ny+kkk,nz+kk)
    f18sub(i,ny+kkk,nz+kk)=buffservice3d(i,ny+kkk,nz+kk)
  end forall
  
#else
  
  call warning(9,ZERO,latt_name)
  call error(12)
  
#endif
  return
  
 end subroutine apply_bounceback_edge_north_rear
 
 subroutine driver_bounceback_edge_north_west(sy,ey)
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the bounceback boundary 
!     condition at the north west edge along the space 
!     defined in input sy,ey
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: sy,ey
  
  call apply_bounceback_edge_north_west(sy,ey,f00R,f01R,f02R,f03R,f04R,&
   f05R,f06R,f07R,f08R,f09R,f10R,f11R,f12R,f13R,f14R,f15R,f16R,f17R, &
   f18R)
  
  if(lsingle_fluid)return
 
  call apply_bounceback_edge_north_west(sy,ey,f00B,f01B,f02B,f03B,f04B,&
   f05B,f06B,f07B,f08B,f09B,f10B,f11B,f12B,f13B,f14B,f15B,f16B,f17B, &
   f18B)
   
  return
  
 end subroutine driver_bounceback_edge_north_west
 
 subroutine apply_bounceback_edge_north_west(sy,ey,f00sub,f01sub, &
  f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub,f09sub,f10sub, &
  f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the bounceback boundary 
!     condition at the north west edge along the space 
!     defined in input sy,ey
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: sy,ey
  real(kind=PRC), allocatable, dimension(:,:,:) :: &
   f00sub,f01sub,f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub, &
   f09sub,f10sub,f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub
  
  integer :: i,j,k
  integer, parameter :: kk = 1
  integer, parameter :: kkk = 1
  
#if LATTICE==319

  ! swap pop1 and pop2 x
  forall(j=sy:ey)
    buffservice3d(1-kkk,j,nz+kk)=f01sub(1-kkk,j,nz+kk)
    f01sub(1-kkk,j,nz+kk)=f02sub(1-kkk,j,nz+kk)
    f02sub(1-kkk,j,nz+kk)=buffservice3d(1-kkk,j,nz+kk)
  end forall
  
  ! swap pop3 and pop4 y
  forall(j=sy:ey)
    buffservice3d(1-kkk,j,nz+kk)=f03sub(1-kkk,j,nz+kk)
    f03sub(1-kkk,j,nz+kk)=f04sub(1-kkk,j,nz+kk)
    f04sub(1-kkk,j,nz+kk)=buffservice3d(1-kkk,j,nz+kk)
  end forall
  
  ! swap pop5 and pop6 z
  forall(j=sy:ey)
    buffservice3d(1-kkk,j,nz+kk)=f05sub(1-kkk,j,nz+kk)
    f05sub(1-kkk,j,nz+kk)=f06sub(1-kkk,j,nz+kk)
    f06sub(1-kkk,j,nz+kk)=buffservice3d(1-kkk,j,nz+kk)
  end forall
  
  ! swap pop7 and pop8
  forall(j=sy:ey)
    buffservice3d(1-kkk,j,nz+kk)=f07sub(1-kkk,j,nz+kk)
    f07sub(1-kkk,j,nz+kk)=f08sub(1-kkk,j,nz+kk)
    f08sub(1-kkk,j,nz+kk)=buffservice3d(1-kkk,j,nz+kk)
  end forall
  
  ! swap pop9 and pop10
  forall(j=sy:ey)
    buffservice3d(1-kkk,j,nz+kk)=f09sub(1-kkk,j,nz+kk)
    f09sub(1-kkk,j,nz+kk)=f10sub(1-kkk,j,nz+kk)
    f10sub(1-kkk,j,nz+kk)=buffservice3d(1-kkk,j,nz+kk)
  end forall
  
  ! swap pop11 and pop12
  forall(j=sy:ey)
    buffservice3d(1-kkk,j,nz+kk)=f11sub(1-kkk,j,nz+kk)
    f11sub(1-kkk,j,nz+kk)=f12sub(1-kkk,j,nz+kk)
    f12sub(1-kkk,j,nz+kk)=buffservice3d(1-kkk,j,nz+kk)
  end forall
  
  ! swap pop13 and pop14
  forall(j=sy:ey)
    buffservice3d(1-kkk,j,nz+kk)=f13sub(1-kkk,j,nz+kk)
    f13sub(1-kkk,j,nz+kk)=f14sub(1-kkk,j,nz+kk)
    f14sub(1-kkk,j,nz+kk)=buffservice3d(1-kkk,j,nz+kk)
  end forall
  
  ! swap pop15 and pop16
  forall(j=sy:ey)
    buffservice3d(1-kkk,j,nz+kk)=f15sub(1-kkk,j,nz+kk)
    f15sub(1-kkk,j,nz+kk)=f16sub(1-kkk,j,nz+kk)
    f16sub(1-kkk,j,nz+kk)=buffservice3d(1-kkk,j,nz+kk)
  end forall
  
  ! swap pop17 and pop18
  forall(j=sy:ey)
    buffservice3d(1-kkk,j,nz+kk)=f17sub(1-kkk,j,nz+kk)
    f17sub(1-kkk,j,nz+kk)=f18sub(1-kkk,j,nz+kk)
    f18sub(1-kkk,j,nz+kk)=buffservice3d(1-kkk,j,nz+kk)
  end forall
  
#else
  
  call warning(9,ZERO,latt_name)
  call error(12)
  
#endif
  return
  
 end subroutine apply_bounceback_edge_north_west
 
 subroutine driver_bounceback_edge_rear_east(sz,ez)
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the bounceback boundary 
!     condition at the rear east edge along the space 
!     defined in input sz,ez
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: sz,ez
  
  call apply_bounceback_edge_rear_east(sz,ez,f00R,f01R,f02R,f03R,f04R,&
   f05R,f06R,f07R,f08R,f09R,f10R,f11R,f12R,f13R,f14R,f15R,f16R,f17R, &
   f18R)
  
  if(lsingle_fluid)return
 
  call apply_bounceback_edge_rear_east(sz,ez,f00B,f01B,f02B,f03B,f04B,&
   f05B,f06B,f07B,f08B,f09B,f10B,f11B,f12B,f13B,f14B,f15B,f16B,f17B, &
   f18B)
   
  return
  
 end subroutine driver_bounceback_edge_rear_east
 
 subroutine apply_bounceback_edge_rear_east(sz,ez,f00sub,f01sub, &
  f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub,f09sub,f10sub, &
  f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the bounceback boundary 
!     condition at the rear east edge along the space 
!     defined in input sz,ez
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: sz,ez
  real(kind=PRC), allocatable, dimension(:,:,:) :: &
   f00sub,f01sub,f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub, &
   f09sub,f10sub,f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub
  
  integer :: i,j,k
  integer, parameter :: kk = 1
  integer, parameter :: kkk = 1
  
#if LATTICE==319

  ! swap pop1 and pop2 x
  forall(k=sz:ez)
    buffservice3d(nx+kkk,ny+kk,k)=f01sub(nx+kkk,ny+kk,k)
    f01sub(nx+kkk,ny+kk,k)=f02sub(nx+kkk,ny+kk,k)
    f02sub(nx+kkk,ny+kk,k)=buffservice3d(nx+kkk,ny+kk,k)
  end forall
  
  ! swap pop3 and pop4 y
  forall(k=sz:ez)
    buffservice3d(nx+kkk,ny+kk,k)=f03sub(nx+kkk,ny+kk,k)
    f03sub(nx+kkk,ny+kk,k)=f04sub(nx+kkk,ny+kk,k)
    f04sub(nx+kkk,ny+kk,k)=buffservice3d(nx+kkk,ny+kk,k)
  end forall
  
  ! swap pop5 and pop6 z
  forall(k=sz:ez)
    buffservice3d(nx+kkk,ny+kk,k)=f05sub(nx+kkk,ny+kk,k)
    f05sub(nx+kkk,ny+kk,k)=f06sub(nx+kkk,ny+kk,k)
    f06sub(nx+kkk,ny+kk,k)=buffservice3d(nx+kkk,ny+kk,k)
  end forall
  
  ! swap pop7 and pop8 y
  forall(k=sz:ez)
    buffservice3d(nx+kkk,ny+kk,k)=f07sub(nx+kkk,ny+kk,k)
    f07sub(nx+kkk,ny+kk,k)=f08sub(nx+kkk,ny+kk,k)
    f08sub(nx+kkk,ny+kk,k)=buffservice3d(nx+kkk,ny+kk,k)
  end forall
  
  ! swap pop9 and pop10
  forall(k=sz:ez)
    buffservice3d(nx+kkk,ny+kk,k)=f09sub(nx+kkk,ny+kk,k)
    f09sub(nx+kkk,ny+kk,k)=f10sub(nx+kkk,ny+kk,k)
    f10sub(nx+kkk,ny+kk,k)=buffservice3d(nx+kkk,ny+kk,k)
  end forall
  
  ! swap pop11 and pop12
  forall(k=sz:ez)
    buffservice3d(nx+kkk,ny+kk,k)=f11sub(nx+kkk,ny+kk,k)
    f11sub(nx+kkk,ny+kk,k)=f12sub(nx+kkk,ny+kk,k)
    f12sub(nx+kkk,ny+kk,k)=buffservice3d(nx+kkk,ny+kk,k)
  end forall
  
  ! swap pop13 and pop14
  forall(k=sz:ez)
    buffservice3d(nx+kkk,ny+kk,k)=f13sub(nx+kkk,ny+kk,k)
    f13sub(nx+kkk,ny+kk,k)=f14sub(nx+kkk,ny+kk,k)
    f14sub(nx+kkk,ny+kk,k)=buffservice3d(nx+kkk,ny+kk,k)
  end forall
  
  ! swap pop15 and pop16
  forall(k=sz:ez)
    buffservice3d(nx+kkk,ny+kk,k)=f15sub(nx+kkk,ny+kk,k)
    f15sub(nx+kkk,ny+kk,k)=f16sub(nx+kkk,ny+kk,k)
    f16sub(nx+kkk,ny+kk,k)=buffservice3d(nx+kkk,ny+kk,k)
  end forall
  
  ! swap pop17 and pop18
  forall(k=sz:ez)
    buffservice3d(nx+kkk,ny+kk,k)=f17sub(nx+kkk,ny+kk,k)
    f17sub(nx+kkk,ny+kk,k)=f18sub(nx+kkk,ny+kk,k)
    f18sub(nx+kkk,ny+kk,k)=buffservice3d(nx+kkk,ny+kk,k)
  end forall
  
#else
  
  call warning(9,ZERO,latt_name)
  call error(12)
  
#endif
  return
  
 end subroutine apply_bounceback_edge_rear_east
 
 subroutine driver_bounceback_edge_rear_west(sz,ez)
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the bounceback boundary 
!     condition at the rear west edge along the space 
!     defined in input sz,ez
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: sz,ez
  
  call apply_bounceback_edge_rear_west(sz,ez,f00R,f01R,f02R,f03R,f04R,&
   f05R,f06R,f07R,f08R,f09R,f10R,f11R,f12R,f13R,f14R,f15R,f16R,f17R, &
   f18R)
  
  if(lsingle_fluid)return
 
  call apply_bounceback_edge_rear_west(sz,ez,f00B,f01B,f02B,f03B,f04B,&
   f05B,f06B,f07B,f08B,f09B,f10B,f11B,f12B,f13B,f14B,f15B,f16B,f17B, &
   f18B)
   
  return
  
 end subroutine driver_bounceback_edge_rear_west
 
 subroutine apply_bounceback_edge_rear_west(sz,ez,f00sub,f01sub, &
  f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub,f09sub,f10sub, &
  f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the bounceback boundary 
!     condition at the rear west edge along the space 
!     defined in input sz,ez
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: sz,ez
  real(kind=PRC), allocatable, dimension(:,:,:) :: &
   f00sub,f01sub,f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub, &
   f09sub,f10sub,f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub
  
  integer :: i,j,k
  integer, parameter :: kk = 1
  integer, parameter :: kkk = 1
  
#if LATTICE==319

  ! swap pop1 and pop2 x
  forall(k=sz:ez)
    buffservice3d(1-kkk,ny+kk,k)=f01sub(1-kkk,ny+kk,k)
    f01sub(1-kkk,ny+kk,k)=f02sub(1-kkk,ny+kk,k)
    f02sub(1-kkk,ny+kk,k)=buffservice3d(1-kkk,ny+kk,k)
  end forall
  
  ! swap pop3 and pop4 y
  forall(k=sz:ez)
    buffservice3d(1-kkk,ny+kk,k)=f03sub(1-kkk,ny+kk,k)
    f03sub(1-kkk,ny+kk,k)=f04sub(1-kkk,ny+kk,k)
    f04sub(1-kkk,ny+kk,k)=buffservice3d(1-kkk,ny+kk,k)
  end forall
  
  ! swap pop5 and pop6 z
  forall(k=sz:ez)
    buffservice3d(1-kkk,ny+kk,k)=f05sub(1-kkk,ny+kk,k)
    f05sub(1-kkk,ny+kk,k)=f06sub(1-kkk,ny+kk,k)
    f06sub(1-kkk,ny+kk,k)=buffservice3d(1-kkk,ny+kk,k)
  end forall
  
  ! swap pop7 and pop8 y
  forall(k=sz:ez)
    buffservice3d(1-kkk,ny+kk,k)=f07sub(1-kkk,ny+kk,k)
    f07sub(1-kkk,ny+kk,k)=f08sub(1-kkk,ny+kk,k)
    f08sub(1-kkk,ny+kk,k)=buffservice3d(1-kkk,ny+kk,k)
  end forall
  
  ! swap pop9 and pop10
  forall(k=sz:ez)
    buffservice3d(1-kkk,ny+kk,k)=f09sub(1-kkk,ny+kk,k)
    f09sub(1-kkk,ny+kk,k)=f10sub(1-kkk,ny+kk,k)
    f10sub(1-kkk,ny+kk,k)=buffservice3d(1-kkk,ny+kk,k)
  end forall
  
  ! swap pop11 and pop12
  forall(k=sz:ez)
    buffservice3d(1-kkk,ny+kk,k)=f11sub(1-kkk,ny+kk,k)
    f11sub(1-kkk,ny+kk,k)=f12sub(1-kkk,ny+kk,k)
    f12sub(1-kkk,ny+kk,k)=buffservice3d(1-kkk,ny+kk,k)
  end forall
  
  ! swap pop13 and pop14
  forall(k=sz:ez)
    buffservice3d(1-kkk,ny+kk,k)=f13sub(1-kkk,ny+kk,k)
    f13sub(1-kkk,ny+kk,k)=f14sub(1-kkk,ny+kk,k)
    f14sub(1-kkk,ny+kk,k)=buffservice3d(1-kkk,ny+kk,k)
  end forall
  
  ! swap pop15 and pop16
  forall(k=sz:ez)
    buffservice3d(1-kkk,ny+kk,k)=f15sub(1-kkk,ny+kk,k)
    f15sub(1-kkk,ny+kk,k)=f16sub(1-kkk,ny+kk,k)
    f16sub(1-kkk,ny+kk,k)=buffservice3d(1-kkk,ny+kk,k)
  end forall
  
  ! swap pop17 and pop18
  forall(k=sz:ez)
    buffservice3d(1-kkk,ny+kk,k)=f17sub(1-kkk,ny+kk,k)
    f17sub(1-kkk,ny+kk,k)=f18sub(1-kkk,ny+kk,k)
    f18sub(1-kkk,ny+kk,k)=buffservice3d(1-kkk,ny+kk,k)
  end forall
  
#else
  
  call warning(9,ZERO,latt_name)
  call error(12)
  
#endif
  return
  
 end subroutine apply_bounceback_edge_rear_west
 
 subroutine driver_bounceback_edge_south_east(sy,ey)
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the bounceback boundary 
!     condition at the south east edge along the space 
!     defined in input sy,ey
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: sy,ey
  
  call apply_bounceback_edge_south_east(sy,ey,f00R,f01R,f02R,f03R,f04R,&
   f05R,f06R,f07R,f08R,f09R,f10R,f11R,f12R,f13R,f14R,f15R,f16R,f17R, &
   f18R)
  
  if(lsingle_fluid)return
 
  call apply_bounceback_edge_south_east(sy,ey,f00B,f01B,f02B,f03B,f04B,&
   f05B,f06B,f07B,f08B,f09B,f10B,f11B,f12B,f13B,f14B,f15B,f16B,f17B, &
   f18B)
   
  return
  
 end subroutine driver_bounceback_edge_south_east
 
 subroutine apply_bounceback_edge_south_east(sy,ey,f00sub,f01sub, &
  f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub,f09sub,f10sub, &
  f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the bounceback boundary 
!     condition at the south east edge along the space 
!     defined in input sy,ey
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: sy,ey
  real(kind=PRC), allocatable, dimension(:,:,:) :: &
   f00sub,f01sub,f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub, &
   f09sub,f10sub,f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub
  
  integer :: i,j,k
  integer, parameter :: kk = 1
  integer, parameter :: kkk = 1
  
#if LATTICE==319
  
  ! swap pop1 and pop2 x
  forall(j=sy:ey)
    buffservice3d(nx+kkk,j,1-kk)=f01sub(nx+kkk,j,1-kk)
    f01sub(nx+kkk,j,1-kk)=f02sub(nx+kkk,j,1-kk)
    f02sub(nx+kkk,j,1-kk)=buffservice3d(nx+kkk,j,1-kk)
  end forall
  
  ! swap pop3 and pop4 y
  forall(j=sy:ey)
    buffservice3d(nx+kkk,j,1-kk)=f03sub(nx+kkk,j,1-kk)
    f03sub(nx+kkk,j,1-kk)=f04sub(nx+kkk,j,1-kk)
    f04sub(nx+kkk,j,1-kk)=buffservice3d(nx+kkk,j,1-kk)
  end forall
  
  ! swap pop5 and pop6 z
  forall(j=sy:ey)
    buffservice3d(nx+kkk,j,1-kk)=f05sub(nx+kkk,j,1-kk)
    f05sub(nx+kkk,j,1-kk)=f06sub(nx+kkk,j,1-kk)
    f06sub(nx+kkk,j,1-kk)=buffservice3d(nx+kkk,j,1-kk)
  end forall
  
  ! swap pop7 and pop8
  forall(j=sy:ey)
    buffservice3d(nx+kkk,j,1-kk)=f07sub(nx+kkk,j,1-kk)
    f07sub(nx+kkk,j,1-kk)=f08sub(nx+kkk,j,1-kk)
    f08sub(nx+kkk,j,1-kk)=buffservice3d(nx+kkk,j,1-kk)
  end forall
  
  ! swap pop9 and pop10
  forall(j=sy:ey)
    buffservice3d(nx+kkk,j,1-kk)=f09sub(nx+kkk,j,1-kk)
    f09sub(nx+kkk,j,1-kk)=f10sub(nx+kkk,j,1-kk)
    f10sub(nx+kkk,j,1-kk)=buffservice3d(nx+kkk,j,1-kk)
  end forall
  
  ! swap pop11 and pop12
  forall(j=sy:ey)
    buffservice3d(nx+kkk,j,1-kk)=f11sub(nx+kkk,j,1-kk)
    f11sub(nx+kkk,j,1-kk)=f12sub(nx+kkk,j,1-kk)
    f12sub(nx+kkk,j,1-kk)=buffservice3d(nx+kkk,j,1-kk)
  end forall
  
  ! swap pop13 and pop14
  forall(j=sy:ey)
    buffservice3d(nx+kkk,j,1-kk)=f13sub(nx+kkk,j,1-kk)
    f13sub(nx+kkk,j,1-kk)=f14sub(nx+kkk,j,1-kk)
    f14sub(nx+kkk,j,1-kk)=buffservice3d(nx+kkk,j,1-kk)
  end forall
  
  ! swap pop15 and pop16
  forall(j=sy:ey)
    buffservice3d(nx+kkk,j,1-kk)=f15sub(nx+kkk,j,1-kk)
    f15sub(nx+kkk,j,1-kk)=f16sub(nx+kkk,j,1-kk)
    f16sub(nx+kkk,j,1-kk)=buffservice3d(nx+kkk,j,1-kk)
  end forall
  
  ! swap pop17 and pop18
  forall(j=sy:ey)
    buffservice3d(nx+kkk,j,1-kk)=f17sub(nx+kkk,j,1-kk)
    f17sub(nx+kkk,j,1-kk)=f18sub(nx+kkk,j,1-kk)
    f18sub(nx+kkk,j,1-kk)=buffservice3d(nx+kkk,j,1-kk)
  end forall
  
#else
  
  call warning(9,ZERO,latt_name)
  call error(12)
  
#endif
  return
  
 end subroutine apply_bounceback_edge_south_east
 
 subroutine driver_bounceback_edge_south_front(sx,ex)
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the bounceback boundary 
!     condition at the south front edge along the space 
!     defined in input sx,ex
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: sx,ex
  
  call apply_bounceback_edge_south_front(sx,ex,f00R,f01R,f02R,f03R,f04R,&
   f05R,f06R,f07R,f08R,f09R,f10R,f11R,f12R,f13R,f14R,f15R,f16R,f17R, &
   f18R)
  
  if(lsingle_fluid)return
 
  call apply_bounceback_edge_south_front(sx,ex,f00B,f01B,f02B,f03B,f04B,&
   f05B,f06B,f07B,f08B,f09B,f10B,f11B,f12B,f13B,f14B,f15B,f16B,f17B, &
   f18B)
   
  return
  
 end subroutine driver_bounceback_edge_south_front
 
 subroutine apply_bounceback_edge_south_front(sx,ex,f00sub,f01sub, &
  f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub,f09sub,f10sub, &
  f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the bounceback boundary 
!     condition at the south east edge along the space 
!     defined in input sx,ex
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: sx,ex
  real(kind=PRC), allocatable, dimension(:,:,:) :: &
   f00sub,f01sub,f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub, &
   f09sub,f10sub,f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub
  
  integer :: i,j,k
  integer, parameter :: kk = 1
  integer, parameter :: kkk = 1
  
#if LATTICE==319
  
  ! swap pop1 and pop2 x
  forall(i=sx:ex)
    buffservice3d(i,1-kkk,1-kk)=f01sub(i,1-kkk,1-kk)
    f01sub(i,1-kkk,1-kk)=f02sub(i,1-kkk,1-kk)
    f02sub(i,1-kkk,1-kk)=buffservice3d(i,1-kkk,1-kk)
  end forall
  
  ! swap pop3 and pop4 y
  forall(i=sx:ex)
    buffservice3d(i,1-kkk,1-kk)=f03sub(i,1-kkk,1-kk)
    f03sub(i,1-kkk,1-kk)=f04sub(i,1-kkk,1-kk)
    f04sub(i,1-kkk,1-kk)=buffservice3d(i,1-kkk,1-kk)
  end forall
  
  ! swap pop5 and pop6 z
  forall(i=sx:ex)
    buffservice3d(i,1-kkk,1-kk)=f05sub(i,1-kkk,1-kk)
    f05sub(i,1-kkk,1-kk)=f06sub(i,1-kkk,1-kk)
    f06sub(i,1-kkk,1-kk)=buffservice3d(i,1-kkk,1-kk)
  end forall
  
  ! swap pop7 and pop8
  forall(i=sx:ex)
    buffservice3d(i,1-kkk,1-kk)=f07sub(i,1-kkk,1-kk)
    f07sub(i,1-kkk,1-kk)=f08sub(i,1-kkk,1-kk)
    f08sub(i,1-kkk,1-kk)=buffservice3d(i,1-kkk,1-kk)
  end forall
  
  ! swap pop9 and pop10
  forall(i=sx:ex)
    buffservice3d(i,1-kkk,1-kk)=f09sub(i,1-kkk,1-kk)
    f09sub(i,1-kkk,1-kk)=f10sub(i,1-kkk,1-kk)
    f10sub(i,1-kkk,1-kk)=buffservice3d(i,1-kkk,1-kk)
  end forall
  
  ! swap pop11 and pop12
  forall(i=sx:ex)
    buffservice3d(i,1-kkk,1-kk)=f11sub(i,1-kkk,1-kk)
    f11sub(i,1-kkk,1-kk)=f12sub(i,1-kkk,1-kk)
    f12sub(i,1-kkk,1-kk)=buffservice3d(i,1-kkk,1-kk)
  end forall
  
  ! swap pop13 and pop14
  forall(i=sx:ex)
    buffservice3d(i,1-kkk,1-kk)=f13sub(i,1-kkk,1-kk)
    f13sub(i,1-kkk,1-kk)=f14sub(i,1-kkk,1-kk)
    f14sub(i,1-kkk,1-kk)=buffservice3d(i,1-kkk,1-kk)
  end forall
  
  ! swap pop15 and pop16
  forall(i=sx:ex)
    buffservice3d(i,1-kkk,1-kk)=f15sub(i,1-kkk,1-kk)
    f15sub(i,1-kkk,1-kk)=f16sub(i,1-kkk,1-kk)
    f16sub(i,1-kkk,1-kk)=buffservice3d(i,1-kkk,1-kk)
  end forall
  
  ! swap pop17 and pop18
  forall(i=sx:ex)
    buffservice3d(i,1-kkk,1-kk)=f17sub(i,1-kkk,1-kk)
    f17sub(i,1-kkk,1-kk)=f18sub(i,1-kkk,1-kk)
    f18sub(i,1-kkk,1-kk)=buffservice3d(i,1-kkk,1-kk)
  end forall
  
#else
  
  call warning(9,ZERO,latt_name)
  call error(12)
  
#endif
  return
  
 end subroutine apply_bounceback_edge_south_front
 
 subroutine driver_bounceback_edge_south_rear(sx,ex)
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the bounceback boundary 
!     condition at the south rear edge along the space 
!     defined in input sx,ex
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: sx,ex
  
  call apply_bounceback_edge_south_rear(sx,ex,f00R,f01R,f02R,f03R,f04R,&
   f05R,f06R,f07R,f08R,f09R,f10R,f11R,f12R,f13R,f14R,f15R,f16R,f17R, &
   f18R)
  
  if(lsingle_fluid)return
 
  call apply_bounceback_edge_south_rear(sx,ex,f00B,f01B,f02B,f03B,f04B,&
   f05B,f06B,f07B,f08B,f09B,f10B,f11B,f12B,f13B,f14B,f15B,f16B,f17B, &
   f18B)
   
  return
  
 end subroutine driver_bounceback_edge_south_rear
 
 subroutine apply_bounceback_edge_south_rear(sx,ex,f00sub,f01sub, &
  f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub,f09sub,f10sub, &
  f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the bounceback boundary 
!     condition at the south rear edge along the space 
!     defined in input sx,ex
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: sx,ex
  real(kind=PRC), allocatable, dimension(:,:,:) :: &
   f00sub,f01sub,f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub, &
   f09sub,f10sub,f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub
  
  integer :: i,j,k
  integer, parameter :: kk = 1
  integer, parameter :: kkk = 1
  
#if LATTICE==319
  
  ! swap pop1 and pop2 x
  forall(i=sx:ex)
    buffservice3d(i,ny+kkk,1-kk)=f01sub(i,ny+kkk,1-kk)
    f01sub(i,ny+kkk,1-kk)=f02sub(i,ny+kkk,1-kk)
    f02sub(i,ny+kkk,1-kk)=buffservice3d(i,ny+kkk,1-kk)
  end forall
  
  ! swap pop3 and pop4 y
  forall(i=sx:ex)
    buffservice3d(i,ny+kkk,1-kk)=f03sub(i,ny+kkk,1-kk)
    f03sub(i,ny+kkk,1-kk)=f04sub(i,ny+kkk,1-kk)
    f04sub(i,ny+kkk,1-kk)=buffservice3d(i,ny+kkk,1-kk)
  end forall
  
  ! swap pop5 and pop6 z
  forall(i=sx:ex)
    buffservice3d(i,ny+kkk,1-kk)=f05sub(i,ny+kkk,1-kk)
    f05sub(i,ny+kkk,1-kk)=f06sub(i,ny+kkk,1-kk)
    f06sub(i,ny+kkk,1-kk)=buffservice3d(i,ny+kkk,1-kk)
  end forall
  
  ! swap pop7 and pop8
  forall(i=sx:ex)
    buffservice3d(i,ny+kkk,1-kk)=f07sub(i,ny+kkk,1-kk)
    f07sub(i,ny+kkk,1-kk)=f08sub(i,ny+kkk,1-kk)
    f08sub(i,ny+kkk,1-kk)=buffservice3d(i,ny+kkk,1-kk)
  end forall
  
  ! swap pop9 and pop10
  forall(i=sx:ex)
    buffservice3d(i,ny+kkk,1-kk)=f09sub(i,ny+kkk,1-kk)
    f09sub(i,ny+kkk,1-kk)=f10sub(i,ny+kkk,1-kk)
    f10sub(i,ny+kkk,1-kk)=buffservice3d(i,ny+kkk,1-kk)
  end forall
  
  ! swap pop11 and pop12
  forall(i=sx:ex)
    buffservice3d(i,ny+kkk,1-kk)=f11sub(i,ny+kkk,1-kk)
    f11sub(i,ny+kkk,1-kk)=f12sub(i,ny+kkk,1-kk)
    f12sub(i,ny+kkk,1-kk)=buffservice3d(i,ny+kkk,1-kk)
  end forall
  
  ! swap pop13 and pop14
  forall(i=sx:ex)
    buffservice3d(i,ny+kkk,1-kk)=f13sub(i,ny+kkk,1-kk)
    f13sub(i,ny+kkk,1-kk)=f14sub(i,ny+kkk,1-kk)
    f14sub(i,ny+kkk,1-kk)=buffservice3d(i,ny+kkk,1-kk)
  end forall
  
  ! swap pop15 and pop16
  forall(i=sx:ex)
    buffservice3d(i,ny+kkk,1-kk)=f15sub(i,ny+kkk,1-kk)
    f15sub(i,ny+kkk,1-kk)=f16sub(i,ny+kkk,1-kk)
    f16sub(i,ny+kkk,1-kk)=buffservice3d(i,ny+kkk,1-kk)
  end forall
  
  ! swap pop17 and pop18
  forall(i=sx:ex)
    buffservice3d(i,ny+kkk,1-kk)=f17sub(i,ny+kkk,1-kk)
    f17sub(i,ny+kkk,1-kk)=f18sub(i,ny+kkk,1-kk)
    f18sub(i,ny+kkk,1-kk)=buffservice3d(i,ny+kkk,1-kk)
  end forall
  
#else
  
  call warning(9,ZERO,latt_name)
  call error(12)
  
#endif
  return
  
 end subroutine apply_bounceback_edge_south_rear
 
 subroutine driver_bounceback_edge_south_west(sy,ey)
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the bounceback boundary 
!     condition at the south west edge along the space 
!     defined in input sy,ey
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: sy,ey
  
  call apply_bounceback_edge_south_west(sy,ey,f00R,f01R,f02R,f03R,f04R,&
   f05R,f06R,f07R,f08R,f09R,f10R,f11R,f12R,f13R,f14R,f15R,f16R,f17R, &
   f18R)
  
  if(lsingle_fluid)return
 
  call apply_bounceback_edge_south_west(sy,ey,f00B,f01B,f02B,f03B,f04B,&
   f05B,f06B,f07B,f08B,f09B,f10B,f11B,f12B,f13B,f14B,f15B,f16B,f17B, &
   f18B)
   
  return
  
 end subroutine driver_bounceback_edge_south_west
 
 subroutine apply_bounceback_edge_south_west(sy,ey,f00sub,f01sub, &
  f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub,f09sub,f10sub, &
  f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the bounceback boundary 
!     condition at the south rear edge along the space 
!     defined in input sy,ey
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: sy,ey
  real(kind=PRC), allocatable, dimension(:,:,:) :: &
   f00sub,f01sub,f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub, &
   f09sub,f10sub,f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub
  
  integer :: i,j,k
  integer, parameter :: kk = 1
  integer, parameter :: kkk = 1
  
#if LATTICE==319
  
  ! swap pop1 and pop2 x
  forall(j=sy:ey)
    buffservice3d(1-kkk,j,1-kk)=f01sub(1-kkk,j,1-kk)
    f01sub(1-kkk,j,1-kk)=f02sub(1-kkk,j,1-kk)
    f02sub(1-kkk,j,1-kk)=buffservice3d(1-kkk,j,1-kk)
  end forall
  
  ! swap pop3 and pop4 y
  forall(j=sy:ey)
    buffservice3d(1-kkk,j,1-kk)=f03sub(1-kkk,j,1-kk)
    f03sub(1-kkk,j,1-kk)=f04sub(1-kkk,j,1-kk)
    f04sub(1-kkk,j,1-kk)=buffservice3d(1-kkk,j,1-kk)
  end forall
  
  ! swap pop5 and pop6 z
  forall(j=sy:ey)
    buffservice3d(1-kkk,j,1-kk)=f05sub(1-kkk,j,1-kk)
    f05sub(1-kkk,j,1-kk)=f06sub(1-kkk,j,1-kk)
    f06sub(1-kkk,j,1-kk)=buffservice3d(1-kkk,j,1-kk)
  end forall
  
  ! swap pop7 and pop8
  forall(j=sy:ey)
    buffservice3d(1-kkk,j,1-kk)=f07sub(1-kkk,j,1-kk)
    f07sub(1-kkk,j,1-kk)=f08sub(1-kkk,j,1-kk)
    f08sub(1-kkk,j,1-kk)=buffservice3d(1-kkk,j,1-kk)
  end forall
  
  ! swap pop9 and pop10
  forall(j=sy:ey)
    buffservice3d(1-kkk,j,1-kk)=f09sub(1-kkk,j,1-kk)
    f09sub(1-kkk,j,1-kk)=f10sub(1-kkk,j,1-kk)
    f10sub(1-kkk,j,1-kk)=buffservice3d(1-kkk,j,1-kk)
  end forall
  
  ! swap pop11 and pop12
  forall(j=sy:ey)
    buffservice3d(1-kkk,j,1-kk)=f11sub(1-kkk,j,1-kk)
    f11sub(1-kkk,j,1-kk)=f12sub(1-kkk,j,1-kk)
    f12sub(1-kkk,j,1-kk)=buffservice3d(1-kkk,j,1-kk)
  end forall
  
  ! swap pop13 and pop14
  forall(j=sy:ey)
    buffservice3d(1-kkk,j,1-kk)=f13sub(1-kkk,j,1-kk)
    f13sub(1-kkk,j,1-kk)=f14sub(1-kkk,j,1-kk)
    f14sub(1-kkk,j,1-kk)=buffservice3d(1-kkk,j,1-kk)
  end forall
  
  ! swap pop15 and pop16
  forall(j=sy:ey)
    buffservice3d(1-kkk,j,1-kk)=f15sub(1-kkk,j,1-kk)
    f15sub(1-kkk,j,1-kk)=f16sub(1-kkk,j,1-kk)
    f16sub(1-kkk,j,1-kk)=buffservice3d(1-kkk,j,1-kk)
  end forall
  
  ! swap pop17 and pop18
  forall(j=sy:ey)
    buffservice3d(1-kkk,j,1-kk)=f17sub(1-kkk,j,1-kk)
    f17sub(1-kkk,j,1-kk)=f18sub(1-kkk,j,1-kk)
    f18sub(1-kkk,j,1-kk)=buffservice3d(1-kkk,j,1-kk)
  end forall
  
#else
  
  call warning(9,ZERO,latt_name)
  call error(12)
  
#endif
  return
  
 end subroutine apply_bounceback_edge_south_west
 
 subroutine driver_bounceback_corner_north_east_front
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the bounceback boundary 
!     condition at the north east front corner
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  call apply_bounceback_corner_north_east_front(f00R,f01R,f02R,f03R,f04R,&
   f05R,f06R,f07R,f08R,f09R,f10R,f11R,f12R,f13R,f14R,f15R,f16R,f17R, &
   f18R)
  
  if(lsingle_fluid)return
 
  call apply_bounceback_corner_north_east_front(f00B,f01B,f02B,f03B,f04B,&
   f05B,f06B,f07B,f08B,f09B,f10B,f11B,f12B,f13B,f14B,f15B,f16B,f17B, &
   f18B)
   
  return
  
 end subroutine driver_bounceback_corner_north_east_front
 
 subroutine apply_bounceback_corner_north_east_front(f00sub,f01sub, &
  f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub,f09sub,f10sub, &
  f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the bounceback boundary 
!     condition at the north east front corner
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), allocatable, dimension(:,:,:) :: &
   f00sub,f01sub,f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub, &
   f09sub,f10sub,f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub
  
  integer :: i,j,k
  integer, parameter :: kk = 1
  integer, parameter :: kkk = 1
  integer, parameter :: kkkk = 1
  
#if LATTICE==319
  

  
  ! swap pop7 and pop8

  buffservice3d(nx+kk,1-kkk,nz+kkkk)=f07sub(nx+kk,1-kkk,nz+kkkk)
  f07sub(nx+kk,1-kkk,nz+kkkk)=f08sub(nx+kk,1-kkk,nz+kkkk)
  f08sub(nx+kk,1-kkk,nz+kkkk)=buffservice3d(nx+kk,1-kkk,nz+kkkk)

  
  ! swap pop9 and pop10

  buffservice3d(nx+kk,1-kkk,nz+kkkk)=f09sub(nx+kk,1-kkk,nz+kkkk)
  f09sub(nx+kk,1-kkk,nz+kkkk)=f10sub(nx+kk,1-kkk,nz+kkkk)
  f10sub(nx+kk,1-kkk,nz+kkkk)=buffservice3d(nx+kk,1-kkk,nz+kkkk)

  
  ! swap pop11 and pop12

  buffservice3d(nx+kk,1-kkk,nz+kkkk)=f11sub(nx+kk,1-kkk,nz+kkkk)
  f11sub(nx+kk,1-kkk,nz+kkkk)=f12sub(nx+kk,1-kkk,nz+kkkk)
  f12sub(nx+kk,1-kkk,nz+kkkk)=buffservice3d(nx+kk,1-kkk,nz+kkkk)

  
  ! swap pop13 and pop14

  buffservice3d(nx+kk,1-kkk,nz+kkkk)=f13sub(nx+kk,1-kkk,nz+kkkk)
  f13sub(nx+kk,1-kkk,nz+kkkk)=f14sub(nx+kk,1-kkk,nz+kkkk)
  f14sub(nx+kk,1-kkk,nz+kkkk)=buffservice3d(nx+kk,1-kkk,nz+kkkk)

  
  ! swap pop15 and pop16

  buffservice3d(nx+kk,1-kkk,nz+kkkk)=f15sub(nx+kk,1-kkk,nz+kkkk)
  f15sub(nx+kk,1-kkk,nz+kkkk)=f16sub(nx+kk,1-kkk,nz+kkkk)
  f16sub(nx+kk,1-kkk,nz+kkkk)=buffservice3d(nx+kk,1-kkk,nz+kkkk)

  
  ! swap pop17 and pop18

  buffservice3d(nx+kk,1-kkk,nz+kkkk)=f17sub(nx+kk,1-kkk,nz+kkkk)
  f17sub(nx+kk,1-kkk,nz+kkkk)=f18sub(nx+kk,1-kkk,nz+kkkk)
  f18sub(nx+kk,1-kkk,nz+kkkk)=buffservice3d(nx+kk,1-kkk,nz+kkkk)
  
#else
  
  call warning(9,ZERO,latt_name)
  call error(12)
  
#endif
  return
  
 end subroutine apply_bounceback_corner_north_east_front
 
 subroutine driver_bounceback_corner_north_east_rear
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the bounceback boundary 
!     condition at the north east rear corner
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  call apply_bounceback_corner_north_east_rear(f00R,f01R,f02R,f03R,f04R,&
   f05R,f06R,f07R,f08R,f09R,f10R,f11R,f12R,f13R,f14R,f15R,f16R,f17R, &
   f18R)
  
  if(lsingle_fluid)return
 
  call apply_bounceback_corner_north_east_rear(f00B,f01B,f02B,f03B,f04B,&
   f05B,f06B,f07B,f08B,f09B,f10B,f11B,f12B,f13B,f14B,f15B,f16B,f17B, &
   f18B)
   
  return
  
 end subroutine driver_bounceback_corner_north_east_rear
 
 subroutine apply_bounceback_corner_north_east_rear(f00sub,f01sub, &
  f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub,f09sub,f10sub, &
  f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the bounceback boundary 
!     condition at the north east rear corner
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), allocatable, dimension(:,:,:) :: &
   f00sub,f01sub,f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub, &
   f09sub,f10sub,f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub
  
  integer :: i,j,k
  integer, parameter :: kk = 1
  integer, parameter :: kkk = 1
  integer, parameter :: kkkk = 1
  
#if LATTICE==319
  
  ! swap pop7 and pop8

  buffservice3d(nx+kk,ny+kkk,nz+kkkk)=f07sub(nx+kk,ny+kkk,nz+kkkk)
  f07sub(nx+kk,ny+kkk,nz+kkkk)=f08sub(nx+kk,ny+kkk,nz+kkkk)
  f08sub(nx+kk,ny+kkk,nz+kkkk)=buffservice3d(nx+kk,ny+kkk,nz+kkkk)

  
  ! swap pop9 and pop10

  buffservice3d(nx+kk,ny+kkk,nz+kkkk)=f09sub(nx+kk,ny+kkk,nz+kkkk)
  f09sub(nx+kk,ny+kkk,nz+kkkk)=f10sub(nx+kk,ny+kkk,nz+kkkk)
  f10sub(nx+kk,ny+kkk,nz+kkkk)=buffservice3d(nx+kk,ny+kkk,nz+kkkk)

  
  ! swap pop11 and pop12

  buffservice3d(nx+kk,ny+kkk,nz+kkkk)=f11sub(nx+kk,ny+kkk,nz+kkkk)
  f11sub(nx+kk,ny+kkk,nz+kkkk)=f12sub(nx+kk,ny+kkk,nz+kkkk)
  f12sub(nx+kk,ny+kkk,nz+kkkk)=buffservice3d(nx+kk,ny+kkk,nz+kkkk)

  
  ! swap pop13 and pop14

  buffservice3d(nx+kk,ny+kkk,nz+kkkk)=f13sub(nx+kk,ny+kkk,nz+kkkk)
  f13sub(nx+kk,ny+kkk,nz+kkkk)=f14sub(nx+kk,ny+kkk,nz+kkkk)
  f14sub(nx+kk,ny+kkk,nz+kkkk)=buffservice3d(nx+kk,ny+kkk,nz+kkkk)

  
  ! swap pop15 and pop16

  buffservice3d(nx+kk,ny+kkk,nz+kkkk)=f15sub(nx+kk,ny+kkk,nz+kkkk)
  f15sub(nx+kk,ny+kkk,nz+kkkk)=f16sub(nx+kk,ny+kkk,nz+kkkk)
  f16sub(nx+kk,ny+kkk,nz+kkkk)=buffservice3d(nx+kk,ny+kkk,nz+kkkk)

  
  ! swap pop17 and pop18

  buffservice3d(nx+kk,ny+kkk,nz+kkkk)=f17sub(nx+kk,ny+kkk,nz+kkkk)
  f17sub(nx+kk,ny+kkk,nz+kkkk)=f18sub(nx+kk,ny+kkk,nz+kkkk)
  f18sub(nx+kk,ny+kkk,nz+kkkk)=buffservice3d(nx+kk,ny+kkk,nz+kkkk)
  
#else
  
  call warning(9,ZERO,latt_name)
  call error(12)
  
#endif
  return
  
 end subroutine apply_bounceback_corner_north_east_rear
 
 subroutine driver_bounceback_corner_north_west_rear
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the bounceback boundary 
!     condition at the north west rear corner
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  call apply_bounceback_corner_north_west_rear(f00R,f01R,f02R,f03R,f04R,&
   f05R,f06R,f07R,f08R,f09R,f10R,f11R,f12R,f13R,f14R,f15R,f16R,f17R, &
   f18R)
  
  if(lsingle_fluid)return
 
  call apply_bounceback_corner_north_west_rear(f00B,f01B,f02B,f03B,f04B,&
   f05B,f06B,f07B,f08B,f09B,f10B,f11B,f12B,f13B,f14B,f15B,f16B,f17B, &
   f18B)
   
  return
  
 end subroutine driver_bounceback_corner_north_west_rear
 
 subroutine apply_bounceback_corner_north_west_rear(f00sub,f01sub, &
  f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub,f09sub,f10sub, &
  f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the bounceback boundary 
!     condition at the north west rear corner
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), allocatable, dimension(:,:,:) :: &
   f00sub,f01sub,f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub, &
   f09sub,f10sub,f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub
  
  integer :: i,j,k
  integer, parameter :: kk = 1
  integer, parameter :: kkk = 1
  integer, parameter :: kkkk = 1
  
#if LATTICE==319
  
  ! swap pop7 and pop8

  buffservice3d(1-kk,ny+kkk,nz+kkkk)=f07sub(1-kk,ny+kkk,nz+kkkk)
  f07sub(1-kk,ny+kkk,nz+kkkk)=f08sub(1-kk,ny+kkk,nz+kkkk)
  f08sub(1-kk,ny+kkk,nz+kkkk)=buffservice3d(1-kk,ny+kkk,nz+kkkk)

  
  ! swap pop9 and pop10

  buffservice3d(1-kk,ny+kkk,nz+kkkk)=f09sub(1-kk,ny+kkk,nz+kkkk)
  f09sub(1-kk,ny+kkk,nz+kkkk)=f10sub(1-kk,ny+kkk,nz+kkkk)
  f10sub(1-kk,ny+kkk,nz+kkkk)=buffservice3d(1-kk,ny+kkk,nz+kkkk)

  
  ! swap pop11 and pop12

  buffservice3d(1-kk,ny+kkk,nz+kkkk)=f11sub(1-kk,ny+kkk,nz+kkkk)
  f11sub(1-kk,ny+kkk,nz+kkkk)=f12sub(1-kk,ny+kkk,nz+kkkk)
  f12sub(1-kk,ny+kkk,nz+kkkk)=buffservice3d(1-kk,ny+kkk,nz+kkkk)

  
  ! swap pop13 and pop14

  buffservice3d(1-kk,ny+kkk,nz+kkkk)=f13sub(1-kk,ny+kkk,nz+kkkk)
  f13sub(1-kk,ny+kkk,nz+kkkk)=f14sub(1-kk,ny+kkk,nz+kkkk)
  f14sub(1-kk,ny+kkk,nz+kkkk)=buffservice3d(1-kk,ny+kkk,nz+kkkk)

  
  ! swap pop15 and pop16

  buffservice3d(1-kk,ny+kkk,nz+kkkk)=f15sub(1-kk,ny+kkk,nz+kkkk)
  f15sub(1-kk,ny+kkk,nz+kkkk)=f16sub(1-kk,ny+kkk,nz+kkkk)
  f16sub(1-kk,ny+kkk,nz+kkkk)=buffservice3d(1-kk,ny+kkk,nz+kkkk)

  
  ! swap pop17 and pop18

  buffservice3d(1-kk,ny+kkk,nz+kkkk)=f17sub(1-kk,ny+kkk,nz+kkkk)
  f17sub(1-kk,ny+kkk,nz+kkkk)=f18sub(1-kk,ny+kkk,nz+kkkk)
  f18sub(1-kk,ny+kkk,nz+kkkk)=buffservice3d(1-kk,ny+kkk,nz+kkkk)
  
#else
  
  call warning(9,ZERO,latt_name)
  call error(12)
  
#endif
  return
  
 end subroutine apply_bounceback_corner_north_west_rear
 
 subroutine driver_bounceback_corner_north_west_front
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the bounceback boundary 
!     condition at the north west front corner
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  call apply_bounceback_corner_north_west_front(f00R,f01R,f02R,f03R,f04R,&
   f05R,f06R,f07R,f08R,f09R,f10R,f11R,f12R,f13R,f14R,f15R,f16R,f17R, &
   f18R)
  
  if(lsingle_fluid)return
 
  call apply_bounceback_corner_north_west_front(f00B,f01B,f02B,f03B,f04B,&
   f05B,f06B,f07B,f08B,f09B,f10B,f11B,f12B,f13B,f14B,f15B,f16B,f17B, &
   f18B)
   
  return
  
 end subroutine driver_bounceback_corner_north_west_front
 
 subroutine apply_bounceback_corner_north_west_front(f00sub,f01sub, &
  f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub,f09sub,f10sub, &
  f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the bounceback boundary 
!     condition at the north west front corner
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), allocatable, dimension(:,:,:) :: &
   f00sub,f01sub,f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub, &
   f09sub,f10sub,f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub
  
  integer :: i,j,k
  integer, parameter :: kk = 1
  integer, parameter :: kkk = 1
  integer, parameter :: kkkk = 1
  
#if LATTICE==319
  
  ! swap pop7 and pop8

  buffservice3d(1-kk,1-kkk,nz+kkkk)=f07sub(1-kk,1-kkk,nz+kkkk)
  f07sub(1-kk,1-kkk,nz+kkkk)=f08sub(1-kk,1-kkk,nz+kkkk)
  f08sub(1-kk,1-kkk,nz+kkkk)=buffservice3d(1-kk,1-kkk,nz+kkkk)

  
  ! swap pop9 and pop10

  buffservice3d(1-kk,1-kkk,nz+kkkk)=f09sub(1-kk,1-kkk,nz+kkkk)
  f09sub(1-kk,1-kkk,nz+kkkk)=f10sub(1-kk,1-kkk,nz+kkkk)
  f10sub(1-kk,1-kkk,nz+kkkk)=buffservice3d(1-kk,1-kkk,nz+kkkk)

  
  ! swap pop11 and pop12

  buffservice3d(1-kk,1-kkk,nz+kkkk)=f11sub(1-kk,1-kkk,nz+kkkk)
  f11sub(1-kk,1-kkk,nz+kkkk)=f12sub(1-kk,1-kkk,nz+kkkk)
  f12sub(1-kk,1-kkk,nz+kkkk)=buffservice3d(1-kk,1-kkk,nz+kkkk)

  
  ! swap pop13 and pop14

  buffservice3d(1-kk,1-kkk,nz+kkkk)=f13sub(1-kk,1-kkk,nz+kkkk)
  f13sub(1-kk,1-kkk,nz+kkkk)=f14sub(1-kk,1-kkk,nz+kkkk)
  f14sub(1-kk,1-kkk,nz+kkkk)=buffservice3d(1-kk,1-kkk,nz+kkkk)

  
  ! swap pop15 and pop16

  buffservice3d(1-kk,1-kkk,nz+kkkk)=f15sub(1-kk,1-kkk,nz+kkkk)
  f15sub(1-kk,1-kkk,nz+kkkk)=f16sub(1-kk,1-kkk,nz+kkkk)
  f16sub(1-kk,1-kkk,nz+kkkk)=buffservice3d(1-kk,1-kkk,nz+kkkk)

  
  ! swap pop17 and pop18

  buffservice3d(1-kk,1-kkk,nz+kkkk)=f17sub(1-kk,1-kkk,nz+kkkk)
  f17sub(1-kk,1-kkk,nz+kkkk)=f18sub(1-kk,1-kkk,nz+kkkk)
  f18sub(1-kk,1-kkk,nz+kkkk)=buffservice3d(1-kk,1-kkk,nz+kkkk)
  
#else
  
  call warning(9,ZERO,latt_name)
  call error(12)
  
#endif
  return
  
 end subroutine apply_bounceback_corner_north_west_front
 
 subroutine driver_bounceback_corner_south_east_front
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the bounceback boundary 
!     condition at the south east front corner
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  call apply_bounceback_corner_south_east_front(f00R,f01R,f02R,f03R,f04R,&
   f05R,f06R,f07R,f08R,f09R,f10R,f11R,f12R,f13R,f14R,f15R,f16R,f17R, &
   f18R)
  
  if(lsingle_fluid)return
 
  call apply_bounceback_corner_south_east_front(f00B,f01B,f02B,f03B,f04B,&
   f05B,f06B,f07B,f08B,f09B,f10B,f11B,f12B,f13B,f14B,f15B,f16B,f17B, &
   f18B)
   
  return
  
 end subroutine driver_bounceback_corner_south_east_front
 
 subroutine apply_bounceback_corner_south_east_front(f00sub,f01sub, &
  f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub,f09sub,f10sub, &
  f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the bounceback boundary 
!     condition at the south east front corner
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), allocatable, dimension(:,:,:) :: &
   f00sub,f01sub,f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub, &
   f09sub,f10sub,f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub
  
  integer :: i,j,k
  integer, parameter :: kk = 1
  integer, parameter :: kkk = 1
  integer, parameter :: kkkk = 1
  
#if LATTICE==319
  
  
  ! swap pop7 and pop8

  buffservice3d(nx+kk,1-kkk,1-kkkk)=f07sub(nx+kk,1-kkk,1-kkkk)
  f07sub(nx+kk,1-kkk,1-kkkk)=f08sub(nx+kk,1-kkk,1-kkkk)
  f08sub(nx+kk,1-kkk,1-kkkk)=buffservice3d(nx+kk,1-kkk,1-kkkk)

  
  ! swap pop9 and pop10

  buffservice3d(nx+kk,1-kkk,1-kkkk)=f09sub(nx+kk,1-kkk,1-kkkk)
  f09sub(nx+kk,1-kkk,1-kkkk)=f10sub(nx+kk,1-kkk,1-kkkk)
  f10sub(nx+kk,1-kkk,1-kkkk)=buffservice3d(nx+kk,1-kkk,1-kkkk)

  
  ! swap pop11 and pop12

  buffservice3d(nx+kk,1-kkk,1-kkkk)=f11sub(nx+kk,1-kkk,1-kkkk)
  f11sub(nx+kk,1-kkk,1-kkkk)=f12sub(nx+kk,1-kkk,1-kkkk)
  f12sub(nx+kk,1-kkk,1-kkkk)=buffservice3d(nx+kk,1-kkk,1-kkkk)

  
  ! swap pop13 and pop14

  buffservice3d(nx+kk,1-kkk,1-kkkk)=f13sub(nx+kk,1-kkk,1-kkkk)
  f13sub(nx+kk,1-kkk,1-kkkk)=f14sub(nx+kk,1-kkk,1-kkkk)
  f14sub(nx+kk,1-kkk,1-kkkk)=buffservice3d(nx+kk,1-kkk,1-kkkk)

  
  ! swap pop15 and pop16

  buffservice3d(nx+kk,1-kkk,1-kkkk)=f15sub(nx+kk,1-kkk,1-kkkk)
  f15sub(nx+kk,1-kkk,1-kkkk)=f16sub(nx+kk,1-kkk,1-kkkk)
  f16sub(nx+kk,1-kkk,1-kkkk)=buffservice3d(nx+kk,1-kkk,1-kkkk)

  
  ! swap pop17 and pop18

  buffservice3d(nx+kk,1-kkk,1-kkkk)=f17sub(nx+kk,1-kkk,1-kkkk)
  f17sub(nx+kk,1-kkk,1-kkkk)=f18sub(nx+kk,1-kkk,1-kkkk)
  f18sub(nx+kk,1-kkk,1-kkkk)=buffservice3d(nx+kk,1-kkk,1-kkkk)
  
#else
  
  call warning(9,ZERO,latt_name)
  call error(12)
  
#endif
  return
  
 end subroutine apply_bounceback_corner_south_east_front
 
 subroutine driver_bounceback_corner_south_west_front
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the bounceback boundary 
!     condition at the south west front corner
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  call apply_bounceback_corner_south_west_front(f00R,f01R,f02R,f03R,f04R,&
   f05R,f06R,f07R,f08R,f09R,f10R,f11R,f12R,f13R,f14R,f15R,f16R,f17R, &
   f18R)
  
  if(lsingle_fluid)return
 
  call apply_bounceback_corner_south_west_front(f00B,f01B,f02B,f03B,f04B,&
   f05B,f06B,f07B,f08B,f09B,f10B,f11B,f12B,f13B,f14B,f15B,f16B,f17B, &
   f18B)
   
  return
  
 end subroutine driver_bounceback_corner_south_west_front
 
 subroutine apply_bounceback_corner_south_west_front(f00sub,f01sub, &
  f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub,f09sub,f10sub, &
  f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the bounceback boundary 
!     condition at the south west front corner
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), allocatable, dimension(:,:,:) :: &
   f00sub,f01sub,f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub, &
   f09sub,f10sub,f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub
  
  integer :: i,j,k
  integer, parameter :: kk = 1
  integer, parameter :: kkk = 1
  integer, parameter :: kkkk = 1
  
#if LATTICE==319
  
  ! swap pop7 and pop8

  buffservice3d(1-kk,1-kkk,1-kkkk)=f07sub(1-kk,1-kkk,1-kkkk)
  f07sub(1-kk,1-kkk,1-kkkk)=f08sub(1-kk,1-kkk,1-kkkk)
  f08sub(1-kk,1-kkk,1-kkkk)=buffservice3d(1-kk,1-kkk,1-kkkk)

  
  ! swap pop9 and pop10

  buffservice3d(1-kk,1-kkk,1-kkkk)=f09sub(1-kk,1-kkk,1-kkkk)
  f09sub(1-kk,1-kkk,1-kkkk)=f10sub(1-kk,1-kkk,1-kkkk)
  f10sub(1-kk,1-kkk,1-kkkk)=buffservice3d(1-kk,1-kkk,1-kkkk)

  
  ! swap pop11 and pop12

  buffservice3d(1-kk,1-kkk,1-kkkk)=f11sub(1-kk,1-kkk,1-kkkk)
  f11sub(1-kk,1-kkk,1-kkkk)=f12sub(1-kk,1-kkk,1-kkkk)
  f12sub(1-kk,1-kkk,1-kkkk)=buffservice3d(1-kk,1-kkk,1-kkkk)

  
  ! swap pop13 and pop14

  buffservice3d(1-kk,1-kkk,1-kkkk)=f13sub(1-kk,1-kkk,1-kkkk)
  f13sub(1-kk,1-kkk,1-kkkk)=f14sub(1-kk,1-kkk,1-kkkk)
  f14sub(1-kk,1-kkk,1-kkkk)=buffservice3d(1-kk,1-kkk,1-kkkk)

  
  ! swap pop15 and pop16

  buffservice3d(1-kk,1-kkk,1-kkkk)=f15sub(1-kk,1-kkk,1-kkkk)
  f15sub(1-kk,1-kkk,1-kkkk)=f16sub(1-kk,1-kkk,1-kkkk)
  f16sub(1-kk,1-kkk,1-kkkk)=buffservice3d(1-kk,1-kkk,1-kkkk)

  
  ! swap pop17 and pop18

  buffservice3d(1-kk,1-kkk,1-kkkk)=f17sub(1-kk,1-kkk,1-kkkk)
  f17sub(1-kk,1-kkk,1-kkkk)=f18sub(1-kk,1-kkk,1-kkkk)
  f18sub(1-kk,1-kkk,1-kkkk)=buffservice3d(1-kk,1-kkk,1-kkkk)
  
#else
  
  call warning(9,ZERO,latt_name)
  call error(12)
  
#endif
  return
  
 end subroutine apply_bounceback_corner_south_west_front
 
 subroutine driver_bounceback_corner_south_west_rear
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the bounceback boundary 
!     condition at the south west rear corner
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  call apply_bounceback_corner_south_west_rear(f00R,f01R,f02R,f03R,f04R,&
   f05R,f06R,f07R,f08R,f09R,f10R,f11R,f12R,f13R,f14R,f15R,f16R,f17R, &
   f18R)
  
  if(lsingle_fluid)return
 
  call apply_bounceback_corner_south_west_rear(f00B,f01B,f02B,f03B,f04B,&
   f05B,f06B,f07B,f08B,f09B,f10B,f11B,f12B,f13B,f14B,f15B,f16B,f17B, &
   f18B)
   
  return
  
 end subroutine driver_bounceback_corner_south_west_rear
 
 subroutine apply_bounceback_corner_south_west_rear(f00sub,f01sub, &
  f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub,f09sub,f10sub, &
  f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the bounceback boundary 
!     condition at the south west rear corner
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), allocatable, dimension(:,:,:) :: &
   f00sub,f01sub,f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub, &
   f09sub,f10sub,f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub
  
  integer :: i,j,k
  integer, parameter :: kk = 1
  integer, parameter :: kkk = 1
  integer, parameter :: kkkk = 1
  
#if LATTICE==319
  
  ! swap pop7 and pop8

  buffservice3d(1-kk,ny+kkk,1-kkkk)=f07sub(1-kk,ny+kkk,1-kkkk)
  f07sub(1-kk,ny+kkk,1-kkkk)=f08sub(1-kk,ny+kkk,1-kkkk)
  f08sub(1-kk,ny+kkk,1-kkkk)=buffservice3d(1-kk,ny+kkk,1-kkkk)

  
  ! swap pop9 and pop10

  buffservice3d(1-kk,ny+kkk,1-kkkk)=f09sub(1-kk,ny+kkk,1-kkkk)
  f09sub(1-kk,ny+kkk,1-kkkk)=f10sub(1-kk,ny+kkk,1-kkkk)
  f10sub(1-kk,ny+kkk,1-kkkk)=buffservice3d(1-kk,ny+kkk,1-kkkk)

  
  ! swap pop11 and pop12

  buffservice3d(1-kk,ny+kkk,1-kkkk)=f11sub(1-kk,ny+kkk,1-kkkk)
  f11sub(1-kk,ny+kkk,1-kkkk)=f12sub(1-kk,ny+kkk,1-kkkk)
  f12sub(1-kk,ny+kkk,1-kkkk)=buffservice3d(1-kk,ny+kkk,1-kkkk)

  
  ! swap pop13 and pop14

  buffservice3d(1-kk,ny+kkk,1-kkkk)=f13sub(1-kk,ny+kkk,1-kkkk)
  f13sub(1-kk,ny+kkk,1-kkkk)=f14sub(1-kk,ny+kkk,1-kkkk)
  f14sub(1-kk,ny+kkk,1-kkkk)=buffservice3d(1-kk,ny+kkk,1-kkkk)

  
  ! swap pop15 and pop16

  buffservice3d(1-kk,ny+kkk,1-kkkk)=f15sub(1-kk,ny+kkk,1-kkkk)
  f15sub(1-kk,ny+kkk,1-kkkk)=f16sub(1-kk,ny+kkk,1-kkkk)
  f16sub(1-kk,ny+kkk,1-kkkk)=buffservice3d(1-kk,ny+kkk,1-kkkk)

  
  ! swap pop17 and pop18

  buffservice3d(1-kk,ny+kkk,1-kkkk)=f17sub(1-kk,ny+kkk,1-kkkk)
  f17sub(1-kk,ny+kkk,1-kkkk)=f18sub(1-kk,ny+kkk,1-kkkk)
  f18sub(1-kk,ny+kkk,1-kkkk)=buffservice3d(1-kk,ny+kkk,1-kkkk)
  
#else
  
  call warning(9,ZERO,latt_name)
  call error(12)
  
#endif
  return
  
 end subroutine apply_bounceback_corner_south_west_rear
 
 subroutine driver_bounceback_corner_south_east_rear
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the bounceback boundary 
!     condition at the south east rear corner
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  call apply_bounceback_corner_south_east_rear(f00R,f01R,f02R,f03R,f04R,&
   f05R,f06R,f07R,f08R,f09R,f10R,f11R,f12R,f13R,f14R,f15R,f16R,f17R, &
   f18R)
  
  if(lsingle_fluid)return
 
  call apply_bounceback_corner_south_east_rear(f00B,f01B,f02B,f03B,f04B,&
   f05B,f06B,f07B,f08B,f09B,f10B,f11B,f12B,f13B,f14B,f15B,f16B,f17B, &
   f18B)
   
  return
  
 end subroutine driver_bounceback_corner_south_east_rear
 
 subroutine apply_bounceback_corner_south_east_rear(f00sub,f01sub, &
  f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub,f09sub,f10sub, &
  f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the bounceback boundary 
!     condition at the south east rear corner
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), allocatable, dimension(:,:,:) :: &
   f00sub,f01sub,f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub, &
   f09sub,f10sub,f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub
  
  integer :: i,j,k
  integer, parameter :: kk = 1
  integer, parameter :: kkk = 1
  integer, parameter :: kkkk = 1
  
#if LATTICE==319
  
  ! swap pop7 and pop8

  buffservice3d(nx+kk,ny+kkk,1-kkkk)=f07sub(nx+kk,ny+kkk,1-kkkk)
  f07sub(nx+kk,ny+kkk,1-kkkk)=f08sub(nx+kk,ny+kkk,1-kkkk)
  f08sub(nx+kk,ny+kkk,1-kkkk)=buffservice3d(nx+kk,ny+kkk,1-kkkk)

  
  ! swap pop9 and pop10

  buffservice3d(nx+kk,ny+kkk,1-kkkk)=f09sub(nx+kk,ny+kkk,1-kkkk)
  f09sub(nx+kk,ny+kkk,1-kkkk)=f10sub(nx+kk,ny+kkk,1-kkkk)
  f10sub(nx+kk,ny+kkk,1-kkkk)=buffservice3d(nx+kk,ny+kkk,1-kkkk)

  
  ! swap pop11 and pop12

  buffservice3d(nx+kk,ny+kkk,1-kkkk)=f11sub(nx+kk,ny+kkk,1-kkkk)
  f11sub(nx+kk,ny+kkk,1-kkkk)=f12sub(nx+kk,ny+kkk,1-kkkk)
  f12sub(nx+kk,ny+kkk,1-kkkk)=buffservice3d(nx+kk,ny+kkk,1-kkkk)

  
  ! swap pop13 and pop14

  buffservice3d(nx+kk,ny+kkk,1-kkkk)=f13sub(nx+kk,ny+kkk,1-kkkk)
  f13sub(nx+kk,ny+kkk,1-kkkk)=f14sub(nx+kk,ny+kkk,1-kkkk)
  f14sub(nx+kk,ny+kkk,1-kkkk)=buffservice3d(nx+kk,ny+kkk,1-kkkk)

  
  ! swap pop15 and pop16

  buffservice3d(nx+kk,ny+kkk,1-kkkk)=f15sub(nx+kk,ny+kkk,1-kkkk)
  f15sub(nx+kk,ny+kkk,1-kkkk)=f16sub(nx+kk,ny+kkk,1-kkkk)
  f16sub(nx+kk,ny+kkk,1-kkkk)=buffservice3d(nx+kk,ny+kkk,1-kkkk)

  
  ! swap pop17 and pop18

  buffservice3d(nx+kk,ny+kkk,1-kkkk)=f17sub(nx+kk,ny+kkk,1-kkkk)
  f17sub(nx+kk,ny+kkk,1-kkkk)=f18sub(nx+kk,ny+kkk,1-kkkk)
  f18sub(nx+kk,ny+kkk,1-kkkk)=buffservice3d(nx+kk,ny+kkk,1-kkkk)
  
#else
  
  call warning(9,ZERO,latt_name)
  call error(12)
  
#endif
  return
  
 end subroutine apply_bounceback_corner_south_east_rear
 
!******************END PART TO MANAGE THE BOUNCEBACK********************

!*****************START PART TO MANAGE THE REFLECTION*******************
 
 subroutine driver_reflect_densities
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the reflection of the density 
!     variables if requested from the boundary conditions
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  select case(ibctype)
  case(0) ! 0 0 0
    call apply_reflection_densities_all
  case(1) ! 1 0 0 
    call apply_reflection_densities_along_yz
  case(2) ! 0 1 0 
    call apply_reflection_densities_along_xz
  case(3) ! 1 1 0 
    call apply_reflection_densities_along_z
  case(4) ! 0 0 1 
    call apply_reflection_densities_along_xy
  case(5) ! 1 0 1 
    call apply_reflection_densities_along_y
  case(6) ! 0 1 1 
    call apply_reflection_densities_along_x
  case(7) ! 1 1 1
    return
  case default
    call error(12)
  end select
  
  return
  
 end subroutine driver_reflect_densities
 
 subroutine apply_reflection_densities_all
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the complete reflection 
!     to density variables
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  !sides 6
  
  !apply reflection at north 1   !z
  call apply_reflection_north(1,nx,1,ny,rhoR)
  
  !apply reflection at south 2   !z
  call apply_reflection_south(rhoR)
  
  !apply reflection at east 3    !x
  call apply_reflection_east(1,ny,1,nz,rhoR)
  
  !apply reflection at west 4    !x
  call apply_reflection_west(1,ny,1,nz,rhoR)
  
  !apply reflection at front 5   !y
  call apply_reflection_front(rhoR)
  
  !apply reflection at rear 6   !y
  call apply_reflection_rear(rhoR)
  
  !edges 12
  
  !apply reflection at front east 1   !xy
  call apply_reflection_edge_front_east(rhoR)
  
  !apply reflection at front west 2   !xy
  call apply_reflection_edge_front_west(rhoR)
  
  !apply reflection at north east 3   !xz
  call apply_reflection_edge_north_east(rhoR)
  
  !apply reflection at north front 4   !yz
  call apply_reflection_edge_north_front(rhoR)
  
  !apply reflection at north rear 5   !yz
  call apply_reflection_edge_north_rear(rhoR)
  
  !apply reflection at north west 6   !xz
  call apply_reflection_edge_north_west(rhoR)
  
  !apply reflection at rear east 7   !xy
  call apply_reflection_edge_rear_east(rhoR)
  
  !apply reflection at rear west 8   !xy
  call apply_reflection_edge_rear_west(rhoR)
  
  !apply reflection at south east 9  !xz
  call apply_reflection_edge_south_east(rhoR)
  
  !apply reflection at south front 10  !yz
  call apply_reflection_edge_south_front(rhoR)
  
  !apply reflection at south rear 11  !yz
  call apply_reflection_edge_south_rear(rhoR)
  
  !apply reflection at south west 12  !xz
  call apply_reflection_edge_south_west(rhoR)
  
  !corner 8
  
  !apply reflection at north east front 1
  call apply_reflection_corner_north_east_front(rhoR)
  
  !apply reflection at north east rear 2
  call apply_reflection_corner_north_east_rear(rhoR)
  
  !apply reflection at north west rear 3
  call apply_reflection_corner_north_west_rear(rhoR)
  
  !apply reflection at north west front 4
  call apply_reflection_corner_north_west_front(rhoR)
  
  !apply reflection at south east front 5
  call apply_reflection_corner_south_east_front(rhoR)
  
  !apply reflection at south west front 6
  call apply_reflection_corner_south_west_front(rhoR)
  
  !apply reflection at south west rear 7
  call apply_reflection_corner_south_west_rear(rhoR)
  
  !apply reflection at south east rear 8
  call apply_reflection_corner_south_east_rear(rhoR)
  
  if(lsingle_fluid)return
  
  !sides 6
  
  !apply reflection at north 1   !z
  call apply_reflection_north(1,nx,1,ny,rhoB)
  
  !apply reflection at south 2   !z
  call apply_reflection_south(rhoB)
  
  !apply reflection at east 3    !x
  call apply_reflection_east(1,ny,1,nz,rhoB)
  
  !apply reflection at west 4    !x
  call apply_reflection_west(1,ny,1,nz,rhoB)
  
  !apply reflection at front 5   !y
  call apply_reflection_front(rhoB)
  
  !apply reflection at rear 6   !y
  call apply_reflection_rear(rhoB)
  
  !edges 12
  
  !apply reflection at front east 1   !xy
  call apply_reflection_edge_front_east(rhoB)
  
  !apply reflection at front west 2   !xy
  call apply_reflection_edge_front_west(rhoB)
  
  !apply reflection at north east 3   !xz
  call apply_reflection_edge_north_east(rhoB)
  
  !apply reflection at north front 4   !yz
  call apply_reflection_edge_north_front(rhoB)
  
  !apply reflection at north rear 5   !yz
  call apply_reflection_edge_north_rear(rhoB)
  
  !apply reflection at north west 6   !xz
  call apply_reflection_edge_north_west(rhoB)
  
  !apply reflection at rear east 7   !xy
  call apply_reflection_edge_rear_east(rhoB)
  
  !apply reflection at rear west 8   !xy
  call apply_reflection_edge_rear_west(rhoB)
  
  !apply reflection at south east 9  !xz
  call apply_reflection_edge_south_east(rhoB)
  
  !apply reflection at south front 10  !yz
  call apply_reflection_edge_south_front(rhoB)
  
  !apply reflection at south rear 11  !yz
  call apply_reflection_edge_south_rear(rhoB)
  
  !apply reflection at south west 12  !xz
  call apply_reflection_edge_south_west(rhoB)
  
  !corner 8
  
  !apply reflection at north east front 1
  call apply_reflection_corner_north_east_front(rhoB)
  
  !apply reflection at north east rear 2
  call apply_reflection_corner_north_east_rear(rhoB)
  
  !apply reflection at north west rear 3
  call apply_reflection_corner_north_west_rear(rhoB)
  
  !apply reflection at north west front 4
  call apply_reflection_corner_north_west_front(rhoB)
  
  !apply reflection at south east front 5
  call apply_reflection_corner_south_east_front(rhoB)
  
  !apply reflection at south west front 6
  call apply_reflection_corner_south_west_front(rhoB)
  
  !apply reflection at south west rear 7
  call apply_reflection_corner_south_west_rear(rhoB)
  
  !apply reflection at south east rear 8
  call apply_reflection_corner_south_east_rear(rhoB)
  
  return
  
 end subroutine apply_reflection_densities_all
 
 subroutine apply_reflection_densities_along_yz
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     to density variables along y and z
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  !sides 6
  
  !apply pbc at front 5 !y
  !red fluid
  call apply_reflection_front_frame_x(rhoR)
  
  !apply pbc at rear 6 !y
  !red fluid
  call apply_reflection_rear_frame_x(rhoR)
  
  !apply pbc at north 1 !z
  !red fluid
  !frame x
  call apply_reflection_north(1-nbuff,nx+nbuff,1,ny,rhoR)
  
  !apply pbc at south 2 !z
  !red fluid
  call apply_reflection_south_frame_x(rhoR)
  
  !edges 12
  
  !apply reflection at front east 1   !xy
  !call apply_reflection_edge_front_east_frame(rhoR)
  
  !apply reflection at front west 2   !xy
  !call apply_reflection_edge_front_west_frame(rhoR)
  
  !apply reflection at north east 3   !xz
  !call apply_reflection_edge_north_east_frame(rhoR)
  
  !apply reflection at north front 4   !yz
  call apply_reflection_edge_north_front_frame(rhoR)
  
  !apply reflection at north rear 5   !yz
  call apply_reflection_edge_north_rear_frame(rhoR)
  
  !apply reflection at north west 6   !xz
  !call apply_reflection_edge_north_west_frame(rhoR)
  
  !apply reflection at rear east 7   !xy
  !call apply_reflection_edge_rear_east_frame(rhoR)
  
  !apply reflection at rear west 8   !xy
  !call apply_reflection_edge_rear_west_frame(rhoR)
  
  !apply reflection at south east 9  !xz
  !call apply_reflection_edge_south_east_frame(rhoR)
  
  !apply reflection at south front 10  !yz
  call apply_reflection_edge_south_front_frame(rhoR)
  
  !apply reflection at south rear 11  !yz
  call apply_reflection_edge_south_rear_frame(rhoR)
  
  !apply reflection at south west 12  !xz
  !call apply_reflection_edge_south_west_frame(rhoR)

  
  if(lsingle_fluid)return
  
  !sides 6
  
  !apply pbc at front 5 !y
  !blue fluid
  call apply_reflection_front_frame_x(rhoB)
  
  !apply pbc at rear 6 !y
  !blue fluid
  call apply_reflection_rear_frame_x(rhoB)
  
  !apply pbc at north 1 !z
  !blue fluid
  !frame x
  call apply_reflection_north(1-nbuff,nx+nbuff,1,ny,rhoB)
  
  !apply pbc at south 2 !z
  !blue fluid
  call apply_reflection_south_frame_x(rhoB)
  
  !edges 12
  
  !apply reflection at front east 1   !xy
  !call apply_reflection_edge_front_east_frame(rhoB)
  
  !apply reflection at front west 2   !xy
  !call apply_reflection_edge_front_west_frame(rhoB)
  
  !apply reflection at north east 3   !xz
  !call apply_reflection_edge_north_east_frame(rhoB)
  
  !apply reflection at north front 4   !yz
  call apply_reflection_edge_north_front_frame(rhoB)
  
  !apply reflection at north rear 5   !yz
  call apply_reflection_edge_north_rear_frame(rhoB)
  
  !apply reflection at north west 6   !xz
  !call apply_reflection_edge_north_west_frame(rhoB)
  
  !apply reflection at rear east 7   !xy
  !call apply_reflection_edge_rear_east_frame(rhoB)
  
  !apply reflection at rear west 8   !xy
  !call apply_reflection_edge_rear_west_frame(rhoB)
  
  !apply reflection at south east 9  !xz
  !call apply_reflection_edge_south_east_frame(rhoB)
  
  !apply reflection at south front 10  !yz
  call apply_reflection_edge_south_front_frame(rhoB)
  
  !apply reflection at south rear 11  !yz
  call apply_reflection_edge_south_rear_frame(rhoB)
  
  !apply reflection at south west 12  !xz
  !call apply_reflection_edge_south_west_frame(rhoB)
  
  return
  
 end subroutine apply_reflection_densities_along_yz
 
 subroutine apply_reflection_densities_along_xz
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     to density variables along x and z
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  !sides 6
  !frame y
  !apply pbc at east 3 !x
  !red fluid
  call apply_reflection_east(1-nbuff,ny+nbuff,1,nz,rhoR)
  
  !apply pbc at west 4 !x
  !red fluid
  call apply_reflection_west(1-nbuff,ny+nbuff,1,nz,rhoR)
  
  !apply pbc at north 1 !z
  !red fluid
  !frame y
  call apply_reflection_north(1,nx,1-nbuff,ny+nbuff,rhoR)
  
  !apply pbc at south 2 !z
  !red fluid
  call apply_reflection_south_frame_y(rhoR)
  
  !edges 12
  
  !apply reflection at front east 1   !xy
  !call apply_reflection_edge_front_east_frame(rhoR)
  
  !apply reflection at front west 2   !xy
  !call apply_reflection_edge_front_west_frame(rhoR)
  
  !apply reflection at north east 3   !xz
  call apply_reflection_edge_north_east_frame(rhoR)
  
  !apply reflection at north front 4   !yz
  !call apply_reflection_edge_north_front_frame(rhoR)
  
  !apply reflection at north rear 5   !yz
  !call apply_reflection_edge_north_rear_frame(rhoR)
  
  !apply reflection at north west 6   !xz
  call apply_reflection_edge_north_west_frame(rhoR)
  
  !apply reflection at rear east 7   !xy
  !call apply_reflection_edge_rear_east_frame(rhoR)
  
  !apply reflection at rear west 8   !xy
  !call apply_reflection_edge_rear_west_frame(rhoR)
  
  !apply reflection at south east 9  !xz
  call apply_reflection_edge_south_east_frame(rhoR)
  
  !apply reflection at south front 10  !yz
  !call apply_reflection_edge_south_front_frame(rhoR)
  
  !apply reflection at south rear 11  !yz
  !call apply_reflection_edge_south_rear_frame(rhoR)
  
  !apply reflection at south west 12  !xz
  call apply_reflection_edge_south_west_frame(rhoR)
  
  if(lsingle_fluid)return
  
  !sides 6
  
  !apply pbc at east 3 !x
  !blue fluid
  call apply_reflection_east(1-nbuff,ny+nbuff,1,nz,rhoB)
  
  !apply pbc at west 4 !x
  !blue fluid
  call apply_reflection_west(1-nbuff,ny+nbuff,1,nz,rhoB)
  
  !apply pbc at north 1 !z
  !blue fluid
  !frame y
  call apply_reflection_north(1,nx,1-nbuff,ny+nbuff,rhoB)
  
  !apply pbc at south 2 !z
  !blue fluid
  call apply_reflection_south_frame_y(rhoB)
  
  !edges 12
  
  !apply reflection at front east 1   !xy
  !call apply_reflection_edge_front_east_frame(rhoB)
  
  !apply reflection at front west 2   !xy
  !call apply_reflection_edge_front_west_frame(rhoB)
  
  !apply reflection at north east 3   !xz
  call apply_reflection_edge_north_east_frame(rhoB)
  
  !apply reflection at north front 4   !yz
  !call apply_reflection_edge_north_front_frame(rhoB)
  
  !apply reflection at north rear 5   !yz
  !call apply_reflection_edge_north_rear_frame(rhoB)
  
  !apply reflection at north west 6   !xz
  call apply_reflection_edge_north_west_frame(rhoB)
  
  !apply reflection at rear east 7   !xy
  !call apply_reflection_edge_rear_east_frame(rhoB)
  
  !apply reflection at rear west 8   !xy
  !call apply_reflection_edge_rear_west_frame(rhoB)
  
  !apply reflection at south east 9  !xz
  call apply_reflection_edge_south_east_frame(rhoB)
  
  !apply reflection at south front 10  !yz
  !call apply_reflection_edge_south_front_frame(rhoB)
  
  !apply reflection at south rear 11  !yz
  !call apply_reflection_edge_south_rear_frame(rhoB)
  
  !apply reflection at south west 12  !xz
  call apply_reflection_edge_south_west_frame(rhoB)
  
  return
  
 end subroutine apply_reflection_densities_along_xz
 
 subroutine apply_reflection_densities_along_xy
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     to density variables along x and y
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  !sides 6
  !frame z
  !apply pbc at east 3 !x
  !red fluid
  call apply_reflection_east(1,ny,1-nbuff,nz+nbuff,rhoR)
  
  !apply pbc at west 4 !x
  !red fluid
  call apply_reflection_west(1,ny,1-nbuff,nz+nbuff,rhoR)
  
  !apply pbc at front 5 !y
  !red fluid
  call apply_reflection_front_frame_z(rhoR)
  
  !apply pbc at rear 6 !y
  !red fluid
  call apply_reflection_rear_frame_z(rhoR)
  
  !edges 12
  
  !apply reflection at front east 1   !xy
  call apply_reflection_edge_front_east_frame(rhoR)
  
  !apply reflection at front west 2   !xy
  call apply_reflection_edge_front_west_frame(rhoR)
  
  !apply reflection at north east 3   !xz
  !call apply_reflection_edge_north_east_frame(rhoR)
  
  !apply reflection at north front 4   !yz
  !call apply_reflection_edge_north_front_frame(rhoR)
  
  !apply reflection at north rear 5   !yz
  !call apply_reflection_edge_north_rear_frame(rhoR)
  
  !apply reflection at north west 6   !xz
  !call apply_reflection_edge_north_west_frame(rhoR)
  
  !apply reflection at rear east 7   !xy
  call apply_reflection_edge_rear_east_frame(rhoR)
  
  !apply reflection at rear west 8   !xy
  call apply_reflection_edge_rear_west_frame(rhoR)
  
  !apply reflection at south east 9  !xz
  !call apply_reflection_edge_south_east_frame(rhoR)
  
  !apply reflection at south front 10  !yz
  !call apply_reflection_edge_south_front_frame(rhoR)
  
  !apply reflection at south rear 11  !yz
  !call apply_reflection_edge_south_rear_frame(rhoR)
  
  !apply reflection at south west 12  !xz
  !call apply_reflection_edge_south_west_frame(rhoR)
  
  if(lsingle_fluid)return
  
  !sides 6
  
  !apply pbc at east 3 !x
  !blue fluid
  call apply_reflection_east(1,ny,1-nbuff,nz+nbuff,rhoB)
  
  !apply pbc at west 4 !x
  !blue fluid
  call apply_reflection_west(1,ny,1-nbuff,nz+nbuff,rhoB)
  
  !apply pbc at front 5 !y
  !blue fluid
  call apply_reflection_front_frame_z(rhoB)
  
  !apply pbc at rear 6 !y
  !blue fluid
  call apply_reflection_rear_frame_z(rhoB)
  
  !edges 12
  
  !apply reflection at front east 1   !xy
  call apply_reflection_edge_front_east_frame(rhoB)
  
  !apply reflection at front west 2   !xy
  call apply_reflection_edge_front_west_frame(rhoB)
  
  !apply reflection at north east 3   !xz
  !call apply_reflection_edge_north_east_frame(rhoB)
  
  !apply reflection at north front 4   !yz
  !call apply_reflection_edge_north_front_frame(rhoB)
  
  !apply reflection at north rear 5   !yz
  !call apply_reflection_edge_north_rear_frame(rhoB)
  
  !apply reflection at north west 6   !xz
  !call apply_reflection_edge_north_west_frame(rhoB)
  
  !apply reflection at rear east 7   !xy
  call apply_reflection_edge_rear_east_frame(rhoB)
  
  !apply reflection at rear west 8   !xy
  call apply_reflection_edge_rear_west_frame(rhoB)
  
  !apply reflection at south east 9  !xz
  !call apply_reflection_edge_south_east_frame(rhoB)
  
  !apply reflection at south front 10  !yz
  !call apply_reflection_edge_south_front_frame(rhoB)
  
  !apply reflection at south rear 11  !yz
  !call apply_reflection_edge_south_rear_frame(rhoB)
  
  !apply reflection at south west 12  !xz
  !call apply_reflection_edge_south_west_frame(rhoB)
  
  return
  
 end subroutine apply_reflection_densities_along_xy
 
 subroutine apply_reflection_densities_along_z
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     to density variables
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  !sides 6
  
  !apply pbc at north 1 !z
  !red fluid
  !frame
  call apply_reflection_north(1-nbuff,nx+nbuff,1-nbuff,ny+nbuff,rhoR)
  
  !apply pbc at south 2 !z
  !red fluid
  call apply_reflection_south_frame(rhoR)
  
  if(lsingle_fluid)return
  
  !sides 6
  
  !apply pbc at north 1 !z
  !blue fluid
  !frame
  call apply_reflection_north(1-nbuff,nx+nbuff,1-nbuff,ny+nbuff,rhoB)
  
  !apply pbc at south 2 !z
  !blue fluid
  call apply_reflection_south_frame(rhoB)
  
  return
  
 end subroutine apply_reflection_densities_along_z
 
 subroutine apply_reflection_densities_along_y
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     to density variables
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  !sides 6
  
  !apply pbc at front 5 !y
  !red fluid
  call apply_reflection_front_frame(rhoR)
  
  !apply pbc at rear 6 !y
  !red fluid
  call apply_reflection_rear_frame(rhoR)
  
  if(lsingle_fluid)return
  
  !sides 6
  
  !apply pbc at front 5 !y
  !blue fluid
  call apply_reflection_front_frame(rhoB)
  
  !apply pbc at rear 6 !y
  !blue fluid
  call apply_reflection_rear_frame(rhoB)
  
  return
  
 end subroutine apply_reflection_densities_along_y
 
 subroutine apply_reflection_densities_along_x
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     to density variables
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  !sides 6
  !frame complete
  !apply pbc at east 3 !x
  !red fluid
  call apply_reflection_east(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,rhoR)
  
  !apply pbc at west 4 !x
  !red fluid
  call apply_reflection_west(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,rhoR)
  
  if(lsingle_fluid)return
  
  !sides 6
  
  !apply pbc at east 3 !x
  !blue fluid
  call apply_reflection_east(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,rhoB)
  
  !apply pbc at west 4 !x
  !blue fluid
  call apply_reflection_west(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,rhoB)
  
  return
  
 end subroutine apply_reflection_densities_along_x
 
 subroutine driver_reflect_pops
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the reflection of the populations
!     if requested from the boundary conditions
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  select case(ibctype)
  case(0) ! 0 0 0 
    call apply_reflection_pops_all
  case(1) ! 1 0 0 
    call apply_reflection_pops_along_yz
  case(2) ! 0 1 0 
    call apply_reflection_pops_along_xz
  case(3) ! 1 1 0 
    call apply_reflection_pops_along_z
  case(4) ! 0 0 1 
    call apply_reflection_pops_along_xy
  case(5) ! 1 0 1 
    call apply_reflection_pops_along_y
  case(6) ! 0 1 1 
    call apply_reflection_pops_along_x
  case(7) ! 1 1 1
    return
  case default
    call error(12)
  end select
  
  return
  
 end subroutine driver_reflect_pops
 
 subroutine apply_reflection_pops_all
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the complete reflection 
!     to fluid populations
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  !sides 6
  
  !apply reflection at north 1   !z
  call apply_reflection_north(1,nx,1,ny,f00R)
  call apply_reflection_north(1,nx,1,ny,f01R)
  call apply_reflection_north(1,nx,1,ny,f02R)
  call apply_reflection_north(1,nx,1,ny,f03R)
  call apply_reflection_north(1,nx,1,ny,f04R)
  call apply_reflection_north(1,nx,1,ny,f05R)
  call apply_reflection_north(1,nx,1,ny,f06R)
  call apply_reflection_north(1,nx,1,ny,f07R)
  call apply_reflection_north(1,nx,1,ny,f08R)
  call apply_reflection_north(1,nx,1,ny,f09R)
  call apply_reflection_north(1,nx,1,ny,f10R)
  call apply_reflection_north(1,nx,1,ny,f11R)
  call apply_reflection_north(1,nx,1,ny,f12R)
  call apply_reflection_north(1,nx,1,ny,f13R)
  call apply_reflection_north(1,nx,1,ny,f14R)
  call apply_reflection_north(1,nx,1,ny,f15R)
  call apply_reflection_north(1,nx,1,ny,f16R)
  call apply_reflection_north(1,nx,1,ny,f17R)
  call apply_reflection_north(1,nx,1,ny,f18R)
  
  !apply reflection at south 2   !z
  call apply_reflection_south(f00R)
  call apply_reflection_south(f01R)
  call apply_reflection_south(f02R)
  call apply_reflection_south(f03R)
  call apply_reflection_south(f04R)
  call apply_reflection_south(f05R)
  call apply_reflection_south(f06R)
  call apply_reflection_south(f07R)
  call apply_reflection_south(f08R)
  call apply_reflection_south(f09R)
  call apply_reflection_south(f10R)
  call apply_reflection_south(f11R)
  call apply_reflection_south(f12R)
  call apply_reflection_south(f13R)
  call apply_reflection_south(f14R)
  call apply_reflection_south(f15R)
  call apply_reflection_south(f16R)
  call apply_reflection_south(f17R)
  call apply_reflection_south(f18R)
  
  !apply reflection at east 3    !x
  call apply_reflection_east(1,ny,1,nz,f00R)
  call apply_reflection_east(1,ny,1,nz,f01R)
  call apply_reflection_east(1,ny,1,nz,f02R)
  call apply_reflection_east(1,ny,1,nz,f03R)
  call apply_reflection_east(1,ny,1,nz,f04R)
  call apply_reflection_east(1,ny,1,nz,f05R)
  call apply_reflection_east(1,ny,1,nz,f06R)
  call apply_reflection_east(1,ny,1,nz,f07R)
  call apply_reflection_east(1,ny,1,nz,f08R)
  call apply_reflection_east(1,ny,1,nz,f09R)
  call apply_reflection_east(1,ny,1,nz,f10R)
  call apply_reflection_east(1,ny,1,nz,f11R)
  call apply_reflection_east(1,ny,1,nz,f12R)
  call apply_reflection_east(1,ny,1,nz,f13R)
  call apply_reflection_east(1,ny,1,nz,f14R)
  call apply_reflection_east(1,ny,1,nz,f15R)
  call apply_reflection_east(1,ny,1,nz,f16R)
  call apply_reflection_east(1,ny,1,nz,f17R)
  call apply_reflection_east(1,ny,1,nz,f18R)
  
  !apply reflection at west 4    !x
  call apply_reflection_west(1,ny,1,nz,f00R)
  call apply_reflection_west(1,ny,1,nz,f01R)
  call apply_reflection_west(1,ny,1,nz,f02R)
  call apply_reflection_west(1,ny,1,nz,f03R)
  call apply_reflection_west(1,ny,1,nz,f04R)
  call apply_reflection_west(1,ny,1,nz,f05R)
  call apply_reflection_west(1,ny,1,nz,f06R)
  call apply_reflection_west(1,ny,1,nz,f07R)
  call apply_reflection_west(1,ny,1,nz,f08R)
  call apply_reflection_west(1,ny,1,nz,f09R)
  call apply_reflection_west(1,ny,1,nz,f10R)
  call apply_reflection_west(1,ny,1,nz,f11R)
  call apply_reflection_west(1,ny,1,nz,f12R)
  call apply_reflection_west(1,ny,1,nz,f13R)
  call apply_reflection_west(1,ny,1,nz,f14R)
  call apply_reflection_west(1,ny,1,nz,f15R)
  call apply_reflection_west(1,ny,1,nz,f16R)
  call apply_reflection_west(1,ny,1,nz,f17R)
  call apply_reflection_west(1,ny,1,nz,f18R)
  
  !apply reflection at front 5   !y
  call apply_reflection_front(f00R)
  call apply_reflection_front(f01R)
  call apply_reflection_front(f02R)
  call apply_reflection_front(f03R)
  call apply_reflection_front(f04R)
  call apply_reflection_front(f05R)
  call apply_reflection_front(f06R)
  call apply_reflection_front(f07R)
  call apply_reflection_front(f08R)
  call apply_reflection_front(f09R)
  call apply_reflection_front(f10R)
  call apply_reflection_front(f11R)
  call apply_reflection_front(f12R)
  call apply_reflection_front(f13R)
  call apply_reflection_front(f14R)
  call apply_reflection_front(f15R)
  call apply_reflection_front(f16R)
  call apply_reflection_front(f17R)
  call apply_reflection_front(f18R)
  
  !apply reflection at rear 6   !y
  call apply_reflection_rear(f00R)
  call apply_reflection_rear(f01R)
  call apply_reflection_rear(f02R)
  call apply_reflection_rear(f03R)
  call apply_reflection_rear(f04R)
  call apply_reflection_rear(f05R)
  call apply_reflection_rear(f06R)
  call apply_reflection_rear(f07R)
  call apply_reflection_rear(f08R)
  call apply_reflection_rear(f09R)
  call apply_reflection_rear(f10R)
  call apply_reflection_rear(f11R)
  call apply_reflection_rear(f12R)
  call apply_reflection_rear(f13R)
  call apply_reflection_rear(f14R)
  call apply_reflection_rear(f15R)
  call apply_reflection_rear(f16R)
  call apply_reflection_rear(f17R)
  call apply_reflection_rear(f18R)
  
  !edges 12
  
  !apply reflection at front east 1   !xy
  call apply_reflection_edge_front_east(f00R)
  call apply_reflection_edge_front_east(f01R)
  call apply_reflection_edge_front_east(f02R)
  call apply_reflection_edge_front_east(f03R)
  call apply_reflection_edge_front_east(f04R)
  call apply_reflection_edge_front_east(f05R)
  call apply_reflection_edge_front_east(f06R)
  call apply_reflection_edge_front_east(f07R)
  call apply_reflection_edge_front_east(f08R)
  call apply_reflection_edge_front_east(f09R)
  call apply_reflection_edge_front_east(f10R)
  call apply_reflection_edge_front_east(f11R)
  call apply_reflection_edge_front_east(f12R)
  call apply_reflection_edge_front_east(f13R)
  call apply_reflection_edge_front_east(f14R)
  call apply_reflection_edge_front_east(f15R)
  call apply_reflection_edge_front_east(f16R)
  call apply_reflection_edge_front_east(f17R)
  call apply_reflection_edge_front_east(f18R)
  
  !apply reflection at front west 2   !xy
  call apply_reflection_edge_front_west(f00R)
  call apply_reflection_edge_front_west(f01R)
  call apply_reflection_edge_front_west(f02R)
  call apply_reflection_edge_front_west(f03R)
  call apply_reflection_edge_front_west(f04R)
  call apply_reflection_edge_front_west(f05R)
  call apply_reflection_edge_front_west(f06R)
  call apply_reflection_edge_front_west(f07R)
  call apply_reflection_edge_front_west(f08R)
  call apply_reflection_edge_front_west(f09R)
  call apply_reflection_edge_front_west(f10R)
  call apply_reflection_edge_front_west(f11R)
  call apply_reflection_edge_front_west(f12R)
  call apply_reflection_edge_front_west(f13R)
  call apply_reflection_edge_front_west(f14R)
  call apply_reflection_edge_front_west(f15R)
  call apply_reflection_edge_front_west(f16R)
  call apply_reflection_edge_front_west(f17R)
  call apply_reflection_edge_front_west(f18R)
  
  !apply reflection at north east 3   !xz
  call apply_reflection_edge_north_east(f00R)
  call apply_reflection_edge_north_east(f01R)
  call apply_reflection_edge_north_east(f02R)
  call apply_reflection_edge_north_east(f03R)
  call apply_reflection_edge_north_east(f04R)
  call apply_reflection_edge_north_east(f05R)
  call apply_reflection_edge_north_east(f06R)
  call apply_reflection_edge_north_east(f07R)
  call apply_reflection_edge_north_east(f08R)
  call apply_reflection_edge_north_east(f09R)
  call apply_reflection_edge_north_east(f10R)
  call apply_reflection_edge_north_east(f11R)
  call apply_reflection_edge_north_east(f12R)
  call apply_reflection_edge_north_east(f13R)
  call apply_reflection_edge_north_east(f14R)
  call apply_reflection_edge_north_east(f15R)
  call apply_reflection_edge_north_east(f16R)
  call apply_reflection_edge_north_east(f17R)
  call apply_reflection_edge_north_east(f18R)
  
  !apply reflection at north front 4   !yz
  call apply_reflection_edge_north_front(f00R)
  call apply_reflection_edge_north_front(f01R)
  call apply_reflection_edge_north_front(f02R)
  call apply_reflection_edge_north_front(f03R)
  call apply_reflection_edge_north_front(f04R)
  call apply_reflection_edge_north_front(f05R)
  call apply_reflection_edge_north_front(f06R)
  call apply_reflection_edge_north_front(f07R)
  call apply_reflection_edge_north_front(f08R)
  call apply_reflection_edge_north_front(f09R)
  call apply_reflection_edge_north_front(f10R)
  call apply_reflection_edge_north_front(f11R)
  call apply_reflection_edge_north_front(f12R)
  call apply_reflection_edge_north_front(f13R)
  call apply_reflection_edge_north_front(f14R)
  call apply_reflection_edge_north_front(f15R)
  call apply_reflection_edge_north_front(f16R)
  call apply_reflection_edge_north_front(f17R)
  call apply_reflection_edge_north_front(f18R)
  
  !apply reflection at north rear 5   !yz
  call apply_reflection_edge_north_rear(f00R)
  call apply_reflection_edge_north_rear(f01R)
  call apply_reflection_edge_north_rear(f02R)
  call apply_reflection_edge_north_rear(f03R)
  call apply_reflection_edge_north_rear(f04R)
  call apply_reflection_edge_north_rear(f05R)
  call apply_reflection_edge_north_rear(f06R)
  call apply_reflection_edge_north_rear(f07R)
  call apply_reflection_edge_north_rear(f08R)
  call apply_reflection_edge_north_rear(f09R)
  call apply_reflection_edge_north_rear(f10R)
  call apply_reflection_edge_north_rear(f11R)
  call apply_reflection_edge_north_rear(f12R)
  call apply_reflection_edge_north_rear(f13R)
  call apply_reflection_edge_north_rear(f14R)
  call apply_reflection_edge_north_rear(f15R)
  call apply_reflection_edge_north_rear(f16R)
  call apply_reflection_edge_north_rear(f17R)
  call apply_reflection_edge_north_rear(f18R)
  
  !apply reflection at north west 6   !xz
  call apply_reflection_edge_north_west(f00R)
  call apply_reflection_edge_north_west(f01R)
  call apply_reflection_edge_north_west(f02R)
  call apply_reflection_edge_north_west(f03R)
  call apply_reflection_edge_north_west(f04R)
  call apply_reflection_edge_north_west(f05R)
  call apply_reflection_edge_north_west(f06R)
  call apply_reflection_edge_north_west(f07R)
  call apply_reflection_edge_north_west(f08R)
  call apply_reflection_edge_north_west(f09R)
  call apply_reflection_edge_north_west(f10R)
  call apply_reflection_edge_north_west(f11R)
  call apply_reflection_edge_north_west(f12R)
  call apply_reflection_edge_north_west(f13R)
  call apply_reflection_edge_north_west(f14R)
  call apply_reflection_edge_north_west(f15R)
  call apply_reflection_edge_north_west(f16R)
  call apply_reflection_edge_north_west(f17R)
  call apply_reflection_edge_north_west(f18R)
  
  !apply reflection at rear east 7   !xy
  call apply_reflection_edge_rear_east(f00R)
  call apply_reflection_edge_rear_east(f01R)
  call apply_reflection_edge_rear_east(f02R)
  call apply_reflection_edge_rear_east(f03R)
  call apply_reflection_edge_rear_east(f04R)
  call apply_reflection_edge_rear_east(f05R)
  call apply_reflection_edge_rear_east(f06R)
  call apply_reflection_edge_rear_east(f07R)
  call apply_reflection_edge_rear_east(f08R)
  call apply_reflection_edge_rear_east(f09R)
  call apply_reflection_edge_rear_east(f10R)
  call apply_reflection_edge_rear_east(f11R)
  call apply_reflection_edge_rear_east(f12R)
  call apply_reflection_edge_rear_east(f13R)
  call apply_reflection_edge_rear_east(f14R)
  call apply_reflection_edge_rear_east(f15R)
  call apply_reflection_edge_rear_east(f16R)
  call apply_reflection_edge_rear_east(f17R)
  call apply_reflection_edge_rear_east(f18R)
  
  !apply reflection at rear west 8   !xy
  call apply_reflection_edge_rear_west(f00R)
  call apply_reflection_edge_rear_west(f01R)
  call apply_reflection_edge_rear_west(f02R)
  call apply_reflection_edge_rear_west(f03R)
  call apply_reflection_edge_rear_west(f04R)
  call apply_reflection_edge_rear_west(f05R)
  call apply_reflection_edge_rear_west(f06R)
  call apply_reflection_edge_rear_west(f07R)
  call apply_reflection_edge_rear_west(f08R)
  call apply_reflection_edge_rear_west(f09R)
  call apply_reflection_edge_rear_west(f10R)
  call apply_reflection_edge_rear_west(f11R)
  call apply_reflection_edge_rear_west(f12R)
  call apply_reflection_edge_rear_west(f13R)
  call apply_reflection_edge_rear_west(f14R)
  call apply_reflection_edge_rear_west(f15R)
  call apply_reflection_edge_rear_west(f16R)
  call apply_reflection_edge_rear_west(f17R)
  call apply_reflection_edge_rear_west(f18R)
  
  !apply reflection at south east 9  !xz
  call apply_reflection_edge_south_east(f00R)
  call apply_reflection_edge_south_east(f01R)
  call apply_reflection_edge_south_east(f02R)
  call apply_reflection_edge_south_east(f03R)
  call apply_reflection_edge_south_east(f04R)
  call apply_reflection_edge_south_east(f05R)
  call apply_reflection_edge_south_east(f06R)
  call apply_reflection_edge_south_east(f07R)
  call apply_reflection_edge_south_east(f08R)
  call apply_reflection_edge_south_east(f09R)
  call apply_reflection_edge_south_east(f10R)
  call apply_reflection_edge_south_east(f11R)
  call apply_reflection_edge_south_east(f12R)
  call apply_reflection_edge_south_east(f13R)
  call apply_reflection_edge_south_east(f14R)
  call apply_reflection_edge_south_east(f15R)
  call apply_reflection_edge_south_east(f16R)
  call apply_reflection_edge_south_east(f17R)
  call apply_reflection_edge_south_east(f18R)
  
  !apply reflection at south front 10  !yz
  call apply_reflection_edge_south_front(f00R)
  call apply_reflection_edge_south_front(f01R)
  call apply_reflection_edge_south_front(f02R)
  call apply_reflection_edge_south_front(f03R)
  call apply_reflection_edge_south_front(f04R)
  call apply_reflection_edge_south_front(f05R)
  call apply_reflection_edge_south_front(f06R)
  call apply_reflection_edge_south_front(f07R)
  call apply_reflection_edge_south_front(f08R)
  call apply_reflection_edge_south_front(f09R)
  call apply_reflection_edge_south_front(f10R)
  call apply_reflection_edge_south_front(f11R)
  call apply_reflection_edge_south_front(f12R)
  call apply_reflection_edge_south_front(f13R)
  call apply_reflection_edge_south_front(f14R)
  call apply_reflection_edge_south_front(f15R)
  call apply_reflection_edge_south_front(f16R)
  call apply_reflection_edge_south_front(f17R)
  call apply_reflection_edge_south_front(f18R)
  
  !apply reflection at south rear 11  !yz
  call apply_reflection_edge_south_rear(f00R)
  call apply_reflection_edge_south_rear(f01R)
  call apply_reflection_edge_south_rear(f02R)
  call apply_reflection_edge_south_rear(f03R)
  call apply_reflection_edge_south_rear(f04R)
  call apply_reflection_edge_south_rear(f05R)
  call apply_reflection_edge_south_rear(f06R)
  call apply_reflection_edge_south_rear(f07R)
  call apply_reflection_edge_south_rear(f08R)
  call apply_reflection_edge_south_rear(f09R)
  call apply_reflection_edge_south_rear(f10R)
  call apply_reflection_edge_south_rear(f11R)
  call apply_reflection_edge_south_rear(f12R)
  call apply_reflection_edge_south_rear(f13R)
  call apply_reflection_edge_south_rear(f14R)
  call apply_reflection_edge_south_rear(f15R)
  call apply_reflection_edge_south_rear(f16R)
  call apply_reflection_edge_south_rear(f17R)
  call apply_reflection_edge_south_rear(f18R)
  
  !apply reflection at south west 12  !xz
  call apply_reflection_edge_south_west(f00R)
  call apply_reflection_edge_south_west(f01R)
  call apply_reflection_edge_south_west(f02R)
  call apply_reflection_edge_south_west(f03R)
  call apply_reflection_edge_south_west(f04R)
  call apply_reflection_edge_south_west(f05R)
  call apply_reflection_edge_south_west(f06R)
  call apply_reflection_edge_south_west(f07R)
  call apply_reflection_edge_south_west(f08R)
  call apply_reflection_edge_south_west(f09R)
  call apply_reflection_edge_south_west(f10R)
  call apply_reflection_edge_south_west(f11R)
  call apply_reflection_edge_south_west(f12R)
  call apply_reflection_edge_south_west(f13R)
  call apply_reflection_edge_south_west(f14R)
  call apply_reflection_edge_south_west(f15R)
  call apply_reflection_edge_south_west(f16R)
  call apply_reflection_edge_south_west(f17R)
  call apply_reflection_edge_south_west(f18R)
  
  !corner 8
  
  !apply reflection at north east front 1
  call apply_reflection_corner_north_east_front(f00R)
  call apply_reflection_corner_north_east_front(f01R)
  call apply_reflection_corner_north_east_front(f02R)
  call apply_reflection_corner_north_east_front(f03R)
  call apply_reflection_corner_north_east_front(f04R)
  call apply_reflection_corner_north_east_front(f05R)
  call apply_reflection_corner_north_east_front(f06R)
  call apply_reflection_corner_north_east_front(f07R)
  call apply_reflection_corner_north_east_front(f08R)
  call apply_reflection_corner_north_east_front(f09R)
  call apply_reflection_corner_north_east_front(f10R)
  call apply_reflection_corner_north_east_front(f11R)
  call apply_reflection_corner_north_east_front(f12R)
  call apply_reflection_corner_north_east_front(f13R)
  call apply_reflection_corner_north_east_front(f14R)
  call apply_reflection_corner_north_east_front(f15R)
  call apply_reflection_corner_north_east_front(f16R)
  call apply_reflection_corner_north_east_front(f17R)
  call apply_reflection_corner_north_east_front(f18R)
  
  !apply reflection at north east rear 2
  call apply_reflection_corner_north_east_rear(f00R)
  call apply_reflection_corner_north_east_rear(f01R)
  call apply_reflection_corner_north_east_rear(f02R)
  call apply_reflection_corner_north_east_rear(f03R)
  call apply_reflection_corner_north_east_rear(f04R)
  call apply_reflection_corner_north_east_rear(f05R)
  call apply_reflection_corner_north_east_rear(f06R)
  call apply_reflection_corner_north_east_rear(f07R)
  call apply_reflection_corner_north_east_rear(f08R)
  call apply_reflection_corner_north_east_rear(f09R)
  call apply_reflection_corner_north_east_rear(f10R)
  call apply_reflection_corner_north_east_rear(f11R)
  call apply_reflection_corner_north_east_rear(f12R)
  call apply_reflection_corner_north_east_rear(f13R)
  call apply_reflection_corner_north_east_rear(f14R)
  call apply_reflection_corner_north_east_rear(f15R)
  call apply_reflection_corner_north_east_rear(f16R)
  call apply_reflection_corner_north_east_rear(f17R)
  call apply_reflection_corner_north_east_rear(f18R)
  
  !apply reflection at north west rear 3
  call apply_reflection_corner_north_west_rear(f00R)
  call apply_reflection_corner_north_west_rear(f01R)
  call apply_reflection_corner_north_west_rear(f02R)
  call apply_reflection_corner_north_west_rear(f03R)
  call apply_reflection_corner_north_west_rear(f04R)
  call apply_reflection_corner_north_west_rear(f05R)
  call apply_reflection_corner_north_west_rear(f06R)
  call apply_reflection_corner_north_west_rear(f07R)
  call apply_reflection_corner_north_west_rear(f08R)
  call apply_reflection_corner_north_west_rear(f09R)
  call apply_reflection_corner_north_west_rear(f10R)
  call apply_reflection_corner_north_west_rear(f11R)
  call apply_reflection_corner_north_west_rear(f12R)
  call apply_reflection_corner_north_west_rear(f13R)
  call apply_reflection_corner_north_west_rear(f14R)
  call apply_reflection_corner_north_west_rear(f15R)
  call apply_reflection_corner_north_west_rear(f16R)
  call apply_reflection_corner_north_west_rear(f17R)
  call apply_reflection_corner_north_west_rear(f18R)
  
  !apply reflection at north west front 4
  call apply_reflection_corner_north_west_front(f00R)
  call apply_reflection_corner_north_west_front(f01R)
  call apply_reflection_corner_north_west_front(f02R)
  call apply_reflection_corner_north_west_front(f03R)
  call apply_reflection_corner_north_west_front(f04R)
  call apply_reflection_corner_north_west_front(f05R)
  call apply_reflection_corner_north_west_front(f06R)
  call apply_reflection_corner_north_west_front(f07R)
  call apply_reflection_corner_north_west_front(f08R)
  call apply_reflection_corner_north_west_front(f09R)
  call apply_reflection_corner_north_west_front(f10R)
  call apply_reflection_corner_north_west_front(f11R)
  call apply_reflection_corner_north_west_front(f12R)
  call apply_reflection_corner_north_west_front(f13R)
  call apply_reflection_corner_north_west_front(f14R)
  call apply_reflection_corner_north_west_front(f15R)
  call apply_reflection_corner_north_west_front(f16R)
  call apply_reflection_corner_north_west_front(f17R)
  call apply_reflection_corner_north_west_front(f18R)
  
  !apply reflection at south east front 5
  call apply_reflection_corner_south_east_front(f00R)
  call apply_reflection_corner_south_east_front(f01R)
  call apply_reflection_corner_south_east_front(f02R)
  call apply_reflection_corner_south_east_front(f03R)
  call apply_reflection_corner_south_east_front(f04R)
  call apply_reflection_corner_south_east_front(f05R)
  call apply_reflection_corner_south_east_front(f06R)
  call apply_reflection_corner_south_east_front(f07R)
  call apply_reflection_corner_south_east_front(f08R)
  call apply_reflection_corner_south_east_front(f09R)
  call apply_reflection_corner_south_east_front(f10R)
  call apply_reflection_corner_south_east_front(f11R)
  call apply_reflection_corner_south_east_front(f12R)
  call apply_reflection_corner_south_east_front(f13R)
  call apply_reflection_corner_south_east_front(f14R)
  call apply_reflection_corner_south_east_front(f15R)
  call apply_reflection_corner_south_east_front(f16R)
  call apply_reflection_corner_south_east_front(f17R)
  call apply_reflection_corner_south_east_front(f18R)
  
  !apply reflection at south west front 6
  call apply_reflection_corner_south_west_front(f00R)
  call apply_reflection_corner_south_west_front(f01R)
  call apply_reflection_corner_south_west_front(f02R)
  call apply_reflection_corner_south_west_front(f03R)
  call apply_reflection_corner_south_west_front(f04R)
  call apply_reflection_corner_south_west_front(f05R)
  call apply_reflection_corner_south_west_front(f06R)
  call apply_reflection_corner_south_west_front(f07R)
  call apply_reflection_corner_south_west_front(f08R)
  call apply_reflection_corner_south_west_front(f09R)
  call apply_reflection_corner_south_west_front(f10R)
  call apply_reflection_corner_south_west_front(f11R)
  call apply_reflection_corner_south_west_front(f12R)
  call apply_reflection_corner_south_west_front(f13R)
  call apply_reflection_corner_south_west_front(f14R)
  call apply_reflection_corner_south_west_front(f15R)
  call apply_reflection_corner_south_west_front(f16R)
  call apply_reflection_corner_south_west_front(f17R)
  call apply_reflection_corner_south_west_front(f18R)
  
  !apply reflection at south west rear 7
  call apply_reflection_corner_south_west_rear(f00R)
  call apply_reflection_corner_south_west_rear(f01R)
  call apply_reflection_corner_south_west_rear(f02R)
  call apply_reflection_corner_south_west_rear(f03R)
  call apply_reflection_corner_south_west_rear(f04R)
  call apply_reflection_corner_south_west_rear(f05R)
  call apply_reflection_corner_south_west_rear(f06R)
  call apply_reflection_corner_south_west_rear(f07R)
  call apply_reflection_corner_south_west_rear(f08R)
  call apply_reflection_corner_south_west_rear(f09R)
  call apply_reflection_corner_south_west_rear(f10R)
  call apply_reflection_corner_south_west_rear(f11R)
  call apply_reflection_corner_south_west_rear(f12R)
  call apply_reflection_corner_south_west_rear(f13R)
  call apply_reflection_corner_south_west_rear(f14R)
  call apply_reflection_corner_south_west_rear(f15R)
  call apply_reflection_corner_south_west_rear(f16R)
  call apply_reflection_corner_south_west_rear(f17R)
  call apply_reflection_corner_south_west_rear(f18R)
  
  !apply reflection at south east rear 8
  call apply_reflection_corner_south_east_rear(f00R)
  call apply_reflection_corner_south_east_rear(f01R)
  call apply_reflection_corner_south_east_rear(f02R)
  call apply_reflection_corner_south_east_rear(f03R)
  call apply_reflection_corner_south_east_rear(f04R)
  call apply_reflection_corner_south_east_rear(f05R)
  call apply_reflection_corner_south_east_rear(f06R)
  call apply_reflection_corner_south_east_rear(f07R)
  call apply_reflection_corner_south_east_rear(f08R)
  call apply_reflection_corner_south_east_rear(f09R)
  call apply_reflection_corner_south_east_rear(f10R)
  call apply_reflection_corner_south_east_rear(f11R)
  call apply_reflection_corner_south_east_rear(f12R)
  call apply_reflection_corner_south_east_rear(f13R)
  call apply_reflection_corner_south_east_rear(f14R)
  call apply_reflection_corner_south_east_rear(f15R)
  call apply_reflection_corner_south_east_rear(f16R)
  call apply_reflection_corner_south_east_rear(f17R)
  call apply_reflection_corner_south_east_rear(f18R)
  
  if(lsingle_fluid)return
  
  !sides 6
  
  !apply reflection at north 1   !z
  call apply_reflection_north(1,nx,1,ny,f00B)
  call apply_reflection_north(1,nx,1,ny,f01B)
  call apply_reflection_north(1,nx,1,ny,f02B)
  call apply_reflection_north(1,nx,1,ny,f03B)
  call apply_reflection_north(1,nx,1,ny,f04B)
  call apply_reflection_north(1,nx,1,ny,f05B)
  call apply_reflection_north(1,nx,1,ny,f06B)
  call apply_reflection_north(1,nx,1,ny,f07B)
  call apply_reflection_north(1,nx,1,ny,f08B)
  call apply_reflection_north(1,nx,1,ny,f09B)
  call apply_reflection_north(1,nx,1,ny,f10B)
  call apply_reflection_north(1,nx,1,ny,f11B)
  call apply_reflection_north(1,nx,1,ny,f12B)
  call apply_reflection_north(1,nx,1,ny,f13B)
  call apply_reflection_north(1,nx,1,ny,f14B)
  call apply_reflection_north(1,nx,1,ny,f15B)
  call apply_reflection_north(1,nx,1,ny,f16B)
  call apply_reflection_north(1,nx,1,ny,f17B)
  call apply_reflection_north(1,nx,1,ny,f18B)
  
  !apply reflection at south 2   !z
  call apply_reflection_south(f00B)
  call apply_reflection_south(f01B)
  call apply_reflection_south(f02B)
  call apply_reflection_south(f03B)
  call apply_reflection_south(f04B)
  call apply_reflection_south(f05B)
  call apply_reflection_south(f06B)
  call apply_reflection_south(f07B)
  call apply_reflection_south(f08B)
  call apply_reflection_south(f09B)
  call apply_reflection_south(f10B)
  call apply_reflection_south(f11B)
  call apply_reflection_south(f12B)
  call apply_reflection_south(f13B)
  call apply_reflection_south(f14B)
  call apply_reflection_south(f15B)
  call apply_reflection_south(f16B)
  call apply_reflection_south(f17B)
  call apply_reflection_south(f18B)
  
  !apply reflection at east 3    !x
  call apply_reflection_east(1,ny,1,nz,f00B)
  call apply_reflection_east(1,ny,1,nz,f01B)
  call apply_reflection_east(1,ny,1,nz,f02B)
  call apply_reflection_east(1,ny,1,nz,f03B)
  call apply_reflection_east(1,ny,1,nz,f04B)
  call apply_reflection_east(1,ny,1,nz,f05B)
  call apply_reflection_east(1,ny,1,nz,f06B)
  call apply_reflection_east(1,ny,1,nz,f07B)
  call apply_reflection_east(1,ny,1,nz,f08B)
  call apply_reflection_east(1,ny,1,nz,f09B)
  call apply_reflection_east(1,ny,1,nz,f10B)
  call apply_reflection_east(1,ny,1,nz,f11B)
  call apply_reflection_east(1,ny,1,nz,f12B)
  call apply_reflection_east(1,ny,1,nz,f13B)
  call apply_reflection_east(1,ny,1,nz,f14B)
  call apply_reflection_east(1,ny,1,nz,f15B)
  call apply_reflection_east(1,ny,1,nz,f16B)
  call apply_reflection_east(1,ny,1,nz,f17B)
  call apply_reflection_east(1,ny,1,nz,f18B)
  
  !apply reflection at west 4    !x
  call apply_reflection_west(1,ny,1,nz,f00B)
  call apply_reflection_west(1,ny,1,nz,f01B)
  call apply_reflection_west(1,ny,1,nz,f02B)
  call apply_reflection_west(1,ny,1,nz,f03B)
  call apply_reflection_west(1,ny,1,nz,f04B)
  call apply_reflection_west(1,ny,1,nz,f05B)
  call apply_reflection_west(1,ny,1,nz,f06B)
  call apply_reflection_west(1,ny,1,nz,f07B)
  call apply_reflection_west(1,ny,1,nz,f08B)
  call apply_reflection_west(1,ny,1,nz,f09B)
  call apply_reflection_west(1,ny,1,nz,f10B)
  call apply_reflection_west(1,ny,1,nz,f11B)
  call apply_reflection_west(1,ny,1,nz,f12B)
  call apply_reflection_west(1,ny,1,nz,f13B)
  call apply_reflection_west(1,ny,1,nz,f14B)
  call apply_reflection_west(1,ny,1,nz,f15B)
  call apply_reflection_west(1,ny,1,nz,f16B)
  call apply_reflection_west(1,ny,1,nz,f17B)
  call apply_reflection_west(1,ny,1,nz,f18B)
  
  !apply reflection at front 5   !y
  call apply_reflection_front(f00B)
  call apply_reflection_front(f01B)
  call apply_reflection_front(f02B)
  call apply_reflection_front(f03B)
  call apply_reflection_front(f04B)
  call apply_reflection_front(f05B)
  call apply_reflection_front(f06B)
  call apply_reflection_front(f07B)
  call apply_reflection_front(f08B)
  call apply_reflection_front(f09B)
  call apply_reflection_front(f10B)
  call apply_reflection_front(f11B)
  call apply_reflection_front(f12B)
  call apply_reflection_front(f13B)
  call apply_reflection_front(f14B)
  call apply_reflection_front(f15B)
  call apply_reflection_front(f16B)
  call apply_reflection_front(f17B)
  call apply_reflection_front(f18B)
  
  !apply reflection at rear 6   !y
  call apply_reflection_rear(f00B)
  call apply_reflection_rear(f01B)
  call apply_reflection_rear(f02B)
  call apply_reflection_rear(f03B)
  call apply_reflection_rear(f04B)
  call apply_reflection_rear(f05B)
  call apply_reflection_rear(f06B)
  call apply_reflection_rear(f07B)
  call apply_reflection_rear(f08B)
  call apply_reflection_rear(f09B)
  call apply_reflection_rear(f10B)
  call apply_reflection_rear(f11B)
  call apply_reflection_rear(f12B)
  call apply_reflection_rear(f13B)
  call apply_reflection_rear(f14B)
  call apply_reflection_rear(f15B)
  call apply_reflection_rear(f16B)
  call apply_reflection_rear(f17B)
  call apply_reflection_rear(f18B)
  
  !edges 12
  
  !apply reflection at front east 1   !xy
  call apply_reflection_edge_front_east(f00B)
  call apply_reflection_edge_front_east(f01B)
  call apply_reflection_edge_front_east(f02B)
  call apply_reflection_edge_front_east(f03B)
  call apply_reflection_edge_front_east(f04B)
  call apply_reflection_edge_front_east(f05B)
  call apply_reflection_edge_front_east(f06B)
  call apply_reflection_edge_front_east(f07B)
  call apply_reflection_edge_front_east(f08B)
  call apply_reflection_edge_front_east(f09B)
  call apply_reflection_edge_front_east(f10B)
  call apply_reflection_edge_front_east(f11B)
  call apply_reflection_edge_front_east(f12B)
  call apply_reflection_edge_front_east(f13B)
  call apply_reflection_edge_front_east(f14B)
  call apply_reflection_edge_front_east(f15B)
  call apply_reflection_edge_front_east(f16B)
  call apply_reflection_edge_front_east(f17B)
  call apply_reflection_edge_front_east(f18B)
  
  !apply reflection at front west 2   !xy
  call apply_reflection_edge_front_west(f00B)
  call apply_reflection_edge_front_west(f01B)
  call apply_reflection_edge_front_west(f02B)
  call apply_reflection_edge_front_west(f03B)
  call apply_reflection_edge_front_west(f04B)
  call apply_reflection_edge_front_west(f05B)
  call apply_reflection_edge_front_west(f06B)
  call apply_reflection_edge_front_west(f07B)
  call apply_reflection_edge_front_west(f08B)
  call apply_reflection_edge_front_west(f09B)
  call apply_reflection_edge_front_west(f10B)
  call apply_reflection_edge_front_west(f11B)
  call apply_reflection_edge_front_west(f12B)
  call apply_reflection_edge_front_west(f13B)
  call apply_reflection_edge_front_west(f14B)
  call apply_reflection_edge_front_west(f15B)
  call apply_reflection_edge_front_west(f16B)
  call apply_reflection_edge_front_west(f17B)
  call apply_reflection_edge_front_west(f18B)
  
  !apply reflection at north east 3   !xz
  call apply_reflection_edge_north_east(f00B)
  call apply_reflection_edge_north_east(f01B)
  call apply_reflection_edge_north_east(f02B)
  call apply_reflection_edge_north_east(f03B)
  call apply_reflection_edge_north_east(f04B)
  call apply_reflection_edge_north_east(f05B)
  call apply_reflection_edge_north_east(f06B)
  call apply_reflection_edge_north_east(f07B)
  call apply_reflection_edge_north_east(f08B)
  call apply_reflection_edge_north_east(f09B)
  call apply_reflection_edge_north_east(f10B)
  call apply_reflection_edge_north_east(f11B)
  call apply_reflection_edge_north_east(f12B)
  call apply_reflection_edge_north_east(f13B)
  call apply_reflection_edge_north_east(f14B)
  call apply_reflection_edge_north_east(f15B)
  call apply_reflection_edge_north_east(f16B)
  call apply_reflection_edge_north_east(f17B)
  call apply_reflection_edge_north_east(f18B)
  
  !apply reflection at north front 4   !yz
  call apply_reflection_edge_north_front(f00B)
  call apply_reflection_edge_north_front(f01B)
  call apply_reflection_edge_north_front(f02B)
  call apply_reflection_edge_north_front(f03B)
  call apply_reflection_edge_north_front(f04B)
  call apply_reflection_edge_north_front(f05B)
  call apply_reflection_edge_north_front(f06B)
  call apply_reflection_edge_north_front(f07B)
  call apply_reflection_edge_north_front(f08B)
  call apply_reflection_edge_north_front(f09B)
  call apply_reflection_edge_north_front(f10B)
  call apply_reflection_edge_north_front(f11B)
  call apply_reflection_edge_north_front(f12B)
  call apply_reflection_edge_north_front(f13B)
  call apply_reflection_edge_north_front(f14B)
  call apply_reflection_edge_north_front(f15B)
  call apply_reflection_edge_north_front(f16B)
  call apply_reflection_edge_north_front(f17B)
  call apply_reflection_edge_north_front(f18B)
  
  !apply reflection at north rear 5   !yz
  call apply_reflection_edge_north_rear(f00B)
  call apply_reflection_edge_north_rear(f01B)
  call apply_reflection_edge_north_rear(f02B)
  call apply_reflection_edge_north_rear(f03B)
  call apply_reflection_edge_north_rear(f04B)
  call apply_reflection_edge_north_rear(f05B)
  call apply_reflection_edge_north_rear(f06B)
  call apply_reflection_edge_north_rear(f07B)
  call apply_reflection_edge_north_rear(f08B)
  call apply_reflection_edge_north_rear(f09B)
  call apply_reflection_edge_north_rear(f10B)
  call apply_reflection_edge_north_rear(f11B)
  call apply_reflection_edge_north_rear(f12B)
  call apply_reflection_edge_north_rear(f13B)
  call apply_reflection_edge_north_rear(f14B)
  call apply_reflection_edge_north_rear(f15B)
  call apply_reflection_edge_north_rear(f16B)
  call apply_reflection_edge_north_rear(f17B)
  call apply_reflection_edge_north_rear(f18B)
  
  !apply reflection at north west 6   !xz
  call apply_reflection_edge_north_west(f00B)
  call apply_reflection_edge_north_west(f01B)
  call apply_reflection_edge_north_west(f02B)
  call apply_reflection_edge_north_west(f03B)
  call apply_reflection_edge_north_west(f04B)
  call apply_reflection_edge_north_west(f05B)
  call apply_reflection_edge_north_west(f06B)
  call apply_reflection_edge_north_west(f07B)
  call apply_reflection_edge_north_west(f08B)
  call apply_reflection_edge_north_west(f09B)
  call apply_reflection_edge_north_west(f10B)
  call apply_reflection_edge_north_west(f11B)
  call apply_reflection_edge_north_west(f12B)
  call apply_reflection_edge_north_west(f13B)
  call apply_reflection_edge_north_west(f14B)
  call apply_reflection_edge_north_west(f15B)
  call apply_reflection_edge_north_west(f16B)
  call apply_reflection_edge_north_west(f17B)
  call apply_reflection_edge_north_west(f18B)
  
  !apply reflection at rear east 7   !xy
  call apply_reflection_edge_rear_east(f00B)
  call apply_reflection_edge_rear_east(f01B)
  call apply_reflection_edge_rear_east(f02B)
  call apply_reflection_edge_rear_east(f03B)
  call apply_reflection_edge_rear_east(f04B)
  call apply_reflection_edge_rear_east(f05B)
  call apply_reflection_edge_rear_east(f06B)
  call apply_reflection_edge_rear_east(f07B)
  call apply_reflection_edge_rear_east(f08B)
  call apply_reflection_edge_rear_east(f09B)
  call apply_reflection_edge_rear_east(f10B)
  call apply_reflection_edge_rear_east(f11B)
  call apply_reflection_edge_rear_east(f12B)
  call apply_reflection_edge_rear_east(f13B)
  call apply_reflection_edge_rear_east(f14B)
  call apply_reflection_edge_rear_east(f15B)
  call apply_reflection_edge_rear_east(f16B)
  call apply_reflection_edge_rear_east(f17B)
  call apply_reflection_edge_rear_east(f18B)
  
  !apply reflection at rear west 8   !xy
  call apply_reflection_edge_rear_west(f00B)
  call apply_reflection_edge_rear_west(f01B)
  call apply_reflection_edge_rear_west(f02B)
  call apply_reflection_edge_rear_west(f03B)
  call apply_reflection_edge_rear_west(f04B)
  call apply_reflection_edge_rear_west(f05B)
  call apply_reflection_edge_rear_west(f06B)
  call apply_reflection_edge_rear_west(f07B)
  call apply_reflection_edge_rear_west(f08B)
  call apply_reflection_edge_rear_west(f09B)
  call apply_reflection_edge_rear_west(f10B)
  call apply_reflection_edge_rear_west(f11B)
  call apply_reflection_edge_rear_west(f12B)
  call apply_reflection_edge_rear_west(f13B)
  call apply_reflection_edge_rear_west(f14B)
  call apply_reflection_edge_rear_west(f15B)
  call apply_reflection_edge_rear_west(f16B)
  call apply_reflection_edge_rear_west(f17B)
  call apply_reflection_edge_rear_west(f18B)
  
  !apply reflection at south east 9  !xz
  call apply_reflection_edge_south_east(f00B)
  call apply_reflection_edge_south_east(f01B)
  call apply_reflection_edge_south_east(f02B)
  call apply_reflection_edge_south_east(f03B)
  call apply_reflection_edge_south_east(f04B)
  call apply_reflection_edge_south_east(f05B)
  call apply_reflection_edge_south_east(f06B)
  call apply_reflection_edge_south_east(f07B)
  call apply_reflection_edge_south_east(f08B)
  call apply_reflection_edge_south_east(f09B)
  call apply_reflection_edge_south_east(f10B)
  call apply_reflection_edge_south_east(f11B)
  call apply_reflection_edge_south_east(f12B)
  call apply_reflection_edge_south_east(f13B)
  call apply_reflection_edge_south_east(f14B)
  call apply_reflection_edge_south_east(f15B)
  call apply_reflection_edge_south_east(f16B)
  call apply_reflection_edge_south_east(f17B)
  call apply_reflection_edge_south_east(f18B)
  
  !apply reflection at south front 10  !yz
  call apply_reflection_edge_south_front(f00B)
  call apply_reflection_edge_south_front(f01B)
  call apply_reflection_edge_south_front(f02B)
  call apply_reflection_edge_south_front(f03B)
  call apply_reflection_edge_south_front(f04B)
  call apply_reflection_edge_south_front(f05B)
  call apply_reflection_edge_south_front(f06B)
  call apply_reflection_edge_south_front(f07B)
  call apply_reflection_edge_south_front(f08B)
  call apply_reflection_edge_south_front(f09B)
  call apply_reflection_edge_south_front(f10B)
  call apply_reflection_edge_south_front(f11B)
  call apply_reflection_edge_south_front(f12B)
  call apply_reflection_edge_south_front(f13B)
  call apply_reflection_edge_south_front(f14B)
  call apply_reflection_edge_south_front(f15B)
  call apply_reflection_edge_south_front(f16B)
  call apply_reflection_edge_south_front(f17B)
  call apply_reflection_edge_south_front(f18B)
  
  !apply reflection at south rear 11  !yz
  call apply_reflection_edge_south_rear(f00B)
  call apply_reflection_edge_south_rear(f01B)
  call apply_reflection_edge_south_rear(f02B)
  call apply_reflection_edge_south_rear(f03B)
  call apply_reflection_edge_south_rear(f04B)
  call apply_reflection_edge_south_rear(f05B)
  call apply_reflection_edge_south_rear(f06B)
  call apply_reflection_edge_south_rear(f07B)
  call apply_reflection_edge_south_rear(f08B)
  call apply_reflection_edge_south_rear(f09B)
  call apply_reflection_edge_south_rear(f10B)
  call apply_reflection_edge_south_rear(f11B)
  call apply_reflection_edge_south_rear(f12B)
  call apply_reflection_edge_south_rear(f13B)
  call apply_reflection_edge_south_rear(f14B)
  call apply_reflection_edge_south_rear(f15B)
  call apply_reflection_edge_south_rear(f16B)
  call apply_reflection_edge_south_rear(f17B)
  call apply_reflection_edge_south_rear(f18B)
  
  !apply reflection at south west 12  !xz
  call apply_reflection_edge_south_west(f00B)
  call apply_reflection_edge_south_west(f01B)
  call apply_reflection_edge_south_west(f02B)
  call apply_reflection_edge_south_west(f03B)
  call apply_reflection_edge_south_west(f04B)
  call apply_reflection_edge_south_west(f05B)
  call apply_reflection_edge_south_west(f06B)
  call apply_reflection_edge_south_west(f07B)
  call apply_reflection_edge_south_west(f08B)
  call apply_reflection_edge_south_west(f09B)
  call apply_reflection_edge_south_west(f10B)
  call apply_reflection_edge_south_west(f11B)
  call apply_reflection_edge_south_west(f12B)
  call apply_reflection_edge_south_west(f13B)
  call apply_reflection_edge_south_west(f14B)
  call apply_reflection_edge_south_west(f15B)
  call apply_reflection_edge_south_west(f16B)
  call apply_reflection_edge_south_west(f17B)
  call apply_reflection_edge_south_west(f18B)
  
  !corner 8
  
  !apply reflection at north east front 1
  call apply_reflection_corner_north_east_front(f00B)
  call apply_reflection_corner_north_east_front(f01B)
  call apply_reflection_corner_north_east_front(f02B)
  call apply_reflection_corner_north_east_front(f03B)
  call apply_reflection_corner_north_east_front(f04B)
  call apply_reflection_corner_north_east_front(f05B)
  call apply_reflection_corner_north_east_front(f06B)
  call apply_reflection_corner_north_east_front(f07B)
  call apply_reflection_corner_north_east_front(f08B)
  call apply_reflection_corner_north_east_front(f09B)
  call apply_reflection_corner_north_east_front(f10B)
  call apply_reflection_corner_north_east_front(f11B)
  call apply_reflection_corner_north_east_front(f12B)
  call apply_reflection_corner_north_east_front(f13B)
  call apply_reflection_corner_north_east_front(f14B)
  call apply_reflection_corner_north_east_front(f15B)
  call apply_reflection_corner_north_east_front(f16B)
  call apply_reflection_corner_north_east_front(f17B)
  call apply_reflection_corner_north_east_front(f18B)
  
  !apply reflection at north east rear 2
  call apply_reflection_corner_north_east_rear(f00B)
  call apply_reflection_corner_north_east_rear(f01B)
  call apply_reflection_corner_north_east_rear(f02B)
  call apply_reflection_corner_north_east_rear(f03B)
  call apply_reflection_corner_north_east_rear(f04B)
  call apply_reflection_corner_north_east_rear(f05B)
  call apply_reflection_corner_north_east_rear(f06B)
  call apply_reflection_corner_north_east_rear(f07B)
  call apply_reflection_corner_north_east_rear(f08B)
  call apply_reflection_corner_north_east_rear(f09B)
  call apply_reflection_corner_north_east_rear(f10B)
  call apply_reflection_corner_north_east_rear(f11B)
  call apply_reflection_corner_north_east_rear(f12B)
  call apply_reflection_corner_north_east_rear(f13B)
  call apply_reflection_corner_north_east_rear(f14B)
  call apply_reflection_corner_north_east_rear(f15B)
  call apply_reflection_corner_north_east_rear(f16B)
  call apply_reflection_corner_north_east_rear(f17B)
  call apply_reflection_corner_north_east_rear(f18B)
  
  !apply reflection at north west rear 3
  call apply_reflection_corner_north_west_rear(f00B)
  call apply_reflection_corner_north_west_rear(f01B)
  call apply_reflection_corner_north_west_rear(f02B)
  call apply_reflection_corner_north_west_rear(f03B)
  call apply_reflection_corner_north_west_rear(f04B)
  call apply_reflection_corner_north_west_rear(f05B)
  call apply_reflection_corner_north_west_rear(f06B)
  call apply_reflection_corner_north_west_rear(f07B)
  call apply_reflection_corner_north_west_rear(f08B)
  call apply_reflection_corner_north_west_rear(f09B)
  call apply_reflection_corner_north_west_rear(f10B)
  call apply_reflection_corner_north_west_rear(f11B)
  call apply_reflection_corner_north_west_rear(f12B)
  call apply_reflection_corner_north_west_rear(f13B)
  call apply_reflection_corner_north_west_rear(f14B)
  call apply_reflection_corner_north_west_rear(f15B)
  call apply_reflection_corner_north_west_rear(f16B)
  call apply_reflection_corner_north_west_rear(f17B)
  call apply_reflection_corner_north_west_rear(f18B)
  
  !apply reflection at north west front 4
  call apply_reflection_corner_north_west_front(f00B)
  call apply_reflection_corner_north_west_front(f01B)
  call apply_reflection_corner_north_west_front(f02B)
  call apply_reflection_corner_north_west_front(f03B)
  call apply_reflection_corner_north_west_front(f04B)
  call apply_reflection_corner_north_west_front(f05B)
  call apply_reflection_corner_north_west_front(f06B)
  call apply_reflection_corner_north_west_front(f07B)
  call apply_reflection_corner_north_west_front(f08B)
  call apply_reflection_corner_north_west_front(f09B)
  call apply_reflection_corner_north_west_front(f10B)
  call apply_reflection_corner_north_west_front(f11B)
  call apply_reflection_corner_north_west_front(f12B)
  call apply_reflection_corner_north_west_front(f13B)
  call apply_reflection_corner_north_west_front(f14B)
  call apply_reflection_corner_north_west_front(f15B)
  call apply_reflection_corner_north_west_front(f16B)
  call apply_reflection_corner_north_west_front(f17B)
  call apply_reflection_corner_north_west_front(f18B)
  
  !apply reflection at south east front 5
  call apply_reflection_corner_south_east_front(f00B)
  call apply_reflection_corner_south_east_front(f01B)
  call apply_reflection_corner_south_east_front(f02B)
  call apply_reflection_corner_south_east_front(f03B)
  call apply_reflection_corner_south_east_front(f04B)
  call apply_reflection_corner_south_east_front(f05B)
  call apply_reflection_corner_south_east_front(f06B)
  call apply_reflection_corner_south_east_front(f07B)
  call apply_reflection_corner_south_east_front(f08B)
  call apply_reflection_corner_south_east_front(f09B)
  call apply_reflection_corner_south_east_front(f10B)
  call apply_reflection_corner_south_east_front(f11B)
  call apply_reflection_corner_south_east_front(f12B)
  call apply_reflection_corner_south_east_front(f13B)
  call apply_reflection_corner_south_east_front(f14B)
  call apply_reflection_corner_south_east_front(f15B)
  call apply_reflection_corner_south_east_front(f16B)
  call apply_reflection_corner_south_east_front(f17B)
  call apply_reflection_corner_south_east_front(f18B)
  
  !apply reflection at south west front 6
  call apply_reflection_corner_south_west_front(f00B)
  call apply_reflection_corner_south_west_front(f01B)
  call apply_reflection_corner_south_west_front(f02B)
  call apply_reflection_corner_south_west_front(f03B)
  call apply_reflection_corner_south_west_front(f04B)
  call apply_reflection_corner_south_west_front(f05B)
  call apply_reflection_corner_south_west_front(f06B)
  call apply_reflection_corner_south_west_front(f07B)
  call apply_reflection_corner_south_west_front(f08B)
  call apply_reflection_corner_south_west_front(f09B)
  call apply_reflection_corner_south_west_front(f10B)
  call apply_reflection_corner_south_west_front(f11B)
  call apply_reflection_corner_south_west_front(f12B)
  call apply_reflection_corner_south_west_front(f13B)
  call apply_reflection_corner_south_west_front(f14B)
  call apply_reflection_corner_south_west_front(f15B)
  call apply_reflection_corner_south_west_front(f16B)
  call apply_reflection_corner_south_west_front(f17B)
  call apply_reflection_corner_south_west_front(f18B)
  
  !apply reflection at south west rear 7
  call apply_reflection_corner_south_west_rear(f00B)
  call apply_reflection_corner_south_west_rear(f01B)
  call apply_reflection_corner_south_west_rear(f02B)
  call apply_reflection_corner_south_west_rear(f03B)
  call apply_reflection_corner_south_west_rear(f04B)
  call apply_reflection_corner_south_west_rear(f05B)
  call apply_reflection_corner_south_west_rear(f06B)
  call apply_reflection_corner_south_west_rear(f07B)
  call apply_reflection_corner_south_west_rear(f08B)
  call apply_reflection_corner_south_west_rear(f09B)
  call apply_reflection_corner_south_west_rear(f10B)
  call apply_reflection_corner_south_west_rear(f11B)
  call apply_reflection_corner_south_west_rear(f12B)
  call apply_reflection_corner_south_west_rear(f13B)
  call apply_reflection_corner_south_west_rear(f14B)
  call apply_reflection_corner_south_west_rear(f15B)
  call apply_reflection_corner_south_west_rear(f16B)
  call apply_reflection_corner_south_west_rear(f17B)
  call apply_reflection_corner_south_west_rear(f18B)
  
  !apply reflection at south east rear 8
  call apply_reflection_corner_south_east_rear(f00B)
  call apply_reflection_corner_south_east_rear(f01B)
  call apply_reflection_corner_south_east_rear(f02B)
  call apply_reflection_corner_south_east_rear(f03B)
  call apply_reflection_corner_south_east_rear(f04B)
  call apply_reflection_corner_south_east_rear(f05B)
  call apply_reflection_corner_south_east_rear(f06B)
  call apply_reflection_corner_south_east_rear(f07B)
  call apply_reflection_corner_south_east_rear(f08B)
  call apply_reflection_corner_south_east_rear(f09B)
  call apply_reflection_corner_south_east_rear(f10B)
  call apply_reflection_corner_south_east_rear(f11B)
  call apply_reflection_corner_south_east_rear(f12B)
  call apply_reflection_corner_south_east_rear(f13B)
  call apply_reflection_corner_south_east_rear(f14B)
  call apply_reflection_corner_south_east_rear(f15B)
  call apply_reflection_corner_south_east_rear(f16B)
  call apply_reflection_corner_south_east_rear(f17B)
  call apply_reflection_corner_south_east_rear(f18B)
  
  return
  
 end subroutine apply_reflection_pops_all
 
 subroutine apply_reflection_pops_along_yz
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     to fluid populations along y and z
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  !sides 4 (6)
  
  !apply pbc at front 5 !y
  !red fluid
  call apply_reflection_front_frame_x(f00R)
  call apply_reflection_front_frame_x(f01R)
  call apply_reflection_front_frame_x(f02R)
  call apply_reflection_front_frame_x(f03R)
  call apply_reflection_front_frame_x(f04R)
  call apply_reflection_front_frame_x(f05R)
  call apply_reflection_front_frame_x(f06R)
  call apply_reflection_front_frame_x(f07R)
  call apply_reflection_front_frame_x(f08R)
  call apply_reflection_front_frame_x(f09R)
  call apply_reflection_front_frame_x(f10R)
  call apply_reflection_front_frame_x(f11R)
  call apply_reflection_front_frame_x(f12R)
  call apply_reflection_front_frame_x(f13R)
  call apply_reflection_front_frame_x(f14R)
  call apply_reflection_front_frame_x(f15R)
  call apply_reflection_front_frame_x(f16R)
  call apply_reflection_front_frame_x(f17R)
  call apply_reflection_front_frame_x(f18R)
  
  !apply pbc at rear 6 !y
  !red fluid
  call apply_reflection_rear_frame_x(f00R)
  call apply_reflection_rear_frame_x(f01R)
  call apply_reflection_rear_frame_x(f02R)
  call apply_reflection_rear_frame_x(f03R)
  call apply_reflection_rear_frame_x(f04R)
  call apply_reflection_rear_frame_x(f05R)
  call apply_reflection_rear_frame_x(f06R)
  call apply_reflection_rear_frame_x(f07R)
  call apply_reflection_rear_frame_x(f08R)
  call apply_reflection_rear_frame_x(f09R)
  call apply_reflection_rear_frame_x(f10R)
  call apply_reflection_rear_frame_x(f11R)
  call apply_reflection_rear_frame_x(f12R)
  call apply_reflection_rear_frame_x(f13R)
  call apply_reflection_rear_frame_x(f14R)
  call apply_reflection_rear_frame_x(f15R)
  call apply_reflection_rear_frame_x(f16R)
  call apply_reflection_rear_frame_x(f17R)
  call apply_reflection_rear_frame_x(f18R)
  
  !apply pbc at north 1 !z
  !red fluid
  !frame x
  call apply_reflection_north(1-nbuff,nx+nbuff,1,ny,f00R)
  call apply_reflection_north(1-nbuff,nx+nbuff,1,ny,f01R)
  call apply_reflection_north(1-nbuff,nx+nbuff,1,ny,f02R)
  call apply_reflection_north(1-nbuff,nx+nbuff,1,ny,f03R)
  call apply_reflection_north(1-nbuff,nx+nbuff,1,ny,f04R)
  call apply_reflection_north(1-nbuff,nx+nbuff,1,ny,f05R)
  call apply_reflection_north(1-nbuff,nx+nbuff,1,ny,f06R)
  call apply_reflection_north(1-nbuff,nx+nbuff,1,ny,f07R)
  call apply_reflection_north(1-nbuff,nx+nbuff,1,ny,f08R)
  call apply_reflection_north(1-nbuff,nx+nbuff,1,ny,f09R)
  call apply_reflection_north(1-nbuff,nx+nbuff,1,ny,f10R)
  call apply_reflection_north(1-nbuff,nx+nbuff,1,ny,f11R)
  call apply_reflection_north(1-nbuff,nx+nbuff,1,ny,f12R)
  call apply_reflection_north(1-nbuff,nx+nbuff,1,ny,f13R)
  call apply_reflection_north(1-nbuff,nx+nbuff,1,ny,f14R)
  call apply_reflection_north(1-nbuff,nx+nbuff,1,ny,f15R)
  call apply_reflection_north(1-nbuff,nx+nbuff,1,ny,f16R)
  call apply_reflection_north(1-nbuff,nx+nbuff,1,ny,f17R)
  call apply_reflection_north(1-nbuff,nx+nbuff,1,ny,f18R)
  
  !apply pbc at south 2 !z
  !red fluid
  call apply_reflection_south_frame_x(f00R)
  call apply_reflection_south_frame_x(f01R)
  call apply_reflection_south_frame_x(f02R)
  call apply_reflection_south_frame_x(f03R)
  call apply_reflection_south_frame_x(f04R)
  call apply_reflection_south_frame_x(f05R)
  call apply_reflection_south_frame_x(f06R)
  call apply_reflection_south_frame_x(f07R)
  call apply_reflection_south_frame_x(f08R)
  call apply_reflection_south_frame_x(f09R)
  call apply_reflection_south_frame_x(f10R)
  call apply_reflection_south_frame_x(f11R)
  call apply_reflection_south_frame_x(f12R)
  call apply_reflection_south_frame_x(f13R)
  call apply_reflection_south_frame_x(f14R)
  call apply_reflection_south_frame_x(f15R)
  call apply_reflection_south_frame_x(f16R)
  call apply_reflection_south_frame_x(f17R)
  call apply_reflection_south_frame_x(f18R)
  
  !edges 4 (12)
  
  !apply reflection at north front 4   !yz
  !red fluid
  call apply_reflection_edge_north_front_frame(f00R)
  call apply_reflection_edge_north_front_frame(f01R)
  call apply_reflection_edge_north_front_frame(f02R)
  call apply_reflection_edge_north_front_frame(f03R)
  call apply_reflection_edge_north_front_frame(f04R)
  call apply_reflection_edge_north_front_frame(f05R)
  call apply_reflection_edge_north_front_frame(f06R)
  call apply_reflection_edge_north_front_frame(f07R)
  call apply_reflection_edge_north_front_frame(f08R)
  call apply_reflection_edge_north_front_frame(f09R)
  call apply_reflection_edge_north_front_frame(f10R)
  call apply_reflection_edge_north_front_frame(f11R)
  call apply_reflection_edge_north_front_frame(f12R)
  call apply_reflection_edge_north_front_frame(f13R)
  call apply_reflection_edge_north_front_frame(f14R)
  call apply_reflection_edge_north_front_frame(f15R)
  call apply_reflection_edge_north_front_frame(f16R)
  call apply_reflection_edge_north_front_frame(f17R)
  call apply_reflection_edge_north_front_frame(f18R)
  
  !apply reflection at north rear 5   !yz
  !red fluid
  call apply_reflection_edge_north_rear_frame(f00R)
  call apply_reflection_edge_north_rear_frame(f01R)
  call apply_reflection_edge_north_rear_frame(f02R)
  call apply_reflection_edge_north_rear_frame(f03R)
  call apply_reflection_edge_north_rear_frame(f04R)
  call apply_reflection_edge_north_rear_frame(f05R)
  call apply_reflection_edge_north_rear_frame(f06R)
  call apply_reflection_edge_north_rear_frame(f07R)
  call apply_reflection_edge_north_rear_frame(f08R)
  call apply_reflection_edge_north_rear_frame(f09R)
  call apply_reflection_edge_north_rear_frame(f10R)
  call apply_reflection_edge_north_rear_frame(f11R)
  call apply_reflection_edge_north_rear_frame(f12R)
  call apply_reflection_edge_north_rear_frame(f13R)
  call apply_reflection_edge_north_rear_frame(f14R)
  call apply_reflection_edge_north_rear_frame(f15R)
  call apply_reflection_edge_north_rear_frame(f16R)
  call apply_reflection_edge_north_rear_frame(f17R)
  call apply_reflection_edge_north_rear_frame(f18R)
  
  !apply reflection at south front 10  !yz
  !red fluid
  call apply_reflection_edge_south_front_frame(f00R)
  call apply_reflection_edge_south_front_frame(f01R)
  call apply_reflection_edge_south_front_frame(f02R)
  call apply_reflection_edge_south_front_frame(f03R)
  call apply_reflection_edge_south_front_frame(f04R)
  call apply_reflection_edge_south_front_frame(f05R)
  call apply_reflection_edge_south_front_frame(f06R)
  call apply_reflection_edge_south_front_frame(f07R)
  call apply_reflection_edge_south_front_frame(f08R)
  call apply_reflection_edge_south_front_frame(f09R)
  call apply_reflection_edge_south_front_frame(f10R)
  call apply_reflection_edge_south_front_frame(f11R)
  call apply_reflection_edge_south_front_frame(f12R)
  call apply_reflection_edge_south_front_frame(f13R)
  call apply_reflection_edge_south_front_frame(f14R)
  call apply_reflection_edge_south_front_frame(f15R)
  call apply_reflection_edge_south_front_frame(f16R)
  call apply_reflection_edge_south_front_frame(f17R)
  call apply_reflection_edge_south_front_frame(f18R)
  
  !apply reflection at south rear 11  !yz
  !red fluid
  call apply_reflection_edge_south_rear_frame(f00R)
  call apply_reflection_edge_south_rear_frame(f01R)
  call apply_reflection_edge_south_rear_frame(f02R)
  call apply_reflection_edge_south_rear_frame(f03R)
  call apply_reflection_edge_south_rear_frame(f04R)
  call apply_reflection_edge_south_rear_frame(f05R)
  call apply_reflection_edge_south_rear_frame(f06R)
  call apply_reflection_edge_south_rear_frame(f07R)
  call apply_reflection_edge_south_rear_frame(f08R)
  call apply_reflection_edge_south_rear_frame(f09R)
  call apply_reflection_edge_south_rear_frame(f10R)
  call apply_reflection_edge_south_rear_frame(f11R)
  call apply_reflection_edge_south_rear_frame(f12R)
  call apply_reflection_edge_south_rear_frame(f13R)
  call apply_reflection_edge_south_rear_frame(f14R)
  call apply_reflection_edge_south_rear_frame(f15R)
  call apply_reflection_edge_south_rear_frame(f16R)
  call apply_reflection_edge_south_rear_frame(f17R)
  call apply_reflection_edge_south_rear_frame(f18R)
  
  if(lsingle_fluid)return
  
  !sides 4 (6)
  
  !apply pbc at front 5 !y
  !blue fluid
  call apply_reflection_front_frame_x(f00B)
  call apply_reflection_front_frame_x(f01B)
  call apply_reflection_front_frame_x(f02B)
  call apply_reflection_front_frame_x(f03B)
  call apply_reflection_front_frame_x(f04B)
  call apply_reflection_front_frame_x(f05B)
  call apply_reflection_front_frame_x(f06B)
  call apply_reflection_front_frame_x(f07B)
  call apply_reflection_front_frame_x(f08B)
  call apply_reflection_front_frame_x(f09B)
  call apply_reflection_front_frame_x(f10B)
  call apply_reflection_front_frame_x(f11B)
  call apply_reflection_front_frame_x(f12B)
  call apply_reflection_front_frame_x(f13B)
  call apply_reflection_front_frame_x(f14B)
  call apply_reflection_front_frame_x(f15B)
  call apply_reflection_front_frame_x(f16B)
  call apply_reflection_front_frame_x(f17B)
  call apply_reflection_front_frame_x(f18B)
  
  !apply pbc at rear 6 !y
  !blue fluid
  call apply_reflection_rear_frame_x(f00B)
  call apply_reflection_rear_frame_x(f01B)
  call apply_reflection_rear_frame_x(f02B)
  call apply_reflection_rear_frame_x(f03B)
  call apply_reflection_rear_frame_x(f04B)
  call apply_reflection_rear_frame_x(f05B)
  call apply_reflection_rear_frame_x(f06B)
  call apply_reflection_rear_frame_x(f07B)
  call apply_reflection_rear_frame_x(f08B)
  call apply_reflection_rear_frame_x(f09B)
  call apply_reflection_rear_frame_x(f10B)
  call apply_reflection_rear_frame_x(f11B)
  call apply_reflection_rear_frame_x(f12B)
  call apply_reflection_rear_frame_x(f13B)
  call apply_reflection_rear_frame_x(f14B)
  call apply_reflection_rear_frame_x(f15B)
  call apply_reflection_rear_frame_x(f16B)
  call apply_reflection_rear_frame_x(f17B)
  call apply_reflection_rear_frame_x(f18B)
  
  !apply pbc at north 1 !z
  !blue fluid
  !frame x
  call apply_reflection_north(1-nbuff,nx+nbuff,1,ny,f00B)
  call apply_reflection_north(1-nbuff,nx+nbuff,1,ny,f01B)
  call apply_reflection_north(1-nbuff,nx+nbuff,1,ny,f02B)
  call apply_reflection_north(1-nbuff,nx+nbuff,1,ny,f03B)
  call apply_reflection_north(1-nbuff,nx+nbuff,1,ny,f04B)
  call apply_reflection_north(1-nbuff,nx+nbuff,1,ny,f05B)
  call apply_reflection_north(1-nbuff,nx+nbuff,1,ny,f06B)
  call apply_reflection_north(1-nbuff,nx+nbuff,1,ny,f07B)
  call apply_reflection_north(1-nbuff,nx+nbuff,1,ny,f08B)
  call apply_reflection_north(1-nbuff,nx+nbuff,1,ny,f09B)
  call apply_reflection_north(1-nbuff,nx+nbuff,1,ny,f10B)
  call apply_reflection_north(1-nbuff,nx+nbuff,1,ny,f11B)
  call apply_reflection_north(1-nbuff,nx+nbuff,1,ny,f12B)
  call apply_reflection_north(1-nbuff,nx+nbuff,1,ny,f13B)
  call apply_reflection_north(1-nbuff,nx+nbuff,1,ny,f14B)
  call apply_reflection_north(1-nbuff,nx+nbuff,1,ny,f15B)
  call apply_reflection_north(1-nbuff,nx+nbuff,1,ny,f16B)
  call apply_reflection_north(1-nbuff,nx+nbuff,1,ny,f17B)
  call apply_reflection_north(1-nbuff,nx+nbuff,1,ny,f18B)
  
  !apply pbc at south 2 !z
  !blue fluid
  call apply_reflection_south_frame_x(f00B)
  call apply_reflection_south_frame_x(f01B)
  call apply_reflection_south_frame_x(f02B)
  call apply_reflection_south_frame_x(f03B)
  call apply_reflection_south_frame_x(f04B)
  call apply_reflection_south_frame_x(f05B)
  call apply_reflection_south_frame_x(f06B)
  call apply_reflection_south_frame_x(f07B)
  call apply_reflection_south_frame_x(f08B)
  call apply_reflection_south_frame_x(f09B)
  call apply_reflection_south_frame_x(f10B)
  call apply_reflection_south_frame_x(f11B)
  call apply_reflection_south_frame_x(f12B)
  call apply_reflection_south_frame_x(f13B)
  call apply_reflection_south_frame_x(f14B)
  call apply_reflection_south_frame_x(f15B)
  call apply_reflection_south_frame_x(f16B)
  call apply_reflection_south_frame_x(f17B)
  call apply_reflection_south_frame_x(f18B)
  
  !edges 4 (12)
  
  !apply reflection at north front 4   !yz
  !blue fluid
  call apply_reflection_edge_north_front_frame(f00B)
  call apply_reflection_edge_north_front_frame(f01B)
  call apply_reflection_edge_north_front_frame(f02B)
  call apply_reflection_edge_north_front_frame(f03B)
  call apply_reflection_edge_north_front_frame(f04B)
  call apply_reflection_edge_north_front_frame(f05B)
  call apply_reflection_edge_north_front_frame(f06B)
  call apply_reflection_edge_north_front_frame(f07B)
  call apply_reflection_edge_north_front_frame(f08B)
  call apply_reflection_edge_north_front_frame(f09B)
  call apply_reflection_edge_north_front_frame(f10B)
  call apply_reflection_edge_north_front_frame(f11B)
  call apply_reflection_edge_north_front_frame(f12B)
  call apply_reflection_edge_north_front_frame(f13B)
  call apply_reflection_edge_north_front_frame(f14B)
  call apply_reflection_edge_north_front_frame(f15B)
  call apply_reflection_edge_north_front_frame(f16B)
  call apply_reflection_edge_north_front_frame(f17B)
  call apply_reflection_edge_north_front_frame(f18B)
  
  !apply reflection at north rear 5   !yz
  !blue fluid
  call apply_reflection_edge_north_rear_frame(f00B)
  call apply_reflection_edge_north_rear_frame(f01B)
  call apply_reflection_edge_north_rear_frame(f02B)
  call apply_reflection_edge_north_rear_frame(f03B)
  call apply_reflection_edge_north_rear_frame(f04B)
  call apply_reflection_edge_north_rear_frame(f05B)
  call apply_reflection_edge_north_rear_frame(f06B)
  call apply_reflection_edge_north_rear_frame(f07B)
  call apply_reflection_edge_north_rear_frame(f08B)
  call apply_reflection_edge_north_rear_frame(f09B)
  call apply_reflection_edge_north_rear_frame(f10B)
  call apply_reflection_edge_north_rear_frame(f11B)
  call apply_reflection_edge_north_rear_frame(f12B)
  call apply_reflection_edge_north_rear_frame(f13B)
  call apply_reflection_edge_north_rear_frame(f14B)
  call apply_reflection_edge_north_rear_frame(f15B)
  call apply_reflection_edge_north_rear_frame(f16B)
  call apply_reflection_edge_north_rear_frame(f17B)
  call apply_reflection_edge_north_rear_frame(f18B)
  
  !apply reflection at south front 10  !yz
  !blue fluid
  call apply_reflection_edge_south_front_frame(f00B)
  call apply_reflection_edge_south_front_frame(f01B)
  call apply_reflection_edge_south_front_frame(f02B)
  call apply_reflection_edge_south_front_frame(f03B)
  call apply_reflection_edge_south_front_frame(f04B)
  call apply_reflection_edge_south_front_frame(f05B)
  call apply_reflection_edge_south_front_frame(f06B)
  call apply_reflection_edge_south_front_frame(f07B)
  call apply_reflection_edge_south_front_frame(f08B)
  call apply_reflection_edge_south_front_frame(f09B)
  call apply_reflection_edge_south_front_frame(f10B)
  call apply_reflection_edge_south_front_frame(f11B)
  call apply_reflection_edge_south_front_frame(f12B)
  call apply_reflection_edge_south_front_frame(f13B)
  call apply_reflection_edge_south_front_frame(f14B)
  call apply_reflection_edge_south_front_frame(f15B)
  call apply_reflection_edge_south_front_frame(f16B)
  call apply_reflection_edge_south_front_frame(f17B)
  call apply_reflection_edge_south_front_frame(f18B)
  
  !apply reflection at south rear 11  !yz
  !blue fluid
  call apply_reflection_edge_south_rear_frame(f00B)
  call apply_reflection_edge_south_rear_frame(f01B)
  call apply_reflection_edge_south_rear_frame(f02B)
  call apply_reflection_edge_south_rear_frame(f03B)
  call apply_reflection_edge_south_rear_frame(f04B)
  call apply_reflection_edge_south_rear_frame(f05B)
  call apply_reflection_edge_south_rear_frame(f06B)
  call apply_reflection_edge_south_rear_frame(f07B)
  call apply_reflection_edge_south_rear_frame(f08B)
  call apply_reflection_edge_south_rear_frame(f09B)
  call apply_reflection_edge_south_rear_frame(f10B)
  call apply_reflection_edge_south_rear_frame(f11B)
  call apply_reflection_edge_south_rear_frame(f12B)
  call apply_reflection_edge_south_rear_frame(f13B)
  call apply_reflection_edge_south_rear_frame(f14B)
  call apply_reflection_edge_south_rear_frame(f15B)
  call apply_reflection_edge_south_rear_frame(f16B)
  call apply_reflection_edge_south_rear_frame(f17B)
  call apply_reflection_edge_south_rear_frame(f18B)
  
  return
  
 end subroutine apply_reflection_pops_along_yz
 
 subroutine apply_reflection_pops_along_xz
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     to fluid populations along x and z
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  !sides 6
  
  !apply pbc at east 3 !x
  !red fluid
  call apply_reflection_east(1-nbuff,ny+nbuff,1,nz,f00R)
  call apply_reflection_east(1-nbuff,ny+nbuff,1,nz,f01R)
  call apply_reflection_east(1-nbuff,ny+nbuff,1,nz,f02R)
  call apply_reflection_east(1-nbuff,ny+nbuff,1,nz,f03R)
  call apply_reflection_east(1-nbuff,ny+nbuff,1,nz,f04R)
  call apply_reflection_east(1-nbuff,ny+nbuff,1,nz,f05R)
  call apply_reflection_east(1-nbuff,ny+nbuff,1,nz,f06R)
  call apply_reflection_east(1-nbuff,ny+nbuff,1,nz,f07R)
  call apply_reflection_east(1-nbuff,ny+nbuff,1,nz,f08R)
  call apply_reflection_east(1-nbuff,ny+nbuff,1,nz,f09R)
  call apply_reflection_east(1-nbuff,ny+nbuff,1,nz,f10R)
  call apply_reflection_east(1-nbuff,ny+nbuff,1,nz,f11R)
  call apply_reflection_east(1-nbuff,ny+nbuff,1,nz,f12R)
  call apply_reflection_east(1-nbuff,ny+nbuff,1,nz,f13R)
  call apply_reflection_east(1-nbuff,ny+nbuff,1,nz,f14R)
  call apply_reflection_east(1-nbuff,ny+nbuff,1,nz,f15R)
  call apply_reflection_east(1-nbuff,ny+nbuff,1,nz,f16R)
  call apply_reflection_east(1-nbuff,ny+nbuff,1,nz,f17R)
  call apply_reflection_east(1-nbuff,ny+nbuff,1,nz,f18R)
  
  !apply pbc at west 4 !x
  !red fluid
  call apply_reflection_west(1-nbuff,ny+nbuff,1,nz,f00R)
  call apply_reflection_west(1-nbuff,ny+nbuff,1,nz,f01R)
  call apply_reflection_west(1-nbuff,ny+nbuff,1,nz,f02R)
  call apply_reflection_west(1-nbuff,ny+nbuff,1,nz,f03R)
  call apply_reflection_west(1-nbuff,ny+nbuff,1,nz,f04R)
  call apply_reflection_west(1-nbuff,ny+nbuff,1,nz,f05R)
  call apply_reflection_west(1-nbuff,ny+nbuff,1,nz,f06R)
  call apply_reflection_west(1-nbuff,ny+nbuff,1,nz,f07R)
  call apply_reflection_west(1-nbuff,ny+nbuff,1,nz,f08R)
  call apply_reflection_west(1-nbuff,ny+nbuff,1,nz,f09R)
  call apply_reflection_west(1-nbuff,ny+nbuff,1,nz,f10R)
  call apply_reflection_west(1-nbuff,ny+nbuff,1,nz,f11R)
  call apply_reflection_west(1-nbuff,ny+nbuff,1,nz,f12R)
  call apply_reflection_west(1-nbuff,ny+nbuff,1,nz,f13R)
  call apply_reflection_west(1-nbuff,ny+nbuff,1,nz,f14R)
  call apply_reflection_west(1-nbuff,ny+nbuff,1,nz,f15R)
  call apply_reflection_west(1-nbuff,ny+nbuff,1,nz,f16R)
  call apply_reflection_west(1-nbuff,ny+nbuff,1,nz,f17R)
  call apply_reflection_west(1-nbuff,ny+nbuff,1,nz,f18R)
  
  !apply pbc at north 1 !z
  !red fluid
  !frame y
  call apply_reflection_north(1,nx,1-nbuff,ny+nbuff,f00R)
  call apply_reflection_north(1,nx,1-nbuff,ny+nbuff,f01R)
  call apply_reflection_north(1,nx,1-nbuff,ny+nbuff,f02R)
  call apply_reflection_north(1,nx,1-nbuff,ny+nbuff,f03R)
  call apply_reflection_north(1,nx,1-nbuff,ny+nbuff,f04R)
  call apply_reflection_north(1,nx,1-nbuff,ny+nbuff,f05R)
  call apply_reflection_north(1,nx,1-nbuff,ny+nbuff,f06R)
  call apply_reflection_north(1,nx,1-nbuff,ny+nbuff,f07R)
  call apply_reflection_north(1,nx,1-nbuff,ny+nbuff,f08R)
  call apply_reflection_north(1,nx,1-nbuff,ny+nbuff,f09R)
  call apply_reflection_north(1,nx,1-nbuff,ny+nbuff,f10R)
  call apply_reflection_north(1,nx,1-nbuff,ny+nbuff,f11R)
  call apply_reflection_north(1,nx,1-nbuff,ny+nbuff,f12R)
  call apply_reflection_north(1,nx,1-nbuff,ny+nbuff,f13R)
  call apply_reflection_north(1,nx,1-nbuff,ny+nbuff,f14R)
  call apply_reflection_north(1,nx,1-nbuff,ny+nbuff,f15R)
  call apply_reflection_north(1,nx,1-nbuff,ny+nbuff,f16R)
  call apply_reflection_north(1,nx,1-nbuff,ny+nbuff,f17R)
  call apply_reflection_north(1,nx,1-nbuff,ny+nbuff,f18R)
  
  !apply pbc at south 2 !z
  !red fluid
  call apply_reflection_south_frame_y(f00R)
  call apply_reflection_south_frame_y(f01R)
  call apply_reflection_south_frame_y(f02R)
  call apply_reflection_south_frame_y(f03R)
  call apply_reflection_south_frame_y(f04R)
  call apply_reflection_south_frame_y(f05R)
  call apply_reflection_south_frame_y(f06R)
  call apply_reflection_south_frame_y(f07R)
  call apply_reflection_south_frame_y(f08R)
  call apply_reflection_south_frame_y(f09R)
  call apply_reflection_south_frame_y(f10R)
  call apply_reflection_south_frame_y(f11R)
  call apply_reflection_south_frame_y(f12R)
  call apply_reflection_south_frame_y(f13R)
  call apply_reflection_south_frame_y(f14R)
  call apply_reflection_south_frame_y(f15R)
  call apply_reflection_south_frame_y(f16R)
  call apply_reflection_south_frame_y(f17R)
  call apply_reflection_south_frame_y(f18R)
  
  !edges 4 (12)
  
  !apply reflection at north east 3   !xz
  !red fluid
  call apply_reflection_edge_north_east_frame(f00R)
  call apply_reflection_edge_north_east_frame(f01R)
  call apply_reflection_edge_north_east_frame(f02R)
  call apply_reflection_edge_north_east_frame(f03R)
  call apply_reflection_edge_north_east_frame(f04R)
  call apply_reflection_edge_north_east_frame(f05R)
  call apply_reflection_edge_north_east_frame(f06R)
  call apply_reflection_edge_north_east_frame(f07R)
  call apply_reflection_edge_north_east_frame(f08R)
  call apply_reflection_edge_north_east_frame(f09R)
  call apply_reflection_edge_north_east_frame(f10R)
  call apply_reflection_edge_north_east_frame(f11R)
  call apply_reflection_edge_north_east_frame(f12R)
  call apply_reflection_edge_north_east_frame(f13R)
  call apply_reflection_edge_north_east_frame(f14R)
  call apply_reflection_edge_north_east_frame(f15R)
  call apply_reflection_edge_north_east_frame(f16R)
  call apply_reflection_edge_north_east_frame(f17R)
  call apply_reflection_edge_north_east_frame(f18R)
  
  !apply reflection at north west 6   !xz
  !red fluid
  call apply_reflection_edge_north_west_frame(f00R)
  call apply_reflection_edge_north_west_frame(f01R)
  call apply_reflection_edge_north_west_frame(f02R)
  call apply_reflection_edge_north_west_frame(f03R)
  call apply_reflection_edge_north_west_frame(f04R)
  call apply_reflection_edge_north_west_frame(f05R)
  call apply_reflection_edge_north_west_frame(f06R)
  call apply_reflection_edge_north_west_frame(f07R)
  call apply_reflection_edge_north_west_frame(f08R)
  call apply_reflection_edge_north_west_frame(f09R)
  call apply_reflection_edge_north_west_frame(f10R)
  call apply_reflection_edge_north_west_frame(f11R)
  call apply_reflection_edge_north_west_frame(f12R)
  call apply_reflection_edge_north_west_frame(f13R)
  call apply_reflection_edge_north_west_frame(f14R)
  call apply_reflection_edge_north_west_frame(f15R)
  call apply_reflection_edge_north_west_frame(f16R)
  call apply_reflection_edge_north_west_frame(f17R)
  call apply_reflection_edge_north_west_frame(f18R)
  
  !apply reflection at south east 9  !xz
  !red fluid
  call apply_reflection_edge_south_east_frame(f00R)
  call apply_reflection_edge_south_east_frame(f01R)
  call apply_reflection_edge_south_east_frame(f02R)
  call apply_reflection_edge_south_east_frame(f03R)
  call apply_reflection_edge_south_east_frame(f04R)
  call apply_reflection_edge_south_east_frame(f05R)
  call apply_reflection_edge_south_east_frame(f06R)
  call apply_reflection_edge_south_east_frame(f07R)
  call apply_reflection_edge_south_east_frame(f08R)
  call apply_reflection_edge_south_east_frame(f09R)
  call apply_reflection_edge_south_east_frame(f10R)
  call apply_reflection_edge_south_east_frame(f11R)
  call apply_reflection_edge_south_east_frame(f12R)
  call apply_reflection_edge_south_east_frame(f13R)
  call apply_reflection_edge_south_east_frame(f14R)
  call apply_reflection_edge_south_east_frame(f15R)
  call apply_reflection_edge_south_east_frame(f16R)
  call apply_reflection_edge_south_east_frame(f17R)
  call apply_reflection_edge_south_east_frame(f18R)
  
  !apply reflection at south west 12  !xz
  !red fluid
  call apply_reflection_edge_south_west_frame(f00R)
  call apply_reflection_edge_south_west_frame(f01R)
  call apply_reflection_edge_south_west_frame(f02R)
  call apply_reflection_edge_south_west_frame(f03R)
  call apply_reflection_edge_south_west_frame(f04R)
  call apply_reflection_edge_south_west_frame(f05R)
  call apply_reflection_edge_south_west_frame(f06R)
  call apply_reflection_edge_south_west_frame(f07R)
  call apply_reflection_edge_south_west_frame(f08R)
  call apply_reflection_edge_south_west_frame(f09R)
  call apply_reflection_edge_south_west_frame(f10R)
  call apply_reflection_edge_south_west_frame(f11R)
  call apply_reflection_edge_south_west_frame(f12R)
  call apply_reflection_edge_south_west_frame(f13R)
  call apply_reflection_edge_south_west_frame(f14R)
  call apply_reflection_edge_south_west_frame(f15R)
  call apply_reflection_edge_south_west_frame(f16R)
  call apply_reflection_edge_south_west_frame(f17R)
  call apply_reflection_edge_south_west_frame(f18R)
  
  if(lsingle_fluid)return
  
  !sides 6
  
  !apply pbc at east 3 !x
  !blue fluid
  call apply_reflection_east(1-nbuff,ny+nbuff,1,nz,f00B)
  call apply_reflection_east(1-nbuff,ny+nbuff,1,nz,f01B)
  call apply_reflection_east(1-nbuff,ny+nbuff,1,nz,f02B)
  call apply_reflection_east(1-nbuff,ny+nbuff,1,nz,f03B)
  call apply_reflection_east(1-nbuff,ny+nbuff,1,nz,f04B)
  call apply_reflection_east(1-nbuff,ny+nbuff,1,nz,f05B)
  call apply_reflection_east(1-nbuff,ny+nbuff,1,nz,f06B)
  call apply_reflection_east(1-nbuff,ny+nbuff,1,nz,f07B)
  call apply_reflection_east(1-nbuff,ny+nbuff,1,nz,f08B)
  call apply_reflection_east(1-nbuff,ny+nbuff,1,nz,f09B)
  call apply_reflection_east(1-nbuff,ny+nbuff,1,nz,f10B)
  call apply_reflection_east(1-nbuff,ny+nbuff,1,nz,f11B)
  call apply_reflection_east(1-nbuff,ny+nbuff,1,nz,f12B)
  call apply_reflection_east(1-nbuff,ny+nbuff,1,nz,f13B)
  call apply_reflection_east(1-nbuff,ny+nbuff,1,nz,f14B)
  call apply_reflection_east(1-nbuff,ny+nbuff,1,nz,f15B)
  call apply_reflection_east(1-nbuff,ny+nbuff,1,nz,f16B)
  call apply_reflection_east(1-nbuff,ny+nbuff,1,nz,f17B)
  call apply_reflection_east(1-nbuff,ny+nbuff,1,nz,f18B)
  
  !apply pbc at west 4 !x
  !blue fluid
  call apply_reflection_west(1-nbuff,ny+nbuff,1,nz,f00B)
  call apply_reflection_west(1-nbuff,ny+nbuff,1,nz,f01B)
  call apply_reflection_west(1-nbuff,ny+nbuff,1,nz,f02B)
  call apply_reflection_west(1-nbuff,ny+nbuff,1,nz,f03B)
  call apply_reflection_west(1-nbuff,ny+nbuff,1,nz,f04B)
  call apply_reflection_west(1-nbuff,ny+nbuff,1,nz,f05B)
  call apply_reflection_west(1-nbuff,ny+nbuff,1,nz,f06B)
  call apply_reflection_west(1-nbuff,ny+nbuff,1,nz,f07B)
  call apply_reflection_west(1-nbuff,ny+nbuff,1,nz,f08B)
  call apply_reflection_west(1-nbuff,ny+nbuff,1,nz,f09B)
  call apply_reflection_west(1-nbuff,ny+nbuff,1,nz,f10B)
  call apply_reflection_west(1-nbuff,ny+nbuff,1,nz,f11B)
  call apply_reflection_west(1-nbuff,ny+nbuff,1,nz,f12B)
  call apply_reflection_west(1-nbuff,ny+nbuff,1,nz,f13B)
  call apply_reflection_west(1-nbuff,ny+nbuff,1,nz,f14B)
  call apply_reflection_west(1-nbuff,ny+nbuff,1,nz,f15B)
  call apply_reflection_west(1-nbuff,ny+nbuff,1,nz,f16B)
  call apply_reflection_west(1-nbuff,ny+nbuff,1,nz,f17B)
  call apply_reflection_west(1-nbuff,ny+nbuff,1,nz,f18B)
  
  !apply pbc at north 1 !z
  !blue fluid
  !frame y
  call apply_reflection_north(1,nx,1-nbuff,ny+nbuff,f00B)
  call apply_reflection_north(1,nx,1-nbuff,ny+nbuff,f01B)
  call apply_reflection_north(1,nx,1-nbuff,ny+nbuff,f02B)
  call apply_reflection_north(1,nx,1-nbuff,ny+nbuff,f03B)
  call apply_reflection_north(1,nx,1-nbuff,ny+nbuff,f04B)
  call apply_reflection_north(1,nx,1-nbuff,ny+nbuff,f05B)
  call apply_reflection_north(1,nx,1-nbuff,ny+nbuff,f06B)
  call apply_reflection_north(1,nx,1-nbuff,ny+nbuff,f07B)
  call apply_reflection_north(1,nx,1-nbuff,ny+nbuff,f08B)
  call apply_reflection_north(1,nx,1-nbuff,ny+nbuff,f09B)
  call apply_reflection_north(1,nx,1-nbuff,ny+nbuff,f10B)
  call apply_reflection_north(1,nx,1-nbuff,ny+nbuff,f11B)
  call apply_reflection_north(1,nx,1-nbuff,ny+nbuff,f12B)
  call apply_reflection_north(1,nx,1-nbuff,ny+nbuff,f13B)
  call apply_reflection_north(1,nx,1-nbuff,ny+nbuff,f14B)
  call apply_reflection_north(1,nx,1-nbuff,ny+nbuff,f15B)
  call apply_reflection_north(1,nx,1-nbuff,ny+nbuff,f16B)
  call apply_reflection_north(1,nx,1-nbuff,ny+nbuff,f17B)
  call apply_reflection_north(1,nx,1-nbuff,ny+nbuff,f18B)
  
  !apply pbc at south 2 !z
  !blue fluid
  call apply_reflection_south_frame_y(f00B)
  call apply_reflection_south_frame_y(f01B)
  call apply_reflection_south_frame_y(f02B)
  call apply_reflection_south_frame_y(f03B)
  call apply_reflection_south_frame_y(f04B)
  call apply_reflection_south_frame_y(f05B)
  call apply_reflection_south_frame_y(f06B)
  call apply_reflection_south_frame_y(f07B)
  call apply_reflection_south_frame_y(f08B)
  call apply_reflection_south_frame_y(f09B)
  call apply_reflection_south_frame_y(f10B)
  call apply_reflection_south_frame_y(f11B)
  call apply_reflection_south_frame_y(f12B)
  call apply_reflection_south_frame_y(f13B)
  call apply_reflection_south_frame_y(f14B)
  call apply_reflection_south_frame_y(f15B)
  call apply_reflection_south_frame_y(f16B)
  call apply_reflection_south_frame_y(f17B)
  call apply_reflection_south_frame_y(f18B)
  
  !edges 4 (12)
  
  !apply reflection at north east 3   !xz
  !blue fluid
  call apply_reflection_edge_north_east_frame(f00B)
  call apply_reflection_edge_north_east_frame(f01B)
  call apply_reflection_edge_north_east_frame(f02B)
  call apply_reflection_edge_north_east_frame(f03B)
  call apply_reflection_edge_north_east_frame(f04B)
  call apply_reflection_edge_north_east_frame(f05B)
  call apply_reflection_edge_north_east_frame(f06B)
  call apply_reflection_edge_north_east_frame(f07B)
  call apply_reflection_edge_north_east_frame(f08B)
  call apply_reflection_edge_north_east_frame(f09B)
  call apply_reflection_edge_north_east_frame(f10B)
  call apply_reflection_edge_north_east_frame(f11B)
  call apply_reflection_edge_north_east_frame(f12B)
  call apply_reflection_edge_north_east_frame(f13B)
  call apply_reflection_edge_north_east_frame(f14B)
  call apply_reflection_edge_north_east_frame(f15B)
  call apply_reflection_edge_north_east_frame(f16B)
  call apply_reflection_edge_north_east_frame(f17B)
  call apply_reflection_edge_north_east_frame(f18B)
  
  !apply reflection at north west 6   !xz
  !blue fluid
  call apply_reflection_edge_north_west_frame(f00B)
  call apply_reflection_edge_north_west_frame(f01B)
  call apply_reflection_edge_north_west_frame(f02B)
  call apply_reflection_edge_north_west_frame(f03B)
  call apply_reflection_edge_north_west_frame(f04B)
  call apply_reflection_edge_north_west_frame(f05B)
  call apply_reflection_edge_north_west_frame(f06B)
  call apply_reflection_edge_north_west_frame(f07B)
  call apply_reflection_edge_north_west_frame(f08B)
  call apply_reflection_edge_north_west_frame(f09B)
  call apply_reflection_edge_north_west_frame(f10B)
  call apply_reflection_edge_north_west_frame(f11B)
  call apply_reflection_edge_north_west_frame(f12B)
  call apply_reflection_edge_north_west_frame(f13B)
  call apply_reflection_edge_north_west_frame(f14B)
  call apply_reflection_edge_north_west_frame(f15B)
  call apply_reflection_edge_north_west_frame(f16B)
  call apply_reflection_edge_north_west_frame(f17B)
  call apply_reflection_edge_north_west_frame(f18B)
  
  !apply reflection at south east 9  !xz
  !blue fluid
  call apply_reflection_edge_south_east_frame(f00B)
  call apply_reflection_edge_south_east_frame(f01B)
  call apply_reflection_edge_south_east_frame(f02B)
  call apply_reflection_edge_south_east_frame(f03B)
  call apply_reflection_edge_south_east_frame(f04B)
  call apply_reflection_edge_south_east_frame(f05B)
  call apply_reflection_edge_south_east_frame(f06B)
  call apply_reflection_edge_south_east_frame(f07B)
  call apply_reflection_edge_south_east_frame(f08B)
  call apply_reflection_edge_south_east_frame(f09B)
  call apply_reflection_edge_south_east_frame(f10B)
  call apply_reflection_edge_south_east_frame(f11B)
  call apply_reflection_edge_south_east_frame(f12B)
  call apply_reflection_edge_south_east_frame(f13B)
  call apply_reflection_edge_south_east_frame(f14B)
  call apply_reflection_edge_south_east_frame(f15B)
  call apply_reflection_edge_south_east_frame(f16B)
  call apply_reflection_edge_south_east_frame(f17B)
  call apply_reflection_edge_south_east_frame(f18B)
  
  !apply reflection at south west 12  !xz
  !blue fluid
  call apply_reflection_edge_south_west_frame(f00B)
  call apply_reflection_edge_south_west_frame(f01B)
  call apply_reflection_edge_south_west_frame(f02B)
  call apply_reflection_edge_south_west_frame(f03B)
  call apply_reflection_edge_south_west_frame(f04B)
  call apply_reflection_edge_south_west_frame(f05B)
  call apply_reflection_edge_south_west_frame(f06B)
  call apply_reflection_edge_south_west_frame(f07B)
  call apply_reflection_edge_south_west_frame(f08B)
  call apply_reflection_edge_south_west_frame(f09B)
  call apply_reflection_edge_south_west_frame(f10B)
  call apply_reflection_edge_south_west_frame(f11B)
  call apply_reflection_edge_south_west_frame(f12B)
  call apply_reflection_edge_south_west_frame(f13B)
  call apply_reflection_edge_south_west_frame(f14B)
  call apply_reflection_edge_south_west_frame(f15B)
  call apply_reflection_edge_south_west_frame(f16B)
  call apply_reflection_edge_south_west_frame(f17B)
  call apply_reflection_edge_south_west_frame(f18B)
  
  return
  
 end subroutine apply_reflection_pops_along_xz
 
 subroutine apply_reflection_pops_along_xy
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     to fluid populations along x and y
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  !sides 4 (6)
  
  !apply pbc at east 3 !x
  !red fluid
  call apply_reflection_east(1,ny,1-nbuff,nz+nbuff,f00R)
  call apply_reflection_east(1,ny,1-nbuff,nz+nbuff,f01R)
  call apply_reflection_east(1,ny,1-nbuff,nz+nbuff,f02R)
  call apply_reflection_east(1,ny,1-nbuff,nz+nbuff,f03R)
  call apply_reflection_east(1,ny,1-nbuff,nz+nbuff,f04R)
  call apply_reflection_east(1,ny,1-nbuff,nz+nbuff,f05R)
  call apply_reflection_east(1,ny,1-nbuff,nz+nbuff,f06R)
  call apply_reflection_east(1,ny,1-nbuff,nz+nbuff,f07R)
  call apply_reflection_east(1,ny,1-nbuff,nz+nbuff,f08R)
  call apply_reflection_east(1,ny,1-nbuff,nz+nbuff,f09R)
  call apply_reflection_east(1,ny,1-nbuff,nz+nbuff,f10R)
  call apply_reflection_east(1,ny,1-nbuff,nz+nbuff,f11R)
  call apply_reflection_east(1,ny,1-nbuff,nz+nbuff,f12R)
  call apply_reflection_east(1,ny,1-nbuff,nz+nbuff,f13R)
  call apply_reflection_east(1,ny,1-nbuff,nz+nbuff,f14R)
  call apply_reflection_east(1,ny,1-nbuff,nz+nbuff,f15R)
  call apply_reflection_east(1,ny,1-nbuff,nz+nbuff,f16R)
  call apply_reflection_east(1,ny,1-nbuff,nz+nbuff,f17R)
  call apply_reflection_east(1,ny,1-nbuff,nz+nbuff,f18R)
  
  !apply pbc at west 4 !x
  !red fluid
  call apply_reflection_west(1,ny,1-nbuff,nz+nbuff,f00R)
  call apply_reflection_west(1,ny,1-nbuff,nz+nbuff,f01R)
  call apply_reflection_west(1,ny,1-nbuff,nz+nbuff,f02R)
  call apply_reflection_west(1,ny,1-nbuff,nz+nbuff,f03R)
  call apply_reflection_west(1,ny,1-nbuff,nz+nbuff,f04R)
  call apply_reflection_west(1,ny,1-nbuff,nz+nbuff,f05R)
  call apply_reflection_west(1,ny,1-nbuff,nz+nbuff,f06R)
  call apply_reflection_west(1,ny,1-nbuff,nz+nbuff,f07R)
  call apply_reflection_west(1,ny,1-nbuff,nz+nbuff,f08R)
  call apply_reflection_west(1,ny,1-nbuff,nz+nbuff,f09R)
  call apply_reflection_west(1,ny,1-nbuff,nz+nbuff,f10R)
  call apply_reflection_west(1,ny,1-nbuff,nz+nbuff,f11R)
  call apply_reflection_west(1,ny,1-nbuff,nz+nbuff,f12R)
  call apply_reflection_west(1,ny,1-nbuff,nz+nbuff,f13R)
  call apply_reflection_west(1,ny,1-nbuff,nz+nbuff,f14R)
  call apply_reflection_west(1,ny,1-nbuff,nz+nbuff,f15R)
  call apply_reflection_west(1,ny,1-nbuff,nz+nbuff,f16R)
  call apply_reflection_west(1,ny,1-nbuff,nz+nbuff,f17R)
  call apply_reflection_west(1,ny,1-nbuff,nz+nbuff,f18R)
  
  !apply pbc at front 5 !y
  !red fluid
  call apply_reflection_front_frame_z(f00R)
  call apply_reflection_front_frame_z(f01R)
  call apply_reflection_front_frame_z(f02R)
  call apply_reflection_front_frame_z(f03R)
  call apply_reflection_front_frame_z(f04R)
  call apply_reflection_front_frame_z(f05R)
  call apply_reflection_front_frame_z(f06R)
  call apply_reflection_front_frame_z(f07R)
  call apply_reflection_front_frame_z(f08R)
  call apply_reflection_front_frame_z(f09R)
  call apply_reflection_front_frame_z(f10R)
  call apply_reflection_front_frame_z(f11R)
  call apply_reflection_front_frame_z(f12R)
  call apply_reflection_front_frame_z(f13R)
  call apply_reflection_front_frame_z(f14R)
  call apply_reflection_front_frame_z(f15R)
  call apply_reflection_front_frame_z(f16R)
  call apply_reflection_front_frame_z(f17R)
  call apply_reflection_front_frame_z(f18R)
  
  !apply pbc at rear 6 !y
  !red fluid
  call apply_reflection_rear_frame_z(f00R)
  call apply_reflection_rear_frame_z(f01R)
  call apply_reflection_rear_frame_z(f02R)
  call apply_reflection_rear_frame_z(f03R)
  call apply_reflection_rear_frame_z(f04R)
  call apply_reflection_rear_frame_z(f05R)
  call apply_reflection_rear_frame_z(f06R)
  call apply_reflection_rear_frame_z(f07R)
  call apply_reflection_rear_frame_z(f08R)
  call apply_reflection_rear_frame_z(f09R)
  call apply_reflection_rear_frame_z(f10R)
  call apply_reflection_rear_frame_z(f11R)
  call apply_reflection_rear_frame_z(f12R)
  call apply_reflection_rear_frame_z(f13R)
  call apply_reflection_rear_frame_z(f14R)
  call apply_reflection_rear_frame_z(f15R)
  call apply_reflection_rear_frame_z(f16R)
  call apply_reflection_rear_frame_z(f17R)
  call apply_reflection_rear_frame_z(f18R)
  
  !edges 4 (12)
  
  !apply reflection at front east 1   !xy
  !red fluid
  call apply_reflection_edge_front_east_frame(f00R)
  call apply_reflection_edge_front_east_frame(f01R)
  call apply_reflection_edge_front_east_frame(f02R)
  call apply_reflection_edge_front_east_frame(f03R)
  call apply_reflection_edge_front_east_frame(f04R)
  call apply_reflection_edge_front_east_frame(f05R)
  call apply_reflection_edge_front_east_frame(f06R)
  call apply_reflection_edge_front_east_frame(f07R)
  call apply_reflection_edge_front_east_frame(f08R)
  call apply_reflection_edge_front_east_frame(f09R)
  call apply_reflection_edge_front_east_frame(f10R)
  call apply_reflection_edge_front_east_frame(f11R)
  call apply_reflection_edge_front_east_frame(f12R)
  call apply_reflection_edge_front_east_frame(f13R)
  call apply_reflection_edge_front_east_frame(f14R)
  call apply_reflection_edge_front_east_frame(f15R)
  call apply_reflection_edge_front_east_frame(f16R)
  call apply_reflection_edge_front_east_frame(f17R)
  call apply_reflection_edge_front_east_frame(f18R)
  
  !apply reflection at front west 2   !xy
  !red fluid
  call apply_reflection_edge_front_west_frame(f00R)
  call apply_reflection_edge_front_west_frame(f01R)
  call apply_reflection_edge_front_west_frame(f02R)
  call apply_reflection_edge_front_west_frame(f03R)
  call apply_reflection_edge_front_west_frame(f04R)
  call apply_reflection_edge_front_west_frame(f05R)
  call apply_reflection_edge_front_west_frame(f06R)
  call apply_reflection_edge_front_west_frame(f07R)
  call apply_reflection_edge_front_west_frame(f08R)
  call apply_reflection_edge_front_west_frame(f09R)
  call apply_reflection_edge_front_west_frame(f10R)
  call apply_reflection_edge_front_west_frame(f11R)
  call apply_reflection_edge_front_west_frame(f12R)
  call apply_reflection_edge_front_west_frame(f13R)
  call apply_reflection_edge_front_west_frame(f14R)
  call apply_reflection_edge_front_west_frame(f15R)
  call apply_reflection_edge_front_west_frame(f16R)
  call apply_reflection_edge_front_west_frame(f17R)
  call apply_reflection_edge_front_west_frame(f18R)
  
  !apply reflection at rear east 7   !xy
  !red fluid
  call apply_reflection_edge_rear_east_frame(f00R)
  call apply_reflection_edge_rear_east_frame(f01R)
  call apply_reflection_edge_rear_east_frame(f02R)
  call apply_reflection_edge_rear_east_frame(f03R)
  call apply_reflection_edge_rear_east_frame(f04R)
  call apply_reflection_edge_rear_east_frame(f05R)
  call apply_reflection_edge_rear_east_frame(f06R)
  call apply_reflection_edge_rear_east_frame(f07R)
  call apply_reflection_edge_rear_east_frame(f08R)
  call apply_reflection_edge_rear_east_frame(f09R)
  call apply_reflection_edge_rear_east_frame(f10R)
  call apply_reflection_edge_rear_east_frame(f11R)
  call apply_reflection_edge_rear_east_frame(f12R)
  call apply_reflection_edge_rear_east_frame(f13R)
  call apply_reflection_edge_rear_east_frame(f14R)
  call apply_reflection_edge_rear_east_frame(f15R)
  call apply_reflection_edge_rear_east_frame(f16R)
  call apply_reflection_edge_rear_east_frame(f17R)
  call apply_reflection_edge_rear_east_frame(f18R)
  
  !apply reflection at rear west 8   !xy
  !red fluid
  call apply_reflection_edge_rear_west_frame(f00R)
  call apply_reflection_edge_rear_west_frame(f01R)
  call apply_reflection_edge_rear_west_frame(f02R)
  call apply_reflection_edge_rear_west_frame(f03R)
  call apply_reflection_edge_rear_west_frame(f04R)
  call apply_reflection_edge_rear_west_frame(f05R)
  call apply_reflection_edge_rear_west_frame(f06R)
  call apply_reflection_edge_rear_west_frame(f07R)
  call apply_reflection_edge_rear_west_frame(f08R)
  call apply_reflection_edge_rear_west_frame(f09R)
  call apply_reflection_edge_rear_west_frame(f10R)
  call apply_reflection_edge_rear_west_frame(f11R)
  call apply_reflection_edge_rear_west_frame(f12R)
  call apply_reflection_edge_rear_west_frame(f13R)
  call apply_reflection_edge_rear_west_frame(f14R)
  call apply_reflection_edge_rear_west_frame(f15R)
  call apply_reflection_edge_rear_west_frame(f16R)
  call apply_reflection_edge_rear_west_frame(f17R)
  call apply_reflection_edge_rear_west_frame(f18R)
  
  
  if(lsingle_fluid)return
  
  !sides 4 (6)
  
  !apply pbc at east 3 !x
  !blue fluid
  call apply_reflection_east(1,ny,1-nbuff,nz+nbuff,f00B)
  call apply_reflection_east(1,ny,1-nbuff,nz+nbuff,f01B)
  call apply_reflection_east(1,ny,1-nbuff,nz+nbuff,f02B)
  call apply_reflection_east(1,ny,1-nbuff,nz+nbuff,f03B)
  call apply_reflection_east(1,ny,1-nbuff,nz+nbuff,f04B)
  call apply_reflection_east(1,ny,1-nbuff,nz+nbuff,f05B)
  call apply_reflection_east(1,ny,1-nbuff,nz+nbuff,f06B)
  call apply_reflection_east(1,ny,1-nbuff,nz+nbuff,f07B)
  call apply_reflection_east(1,ny,1-nbuff,nz+nbuff,f08B)
  call apply_reflection_east(1,ny,1-nbuff,nz+nbuff,f09B)
  call apply_reflection_east(1,ny,1-nbuff,nz+nbuff,f10B)
  call apply_reflection_east(1,ny,1-nbuff,nz+nbuff,f11B)
  call apply_reflection_east(1,ny,1-nbuff,nz+nbuff,f12B)
  call apply_reflection_east(1,ny,1-nbuff,nz+nbuff,f13B)
  call apply_reflection_east(1,ny,1-nbuff,nz+nbuff,f14B)
  call apply_reflection_east(1,ny,1-nbuff,nz+nbuff,f15B)
  call apply_reflection_east(1,ny,1-nbuff,nz+nbuff,f16B)
  call apply_reflection_east(1,ny,1-nbuff,nz+nbuff,f17B)
  call apply_reflection_east(1,ny,1-nbuff,nz+nbuff,f18B)
  
  !apply pbc at west 4 !x
  !blue fluid
  call apply_reflection_west(1,ny,1-nbuff,nz+nbuff,f00B)
  call apply_reflection_west(1,ny,1-nbuff,nz+nbuff,f01B)
  call apply_reflection_west(1,ny,1-nbuff,nz+nbuff,f02B)
  call apply_reflection_west(1,ny,1-nbuff,nz+nbuff,f03B)
  call apply_reflection_west(1,ny,1-nbuff,nz+nbuff,f04B)
  call apply_reflection_west(1,ny,1-nbuff,nz+nbuff,f05B)
  call apply_reflection_west(1,ny,1-nbuff,nz+nbuff,f06B)
  call apply_reflection_west(1,ny,1-nbuff,nz+nbuff,f07B)
  call apply_reflection_west(1,ny,1-nbuff,nz+nbuff,f08B)
  call apply_reflection_west(1,ny,1-nbuff,nz+nbuff,f09B)
  call apply_reflection_west(1,ny,1-nbuff,nz+nbuff,f10B)
  call apply_reflection_west(1,ny,1-nbuff,nz+nbuff,f11B)
  call apply_reflection_west(1,ny,1-nbuff,nz+nbuff,f12B)
  call apply_reflection_west(1,ny,1-nbuff,nz+nbuff,f13B)
  call apply_reflection_west(1,ny,1-nbuff,nz+nbuff,f14B)
  call apply_reflection_west(1,ny,1-nbuff,nz+nbuff,f15B)
  call apply_reflection_west(1,ny,1-nbuff,nz+nbuff,f16B)
  call apply_reflection_west(1,ny,1-nbuff,nz+nbuff,f17B)
  call apply_reflection_west(1,ny,1-nbuff,nz+nbuff,f18B)
  
  !apply pbc at front 5 !y
  !blue fluid
  call apply_reflection_front_frame_z(f00B)
  call apply_reflection_front_frame_z(f01B)
  call apply_reflection_front_frame_z(f02B)
  call apply_reflection_front_frame_z(f03B)
  call apply_reflection_front_frame_z(f04B)
  call apply_reflection_front_frame_z(f05B)
  call apply_reflection_front_frame_z(f06B)
  call apply_reflection_front_frame_z(f07B)
  call apply_reflection_front_frame_z(f08B)
  call apply_reflection_front_frame_z(f09B)
  call apply_reflection_front_frame_z(f10B)
  call apply_reflection_front_frame_z(f11B)
  call apply_reflection_front_frame_z(f12B)
  call apply_reflection_front_frame_z(f13B)
  call apply_reflection_front_frame_z(f14B)
  call apply_reflection_front_frame_z(f15B)
  call apply_reflection_front_frame_z(f16B)
  call apply_reflection_front_frame_z(f17B)
  call apply_reflection_front_frame_z(f18B)
  
  !apply pbc at rear 6 !y
  !blue fluid
  call apply_reflection_rear_frame_z(f00B)
  call apply_reflection_rear_frame_z(f01B)
  call apply_reflection_rear_frame_z(f02B)
  call apply_reflection_rear_frame_z(f03B)
  call apply_reflection_rear_frame_z(f04B)
  call apply_reflection_rear_frame_z(f05B)
  call apply_reflection_rear_frame_z(f06B)
  call apply_reflection_rear_frame_z(f07B)
  call apply_reflection_rear_frame_z(f08B)
  call apply_reflection_rear_frame_z(f09B)
  call apply_reflection_rear_frame_z(f10B)
  call apply_reflection_rear_frame_z(f11B)
  call apply_reflection_rear_frame_z(f12B)
  call apply_reflection_rear_frame_z(f13B)
  call apply_reflection_rear_frame_z(f14B)
  call apply_reflection_rear_frame_z(f15B)
  call apply_reflection_rear_frame_z(f16B)
  call apply_reflection_rear_frame_z(f17B)
  call apply_reflection_rear_frame_z(f18B)
  
  !edges 4 (12)
  
  !apply reflection at front east 1   !xy
  !blue fluid
  call apply_reflection_edge_front_east_frame(f00B)
  call apply_reflection_edge_front_east_frame(f01B)
  call apply_reflection_edge_front_east_frame(f02B)
  call apply_reflection_edge_front_east_frame(f03B)
  call apply_reflection_edge_front_east_frame(f04B)
  call apply_reflection_edge_front_east_frame(f05B)
  call apply_reflection_edge_front_east_frame(f06B)
  call apply_reflection_edge_front_east_frame(f07B)
  call apply_reflection_edge_front_east_frame(f08B)
  call apply_reflection_edge_front_east_frame(f09B)
  call apply_reflection_edge_front_east_frame(f10B)
  call apply_reflection_edge_front_east_frame(f11B)
  call apply_reflection_edge_front_east_frame(f12B)
  call apply_reflection_edge_front_east_frame(f13B)
  call apply_reflection_edge_front_east_frame(f14B)
  call apply_reflection_edge_front_east_frame(f15B)
  call apply_reflection_edge_front_east_frame(f16B)
  call apply_reflection_edge_front_east_frame(f17B)
  call apply_reflection_edge_front_east_frame(f18B)
  
  !apply reflection at front west 2   !xy
  !blue fluid
  call apply_reflection_edge_front_west_frame(f00B)
  call apply_reflection_edge_front_west_frame(f01B)
  call apply_reflection_edge_front_west_frame(f02B)
  call apply_reflection_edge_front_west_frame(f03B)
  call apply_reflection_edge_front_west_frame(f04B)
  call apply_reflection_edge_front_west_frame(f05B)
  call apply_reflection_edge_front_west_frame(f06B)
  call apply_reflection_edge_front_west_frame(f07B)
  call apply_reflection_edge_front_west_frame(f08B)
  call apply_reflection_edge_front_west_frame(f09B)
  call apply_reflection_edge_front_west_frame(f10B)
  call apply_reflection_edge_front_west_frame(f11B)
  call apply_reflection_edge_front_west_frame(f12B)
  call apply_reflection_edge_front_west_frame(f13B)
  call apply_reflection_edge_front_west_frame(f14B)
  call apply_reflection_edge_front_west_frame(f15B)
  call apply_reflection_edge_front_west_frame(f16B)
  call apply_reflection_edge_front_west_frame(f17B)
  call apply_reflection_edge_front_west_frame(f18B)
  
  !apply reflection at rear east 7   !xy
  !blue fluid
  call apply_reflection_edge_rear_east_frame(f00B)
  call apply_reflection_edge_rear_east_frame(f01B)
  call apply_reflection_edge_rear_east_frame(f02B)
  call apply_reflection_edge_rear_east_frame(f03B)
  call apply_reflection_edge_rear_east_frame(f04B)
  call apply_reflection_edge_rear_east_frame(f05B)
  call apply_reflection_edge_rear_east_frame(f06B)
  call apply_reflection_edge_rear_east_frame(f07B)
  call apply_reflection_edge_rear_east_frame(f08B)
  call apply_reflection_edge_rear_east_frame(f09B)
  call apply_reflection_edge_rear_east_frame(f10B)
  call apply_reflection_edge_rear_east_frame(f11B)
  call apply_reflection_edge_rear_east_frame(f12B)
  call apply_reflection_edge_rear_east_frame(f13B)
  call apply_reflection_edge_rear_east_frame(f14B)
  call apply_reflection_edge_rear_east_frame(f15B)
  call apply_reflection_edge_rear_east_frame(f16B)
  call apply_reflection_edge_rear_east_frame(f17B)
  call apply_reflection_edge_rear_east_frame(f18B)
  
  !apply reflection at rear west 8   !xy
  !blue fluid
  call apply_reflection_edge_rear_west_frame(f00B)
  call apply_reflection_edge_rear_west_frame(f01B)
  call apply_reflection_edge_rear_west_frame(f02B)
  call apply_reflection_edge_rear_west_frame(f03B)
  call apply_reflection_edge_rear_west_frame(f04B)
  call apply_reflection_edge_rear_west_frame(f05B)
  call apply_reflection_edge_rear_west_frame(f06B)
  call apply_reflection_edge_rear_west_frame(f07B)
  call apply_reflection_edge_rear_west_frame(f08B)
  call apply_reflection_edge_rear_west_frame(f09B)
  call apply_reflection_edge_rear_west_frame(f10B)
  call apply_reflection_edge_rear_west_frame(f11B)
  call apply_reflection_edge_rear_west_frame(f12B)
  call apply_reflection_edge_rear_west_frame(f13B)
  call apply_reflection_edge_rear_west_frame(f14B)
  call apply_reflection_edge_rear_west_frame(f15B)
  call apply_reflection_edge_rear_west_frame(f16B)
  call apply_reflection_edge_rear_west_frame(f17B)
  call apply_reflection_edge_rear_west_frame(f18B)
  
  return
  
 end subroutine apply_reflection_pops_along_xy
 
 subroutine apply_reflection_pops_along_z
 
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
  
  !apply pbc at north 1 !z
  !red fluid
  !frame
  call apply_reflection_north(1-nbuff,nx+nbuff,1-nbuff,ny+nbuff,f00R)
  call apply_reflection_north(1-nbuff,nx+nbuff,1-nbuff,ny+nbuff,f01R)
  call apply_reflection_north(1-nbuff,nx+nbuff,1-nbuff,ny+nbuff,f02R)
  call apply_reflection_north(1-nbuff,nx+nbuff,1-nbuff,ny+nbuff,f03R)
  call apply_reflection_north(1-nbuff,nx+nbuff,1-nbuff,ny+nbuff,f04R)
  call apply_reflection_north(1-nbuff,nx+nbuff,1-nbuff,ny+nbuff,f05R)
  call apply_reflection_north(1-nbuff,nx+nbuff,1-nbuff,ny+nbuff,f06R)
  call apply_reflection_north(1-nbuff,nx+nbuff,1-nbuff,ny+nbuff,f07R)
  call apply_reflection_north(1-nbuff,nx+nbuff,1-nbuff,ny+nbuff,f08R)
  call apply_reflection_north(1-nbuff,nx+nbuff,1-nbuff,ny+nbuff,f09R)
  call apply_reflection_north(1-nbuff,nx+nbuff,1-nbuff,ny+nbuff,f10R)
  call apply_reflection_north(1-nbuff,nx+nbuff,1-nbuff,ny+nbuff,f11R)
  call apply_reflection_north(1-nbuff,nx+nbuff,1-nbuff,ny+nbuff,f12R)
  call apply_reflection_north(1-nbuff,nx+nbuff,1-nbuff,ny+nbuff,f13R)
  call apply_reflection_north(1-nbuff,nx+nbuff,1-nbuff,ny+nbuff,f14R)
  call apply_reflection_north(1-nbuff,nx+nbuff,1-nbuff,ny+nbuff,f15R)
  call apply_reflection_north(1-nbuff,nx+nbuff,1-nbuff,ny+nbuff,f16R)
  call apply_reflection_north(1-nbuff,nx+nbuff,1-nbuff,ny+nbuff,f17R)
  call apply_reflection_north(1-nbuff,nx+nbuff,1-nbuff,ny+nbuff,f18R)
  
  !apply pbc at south 2 !z
  !red fluid
  call apply_reflection_south_frame(f00R)
  call apply_reflection_south_frame(f01R)
  call apply_reflection_south_frame(f02R)
  call apply_reflection_south_frame(f03R)
  call apply_reflection_south_frame(f04R)
  call apply_reflection_south_frame(f05R)
  call apply_reflection_south_frame(f06R)
  call apply_reflection_south_frame(f07R)
  call apply_reflection_south_frame(f08R)
  call apply_reflection_south_frame(f09R)
  call apply_reflection_south_frame(f10R)
  call apply_reflection_south_frame(f11R)
  call apply_reflection_south_frame(f12R)
  call apply_reflection_south_frame(f13R)
  call apply_reflection_south_frame(f14R)
  call apply_reflection_south_frame(f15R)
  call apply_reflection_south_frame(f16R)
  call apply_reflection_south_frame(f17R)
  call apply_reflection_south_frame(f18R)
  
  if(lsingle_fluid)return
  
  !sides 6
  
  !apply pbc at north 1 !z
  !blue fluid
  !frame
  call apply_reflection_north(1-nbuff,nx+nbuff,1-nbuff,ny+nbuff,f00B)
  call apply_reflection_north(1-nbuff,nx+nbuff,1-nbuff,ny+nbuff,f01B)
  call apply_reflection_north(1-nbuff,nx+nbuff,1-nbuff,ny+nbuff,f02B)
  call apply_reflection_north(1-nbuff,nx+nbuff,1-nbuff,ny+nbuff,f03B)
  call apply_reflection_north(1-nbuff,nx+nbuff,1-nbuff,ny+nbuff,f04B)
  call apply_reflection_north(1-nbuff,nx+nbuff,1-nbuff,ny+nbuff,f05B)
  call apply_reflection_north(1-nbuff,nx+nbuff,1-nbuff,ny+nbuff,f06B)
  call apply_reflection_north(1-nbuff,nx+nbuff,1-nbuff,ny+nbuff,f07B)
  call apply_reflection_north(1-nbuff,nx+nbuff,1-nbuff,ny+nbuff,f08B)
  call apply_reflection_north(1-nbuff,nx+nbuff,1-nbuff,ny+nbuff,f09B)
  call apply_reflection_north(1-nbuff,nx+nbuff,1-nbuff,ny+nbuff,f10B)
  call apply_reflection_north(1-nbuff,nx+nbuff,1-nbuff,ny+nbuff,f11B)
  call apply_reflection_north(1-nbuff,nx+nbuff,1-nbuff,ny+nbuff,f12B)
  call apply_reflection_north(1-nbuff,nx+nbuff,1-nbuff,ny+nbuff,f13B)
  call apply_reflection_north(1-nbuff,nx+nbuff,1-nbuff,ny+nbuff,f14B)
  call apply_reflection_north(1-nbuff,nx+nbuff,1-nbuff,ny+nbuff,f15B)
  call apply_reflection_north(1-nbuff,nx+nbuff,1-nbuff,ny+nbuff,f16B)
  call apply_reflection_north(1-nbuff,nx+nbuff,1-nbuff,ny+nbuff,f17B)
  call apply_reflection_north(1-nbuff,nx+nbuff,1-nbuff,ny+nbuff,f18B)
  
  !apply pbc at south 2 !z
  !blue fluid
  call apply_reflection_south_frame(f00B)
  call apply_reflection_south_frame(f01B)
  call apply_reflection_south_frame(f02B)
  call apply_reflection_south_frame(f03B)
  call apply_reflection_south_frame(f04B)
  call apply_reflection_south_frame(f05B)
  call apply_reflection_south_frame(f06B)
  call apply_reflection_south_frame(f07B)
  call apply_reflection_south_frame(f08B)
  call apply_reflection_south_frame(f09B)
  call apply_reflection_south_frame(f10B)
  call apply_reflection_south_frame(f11B)
  call apply_reflection_south_frame(f12B)
  call apply_reflection_south_frame(f13B)
  call apply_reflection_south_frame(f14B)
  call apply_reflection_south_frame(f15B)
  call apply_reflection_south_frame(f16B)
  call apply_reflection_south_frame(f17B)
  call apply_reflection_south_frame(f18B)
  
  return
  
 end subroutine apply_reflection_pops_along_z
 
 subroutine apply_reflection_pops_along_y
 
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
  
  !apply pbc at front 5 !y
  !red fluid
  call apply_reflection_front_frame(f00R)
  call apply_reflection_front_frame(f01R)
  call apply_reflection_front_frame(f02R)
  call apply_reflection_front_frame(f03R)
  call apply_reflection_front_frame(f04R)
  call apply_reflection_front_frame(f05R)
  call apply_reflection_front_frame(f06R)
  call apply_reflection_front_frame(f07R)
  call apply_reflection_front_frame(f08R)
  call apply_reflection_front_frame(f09R)
  call apply_reflection_front_frame(f10R)
  call apply_reflection_front_frame(f11R)
  call apply_reflection_front_frame(f12R)
  call apply_reflection_front_frame(f13R)
  call apply_reflection_front_frame(f14R)
  call apply_reflection_front_frame(f15R)
  call apply_reflection_front_frame(f16R)
  call apply_reflection_front_frame(f17R)
  call apply_reflection_front_frame(f18R)
  
  !apply pbc at rear 6 !y
  !red fluid
  call apply_reflection_rear_frame(f00R)
  call apply_reflection_rear_frame(f01R)
  call apply_reflection_rear_frame(f02R)
  call apply_reflection_rear_frame(f03R)
  call apply_reflection_rear_frame(f04R)
  call apply_reflection_rear_frame(f05R)
  call apply_reflection_rear_frame(f06R)
  call apply_reflection_rear_frame(f07R)
  call apply_reflection_rear_frame(f08R)
  call apply_reflection_rear_frame(f09R)
  call apply_reflection_rear_frame(f10R)
  call apply_reflection_rear_frame(f11R)
  call apply_reflection_rear_frame(f12R)
  call apply_reflection_rear_frame(f13R)
  call apply_reflection_rear_frame(f14R)
  call apply_reflection_rear_frame(f15R)
  call apply_reflection_rear_frame(f16R)
  call apply_reflection_rear_frame(f17R)
  call apply_reflection_rear_frame(f18R)
  
  if(lsingle_fluid)return
  
  !sides 6
  
  !apply pbc at front 5 !y
  !blue fluid
  call apply_reflection_front_frame(f00B)
  call apply_reflection_front_frame(f01B)
  call apply_reflection_front_frame(f02B)
  call apply_reflection_front_frame(f03B)
  call apply_reflection_front_frame(f04B)
  call apply_reflection_front_frame(f05B)
  call apply_reflection_front_frame(f06B)
  call apply_reflection_front_frame(f07B)
  call apply_reflection_front_frame(f08B)
  call apply_reflection_front_frame(f09B)
  call apply_reflection_front_frame(f10B)
  call apply_reflection_front_frame(f11B)
  call apply_reflection_front_frame(f12B)
  call apply_reflection_front_frame(f13B)
  call apply_reflection_front_frame(f14B)
  call apply_reflection_front_frame(f15B)
  call apply_reflection_front_frame(f16B)
  call apply_reflection_front_frame(f17B)
  call apply_reflection_front_frame(f18B)
  
  !apply pbc at rear 6 !y
  !blue fluid
  call apply_reflection_rear_frame(f00B)
  call apply_reflection_rear_frame(f01B)
  call apply_reflection_rear_frame(f02B)
  call apply_reflection_rear_frame(f03B)
  call apply_reflection_rear_frame(f04B)
  call apply_reflection_rear_frame(f05B)
  call apply_reflection_rear_frame(f06B)
  call apply_reflection_rear_frame(f07B)
  call apply_reflection_rear_frame(f08B)
  call apply_reflection_rear_frame(f09B)
  call apply_reflection_rear_frame(f10B)
  call apply_reflection_rear_frame(f11B)
  call apply_reflection_rear_frame(f12B)
  call apply_reflection_rear_frame(f13B)
  call apply_reflection_rear_frame(f14B)
  call apply_reflection_rear_frame(f15B)
  call apply_reflection_rear_frame(f16B)
  call apply_reflection_rear_frame(f17B)
  call apply_reflection_rear_frame(f18B)
  
  return
  
 end subroutine apply_reflection_pops_along_y
 
 subroutine apply_reflection_pops_along_x
 
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
  
  !apply pbc at east 3 !x
  !red fluid
  call apply_reflection_east(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f00R)
  call apply_reflection_east(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f01R)
  call apply_reflection_east(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f02R)
  call apply_reflection_east(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f03R)
  call apply_reflection_east(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f04R)
  call apply_reflection_east(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f05R)
  call apply_reflection_east(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f06R)
  call apply_reflection_east(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f07R)
  call apply_reflection_east(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f08R)
  call apply_reflection_east(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f09R)
  call apply_reflection_east(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f10R)
  call apply_reflection_east(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f11R)
  call apply_reflection_east(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f12R)
  call apply_reflection_east(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f13R)
  call apply_reflection_east(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f14R)
  call apply_reflection_east(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f15R)
  call apply_reflection_east(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f16R)
  call apply_reflection_east(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f17R)
  call apply_reflection_east(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f18R)
  
  !apply pbc at west 4 !x
  !red fluid
  call apply_reflection_west(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f00R)
  call apply_reflection_west(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f01R)
  call apply_reflection_west(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f02R)
  call apply_reflection_west(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f03R)
  call apply_reflection_west(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f04R)
  call apply_reflection_west(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f05R)
  call apply_reflection_west(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f06R)
  call apply_reflection_west(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f07R)
  call apply_reflection_west(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f08R)
  call apply_reflection_west(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f09R)
  call apply_reflection_west(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f10R)
  call apply_reflection_west(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f11R)
  call apply_reflection_west(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f12R)
  call apply_reflection_west(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f13R)
  call apply_reflection_west(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f14R)
  call apply_reflection_west(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f15R)
  call apply_reflection_west(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f16R)
  call apply_reflection_west(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f17R)
  call apply_reflection_west(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f18R)
  
  if(lsingle_fluid)return
  
  !sides 6
  
  !apply pbc at east 3 !x
  !blue fluid
  call apply_reflection_east(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f00B)
  call apply_reflection_east(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f01B)
  call apply_reflection_east(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f02B)
  call apply_reflection_east(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f03B)
  call apply_reflection_east(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f04B)
  call apply_reflection_east(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f05B)
  call apply_reflection_east(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f06B)
  call apply_reflection_east(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f07B)
  call apply_reflection_east(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f08B)
  call apply_reflection_east(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f09B)
  call apply_reflection_east(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f10B)
  call apply_reflection_east(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f11B)
  call apply_reflection_east(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f12B)
  call apply_reflection_east(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f13B)
  call apply_reflection_east(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f14B)
  call apply_reflection_east(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f15B)
  call apply_reflection_east(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f16B)
  call apply_reflection_east(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f17B)
  call apply_reflection_east(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f18B)
  
  !apply pbc at west 4 !x
  !blue fluid
  call apply_reflection_west(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f00B)
  call apply_reflection_west(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f01B)
  call apply_reflection_west(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f02B)
  call apply_reflection_west(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f03B)
  call apply_reflection_west(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f04B)
  call apply_reflection_west(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f05B)
  call apply_reflection_west(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f06B)
  call apply_reflection_west(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f07B)
  call apply_reflection_west(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f08B)
  call apply_reflection_west(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f09B)
  call apply_reflection_west(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f10B)
  call apply_reflection_west(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f11B)
  call apply_reflection_west(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f12B)
  call apply_reflection_west(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f13B)
  call apply_reflection_west(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f14B)
  call apply_reflection_west(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f15B)
  call apply_reflection_west(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f16B)
  call apply_reflection_west(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f17B)
  call apply_reflection_west(1-nbuff,ny+nbuff,1-nbuff,nz+nbuff,f18B)
  
  return
  
 end subroutine apply_reflection_pops_along_x
 
 subroutine apply_reflection_east(sy,ey,sz,ez,dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection
!     at the east side along the space defined in input sy,ey,sz,ez
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: sy,ey,sz,ez
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk
  
  forall(kk=1:nbuff,j=sy:ey,k=sz:ez)
    dtemp(nx+kk,j,k)=dtemp(nx-kk+1,j,k)
  end forall
  
  return
  
 end subroutine apply_reflection_east
 
 subroutine apply_reflection_west(sy,ey,sz,ez,dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection
!     at the west side along the space defined in input sy,ey,sz,ez
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: sy,ey,sz,ez
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk

  forall(kk=1:nbuff,j=sy:ey,k=sz:ez)
    dtemp(1-kk,j,k)=dtemp(kk,j,k)
  end forall
  
  return
  
 end subroutine apply_reflection_west
 
 subroutine apply_reflection_north(sx,ex,sy,ey,dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection
!     at the north side along the space defined in input sx,ex,sy,ey
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  integer, intent(in) :: sx,ex,sy,ey
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk
  
  forall(i=sx:ex,j=sy:ey,kk=1:nbuff)
    dtemp(i,j,nz+kk)= dtemp(i,j,nz-kk+1)
  end forall
  
  return
  
 end subroutine apply_reflection_north
 
 subroutine apply_reflection_south_frame(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     at the south side alongside with its frame
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk
  
  forall(i=1-nbuff:nx+nbuff,j=1-nbuff:ny+nbuff,kk=1:nbuff)
    dtemp(i,j,1-kk)=dtemp(i,j,kk)
  end forall
  
  return
  
 end subroutine apply_reflection_south_frame
 
 subroutine apply_reflection_south_frame_x(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     at the south side alongside with its frame only along x
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk
  
  forall(i=1-nbuff:nx+nbuff,j=1:ny,kk=1:nbuff)
    dtemp(i,j,1-kk)=dtemp(i,j,kk)
  end forall
  
  return
  
 end subroutine apply_reflection_south_frame_x
 
 subroutine apply_reflection_south_frame_y(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     at the south side alongside with its frame only along y
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk
  
  forall(i=1:nx,j=1-nbuff:ny+nbuff,kk=1:nbuff)
    dtemp(i,j,1-kk)=dtemp(i,j,kk)
  end forall
  
  return
  
 end subroutine apply_reflection_south_frame_y
 
 subroutine apply_reflection_front_frame(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection
!     at the front side alongside with its frame
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk
  
  forall(i=1-nbuff:nx+nbuff,kk=1:nbuff,k=1-nbuff:nz+nbuff)
    dtemp(i,1-kk,k)= dtemp(i,kk,k)
  end forall

  return
  
 end subroutine apply_reflection_front_frame
 
 subroutine apply_reflection_front_frame_x(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection
!     at the front side alongside with its frame only along x
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk
  
  forall(i=1-nbuff:nx+nbuff,kk=1:nbuff,k=1:nz)
    dtemp(i,1-kk,k)= dtemp(i,kk,k)
  end forall

  return
  
 end subroutine apply_reflection_front_frame_x
 
 subroutine apply_reflection_front_frame_z(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection
!     at the front side alongside with its frame only along z
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk
  
  forall(i=1:nx,kk=1:nbuff,k=1-nbuff:nz+nbuff)
    dtemp(i,1-kk,k)= dtemp(i,kk,k)
  end forall

  return
  
 end subroutine apply_reflection_front_frame_z
 
 subroutine apply_reflection_rear_frame(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     at the rear side alongside with its frame
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk
  
  forall(i=1-nbuff:nx+nbuff,kk=1:nbuff,k=1-nbuff:nz+nbuff)
    dtemp(i,ny+kk,k)= dtemp(i,ny-kk+1,k)
  end forall

  return
  
 end subroutine apply_reflection_rear_frame
 
 subroutine apply_reflection_rear_frame_x(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     at the rear side alongside with its frame only along x
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk
  
  forall(i=1-nbuff:nx+nbuff,kk=1:nbuff,k=1:nz)
    dtemp(i,ny+kk,k)= dtemp(i,ny-kk+1,k)
  end forall

  return
  
 end subroutine apply_reflection_rear_frame_x
 
 subroutine apply_reflection_rear_frame_z(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     at the rear side alongside with its frame only along z
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk
  
  forall(i=1:nx,kk=1:nbuff,k=1-nbuff:nz+nbuff)
    dtemp(i,ny+kk,k)= dtemp(i,ny-kk+1,k)
  end forall

  return
  
 end subroutine apply_reflection_rear_frame_z
 
 subroutine apply_reflection_south(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     at the south side
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
    dtemp(i,j,1-kk)=dtemp(i,j,kk)
  end forall
  
  return
  
 end subroutine apply_reflection_south
 
 subroutine apply_reflection_front(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection
!     at the front side
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
    dtemp(i,1-kk,k)= dtemp(i,kk,k)
  end forall

  return
  
 end subroutine apply_reflection_front
 
 subroutine apply_reflection_rear(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     at the rear side
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
    dtemp(i,ny+kk,k)= dtemp(i,ny-kk+1,k)
  end forall

  return
  
 end subroutine apply_reflection_rear
 
 subroutine apply_reflection_edge_front_east(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     at the front east edge
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
    dtemp(nx+kkk,1-kk,k)=dtemp(nx-kkk+1,kk,k)
  end forall

  return
  
 end subroutine apply_reflection_edge_front_east
 
 subroutine apply_reflection_edge_front_east_frame(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     at the front east edge with its frame
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk,kkk
  
  forall(kkk=1:nbuff,kk=1:nbuff,k=1-nbuff:nz+nbuff)
    dtemp(nx+kkk,1-kk,k)=dtemp(nx-kkk+1,kk,k)
  end forall

  return
  
 end subroutine apply_reflection_edge_front_east_frame
 
 subroutine apply_reflection_edge_front_west(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     at the front west edge
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
    dtemp(1-kkk,1-kk,k)=dtemp(kkk,kk,k)
  end forall

  return
  
 end subroutine apply_reflection_edge_front_west
 
 subroutine apply_reflection_edge_front_west_frame(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     at the front west edge with its frame
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk,kkk
  
  forall(kkk=1:nbuff,kk=1:nbuff,k=1-nbuff:nz+nbuff)
    dtemp(1-kkk,1-kk,k)=dtemp(kkk,kk,k)
  end forall

  return
  
 end subroutine apply_reflection_edge_front_west_frame
 
 subroutine apply_reflection_edge_north_east(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     at the north east edge
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
    dtemp(nx+kkk,j,nz+kk)=dtemp(nx-kkk+1,j,nz-kk+1)
  end forall

  return
  
 end subroutine apply_reflection_edge_north_east
 
 subroutine apply_reflection_edge_north_east_frame(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     at the north east edge with its frame
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk,kkk
  
  forall(kkk=1:nbuff,j=1-nbuff:ny+nbuff,kk=1:nbuff)
    dtemp(nx+kkk,j,nz+kk)=dtemp(nx-kkk+1,j,nz-kk+1)
  end forall

  return
  
 end subroutine apply_reflection_edge_north_east_frame
 
 subroutine apply_reflection_edge_north_front(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     at the north front edge
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
    dtemp(i,1-kkk,nz+kk)=dtemp(i,kkk,nz-kk+1)
  end forall

  return
  
 end subroutine apply_reflection_edge_north_front
 
 subroutine apply_reflection_edge_north_front_frame(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     at the north front edge with its frame
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk,kkk
  
  forall(i=1-nbuff:nx+nbuff,kkk=1:nbuff,kk=1:nbuff)
    dtemp(i,1-kkk,nz+kk)=dtemp(i,kkk,nz-kk+1)
  end forall

  return
  
 end subroutine apply_reflection_edge_north_front_frame
 
 subroutine apply_reflection_edge_north_rear(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     at the north rear edge
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
    dtemp(i,ny+kkk,nz+kk)= dtemp(i,ny-kkk+1,nz-kk+1)
  end forall

  return
  
 end subroutine apply_reflection_edge_north_rear
 
 subroutine apply_reflection_edge_north_rear_frame(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     at the north rear edge with its frame
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk,kkk
  
  forall(i=1-nbuff:nx+nbuff,kkk=1:nbuff,kk=1:nbuff)
    dtemp(i,ny+kkk,nz+kk)= dtemp(i,ny-kkk+1,nz-kk+1)
  end forall

  return
  
 end subroutine apply_reflection_edge_north_rear_frame
 
 subroutine apply_reflection_edge_north_west(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     at the north west edge
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
    dtemp(1-kkk,j,nz+kk)= dtemp(kkk,j,nz-kk+1)
  end forall
  
  return
  
 end subroutine apply_reflection_edge_north_west
 
 subroutine apply_reflection_edge_north_west_frame(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     at the north west edge with its frame
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk,kkk
  
  forall(kkk=1:nbuff,j=1-nbuff:ny+nbuff,kk=1:nbuff)
    dtemp(1-kkk,j,nz+kk)= dtemp(kkk,j,nz-kk+1)
  end forall
  
  return
  
 end subroutine apply_reflection_edge_north_west_frame
 
 subroutine apply_reflection_edge_rear_east(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     at the rear east edge
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
    dtemp(nx+kkk,ny+kk,k)= dtemp(nx-kkk+1,ny-kk+1,k)
  end forall
  
  return
  
 end subroutine apply_reflection_edge_rear_east
 
 subroutine apply_reflection_edge_rear_east_frame(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     at the rear east edge with its frame
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk,kkk
  
  forall(kkk=1:nbuff,kk=1:nbuff,k=1-nbuff:nz+nbuff)
    dtemp(nx+kkk,ny+kk,k)= dtemp(nx-kkk+1,ny-kk+1,k)
  end forall
  
  return
  
 end subroutine apply_reflection_edge_rear_east_frame
 
 subroutine apply_reflection_edge_rear_west(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     at the rear west edge
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
    dtemp(1-kkk,ny+kk,k)=dtemp(kkk,ny-kk+1,k)
  end forall
  
  return
  
 end subroutine apply_reflection_edge_rear_west
 
 subroutine apply_reflection_edge_rear_west_frame(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     at the rear west edge with its frame
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk,kkk
  
  forall(kkk=1:nbuff,kk=1:nbuff,k=1-nbuff:nz+nbuff)
    dtemp(1-kkk,ny+kk,k)=dtemp(kkk,ny-kk+1,k)
  end forall
  
  return
  
 end subroutine apply_reflection_edge_rear_west_frame
 
 subroutine apply_reflection_edge_south_east(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     at the south east edge
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
    dtemp(nx+kkk,j,1-kk)=dtemp(nx-kkk+1,j,kk)
  end forall
  
  return
  
 end subroutine apply_reflection_edge_south_east
 
 subroutine apply_reflection_edge_south_east_frame(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     at the south east edge with its frame
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk,kkk
  
  forall(kkk=1:nbuff,j=1-nbuff:ny+nbuff,kk=1:nbuff)
    dtemp(nx+kkk,j,1-kk)=dtemp(nx-kkk+1,j,kk)
  end forall
  
  return
  
 end subroutine apply_reflection_edge_south_east_frame
 
 subroutine apply_reflection_edge_south_front(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     at the south east edge
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
    dtemp(i,1-kkk,1-kk)=dtemp(i,kkk,kk)
  end forall
  
  return
  
 end subroutine apply_reflection_edge_south_front
 
 subroutine apply_reflection_edge_south_front_frame(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     at the south east edge with its frame
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk,kkk
  
  forall(i=1-nbuff:nx+nbuff,kkk=1:nbuff,kk=1:nbuff)
    dtemp(i,1-kkk,1-kk)=dtemp(i,kkk,kk)
  end forall
  
  return
  
 end subroutine apply_reflection_edge_south_front_frame
 
 subroutine apply_reflection_edge_south_rear(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     at the south rear edge
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
    dtemp(i,ny+kkk,1-kk)=dtemp(i,ny-kkk+1,kk)
  end forall
  
  return
  
 end subroutine apply_reflection_edge_south_rear
 
 subroutine apply_reflection_edge_south_rear_frame(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     at the south rear edge with its frame
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk,kkk
  
  forall(i=1-nbuff:nx+nbuff,kkk=1:nbuff,kk=1:nbuff)
    dtemp(i,ny+kkk,1-kk)=dtemp(i,ny-kkk+1,kk)
  end forall
  
  return
  
 end subroutine apply_reflection_edge_south_rear_frame
 
 subroutine apply_reflection_edge_south_west(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     at the south rear edge
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
    dtemp(1-kkk,j,1-kk)=dtemp(kkk,j,kk)
  end forall
  
  return
  
 end subroutine apply_reflection_edge_south_west
 
 subroutine apply_reflection_edge_south_west_frame(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     at the south rear edge with its frame
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i,j,k,kk,kkk
  
  forall(kkk=1:nbuff,j=1-nbuff:ny+nbuff,kk=1:nbuff)
    dtemp(1-kkk,j,1-kk)=dtemp(kkk,j,kk)
  end forall
  
  return
  
 end subroutine apply_reflection_edge_south_west_frame
 
 subroutine apply_reflection_corner_north_east_front(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     at the north east front corner
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
    dtemp(nx+kk,1-kkk,nz+kkkk)=dtemp(nx-kk+1,kkk,nz-kkkk+1)
  end forall
  
  return
  
 end subroutine apply_reflection_corner_north_east_front
 
 subroutine apply_reflection_corner_north_east_rear(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     at the north east rear corner
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
    dtemp(nx+kk,ny+kkk,nz+kkkk)=dtemp(nx-kk+1,ny-kkk+1,nz-kkkk+1)
  end forall
  
  return
  
 end subroutine apply_reflection_corner_north_east_rear
 
 subroutine apply_reflection_corner_north_west_rear(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     at the north west rear corner
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
    dtemp(1-kk,ny+kkk,nz+kkkk)=dtemp(kk,ny-kkk+1,nz-kkkk+1)
  end forall
  
  return
  
 end subroutine apply_reflection_corner_north_west_rear
 
 subroutine apply_reflection_corner_north_west_front(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     at the north west front corner
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
    dtemp(1-kk,1-kkk,nz+kkkk)=dtemp(kk,kkk,nz-kkkk+1)
  end forall
  
  return
  
 end subroutine apply_reflection_corner_north_west_front
 
 subroutine apply_reflection_corner_south_east_front(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     at the south east front corner
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
    dtemp(nx+kk,1-kkk,1-kkkk)=dtemp(nx-kk+1,kkk,kkkk)
  end forall
  
  return
  
 end subroutine apply_reflection_corner_south_east_front
 
 subroutine apply_reflection_corner_south_west_front(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     at the south west front corner
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
    dtemp(1-kk,1-kkk,1-kkkk)=dtemp(kk,kkk,kkkk)
  end forall
  
  return
  
 end subroutine apply_reflection_corner_south_west_front
 
 subroutine apply_reflection_corner_south_west_rear(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     at the south west rear corner
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
    dtemp(1-kk,ny+kkk,1-kkkk)=dtemp(kk,ny-kkk+1,kkkk)
  end forall
  
  return
  
 end subroutine apply_reflection_corner_south_west_rear
 
 subroutine apply_reflection_corner_south_east_rear(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the reflection 
!     at the south east rear corner
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
    dtemp(nx+kk,ny+kkk,1-kkkk)=dtemp(nx-kk+1,ny-kkk+1,kkkk)
  end forall
  
  return
  
 end subroutine apply_reflection_corner_south_east_rear
 
!******************END PART TO MANAGE THE REFLECTION********************
 
 end module fluids_mod
