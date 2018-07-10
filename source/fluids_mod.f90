
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
 
 use version_mod,    only : idrank
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
 
 logical, save, protected, public :: lalloc_hvar=.false.
 logical, save, protected, public :: lalloc_pop=.false.
 
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
 
 
 integer, save, protected, public, allocatable, dimension(:,:,:) :: isfluid
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: rhoR
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: rhoB
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: u
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: v
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: w
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: fuR
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: fvR
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: fwR
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: fuB
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: fvB
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: fwB
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: psixR
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: psiyR
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: psizR
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: psixB
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: psiyB
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: psizB
 
#if LATTICE==319
 
 integer, parameter, public :: latt_dim=3
 integer, parameter, public :: links=18
 
 integer, parameter, public :: nbuff=1
 
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
  
  allocate(psixR(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(12))
  allocate(psiyR(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(13))
  allocate(psizR(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(14))
  
  allocate(psixB(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(15))
  allocate(psiyB(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(16))
  allocate(psizB(myzero:mynx,myzero:myny,myzero:mynz),stat=istat(17))
  
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
  
  if(any(istat.ne.0))then
    do i=1,nistatmax
      if(istat(i).ne.0)exit
    enddo
    call warning(2,dble(i))
    call error(4)
  endif
  
  call allocate_array_buffservice3d(myzero,mynx,myzero,myny,myzero,mynz)
  
  return
  
 end subroutine allocate_fluids
 
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
  
  return
  
 end subroutine initialize_fluids
 
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
 
 subroutine apply_pbc_hvar
 
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
  
  !apply pbc at north west front 5
  call apply_pbc_corner_north_west_front(rhoR)
  call apply_pbc_corner_north_west_front(rhoB)
  
  call apply_pbc_corner_north_west_front(u)
  call apply_pbc_corner_north_west_front(v)
  call apply_pbc_corner_north_west_front(w)
  
  !apply pbc at south east front 6
  call apply_pbc_corner_south_east_front(rhoR)
  call apply_pbc_corner_south_east_front(rhoB)
  
  call apply_pbc_corner_south_east_front(u)
  call apply_pbc_corner_south_east_front(v)
  call apply_pbc_corner_south_east_front(w)
  
  !apply pbc at south west front 7
  call apply_pbc_corner_south_west_front(rhoR)
  call apply_pbc_corner_south_west_front(rhoB)
  
  call apply_pbc_corner_south_west_front(u)
  call apply_pbc_corner_south_west_front(v)
  call apply_pbc_corner_south_west_front(w)
  
  !apply pbc at south west rear 8
  call apply_pbc_corner_south_west_rear(rhoR)
  call apply_pbc_corner_south_west_rear(rhoB)
  
  call apply_pbc_corner_south_west_rear(u)
  call apply_pbc_corner_south_west_rear(v)
  call apply_pbc_corner_south_west_rear(w)
  
  return
  
 end subroutine apply_pbc_hvar
 
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
 
 end module fluids_mod
