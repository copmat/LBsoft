
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
   !    0,1, 2,3,4,5,6,7, 8, 9,10,11,12,13,14,15,16,17,18
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
 
 end module fluids_mod
