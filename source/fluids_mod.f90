
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
 
 implicit none
 
 private
 
 integer, save, protected, public :: nx=10
 integer, save, protected, public :: ny=10
 integer, save, protected, public :: nz=10
 
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
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: psix
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: psiy
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: psiz
 
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:,:) :: fR
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:,:) :: fB
 
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:,:) :: fR_b
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:,:) :: fB_b
 
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
 
 public :: set_random_dens_fluids
 public :: set_initial_dens_fluids
 public :: set_initial_vel_fluids
 public :: set_mean_value_dens_fluids
 public :: set_stdev_value_dens_fluids
 public :: set_initial_dist_type
 
 contains
 
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
 
 end module fluids_mod
