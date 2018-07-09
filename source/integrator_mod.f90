
#include <default_macro.h>
 module integrator_mod
 
!***********************************************************************
!     
!     LBsoft module containing integrators data routines
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification july 2018
!     
!***********************************************************************


 implicit none

 private
 
 integer, public, protected, save :: integrator
 integer, public, protected, save :: nstep=0
 logical, public, protected, save :: lintegrator=.false.
 
 real(kind=PRC), public, save :: tstep=0.d0
 real(kind=PRC), public, save :: initime = ZERO
 real(kind=PRC), public, save :: endtime = ZERO
 logical, public, save :: lendtime
 
 public :: set_nstep
 public :: update_nstep
 
 contains
 
 subroutine set_nstep(itemp)
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the nstep counter
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification july 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: itemp
  
  nstep = itemp
  
  return
  
 end subroutine set_nstep
 
 subroutine update_nstep
 
!***********************************************************************
!     
!     LBsoft subroutine for updating the nstep counter
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification july 2018
!     
!***********************************************************************
 
  implicit none
  
  nstep = nstep + 1
  
  return
  
 end subroutine update_nstep

 end module integrator_mod

