
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
 
 integer, public, save :: integrator
 logical, public, save :: lintegrator=.false.
 
 real(kind=PRC), public, save :: tstep=0.d0
 real(kind=PRC), public, save :: initime = ZERO
 real(kind=PRC), public, save :: endtime = ZERO
 logical, public, save :: lendtime
 

 
 end module integrator_mod

