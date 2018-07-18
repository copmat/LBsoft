
#include <default_macro.h>
 module statistic_mod
 
!***********************************************************************
!     
!     LBsoft module containing statistical data routines
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
 use version_mod,           only : time_world
 use fluids_mod,            only : nx,ny,nz,rhoR,rhoB,u,v,w, &
  lsingle_fluid
 
 implicit none
 
 private
 
 integer, public, parameter :: nmaxstatdata=16
 
 real(kind=PRC), public, save, dimension(nmaxstatdata) :: statdata
 real(kind=PRC), public, save :: meancputime=0.d0
 real(kind=PRC), public, save :: elapsedcputime=0.d0
 
 integer, public, save :: reprinttime=0
 
 
 public :: statistic_driver
 public :: compute_statistic
 
 contains
 
 subroutine statistic_driver(nstepsub,timesub)
 
!***********************************************************************
!     
!     LBsoft subroutine for controlling calls to
!     subroutines which store statistical data 
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: nstepsub
  real(kind=PRC), intent(in) :: timesub
  
! update the counter for all the observables

  
  return
  
 end subroutine statistic_driver
 
 subroutine compute_statistic(tempint,timesub,nstepsub)
 
!***********************************************************************
!     
!     LBsoft subroutine for controlling subroutine calls to
!     subroutines which compute statistical data 
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC), intent(in) :: tempint,timesub
  integer, intent(in) :: nstepsub
  
  integer, save :: icount=0
  integer :: i
  
  integer, save :: nstepsubold=0
  integer, save :: nmulstepdoneold=0
  
  real(kind=PRC) :: dnorm,dsum
  
  
  call compute_elapsed_cpu_time()
  call compute_cpu_time()
  
  dnorm=real(nx*ny*nz,kind=PRC)
  
  forall(i=1:nmaxstatdata)statdata(i)=ZERO
  
! store all the observables in the statdata array to be printed
  
  statdata(1)=sum(rhoR(1:nx,1:ny,1:nz))/dnorm
  statdata(3)=maxval(rhoR(1:nx,1:ny,1:nz))
  statdata(5)=minval(rhoR(1:nx,1:ny,1:nz))
  
  if(.not. lsingle_fluid)then
    statdata(2)=sum(rhoB(1:nx,1:ny,1:nz))/dnorm
    statdata(4)=maxval(rhoB(1:nx,1:ny,1:nz))
    statdata(6)=minval(rhoB(1:nx,1:ny,1:nz))
  endif
  
  statdata(7)=maxval(u(1:nx,1:ny,1:nz))
  statdata(8)=minval(u(1:nx,1:ny,1:nz))
  statdata(9)=maxval(v(1:nx,1:ny,1:nz))
  statdata(10)=minval(v(1:nx,1:ny,1:nz))
  statdata(11)=maxval(w(1:nx,1:ny,1:nz))
  statdata(12)=minval(w(1:nx,1:ny,1:nz))
! estimate the remaining CPU time
  statdata(13)=meancputime*(reprinttime-icount)
  statdata(14)=elapsedcputime
  statdata(15)=meancputime
  statdata(16)=timesub
  
! update the counter
  icount=icount+1
  
  return
  
 end subroutine compute_statistic
 
 subroutine compute_cpu_time()
 
!***********************************************************************
!     
!     LBsoft subroutine for registering the actual CPU time
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC), save :: timeold=0.d0
  real(kind=PRC) :: timenew
  
  call time_world(timenew)
  meancputime=timenew-timeold
  timeold=timenew
  
  return
  
 end subroutine compute_cpu_time
 
 subroutine compute_elapsed_cpu_time()
 
!***********************************************************************
!     
!     LBsoft subroutine for computing the elapsed CPU time
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC), save :: initialtime=0.d0
  logical, save :: lfirst=.true.
  
  real(kind=PRC) :: timenew
  
  if(lfirst)then
    lfirst=.false.
    call time_world(initialtime)
  endif
  
  call time_world(timenew)
  elapsedcputime=timenew-initialtime
  
  return
  
 end subroutine compute_elapsed_cpu_time
 
 end module statistic_mod


