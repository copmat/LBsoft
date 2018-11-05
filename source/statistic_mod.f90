
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
 
 use version_mod,           only : time_world,idrank,mxrank,sum_world_farr, &
  max_world_farr,min_world_farr,get_sync_world,finalize_world
 use fluids_mod,            only : nx,ny,nz,rhoR,rhoB,u,v,w, &
  lsingle_fluid, minx, maxx, miny, maxy, minz, maxz
 use particles_mod,         only : lparticles,engke,engcfg,engtot, &
  engrot,tempboltz,degfre
 
 implicit none
 
 private
 
 integer, public, parameter :: nmaxstatdata=20
 
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
  
  real(kind=PRC) :: dnorm,dsum,dtemp(10)
  
  
  call compute_elapsed_cpu_time()
  call compute_cpu_time()
  
  dnorm=real(nx*ny*nz,kind=PRC)
  
  forall(i=1:nmaxstatdata)statdata(i)=ZERO
  
! store all the observables in the statdata array to be printed
  
  statdata(1)=sum(rhoR(minx:maxx,miny:maxy,minz:maxz))
  
  statdata(3)=maxval(rhoR(minx:maxx,miny:maxy,minz:maxz))
  statdata(5)=minval(rhoR(minx:maxx,miny:maxy,minz:maxz))
  
  if(.not. lsingle_fluid)then
    statdata(2)=sum(rhoB(minx:maxx,miny:maxy,minz:maxz))
    statdata(4)=maxval(rhoB(minx:maxx,miny:maxy,minz:maxz))
    statdata(6)=minval(rhoB(minx:maxx,miny:maxy,minz:maxz))
  endif
  
  statdata(7)=maxval(u(minx:maxx,miny:maxy,minz:maxz))
  statdata(8)=minval(u(minx:maxx,miny:maxy,minz:maxz))
  statdata(9)=maxval(v(minx:maxx,miny:maxy,minz:maxz))
  statdata(10)=minval(v(minx:maxx,miny:maxy,minz:maxz))
  statdata(11)=maxval(w(minx:maxx,miny:maxy,minz:maxz))
  statdata(12)=minval(w(minx:maxx,miny:maxy,minz:maxz))
! estimate the remaining CPU time
  statdata(13)=meancputime*(reprinttime-icount)
  statdata(14)=elapsedcputime
  statdata(15)=meancputime
  statdata(16)=timesub
  
! these were already sum all reduced in particle_mod.f90 
  statdata(17)=engke
  statdata(18)=engcfg
  statdata(19)=engtot
! particle temperature as ratio of KbT
  statdata(20)=TWO*(engke+engrot)/(tempboltz*degfre)
  
  if(mxrank>1)then
    dtemp(1)=statdata(1)
    dtemp(2)=statdata(2)
    call sum_world_farr(dtemp,2)
    statdata(1)=dtemp(1)
    if(.not. lsingle_fluid)statdata(2)=dtemp(2)
    
    dtemp(1)=statdata(3)
    dtemp(2)=statdata(4)
    dtemp(3)=statdata(7)
    dtemp(4)=statdata(9)
    dtemp(5)=statdata(11)
    call max_world_farr(dtemp,5)
    statdata(3)=dtemp(1)
    if(.not. lsingle_fluid)statdata(4)=dtemp(2)
    statdata(7)=dtemp(3)
    statdata(9)=dtemp(4)
    statdata(11)=dtemp(5)
    
    dtemp(1)=statdata(5)
    dtemp(2)=statdata(6)
    dtemp(3)=statdata(8)
    dtemp(4)=statdata(10)
    dtemp(5)=statdata(12)
    call min_world_farr(dtemp,5)
    statdata(5)=dtemp(1)
    if(.not. lsingle_fluid)statdata(6)=dtemp(2)
    statdata(8)=dtemp(3)
    statdata(10)=dtemp(4)
    statdata(12)=dtemp(5)
  endif
  
  statdata(1)=statdata(1)/dnorm
  if(.not. lsingle_fluid)statdata(2)=statdata(2)/dnorm
  
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


