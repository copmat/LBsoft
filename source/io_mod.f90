
#include <default_macro.h>
 module io_mod
 
!***********************************************************************
!     
!     JETSPIN module for input/output routines
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification September 2017
!     
!***********************************************************************
 
 use version_mod
 use parse_module
 use error_mod
 use utility_mod,           only : write_fmtnumb,pi,get_prntime
 use fluids_mod,            only : set_initial_dist_type, &
  set_mean_value_dens_fluids,set_stdev_value_dens_fluids,idistselect, &
  meanR,meanB,stdevR,stdevB
  
 implicit none

 private
 
 integer, parameter :: maxlen=150
 
 
 
 real(kind=PRC), public, save :: timcls=0.d0
 real(kind=PRC), public, save :: timjob=0.d0
 
 public :: print_logo
 public :: read_input
 
 contains
 
 subroutine print_logo(iu)
 
!***********************************************************************
!     
!     JETSPIN subroutine for printing the logo
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification January 2017
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: iu
  character(len=*),parameter :: of='(a)'
  character(len=50) :: stringversion
  
  if(idrank/=0)return
  
  call print_version(stringversion)
  
  write(iu,*)
  write(iu,*)
  write(iu,of)"*******************************************************************************"
  write(iu,of)"*                                                                             *"
  write(iu,of)"*                                 O U T P U T                                 *"
  write(iu,of)"*                                    from                                     *"
  write(iu,of)"*         __      .______        _______.  ______    _______ .___________.    *"
  write(iu,of)"*        |  |     |   _  \      /       | /  __  \  |   ____||           |    *"
  write(iu,of)"*        |  |     |  |_)  |    |   (----`|  |  |  | |  |__   `---|  |----`    *"
  write(iu,of)"*        |  |     |   _  <      \   \    |  |  |  | |   __|      |  |         *"
  write(iu,of)"*        |  `----.|  |_)  | .----)   |   |  `--'  | |  |         |  |         *"
  write(iu,of)"*        |_______||______/  |_______/     \______/  |__|         |__|         *"
  write(iu,of)"*                                                                             *"
  write(iu,of)"*                                                                             *"
  write(iu,of)"*******************************************************************************"
  write(iu,of)"*                                                                             *"
  write(iu,of)"*         =========================================================           *"
  write(iu,of)"*              LBsoft is a specific-purpose open-source software              *"
  write(iu,of)"*                    for soft glassy emulsion simulations                     *"
  write(iu,of)"*         =========================================================           *"
  write(iu,of)"*                                                                             *"
  write(iu,of)"*                                                                             *"
  write(iu,of)"*    Version 0.01 (July 2018)                                                 *"
!  if(ldevelopers) &
!  write(iu,of)"*    Compiled in developer mode                                               *"
  write(iu,of)"*                                                                             *"
  write(iu,of)"*                                                                             *"
  write(iu,of)"*    The code was written by:                                                 *"
  write(iu,of)"*                                                                             *"
  write(iu,of)"*    Fabio Bonaccorso         IIT-CLNS, Rome                    Italy         *"
  write(iu,of)"*    Marco Lauricella         IAC-CNR, Rome                     Italy         *"
  write(iu,of)"*    Andrea Montessori        IAC-CNR, Rome                     Italy         *"
  write(iu,of)"*                                                                             *"
  write(iu,of)"*    with contributions from:                                                 *"
  write(iu,of)"*                                                                             *"
  write(iu,of)"*    Massimo Bernaschi        IAC-CNR, Rome                     Italy         *"
  write(iu,of)"*    Sauro Succi              IAC-CNR, Rome                     Italy         *"
  write(iu,of)"*                                                                             *"
  write(iu,of)"*                                                                             *"
  write(iu,of)"*    This is an experimental code. The authors accept no responsibility       *"
  write(iu,of)"*    for the performance of the code or for the correctness of the results.   *"
  write(iu,of)"*                                                                             *"
  write(iu,of)"*    The code is licensed under Open Software License v. 3.0 (OSL-3.0).       *"
  write(iu,of)"*    The full text of the licence can be found on the website:                *"
  write(iu,of)"*    http://opensource.org/licenses/OSL-3.0                                   *"
  write(iu,of)"*                                                                             *"
  write(iu,of)"*    A brief explanation of this license is available on the website:         *"
  write(iu,of)"*    http://rosenlaw.com/OSL3.0-explained.htm                                 *"
  write(iu,of)"*                                                                             *"
  write(iu,of)"*    The software development process has received funding from the           *"
  write(iu,of)"*    European Research Council under the Horizon 2020 Programme               *"
  write(iu,of)"*    Grant Agreement n. 739964 (COPMAT).                                      *"
  write(iu,of)"*                                                                             *"
  write(iu,of)"*                                                                             *"
  write(iu,of)"*                                                                             *"
  write(iu,of)"*******************************************************************************"
  write(iu,of)"*                                                                             *"
  write(iu,'(a5,a50,a24)')"*    ",stringversion,"                       *"
  write(iu,of)"*                                                                             *"
  write(iu,of)"*******************************************************************************"
  write(iu,of)"                                                                               "
  
  return
  
 end subroutine print_logo
 
 subroutine read_input(inputunit,inputname)
 
!***********************************************************************
!     
!     LBsoft subroutine for reading the input file
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: inputunit
  character(len=*), intent(in) :: inputname
  
  character(len=maxlen) :: redstring,directive
  logical :: safe,lredo,lredo2,ltest,lprintlisterror,lprintlisterror2
  logical :: ltestread,lexists,lfoundprint
  integer :: inumchar,i,nwords,iline,itest
  character(len=maxlen) ,allocatable :: outwords(:),outwords2(:)
  
  character(len=1) :: hms=' '
  logical :: lbox=.false.
  logical :: ltimjob=.false.
  integer :: temp_idistselect=0
  integer :: temp_nx=0
  integer :: temp_ny=0
  integer :: temp_nz=0
  real(kind=PRC) :: dtemp_meanR = ZERO
  real(kind=PRC) :: dtemp_meanB = ZERO
  real(kind=PRC) :: dtemp_stdevR = ZERO
  real(kind=PRC) :: dtemp_stdevB = ZERO
  real(kind=PRC) :: prntim = ZERO
  
! initialize parameters  

  
! note the parameters are read only by the zero node
  if(idrank==0)then
    
!   check if the input file exist
    inquire(file=inputname,exist=lexists)
    if(.not.lexists)call error(1)
    
!   open the inout file
    open(unit=inputunit,file=inputname,status='old',action='read', &
    iostat=itest)
    if(itest/=0)call error(2)
    
!   initialize the lredo condition
    lredo=.true.
    
!   counter the line which are read in the input file
    iline=0
    
!   read the input file as long as the finish directive is read
    do while(lredo)
      call getline(safe,inputunit,maxlen,redstring)
      if(.not.safe)call error(5)
      iline=iline+1
      call strip(redstring,maxlen)
      call lowcase(redstring,maxlen)
      call copystring(redstring,directive,maxlen)
 
      if(redstring(1:1)=='#')then
        cycle
      elseif(redstring(1:1)=='!')then
        cycle
      elseif(redstring(1:1)==' ')then
        cycle
      elseif(findstring('[room',directive,inumchar,maxlen))then
        if(findstring('system',directive,inumchar,maxlen))then
          lredo2=.true.
          do while(lredo2)
            call getline(safe,inputunit,maxlen,redstring)
            if(.not.safe)call error(5)
            iline=iline+1
            call strip(redstring,maxlen)
            call lowcase(redstring,maxlen)
            call copystring(redstring,directive,maxlen)
            if(redstring(1:1)=='#')then
              cycle
            elseif(redstring(1:1)=='!')then
              cycle
            elseif(redstring(1:1)==' ')then
              cycle
            elseif(findstring('job time',directive,inumchar,maxlen))then
              ltimjob=.true.
              if(findstring('indef',directive,inumchar,maxlen))then
                timjob=1.0d6*365.25d0*24.d0*60.d0*60.d0
              else
                timjob=dblstr(directive,maxlen,inumchar)
                if(findstring('m',directive,inumchar,maxlen))then
                  timjob=6.0d1*timjob
                elseif(findstring('h',directive,inumchar,maxlen))then
                  timjob=3.6d3*timjob
                elseif(findstring('d',directive,inumchar,maxlen))then
                  timjob=8.64d4*timjob
                endif
              endif
            elseif(findstring('close time',directive,inumchar,maxlen))then
              timcls=dblstr(directive,maxlen,inumchar)
            elseif(findstring('box',directive,inumchar,maxlen))then
              temp_nx=intstr(directive,maxlen,inumchar)
              temp_ny=intstr(directive,maxlen,inumchar)
              temp_nz=intstr(directive,maxlen,inumchar)
              lbox=.true.
            elseif(findstring('[end room',directive,inumchar,maxlen))then
              lredo2=.false.
            elseif(findstring('[end',directive,inumchar,maxlen))then
              lredo=.false.
              lredo2=.false.
            endif
          enddo
        elseif(findstring('hydro_vals',directive,inumchar,maxlen))then
          lredo2=.true.
          do while(lredo2)
            call getline(safe,inputunit,maxlen,redstring)
            if(.not.safe)call error(5)
            iline=iline+1
            call strip(redstring,maxlen)
            call lowcase(redstring,maxlen)
            call copystring(redstring,directive,maxlen)
            if(redstring(1:1)=='#')then
              cycle
            elseif(redstring(1:1)=='!')then
              cycle
            elseif(redstring(1:1)==' ')then
              cycle
            elseif(findstring('dens',directive,inumchar,maxlen))then
              if(findstring('mean',directive,inumchar,maxlen))then
                dtemp_meanR=dblstr(directive,maxlen,inumchar)
                dtemp_meanB=dblstr(directive,maxlen,inumchar)
              elseif(findstring('stdev',directive,inumchar,maxlen))then
                dtemp_stdevR=dblstr(directive,maxlen,inumchar)
                dtemp_stdevB=dblstr(directive,maxlen,inumchar)
              elseif(findstring('gaussian',directive,inumchar,maxlen))then
                temp_idistselect=1
              endif
            elseif(findstring('[end room',directive,inumchar,maxlen))then
              lredo2=.false.
            elseif(findstring('[end',directive,inumchar,maxlen))then
              lredo=.false.
              lredo2=.false.
            endif
          enddo
        endif
      elseif(findstring('[end',directive,inumchar,maxlen))then
         lredo=.false.
      else
        call warning(37,1.d0,redstring)
        call error(6)
      endif
    enddo
    
    close(inputunit)
    
  endif
  
  if(idrank==0)then
    write(6,'(/,3a,/)')repeat('*',10),"parameters from input file",repeat('*',10)
  endif
  
  call bcast_world(lbox)
  
! send the read parameters to all the nodes
  if(temp_idistselect>0)then
    call set_initial_dist_type(temp_idistselect)
    if(idrank==0)then
      select case(idistselect)
      case(1)
        write(6,'(a)')"initial density distribution: gaussian"
      end select
    endif
  endif
  
  if(dtemp_meanR>ZERO .or. dtemp_meanB>ZERO)then
    call set_mean_value_dens_fluids(dtemp_meanR,dtemp_meanB)
    if(idrank==0)then
      write(6,'(a,2f16.8)')"initial density mean values: ",meanR,meanB
    endif
  endif
  
  if(dtemp_stdevR>ZERO .or. dtemp_stdevB>ZERO)then
    call set_stdev_value_dens_fluids(dtemp_stdevR,dtemp_stdevB)
    if(idrank==0)then
      write(6,'(a,2f16.8)')"initial density stdev values: ",stdevR,stdevB
    endif
  endif
  
  if(.not. ltimjob)then
    timjob=1.0d6*365.25d0*24.d0*60.d0*60.d0
  else
    if(idrank==0)then
      call get_prntime(hms,timjob,prntim)
      write(6,'(a,a1,a,3x,f8.4)')"user allocated job time (",hms,") ",prntim
      write(6,'(a,a1,a,3x,f8.4)')"job closure time (","s",") ",timcls
    endif
  endif
  
  if(idrank==0)then
    write(6,'(/,3a,/)')repeat('*',16),"end input file",repeat('*',16)
  endif
stop
! print warning and check error
  
  return
  
 end subroutine read_input
 
 end module io_mod


