
#include <default_macro.h>
 module io_mod
 
!***********************************************************************
!     
!     LBsoft module for input/output routines
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
 use fluids_mod,            only : nx,ny,nz,set_initial_dist_type, &
  set_mean_value_dens_fluids,set_stdev_value_dens_fluids,idistselect, &
  meanR,meanB,stdevR,stdevB,set_initial_dim_box,initial_u,initial_v, &
  initial_w,set_mean_value_vel_fluids,set_boundary_conditions_type, &
  ibctype,set_value_ext_force_fluids,ext_fu,ext_fv,ext_fw
  
 implicit none

 private
 
 integer, parameter :: maxlen=150
 
 integer, public, protected, save :: init_seed=317
 
 real(kind=PRC), public, save :: timcls=0.d0
 real(kind=PRC), public, save :: timjob=0.d0
 
 public :: print_logo
 public :: read_input
 public :: print_memory_registration
 
 contains
 
 subroutine print_logo(iu)
 
!***********************************************************************
!     
!     LBsoft subroutine for printing the logo
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
 
 subroutine print_memory_registration(iu,mybanner,mymemory)
 
!***********************************************************************
!     
!     LBsoft subroutine for printing the memory registration
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: iu
  character(len=*), intent(in) :: mybanner
  real(kind=PRC), intent(in) :: mymemory
  real(kind=PRC) :: mymemoryt,mymega,mykilo
  
  character(len=12) :: r_char,r_char2
  
  character(len=*),parameter :: of='(a)'
  
  real(kind=PRC), parameter :: convert1=real(1024.d0,kind=PRC)
  real(kind=PRC), parameter :: convert2=real(1000.d0,kind=PRC)
  
  mymemoryt=mymemory/convert1
  mymega=floor(mymemoryt)
  mykilo=(mymemoryt-mymega)*convert2
  
  if(idrank/=0)return
  write(iu,of)"                                                                               "
  write(iu,of)"********************************MEMORY MONITOR*********************************"
  write(iu,of)"                                                                               "
  write (r_char,'(i12)')nint(mymemoryt)
  write (r_char2,'(i12)')nint(mykilo)
  write(iu,'(6a)')trim(mybanner)," = ",trim(adjustl(r_char)),".",trim(adjustl(r_char2))," (mb)"
  write(iu,of)"                                                                               "
  write(iu,of)"*******************************************************************************"
  write(iu,of)"                                                                               "
  
  return
  
 end subroutine print_memory_registration
 
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
  integer :: temp_ibcx=0
  integer :: temp_ibcy=0
  integer :: temp_ibcz=0
  logical :: temp_ibc=.false.
  logical :: linit_seed=.false.
  real(kind=PRC) :: dtemp_meanR = ZERO
  real(kind=PRC) :: dtemp_meanB = ZERO
  real(kind=PRC) :: dtemp_stdevR = ZERO
  real(kind=PRC) :: dtemp_stdevB = ZERO
  real(kind=PRC) :: dtemp_initial_u = ZERO
  real(kind=PRC) :: dtemp_initial_v = ZERO
  real(kind=PRC) :: dtemp_initial_w = ZERO
  real(kind=PRC) :: dtemp_ext_fu = ZERO
  real(kind=PRC) :: dtemp_ext_fv = ZERO
  real(kind=PRC) :: dtemp_ext_fw = ZERO
  real(kind=PRC) :: prntim = ZERO
  
  integer, parameter :: dimprint=28
  character(len=dimprint) :: mystring
    
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
      if(.not.safe)then
        call warning(1,dble(iline),redstring)
        call error(5)
      endif
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
            if(.not.safe)then
              call warning(1,dble(iline),redstring)
              call error(5)
            endif
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
                else
                  call warning(1,dble(iline),redstring)
                  call error(6)
                endif
              endif
            elseif(findstring('bound',directive,inumchar,maxlen))then
              if(findstring('cond',directive,inumchar,maxlen))then
                temp_ibc=.true.
                temp_ibcx=intstr(directive,maxlen,inumchar)
                temp_ibcy=intstr(directive,maxlen,inumchar)
                temp_ibcz=intstr(directive,maxlen,inumchar)
              endif
            elseif(findstring('seed',directive,inumchar,maxlen))then
              init_seed=intstr(directive,maxlen,inumchar)
              linit_seed=.true.
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
            else
              call warning(1,dble(iline),redstring)
              call error(6)
            endif
          enddo
        elseif(findstring('lb',directive,inumchar,maxlen))then
          lredo2=.true.
          do while(lredo2)
            call getline(safe,inputunit,maxlen,redstring)
            if(.not.safe)then
             call warning(1,dble(iline),redstring)
             call error(5)
            endif
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
              elseif(findstring('gauss',directive,inumchar,maxlen))then
                temp_idistselect=1
              else
                call warning(1,dble(iline),redstring)
                call error(6)
              endif
            elseif(findstring('veloc',directive,inumchar,maxlen))then
              if(findstring('mean',directive,inumchar,maxlen))then
                dtemp_initial_u=dblstr(directive,maxlen,inumchar)
                dtemp_initial_v=dblstr(directive,maxlen,inumchar)
                dtemp_initial_w=dblstr(directive,maxlen,inumchar)
              else
                call warning(1,dble(iline),redstring)
                call error(6)
              endif
            elseif(findstring('ext',directive,inumchar,maxlen))then
              if(findstring('force',directive,inumchar,maxlen))then
                dtemp_ext_fu=dblstr(directive,maxlen,inumchar)
                dtemp_ext_fv=dblstr(directive,maxlen,inumchar)
                dtemp_ext_fw=dblstr(directive,maxlen,inumchar)
              else
                call warning(1,dble(iline),redstring)
                call error(6)
              endif
            elseif(findstring('[end room',directive,inumchar,maxlen))then
              lredo2=.false.
            elseif(findstring('[end',directive,inumchar,maxlen))then
              lredo=.false.
              lredo2=.false.
            else
              call warning(1,dble(iline),redstring)
              call error(6)
            endif
          enddo
        endif
      elseif(findstring('[end',directive,inumchar,maxlen))then
         lredo=.false.
      else
        call warning(1,dble(iline),redstring)
        call error(4)
      endif
    enddo
    
    close(inputunit)
    
  endif
  
! send the read parameters to all the nodes and print them on terminal
  if(idrank==0)then
    write(6,'(/,3a,/)')repeat('*',26),"parameters from input file",repeat('*',27)
  endif
  
  if(idrank==0)then
    write(6,'(/,3a,/)')repeat('*',27),"parameters of system room",repeat('*',27)
  endif
  
  call bcast_world(lbox)
  if(lbox)then
    call bcast_world(temp_nx)
    call bcast_world(temp_ny)
    call bcast_world(temp_nz)
    call set_initial_dim_box(temp_nx,temp_ny,temp_nz)
    if(idrank==0)then
      mystring=repeat(' ',dimprint)
      mystring='system dimensions'
      write(6,'(2a,3i8)')mystring,": ",nx,ny,nz
    endif
  else
    call error(3)
  endif
  
  call bcast_world(temp_ibc)
  call bcast_world(temp_ibcx)
  call bcast_world(temp_ibcy)
  call bcast_world(temp_ibcz)
  if(.not. temp_ibc)call error(10)
  call set_boundary_conditions_type(temp_ibcx,temp_ibcy,temp_ibcz)
  ! check only full periodic
  if(ibctype.ne.7)then
    call warning(3)
    call error(9)
  endif
  if(idrank==0)then
    mystring=repeat(' ',dimprint)
    mystring='boundary conditions'
    write(6,'(2a,3i8)')mystring,": ",temp_ibcx,temp_ibcy,temp_ibcz
  endif
  
  call bcast_world(init_seed)
  if(linit_seed)then
    mystring=repeat(' ',dimprint)
    mystring='initial random seed'
    write(6,'(2a,i8)')mystring,": ",init_seed
  endif
  
  if(.not. ltimjob)then
    timjob=1.0d6*365.25d0*24.d0*60.d0*60.d0
  else
    call bcast_world(timjob)
    if(idrank==0)then
      call get_prntime(hms,timjob,prntim)
      mystring=repeat(' ',dimprint)
      mystring="user allocated job time ("//hms//") "
      write(6,'(2a,f8.4)')mystring,": ",prntim
      mystring=repeat(' ',dimprint)
      mystring='job closure time (s) '
      write(6,'(2a,f8.4)')mystring,": ",timcls
    endif
  endif
  
  if(idrank==0)then
    write(6,'(/,3a,/)')repeat('*',29),"parameters of LB room",repeat('*',29)
  endif
  
  call bcast_world(temp_idistselect)
  if(temp_idistselect>0)then
    call set_initial_dist_type(temp_idistselect)
    if(idrank==0)then
      select case(idistselect)
      case(1)
        write(6,'(a)')"initial density distribution: gaussian"
      end select
    endif
  endif
  
  call bcast_world(dtemp_meanR)
  call bcast_world(dtemp_meanB)
  if(dtemp_meanR>ZERO .or. dtemp_meanB>ZERO)then
    call set_mean_value_dens_fluids(dtemp_meanR,dtemp_meanB)
    if(idrank==0)then
      mystring=repeat(' ',dimprint)
      mystring='initial density mean values'
      write(6,'(2a,2f16.8)')mystring,": ",meanR,meanB
    endif
  endif
  
  call bcast_world(dtemp_stdevR)
  call bcast_world(dtemp_stdevB)
  if(dtemp_stdevR>ZERO .or. dtemp_stdevB>ZERO)then
    call set_stdev_value_dens_fluids(dtemp_stdevR,dtemp_stdevB)
    if(idrank==0)then
      mystring=repeat(' ',dimprint)
      mystring='initial density stdev values'
      write(6,'(2a,2f16.8)')mystring,": ",stdevR,stdevB
    endif
  endif
  
  call bcast_world(dtemp_initial_u)
  call bcast_world(dtemp_initial_v)
  call bcast_world(dtemp_initial_w)
  if(dtemp_initial_u>ZERO .or. dtemp_initial_v>ZERO .or. &
   dtemp_initial_w>ZERO)then
    call set_mean_value_vel_fluids(dtemp_initial_u,dtemp_initial_v, &
     dtemp_initial_w)
    if(idrank==0)then
      mystring=repeat(' ',dimprint)
      mystring='initial velocity mean values'
      write(6,'(2a,3f16.8)')mystring,": ",initial_u,initial_v,initial_w
    endif
  endif
  
  call bcast_world(dtemp_ext_fu)
  call bcast_world(dtemp_ext_fv)
  call bcast_world(dtemp_ext_fw)
  if(dtemp_ext_fu/=ZERO .or. dtemp_ext_fv/=ZERO .or. &
   dtemp_ext_fw/=ZERO)then
    call set_value_ext_force_fluids(dtemp_ext_fu,dtemp_ext_fv, &
     dtemp_ext_fw)
    if(idrank==0)then
      mystring=repeat(' ',dimprint)
      mystring='external force on fluids'
      write(6,'(2a,3f16.8)')mystring,": ",ext_fu,ext_fv,ext_fw
    endif
  endif
  
  if(idrank==0)then
    write(6,'(/,3a,/)')repeat('*',32),"end input file",repeat('*',33)
  endif

! print warning and check error
  
  return
  
 end subroutine read_input
 
 end module io_mod


