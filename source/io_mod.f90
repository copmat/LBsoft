
#include <default_macro.h>
 module io_mod
 
!***********************************************************************
!     
!     LBsoft module for input/output routines
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
 use version_mod
 use parse_module
 use error_mod
 use profiling_mod,         only : set_value_idiagnostic, &
  set_value_ldiagnostic,idiagnostic,ldiagnostic
 use utility_mod,           only : write_fmtnumb,pi,get_prntime, &
  linit_seed,ltest_mode,conv_rad
 use lbempi_mod,            only : set_domdec,set_domain,domdec, &
  nprocx,nprocy,nprocz

 use fluids_mod,            only : nx,ny,nz,set_initial_dist_type, &
  set_mean_value_dens_fluids,set_stdev_value_dens_fluids,idistselect, &
  meanR,meanB,stdevR,stdevB,set_initial_dim_box,initial_u,initial_v, &
  initial_w,set_mean_value_vel_fluids,set_boundary_conditions_type, &
  ibctype,set_value_ext_force_fluids,ext_fu,ext_fv,ext_fw,lpair_SC, &
  pair_SC,set_fluid_force_sc,set_value_viscosity,set_value_tau, &
  viscR,viscB,tauR,tauB,lunique_omega,lforce_add,set_lsingle_fluid, &
  lsingle_fluid,latt_name,set_value_bc_east,set_value_bc_west, &
  set_value_bc_front,set_value_bc_rear,set_value_bc_north, &
  set_value_bc_south,bc_rhoR_east,bc_rhoB_east,bc_u_east,bc_v_east,bc_w_east,&
  bc_type_east,bc_rhoR_west,bc_rhoB_west,bc_u_west,bc_v_west,bc_w_west,&
  bc_type_west,bc_rhoR_north,bc_rhoB_north,bc_u_north,bc_v_north,bc_w_north,&
  bc_type_north,bc_rhoR_south,bc_rhoB_south,bc_u_south,bc_v_south,bc_w_south,&
  bc_type_south,bc_rhoR_front,bc_rhoB_front,bc_u_front,bc_v_front,bc_w_front,&
  bc_type_front,bc_rhoR_rear,bc_rhoB_rear,bc_u_rear,bc_v_rear,bc_w_rear,&
  bc_type_rear,set_fluid_wall_sc,wallR_SC,wallB_SC,LBintegrator, &
  set_LBintegrator_type,lbc_halfway,set_lbc_halfway,set_lbc_fullway, &
  lbc_fullway,partR_SC,partB_SC,set_fluid_particle_sc,theta_SC, &
  devtheta_SC,set_theta_particle_sc,set_back_value_dens_fluids, &
  backR,backB,set_init_area_dens_fluids,areaR,set_cap_force,cap_force, &
  imass_rescale,set_mass_rescale,dmass_rescale,iselwetting, &
  densR_wetting,densB_wetting,set_fluid_wall_mode, set_restore

 use particles_mod,         only : set_natms_tot,natms_tot,lparticles, &
  set_ishape,set_densvar,densvar,ishape,set_rcut,rcut,delr, &
  set_delr,rotmat_2_quat,lrotate,set_lrotate,allocate_field_array, &
  set_ntpvdw,ntpvdw,set_field_array,mxpvdw,mxvdw,ltpvdw,prmvdw, &
  set_umass,lumass,linit_temp,init_temp,set_init_temp,lvv, &
  set_lvv,set_rdim,lurdim,set_lparticles,set_atmnamtype, &
  set_ntype,ntype,mxntype,rdimx,rdimy,rdimz,atmnamtype,weight, &
  find_type,allocate_particle_features,mskvdw,set_lsidewall_md,lsidewall_md, &
  set_sidewall_md_k,sidewall_md_k,sidewall_md_k,sidewall_md_rdist, &
  llubrication,set_lubrication,set_sidewall_md_rdist, ext_fxx, &
  ext_fyy,ext_fzz,set_value_ext_force_particles,ext_tqx,ext_tqy, &
  ext_tqz,set_value_ext_torque_particles,set_cap_force_part, &
  cap_force_part,ltype,lfix_moment,ifix_moment,set_fix_moment, &
  lubricrparcut,lubricrparcap,lubricconst
 use write_output_mod,      only: set_value_ivtkevery,ivtkevery, &
  lvtkfile,set_value_ixyzevery,lxyzfile,ixyzevery, &
  set_value_ibinevery,set_value_dumpevery,idumpevery
 use integrator_mod,        only : set_nstepmax,nstepmax,tstep,endtime
 use statistic_mod,         only : reprinttime,compute_statistic, &
  statdata
  
  
 implicit none

 private
 
 integer, parameter :: maxlen=150
 
 integer, public, protected, save :: init_seed=317
 integer, public, protected, save :: nprintlist=0
 integer, public, protected, save :: iprinttime=0
 integer, public, protected, save :: eprinttime=0
 integer, public, protected, save :: nxyzlist=0
 integer, public, protected, save :: nprintlistvtk=0
 logical, public, protected, save :: lprintlist=.false.
 logical, public, protected, save :: lprintlistvtk=.false.
 logical, public, protected, save :: lprinttime=.false.
 real(kind=PRC), public, protected, save :: printtime=FIFTY
 real(kind=PRC), public, protected, save :: timcls=ZERO
 real(kind=PRC), public, protected, save :: timjob=ZERO
 
 
 character(len=20), public, protected, save, allocatable, dimension(:) :: printarg(:)
 integer, public, protected, save, allocatable, dimension(:) :: printlist
 integer, public, protected, save, allocatable, dimension(:) :: printlistvtk
 real(kind=PRC), public, protected, allocatable, save :: xprint(:)
 
 integer, public, protected, allocatable, dimension(:), save :: xyzlist
 
 public :: print_logo
 public :: read_input
 public :: print_memory_registration
 public :: outprint_driver
 public :: allocate_print
 public :: read_input_atom
 
 contains
 
 subroutine allocate_print()
  
  implicit none
  
  allocate(xprint(1:nprintlist))
  
  return
  
 end subroutine allocate_print
 
 
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
  write(iu,of)"*    Version 0.02 (January 2019)                                              *"
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
  write(iu,'(a48,a6,a25)')"*    the code was compiled with the LB scheme : ",latt_name, &
   "                        *"
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
 
 subroutine print_legend_observables(iu)
 
!***********************************************************************
!     
!     LBsoft subroutine for printing the legend of the observables
!     requested in input file
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: iu
  
  integer :: i
  
  character(len=*),parameter :: of='(a)'
  
  if(idrank/=0)return
  
  write(iu,of)"                                                                               "
  write(iu,of)"*********************PRINT LIST SPECIFIED IN INPUT FILE************************"
  write(iu,of)"                                                                               "
  do i=1,nprintlist
    write(iu,'(a)')legendobs(i,printlist)
  enddo
  write(iu,of)"                                                                               "
  write(iu,of)"*******************************************************************************"
  write(iu,of)"                                                                               "
  
  return
  
 end subroutine print_legend_observables
 
 
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
  logical :: safe,lredo,lredo2,ltest
  logical :: ltestread,lexists,lfoundprint,lfoundprintvtk
  integer :: inumchar,i,j,k,nwords,iline,itest,itype,jtype
  character(len=maxlen) ,allocatable :: outwords(:),outwords2(:)
  
  character(len=1) :: hms=' '
  character(len=8), dimension(mxntype) :: temp_atmnamtype
  logical :: lbox=.false.
  logical :: ltimjob=.false.
  integer :: temp_idistselect=0
  integer :: temp_nx=0
  integer :: temp_ny=0
  integer :: temp_nz=0
  integer :: temp_ibcx=0
  integer :: temp_ibcy=0
  integer :: temp_ibcz=0
  integer :: temp_nstepmax=0
  integer :: temp_idiagnostic=1
  integer :: temp_ivtkevery=1

  integer :: temp_ibinevery   =10000000
  integer :: temp_idumpevery  =10000

  integer :: temp_ixyzevery=1
  integer :: temp_nfluid=0
  integer :: temp_bc_type_east=0
  integer :: temp_bc_type_west=0
  integer :: temp_bc_type_front=0
  integer :: temp_bc_type_rear=0
  integer :: temp_bc_type_north=0
  integer :: temp_bc_type_south=0
  integer :: temp_domdec=7
  integer :: temp_nprocx=1
  integer :: temp_nprocy=1
  integer :: temp_nprocz=1
  integer :: ifield_pair=0
  integer, dimension(mxntype) :: temp_ishape(1:mxntype)=0
  integer, dimension(mxntype,mxntype) :: temp_mskvdw(1:mxntype,1:mxntype)=0
  integer :: temp_ntype=0
  integer :: temp_areaR(1:2,1:3)=0
  integer :: temp_imass_rescale=100
  integer :: temp_ifix_moment=1
  integer :: temp_iselwetting=0
  logical :: temp_ibc=.false.
  logical :: temp_lpair_SC=.false.
  logical :: temp_ldomdec=.false.
  logical :: temp_lnprocn=.false.
  logical :: lvisc=.false.
  logical :: ltau=.false.
  logical :: lnstepmax=.false.
  logical :: lprintlisterror=.false.
  logical :: lprintlisterrorvtk=.false.
  logical :: temp_ldiagnostic=.false.
  logical :: temp_lvtkfile=.false.
  logical :: temp_livtkfile=.false.
  logical :: temp_lxyzfile=.false.
  logical :: temp_lbinary=.false.
  logical :: lidiagnostic=.false.
  logical :: temp_lnfluid=.false.
  logical :: temp_wall_SC=.false.
  logical :: temp_lbc_halfway=.false.
  logical :: temp_lbc_fullway=.false.
  logical :: temp_lparticles=.false.
  logical :: lerror1=.false.
  logical :: lerror2=.false.
  logical :: lerror3=.false.
  logical :: lerror4=.false.
  logical :: lerror5=.false.
  logical :: lerror6=.false.
  logical :: lerror7=.false.
  logical :: lmd=.false.
  logical :: temp_densvar=.false.
  logical :: temp_rcut=.false.
  logical :: temp_delr=.false.
  logical :: temp_lrotate=.false.
  logical :: temp_field_pair=.false.
  logical, dimension(1:mxntype) :: temp_lumass(1:mxntype)=.false.
  logical :: temp_linit_temp=.false.
  logical :: temp_lvv=.false.
  logical, dimension(mxntype) :: temp_lurdim(1:mxntype)=.false.
  logical :: temp_lparticlentype=.false.
  logical :: temp_lsidewall_md=.false.
  logical :: temp_llubrication=.false.
  logical :: temp_sidewall_md_k=.false.
  logical :: temp_sidewall_md_rdist=.false.
  logical :: temp_lpart_SC=.false.
  logical :: temp_ltheta_SC=.false.
  logical :: temp_lcap_force_part=.false.
  logical :: temp_lcap_force=.false.
  logical :: temp_lmass_rescale=.false.
  logical :: temp_lfix_moment=.false.
  logical :: temp_lselwetting=.false.
  logical :: temp_lrestore  = .false.
  logical :: temp_lidumpevery= .false.

  real(kind=PRC) :: dtemp_meanR = ZERO
  real(kind=PRC) :: dtemp_meanB = ZERO
  real(kind=PRC) :: dtemp_stdevR = ZERO
  real(kind=PRC) :: dtemp_stdevB = ZERO
  real(kind=PRC) :: dtemp_backR = ZERO
  real(kind=PRC) :: dtemp_backB = ZERO
  real(kind=PRC) :: dtemp_initial_u = ZERO
  real(kind=PRC) :: dtemp_initial_v = ZERO
  real(kind=PRC) :: dtemp_initial_w = ZERO
  real(kind=PRC) :: dtemp_ext_fu = ZERO
  real(kind=PRC) :: dtemp_ext_fv = ZERO
  real(kind=PRC) :: dtemp_ext_fw = ZERO
  real(kind=PRC) :: prntim = ZERO
  real(kind=PRC) :: temp_pair_SC = ZERO
  real(kind=PRC) :: dtemp_viscR = ZERO
  real(kind=PRC) :: dtemp_viscB = ZERO
  real(kind=PRC) :: dtemp_tauR = ZERO
  real(kind=PRC) :: dtemp_tauB = ZERO
  
  real(kind=PRC) :: dtemp_wallR_SC = ONE
  real(kind=PRC) :: dtemp_wallB_SC = ONE
  
  real(kind=PRC) :: dtemp_bc_rhoR_east = ZERO
  real(kind=PRC) :: dtemp_bc_rhoR_west = ZERO
  real(kind=PRC) :: dtemp_bc_rhoR_front= ZERO
  real(kind=PRC) :: dtemp_bc_rhoR_rear = ZERO
  real(kind=PRC) :: dtemp_bc_rhoR_north= ZERO
  real(kind=PRC) :: dtemp_bc_rhoR_south= ZERO
 
  real(kind=PRC) :: dtemp_bc_rhoB_east = ZERO
  real(kind=PRC) :: dtemp_bc_rhoB_west = ZERO
  real(kind=PRC) :: dtemp_bc_rhoB_front= ZERO
  real(kind=PRC) :: dtemp_bc_rhoB_rear = ZERO
  real(kind=PRC) :: dtemp_bc_rhoB_north= ZERO
  real(kind=PRC) :: dtemp_bc_rhoB_south= ZERO
 
  real(kind=PRC) :: dtemp_bc_u_east = ZERO
  real(kind=PRC) :: dtemp_bc_u_west = ZERO
  real(kind=PRC) :: dtemp_bc_u_front= ZERO
  real(kind=PRC) :: dtemp_bc_u_rear = ZERO
  real(kind=PRC) :: dtemp_bc_u_north= ZERO
  real(kind=PRC) :: dtemp_bc_u_south= ZERO
 
  real(kind=PRC) :: dtemp_bc_v_east = ZERO
  real(kind=PRC) :: dtemp_bc_v_west = ZERO
  real(kind=PRC) :: dtemp_bc_v_front= ZERO
  real(kind=PRC) :: dtemp_bc_v_rear = ZERO
  real(kind=PRC) :: dtemp_bc_v_north= ZERO
  real(kind=PRC) :: dtemp_bc_v_south= ZERO
 
  real(kind=PRC) :: dtemp_bc_w_east = ZERO
  real(kind=PRC) :: dtemp_bc_w_west = ZERO
  real(kind=PRC) :: dtemp_bc_w_front= ZERO
  real(kind=PRC) :: dtemp_bc_w_rear = ZERO
  real(kind=PRC) :: dtemp_bc_w_north= ZERO
  real(kind=PRC) :: dtemp_bc_w_south= ZERO
  
  real(kind=PRC) :: dtemp_densvar= ZERO
  real(kind=PRC) :: dtemp_rcut= ZERO
  real(kind=PRC) :: dtemp_delr= ZERO
  
  real(kind=PRC) :: dtemp_theta_SC= ZERO
  real(kind=PRC) :: dtemp_devtheta_SC= ZERO
  real(kind=PRC) :: dtemp_partR_SC= ZERO
  real(kind=PRC) :: dtemp_partB_SC= ZERO
  
  real(kind=PRC) :: dtemp_cap_force_part= ZERO
  real(kind=PRC) :: dtemp_cap_force= ZERO
  real(kind=PRC), dimension(1:mxntype) :: dtemp_umass(1:mxntype)= ZERO
  
  real(kind=PRC), dimension(mxntype) :: dtemp_urdim(1:mxntype)=ZERO
  real(kind=PRC), dimension(mxntype) :: dtemp_urdimy(1:mxntype)=ZERO
  real(kind=PRC), dimension(mxntype) :: dtemp_urdimz(1:mxntype)=ZERO
  
  real(kind=PRC) :: dtemp_init_temp= ZERO
  
  real(kind=PRC) :: dtemp_sidewall_md_k=ZERO
  real(kind=PRC) :: dtemp_sidewall_md_rdist=ONE
  
  real(kind=PRC) :: dtemp_ext_fxx = ZERO
  real(kind=PRC) :: dtemp_ext_fyy = ZERO
  real(kind=PRC) :: dtemp_ext_fzz = ZERO
  
  real(kind=PRC) :: dtemp_ext_tqx = ZERO
  real(kind=PRC) :: dtemp_ext_tqy = ZERO
  real(kind=PRC) :: dtemp_ext_tqz = ZERO
  
  real(kind=PRC) :: temp_dmass_rescale = ZERO
  
  real(kind=PRC) :: temp_lubricrparcut = ONE/TEN
  real(kind=PRC) :: temp_lubricrparcap = TWO/THREE
  real(kind=PRC) :: temp_lubricconst = ONE
  
  real(kind=PRC) :: temp_densR_wetting = ZERO
  real(kind=PRC) :: temp_densB_wetting = ZERO
  
  integer, dimension(mxvdw) :: temp_ltpvdw
  real(kind=PRC), dimension(mxpvdw,mxvdw) :: dtemp_prmvdw
  
  integer, parameter :: dimprint=36
  integer, parameter :: dimprint2=12
  character(len=dimprint) :: mystring
  character(len=dimprint) :: mystring36
  character(len=dimprint2) :: mystring12
  character(len=1) :: r_char
  
    
! initialize parameters  

  
! note the parameters are read only by the zero node
  if(idrank==0)then
  
!   initialize the lredo condition
    lredo=.true.
    
!   check if the input file exist
    inquire(file=inputname,exist=lexists)
    if(.not.lexists)then
      lerror1=.true.
      lredo=.false.
    endif
    
!   open the inout file
    open(unit=inputunit,file=inputname,status='old',action='read', &
    iostat=itest)
    if(itest/=0)then
      lerror2=.true.
      lredo=.false.
    endif
    
!   counter the line which are read in the input file
    iline=0
    
!   read the input file as long as the finish directive is read
    do while(lredo)
      call getline(safe,inputunit,maxlen,redstring)
      if(.not.safe)then
        call warning(1,dble(iline),redstring)
        lerror5=.true.
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
              lerror5=.true.
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
                  lerror6=.true.
                endif
              endif
            elseif(findstring('diagnostic',directive,inumchar,maxlen))then
              if(findstring('yes',directive,inumchar,maxlen))then
                temp_ldiagnostic=.true.
              elseif(findstring('every',directive,inumchar,maxlen))then
                temp_idiagnostic=intstr(directive,maxlen,inumchar)
                lidiagnostic=.true.
              endif
            elseif(findstring('decompos',directive,inumchar,maxlen))then
              if(findstring('type',directive,inumchar,maxlen))then
                temp_ldomdec=.true.
                temp_domdec=intstr(directive,maxlen,inumchar)
              elseif(findstring('dimen',directive,inumchar,maxlen))then
                temp_lnprocn=.true.
                temp_nprocx=intstr(directive,maxlen,inumchar)
                temp_nprocy=intstr(directive,maxlen,inumchar)   
                temp_nprocz=intstr(directive,maxlen,inumchar)   
              endif
            elseif(findstring('restart',directive,inumchar,maxlen))then
              if(findstring('yes',directive,inumchar,maxlen))then
                temp_lrestore=.true.
              elseif(findstring('every',directive,inumchar,maxlen))then
                temp_idumpevery=intstr(directive,maxlen,inumchar)
                temp_lidumpevery=.true.
              endif
            elseif(findstring('print',directive,inumchar,maxlen))then
              if(findstring('vtk',directive,inumchar,maxlen))then
                if(findstring('every',directive,inumchar,maxlen))then
                  temp_livtkfile=.true.
                  temp_ivtkevery=intstr(directive,maxlen,inumchar)
                elseif(findstring('list',directive,inumchar,maxlen))then
                  lprintlistvtk=.true.
                  lprintlisterrorvtk=.false.
                  call findwords(nwords,outwords,directive,maxlen)
                  nprintlistvtk=max(0,nwords-3)
                  directive=outwords(1)
                  if(.not.findstring('print',directive,inumchar,maxlen))then
                    lprintlisterrorvtk=.true.
                  endif
                  directive=outwords(2)
                  if(.not.findstring('vtk',directive,inumchar,maxlen))then
                    lprintlisterrorvtk=.true.
                  endif
                  directive=outwords(3)
                  if(.not.findstring('list',directive,inumchar,maxlen))then
                    lprintlisterrorvtk=.true.
                  endif
                  if(nprintlistvtk<=0)lprintlisterrorvtk=.true.
                  if(nprintlistvtk>0)allocate(printlistvtk(nprintlist))
                  if(allocated(outwords2))deallocate(outwords2)
                  if(nprintlistvtk>0)allocate(outwords2(nprintlist))
                  do i=1,nprintlistvtk
                    outwords2(i)=outwords(i+3)
                    call identify_argument_vtk(i,outwords2,printlistvtk,maxlen, &
                     lfoundprintvtk)
                    lprintlisterrorvtk=(lprintlisterrorvtk .or. (.not.lfoundprintvtk))
                  enddo
                elseif(findstring('yes',directive,inumchar,maxlen))then
                  temp_lvtkfile=.true.
                else
                  call warning(1,dble(iline),redstring)
                  lerror6=.true.
                endif
              elseif(findstring('list',directive,inumchar,maxlen))then
                lprintlist=.true.
                lprintlisterror=.false.
                call findwords(nwords,outwords,directive,maxlen)
                nprintlist=max(0,nwords-2)
                directive=outwords(1)
                if(.not.findstring('print',directive,inumchar,maxlen))then
                  lprintlisterror=.true.
                endif
                directive=outwords(2)
                if(.not.findstring('list',directive,inumchar,maxlen))then
                  lprintlisterror=.true.
                endif
                if(nprintlist>0)allocate(printlist(nprintlist))
                if(allocated(outwords2))deallocate(outwords2)
                if(nprintlist>0)allocate(outwords2(nprintlist))
                do i=1,nprintlist
                  outwords2(i)=outwords(i+2)
                  call identify_argument(i,outwords2,printlist,maxlen, &
                   lfoundprint)
                  lprintlisterror=(lprintlisterror .or. (.not.lfoundprint))
                enddo
              elseif(findstring('every',directive,inumchar,maxlen))then
                printtime=real(intstr(directive,maxlen,inumchar),kind=PRC)
                lprinttime=.true.
              elseif(findstring('xyz',directive,inumchar,maxlen))then
                temp_ixyzevery=intstr(directive,maxlen,inumchar)
                temp_lxyzfile=.true.
              elseif(findstring('binary',directive,inumchar,maxlen))then
                temp_ibinevery=intstr(directive,maxlen,inumchar)
                temp_lbinary = .true.
              else
                call warning(1,dble(iline),redstring)
                lerror6=.true.
              endif
            elseif(findstring('bound',directive,inumchar,maxlen))then
              if(findstring('cond',directive,inumchar,maxlen))then
                temp_ibc=.true.
                temp_ibcx=intstr(directive,maxlen,inumchar)
                temp_ibcy=intstr(directive,maxlen,inumchar)
                temp_ibcz=intstr(directive,maxlen,inumchar)
              else
                call warning(1,dble(iline),redstring)
                lerror6=.true.
              endif
            elseif(findstring('test',directive,inumchar,maxlen))then
              if(findstring('yes',directive,inumchar,maxlen))then
                ltest_mode=.true.
              elseif(findstring('no',directive,inumchar,maxlen))then
                ltest_mode=.false.
              else
                call warning(1,dble(iline),redstring)
                lerror6=.true.
              endif
            elseif(findstring('steps',directive,inumchar,maxlen))then
              lnstepmax=.true.
              temp_nstepmax=intstr(directive,maxlen,inumchar)
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
              lerror6=.true.
            endif
          enddo
        elseif(findstring('lb',directive,inumchar,maxlen))then
          lredo2=.true.
          do while(lredo2)
            call getline(safe,inputunit,maxlen,redstring)
            if(.not.safe)then
             call warning(1,dble(iline),redstring)
             lerror5=.true.
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
            elseif(findstring('compon',directive,inumchar,maxlen))then
              temp_nfluid=intstr(directive,maxlen,inumchar)
              temp_lnfluid=.true.
            elseif(findstring('bounc',directive,inumchar,maxlen))then
              if(findstring('half',directive,inumchar,maxlen))then
                if(findstring('yes',directive,inumchar,maxlen))then
                  temp_lbc_halfway=.true.
                elseif(findstring('no',directive,inumchar,maxlen))then
                  temp_lbc_halfway=.false.
                else
                  call warning(1,dble(iline),redstring)
                  lerror6=.true.
                endif
              elseif(findstring('full',directive,inumchar,maxlen))then
                if(findstring('yes',directive,inumchar,maxlen))then
                  temp_lbc_fullway=.true.
                elseif(findstring('no',directive,inumchar,maxlen))then
                  temp_lbc_fullway=.false.
                else
                  call warning(1,dble(iline),redstring)
                  lerror6=.true.
                endif
              else
                call warning(1,dble(iline),redstring)
                lerror6=.true.
              endif
            elseif(findstring('bound',directive,inumchar,maxlen))then
              if(findstring('open',directive,inumchar,maxlen))then
                if(findstring('east',directive,inumchar,maxlen))then
                  if(findstring('type',directive,inumchar,maxlen))then
                    temp_bc_type_east=intstr(directive,maxlen,inumchar)
                  elseif(findstring('dens',directive,inumchar,maxlen))then
                    dtemp_bc_rhoR_east=dblstr(directive,maxlen,inumchar)
                    dtemp_bc_rhoB_east=dblstr(directive,maxlen,inumchar)
                  elseif(findstring('veloc',directive,inumchar,maxlen))then
                    dtemp_bc_u_east=dblstr(directive,maxlen,inumchar)
                    dtemp_bc_v_east=dblstr(directive,maxlen,inumchar)
                    dtemp_bc_w_east=dblstr(directive,maxlen,inumchar)
                  else
                    call warning(1,dble(iline),redstring)
                    lerror6=.true.
                  endif
                elseif(findstring('west',directive,inumchar,maxlen))then
                  if(findstring('type',directive,inumchar,maxlen))then
                    temp_bc_type_west=intstr(directive,maxlen,inumchar)
                  elseif(findstring('dens',directive,inumchar,maxlen))then
                    dtemp_bc_rhoR_west=dblstr(directive,maxlen,inumchar)
                    dtemp_bc_rhoB_west=dblstr(directive,maxlen,inumchar)
                  elseif(findstring('veloc',directive,inumchar,maxlen))then
                    dtemp_bc_u_west=dblstr(directive,maxlen,inumchar)
                    dtemp_bc_v_west=dblstr(directive,maxlen,inumchar)
                    dtemp_bc_w_west=dblstr(directive,maxlen,inumchar)
                  else
                    call warning(1,dble(iline),redstring)
                    lerror6=.true.
                  endif
                elseif(findstring('front',directive,inumchar,maxlen))then
                  if(findstring('type',directive,inumchar,maxlen))then
                    temp_bc_type_front=intstr(directive,maxlen,inumchar)
                  elseif(findstring('dens',directive,inumchar,maxlen))then
                    dtemp_bc_rhoR_front=dblstr(directive,maxlen,inumchar)
                    dtemp_bc_rhoB_front=dblstr(directive,maxlen,inumchar)
                  elseif(findstring('veloc',directive,inumchar,maxlen))then
                    dtemp_bc_u_front=dblstr(directive,maxlen,inumchar)
                    dtemp_bc_v_front=dblstr(directive,maxlen,inumchar)
                    dtemp_bc_w_front=dblstr(directive,maxlen,inumchar)
                  else
                    call warning(1,dble(iline),redstring)
                    lerror6=.true.
                  endif
                elseif(findstring('rear',directive,inumchar,maxlen))then
                  if(findstring('type',directive,inumchar,maxlen))then
                    temp_bc_type_rear=intstr(directive,maxlen,inumchar)
                  elseif(findstring('dens',directive,inumchar,maxlen))then
                    dtemp_bc_rhoR_rear=dblstr(directive,maxlen,inumchar)
                    dtemp_bc_rhoB_rear=dblstr(directive,maxlen,inumchar)
                  elseif(findstring('veloc',directive,inumchar,maxlen))then
                    dtemp_bc_u_rear=dblstr(directive,maxlen,inumchar)
                    dtemp_bc_v_rear=dblstr(directive,maxlen,inumchar)
                    dtemp_bc_w_rear=dblstr(directive,maxlen,inumchar)
                  else
                    call warning(1,dble(iline),redstring)
                    lerror6=.true.
                  endif
                elseif(findstring('north',directive,inumchar,maxlen))then
                  if(findstring('type',directive,inumchar,maxlen))then
                    temp_bc_type_north=intstr(directive,maxlen,inumchar)
                  elseif(findstring('dens',directive,inumchar,maxlen))then
                    dtemp_bc_rhoR_north=dblstr(directive,maxlen,inumchar)
                    dtemp_bc_rhoB_north=dblstr(directive,maxlen,inumchar)
                  elseif(findstring('veloc',directive,inumchar,maxlen))then
                    dtemp_bc_u_north=dblstr(directive,maxlen,inumchar)
                    dtemp_bc_v_north=dblstr(directive,maxlen,inumchar)
                    dtemp_bc_w_north=dblstr(directive,maxlen,inumchar)
                  else
                    call warning(1,dble(iline),redstring)
                    lerror6=.true.
                  endif
                elseif(findstring('south',directive,inumchar,maxlen))then
                  if(findstring('type',directive,inumchar,maxlen))then
                    temp_bc_type_south=intstr(directive,maxlen,inumchar)
                  elseif(findstring('dens',directive,inumchar,maxlen))then
                    dtemp_bc_rhoR_south=dblstr(directive,maxlen,inumchar)
                    dtemp_bc_rhoB_south=dblstr(directive,maxlen,inumchar)
                  elseif(findstring('veloc',directive,inumchar,maxlen))then
                    dtemp_bc_u_south=dblstr(directive,maxlen,inumchar)
                    dtemp_bc_v_south=dblstr(directive,maxlen,inumchar)
                    dtemp_bc_w_south=dblstr(directive,maxlen,inumchar)
                  else
                    call warning(1,dble(iline),redstring)
                    lerror6=.true.
                  endif
                else
                  call warning(1,dble(iline),redstring)
                  lerror6=.true.
                endif
              else
                call warning(1,dble(iline),redstring)
                lerror6=.true.
              endif
            elseif(findstring('dens',directive,inumchar,maxlen))then
              if(findstring('mean',directive,inumchar,maxlen))then
                dtemp_meanR=dblstr(directive,maxlen,inumchar)
                dtemp_meanB=dblstr(directive,maxlen,inumchar)
              elseif(findstring('stdev',directive,inumchar,maxlen))then
                dtemp_stdevR=dblstr(directive,maxlen,inumchar)
                dtemp_stdevB=dblstr(directive,maxlen,inumchar)
              elseif(findstring('back',directive,inumchar,maxlen))then
                dtemp_backR=dblstr(directive,maxlen,inumchar)
                dtemp_backB=dblstr(directive,maxlen,inumchar)
              elseif(findstring('selec',directive,inumchar,maxlen))then
                temp_areaR(1,1)=intstr(directive,maxlen,inumchar)
                temp_areaR(2,1)=intstr(directive,maxlen,inumchar)
                temp_areaR(1,2)=intstr(directive,maxlen,inumchar)
                temp_areaR(2,2)=intstr(directive,maxlen,inumchar)
                temp_areaR(1,3)=intstr(directive,maxlen,inumchar)
                temp_areaR(2,3)=intstr(directive,maxlen,inumchar)
              elseif(findstring('rescal',directive,inumchar,maxlen))then
                temp_imass_rescale=intstr(directive,maxlen,inumchar)
                temp_dmass_rescale=dblstr(directive,maxlen,inumchar)
                temp_lmass_rescale=.true.
              elseif(findstring('gauss',directive,inumchar,maxlen))then
                temp_idistselect=1
              elseif(findstring('unifo',directive,inumchar,maxlen))then
                temp_idistselect=2
              elseif(findstring('fake',directive,inumchar,maxlen))then
                temp_idistselect=3
              elseif(findstring('area',directive,inumchar,maxlen))then
                temp_idistselect=4
              else
                call warning(1,dble(iline),redstring)
                lerror6=.true.
              endif
            elseif(findstring('veloc',directive,inumchar,maxlen))then
              if(findstring('mean',directive,inumchar,maxlen))then
                dtemp_initial_u=dblstr(directive,maxlen,inumchar)
                dtemp_initial_v=dblstr(directive,maxlen,inumchar)
                dtemp_initial_w=dblstr(directive,maxlen,inumchar)
              else
                call warning(1,dble(iline),redstring)
                lerror6=.true.
              endif
            elseif(findstring('force',directive,inumchar,maxlen))then
              if(findstring('ext',directive,inumchar,maxlen))then
                dtemp_ext_fu=dblstr(directive,maxlen,inumchar)
                dtemp_ext_fv=dblstr(directive,maxlen,inumchar)
                dtemp_ext_fw=dblstr(directive,maxlen,inumchar)
              elseif(findstring('shanc',directive,inumchar,maxlen))then
                if(findstring('pair',directive,inumchar,maxlen))then
                  temp_lpair_SC=.true.
                  temp_pair_SC=dblstr(directive,maxlen,inumchar)
                elseif(findstring('wall',directive,inumchar,maxlen))then
                  temp_wall_SC = .true.
                  dtemp_wallR_SC = dblstr(directive,maxlen,inumchar)
                  dtemp_wallB_SC = dblstr(directive,maxlen,inumchar)
                else
                  call warning(1,dble(iline),redstring)
                  lerror6=.true.
                endif
              elseif(findstring('cap',directive,inumchar,maxlen))then
                dtemp_cap_force=dblstr(directive,maxlen,inumchar)
                temp_lcap_force=.true.
              else
                call warning(1,dble(iline),redstring)
                lerror6=.true.
              endif
            elseif(findstring('visc',directive,inumchar,maxlen))then
              lvisc=.true.
              dtemp_viscR=dblstr(directive,maxlen,inumchar)
              dtemp_viscB=dblstr(directive,maxlen,inumchar)
            elseif(findstring('tau',directive,inumchar,maxlen))then
              ltau=.true.
              dtemp_tauR=dblstr(directive,maxlen,inumchar)
              dtemp_tauB=dblstr(directive,maxlen,inumchar)
            elseif(findstring('wetting',directive,inumchar,maxlen))then
              if(findstring('mode',directive,inumchar,maxlen))then
                temp_lselwetting = .true.
                temp_iselwetting = intstr(directive,maxlen,inumchar)
                temp_densR_wetting = dblstr(directive,maxlen,inumchar)
                temp_densB_wetting = dblstr(directive,maxlen,inumchar)
              else
                call warning(1,dble(iline),redstring)
                lerror6=.true.
              endif
            elseif(findstring('[end room',directive,inumchar,maxlen))then
              lredo2=.false.
            elseif(findstring('[end',directive,inumchar,maxlen))then
              lredo=.false.
              lredo2=.false.
            else
              call warning(1,dble(iline),redstring)
              lerror6=.true.
            endif
          enddo
        elseif(findstring('md',directive,inumchar,maxlen))then
          lredo2=.true.
          lmd=.true.
          do while(lredo2)
            call getline(safe,inputunit,maxlen,redstring)
            if(.not.safe)then
             call warning(1,dble(iline),redstring)
             lerror5=.true.
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
            elseif(findstring('particle',directive,inumchar,maxlen))then
              if(findstring('type',directive,inumchar,maxlen))then
                temp_ntype=intstr(directive,maxlen,inumchar)
                if(temp_ntype<=0)cycle
                if(temp_ntype>mxntype)then
                  lerror5=.true.
                  lredo=.false.
                  call warning(39,dble(mxntype))
                endif
                temp_lparticlentype=.true.
                do i=1,temp_ntype
                  call getline(safe,inputunit,maxlen,redstring)
                  if(.not.safe)then
                    call warning(1,dble(iline),redstring)
                    lerror5=.true.
                    lredo=.false.
                    exit
                  endif
                  iline=iline+1
                  call strip(redstring,maxlen)
                  call copystring(redstring,directive,maxlen)
                  if(redstring(1:1)=='#')then
                    lerror7=.true.
                    lredo=.false.
                    call warning(37,dble(iline))
                    exit
                  elseif(redstring(1:1)=='!')then
                    lerror7=.true.
                    lredo=.false.
                    call warning(37,dble(iline))
                    exit
                  elseif(redstring(1:1)==' ')then
                    lerror7=.true.
                    lredo=.false.
                    call warning(37,dble(iline))
                    exit
                  else
                    call findwords(nwords,outwords,directive,maxlen)
                    if(nwords>0)then
                      directive=outwords(1)
                      temp_atmnamtype(i)=directive(1:8)
                    else
                      lerror7=.true.
                      lredo=.false.
                      call warning(37,dble(iline))
                      exit
                    endif
                  endif
                enddo
              elseif(findstring('yes',directive,inumchar,maxlen))then
                temp_lparticles=.true.
              elseif(findstring('no',directive,inumchar,maxlen))then
                temp_lparticles=.false.
              else
                call warning(1,dble(iline),redstring)
                lerror6=.true.
              endif
            elseif(findstring('shape',directive,inumchar,maxlen))then
              if(findstring('spher',directive,inumchar,maxlen))then
                itype=intstr(directive,maxlen,inumchar)
                if(itype==0)then
                  call warning(43)
                  call warning(1,dble(iline),redstring)
                  lerror5=.true.
                elseif(itype>mxntype)then
                  call warning(42,dble(mxntype))
                  call warning(1,dble(iline),redstring)
                  lerror5=.true.
                else
                  temp_ishape(itype)=0
                  temp_lurdim(itype)=.true.
                  dtemp_urdim(itype)=dblstr(directive,maxlen,inumchar)
                endif
                if(mod(dtemp_urdim(itype)-HALF,ONE)/=ZERO)then
                  call warning(30)
                  call warning(31,dble(iline),redstring)
                  lerror5=.true.
                endif
              else
                call warning(1,dble(iline),redstring)
                lerror6=.true.
              endif
            elseif(findstring('init',directive,inumchar,maxlen))then
              if(findstring('temperat',directive,inumchar,maxlen))then
                temp_linit_temp=.true.
                dtemp_init_temp=dblstr(directive,maxlen,inumchar)
              else
                call warning(1,dble(iline),redstring)
                lerror6=.true.
              endif
            elseif(findstring('densvar',directive,inumchar,maxlen))then
              dtemp_densvar=dblstr(directive,maxlen,inumchar)
              temp_densvar=.true.
            elseif(findstring('rcut',directive,inumchar,maxlen))then
              dtemp_rcut=dblstr(directive,maxlen,inumchar)
              temp_rcut=.true.
            elseif(findstring('delr',directive,inumchar,maxlen))then
              dtemp_delr=dblstr(directive,maxlen,inumchar)
              temp_delr=.true.
            elseif(findstring('rotat',directive,inumchar,maxlen))then
              if(findstring('yes',directive,inumchar,maxlen))then
                temp_lrotate=.true.
              elseif(findstring('no',directive,inumchar,maxlen))then
                temp_lrotate=.false.
              else
                call warning(1,dble(iline),redstring)
                lerror6=.true.
              endif
!            elseif(findstring('velocity',directive,inumchar,maxlen))then
!              if(findstring('verlet',directive,inumchar,maxlen))then
!                if(findstring('yes',directive,inumchar,maxlen))then
!                  temp_lvv=.true.
!                elseif(findstring('no',directive,inumchar,maxlen))then
!                  temp_lvv=.false.
!                else
!                  call warning(1,dble(iline),redstring)
!                  lerror6=.true.
!                endif
!              else
!                call warning(1,dble(iline),redstring)
!                lerror6=.true.
!              endif
            elseif(findstring('mass',directive,inumchar,maxlen))then
              itype=intstr(directive,maxlen,inumchar)
              if(itype<=0)then
                call warning(43)
                call warning(1,dble(iline),redstring)
                lerror5=.true.
              elseif(itype>mxntype)then
                call warning(42,dble(mxntype))
                call warning(1,dble(iline),redstring)
                lerror5=.true.
              else
                temp_lumass(itype)=.true.
                dtemp_umass(itype)=dblstr(directive,maxlen,inumchar)
              endif
            elseif(findstring('force',directive,inumchar,maxlen))then
              if(findstring('ext',directive,inumchar,maxlen))then
                dtemp_ext_fxx=dblstr(directive,maxlen,inumchar)
                dtemp_ext_fyy=dblstr(directive,maxlen,inumchar)
                dtemp_ext_fzz=dblstr(directive,maxlen,inumchar)
              elseif(findstring('shanc',directive,inumchar,maxlen))then
                if(findstring('angle',directive,inumchar,maxlen))then
                  temp_ltheta_SC=.true.
                  dtemp_theta_SC=dblstr(directive,maxlen,inumchar)
                  dtemp_devtheta_SC=dblstr(directive,maxlen,inumchar)
                elseif(findstring('part',directive,inumchar,maxlen))then
                  temp_lpart_SC = .true.
                  dtemp_partR_SC = dblstr(directive,maxlen,inumchar)
                  dtemp_partB_SC = dblstr(directive,maxlen,inumchar)
                else
                  call warning(1,dble(iline),redstring)
                  lerror6=.true.
                endif
              elseif(findstring('cap',directive,inumchar,maxlen))then
                dtemp_cap_force_part=dblstr(directive,maxlen,inumchar)
                temp_lcap_force_part=.true.
              else
                call warning(1,dble(iline),redstring)
                lerror6=.true.
              endif
            elseif(findstring('torque',directive,inumchar,maxlen))then
              if(findstring('ext',directive,inumchar,maxlen))then
                dtemp_ext_tqx=dblstr(directive,maxlen,inumchar)
                dtemp_ext_tqy=dblstr(directive,maxlen,inumchar)
                dtemp_ext_tqz=dblstr(directive,maxlen,inumchar)
              else
                call warning(1,dble(iline),redstring)
                lerror6=.true.
              endif
            elseif(findstring('side',directive,inumchar,maxlen))then
              if(findstring('wall',directive,inumchar,maxlen))then
                if(findstring('const',directive,inumchar,maxlen))then
                  temp_sidewall_md_k=.true.
                  dtemp_sidewall_md_k=dblstr(directive,maxlen,inumchar)
                elseif(findstring('dist',directive,inumchar,maxlen))then
                  temp_sidewall_md_rdist=.true.
                  dtemp_sidewall_md_rdist= &
                   dblstr(directive,maxlen,inumchar)
                else
                  call warning(1,dble(iline),redstring)
                  lerror6=.true.
                endif
              else
                call warning(1,dble(iline),redstring)
                lerror6=.true.
              endif
            elseif(findstring('field',directive,inumchar,maxlen))then
              if(findstring('pair',directive,inumchar,maxlen))then
                if(findstring('wca',directive,inumchar,maxlen))then
                  temp_field_pair=.true.
                  ifield_pair=ifield_pair+1
                  if(ifield_pair>mxvdw)then
                    call warning(24,dble(mxvdw))
                    call warning(1,dble(iline),redstring)
                    lerror6=.true.
                  else
                    temp_ltpvdw(ifield_pair)=1
                    itype=intstr(directive,maxlen,inumchar)
                    jtype=intstr(directive,maxlen,inumchar)
                    if(itype<=0 .or. jtype<=0 )then
                      call warning(46)
                      call warning(1,dble(iline),redstring)
                      lerror6=.true.
                    elseif(itype>mxntype .or. jtype>mxntype)then
                      call warning(42,dble(mxntype))
                      call warning(1,dble(iline),redstring)
                      lerror6=.true.
                    else
                      if(temp_mskvdw(itype,jtype)/=0 .or. &
                       temp_mskvdw(jtype,itype)/=0)then
                        call warning(45)
                        call warning(1,dble(iline),redstring)
                        lerror6=.true.
                      else
                        temp_mskvdw(itype,jtype)=ifield_pair
                        temp_mskvdw(jtype,itype)=ifield_pair
                        dtemp_prmvdw(1,ifield_pair)=dblstr(directive,maxlen,inumchar)
                        dtemp_prmvdw(2,ifield_pair)=dblstr(directive,maxlen,inumchar)
                      endif
                    endif
                  endif
                elseif(findstring('lj',directive,inumchar,maxlen))then
                  temp_field_pair=.true.
                  ifield_pair=ifield_pair+1
                  if(ifield_pair>mxvdw)then
                    call warning(24,dble(mxvdw))
                    call warning(1,dble(iline),redstring)
                    lerror6=.true.
                  else
                    temp_ltpvdw(ifield_pair)=2
                    itype=intstr(directive,maxlen,inumchar)
                    jtype=intstr(directive,maxlen,inumchar)
                    if(itype<=0 .or. jtype<=0 )then
                      call warning(46)
                      call warning(1,dble(iline),redstring)
                      lerror6=.true.
                    elseif(itype>mxntype .or. jtype>mxntype)then
                      call warning(42,dble(mxntype))
                      call warning(1,dble(iline),redstring)
                      lerror6=.true.
                    else
                      if(temp_mskvdw(itype,jtype)/=0 .or. &
                       temp_mskvdw(jtype,itype)/=0)then
                        call warning(45)
                        call warning(1,dble(iline),redstring)
                        lerror6=.true.
                      else
                        temp_mskvdw(itype,jtype)=ifield_pair
                        temp_mskvdw(jtype,itype)=ifield_pair
                        dtemp_prmvdw(1,ifield_pair)=dblstr(directive,maxlen,inumchar)
                        dtemp_prmvdw(2,ifield_pair)=dblstr(directive,maxlen,inumchar)
                      endif
                    endif
                  endif
                elseif(findstring('hz',directive,inumchar,maxlen))then
                  temp_field_pair=.true.
                  ifield_pair=ifield_pair+1
                  if(ifield_pair>mxvdw)then
                    call warning(24,dble(mxvdw))
                    call warning(1,dble(iline),redstring)
                    lerror6=.true.
                  else
                    temp_ltpvdw(ifield_pair)=3
                    itype=intstr(directive,maxlen,inumchar)
                    jtype=intstr(directive,maxlen,inumchar)
                    if(itype<=0 .or. jtype<=0 )then
                      call warning(46)
                      call warning(1,dble(iline),redstring)
                      lerror6=.true.
                    elseif(itype>mxntype .or. jtype>mxntype)then
                      call warning(42,dble(mxntype))
                      call warning(1,dble(iline),redstring)
                      lerror6=.true.
                    else
                      if(temp_mskvdw(itype,jtype)/=0 .or. &
                       temp_mskvdw(jtype,itype)/=0)then
                        call warning(45)
                        call warning(1,dble(iline),redstring)
                        lerror6=.true.
                      else
                        temp_mskvdw(itype,jtype)=ifield_pair
                        temp_mskvdw(jtype,itype)=ifield_pair
                        dtemp_prmvdw(1,ifield_pair)=dblstr(directive,maxlen,inumchar)
                        dtemp_prmvdw(2,ifield_pair)=dblstr(directive,maxlen,inumchar)
                        dtemp_prmvdw(3,ifield_pair)=dblstr(directive,maxlen,inumchar)
                      endif
                    endif
                  endif
                else
                  call warning(1,dble(iline),redstring)
                  lerror6=.true.
                endif
              else
                call warning(1,dble(iline),redstring)
                lerror6=.true.
              endif
            elseif(findstring('lubric',directive,inumchar,maxlen))then
              if(findstring('yes',directive,inumchar,maxlen))then
                temp_lubricconst=dblstr(directive,maxlen,inumchar)
                temp_lubricrparcut=dblstr(directive,maxlen,inumchar)
                temp_lubricrparcap=dblstr(directive,maxlen,inumchar)
                temp_llubrication=.true.
              elseif(findstring('no',directive,inumchar,maxlen))then
                temp_llubrication=.false.
              else
                call warning(1,dble(iline),redstring)
                lerror6=.true.
              endif
            elseif(findstring('moment',directive,inumchar,maxlen))then
              if(findstring('fix',directive,inumchar,maxlen))then
                temp_ifix_moment=intstr(directive,maxlen,inumchar)
                temp_lfix_moment=.true.
              else
                call warning(1,dble(iline),redstring)
                lerror6=.true.
              endif
            elseif(findstring('[end room',directive,inumchar,maxlen))then
              lredo2=.false.
            elseif(findstring('[end',directive,inumchar,maxlen))then
              lredo=.false.
              lredo2=.false.
            else
              call warning(1,dble(iline),redstring)
              lerror6=.true.
            endif
          enddo
        endif
      elseif(findstring('[end',directive,inumchar,maxlen))then
        lredo=.false.
      else
        call warning(1,dble(iline),redstring)
        lerror4=.true.
      endif
    enddo
    
    close(inputunit)
    
  endif
  
  call bcast_world_l(lerror1)
  call bcast_world_l(lerror2)
  call bcast_world_l(lerror4)
  call bcast_world_l(lerror5)
  call bcast_world_l(lerror6)
  call bcast_world_l(lerror7)
  
  if(lerror1)call error(1)
  if(lerror2)call error(2)
  if(lerror4)call error(4)
  if(lerror5)call error(5)
  if(lerror6)call error(6)
  if(lerror7)call error(7)
  
! send the read parameters to all the nodes and print them on terminal
  if(idrank==0)then
    write(6,'(/,3a,/)')repeat('*',26),"parameters from input file",repeat('*',27)
  endif
  
  if(idrank==0)then
    write(6,'(/,3a,/)')repeat('*',27),"parameters of system room",repeat('*',27)
  endif
  
  
  call bcast_world_l(lbox)
  if(lbox)then
    call bcast_world_i(temp_nx)
    call bcast_world_i(temp_ny)
    call bcast_world_i(temp_nz)
    call set_initial_dim_box(temp_nx,temp_ny,temp_nz)
    if(idrank==0)then
      mystring=repeat(' ',dimprint)
      mystring='system dimensions'
      write(6,'(2a,3i12)')mystring,": ",nx,ny,nz
    endif
  else
    call error(3)
  endif
  
  call bcast_world_l(temp_ldomdec)
  if(temp_ldomdec)then
    call bcast_world_i(temp_domdec)
    call set_domdec(temp_domdec)
    if(idrank==0)then
      mystring=repeat(' ',dimprint)
      mystring='decomposition type'
      write(6,'(2a,i12)')mystring,": ",domdec
    endif
  endif
  
  call bcast_world_l(temp_lnprocn)
  if(temp_lnprocn)then
    call bcast_world_i(temp_nprocx)
    call bcast_world_i(temp_nprocy)
    call bcast_world_i(temp_nprocz)
    call set_domain(temp_nprocx,temp_nprocy,temp_nprocz)
    if(idrank==0)then
      mystring=repeat(' ',dimprint)
      mystring='decomposition dimensions'
      write(6,'(2a,3i12)')mystring,": ",nprocx,nprocy,nprocz
    endif
  endif
  
  call bcast_world_l(lnstepmax)
  call bcast_world_i(temp_nstepmax)
  if(lnstepmax)then
    call set_nstepmax(temp_nstepmax)
    if(idrank==0)then
      mystring=repeat(' ',dimprint)
      mystring='steps'
      write(6,'(2a,i12)')mystring,": ",nstepmax
    endif
  else
    call warning(5)
    call error(7)
  endif
  
  call bcast_world_l(temp_ibc)
  call bcast_world_i(temp_ibcx)
  call bcast_world_i(temp_ibcy)
  call bcast_world_i(temp_ibcz)
  if(.not. temp_ibc)call error(10)
  call set_boundary_conditions_type(temp_ibcx,temp_ibcy,temp_ibcz)
  ! check boundary conditions if are supported/implemented
  ! 0 F F F
  ! 1 T F F
  ! 2 F T F
  ! 3 T T F
  ! 4 F F T
  ! 5 T F T
  ! 6 F T T
  ! 7 T T T
  if(ibctype<0 .or. ibctype>7)then
    call warning(3)
    call error(9)
  endif
  if(idrank==0)then
    mystring=repeat(' ',dimprint)
    mystring='boundary conditions'
    write(6,'(2a,3i12)')mystring,": ",temp_ibcx,temp_ibcy,temp_ibcz
  endif
  
  call bcast_world_i(nprintlist)
  call bcast_world_l(lprintlist)
  call bcast_world_l(lprintlisterror)
  if(lprintlist)then
    if(idrank/=0)allocate(printlist(nprintlist))
    call bcast_world_iarr(printlist,nprintlist)
    call label_argument(nprintlist,printlist,printarg)
    if(idrank==0)then
      mystring=repeat(' ',dimprint)
      mystring='explicit print list'
      mystring12=repeat(' ',dimprint2)
      mystring12='yes'
      mystring12=adjustr(mystring12)
      write(6,'(3a)')mystring,": ",mystring12
    endif
  endif
  if(lprintlisterror)then
    call warning(7)
    call warning(6)
    call error(7)
  endif
  call bcast_world_l(lprinttime)
  call bcast_world_f(printtime)
  if(lprinttime)then
   if(idrank==0)then
     mystring=repeat(' ',dimprint)
     mystring='print on terminal every'
     write(6,'(2a,i12)')mystring,": ",nint(printtime)
   endif
  endif
  
  call bcast_world_l(temp_lvtkfile)
  call bcast_world_l(temp_livtkfile)
  if(temp_lvtkfile)then
    if(.not. temp_livtkfile)then
      call warning(55)
      call error(7)
    endif
    call bcast_world_i(nprintlistvtk)
    call bcast_world_l(lprintlistvtk)
    call bcast_world_l(lprintlisterrorvtk)
    if(lprintlisterrorvtk)then
      call warning(57)
      call warning(56)
      call error(7)
    endif
    if(lprintlistvtk)then
      if(idrank/=0)allocate(printlistvtk(nprintlistvtk))
      call bcast_world_iarr(printlistvtk,nprintlistvtk)
      call bcast_world_i(temp_ivtkevery)
      call set_value_ivtkevery(temp_lvtkfile, &
       temp_ivtkevery)
      if(idrank==0)then
        mystring=repeat(' ',dimprint)
        mystring='print VTK file every'
        write(6,'(2a,i12)')mystring,": ",ivtkevery
        mystring=repeat(' ',dimprint)
        mystring='print VTK list'
        mystring12=repeat(' ',dimprint2)
        mystring12='yes'
        mystring12=adjustr(mystring12)
        write(6,'(3a)')mystring,": ",mystring12
        mystring=repeat(' ',dimprint)
        mystring='list specified in VTK file'
        mystring12=repeat(' ',dimprint2)
        write(6,'(3a)')mystring,": ",mystring12
        do i=1,nprintlistvtk
          write(6,'(36x,2a)')": ",legendobs_vtk(i,printlistvtk)
        enddo
      endif
    else
      call warning(58)
      call warning(56)
      call error(7)
    endif
  endif

  call bcast_world_i(temp_ibinevery)
  call bcast_world_l(temp_lbinary)
  if (temp_lbinary) call set_value_ibinevery(temp_lbinary,temp_ibinevery)
  if(idrank==0) then
     mystring=repeat(' ',dimprint)
     mystring='print binary every'
     write(6,'(2a,i12)')mystring,": ",temp_ibinevery
  endif

  call bcast_world_l(temp_lrestore)
  if (temp_lrestore) call set_restore(.true.)
  
  call bcast_world_l(temp_lidumpevery)
  if(temp_lidumpevery)then
    call bcast_world_i(temp_idumpevery)
    call set_value_dumpevery(temp_idumpevery)
    if(idrank==0)then
     mystring=repeat(' ',dimprint)
     mystring='print restart file every'
     write(6,'(2a,i12)')mystring,": ",idumpevery
   endif
  endif
  
  call bcast_world_l(temp_lxyzfile)
  if(temp_lxyzfile)then
    call bcast_world_i(temp_ixyzevery)
    call set_value_ixyzevery(temp_lxyzfile,temp_ixyzevery)
    if(idrank==0)then
     mystring=repeat(' ',dimprint)
     mystring='print xyz file every'
     write(6,'(2a,i12)')mystring,": ",ixyzevery
   endif
  endif
  
  call bcast_world_l(temp_ldiagnostic)
  if(temp_ldiagnostic)then
    call set_value_ldiagnostic(temp_ldiagnostic)
    if(idrank==0)then
      mystring=repeat(' ',dimprint)
      mystring='diagnostic profile'
      mystring12=repeat(' ',dimprint2)
      mystring12='yes'
      mystring12=adjustr(mystring12)
      write(6,'(3a)')mystring,": ",mystring12
    endif
  endif
  
  call bcast_world_l(lidiagnostic)
  if(lidiagnostic)then
    call bcast_world_i(temp_idiagnostic)
    call set_value_idiagnostic(temp_idiagnostic)
    if(idrank==0)then
      mystring=repeat(' ',dimprint)
      mystring='diagnostic profile every'
      write(6,'(2a,i12)')mystring,": ",idiagnostic
    endif
  endif
  
  call bcast_world_i(init_seed)
  call bcast_world_l(linit_seed)
  if(linit_seed)then
    if(idrank==0)then
      mystring=repeat(' ',dimprint)
      mystring='initial random seed'
      write(6,'(2a,i12)')mystring,": ",init_seed
    endif
  else
    if(idrank==0)then
      mystring=repeat(' ',dimprint)
      mystring='no seed mode'
      mystring12=repeat(' ',dimprint2)
      mystring12='yes'
      mystring12=adjustr(mystring12)
      write(6,'(3a)')mystring,": ",mystring12
    endif
  endif
  
  call bcast_world_l(ltest_mode)
  if(ltest_mode)then
    if(idrank==0)then
      mystring=repeat(' ',dimprint)
      mystring='test mode'
      mystring12=repeat(' ',dimprint2)
      mystring12='yes'
      mystring12=adjustr(mystring12)
      write(6,'(3a)')mystring,": ",mystring12
    endif
  endif
  
  call bcast_world_l(ltimjob)
  if(.not. ltimjob)then
    timjob=1.0d6*365.25d0*24.d0*60.d0*60.d0
  else
    call bcast_world_f(timjob)
    call get_prntime(hms,timjob,prntim)
    if(idrank==0)then
      mystring=repeat(' ',dimprint)
      mystring="user allocated job time ("//hms//") "
      write(6,'(2a,f12.6)')mystring,": ",prntim
      mystring=repeat(' ',dimprint)
      mystring='job closure time (s) '
      write(6,'(2a,f12.6)')mystring,": ",timcls
    endif
  endif
  
  if(idrank==0)then
    write(6,'(/,3a,/)')repeat('*',29),"parameters of LB room",repeat('*',29)
  endif
  
  call bcast_world_l(temp_lnfluid)
  if(temp_lnfluid)then
    call bcast_world_i(temp_nfluid)
    if(temp_nfluid/=1 .and. temp_nfluid/=2)then
      call warning(8,real(temp_nfluid,kind=PRC))
      call error(5)
    endif
    if(temp_nfluid==1)call set_lsingle_fluid(.true.)
    if(idrank==0)then
      mystring=repeat(' ',dimprint)
      mystring='fluid components'
      write(6,'(2a,i12)')mystring,": ",temp_nfluid
    endif
  endif
  
  call bcast_world_l(temp_lbc_halfway)
  call bcast_world_l(temp_lbc_fullway)
  if(temp_lbc_halfway.and.temp_lbc_fullway)then
    call warning(29)
    call error(5)
  endif
  
  if(temp_lbc_fullway)then
    call set_lbc_fullway(temp_lbc_fullway)
    if(idrank==0)then
      mystring=repeat(' ',dimprint)
      mystring='fullway bounce back'
      mystring12=repeat(' ',dimprint2)
      mystring12='yes'
      mystring12=adjustr(mystring12)
      write(6,'(3a)')mystring,": ",mystring12
    endif
  endif
  
  if(temp_lbc_halfway)then
    call set_lbc_halfway(temp_lbc_halfway)
    if(idrank==0)then
      mystring=repeat(' ',dimprint)
      mystring='halfway bounce back'
      mystring12=repeat(' ',dimprint2)
      mystring12='yes'
      mystring12=adjustr(mystring12)
      write(6,'(3a)')mystring,": ",mystring12
    endif
  endif
  
  call bcast_world_i(temp_bc_type_east)
  call bcast_world_i(temp_bc_type_west)
  call bcast_world_i(temp_bc_type_front)
  call bcast_world_i(temp_bc_type_rear)
  call bcast_world_i(temp_bc_type_north)
  call bcast_world_i(temp_bc_type_south)
  
  if(temp_bc_type_east/=0)then
    call bcast_world_f(dtemp_bc_rhoR_east)
    call bcast_world_f(dtemp_bc_rhoB_east)
    call bcast_world_f(dtemp_bc_u_east)
    call bcast_world_f(dtemp_bc_v_east)
    call bcast_world_f(dtemp_bc_w_east)
    call set_value_bc_east(temp_bc_type_east,dtemp_bc_rhoR_east, &
     dtemp_bc_rhoB_east,dtemp_bc_u_east,dtemp_bc_v_east,dtemp_bc_w_east)
    if(idrank==0)then
      mystring=repeat(' ',dimprint)
      mystring='open boundary east'
      write(6,'(2a,i12)')mystring,": ",bc_type_east
      mystring=repeat(' ',dimprint)
      mystring='open boundary east density'
      write(6,'(2a,2f12.6)')mystring,": ",bc_rhoR_east,bc_rhoB_east
      mystring=repeat(' ',dimprint)
      mystring='open boundary east vel'
      write(6,'(2a,3f12.6)')mystring,": ",bc_u_east,bc_v_east,bc_w_east
    endif
  endif
  
  if(temp_bc_type_west/=0)then
    call bcast_world_f(dtemp_bc_rhoR_west)
    call bcast_world_f(dtemp_bc_rhoB_west)
    call bcast_world_f(dtemp_bc_u_west)
    call bcast_world_f(dtemp_bc_v_west)
    call bcast_world_f(dtemp_bc_w_west)
    call set_value_bc_west(temp_bc_type_west,dtemp_bc_rhoR_west, &
     dtemp_bc_rhoB_west,dtemp_bc_u_west,dtemp_bc_v_west,dtemp_bc_w_west)
    if(idrank==0)then
      mystring=repeat(' ',dimprint)
      mystring='open boundary west'
      write(6,'(2a,i12)')mystring,": ",bc_type_west
      mystring=repeat(' ',dimprint)
      mystring='open boundary west density'
      write(6,'(2a,2f12.6)')mystring,": ",bc_rhoR_west,bc_rhoB_west
      mystring=repeat(' ',dimprint)
      mystring='open boundary west vel'
      write(6,'(2a,3f12.6)')mystring,": ",bc_u_west,bc_v_west,bc_w_west
    endif
  endif
  
  if(temp_bc_type_front/=0)then
    call bcast_world_f(dtemp_bc_rhoR_front)
    call bcast_world_f(dtemp_bc_rhoB_front)
    call bcast_world_f(dtemp_bc_u_front)
    call bcast_world_f(dtemp_bc_v_front)
    call bcast_world_f(dtemp_bc_w_front)
    call set_value_bc_front(temp_bc_type_front,dtemp_bc_rhoR_front, &
     dtemp_bc_rhoB_front,dtemp_bc_u_front,dtemp_bc_v_front,dtemp_bc_w_front)
    if(idrank==0)then
      mystring=repeat(' ',dimprint)
      mystring='open boundary front'
      write(6,'(2a,i12)')mystring,": ",bc_type_front
      mystring=repeat(' ',dimprint)
      mystring='open boundary front density'
      write(6,'(2a,2f12.6)')mystring,": ",bc_rhoR_front,bc_rhoB_front
      mystring=repeat(' ',dimprint)
      mystring='open boundary front vel'
      write(6,'(2a,3f12.6)')mystring,": ",bc_u_front,bc_v_front,bc_w_front
    endif
  endif
  
  if(temp_bc_type_rear/=0)then
    call bcast_world_f(dtemp_bc_rhoR_rear)
    call bcast_world_f(dtemp_bc_rhoB_rear)
    call bcast_world_f(dtemp_bc_u_rear)
    call bcast_world_f(dtemp_bc_v_rear)
    call bcast_world_f(dtemp_bc_w_rear)
    call set_value_bc_rear(temp_bc_type_rear,dtemp_bc_rhoR_rear, &
     dtemp_bc_rhoB_rear,dtemp_bc_u_rear,dtemp_bc_v_rear,dtemp_bc_w_rear)
    if(idrank==0)then
      mystring=repeat(' ',dimprint)
      mystring='open boundary rear'
      write(6,'(2a,i12)')mystring,": ",bc_type_rear
      mystring=repeat(' ',dimprint)
      mystring='open boundary rear density'
      write(6,'(2a,2f12.6)')mystring,": ",bc_rhoR_rear,bc_rhoB_rear
      mystring=repeat(' ',dimprint)
      mystring='open boundary rear vel'
      write(6,'(2a,3f12.6)')mystring,": ",bc_u_rear,bc_v_rear,bc_w_rear
    endif
  endif
  
  if(temp_bc_type_north/=0)then
    call bcast_world_f(dtemp_bc_rhoR_north)
    call bcast_world_f(dtemp_bc_rhoB_north)
    call bcast_world_f(dtemp_bc_u_north)
    call bcast_world_f(dtemp_bc_v_north)
    call bcast_world_f(dtemp_bc_w_north)
    call set_value_bc_north(temp_bc_type_north,dtemp_bc_rhoR_north, &
     dtemp_bc_rhoB_north,dtemp_bc_u_north,dtemp_bc_v_north,dtemp_bc_w_north)
    if(idrank==0)then
      mystring=repeat(' ',dimprint)
      mystring='open boundary north'
      write(6,'(2a,i12)')mystring,": ",bc_type_north
      mystring=repeat(' ',dimprint)
      mystring='open boundary north density'
      write(6,'(2a,2f12.6)')mystring,": ",bc_rhoR_north,bc_rhoB_north
      mystring=repeat(' ',dimprint)
      mystring='open boundary north vel'
      write(6,'(2a,3f12.6)')mystring,": ",bc_u_north,bc_v_north,bc_w_north
    endif
  endif
  
  if(temp_bc_type_south/=0)then
    call bcast_world_f(dtemp_bc_rhoR_south)
    call bcast_world_f(dtemp_bc_rhoB_south)
    call bcast_world_f(dtemp_bc_u_south)
    call bcast_world_f(dtemp_bc_v_south)
    call bcast_world_f(dtemp_bc_w_south)
    call set_value_bc_south(temp_bc_type_south,dtemp_bc_rhoR_south, &
     dtemp_bc_rhoB_south,dtemp_bc_u_south,dtemp_bc_v_south,dtemp_bc_w_south)
    if(idrank==0)then
      mystring=repeat(' ',dimprint)
      mystring='open boundary south'
      write(6,'(2a,i12)')mystring,": ",bc_type_south
      mystring=repeat(' ',dimprint)
      mystring='open boundary south density'
      write(6,'(2a,2f12.6)')mystring,": ",bc_rhoR_south,bc_rhoB_south
      mystring=repeat(' ',dimprint)
      mystring='open boundary south vel'
      write(6,'(2a,3f12.6)')mystring,": ",bc_u_south,bc_v_south,bc_w_south
    endif
  endif
  
  call bcast_world_l(lvisc)
  call bcast_world_l(ltau)
  if(.not. (lvisc.or.ltau))then
    call warning(4)
    call error(7)
  endif
  if(lvisc)then
    call bcast_world_f(dtemp_viscR)
    call bcast_world_f(dtemp_viscB)
    call set_value_viscosity(dtemp_viscR,dtemp_viscB,lsingle_fluid)
  else
    call bcast_world_f(dtemp_tauR)
    call bcast_world_f(dtemp_tauB)
    call set_value_tau(dtemp_tauR,dtemp_tauB,lsingle_fluid)
  endif
  
  if(idrank==0)then
    mystring=repeat(' ',dimprint)
    mystring='fluid viscosity'
    write(6,'(2a,2f12.6)')mystring,": ",viscR,viscB
    mystring=repeat(' ',dimprint)
    mystring='fluid tau'
    write(6,'(2a,2f12.6)')mystring,": ",tauR,tauB
    if(lunique_omega)then
      mystring=repeat(' ',dimprint)
      mystring='constant omega mode'
      mystring12=repeat(' ',dimprint2)
      mystring12='yes'
      mystring12=adjustr(mystring12)
      write(6,'(3a)')mystring,": ",mystring12
    endif
  endif
  
  call bcast_world_i(temp_idistselect)
  if(temp_idistselect>0)then
    call set_initial_dist_type(temp_idistselect)
    if(idrank==0)then
      select case(idistselect)
      case(1)
        mystring=repeat(' ',dimprint)
        mystring='initial density distribution'
        mystring12=repeat(' ',dimprint2)
        mystring12='gaussian'
        mystring12=adjustr(mystring12)
        write(6,'(3a)')mystring,": ",mystring12
      case(2)
        mystring=repeat(' ',dimprint)
        mystring='initial density distribution'
        mystring12=repeat(' ',dimprint2)
        mystring12='uniform'
        mystring12=adjustr(mystring12)
        write(6,'(3a)')mystring,": ",mystring12
      case(3)
        mystring=repeat(' ',dimprint)
        mystring='initial density distribution'
        mystring12=repeat(' ',dimprint2)
        mystring12='fake from id fluid node'
        mystring12=adjustr(mystring12)
        write(6,'(3a)')mystring,": ",mystring12
      case(4)
        mystring=repeat(' ',dimprint)
        mystring='initial density distribution'
        mystring12=repeat(' ',dimprint2)
        mystring12='uniform areas'
        mystring12=adjustr(mystring12)
        write(6,'(3a)')mystring,": ",mystring12
      end select
    endif
  endif
  
  call bcast_world_f(dtemp_meanR)
  call bcast_world_f(dtemp_meanB)
  if(dtemp_meanR>ZERO .or. dtemp_meanB>ZERO)then
    call set_mean_value_dens_fluids(dtemp_meanR,dtemp_meanB)
    if(idrank==0)then
      mystring=repeat(' ',dimprint)
      mystring='initial density mean values'
      write(6,'(2a,2f12.6)')mystring,": ",meanR,meanB
    endif
  endif
  
  call bcast_world_f(dtemp_stdevR)
  call bcast_world_f(dtemp_stdevB)
  if(dtemp_stdevR>ZERO .or. dtemp_stdevB>ZERO)then
    call set_stdev_value_dens_fluids(dtemp_stdevR,dtemp_stdevB)
    if(idrank==0)then
      mystring=repeat(' ',dimprint)
      mystring='initial density stdev values'
      write(6,'(2a,2f12.6)')mystring,": ",stdevR,stdevB
    endif
  endif
  
  if(idistselect==4)then
    call bcast_world_iarr(temp_areaR,6)
    if(any(temp_areaR/=0))then
      call bcast_world_f(dtemp_backR)
      call bcast_world_f(dtemp_backB)
      if(lsingle_fluid)then
        if(dtemp_backR<=ZERO)then
          call warning(52)
          call error(7)
        endif
      else
        if(dtemp_backR<=ZERO .or. dtemp_backB<=ZERO)then
          call warning(52)
          call error(7)
        endif
      endif
      call set_back_value_dens_fluids(dtemp_backR,dtemp_backB)
      if(idrank==0)then
        mystring=repeat(' ',dimprint)
        mystring='initial background density values'
        write(6,'(2a,2f12.6)')mystring,": ",backR,backB
      endif
      call set_init_area_dens_fluids(temp_areaR)
      if(idrank==0)then
        mystring=repeat(' ',dimprint)
        mystring='initial density area extremes in x'
        write(6,'(2a,2i12)')mystring,": ",areaR(1,1),areaR(2,1)
        mystring='initial density area extremes in y'
        write(6,'(2a,2i12)')mystring,": ",areaR(1,2),areaR(2,2)
        mystring='initial density area extremes in z'
        write(6,'(2a,2i12)')mystring,": ",areaR(1,3),areaR(2,3)
      endif
    endif
  endif
  
  call bcast_world_l(temp_lmass_rescale)
  if(temp_lmass_rescale)then
    call bcast_world_i(temp_imass_rescale)
    !remove the percentage
    temp_dmass_rescale=temp_dmass_rescale/real(100.d0,kind=PRC)
    call bcast_world_f(temp_dmass_rescale)
    call set_mass_rescale(temp_lmass_rescale,temp_imass_rescale, &
     temp_dmass_rescale)
    if(idrank==0)then
      mystring=repeat(' ',dimprint)
      mystring='fluid mass rescaled every'
      write(6,'(2a,i12)')mystring,": ",imass_rescale
      mystring=repeat(' ',dimprint)
      mystring='fluid mass rescaled tolerance %'
      write(6,'(2a,f12.6)')mystring,": ", &
       dmass_rescale*real(100.d0,kind=PRC)
    endif
  endif
  
  call bcast_world_f(dtemp_initial_u)
  call bcast_world_f(dtemp_initial_v)
  call bcast_world_f(dtemp_initial_w)
  if(dtemp_initial_u>ZERO .or. dtemp_initial_v>ZERO .or. &
   dtemp_initial_w>ZERO)then
    call set_mean_value_vel_fluids(dtemp_initial_u,dtemp_initial_v, &
     dtemp_initial_w)
    if(idrank==0)then
      mystring=repeat(' ',dimprint)
      mystring='initial velocity mean values'
      write(6,'(2a,3f12.6)')mystring,": ",initial_u,initial_v,initial_w
    endif
  endif
  
  call bcast_world_f(dtemp_ext_fu)
  call bcast_world_f(dtemp_ext_fv)
  call bcast_world_f(dtemp_ext_fw)
  if(dtemp_ext_fu/=ZERO .or. dtemp_ext_fv/=ZERO .or. &
   dtemp_ext_fw/=ZERO)then
    call set_value_ext_force_fluids(dtemp_ext_fu,dtemp_ext_fv, &
     dtemp_ext_fw)
    if(idrank==0)then
      mystring=repeat(' ',dimprint)
      mystring='external force on fluids'
      write(6,'(2a,3f12.6)')mystring,": ",ext_fu,ext_fv,ext_fw
    endif
  endif
  
  call bcast_world_l(temp_lpair_SC)
  call bcast_world_f(temp_pair_SC)
  if(temp_lpair_SC)then
    if(temp_lpair_SC .and. lsingle_fluid)then
      call warning(10)
      call error(5)
    endif
    call set_fluid_force_sc(temp_lpair_SC,temp_pair_SC)
    if(idrank==0)then
      mystring=repeat(' ',dimprint)
      mystring='pair SC force constant on fluids'
      write(6,'(2a,f12.6)')mystring,": ",pair_SC
    endif
  endif
  
  call bcast_world_l(temp_wall_SC)
  if(temp_lpair_SC)then
    call bcast_world_f(dtemp_wallR_SC)
    call bcast_world_f(dtemp_wallB_SC)
    call set_fluid_wall_sc(dtemp_wallR_SC,dtemp_wallB_SC)
    if(idrank==0)then
      mystring=repeat(' ',dimprint)
      mystring='pair SC wall constant on fluids'
      write(6,'(2a,2f12.6)')mystring,": ",wallR_SC,wallB_SC
    endif
  endif
  
  if(lforce_add)then
    call set_LBintegrator_type(1)
    if(idrank==0)then
      mystring=repeat(' ',dimprint)
      mystring='extra force term'
      mystring12=repeat(' ',dimprint2)
      mystring12='yes'
      mystring12=adjustr(mystring12)
      write(6,'(3a)')mystring,": ",mystring12
    endif
  else
    call set_LBintegrator_type(0)
    if(idrank==0)then
      mystring=repeat(' ',dimprint)
      mystring='extra force term'
      mystring12=repeat(' ',dimprint2)
      mystring12='no'
      mystring12=adjustr(mystring12)
      write(6,'(3a)')mystring,": ",mystring12
    endif
  endif
  
  
  call bcast_world_l(temp_lcap_force)
  if(temp_lcap_force)then
    call bcast_world_f(dtemp_cap_force)
    call set_cap_force(temp_lcap_force,dtemp_cap_force)
    if(idrank==0)then
      mystring=repeat(' ',dimprint)
      mystring='capping force on fluids'
      write(6,'(2a,f12.6)')mystring,": ",cap_force
    endif
  endif
  
  call bcast_world_l(temp_lselwetting)
  if(temp_lselwetting)then
    call bcast_world_i(temp_iselwetting)
    call bcast_world_f(temp_densR_wetting)
    call bcast_world_f(temp_densB_wetting)
    call set_fluid_wall_mode(temp_iselwetting,temp_densR_wetting, &
     temp_densB_wetting)
    if(idrank==0)then
      select case(iselwetting)
      case(1)
        mystring=repeat(' ',dimprint)
        mystring='wetting wall mode'
        mystring12=repeat(' ',dimprint2)
        mystring12='fixed density'
        mystring12=adjustr(mystring12)
        write(6,'(3a)')mystring,": ",mystring12
        mystring=repeat(' ',dimprint)
        mystring='fixed fluid density'
        write(6,'(2a,2f12.6)')mystring,": ",densR_wetting,densB_wetting
      case default
        mystring=repeat(' ',dimprint)
        mystring='wetting wall mode'
        mystring12=repeat(' ',dimprint2)
        mystring12='averaged'
        mystring12=adjustr(mystring12)
        write(6,'(3a)')mystring,": ",mystring12
      end select
    endif
  endif
  
  if(LBintegrator<0 .or. LBintegrator>1)call error(14)
  
  if(idrank==0)then
    mystring=repeat(' ',dimprint)
    mystring='LB integrator'
    mystring12=repeat(' ',dimprint2)
    select case(LBintegrator)
    case(0)
      mystring12='BGK'
      mystring12=adjustr(mystring12)
      write(6,'(3a)')mystring,": ",mystring12
    case(1)
      mystring12='EDM'
      mystring12=adjustr(mystring12)
      write(6,'(3a)')mystring,": ",mystring12
    end select
  endif
  
  call bcast_world_l(lmd)
  if(lmd)then
    if(idrank==0)write(6,'(/,3a,/)')repeat('*',29), &
     "parameters of MD room",repeat('*',29)
    
    call bcast_world_l(temp_lparticles)
    call set_lparticles(temp_lparticles)
    
    call bcast_world_l(temp_lrotate)
    call set_lrotate(temp_lrotate)
    
    if(lparticles)call allocate_particle_features
    
    if(idrank==0)then
      mystring=repeat(' ',dimprint)
      mystring='particles'
      mystring12=repeat(' ',dimprint2)
      if(lparticles)then
        mystring12='yes'
      else
        mystring12='no'
      endif
      mystring12=adjustr(mystring12)
      write(6,'(3a)')mystring,": ",mystring12
    endif
    
    
    if(lrotate)then
      if(idrank==0 .and. lparticles)then
        mystring=repeat(' ',dimprint)
        mystring='rotation'
        mystring12=repeat(' ',dimprint2)
        if(lrotate)then
          mystring12='yes'
        else
          mystring12='no'
        endif
        mystring12=adjustr(mystring12)
        write(6,'(3a)')mystring,": ",mystring12
      endif
    endif
    
    call bcast_world_l(temp_lparticlentype)
    if(lparticles)then
      if(.not. temp_lparticlentype)then
        call warning(38)
        call error(7)
      endif
      call bcast_world_i(temp_ntype)
      call set_ntype(temp_ntype)
      call bcast_world_carr(8,temp_atmnamtype,mxntype)
      call set_atmnamtype(temp_atmnamtype)
      if(idrank==0 .and. lparticles)then
        mystring=repeat(' ',dimprint)
        mystring='particle types'
        write(6,'(2a,i12)')mystring,": ",ntype
        do itype=1,ntype
          mystring=repeat(' ',dimprint)
          write (r_char,'(i1)')itype
          mystring=repeat(' ',4)//r_char//'°type label'
          mystring12=repeat(' ',12)
          mystring12=atmnamtype(itype)//repeat(' ',4)
          mystring12=trim(adjustr(mystring12))
          write(6,'(3a)')mystring," : ",mystring12
        enddo
      endif
    endif
    
    if(lparticles)then
      call bcast_world_iarr(temp_ishape,mxntype)
      if(any(temp_ishape>=0))then
        call set_ishape(temp_ishape)
        if(idrank==0)then
          mystring=repeat(' ',dimprint)
          mystring='particle shape'
          write(6,'(2a,i12)')mystring,": ",ntype
          do itype=1,ntype
            mystring=repeat(' ',dimprint)
            write (r_char,'(i1)')itype
            mystring=repeat(' ',4)//r_char//'°type shape '
            select case(ishape(itype))
            case(0)
              mystring12='spherical'
            case(1)
              mystring12='ellipsoid'
            end select
            mystring12=adjustr(mystring12)
            write(6,'(3a)')mystring," : ",mystring12
          enddo
        endif
      else
        call warning(34)
        call error(7)
      endif
    endif
    
    
    
    call bcast_world_larr(temp_lurdim,mxntype)
    if(.not. all(temp_lurdim(1:ntype)))then
      call warning(44)
      call error(7)
    endif
    if(any(temp_lurdim) .and. lparticles)then
      call bcast_world_farr(dtemp_urdim,mxntype)
      call bcast_world_farr(dtemp_urdimy,mxntype)
      call bcast_world_farr(dtemp_urdimz,mxntype)
      call set_rdim(dtemp_urdim,dtemp_urdimy,dtemp_urdimz)
      if(idrank==0 .and. lparticles)then
        mystring=repeat(' ',dimprint)
        mystring='particle radius'
        write(6,'(2a,i12)')mystring,": ",ntype
        do itype=1,ntype
          mystring=repeat(' ',dimprint)
          write (r_char,'(i1)')itype
          mystring=repeat(' ',4)//r_char//'°type radius '
          select case(ishape(itype))
          case(0)
            write(6,'(2a,f12.6)')mystring," : ",rdimx(itype)
          case(1)
            write(6,'(2a,3f12.6)')mystring," : ", &
             rdimx(itype),rdimy(itype),rdimz(itype)
          end select
        enddo
      endif
    endif
    
    
    call bcast_world_l(temp_linit_temp)
    if(temp_linit_temp)then
      call bcast_world_f(dtemp_init_temp)
      call set_init_temp(dtemp_init_temp)
      if(idrank==0)then
        mystring=repeat(' ',dimprint)
        mystring='initial temperature'
        write(6,'(2a,f12.6)')mystring,": ",init_temp
      endif
    endif
    
    call bcast_world_l(temp_lvv)
    if(temp_lvv)then
      call set_lvv(temp_lvv)
      if(idrank==0 .and. lparticles)then
        mystring=repeat(' ',dimprint)
        mystring='velocity Verlet'
        mystring12=repeat(' ',dimprint2)
        if(lvv)then
          mystring12='yes'
        else
          mystring12='no'
        endif
        mystring12=adjustr(mystring12)
        write(6,'(3a)')mystring,": ",mystring12
      endif
    endif
    
    call bcast_world_f(dtemp_ext_fxx)
    call bcast_world_f(dtemp_ext_fyy)
    call bcast_world_f(dtemp_ext_fzz)
    if(dtemp_ext_fxx/=ZERO .or. dtemp_ext_fyy/=ZERO .or. &
     dtemp_ext_fzz/=ZERO)then
      call set_value_ext_force_particles(dtemp_ext_fxx,dtemp_ext_fyy, &
       dtemp_ext_fzz)
      if(idrank==0)then
        mystring=repeat(' ',dimprint)
        mystring='external force on particles'
        write(6,'(2a,3f12.6)')mystring,": ",ext_fxx,ext_fyy,ext_fzz
      endif
    endif
    
    if(lrotate)then
      call bcast_world_f(dtemp_ext_tqx)
      call bcast_world_f(dtemp_ext_tqy)
      call bcast_world_f(dtemp_ext_tqz)
      if(dtemp_ext_tqx/=ZERO .or. dtemp_ext_tqy/=ZERO .or. &
       dtemp_ext_tqz/=ZERO)then
        call set_value_ext_torque_particles(dtemp_ext_tqx, &
         dtemp_ext_tqy,dtemp_ext_tqz)
        if(idrank==0)then
          mystring=repeat(' ',dimprint)
          mystring='external torque on particles'
          write(6,'(2a,3f12.6)')mystring,": ",ext_tqx,ext_tqy,ext_tqz
        endif
      endif
    endif
    
    call bcast_world_l(temp_lpart_SC)
    call bcast_world_f(dtemp_partR_SC)
    call bcast_world_f(dtemp_partB_SC)
    if(temp_lpart_SC)then
      call set_fluid_particle_sc(dtemp_partR_SC,dtemp_partB_SC)
      if(idrank==0)then
        mystring=repeat(' ',dimprint)
        mystring='pair SC color particle on fluids'
        write(6,'(2a,2f12.6)')mystring,": ",partR_SC,partB_SC
      endif
    endif
    
    if(lrotate)then
      call bcast_world_l(temp_ltheta_SC)
      call bcast_world_f(dtemp_theta_SC)
      call bcast_world_f(dtemp_devtheta_SC)
      if(temp_lpart_SC)then
        dtemp_theta_SC=dtemp_theta_SC*conv_rad
        dtemp_devtheta_SC=dtemp_devtheta_SC*conv_rad
        call set_theta_particle_sc(dtemp_theta_SC,dtemp_devtheta_SC)
        if(idrank==0)then
          mystring=repeat(' ',dimprint)
          mystring='pair SC contact angle of particles'
          write(6,'(2a,2f12.6)')mystring,": ",theta_SC/conv_rad, &
           devtheta_SC/conv_rad
        endif
      endif
    endif
    
    call bcast_world_l(temp_lcap_force_part)
    if(temp_lcap_force_part)then
      call bcast_world_f(dtemp_cap_force_part)
      call set_cap_force_part(temp_lcap_force_part,dtemp_cap_force_part)
      if(idrank==0)then
        mystring=repeat(' ',dimprint)
        mystring='capping force on particles'
        write(6,'(2a,f12.6)')mystring,": ",cap_force_part
      endif
    endif
    
    call bcast_world_l(temp_lfix_moment)
    if(temp_lfix_moment)then
      call bcast_world_i(temp_ifix_moment)
      call set_fix_moment(temp_lfix_moment,temp_ifix_moment)
      if(idrank==0)then
        mystring=repeat(' ',dimprint)
        mystring='zero CoM moment on particles every'
        write(6,'(2a,i12)')mystring,": ",ifix_moment
      endif
    endif
    
    call bcast_world_l(temp_rcut)
    if(temp_rcut)then
      call bcast_world_f(dtemp_rcut)
      call set_rcut(dtemp_rcut)
      if(idrank==0 .and. lparticles)then
        mystring=repeat(' ',dimprint)
        mystring='potential cut-off'
        write(6,'(2a,f12.6)')mystring,": ",rcut
      endif
    else
      call warning(17)
      call error(7)
    endif
    
    call bcast_world_larr(temp_lumass,mxntype)
    if(.not. all(temp_lumass(1:ntype)))then
      call warning(40)
      call error(7)
    endif
    if(any(temp_lumass))then
      call bcast_world_farr(dtemp_umass,mxntype)
      call set_umass(dtemp_umass)
      if(idrank==0 .and. lparticles)then
        mystring=repeat(' ',dimprint)
        mystring='particle mass'
        write(6,'(2a,i12)')mystring,": ",ntype
        do itype=1,ntype
          mystring=repeat(' ',dimprint)
          write (r_char,'(i1)')itype
          mystring=repeat(' ',4)//r_char//'°type radius '
          write(6,'(2a,f12.6)')mystring," : ",weight(itype)
        enddo
      endif
    else
      if(lparticles)then
        call warning(35)
        call error(7)
      endif
    endif
    
    if(lparticles .and. ibctype/=7)then
      temp_lsidewall_md=.true.
      call set_lsidewall_md(temp_lsidewall_md)
      if(idrank==0)then
        mystring=repeat(' ',dimprint)
        mystring='side wall'
        mystring12=repeat(' ',dimprint2)
        if(lsidewall_md)then
          mystring12='yes'
        else
          mystring12='no'
        endif
        mystring12=adjustr(mystring12)
        write(6,'(3a)')mystring,": ",mystring12
      endif
      call bcast_world_l(temp_sidewall_md_k)
      if(temp_sidewall_md_k)then
        call bcast_world_f(dtemp_sidewall_md_k)
        if(dtemp_sidewall_md_k>ZERO)then
          call set_sidewall_md_k(dtemp_sidewall_md_k)
          if(idrank==0)then
            mystring=repeat(' ',dimprint)
            mystring='side wall constant'
            write(6,'(2a,f12.6)')mystring,": ",sidewall_md_k
          endif
        else
          call warning(48)
          call error(7)
        endif
      else
        call warning(47)
        call error(7)
      endif
      call bcast_world_l(temp_sidewall_md_rdist)
      if(temp_sidewall_md_rdist)then
        call bcast_world_f(dtemp_sidewall_md_rdist)
        if(dtemp_sidewall_md_rdist>ZERO)then
          call set_sidewall_md_rdist(dtemp_sidewall_md_rdist)
          if(idrank==0)then
            mystring=repeat(' ',dimprint)
            mystring='side wall distance'
            write(6,'(2a,f12.6)')mystring,": ",sidewall_md_rdist
          endif
        else
          call warning(50)
          call error(7)
        endif
      else
        call warning(51)
        call error(7)
      endif
    endif
    
    call bcast_world_l(temp_field_pair)
    if(temp_field_pair)then
      call bcast_world_i(ifield_pair)
      call set_ntpvdw(ifield_pair)
      call allocate_field_array
      call bcast_world_iarr(temp_ltpvdw,mxvdw)
      call bcast_world_iarr(temp_mskvdw,mxntype*mxntype)
      call bcast_world_farr(dtemp_prmvdw,mxpvdw*mxvdw)
      call set_field_array(ntpvdw,mxpvdw,temp_ltpvdw,dtemp_prmvdw, &
       temp_mskvdw)
       
      if(idrank==0 .and. lparticles)then
        mystring=repeat(' ',dimprint)
        mystring='number of field pairs'
        write(6,'(2a,i12)')mystring,": ",ntpvdw
        mystring=repeat(' ',dimprint)
        mystring='list of force fields'
        write(6,'(2a)')mystring,": "
        do i=1,ntpvdw
          mystring=repeat(' ',dimprint)
          write(mystring,'(i2)')i
          mystring=adjustr(mystring)
          mystring=mystring(35:36)//'° type of force field'
          mystring36=repeat(' ',dimprint)
          select case(ltpvdw(i))
          case(1)
            mystring36='Weeks-Chandler-Andersen'
          case(2)
            mystring36='Lennard-Jones'
          case(3)
            mystring36='Hertzian'
          end select
          mystring36=adjustl(mystring36)
          write(6,'(3a)')mystring," : ",mystring36
          mystring=repeat(' ',dimprint)
          mystring='    itype jtype pair'
          do j=1,mxntype
            do k=1,mxntype
              if(mskvdw(j,k)==i)then
                itype=j
                jtype=k
              endif
            enddo
          enddo
          write(6,'(2a,2a12)')mystring,": ",atmnamtype(itype), &
           atmnamtype(jtype)
          mystring=repeat(' ',dimprint)
          mystring='    kappa'
          write(6,'(2a,f12.6)')mystring,": ",prmvdw(1,i)
          mystring=repeat(' ',dimprint)
          mystring='    r min'
          write(6,'(2a,f12.6)')mystring,": ",prmvdw(2,i)
          mystring=repeat(' ',dimprint)
          mystring='    r cap'
          write(6,'(2a,f12.6)')mystring,": ",prmvdw(3,i)
        enddo
      endif
    endif
    
    call bcast_world_l(temp_llubrication)
    if(temp_llubrication)then
      call bcast_world_f(temp_lubricconst)
      call bcast_world_f(temp_lubricrparcut)
      call bcast_world_f(temp_lubricrparcap)
      call set_lubrication(temp_llubrication,temp_lubricrparcut, &
       temp_lubricrparcap,temp_lubricconst)
      if(idrank==0 .and. llubrication)then
        mystring=repeat(' ',dimprint)
        mystring='lubrication'
        mystring12=repeat(' ',dimprint2)
        if(llubrication)then
          mystring12='yes'
        else
          mystring12='no'
        endif
        mystring12=adjustr(mystring12)
        write(6,'(3a)')mystring,": ",mystring12
        mystring=repeat(' ',dimprint)
        mystring='    kappa'
        write(6,'(2a,f12.6)')mystring,": ",lubricconst
        mystring=repeat(' ',dimprint)
        mystring='    r min'
        write(6,'(2a,f12.6)')mystring,": ",lubricrparcut
        mystring=repeat(' ',dimprint)
        mystring='    r cap'
        write(6,'(2a,f12.6)')mystring,": ",lubricrparcap
      endif
    endif
    
    call bcast_world_l(temp_delr)
    if(temp_delr)then
      call bcast_world_f(dtemp_delr)
      call set_delr(dtemp_delr)
      if(idrank==0)then
        mystring=repeat(' ',dimprint)
        mystring='width of verlet border'
        write(6,'(2a,f12.6)')mystring,": ",delr
      endif
    else
      call warning(18)
      call error(7)
    endif
    
    call bcast_world_l(temp_densvar)
    if(temp_densvar)then
      call bcast_world_f(dtemp_densvar)
      if(dtemp_densvar<ONE)then
        call warning(14,dtemp_densvar)
        call error(5)
      endif
      call set_densvar(dtemp_densvar)
      if(idrank==0)then
        mystring=repeat(' ',dimprint)
        mystring='density variance (densvar)'
        write(6,'(2a,f12.6)')mystring,": ",densvar
      endif
    endif
  endif
  
  if(idrank==0)then
    write(6,'(/,3a,/)')repeat('*',32),"end input file",repeat('*',33)
  endif
  
  endtime = tstep*REAL(nstepmax,kind=PRC)
  iprinttime=nint(printtime/tstep)
  eprinttime=nint(endtime/tstep)
  
  reprinttime=nint(REAL(eprinttime,kind=PRC)/REAL(iprinttime,kind=PRC))
  
  call get_sync_world
  
  call print_legend_observables(6)

! print warning and check error
  
  return
  
 end subroutine read_input
 
 function legendobs(iarg,printcodsub)
 
!***********************************************************************
!     
!     LBsoft function for returning the legend which is associated to
!     the integer contained in the printcod array
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: iarg
  integer, intent(in), allocatable, dimension(:) :: printcodsub
  
  character(len=52) :: legendobs
  
  legendobs=repeat(' ',52)
  if(printcodsub(iarg)==1)then
    legendobs='dens1 =  mean fluid density of first component    '
  elseif(printcodsub(iarg)==2)then
    legendobs='dens2 =  mean fluid density of second component   '
  elseif(printcodsub(iarg)==3)then
    legendobs='maxd1 =  max fluid density of first component     '
  elseif(printcodsub(iarg)==4)then
    legendobs='maxd2 =  max fluid density of second component    '
  elseif(printcodsub(iarg)==5)then
    legendobs='mind1 =  min fluid density of first component     '
  elseif(printcodsub(iarg)==6)then
    legendobs='mind2 =  min fluid density of second component    '
  elseif(printcodsub(iarg)==7)then
    legendobs='maxvx =  max fluid velocity of fluid mix along x  '
  elseif(printcodsub(iarg)==8)then
    legendobs='minvx =  min fluid velocity of fluid mix along x  '
  elseif(printcodsub(iarg)==9)then
    legendobs='maxvy =  max fluid velocity of fluid mix along y  '
  elseif(printcodsub(iarg)==10)then
    legendobs='minvy =  min fluid velocity of fluid mix along y  '
  elseif(printcodsub(iarg)==11)then
    legendobs='maxvz =  max fluid velocity of fluid mix along z  '
  elseif(printcodsub(iarg)==12)then
    legendobs='minvz =  min fluid velocity of fluid mix along z  '
  elseif(printcodsub(iarg)==13)then
    legendobs='cpur  =  remaining time to the end                '
  elseif(printcodsub(iarg)==14)then
    legendobs='cpue  =  elapsed time                             '
  elseif(printcodsub(iarg)==15)then
    legendobs='cpu   =  time for every print interval            '
  elseif(printcodsub(iarg)==16)then
    legendobs='t     =  unscaled time                            '
  elseif(printcodsub(iarg)==17)then
    legendobs='engke =  particle kinetic energy                  '
  elseif(printcodsub(iarg)==18)then
    legendobs='engcf =  configurational energy                   '
  elseif(printcodsub(iarg)==19)then
    legendobs='engto =  total energy                             '
  elseif(printcodsub(iarg)==20)then
    legendobs='tempp =  particle temperature                     '
  elseif(printcodsub(iarg)==21)then
    legendobs='maxpv =  max particle velocity                    '
  elseif(printcodsub(iarg)==22)then
    legendobs='pvx   =  mean particle velocity along x           '
  elseif(printcodsub(iarg)==23)then
    legendobs='pvy   =  mean particle velocity along y           '
  elseif(printcodsub(iarg)==24)then
    legendobs='pvz   =  mean particle velocity along z           '
  elseif(printcodsub(iarg)==25)then
    legendobs='fvx   =  mean fluid velocity along x              '
  elseif(printcodsub(iarg)==26)then
    legendobs='fvy   =  mean fluid velocity along y              '
  elseif(printcodsub(iarg)==27)then
    legendobs='fvz   =  mean fluid velocity along z              '
  elseif(printcodsub(iarg)==28)then
    legendobs='pfx   =  mean particle force along x              '
  elseif(printcodsub(iarg)==29)then
    legendobs='pfy   =  mean particle force along y              '
  elseif(printcodsub(iarg)==30)then
    legendobs='pfz   =  mean particle force along z              '
  elseif(printcodsub(iarg)==31)then
    legendobs='engrt =  rotational energy                        '
  elseif(printcodsub(iarg)==32)then
    legendobs='engkf =  fluid kinetic energy                     '
  elseif(printcodsub(iarg)==33)then
    legendobs='intph =  fluid interphase volume fraction         '
  elseif(printcodsub(iarg)==34)then
    legendobs='rminp =  min pair distance between particles      '
  elseif(printcodsub(iarg)==35)then
    legendobs='pvm   =  mean particle velocity module            '
  elseif(printcodsub(iarg)==36)then
    legendobs='pfm   =  mean particle force module               '
  endif
  legendobs=adjustl(legendobs)
  
  return
  
 end function legendobs
 
 subroutine identify_argument(iarg,arg,printcodsub,lenstring,lfound)
 
!***********************************************************************
!     
!     LBsoft subroutine for identifying the symbolic string of the 
!     output observables
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: lenstring
  integer ,intent(in) :: iarg
  character(len=lenstring) ,allocatable, intent(in) :: arg(:)
  logical, intent(out) :: lfound
  integer, allocatable, intent(inout) :: printcodsub(:)
  
  integer :: inumchar
  character(len=lenstring) :: temps
  
  lfound=.false.
  temps=arg(iarg)
  
! for each symbolic string we associate an integer defined in 
!  printcodsub

  if(findstring('dens1',temps,inumchar,lenstring))then
    printcodsub(iarg)=1
    lfound=.true.
  elseif(findstring('dens2',temps,inumchar,lenstring))then
    printcodsub(iarg)=2
    lfound=.true.
  elseif(findstring('intph',temps,inumchar,lenstring))then
    printcodsub(iarg)=33
    lfound=.true.
  elseif(findstring('maxd1',temps,inumchar,lenstring))then
    printcodsub(iarg)=3
    lfound=.true.
  elseif(findstring('maxd2',temps,inumchar,lenstring))then
    printcodsub(iarg)=4
    lfound=.true.
  elseif(findstring('mind1',temps,inumchar,lenstring))then
    printcodsub(iarg)=5
    lfound=.true.
  elseif(findstring('mind2',temps,inumchar,lenstring))then
    printcodsub(iarg)=6
    lfound=.true.
  elseif(findstring('maxvx',temps,inumchar,lenstring))then
    printcodsub(iarg)=7
    lfound=.true.
  elseif(findstring('minvx',temps,inumchar,lenstring))then
    printcodsub(iarg)=8
    lfound=.true.
  elseif(findstring('maxvy',temps,inumchar,lenstring))then
    printcodsub(iarg)=9
    lfound=.true.
  elseif(findstring('minvy',temps,inumchar,lenstring))then
    printcodsub(iarg)=10
    lfound=.true.
  elseif(findstring('maxvz',temps,inumchar,lenstring))then
    printcodsub(iarg)=11
    lfound=.true.
  elseif(findstring('minvz',temps,inumchar,lenstring))then
    printcodsub(iarg)=12
    lfound=.true.
  elseif(findstring('rminp',temps,inumchar,lenstring))then
    printcodsub(iarg)=34
    lfound=.true.
  elseif(findstring('engke',temps,inumchar,lenstring))then
    printcodsub(iarg)=17
    lfound=.true.
  elseif(findstring('engkf',temps,inumchar,lenstring))then
    printcodsub(iarg)=32
    lfound=.true.
  elseif(findstring('engcf',temps,inumchar,lenstring))then
    printcodsub(iarg)=18
    lfound=.true.
  elseif(findstring('engto',temps,inumchar,lenstring))then
    printcodsub(iarg)=19
    lfound=.true.
  elseif(findstring('engrt',temps,inumchar,lenstring))then
    printcodsub(iarg)=31
    lfound=.true.
  elseif(findstring('tempp',temps,inumchar,lenstring))then
    printcodsub(iarg)=20
    lfound=.true.
  elseif(findstring('maxpv',temps,inumchar,lenstring))then
    printcodsub(iarg)=21
    lfound=.true.
  elseif(findstring('cpur',temps,inumchar,lenstring))then
    printcodsub(iarg)=13
    lfound=.true.
  elseif(findstring('cpue',temps,inumchar,lenstring))then
    printcodsub(iarg)=14
    lfound=.true.
  elseif(findstring('cpu',temps,inumchar,lenstring))then
    printcodsub(iarg)=15
    lfound=.true.
  elseif(findstring('pvm',temps,inumchar,lenstring))then
    printcodsub(iarg)=35
    lfound=.true.
  elseif(findstring('pvx',temps,inumchar,lenstring))then
    printcodsub(iarg)=22
    lfound=.true.
  elseif(findstring('pvy',temps,inumchar,lenstring))then
    printcodsub(iarg)=23
    lfound=.true.
  elseif(findstring('pvz',temps,inumchar,lenstring))then
    printcodsub(iarg)=24
    lfound=.true.
  elseif(findstring('fvx',temps,inumchar,lenstring))then
    printcodsub(iarg)=25
    lfound=.true.
  elseif(findstring('fvy',temps,inumchar,lenstring))then
    printcodsub(iarg)=26
    lfound=.true.
  elseif(findstring('fvz',temps,inumchar,lenstring))then
    printcodsub(iarg)=27
    lfound=.true.
  elseif(findstring('pfm',temps,inumchar,lenstring))then
    printcodsub(iarg)=36
    lfound=.true.
  elseif(findstring('pfx',temps,inumchar,lenstring))then
    printcodsub(iarg)=28
    lfound=.true.
  elseif(findstring('pfy',temps,inumchar,lenstring))then
    printcodsub(iarg)=29
    lfound=.true.
  elseif(findstring('pfz',temps,inumchar,lenstring))then
    printcodsub(iarg)=30
    lfound=.true.
  elseif(findstring('t',temps,inumchar,lenstring))then
    printcodsub(iarg)=16
    lfound=.true.
  else
    lfound=.false.
  endif
  
  
  return
  
 end subroutine identify_argument
 
 subroutine label_argument(narg,printcodsub,printlisub)
 
!***********************************************************************
!     
!     LBsoft subroutine for associating a description string 
!     for each output observables requested in input file
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer ,intent(in) :: narg
  integer, allocatable, intent(in) :: printcodsub(:)
  character(len=20), allocatable, intent(out) :: printlisub(:)
  
  integer :: iarg
  
  if(allocated(printlisub))deallocate(printlisub)
  allocate(printlisub(narg))
  
  do iarg=1,narg
  printlisub(iarg)=repeat(' ',20)
  if(printcodsub(iarg)==1)then
    printlisub(iarg)='dens1 (lu)'
  elseif(printcodsub(iarg)==2)then
    printlisub(iarg)='dens2 (lu)'
  elseif(printcodsub(iarg)==3)then
    printlisub(iarg)='maxd1 (lu)'
  elseif(printcodsub(iarg)==4)then
    printlisub(iarg)='maxd2 (lu)'
  elseif(printcodsub(iarg)==5)then
    printlisub(iarg)='mind1 (lu)'
  elseif(printcodsub(iarg)==6)then
    printlisub(iarg)='mind2 (lu)'
  elseif(printcodsub(iarg)==7)then
    printlisub(iarg)='maxvx (lu)'
  elseif(printcodsub(iarg)==8)then
    printlisub(iarg)='minvx (lu)'
  elseif(printcodsub(iarg)==9)then
    printlisub(iarg)='maxvy (lu)'
  elseif(printcodsub(iarg)==10)then
    printlisub(iarg)='minvy (lu)'
  elseif(printcodsub(iarg)==11)then
    printlisub(iarg)='maxvz (lu)'
  elseif(printcodsub(iarg)==12)then
    printlisub(iarg)='minvz (lu)'
  elseif(printcodsub(iarg)==13)then
    printlisub(iarg)='cpur (s)'
  elseif(printcodsub(iarg)==14)then
    printlisub(iarg)='cpue (s)'
  elseif(printcodsub(iarg)==15)then
    printlisub(iarg)='cpu (s)'
  elseif(printcodsub(iarg)==16)then
    printlisub(iarg)='t (lu)'
  elseif(printcodsub(iarg)==17)then
    printlisub(iarg)='engke (lu)'
  elseif(printcodsub(iarg)==18)then
    printlisub(iarg)='engcf (lu)'
  elseif(printcodsub(iarg)==19)then
    printlisub(iarg)='engto (lu)'
  elseif(printcodsub(iarg)==20)then
    printlisub(iarg)='tempp (KbT)'
  elseif(printcodsub(iarg)==21)then
    printlisub(iarg)='maxpv (lu)'
  elseif(printcodsub(iarg)==22)then
    printlisub(iarg)='pvx (lu)'
  elseif(printcodsub(iarg)==23)then
    printlisub(iarg)='pvy (lu)'
  elseif(printcodsub(iarg)==24)then
    printlisub(iarg)='pvz (lu)'
  elseif(printcodsub(iarg)==25)then
    printlisub(iarg)='fvx (lu)'
  elseif(printcodsub(iarg)==26)then
    printlisub(iarg)='fvy (lu)'
  elseif(printcodsub(iarg)==27)then
    printlisub(iarg)='fvz (lu)'
  elseif(printcodsub(iarg)==28)then
    printlisub(iarg)='pfx (lu)'
  elseif(printcodsub(iarg)==29)then
    printlisub(iarg)='pfy (lu)'
  elseif(printcodsub(iarg)==30)then
    printlisub(iarg)='pfz (lu)'
  elseif(printcodsub(iarg)==31)then
    printlisub(iarg)='engrt (lu)'
  elseif(printcodsub(iarg)==32)then
    printlisub(iarg)='engkf (lu)'
  elseif(printcodsub(iarg)==33)then
    printlisub(iarg)='intph (#)'
  elseif(printcodsub(iarg)==34)then
    printlisub(iarg)='rminp (lu)'
  elseif(printcodsub(iarg)==35)then
    printlisub(iarg)='pvm (lu)'
  elseif(printcodsub(iarg)==36)then
    printlisub(iarg)='pfm (lu)'
  endif
  printlisub(iarg)=adjustr(printlisub(iarg))
  enddo
  
  return
  
 end subroutine label_argument
 
 function legendobs_vtk(iarg,printcodsub)
 
!***********************************************************************
!     
!     LBsoft function for returning the legend_vtk which is associated
!      tothe integer contained in the printcod array for vtk files
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: iarg
  integer, intent(in), allocatable, dimension(:) :: printcodsub
  
  character(len=52) :: legendobs_vtk
  
  legendobs_vtk=repeat(' ',52)
  if(printcodsub(iarg)==1)then
    legendobs_vtk='rho1  =  fluid density of first component         '
  elseif(printcodsub(iarg)==2)then
    legendobs_vtk='rho2  =  fluid density of second component        '
  elseif(printcodsub(iarg)==3)then
    legendobs_vtk='phase =  phase field of the two fluid components  '
  elseif(printcodsub(iarg)==4)then
    legendobs_vtk='vel   =  vector velocity field of the fluid       '
  elseif(printcodsub(iarg)==5)then
    legendobs_vtk='part  =  particle positions and orientations      '
  endif
  legendobs_vtk=adjustl(legendobs_vtk)
  
  return
  
 end function legendobs_vtk
 
 subroutine identify_argument_vtk(iarg,arg,printcodsub,lenstring,lfound)
 
!***********************************************************************
!     
!     LBsoft subroutine for identifying the symbolic string of the 
!     output observables for the vtk files
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: lenstring
  integer ,intent(in) :: iarg
  character(len=lenstring) ,allocatable, intent(in) :: arg(:)
  logical, intent(out) :: lfound
  integer, allocatable, intent(inout) :: printcodsub(:)
  
  integer :: inumchar
  character(len=lenstring) :: temps
  
  lfound=.false.
  temps=arg(iarg)
  
! for each symbolic string we associate an integer defined in 
!  printcodsub

  if(findstring('phase',temps,inumchar,lenstring))then
    printcodsub(iarg)=3
    lfound=.true.
  elseif(findstring('part',temps,inumchar,lenstring))then
    printcodsub(iarg)=5
    lfound=.true.
  elseif(findstring('rho1',temps,inumchar,lenstring))then
    printcodsub(iarg)=1
    lfound=.true.
  elseif(findstring('rho2',temps,inumchar,lenstring))then
    printcodsub(iarg)=2
    lfound=.true.
  elseif(findstring('vel',temps,inumchar,lenstring))then
    printcodsub(iarg)=4
    lfound=.true.
  else
    lfound=.false.
  endif
  
  
  return
  
 end subroutine identify_argument_vtk
 
 subroutine identify_argument_xyz(iarg,arg,printcodsub,lenstring,lfound)
 
!***********************************************************************
!     
!     LBsoft subroutine for identifying the symbolic string of the 
!     input xyz variables that have to be read
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: lenstring
  integer ,intent(in) :: iarg
  character(len=lenstring) ,allocatable, intent(in) :: arg(:)
  logical, intent(out) :: lfound
  integer, allocatable, intent(inout) :: printcodsub(:)
  
  integer :: inumchar
  character(len=lenstring) :: temps
  
  lfound=.false.
  temps=arg(iarg)
  
! for each symbolic string we associate an integer defined in 
!  printcodsub

!  if(findstring('char',temps,inumchar,lenstring))then
!    printcodsub(iarg)=1
!    lfound=.true.
!  else
  if(findstring('mass',temps,inumchar,lenstring))then
    printcodsub(iarg)=2
    lfound=.true.
  elseif(findstring('radx',temps,inumchar,lenstring))then
    printcodsub(iarg)=3
    lfound=.true.
  elseif(findstring('rady',temps,inumchar,lenstring))then
    printcodsub(iarg)=4
    lfound=.true.
  elseif(findstring('radz',temps,inumchar,lenstring))then
    printcodsub(iarg)=5
    lfound=.true.
  elseif(findstring('rad',temps,inumchar,lenstring))then
    printcodsub(iarg)=6
    lfound=.true.
  elseif(findstring('vx',temps,inumchar,lenstring))then
    printcodsub(iarg)=7
    lfound=.true.
  elseif(findstring('vy',temps,inumchar,lenstring))then
    printcodsub(iarg)=8
    lfound=.true.
  elseif(findstring('vz',temps,inumchar,lenstring))then
    printcodsub(iarg)=9
    lfound=.true.
  elseif(findstring('ox',temps,inumchar,lenstring))then
    printcodsub(iarg)=10
    lfound=.true.
  elseif(findstring('oy',temps,inumchar,lenstring))then
    printcodsub(iarg)=11
    lfound=.true.
  elseif(findstring('oz',temps,inumchar,lenstring))then
    printcodsub(iarg)=12
    lfound=.true.
  elseif(findstring('phi',temps,inumchar,lenstring))then
    printcodsub(iarg)=13
    lfound=.true.
  elseif(findstring('theta',temps,inumchar,lenstring))then
    printcodsub(iarg)=14
    lfound=.true.
  elseif(findstring('psi',temps,inumchar,lenstring))then
    printcodsub(iarg)=15
    lfound=.true.
  else
    lfound=.false.
  endif
  
  
  return
  
 end subroutine identify_argument_xyz
 
 function legendobs_xyz(iarg,xyzcod)
 
!***********************************************************************
!     
!     LBsoft function for returning the legend which is associated to
!     the integer contained in the xyzcod array
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: iarg
  integer, intent(in), allocatable, dimension(:) :: xyzcod
  
  character(len=52) :: legendobs_xyz
  
  legendobs_xyz=repeat(' ',52)
!  if(xyzcod(iarg)==1)then
!    legendobs_xyz='char =  charge of the particle                    '
!  else
  if(xyzcod(iarg)==2)then
    legendobs_xyz='mass =  mass of the particle                      '
  elseif(xyzcod(iarg)==3)then
    legendobs_xyz='radx =  radius of the particle along x            '
  elseif(xyzcod(iarg)==4)then
    legendobs_xyz='rady =  radius of the particle along y            '
  elseif(xyzcod(iarg)==5)then
    legendobs_xyz='radz =  radius of the particle along z            '
  elseif(xyzcod(iarg)==6)then
    legendobs_xyz='rad  =  radius of the spherical particle          '
  elseif(xyzcod(iarg)==7)then
    legendobs_xyz='vx   =  particle velocity along x                 '
  elseif(xyzcod(iarg)==8)then
    legendobs_xyz='vy   =  particle velocity along y                 '
  elseif(xyzcod(iarg)==9)then
    legendobs_xyz='vz   =  particle velocity along z                 '
  elseif(xyzcod(iarg)==10)then
    legendobs_xyz='ox   =  particle angular velocity along x         '
  elseif(xyzcod(iarg)==11)then
    legendobs_xyz='oy   =  particle angular velocity along y         '
  elseif(xyzcod(iarg)==12)then
    legendobs_xyz='oz   =  particle angular velocity along z         '
  elseif(xyzcod(iarg)==13)then
    legendobs_xyz='phi  =  first Euler angle [−Pi,Pi]                '
  elseif(xyzcod(iarg)==14)then
    legendobs_xyz='theta=  second Euler angle [−Pi/2,Pi/2]           '
  elseif(xyzcod(iarg)==15)then
    legendobs_xyz='psi  =  third Euler angle [−Pi,Pi]                '
  endif
  legendobs_xyz=adjustl(legendobs_xyz)
  
  return
  
 end function legendobs_xyz
 
 subroutine outprint_driver(k,timesub)
 
!***********************************************************************
!     
!     LBsoft subroutine for managing the output subroutines
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: k
  double precision, intent(in) :: timesub
  
  integer :: i,j
  double precision :: tempint
  
  if(.not.lprintlist)return
  if(mod(k,iprinttime)/=0)return
  
  tempint=dble(iprinttime)*tstep
  call compute_statistic(tempint,timesub,k)
  
  if(idrank/=0)return
  
  do i=1,nprintlist
    j=printlist(i)
    xprint(i)=statdata(j)
  enddo

  call outprint_term(k,nprintlist,xprint,printarg)
  
  return
  
 end subroutine outprint_driver
 
 subroutine outprint_term(k,nprint1,dprint,printargsub)
 
!***********************************************************************
!     
!     LBsoft subroutine for printing records on the terminal
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  integer,intent(in) :: k,nprint1
  double precision,intent(in),dimension(1:nprint1) :: dprint
  character(len=20),intent(in),dimension(1:nprint1) :: printargsub
  
  integer :: i
  
  logical, save :: lfirst=.true.
  integer,save :: icount=0
  
  if(idrank/=0)return
  
  if(lfirst .or. (mod(icount,100)==0))then
    lfirst=.false.    
    write(6,'(a10)',advance="no")'#*********'
    do i=1,nprint1
      write(6,'(a20)',advance="no")' *******************'
    enddo
    write(6,'(1x)')           
    write(6,'(a10)',advance="no")'#    nstep'
    do i=1,nprint1
      write(6,'(a20)',advance="no")printargsub(i)
    enddo
    write(6,'(1x)')
    write(6,'(a10)',advance="no")'#*********'
    do i=1,nprint1
      write(6,'(a20)',advance="no")' *******************'
    enddo
    write(6,'(1x)')
  endif
  
  write(6,'(i10)',advance="no")k
  do i=1,nprint1
    write(6,'(g20.10)',advance="no")dprint(i)
  enddo
  
  write(6,'(1x)')
  call flush(6)
  
  icount=icount+1
  
  return
  
 end subroutine outprint_term
 
 subroutine read_input_atom(inputunit,inputname,ltypes,xxs,yys,zzs,others, &
  lvelocity)
 
!***********************************************************************
!     
!     LBsoft subroutine for reading the input atomic file
!
!     note: the input file is written in a modified xyz format
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
  integer, intent(in) :: inputunit
  character(len=*), intent(in) :: inputname
  integer, allocatable, dimension(:) :: ltypes
  real(kind=PRC), allocatable, dimension(:) :: xxs,yys,zzs
  real(kind=PRC), allocatable, dimension(:,:) :: others
  logical :: lvelocity
  
  character(len=maxlen) :: redstring,directive
  integer :: inumchar,i,j,nwords,iline,itest,itype
  character(len=maxlen) ,allocatable :: outwords(:),outwords2(:)
  logical :: lerror1,lerror2,lerror3,lerror4,safe,lexists
  integer :: temp_natms_tot
  logical :: lxyzlisterror,lfoundprintxyz,lxyzlist,lxyzlist1,lxyzlist2, &
   lqinput,lerror5,lmass,lvel(3),lrdim
  
  character(len=*),parameter :: of='(a)'
  integer, parameter :: natmname=4
  real(kind=PRC) :: volatms = ZERO
  real(kind=PRC) :: volfrac = ZERO
  
  
  temp_natms_tot=0
  lerror1=.false.
  lerror2=.false.
  lerror3=.false.
  lerror4=.false.
  lerror5=.false.
  lxyzlist1=.false.
  lxyzlist2=.false.
  lxyzlisterror=.false.
  
  if(idrank==0)then
  
!   check if the input file exist
    inquire(file=inputname,exist=lexists)
    if(.not.lexists)then
      lerror1=.true.
      goto 120
    endif
    
!   open the inout file
    open(unit=inputunit,file=inputname,status='old',action='read', &
    iostat=itest)
    if(itest/=0)then
      lerror2=.true.
      goto 120
    endif
    
    write(6,of)"                                                                               "
    write(6,of)"*****************************READ XYZ INPUT FILE*******************************"
    write(6,of)"                                                                               "
    
    write(6,'(/,3a)')'start reading ',trim(inputname),' file'
    
!   counter the line which are read in the input file
    iline=0
    
!   read the input file

    call getline(safe,inputunit,maxlen,redstring)
    iline=iline+1
    if(.not.safe)then
      call warning(11,dble(iline))
      lerror3=.true.
      goto 120
    endif
    call strip(redstring,maxlen)
    call copystring(redstring,directive,maxlen)
    temp_natms_tot=intstr(directive,maxlen,inumchar)
    if(temp_natms_tot<1)then
      call warning(11,dble(iline))
      lerror3=.true.
      goto 120
    endif
    
    call getline(safe,inputunit,maxlen,redstring)
    iline=iline+1
    if(.not.safe)then
      call warning(11,dble(iline))
      lerror3=.true.
      goto 120
    endif
    call strip(redstring,maxlen)
    call copystring(redstring,directive,maxlen)
    call findwords(nwords,outwords,directive,maxlen)
    nxyzlist=max(0,nwords-2)
    if(nwords>=2)then
      directive=outwords(1)
      if(findstring('read',directive,inumchar,maxlen))then
        lxyzlist1=.true.
      endif
      directive=outwords(2)
      if(findstring('list',directive,inumchar,maxlen))then
        lxyzlist2=.true.
      endif
    endif
    lxyzlist=(lxyzlist1.and.lxyzlist2)
    if(lxyzlist)then
      write(6,'(/,a)')'read list invoked'
      if(nxyzlist>0)allocate(xyzlist(nxyzlist))
      if(allocated(outwords2))deallocate(outwords2)
      if(nxyzlist>0)allocate(outwords2(nxyzlist))
      do i=1,nxyzlist
        outwords2(i)=outwords(i+2)
        call identify_argument_xyz(i,outwords2,xyzlist,maxlen, &
         lfoundprintxyz)
        lxyzlisterror=(lxyzlisterror .or. (.not.lfoundprintxyz))
        if(lxyzlisterror)then
          call warning(11,dble(iline))
          call warning(12)
          call warning(13)
          goto 120
        endif
      enddo
      call print_legend_xyz(6)
      lqinput=.false.
      lmass=.false.
      lrdim=.false.
      lvel(1:3)=.false.
      do i=1,nxyzlist
        if(xyzlist(i)>=3 .and. xyzlist(i)<=5)then
          if(all(ishape==0))then
            lerror4=.true.
            call warning(11,dble(iline))
            call warning(16)
            goto 120
          endif
        endif
        if(xyzlist(i)>=13 .and. xyzlist(i)<=15)lqinput=.true.
        if(xyzlist(i)==2)lmass=.true.
        if(xyzlist(i)==6)lrdim=.true.
        if(xyzlist(i)==7)lvel(1)=.true.
        if(xyzlist(i)==8)lvel(2)=.true.
        if(xyzlist(i)==9)lvel(3)=.true.
      enddo
      if((.not. all(lumass(1:ntype))).and.(.not. lmass))then
        lerror5=.true.
        call warning(11,dble(iline))
        call warning(25)
        goto 120
      endif
      if(all(ishape==0))then
        if((.not. lrdim).and.(.not. all(lurdim(1:ntype))))then
          lerror5=.true.
          call warning(11,dble(iline))
          call warning(28)
          goto 120
        endif
      endif
      if(all(lvel))then
        lvelocity=.true.
      else
        lvelocity=.false.
      endif
      if((.not. lvelocity).and.(.not. linit_temp))then
        lerror5=.true.
        call warning(11,dble(iline))
        call warning(27)
        goto 120
      endif
    endif
    if(allocated(ltypes))deallocate(ltypes)
    if(allocated(xxs))deallocate(xxs)
    if(allocated(yys))deallocate(yys)
    if(allocated(zzs))deallocate(zzs)
    allocate(ltypes(temp_natms_tot))
    allocate(xxs(temp_natms_tot))
    allocate(yys(temp_natms_tot))
    allocate(zzs(temp_natms_tot))
    if(nxyzlist>0)then
      if(allocated(others))deallocate(others)
      allocate(others(nxyzlist,temp_natms_tot))
      do i=1,temp_natms_tot
        call getline(safe,inputunit,maxlen,redstring)
        iline=iline+1
        if(.not.safe)then
          call warning(11,dble(iline))
          lerror3=.true.
          goto 120
        endif
        call strip(redstring,maxlen)
        call find_type(redstring(1:natmname),natmname,itype)
        if(itype==0)then
          !error: itype not found
          lerror5=.true.
          call warning(11,dble(iline))
          call warning(41,dble(iline),redstring(1:natmname))
          goto 120
        endif
        ltypes(i)=itype
        directive(1:maxlen)=redstring(1+natmname:maxlen)// &
         repeat(' ',natmname)
        xxs(i)=dblstr(directive,maxlen,inumchar)
        yys(i)=dblstr(directive,maxlen,inumchar)
        zzs(i)=dblstr(directive,maxlen,inumchar)
        do j=1,nxyzlist
          others(j,i)=dblstr(directive,maxlen,inumchar)
          if(xyzlist(j)==6)then
            if(mod(others(j,i)-HALF,ONE)/=ZERO)then
              call warning(32)
              call warning(11,dble(iline))
              lerror3=.true.
            endif
          endif
        enddo
      enddo
    else
      do i=1,temp_natms_tot
        call getline(safe,inputunit,maxlen,redstring)
        iline=iline+1
        if(.not.safe)then
          call warning(11,dble(iline))
          lerror3=.true.
          goto 120
        endif
        call strip(redstring,maxlen)
        call find_type(redstring(1:natmname),natmname,itype)
        if(itype==0)then
          !error: itype not found
          lerror5=.true.
          call warning(11,dble(iline))
          call warning(41,dble(iline),redstring(1:natmname))
          goto 120
        endif
        ltypes(i)=itype
        directive(1:maxlen)=redstring(1+natmname:maxlen)// &
         repeat(' ',natmname)
        xxs(i)=dblstr(directive,maxlen,inumchar)
        yys(i)=dblstr(directive,maxlen,inumchar)
        zzs(i)=dblstr(directive,maxlen,inumchar)
      enddo
    endif
    
  endif
  
120 continue
  
  call bcast_world_l(lerror1)
  call bcast_world_l(lerror2)
  call bcast_world_l(lerror3)
  call bcast_world_l(lerror4)
  call bcast_world_l(lerror5)
  call bcast_world_l(lxyzlisterror)
  
  if(lerror1)call error(17)
  if(lerror2)call error(18)
  if(lerror3)call error(19)
  if(lerror4)call error(19)
  if(lerror5)call error(19)
  if(lxyzlisterror)call error(19)
  
  call bcast_world_i(nxyzlist)
  if(nxyzlist>0)then
    if(idrank/=0)then
      allocate(xyzlist(nxyzlist))
    endif
    call bcast_world_iarr(xyzlist,nxyzlist)
  endif
  
  if(idrank==0)then
    !compute total volume of particles
    volatms = ZERO
    do i=1,temp_natms_tot
      itype=ltypes(i)
      if(ishape(itype)==1)then
        !oblate to be added for now is a sphere
        volatms = volatms + FOUR/THREE*Pi*(rdimx(itype))**THREE
      else
        volatms = volatms + FOUR/THREE*Pi*(rdimx(itype))**THREE
      endif
    enddo
    volfrac=volatms/real(nx*ny*nz,kind=PRC)
    
    write(6,'(3a,/)')'file ',trim(inputname),' correctly closed'
    write(6,'(a,i12)') 'Total number of atoms    : ', temp_natms_tot
    write(6,'(a,f12.6)') 'Particle volume fraction : ', volfrac
    write(6,of)"                                                                               "
    write(6,of)"*******************************************************************************"
    write(6,of)"                                                                               "
  endif
  
  call bcast_world_i(temp_natms_tot)
  call set_natms_tot(temp_natms_tot)
  
  return
  
 end subroutine read_input_atom

 
 subroutine print_legend_xyz(iu)
 
!***********************************************************************
!     
!     LBsoft subroutine for printing the legend of the observables
!     written in input xyz file
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: iu
  
  integer :: i
  
  character(len=*),parameter :: of='(a)'
  
  if(idrank/=0)return
  
  write(iu,of)"                                                                               "
  write(iu,of)"LIST SPECIFIED IN XYZ FILE:"
  do i=1,nxyzlist
    write(iu,'(a)')legendobs_xyz(i,xyzlist)
  enddo
  write(iu,of)"                                                                               "
  
  return
  
 end subroutine print_legend_xyz
 
 end module io_mod


