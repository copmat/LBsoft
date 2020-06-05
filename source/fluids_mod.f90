
#include <default_macro.h>

 module fluids_mod
 
!***********************************************************************
!     
!     LBsoft module containing variable and subroutines of the
!     fluid components
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
 use fluids_lattices_mod
 use fluids_equilibrium_mod
 use version_mod,    only : idrank,mxrank,or_world_larr,finalize_world, &
                   get_sync_world,sum_world_farr,commspop,commrpop, &
                   i4find,i4back,ownern,deallocate_ownern,commexch_dens, &
                   commwait_dens,comm_init_isfluid,commexch_vel_component, &
                   commwait_vel_component,commexch_isfluid, &
                   commwait_isfluid,commexch_single_halo, &
                   commwait_single_halo,collective_readFile_pops, &
                   collective_writeFile_pops,collective_readFile, &
                   collective_readFile_int,collective_writeFile, &
                   collective_writeFile_int
 use error_mod
 use utility_mod, only : Pi,modulvec,cross,dot,gauss,ibuffservice, &
                   allocate_array_ibuffservice,buffservice,conv_rad, &
                   allocate_array_buffservice,lbuffservice, &
                   allocate_array_lbuffservice,xcross,ycross,zcross, &
                   nbuffservice3d,buffservice3d,int_cube_sphere,fcut, &
                   rand_noseeded,linit_seed,gauss_noseeded,write_fmtnumb, &
                   openLogFile
 use fluids_bc_mod,only: ibctype,set_boundary_conditions_type, &
                   ixpbc,iypbc,izpbc,isfluid,new_isfluid,bcfluid, &
                   allocate_isfluid,nbcdir,bc_rhoR,bc_rhoB,bc_u,bc_v, &
                   bc_w,bc_flow,bc_omega,set_isfluid_bcfluid, &
                   set_lread_isfluid,lread_isfluid,pimage,ndouble, &
                   exdouble,eydouble,ezdouble,manage_bc_hvar_selfcomm, &
                   set_bc_variable_hvar,apply_bounceback_pop_halfway, &
                   manage_bc_singlehalo_selfcomm,manage_bc_pop_selfcomm, &
                   compute_densities_wall,driver_bc_isfluid, &
                   initialiaze_manage_bc_pop_selfcomm,set_double_range

 
 implicit none
 
 private
 
 integer, save, protected, public :: nx=0
 integer, save, protected, public :: ny=0
 integer, save, protected, public :: nz=0
 
 !max
 integer, save, public :: minx, maxx, miny, maxy, minz, maxz
 integer, save, public :: wminx, wmaxx, wminy, wmaxy, wminz, wmaxz
 integer, save, public :: ix,iy,iz,mynx,myny,mynz
 
 !max
 integer, save, protected, public :: LBintegrator=0
 
 integer, save, protected, public :: idistselect=0
 
 !timestep interval to rescale the total fluid mass
 integer, save, protected, public :: imass_rescale=100
 
 !interpolation mode of omega (0=phase field interp;1=Grunau interp)
 integer, parameter :: iselomegaint=0
  
 integer, save, protected, public, allocatable, dimension(:) :: &
  typeobjectliq
 
 !integer denoting the LOCAL fluid nodes at actual time 
 !(it changes if particles are presents)
 integer, save :: nfluidnodes=0
 !integer denoting the total fluid nodes at initial time
 integer, save :: totnfluidnodes_init=0
 
 !number of object drawn in fluid initialization
 integer, save, protected, public :: nobjectliq=0
 
 logical, save, protected, public :: lalloc_hvars=.false.
 logical, save, protected, public :: lalloc_pops=.false.
 
 logical, save, protected, public :: lforce_add=.false.
 logical, save, protected, public :: lShanChen=.false.
 logical, save, protected, public :: lsing_SC=.false.
 logical, save, protected, public :: lpair_SC=.false.
 
 logical, save, protected, public :: lColourG=.false.
 
 logical, save, protected, public :: lunique_omega=.false.
 
 logical, save, protected, public :: lcap_force=.false.
 
 !rescale the total fluid mass
 logical, save, protected, public :: lmass_rescale=.false.
 
 logical, save, protected, public :: lsingle_fluid=.false.
 
 logical, save, protected, public :: lvorticity=.false.
 
 logical, save, protected, public :: lwall_SC=.false.
 
 logical, save, protected, public :: lpart_SC=.false.
 
 logical, save, protected, public :: lexch_dens=.false.
 logical, save, protected, public :: lexch_u=.false.
 logical, save, protected, public :: lexch_v=.false.
 logical, save, protected, public :: lexch_w=.false.
 
 !bounceback default is halfway
 logical, save, protected, public :: lbc_halfway=.true.
 logical, save, protected, public :: lbc_fullway=.false.
 
 !are there specified object in fuild initialization
 logical, save, protected, public :: lobjectliq=.false.
 
 real(kind=PRC), save, protected, public :: t_LB = ONE
 
 real(kind=PRC), save, protected, public :: beta         = ZERO
 real(kind=PRC), save, protected, public :: akl          = ZERO
 real(kind=PRC), save, protected, public :: awall        = ZERO
 real(kind=PRC), save, protected, public :: tauR         = ZERO
 real(kind=PRC), save, protected, public :: tauB         = ZERO
 real(kind=PRC), save, protected, public :: unique_omega = ZERO
 real(kind=PRC), save, protected, public :: viscR        = ZERO
 real(kind=PRC), save, protected, public :: viscB        = ZERO
 real(kind=PRC), save, protected, public :: meanR        = ZERO
 real(kind=PRC), save, protected, public :: meanB        = ZERO
 real(kind=PRC), save, protected, public :: stdevR       = ZERO
 real(kind=PRC), save, protected, public :: stdevB       = ZERO
 real(kind=PRC), save, protected, public :: backR        = ZERO
 real(kind=PRC), save, protected, public :: backB        = ZERO
 real(kind=PRC), save, protected, public :: initial_u    = ZERO
 real(kind=PRC), save, protected, public :: initial_v    = ZERO
 real(kind=PRC), save, protected, public :: initial_w    = ZERO
 real(kind=PRC), save, protected, public :: ext_fu       = ZERO
 real(kind=PRC), save, protected, public :: ext_fv       = ZERO
 real(kind=PRC), save, protected, public :: ext_fw       = ZERO
 real(kind=PRC), save, protected, public :: cap_force    = ZERO
 
 !Shan Chen coupling constant 
 real(kind=PRC), save, protected, public :: pair_SC = ZERO
 !Shan Chen wall coupling constant
 real(kind=PRC), save, protected, public :: wallR_SC = ONE
 real(kind=PRC), save, protected, public :: wallB_SC = ONE
 !contact angle of shan chen force on particle (in radiant)
 real(kind=PRC), save, protected, public :: theta_SC = ZERO
 !deviation of contact angle of shan chen force on particle (in radiant)
 real(kind=PRC), save, protected, public :: devtheta_SC = ZERO
 
 !Shan Chen particle color factor(see definition at page 4 of PRE 83,046707)
 real(kind=PRC), save, protected, public :: partR_SC = ZERO
 real(kind=PRC), save, protected, public :: partB_SC = ZERO
 
 !Colour gradient coupling constant 
 real(kind=PRC), save, protected, public :: sigma_CG = ZERO
 real(kind=PRC), save, protected, public :: alphaR_CG = ZERO
 real(kind=PRC), save, protected, public :: alphaB_CG = ZERO
 real(kind=PRC), save, protected, public :: beta_CG = ZERO
 
 !Grunau mode of omega interpolation
 real(kind=PRC), parameter :: Delta_Grunau = real(0.98d0,kind=PRC)
 real(kind=PRC), save :: s1_R = ZERO
 real(kind=PRC), save :: s2_R = ZERO
 real(kind=PRC), save :: s3_R = ZERO
 real(kind=PRC), save :: t1_B = ZERO
 real(kind=PRC), save :: t2_B = ZERO
 real(kind=PRC), save :: t3_B = ZERO
 
 real(kind=PRC), save, protected, public :: dmass_rescale = ZERO
 
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: wwfluid
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:),target :: rhoR
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:),target :: rhoB
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:),target  :: u
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:),target  :: v
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:),target  :: w
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: omega
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: fuR
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: fvR
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: fwR
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: fuB
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: fvB
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: fwB
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:),target :: psiR
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:),target :: psiB
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: gradpsixR
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: gradpsiyR
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: gradpsizR
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: gradpsixB
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: gradpsiyB
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:,:) :: gradpsizB
 
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:) :: &
  objectdata
 
!real denoting the LOCAL fluid mass at actual time for red and blue
!(it changes if particles are presents)
 real(kind=PRC), save :: nfluidmassR = ZERO
 real(kind=PRC), save :: nfluidmassB = ZERO
!real denoting the TOTAL fluid mass at initial time for red and blue
 real(kind=PRC), save :: totnfluidmassR_init = ZERO
 real(kind=PRC), save :: totnfluidmassB_init = ZERO
 
! coeff for Colour Gradient
 real(kind=PRC), dimension(0:links), protected, public :: phiR_CG
 real(kind=PRC), dimension(0:links), protected, public :: phiB_CG
 
 real(kind=PRC),save,public,allocatable,dimension(:,:,:,:) :: &
  fluidR,fluidB
 

 integer,save, private :: iounit(3)

 
 public :: set_random_dens_fluids
 public :: set_initial_dens_fluids
 public :: set_initial_vel_fluids
 public :: set_mean_value_dens_fluids
 public :: set_stdev_value_dens_fluids
 public :: set_back_value_dens_fluids
 public :: set_initial_dist_type
 public :: set_initial_dim_box
 public :: set_mean_value_vel_fluids
 public :: allocate_fluids
 public :: initialize_fluids
 public :: initialize_fluid_force
 public :: set_value_ext_force_fluids
 public :: compute_fluid_force_sc
 public :: driver_collision_fluids
 public :: driver_bc_densities
 public :: driver_bc_velocities
 public :: set_fluid_force_sc
 public :: set_value_viscosity
 public :: set_value_tau
 public :: compute_omega
 public :: driver_streaming_fluids
 public :: moments_fluids
 public :: compute_psi_sc
 public :: set_lsingle_fluid
 public :: probe_red_moments_in_node
 public :: probe_blue_moments_in_node
 public :: set_fluid_wall_sc
 public :: set_fluid_particle_sc
 public :: set_theta_particle_sc
 public :: set_LBintegrator_type
 public :: initialize_isfluid_bcfluid
 public :: test_fake_pops
 public :: probe_pops_in_node
 public :: driver_bc_pop_selfcomm
 public :: driver_initialiaze_manage_bc_selfcomm
 public :: set_lbc_halfway
 public :: set_lbc_fullway
 public :: driver_apply_bounceback_halfway_pop
 public :: print_all_hvar
 public :: init_particle_2_isfluid
 public :: push_comm_isfluid
 public :: particle_bounce_back
 public :: particle_bounce_back_phase2
 public :: particle_delete_fluids
 public :: particle_create_fluids
 public :: erase_fluids_in_particles
 public :: print_all_pops
 public :: print_all_pops_center
 public :: print_all_pops_area_shpere
 public :: omega_to_viscosity
 public :: compute_sc_particle_interact_phase1
 public :: compute_sc_particle_interact_phase2
 public :: print_all_pops2
 public :: driver_bc_all_pops
 public :: dumpHvar, restorePops, restoreHvar
 public :: driver_bc_psi
 public :: fixPops
 public :: set_cap_force
 public :: set_mass_rescale
 public :: rescale_fluid_mass
 public :: bin_oneFile
 public :: dump_oneFile
 public :: restore_oneFile
 public :: set_objectliq
 public :: set_objectdata
 public :: set_fluid_force_cg
 
 contains
 
 subroutine allocate_fluids(lparticles)

!***********************************************************************
!     
!     LBsoft subroutine for allocating arrays which 
!     describe the fluids
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************

     implicit none
     
     logical, intent(in) :: lparticles

     integer, parameter :: nistatmax=100
     integer, dimension(nistatmax) :: istat
!max     integer :: ix,iy,iz,mynx,myny,mynz

     integer :: i,j,k,l
     logical, dimension(1) :: ltest=.false.
     
     logical, parameter :: lverbose=.true.

     ix=minx-nbuff
     iy=miny-nbuff
     iz=minz-nbuff
     mynx=maxx+nbuff
     myny=maxy+nbuff
     mynz=maxz+nbuff

     istat=0
     
     call allocate_isfluid(lparticles)
     
     allocate(rhoR(ix:mynx,iy:myny,iz:mynz),stat=istat(2))

     allocate(u(ix:mynx,iy:myny,iz:mynz),stat=istat(3))
     allocate(v(ix:mynx,iy:myny,iz:mynz),stat=istat(4))
     allocate(w(ix:mynx,iy:myny,iz:mynz),stat=istat(5))

     allocate(fuR(ix:mynx,iy:myny,iz:mynz),stat=istat(6))
     allocate(fvR(ix:mynx,iy:myny,iz:mynz),stat=istat(7))
     allocate(fwR(ix:mynx,iy:myny,iz:mynz),stat=istat(8))

     if(lShanChen)then
        allocate(gradpsixR(ix:mynx,iy:myny,iz:mynz),stat=istat(12))
        allocate(gradpsiyR(ix:mynx,iy:myny,iz:mynz),stat=istat(13))
        allocate(gradpsizR(ix:mynx,iy:myny,iz:mynz),stat=istat(14))
        allocate(psiR(ix:mynx,iy:myny,iz:mynz),stat=istat(19))
     endif
     
     allocate(wwfluid(ix:mynx,iy:myny,iz:mynz),stat=istat(20))

     allocate(omega(ix:mynx,iy:myny,iz:mynz),stat=istat(18))
     
     allocate(fluidR(ix:mynx,iy:myny,iz:mynz,0:links),stat=istat(30))


     
     

     if(.not. lsingle_fluid)then
        
        allocate(rhoB(ix:mynx,iy:myny,iz:mynz),stat=istat(2))

        allocate(fuB(ix:mynx,iy:myny,iz:mynz),stat=istat(9))
        allocate(fvB(ix:mynx,iy:myny,iz:mynz),stat=istat(10))
        allocate(fwB(ix:mynx,iy:myny,iz:mynz),stat=istat(11))

        if(lShanChen)then
           allocate(gradpsixB(ix:mynx,iy:myny,iz:mynz),stat=istat(15))
           allocate(gradpsiyB(ix:mynx,iy:myny,iz:mynz),stat=istat(16))
           allocate(gradpsizB(ix:mynx,iy:myny,iz:mynz),stat=istat(17))
           allocate(psiB(ix:mynx,iy:myny,iz:mynz),stat=istat(20))
        endif
        
        allocate(fluidB(ix:mynx,iy:myny,iz:mynz,0:links),stat=istat(60))

        
     endif

     ltest=.false.
     if(any(istat.ne.0))then
        do i=1,nistatmax
           if(istat(i).ne.0)exit
        enddo
        call warning(2,dble(i))
        ltest=.true.
     endif

     call or_world_larr(ltest,1)
     if(ltest(1))call error(4)
     
     wminx=minx
     wminy=miny
     wminz=minz
     wmaxx=maxx
     wmaxy=maxy
     wmaxz=maxz

     if(minx.le.1) then
        if(ixpbc.eq.1) then
           wminx=1
        else
           wminx=0
        endif
        minx=1
     else 
        wminx=minx-1
     endif
     if(miny.le.1) then
        if(iypbc.eq.1) then
           wminy=1
        else
           wminy=0
        endif
        miny=1
     else 
        wminy=miny-1
     endif
     if(minz.le.1) then
        if(izpbc.eq.1) then
           wminz=1
        else
           wminz=0
        endif
        minz=1
     else 
        wminz=minz-1
     endif
     if(maxx.ge.nx) then
        if(ixpbc.eq.1) then
           wmaxx=nx
        else
           wmaxx=nx+1
        endif
        maxx=nx
     else 
        wmaxx=maxx+1
     endif
     if(maxy.ge.ny) then
        if(iypbc.eq.1) then
           wmaxy=ny
        else
           wmaxy=ny+1
        endif
        maxy=ny
     else 
        wmaxy=maxy+1
     endif
     if(maxz.ge.nz) then
        if(izpbc.eq.1) then
           wmaxz=nz
        else
           wmaxz=nz+1
        endif
        maxz=nz
     else 
        wmaxz=maxz+1
     endif
      
   if(lverbose)then
     if(idrank==0)write(6,*)'ixpbc=',ixpbc,'iypbc=',iypbc,'izpbc=',izpbc
     call flush(6)
     call get_sync_world
     ! Only for few procs
     if (10 >= mxrank) then
     do i=0,mxrank-1
       if(i==idrank)then
         write(6,*)'fluids id=',idrank,'minx=',minx,'maxx=',maxx, &
          'miny=',miny,'maxy=',maxy,'minz=',minz,'maxz=',maxz
       endif
       call flush(6)
       call get_sync_world
     enddo
     endif
   endif
   
   rhoR(ix:mynx,iy:myny,iz:mynz)=MINDENS

   u(ix:mynx,iy:myny,iz:mynz)=ZERO
   v(ix:mynx,iy:myny,iz:mynz)=ZERO
   w(ix:mynx,iy:myny,iz:mynz)=ZERO

   fuR(ix:mynx,iy:myny,iz:mynz)=ZERO
   fvR(ix:mynx,iy:myny,iz:mynz)=ZERO
   fwR(ix:mynx,iy:myny,iz:mynz)=ZERO

   if(lShanChen)then
     gradpsixR(ix:mynx,iy:myny,iz:mynz)=ZERO
     gradpsiyR(ix:mynx,iy:myny,iz:mynz)=ZERO
     gradpsizR(ix:mynx,iy:myny,iz:mynz)=ZERO
     psiR(ix:mynx,iy:myny,iz:mynz)=ZERO
     wwfluid(ix:mynx,iy:myny,iz:mynz)=0
   endif

   omega(ix:mynx,iy:myny,iz:mynz)=ZERO

   fluidR(ix:mynx,iy:myny,iz:mynz,0:links)=ZERO
   
   if(.not. lsingle_fluid)then
        
     rhoB(ix:mynx,iy:myny,iz:mynz)=MINDENS

      fuB(ix:mynx,iy:myny,iz:mynz)=ZERO
      fvB(ix:mynx,iy:myny,iz:mynz)=ZERO
      fwB(ix:mynx,iy:myny,iz:mynz)=ZERO

        if(lShanChen)then
         gradpsixB(ix:mynx,iy:myny,iz:mynz)=ZERO
         gradpsiyB(ix:mynx,iy:myny,iz:mynz)=ZERO
         gradpsizB(ix:mynx,iy:myny,iz:mynz)=ZERO
         psiB(ix:mynx,iy:myny,iz:mynz)=ZERO
        endif

      fluidB(ix:mynx,iy:myny,iz:mynz,0:links)=ZERO
   endif
   
   call set_double_range
   
   return
   
 end subroutine allocate_fluids
 
 subroutine driver_initialiaze_manage_bc_selfcomm()
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the initialization of 
!     the communication within the same process bewteen nodes linked
!     because of the periodic bc
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
  
  implicit none
  
  call initialiaze_manage_bc_pop_selfcomm
  
  return
  
 end subroutine driver_initialiaze_manage_bc_selfcomm
 
 subroutine initialize_fluids
 
!***********************************************************************
!     
!     LBsoft subroutine for initialize the fluids and read the
!     restart file if requested
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer :: i,j,k,l,imio3(3)
  logical :: ltestout
  character(len=32) :: itesto
  real(kind=PRC) :: dtemp(3)
  
! initialize isfluid
  
  select case(idistselect)
  case(1)
    call set_random_dens_fluids
  case(2)
    call set_uniform_dens_fluids
  case(3)
    call set_fake_dens_fluids
  case(4)
    call set_special_objects_dens_fluids
  case default
    call set_initial_dens_fluids
  end select
  
  call set_initial_vel_fluids
  
#ifdef DIAGNINIT
  call print_all_hvar(100,'inithvar_prima',0,rhoR,u,v,w)
#endif
  

  call commexch_dens(rhoR,rhoB)
  
  call manage_bc_hvar_selfcomm(rhoR)
  if(.not. lsingle_fluid)then
    call manage_bc_hvar_selfcomm(rhoB)
  endif
  
  call commwait_dens(rhoR,rhoB)
  
  
#ifdef DIAGNINIT
  call print_all_hvar(100,'inithvar_dopo',0,rhoR,u,v,w)
#endif
  
  call compute_densities_wall(0,lsingle_fluid,rhoR,rhoB)
  
#ifdef DIAGNINIT
  call print_all_hvar(100,'inithvar_fine',0,rhoR,u,v,w)
#endif
  
  !save initial local fluid mass
  !red fluid
  if(lsingle_fluid)then
    totnfluidnodes_init=0
    totnfluidmassR_init=ZERO
    ! forall(i=minx:maxx,j=miny:maxy,k=minz:maxz)
    do k=minz,maxz
      do j=miny,maxy
        do i=minx,maxx
          if(isfluid(i,j,k)==1) then
            totnfluidnodes_init=totnfluidnodes_init+1
            totnfluidmassR_init=totnfluidmassR_init+rhoR(i,j,k)
          endif
        enddo
      enddo
    enddo
    if(mxrank>1)then
      dtemp(1)=real(totnfluidnodes_init,kind=PRC)
      dtemp(2)=totnfluidmassR_init
      call sum_world_farr(dtemp,2)
      totnfluidnodes_init=nint(dtemp(1))
      totnfluidmassR_init=dtemp(2)
    endif
  else
  !both fluids
    totnfluidnodes_init=0
    totnfluidmassR_init=ZERO
    totnfluidmassB_init=ZERO
    ! forall(i=minx:maxx,j=miny:maxy,k=minz:maxz)
    do k=minz,maxz
      do j=miny,maxy
        do i=minx,maxx
          if(isfluid(i,j,k)==1) then
            totnfluidnodes_init=totnfluidnodes_init+1
            totnfluidmassR_init=totnfluidmassR_init+rhoR(i,j,k)
            totnfluidmassB_init=totnfluidmassB_init+rhoB(i,j,k)
          endif
        enddo
      enddo
    enddo
    if(mxrank>1)then
      dtemp(1)=real(totnfluidnodes_init,kind=PRC)
      dtemp(2)=totnfluidmassR_init
      dtemp(3)=totnfluidmassB_init
      call sum_world_farr(dtemp,3)
      totnfluidnodes_init=nint(dtemp(1))
      totnfluidmassR_init=dtemp(2)
      totnfluidmassB_init=dtemp(3)
    endif
  endif
  
  if(idistselect==3)then
    !perform a fake test with the density set as the real(i4) value
    call test_fake_dens(rhoB,ltestout)
    if(idrank==0)then
      if(ltestout)then
        write(6,*)'TEST: NO'
      else
        write(6,*)'TEST: OK'
      endif
    endif
    call finalize_world
    stop
  endif
  
  if(lvorticity)call driver_bc_velocities
  
  if(lunique_omega)then
    call set_unique_omega
  else
    call compute_omega
  endif
  
  call driver_set_initial_pop_fluids
  
  if(lShanChen)lexch_dens=.true.
  
  if(lColourG)lexch_dens=.true.
  
#ifdef DIAGNINIT
  call print_all_pops(100,'initpop',0,fluidR)
#endif
  
  !restart to be added
  
  return
  
 end subroutine initialize_fluids
 
 subroutine push_comm_isfluid()
 
!***********************************************************************
!     
!     LBsoft subroutine for pushing the isfluid communication
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
 
  implicit none
  
  !communicate isfluid over the processes applying the bc if necessary
  call comm_init_isfluid(isfluid)
  
  return
 
 end subroutine push_comm_isfluid
 
 subroutine initialize_fluid_force
 
!***********************************************************************
!     
!     LBsoft subroutine for initialize the force fluids
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer :: i,j,k
  
  if(.not. lforce_add)return
  
  forall(i=minx:maxx,j=miny:maxy,k=minz:maxz)fuR(i,j,k)=ext_fu
  forall(i=minx:maxx,j=miny:maxy,k=minz:maxz)fvR(i,j,k)=ext_fv
  forall(i=minx:maxx,j=miny:maxy,k=minz:maxz)fwR(i,j,k)=ext_fw
  
  if(lsingle_fluid)return
  
  forall(i=minx:maxx,j=miny:maxy,k=minz:maxz)fuB(i,j,k)=ext_fu
  forall(i=minx:maxx,j=miny:maxy,k=minz:maxz)fvB(i,j,k)=ext_fv
  forall(i=minx:maxx,j=miny:maxy,k=minz:maxz)fwB(i,j,k)=ext_fw
  
  return
 
 end subroutine initialize_fluid_force
 
 subroutine set_random_dens_fluids
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the initial density of fluids 
!     following a random gaussian distribution
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  integer :: i,j,k
  
  if(linit_seed)then
  
    do k=minz,maxz
      do j=miny,maxy
        do i=minx,maxx
          rhoR(i,j,k)=meanR+stdevR*gauss()
        enddo
      enddo
    enddo
    
    if(lsingle_fluid)return
    
    do k=minz,maxz
      do j=miny,maxy
        do i=minx,maxx
          rhoB(i,j,k)=meanB+stdevB*gauss()
        enddo
      enddo
    enddo
  
  else
  
    do k=minz,maxz
      do j=miny,maxy
        do i=minx,maxx
          rhoR(i,j,k)=meanR+stdevR*gauss_noseeded(i,j,k,1)
        enddo
      enddo
    enddo
    
    if(lsingle_fluid)return
    
    do k=minz,maxz
      do j=miny,maxy
        do i=minx,maxx
          rhoB(i,j,k)=meanB+stdevB*gauss_noseeded(i,j,k,100)
        enddo
      enddo
    enddo
  
  endif
  
  return
  
 end subroutine set_random_dens_fluids
 
 subroutine set_uniform_dens_fluids
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the initial density of fluids 
!     following a random uniform distribution
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  integer :: i,j,k
  
  real(kind=PRC) :: dtemp1
  
  if(linit_seed)then
  
    do k=1,nz
      do j=1,ny
        do i=1,nx
          call random_number(dtemp1)
          rhoR(i,j,k)=meanR+stdevR*dtemp1
        enddo
      enddo
    enddo
    
    if(lsingle_fluid)return
    
    do k=1,nz
      do j=1,ny
        do i=1,nx
          call random_number(dtemp1)
          rhoB(i,j,k)=meanB+stdevB*dtemp1
        enddo
      enddo
    enddo
  
  else
  
    forall (i=minx:maxx,j=miny:maxy,k=minz:maxz)
      rhoR(i,j,k)=meanR+stdevR*rand_noseeded(i,j,k,1)
    end forall
    
    if(lsingle_fluid)return
    
    forall (i=minx:maxx,j=miny:maxy,k=minz:maxz)
      rhoB(i,j,k)=meanB+stdevB*rand_noseeded(i,j,k,100)
    end forall
  
  endif
  
  return
  
 end subroutine set_uniform_dens_fluids
 
 subroutine set_special_objects_dens_fluids
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the initial density of fluids 
!     following a set of given objects
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
  
  implicit none
  
  integer :: i,j,k,l,iframe
  real(kind=PRC), dimension(3) :: dists,dbox
  real(kind=PRC) :: rdist,dtemp,xtemp,ytemp,ztemp
  
  dbox=[real(nx,kind=PRC),real(ny,kind=PRC),real(nz,kind=PRC)]
  
  if(lsingle_fluid)then
    do k=minz,maxz
      do j=miny,maxy
        do i=minx,maxx
          rhoR(i,j,k)=backR
        enddo
      enddo
    enddo
    do k=minz,maxz
      do j=miny,maxy
        do i=minx,maxx
          do l=1,nobjectliq
            select case(typeobjectliq(l))
            case(1)
              dists(1)=real(i,kind=PRC)-objectdata(1,l)
              dists(2)=real(j,kind=PRC)-objectdata(2,l)
              dists(3)=real(k,kind=PRC)-objectdata(3,l)
              call minimage_onevec(ibctype,dbox,dists)
              rdist=sqrt(dists(1)**TWO+dists(2)**TWO+dists(3)**TWO)
              dtemp=fcut(rdist,objectdata(4,l),objectdata(4,l)+objectdata(5,l))
              rhoR(i,j,k)=backR+(objectdata(6,l)-backR)*dtemp
            case(2)
              iframe=nint(objectdata(9,l))
              if((nint(objectdata(1,l))-iframe)<=i .and. (nint(objectdata(2,l))+iframe)>=i .and. &
               (nint(objectdata(3,l))-iframe)<=j .and. (nint(objectdata(4,l))+iframe)>=j .and. &
               (nint(objectdata(5,l))-iframe)<=k .and. (nint(objectdata(6,l))+iframe)>=k ) then
                if(nint(objectdata(1,l))<=i .and. nint(objectdata(2,l))>=i .and. &
                 nint(objectdata(3,l))<=j .and. nint(objectdata(4,l))>=j .and. &
                 nint(objectdata(5,l))<=k .and. nint(objectdata(6,l))>=k) then
                  rhoR(i,j,k)=objectdata(7,l)
                else
                  xtemp=ONE
                  ytemp=ONE
                  ztemp=ONE
                  if((nint(objectdata(1,l))-iframe)<=i .and. &
                   nint(objectdata(2,l))>i)then
                    rdist=objectdata(1,l)-real(i,kind=PRC)
                    xtemp=fcut(rdist,ZERO,objectdata(9,l))
                  endif
                  if((nint(objectdata(2,l))+iframe)>=i .and. &
                   i>nint(objectdata(1,l)))then
                    rdist=real(i,kind=PRC)-objectdata(2,l)
                    xtemp=fcut(rdist,ZERO,objectdata(9,l))
                  endif
                  if((nint(objectdata(3,l))-iframe)<=j .and. &
                   nint(objectdata(4,l))>j)then
                    rdist=objectdata(3,l)-real(j,kind=PRC)
                    ytemp=fcut(rdist,ZERO,objectdata(9,l))
                  endif
                  if((nint(objectdata(4,l))+iframe)>=j .and. &
                   j>nint(objectdata(3,l)))then
                    rdist=real(j,kind=PRC)-objectdata(4,l)
                    ytemp=fcut(rdist,ZERO,objectdata(9,l))
                  endif
                  if((nint(objectdata(5,l))-iframe)<=k .and. &
                   nint(objectdata(6,l))>k)then
                    rdist=objectdata(5,l)-real(k,kind=PRC)
                    ztemp=fcut(rdist,ZERO,objectdata(9,l))
                  endif
                  if((nint(objectdata(6,l))+iframe)>=k .and. &
                   k>nint(objectdata(5,l)))then
                    rdist=real(k,kind=PRC)-objectdata(6,l)
                    ztemp=fcut(rdist,ZERO,objectdata(9,l))
                  endif
                  dtemp=xtemp*ytemp*ztemp
                  rhoR(i,j,k)=backR+(objectdata(7,l)-backR)*dtemp
                endif
              endif
            end select
          enddo
        enddo
      enddo
    enddo
  else
    do k=minz,maxz
      do j=miny,maxy
        do i=minx,maxx
          rhoR(i,j,k)=backR
          rhoB(i,j,k)=backB
        enddo
      enddo
    enddo
    do k=minz,maxz
      do j=miny,maxy
        do i=minx,maxx
          do l=1,nobjectliq
            select case(typeobjectliq(l))
            case(1)
              dists(1)=real(i,kind=PRC)-objectdata(1,l)
              dists(2)=real(j,kind=PRC)-objectdata(2,l)
              dists(3)=real(k,kind=PRC)-objectdata(3,l)
              call minimage_onevec(ibctype,dbox,dists)
              rdist=sqrt(dists(1)**TWO+dists(2)**TWO+dists(3)**TWO)
              dtemp=fcut(rdist,objectdata(4,l),objectdata(4,l)+objectdata(5,l))
              rhoR(i,j,k)=backR+(objectdata(6,l)-backR)*dtemp
              rhoB(i,j,k)=backB+(objectdata(7,l)-backB)*dtemp
            case(2)
              iframe=nint(objectdata(9,l))
              if((nint(objectdata(1,l))-iframe)<=i .and. (nint(objectdata(2,l))+iframe)>=i .and. &
               (nint(objectdata(3,l))-iframe)<=j .and. (nint(objectdata(4,l))+iframe)>=j .and. &
               (nint(objectdata(5,l))-iframe)<=k .and. (nint(objectdata(6,l))+iframe)>=k ) then
                if(nint(objectdata(1,l))<=i .and. nint(objectdata(2,l))>=i .and. &
                 nint(objectdata(3,l))<=j .and. nint(objectdata(4,l))>=j .and. &
                 nint(objectdata(5,l))<=k .and. nint(objectdata(6,l))>=k) then
                  rhoR(i,j,k)=objectdata(7,l)
                  rhoB(i,j,k)=objectdata(8,l)
                else
                  xtemp=ONE
                  ytemp=ONE
                  ztemp=ONE
                  if((nint(objectdata(1,l))-iframe)<=i .and. &
                   nint(objectdata(2,l))>i)then
                    rdist=objectdata(1,l)-real(i,kind=PRC)
                    xtemp=fcut(rdist,ZERO,objectdata(9,l))
                  endif
                  if((nint(objectdata(2,l))+iframe)>=i .and. &
                   i>nint(objectdata(1,l)))then
                    rdist=real(i,kind=PRC)-objectdata(2,l)
                    xtemp=fcut(rdist,ZERO,objectdata(9,l))
                  endif
                  if((nint(objectdata(3,l))-iframe)<=j .and. &
                   nint(objectdata(4,l))>j)then
                    rdist=objectdata(3,l)-real(j,kind=PRC)
                    ytemp=fcut(rdist,ZERO,objectdata(9,l))
                  endif
                  if((nint(objectdata(4,l))+iframe)>=j .and. &
                   j>nint(objectdata(3,l)))then
                    rdist=real(j,kind=PRC)-objectdata(4,l)
                    ytemp=fcut(rdist,ZERO,objectdata(9,l))
                  endif
                  if((nint(objectdata(5,l))-iframe)<=k .and. &
                   nint(objectdata(6,l))>k)then
                    rdist=objectdata(5,l)-real(k,kind=PRC)
                    ztemp=fcut(rdist,ZERO,objectdata(9,l))
                  endif
                  if((nint(objectdata(6,l))+iframe)>=k .and. &
                   k>nint(objectdata(5,l)))then
                    rdist=real(k,kind=PRC)-objectdata(6,l)
                    ztemp=fcut(rdist,ZERO,objectdata(9,l))
                  endif
                  dtemp=xtemp*ytemp*ztemp
                  rhoR(i,j,k)=backR+(objectdata(7,l)-backR)*dtemp
                  rhoB(i,j,k)=backB+(objectdata(8,l)-backB)*dtemp
                endif
              endif
            end select
          enddo
        enddo
      enddo
    enddo
  endif
  
  return
  
 end subroutine set_special_objects_dens_fluids
 
 subroutine set_fake_dens_fluids
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the initial density of fluids 
!     following a fake distribution equal to the id fluid node
!     ONLY FOR DIAGNOSTIC PURPOSE
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification September 2018
!     
!***********************************************************************
  
  implicit none
  
  integer :: i,j,k,l,imio3(3)
  
  
  
  !check if the function i4back follows the natural order and if the
  !function i4find works correctly
  l=0
  do k=1-nbuff,nz+nbuff
    do j=1-nbuff,ny+nbuff
      do i=1-nbuff,nx+nbuff
        l=l+1
        !write(6,*)i4back(i,j,k),i,j,k,i4find(i4back(i,j,k))
        imio3=i4find(i4back(i,j,k))
        if(imio3(1)/=i.or.imio3(2)/=j.or.imio3(3)/=k.or.l/=i4back(i,j,k))then
          write(6,'(a,i8,a,3i8,a,3i8,a,i8)')'error in order ',l, &
          ' cord ',i,j,k,' i4find ',i4find(i4back(i,j,k)),' i4back ', &
          i4back(i,j,k)
          stop
        endif
      enddo
    enddo
  enddo
  
  !forall (i=minx:maxx,j=miny:maxy,k=minz:maxz)
  do k=minz,maxz
    do j=miny,maxy
      do i=minx,maxx
        rhoR(i,j,k)=real(i4back(i,j,k),kind=PRC)
      enddo
    enddo
  enddo
  !end forall
    
  if(lsingle_fluid)return
  
  !forall (i=minx:maxx,j=miny:maxy,k=minz:maxz)
  do k=minz,maxz
    do j=miny,maxy
      do i=minx,maxx
        rhoB(i,j,k)=real(i4back(i,j,k),kind=PRC)
      enddo
    enddo
  enddo
  !end forall
  
  return
  
 end subroutine set_fake_dens_fluids
 
 subroutine set_initial_dens_fluids
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the initial density of fluids 
!     equal to a given value
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  integer :: i,j,k
  
  forall(i=minx:maxx,j=miny:maxy,k=minz:maxz)rhoR(i,j,k)=meanR
  if(lsingle_fluid)return
  forall(i=minx:maxx,j=miny:maxy,k=minz:maxz)rhoB(i,j,k)=meanB
  
  return
  
 end subroutine set_initial_dens_fluids
 
 subroutine set_initial_vel_fluids
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the initial velocity of fluids 
!     equal to given values
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  integer :: i,j,k
  
  forall(i=minx:maxx,j=miny:maxy,k=minz:maxz)u(i,j,k)=initial_u
  forall(i=minx:maxx,j=miny:maxy,k=minz:maxz)v(i,j,k)=initial_v
  forall(i=minx:maxx,j=miny:maxy,k=minz:maxz)w(i,j,k)=initial_w
  
  return
  
 end subroutine set_initial_vel_fluids
 
 subroutine driver_set_initial_pop_fluids
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the setting of the initial 
!     fluid populations
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  integer :: i,j,k,l
  real*8 :: alphamean,uv,feqR,feqB
  real*8 :: locrhoR,locrhoB,locu,locv,locw
  
  if(lcolourG)then
    if(lunique_omega)then
      call set_initial_pop_fluids_CG_unique
    else
      call set_initial_pop_fluids_CG
    endif
    
  else
    
    call set_initial_pop_fluids(rhoR,u,v,w,fluidR)
    
    if(lsingle_fluid)return
    
    call set_initial_pop_fluids(rhoB,u,v,w,fluidB)
    
  endif
  
  return
  
 end subroutine driver_set_initial_pop_fluids
 
 subroutine set_initial_pop_fluids(rhosub,usub,vsub,wsub,fluidsub)
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the initial populations of fluids 
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC), allocatable, dimension(:,:,:)  :: rhosub,usub,vsub, &
   wsub
  real(kind=PRC), allocatable, dimension(:,:,:,:)  :: fluidsub
  
  integer :: i,j,k,l
  
  
  
  
  do k=minz-1,maxz+1
    do j=miny-1,maxy+1
      do i=minx-1,maxx+1
        fluidsub(i,j,k,0)=equil_pop00(rhosub(i,j,k),usub(i,j,k), &
         vsub(i,j,k),wsub(i,j,k))
        
        fluidsub(i,j,k,1)=equil_pop01(rhosub(i,j,k),usub(i,j,k), &
         vsub(i,j,k),wsub(i,j,k))
        
        fluidsub(i,j,k,2)=equil_pop02(rhosub(i,j,k),usub(i,j,k), &
         vsub(i,j,k),wsub(i,j,k))
        
        fluidsub(i,j,k,3)=equil_pop03(rhosub(i,j,k),usub(i,j,k), &
         vsub(i,j,k),wsub(i,j,k))
        
        fluidsub(i,j,k,4)=equil_pop04(rhosub(i,j,k),usub(i,j,k), &
         vsub(i,j,k),wsub(i,j,k))
        
        fluidsub(i,j,k,5)=equil_pop05(rhosub(i,j,k),usub(i,j,k), &
         vsub(i,j,k),wsub(i,j,k))
        
        fluidsub(i,j,k,6)=equil_pop06(rhosub(i,j,k),usub(i,j,k), &
         vsub(i,j,k),wsub(i,j,k))
        
        fluidsub(i,j,k,7)=equil_pop07(rhosub(i,j,k),usub(i,j,k), &
         vsub(i,j,k),wsub(i,j,k))
        
        fluidsub(i,j,k,8)=equil_pop08(rhosub(i,j,k),usub(i,j,k), &
         vsub(i,j,k),wsub(i,j,k))
        
        fluidsub(i,j,k,9)=equil_pop09(rhosub(i,j,k),usub(i,j,k), &
         vsub(i,j,k),wsub(i,j,k))
        
        fluidsub(i,j,k,10)=equil_pop10(rhosub(i,j,k),usub(i,j,k), &
         vsub(i,j,k),wsub(i,j,k))
        
        fluidsub(i,j,k,11)=equil_pop11(rhosub(i,j,k),usub(i,j,k), &
         vsub(i,j,k),wsub(i,j,k))
         
        fluidsub(i,j,k,12)=equil_pop12(rhosub(i,j,k),usub(i,j,k), &
         vsub(i,j,k),wsub(i,j,k))
         
        fluidsub(i,j,k,13)=equil_pop13(rhosub(i,j,k),usub(i,j,k), &
         vsub(i,j,k),wsub(i,j,k))
         
        fluidsub(i,j,k,14)=equil_pop14(rhosub(i,j,k),usub(i,j,k), &
         vsub(i,j,k),wsub(i,j,k))
         
        fluidsub(i,j,k,15)=equil_pop15(rhosub(i,j,k),usub(i,j,k), &
         vsub(i,j,k),wsub(i,j,k))
         
        fluidsub(i,j,k,16)=equil_pop16(rhosub(i,j,k),usub(i,j,k), &
         vsub(i,j,k),wsub(i,j,k))
         
        fluidsub(i,j,k,17)=equil_pop17(rhosub(i,j,k),usub(i,j,k), &
         vsub(i,j,k),wsub(i,j,k))
         
        fluidsub(i,j,k,18)=equil_pop18(rhosub(i,j,k),usub(i,j,k), &
         vsub(i,j,k),wsub(i,j,k))
        
      enddo
    enddo
  enddo
  
  return
  
 end subroutine set_initial_pop_fluids
 
 subroutine set_initial_pop_fluids_CG_unique()
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the initial populations of fluids
!     in the CG model with a unique omega
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification November 2019
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC) :: temp_omega
  real(kind=PRC) :: locrhoR,locrhoB,locu,locv,locw
  real(kind=PRC) :: grad_rhoRx,grad_rhoRy,grad_rhoRz
  real(kind=PRC) :: grad_rhoBx,grad_rhoBy,grad_rhoBz
  integer :: i,j,k,l
  
  
  temp_omega=unique_omega
  
   do k=minz,maxz
     do j=miny,maxy
       do i=minx,maxx
         if(isfluid(i,j,k)/=1)cycle
         locrhoR = rhoR(i,j,k)
         locrhoB = rhoB(i,j,k)
         locu = u(i,j,k)
         locv = v(i,j,k)
         locw = w(i,j,k)
         
#ifdef CG_CORRECTION
         
         grad_rhoRx = ZERO
         grad_rhoRy = ZERO
         grad_rhoRz = ZERO
           
         grad_rhoBx = ZERO
         grad_rhoBy = ZERO
         grad_rhoBz = ZERO
  
#ifdef GRADIENTD3Q27
         do l=1,linksd3q27
           grad_rhoRx=grad_rhoRx + ad3q27(l)*dexd3q27(l)* &
            (rhoR(i+exd3q27(l),j+eyd3q27(l),k+ezd3q27(l)))
           grad_rhoRy=grad_rhoRy + ad3q27(l)*deyd3q27(l)* &
            (rhoR(i+exd3q27(l),j+eyd3q27(l),k+ezd3q27(l)))
           grad_rhoRz=grad_rhoRz + ad3q27(l)*dezd3q27(l)* &
            (rhoR(i+exd3q27(l),j+eyd3q27(l),k+ezd3q27(l)))
              
          grad_rhoBx=grad_rhoBx + ad3q27(l)*dexd3q27(l)* &
            (rhoB(i+exd3q27(l),j+eyd3q27(l),k+ezd3q27(l)))
           grad_rhoBy=grad_rhoBy + ad3q27(l)*deyd3q27(l)* &
            (rhoB(i+exd3q27(l),j+eyd3q27(l),k+ezd3q27(l)))
           grad_rhoBz=grad_rhoBz + ad3q27(l)*dezd3q27(l)* &
            (rhoB(i+exd3q27(l),j+eyd3q27(l),k+ezd3q27(l)))
         enddo
#else             
         do l=1,links
           grad_rhoRx=grad_rhoRx + a(l)*dex(l)* &
            (rhoR(i+ex(l),j+ey(l),k+ez(l)))
           grad_rhoRy=grad_rhoRy + a(l)*dey(l)* &
            (rhoR(i+ex(l),j+ey(l),k+ez(l)))
           grad_rhoRz=grad_rhoRz + a(l)*dez(l)* &
            (rhoR(i+ex(l),j+ey(l),k+ez(l)))
             
           grad_rhoBx=grad_rhoBx + a(l)*dex(l)* &
            (rhoB(i+ex(l),j+ey(l),k+ez(l)))
           grad_rhoBy=grad_rhoBy + a(l)*dey(l)* &
            (rhoB(i+ex(l),j+ey(l),k+ez(l)))
           grad_rhoBz=grad_rhoBz + a(l)*dez(l)* &
            (rhoB(i+ex(l),j+ey(l),k+ez(l)))
         enddo
#endif    
         

         do l=0,links
           fluidR(i,j,k,l)= &
            equil_popCG_corr(l,temp_omega,alphaR_CG,locrhoR, &
            locu,locv,locw,grad_rhoRx,grad_rhoRy,grad_rhoRz)
           fluidB(i,j,k,l)= &
            equil_popCG_corr(l,temp_omega,alphaB_CG,locrhoB, &
            locu,locv,locw,grad_rhoBx,grad_rhoBy,grad_rhoBz)
         enddo
#else
         do l=0,links
           fluidR(i,j,k,l)= &
            equil_popCG(l,temp_omega,alphaR_CG,locrhoR,locu,locv,locw)
           fluidB(i,j,k,l)= &
            equil_popCG(l,temp_omega,alphaB_CG,locrhoB,locu,locv,locw)
         enddo
#endif
       enddo
     enddo
   enddo
  
  return
  
 end subroutine set_initial_pop_fluids_CG_unique
 
 subroutine set_initial_pop_fluids_CG()
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the initial populations of fluids
!     in the CG model 
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification November 2019
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC) :: temp_omega
  real(kind=PRC) :: locrhoR,locrhoB,locu,locv,locw
  real(kind=PRC) :: grad_rhoRx,grad_rhoRy,grad_rhoRz
  real(kind=PRC) :: grad_rhoBx,grad_rhoBy,grad_rhoBz
  integer :: i,j,k,l
  
  
   do k=minz,maxz
     do j=miny,maxy
       do i=minx,maxx
         if(isfluid(i,j,k)/=1)cycle
         temp_omega=omega(i,j,k)
         locrhoR = rhoR(i,j,k)
         locrhoB = rhoB(i,j,k)
         locu = u(i,j,k)
         locv = v(i,j,k)
         locw = w(i,j,k)
         
#ifdef CG_CORRECTION
         
           grad_rhoRx = ZERO
           grad_rhoRy = ZERO
           grad_rhoRz = ZERO
           
           grad_rhoBx = ZERO
           grad_rhoBy = ZERO
           grad_rhoBz = ZERO
  
#ifdef GRADIENTD3Q27
           do l=1,linksd3q27
             grad_rhoRx=grad_rhoRx + ad3q27(l)*dexd3q27(l)* &
              (rhoR(i+exd3q27(l),j+eyd3q27(l),k+ezd3q27(l)))
             grad_rhoRy=grad_rhoRy + ad3q27(l)*deyd3q27(l)* &
              (rhoR(i+exd3q27(l),j+eyd3q27(l),k+ezd3q27(l)))
             grad_rhoRz=grad_rhoRz + ad3q27(l)*dezd3q27(l)* &
              (rhoR(i+exd3q27(l),j+eyd3q27(l),k+ezd3q27(l)))
              
            grad_rhoBx=grad_rhoBx + ad3q27(l)*dexd3q27(l)* &
              (rhoB(i+exd3q27(l),j+eyd3q27(l),k+ezd3q27(l)))
             grad_rhoBy=grad_rhoBy + ad3q27(l)*deyd3q27(l)* &
              (rhoB(i+exd3q27(l),j+eyd3q27(l),k+ezd3q27(l)))
             grad_rhoBz=grad_rhoBz + ad3q27(l)*dezd3q27(l)* &
              (rhoB(i+exd3q27(l),j+eyd3q27(l),k+ezd3q27(l)))
           enddo
#else             
           do l=1,links
             grad_rhoRx=grad_rhoRx + a(l)*dex(l)* &
              (rhoR(i+ex(l),j+ey(l),k+ez(l)))
             grad_rhoRy=grad_rhoRy + a(l)*dey(l)* &
              (rhoR(i+ex(l),j+ey(l),k+ez(l)))
             grad_rhoRz=grad_rhoRz + a(l)*dez(l)* &
              (rhoR(i+ex(l),j+ey(l),k+ez(l)))
             
             grad_rhoBx=grad_rhoBx + a(l)*dex(l)* &
              (rhoB(i+ex(l),j+ey(l),k+ez(l)))
             grad_rhoBy=grad_rhoBy + a(l)*dey(l)* &
              (rhoB(i+ex(l),j+ey(l),k+ez(l)))
             grad_rhoBz=grad_rhoBz + a(l)*dez(l)* &
              (rhoB(i+ex(l),j+ey(l),k+ez(l)))
           enddo
#endif    


          do l=0,links
            fluidR(i,j,k,l)= &
             equil_popCG_corr(l,temp_omega,alphaR_CG,locrhoR, &
             locu,locv,locw,grad_rhoRx,grad_rhoRy,grad_rhoRz)
            fluidB(i,j,k,l)= &
             equil_popCG_corr(l,temp_omega,alphaB_CG,locrhoB, &
             locu,locv,locw,grad_rhoBx,grad_rhoBy,grad_rhoBz)
          enddo
#else
          do l=0,links
            fluidR(i,j,k,l)= &
             equil_popCG(l,temp_omega,alphaR_CG,locrhoR,locu,locv,locw)
            fluidB(i,j,k,l)= &
             equil_popCG(l,temp_omega,alphaB_CG,locrhoB,locu,locv,locw)
          enddo
#endif
       enddo
     enddo
   enddo
  
  return
  
 end subroutine set_initial_pop_fluids_CG
 
 subroutine set_mean_value_dens_fluids(dtemp1,dtemp2)
 
!***********************************************************************
!     
!     LBsoft subroutine for change the mean density of fluids 
!     inside this module
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
 
 subroutine set_back_value_dens_fluids(dtemp1,dtemp2)
 
!***********************************************************************
!     
!     LBsoft subroutine for change the background density of fluids 
!     inside this module
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification January 2019
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC), intent(in) :: dtemp1,dtemp2
  
  backR = dtemp1
  backB = dtemp2
  
  return
  
 end subroutine set_back_value_dens_fluids
 
 subroutine set_initial_value_vel_fluids(dtemp1,dtemp2,dtemp3)
 
!***********************************************************************
!     
!     LBsoft subroutine for change the initial velocity values of fluids 
!     inside this module
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
 
 subroutine set_value_ext_force_fluids(dtemp1,dtemp2,dtemp3)
 
!***********************************************************************
!     
!     LBsoft subroutine for set the external force values acting
!     on the fluids 
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC), intent(in) :: dtemp1,dtemp2,dtemp3
  
  lforce_add=.true.
  ext_fu = dtemp1
  ext_fv = dtemp2
  ext_fw = dtemp3
  
  return
  
 end subroutine set_value_ext_force_fluids
 
 subroutine set_LBintegrator_type(itemp1)
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the LB integrator type for
!     the collisional step
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  integer, intent(in) :: itemp1
  
  LBintegrator = itemp1
 end subroutine set_LBintegrator_type

 
 subroutine set_initial_dist_type(itemp1)
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the initial distribution type for
!     the density fluid initialization
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
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
 
 subroutine set_lbc_halfway(ltemp1)
  
!***********************************************************************
!     
!     LBsoft subroutine for set the value of lbc_halfway for select 
!     the halfway bounce back mode
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  logical, intent(in) :: ltemp1
  
  lbc_halfway=ltemp1
  lbc_fullway=(.not.ltemp1)
  
  return
  
 end subroutine set_lbc_halfway
 
 subroutine set_lbc_fullway(ltemp1)
  
!***********************************************************************
!     
!     LBsoft subroutine for set the value of lbc_fullway for select 
!     the fullway bounce back mode
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
  
  implicit none
  
  logical, intent(in) :: ltemp1
  
  lbc_fullway=ltemp1
  lbc_halfway=(.not.ltemp1)
  
  call error(41)
  
  return
  
 end subroutine set_lbc_fullway
 
 subroutine set_lsingle_fluid(ltemp1)
 
!***********************************************************************
!     
!     LBsoft subroutine for set the value of lsingle_fluid for select 
!     the single fluid mode
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  logical, intent(in) :: ltemp1
  
  lsingle_fluid=ltemp1
  
  return
  
 end subroutine set_lsingle_fluid
 
 subroutine set_fluid_force_sc(ltemp1,dtemp1)
 
!***********************************************************************
!     
!     LBsoft subroutine for set the value of the pair ShanChen 
!     force constant
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  logical, intent(in) :: ltemp1
  real(kind=PRC), intent(in) :: dtemp1
  
  lforce_add=.true.
  lpair_SC = ltemp1
  lShanChen = (lShanChen .or. lpair_SC)
  pair_SC = dtemp1
  
  return
  
 end subroutine set_fluid_force_sc
 
 subroutine set_fluid_force_cg(ltemp1,dtemp1,dtemp2,dtemp3,dtemp4)
 
!***********************************************************************
!     
!     LBsoft subroutine for set the value of the surface tension 
!     force in the colour gradient model
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification November 2019
!     
!***********************************************************************
  
  implicit none
  
  logical, intent(in) :: ltemp1
  real(kind=PRC), intent(in) :: dtemp1,dtemp2,dtemp3,dtemp4
  
  lColourG = ltemp1
  sigma_CG = dtemp1
  alphaR_CG = dtemp2
  alphaB_CG = dtemp3
  beta_CG = dtemp4
  
  call initialize_phi_coeff_cg(alphaR_CG,alphaB_CG)
  
  return
  
 end subroutine set_fluid_force_cg
 
 subroutine initialize_phi_coeff_cg(dtemp1,dtemp2)
  
!***********************************************************************
!     
!     LBsoft subroutine for set the value of the phi coefficients
!     of the lattice in the colour gradient model
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification November 2019
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC), intent(in) :: dtemp1,dtemp2
  
  phiR_cg(0:links)=phi_CG(0:links)+dtemp1*varphi_CG(0:links)
  phiB_cg(0:links)=phi_CG(0:links)+dtemp2*varphi_CG(0:links)
  
  return
  
 end subroutine initialize_phi_coeff_cg
 
 subroutine set_fluid_wall_sc(dtemp1,dtemp2)
 
!***********************************************************************
!     
!     LBsoft subroutine for set the value of the pair ShanChen 
!     wall coupling constant
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC), intent(in) :: dtemp1,dtemp2
  
  wallR_SC = dtemp1
  wallB_SC = dtemp2
  lwall_SC = (wallR_SC/=ONE .or. wallB_SC/=ONE )
  
  return
  
 end subroutine set_fluid_wall_sc
 
 subroutine set_fluid_particle_sc(dtemp1,dtemp2)
 
!***********************************************************************
!     
!     LBsoft subroutine for set the value of the pair ShanChen 
!     particle coupling constant
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification January 2019
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC), intent(in) :: dtemp1,dtemp2
  
  partR_SC = dtemp1
  partB_SC = dtemp2
  lpart_SC = (partR_SC/=ONE .or. partB_SC/=ONE )
  
  return
  
 end subroutine set_fluid_particle_sc
 
 subroutine set_theta_particle_sc(dtemp1,dtemp2)
 
!***********************************************************************
!     
!     LBsoft subroutine for set the value of the contact angle 
!     in the particle Shan Chen coupling
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification January 2019
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC), intent(in) :: dtemp1,dtemp2
  
  theta_SC = dtemp1
  devtheta_SC = dtemp2
  
  return
  
 end subroutine set_theta_particle_sc
 
 subroutine set_value_viscosity(dtemp1,dtemp2,temp_lsingle_fluid)
 
!***********************************************************************
!     
!     LBsoft subroutine for set the viscosity value and the relaxation 
!     time tau
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC), intent(in) :: dtemp1,dtemp2
  logical, intent(in) :: temp_lsingle_fluid
  
  if(dtemp1==dtemp2 .or. temp_lsingle_fluid)then
    lunique_omega=.true.
    viscR=dtemp1
    viscB=dtemp2
    tauR=viscR/cssq+HALF
    tauB=viscB/cssq+HALF
    unique_omega = viscosity_to_omega(viscR)
  else
    lunique_omega=.false.
    viscR=dtemp1
    viscB=dtemp2
    tauR=viscR/cssq+HALF
    tauB=viscB/cssq+HALF
    
  endif
  
  return
  
 end subroutine set_value_viscosity
 
 subroutine set_value_tau(dtemp1,dtemp2,temp_lsingle_fluid)
 
!***********************************************************************
!     
!     LBsoft subroutine for set the relaxation time tau value and 
!     the viscosity
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC), intent(in) :: dtemp1,dtemp2
  logical, intent(in) :: temp_lsingle_fluid
  
  if(dtemp1==dtemp2 .or. temp_lsingle_fluid)then
    lunique_omega=.true.
    tauR=dtemp1
    tauB=dtemp2
    viscR=cssq*(tauR-HALF)
    viscB=cssq*(tauB-HALF)
    unique_omega = viscosity_to_omega(viscR)
  else
    lunique_omega=.false.
    tauR=dtemp1
    tauB=dtemp2
    viscR=cssq*(tauR-HALF)
    viscB=cssq*(tauB-HALF)
    if(iselomegaint==1)then
      s1_R=TWO*(tauR*tauB)/(tauR+tauB)
      s2_R=TWO*(tauR-s1_R)/Delta_Grunau
      s3_R=-s2_R/(TWO*Delta_Grunau)
      t1_B=TWO*(tauR*tauB)/(tauR+tauB)
      t2_B=TWO*(t1_B-tauB)/Delta_Grunau
      t3_B=t2_B/(TWO*Delta_Grunau)
    endif
  endif
  
  return
  
 end subroutine set_value_tau
 
 subroutine set_unique_omega
 
!***********************************************************************
!     
!     LBsoft subroutine for set a unique omega value
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  integer :: i,j,k
  
  forall(i=minx:maxx,j=miny:maxy,k=minz:maxz)
    omega(i,j,k)=unique_omega
  end forall
  
  return
  
 end subroutine set_unique_omega
 
 subroutine set_cap_force(ltemp1,dtemp1)
 
!***********************************************************************
!     
!     LBsoft subroutine for set the capping force value acting
!     on fluids 
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2019
!     
!***********************************************************************
  
  implicit none
  
  logical, intent(in) :: ltemp1
  real(kind=PRC), intent(in) :: dtemp1
  
  lcap_force = ltemp1
  cap_force = dtemp1
  
  return
  
 end subroutine set_cap_force
 
 subroutine set_mass_rescale(ltemp1,itemp1,dtemp1)
 
!***********************************************************************
!     
!     LBsoft subroutine for set the rescaling of fluid mass
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2019
!     
!***********************************************************************
  
  implicit none
  
  logical, intent(in) :: ltemp1
  integer, intent(in) :: itemp1
  real(kind=PRC), intent(in) :: dtemp1
  
  lmass_rescale = ltemp1
  imass_rescale = itemp1
  dmass_rescale = dtemp1
  
  return
  
 end subroutine set_mass_rescale
 
 subroutine set_objectliq(ltemp1,itemp1)
 
!***********************************************************************
!     
!     LBsoft subroutine for set the number and type of objects in  
!     the special fluid initialization
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
  
  implicit none
  
  logical, intent(in) :: ltemp1
  integer, intent(in) :: itemp1
  
  lobjectliq=ltemp1
  nobjectliq=itemp1
  
  if(allocated(typeobjectliq))deallocate(typeobjectliq)
  if(allocated(objectdata))deallocate(objectdata)
  allocate(typeobjectliq(nobjectliq))
  allocate(objectdata(9,nobjectliq))
  
  return
  
 end subroutine set_objectliq
 
 subroutine set_objectdata(itemp1,itemp2,dtemp1)
 
!***********************************************************************
!     
!     LBsoft subroutine for set the value of objectdata for special 
!     fluid initialization
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: itemp1
  integer, intent(in), dimension(itemp1) :: itemp2
  real(kind=PRC), intent(in), dimension(9,itemp1) :: dtemp1
    
  typeobjectliq(1:itemp1)=itemp2(1:itemp1)
  
  objectdata(1:9,1:itemp1)=dtemp1(1:9,1:itemp1)
  
  return
  
 end subroutine set_objectdata
 
 pure function viscosity_to_omega(dtemp1)
 
!***********************************************************************
!     
!     LBsoft function to convert the viscosity to omega
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: dtemp1
  
  real(kind=PRC) :: viscosity_to_omega
  
  viscosity_to_omega = ONE / ( dtemp1 /cssq  + HALF)
  
  return
 
 end function viscosity_to_omega
 
 pure function omega_to_viscosity(dtemp1)
 
!***********************************************************************
!     
!     LBsoft function to convert the omega to viscosity
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification December 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: dtemp1
  
  real(kind=PRC) :: omega_to_viscosity
  
  omega_to_viscosity = cssq * ( ONE / dtemp1 - HALF )
  
  return
 
 end function omega_to_viscosity
 
 subroutine rescale_fluid_mass(nstep)
 
!***********************************************************************
!     
!     LBsoft subroutine for rescaling total fluid mass
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2019
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: nstep
  
  real(kind=PRC) :: dtemp(3),rescalef
  real(kind=PRC) :: totnfluidmassR,totnfluidmassB,totnfluidnodes
  
  if(mod(nstep,imass_rescale)/=0)return
  
  if(lsingle_fluid)then
    if(mxrank>1)then
      dtemp(1)=real(nfluidnodes,kind=PRC)
      dtemp(2)=nfluidmassR
      call sum_world_farr(dtemp,2)
      totnfluidnodes=dtemp(1)
      totnfluidmassR=dtemp(2)
    else
      totnfluidnodes=real(nfluidnodes,kind=PRC)
      totnfluidmassR=nfluidmassR
    endif
    !check the deviation
    if((abs(totnfluidmassR-totnfluidmassR_init)/totnfluidmassR_init)> &
     dmass_rescale)then
      rescalef=(totnfluidmassR_init*totnfluidnodes)/ &
       (totnfluidmassR*real(totnfluidnodes_init,kind=PRC))
      call rescale_pops(rescalef,fluidR)
    endif
  else
    if(mxrank>1)then
      dtemp(1)=real(nfluidnodes,kind=PRC)
      dtemp(2)=nfluidmassR
      dtemp(3)=nfluidmassB
      call sum_world_farr(dtemp,3)
      totnfluidnodes=dtemp(1)
      totnfluidmassR=dtemp(2)
      totnfluidmassB=dtemp(3)
    else
      totnfluidnodes=real(nfluidnodes,kind=PRC)
      totnfluidmassR=nfluidmassR
      totnfluidmassB=nfluidmassB
    endif
    !check the deviation
    if((abs(totnfluidmassR-totnfluidmassR_init)/totnfluidmassR_init)> &
     dmass_rescale)then
      rescalef=(totnfluidmassR_init*totnfluidnodes)/ &
       (totnfluidmassR*real(totnfluidnodes_init,kind=PRC))
      call rescale_pops(rescalef,fluidR)
    endif
    !check the deviation
    if((abs(totnfluidmassB-totnfluidmassB_init)/totnfluidmassB_init)> &
     dmass_rescale)then
      rescalef=(totnfluidmassB_init*totnfluidnodes)/ &
       (totnfluidmassB*real(totnfluidnodes_init,kind=PRC))
      call rescale_pops(rescalef,fluidB)
    endif
    
  endif
  
  return
  
 end subroutine rescale_fluid_mass
 
 subroutine rescale_pops(rescalesub,fluidsub)
 
!***********************************************************************
!     
!     LBsoft subroutine for scaling all the populations of a factor
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2019
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC), intent(in) :: rescalesub
  real(kind=PRC), allocatable, dimension(:,:,:,:)  :: fluidsub
  
  integer :: i,j,k,l
  
  do l=0,links
    do k=minz,maxz
      do j=miny,maxy
        do i=minx,maxx
          if(isfluid(i,j,k)/=1)cycle
          fluidsub(i,j,k,l)=fluidsub(i,j,k,l)*rescalesub
        enddo
      enddo
    enddo
  enddo

  return
  
 end subroutine rescale_pops
 
 subroutine convert_fluid_force_to_velshifted
 
!***********************************************************************
!     
!     LBsoft subroutine for converting the forces to the correspondent
!     velocity shifted for the forces
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer :: i,j,k
  
  !red fluid
  !forall(i=minx:maxx,j=miny:maxy,k=minz:maxz,isfluid(i,j,k)==1)
  do k=minz,maxz
   do j=miny,maxy
    do i=minx,maxx

    if (isfluid(i,j,k)/=1) cycle
    fuR(i,j,k) = fuR(i,j,k)*t_LB / rhoR(i,j,k) + u(i,j,k)
    fvR(i,j,k) = fvR(i,j,k)*t_LB / rhoR(i,j,k) + v(i,j,k)
    fwR(i,j,k) = fwR(i,j,k)*t_LB / rhoR(i,j,k) + w(i,j,k)
  enddo
    enddo
     enddo
  !end forall
  
  if(lsingle_fluid)return
  
  !blue fluid
  !forall(i=minx:maxx,j=miny:maxy,k=minz:maxz,isfluid(i,j,k)==1)
  do k=minz,maxz
   do j=miny,maxy
    do i=minx,maxx

    if (isfluid(i,j,k)/=1) cycle

    fuB(i,j,k) = fuB(i,j,k)*t_LB / rhoB(i,j,k) + u(i,j,k)
    fvB(i,j,k) = fvB(i,j,k)*t_LB / rhoB(i,j,k) + v(i,j,k)
    fwB(i,j,k) = fwB(i,j,k)*t_LB / rhoB(i,j,k) + w(i,j,k)
  enddo
    enddo
     enddo
  !end forall
  
  return
  
 end subroutine convert_fluid_force_to_velshifted
 
 subroutine driver_collision_fluids(nstep)
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the collisional step
!     on the Boltzmann populations 
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  integer, intent(in) :: nstep
  
  
#ifdef DIAGNHVAR
  integer, save :: iter=0
  
  iter=iter+1
  if(iter.eq.NDIAGNHVAR)call print_all_hvar(1320,'testhvar',iter, &
   rhoR,u,v,w)
#endif
  
  if(.not. lColourG)then
  
    select case(LBintegrator)
    case (0)
      call collision_fluids_BGK(rhoR,u,v,w,omega,fluidR,nstep)
      if(lsingle_fluid)return
      call collision_fluids_BGK(rhoB,u,v,w,omega,fluidB,nstep)

    case (1)
      ! done in collision_EDM) call convert_fluid_force_to_velshifted
      call collision_fluids_EDM(rhoR,u,v,w,fuR,fvR,fwR,omega,fluidR,nstep)
      if(lsingle_fluid)return
      call collision_fluids_EDM(rhoB,u,v,w,fuB,fvB,fwB,omega,fluidB,nstep)
      call update_velocity_EDM(rhoR,rhoB,fuR,fvR,fwR,fuB,fvB,fwB, &
       u,v,w,nstep)
    case default
      call error(14)
    end select
  
  else
    
    call driver_collision_fluids_CG(nstep)
    
  endif
  
  
 end subroutine driver_collision_fluids

 
 subroutine collision_fluids_BGK(rhosub,usub,vsub,wsub,omegas, &
  fluidsub,nstep)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the collisional step
!     on the Boltzmann populations 
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: nstep
  real(kind=PRC), allocatable, dimension(:,:,:)  :: rhosub,usub,vsub, &
   wsub,omegas
  real(kind=PRC), allocatable, dimension(:,:,:,:)  :: fluidsub
  integer :: l
  real(kind=PRC) :: locrho,locu,locv,locw, oneminusomega, temp_omega
  real(kind=PRC), dimension(0:links) :: feq
#ifdef REGULARIZED
  real(kind=PRC) :: pxx,pyy,pzz,pxy,pxz,pyz
  real(kind=PRC) :: fneq
#endif
  
  integer :: i,j,k
  
  ! forall(i=minx:maxx,j=miny:maxy,k=minz:maxz,isfluid(i,j,k)==1)
  if (lunique_omega) then
      temp_omega=unique_omega
      oneminusomega = ONE-temp_omega
      do k=minz,maxz
       do j=miny,maxy
        do i=minx,maxx

        if (isfluid(i,j,k)/=1) cycle

        locrho = rhosub(i,j,k)
        locu = usub(i,j,k)
        locv = vsub(i,j,k)
        locw = wsub(i,j,k)
        
#ifdef REGULARIZED
        pxx = ZERO
        pyy = ZERO
        pzz = ZERO
		pxy = ZERO
        pxz = ZERO
        pyz = ZERO
          
        !compute the equilibrium
        feq(0) = equil_pop00(locrho,locu,locv,locw)
        feq(1) = equil_pop01(locrho,locu,locv,locw)
        feq(2) = equil_pop02(locrho,locu,locv,locw)
        feq(3) = equil_pop03(locrho,locu,locv,locw)
        feq(4) = equil_pop04(locrho,locu,locv,locw)
        feq(5) = equil_pop05(locrho,locu,locv,locw)
        feq(6) = equil_pop06(locrho,locu,locv,locw)
        feq(7) = equil_pop07(locrho,locu,locv,locw)
        feq(8) = equil_pop08(locrho,locu,locv,locw)
        feq(9) = equil_pop09(locrho,locu,locv,locw)
        feq(10) = equil_pop10(locrho,locu,locv,locw)
        feq(11) = equil_pop11(locrho,locu,locv,locw)
        feq(12) = equil_pop12(locrho,locu,locv,locw)
        feq(13) = equil_pop13(locrho,locu,locv,locw)
        feq(14) = equil_pop14(locrho,locu,locv,locw)
        feq(15) = equil_pop15(locrho,locu,locv,locw)
        feq(16) = equil_pop16(locrho,locu,locv,locw)
        feq(17) = equil_pop17(locrho,locu,locv,locw)
        feq(18) = equil_pop18(locrho,locu,locv,locw)

        do l=0,links
          !compute the equilibrium for the two components
          fneq = fluidsub(i,j,k,l) - feq(l)
			
	      !non equilibrium part of the momentum flux tensor
          pxx=pxx + (dex(l)*dex(l) - cssq)*fneq
          pyy=pyy + (dey(l)*dey(l) - cssq)*fneq
          pzz=pzz + (dez(l)*dez(l) - cssq)*fneq
          pxy=pxy +  dex(l)*dey(l)*fneq
          pxz=pxz +  dex(l)*dez(l)*fneq
          pyz=pyz +  dey(l)*dez(l)*fneq
        enddo
        
        do l=0,links
          fluidsub(i,j,k,l)= feq(l) + &
           ((HALF*p(l))/(cssq**TWO))*((dex(l)*dex(l) - cssq)*pxx + &
           (dey(l)*dey(l) - cssq)*pyy + &
           (dez(l)*dez(l) - cssq)*pzz + &
           TWO*dex(l)*dey(l)*pxy + &
           TWO*dex(l)*dez(l)*pxz + &
           TWO*dey(l)*dez(l)*pyz)           
        enddo
        
#else
        
        !compute the equilibrium
        feq(0) = equil_pop00(locrho,locu,locv,locw)
        feq(1) = equil_pop01(locrho,locu,locv,locw)
        feq(2) = equil_pop02(locrho,locu,locv,locw)
        feq(3) = equil_pop03(locrho,locu,locv,locw)
        feq(4) = equil_pop04(locrho,locu,locv,locw)
        feq(5) = equil_pop05(locrho,locu,locv,locw)
        feq(6) = equil_pop06(locrho,locu,locv,locw)
        feq(7) = equil_pop07(locrho,locu,locv,locw)
        feq(8) = equil_pop08(locrho,locu,locv,locw)
        feq(9) = equil_pop09(locrho,locu,locv,locw)
        feq(10) = equil_pop10(locrho,locu,locv,locw)
        feq(11) = equil_pop11(locrho,locu,locv,locw)
        feq(12) = equil_pop12(locrho,locu,locv,locw)
        feq(13) = equil_pop13(locrho,locu,locv,locw)
        feq(14) = equil_pop14(locrho,locu,locv,locw)
        feq(15) = equil_pop15(locrho,locu,locv,locw)
        feq(16) = equil_pop16(locrho,locu,locv,locw)
        feq(17) = equil_pop17(locrho,locu,locv,locw)
        feq(18) = equil_pop18(locrho,locu,locv,locw)
        
#endif        
        
        do l=0,links
          fluidsub(i,j,k,l)=fluidsub(i,j,k,l)*oneminusomega+feq(l)*temp_omega
        enddo

      enddo
     enddo
    enddo

    return
  endif

  ! Variable omega

  do k=minz,maxz
   do j=miny,maxy
    do i=minx,maxx

    if (isfluid(i,j,k)/=1) cycle
    
    temp_omega=omegas(i,j,k)
    oneminusomega = ONE-temp_omega
    
#ifdef REGULARIZED
        pxx = ZERO
        pyy = ZERO
        pzz = ZERO
		pxy = ZERO
        pxz = ZERO
        pyz = ZERO
          
        !compute the equilibrium
        feq(0) = equil_pop00(locrho,locu,locv,locw)
        feq(1) = equil_pop01(locrho,locu,locv,locw)
        feq(2) = equil_pop02(locrho,locu,locv,locw)
        feq(3) = equil_pop03(locrho,locu,locv,locw)
        feq(4) = equil_pop04(locrho,locu,locv,locw)
        feq(5) = equil_pop05(locrho,locu,locv,locw)
        feq(6) = equil_pop06(locrho,locu,locv,locw)
        feq(7) = equil_pop07(locrho,locu,locv,locw)
        feq(8) = equil_pop08(locrho,locu,locv,locw)
        feq(9) = equil_pop09(locrho,locu,locv,locw)
        feq(10) = equil_pop10(locrho,locu,locv,locw)
        feq(11) = equil_pop11(locrho,locu,locv,locw)
        feq(12) = equil_pop12(locrho,locu,locv,locw)
        feq(13) = equil_pop13(locrho,locu,locv,locw)
        feq(14) = equil_pop14(locrho,locu,locv,locw)
        feq(15) = equil_pop15(locrho,locu,locv,locw)
        feq(16) = equil_pop16(locrho,locu,locv,locw)
        feq(17) = equil_pop17(locrho,locu,locv,locw)
        feq(18) = equil_pop18(locrho,locu,locv,locw)

        do l=0,links
          !compute the equilibrium for the two components
          fneq = fluidsub(i,j,k,l) - feq(l)
			
	      !non equilibrium part of the momentum flux tensor
          pxx=pxx + (dex(l)*dex(l) - cssq)*fneq
          pyy=pyy + (dey(l)*dey(l) - cssq)*fneq
          pzz=pzz + (dez(l)*dez(l) - cssq)*fneq
          pxy=pxy +  dex(l)*dey(l)*fneq
          pxz=pxz +  dex(l)*dez(l)*fneq
          pyz=pyz +  dey(l)*dez(l)*fneq
        enddo
        
        do l=0,links
          fluidsub(i,j,k,l)= feq(l) + &
           ((HALF*p(l))/(cssq**TWO))*((dex(l)*dex(l) - cssq)*pxx + &
           (dey(l)*dey(l) - cssq)*pyy + &
           (dez(l)*dez(l) - cssq)*pzz + &
           TWO*dex(l)*dey(l)*pxy + &
           TWO*dex(l)*dez(l)*pxz + &
           TWO*dey(l)*dez(l)*pyz)           
        enddo
        
#else
        
        !compute the equilibrium
        feq(0) = equil_pop00(locrho,locu,locv,locw)
        feq(1) = equil_pop01(locrho,locu,locv,locw)
        feq(2) = equil_pop02(locrho,locu,locv,locw)
        feq(3) = equil_pop03(locrho,locu,locv,locw)
        feq(4) = equil_pop04(locrho,locu,locv,locw)
        feq(5) = equil_pop05(locrho,locu,locv,locw)
        feq(6) = equil_pop06(locrho,locu,locv,locw)
        feq(7) = equil_pop07(locrho,locu,locv,locw)
        feq(8) = equil_pop08(locrho,locu,locv,locw)
        feq(9) = equil_pop09(locrho,locu,locv,locw)
        feq(10) = equil_pop10(locrho,locu,locv,locw)
        feq(11) = equil_pop11(locrho,locu,locv,locw)
        feq(12) = equil_pop12(locrho,locu,locv,locw)
        feq(13) = equil_pop13(locrho,locu,locv,locw)
        feq(14) = equil_pop14(locrho,locu,locv,locw)
        feq(15) = equil_pop15(locrho,locu,locv,locw)
        feq(16) = equil_pop16(locrho,locu,locv,locw)
        feq(17) = equil_pop17(locrho,locu,locv,locw)
        feq(18) = equil_pop18(locrho,locu,locv,locw)
        
#endif
        
        do l=0,links
          fluidsub(i,j,k,l)=fluidsub(i,j,k,l)*oneminusomega+feq(l)*temp_omega
        enddo

      enddo
    enddo
  enddo

  ! end forall
 end subroutine collision_fluids_BGK
 

 subroutine collision_fluids_EDM(rhosub,usub,vsub,wsub,fusub,fvsub, &
  fwsub,omegas,fluidsub,nstep)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the collisional step
!     on the Boltzmann populations 
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  integer, intent(in) :: nstep
  
  real(kind=PRC), allocatable, dimension(:,:,:)  :: rhosub,usub,vsub, &
   wsub,fusub,fvsub,fwsub,omegas
  real(kind=PRC), allocatable, dimension(:,:,:,:)  :: fluidsub
  real(kind=PRC) :: locrho,locu,locv,locw, oneminusomega, &
        locfu,locfv,locfw,temp_omega,omegaminusone
  
  integer :: i,j,k,l
  
#ifdef REGULARIZED
  real(kind=PRC) :: pxx,pyy,pzz,pxy,pxz,pyz
  real(kind=PRC) :: fneq
#endif
  real(kind=PRC), dimension(0:links) :: feq,feqshift
  
  
  if (lunique_omega) then
    temp_omega=unique_omega
    oneminusomega = ONE-temp_omega
    omegaminusone = temp_omega-ONE

    do k=minz,maxz
      do j=miny,maxy
        do i=minx,maxx

          if (isfluid(i,j,k)/=1) cycle

          locrho = rhosub(i,j,k)
          locu = usub(i,j,k)
          locv = vsub(i,j,k)
          locw = wsub(i,j,k)
        

          locfu = fusub(i,j,k)*t_LB / locrho + locu
          locfv = fvsub(i,j,k)*t_LB / locrho + locv
          locfw = fwsub(i,j,k)*t_LB / locrho + locw

          
          
#ifdef REGULARIZED
          pxx = ZERO
          pyy = ZERO
          pzz = ZERO
 	      pxy = ZERO
          pxz = ZERO
          pyz = ZERO
           
          !compute the equilibrium
          feq(0) = equil_pop00(locrho,locu,locv,locw)
          feq(1) = equil_pop01(locrho,locu,locv,locw)
          feq(2) = equil_pop02(locrho,locu,locv,locw)
          feq(3) = equil_pop03(locrho,locu,locv,locw)
          feq(4) = equil_pop04(locrho,locu,locv,locw)
          feq(5) = equil_pop05(locrho,locu,locv,locw)
          feq(6) = equil_pop06(locrho,locu,locv,locw)
          feq(7) = equil_pop07(locrho,locu,locv,locw)
          feq(8) = equil_pop08(locrho,locu,locv,locw)
          feq(9) = equil_pop09(locrho,locu,locv,locw)
          feq(10) = equil_pop10(locrho,locu,locv,locw)
          feq(11) = equil_pop11(locrho,locu,locv,locw)
          feq(12) = equil_pop12(locrho,locu,locv,locw)
          feq(13) = equil_pop13(locrho,locu,locv,locw)
          feq(14) = equil_pop14(locrho,locu,locv,locw)
          feq(15) = equil_pop15(locrho,locu,locv,locw)
          feq(16) = equil_pop16(locrho,locu,locv,locw)
          feq(17) = equil_pop17(locrho,locu,locv,locw)
          feq(18) = equil_pop18(locrho,locu,locv,locw)
          
          do l=0,links
            !compute the equilibrium for the two components
            fneq = fluidsub(i,j,k,l) - feq(l)
		    
	        !non equilibrium part of the momentum flux tensor
            pxx=pxx + (dex(l)*dex(l) - cssq)*fneq
            pyy=pyy + (dey(l)*dey(l) - cssq)*fneq
            pzz=pzz + (dez(l)*dez(l) - cssq)*fneq
            pxy=pxy +  dex(l)*dey(l)*fneq
            pxz=pxz +  dex(l)*dez(l)*fneq
            pyz=pyz +  dey(l)*dez(l)*fneq
          enddo
        
          do l=0,links
            fluidsub(i,j,k,l)= feq(l) + &
             ((HALF*p(l))/(cssq**TWO))*((dex(l)*dex(l) - cssq)*pxx + &
             (dey(l)*dey(l) - cssq)*pyy + &
             (dez(l)*dez(l) - cssq)*pzz + &
             TWO*dex(l)*dey(l)*pxy + &
             TWO*dex(l)*dez(l)*pxz + &
             TWO*dey(l)*dez(l)*pyz)           
          enddo
          
#else
          
          !compute the equilibrium
          feq(0) = equil_pop00(locrho,locu,locv,locw)
          feq(1) = equil_pop01(locrho,locu,locv,locw)
          feq(2) = equil_pop02(locrho,locu,locv,locw)
          feq(3) = equil_pop03(locrho,locu,locv,locw)
          feq(4) = equil_pop04(locrho,locu,locv,locw)
          feq(5) = equil_pop05(locrho,locu,locv,locw)
          feq(6) = equil_pop06(locrho,locu,locv,locw)
          feq(7) = equil_pop07(locrho,locu,locv,locw)
          feq(8) = equil_pop08(locrho,locu,locv,locw)
          feq(9) = equil_pop09(locrho,locu,locv,locw)
          feq(10) = equil_pop10(locrho,locu,locv,locw)
          feq(11) = equil_pop11(locrho,locu,locv,locw)
          feq(12) = equil_pop12(locrho,locu,locv,locw)
          feq(13) = equil_pop13(locrho,locu,locv,locw)
          feq(14) = equil_pop14(locrho,locu,locv,locw)
          feq(15) = equil_pop15(locrho,locu,locv,locw)
          feq(16) = equil_pop16(locrho,locu,locv,locw)
          feq(17) = equil_pop17(locrho,locu,locv,locw)
          feq(18) = equil_pop18(locrho,locu,locv,locw)
          
#endif          
          
          
          !compute the shifted equilibrium
          feqshift(0) = equil_pop00(locrho,locfu,locfv,locfw)
          feqshift(1) = equil_pop01(locrho,locfu,locfv,locfw)
          feqshift(2) = equil_pop02(locrho,locfu,locfv,locfw)
          feqshift(3) = equil_pop03(locrho,locfu,locfv,locfw)
          feqshift(4) = equil_pop04(locrho,locfu,locfv,locfw)
          feqshift(5) = equil_pop05(locrho,locfu,locfv,locfw)
          feqshift(6) = equil_pop06(locrho,locfu,locfv,locfw)
          feqshift(7) = equil_pop07(locrho,locfu,locfv,locfw)
          feqshift(8) = equil_pop08(locrho,locfu,locfv,locfw)
          feqshift(9) = equil_pop09(locrho,locfu,locfv,locfw)
          feqshift(10) = equil_pop10(locrho,locfu,locfv,locfw)
          feqshift(11) = equil_pop11(locrho,locfu,locfv,locfw)
          feqshift(12) = equil_pop12(locrho,locfu,locfv,locfw)
          feqshift(13) = equil_pop13(locrho,locfu,locfv,locfw)
          feqshift(14) = equil_pop14(locrho,locfu,locfv,locfw)
          feqshift(15) = equil_pop15(locrho,locfu,locfv,locfw)
          feqshift(16) = equil_pop16(locrho,locfu,locfv,locfw)
          feqshift(17) = equil_pop17(locrho,locfu,locfv,locfw)
          feqshift(18) = equil_pop18(locrho,locfu,locfv,locfw)
          
          
          do l=0,links
            fluidsub(i,j,k,l)=oneminusomega*fluidsub(i,j,k,l)+ &
            omegaminusone*feq(l)+feqshift(l)
          enddo
           
        enddo
      enddo
    enddo
    
    return
  endif

  ! Variable omega

  ! forall(i=minx:maxx,j=miny:maxy,k=minz:maxz,isfluid(i,j,k)==1)
  do k=minz,maxz
    do j=miny,maxy
      do i=minx,maxx

        if (isfluid(i,j,k)/=1) cycle
        
        temp_omega=omegas(i,j,k)
        oneminusomega = ONE-temp_omega
        omegaminusone = temp_omega-ONE
        
        locrho = rhosub(i,j,k)
        locu = usub(i,j,k)
        locv = vsub(i,j,k)
        locw = wsub(i,j,k)
        

        locfu = fusub(i,j,k)*t_LB / locrho + locu
        locfv = fvsub(i,j,k)*t_LB / locrho + locv
        locfw = fwsub(i,j,k)*t_LB / locrho + locw

        
        
#ifdef REGULARIZED
        pxx = ZERO
        pyy = ZERO
        pzz = ZERO
 	    pxy = ZERO
        pxz = ZERO
        pyz = ZERO
         
        !compute the equilibrium
        feq(0) = equil_pop00(locrho,locu,locv,locw)
        feq(1) = equil_pop01(locrho,locu,locv,locw)
        feq(2) = equil_pop02(locrho,locu,locv,locw)
        feq(3) = equil_pop03(locrho,locu,locv,locw)
        feq(4) = equil_pop04(locrho,locu,locv,locw)
        feq(5) = equil_pop05(locrho,locu,locv,locw)
        feq(6) = equil_pop06(locrho,locu,locv,locw)
        feq(7) = equil_pop07(locrho,locu,locv,locw)
        feq(8) = equil_pop08(locrho,locu,locv,locw)
        feq(9) = equil_pop09(locrho,locu,locv,locw)
        feq(10) = equil_pop10(locrho,locu,locv,locw)
        feq(11) = equil_pop11(locrho,locu,locv,locw)
        feq(12) = equil_pop12(locrho,locu,locv,locw)
        feq(13) = equil_pop13(locrho,locu,locv,locw)
        feq(14) = equil_pop14(locrho,locu,locv,locw)
        feq(15) = equil_pop15(locrho,locu,locv,locw)
        feq(16) = equil_pop16(locrho,locu,locv,locw)
        feq(17) = equil_pop17(locrho,locu,locv,locw)
        feq(18) = equil_pop18(locrho,locu,locv,locw)
        
        do l=0,links
          !compute the equilibrium for the two components
          fneq = fluidsub(i,j,k,l) - feq(l)
          
	      !non equilibrium part of the momentum flux tensor
          pxx=pxx + (dex(l)*dex(l) - cssq)*fneq
          pyy=pyy + (dey(l)*dey(l) - cssq)*fneq
          pzz=pzz + (dez(l)*dez(l) - cssq)*fneq
          pxy=pxy +  dex(l)*dey(l)*fneq
          pxz=pxz +  dex(l)*dez(l)*fneq
          pyz=pyz +  dey(l)*dez(l)*fneq
        enddo
        
        do l=0,links
          fluidsub(i,j,k,l)= feq(l) + &
           ((HALF*p(l))/(cssq**TWO))*((dex(l)*dex(l) - cssq)*pxx + &
           (dey(l)*dey(l) - cssq)*pyy + &
           (dez(l)*dez(l) - cssq)*pzz + &
           TWO*dex(l)*dey(l)*pxy + &
           TWO*dex(l)*dez(l)*pxz + &
           TWO*dey(l)*dez(l)*pyz)           
        enddo
           
#else
          
        !compute the equilibrium
        feq(0) = equil_pop00(locrho,locu,locv,locw)
        feq(1) = equil_pop01(locrho,locu,locv,locw)
        feq(2) = equil_pop02(locrho,locu,locv,locw)
        feq(3) = equil_pop03(locrho,locu,locv,locw)
        feq(4) = equil_pop04(locrho,locu,locv,locw)
        feq(5) = equil_pop05(locrho,locu,locv,locw)
        feq(6) = equil_pop06(locrho,locu,locv,locw)
        feq(7) = equil_pop07(locrho,locu,locv,locw)
        feq(8) = equil_pop08(locrho,locu,locv,locw)
        feq(9) = equil_pop09(locrho,locu,locv,locw)
        feq(10) = equil_pop10(locrho,locu,locv,locw)
        feq(11) = equil_pop11(locrho,locu,locv,locw)
        feq(12) = equil_pop12(locrho,locu,locv,locw)
        feq(13) = equil_pop13(locrho,locu,locv,locw)
        feq(14) = equil_pop14(locrho,locu,locv,locw)
        feq(15) = equil_pop15(locrho,locu,locv,locw)
        feq(16) = equil_pop16(locrho,locu,locv,locw)
        feq(17) = equil_pop17(locrho,locu,locv,locw)
        feq(18) = equil_pop18(locrho,locu,locv,locw)
          
#endif          
          
          
        !compute the shifted equilibrium
        feqshift(0) = equil_pop00(locrho,locfu,locfv,locfw)
        feqshift(1) = equil_pop01(locrho,locfu,locfv,locfw)
        feqshift(2) = equil_pop02(locrho,locfu,locfv,locfw)
        feqshift(3) = equil_pop03(locrho,locfu,locfv,locfw)
        feqshift(4) = equil_pop04(locrho,locfu,locfv,locfw)
        feqshift(5) = equil_pop05(locrho,locfu,locfv,locfw)
        feqshift(6) = equil_pop06(locrho,locfu,locfv,locfw)
        feqshift(7) = equil_pop07(locrho,locfu,locfv,locfw)
        feqshift(8) = equil_pop08(locrho,locfu,locfv,locfw)
        feqshift(9) = equil_pop09(locrho,locfu,locfv,locfw)
        feqshift(10) = equil_pop10(locrho,locfu,locfv,locfw)
        feqshift(11) = equil_pop11(locrho,locfu,locfv,locfw)
        feqshift(12) = equil_pop12(locrho,locfu,locfv,locfw)
        feqshift(13) = equil_pop13(locrho,locfu,locfv,locfw)
        feqshift(14) = equil_pop14(locrho,locfu,locfv,locfw)
        feqshift(15) = equil_pop15(locrho,locfu,locfv,locfw)
        feqshift(16) = equil_pop16(locrho,locfu,locfv,locfw)
        feqshift(17) = equil_pop17(locrho,locfu,locfv,locfw)
        feqshift(18) = equil_pop18(locrho,locfu,locfv,locfw)
          
          
        do l=0,links
          fluidsub(i,j,k,l)=oneminusomega*fluidsub(i,j,k,l)+ &
          omegaminusone*feq(l)+feqshift(l)
        enddo

      enddo
    enddo
  enddo
  
  return
  
 end subroutine collision_fluids_EDM
 
 subroutine update_velocity_EDM(rhoRsub,rhoBsub,fuRsub,fvRsub, &
  fwRsub,fuBsub,fvBsub,fwBsub,usub,vsub,wsub,nstep)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the update to the velocity at
!     half time step in case of EDM forcing scheme
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification January 2020
!     
!***********************************************************************
  
  implicit none
  integer, intent(in) :: nstep
  
  real(kind=PRC), allocatable, dimension(:,:,:)  :: rhoRsub, &
   rhoBsub,fuRsub,fvRsub,fwRsub,fuBsub,fvBsub,fwBsub,usub,vsub,wsub
  real(kind=PRC) :: locrho,locu,locv,locw
  
  integer :: i,j,k,l
  
  
  do k=minz,maxz
    do j=miny,maxy
      do i=minx,maxx

        if (isfluid(i,j,k)/=1) cycle

        locrho = rhoRsub(i,j,k) + rhoBsub(i,j,k)
        locu = usub(i,j,k)
        locv = vsub(i,j,k)
        locw = wsub(i,j,k)
        
        usub(i,j,k) = HALF*(fuRsub(i,j,k)+fuBsub(i,j,k))*t_LB/locrho+locu
        vsub(i,j,k) = HALF*(fvRsub(i,j,k)+fvBsub(i,j,k))*t_LB/locrho+locv
        wsub(i,j,k) = HALF*(fwRsub(i,j,k)+fwBsub(i,j,k))*t_LB/locrho+locw
        
      enddo
    enddo
  enddo
  
  return
  
 end subroutine update_velocity_EDM
 
 subroutine compute_omega
 
!***********************************************************************
!     
!     LBsoft subroutine for computing the relative omega of the
!     fluid system
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  integer :: i,j,k
  real(kind=PRC) :: phis,locrhoR,locrhoB
  
  
  if(lunique_omega)return
  
  select case(iselomegaint)
  case(1)
    do k=minz,maxz
      do j=miny,maxy
        do i=minx,maxx
          phis=(rhoR(i,j,k)/meanR-rhoB(i,j,k)/meanB)/ &
           (rhoR(i,j,k)/meanR+rhoB(i,j,k)/meanB)
          if(abs(phis)>Delta_Grunau)then
            if(phis>ZERO)then
              omega(i,j,k)=ONE/tauR
            else
              omega(i,j,k)=ONE/TauB
            endif
          else
            if(phis>ZERO)then
              omega(i,j,k)=ONE/(s1_R+s2_R*phis+s3_r*phis*phis)
            else
              omega(i,j,k)=ONE/(t1_B+t2_B*phis+t3_B*phis*phis)
            endif
          endif
        enddo
      enddo
    enddo
  case default
    do k=minz,maxz
      do j=miny,maxy
        do i=minx,maxx
          omega(i,j,k)= &
           viscosity_to_omega(qaverage1(rhoR(i,j,k),rhoB(i,j,k),viscR,viscB))
        enddo
      enddo
    enddo
  end select
  
  return
  
 end subroutine compute_omega
 
 subroutine driver_collision_fluids_CG(nstep)

!***********************************************************************
!     
!     LBsoft subroutine for driving the collisional step of CG model
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification November 2019
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: nstep
  
  integer :: i,j,k
    
  
  if(lunique_omega)then
    call collision_fluids_CG_unique(nstep)
  else    
    call collision_fluids_CG(nstep)
  endif
  
  return
  
 end subroutine driver_collision_fluids_CG
    
 pure function qaverage(q,rhoa,rhob,mya,myb)
  
  implicit none
  
  integer, intent(in) :: q 
  
  real*8, intent(in) :: rhoa,rhob,mya,myb
  
  real*8 :: qaverage,rhosum
  
  rhosum=rhoa+rhob
  
  select case(q)
  case(1)
    qaverage=rhoa/rhosum*mya+rhob/rhosum*myb
  case(-1)
    qaverage=1.d0/(rhoa/rhosum*(1.d0/mya)+rhob/rhosum*(1.d0/myb))
  case(0)
    qaverage=(mya**(rhoa)+myb**(rhob))**(1.d0/rhosum)
  case default
    qaverage=0.d0
  end select
  
  return
  
 end function
 
 pure function qaverage1(rhoa,rhob,mya,myb)
  
  implicit none
  
  real*8, intent(in) :: rhoa,rhob,mya,myb
  
  real*8 :: qaverage1,rhosum
  
  rhosum=rhoa+rhob
  
  qaverage1=rhoa/rhosum*mya+rhob/rhosum*myb
  
  return
  
 end function
 
 pure function qaverage0(rhoa,rhob,mya,myb)
  
  implicit none
  
  real*8, intent(in) :: rhoa,rhob,mya,myb
  
  real*8 :: qaverage0,rhosum
  
  rhosum=rhoa+rhob
  
  qaverage0=(mya**(rhoa)+myb**(rhob))**(1.d0/rhosum)
  
  return
  
 end function
 
 pure function qaveragem1(rhoa,rhob,mya,myb)
  
  implicit none
  
  real*8, intent(in) :: rhoa,rhob,mya,myb
  
  real*8 :: qaveragem1,rhosum
  
  rhosum=rhoa+rhob
  
  qaveragem1=1.d0/(rhoa/rhosum*(1.d0/mya)+rhob/rhosum*(1.d0/myb))
  
  return
  
 end function
 
 subroutine collision_fluids_CG_unique(nstep)
 
!***********************************************************************
!     
!     LBsoft subroutine for adding the surface tension and the
!     recolouring step in the collisional step
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification November 2019
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: nstep
  
  real(kind=PRC), parameter :: gradlim = real(1.d-16,kind=PRC)
  real(kind=PRC), parameter :: maxvel = real(0.2d0,kind=PRC)
  
  real(kind=PRC) :: psinorm,psinorm_sq,e_dot_psi,temp,rhosum
  real(kind=PRC) :: psix,psiy,psiz,cosphi,fsum,acoeff
  real(kind=PRC), dimension(0:links) :: feq,feqR,feqB,rhosumtemp,rhodiff
  real(kind=PRC) :: oneminusomega
  real(kind=PRC) :: locrhoR,locrhoB,locu,locv,locw
  real(kind=PRC) :: locrhoR_shifted,locrhoB_shifted,locrhosum_shifted
  integer :: i,j,k,l
  
  real(kind=PRC) :: temp_omega,grad_rhoRx,grad_rhoRy,grad_rhoRz, &
   grad_rhoBx,grad_rhoBy,grad_rhoBz
   
#if PRC==8
  real(kind=PRC), parameter :: mylimit=1.d-100
#elif PRC==4
  real(kind=PRC), parameter :: mylimit=1.e-50
#endif
  
#ifdef REGULARIZED
  real(kind=PRC) :: pxxR,pyyR,pzzR,pxyR,pxzR,pyzR
  real(kind=PRC) :: pxxB,pyyB,pzzB,pxyB,pxzB,pyzB
  real(kind=PRC), dimension(0:links) :: feqR,feqB
  real(kind=PRC) :: fneqR,fneqB
#endif
  
  temp_omega=unique_omega
  oneminusomega = ONE-temp_omega
  
  do k=minz,maxz
    do j=miny,maxy
      do i=minx,maxx
  
        if(isfluid(i,j,k).eq.1)then
          
          locrhoR = rhoR(i,j,k)
          locrhoB = rhoB(i,j,k)
          rhosum = locrhoR+locrhoB
          locu = u(i,j,k)
          locv = v(i,j,k)
          locw = w(i,j,k)
          
#ifdef STABILIZED
          if(abs(locu)>MAXVEL)locu=sign(ONE,locu)*MAXVEL
          if(abs(locv)>MAXVEL)locv=sign(ONE,locv)*MAXVEL
          if(abs(locw)>MAXVEL)locw=sign(ONE,locw)*MAXVEL
#endif
         
         grad_rhoRx = ZERO
         grad_rhoRy = ZERO
         grad_rhoRz = ZERO
         
         grad_rhoBx = ZERO
         grad_rhoBy = ZERO
         grad_rhoBz = ZERO

#ifdef GRADIENTD3Q27
         do l=1,linksd3q27
           locrhoR_shifted=rhoR(i+exd3q27(l),j+eyd3q27(l),k+ezd3q27(l))
           locrhoB_shifted=rhoB(i+exd3q27(l),j+eyd3q27(l),k+ezd3q27(l))
           locrhosum_shifted=locrhoR_shifted+locrhoB_shifted
           grad_rhoRx=grad_rhoRx + ad3q27(l)*dexd3q27(l)* &
            (locrhoR_shifted/locrhosum_shifted)
           grad_rhoRy=grad_rhoRy + ad3q27(l)*deyd3q27(l)* &
            (locrhoR_shifted/locrhosum_shifted)
           grad_rhoRz=grad_rhoRz + ad3q27(l)*dezd3q27(l)* &
            (locrhoR_shifted/locrhosum_shifted)
            
           grad_rhoBx=grad_rhoBx + ad3q27(l)*dexd3q27(l)* &
            (locrhoB_shifted/locrhosum_shifted)
           grad_rhoBy=grad_rhoBy + ad3q27(l)*deyd3q27(l)* &
            (locrhoB_shifted/locrhosum_shifted)
           grad_rhoBz=grad_rhoBz + ad3q27(l)*dezd3q27(l)* &
            (locrhoB_shifted/locrhosum_shifted)
         enddo
#else             
         do l=1,links
           locrhoR_shifted=rhoR(i+ex(l),j+ey(l),k+ez(l))
           locrhoB_shifted=rhoB(i+ex(l),j+ey(l),k+ez(l))
           locrhosum_shifted=locrhoR_shifted+locrhoB_shifted
           grad_rhoRx=grad_rhoRx + a(l)*dex(l)* &
            (locrhoR_shifted/locrhosum_shifted)
           grad_rhoRy=grad_rhoRy + a(l)*dey(l)* &
            (locrhoR_shifted/locrhosum_shifted)
           grad_rhoRz=grad_rhoRz + a(l)*dez(l)* &
            (locrhoR_shifted/locrhosum_shifted)
           
           grad_rhoBx=grad_rhoBx + a(l)*dex(l)* &
            (locrhoB_shifted/locrhosum_shifted)
           grad_rhoBy=grad_rhoBy + a(l)*dey(l)* &
            (locrhoB_shifted/locrhosum_shifted)
           grad_rhoBz=grad_rhoBz + a(l)*dez(l)* &
            (locrhoB_shifted/locrhosum_shifted)
         enddo
#endif    

#ifdef CG_CORRECTION
         do l=0,links
           feqR(l)= &
            equil_popCG_corr(l,temp_omega,alphaR_CG,locrhoR, &
            locu,locv,locw,grad_rhoRx,grad_rhoRy,grad_rhoRz)
           feqB(l)= &
            equil_popCG_corr(l,temp_omega,alphaB_CG,locrhoB, &
            locu,locv,locw,grad_rhoBx,grad_rhoBy,grad_rhoBz)
          enddo
#else
          do l=0,links
           feqR(l)= &
            equil_popCG(l,temp_omega,alphaR_CG,locrhoR,locu,locv,locw)
           feqB(l)= &
            equil_popCG(l,temp_omega,alphaB_CG,locrhoB,locu,locv,locw)
          enddo
#endif
                    
#ifdef REGULARIZED
          pxxR = ZERO
		  pyyR = ZERO
          pzzR = ZERO
		  pxyR = ZERO
          pxzR = ZERO
          pyzR = ZERO
		  pxxB = ZERO
		  pyyB = ZERO
          pzzB = ZERO
		  pxyB = ZERO
          pxzB = ZERO
          pyzB = ZERO
          do l=0,links
            !compute the equilibrium for the two components
             
            fneqR = fluidR(i,j,k,l) - feqR(l)
            fneqB = fluidB(i,j,k,l) - feqB(l)
			
			!non equilibrium part of the momentum flux tensor
			pxxR=pxxR + (dex(l)*dex(l) - cssq)*fneqR
			pyyR=pyyR + (dey(l)*dey(l) - cssq)*fneqR
			pzzR=pzzR + (dez(l)*dez(l) - cssq)*fneqR
			pxyR=pxyR +  dex(l)*dey(l)*fneqR
			pxzR=pxzR +  dex(l)*dez(l)*fneqR
            pyzR=pyzR +  dey(l)*dez(l)*fneqR
			
			pxxB=pxxB + (dex(l)*dex(l) - cssq)*fneqB
			pyyB=pyyB + (dey(l)*dey(l) - cssq)*fneqB
			pzzB=pzzB + (dez(l)*dez(l) - cssq)*fneqB
			pxyB=pxyB +  dex(l)*dey(l)*fneqB
			pxzB=pxzB +  dex(l)*dez(l)*fneqB
            pyzB=pyzB +  dey(l)*dez(l)*fneqB
          enddo
          
          do l=0,links
            fluidR(i,j,k,l)= feqR(l) + &
             ((HALF*p(l))/(cssq**TWO))*((dex(l)*dex(l) - cssq)*pxxR + &
             (dey(l)*dey(l) - cssq)*pyyR + &
             (dez(l)*dez(l) - cssq)*pzzR + &
             TWO*dex(l)*dey(l)*pxyR + &
             TWO*dex(l)*dez(l)*pxzR + &
             TWO*dey(l)*dez(l)*pyzR)
            
            fluidB(i,j,k,l)= feqB(l) + &
             ((HALF*p(l))/(cssq**TWO))*((dex(l)*dex(l) - cssq)*pxxB + &
             (dey(l)*dey(l) - cssq)*pyyB + &
             (dez(l)*dez(l) - cssq)*pzzB + &
             TWO*dex(l)*dey(l)*pxyB + &
             TWO*dex(l)*dez(l)*pxzB + &
             TWO*dey(l)*dez(l)*pyzB)
          enddo
          
#endif
          !bgk step
          do l=0,links
            fluidR(i,j,k,l)=fluidR(i,j,k,l)*oneminusomega + &
             feqR(l)* temp_omega
            fluidB(i,j,k,l)=fluidB(i,j,k,l)*oneminusomega + &
             feqB(l)* temp_omega 
          enddo
          
          !perturbation step
          psix=(locrhoB/rhosum)*grad_rhoRx-(locrhoR/rhosum)*grad_rhoBx
          psiy=(locrhoB/rhosum)*grad_rhoRy-(locrhoR/rhosum)*grad_rhoBy
          psiz=(locrhoB/rhosum)*grad_rhoRz-(locrhoR/rhosum)*grad_rhoBz
          
          psinorm_sq=psix**TWO + psiy**TWO + psiz**TWO
          if(psinorm_sq<=mylimit)cycle
          psinorm=dsqrt(psinorm_sq)
          acoeff=( NINE / FOUR )*temp_omega*sigma_cg
          do l=0,links
            e_dot_psi=dex(l)*psix + dey(l)*psiy + dez(l)*psiz
            temp=psinorm*(p(l)*(e_dot_psi**TWO)/psinorm_sq - b_l(l))
            if(isnan(temp)) temp=ZERO
            fluidR(i,j,k,l)=fluidR(i,j,k,l) + (HALF*acoeff)*temp
            fluidB(i,j,k,l)=fluidB(i,j,k,l) + (HALF*acoeff)*temp
          enddo
          
          !recoloring step
          do l=0,links
            feq(l)= &
             equil_popCG(l,temp_omega,alphaR_CG,locrhoR,ZERO,ZERO,ZERO) + &
             equil_popCG(l,temp_omega,alphaB_CG,locrhoB,ZERO,ZERO,ZERO)
          enddo
          
          do l=0,links
            e_dot_psi=dex(l)*psix + dey(l)*psiy + dez(l)*psiz
            temp=sqrt(dex(l)**TWO + dey(l)**TWO+ dez(l)**TWO)*psinorm
            if(temp<=mylimit)then
              cosphi=ZERO
            else
              cosphi=e_dot_psi/temp
            endif
            temp=beta_CG*locrhoR*locrhoB*cosphi/(rhosum**TWO)
            fsum=fluidR(i,j,k,l) + fluidB(i,j,k,l)
            fluidR(i,j,k,l)=fsum*locrhoR/rhosum + temp*feq(l)
            fluidB(i,j,k,l)=fsum*locrhoB/rhosum - temp*feq(l)
          enddo
        endif
      enddo
    enddo
  enddo
  
  return
  
 end subroutine collision_fluids_CG_unique
 
 subroutine collision_fluids_CG(nstep)
 
!***********************************************************************
!     
!     LBsoft subroutine for adding the surface tension and the
!     recolouring step in the collisional step
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification November 2019
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: nstep
  
  real(kind=PRC), parameter :: gradlim = real(1.d-16,kind=PRC)
  real(kind=PRC), parameter :: maxvel = real(0.2d0,kind=PRC)
  
  real(kind=PRC) :: psinorm,psinorm_sq,e_dot_psi,temp,rhosum
  real(kind=PRC) :: psix,psiy,psiz,cosphi,fsum,acoeff,omega_rel
  real(kind=PRC), dimension(0:links) :: feq
  real(kind=PRC), dimension(0:links) :: feqR,feqB
  real(kind=PRC) :: oneminusomega,uv
  real(kind=PRC) :: locrhoR,locrhoB,locu,locv,locw
  real(kind=PRC) :: locrhoR_shifted,locrhoB_shifted,locrhosum_shifted
  integer :: i,j,k,l
  
  real(kind=PRC) :: temp_omega,grad_rhoRx,grad_rhoRy,grad_rhoRz, &
   grad_rhoBx,grad_rhoBy,grad_rhoBz
  
#if PRC==8
  real(kind=PRC), parameter :: mylimit=1.d-100
#elif PRC==4
  real(kind=PRC), parameter :: mylimit=1.e-50
#endif
  
#ifdef REGULARIZED
  real(kind=PRC) :: pxxR,pyyR,pzzR,pxyR,pxzR,pyzR
  real(kind=PRC) :: pxxB,pyyB,pzzB,pxyB,pxzB,pyzB
  real(kind=PRC) :: fneqR,fneqB
#endif
  
  do k=minz,maxz
    do j=miny,maxy
      do i=minx,maxx
        if(isfluid(i,j,k).eq.1)then
          
          temp_omega = omega(i,j,k)
          oneminusomega = ONE-temp_omega
          
          locrhoR = rhoR(i,j,k)
          locrhoB = rhoB(i,j,k)
          rhosum = locrhoR+locrhoB
          locu = u(i,j,k)
          locv = v(i,j,k)
          locw = w(i,j,k)
          
#ifdef STABILIZED
          if(abs(locu)>MAXVEL)locu=sign(ONE,locu)*MAXVEL
          if(abs(locv)>MAXVEL)locv=sign(ONE,locv)*MAXVEL
          if(abs(locw)>MAXVEL)locw=sign(ONE,locw)*MAXVEL
#endif
         
         grad_rhoRx = ZERO
         grad_rhoRy = ZERO
         grad_rhoRz = ZERO
         
         grad_rhoBx = ZERO
         grad_rhoBy = ZERO
         grad_rhoBz = ZERO

#ifdef GRADIENTD3Q27
         do l=1,linksd3q27
           locrhoR_shifted=rhoR(i+exd3q27(l),j+eyd3q27(l),k+ezd3q27(l))
           locrhoB_shifted=rhoB(i+exd3q27(l),j+eyd3q27(l),k+ezd3q27(l))
           locrhosum_shifted=locrhoR_shifted+locrhoB_shifted
           grad_rhoRx=grad_rhoRx + ad3q27(l)*dexd3q27(l)* &
            (locrhoR_shifted/locrhosum_shifted)
           grad_rhoRy=grad_rhoRy + ad3q27(l)*deyd3q27(l)* &
            (locrhoR_shifted/locrhosum_shifted)
           grad_rhoRz=grad_rhoRz + ad3q27(l)*dezd3q27(l)* &
            (locrhoR_shifted/locrhosum_shifted)
            
           grad_rhoBx=grad_rhoBx + ad3q27(l)*dexd3q27(l)* &
            (locrhoB_shifted/locrhosum_shifted)
           grad_rhoBy=grad_rhoBy + ad3q27(l)*deyd3q27(l)* &
            (locrhoB_shifted/locrhosum_shifted)
           grad_rhoBz=grad_rhoBz + ad3q27(l)*dezd3q27(l)* &
            (locrhoB_shifted/locrhosum_shifted)
         enddo
#else             
         do l=1,links
           locrhoR_shifted=rhoR(i+ex(l),j+ey(l),k+ez(l))
           locrhoB_shifted=rhoB(i+ex(l),j+ey(l),k+ez(l))
           locrhosum_shifted=locrhoR_shifted+locrhoB_shifted
           grad_rhoRx=grad_rhoRx + a(l)*dex(l)* &
            (locrhoR_shifted/locrhosum_shifted)
           grad_rhoRy=grad_rhoRy + a(l)*dey(l)* &
            (locrhoR_shifted/locrhosum_shifted)
           grad_rhoRz=grad_rhoRz + a(l)*dez(l)* &
            (locrhoR_shifted/locrhosum_shifted)
           
           grad_rhoBx=grad_rhoBx + a(l)*dex(l)* &
            (locrhoB_shifted/locrhosum_shifted)
           grad_rhoBy=grad_rhoBy + a(l)*dey(l)* &
            (locrhoB_shifted/locrhosum_shifted)
           grad_rhoBz=grad_rhoBz + a(l)*dez(l)* &
            (locrhoB_shifted/locrhosum_shifted)
         enddo
#endif     

          
#ifdef CG_CORRECTION
         do l=0,links
           feqR(l)= &
            equil_popCG_corr(l,temp_omega,alphaR_CG,locrhoR, &
            locu,locv,locw,grad_rhoRx,grad_rhoRy,grad_rhoRz)
           feqB(l)= &
            equil_popCG_corr(l,temp_omega,alphaB_CG,locrhoB, &
            locu,locv,locw,grad_rhoBx,grad_rhoBy,grad_rhoBz)
          enddo
#else
          do l=0,links
           feqR(l)= &
            equil_popCG(l,temp_omega,alphaR_CG,locrhoR,locu,locv,locw)
           feqB(l)= &
            equil_popCG(l,temp_omega,alphaB_CG,locrhoB,locu,locv,locw)
          enddo
#endif
          
          
#ifdef REGULARIZED
          pxxR = ZERO
		  pyyR = ZERO
          pzzR = ZERO
		  pxyR = ZERO
          pxzR = ZERO
          pyzR = ZERO
		  pxxB = ZERO
		  pyyB = ZERO
          pzzB = ZERO
		  pxyB = ZERO
          pxzB = ZERO
          pyzB = ZERO
          do l=0,links
            !compute the equilibrium for the two components
             
            fneqR = fluidR(i,j,k,l) - feqR(l)
            fneqB = fluidB(i,j,k,l) - feqB(l)
			
			!non equilibrium part of the momentum flux tensor
			pxxR=pxxR + (dex(l)*dex(l) - cssq)*fneqR
			pyyR=pyyR + (dey(l)*dey(l) - cssq)*fneqR
			pzzR=pzzR + (dez(l)*dez(l) - cssq)*fneqR
			pxyR=pxyR +  dex(l)*dey(l)*fneqR
			pxzR=pxzR +  dex(l)*dez(l)*fneqR
            pyzR=pyzR +  dey(l)*dez(l)*fneqR
			
			pxxB=pxxB + (dex(l)*dex(l) - cssq)*fneqB
			pyyB=pyyB + (dey(l)*dey(l) - cssq)*fneqB
			pzzB=pzzB + (dez(l)*dez(l) - cssq)*fneqB
			pxyB=pxyB +  dex(l)*dey(l)*fneqB
			pxzB=pxzB +  dex(l)*dez(l)*fneqB
            pyzB=pyzB +  dey(l)*dez(l)*fneqB
          enddo
          
          do l=0,links
            fluidR(i,j,k,l)= feqR(l) + &
             ((HALF*p(l))/(cssq**TWO))*((dex(l)*dex(l) - cssq)*pxxR + &
             (dey(l)*dey(l) - cssq)*pyyR + &
             (dez(l)*dez(l) - cssq)*pzzR + &
             TWO*dex(l)*dey(l)*pxyR + &
             TWO*dex(l)*dez(l)*pxzR + &
             TWO*dey(l)*dez(l)*pyzR)
            
            fluidB(i,j,k,l)= feqB(l) + &
             ((HALF*p(l))/(cssq**TWO))*((dex(l)*dex(l) - cssq)*pxxB + &
             (dey(l)*dey(l) - cssq)*pyyB + &
             (dez(l)*dez(l) - cssq)*pzzB + &
             TWO*dex(l)*dey(l)*pxyB + &
             TWO*dex(l)*dez(l)*pxzB + &
             TWO*dey(l)*dez(l)*pyzB)
          enddo
          
#endif
          !bgk step
          do l=0,links
            fluidR(i,j,k,l)=fluidR(i,j,k,l)*oneminusomega + &
             feqR(l)* temp_omega
            fluidB(i,j,k,l)=fluidB(i,j,k,l)*oneminusomega + &
             feqB(l)* temp_omega 
          enddo
          
          !perturbation step
          
          psix=(locrhoB/rhosum)*grad_rhoRx-(locrhoR/rhosum)*grad_rhoBx
          psiy=(locrhoB/rhosum)*grad_rhoRy-(locrhoR/rhosum)*grad_rhoBy
          psiz=(locrhoB/rhosum)*grad_rhoRz-(locrhoR/rhosum)*grad_rhoBz
          
          psinorm_sq=psix**TWO + psiy**TWO + psiz**TWO
          if(psinorm_sq<=mylimit)cycle
          psinorm=dsqrt(psinorm_sq)
          acoeff=( NINE / FOUR )*temp_omega*sigma_cg
          do l=0,links
            e_dot_psi=dex(l)*psix + dey(l)*psiy + dez(l)*psiz
            temp=psinorm*(p(l)*(e_dot_psi**TWO)/psinorm_sq - b_l(l))
            if(isnan(temp)) temp=ZERO
            fluidR(i,j,k,l)=fluidR(i,j,k,l) + (HALF*acoeff)*temp
            fluidB(i,j,k,l)=fluidB(i,j,k,l) + (HALF*acoeff)*temp
          enddo
          
          !recoloring step
          do l=0,links
            feq(l)= &
             equil_popCG(l,temp_omega,alphaR_CG,locrhoR,ZERO,ZERO,ZERO) + &
             equil_popCG(l,temp_omega,alphaB_CG,locrhoB,ZERO,ZERO,ZERO)
          enddo
          
          do l=0,links
            e_dot_psi=dex(l)*psix + dey(l)*psiy + dez(l)*psiz
            temp=sqrt(dex(l)**TWO + dey(l)**TWO+ dez(l)**TWO)*psinorm
            if(temp<=mylimit)then
              cosphi=ZERO
            else
              cosphi=e_dot_psi/temp
            endif
            temp=beta_CG*locrhoR*locrhoB*cosphi/(rhosum**TWO)
            fsum=fluidR(i,j,k,l) + fluidB(i,j,k,l)
            fluidR(i,j,k,l)=fsum*locrhoR/rhosum + temp*feq(l)
            fluidB(i,j,k,l)=fsum*locrhoB/rhosum - temp*feq(l)
          enddo
        endif
      enddo
    enddo
  enddo
  
  return
  
 end subroutine collision_fluids_CG
 
 subroutine driver_bc_pop_selfcomm(lparticles)
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the pbc within the same process
!     on the Boltzmann populations
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************

  implicit none
  
  logical, intent(in) :: lparticles
  
  call manage_bc_pop_selfcomm(fluidR,lparticles)
  
  if(lsingle_fluid)return
  
  call manage_bc_pop_selfcomm(fluidB,lparticles)
  
 
 end subroutine driver_bc_pop_selfcomm
 


  subroutine stream_nocopy(fluidsub)
    
!***********************************************************************
!     
!     LBsoft subroutine for streaming the populations without
!     a buffer copy
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: F. Bonaccorso
!     last modification April 2019
!     
!***********************************************************************

  implicit none
  real(kind=PRC), allocatable, dimension(:,:,:,:)  :: fluidsub
  integer :: i,j,k
  logical, save :: isFirst = .true.


  if (isFirst) then
        isFirst = .false.
        if (idrank == 0) write(6,*) "STREAM: Using stream_nocopy"
  endif


      do i=maxx+1, minx, -1
         fluidsub(i,:,:,1) =  fluidsub(i-1,:,:,1)
      enddo
      do i=minx-1, maxx
         fluidsub(i,:,:,2) =  fluidsub(i+1,:,:,2)
      enddo

      do j=maxy+1, miny, -1
         fluidsub(:,j,:,3) =  fluidsub(:,j-1,:,3)
      enddo
      do j=miny-1, maxy
         fluidsub(:,j,:,4) =  fluidsub(:,j+1,:,4)
      enddo

      do k=maxz+1, minz, -1
         fluidsub(:,:,k,5) =  fluidsub(:,:,k-1,5)
      enddo
      do k=minz-1, maxz
         fluidsub(:,:,k,6) =  fluidsub(:,:,k+1,6)
      enddo

      do i=maxx+1, minx, -1
       do j=maxy+1, miny, -1
         fluidsub(i,j,:,7) =  fluidsub(i-1,j-1,:,7)
       enddo
      enddo
      do i=minx-1, maxx
       do j=miny-1, maxy
         fluidsub(i,j,:,8) =  fluidsub(i+1,j+1,:,8)
       enddo
      enddo

      do i=minx-1, maxx
       do j=maxy+1, miny, -1
         fluidsub(i,j,:,9) =  fluidsub(i+1,j-1,:,9)
       enddo
      enddo
      do i=maxx+1, minx, -1
       do j=miny-1, maxy
         fluidsub(i,j,:,10) =  fluidsub(i-1,j+1,:,10)
       enddo
      enddo

      do i=maxx+1, minx, -1
       do k=maxz+1, minz, -1
         fluidsub(i,:,k,11) =  fluidsub(i-1,:,k-1,11)
       enddo
      enddo
      do i=minx-1, maxx
       do k=minz-1, maxz
         fluidsub(i,:,k,12) =  fluidsub(i+1,:,k+1,12)
       enddo
      enddo

      do i=minx-1, maxx
       do k=maxz+1, minz, -1
         fluidsub(i,:,k,13) =  fluidsub(i+1,:,k-1,13)
       enddo
      enddo
      do i=maxx+1, minx, -1
       do k=minz-1, maxz
         fluidsub(i,:,k,14) =  fluidsub(i-1,:,k+1,14)
       enddo
      enddo

      do j=maxy+1, miny, -1
       do k=maxz+1, minz, -1
         fluidsub(:,j,k,15) =  fluidsub(:,j-1,k-1,15)
       enddo
      enddo
      do j=miny-1, maxy
       do k=minz-1, maxz
         fluidsub(:,j,k,16) =  fluidsub(:,j+1,k+1,16)
       enddo
      enddo

      do j=miny-1, maxy
       do k=maxz+1, minz, -1
         fluidsub(:,j,k,17) =  fluidsub(:,j+1,k-1,17)
       enddo
      enddo
      do j=maxy+1, miny, -1
       do k=minz-1, maxz
         fluidsub(:,j,k,18) =  fluidsub(:,j-1,k+1,18)
       enddo
      enddo

     end subroutine stream_nocopy
     
  subroutine stream_nocopy_isfluid(fluidsub)
    
!***********************************************************************
!     
!     LBsoft subroutine for streaming the populations without
!     a buffer copy
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: F. Bonaccorso
!     last modification April 2019
!     
!***********************************************************************

  implicit none
  
  real(kind=PRC), allocatable, dimension(:,:,:,:)  :: fluidsub
  
  integer :: i,j,k
  logical, save :: isFirst = .true.


  if (isFirst) then
        isFirst = .false.
        if (idrank == 0) write(6,*) "STREAM: Using stream_nocopy_isfluid"
  endif
     
  do k=minz-1, maxz+1
    do j=miny-1, maxy+1
      do i=maxx+1, minx, -1
        if ( isfluid(i-1,j,k)<3 .or. isfluid(i-1,j,k)>4) then
          fluidsub(i,j,k,1) =  fluidsub(i-1,j,k,1)
        endif
      enddo
    enddo
  enddo
  
  do k=minz-1, maxz+1
    do j=miny-1, maxy+1
      do i=minx-1, maxx
        if ( isfluid(i+1,j,k)<3 .or. isfluid(i+1,j,k)>4) then
          fluidsub(i,j,k,2) =  fluidsub(i+1,j,k,2)
        endif
      enddo
    enddo
  enddo
  
  do k=minz-1, maxz+1
    do j=maxy+1, miny, -1
      do i=minx-1, maxx+1
        if ( isfluid(i,j-1,k)<3 .or. isfluid(i,j-1,k)>4) then
          fluidsub(i,j,k,3) =  fluidsub(i,j-1,k,3)
        endif
      enddo
    enddo
  enddo
  
  do k=minz-1, maxz+1
    do j=miny-1, maxy
      do i=minx-1, maxx+1
        if ( isfluid(i,j+1,k)<3 .or. isfluid(i,j+1,k)>4) then
          fluidsub(i,j,k,4) =  fluidsub(i,j+1,k,4)
        endif
      enddo
    enddo
  enddo
  
  do k=maxz+1, minz, -1
    do j=miny-1, maxy+1
      do i=minx-1, maxx+1
        if ( isfluid(i,j,k-1)<3 .or. isfluid(i,j,k-1)>4) then
          fluidsub(i,j,k,5) =  fluidsub(i,j,k-1,5)
        endif
      enddo
    enddo
  enddo
  
  do k=minz-1, maxz
    do j=miny-1, maxy+1
      do i=minx-1, maxx+1
        if ( isfluid(i,j,k+1)<3 .or. isfluid(i,j,k+1)>4) then
          fluidsub(i,j,k,6) =  fluidsub(i,j,k+1,6)
        endif
      enddo
    enddo
  enddo
  
  do k=minz-1, maxz+1
    do i=maxx+1, minx, -1
      do j=maxy+1, miny, -1
        if ( isfluid(i-1,j-1,k)<3 .or. isfluid(i-1,j-1,k)>4) then
          fluidsub(i,j,k,7) =  fluidsub(i-1,j-1,k,7)
        endif
      enddo
    enddo
  enddo
  
  do k=minz-1, maxz+1
    do i=minx-1, maxx
      do j=miny-1, maxy
        if ( isfluid(i+1,j+1,k)<3 .or. isfluid(i+1,j+1,k)>4) then
          fluidsub(i,j,k,8) =  fluidsub(i+1,j+1,k,8)
        endif
      enddo
    enddo
  enddo
  
  do k=minz-1, maxz+1
    do i=minx-1, maxx
      do j=maxy+1, miny, -1
        if ( isfluid(i+1,j-1,k)<3 .or. isfluid(i+1,j-1,k)>4) then
          fluidsub(i,j,k,9) =  fluidsub(i+1,j-1,k,9)
        endif
      enddo
    enddo
  enddo
   
  do k=minz-1, maxz+1
    do i=maxx+1, minx, -1
      do j=miny-1, maxy
        if ( isfluid(i-1,j+1,k)<3 .or. isfluid(i-1,j+1,k)>4) then
          fluidsub(i,j,k,10) =  fluidsub(i-1,j+1,k,10)
        endif
      enddo
    enddo
  enddo
   
  do k=maxz+1, minz, -1
    do j=miny-1, maxy+1
      do i=maxx+1, minx, -1
        if ( isfluid(i-1,j,k-1)<3 .or. isfluid(i-1,j,k-1)>4) then
          fluidsub(i,j,k,11) =  fluidsub(i-1,j,k-1,11)
        endif
      enddo
    enddo
  enddo
    
  do k=minz-1, maxz
    do j=miny-1, maxy+1
      do i=minx-1, maxx
        if ( isfluid(i+1,j,k+1)<3 .or. isfluid(i+1,j,k+1)>4) then
          fluidsub(i,j,k,12) =  fluidsub(i+1,j,k+1,12)
        endif
      enddo
    enddo
  enddo
    
  do k=maxz+1, minz, -1
    do j=miny-1, maxy+1
      do i=minx-1, maxx
        if ( isfluid(i+1,j,k-1)<3 .or. isfluid(i+1,j,k-1)>4) then
          fluidsub(i,j,k,13) =  fluidsub(i+1,j,k-1,13)
        endif
       enddo
     enddo
  enddo
  
  do k=minz-1, maxz
    do j=miny-1, maxy+1
      do i=maxx+1, minx, -1
        if ( isfluid(i-1,j,k+1)<3 .or. isfluid(i-1,j,k+1)>4) then
          fluidsub(i,j,k,14) =  fluidsub(i-1,j,k+1,14)
        endif
      enddo
    enddo
  enddo
  
  do k=maxz+1, minz, -1
    do j=maxy+1, miny, -1
      do i=minx-1, maxx+1
        if ( isfluid(i,j-1,k-1)<3 .or. isfluid(i,j-1,k-1)>4) then
          fluidsub(i,j,k,15) =  fluidsub(i,j-1,k-1,15)
        endif
      enddo
    enddo
  enddo
      
  do k=minz-1, maxz
    do j=miny-1, maxy
      do i=minx-1, maxx+1
        if ( isfluid(i,j+1,k+1)<3 .or. isfluid(i,j+1,k+1)>4) then
          fluidsub(i,j,k,16) =  fluidsub(i,j+1,k+1,16)
        endif
      enddo
    enddo
  enddo
   
  do k=maxz+1, minz, -1
    do j=miny-1, maxy
      do i=minx-1, maxx+1
        if ( isfluid(i,j+1,k-1)<3 .or. isfluid(i,j+1,k-1)>4) then
          fluidsub(i,j,k,17) =  fluidsub(i,j+1,k-1,17)
        endif
      enddo
    enddo
  enddo
      
  do k=minz-1, maxz
    do j=maxy+1, miny, -1
      do i=minx-1, maxx+1
        if ( isfluid(i,j-1,k+1)<3 .or. isfluid(i,j-1,k+1)>4) then
          fluidsub(i,j,k,18) =  fluidsub(i,j-1,k+1,18)
        endif
      enddo
    enddo
  enddo
  
  return
  
  end subroutine stream_nocopy_isfluid

  subroutine stream_copy(fluidsub)
  
!***********************************************************************
!     
!     LBsoft subroutine for streaming the populations without
!     a buffer copy
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: F. Bonaccorso
!     last modification May 2020
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC), allocatable, dimension(:,:,:,:)  :: fluidsub
  integer :: i,j,k,l,ishift,jshift,kshift
  logical, save :: isFirst = .true.


  if (isFirst) then
        isFirst = .false.
        if (idrank == 0) write(6,*) "STREAM: Using stream_copy"
  endif
  
  do l=1,links
    ishift=ex(l)
    jshift=ey(l)
    kshift=ez(l)
    do k=minz-1,maxz+1
     do j=miny-1,maxy+1
      do i=minx-1,maxx+1
        if ( isfluid(i,j,k)<3 .or. isfluid(i,j,k)>4) then
         buffservice3d(i+ishift,j+jshift,k+kshift) = fluidsub(i,j,k,l)
        endif
      enddo
     enddo
    enddo

    do k=minz-1,maxz+1
     do j=miny-1,maxy+1
      do i=minx-1,maxx+1
        if ( isfluid(i,j,k)<3 .or. isfluid(i,j,k)>4) then
         fluidsub(i,j,k,l) = buffservice3d(i,j,k)
        endif
      enddo
     enddo
    enddo
  enddo
  
  return
  
 end subroutine stream_copy

 subroutine driver_streaming_fluids(lparticles)
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the streaming step
!     on the Boltzmann populations
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  logical, intent(in) :: lparticles
  integer :: i,j,k,l,ishift,jshift,kshift
  integer, save :: iter=0
  
  iter=iter+1


  



#ifdef DIAGNSTREAM
   if(iter.eq.NDIAGNSTREAM) then
     if(allocated(ownern))then
       call print_all_pops(100,'mioprima',iter,fluidR)
     endif
   endif
#endif
    

  call commspop(fluidR)



#ifdef STREAM_NOCOPY

  if(lparticles)then
  
  call stream_nocopy_isfluid(fluidR)
  
  else
  
    if(lread_isfluid)then
      call stream_nocopy_isfluid(fluidR)
    else
      call stream_nocopy(fluidR)
    endif
  endif
#else
  call stream_copy(fluidR)
#endif /* STREAM_NOCOPY */





  call manage_bc_pop_selfcomm(fluidR,lparticles)

  call commrpop(fluidR,lparticles,isfluid)



#ifdef DIAGNSTREAM
  if(iter.eq.NDIAGNSTREAM) then
    if(allocated(ownern))then
       call print_all_pops(300,'miodopo',iter,fluidR)
    endif
  endif
#endif


  
  if(lsingle_fluid)return


  call commspop(fluidB)




#ifdef STREAM_NOCOPY
  if(lparticles)then
  
  call stream_nocopy_isfluid(fluidB)
  
  else
  if(lread_isfluid)then
      call stream_nocopy_isfluid(fluidB)
    else
      call stream_nocopy(fluidB)
    endif
  endif
#else
  call stream_copy(fluidB)
#endif /* STREAM_NOCOPY */
  
  call manage_bc_pop_selfcomm(fluidB,lparticles)

  call commrpop(fluidB,lparticles,isfluid)
  
  return
  
 end subroutine driver_streaming_fluids
 
 subroutine moments_fluids(nstep)
 
!***********************************************************************
!     
!     LBsoft subroutine for computing the moments on the 
!     Boltzmann populations and estimate the hydrodynamic variables
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  integer, intent(in) :: nstep
  integer :: i,j,k,l
  real(kind=PRC) :: factR, factB
  real(kind=PRC) :: locrho,locu,locv,locw, invrho
  real(kind=PRC) :: locrhor, locrhob, weight_RB,rhosum
  real(kind=PRC) :: ddx,ddy,ddz,ddxB,ddyB,ddzB
  real(kind=8), parameter :: RHOEPS = 9.0E-15
  


  !compute density and accumulate mass flux
  
  !red fluid
  
  if(lsingle_fluid)then
    nfluidnodes=0
    nfluidmassR=ZERO
    
    do k=minz,maxz
      do j=miny,maxy
        do i=minx,maxx
          locrho = &
           fluidR(i,j,k,0) + fluidR(i,j,k,1) + fluidR(i,j,k,2) + &
           fluidR(i,j,k,3) + fluidR(i,j,k,4) + &
           fluidR(i,j,k,5) + fluidR(i,j,k,6) + fluidR(i,j,k,7) + &
           fluidR(i,j,k,8) + fluidR(i,j,k,9) + &
           fluidR(i,j,k,10) + fluidR(i,j,k,11) + fluidR(i,j,k,12) + &
           fluidR(i,j,k,13) + fluidR(i,j,k,14) + &
           fluidR(i,j,k,15) + fluidR(i,j,k,16) + &
           fluidR(i,j,k,17) + fluidR(i,j,k,18)

          invrho = ONE / locrho

          locu   = invrho * ( &
           fluidR(i,j,k,1) - fluidR(i,j,k,2) + fluidR(i,j,k,7) - &
           fluidR(i,j,k,8) - fluidR(i,j,k,9) + &
           fluidR(i,j,k,10) + fluidR(i,j,k,11) - fluidR(i,j,k,12) - &
           fluidR(i,j,k,13) + fluidR(i,j,k,14) )

          locv    = invrho * ( &
           fluidR(i,j,k,3) - fluidR(i,j,k,4) + fluidR(i,j,k,7) - &
           fluidR(i,j,k,8) + fluidR(i,j,k,9) - &
           fluidR(i,j,k,10) + fluidR(i,j,k,15) - fluidR(i,j,k,16) - &
           fluidR(i,j,k,17) + fluidR(i,j,k,18) )

          locw    = invrho * ( &
           fluidR(i,j,k,5) - fluidR(i,j,k,6) + fluidR(i,j,k,11) - &
           fluidR(i,j,k,12) + fluidR(i,j,k,13) - &
           fluidR(i,j,k,14) + fluidR(i,j,k,15) - fluidR(i,j,k,16) + &
           fluidR(i,j,k,17) - fluidR(i,j,k,18) )

          if (isfluid(i,j,k)==1) then
            if (locrho < MINDENS) locrho = MINDENS
            rhoR(i,j,k) = locrho
            u(i,j,k) = locu
            v(i,j,k) = locv
            w(i,j,k) = locw
            nfluidnodes=nfluidnodes+1
            nfluidmassR=nfluidmassR+locrho
          else
            rhoR(i,j,k) = MINDENS
            u(i,j,k) = ZERO
            v(i,j,k) = ZERO
            w(i,j,k) = ZERO
          endif

        enddo
      enddo
    enddo
    

    return
    
  endif
  

  nfluidnodes=0
  nfluidmassR=ZERO
  nfluidmassB=ZERO
  do k=minz,maxz
    do j=miny,maxy
      do i=minx,maxx
        locrhor = &
         fluidR(i,j,k,0) + fluidR(i,j,k,1) + fluidR(i,j,k,2) + &
         fluidR(i,j,k,3) + fluidR(i,j,k,4) + &
         fluidR(i,j,k,5) + fluidR(i,j,k,6) + fluidR(i,j,k,7) + &
         fluidR(i,j,k,8) + fluidR(i,j,k,9) + &
         fluidR(i,j,k,10) + fluidR(i,j,k,11) + fluidR(i,j,k,12) + &
         fluidR(i,j,k,13) + fluidR(i,j,k,14) + &
         fluidR(i,j,k,15) + fluidR(i,j,k,16) + &
         fluidR(i,j,k,17) + fluidR(i,j,k,18)

        locrhob = &
         fluidB(i,j,k,0) + fluidB(i,j,k,1) + fluidB(i,j,k,2) + &
         fluidB(i,j,k,3) + fluidB(i,j,k,4) + &
         fluidB(i,j,k,5) + fluidB(i,j,k,6) + fluidB(i,j,k,7) + &
         fluidB(i,j,k,8) + fluidB(i,j,k,9) + &
         fluidB(i,j,k,10) + fluidB(i,j,k,11) + fluidB(i,j,k,12) + &
         fluidB(i,j,k,13) + fluidB(i,j,k,14) + &
         fluidB(i,j,k,15) + fluidB(i,j,k,16) + &
         fluidB(i,j,k,17) + fluidB(i,j,k,18)
        
        locu   = ( &
         fluidR(i,j,k,1) - fluidR(i,j,k,2) + fluidR(i,j,k,7) - &
         fluidR(i,j,k,8) - fluidR(i,j,k,9) + &
         fluidR(i,j,k,10) + fluidR(i,j,k,11) - fluidR(i,j,k,12) - &
         fluidR(i,j,k,13) + fluidR(i,j,k,14) ) + ( &
         fluidB(i,j,k,1) - fluidB(i,j,k,2) + fluidB(i,j,k,7) - &
         fluidB(i,j,k,8) - fluidB(i,j,k,9) + &
         fluidB(i,j,k,10) + fluidB(i,j,k,11) - fluidB(i,j,k,12) - &
         fluidB(i,j,k,13) + fluidB(i,j,k,14) )

        locv    = ( &
         fluidR(i,j,k,3) - fluidR(i,j,k,4) + fluidR(i,j,k,7) - &
         fluidR(i,j,k,8) + fluidR(i,j,k,9) - &
         fluidR(i,j,k,10) + fluidR(i,j,k,15) - fluidR(i,j,k,16) - &
         fluidR(i,j,k,17) + fluidR(i,j,k,18) ) + ( &
         fluidB(i,j,k,3) - fluidB(i,j,k,4) + fluidB(i,j,k,7) - &
         fluidB(i,j,k,8) + fluidB(i,j,k,9) - &
         fluidB(i,j,k,10) + fluidB(i,j,k,15) - fluidB(i,j,k,16) - &
         fluidB(i,j,k,17) + fluidB(i,j,k,18) ) 

        locw    = ( &
         fluidR(i,j,k,5) - fluidR(i,j,k,6) + fluidR(i,j,k,11) - &
         fluidR(i,j,k,12) + fluidR(i,j,k,13) - &
         fluidR(i,j,k,14) + fluidR(i,j,k,15) - fluidR(i,j,k,16) + &
         fluidR(i,j,k,17) - fluidR(i,j,k,18) ) + ( &
         fluidB(i,j,k,5) - fluidB(i,j,k,6) + fluidB(i,j,k,11) - &
         fluidB(i,j,k,12) + fluidB(i,j,k,13) - &
         fluidB(i,j,k,14) + fluidB(i,j,k,15) - fluidB(i,j,k,16) + &
         fluidB(i,j,k,17) - fluidB(i,j,k,18) )

         if (isfluid(i,j,k)==1) then
           if (locrhor < MINDENS) locrhor = MINDENS
           if (locrhob < MINDENS) locrhob = MINDENS
           rhoR(i,j,k) = locrhor
           rhoB(i,j,k) = locrhob
           weight_RB = ONE / (locrhor + locrhob)
           u(i,j,k) = locu * weight_RB
           v(i,j,k) = locv * weight_RB
           w(i,j,k) = locw * weight_RB
           nfluidnodes=nfluidnodes+1
           nfluidmassR=nfluidmassR+locrhor
           nfluidmassB=nfluidmassB+locrhob
         else
           rhoR(i,j,k) = MINDENS
           rhoB(i,j,k) = MINDENS
           u(i,j,k) = ZERO
           v(i,j,k) = ZERO
           w(i,j,k) = ZERO
         endif
       enddo
     enddo
   enddo
    
   return
    
 end subroutine moments_fluids


 subroutine probe_red_moments_in_node(i,j,k,dtemp1,dtemp2,dtemp3,dtemp4)
 
!***********************************************************************
!     
!     LBsoft subroutine for probing the Red hydrodynamic moments
!     at the point i j k
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: i,j,k
  real(kind=PRC), intent(out) :: dtemp1,dtemp2,dtemp3,dtemp4
  
  integer :: l
  
  dtemp1=ZERO
  dtemp2=ZERO
  dtemp3=ZERO
  dtemp4=ZERO
  
  do l=0,links
    dtemp1 = fluidR(i,j,k,l) + dtemp1 
    dtemp2    = fluidR(i,j,k,l)*dex(l) + dtemp2
    dtemp3    = fluidR(i,j,k,l)*dey(l) + dtemp3
    dtemp4    = fluidR(i,j,k,l)*dez(l) + dtemp4
  enddo
  
  return
  
 end subroutine probe_red_moments_in_node
 
 subroutine probe_blue_moments_in_node(i,j,k,dtemp1,dtemp2,dtemp3,dtemp4)
 
!***********************************************************************
!     
!     LBsoft subroutine for probing the Blue hydrodynamic moments
!     at the point i j k
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: i,j,k
  real(kind=PRC), intent(out) :: dtemp1,dtemp2,dtemp3,dtemp4
  
  integer :: l
  
  dtemp1=ZERO
  dtemp2=ZERO
  dtemp3=ZERO
  dtemp4=ZERO
  
  do l=0,links
    dtemp1 = fluidB(i,j,k,l) + dtemp1 
    dtemp2    = fluidB(i,j,k,l)*dex(l) + dtemp2
    dtemp3    = fluidB(i,j,k,l)*dey(l) + dtemp3
    dtemp4    = fluidB(i,j,k,l)*dez(l) + dtemp4
  enddo
  
  return
  
 end subroutine probe_blue_moments_in_node
 
 subroutine probe_pops_in_node(itemp,jtemp,ktemp,nstepsub,fluidsub,rhosub)
  
!***********************************************************************
!     
!     LBsoft subroutine for probing the population values
!     at the point i j k
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification September 2018
!     
!***********************************************************************
  
  integer, intent(in) :: itemp,jtemp,ktemp,nstepsub
  real(kind=PRC), allocatable :: fluidsub(:,:,:,:)
  real(kind=PRC), allocatable :: rhosub(:,:,:)
  integer(kind=IPRC) :: i4
  
  integer :: l
  
  i4=i4back(itemp,jtemp,ktemp)
  if(ownern(i4)==idrank)then
    write(6,*)'------------------------------------------------','id =',idrank, &
     rhosub(itemp,jtemp,ktemp)
    do l=0,18
    write(6,*)'l =',l,nstepsub,fluidsub(itemp,jtemp,ktemp,l)
    call flush(6)
    enddo
  endif
  return
  
 end subroutine probe_pops_in_node
 

!*****************START PART TO MANAGE THE PERIODIC BC******************

 subroutine initialize_isfluid_bcfluid(lvtkfilesub,lparticles)
 
!***********************************************************************
!     
!     LBsoft subroutine for initialize the isfluid and bcfluid
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  logical, intent(in) :: lvtkfilesub, lparticles
  
  call set_isfluid_bcfluid(lvtkfilesub,lsingle_fluid,lparticles, &
   lunique_omega,lexch_dens,lexch_u,lexch_v,lexch_w,unique_omega,omega, &
   rhoR,rhoB,u,v,w)
  
  return
  
 end subroutine initialize_isfluid_bcfluid
 
 subroutine driver_bc_densities
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the boundary conditions
!     to fluid densities
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  Logical, save   :: isFirst = .true.
  
  if(lexch_dens)then
  
    call commexch_dens(rhoR,rhoB)
  
    call manage_bc_hvar_selfcomm(rhoR)
    if(.not. lsingle_fluid)then
      call manage_bc_hvar_selfcomm(rhoB)
    endif
  
    call commwait_dens(rhoR,rhoB)
  
  endif
  
 end subroutine driver_bc_densities

 subroutine driver_bc_psi
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the boundary conditions
!     to fluid densities
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  implicit none
  
  
  call commexch_dens(psiR,psiB)
  

  call manage_bc_hvar_selfcomm(psiR)
  if(.not. lsingle_fluid)then
    call manage_bc_hvar_selfcomm(psiB)
  endif

  
  call commwait_dens(psiR,psiB)
  
  
 end subroutine driver_bc_psi
 
 subroutine driver_bc_velocities
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the boundary conditions
!     to fluid densities
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  Logical, save   :: isFirst = .true.

  if(lexch_u)then
  

    call commexch_vel_component(u)

    call manage_bc_hvar_selfcomm(u)

    call commwait_vel_component(u)

  
  endif
  
  if(lexch_v)then
  

    call commexch_vel_component(v)

    call manage_bc_hvar_selfcomm(v)

    call commwait_vel_component(v)
  
  endif
  
  if(lexch_w)then
  

    call commexch_vel_component(w)

    call manage_bc_hvar_selfcomm(w)

    call commwait_vel_component(w)

  
  endif
  
  isFirst = .false.

  return

  
 end subroutine driver_bc_velocities
 
 subroutine driver_bc_all_pops()
 
!***********************************************************************
!     
!     LBsoft subroutine to drive the population buffer fluid nodes
!     among MPI processing and within the same process for all the
!     populations
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification May 2020
!     
!*********************************************************************** 
 
  implicit none
  
  integer :: l
  
  do l=0,links
    call commexch_single_halo(fluidR(:,:,:,l), &
     nbuff,minx,maxx,miny,maxy,minz,maxz)
    call manage_bc_singlehalo_selfcomm(fluidR(:,:,:,l), &
     nbuff,minx,maxx,miny,maxy,minz,maxz)
    call commwait_single_halo(fluidR(:,:,:,l), &
     nbuff,minx,maxx,miny,maxy,minz,maxz)
  enddo
  
  if(lsingle_fluid)return
 
  do l=0,links
    call commexch_single_halo(fluidB(:,:,:,l), &
     nbuff,minx,maxx,miny,maxy,minz,maxz)
    call manage_bc_singlehalo_selfcomm(fluidB(:,:,:,l), &
     nbuff,minx,maxx,miny,maxy,minz,maxz)
    call commwait_single_halo(fluidB(:,:,:,l), &
     nbuff,minx,maxx,miny,maxy,minz,maxz)
  enddo
 
  return
  
 end subroutine driver_bc_all_pops

!******************END PART TO MANAGE THE PERIODIC BC*******************
 
!*****************START PART TO MANAGE THE BOUNCEBACK*******************
 
 subroutine driver_apply_bounceback_halfway_pop(nstep)
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the bounce back of the fluid
!     populations if requested from the boundary conditions 
!     in halfway mode.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: nstep
  
  integer ::i,j,k,l
  integer, save :: iter=0
  integer :: idir,inits,ends
  
  iter=iter+1
  
  call set_bc_variable_hvar(lsingle_fluid,lunique_omega,unique_omega, &
   omega,rhoR,rhoB,u,v,w)
  
  
  call apply_bounceback_pop_halfway(nstep,lColourG,bc_rhoR,bc_u,bc_v,bc_w, &
   rhoR,rhoB,fluidR,alphaR_CG)
  
#ifdef DIAGNSTREAM
  if(iter==NDIAGNSTREAM)call print_all_pops(100,'miodopobounce',iter,fluidR)
#endif
#if 0
if(mod(nstep,1)==0)then
     do idir=1,nbcdir
      inits=nbounce7dir(idir-1)+1
      ends=nbounce7dir(idir)
      if(inits>ends)cycle
      do i=inits,ends
        if(ibounce(1,i)==30 .and. ibounce(2,i)==0 .and. &
         ibounce(3,i)==1)then
           write(6,'(a,2i8,6f12.6)')'red ',nstep,idir,bc_flow(i),bc_rhoR(i), &
            bc_u(i),f03R(30,0,1),f04R(30,1,1)
        endif
!        if(ibounce(1,i)==16 .and. ibounce(2,i)==16 .and. &
!         ibounce(3,i)==0)then
!           write(6,'(a,2i8,6f12.6)')'blu ',nstep,idir,bc_flow(i),bc_rhoR(i), &
!            bc_rhoB(i),bc_w(i),f05B(16,16,0),f06B(16,16,0)
!        endif
      enddo
    enddo
    endif
#endif
  if(lsingle_fluid)return
  
  call apply_bounceback_pop_halfway(nstep,lColourG,bc_rhoB,bc_u,bc_v,bc_w, &
   rhoR,rhoB,fluidB,alphaB_CG)
  
#if 0
  if(mod(nstep,50)==0)then
     do idir=1,nbcdir
      inits=nbounce9dir(idir-1)+1
      ends=nbounce9dir(idir)
      if(inits>ends)cycle
      do i=inits,ends
        if(ibounce(1,i)==16 .and. ibounce(2,i)==16 .and. &
         ibounce(3,i)==33)then
           write(6,'(a,2i8,6f12.6)')'red ',nstep,idir,bc_flow(i),bc_rhoR(i), &
            bc_rhoB(i),bc_w(i),f06R(16,16,33),f05R(16,16,32)
        endif
        if(ibounce(1,i)==16 .and. ibounce(2,i)==16 .and. &
         ibounce(3,i)==0)then
           write(6,'(a,2i8,6f12.6)')'blu ',nstep,idir,bc_flow(i),bc_rhoR(i), &
            bc_rhoB(i),bc_w(i),f05B(16,16,0),f06B(16,16,0)
        endif
      enddo
    enddo
    endif
#endif
  
  return
  
 end subroutine driver_apply_bounceback_halfway_pop
 
!******************END PART TO MANAGE THE BOUNCEBACK********************



!************START PART TO MANAGE THE SHAN CHEN INTERACTION*************

 subroutine compute_psi_sc
 
!***********************************************************************
!     
!     LBsoft subroutine for computing the Shan Chen pseudo 
!     potential values for the fluid part
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification January 2019
!     
!***********************************************************************
 
  implicit none
  
  integer :: i,j,k
  integer, save :: io,jo,ko,ie,je,ke
  logical, save :: lfirst=.true.
  
  
  if(lfirst)then
    lfirst=.false.
    io=minx-nbuff
    ie=maxx+nbuff
    jo=miny-nbuff
    je=maxy+nbuff
    ko=minz-nbuff
    ke=maxz+nbuff
  endif
  
  if(lsingle_fluid)then

    where(isfluid(io:ie,jo:je,ko:ke)<3 .or. isfluid(io:ie,jo:je,ko:ke)>5)
      psiR(io:ie,jo:je,ko:ke)=rhoR(io:ie,jo:je,ko:ke)
    elsewhere
      psiR(io:ie,jo:je,ko:ke)=ZERO
    endwhere
    
  else

    where(isfluid(io:ie,jo:je,ko:ke)<3 .or. isfluid(io:ie,jo:je,ko:ke)>5)
      psiR(io:ie,jo:je,ko:ke)=rhoR(io:ie,jo:je,ko:ke)
      psiB(io:ie,jo:je,ko:ke)=rhoB(io:ie,jo:je,ko:ke)
    elsewhere
      psiR(io:ie,jo:je,ko:ke)=ZERO
      psiB(io:ie,jo:je,ko:ke)=ZERO
    endwhere
    
  endif
 end subroutine compute_psi_sc

 
 subroutine compute_sc_particle_interact_old(debug,nstep,iatm,lown, &
  lrotate,isub,jsub,ksub,nspheres,spherelists,spheredists,rdimx,rdimy, &
  rdimz,xx,yy,zz,vx,vy,vz,fx,fy,fz, A, ux,uy,uz,tx,ty,tz)
 
 
  implicit none
  
  logical, intent(in) :: lown,lrotate, debug
  integer, intent(in) :: nstep,iatm,isub,jsub,ksub,nspheres
  integer, allocatable, dimension(:,:), intent(in) :: spherelists
  real(kind=PRC), allocatable, dimension(:), intent(in) :: spheredists
  real(kind=PRC), intent(in) :: rdimx,rdimy,rdimz
  real(kind=PRC), intent(in) :: xx,yy,zz
  real(kind=PRC), intent(in) :: vx,vy,vz

#ifdef QUAD_FORCEINT
  real(kind=PRC*2), intent(inout) :: fx,fy,fz
  real(kind=PRC*2), intent(inout), optional :: tx,ty,tz
#else
  real(kind=PRC), intent(inout) :: fx,fy,fz
  real(kind=PRC), intent(inout), optional :: tx,ty,tz
#endif

#ifdef QUAD_FORCEINT
  real(kind=PRC*2), allocatable, dimension(:,:,:) :: A
#else
  real(kind=PRC), allocatable, dimension(:,:,:) :: A
#endif
  
  real(kind=PRC), intent(in), optional :: ux,uy,uz
  
  integer :: i,j,k,l,ii,jj,kk
  integer, save :: imin,imax,jmin,jmax,kmin,kmax
  logical, save :: lfirst=.true.
  real(kind=PRC) :: vxt,vyt,vzt,modr,ftx,fty,ftz,mytheta,factR,factB
  real(kind=PRC), dimension(3) :: rtemp,otemp,ftemp,urtemp,initf,ttemp
  character(len=120) :: mynamefile

  
  if(lfirst)then
    lfirst=.false.
    imin=minx-1
    imax=maxx+1
    jmin=miny-1
    jmax=maxy+1
    kmin=minz-1
    kmax=maxz+1
  endif
  
  if(lrotate)then
    otemp(1)=ux
    otemp(2)=uy
    otemp(3)=uz
    modr=modulvec(otemp)
    otemp=otemp/modr
  endif
  
  ttemp = ZERO

  if (nstep==1) then
    mynamefile=repeat(' ',120)
    mynamefile="forceInt_SC.atom"//write_fmtnumb(iatm)
    iounit(1) = 114
    call OpenLogFile(nstep, mynamefile, iounit(1))
  endif

  initf(1)=fx
  initf(2)=fy
  initf(3)=fz

  if(lsingle_fluid)then
    do l=1,nspheres
      i=isub+spherelists(1,l)
      j=jsub+spherelists(2,l)
      k=ksub+spherelists(3,l)
      ii=i
      jj=j
      kk=k
      !apply periodic conditions if necessary
      i=pimage(ixpbc,i,nx)
      j=pimage(iypbc,j,ny)
      k=pimage(izpbc,k,nz)
      if(lrotate)then
        rtemp(1)=real(ii,kind=PRC)-xx
        rtemp(2)=real(jj,kind=PRC)-yy
        rtemp(3)=real(kk,kind=PRC)-zz
        modr=modulvec(rtemp)
        rtemp(1:3)=rdimx/modr*rtemp(1:3)
        modr=modulvec(rtemp)
        urtemp(1:3)=rtemp(1:3)/modr
        mytheta=acos(dot(otemp,urtemp))
        factR= partR_SC*(TWO* &
         fcut(mytheta,theta_SC-devtheta_SC,theta_SC+devtheta_SC)-ONE)
      else
        factR=partR_SC*ONE
      endif


      !the fluid bounce back is local so I have to do it
      if(i==ii.and.j==jj.and.k==kk)then
        if(i<imin .or. i>imax)cycle
        if(j<jmin .or. j>jmax)cycle
        if(k<kmin .or. k>kmax)cycle
        psiR(i,j,k)=rhoR(i,j,k)*(ONE+factR)
      else
        if(i>=imin .and. i<=imax .and. j>=jmin .and. j<=jmax .and. &
         k>=kmin .and. k<=kmax)then
          psiR(i,j,k)=rhoR(i,j,k)*(ONE+factR)
        endif
        if(ii>=imin .and. ii<=imax .and. jj>=jmin .and. jj<=jmax .and. &
         kk>=kmin .and. kk<=kmax)then
          psiR(ii,jj,kk)=rhoR(ii,jj,kk)*(ONE+factR)
        endif
      endif

      if(i>=imin .and. i<=imax .and. j>=jmin .and. j<=jmax .and. &
       k>=kmin .and. k<=kmax)then
        call compute_grad_on_particle_old(debug,nstep,i,j,k,psiR,gradpsixR,gradpsiyR,gradpsizR)
        ftemp(1)=- pair_SC*psiR(i,j,k)*gradpsixR(i,j,k)
        ftemp(2)=- pair_SC*psiR(i,j,k)*gradpsiyR(i,j,k)
        ftemp(3)=- pair_SC*psiR(i,j,k)*gradpsizR(i,j,k)

#ifndef DEBUG_FORCEINT
        fx = fx + ftemp(1)
        fy = fy + ftemp(2)
        fz = fz + ftemp(3)
#else
        A(iatm, l, 1) = A(iatm, l, 1) + ftemp(1)
        A(iatm, l, 2) = A(iatm, l, 2) + ftemp(2)
        A(iatm, l, 3) = A(iatm, l, 3) + ftemp(3)
#endif

        if(lrotate)then
          ttemp(1) = xcross(rtemp,ftemp)
          ttemp(2) = ycross(rtemp,ftemp)
          ttemp(3) = zcross(rtemp,ftemp)
#ifndef DEBUG_FORCEINT
          tx = tx + ttemp(1)
          ty = ty + ttemp(2)
          tz = tz + ttemp(3)
#else
          A(iatm, l, 4) = A(iatm, l, 4) + ttemp(1)
          A(iatm, l, 5) = A(iatm, l, 5) + ttemp(2)
          A(iatm, l, 6) = A(iatm, l, 6) + ttemp(3)
#endif
        endif

        ! if (debug) then
        !  write (iounit(1),*) __FILE__,8457, i,j,k, &
        !     "ftemp=", ftemp, "ttemp=", ttemp
        ! endif

      endif
    enddo
    
  else
    do l=1,nspheres
      i=isub+spherelists(1,l)
      j=jsub+spherelists(2,l)
      k=ksub+spherelists(3,l)
      ii=i
      jj=j
      kk=k
      !apply periodic conditions if necessary
      i=pimage(ixpbc,i,nx)
      j=pimage(iypbc,j,ny)
      k=pimage(izpbc,k,nz)
      if(lrotate)then
        rtemp(1)=real(ii,kind=PRC)-xx
        rtemp(2)=real(jj,kind=PRC)-yy
        rtemp(3)=real(kk,kind=PRC)-zz
        modr=modulvec(rtemp)
        rtemp(1:3)=rdimx/modr*rtemp(1:3)
        modr=modulvec(rtemp)
        urtemp(1:3)=rtemp(1:3)/modr
        mytheta=acos(dot(otemp,urtemp))
        factR= partR_SC*(TWO* &
         fcut(mytheta,theta_SC-devtheta_SC,theta_SC+devtheta_SC)-ONE)
        factB= partB_SC*(TWO*(ONE- &
         fcut(mytheta,theta_SC-devtheta_SC,theta_SC+devtheta_SC))-ONE)
      else
        factR=partR_SC*ONE
        factB=partB_SC*ONE
      endif

      !the fluid bounce back is local so I have to do it
      if(i==ii.and.j==jj.and.k==kk)then
        if(i<imin .or. i>imax)cycle
        if(j<jmin .or. j>jmax)cycle
        if(k<kmin .or. k>kmax)cycle
        psiR(i,j,k)=rhoR(i,j,k)*(ONE+factR)
        psiB(i,j,k)=rhoB(i,j,k)*(ONE+factB)
      else
        if(i>=imin .and. i<=imax .and. j>=jmin .and. j<=jmax .and. &
         k>=kmin .and. k<=kmax)then
          psiR(i,j,k)=rhoR(i,j,k)*(ONE+factR)
          psiB(i,j,k)=rhoB(i,j,k)*(ONE+factB)
        endif
        if(ii>=imin .and. ii<=imax .and. jj>=jmin .and. jj<=jmax .and. &
         kk>=kmin .and. kk<=kmax)then
          psiR(ii,jj,kk)=rhoR(ii,jj,kk)*(ONE+factR)
          psiB(ii,jj,kk)=rhoB(ii,jj,kk)*(ONE+factB)
        endif
      endif


      if(i>=imin .and. i<=imax .and. j>=jmin .and. j<=jmax .and. &
       k>=kmin .and. k<=kmax)then
        call compute_grad_on_particle_old(debug,nstep,i,j,k,psiR,gradpsixR,gradpsiyR,gradpsizR)
        call compute_grad_on_particle_old(debug,nstep,i,j,k,psiB,gradpsixB,gradpsiyB,gradpsizB)
        ftemp(1)=- pair_SC*psiR(i,j,k)*gradpsixB(i,j,k) - &
         pair_SC*psiB(i,j,k)*gradpsixR(i,j,k)
        ftemp(2)=- pair_SC*psiR(i,j,k)*gradpsiyB(i,j,k) - &
         pair_SC*psiB(i,j,k)*gradpsiyR(i,j,k)
        ftemp(3)=- pair_SC*psiR(i,j,k)*gradpsizB(i,j,k) - &
         pair_SC*psiB(i,j,k)*gradpsizR(i,j,k)
#ifndef DEBUG_FORCEINT
        fx = fx + ftemp(1)
        fy = fy + ftemp(2)
        fz = fz + ftemp(3)
#else
        A(iatm, l, 1) = A(iatm, l, 1) + ftemp(1)
        A(iatm, l, 2) = A(iatm, l, 2) + ftemp(2)
        A(iatm, l, 3) = A(iatm, l, 3) + ftemp(3)
#endif
        if(lrotate)then
          ttemp(1) = xcross(rtemp,ftemp)
          ttemp(2) = ycross(rtemp,ftemp)
          ttemp(3) = zcross(rtemp,ftemp)
#ifndef DEBUG_FORCEINT
          tx = tx + ttemp(1)
          ty = ty + ttemp(2)
          tz = tz + ttemp(3)
#else
          A(iatm, l, 4) = A(iatm, l, 4) + ttemp(1)
          A(iatm, l, 5) = A(iatm, l, 5) + ttemp(2)
          A(iatm, l, 6) = A(iatm, l, 6) + ttemp(3)
#endif
        endif

        if (nstep==1) then
         write (iounit(1),*) __FILE__,0, i,j,k, &
            "ftemp=", ftemp, "ttemp=", ttemp, &
            "psiR",psiR(i,j,k),"gradR",gradpsixR(i,j,k),gradpsiyR(i,j,k),gradpsizR(i,j,k), &
            "psiB",psiB(i,j,k),"gradB",gradpsixB(i,j,k),gradpsiyB(i,j,k),gradpsizB(i,j,k), &
            "l",l
        endif

      endif

    enddo
  endif
  
  if (nstep==1) then
    close(iounit(1))
  endif
 end subroutine compute_sc_particle_interact_old
  
 subroutine compute_sc_particle_interact_phase1(debug,nstep,iatm,lown, &
  lrotate,isub,jsub,ksub,nspheres,spherelists,spheredists,rdimx,rdimy, &
  rdimz,xx,yy,zz,vx,vy,vz,fx,fy,fz, A, ux,uy,uz,tx,ty,tz)
 
!***********************************************************************
!     
!     LBsoft subroutine for computing the Shan Chen force acting on 
!     particles and the pseudo potential values for the fluid part
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification January 2019
!     
!***********************************************************************
 
  implicit none
  
  logical, intent(in) :: lown,lrotate, debug
  integer, intent(in) :: nstep,iatm,isub,jsub,ksub,nspheres
  integer, allocatable, dimension(:,:), intent(in) :: spherelists
  real(kind=PRC), allocatable, dimension(:), intent(in) :: spheredists
  real(kind=PRC), intent(in) :: rdimx,rdimy,rdimz
  real(kind=PRC), intent(in) :: xx,yy,zz
  real(kind=PRC), intent(in) :: vx,vy,vz

#ifdef QUAD_FORCEINT
  real(kind=PRC*2), intent(inout) :: fx,fy,fz
  real(kind=PRC*2), intent(inout), optional :: tx,ty,tz
#else
  real(kind=PRC), intent(inout) :: fx,fy,fz
  real(kind=PRC), intent(inout), optional :: tx,ty,tz
#endif

#ifdef QUAD_FORCEINT
  real(kind=PRC*2), allocatable, dimension(:,:,:) :: A
#else
  real(kind=PRC), allocatable, dimension(:,:,:) :: A
#endif
  
  real(kind=PRC), intent(in), optional :: ux,uy,uz
  
  integer :: i,j,k,l,ii,jj,kk
  integer, save :: imin,imax,jmin,jmax,kmin,kmax
  logical, save :: lfirst=.true.
  real(kind=PRC) :: vxt,vyt,vzt,modr,ftx,fty,ftz,mytheta,factR,factB
  real(kind=PRC), dimension(3) :: rtemp,otemp,ftemp,urtemp,initf,ttemp
  character(len=120) :: mynamefile

  
  if(lfirst)then
    lfirst=.false.
    imin=minx-1
    imax=maxx+1
    jmin=miny-1
    jmax=maxy+1
    kmin=minz-1
    kmax=maxz+1
  endif
  
  if(lrotate)then
    otemp(1)=ux
    otemp(2)=uy
    otemp(3)=uz
    modr=modulvec(otemp)
    otemp=otemp/modr
  endif
  
  ttemp = ZERO

  initf(1)=fx
  initf(2)=fy
  initf(3)=fz

  if(lsingle_fluid)then
    do l=1,nspheres
      i=isub+spherelists(1,l)
      j=jsub+spherelists(2,l)
      k=ksub+spherelists(3,l)
      ii=i
      jj=j
      kk=k
      !apply periodic conditions if necessary
      i=pimage(ixpbc,i,nx)
      j=pimage(iypbc,j,ny)
      k=pimage(izpbc,k,nz)
      if(lrotate)then
        rtemp(1)=real(ii,kind=PRC)-xx
        rtemp(2)=real(jj,kind=PRC)-yy
        rtemp(3)=real(kk,kind=PRC)-zz
        modr=modulvec(rtemp)
        rtemp(1:3)=rdimx/modr*rtemp(1:3)
        modr=modulvec(rtemp)
        urtemp(1:3)=rtemp(1:3)/modr
        mytheta=acos(dot(otemp,urtemp))
        factR= partR_SC*(TWO* &
         fcut(mytheta,theta_SC-devtheta_SC,theta_SC+devtheta_SC)-ONE)
      else
        factR=partR_SC*ONE
      endif


      !the fluid bounce back is local so I have to do it
      if(i==ii.and.j==jj.and.k==kk)then
        if(i<imin .or. i>imax)cycle
        if(j<jmin .or. j>jmax)cycle
        if(k<kmin .or. k>kmax)cycle
        psiR(i,j,k) = psiR(i,j,k) + psiR(i,j,k)*factR
      else
        if(i>=imin .and. i<=imax .and. j>=jmin .and. j<=jmax .and. &
         k>=kmin .and. k<=kmax)then
          psiR(i,j,k) = psiR(i,j,k) + psiR(i,j,k)*factR
        endif
        if(ii>=imin .and. ii<=imax .and. jj>=jmin .and. jj<=jmax .and. &
         kk>=kmin .and. kk<=kmax)then
          psiR(ii,jj,kk) = psiR(ii,jj,kk) + psiR(ii,jj,kk)*factR
        endif
      endif
    enddo
    
  else
    do l=1,nspheres
      i=isub+spherelists(1,l)
      j=jsub+spherelists(2,l)
      k=ksub+spherelists(3,l)
      ii=i
      jj=j
      kk=k
      !apply periodic conditions if necessary
      i=pimage(ixpbc,i,nx)
      j=pimage(iypbc,j,ny)
      k=pimage(izpbc,k,nz)
      if(lrotate)then
        rtemp(1)=real(ii,kind=PRC)-xx
        rtemp(2)=real(jj,kind=PRC)-yy
        rtemp(3)=real(kk,kind=PRC)-zz
        modr=modulvec(rtemp)
        rtemp(1:3)=rdimx/modr*rtemp(1:3)
        modr=modulvec(rtemp)
        urtemp(1:3)=rtemp(1:3)/modr
        mytheta=acos(dot(otemp,urtemp))
        factR= partR_SC*(TWO* &
         fcut(mytheta,theta_SC-devtheta_SC,theta_SC+devtheta_SC)-ONE)
        factB= partB_SC*(TWO*(ONE- &
         fcut(mytheta,theta_SC-devtheta_SC,theta_SC+devtheta_SC))-ONE)
      else
        factR=partR_SC*ONE
        factB=partB_SC*ONE
      endif

      !the fluid bounce back is local so I have to do it
      if(i==ii.and.j==jj.and.k==kk)then
        if(i<imin .or. i>imax)cycle
        if(j<jmin .or. j>jmax)cycle
        if(k<kmin .or. k>kmax)cycle
        psiR(i,j,k) = psiR(i,j,k) + psiR(i,j,k)*factR
        psiB(i,j,k) = psiB(i,j,k) + psiB(i,j,k)*factB
      else
        if(i>=imin .and. i<=imax .and. j>=jmin .and. j<=jmax .and. &
         k>=kmin .and. k<=kmax)then
          psiR(i,j,k) = psiR(i,j,k) + psiR(i,j,k)*factR
          psiB(i,j,k) = psiB(i,j,k) + psiB(i,j,k)*factB
        endif
        if(ii>=imin .and. ii<=imax .and. jj>=jmin .and. jj<=jmax .and. &
         kk>=kmin .and. kk<=kmax)then
          psiR(ii,jj,kk) = psiR(ii,jj,kk) + psiR(ii,jj,kk)*factR
          psiB(ii,jj,kk) = psiB(ii,jj,kk) + psiB(ii,jj,kk)*factB
        endif
      endif
    enddo
  endif
  
 end subroutine compute_sc_particle_interact_phase1

 subroutine compute_sc_particle_interact_phase2(debug,nstep,iatm,lown, &
  lrotate,isub,jsub,ksub,nspheres,spherelists,spheredists,rdimx,rdimy, &
  rdimz,xx,yy,zz,vx,vy,vz,fx,fy,fz, A, ux,uy,uz,tx,ty,tz)
 
!***********************************************************************
!     
!     LBsoft subroutine for computing the Shan Chen force acting on 
!     particles and the pseudo potential values for the fluid part
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification January 2019
!     
!***********************************************************************
 
  implicit none
  
  logical, intent(in) :: lown,lrotate, debug
  integer, intent(in) :: nstep,iatm,isub,jsub,ksub,nspheres
  integer, allocatable, dimension(:,:), intent(in) :: spherelists
  real(kind=PRC), allocatable, dimension(:), intent(in) :: spheredists
  real(kind=PRC), intent(in) :: rdimx,rdimy,rdimz
  real(kind=PRC), intent(in) :: xx,yy,zz
  real(kind=PRC), intent(in) :: vx,vy,vz

#ifdef QUAD_FORCEINT
  real(kind=PRC*2), intent(inout) :: fx,fy,fz
  real(kind=PRC*2), intent(inout), optional :: tx,ty,tz
#else
  real(kind=PRC), intent(inout) :: fx,fy,fz
  real(kind=PRC), intent(inout), optional :: tx,ty,tz
#endif

#ifdef QUAD_FORCEINT
  real(kind=PRC*2), allocatable, dimension(:,:,:) :: A
#else
  real(kind=PRC), allocatable, dimension(:,:,:) :: A
#endif
  
  real(kind=PRC), intent(in), optional :: ux,uy,uz
  
  integer :: i,j,k,l,ii,jj,kk
  integer, save :: imin,imax,jmin,jmax,kmin,kmax
  logical, save :: lfirst=.true.
  real(kind=PRC) :: vxt,vyt,vzt,modr,ftx,fty,ftz,mytheta,factR,factB
  real(kind=PRC), dimension(3) :: rtemp,otemp,ftemp,urtemp,initf,ttemp
  character(len=120) :: mynamefile

  
  if(lfirst)then
    lfirst=.false.
    imin=minx-1
    imax=maxx+1
    jmin=miny-1
    jmax=maxy+1
    kmin=minz-1
    kmax=maxz+1
  endif
  
  if(lrotate)then
    otemp(1)=ux
    otemp(2)=uy
    otemp(3)=uz
    modr=modulvec(otemp)
    otemp=otemp/modr
  endif
  
  ttemp = ZERO

  ! if (debug) then
  !   mynamefile=repeat(' ',120)
  !   mynamefile="forceInt_SC.atom"//write_fmtnumb(iatm)
  !   iounit(1) = 114
  !   call OpenLogFile(nstep, mynamefile, iounit(1))

  !   mynamefile=repeat(' ',120)
  !   mynamefile="forceInt_SC_grad.atom"//write_fmtnumb(iatm)
  !   iounit(2) = 115
  !   call OpenLogFile(nstep, mynamefile, iounit(2))
  ! endif

  initf(1)=fx
  initf(2)=fy
  initf(3)=fz

  if(lsingle_fluid)then
    do l=1,nspheres
      i=isub+spherelists(1,l)
      j=jsub+spherelists(2,l)
      k=ksub+spherelists(3,l)
      ii=i
      jj=j
      kk=k
      !apply periodic conditions if necessary
      i=pimage(ixpbc,i,nx)
      j=pimage(iypbc,j,ny)
      k=pimage(izpbc,k,nz)
      if(lrotate)then
        rtemp(1)=real(ii,kind=PRC)-xx
        rtemp(2)=real(jj,kind=PRC)-yy
        rtemp(3)=real(kk,kind=PRC)-zz
        modr=modulvec(rtemp)
        rtemp(1:3)=rdimx/modr*rtemp(1:3)
      endif

      if(i>=minx .and. i<=maxx .and. j>=miny .and. j<=maxy .and. k>=minz .and. k<=maxz)then
        call compute_grad_on_particle1Fl(lrotate, i,j,k,psiR,rtemp, debug,iatm, A, l, &
                fx,fy,fz, tx,ty,tz)
      endif
    enddo
    
  else
    do l=1,nspheres
      i=isub+spherelists(1,l)
      j=jsub+spherelists(2,l)
      k=ksub+spherelists(3,l)
      ii=i
      jj=j
      kk=k
      !apply periodic conditions if necessary
      i=pimage(ixpbc,i,nx)
      j=pimage(iypbc,j,ny)
      k=pimage(izpbc,k,nz)
      if(lrotate)then
        rtemp(1)=real(ii,kind=PRC)-xx
        rtemp(2)=real(jj,kind=PRC)-yy
        rtemp(3)=real(kk,kind=PRC)-zz
        modr=modulvec(rtemp)
        rtemp(1:3)=rdimx/modr*rtemp(1:3)
      endif

      if(i>=minx .and. i<=maxx .and. j>=miny .and. j<=maxy .and. k>=minz .and. k<=maxz)then
        call compute_grad_on_particle2Fl(lrotate, i,j,k,psiR,psiB, rtemp, debug,iatm, A, l, &
                fx,fy,fz, tx,ty,tz)
      endif

    enddo
  endif
  
  ! if (debug) then
  !   close(iounit(1))
  !   close(iounit(2))
  ! endif
  
 end subroutine compute_sc_particle_interact_phase2
 
 subroutine compute_grad_on_particle1Fl(lrotate, i,j,k,psiR,rtemp, debug,iatm, A, l, &
                 fx,fy,fz,tx,ty,tz)

!***********************************************************************
!
!     LBsoft subroutine to computing the Shan Chen force acting
!     on a particle surface node including the pbc
!
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification January 2019
!
!***********************************************************************

  implicit none
  logical, intent(in) :: debug, lrotate
  integer, intent(in) :: iatm,l, i,j,k
  real(kind=PRC), allocatable, dimension(:,:,:)  :: psiR
  real(kind=PRC), intent(in)  :: rtemp(3)
#ifdef QUAD_FORCEINT
  real(kind=PRC*2), allocatable, dimension(:,:,:) :: A
#else
  real(kind=PRC), allocatable, dimension(:,:,:) :: A
#endif
#ifdef QUAD_FORCEINT
  real(kind=PRC*2), intent(inout) :: fx,fy,fz
  real(kind=PRC*2), intent(inout), optional :: tx,ty,tz
#else
  real(kind=PRC), intent(inout) :: fx,fy,fz
  real(kind=PRC), intent(inout), optional :: tx,ty,tz
#endif

  real(kind=PRC)  :: gradpsixR,gradpsiyR,gradpsizR
  real(kind=PRC)  :: ftemp(3), ttemp(3)
  integer :: ii,jj,kk,io,jo,ko, iloop,indlow,indhig


  ttemp = ZERO

  do iloop = 1, 9
  indlow = iloop*2 - 1
  indhig = iloop*2

  !force on particle fx fy fz
  !eq. 11.2 from page 437 Kruger's book "the lattice boltzmann method"

  ii=i+ex(indlow)
  jj=j+ey(indlow)
  kk=k+ez(indlow)
!  ii=pimage(ixpbc,ii,nx)
!  jj=pimage(iypbc,jj,ny)
!  kk=pimage(izpbc,kk,nz)
  io=i+ex(indhig)
  jo=j+ey(indhig)
  ko=k+ez(indhig)
!  io=pimage(ixpbc,io,nx)
!  jo=pimage(iypbc,jo,ny)
!  ko=pimage(izpbc,ko,nz)
  if(ii>=minx-1 .and. ii<=maxx+1 .and. jj>=miny-1 .and. jj<=maxy+1 .and. &
     kk>=minz-1 .and. kk<=maxz+1)then
   if(io>=minx-1 .and. io<=maxx+1 .and. jo>=miny-1 .and. jo<=maxy+1 .and. &
     ko>=minz-1 .and. ko<=maxz+1)then
    if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
      gradpsixR = psiR(ii,jj,kk)*p(indlow)*ex(indlow)
      gradpsiyR = psiR(ii,jj,kk)*p(indlow)*ey(indlow)
      gradpsizR = psiR(ii,jj,kk)*p(indlow)*ez(indlow)

      ftemp(1)=- pair_SC*psiR(i,j,k)*gradpsixR
      ftemp(2)=- pair_SC*psiR(i,j,k)*gradpsiyR
      ftemp(3)=- pair_SC*psiR(i,j,k)*gradpsizR
#ifndef DEBUG_FORCEINT
      fx = fx + ftemp(1)
      fy = fy + ftemp(2)
      fz = fz + ftemp(3)
#else
      A(iatm, l, 1) = A(iatm, l, 1) + ftemp(1)
      A(iatm, l, 2) = A(iatm, l, 2) + ftemp(2)
      A(iatm, l, 3) = A(iatm, l, 3) + ftemp(3)
#endif
      if(lrotate)then
          ttemp(1) = xcross(rtemp,ftemp)
          ttemp(2) = ycross(rtemp,ftemp)
          ttemp(3) = zcross(rtemp,ftemp)
#ifndef DEBUG_FORCEINT
          tx = tx + ttemp(1)
          ty = ty + ttemp(2)
          tz = tz + ttemp(3)
#else
          A(iatm, l, 4) = A(iatm, l, 4) + ttemp(1)
          A(iatm, l, 5) = A(iatm, l, 5) + ttemp(2)
          A(iatm, l, 6) = A(iatm, l, 6) + ttemp(3)
#endif
      endif

      ! if (debug) then
      !   write (iounit(1),*) __FILE__,indlow, i,j,k, &
      !       "ftemp=", ftemp, "ttemp=", ttemp, &
      !       "psiR",psiR(i,j,k),"gradR",gradpsixR,gradpsiyR,gradpsizR
      !   write (iounit(2),*) indlow, i,j,k, gradpsixR,gradpsiyR,gradpsizR
      ! endif
    endif
  endif
  endif

  ii=i+ex(indhig)
  jj=j+ey(indhig)
  kk=k+ez(indhig)
!  ii=pimage(ixpbc,ii,nx)
!  jj=pimage(iypbc,jj,ny)
!  kk=pimage(izpbc,kk,nz)
  io=i+ex(indlow)
  jo=j+ey(indlow)
  ko=k+ez(indlow)
!  io=pimage(ixpbc,io,nx)
!  jo=pimage(iypbc,jo,ny)
!  ko=pimage(izpbc,ko,nz)
  if(ii>=minx-1 .and. ii<=maxx+1 .and. jj>=miny-1 .and. jj<=maxy+1 .and. &
     kk>=minz-1 .and. kk<=maxz+1)then
   if(io>=minx-1 .and. io<=maxx+1 .and. jo>=miny-1 .and. jo<=maxy+1 .and. &
     ko>=minz-1 .and. ko<=maxz+1)then
    if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
      gradpsixR = psiR(ii,jj,kk)*p(indhig)*ex(indhig)
      gradpsiyR = psiR(ii,jj,kk)*p(indhig)*ey(indhig)
      gradpsizR = psiR(ii,jj,kk)*p(indhig)*ez(indhig)

      ftemp(1)=- pair_SC*psiR(i,j,k)*gradpsixR
      ftemp(2)=- pair_SC*psiR(i,j,k)*gradpsiyR
      ftemp(3)=- pair_SC*psiR(i,j,k)*gradpsizR
#ifndef DEBUG_FORCEINT
      fx = fx + ftemp(1)
      fy = fy + ftemp(2)
      fz = fz + ftemp(3)
#else
      A(iatm, l, 1) = A(iatm, l, 1) + ftemp(1)
      A(iatm, l, 2) = A(iatm, l, 2) + ftemp(2)
      A(iatm, l, 3) = A(iatm, l, 3) + ftemp(3)
#endif
      if(lrotate)then
          ttemp(1) = xcross(rtemp,ftemp)
          ttemp(2) = ycross(rtemp,ftemp)
          ttemp(3) = zcross(rtemp,ftemp)
#ifndef DEBUG_FORCEINT
          tx = tx + ttemp(1)
          ty = ty + ttemp(2)
          tz = tz + ttemp(3)
#else
          A(iatm, l, 4) = A(iatm, l, 4) + ttemp(1)
          A(iatm, l, 5) = A(iatm, l, 5) + ttemp(2)
          A(iatm, l, 6) = A(iatm, l, 6) + ttemp(3)
#endif
      endif

      ! if (debug) then
      !   write (iounit(1),*) __FILE__,indhig, i,j,k, &
      !       "ftemp=", ftemp, "ttemp=", ttemp, &
      !       "psiR",psiR(i,j,k),"gradR",gradpsixR,gradpsiyR,gradpsizR
      !   write (iounit(2),*) indhig, i,j,k, gradpsixR,gradpsiyR,gradpsizR
      ! endif
    endif
  endif
  endif

 enddo

 end subroutine compute_grad_on_particle1Fl

 subroutine compute_grad_on_particle2Fl(lrotate, i,j,k,psiR,psiB, rtemp, debug,iatm, A, l, &
                 fx,fy,fz,tx,ty,tz)
 
!***********************************************************************
!     
!     LBsoft subroutine to computing the Shan Chen force acting
!     on a particle surface node including the pbc
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification January 2019
!     
!***********************************************************************
 
  implicit none
  logical, intent(in) :: debug, lrotate
  integer, intent(in) :: iatm,l, i,j,k
  real(kind=PRC), allocatable, dimension(:,:,:)  :: psiR,psiB
  real(kind=PRC), intent(in)  :: rtemp(3)
#ifdef QUAD_FORCEINT
  real(kind=PRC*2), allocatable, dimension(:,:,:) :: A
#else
  real(kind=PRC), allocatable, dimension(:,:,:) :: A
#endif
#ifdef QUAD_FORCEINT
  real(kind=PRC*2), intent(inout) :: fx,fy,fz
  real(kind=PRC*2), intent(inout), optional :: tx,ty,tz
#else
  real(kind=PRC), intent(inout) :: fx,fy,fz
  real(kind=PRC), intent(inout), optional :: tx,ty,tz
#endif

  real(kind=PRC)  :: gradpsixR,gradpsiyR,gradpsizR, gradpsixB,gradpsiyB,gradpsizB
  real(kind=PRC)  :: ftemp(3), ttemp(3)
  integer :: ii,jj,kk,io,jo,ko, iloop,indlow,indhig
  

  ttemp = ZERO

  do iloop = 1, 9
  indlow = iloop*2 - 1
  indhig = iloop*2
  
  !force on particle fx fy fz
  !eq. 11.2 from page 437 Kruger's book "the lattice boltzmann method"
    
  ii=i+ex(indlow)
  jj=j+ey(indlow)
  kk=k+ez(indlow)
!  ii=pimage(ixpbc,ii,nx)
!  jj=pimage(iypbc,jj,ny)
!  kk=pimage(izpbc,kk,nz)
  io=i+ex(indhig)
  jo=j+ey(indhig)
  ko=k+ez(indhig)
!  io=pimage(ixpbc,io,nx)
!  jo=pimage(iypbc,jo,ny)
!  ko=pimage(izpbc,ko,nz)
  if(ii>=minx-1 .and. ii<=maxx+1 .and. jj>=miny-1 .and. jj<=maxy+1 .and. &
     kk>=minz-1 .and. kk<=maxz+1)then
   if(io>=minx-1 .and. io<=maxx+1 .and. jo>=miny-1 .and. jo<=maxy+1 .and. &
     ko>=minz-1 .and. ko<=maxz+1)then
    if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
      gradpsixR = psiR(ii,jj,kk)*p(indlow)*ex(indlow)
      gradpsiyR = psiR(ii,jj,kk)*p(indlow)*ey(indlow)
      gradpsizR = psiR(ii,jj,kk)*p(indlow)*ez(indlow)
      gradpsixB = psiB(ii,jj,kk)*p(indlow)*ex(indlow)
      gradpsiyB = psiB(ii,jj,kk)*p(indlow)*ey(indlow)
      gradpsizB = psiB(ii,jj,kk)*p(indlow)*ez(indlow)

      ftemp(1)=- pair_SC*(psiR(i,j,k)*gradpsixB + psiB(i,j,k)*gradpsixR)
      ftemp(2)=- pair_SC*(psiR(i,j,k)*gradpsiyB + psiB(i,j,k)*gradpsiyR)
      ftemp(3)=- pair_SC*(psiR(i,j,k)*gradpsizB + psiB(i,j,k)*gradpsizR)
#ifndef DEBUG_FORCEINT
      fx = fx + ftemp(1)
      fy = fy + ftemp(2)
      fz = fz + ftemp(3)
#else
      A(iatm, l, 1) = A(iatm, l, 1) + ftemp(1)
      A(iatm, l, 2) = A(iatm, l, 2) + ftemp(2)
      A(iatm, l, 3) = A(iatm, l, 3) + ftemp(3)
#endif
      if(lrotate)then
          ttemp(1) = xcross(rtemp,ftemp)
          ttemp(2) = ycross(rtemp,ftemp)
          ttemp(3) = zcross(rtemp,ftemp)
#ifndef DEBUG_FORCEINT
          tx = tx + ttemp(1)
          ty = ty + ttemp(2)
          tz = tz + ttemp(3)
#else
          A(iatm, l, 4) = A(iatm, l, 4) + ttemp(1)
          A(iatm, l, 5) = A(iatm, l, 5) + ttemp(2)
          A(iatm, l, 6) = A(iatm, l, 6) + ttemp(3)
#endif
      endif

      ! if (debug) then
      !   write (iounit(1),*) __FILE__,indlow, i,j,k, &
      !       "ftemp=", real(ftemp,kind=PRC*2), "ttemp=", real(ttemp,kind=PRC*2), &
      !       "psiR",psiR(i,j,k),"gradR",gradpsixR,gradpsiyR,gradpsizR, &
      !       "psiB",psiB(i,j,k),"gradB",gradpsixB,gradpsiyB,gradpsizB, &
      !       "l",l
      !   write (iounit(2),*) indlow, i,j,k, gradpsixR,gradpsiyR,gradpsizR,rtemp
      !   write (iounit(2),*) indlow, i,j,k, gradpsixB,gradpsiyB,gradpsizB,rtemp
      ! endif
    endif
  endif
  endif
    
  ii=i+ex(indhig)
  jj=j+ey(indhig)
  kk=k+ez(indhig)
!  ii=pimage(ixpbc,ii,nx)
!  jj=pimage(iypbc,jj,ny)
!  kk=pimage(izpbc,kk,nz)
  io=i+ex(indlow)
  jo=j+ey(indlow)
  ko=k+ez(indlow)
!  io=pimage(ixpbc,io,nx)
!  jo=pimage(iypbc,jo,ny)
!  ko=pimage(izpbc,ko,nz)
  if(ii>=minx-1 .and. ii<=maxx+1 .and. jj>=miny-1 .and. jj<=maxy+1 .and. &
     kk>=minz-1 .and. kk<=maxz+1)then
   if(io>=minx-1 .and. io<=maxx+1 .and. jo>=miny-1 .and. jo<=maxy+1 .and. &
     ko>=minz-1 .and. ko<=maxz+1)then
    if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
      gradpsixR = psiR(ii,jj,kk)*p(indhig)*ex(indhig)
      gradpsiyR = psiR(ii,jj,kk)*p(indhig)*ey(indhig)
      gradpsizR = psiR(ii,jj,kk)*p(indhig)*ez(indhig)
      gradpsixB = psiB(ii,jj,kk)*p(indhig)*ex(indhig)
      gradpsiyB = psiB(ii,jj,kk)*p(indhig)*ey(indhig)
      gradpsizB = psiB(ii,jj,kk)*p(indhig)*ez(indhig)

      ftemp(1)=- pair_SC*(psiR(i,j,k)*gradpsixB + psiB(i,j,k)*gradpsixR)
      ftemp(2)=- pair_SC*(psiR(i,j,k)*gradpsiyB + psiB(i,j,k)*gradpsiyR)
      ftemp(3)=- pair_SC*(psiR(i,j,k)*gradpsizB + psiB(i,j,k)*gradpsizR)
#ifndef DEBUG_FORCEINT
      fx = fx + ftemp(1)
      fy = fy + ftemp(2)
      fz = fz + ftemp(3)
#else
      A(iatm, l, 1) = A(iatm, l, 1) + ftemp(1)
      A(iatm, l, 2) = A(iatm, l, 2) + ftemp(2)
      A(iatm, l, 3) = A(iatm, l, 3) + ftemp(3)
#endif
      if(lrotate)then
          ttemp(1) = xcross(rtemp,ftemp)
          ttemp(2) = ycross(rtemp,ftemp)
          ttemp(3) = zcross(rtemp,ftemp)
#ifndef DEBUG_FORCEINT
          tx = tx + ttemp(1)
          ty = ty + ttemp(2)
          tz = tz + ttemp(3)
#else
          A(iatm, l, 4) = A(iatm, l, 4) + ttemp(1)
          A(iatm, l, 5) = A(iatm, l, 5) + ttemp(2)
          A(iatm, l, 6) = A(iatm, l, 6) + ttemp(3)
#endif
      endif

      ! if (debug) then
      !   write (iounit(1),*) __FILE__,indhig, i,j,k, &
      !       "ftemp=", real(ftemp,kind=PRC*2), "ttemp=", real(ttemp,kind=PRC*2), &
      !       "psiR",psiR(i,j,k),"gradR",gradpsixR,gradpsiyR,gradpsizR, &
      !       "psiB",psiB(i,j,k),"gradB",gradpsixB,gradpsiyB,gradpsizB, &
      !       "l",l
      !   write (iounit(2),*) indhig, i,j,k, gradpsixR,gradpsiyR,gradpsizR,rtemp
      !   write (iounit(2),*) indhig, i,j,k, gradpsixB,gradpsiyB,gradpsizB,rtemp
      ! endif
    endif
  endif
  endif

 enddo

 end subroutine compute_grad_on_particle2Fl


 subroutine compute_grad_on_particle_old(debug,nstep,i,j,k,psisub,mygradx,mygrady,mygradz)

!***********************************************************************
!
!     LBsoft subroutine to computing the Shan Chen force acting
!     on a particle surface node including the pbc
!
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification January 2019
!
!***********************************************************************

  implicit none
  logical, intent(in) :: debug
  integer, intent(in) :: nstep,i,j,k
  real(kind=PRC), allocatable, dimension(:,:,:)  :: psisub
  real(kind=PRC), allocatable, dimension(:,:,:)  :: mygradx,mygrady,mygradz

  integer :: ii,jj,kk,io,jo,ko

  mygradx(i,j,k)=ZERO
  mygrady(i,j,k)=ZERO
  mygradz(i,j,k)=ZERO


  !force on particle fx fy fz
  !eq. 11.2 from page 437 Kruger's book "the lattice boltzmann method"

  ii=i+ex(1)
  jj=j
  kk=k
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i+ex(2)
  jo=j
  ko=k
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      mygradx(i,j,k)=mygradx(i,j,k)+psisub(ii,jj,kk)*p(1)*dex(1)
      ! if (debug) then
      !   write (iounit(2),*) 1, i,j,k,psisub(ii,jj,kk)*p(1)*dex(1),0,0
      ! endif
    endif
  endif

  ii=i+ex(2)
  jj=j
  kk=k
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i+ex(1)
  jo=j
  ko=k
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      mygradx(i,j,k)=mygradx(i,j,k)+psisub(ii,jj,kk)*p(2)*dex(2)
      ! if (debug) then
      !   write (iounit(2),*) 2, i,j,k, psisub(ii,jj,kk)*p(2)*dex(2),0,0
      ! endif
    endif
  endif

  ii=i
  jj=j+ey(3)
  kk=k
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i
  jo=j+ey(4)
  ko=k
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      mygrady(i,j,k)=mygrady(i,j,k)+psisub(ii,jj,kk)*p(3)*dey(3)
      ! if (debug) then
      !   write (iounit(2),*) 3, i,j,k, 0,psisub(ii,jj,kk)*p(3)*dey(3),0
      ! endif
    endif
  endif

  ii=i
  jj=j+ey(4)
  kk=k
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i
  jo=j+ey(3)
  ko=k
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      mygrady(i,j,k)=mygrady(i,j,k)+psisub(ii,jj,kk)*p(4)*dey(4)
      ! if (debug) then
      !   write (iounit(2),*) 4, i,j,k, 0,psisub(ii,jj,kk)*p(4)*dey(4),0
      ! endif
    endif
  endif

  ii=i
  jj=j
  kk=k+ez(5)
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i
  jo=j
  ko=k+ez(6)
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      mygradz(i,j,k)=mygradz(i,j,k)+psisub(ii,jj,kk)*p(5)*dez(5)
      ! if (debug) then
      !   write (iounit(2),*) 5, i,j,k, 0,0,psisub(ii,jj,kk)*p(5)*dez(5)
      ! endif
    endif
  endif

  ii=i
  jj=j
  kk=k+ez(6)
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i
  jo=j
  ko=k+ez(5)
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      mygradz(i,j,k)=mygradz(i,j,k)+psisub(ii,jj,kk)*p(6)*dez(6)
      ! if (debug) then
      !   write (iounit(2),*) 6, i,j,k, 0,0,psisub(ii,jj,kk)*p(6)*dez(6)
      ! endif
    endif
  endif

  ii=i+ex(7)
  jj=j+ey(7)
  kk=k
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i+ex(8)
  jo=j+ey(8)
  ko=k
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      mygradx(i,j,k)=mygradx(i,j,k)+psisub(ii,jj,kk)*p(7)*dex(7)
      mygrady(i,j,k)=mygrady(i,j,k)+psisub(ii,jj,kk)*p(7)*dey(7)
      ! if (debug) then
      !   write (iounit(2),*) 7, i,j,k, psisub(ii,jj,kk)*p(7)*dex(7),psisub(ii,jj,kk)*p(7)*dey(7),0
      ! endif
    endif
  endif

  ii=i+ex(8)
  jj=j+ey(8)
  kk=k
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i+ex(7)
  jo=j+ey(7)
  ko=k
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      mygradx(i,j,k)=mygradx(i,j,k)+psisub(ii,jj,kk)*p(8)*dex(8)
      mygrady(i,j,k)=mygrady(i,j,k)+psisub(ii,jj,kk)*p(8)*dey(8)
      ! if (debug) then
      !   write (iounit(2),*) 8, i,j,k, psisub(ii,jj,kk)*p(8)*dex(8),psisub(ii,jj,kk)*p(8)*dey(8),0
      ! endif
    endif
  endif

  ii=i+ex(9)
  jj=j+ey(9)
  kk=k
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i+ex(10)
  jo=j+ey(10)
  ko=k
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      mygradx(i,j,k)=mygradx(i,j,k)+psisub(ii,jj,kk)*p(9)*dex(9)
      mygrady(i,j,k)=mygrady(i,j,k)+psisub(ii,jj,kk)*p(9)*dey(9)
      ! if (debug) then
      !   write (iounit(2),*) 9, i,j,k, psisub(ii,jj,kk)*p(9)*dex(9),psisub(ii,jj,kk)*p(9)*dey(9),0
      ! endif
    endif
  endif

  ii=i+ex(10)
  jj=j+ey(10)
  kk=k
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i+ex(9)
  jo=j+ey(9)
  ko=k
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      mygradx(i,j,k)=mygradx(i,j,k)+psisub(ii,jj,kk)*p(10)*dex(10)
      mygrady(i,j,k)=mygrady(i,j,k)+psisub(ii,jj,kk)*p(10)*dey(10)
      ! if (debug) then
      !   write (iounit(2),*) 10, i,j,k, psisub(ii,jj,kk)*p(10)*dex(10),psisub(ii,jj,kk)*p(10)*dey(10),0
      ! endif
    endif
  endif

  ii=i+ex(11)
  jj=j
  kk=k+ez(11)
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i+ex(12)
  jo=j
  ko=k+ez(12)
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      mygradx(i,j,k)=mygradx(i,j,k)+psisub(ii,jj,kk)*p(11)*dex(11)
      mygradz(i,j,k)=mygradz(i,j,k)+psisub(ii,jj,kk)*p(11)*dez(11)
      ! if (debug) then
      !   write (iounit(2),*) 11, i,j,k, psisub(ii,jj,kk)*p(11)*dex(11),0,psisub(ii,jj,kk)*p(11)*dez(11)
      ! endif
    endif
  endif

  ii=i+ex(12)
  jj=j
  kk=k+ez(12)
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i+ex(11)
  jo=j
  ko=k+ez(11)
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      mygradx(i,j,k)=mygradx(i,j,k)+psisub(ii,jj,kk)*p(12)*dex(12)
      mygradz(i,j,k)=mygradz(i,j,k)+psisub(ii,jj,kk)*p(12)*dez(12)
      ! if (debug) then
      !   write (iounit(2),*) 12, i,j,k, psisub(ii,jj,kk)*p(12)*dex(12),0,psisub(ii,jj,kk)*p(12)*dez(12)
      ! endif
    endif
  endif

  ii=i+ex(13)
  jj=j
  kk=k+ez(13)
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i+ex(14)
  jo=j
  ko=k+ez(14)
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      mygradx(i,j,k)=mygradx(i,j,k)+psisub(ii,jj,kk)*p(13)*dex(13)
      mygradz(i,j,k)=mygradz(i,j,k)+psisub(ii,jj,kk)*p(13)*dez(13)
      ! if (debug) then
      !   write (iounit(2),*) 13, i,j,k, psisub(ii,jj,kk)*p(13)*dex(13),0,psisub(ii,jj,kk)*p(13)*dez(13)
      ! endif
    endif
  endif

  ii=i+ex(14)
  jj=j
  kk=k+ez(14)
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i+ex(13)
  jo=j
  ko=k+ez(13)
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      mygradx(i,j,k)=mygradx(i,j,k)+psisub(ii,jj,kk)*p(14)*dex(14)
      mygradz(i,j,k)=mygradz(i,j,k)+psisub(ii,jj,kk)*p(14)*dez(14)
      ! if (debug) then
      !   write (iounit(2),*) 14, i,j,k, psisub(ii,jj,kk)*p(14)*dex(14),0,psisub(ii,jj,kk)*p(14)*dez(14)
      ! endif
    endif
  endif

  ii=i
  jj=j+ey(15)
  kk=k+ez(15)
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i
  jo=j+ey(16)
  ko=k+ez(16)
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      mygrady(i,j,k)=mygrady(i,j,k)+psisub(ii,jj,kk)*p(15)*dey(15)
      mygradz(i,j,k)=mygradz(i,j,k)+psisub(ii,jj,kk)*p(15)*dez(15)
      ! if (debug) then
      !   write (iounit(2),*) 15, i,j,k, 0,psisub(ii,jj,kk)*p(15)*dey(15),psisub(ii,jj,kk)*p(15)*dez(15)
      ! endif
    endif
  endif

  ii=i
  jj=j+ey(16)
  kk=k+ez(16)
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i
  jo=j+ey(15)
  ko=k+ez(15)
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      mygrady(i,j,k)=mygrady(i,j,k)+psisub(ii,jj,kk)*p(16)*dey(16)
      mygradz(i,j,k)=mygradz(i,j,k)+psisub(ii,jj,kk)*p(16)*dez(16)
      ! if (debug) then
      !   write (iounit(2),*) 16, i,j,k, 0,psisub(ii,jj,kk)*p(16)*dey(16),psisub(ii,jj,kk)*p(16)*dez(16)
      ! endif
    endif
  endif

  ii=i
  jj=j+ey(17)
  kk=k+ez(17)
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i
  jo=j+ey(18)
  ko=k+ez(18)
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      mygrady(i,j,k)=mygrady(i,j,k)+psisub(ii,jj,kk)*p(17)*dey(17)
      mygradz(i,j,k)=mygradz(i,j,k)+psisub(ii,jj,kk)*p(17)*dez(17)
      ! if (debug) then
      !   write (iounit(2),*) 17, i,j,k, 0,psisub(ii,jj,kk)*p(17)*dey(17),psisub(ii,jj,kk)*p(17)*dez(17)
      ! endif
    endif
  endif

  ii=i
  jj=j+ey(18)
  kk=k+ez(18)
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i
  jo=j+ey(17)
  ko=k+ez(17)
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      mygrady(i,j,k)=mygrady(i,j,k)+psisub(ii,jj,kk)*p(18)*dey(18)
      mygradz(i,j,k)=mygradz(i,j,k)+psisub(ii,jj,kk)*p(18)*dez(18)
      ! if (debug) then
      !   write (iounit(2),*) 18, i,j,k, 0,psisub(ii,jj,kk)*p(18)*dey(18),psisub(ii,jj,kk)*p(18)*dez(18)
      ! endif
    endif
  endif

  return

 end subroutine compute_grad_on_particle_old

 
 subroutine compute_fluid_force_sc_old

!***********************************************************************
!
!     LBsoft subroutine for computing the Shan Chen pair interaction
!     forces
!
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!
!***********************************************************************

  implicit none

  integer :: i,j,k

  !red fluid
  call compute_grad_on_lattice(psiR,gradpsixR,gradpsiyR,gradpsizR)
  !blue fluid
  call compute_grad_on_lattice(psiB,gradpsixB,gradpsiyB,gradpsizB)



  !red fluid
  forall(i=minx:maxx,j=miny:maxy,k=minz:maxz,isfluid(i,j,k)==1)
    fuR(i,j,k) = fuR(i,j,k) - pair_SC*psiR(i,j,k)*gradpsixB(i,j,k)
    fvR(i,j,k) = fvR(i,j,k) - pair_SC*psiR(i,j,k)*gradpsiyB(i,j,k)
    fwR(i,j,k) = fwR(i,j,k) - pair_SC*psiR(i,j,k)*gradpsizB(i,j,k)
  end forall
  !blue fluid
  forall(i=minx:maxx,j=miny:maxy,k=minz:maxz,isfluid(i,j,k)==1)
    fuB(i,j,k) = fuB(i,j,k) - pair_SC*psiB(i,j,k)*gradpsixR(i,j,k)
    fvB(i,j,k) = fvB(i,j,k) - pair_SC*psiB(i,j,k)*gradpsiyR(i,j,k)
    fwB(i,j,k) = fwB(i,j,k) - pair_SC*psiB(i,j,k)*gradpsizR(i,j,k)
  end forall

  return

 end subroutine compute_fluid_force_sc_old

 subroutine compute_grad_on_lattice_old(myarr,mygradx,mygrady,mygradz)
 
!***********************************************************************
!     
!     LBsoft subroutine for computing the gradient of a scalar
!     over the lattice
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC), allocatable, dimension(:,:,:), intent(in) :: myarr
  real(kind=PRC), allocatable, dimension(:,:,:), intent(inout) :: &
   mygradx,mygrady,mygradz
   
  integer :: i,j,k,l
  
#if LATTICE==319
  !tolti gli zeri su dex dey dez con pazienza
  forall(i=minx:maxx,j=miny:maxy,k=minz:maxz,isfluid(i,j,k)==1)
    mygradx(i,j,k)= &
     myarr(i+ex(1),j+ey(1),k+ez(1))*p(1)*dex(1)+ & !01
     myarr(i+ex(2),j+ey(2),k+ez(2))*p(2)*dex(2)+ & !02
     myarr(i+ex(7),j+ey(7),k+ez(7))*p(7)*dex(7)+ & !07
     myarr(i+ex(8),j+ey(8),k+ez(8))*p(8)*dex(8)+ & !08
     myarr(i+ex(9),j+ey(9),k+ez(9))*p(9)*dex(9)+ & !09
     myarr(i+ex(10),j+ey(10),k+ez(10))*p(10)*dex(10)+ & !10
     myarr(i+ex(11),j+ey(11),k+ez(11))*p(11)*dex(11)+ & !11
     myarr(i+ex(12),j+ey(12),k+ez(12))*p(12)*dex(12)+ & !12
     myarr(i+ex(13),j+ey(13),k+ez(13))*p(13)*dex(13)+ & !13
     myarr(i+ex(14),j+ey(14),k+ez(14))*p(14)*dex(14)    !14
  end forall
  
  forall(i=minx:maxx,j=miny:maxy,k=minz:maxz,isfluid(i,j,k)==1)
    mygrady(i,j,k)= &
     myarr(i+ex(3),j+ey(3),k+ez(3))*p(3)*dey(3)+ & !03
     myarr(i+ex(4),j+ey(4),k+ez(4))*p(4)*dey(4)+ & !04
     myarr(i+ex(7),j+ey(7),k+ez(7))*p(7)*dey(7)+ & !07
     myarr(i+ex(8),j+ey(8),k+ez(8))*p(8)*dey(8)+ & !08
     myarr(i+ex(9),j+ey(9),k+ez(9))*p(9)*dey(9)+ & !09
     myarr(i+ex(10),j+ey(10),k+ez(10))*p(10)*dey(10)+ & !10
     myarr(i+ex(15),j+ey(15),k+ez(15))*p(15)*dey(15)+ & !15
     myarr(i+ex(16),j+ey(16),k+ez(16))*p(16)*dey(16)+ & !16
     myarr(i+ex(17),j+ey(17),k+ez(17))*p(17)*dey(17)+ & !17
     myarr(i+ex(18),j+ey(18),k+ez(18))*p(18)*dey(18)    !18
  end forall
  
  forall(i=minx:maxx,j=miny:maxy,k=minz:maxz,isfluid(i,j,k)==1)
    mygradz(i,j,k)= &
     myarr(i+ex(5),j+ey(5),k+ez(5))*p(5)*dez(5)+ & !05
     myarr(i+ex(6),j+ey(6),k+ez(6))*p(6)*dez(6)+ & !06
     myarr(i+ex(11),j+ey(11),k+ez(11))*p(11)*dez(11)+ & !11
     myarr(i+ex(12),j+ey(12),k+ez(12))*p(12)*dez(12)+ & !12
     myarr(i+ex(13),j+ey(13),k+ez(13))*p(13)*dez(13)+ & !13
     myarr(i+ex(14),j+ey(14),k+ez(14))*p(14)*dez(14)+ & !14
     myarr(i+ex(15),j+ey(15),k+ez(15))*p(15)*dez(15)+ & !15
     myarr(i+ex(16),j+ey(16),k+ez(16))*p(16)*dez(16)+ & !16
     myarr(i+ex(17),j+ey(17),k+ez(17))*p(17)*dez(17)+ & !17
     myarr(i+ex(18),j+ey(18),k+ez(18))*p(18)*dez(18)    !18
  end forall
#else
  !occhio gli zeri su dex dey dez andrebbero tolti con pazienza
  forall(i=minx:maxx,j=miny:maxy,k=minz:maxz,isfluid(i,j,k)==1)
    mygradx(i,j,k)= &
     myarr(i+ex(1),j+ey(1),k+ez(1))*p(1)*dex(1)+ & !01
     myarr(i+ex(2),j+ey(2),k+ez(2))*p(2)*dex(2)+ & !02
     myarr(i+ex(3),j+ey(3),k+ez(3))*p(3)*dex(3)+ & !03
     myarr(i+ex(4),j+ey(4),k+ez(4))*p(4)*dex(4)+ & !04
     myarr(i+ex(5),j+ey(5),k+ez(5))*p(5)*dex(5)+ & !05
     myarr(i+ex(6),j+ey(6),k+ez(6))*p(6)*dex(6)+ & !06
     myarr(i+ex(7),j+ey(7),k+ez(7))*p(7)*dex(7)+ & !07
     myarr(i+ex(8),j+ey(8),k+ez(8))*p(8)*dex(8)+ & !08
     myarr(i+ex(9),j+ey(9),k+ez(9))*p(9)*dex(9)+ & !09
     myarr(i+ex(10),j+ey(10),k+ez(10))*p(10)*dex(10)+ & !10
     myarr(i+ex(11),j+ey(11),k+ez(11))*p(11)*dex(11)+ & !11
     myarr(i+ex(12),j+ey(12),k+ez(12))*p(12)*dex(12)+ & !12
     myarr(i+ex(13),j+ey(13),k+ez(13))*p(13)*dex(13)+ & !13
     myarr(i+ex(14),j+ey(14),k+ez(14))*p(14)*dex(14)+ & !14
     myarr(i+ex(15),j+ey(15),k+ez(15))*p(15)*dex(15)+ & !15
     myarr(i+ex(16),j+ey(16),k+ez(16))*p(16)*dex(16)+ & !16
     myarr(i+ex(17),j+ey(17),k+ez(17))*p(17)*dex(17)+ & !17
     myarr(i+ex(18),j+ey(18),k+ez(18))*p(18)*dex(18)    !18
  end forall
  
  forall(i=minx:maxx,j=miny:maxy,k=minz:maxz,isfluid(i,j,k)==1)
    mygrady(i,j,k)= &
     myarr(i+ex(1),j+ey(1),k+ez(1))*p(1)*dey(1)+ & !01
     myarr(i+ex(2),j+ey(2),k+ez(2))*p(2)*dey(2)+ & !02
     myarr(i+ex(3),j+ey(3),k+ez(3))*p(3)*dey(3)+ & !03
     myarr(i+ex(4),j+ey(4),k+ez(4))*p(4)*dey(4)+ & !04
     myarr(i+ex(5),j+ey(5),k+ez(5))*p(5)*dey(5)+ & !05
     myarr(i+ex(6),j+ey(6),k+ez(6))*p(6)*dey(6)+ & !06
     myarr(i+ex(7),j+ey(7),k+ez(7))*p(7)*dey(7)+ & !07
     myarr(i+ex(8),j+ey(8),k+ez(8))*p(8)*dey(8)+ & !08
     myarr(i+ex(9),j+ey(9),k+ez(9))*p(9)*dey(9)+ & !09
     myarr(i+ex(10),j+ey(10),k+ez(10))*p(10)*dey(10)+ & !10
     myarr(i+ex(11),j+ey(11),k+ez(11))*p(11)*dey(11)+ & !11
     myarr(i+ex(12),j+ey(12),k+ez(12))*p(12)*dey(12)+ & !12
     myarr(i+ex(13),j+ey(13),k+ez(13))*p(13)*dey(13)+ & !13
     myarr(i+ex(14),j+ey(14),k+ez(14))*p(14)*dey(14)+ & !14
     myarr(i+ex(15),j+ey(15),k+ez(15))*p(15)*dey(15)+ & !15
     myarr(i+ex(16),j+ey(16),k+ez(16))*p(16)*dey(16)+ & !16
     myarr(i+ex(17),j+ey(17),k+ez(17))*p(17)*dey(17)+ & !17
     myarr(i+ex(18),j+ey(18),k+ez(18))*p(18)*dey(18)    !18
  end forall
  
  forall(i=minx:maxx,j=miny:maxy,k=minz:maxz,isfluid(i,j,k)==1)
    mygradz(i,j,k)= &
     myarr(i+ex(1),j+ey(1),k+ez(1))*p(1)*dez(1)+ & !01
     myarr(i+ex(2),j+ey(2),k+ez(2))*p(2)*dez(2)+ & !02
     myarr(i+ex(3),j+ey(3),k+ez(3))*p(3)*dez(3)+ & !03
     myarr(i+ex(4),j+ey(4),k+ez(4))*p(4)*dez(4)+ & !04
     myarr(i+ex(5),j+ey(5),k+ez(5))*p(5)*dez(5)+ & !05
     myarr(i+ex(6),j+ey(6),k+ez(6))*p(6)*dez(6)+ & !06
     myarr(i+ex(7),j+ey(7),k+ez(7))*p(7)*dez(7)+ & !07
     myarr(i+ex(8),j+ey(8),k+ez(8))*p(8)*dez(8)+ & !08
     myarr(i+ex(9),j+ey(9),k+ez(9))*p(9)*dez(9)+ & !09
     myarr(i+ex(10),j+ey(10),k+ez(10))*p(10)*dez(10)+ & !10
     myarr(i+ex(11),j+ey(11),k+ez(11))*p(11)*dez(11)+ & !11
     myarr(i+ex(12),j+ey(12),k+ez(12))*p(12)*dez(12)+ & !12
     myarr(i+ex(13),j+ey(13),k+ez(13))*p(13)*dez(13)+ & !13
     myarr(i+ex(14),j+ey(14),k+ez(14))*p(14)*dez(14)+ & !14
     myarr(i+ex(15),j+ey(15),k+ez(15))*p(15)*dez(15)+ & !15
     myarr(i+ex(16),j+ey(16),k+ez(16))*p(16)*dez(16)+ & !16
     myarr(i+ex(17),j+ey(17),k+ez(17))*p(17)*dez(17)+ & !17
     myarr(i+ex(18),j+ey(18),k+ez(18))*p(18)*dez(18)    !18
  end forall
#endif
  
  return
  
 end subroutine compute_grad_on_lattice_old

 subroutine compute_fluid_force_sc
 
!***********************************************************************
!     
!     LBsoft subroutine for computing the Shan Chen pair interaction
!     forces
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer :: i,j,k
  

  
  !red fluid
  call compute_grad_on_lattice(psiR,gradpsixR,gradpsiyR,gradpsizR)
  !blue fluid
  call compute_grad_on_lattice(psiB,gradpsixB,gradpsiyB,gradpsizB)
  
  !red fluid
  do k=minz,maxz
   do j=miny,maxy
    do i=minx,maxx

        if (isfluid(i,j,k)/=1) cycle
    fuR(i,j,k) = fuR(i,j,k) - pair_SC*psiR(i,j,k)*gradpsixB(i,j,k)
    fvR(i,j,k) = fvR(i,j,k) - pair_SC*psiR(i,j,k)*gradpsiyB(i,j,k)
    fwR(i,j,k) = fwR(i,j,k) - pair_SC*psiR(i,j,k)*gradpsizB(i,j,k)
      enddo
   enddo
  enddo
   
  !blue fluid
  do k=minz,maxz
   do j=miny,maxy
    do i=minx,maxx

        if (isfluid(i,j,k)/=1) cycle
    fuB(i,j,k) = fuB(i,j,k) - pair_SC*psiB(i,j,k)*gradpsixR(i,j,k)
    fvB(i,j,k) = fvB(i,j,k) - pair_SC*psiB(i,j,k)*gradpsiyR(i,j,k)
    fwB(i,j,k) = fwB(i,j,k) - pair_SC*psiB(i,j,k)*gradpsizR(i,j,k)
    enddo
   enddo
  enddo

  return
  
 end subroutine compute_fluid_force_sc
 
 subroutine compute_grad_on_lattice(myarr,mygradx,mygrady,mygradz)
 
!***********************************************************************
!     
!     LBsoft subroutine for computing the gradient of a scalar
!     over the lattice
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC), allocatable, dimension(:,:,:), intent(in) :: myarr
  real(kind=PRC), allocatable, dimension(:,:,:), intent(inout) :: &
   mygradx,mygrady,mygradz
   
  integer :: i,j,k,l
  
#if LATTICE==319
  !tolti gli zeri su dex dey dez con pazienza
  !forall(i=minx:maxx,j=miny:maxy,k=minz:maxz,isfluid(i,j,k)==1)
  do k=minz,maxz
   do j=miny,maxy
    do i=minx,maxx

        if (isfluid(i,j,k)/=1) cycle
    mygradx(i,j,k)= &
     myarr(i+ex(1),j+ey(1),k+ez(1))*p(1)*dex(1)+ & !01
     myarr(i+ex(2),j+ey(2),k+ez(2))*p(2)*dex(2)+ & !02
     myarr(i+ex(7),j+ey(7),k+ez(7))*p(7)*dex(7)+ & !07
     myarr(i+ex(8),j+ey(8),k+ez(8))*p(8)*dex(8)+ & !08
     myarr(i+ex(9),j+ey(9),k+ez(9))*p(9)*dex(9)+ & !09
     myarr(i+ex(10),j+ey(10),k+ez(10))*p(10)*dex(10)+ & !10
     myarr(i+ex(11),j+ey(11),k+ez(11))*p(11)*dex(11)+ & !11
     myarr(i+ex(12),j+ey(12),k+ez(12))*p(12)*dex(12)+ & !12
     myarr(i+ex(13),j+ey(13),k+ez(13))*p(13)*dex(13)+ & !13
     myarr(i+ex(14),j+ey(14),k+ez(14))*p(14)*dex(14)    !14
    enddo
   enddo
  enddo
  !end forall
  
  !forall(i=minx:maxx,j=miny:maxy,k=minz:maxz,isfluid(i,j,k)==1)
  do k=minz,maxz
   do j=miny,maxy
    do i=minx,maxx

        if (isfluid(i,j,k)/=1) cycle

    mygrady(i,j,k)= &
     myarr(i+ex(3),j+ey(3),k+ez(3))*p(3)*dey(3)+ & !03
     myarr(i+ex(4),j+ey(4),k+ez(4))*p(4)*dey(4)+ & !04
     myarr(i+ex(7),j+ey(7),k+ez(7))*p(7)*dey(7)+ & !07
     myarr(i+ex(8),j+ey(8),k+ez(8))*p(8)*dey(8)+ & !08
     myarr(i+ex(9),j+ey(9),k+ez(9))*p(9)*dey(9)+ & !09
     myarr(i+ex(10),j+ey(10),k+ez(10))*p(10)*dey(10)+ & !10
     myarr(i+ex(15),j+ey(15),k+ez(15))*p(15)*dey(15)+ & !15
     myarr(i+ex(16),j+ey(16),k+ez(16))*p(16)*dey(16)+ & !16
     myarr(i+ex(17),j+ey(17),k+ez(17))*p(17)*dey(17)+ & !17
     myarr(i+ex(18),j+ey(18),k+ez(18))*p(18)*dey(18)    !18
    enddo
   enddo
  enddo
  !end forall
  
  !forall(i=minx:maxx,j=miny:maxy,k=minz:maxz,isfluid(i,j,k)==1)
  do k=minz,maxz
   do j=miny,maxy
    do i=minx,maxx

        if (isfluid(i,j,k)/=1) cycle

    mygradz(i,j,k)= &
     myarr(i+ex(5),j+ey(5),k+ez(5))*p(5)*dez(5)+ & !05
     myarr(i+ex(6),j+ey(6),k+ez(6))*p(6)*dez(6)+ & !06
     myarr(i+ex(11),j+ey(11),k+ez(11))*p(11)*dez(11)+ & !11
     myarr(i+ex(12),j+ey(12),k+ez(12))*p(12)*dez(12)+ & !12
     myarr(i+ex(13),j+ey(13),k+ez(13))*p(13)*dez(13)+ & !13
     myarr(i+ex(14),j+ey(14),k+ez(14))*p(14)*dez(14)+ & !14
     myarr(i+ex(15),j+ey(15),k+ez(15))*p(15)*dez(15)+ & !15
     myarr(i+ex(16),j+ey(16),k+ez(16))*p(16)*dez(16)+ & !16
     myarr(i+ex(17),j+ey(17),k+ez(17))*p(17)*dez(17)+ & !17
     myarr(i+ex(18),j+ey(18),k+ez(18))*p(18)*dez(18)    !18
    enddo
   enddo
  enddo
  !end forall
#else
  !occhio gli zeri su dex dey dez andrebbero tolti con pazienza
  forall(i=minx:maxx,j=miny:maxy,k=minz:maxz,isfluid(i,j,k)==1)
    mygradx(i,j,k)= &
     myarr(i+ex(1),j+ey(1),k+ez(1))*p(1)*dex(1)+ & !01
     myarr(i+ex(2),j+ey(2),k+ez(2))*p(2)*dex(2)+ & !02
     myarr(i+ex(3),j+ey(3),k+ez(3))*p(3)*dex(3)+ & !03
     myarr(i+ex(4),j+ey(4),k+ez(4))*p(4)*dex(4)+ & !04
     myarr(i+ex(5),j+ey(5),k+ez(5))*p(5)*dex(5)+ & !05
     myarr(i+ex(6),j+ey(6),k+ez(6))*p(6)*dex(6)+ & !06
     myarr(i+ex(7),j+ey(7),k+ez(7))*p(7)*dex(7)+ & !07
     myarr(i+ex(8),j+ey(8),k+ez(8))*p(8)*dex(8)+ & !08
     myarr(i+ex(9),j+ey(9),k+ez(9))*p(9)*dex(9)+ & !09
     myarr(i+ex(10),j+ey(10),k+ez(10))*p(10)*dex(10)+ & !10
     myarr(i+ex(11),j+ey(11),k+ez(11))*p(11)*dex(11)+ & !11
     myarr(i+ex(12),j+ey(12),k+ez(12))*p(12)*dex(12)+ & !12
     myarr(i+ex(13),j+ey(13),k+ez(13))*p(13)*dex(13)+ & !13
     myarr(i+ex(14),j+ey(14),k+ez(14))*p(14)*dex(14)+ & !14
     myarr(i+ex(15),j+ey(15),k+ez(15))*p(15)*dex(15)+ & !15
     myarr(i+ex(16),j+ey(16),k+ez(16))*p(16)*dex(16)+ & !16
     myarr(i+ex(17),j+ey(17),k+ez(17))*p(17)*dex(17)+ & !17
     myarr(i+ex(18),j+ey(18),k+ez(18))*p(18)*dex(18)    !18
  end forall
  
  forall(i=minx:maxx,j=miny:maxy,k=minz:maxz,isfluid(i,j,k)==1)
    mygrady(i,j,k)= &
     myarr(i+ex(1),j+ey(1),k+ez(1))*p(1)*dey(1)+ & !01
     myarr(i+ex(2),j+ey(2),k+ez(2))*p(2)*dey(2)+ & !02
     myarr(i+ex(3),j+ey(3),k+ez(3))*p(3)*dey(3)+ & !03
     myarr(i+ex(4),j+ey(4),k+ez(4))*p(4)*dey(4)+ & !04
     myarr(i+ex(5),j+ey(5),k+ez(5))*p(5)*dey(5)+ & !05
     myarr(i+ex(6),j+ey(6),k+ez(6))*p(6)*dey(6)+ & !06
     myarr(i+ex(7),j+ey(7),k+ez(7))*p(7)*dey(7)+ & !07
     myarr(i+ex(8),j+ey(8),k+ez(8))*p(8)*dey(8)+ & !08
     myarr(i+ex(9),j+ey(9),k+ez(9))*p(9)*dey(9)+ & !09
     myarr(i+ex(10),j+ey(10),k+ez(10))*p(10)*dey(10)+ & !10
     myarr(i+ex(11),j+ey(11),k+ez(11))*p(11)*dey(11)+ & !11
     myarr(i+ex(12),j+ey(12),k+ez(12))*p(12)*dey(12)+ & !12
     myarr(i+ex(13),j+ey(13),k+ez(13))*p(13)*dey(13)+ & !13
     myarr(i+ex(14),j+ey(14),k+ez(14))*p(14)*dey(14)+ & !14
     myarr(i+ex(15),j+ey(15),k+ez(15))*p(15)*dey(15)+ & !15
     myarr(i+ex(16),j+ey(16),k+ez(16))*p(16)*dey(16)+ & !16
     myarr(i+ex(17),j+ey(17),k+ez(17))*p(17)*dey(17)+ & !17
     myarr(i+ex(18),j+ey(18),k+ez(18))*p(18)*dey(18)    !18
  end forall
  
  forall(i=minx:maxx,j=miny:maxy,k=minz:maxz,isfluid(i,j,k)==1)
    mygradz(i,j,k)= &
     myarr(i+ex(1),j+ey(1),k+ez(1))*p(1)*dez(1)+ & !01
     myarr(i+ex(2),j+ey(2),k+ez(2))*p(2)*dez(2)+ & !02
     myarr(i+ex(3),j+ey(3),k+ez(3))*p(3)*dez(3)+ & !03
     myarr(i+ex(4),j+ey(4),k+ez(4))*p(4)*dez(4)+ & !04
     myarr(i+ex(5),j+ey(5),k+ez(5))*p(5)*dez(5)+ & !05
     myarr(i+ex(6),j+ey(6),k+ez(6))*p(6)*dez(6)+ & !06
     myarr(i+ex(7),j+ey(7),k+ez(7))*p(7)*dez(7)+ & !07
     myarr(i+ex(8),j+ey(8),k+ez(8))*p(8)*dez(8)+ & !08
     myarr(i+ex(9),j+ey(9),k+ez(9))*p(9)*dez(9)+ & !09
     myarr(i+ex(10),j+ey(10),k+ez(10))*p(10)*dez(10)+ & !10
     myarr(i+ex(11),j+ey(11),k+ez(11))*p(11)*dez(11)+ & !11
     myarr(i+ex(12),j+ey(12),k+ez(12))*p(12)*dez(12)+ & !12
     myarr(i+ex(13),j+ey(13),k+ez(13))*p(13)*dez(13)+ & !13
     myarr(i+ex(14),j+ey(14),k+ez(14))*p(14)*dez(14)+ & !14
     myarr(i+ex(15),j+ey(15),k+ez(15))*p(15)*dez(15)+ & !15
     myarr(i+ex(16),j+ey(16),k+ez(16))*p(16)*dez(16)+ & !16
     myarr(i+ex(17),j+ey(17),k+ez(17))*p(17)*dez(17)+ & !17
     myarr(i+ex(18),j+ey(18),k+ez(18))*p(18)*dez(18)    !18
  end forall
#endif
  
  return
  
 end subroutine compute_grad_on_lattice

!*********START PART TO MANAGE THE INTERACTION WITH PARTICLES***********
 subroutine helper_init_particle_2_isfluid(isub,jsub,ksub, nspheres, spherelists, isInternal)
 implicit none
 integer, intent(in) :: isub,jsub,ksub,nspheres
 integer, allocatable, dimension(:,:), intent(in) :: spherelists
 logical, intent(in) :: isInternal
 integer :: i,j,k,l
 integer :: imin,imax,jmin,jmax,kmin,kmax

 imin=minx-nbuff
 imax=maxx+nbuff
 jmin=miny-nbuff
 jmax=maxy+nbuff
 kmin=minz-nbuff
 kmax=maxz+nbuff


 do l=1,nspheres
  i=isub+spherelists(1,l)
  j=jsub+spherelists(2,l)
  k=ksub+spherelists(3,l)
  !apply periodic conditions if necessary
  i=pimage(ixpbc,i,nx)
  j=pimage(iypbc,j,ny)
  k=pimage(izpbc,k,nz)

  CYCLE_OUT_INTERVAL(i, imin, imax)
  CYCLE_OUT_INTERVAL(j, jmin, jmax)
  CYCLE_OUT_INTERVAL(k, kmin, kmax)

  if (isInternal) then
    if(isfluid(i,j,k)/=4)isfluid(i,j,k)=2
  else
    isfluid(i,j,k)=4
  endif

  rhoR(i,j,k)=MINDENS
  if (.not. lsingle_fluid) then
    rhoB(i,j,k)=MINDENS
  endif
  u(i,j,k)=ZERO
  v(i,j,k)=ZERO
  w(i,j,k)=ZERO
 enddo
 end subroutine helper_init_particle_2_isfluid


 subroutine init_particle_2_isfluid(myi,isub,jsub,ksub,nspheres, &
  spherelists,spheredists,nspheredeads,spherelistdeads, isInternal)
  
!***********************************************************************
!     
!     LBsoft subroutine to initialize isfluid and hydrodynamic
!     variables according to the particle presence
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: myi,isub,jsub,ksub,nspheres,nspheredeads
  integer, allocatable, dimension(:,:), intent(in) :: spherelists
  integer, allocatable, dimension(:,:), intent(in) :: spherelistdeads
  real(kind=PRC), allocatable, dimension(:), intent(in) :: spheredists
  logical, intent(in) :: isInternal
  

  call helper_init_particle_2_isfluid(isub,jsub,ksub,nspheres, spherelists, .true.)
  call helper_init_particle_2_isfluid(isub,jsub,ksub,nspheredeads, spherelistdeads, .false.)

  if (isInternal) then
    isfluid(isub,jsub,ksub)=5
    rhoR(isub,jsub,ksub)=MINDENS
    u(isub,jsub,ksub)=ZERO
    v(isub,jsub,ksub)=ZERO
    w(isub,jsub,ksub)=ZERO

    if(.not. lsingle_fluid)then
     rhoB(isub,jsub,ksub)=MINDENS
    endif
  endif
 
 end subroutine init_particle_2_isfluid

 subroutine particle_bounce_back(debug, nstep,iatm,lown,lrotate,isub,jsub,ksub,nspheres, &
  spherelists,spheredists,rdimx,rdimy,rdimz,xx,yy,zz,vx,vy,vz,&
  fx,fy,fz, A, ox,oy,oz,tx,ty,tz)
  
!***********************************************************************
!     
!     LBsoft subroutine to initialize isfluid and hydrodynamic
!     variables according to the particle presence
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
  
  implicit none
  
  logical, intent(in) :: lown,lrotate, debug
  integer, intent(in) :: nstep,iatm,isub,jsub,ksub,nspheres
  integer, allocatable, dimension(:,:), intent(in) :: spherelists
  real(kind=PRC), allocatable, dimension(:), intent(in) :: spheredists
  real(kind=PRC), intent(in) :: rdimx,rdimy,rdimz
  real(kind=PRC), intent(in) :: xx,yy,zz
  real(kind=PRC), intent(in) :: vx,vy,vz
  real(kind=PRC), intent(in), optional :: ox,oy,oz

#ifdef QUAD_FORCEINT
  real(kind=PRC*2), intent(inout) :: fx,fy,fz
  real(kind=PRC*2), intent(inout), optional :: tx,ty,tz
#else
  real(kind=PRC),   intent(inout) :: fx,fy,fz
  real(kind=PRC),   intent(inout), optional :: tx,ty,tz
#endif
  
  integer :: i,j,k,l,ii,jj,kk
  integer, save :: imin,imax,jmin,jmax,kmin,kmax
  logical, save :: lfirst=.true.
  real(kind=PRC) :: vxt,vyt,vzt,modr,ftx,fty,ftz
  real(kind=PRC), dimension(3) :: rtemp,otemp,ftemp
  
#ifdef QUAD_FORCEINT
  real(kind=PRC*2), allocatable, dimension(:,:,:) :: A
#else
  real(kind=PRC), allocatable, dimension(:,:,:) :: A
#endif
  character(len=120) :: mynamefile




  if(lfirst)then
    lfirst=.false.
    imin=minx-1
    imax=maxx+1
    jmin=miny-1
    jmax=maxy+1
    kmin=minz-1
    kmax=maxz+1
  endif
  
  if(lrotate)then
    otemp(1)=ox
    otemp(2)=oy
    otemp(3)=oz
  endif


  ! if (debug) then
  !  do l=1,3
  !   mynamefile=repeat(' ',120)
  !   mynamefile="forceInt.atom"//write_fmtnumb(iatm)//'.step'//write_fmtnumb(l)
  !   iounit(l) = 113 + l
  !   call OpenLogFile(nstep, mynamefile, iounit(l))
  !  enddo
  ! endif

  
  if(lsingle_fluid)then
    do l=1,nspheres
      i=isub+spherelists(1,l)
      j=jsub+spherelists(2,l)
      k=ksub+spherelists(3,l)
      ii=i
      jj=j
      kk=k

      !apply periodic conditions if necessary
      i=pimage(ixpbc,i,nx)
      j=pimage(iypbc,j,ny)
      k=pimage(izpbc,k,nz)

      if(lrotate)then
        rtemp(1)=real(ii,kind=PRC)-xx
        rtemp(2)=real(jj,kind=PRC)-yy
        rtemp(3)=real(kk,kind=PRC)-zz
        modr=modulvec(rtemp)
        rtemp(1:3)=rdimx/modr*rtemp(1:3)
      endif
      if(i>=minx .and. i<=maxx .and. j>=miny .and. j<=maxy .and. k>=minz .and. k<=maxz)then
        call node_to_particle_bounce_back_bc2(lrotate,nstep,i,j,k,rtemp, &
           otemp,vx,vy,vz,fx,fy,fz,tx,ty,tz,rhoR,fluidR, debug,iatm, A, l)
        ! if (debug) then
        !  if(lrotate) then
        !    write(iounit(1),*) __FILE__,9249, "i,j,k=", i,j,k, "xx,..", xx,yy,zz, &
        !             "vx,..",vx,vy,vz, rtemp,otemp
        !  else
        !    write(iounit(1),*) __FILE__,9249, "i,j,k=", i,j,k, "xx,..", xx,yy,zz, &
        !             "vx,..",vx,vy,vz
        !  endif
        ! endif
      endif

      if(i>=minx .and. i<=maxx .and. j>=miny .and. j<=maxy .and. k>=minz .and. k<=maxz)then
        call particle_to_node_bounce_back_bc2_phase1(lrotate,nstep,i,j,k,fluidR, debug,iatm)
      endif

      !the fluid bounce back is local so I have to do it
      if(i/=ii.or.j/=jj.or.k/=kk)then
       if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. kk>=minz .and. kk<=maxz)then
          call particle_to_node_bounce_back_bc2_phase1(lrotate,nstep,ii,jj,kk,fluidR, debug,iatm)
        endif
      endif

    enddo
  else
    do l=1,nspheres
      i=isub+spherelists(1,l)
      j=jsub+spherelists(2,l)
      k=ksub+spherelists(3,l)
      ii=i
      jj=j
      kk=k
      !apply periodic conditions if necessary
      i=pimage(ixpbc,i,nx)
      j=pimage(iypbc,j,ny)
      k=pimage(izpbc,k,nz)

      if(lrotate)then
        rtemp(1)=real(ii,kind=PRC)-xx
        rtemp(2)=real(jj,kind=PRC)-yy
        rtemp(3)=real(kk,kind=PRC)-zz
        modr=modulvec(rtemp)
        rtemp(1:3)=rdimx/modr*rtemp(1:3)
      endif

      if(i>=minx .and. i<=maxx .and. j>=miny .and. j<=maxy .and. k>=minz .and. k<=maxz)then
        call node_to_particle_bounce_back_bc2(lrotate,nstep,i,j,k,rtemp, &
           otemp,vx,vy,vz,fx,fy,fz,tx,ty,tz,rhoR,fluidR, debug,iatm, A, l)
        call node_to_particle_bounce_back_bc2(lrotate,nstep,i,j,k,rtemp, &
           otemp,vx,vy,vz,fx,fy,fz,tx,ty,tz,rhoB,fluidB, debug,iatm, A, l)

        ! if (debug) then
        !  if(lrotate) then
        !    write(iounit(1),*) __FILE__,9249, "i,j,k=", i,j,k, "xx,..", xx,yy,zz, &
        !             "vx,..",vx,vy,vz, rtemp,otemp
        !  else
        !    write(iounit(1),*) __FILE__,9249, "i,j,k=", i,j,k, "xx,..", xx,yy,zz, &
        !             "vx,..",vx,vy,vz
        !  endif
        ! endif
      endif

      if(i>=minx .and. i<=maxx .and. j>=miny .and. j<=maxy .and. k>=minz .and. k<=maxz)then
        call particle_to_node_bounce_back_bc2_phase1(lrotate,nstep,i,j,k, &
         fluidR, debug,iatm)
        call particle_to_node_bounce_back_bc2_phase1(lrotate,nstep,i,j,k, &
         fluidB, debug,iatm)
      endif

      !the fluid bounce back is local so I have to do it
      if(i/=ii.or.j/=jj.or.k/=kk)then
       if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. kk>=minz .and. kk<=maxz)then
          call particle_to_node_bounce_back_bc2_phase1(lrotate,nstep,ii,jj,kk, &
           fluidR, debug,iatm)
          call particle_to_node_bounce_back_bc2_phase1(lrotate,nstep,ii,jj,kk, &
           fluidB, debug,iatm)
        endif
      endif
      
    enddo
  endif

  ! if (debug) then
  !  do l=1,3
  !   close(iounit(l))
  !  enddo
  ! endif
 end subroutine particle_bounce_back

 subroutine particle_bounce_back_phase2(debug, nstep,iatm,lown,lrotate,isub,jsub,ksub,nspheres, &
  spherelists,spheredists,rdimx,rdimy,rdimz,xx,yy,zz,vx,vy,vz,&
  fx,fy,fz, A, ox,oy,oz,tx,ty,tz)

!***********************************************************************
!
!     LBsoft subroutine to initialize isfluid and hydrodynamic
!     variables according to the particle presence
!
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification November 2018
!
!***********************************************************************

  implicit none

  logical, intent(in) :: lown,lrotate, debug
  integer, intent(in) :: nstep,iatm,isub,jsub,ksub,nspheres
  integer, allocatable, dimension(:,:), intent(in) :: spherelists
  real(kind=PRC), allocatable, dimension(:), intent(in) :: spheredists
  real(kind=PRC), intent(in) :: rdimx,rdimy,rdimz
  real(kind=PRC), intent(in) :: xx,yy,zz
  real(kind=PRC), intent(in) :: vx,vy,vz
  real(kind=PRC), intent(in), optional :: ox,oy,oz

#ifdef QUAD_FORCEINT
  real(kind=PRC*2), intent(inout) :: fx,fy,fz
  real(kind=PRC*2), intent(inout), optional :: tx,ty,tz
#else
  real(kind=PRC),   intent(inout) :: fx,fy,fz
  real(kind=PRC),   intent(inout), optional :: tx,ty,tz
#endif

  integer :: i,j,k,l,ii,jj,kk
  integer, save :: imin,imax,jmin,jmax,kmin,kmax
  logical, save :: lfirst=.true.
  real(kind=PRC) :: vxt,vyt,vzt,modr,ftx,fty,ftz
  real(kind=PRC), dimension(3) :: rtemp,otemp,ftemp

#ifdef QUAD_FORCEINT
  real(kind=PRC*2), allocatable, dimension(:,:,:) :: A
#else
  real(kind=PRC), allocatable, dimension(:,:,:) :: A
#endif
  character(len=120) :: mynamefile




  if(lfirst)then
    lfirst=.false.
    imin=minx-1
    imax=maxx+1
    jmin=miny-1
    jmax=maxy+1
    kmin=minz-1
    kmax=maxz+1
  endif

  if(lrotate)then
    otemp(1)=ox
    otemp(2)=oy
    otemp(3)=oz
  endif


  ! if (debug) then
  !  do l=1,3
  !   mynamefile=repeat(' ',120)
  !   mynamefile="forceInt_phase2.atom"//write_fmtnumb(iatm)//'.step'//write_fmtnumb(l)
  !   iounit(l) = 113 + l
  !   call OpenLogFile(nstep, mynamefile, iounit(l))
  !  enddo
  ! endif


  if(lsingle_fluid)then
    do l=1,nspheres
      i=isub+spherelists(1,l)
      j=jsub+spherelists(2,l)
      k=ksub+spherelists(3,l)
      ii=i
      jj=j
      kk=k

      !apply periodic conditions if necessary
      i=pimage(ixpbc,i,nx)
      j=pimage(iypbc,j,ny)
      k=pimage(izpbc,k,nz)

      if(lrotate)then
        rtemp(1)=real(ii,kind=PRC)-xx
        rtemp(2)=real(jj,kind=PRC)-yy
        rtemp(3)=real(kk,kind=PRC)-zz
        modr=modulvec(rtemp)
        rtemp(1:3)=rdimx/modr*rtemp(1:3)
      endif


      if(i>=minx .and. i<=maxx .and. j>=miny .and. j<=maxy .and. k>=minz .and. k<=maxz)then
        call particle_to_node_bounce_back_bc2_phase2(lrotate,nstep,i,j,k,rtemp, &
         otemp,vx,vy,vz,rhoR,fluidR, debug,iatm)
      endif

      !the fluid bounce back is local so I have to do it
      if(i/=ii.or.j/=jj.or.k/=kk)then
       if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. kk>=minz .and. kk<=maxz)then
          call particle_to_node_bounce_back_bc2_phase2(lrotate,nstep,ii,jj,kk,rtemp, &
           otemp,vx,vy,vz,rhoR,fluidR, debug,iatm)
        endif
      endif

    enddo
  else
    do l=1,nspheres
      i=isub+spherelists(1,l)
      j=jsub+spherelists(2,l)
      k=ksub+spherelists(3,l)
      ii=i
      jj=j
      kk=k
      !apply periodic conditions if necessary
      i=pimage(ixpbc,i,nx)
      j=pimage(iypbc,j,ny)
      k=pimage(izpbc,k,nz)

      if(lrotate)then
        rtemp(1)=real(ii,kind=PRC)-xx
        rtemp(2)=real(jj,kind=PRC)-yy
        rtemp(3)=real(kk,kind=PRC)-zz
        modr=modulvec(rtemp)
        rtemp(1:3)=rdimx/modr*rtemp(1:3)
      endif

      if(i>=minx .and. i<=maxx .and. j>=miny .and. j<=maxy .and. k>=minz .and. k<=maxz)then
        call particle_to_node_bounce_back_bc2_phase2(lrotate,nstep,i,j,k,rtemp, &
         otemp,vx,vy,vz,rhoR,fluidR, debug,iatm)
        call particle_to_node_bounce_back_bc2_phase2(lrotate,nstep,i,j,k,rtemp, &
         otemp,vx,vy,vz,rhoB,fluidB, debug,iatm)
      endif

      !the fluid bounce back is local so I have to do it
      if(i/=ii.or.j/=jj.or.k/=kk)then
       if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. kk>=minz .and. kk<=maxz)then
          call particle_to_node_bounce_back_bc2_phase2(lrotate,nstep,ii,jj,kk,rtemp, &
           otemp,vx,vy,vz,rhoR,fluidR, debug,iatm)
          call particle_to_node_bounce_back_bc2_phase2(lrotate,nstep,ii,jj,kk,rtemp, &
           otemp,vx,vy,vz,rhoB,fluidB, debug,iatm)
        endif
      endif

    enddo
  endif

  ! if (debug) then
  !  do l=1,3
  !   close(iounit(l))
  !  enddo
  ! endif
 end subroutine particle_bounce_back_phase2

 subroutine node_to_particle_bounce_back_bc2(lrotate,nstep,i,j,k,rversor, &
   otemp,vxs,vys,vzs,fx,fy,fz,tx,ty,tz,rhosub,fluidsub, debug,iatm, A, l)
 
!***********************************************************************
!     
!     LBsoft subroutine to apply the bounce back of the fluid 
!     on a particle surface node including the pbc
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
 
  implicit none
  
  logical, intent(in) :: lrotate, debug
  integer, intent(in) :: nstep,i,j,k,iatm, l
  real(kind=PRC), intent(in) :: vxs,vys,vzs
  real(kind=PRC), intent(in), dimension(3) :: rversor,otemp

#ifdef QUAD_FORCEINT
  real(kind=PRC*2), intent(inout) :: fx,fy,fz
  real(kind=PRC*2), intent(inout) :: tx,ty,tz
#else
  real(kind=PRC), intent(inout) :: fx,fy,fz
  real(kind=PRC), intent(inout) :: tx,ty,tz
#endif

  real(kind=PRC), allocatable, dimension(:,:,:)  :: rhosub
  real(kind=PRC), allocatable :: fluidsub(:,:,:,:)
  
  real(kind=PRC), parameter :: onesixth=ONE/SIX
  real(kind=PRC), parameter :: pref_bouzidi=TWO/cssq
  real(kind=PRC) :: f2p,vx,vy,vz
  real(kind=PRC), dimension(3) :: rtemp,ftemp, ttemp
  
  integer :: ii,jj,kk,io,jo,ko, iloop,indlow,indhig

#ifdef QUAD_FORCEINT
  real(kind=PRC*2), allocatable, dimension(:,:,:) :: A
#else
  real(kind=PRC), allocatable, dimension(:,:,:) :: A
#endif


  ttemp = ZERO
  
  vx=vxs
  vy=vys
  vz=vzs
  rtemp=rversor
  
  !force on particle fx fy fz
  !eq. 11.2 from page 437 Kruger's book "the lattice boltzmann method"
    
 do iloop = 1, 9
  
  indlow = iloop*2 - 1
  indhig = iloop*2


  ii=i+ex(indlow)
  jj=j+ey(indlow)
  kk=k+ez(indlow)
!  ii=pimage(ixpbc,ii,nx)
!  jj=pimage(iypbc,jj,ny)
!  kk=pimage(izpbc,kk,nz)

  io=i+ex(indhig)
  jo=j+ey(indhig)
  ko=k+ez(indhig)
!  io=pimage(ixpbc,io,nx)
!  jo=pimage(iypbc,jo,ny)
!  ko=pimage(izpbc,ko,nz)

  if(ii>=minx-1 .and. ii<=maxx+1 .and. jj>=miny-1 .and. jj<=maxy+1 .and. &
   kk>=minz-1 .and. kk<=maxz+1)then
    if(io>=minx-1 .and. io<=maxx+1 .and. jo>=miny-1 .and. jo<=maxy+1 .and. &
     ko>=minz-1 .and. ko<=maxz+1)then

  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
      if(lrotate)then
        rtemp=rversor
        rtemp(1)=rtemp(1)+HALF*dex(indlow)
        rtemp(2)=rtemp(2)+HALF*dey(indlow)
        rtemp(3)=rtemp(3)+HALF*dez(indlow)
        vx=vxs+xcross(otemp,rtemp)
        vy=vys+ycross(otemp,rtemp)
        vz=vzs+zcross(otemp,rtemp)
      endif

      f2p = TWO * fluidsub(ii,jj,kk,indhig)- &
       p(indhig)*pref_bouzidi*rhosub(ii,jj,kk)*(dex(indhig)*vx+ &
       dey(indhig)*vy+dez(indhig)*vz)

      ftemp(1) = f2p*dex(indhig)
      ftemp(2) = f2p*dey(indhig)
      ftemp(3) = f2p*dez(indhig)
#ifndef DEBUG_FORCEINT
      fx = fx + ftemp(1)
      fy = fy + ftemp(2)
      fz = fz + ftemp(3)
#else
      A(iatm, l, 1) = A(iatm, l, 1) + ftemp(1)
      A(iatm, l, 2) = A(iatm, l, 2) + ftemp(2)
      A(iatm, l, 3) = A(iatm, l, 3) + ftemp(3)
#endif

      if(lrotate)then
        ttemp(1) = xcross(rtemp,ftemp) * abs(dex(indhig))
        ttemp(2) = ycross(rtemp,ftemp) * abs(dey(indhig))
        ttemp(3) = zcross(rtemp,ftemp) * abs(dez(indhig))
#ifndef DEBUG_FORCEINT
        tx = tx + ttemp(1)
        ty = ty + ttemp(2)
        tz = tz + ttemp(3)
#else
        A(iatm, l, 4) = A(iatm, l, 4) + ttemp(1)
        A(iatm, l, 5) = A(iatm, l, 5) + ttemp(2)
        A(iatm, l, 6) = A(iatm, l, 6) + ttemp(3)
#endif
      endif

      ! if (debug) then
      !   write (iounit(2),*) __FILE__,9466, pimage(ixpbc,ii,nx),&
      !       pimage(iypbc,jj,ny),pimage(izpbc,kk,nz), "pop=", indlow,indhig, &
      !       "ftemp=", ftemp, "ttemp=", ttemp, "f2p", f2p, &
      !       "pop",fluidsub(ii,jj,kk,indhig), "rho",rhosub(ii,jj,kk), &
      !       "v", vx,vy,vz, vxs,vys,vzs
      ! endif
    endif
  endif
  endif


  ii=i+ex(indhig)
  jj=j+ey(indhig)
  kk=k+ez(indhig)
!  ii=pimage(ixpbc,ii,nx)
!  jj=pimage(iypbc,jj,ny)
!  kk=pimage(izpbc,kk,nz)

  io=i+ex(indlow)
  jo=j+ey(indlow)
  ko=k+ez(indlow)
!  io=pimage(ixpbc,io,nx)
!  jo=pimage(iypbc,jo,ny)
!  ko=pimage(izpbc,ko,nz)

  if(ii>=minx-1 .and. ii<=maxx+1 .and. jj>=miny-1 .and. jj<=maxy+1 .and. &
   kk>=minz-1 .and. kk<=maxz+1)then
    if(io>=minx-1 .and. io<=maxx+1 .and. jo>=miny-1 .and. jo<=maxy+1 .and. &
     ko>=minz-1 .and. ko<=maxz+1)then

  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
      if(lrotate)then
        rtemp=rversor
        rtemp(1)=rtemp(1)+HALF*dex(indhig)
        rtemp(2)=rtemp(2)+HALF*dey(indhig)
        rtemp(3)=rtemp(3)+HALF*dez(indhig)
        vx=vxs+xcross(otemp,rtemp)
        vy=vys+ycross(otemp,rtemp)
        vz=vzs+zcross(otemp,rtemp)
      endif

      f2p = TWO * fluidsub(ii,jj,kk,indlow) - &
       p(indlow)*pref_bouzidi*rhosub(ii,jj,kk)*(dex(indlow)*vx+dey(indlow)*vy+dez(indlow)*vz)

      ftemp(1) = f2p*dex(indlow)
      ftemp(2) = f2p*dey(indlow)
      ftemp(3) = f2p*dez(indlow)
#ifndef DEBUG_FORCEINT
      fx = fx + ftemp(1)
      fy = fy + ftemp(2)
      fz = fz + ftemp(3)
#else
      A(iatm, l, 1) = A(iatm, l, 1) + ftemp(1)
      A(iatm, l, 2) = A(iatm, l, 2) + ftemp(2)
      A(iatm, l, 3) = A(iatm, l, 3) + ftemp(3)
#endif

      if(lrotate)then
        ttemp(1) = xcross(rtemp,ftemp) * abs(dex(indlow))
        ttemp(2) = ycross(rtemp,ftemp) * abs(dey(indlow))
        ttemp(3) = zcross(rtemp,ftemp) * abs(dez(indlow))
#ifndef DEBUG_FORCEINT
        tx = tx + ttemp(1)
        ty = ty + ttemp(2)
        tz = tz + ttemp(3)
#else
        A(iatm, l, 4) = A(iatm, l, 4) + ttemp(1)
        A(iatm, l, 5) = A(iatm, l, 5) + ttemp(2)
        A(iatm, l, 6) = A(iatm, l, 6) + ttemp(3)
#endif
      endif

      ! if (debug) then
      !   write (iounit(2),*) __FILE__,9538, pimage(ixpbc,ii,nx),&
      !       pimage(iypbc,jj,ny),pimage(izpbc,kk,nz), "pop=", indhig,indlow, &
      !       "ftemp=", ftemp, "ttemp=", ttemp, "f2p", f2p, &
      !       "pop",fluidsub(ii,jj,kk,indlow), "rho",rhosub(ii,jj,kk), &
      !       "v", vx,vy,vz, vxs,vys,vzs
      ! endif
    endif
  endif
  endif
 enddo
  
 end subroutine node_to_particle_bounce_back_bc2
 
 subroutine node_to_particle_bounce_back_bc(lrotate,nstep,i,j,k,rversor, &
   otemp,vxs,vys,vzs,fx,fy,fz,tx,ty,tz,rhosub,fluidsub)
 
!***********************************************************************
!     
!     LBsoft subroutine to apply the bounce back of the fluid 
!     on a particle surface node including the pbc
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
 
  implicit none
  
  logical, intent(in) :: lrotate
  integer, intent(in) :: nstep,i,j,k
  real(kind=PRC), intent(in) :: vxs,vys,vzs
  real(kind=PRC), intent(in), dimension(3) :: rversor,otemp
  real(kind=PRC), intent(inout) :: fx,fy,fz
  real(kind=PRC), intent(inout) :: tx,ty,tz
  real(kind=PRC), allocatable, dimension(:,:,:)  :: rhosub
  real(kind=PRC), allocatable :: fluidsub(:,:,:,:)
  
  real(kind=PRC), parameter :: onesixth=ONE/SIX
  real(kind=PRC), parameter :: pref_bouzidi=TWO/cssq
  real(kind=PRC) :: f2p,vx,vy,vz
  real(kind=PRC), dimension(3) :: rtemp,ftemp
  
  integer :: ii,jj,kk,io,jo,ko, iloop,indlow,indhig
  
  vx=vxs
  vy=vys
  vz=vzs
  rtemp=rversor
  
  !force on particle fx fy fz
  !eq. 11.2 from page 437 Kruger's book "the lattice boltzmann method"


  ii=i+ex(1)
  jj=j
  kk=k
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i+ex(2)
  jo=j
  ko=k
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      if(lrotate)then
        rtemp=rversor
        rtemp(1)=rtemp(1)+HALF*dex(1)
        vx=vxs+xcross(otemp,rtemp)
      endif
      f2p=TWO*fluidsub(ii,jj,kk,2)- &
       p(2)*pref_bouzidi*rhosub(ii,jj,kk)*(dex(2)*vx)
      fx=fx+f2p*dex(2)
      if(lrotate)then
        ftemp(1:3)=ZERO
        ftemp(1)=f2p*dex(2)
        tx=tx+xcross(rtemp,ftemp)
      endif
    endif
  endif
    
  ii=i+ex(2)
  jj=j
  kk=k
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i+ex(1)
  jo=j
  ko=k
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      if(lrotate)then
        rtemp=rversor
        rtemp(1)=rtemp(1)+HALF*dex(2)
        vx=vxs+xcross(otemp,rtemp)
      endif
      f2p=TWO*fluidsub(ii,jj,kk,1)- &
       p(1)*pref_bouzidi*rhosub(ii,jj,kk)*(dex(1)*vx)
      fx=fx+f2p*dex(1)
      if(lrotate)then
        ftemp(1:3)=ZERO
        ftemp(1)=f2p*dex(1)
        tx=tx+xcross(rtemp,ftemp)
      endif
    endif
  endif
    
  ii=i
  jj=j+ey(3)
  kk=k
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i
  jo=j+ey(4)
  ko=k
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      if(lrotate)then
        rtemp=rversor
        rtemp(2)=rtemp(2)+HALF*dey(3)
        vy=vys+ycross(otemp,rtemp)
      endif
      f2p=TWO*fluidsub(ii,jj,kk,4)- &
       p(4)*pref_bouzidi*rhosub(ii,jj,kk)*(dey(4)*vy)
      fy=fy+f2p*dey(4)
      if(lrotate)then
        ftemp(1:3)=ZERO
        ftemp(2)=f2p*dey(4)
        ty=ty+ycross(rtemp,ftemp)
      endif
    endif
  endif
    
  ii=i
  jj=j+ey(4)
  kk=k
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i
  jo=j+ey(3)
  ko=k
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      if(lrotate)then
        rtemp=rversor
        rtemp(2)=rtemp(2)+HALF*dey(4)
        vy=vys+ycross(otemp,rtemp)
      endif
      f2p=TWO*fluidsub(ii,jj,kk,3)- &
       p(3)*pref_bouzidi*rhosub(ii,jj,kk)*(dey(3)*vy)
      fy=fy+f2p*dey(3)
      if(lrotate)then
        ftemp(1:3)=ZERO
        ftemp(2)=f2p*dey(3)
        ty=ty+ycross(rtemp,ftemp)
      endif
    endif
  endif
    
  ii=i
  jj=j
  kk=k+ez(5)
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i
  jo=j
  ko=k+ez(6)
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      if(lrotate)then
        rtemp=rversor
        rtemp(3)=rtemp(3)+HALF*dez(5)
        vz=vzs+zcross(otemp,rtemp)
      endif
      f2p=TWO*fluidsub(ii,jj,kk,6)- &
       p(6)*pref_bouzidi*rhosub(ii,jj,kk)*(dez(6)*vz)
      fz=fz+f2p*dez(6)
      if(lrotate)then
        ftemp(1:3)=ZERO
        ftemp(3)=f2p*dez(6)
        tz=tz+zcross(rtemp,ftemp)
      endif
    endif
  endif
    
  ii=i
  jj=j
  kk=k+ez(6)
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i
  jo=j
  ko=k+ez(5)
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      if(lrotate)then
        rtemp=rversor
        rtemp(3)=rtemp(3)+HALF*dez(6)
        vz=vzs+zcross(otemp,rtemp)
      endif
      f2p=TWO*fluidsub(ii,jj,kk,5)- &
       p(5)*pref_bouzidi*rhosub(ii,jj,kk)*(dez(5)*vz)
      fz=fz+f2p*dez(5)
      if(lrotate)then
        ftemp(1:3)=ZERO
        ftemp(3)=f2p*dez(5)
        tz=tz+zcross(rtemp,ftemp)
      endif
    endif
  endif
    
  ii=i+ex(7)
  jj=j+ey(7)
  kk=k
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i+ex(8)
  jo=j+ey(8)
  ko=k
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      if(lrotate)then
        rtemp=rversor
        rtemp(1)=rtemp(1)+HALF*dex(7)
        rtemp(2)=rtemp(2)+HALF*dey(7)
        vx=vxs+xcross(otemp,rtemp)
        vy=vys+ycross(otemp,rtemp)
      endif
      f2p=TWO*fluidsub(ii,jj,kk,8)- &
       p(8)*pref_bouzidi*rhosub(ii,jj,kk)*(dex(8)*vx+dey(8)*vy)
      fx=fx+f2p*dex(8)
      fy=fy+f2p*dey(8)
      if(lrotate)then
        ftemp(1:3)=ZERO
        ftemp(1)=f2p*dex(8)
        ftemp(2)=f2p*dey(8)
        tx=tx+xcross(rtemp,ftemp)
        ty=ty+ycross(rtemp,ftemp)
      endif
    endif
  endif
    
  ii=i+ex(8)
  jj=j+ey(8)
  kk=k
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i+ex(7)
  jo=j+ey(7)
  ko=k
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      if(lrotate)then
        rtemp=rversor
        rtemp(1)=rtemp(1)+HALF*dex(8)
        rtemp(2)=rtemp(2)+HALF*dey(8)
        vx=vxs+xcross(otemp,rtemp)
        vy=vys+ycross(otemp,rtemp)
      endif
      f2p=TWO*fluidsub(ii,jj,kk,7)- &
       p(7)*pref_bouzidi*rhosub(ii,jj,kk)*(dex(7)*vx+dey(7)*vy)
      fx=fx+f2p*dex(7)
      fy=fy+f2p*dey(7)
      if(lrotate)then
        ftemp(1:3)=ZERO
        ftemp(1)=f2p*dex(7)
        ftemp(2)=f2p*dey(7)
        tx=tx+xcross(rtemp,ftemp)
        ty=ty+ycross(rtemp,ftemp)
      endif
    endif
  endif
    
  ii=i+ex(9)
  jj=j+ey(9)
  kk=k
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i+ex(10)
  jo=j+ey(10)
  ko=k
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      if(lrotate)then
        rtemp=rversor
        rtemp(1)=rtemp(1)+HALF*dex(9)
        rtemp(2)=rtemp(2)+HALF*dey(9)
        vx=vxs+xcross(otemp,rtemp)
        vy=vys+ycross(otemp,rtemp)
      endif
      f2p=TWO*fluidsub(ii,jj,kk,10)- &
       p(10)*pref_bouzidi*rhosub(ii,jj,kk)*(dex(10)*vx+dey(10)*vy)
      fx=fx+f2p*dex(10)
      fy=fy+f2p*dey(10)
      if(lrotate)then
        ftemp(1:3)=ZERO
        ftemp(1)=f2p*dex(10)
        ftemp(2)=f2p*dey(10)
        tx=tx+xcross(rtemp,ftemp)
        ty=ty+ycross(rtemp,ftemp)
      endif
    endif
  endif
    
  ii=i+ex(10)
  jj=j+ey(10)
  kk=k
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i+ex(9)
  jo=j+ey(9)
  ko=k
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      if(lrotate)then
        rtemp=rversor
        rtemp(1)=rtemp(1)+HALF*dex(10)
        rtemp(2)=rtemp(2)+HALF*dey(10)
        vx=vxs+xcross(otemp,rtemp)
        vy=vys+ycross(otemp,rtemp)
      endif
      f2p=TWO*fluidsub(ii,jj,kk,9)- &
       p(9)*pref_bouzidi*rhosub(ii,jj,kk)*(dex(9)*vx+dey(9)*vy)
      fx=fx+f2p*dex(9)
      fy=fy+f2p*dey(9)
      if(lrotate)then
        ftemp(1:3)=ZERO
        ftemp(1)=f2p*dex(9)
        ftemp(2)=f2p*dey(9)
        tx=tx+xcross(rtemp,ftemp)
        ty=ty+ycross(rtemp,ftemp)
      endif
    endif
  endif
    
  ii=i+ex(11)
  jj=j
  kk=k+ez(11)
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i+ex(12)
  jo=j
  ko=k+ez(12)
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      if(lrotate)then
        rtemp=rversor
        rtemp(1)=rtemp(1)+HALF*dex(11)
        rtemp(3)=rtemp(3)+HALF*dez(11)
        vx=vxs+xcross(otemp,rtemp)
        vz=vzs+zcross(otemp,rtemp)
      endif
      f2p=TWO*fluidsub(ii,jj,kk,12)- &
       p(12)*pref_bouzidi*rhosub(ii,jj,kk)*(dex(12)*vx+dez(12)*vz)
      fx=fx+f2p*dex(12)
      fz=fz+f2p*dez(12)
      if(lrotate)then
        ftemp(1:3)=ZERO
        ftemp(1)=f2p*dex(12)
        ftemp(3)=f2p*dez(12)
        tx=tx+xcross(rtemp,ftemp)
        tz=tz+zcross(rtemp,ftemp)
      endif
    endif
  endif
    
  ii=i+ex(12)
  jj=j
  kk=k+ez(12)
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i+ex(11)
  jo=j
  ko=k+ez(11)
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      if(lrotate)then
        rtemp=rversor
        rtemp(1)=rtemp(1)+HALF*dex(12)
        rtemp(3)=rtemp(3)+HALF*dez(12)
        vx=vxs+xcross(otemp,rtemp)
        vz=vzs+zcross(otemp,rtemp)
      endif
      f2p=TWO*fluidsub(ii,jj,kk,11)- &
       p(11)*pref_bouzidi*rhosub(ii,jj,kk)*(dex(11)*vx+dez(11)*vz)
      fx=fx+f2p*dex(11)
      fz=fz+f2p*dez(11)
      if(lrotate)then
        ftemp(1:3)=ZERO
        ftemp(1)=f2p*dex(11)
        ftemp(3)=f2p*dez(11)
        tx=tx+xcross(rtemp,ftemp)
        tz=tz+zcross(rtemp,ftemp)
      endif
    endif
  endif
    
  ii=i+ex(13)
  jj=j
  kk=k+ez(13)
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i+ex(14)
  jo=j
  ko=k+ez(14)
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      if(lrotate)then
        rtemp=rversor
        rtemp(1)=rtemp(1)+HALF*dex(13)
        rtemp(3)=rtemp(3)+HALF*dez(13)
        vx=vxs+xcross(otemp,rtemp)
        vz=vzs+zcross(otemp,rtemp)
      endif
      f2p=TWO*fluidsub(ii,jj,kk,14)- &
       p(14)*pref_bouzidi*rhosub(ii,jj,kk)*(dex(14)*vx+dez(14)*vz)
      fx=fx+f2p*dex(14)
      fz=fz+f2p*dez(14)
      if(lrotate)then
        ftemp(1:3)=ZERO
        ftemp(1)=f2p*dex(14)
        ftemp(3)=f2p*dez(14)
        tx=tx+xcross(rtemp,ftemp)
        tz=tz+zcross(rtemp,ftemp)
      endif
    endif
  endif
    
  ii=i+ex(14)
  jj=j
  kk=k+ez(14)
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i+ex(13)
  jo=j
  ko=k+ez(13)
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      if(lrotate)then
        rtemp=rversor
        rtemp(1)=rtemp(1)+HALF*dex(14)
        rtemp(3)=rtemp(3)+HALF*dez(14)
        vx=vxs+xcross(otemp,rtemp)
        vz=vzs+zcross(otemp,rtemp)
      endif
      f2p=TWO*fluidsub(ii,jj,kk,13)- &
       p(13)*pref_bouzidi*rhosub(ii,jj,kk)*(dex(13)*vx+dez(13)*vz)
      fx=fx+f2p*dex(13)
      fz=fz+f2p*dez(13)
      if(lrotate)then
        ftemp(1:3)=ZERO
        ftemp(1)=f2p*dex(13)
        ftemp(3)=f2p*dez(13)
        tx=tx+xcross(rtemp,ftemp)
        tz=tz+zcross(rtemp,ftemp)
      endif
    endif
  endif
    
  ii=i
  jj=j+ey(15)
  kk=k+ez(15)
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i
  jo=j+ey(16)
  ko=k+ez(16)
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      if(lrotate)then
        rtemp=rversor
        rtemp(2)=rtemp(2)+HALF*dey(15)
        rtemp(3)=rtemp(3)+HALF*dez(15)
        vy=vys+ycross(otemp,rtemp)
        vz=vzs+zcross(otemp,rtemp)
      endif
      f2p=TWO*fluidsub(ii,jj,kk,16)- &
       p(16)*pref_bouzidi*rhosub(ii,jj,kk)*(dey(16)*vy+dez(16)*vz)
      fy=fy+f2p*dey(16)
      fz=fz+f2p*dez(16)
      if(lrotate)then
        ftemp(1:3)=ZERO
        ftemp(2)=f2p*dey(16)
        ftemp(3)=f2p*dez(16)
        ty=ty+ycross(rtemp,ftemp)
        tz=tz+zcross(rtemp,ftemp)
      endif
    endif
  endif
    
  ii=i
  jj=j+ey(16)
  kk=k+ez(16)
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i
  jo=j+ey(15)
  ko=k+ez(15)
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      if(lrotate)then
        rtemp=rversor
        rtemp(2)=rtemp(2)+HALF*dey(16)
        rtemp(3)=rtemp(3)+HALF*dez(16)
        vy=vys+ycross(otemp,rtemp)
        vz=vzs+zcross(otemp,rtemp)
      endif
      f2p=TWO*fluidsub(ii,jj,kk,15)- &
       p(15)*pref_bouzidi*rhosub(ii,jj,kk)*(dey(15)*vy+dez(15)*vz)
      fy=fy+f2p*dey(15)
      fz=fz+f2p*dez(15)
      if(lrotate)then
        ftemp(1:3)=ZERO
        ftemp(2)=f2p*dey(15)
        ftemp(3)=f2p*dez(15)
        ty=ty+ycross(rtemp,ftemp)
        tz=tz+zcross(rtemp,ftemp)
      endif
    endif
  endif
    
  ii=i
  jj=j+ey(17)
  kk=k+ez(17)
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i
  jo=j+ey(18)
  ko=k+ez(18)
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      if(lrotate)then
        rtemp=rversor
        rtemp(2)=rtemp(2)+HALF*dey(17)
        rtemp(3)=rtemp(3)+HALF*dez(17)
        vy=vys+ycross(otemp,rtemp)
        vz=vzs+zcross(otemp,rtemp)
      endif
      f2p=TWO*fluidsub(ii,jj,kk,18)- &
       p(18)*pref_bouzidi*rhosub(ii,jj,kk)*(dey(18)*vy+dez(18)*vz)
      fy=fy+f2p*dey(18)
      fz=fz+f2p*dez(18)
      if(lrotate)then
        ftemp(1:3)=ZERO
        ftemp(2)=f2p*dey(18)
        ftemp(3)=f2p*dez(18)
        ty=ty+ycross(rtemp,ftemp)
        tz=tz+zcross(rtemp,ftemp)
      endif
    endif
  endif
    
  ii=i
  jj=j+ey(18)
  kk=k+ez(18)
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i
  jo=j+ey(17)
  ko=k+ez(17)
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      if(lrotate)then
        rtemp=rversor
        rtemp(2)=rtemp(2)+HALF*dey(18)
        rtemp(3)=rtemp(3)+HALF*dez(18)
        vy=vys+ycross(otemp,rtemp)
        vz=vzs+zcross(otemp,rtemp)
      endif
      f2p=TWO*fluidsub(ii,jj,kk,17)- &
       p(17)*pref_bouzidi*rhosub(ii,jj,kk)*(dey(17)*vy+dez(17)*vz)
      fy=fy+f2p*dey(17)
      fz=fz+f2p*dez(17)
      if(lrotate)then
        ftemp(1:3)=ZERO
        ftemp(2)=f2p*dey(17)
        ftemp(3)=f2p*dez(17)
        ty=ty+ycross(rtemp,ftemp)
        tz=tz+zcross(rtemp,ftemp)
      endif
    endif
  endif
 end subroutine node_to_particle_bounce_back_bc

 subroutine particle_to_node_bounce_back_bc2_phase1(lrotate,nstep,i,j,k, fluidsub, debug,iatm)
 
!***********************************************************************
!     
!     LBsoft subroutine to apply the bounce back of the fluid due to 
!     the particle surface
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
 
  implicit none
  
  logical, intent(in) :: lrotate, debug
  integer, intent(in) :: nstep,i,j,k,iatm
  real(kind=PRC), allocatable :: fluidsub(:,:,:,:)
  
  integer :: ii,jj,kk,io,jo,ko, iloop,indlow,indhig
  

  !moving walls bounce-back approach
  !from page 180 Kruger's book "the lattice boltzmann method"
  !NOTE de[x,y,z]=zero eliminated
  
 do iloop = 1, 9
  
  indlow = iloop*2 - 1
  indhig = iloop*2


  ii=i+ex(indlow)
  jj=j+ey(indlow)
  kk=k+ez(indlow)
!  ii=pimage(ixpbc,ii,nx)
!  jj=pimage(iypbc,jj,ny)
!  kk=pimage(izpbc,kk,nz)

  io=i+ex(indhig)
  jo=j+ey(indhig)
  ko=k+ez(indhig)
!  io=pimage(ixpbc,io,nx)
!  jo=pimage(iypbc,jo,ny)
!  ko=pimage(izpbc,ko,nz)

  if(ii>=minx-1 .and. ii<=maxx+1 .and. jj>=miny-1 .and. jj<=maxy+1 .and. &
   kk>=minz-1 .and. kk<=maxz+1)then
    if(io>=minx-1 .and. io<=maxx+1 .and. jo>=miny-1 .and. jo<=maxy+1 .and. &
     ko>=minz-1 .and. ko<=maxz+1)then
       if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
	      fluidsub(i,j,k,indlow) = fluidsub(ii,jj,kk,indhig)
           ! if (debug) then
           !   write (iounit(3),*) __FILE__,__LINE__, ii,jj,kk, "pop=", indlow,indhig, &
           !       "pop now",fluidsub(ii,jj,kk,indhig), "@",i,j,k
           ! endif
       endif
    endif
  endif

  ii=i+ex(indhig)
  jj=j+ey(indhig)
  kk=k+ez(indhig)
!  ii=pimage(ixpbc,ii,nx)
!  jj=pimage(iypbc,jj,ny)
!  kk=pimage(izpbc,kk,nz)

  io=i+ex(indlow)
  jo=j+ey(indlow)
  ko=k+ez(indlow)
!  io=pimage(ixpbc,io,nx)
!  jo=pimage(iypbc,jo,ny)
!  ko=pimage(izpbc,ko,nz)

  if(ii>=minx-1 .and. ii<=maxx+1 .and. jj>=miny-1 .and. jj<=maxy+1 .and. &
   kk>=minz-1 .and. kk<=maxz+1)then
    if(io>=minx-1 .and. io<=maxx+1 .and. jo>=miny-1 .and. jo<=maxy+1 .and. &
     ko>=minz-1 .and. ko<=maxz+1)then
       if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
	      fluidsub(i,j,k,indhig) = fluidsub(ii,jj,kk,indlow)
	      ! if (debug) then
              !   write (iounit(3),*) __FILE__,__LINE__, ii,jj,kk, "pop=", indhig,indlow, &
              !    "pop now",fluidsub(ii,jj,kk,indlow), "@",i,j,k
              ! endif
       endif
    endif
  endif
 enddo
  
 end subroutine particle_to_node_bounce_back_bc2_phase1

 subroutine particle_to_node_bounce_back_bc2_phase2(lrotate,nstep,i,j,k,rversor, &
   otemp,vxs,vys,vzs,rhosub,fluidsub, debug,iatm)

!***********************************************************************
!
!     LBsoft subroutine to apply the bounce back of the fluid due to
!     the particle surface
!
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification November 2018
!
!***********************************************************************

  implicit none

  logical, intent(in) :: lrotate, debug
  integer, intent(in) :: nstep,i,j,k,iatm
  real(kind=PRC), intent(in) :: vxs,vys,vzs
  real(kind=PRC), intent(in), dimension(3) :: rversor,otemp
  real(kind=PRC), allocatable, dimension(:,:,:)  :: rhosub
  real(kind=PRC), allocatable :: fluidsub(:,:,:,:)

  real(kind=PRC), parameter :: onesixth=ONE/SIX
  real(kind=PRC), parameter :: pref_bouzidi=TWO/cssq

  integer :: ii,jj,kk,io,jo,ko, iloop,indlow,indhig
  real(kind=PRC) :: vx,vy,vz
  real(kind=PRC), dimension(3) :: rtemp

  vx=vxs
  vy=vys
  vz=vzs
  rtemp=rversor


  !moving walls bounce-back approach
  !from page 180 Kruger's book "the lattice boltzmann method"
  !NOTE de[x,y,z]=zero eliminated

 do iloop = 1, 9

  indlow = iloop*2 - 1
  indhig = iloop*2


  ii=i+ex(indlow)
  jj=j+ey(indlow)
  kk=k+ez(indlow)
!  ii=pimage(ixpbc,ii,nx)
!  jj=pimage(iypbc,jj,ny)
!  kk=pimage(izpbc,kk,nz)

  io=i+ex(indhig)
  jo=j+ey(indhig)
  ko=k+ez(indhig)
!  io=pimage(ixpbc,io,nx)
!  jo=pimage(iypbc,jo,ny)
!  ko=pimage(izpbc,ko,nz)

  if(ii>=minx-1 .and. ii<=maxx+1 .and. jj>=miny-1 .and. jj<=maxy+1 .and. &
   kk>=minz-1 .and. kk<=maxz+1)then
    if(io>=minx-1 .and. io<=maxx+1 .and. jo>=miny-1 .and. jo<=maxy+1 .and. &
     ko>=minz-1 .and. ko<=maxz+1)then
       if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
          if(lrotate)then
            rtemp=rversor
            rtemp(1)=rtemp(1)+HALF*dex(indlow)
            rtemp(2)=rtemp(2)+HALF*dey(indlow)
            rtemp(3)=rtemp(3)+HALF*dez(indlow)
            vx=vxs+xcross(otemp,rtemp)
            vy=vys+ycross(otemp,rtemp)
            vz=vzs+zcross(otemp,rtemp)
          endif
          fluidsub(i,j,k,indlow)= fluidsub(i,j,k,indlow) - &
           p(indhig)*pref_bouzidi*rhosub(ii,jj,kk)*(dex(indhig)*vx+dey(indhig)*vy+dez(indhig)*vz)
           ! if (debug) then
           !   write (iounit(3),*) __FILE__,__LINE__, ii,jj,kk, "pop=", indlow,indhig, &
           !       "rho",rhosub(ii,jj,kk), &
           !       "->pop low", fluidsub(i,j,k,indlow),i,j,k
           ! endif
       endif
    endif
  endif

  ii=i+ex(indhig)
  jj=j+ey(indhig)
  kk=k+ez(indhig)
!  ii=pimage(ixpbc,ii,nx)
!  jj=pimage(iypbc,jj,ny)
!  kk=pimage(izpbc,kk,nz)

  io=i+ex(indlow)
  jo=j+ey(indlow)
  ko=k+ez(indlow)
!  io=pimage(ixpbc,io,nx)
!  jo=pimage(iypbc,jo,ny)
!  ko=pimage(izpbc,ko,nz)

  if(ii>=minx-1 .and. ii<=maxx+1 .and. jj>=miny-1 .and. jj<=maxy+1 .and. &
   kk>=minz-1 .and. kk<=maxz+1)then
    if(io>=minx-1 .and. io<=maxx+1 .and. jo>=miny-1 .and. jo<=maxy+1 .and. &
     ko>=minz-1 .and. ko<=maxz+1)then
       if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
          if(lrotate)then
            rtemp=rversor
            rtemp(1)=rtemp(1)+HALF*dex(indhig)
            rtemp(2)=rtemp(2)+HALF*dey(indhig)
            rtemp(3)=rtemp(3)+HALF*dez(indhig)
            vx=vxs+xcross(otemp,rtemp)
            vy=vys+ycross(otemp,rtemp)
            vz=vzs+zcross(otemp,rtemp)
          endif
          fluidsub(i,j,k,indhig) = fluidsub(i,j,k,indhig) - &
           p(indlow)*pref_bouzidi*rhosub(ii,jj,kk)*(dex(indlow)*vx+dey(indlow)*vy+dez(indlow)*vz)
          ! if (debug) then
          !    write (iounit(3),*) __FILE__,__LINE__, ii,jj,kk, "pop=", indhig,indlow, &
          !        "rho",rhosub(ii,jj,kk), &
          !        "->pop high", fluidsub(i,j,k,indhig),i,j,k
          !  endif
       endif
    endif
  endif
 enddo

 end subroutine particle_to_node_bounce_back_bc2_phase2

 subroutine particle_to_node_bounce_back_bc(lrotate,nstep,i,j,k,rversor, &
   otemp,vxs,vys,vzs,rhosub,fluidsub)
 
!***********************************************************************
!
!     LBsoft subroutine to apply the bounce back of the fluid due to
!     the particle surface
!
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification November 2018
!
!*********************************************************************** 
 
  implicit none
  
  logical, intent(in) :: lrotate
  integer, intent(in) :: nstep,i,j,k
  real(kind=PRC), intent(in) :: vxs,vys,vzs
  real(kind=PRC), intent(in), dimension(3) :: rversor,otemp
  real(kind=PRC), allocatable, dimension(:,:,:)  :: rhosub
  real(kind=PRC), allocatable :: fluidsub(:,:,:,:)
  
  real(kind=PRC), parameter :: onesixth=ONE/SIX
  real(kind=PRC), parameter :: pref_bouzidi=TWO/cssq
  
  integer :: ii,jj,kk,io,jo,ko
  real(kind=PRC) :: vx,vy,vz
  real(kind=PRC), dimension(3) :: rtemp
  
  vx=vxs
  vy=vys
  vz=vzs
  rtemp=rversor
  

  !moving walls bounce-back approach
  !from page 180 Kruger's book "the lattice boltzmann method"
  !NOTE de[x,y,z]=zero eliminated
  
  ii=i+ex(1)
  jj=j
  kk=k
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i+ex(2)
  jo=j
  ko=k
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
      kk>=minz .and. kk<=maxz)then
      if(lrotate)then
        rtemp=rversor
        rtemp(1)=rtemp(1)+HALF*dex(1)
        vx=vxs+xcross(otemp,rtemp)
      endif
      fluidsub(i,j,k,1)=fluidsub(ii,jj,kk,2)- &
       p(2)*pref_bouzidi*rhosub(ii,jj,kk)*(dex(2)*vx)
    endif
  endif
  
  ii=i+ex(2)
  jj=j
  kk=k
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i+ex(1)
  jo=j
  ko=k
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
      kk>=minz .and. kk<=maxz)then
      if(lrotate)then
        rtemp=rversor
        rtemp(1)=rtemp(1)+HALF*dex(2)
        vx=vxs+xcross(otemp,rtemp)
      endif
      fluidsub(i,j,k,2)=fluidsub(ii,jj,kk,1)- &
       p(1)*pref_bouzidi*rhosub(ii,jj,kk)*(dex(1)*vx)
    endif
  endif
  
  ii=i
  jj=j+ey(3)
  kk=k
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i
  jo=j+ey(4)
  ko=k
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      if(lrotate)then
        rtemp=rversor
        rtemp(2)=rtemp(2)+HALF*dey(3)
        vy=vys+ycross(otemp,rtemp)
      endif
      fluidsub(i,j,k,3)=fluidsub(ii,jj,kk,4)- &
       p(4)*pref_bouzidi*rhosub(ii,jj,kk)*(dey(4)*vy)
    endif
  endif
  
  ii=i
  jj=j+ey(4)
  kk=k
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i
  jo=j+ey(3)
  ko=k
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      if(lrotate)then
        rtemp=rversor
        rtemp(2)=rtemp(2)+HALF*dey(4)
        vy=vys+ycross(otemp,rtemp)
      endif
      fluidsub(i,j,k,4)=fluidsub(ii,jj,kk,3)- &
       p(3)*pref_bouzidi*rhosub(ii,jj,kk)*(dey(3)*vy)
    endif
  endif
  
  ii=i
  jj=j
  kk=k+ez(5)
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i
  jo=j
  ko=k+ez(6)
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      if(lrotate)then
        rtemp=rversor
        rtemp(3)=rtemp(3)+HALF*dez(5)
        vz=vzs+zcross(otemp,rtemp)
      endif
      fluidsub(i,j,k,5)=fluidsub(ii,jj,kk,6)- &
       p(6)*pref_bouzidi*rhosub(ii,jj,kk)*(dez(6)*vz)
    endif
  endif
  
  ii=i
  jj=j
  kk=k+ez(6)
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i
  jo=j
  ko=k+ez(5)
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      if(lrotate)then
        rtemp=rversor
        rtemp(3)=rtemp(3)+HALF*dez(6)
        vz=vzs+zcross(otemp,rtemp)
      endif
      fluidsub(i,j,k,6)=fluidsub(ii,jj,kk,5)- &
       p(5)*pref_bouzidi*rhosub(ii,jj,kk)*(dez(5)*vz)
    endif
  endif
  
  ii=i+ex(7)
  jj=j+ey(7)
  kk=k
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i+ex(8)
  jo=j+ey(8)
  ko=k
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      if(lrotate)then
        rtemp=rversor
        rtemp(1)=rtemp(1)+HALF*dex(7)
        rtemp(2)=rtemp(2)+HALF*dey(7)
        vx=vxs+xcross(otemp,rtemp)
        vy=vys+ycross(otemp,rtemp)
      endif
      fluidsub(i,j,k,7)=fluidsub(ii,jj,kk,8)- &
       p(8)*pref_bouzidi*rhosub(ii,jj,kk)*(dex(8)*vx+dey(8)*vy)
    endif
  endif
  
  ii=i+ex(8)
  jj=j+ey(8)
  kk=k
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i+ex(7)
  jo=j+ey(7)
  ko=k
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      if(lrotate)then
        rtemp=rversor
        rtemp(1)=rtemp(1)+HALF*dex(8)
        rtemp(2)=rtemp(2)+HALF*dey(8)
        vx=vxs+xcross(otemp,rtemp)
        vy=vys+ycross(otemp,rtemp)
      endif
      fluidsub(i,j,k,8)=fluidsub(ii,jj,kk,7)- &
       p(7)*pref_bouzidi*rhosub(ii,jj,kk)*(dex(7)*vx+dey(7)*vy)
    endif
  endif
  
  ii=i+ex(9)
  jj=j+ey(9)
  kk=k
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i+ex(10)
  jo=j+ey(10)
  ko=k
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      if(lrotate)then
        rtemp=rversor
        rtemp(1)=rtemp(1)+HALF*dex(9)
        rtemp(2)=rtemp(2)+HALF*dey(9)
        vx=vxs+xcross(otemp,rtemp)
        vy=vys+ycross(otemp,rtemp)
      endif
      fluidsub(i,j,k,9)=fluidsub(ii,jj,kk,10)- &
       p(10)*pref_bouzidi*rhosub(ii,jj,kk)*(dex(10)*vx+dey(10)*vy)
    endif
  endif
  
  ii=i+ex(10)
  jj=j+ey(10)
  kk=k
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i+ex(9)
  jo=j+ey(9)
  ko=k
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      if(lrotate)then
        rtemp=rversor
        rtemp(1)=rtemp(1)+HALF*dex(10)
        rtemp(2)=rtemp(2)+HALF*dey(10)
        vx=vxs+xcross(otemp,rtemp)
        vy=vys+ycross(otemp,rtemp)
      endif
      fluidsub(i,j,k,10)=fluidsub(ii,jj,kk,9)- &
       p(9)*pref_bouzidi*rhosub(ii,jj,kk)*(dex(9)*vx+dey(9)*vy)
    endif
  endif
  
  ii=i+ex(11)
  jj=j
  kk=k+ez(11)
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i+ex(12)
  jo=j
  ko=k+ez(12)
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      if(lrotate)then
        rtemp=rversor
        rtemp(1)=rtemp(1)+HALF*dex(11)
        rtemp(3)=rtemp(3)+HALF*dez(11)
        vx=vxs+xcross(otemp,rtemp)
        vz=vzs+zcross(otemp,rtemp)
      endif
      fluidsub(i,j,k,11)=fluidsub(ii,jj,kk,12)- &
       p(12)*pref_bouzidi*rhosub(ii,jj,kk)*(dex(12)*vx+dez(12)*vz)
    endif
  endif
  
  ii=i+ex(12)
  jj=j
  kk=k+ez(12)
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i+ex(11)
  jo=j
  ko=k+ez(11)
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      if(lrotate)then
        rtemp=rversor
        rtemp(1)=rtemp(1)+HALF*dex(12)
        rtemp(3)=rtemp(3)+HALF*dez(12)
        vx=vxs+xcross(otemp,rtemp)
        vz=vzs+zcross(otemp,rtemp)
      endif
      fluidsub(i,j,k,12)=fluidsub(ii,jj,kk,11)- &
       p(11)*pref_bouzidi*rhosub(ii,jj,kk)*(dex(11)*vx+dez(11)*vz)
    endif
  endif
  
  ii=i+ex(13)
  jj=j
  kk=k+ez(13)
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i+ex(14)
  jo=j
  ko=k+ez(14)
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      if(lrotate)then
        rtemp=rversor
        rtemp(1)=rtemp(1)+HALF*dex(13)
        rtemp(3)=rtemp(3)+HALF*dez(13)
        vx=vxs+xcross(otemp,rtemp)
        vz=vzs+zcross(otemp,rtemp)
      endif
      fluidsub(i,j,k,13)=fluidsub(ii,jj,kk,14)- &
       p(14)*pref_bouzidi*rhosub(ii,jj,kk)*(dex(14)*vx+dez(14)*vz)
    endif
  endif
  
  ii=i+ex(14)
  jj=j
  kk=k+ez(14)
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i+ex(13)
  jo=j
  ko=k+ez(13)
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      if(lrotate)then
        rtemp=rversor
        rtemp(1)=rtemp(1)+HALF*dex(14)
        rtemp(3)=rtemp(3)+HALF*dez(14)
        vx=vxs+xcross(otemp,rtemp)
        vz=vzs+zcross(otemp,rtemp)
      endif
      fluidsub(i,j,k,14)=fluidsub(ii,jj,kk,13)- &
       p(13)*pref_bouzidi*rhosub(ii,jj,kk)*(dex(13)*vx+dez(13)*vz)
    endif
  endif
  
  ii=i
  jj=j+ey(15)
  kk=k+ez(15)
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i
  jo=j+ey(16)
  ko=k+ez(16)
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      if(lrotate)then
        rtemp=rversor
        rtemp(2)=rtemp(2)+HALF*dey(15)
        rtemp(3)=rtemp(3)+HALF*dez(15)
        vy=vys+ycross(otemp,rtemp)
        vz=vzs+zcross(otemp,rtemp)
      endif
      fluidsub(i,j,k,15)=fluidsub(ii,jj,kk,16)- &
       p(16)*pref_bouzidi*rhosub(ii,jj,kk)*(dey(16)*vy+dez(16)*vz)
    endif
  endif
  
  ii=i
  jj=j+ey(16)
  kk=k+ez(16)
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i
  jo=j+ey(15)
  ko=k+ez(15)
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      if(lrotate)then
        rtemp=rversor
        rtemp(2)=rtemp(2)+HALF*dey(16)
        rtemp(3)=rtemp(3)+HALF*dez(16)
        vy=vys+ycross(otemp,rtemp)
        vz=vzs+zcross(otemp,rtemp)
      endif
      fluidsub(i,j,k,16)=fluidsub(ii,jj,kk,15)- &
       p(15)*pref_bouzidi*rhosub(ii,jj,kk)*(dey(15)*vy+dez(15)*vz)
    endif
  endif
  
  ii=i
  jj=j+ey(17)
  kk=k+ez(17)
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i
  jo=j+ey(18)
  ko=k+ez(18)
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      if(lrotate)then
        rtemp=rversor
        rtemp(2)=rtemp(2)+HALF*dey(17)
        rtemp(3)=rtemp(3)+HALF*dez(17)
        vy=vys+ycross(otemp,rtemp)
        vz=vzs+zcross(otemp,rtemp)
      endif
      fluidsub(i,j,k,17)=fluidsub(ii,jj,kk,18)- &
       p(18)*pref_bouzidi*rhosub(ii,jj,kk)*(dey(18)*vy+dez(18)*vz)
    endif
  endif
  
  ii=i
  jj=j+ey(18)
  kk=k+ez(18)
  ii=pimage(ixpbc,ii,nx)
  jj=pimage(iypbc,jj,ny)
  kk=pimage(izpbc,kk,nz)
  io=i
  jo=j+ey(17)
  ko=k+ez(17)
  io=pimage(ixpbc,io,nx)
  jo=pimage(iypbc,jo,ny)
  ko=pimage(izpbc,ko,nz)
  if(isfluid(ii,jj,kk)==1 .and. isfluid(io,jo,ko)/=1)then
    if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
     kk>=minz .and. kk<=maxz)then
      if(lrotate)then
        rtemp=rversor
        rtemp(2)=rtemp(2)+HALF*dey(18)
        rtemp(3)=rtemp(3)+HALF*dez(18)
        vy=vys+ycross(otemp,rtemp)
        vz=vzs+zcross(otemp,rtemp)
      endif
      fluidsub(i,j,k,18)=fluidsub(ii,jj,kk,17)- &
       p(17)*pref_bouzidi*rhosub(ii,jj,kk)*(dey(17)*vy+dez(17)*vz)
    endif
  endif
  
  
  return
  
 end subroutine particle_to_node_bounce_back_bc

 subroutine fixPops(debug, nstep,isub,jsub,ksub,nspheres,spherelists)
                 
!***********************************************************************
!     
!     LBsoft subroutine to fix the fluid population with a minimum value 
!     if inside the particle template
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     authors: F. Bonaccorso & M. Lauricella
!     last modification November 2018
!     
!*********************************************************************** 


  implicit none
  logical, intent(in) :: debug
  integer, intent(in) :: nstep,isub,jsub,ksub,nspheres
  integer, allocatable, dimension(:,:), intent(in) :: spherelists

  real(kind=PRC), parameter :: MINPOP = 1.0E-15
  integer :: i,j,k,l,ll,ii,jj,kk


  if(lsingle_fluid)then
    do l=1,nspheres
      i=isub+spherelists(1,l)
      j=jsub+spherelists(2,l)
      k=ksub+spherelists(3,l)
      ii=i
      jj=j
      kk=k

      !apply periodic conditions if necessary
      i=pimage(ixpbc,i,nx)
      j=pimage(iypbc,j,ny)
      k=pimage(izpbc,k,nz)

      if(i>=minx .and. i<=maxx .and. j>=miny .and. j<=maxy .and. k>=minz .and. k<=maxz)then
        do ll=1,links
          if( fluidR(i,j,k,ll)<MINPOP ) fluidR(i,j,k,ll) = MINPOP
        enddo
      endif

      !the fluid bounce back is local so I have to do it
      if(i/=ii.or.j/=jj.or.k/=kk)then
       if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. kk>=minz .and. kk<=maxz)then
         do ll=1,links
           if( fluidR(ii,jj,kk,ll)<MINPOP ) fluidR(ii,jj,kk,ll) = MINPOP
         enddo
       endif
      endif

    enddo
  else
    do l=1,nspheres
      i=isub+spherelists(1,l)
      j=jsub+spherelists(2,l)
      k=ksub+spherelists(3,l)
      ii=i
      jj=j
      kk=k
      !apply periodic conditions if necessary
      i=pimage(ixpbc,i,nx)
      j=pimage(iypbc,j,ny)
      k=pimage(izpbc,k,nz)

      if(i>=minx .and. i<=maxx .and. j>=miny .and. j<=maxy .and. k>=minz .and. k<=maxz)then
        do ll=1,links
          if( fluidR(i,j,k,ll)<MINPOP ) fluidR(i,j,k,ll) = MINPOP
        enddo
        do ll=1,links
          if( fluidB(i,j,k,ll)<MINPOP ) fluidB(i,j,k,ll) = MINPOP
        enddo
      endif

      !the fluid bounce back is local so I have to do it
      if(i/=ii.or.j/=jj.or.k/=kk)then
        if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. kk>=minz .and. kk<=maxz)then
          do ll=1,links
            if( fluidR(ii,jj,kk,ll)<MINPOP )fluidR(ii,jj,kk,ll) = MINPOP
          enddo
          do ll=1,links
            if( fluidB(ii,jj,kk,ll)<MINPOP )fluidB(ii,jj,kk,ll) = MINPOP
          enddo
         endif
      endif
    enddo
  endif

  ! if (debug) then
  !  do l=1,3
  !   close(iounit(l))
  !  enddo
  ! endif
 end subroutine fixPops

 subroutine helper_particle_delete_fluids(debug,iatm,itype, isub,jsub,ksub, nspheres,spherelists, lrotate, &
    rdimx,rdimy,rdimz, xx,yy,zz, fx,fy,fz, tx,ty,tz, A)
    
!***********************************************************************
!     
!     LBsoft subroutine to accumulate in A the fluid nodes deleted 
!     if inside the particle template
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     authors: F. Bonaccorso & M. Lauricella
!     last modification November 2018
!     
!***********************************************************************    
    
   implicit none
   integer, intent(in) :: iatm,itype, isub,jsub,ksub,nspheres
   integer, allocatable, dimension(:,:), intent(in) :: spherelists
   logical, intent(in) :: lrotate, debug
   real(kind=PRC), allocatable, dimension(:), intent(in) :: rdimx,rdimy,rdimz
   real(kind=PRC), allocatable, dimension(:), intent(in) :: xx,yy,zz
#ifdef QUAD_FORCEINT
  real(kind=PRC*2), allocatable, dimension(:), intent(inout) :: fx,fy,fz
  real(kind=PRC*2), allocatable, dimension(:), intent(inout) :: tx,ty,tz
#else
  real(kind=PRC), allocatable, dimension(:), intent(inout) :: fx,fy,fz
  real(kind=PRC), allocatable, dimension(:), intent(inout) :: tx,ty,tz
#endif
   real(kind=PRC) :: rtemp(3), ftemp(3),ttemp(3), modr
   integer :: l, i,j,k, ii,jj,kk
   integer, save :: imin,imax,jmin,jmax,kmin,kmax
  logical, save :: lfirst=.true.
#ifdef QUAD_FORCEINT
  real(kind=PRC*2), allocatable, dimension(:,:,:), optional :: A
#else
  real(kind=PRC), allocatable, dimension(:,:,:), optional :: A
#endif


  ttemp = ZERO

  if(lfirst)then
    lfirst=.false.
    imin=minx
    imax=maxx
    jmin=miny
    jmax=maxy
    kmin=minz
    kmax=maxz
  endif

    do l=1,nspheres
        i=isub+spherelists(1,l)
        j=jsub+spherelists(2,l)
        k=ksub+spherelists(3,l)

        ii=i
        jj=j
        kk=k

        !apply periodic conditions if necessary
        i=pimage(ixpbc,i,nx)
        j=pimage(iypbc,j,ny)
        k=pimage(izpbc,k,nz)

        CYCLE_OUT_INTERVAL(i, imin, imax)
        CYCLE_OUT_INTERVAL(j, jmin, jmax)
        CYCLE_OUT_INTERVAL(k, kmin, kmax)

        if(isfluid(i,j,k)==1 .and. new_isfluid(i,j,k)/=1)then
          !fluid node is trasformed to solid
          !formula taken from eq. 18 of PRE 83, 046707 (2011)
          if (lsingle_fluid) then
            ftemp(1)=rhoR(i,j,k)*u(i,j,k)
            ftemp(2)=rhoR(i,j,k)*v(i,j,k)
            ftemp(3)=rhoR(i,j,k)*w(i,j,k)
          else
            ftemp(1)=(rhoR(i,j,k)+rhoB(i,j,k))*u(i,j,k)
            ftemp(2)=(rhoR(i,j,k)+rhoB(i,j,k))*v(i,j,k)
            ftemp(3)=(rhoR(i,j,k)+rhoB(i,j,k))*w(i,j,k)
          endif
#ifndef DEBUG_FORCEINT
          fx(iatm)=fx(iatm) + ftemp(1)
          fy(iatm)=fy(iatm) + ftemp(2)
          fz(iatm)=fz(iatm) + ftemp(3)
#endif
            if(lrotate)then
              rtemp(1)=real(ii,kind=PRC)-xx(iatm)
              rtemp(2)=real(jj,kind=PRC)-yy(iatm)
              rtemp(3)=real(kk,kind=PRC)-zz(iatm)
              modr=modulvec(rtemp)
              rtemp(1:3)=rdimx(itype)/modr*rtemp(1:3)
              ttemp(1) = xcross(rtemp,ftemp)
              ttemp(2) = ycross(rtemp,ftemp)
              ttemp(3) = zcross(rtemp,ftemp)
#ifndef DEBUG_FORCEINT
              tx(iatm)=tx(iatm) + ttemp(1)
              ty(iatm)=ty(iatm) + ttemp(2)
              tz(iatm)=tz(iatm) + ttemp(3)
#endif
            endif

#ifdef DEBUG_FORCEINT
            A(iatm, l, 1) = A(iatm, l, 1) + ftemp(1)
            A(iatm, l, 2) = A(iatm, l, 2) + ftemp(2)
            A(iatm, l, 3) = A(iatm, l, 3) + ftemp(3)
            if(lrotate)then
            A(iatm, l, 4) = A(iatm, l, 4) + ttemp(1)
            A(iatm, l, 5) = A(iatm, l, 5) + ttemp(2)
            A(iatm, l, 6) = A(iatm, l, 6) + ttemp(3)
            endif
#endif

          ! if (debug) then
          !   if (lsingle_fluid) then
          !     write (iounit(1),*) __LINE__,i,j,k, "ftemp", ftemp, "ttemp", ttemp, &
          !        "R", rhoR(i,j,k), "u,v,w=", u(i,j,k),v(i,j,k),w(i,j,k), l
          !   else
          !     write (iounit(1),*) __LINE__,i,j,k, "ftemp", ftemp, "ttemp", ttemp, &
          !        "R", rhoR(i,j,k),"B", rhoB(i,j,k), "u,v,w=", u(i,j,k),v(i,j,k),w(i,j,k), l
          !   endif
          ! endif

        endif
      enddo

 end subroutine helper_particle_delete_fluids


 subroutine particle_delete_fluids(debug,nstep,natmssub,atmbook,nspheres,spherelists, &
   spheredists,nspheredeads,spherelistdeads,lmoved,lrotate,ltype,xx,yy,zz, &
   vx,vy,vz,fx,fy,fz,tx,ty,tz,xo,yo,zo,rdimx,rdimy,rdimz, A)
  
!***********************************************************************
!     
!     LBsoft subroutine to delete fluid nodes according to the 
!     particle presence
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     authors: F. Bonaccorso & M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: nstep,natmssub,nspheres,nspheredeads
  integer, allocatable, intent(in) :: atmbook(:)
  integer, allocatable, dimension(:,:), intent(in) :: spherelists, &
   spherelistdeads
  real(kind=PRC), allocatable, dimension(:), intent(in) :: spheredists
  logical(kind=1), allocatable, dimension(:), intent(in) :: lmoved
  logical, intent(in) :: lrotate, debug
  integer, allocatable, dimension(:), intent(in) :: ltype
  real(kind=PRC), allocatable, dimension(:), intent(in) :: xx,yy,zz
  real(kind=PRC), allocatable, dimension(:), intent(in) :: vx,vy,vz
#ifdef QUAD_FORCEINT
  real(kind=PRC*2), allocatable, dimension(:), intent(inout) :: fx,fy,fz
  real(kind=PRC*2), allocatable, dimension(:), intent(inout) :: tx,ty,tz
#else
  real(kind=PRC), allocatable, dimension(:), intent(inout) :: fx,fy,fz
  real(kind=PRC), allocatable, dimension(:), intent(inout) :: tx,ty,tz
#endif
#ifdef QUAD_FORCEINT
  real(kind=PRC*2), allocatable, dimension(:,:,:), optional :: A
#else
  real(kind=PRC), allocatable, dimension(:,:,:), optional :: A
#endif
  real(kind=PRC), allocatable, dimension(:), intent(in) :: xo,yo,zo
  real(kind=PRC), allocatable, dimension(:), intent(in) :: rdimx,rdimy,rdimz
  integer :: isub,jsub,ksub, iatm,myi, itype
  character(len=120) :: mynamefile

  
  !delete fluid 
  do myi=1,natmssub
    iatm=atmbook(myi)
    if(.not. lmoved(iatm))cycle

    ! if (debug) then
    !     mynamefile=repeat(' ',120)
    !     mynamefile="part_delfluid.atom"//write_fmtnumb(iatm)
    !     iounit(1) = 113
    !     call OpenLogFile(nstep, mynamefile, iounit(1))
    ! endif


    isub=nint(xx(iatm))
    jsub=nint(yy(iatm))
    ksub=nint(zz(iatm))

    itype=ltype(iatm)

    call helper_particle_delete_fluids(debug,iatm,itype, isub,jsub,ksub, nspheres, spherelists, lrotate, &
        rdimx,rdimy,rdimz, xx,yy,zz, fx,fy,fz,tx,ty,tz, A)
    call helper_particle_delete_fluids(debug,iatm,itype, isub,jsub,ksub, nspheredeads, spherelistdeads, lrotate, &
        rdimx,rdimy,rdimz, xx,yy,zz, fx,fy,fz,tx,ty,tz, A)

    ! if (debug) close(iounit(1))
  enddo
  
 end subroutine particle_delete_fluids
 


 subroutine compute_onebelt_density(i,j,k, Rsum, Dsum)
 
!***********************************************************************
!     
!     LBsoft subroutine to compute the average of one fluid by the
!     interpolation over the first belt
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     authors: F. Bonaccorso & M. Lauricella
!     last modification April 2020
!     
!***********************************************************************
 
 
    implicit none
    real(kind=PRC), intent(out) :: Rsum,Dsum
    integer, intent(in) :: i,j,k
    integer :: l, ishift,jshift,kshift

    Rsum=ZERO; Dsum=ZERO
    !compute mean density value eq. 23 of PRE 83, 046707 (2011)
    do l=1,linksd3q27
        ishift=i+exd3q27(l)
        jshift=j+eyd3q27(l)
        kshift=k+ezd3q27(l)
!        ishift=pimage(ixpbc,ishift,nx)
!        jshift=pimage(iypbc,jshift,ny)
!        kshift=pimage(izpbc,kshift,nz)
        CYCLE_OUT_INTERVAL(ishift, minx-nbuff, maxx+nbuff)
        CYCLE_OUT_INTERVAL(jshift, miny-nbuff, maxy+nbuff)
        CYCLE_OUT_INTERVAL(kshift, minz-nbuff, maxz+nbuff)
        if(isfluid(ishift,jshift,kshift)==1 .and. &
            new_isfluid(ishift,jshift,kshift)==1)then
            Rsum=Rsum+rhoR(ishift,jshift,kshift)
            Dsum=Dsum+ONE
        endif
    enddo
    if (Dsum /= ZERO) Rsum=Rsum/Dsum
    
    return
    
 end subroutine compute_onebelt_density

 subroutine compute_onebelt_density_twofluids(i,j,k, Rsum, Bsum, Dsum)
 
!***********************************************************************
!     
!     LBsoft subroutine to compute the average of two fluids by the
!     interpolation over the first belt
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     authors: F. Bonaccorso & M. Lauricella
!     last modification April 2020
!     
!*********************************************************************** 
 
    implicit none
    real(kind=PRC), intent(out) :: Rsum,Bsum,Dsum
    integer, intent(in) :: i,j,k
    integer :: l, ishift,jshift,kshift

    Rsum=ZERO; Bsum=ZERO; Dsum=ZERO
    !compute mean density value
    do l=1,linksd3q27
        ishift=i+exd3q27(l)
        jshift=j+eyd3q27(l)
        kshift=k+ezd3q27(l)
!        ishift=pimage(ixpbc,ishift,nx)
!        jshift=pimage(iypbc,jshift,ny)
!        kshift=pimage(izpbc,kshift,nz)
        CYCLE_OUT_INTERVAL(ishift, minx-nbuff, maxx+nbuff)
        CYCLE_OUT_INTERVAL(jshift, miny-nbuff, maxy+nbuff)
        CYCLE_OUT_INTERVAL(kshift, minz-nbuff, maxz+nbuff)
        if(isfluid(ishift,jshift,kshift)==1 .and. &
            new_isfluid(ishift,jshift,kshift)==1)then
            Rsum=Rsum+rhoR(ishift,jshift,kshift)
            Bsum=Bsum+rhoB(ishift,jshift,kshift)
            Dsum=Dsum+ONE
        endif
    enddo
    if (Dsum /= ZERO) then
        Rsum=Rsum/Dsum
        Bsum=Bsum/Dsum
    endif
    
    return
    
 end subroutine compute_onebelt_density_twofluids

 subroutine compute_secbelt_density(i,j,k, Rsum, Dsum)
 
!***********************************************************************
!     
!     LBsoft subroutine to compute the average of one fluid by the
!     interpolation over the first and second belts
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     authors: F. Bonaccorso & M. Lauricella
!     last modification April 2020
!     
!***********************************************************************  
 
    implicit none
    real(kind=PRC), intent(out) :: Rsum,Dsum
    integer, intent(in) :: i,j,k
    integer :: l, ishift,jshift,kshift

    Rsum=ZERO; Dsum=ZERO
    do l=1,ndouble
        ishift=i+exdouble(l)
        jshift=j+eydouble(l)
        kshift=k+ezdouble(l)
        CYCLE_OUT_INTERVAL(ishift, minx-nbuff, maxx+nbuff)
        CYCLE_OUT_INTERVAL(jshift, miny-nbuff, maxy+nbuff)
        CYCLE_OUT_INTERVAL(kshift, minz-nbuff, maxz+nbuff)
        if(isfluid(ishift,jshift,kshift)==1 .and. &
            new_isfluid(ishift,jshift,kshift)==1)then
            Rsum=Rsum+rhoR(ishift,jshift,kshift)
            Dsum=Dsum+ONE
        endif
    enddo
    
    if (Dsum /= ZERO) Rsum=Rsum/Dsum
    
    return
    
 end subroutine compute_secbelt_density

 subroutine compute_secbelt_density_twofluids(i,j,k, Rsum, Bsum, Dsum)
 
!***********************************************************************
!     
!     LBsoft subroutine to compute the average of two fluids by the
!     interpolation over the first and second belts
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     authors: F. Bonaccorso & M. Lauricella
!     last modification April 2020
!     
!***********************************************************************  
 
    implicit none
    real(kind=PRC), intent(out) :: Rsum,Bsum,Dsum
    integer, intent(in) :: i,j,k
    integer :: l, ishift,jshift,kshift

    Rsum=ZERO
    Bsum=ZERO
    Dsum=ZERO
    do l=1,ndouble
        ishift=i+exdouble(l)
        jshift=j+eydouble(l)
        kshift=k+ezdouble(l)
        CYCLE_OUT_INTERVAL(ishift, minx-nbuff, maxx+nbuff)
        CYCLE_OUT_INTERVAL(jshift, miny-nbuff, maxy+nbuff)
        CYCLE_OUT_INTERVAL(kshift, minz-nbuff, maxz+nbuff)
        if(isfluid(ishift,jshift,kshift)==1 .and. &
            new_isfluid(ishift,jshift,kshift)==1)then
            Rsum=Rsum+rhoR(ishift,jshift,kshift)
            Bsum=Bsum+rhoB(ishift,jshift,kshift)
            Dsum=Dsum+ONE
        endif
    enddo
    if (Dsum /= ZERO) then
        Rsum=Rsum/Dsum
        Bsum=Bsum/Dsum
    endif
    
    return
    
 end subroutine compute_secbelt_density_twofluids


 subroutine fix_onebelt_density(i,j,k, Rsum, Dsum)
 
!***********************************************************************
!     
!     LBsoft subroutine to compute the average of one fluid by the
!     interpolation over the first belt D3q27 stencil and second belt
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     authors: F. Bonaccorso & M. Lauricella
!     last modification April 2020
!     
!*********************************************************************** 
 
    implicit none
    real(kind=PRC), intent(out) :: Rsum,Dsum
    integer, intent(in) :: i,j,k
    integer :: l, ishift,jshift,kshift

    Rsum=ZERO; Dsum=ZERO
    ! We accept old fluid points
    do l=1,linksd3q27
        ishift=i+exd3q27(l)
        jshift=j+eyd3q27(l)
        kshift=k+ezd3q27(l)
        CYCLE_OUT_INTERVAL(ishift, minx-nbuff, maxx+nbuff)
        CYCLE_OUT_INTERVAL(jshift, miny-nbuff, maxy+nbuff)
        CYCLE_OUT_INTERVAL(kshift, minz-nbuff, maxz+nbuff)
        if(isfluid(ishift,jshift,kshift)==1)then
            Rsum=Rsum+rhoR(ishift,jshift,kshift)
            Dsum=Dsum+ONE
        endif
    enddo
    if (Dsum /= ZERO) then
        Rsum=Rsum/Dsum
        return
    endif

    Rsum=ZERO; Dsum=ZERO
    ! We accept old fluid points
    do l=1,ndouble
        ishift=i+exdouble(l)
        jshift=j+eydouble(l)
        kshift=k+ezdouble(l)
        CYCLE_OUT_INTERVAL(ishift, minx-nbuff, maxx+nbuff)
        CYCLE_OUT_INTERVAL(jshift, miny-nbuff, maxy+nbuff)
        CYCLE_OUT_INTERVAL(kshift, minz-nbuff, maxz+nbuff)
        if(isfluid(ishift,jshift,kshift)==1)then
            Rsum=Rsum+rhoR(ishift,jshift,kshift)
            Dsum=Dsum+ONE
        endif
    enddo
    if (Dsum /= ZERO) then
        Rsum=Rsum/Dsum
        return
    endif

    Rsum=ZERO; Dsum=ZERO
    ! We accept also old walls
    do l=1,linksd3q27
        ishift=i+exd3q27(l)
        jshift=j+eyd3q27(l)
        kshift=k+ezd3q27(l)
        CYCLE_OUT_INTERVAL(ishift, minx-nbuff, maxx+nbuff)
        CYCLE_OUT_INTERVAL(jshift, miny-nbuff, maxy+nbuff)
        CYCLE_OUT_INTERVAL(kshift, minz-nbuff, maxz+nbuff)
        if(isfluid(ishift,jshift,kshift)==0 .or. isfluid(ishift,jshift,kshift)==2)then
          if (rhoR(ishift,jshift,kshift) <= MINDENS) cycle
            Rsum=Rsum+rhoR(ishift,jshift,kshift)
            Dsum=Dsum+ONE
        endif
    enddo
    if (Dsum /= ZERO) then
        Rsum=Rsum/Dsum
        return
    endif

    Rsum=ZERO; Dsum=ZERO
    ! We accept also old walls
    do l=1,ndouble
        ishift=i+exdouble(l)
        jshift=j+eydouble(l)
        kshift=k+ezdouble(l)
        CYCLE_OUT_INTERVAL(ishift, minx-nbuff, maxx+nbuff)
        CYCLE_OUT_INTERVAL(jshift, miny-nbuff, maxy+nbuff)
        CYCLE_OUT_INTERVAL(kshift, minz-nbuff, maxz+nbuff)
        if(isfluid(ishift,jshift,kshift)==0 .or. isfluid(ishift,jshift,kshift)==2)then
          if (rhoR(ishift,jshift,kshift) <= MINDENS) cycle
            Rsum=Rsum+rhoR(ishift,jshift,kshift)
            Dsum=Dsum+ONE
        endif
    enddo
    if (Dsum /= ZERO) then
        Rsum=Rsum/Dsum
        return
    endif

    Rsum = MINDENS
    Dsum = -100000
    
    
 end subroutine fix_onebelt_density

 subroutine fix_onebelt_density_twofluids(i,j,k, Rsum, Bsum, Dsum)
 
!***********************************************************************
!     
!     LBsoft subroutine to compute the average of two fluids by the
!     interpolation over the first belt D3q27 stencil and second belt
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification April 2020
!     
!*********************************************************************** 

    implicit none
    real(kind=PRC), intent(out) :: Rsum,Bsum,Dsum
    integer, intent(in) :: i,j,k
    integer :: l, ishift,jshift,kshift

    Rsum=ZERO; Bsum=ZERO; Dsum=ZERO
    ! We accept old fluid points
    do l=1,linksd3q27
        ishift=i+exd3q27(l)
        jshift=j+eyd3q27(l)
        kshift=k+ezd3q27(l)
        CYCLE_OUT_INTERVAL(ishift, minx-nbuff, maxx+nbuff)
        CYCLE_OUT_INTERVAL(jshift, miny-nbuff, maxy+nbuff)
        CYCLE_OUT_INTERVAL(kshift, minz-nbuff, maxz+nbuff)
        if(isfluid(ishift,jshift,kshift)==1)then
            Rsum=Rsum+rhoR(ishift,jshift,kshift)
            Bsum=Bsum+rhoB(ishift,jshift,kshift)
            Dsum=Dsum+ONE
        endif
    enddo
    if (Dsum /= ZERO) then
        Rsum=Rsum/Dsum
        Bsum=Bsum/Dsum
        return
    endif

    Rsum=ZERO; Bsum=ZERO; Dsum=ZERO
    ! We accept old fluid points
    do l=1,ndouble
        ishift=i+exdouble(l)
        jshift=j+eydouble(l)
        kshift=k+ezdouble(l)
        CYCLE_OUT_INTERVAL(ishift, minx-nbuff, maxx+nbuff)
        CYCLE_OUT_INTERVAL(jshift, miny-nbuff, maxy+nbuff)
        CYCLE_OUT_INTERVAL(kshift, minz-nbuff, maxz+nbuff)
        if(isfluid(ishift,jshift,kshift)==1)then
            Rsum=Rsum+rhoR(ishift,jshift,kshift)
            Bsum=Bsum+rhoB(ishift,jshift,kshift)
            Dsum=Dsum+ONE
        endif
    enddo
    if (Dsum /= ZERO) then
        Rsum=Rsum/Dsum
        Bsum=Bsum/Dsum
        return
    endif

    Rsum=ZERO; Bsum=ZERO; Dsum=ZERO
    ! We accept also old walls
    do l=1,linksd3q27
        ishift=i+exd3q27(l)
        jshift=j+eyd3q27(l)
        kshift=k+ezd3q27(l)
        CYCLE_OUT_INTERVAL(ishift, minx-nbuff, maxx+nbuff)
        CYCLE_OUT_INTERVAL(jshift, miny-nbuff, maxy+nbuff)
        CYCLE_OUT_INTERVAL(kshift, minz-nbuff, maxz+nbuff)
        if(isfluid(ishift,jshift,kshift)==0 .or. isfluid(ishift,jshift,kshift)==2)then
          if (rhoR(ishift,jshift,kshift) <= MINDENS) cycle
            Rsum=Rsum+rhoR(ishift,jshift,kshift)
            Bsum=Bsum+rhoB(ishift,jshift,kshift)
            Dsum=Dsum+ONE
        endif
    enddo
    if (Dsum /= ZERO) then
        Rsum=Rsum/Dsum
        Bsum=Bsum/Dsum
        return
    endif

    Rsum=ZERO; Bsum=ZERO; Dsum=ZERO
    ! We accept also old walls
    do l=1,ndouble
        ishift=i+exdouble(l)
        jshift=j+eydouble(l)
        kshift=k+ezdouble(l)
        CYCLE_OUT_INTERVAL(ishift, minx-nbuff, maxx+nbuff)
        CYCLE_OUT_INTERVAL(jshift, miny-nbuff, maxy+nbuff)
        CYCLE_OUT_INTERVAL(kshift, minz-nbuff, maxz+nbuff)
        if(isfluid(ishift,jshift,kshift)==0 .or. isfluid(ishift,jshift,kshift)==2)then
          if (rhoR(ishift,jshift,kshift) <= MINDENS) cycle
            Rsum=Rsum+rhoR(ishift,jshift,kshift)
            Bsum=Bsum+rhoB(ishift,jshift,kshift)
            Dsum=Dsum+ONE
        endif
    enddo
    if (Dsum /= ZERO) then
        Rsum=Rsum/Dsum
        Bsum=Bsum/Dsum
        return
    endif

    Rsum = MINDENS
    Bsum = MINDENS
    Dsum = -100000
    
    return
    
 end subroutine fix_onebelt_density_twofluids


 subroutine particle_create_fluids(debug,nstep,natmssub,atmbook,nspheres,spherelists, &
   spheredists,nspheredeads,spherelistdeads,lmoved,lrotate,ltype,xx,yy,zz, &
   vx,vy,vz,fx,fy,fz,tx,ty,tz,xo,yo,zo,rdimx,rdimy,rdimz, A)
  
!***********************************************************************
!     
!     LBsoft subroutine to create fluid nodes according to the 
!     particle presence
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     authors: F. Bonaccorso & M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
  
  implicit none
  integer, intent(in) :: nstep,natmssub,nspheres,nspheredeads
  integer, allocatable, intent(in) :: atmbook(:)
  integer, allocatable, dimension(:,:), intent(in) :: spherelists, &
   spherelistdeads
  real(kind=PRC), allocatable, dimension(:), intent(in) :: spheredists
  logical(kind=1), allocatable, dimension(:), intent(in) :: lmoved
  logical, intent(in) :: lrotate, debug
  integer, allocatable, dimension(:), intent(in) :: ltype
  real(kind=PRC), allocatable, dimension(:), intent(in) :: xx,yy,zz
  real(kind=PRC), allocatable, dimension(:), intent(in) :: vx,vy,vz
#ifdef QUAD_FORCEINT
  real(kind=PRC*2), allocatable, dimension(:), intent(inout) :: fx,fy,fz
  real(kind=PRC*2), allocatable, dimension(:), intent(inout) :: tx,ty,tz
#else
  real(kind=PRC), allocatable, dimension(:), intent(inout) :: fx,fy,fz
  real(kind=PRC), allocatable, dimension(:), intent(inout) :: tx,ty,tz
#endif
  real(kind=PRC), allocatable, dimension(:), intent(in) :: xo,yo,zo
  real(kind=PRC), allocatable, dimension(:), intent(in) :: rdimx,rdimy,rdimz
#ifdef QUAD_FORCEINT
  real(kind=PRC*2), allocatable, dimension(:,:,:), optional :: A
#else
  real(kind=PRC), allocatable, dimension(:,:,:), optional :: A
#endif
  
  integer :: i,j,k, m,iatm,io,jo,ko,itype,myi
  integer :: ii,jj,kk,ishift,jshift,kshift, errorCount
  integer, save :: imin,imax,jmin,jmax,kmin,kmax
  logical, save :: lfirst=.true.
  logical :: lfind,ltest(1)
  real(kind=PRC) :: Rsum,Bsum,Dsum,myu,myv,myw,ftemp(3),rtemp(3),modr, ttemp(3)
  real(kind=PRC) :: dtemp1,dtemp2,dtemp3,dtemp4
  character(len=120) :: mynamefile

  ttemp = ZERO
  errorCount = 0
  
  if(lfirst)then
    lfirst=.false.
    imin=minx
    imax=maxx
    jmin=miny
    jmax=maxy
    kmin=minz
    kmax=maxz
  endif

  !create fluid 
  ltest(1)=.false.
  do myi=1,natmssub
    iatm=atmbook(myi)
    itype=ltype(iatm)
    io=nint(xo(iatm))
    jo=nint(yo(iatm))
    ko=nint(zo(iatm))

    ! if (debug) then
    !     mynamefile=repeat(' ',120)
    !     mynamefile="part_createfluid.atom"//write_fmtnumb(iatm)
    !     call OpenLogFile(nstep, mynamefile, 117)
    ! endif

    if(.not. lmoved(iatm))cycle

      if(lsingle_fluid)then
        do m=1,nspheres
          i=io+spherelists(1,m)
          j=jo+spherelists(2,m)
          k=ko+spherelists(3,m)
          ii=i
          jj=j
          kk=k
         !apply periodic conditions if necessary
          i=pimage(ixpbc,i,nx)
          j=pimage(iypbc,j,ny)
          k=pimage(izpbc,k,nz)

          CYCLE_OUT_INTERVAL(i, imin, imax)
          CYCLE_OUT_INTERVAL(j, jmin, jmax)
          CYCLE_OUT_INTERVAL(k, kmin, kmax)

          if(new_isfluid(i,j,k)==1 .and. isfluid(i,j,k)/=1)then
            !solid node is trasformed to fluid node
            call compute_onebelt_density(i,j,k, Rsum, Dsum)
            if(Dsum==ZERO)then
              call compute_secbelt_density(i,j,k, Rsum, Dsum)
              if(Dsum==ZERO)then
                call fix_onebelt_density(i,j,k, Rsum, Dsum)
                errorCount = errorCount + nint(Dsum)
                ltest(1)=.true.
              endif
            endif

            myu=vx(iatm)
            myv=vy(iatm)
            myw=vz(iatm)
            rhoR(i,j,k)=Rsum
            u(i,j,k)=myu
            v(i,j,k)=myv
            w(i,j,k)=myw

            !formula taken from eq. 25 of PRE 83, 046707 (2011)
            call initialize_newnode_fluid(i,j,k,Rsum,myu,myv,myw,fluidR)
            !formula taken from eq. 26 of PRE 83, 046707 (2011)
            ftemp(1)=-Rsum*myu
            ftemp(2)=-Rsum*myv
            ftemp(3)=-Rsum*myw
#ifndef DEBUG_FORCEINT
            fx(iatm)=fx(iatm) + ftemp(1)
            fy(iatm)=fy(iatm) + ftemp(2)
            fz(iatm)=fz(iatm) + ftemp(3)
#endif
            if(lrotate)then
              rtemp(1)=real(ii,kind=PRC)-xx(iatm)
              rtemp(2)=real(jj,kind=PRC)-yy(iatm)
              rtemp(3)=real(kk,kind=PRC)-zz(iatm)
              modr=modulvec(rtemp)
              rtemp(1:3)=rdimx(itype)/modr*rtemp(1:3)
              ttemp(1) = xcross(rtemp,ftemp)
              ttemp(2) = ycross(rtemp,ftemp)
              ttemp(3) = zcross(rtemp,ftemp)
#ifndef DEBUG_FORCEINT
              tx(iatm)=tx(iatm) + ttemp(1)
              ty(iatm)=ty(iatm) + ttemp(2)
              tz(iatm)=tz(iatm) + ttemp(3)
#endif
            endif
#ifdef DEBUG_FORCEINT
            A(iatm, m, 1) = A(iatm, m, 1) + ftemp(1)
            A(iatm, m, 2) = A(iatm, m, 2) + ftemp(2)
            A(iatm, m, 3) = A(iatm, m, 3) + ftemp(3)
            if(lrotate)then
            A(iatm, m, 4) = A(iatm, m, 4) + ttemp(1)
            A(iatm, m, 5) = A(iatm, m, 5) + ttemp(2)
            A(iatm, m, 6) = A(iatm, m, 6) + ttemp(3)
            endif
#endif
            ! if (debug) write (117,*) __LINE__, i,j,k, "ftemp", ftemp, "ttemp", ttemp, &
            !      "Rsum",Rsum,"Dsum",Dsum, "v(iatm)=", myu,myv,myw
          endif
        enddo
        do m=1,nspheredeads
          i=io+spherelistdeads(1,m)
          j=jo+spherelistdeads(2,m)
          k=ko+spherelistdeads(3,m)
          ii=i
          jj=j
          kk=k
          !apply periodic conditions if necessary
          i=pimage(ixpbc,i,nx)
          j=pimage(iypbc,j,ny)
          k=pimage(izpbc,k,nz)

          CYCLE_OUT_INTERVAL(i, imin, imax)
          CYCLE_OUT_INTERVAL(j, jmin, jmax)
          CYCLE_OUT_INTERVAL(k, kmin, kmax)
          if(new_isfluid(i,j,k)==1 .and. isfluid(i,j,k)/=1)then
            !solid node is trasformed to fluid node
            call compute_onebelt_density(i,j,k, Rsum, Dsum)
            if(Dsum==ZERO)then
              call compute_secbelt_density(i,j,k, Rsum, Dsum)
              if(Dsum==ZERO)then
                call fix_onebelt_density(i,j,k, Rsum, Dsum)
                errorCount = errorCount + nint(Dsum)
                ltest(1)=.true.
              endif
            endif

            myu=vx(iatm)
            myv=vy(iatm)
            myw=vz(iatm)
            rhoR(i,j,k)=Rsum
            u(i,j,k)=myu
            v(i,j,k)=myv
            w(i,j,k)=myw
            !formula taken from eq. 25 of PRE 83, 046707 (2011)
            call initialize_newnode_fluid(i,j,k,Rsum,myu,myv,myw,fluidR)
            !formula taken from eq. 26 of PRE 83, 046707 (2011)
            ftemp(1)=-Rsum*myu
            ftemp(2)=-Rsum*myv
            ftemp(3)=-Rsum*myw
#ifndef DEBUG_FORCEINT
            fx(iatm)=fx(iatm) + ftemp(1)
            fy(iatm)=fy(iatm) + ftemp(2)
            fz(iatm)=fz(iatm) + ftemp(3)
#endif
            if(lrotate)then
              rtemp(1)=real(ii,kind=PRC)-xx(iatm)
              rtemp(2)=real(jj,kind=PRC)-yy(iatm)
              rtemp(3)=real(kk,kind=PRC)-zz(iatm)
              modr=modulvec(rtemp)
              rtemp(1:3)=rdimx(itype)/modr*rtemp(1:3)
              ttemp(1) = xcross(rtemp,ftemp)
              ttemp(2) = ycross(rtemp,ftemp)
              ttemp(3) = zcross(rtemp,ftemp)
#ifndef DEBUG_FORCEINT
              tx(iatm)=tx(iatm) + ttemp(1)
              ty(iatm)=ty(iatm) + ttemp(2)
              tz(iatm)=tz(iatm) + ttemp(3)
#endif
            endif
#ifdef DEBUG_FORCEINT
            A(iatm, m, 1) = A(iatm, m, 1) + ftemp(1)
            A(iatm, m, 2) = A(iatm, m, 2) + ftemp(2)
            A(iatm, m, 3) = A(iatm, m, 3) + ftemp(3)
            if(lrotate)then
            A(iatm, m, 4) = A(iatm, m, 4) + ttemp(1)
            A(iatm, m, 5) = A(iatm, m, 5) + ttemp(2)
            A(iatm, m, 6) = A(iatm, m, 6) + ttemp(3)
            endif
#endif
            ! if (debug) write (117,*) __LINE__, i,j,k, "ftemp", ftemp, "ttemp", ttemp, &
            !      "Rsum",Rsum,"Dsum",Dsum, "v(iatm)=", myu,myv,myw
          endif
        enddo
        ! We tried to fix the problem
        ! if(ltest(1))exit
      else
        do m=1,nspheres
          i=io+spherelists(1,m)
          j=jo+spherelists(2,m)
          k=ko+spherelists(3,m)
          ii=i
          jj=j
          kk=k
         !apply periodic conditions if necessary
          i=pimage(ixpbc,i,nx)
          j=pimage(iypbc,j,ny)
          k=pimage(izpbc,k,nz)

          CYCLE_OUT_INTERVAL(i, imin, imax)
          CYCLE_OUT_INTERVAL(j, jmin, jmax)
          CYCLE_OUT_INTERVAL(k, kmin, kmax)
          if(new_isfluid(i,j,k)==1 .and. isfluid(i,j,k)/=1)then

            call compute_onebelt_density_twofluids(i,j,k, Rsum, Bsum, Dsum)
            if(Dsum==ZERO)then
              call compute_secbelt_density_twofluids(i,j,k, Rsum, Bsum,Dsum)
              if(Dsum==ZERO)then
                call fix_onebelt_density_twofluids(i,j,k, Rsum, Bsum, Dsum)
                errorCount = errorCount + nint(Dsum)
                ltest(1)=.true.
              endif
            endif

            myu=vx(iatm)
            myv=vy(iatm)
            myw=vz(iatm)
            rhoR(i,j,k)=Rsum
            rhoB(i,j,k)=Bsum
            u(i,j,k)=myu
            v(i,j,k)=myv
            w(i,j,k)=myw
            !formula taken from eq. 25 of PRE 83, 046707 (2011)
            call initialize_newnode_fluid(i,j,k,Rsum,myu,myv,myw,fluidR)
            call initialize_newnode_fluid(i,j,k,Bsum,myu,myv,myw,fluidB)
            !formula taken from eq. 26 of PRE 83, 046707 (2011)
            ftemp(1)= -(Rsum+Bsum)*myu
            ftemp(2)= -(Rsum+Bsum)*myv
            ftemp(3)= -(Rsum+Bsum)*myw
#ifndef DEBUG_FORCEINT
            fx(iatm)=fx(iatm) + ftemp(1)
            fy(iatm)=fy(iatm) + ftemp(2)
            fz(iatm)=fz(iatm) + ftemp(3)
#endif
            if(lrotate)then
              rtemp(1)=real(ii,kind=PRC)-xx(iatm)
              rtemp(2)=real(jj,kind=PRC)-yy(iatm)
              rtemp(3)=real(kk,kind=PRC)-zz(iatm)
              modr=modulvec(rtemp)
              rtemp(1:3)=rdimx(itype)/modr*rtemp(1:3)
              ttemp(1) = xcross(rtemp,ftemp)
              ttemp(2) = ycross(rtemp,ftemp)
              ttemp(3) = zcross(rtemp,ftemp)
#ifndef DEBUG_FORCEINT
              tx(iatm)=tx(iatm) + ttemp(1)
              ty(iatm)=ty(iatm) + ttemp(2)
              tz(iatm)=tz(iatm) + ttemp(3)
#endif
            endif
#ifdef DEBUG_FORCEINT
            A(iatm, m, 1) = A(iatm, m, 1) + ftemp(1)
            A(iatm, m, 2) = A(iatm, m, 2) + ftemp(2)
            A(iatm, m, 3) = A(iatm, m, 3) + ftemp(3)
            if(lrotate)then
            A(iatm, m, 4) = A(iatm, m, 4) + ttemp(1)
            A(iatm, m, 5) = A(iatm, m, 5) + ttemp(2)
            A(iatm, m, 6) = A(iatm, m, 6) + ttemp(3)
            endif
#endif
            ! if (debug) write (117,*) __LINE__, i,j,k, "ftemp", ftemp, "ttemp", ttemp, &
            !      "Rsum",Rsum,"Bsum",Bsum,"Dsum",Dsum,"v(iatm)=", myu,myv,myw
          endif
        enddo
        do m=1,nspheredeads
          i=io+spherelistdeads(1,m)
          j=jo+spherelistdeads(2,m)
          k=ko+spherelistdeads(3,m)
          ii=i
          jj=j
          kk=k
          !apply periodic conditions if necessary
          i=pimage(ixpbc,i,nx)
          j=pimage(iypbc,j,ny)
          k=pimage(izpbc,k,nz)

          CYCLE_OUT_INTERVAL(i, imin, imax)
          CYCLE_OUT_INTERVAL(j, jmin, jmax)
          CYCLE_OUT_INTERVAL(k, kmin, kmax)
          if(new_isfluid(i,j,k)==1 .and. isfluid(i,j,k)/=1)then
            !solid node is trasformed to fluid node
            call compute_onebelt_density_twofluids(i,j,k, Rsum, Bsum, Dsum)
            if(Dsum==ZERO)then
              call compute_secbelt_density_twofluids(i,j,k, Rsum, Bsum,Dsum)
              if(Dsum==ZERO)then
                call fix_onebelt_density_twofluids(i,j,k, Rsum, Bsum, Dsum)
                errorCount = errorCount + nint(Dsum)
                ltest(1)=.true.
              endif
            endif

            myu=vx(iatm)
            myv=vy(iatm)
            myw=vz(iatm)
            rhoR(i,j,k)=Rsum
            rhoB(i,j,k)=Bsum
            u(i,j,k)=myu
            v(i,j,k)=myv
            w(i,j,k)=myw
            !formula taken from eq. 25 of PRE 83, 046707 (2011)
            call initialize_newnode_fluid(i,j,k,Rsum,myu,myv,myw,fluidR)
            call initialize_newnode_fluid(i,j,k,Bsum,myu,myv,myw,fluidB)
            !formula taken from eq. 26 of PRE 83, 046707 (2011)
            ftemp(1)=-(Rsum+Bsum)*myu
            ftemp(2)=-(Rsum+Bsum)*myv
            ftemp(3)=-(Rsum+Bsum)*myw
#ifndef DEBUG_FORCEINT
            fx(iatm)=fx(iatm) + ftemp(1)
            fy(iatm)=fy(iatm) + ftemp(2)
            fz(iatm)=fz(iatm) + ftemp(3)
#endif
            if(lrotate)then
              rtemp(1)=real(ii,kind=PRC)-xx(iatm)
              rtemp(2)=real(jj,kind=PRC)-yy(iatm)
              rtemp(3)=real(kk,kind=PRC)-zz(iatm)
              modr=modulvec(rtemp)
              rtemp(1:3)=rdimx(itype)/modr*rtemp(1:3)
              ttemp(1) = xcross(rtemp,ftemp)
              ttemp(2) = ycross(rtemp,ftemp)
              ttemp(3) = zcross(rtemp,ftemp)
#ifndef DEBUG_FORCEINT
              tx(iatm)=tx(iatm) + ttemp(1)
              ty(iatm)=ty(iatm) + ttemp(2)
              tz(iatm)=tz(iatm) + ttemp(3)
#endif
            endif
#ifdef DEBUG_FORCEINT
            A(iatm, m, 1) = A(iatm, m, 1) + ftemp(1)
            A(iatm, m, 2) = A(iatm, m, 2) + ftemp(2)
            A(iatm, m, 3) = A(iatm, m, 3) + ftemp(3)
            if(lrotate)then
            A(iatm, m, 4) = A(iatm, m, 4) + ttemp(1)
            A(iatm, m, 5) = A(iatm, m, 5) + ttemp(2)
            A(iatm, m, 6) = A(iatm, m, 6) + ttemp(3)
            endif
#endif
            ! if (debug) write (117,*) __LINE__, i,j,k, "ftemp", ftemp, "ttemp", ttemp, &
            !      "Rsum",Rsum,"Bsum",Bsum,"Dsum",Dsum,"v(iatm)=", myu,myv,myw
          endif
        enddo
        ! We tried to fix the problem
        ! if(ltest(1))exit
      endif

      ! if (debug) close(117)
  enddo
  
  !call or_world_larr(ltest,1)
  !if(ltest(1))call error(34)
  if(ltest(1)) write(217,fmt=1004) nstep, idrank, errorCount
1004 format ("step:", I6, " id:",I4,I10," fix in particle_create_fluids", X)

 end subroutine particle_create_fluids
 

 subroutine erase_fluids_in_particles(natmssub,atmbook,nspheres,spherelists, &
   spheredists,nspheredeads,spherelistdeads,ltype,xx,yy,zz,vx,vy,vz, &
   rdimx,rdimy,rdimz)
  
!***********************************************************************
!     
!     LBsoft subroutine to create fluid nodes according to the 
!     particle presence
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     authors: F. Bonaccorso & M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
  
  implicit none
  integer, intent(in) :: natmssub,nspheres,nspheredeads
  integer, allocatable, intent(in) :: atmbook(:)
  integer, allocatable, dimension(:,:), intent(in) :: spherelists, &
   spherelistdeads
  real(kind=PRC), allocatable, dimension(:), intent(in) :: spheredists
  integer, allocatable, dimension(:), intent(in) :: ltype
  real(kind=PRC), allocatable, dimension(:), intent(in) :: xx,yy,zz
  real(kind=PRC), allocatable, dimension(:), intent(in) :: vx,vy,vz
  real(kind=PRC), allocatable, dimension(:), intent(in) :: rdimx,rdimy,rdimz
  
  integer :: i,j,k,l,ll,isub,jsub,ksub,iatm,io,jo,ko,itype
  integer :: ishift,jshift,kshift, myi
  integer, save :: imin,imax,jmin,jmax,kmin,kmax
  logical, save :: lfirst=.true.
  real(kind=PRC) :: Rsum,Bsum,Dsum,ii,jj,kk,eps,Usum,Vsum,Wsum
  
  
  if(lfirst)then
    lfirst=.false.
    imin=minx
    imax=maxx
    jmin=miny
    jmax=maxy
    kmin=minz
    kmax=maxz
  endif
  
  if(lsingle_fluid)then
    do myi=1,natmssub
      iatm=atmbook(myi)
      itype=ltype(iatm)
      isub=nint(xx(iatm))
      jsub=nint(yy(iatm))
      ksub=nint(zz(iatm))

!      CYCLE_OUT_INTERVAL(isub, imin, imax)
!      CYCLE_OUT_INTERVAL(jsub, jmin, jmax)
!      CYCLE_OUT_INTERVAL(ksub, kmin, kmax)

      do l=1,nspheres
        i=isub+spherelists(1,l)
        j=jsub+spherelists(2,l)
        k=ksub+spherelists(3,l)
        ii=real(i,kind=PRC)
        jj=real(j,kind=PRC)
        kk=real(k,kind=PRC)
        !apply periodic conditions if necessary
        i=pimage(ixpbc,i,nx)
        j=pimage(iypbc,j,ny)
        k=pimage(izpbc,k,nz)
        if(i<imin .or. i>imax)cycle
        if(j<jmin .or. j>jmax)cycle
        if(k<kmin .or. k>kmax)cycle
        !solid node is trasformed to fluid node
        Rsum=ZERO
        Dsum=ZERO
        !compute mean density value eq. 23 of PRE 83, 046707 (2011)
        do ll=1,linksd3q27
          ishift=i+exd3q27(ll)
          jshift=j+eyd3q27(ll)
          kshift=k+ezd3q27(ll)
          if(isfluid(ishift,jshift,kshift)==1)then
            Rsum=Rsum+pd3q27(ll)*rhoR(ishift,jshift,kshift)
            Dsum=Dsum+pd3q27(ll)
          endif
        enddo
        if(Dsum/=ZERO)then
          Rsum=Rsum/Dsum
          eps=int_cube_sphere(rdimx(itype),xx(iatm),yy(iatm),zz(iatm), &
           ONE,ii,jj,kk,10)
          rhoR(i,j,k)=eps*Rsum
        else
          rhoR(i,j,k)=MINDENS
        endif
        u(i,j,k)=ZERO
        v(i,j,k)=ZERO
        w(i,j,k)=ZERO
      enddo
      do l=1,nspheredeads
        i=isub+spherelistdeads(1,l)
        j=jsub+spherelistdeads(2,l)
        k=ksub+spherelistdeads(3,l)
        !apply periodic conditions if necessary
        i=pimage(ixpbc,i,nx)
        j=pimage(iypbc,j,ny)
        k=pimage(izpbc,k,nz)
        if(i<imin .or. i>imax)cycle
        if(j<jmin .or. j>jmax)cycle
        if(k<kmin .or. k>kmax)cycle
        rhoR(i,j,k)=MINDENS
        u(i,j,k)=ZERO
        v(i,j,k)=ZERO
        w(i,j,k)=ZERO
      enddo
      rhoR(isub,jsub,ksub)=MINDENS
      u(isub,jsub,ksub)=ZERO
      v(isub,jsub,ksub)=ZERO
      w(isub,jsub,ksub)=ZERO
    enddo
  else
    do myi=1,natmssub
      iatm=atmbook(myi)
      itype=ltype(iatm)
      isub=nint(xx(iatm))
      jsub=nint(yy(iatm))
      ksub=nint(zz(iatm))
      do l=1,nspheres
        i=isub+spherelists(1,l)
        j=jsub+spherelists(2,l)
        k=ksub+spherelists(3,l)
        ii=real(i,kind=PRC)
        jj=real(j,kind=PRC)
        kk=real(k,kind=PRC)
        !apply periodic conditions if necessary
        i=pimage(ixpbc,i,nx)
        j=pimage(iypbc,j,ny)
        k=pimage(izpbc,k,nz)
        if(i<imin .or. i>imax)cycle
        if(j<jmin .or. j>jmax)cycle
        if(k<kmin .or. k>kmax)cycle
        !solid node is trasformed to fluid node
        Rsum=ZERO
        Bsum=ZERO
        Dsum=ZERO
        !compute mean density value eq. 23 of PRE 83, 046707 (2011)
        do ll=1,linksd3q27
          ishift=i+exd3q27(ll)
          jshift=j+eyd3q27(ll)
          kshift=k+ezd3q27(ll)
          if(isfluid(ishift,jshift,kshift)==1)then
            Rsum=Rsum+pd3q27(ll)*rhoR(ishift,jshift,kshift)
            Bsum=Bsum+pd3q27(ll)*rhoB(ishift,jshift,kshift)
            Dsum=Dsum+pd3q27(ll)
          endif
        enddo
        if(Dsum/=ZERO)then
          Rsum=Rsum/Dsum
          Bsum=Bsum/Dsum
          eps=int_cube_sphere(rdimx(itype),xx(iatm),yy(iatm),zz(iatm), &
           ONE,ii,jj,kk,10)
          rhoR(i,j,k)=eps*Rsum
          rhoB(i,j,k)=eps*Bsum
        else
          rhoR(i,j,k)=MINDENS
          rhoB(i,j,k)=MINDENS
        endif
        u(i,j,k)=ZERO
        v(i,j,k)=ZERO
        w(i,j,k)=ZERO
      enddo
      do l=1,nspheredeads
        i=isub+spherelistdeads(1,l)
        j=jsub+spherelistdeads(2,l)
        k=ksub+spherelistdeads(3,l)
        !apply periodic conditions if necessary
        i=pimage(ixpbc,i,nx)
        j=pimage(iypbc,j,ny)
        k=pimage(izpbc,k,nz)
        if(i<imin .or. i>imax)cycle
        if(j<jmin .or. j>jmax)cycle
        if(k<kmin .or. k>kmax)cycle
        rhoR(i,j,k)=MINDENS
        rhoB(i,j,k)=MINDENS
        u(i,j,k)=ZERO
        v(i,j,k)=ZERO
        w(i,j,k)=ZERO
      enddo
      rhoR(isub,jsub,ksub)=MINDENS
      rhoB(isub,jsub,ksub)=MINDENS
      u(isub,jsub,ksub)=ZERO
      v(isub,jsub,ksub)=ZERO
      w(isub,jsub,ksub)=ZERO
    enddo
  endif
 end subroutine erase_fluids_in_particles
 
 subroutine initialize_newnode_fluid(i,j,k,rhosub,usub,vsub,wsub,fluidsub)
 
!***********************************************************************
!     
!     LBsoft subroutine for initialize the new fluid node
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: i,j,k
  real(kind=PRC), intent(in) :: rhosub,usub,vsub,wsub
  real(kind=PRC), allocatable :: fluidsub(:,:,:,:)
  
  !formula taken from eq. 19 of PRE 83, 046707 (2011)
    fluidsub(i,j,k,0)=equil_pop00(rhosub,usub,vsub, &
     wsub)
  
  
  
    fluidsub(i,j,k,1)=equil_pop01(rhosub,usub,vsub, &
     wsub)
  
  
  
    fluidsub(i,j,k,2)=equil_pop02(rhosub,usub,vsub, &
     wsub)
  
  
  
    fluidsub(i,j,k,3)=equil_pop03(rhosub,usub,vsub, &
     wsub)
  
  
  
    fluidsub(i,j,k,4)=equil_pop04(rhosub,usub,vsub, &
     wsub)
  
  
  
    fluidsub(i,j,k,5)=equil_pop05(rhosub,usub,vsub, &
     wsub)
  
  
  
    fluidsub(i,j,k,6)=equil_pop06(rhosub,usub,vsub, &
     wsub)
  
  
  
    fluidsub(i,j,k,7)=equil_pop07(rhosub,usub,vsub, &
     wsub)
  
  
  
    fluidsub(i,j,k,8)=equil_pop08(rhosub,usub,vsub, &
     wsub)
  
  
  
    fluidsub(i,j,k,9)=equil_pop09(rhosub,usub,vsub, &
     wsub)
  
  
  
    fluidsub(i,j,k,10)=equil_pop10(rhosub,usub,vsub, &
     wsub)
  
  
  
    fluidsub(i,j,k,11)=equil_pop11(rhosub,usub,vsub, &
     wsub)
  
  
  
    fluidsub(i,j,k,12)=equil_pop12(rhosub,usub,vsub, &
     wsub)
  
  
  
    fluidsub(i,j,k,13)=equil_pop13(rhosub,usub,vsub, &
     wsub)
  
  
  
    fluidsub(i,j,k,14)=equil_pop14(rhosub,usub,vsub, &
     wsub)
  
  
  
    fluidsub(i,j,k,15)=equil_pop15(rhosub,usub,vsub, &
     wsub)
  
  
  
    fluidsub(i,j,k,16)=equil_pop16(rhosub,usub,vsub, &
     wsub)
  
  
  
    fluidsub(i,j,k,17)=equil_pop17(rhosub,usub,vsub, &
     wsub)
  
  
  
    fluidsub(i,j,k,18)=equil_pop18(rhosub,usub,vsub, &
     wsub)
  
  
  return
  
 end subroutine initialize_newnode_fluid
 
!***********END PART TO MANAGE THE INTERACTION WITH PARTICLES***********

 subroutine write_test_map_pop(fluidsub)
  
!***********************************************************************
!     
!     LBsoft subroutine for writing the population
!     in ASCII format for diagnostic purposes
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification july 2018
!     
!***********************************************************************

  implicit none
  
  real(kind=PRC), allocatable :: fluidsub(:,:,:,:)
  
  integer :: i,j,k,l
  
  integer, parameter :: iomap=250
  
  character(len=256) :: file_loc_proc

  
  if(idrank==0)then 
    open(iomap-1,file='global.map',status='replace',action='write')
    write(iomap-1,'(4i10)')mxrank,nx,ny,nz
    close(iomap-1)
  endif

  
  file_loc_proc = 'output'//trim(write_fmtnumb(idrank))//'.map'
  open(iomap,file=trim(file_loc_proc),status='replace',action='write')
  write(iomap,'(i10)')mxrank
    write(iomap,'(8i10)')minx,maxx,miny,maxy,minz,maxz, PRC, 4
    do k=minz,maxz
      do j=miny,maxy
        do i=minx,maxx
          write(iomap,'(3i8)')i,j,k
          do l=0,links
            write(iomap,'(f20.10)')fluidsub(i,j,k,l) 
          enddo
        enddo
      enddo
    enddo

  
  close(iomap)
  
 end subroutine write_test_map_pop
 
 subroutine test_fake_dens(dtemp,ltestout)
 
!***********************************************************************
!     
!     LBsoft subroutine for testing the density in ASCII format 
!     for diagnostic purposes always ordered also in parallel
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2018
!     
!*********************************************************************** 
 
  implicit none
  
  logical :: lwritesub
  real(kind=PRC), intent(inout), allocatable, dimension(:,:,:) :: dtemp
  logical, intent(out) :: ltestout
  integer, parameter :: iotest=166
  character(len=*), parameter :: filetest='fake_test.dat'
  
  real(kind=PRC), allocatable, dimension(:,:,:) :: dtemp2
  integer, dimension(4) :: itempp
  logical :: ltest(1)
  integer :: i,j,k,imy(3)
  logical :: lexist
  
  ltest(1)=.false.
  
  inquire(file=filetest,exist=lexist)
  lwritesub=(.not. lexist)
  
  if(lwritesub)then
    if(mxrank/=1)then
      if(idrank==0)then
        write(6,*)'error in test_fake_dens'
        write(6,*)'to write ',filetest,' you should run in serial'
        write(6,*)'actual mxrank = ',mxrank
        call flush(6)
      endif
      call finalize_world
      stop
    endif
    open(unit=iotest,file=filetest,status='replace',action='write')
    write(iotest,'(4i10)')mxrank,nx,ny,nz
    do k=1-nbuff,nz+nbuff
      do j=1-nbuff,ny+nbuff
        do i=1-nbuff,nx+nbuff
          imy=i4find(nint(dtemp(i,j,k)))
          write(iotest,'(4i8,g20.10,3i8)')i,j,k,i4back(i,j,k),dtemp(i,j,k),imy(1:3)
        enddo
      enddo
    enddo
    close(iotest)
  else
    open(unit=iotest,file=filetest,status='old',action='read')
    read(iotest,'(4i10)')itempp(1:4)
    if(itempp(2)/=nx .or. itempp(3)/=ny .or. itempp(4)/=nz)then
      write(6,*)'error in test_fake_dens'
      write(6,*)'nx ny nz = ',nx,itempp(2),ny,itempp(3),nz,itempp(4)
      call finalize_world
      stop
    endif
    allocate(dtemp2(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff))
    do k=1-nbuff,nz+nbuff
      do j=1-nbuff,ny+nbuff
        do i=1-nbuff,nx+nbuff
          read(iotest,'(4i8,g20.10,3i8)')itempp(1:4),dtemp2(i,j,k),imy(1:3)
        enddo
      enddo
    enddo
    close(iotest)
    do k=minz-1,maxz+1
      do j=miny-1,maxy+1
        do i=minx-1,maxx+1
          if(dtemp(i,j,k)/=dtemp2(i,j,k))then
            ltest(1)=.true.
            goto 125
          endif
        enddo
      enddo
    enddo
 125 continue
    if(ltest(1))then
      write(6,*)'error test_fake_dens idrank = ',idrank
      write(6,*)'node = ',i,j,k,dtemp(i,j,k),dtemp2(i,j,k)
    endif
    deallocate(dtemp2)
    call or_world_larr(ltest,1)
  endif
  
  ltestout=ltest(1)
  
  return
  
 end subroutine test_fake_dens
 
 subroutine test_fake_pops(fluidsub,ltestout)
 
!***********************************************************************
!     
!     LBsoft subroutine for testing the populations in ASCII format 
!     for diagnostic purposes always ordered also in parallel
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
  
  logical :: lwritesub
  real(kind=PRC), allocatable :: fluidsub(:,:,:,:)
  logical, intent(out) :: ltestout
  integer, parameter :: iotest=166
  character(len=60) :: filetest
  
  real(kind=PRC), allocatable, dimension(:,:,:,:) :: dtemp2
  integer, dimension(4) :: itempp,itempp2
  logical :: ltest(1)
  integer :: i,j,k,l,imy(3)
  logical :: lexist
  integer, save :: iter=0
  logical, parameter :: lbinary=.false. 
  integer(kind=1) :: istemp
  
  ltest(1)=.false.
  
  iter=iter+1
  
  
  filetest=repeat(' ',60)
  filetest='fake_test_pop'//write_fmtnumb(iter)//'.dat'
  
  inquire(file=trim(filetest),exist=lexist)
  lwritesub=(.not. lexist)
  
  if(lwritesub)then
    if(mxrank/=1)then
      if(idrank==0)then
        write(6,*)'error in test_fake_dens'
        write(6,*)'to write ',trim(filetest),' you should run in serial'
        write(6,*)'actual mxrank = ',mxrank
        call flush(6)
      endif
      call finalize_world
      stop
    endif
    if(lbinary)then
      open(unit=iotest,file=trim(filetest),status='replace',action='write',form='unformatted')
      write(iotest)mxrank,nx,ny,nz
      do l=0,links
        do k=1-nbuff,nz+nbuff
          do j=1-nbuff,ny+nbuff
            do i=1-nbuff,nx+nbuff
              write(iotest)i,j,k,l,fluidsub(i,j,k,l), &
               isfluid(i,j,k),i+ex(l),j+ey(l),k+ez(l)
            enddo
          enddo
        enddo
      enddo
    else
      open(unit=iotest,file=trim(filetest),status='replace',action='write')
      write(iotest,'(4i10)')mxrank,nx,ny,nz
      do l=0,links
        do k=1-nbuff,nz+nbuff
          do j=1-nbuff,ny+nbuff
            do i=1-nbuff,nx+nbuff
              write(iotest,'(4i8,g20.10,4i4)')i,j,k,l,fluidsub(i,j,k,l), &
               isfluid(i,j,k),i+ex(l),j+ey(l),k+ez(l)
            enddo
          enddo
        enddo
      enddo
    endif
    close(iotest)
  else
    if(lbinary)then
      open(unit=iotest,file=trim(filetest),status='old',action='read',form='unformatted')
      read(iotest)itempp(1:4)
    else
      open(unit=iotest,file=trim(filetest),status='old',action='read')
      read(iotest,'(4i10)')itempp(1:4)
    endif
    if(itempp(2)/=nx .or. itempp(3)/=ny .or. itempp(4)/=nz)then
      write(6,*)'error in test_fake_dens'
      write(6,*)'nx ny nz = ',nx,itempp(2),ny,itempp(3),nz,itempp(4)
      call finalize_world
      stop
    endif
    allocate(dtemp2(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff,0:links))
    if(lbinary)then
      do l=0,links
        do k=1-nbuff,nz+nbuff
          do j=1-nbuff,ny+nbuff
            do i=1-nbuff,nx+nbuff
              read(iotest)itempp(1:4),dtemp2(i,j,k,l),istemp,itempp2(1:3)
            enddo
          enddo
        enddo
      enddo
    else
      do l=0,links
        do k=1-nbuff,nz+nbuff
          do j=1-nbuff,ny+nbuff
            do i=1-nbuff,nx+nbuff
              read(iotest,'(4i8,g20.10,4i4)')itempp(1:4),dtemp2(i,j,k,l),itempp2(1:4)
            enddo
          enddo
        enddo
      enddo
    endif
    close(iotest)
    if(lbinary)then
      do l=0,links
        do k=minz,maxz
          do j=miny,maxy
            do i=minx,maxx
              if(fluidsub(i,j,k,l)/=dtemp2(i,j,k,l))then
                ltest(1)=.true.
                goto 125
              endif
            enddo
          enddo
        enddo
      enddo
    else
      do l=0,links
        do k=minz,maxz
          do j=miny,maxy
            do i=minx,maxx
              if(dabs(fluidsub(i,j,k,l)-dtemp2(i,j,k,l))>1.d-9)then
                ltest(1)=.true.
                goto 125
              endif
            enddo
          enddo
        enddo
      enddo
    endif
 125 continue
    if(ltest(1))then
      write(6,*)'error test_fake_dens idrank = ',idrank
      write(6,*)'node = ',i,j,k,l,fluidsub(i,j,k,l),dtemp2(i,j,k,l)
      call flush(6)
    endif
    deallocate(dtemp2)
    call or_world_larr(ltest,1)
  endif
  
  ltestout=ltest(1)
  if(idrank==0)then
    write(6,*)'ITER: ',iter
    if(lwritesub)then
      write(6,*)'WRITE FILE: ',trim(filetest)
    else
      if(ltestout)then
        write(6,*)'TEST: NO'
      else
        write(6,*)'TEST: OK'
      endif
    endif
    call flush(6)
  endif
  call get_sync_world
  return
  
 end subroutine test_fake_pops
 
 subroutine print_all_pops(iosub,filenam,itersub,fluidsub)
 
!***********************************************************************
!     
!     LBsoft subroutine for writing the populations in ASCII format 
!     for diagnostic purposes always ordered also in parallel
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: iosub,itersub
  character(len=*), intent(in) :: filenam
  real(kind=PRC), allocatable :: fluidsub(:,:,:,:)
  
  character(len=120) :: mynamefile
  integer :: i,j,k,l
  
  !ownern array is mandatory
  if(.not. allocated(ownern))then
    if(idrank==0)then
      write(6,*)'Error in print_all_pops'
      write(6,*)'ownern not allocated'
    endif
    call error(-1)
  endif
  
  mynamefile=repeat(' ',120)
  mynamefile=trim(filenam)//write_fmtnumb(itersub)//'.dat'
  
  if(idrank==0) then
    open(unit=iosub*idrank+23,file=trim(mynamefile),status='replace')
    close(iosub*idrank+23)
  endif
  do k=1,nz
    do j=1,ny
      do i=1,nx
        if(ownern(i4back(i,j,k))==idrank)then
          open(unit=iosub*idrank+23,file=trim(mynamefile),status='old',position='append')
          do l=1,links
            write(iosub*idrank+23,*)i,j,k,l,fluidsub(i,j,k,l)
          enddo
          close(iosub*idrank+23)
        endif
        call get_sync_world
      enddo
    enddo
  enddo
  
  return
  
 end subroutine print_all_pops
 
 subroutine print_all_pops_center(iosub,filenam,itersub, &
  rdimx,rdimy,rdimz,fluidsub,xc,yc,zc,fmiosss,nspheres, &
  spherelists)
 
!***********************************************************************
!     
!     LBsoft subroutine for writing the populations in ASCII format 
!     for diagnostic purposes always ordered also in parallel
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: iosub,itersub,xc,yc,zc,nspheres
  character(len=*), intent(in) :: filenam
  real(kind=PRC), intent(in) :: rdimx,rdimy,rdimz,fmiosss(3,3)
  real(kind=PRC), allocatable :: fluidsub(:,:,:,:)
  integer, allocatable, dimension(:,:), intent(in) :: spherelists
  
  character(len=120) :: mynamefile
  integer :: i,j,k,l,ii,jj,kk,ir,i2,j2,k2,m
  
  !ownern array is mandatory
  if(.not. allocated(ownern))then
    if(idrank==0)then
      write(6,*)'Error in print_all_pops'
      write(6,*)'ownern not allocated'
    endif
    call error(-1)
  endif
  
  mynamefile=repeat(' ',120)
  mynamefile=trim(filenam)//write_fmtnumb(itersub)//'.dat'
  
  if(idrank==0) then
    open(unit=iosub*idrank+23,file=trim(mynamefile),status='replace')
    write(iosub*idrank+23,*)fmiosss
    !close(iosub*idrank+23)
  endif
  
  ir=ceiling(rdimx)+1
  do k=-ir,ir
    do j=-ir,ir
      do i=-ir,ir
        ii=i+xc
        jj=j+yc
        kk=k+zc
        ii=pimage(ixpbc,ii,nx)
        jj=pimage(iypbc,jj,ny)
        kk=pimage(izpbc,kk,nz)
        !if(ownern(i4back(i,j,k))==idrank)then
          if(isfluid(ii,jj,kk)==1)then
!            open(unit=iosub*idrank+23,file=trim(mynamefile), &
!             status='old',position='append')
            do l=1,links
              write(iosub*idrank+23,*)i,j,k,l,isfluid(ii,jj,kk), &
               fluidsub(ii,jj,kk,l)
            enddo
            !write(iosub*idrank+23,*)i,j,k,isfluid(ii,jj,kk)
            !close(iosub*idrank+23)
            !call get_sync_world
          endif
        !endif
      enddo
    enddo
  enddo
  
  do m=1,nspheres
    ii=spherelists(1,m)
    jj=spherelists(2,m)
    kk=spherelists(3,m)
    i=xc+spherelists(1,m)
    j=yc+spherelists(2,m)
    k=zc+spherelists(3,m)
    i=pimage(ixpbc,i,nx)
    j=pimage(iypbc,j,ny)
    k=pimage(izpbc,k,nz)
    do l=1,links
      i2=i+ex(l)
      j2=j+ey(l)
      k2=k+ez(l)
      i2=pimage(ixpbc,i2,nx)
      j2=pimage(iypbc,j2,ny)
      k2=pimage(izpbc,k2,nz)
      if(isfluid(i2,j2,k2)==1)then
        write(iosub*idrank+23,*)ii,jj,kk,l,isfluid(i,j,k), &
               fluidsub(i2,j2,k2,opp(l))
      endif
    enddo
  enddo
  
  close(iosub*idrank+23)
  
  return
  
 end subroutine print_all_pops_center
 
 subroutine print_all_pops_area_shpere(iosub,filenam,itersub, &
  rdimx,rdimy,rdimz,fluidsub,xc,yc,zc,fmiosss,nspheres, &
  spherelists)
 
!***********************************************************************
!     
!     LBsoft subroutine for writing the populations in ASCII format 
!     for diagnostic purposes always ordered also in parallel
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: iosub,itersub,xc,yc,zc,nspheres
  character(len=*), intent(in) :: filenam
  real(kind=PRC), intent(in) :: rdimx,rdimy,rdimz,fmiosss(3,3)
  real(kind=PRC), allocatable :: fluidsub(:,:,:,:)
  integer, allocatable, dimension(:,:), intent(in) :: spherelists
  
  character(len=120) :: mynamefile
  integer :: i,j,k,l,ii,jj,kk,ir,i2,j2,k2,itar,jtar,ktar,m,i3,j3,k3
  
  integer, allocatable :: nmiapop(:,:)
  logical, allocatable :: lmiapop(:,:)
  real(kind=PRC), allocatable :: miapop(:,:)
  
  !ownern array is mandatory
  if(.not. allocated(ownern))then
    if(idrank==0)then
      write(6,*)'Error in print_all_pops'
      write(6,*)'ownern not allocated'
    endif
    call error(-1)
  endif
  
  mynamefile=repeat(' ',120)
  mynamefile=trim(filenam)//write_fmtnumb(itersub)//'.dat'
  
  if(idrank==0) then
    open(unit=iosub*idrank+23,file=trim(mynamefile),status='replace')
    write(iosub*idrank+23,*)fmiosss
    
  endif
  
  allocate(miapop(18,nspheres))
  allocate(nmiapop(18,nspheres))
  allocate(lmiapop(18,nspheres))
  
  miapop(:,:)=-666.d0
  nmiapop(:,:)=0
  lmiapop(:,:)=.false.
  
  do l=1,nspheres
    i=xc+spherelists(1,l)
    j=yc+spherelists(2,l)
    k=zc+spherelists(3,l)
    ii=i
    jj=j
    kk=k
    
    !apply periodic conditions if necessary
    i=pimage(ixpbc,i,nx)
    j=pimage(iypbc,j,ny)
    k=pimage(izpbc,k,nz)
    if(i==ii.and.j==jj.and.k==kk)then
      itar=i
      jtar=j
      ktar=k
      do m=1,links
        i2=itar+ex(m)
        j2=jtar+ey(m)
        k2=ktar+ez(m)
        i3=pimage(ixpbc,i2,nx)
        j3=pimage(iypbc,j2,ny)
        k3=pimage(izpbc,k2,nz)
        if(l==21 .and. m==7)then
          write(6,*)'A ',itar,jtar,ktar,i2,j2,k2,isfluid(i2,j2,k2),isfluid(i3,j3,k3)
        endif
        !i2=pimage(ixpbc,i2,nx)
        !j2=pimage(iypbc,j2,ny)
        !k2=pimage(izpbc,k2,nz)
        if(isfluid(i2,j2,k2)==1)then
          nmiapop(m,l)=nmiapop(m,l)+1
          lmiapop(m,l)=.true.
          miapop(m,l)=fluidsub(itar,jtar,ktar,m)
        endif
      enddo
    else
      itar=i
      jtar=j
      ktar=k
      
      do m=1,links
        i2=itar+ex(m)
        j2=jtar+ey(m)
        k2=ktar+ez(m)
        if(l==21 .and. m==7)then
          write(6,*)'B ',itar,jtar,ktar,i2,j2,k2
        endif
        !i2=pimage(ixpbc,i2,nx)
        !j2=pimage(iypbc,j2,ny)
        !k2=pimage(izpbc,k2,nz)
        if(isfluid(i2,j2,k2)==1)then
          if(i2<=maxx .and. i2>=minx .and. j2<=maxy .and. j2>=miny .and. &
           k2<=maxz .and. k2>=minz)then
            nmiapop(m,l)=nmiapop(m,l)+1
            lmiapop(m,l)=.true.
            miapop(m,l)=fluidsub(itar,jtar,ktar,m)
          endif
        endif
      enddo
      itar=ii
      jtar=jj
      ktar=kk
      do m=1,links
        i2=itar+ex(m)
        j2=jtar+ey(m)
        k2=ktar+ez(m)
        if(l==21 .and. m==7)then
          write(6,*)'C ',itar,jtar,ktar,i2,j2,k2
        endif
        !i2=pimage(ixpbc,i2,nx)
        !j2=pimage(iypbc,j2,ny)
        !k2=pimage(izpbc,k2,nz)
        if(isfluid(i2,j2,k2)==1)then
          if(i2<=maxx .and. i2>=minx .and. j2<=maxy .and. j2>=miny .and. &
           k2<=maxz .and. k2>=minz)then
            nmiapop(m,l)=nmiapop(m,l)+1
            lmiapop(m,l)=.true.
            miapop(m,l)=fluidsub(itar,jtar,ktar,m)
          endif
        endif
      enddo
    endif
  enddo
  
  
  
  if(any(nmiapop>1))then
    write(6,*)'problem!'
    do l=1,nspheres
      do m=1,links
        write(6,*)m,l,nmiapop(m,l)
      enddo
    enddo
  endif
  
  do l=1,nspheres
    do m=1,links
      if(lmiapop(m,l))then
        write(iosub*idrank+23,*)l,m,miapop(m,l)
      endif
    enddo
  enddo
  
  close(iosub*idrank+23)

  
  return
  
 end subroutine print_all_pops_area_shpere
 
 subroutine print_all_hvar(iosub,filenam,itersub,rhosub,usub,vsub, &
   wsub)
 
!***********************************************************************
!     
!     LBsoft subroutine for writing the hvars in ASCII format 
!     for diagnostic purposes always ordered also in parallel
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: iosub,itersub
  character(len=*), intent(in) :: filenam
  real(kind=PRC), allocatable, dimension(:,:,:)  :: rhosub,usub,vsub, &
   wsub
  
  character(len=120) :: mynamefile
  integer :: i,j,k,l
  
  !ownern array is mandatory
  if(.not. allocated(ownern))then
    if(idrank==0)then
      write(6,*)'Error in print_all_hvar'
      write(6,*)'ownern not allocated'
    endif
    call error(-1)
  endif
  
  mynamefile=repeat(' ',120)
  mynamefile=trim(filenam)//write_fmtnumb(itersub)//'.dat'
  
  if(idrank==0) then
    open(unit=iosub*idrank+23,file=trim(mynamefile),status='replace')
    close(iosub*idrank+23)
  endif
  call get_sync_world

  do k=1,nz
    do j=1,ny
      do i=1,nx
        if(ownern(i4back(i,j,k))==idrank)then
          open(unit=iosub*idrank+23,file=trim(mynamefile),status='old',position='append')
          write(iosub*idrank+23,*)i,j,k,rhosub(i,j,k),usub(i,j,k), &
           vsub(i,j,k),wsub(i,j,k)
          close(iosub*idrank+23)
        endif
        call get_sync_world
      enddo
    enddo
  enddo
  
  return
  
 end subroutine print_all_hvar
 
 subroutine setTest
 
!***********************************************************************
!     
!     LBsoft subroutine to set fake pop values for diagnostic
!     purposes
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
  
  implicit none
  integer :: i,j,k,l

  do l=0, links
    forall(i=minx-1:maxx+1, j=miny-1:maxy+1, k=minz-1:maxz+1)
       fluidR(i,j,k,l) = i*1000000 + j*10000 + k*100 + l
     end forall

     where(isfluid(minx-1:maxx+1,miny-1:maxy+1,minz-1:maxz+1)==0)
       fluidR(minx-1:maxx+1,miny-1:maxy+1,minz-1:maxz+1,l) = &
        - fluidR(minx-1:maxx+1,miny-1:maxy+1,minz-1:maxz+1,l)
    end where
  enddo
  
  return
  
 end subroutine setTest
  
 subroutine checkTest(it)
    
!***********************************************************************
!     
!     LBsoft subroutine to check fake pop values for diagnostic
!     purposes
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************

  implicit none
  
  integer, intent(in) :: it
  integer :: i,j,k,l,ishift,jshift,kshift
  integer :: itemp, jtemp, ktemp, ltemp, errtemp
  real(kind=PRC)    :: destVal, val

  if (it /= 1) return

  do l=0, links
    ishift=ex(l)
    jshift=ey(l)
    kshift=ez(l)
    do k=minz,maxz
      do j=miny,maxy
        do i=minx,maxx

          destVal = fluidR(i,j,k,l)

          ltemp = mod(destVal, 100.)
          itemp =  destVal / 10000
          jtemp = (destVal - itemp*10000) / 1000
          ktemp = (destVal - itemp*10000 - jtemp*1000) / 100


          errtemp = 0
          if ( ltemp /= l) then
            errtemp = errtemp + 1
          endif
          ktemp = ktemp+kshift
          if ( ktemp == 0) ktemp=nz
          if ( ktemp == nz+1) ktemp=1
          if ( ktemp /= k) then
            errtemp = errtemp + 100
          endif
          jtemp = jtemp+jshift
          if ( jtemp == 0) jtemp=ny
          if ( jtemp == ny+1) jtemp=1
          if ( jtemp /= j) then
              errtemp = errtemp + 1000
          endif
          itemp = itemp+ishift
          if ( itemp == 0) itemp=nx
          if ( itemp == nx+1) itemp=1
          if ( itemp /= i) then
            errtemp = errtemp + 10000
          endif
          if (errtemp /= 0) then
            val = itemp*10000 + jtemp*1000 +ktemp*100 + l
            write(6,'("Err in stream:", 3I3, " pop:",I3, F8.0, " vs:", F8.0, " errNr:", I6, " rank=", I3)') &
             i,j,k,l,destVal,val,errtemp, idrank
            write(6,*) i,j,k,l,destVal, "vs", itemp,jtemp,ktemp,ltemp, val
          endif
        enddo
      enddo
    enddo
  enddo
        
  return
        
 end subroutine checkTest

 subroutine print_all_pops2(iosub,filenam,itersub)
 
!***********************************************************************
!     
!     LBsoft subroutine for writing the populations in ASCII format 
!     for diagnostic purposes
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
  integer, intent(in) :: iosub,itersub
  character(len=*), intent(in) :: filenam
  character(len=120) :: mynamefile, mynamefile1
  integer :: i,j,k,l, iosub1


  mynamefile=repeat(' ',120)
  mynamefile=trim(filenam)//write_fmtnumb(itersub)//'.'//write_fmtnumb(idrank)//'.dat'
  open(unit=iosub, file=trim(mynamefile), status='replace')

  mynamefile1=repeat(' ',120)
  mynamefile1=trim(filenam)//write_fmtnumb(itersub)//'.'//write_fmtnumb(idrank)//'.dat1'
  iosub1 = iosub + 1
  open(unit=iosub1,file=trim(mynamefile1),status='replace')

  do k=minz-1,maxz+1
    do j=miny-1,maxy+1
      do i=minx-1,maxx+1
          do l=0,links
!!                        1,2,3,4,  5,                  6
            write(iosub,*)i,j,k,l,fluidR(i,j,k,l),fluidB(i,j,k,l)
          enddo
        
!!                        1,2,3, 4,          5,          6,          
          write(iosub1,*) i,j,k, rhoR(i,j,k),rhoB(i,j,k),u(i,j,k), &
!!              7,       8,       9,                  
                v(i,j,k),w(i,j,k),isfluid(i,j,k), &
!!              10,        11,        12,         13,        14,        15,                        
                fuR(i,j,k),fvR(i,j,k),fwR(i,j,k), fuB(i,j,k),fvB(i,j,k),fwB(i,j,k), &
!!              16,          17,         18                        
                psiR(i,j,k), psiB(i,j,k),new_isfluid(i,j,k)
      enddo
    enddo
  enddo

  close(iosub)
  close(iosub1)

 end subroutine print_all_pops2

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! SUBS TO DO DUMP & RESTORE USING 1 FILE PER PROC
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine dumpHvar(nstep)
  
!***********************************************************************
!     
!     LBsoft subroutine for writing the Hvar in ASCII format 
!     for diagnostic purposes
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
  
     implicit none
     integer, intent(in) :: nstep
     character(len=120) :: mynamefile

     mynamefile=repeat(' ',120)
     mynamefile='dumpStat/dumpHVAR.' // write_fmtnumb(nstep) //'.'//write_fmtnumb(idrank)//'.dat'
     open(unit=133,file=trim(mynamefile),form='unformatted',status='replace')
     write(133) rhoR(minx:maxx,miny:maxy,minz:maxz)
     if(.not. lsingle_fluid) write(133) rhoB(minx:maxx,miny:maxy,minz:maxz)
     write(133) u(minx:maxx,miny:maxy,minz:maxz)
     write(133) v(minx:maxx,miny:maxy,minz:maxz)
     write(133) w(minx:maxx,miny:maxy,minz:maxz)
     close(133)
     
  return
     
 end subroutine dumpHvar

 subroutine restoreHvar(nstep)
  
!***********************************************************************
!     
!     LBsoft subroutine for restoring the Hvar in ASCII format 
!     for diagnostic purposes
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
  
  implicit none
  integer, intent(in) :: nstep
  character(len=120) :: mynamefile

  mynamefile=repeat(' ',120)
  mynamefile='dumpHVAR.' // write_fmtnumb(nstep) //'.'//write_fmtnumb(idrank)//'.dat'
  open(unit=133,file=trim(mynamefile),form='unformatted',status='old')
  read(133) rhoR(minx:maxx,miny:maxy,minz:maxz)
  if(.not. lsingle_fluid) read(133) rhoB(minx:maxx,miny:maxy,minz:maxz)
  read(133) u(minx:maxx,miny:maxy,minz:maxz)
  read(133) v(minx:maxx,miny:maxy,minz:maxz)
  read(133) w(minx:maxx,miny:maxy,minz:maxz)
  close(133)
  
  return
  
 end subroutine restoreHvar

 subroutine restorePops(nstep)
 
!***********************************************************************
!     
!     LBsoft subroutine for restoring the pops in ASCII format 
!     for diagnostic purposes
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
  integer, intent(in) :: nstep
  character(len=120) :: mynamefile
  integer :: l

  call print_all_pops2(1001, "red_bef",nstep)

  mynamefile=repeat(' ',120)
  mynamefile='restart.' // write_fmtnumb(nstep) //'.'//write_fmtnumb(idrank)//'.dat'
  open(unit=133,file=trim(mynamefile),form='unformatted',status='old')
  do l=0, links
    read(133) fluidR(minx:maxx,miny:maxy,minz:maxz,l)
    read(133) fluidB(minx:maxx,miny:maxy,minz:maxz,l)
  enddo
  read(133) rhoR(minx:maxx,miny:maxy,minz:maxz)
  read(133) rhoB(minx:maxx,miny:maxy,minz:maxz)
  read(133) u(minx:maxx,miny:maxy,minz:maxz)
  read(133) v(minx:maxx,miny:maxy,minz:maxz)
  read(133) w(minx:maxx,miny:maxy,minz:maxz)
  read(133) isfluid(minx:maxx,miny:maxy,minz:maxz)
  close(133)

  call print_all_pops2(1001, "red_aft",nstep)
     
  return
     
 end subroutine restorePops

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! SUBS TO DO DUMP & RESTORE USING MPI-IO & SINGLE FILE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine dump_oneFile(nstep)
  
!***********************************************************************
!     
!     LBsoft subroutine for print the restart files in binary
!     format
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
     implicit none
     integer, intent(in) :: nstep
     character(len=120) :: mynamefile
     character(len=*), parameter :: mysave='save'
     
     integer :: ierror,myo,nn
     integer, parameter :: myunit=116
     
#if 0
     if(idrank==0) then
             write (6,*) "Making DUMP file for fluid at step:", nstep
     endif
#endif
     
 if(mxrank==1)then
 
   nn=nx*ny*nz
 
 
   myo=myunit+1
   mynamefile=repeat(' ',120)
   mynamefile='dumpPopsR.' // mysave //'.dat'
   call open_file_binary(myo,trim(mynamefile),ierror)
   call print_binary_pop_array(myo,nn, &
    reshape(fluidR(1:nx,1:ny,1:nz,0),(/nn/)), &
    reshape(fluidR(1:nx,1:ny,1:nz,1),(/nn/)), &
    reshape(fluidR(1:nx,1:ny,1:nz,2),(/nn/)), &
    reshape(fluidR(1:nx,1:ny,1:nz,3),(/nn/)), &
    reshape(fluidR(1:nx,1:ny,1:nz,4),(/nn/)), &
    reshape(fluidR(1:nx,1:ny,1:nz,5),(/nn/)), &
    reshape(fluidR(1:nx,1:ny,1:nz,6),(/nn/)), &
    reshape(fluidR(1:nx,1:ny,1:nz,7),(/nn/)), &
    reshape(fluidR(1:nx,1:ny,1:nz,8),(/nn/)), &
    reshape(fluidR(1:nx,1:ny,1:nz,9),(/nn/)), &
    reshape(fluidR(1:nx,1:ny,1:nz,10),(/nn/)), &
    reshape(fluidR(1:nx,1:ny,1:nz,11),(/nn/)), &
    reshape(fluidR(1:nx,1:ny,1:nz,12),(/nn/)), &
    reshape(fluidR(1:nx,1:ny,1:nz,13),(/nn/)), &
    reshape(fluidR(1:nx,1:ny,1:nz,14),(/nn/)), &
    reshape(fluidR(1:nx,1:ny,1:nz,15),(/nn/)), &
    reshape(fluidR(1:nx,1:ny,1:nz,16),(/nn/)), &
    reshape(fluidR(1:nx,1:ny,1:nz,17),(/nn/)), &
    reshape(fluidR(1:nx,1:ny,1:nz,18),(/nn/)), &
    ierror)
   call close_file_binary(myo,ierror)
     
   if(.not. lsingle_fluid) then
     myo=myunit+2
     mynamefile=repeat(' ',120)
     mynamefile='dumpPopsB.' // mysave //'.dat'
     call open_file_binary(myo,trim(mynamefile),ierror)
     call print_binary_pop_array(myo,nn, &
      reshape(fluidB(1:nx,1:ny,1:nz,0),(/nn/)), &
      reshape(fluidB(1:nx,1:ny,1:nz,1),(/nn/)), &
      reshape(fluidB(1:nx,1:ny,1:nz,2),(/nn/)), &
      reshape(fluidB(1:nx,1:ny,1:nz,3),(/nn/)), &
      reshape(fluidB(1:nx,1:ny,1:nz,4),(/nn/)), &
      reshape(fluidB(1:nx,1:ny,1:nz,5),(/nn/)), &
      reshape(fluidB(1:nx,1:ny,1:nz,6),(/nn/)), &
      reshape(fluidB(1:nx,1:ny,1:nz,7),(/nn/)), &
      reshape(fluidB(1:nx,1:ny,1:nz,8),(/nn/)), &
      reshape(fluidB(1:nx,1:ny,1:nz,9),(/nn/)), &
      reshape(fluidB(1:nx,1:ny,1:nz,10),(/nn/)), &
      reshape(fluidB(1:nx,1:ny,1:nz,11),(/nn/)), &
      reshape(fluidB(1:nx,1:ny,1:nz,12),(/nn/)), &
      reshape(fluidB(1:nx,1:ny,1:nz,13),(/nn/)), &
      reshape(fluidB(1:nx,1:ny,1:nz,14),(/nn/)), &
      reshape(fluidB(1:nx,1:ny,1:nz,15),(/nn/)), &
      reshape(fluidB(1:nx,1:ny,1:nz,16),(/nn/)), &
      reshape(fluidB(1:nx,1:ny,1:nz,17),(/nn/)), &
      reshape(fluidB(1:nx,1:ny,1:nz,18),(/nn/)), &
      ierror)
     call close_file_binary(myo,ierror)
   endif
     
   myo=myunit+3
   mynamefile=repeat(' ',120)
   mynamefile='dumpRhoR.' // mysave //'.dat'
   call open_file_binary(myo,trim(mynamefile),ierror)
   call print_binary_1d_array(myo,nn, &
    reshape(rhoR(1:nx,1:ny,1:nz),(/nn/)),ierror)
   call close_file_binary(myo,ierror)

   if(.not. lsingle_fluid) then
     myo=myunit+4
     mynamefile=repeat(' ',120)
     mynamefile='dumpRhoB.' // mysave //'.dat'
     call open_file_binary(myo,trim(mynamefile),ierror)
     call print_binary_1d_array(myo,nn, &
      reshape(rhoB(1:nx,1:ny,1:nz),(/nn/)),ierror)
     call close_file_binary(myo,ierror)
   endif
     
   myo=myunit+5
   mynamefile=repeat(' ',120)
   mynamefile='dumpU.' // mysave //'.dat'
   call open_file_binary(myo,trim(mynamefile),ierror)
   call print_binary_1d_array(myo,nn, &
    reshape(u(1:nx,1:ny,1:nz),(/nn/)),ierror)
   call close_file_binary(myo,ierror)
   
   myo=myunit+6
   mynamefile=repeat(' ',120)
   mynamefile='dumpV.' // mysave //'.dat'
   call open_file_binary(myo,trim(mynamefile),ierror)
   call print_binary_1d_array(myo,nn, &
    reshape(v(1:nx,1:ny,1:nz),(/nn/)),ierror)
   call close_file_binary(myo,ierror)
   
   myo=myunit+7
   mynamefile=repeat(' ',120)
   mynamefile='dumpW.' // mysave //'.dat'
   call open_file_binary(myo,trim(mynamefile),ierror)
   call print_binary_1d_array(myo,nn, &
    reshape(w(1:nx,1:ny,1:nz),(/nn/)),ierror)
   call close_file_binary(myo,ierror)
   
   myo=myunit+8
   mynamefile=repeat(' ',120)
   mynamefile='dumpisFluid.' // mysave //'.dat'
   call open_file_binary(myo,trim(mynamefile),ierror)
   call print_binary_1d_array_int1(myo,nn, &
    reshape(isfluid(1:nx,1:ny,1:nz),(/nn/)),ierror)
   call close_file_binary(myo,ierror)
     
 else
     

   mynamefile=repeat(' ',120)
   !mynamefile='dumpPopsR.' // write_fmtnumb(nstep) //'.dat'
   mynamefile='dumpPopsR.' // mysave //'.dat'

   call collective_writeFile_pops(nstep,mynamefile, &
     fluidR(:,:,:,0),fluidR(:,:,:,1),fluidR(:,:,:,2),fluidR(:,:,:,3), &
     fluidR(:,:,:,4),fluidR(:,:,:,5),fluidR(:,:,:,6),fluidR(:,:,:,7), &
     fluidR(:,:,:,8),fluidR(:,:,:,9),fluidR(:,:,:,10),fluidR(:,:,:,11), &
     fluidR(:,:,:,12),fluidR(:,:,:,13),fluidR(:,:,:,14), &
     fluidR(:,:,:,15),fluidR(:,:,:,16),fluidR(:,:,:,17), &
     fluidR(:,:,:,18),nbuff,minx,maxx,miny,maxy,minz,maxz)

   if(.not. lsingle_fluid) then
    mynamefile=repeat(' ',120)
    mynamefile='dumpPopsB.' // mysave //'.dat'
    call collective_writeFile_pops(nstep,mynamefile, &
     fluidB(:,:,:,0),fluidB(:,:,:,1),fluidB(:,:,:,2),fluidB(:,:,:,3), &
     fluidB(:,:,:,4),fluidB(:,:,:,5),fluidB(:,:,:,6),fluidB(:,:,:,7), &
     fluidB(:,:,:,8),fluidB(:,:,:,9),fluidB(:,:,:,10),fluidB(:,:,:,11), &
     fluidB(:,:,:,12),fluidB(:,:,:,13),fluidB(:,:,:,14), &
     fluidB(:,:,:,15),fluidB(:,:,:,16),fluidB(:,:,:,17), &
     fluidB(:,:,:,18),nbuff,minx,maxx,miny,maxy,minz,maxz)
   endif

   ! Conservative: dump also hvars & isfluid
   mynamefile=repeat(' ',120)
   mynamefile='dumpRhoR.' // mysave //'.dat'
   call collective_writeFile(nstep,mynamefile,rhoR,nbuff, &
    minx,maxx,miny,maxy,minz,maxz)

   if(.not. lsingle_fluid) then
    mynamefile=repeat(' ',120)
    mynamefile='dumpRhoB.' // mysave //'.dat'
    call collective_writeFile(nstep,mynamefile,rhoB,nbuff, &
     minx,maxx,miny,maxy,minz,maxz)
   endif

   mynamefile=repeat(' ',120)
   mynamefile='dumpU.' // mysave //'.dat'
   call collective_writeFile(nstep,mynamefile,u,nbuff, &
    minx,maxx,miny,maxy,minz,maxz)

   mynamefile=repeat(' ',120)
   mynamefile='dumpV.' // mysave //'.dat'
   call collective_writeFile(nstep,mynamefile,v,nbuff, &
    minx,maxx,miny,maxy,minz,maxz)

   mynamefile=repeat(' ',120)
   mynamefile='dumpW.' // mysave //'.dat'
   call collective_writeFile(nstep,mynamefile,w,nbuff, &
    minx,maxx,miny,maxy,minz,maxz)

   mynamefile=repeat(' ',120)
   mynamefile='dumpisFluid.' // mysave //'.dat'
   call collective_writeFile_int(nstep,mynamefile,isfluid,nbuff, &
    minx,maxx,miny,maxy,minz,maxz)


 endif
   
 return

 end subroutine dump_oneFile

 subroutine restore_oneFile(nstep)
  
!***********************************************************************
!     
!     LBsoft subroutine for restoring all data from one dumped 
!     file
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
  
  implicit none
  integer, intent(in) :: nstep
  character(len=120) :: mynamefile
  logical :: lexist

  if(idrank==0) then
    write (6,*) "Restore from DUMP file for fluid at step:", nstep
  endif

     mynamefile=repeat(' ',120)
     mynamefile='dumpPopsR.restart.dat'
     inquire(file=trim(mynamefile),exist=lexist)
     if(.not. lexist)then
       call warning(54,real(1.0,kind=PRC),trim(mynamefile))
       call error(36)
     endif
     call collective_readFile_pops(nstep,mynamefile, &
      fluidR(:,:,:,0),fluidR(:,:,:,1),fluidR(:,:,:,2),fluidR(:,:,:,3), &
      fluidR(:,:,:,4),fluidR(:,:,:,5),fluidR(:,:,:,6),fluidR(:,:,:,7), &
      fluidR(:,:,:,8),fluidR(:,:,:,9),fluidR(:,:,:,10),fluidR(:,:,:,11), &
      fluidR(:,:,:,12),fluidR(:,:,:,13),fluidR(:,:,:,14), &
      fluidR(:,:,:,15),fluidR(:,:,:,16),fluidR(:,:,:,17), &
      fluidR(:,:,:,18),nbuff,minx,maxx,miny,maxy,minz,maxz)

     if(.not. lsingle_fluid) then
      mynamefile=repeat(' ',120)
      mynamefile='dumpPopsB.restart.dat'
      inquire(file=trim(mynamefile),exist=lexist)
      if(.not. lexist)then
        call warning(54,real(1.0,kind=PRC),trim(mynamefile))
        call error(36)
      endif
      call collective_readFile_pops(nstep,mynamefile, &
       fluidB(:,:,:,0),fluidB(:,:,:,1),fluidB(:,:,:,2),fluidB(:,:,:,3), &
       fluidB(:,:,:,4),fluidB(:,:,:,5),fluidB(:,:,:,6),fluidB(:,:,:,7), &
       fluidB(:,:,:,8),fluidB(:,:,:,9),fluidB(:,:,:,10),fluidB(:,:,:,11), &
       fluidB(:,:,:,12),fluidB(:,:,:,13),fluidB(:,:,:,14), &
       fluidB(:,:,:,15),fluidB(:,:,:,16),fluidB(:,:,:,17), &
       fluidB(:,:,:,18),nbuff,minx,maxx,miny,maxy,minz,maxz)
     endif

     ! Conservative: read also hvars & isfluid
     mynamefile=repeat(' ',120)
     mynamefile='dumpRhoR.restart.dat'
     inquire(file=trim(mynamefile),exist=lexist)
     if(.not. lexist)then
       call warning(54,real(1.0,kind=PRC),trim(mynamefile))
       call error(36)
     endif
     call collective_readFile(nstep,mynamefile,rhoR,nbuff, &
      minx,maxx,miny,maxy,minz,maxz)

     if(.not. lsingle_fluid) then
      mynamefile=repeat(' ',120)
      mynamefile='dumpRhoB.restart.dat'
      inquire(file=trim(mynamefile),exist=lexist)
      if(.not. lexist)then
        call warning(54,real(1.0,kind=PRC),trim(mynamefile))
        call error(36)
      endif
      call collective_readFile(nstep,mynamefile,rhoB,nbuff, &
       minx,maxx,miny,maxy,minz,maxz)
     endif

     mynamefile=repeat(' ',120)
     mynamefile='dumpU.restart.dat'
     inquire(file=trim(mynamefile),exist=lexist)
     if(.not. lexist)then
       call warning(54,real(1.0,kind=PRC),trim(mynamefile))
       call error(36)
     endif
     call collective_readFile(nstep,mynamefile,u,nbuff, &
      minx,maxx,miny,maxy,minz,maxz)

     mynamefile=repeat(' ',120)
     mynamefile='dumpV.restart.dat'
     inquire(file=trim(mynamefile),exist=lexist)
     if(.not. lexist)then
       call warning(54,real(1.0,kind=PRC),trim(mynamefile))
       call error(36)
     endif
     call collective_readFile(nstep,mynamefile,v,nbuff, &
      minx,maxx,miny,maxy,minz,maxz)

     mynamefile=repeat(' ',120)
     mynamefile='dumpW.restart.dat'
     inquire(file=trim(mynamefile),exist=lexist)
     if(.not. lexist)then
       call warning(54,real(1.0,kind=PRC),trim(mynamefile))
       call error(36)
     endif
     call collective_readFile(nstep,mynamefile,w,nbuff, &
      minx,maxx,miny,maxy,minz,maxz)

     mynamefile=repeat(' ',120)
     mynamefile='dumpisFluid.restart.dat'
     inquire(file=trim(mynamefile),exist=lexist)
     if(.not. lexist)then
       call warning(54,real(1.0,kind=PRC),trim(mynamefile))
       call error(36)
     endif
     call collective_readFile_int(nstep,mynamefile,isfluid,nbuff, &
      minx,maxx,miny,maxy,minz,maxz)

     call driver_bc_densities
     call driver_bc_velocities
     call driver_bc_isfluid
  end subroutine restore_oneFile

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! SUBS TO CREATE FILE FOR STATS
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine bin_oneFile(nstep)
 
!***********************************************************************
!     
!     LBsoft subroutine for writing output data with double
!     precision in binary format
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
  
  implicit none
  integer, intent(in) :: nstep
  character(len=120) :: mynamefile
  integer, parameter :: myunit=86
  
  integer :: myo,ierror,nn

#if 0
  if(idrank==0) then
    write (6,'(a)') "Making STAT file at step:", nstep
  endif
#endif

  if(mxrank==1)then
     
     nn=nx*ny*nz
     
     myo=myunit+1
     mynamefile=repeat(' ',120)
     mynamefile='dumpBin/RhoR.' // write_fmtnumb(nstep) //'.raw'
     call open_file_binary(myo,trim(mynamefile),ierror)
     call print_binary_1d_array(myo,nn, &
      reshape(rhoR(1:nx,1:ny,1:nz),(/nn/)),ierror)
     call close_file_binary(myo,ierror)

     if(.not. lsingle_fluid) then
       myo=myunit+2
       mynamefile=repeat(' ',120)
       mynamefile='dumpBin/RhoB.' // write_fmtnumb(nstep) //'.raw'
       call open_file_binary(myo,trim(mynamefile),ierror)
       call print_binary_1d_array(myo,nn, &
        reshape(rhoB(1:nx,1:ny,1:nz),(/nn/)),ierror)
       call close_file_binary(myo,ierror)
     endif
     
     myo=myunit+3
     mynamefile=repeat(' ',120)
     mynamefile='dumpBin/U.' // write_fmtnumb(nstep) //'.raw'
     call open_file_binary(myo,trim(mynamefile),ierror)
     call print_binary_1d_array(myo,nn, &
      reshape(u(1:nx,1:ny,1:nz),(/nn/)),ierror)
     call close_file_binary(myo,ierror)
     
     myo=myunit+4
     mynamefile=repeat(' ',120)
     mynamefile='dumpBin/V.' // write_fmtnumb(nstep) //'.raw'
     call open_file_binary(myo,trim(mynamefile),ierror)
     call print_binary_1d_array(myo,nn, &
      reshape(v(1:nx,1:ny,1:nz),(/nn/)),ierror)
     call close_file_binary(myo,ierror)
     
     myo=myunit+5
     mynamefile=repeat(' ',120)
     mynamefile='dumpBin/W.' // write_fmtnumb(nstep) //'.raw'
     call open_file_binary(myo,trim(mynamefile),ierror)
     call print_binary_1d_array(myo,nn, &
      reshape(w(1:nx,1:ny,1:nz),(/nn/)),ierror)
     call close_file_binary(myo,ierror)
     
  else

     mynamefile=repeat(' ',120)
     mynamefile='dumpBin/RhoR.' // write_fmtnumb(nstep) //'.raw'

     call collective_writeFile(nstep,mynamefile,rhoR,nbuff, &
      minx,maxx,miny,maxy,minz,maxz)

     if(.not. lsingle_fluid) then
      mynamefile=repeat(' ',120)
      mynamefile='dumpBin/RhoB.' // write_fmtnumb(nstep) //'.raw'
      call collective_writeFile(nstep,mynamefile,rhoB,nbuff, &
       minx,maxx,miny,maxy,minz,maxz)
     endif

     mynamefile=repeat(' ',120)
     mynamefile='dumpBin/U.' // write_fmtnumb(nstep) //'.raw'
     call collective_writeFile(nstep,mynamefile,u,nbuff, &
      minx,maxx,miny,maxy,minz,maxz)

     mynamefile=repeat(' ',120)
     mynamefile='dumpBin/V.' // write_fmtnumb(nstep) //'.raw'
     call collective_writeFile(nstep,mynamefile,v,nbuff, &
      minx,maxx,miny,maxy,minz,maxz)

     mynamefile=repeat(' ',120)
     mynamefile='dumpBin/W.' // write_fmtnumb(nstep) //'.raw'
     call collective_writeFile(nstep,mynamefile,w,nbuff, &
      minx,maxx,miny,maxy,minz,maxz)

  endif
  
  return
  
 end subroutine bin_oneFile
 
 subroutine open_file_binary(iotest,myname,E_IO)
 
!***********************************************************************
!     
!     LBsoft subroutine for opening a binary stream file
!     in serial IO
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: iotest
  character(len=*), intent(in) :: myname
  integer, intent(out) :: e_io
  
  open(unit=iotest,file=trim(myname),&
   form='unformatted',access='stream',action='write',status='replace', &
   iostat=e_io)
  
  return
  
 endsubroutine open_file_binary
 
 subroutine print_binary_1d_array(iotest,nn,myvar1d,E_IO)
 
!***********************************************************************
!     
!     LBsoft subroutine for writing a vector with double
!     precision in binary format in serial IO
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: iotest,nn
  double precision, intent(in), dimension(nn) :: myvar1d
  integer, intent(out) :: E_IO
  
  integer :: ii
  
  write(iotest,iostat=E_IO)(myvar1d(ii),ii=1,nn)
  
  return
  
 end subroutine print_binary_1d_array
 
 subroutine print_binary_pop_array(iotest,nn,myvar00,myvar01, &
   myvar02,myvar03,myvar04,myvar05,myvar06,myvar07,myvar08,myvar09, &
   myvar10,myvar11,myvar12,myvar13,myvar14,myvar15,myvar16,myvar17, &
   myvar18,E_IO)
 
!***********************************************************************
!     
!     LBsoft subroutine for writing a vector with double
!     precision in binary format in serial IO
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: iotest,nn
  double precision, intent(in), dimension(nn) :: myvar00,myvar01, &
   myvar02,myvar03,myvar04,myvar05,myvar06,myvar07,myvar08,myvar09, &
   myvar10,myvar11,myvar12,myvar13,myvar14,myvar15,myvar16,myvar17, &
   myvar18
  
  integer, intent(out) :: E_IO
  
  integer :: ii
  
  write(iotest,iostat=E_IO)(myvar00(ii),ii=1,nn)
  write(iotest,iostat=E_IO)(myvar01(ii),ii=1,nn)
  write(iotest,iostat=E_IO)(myvar02(ii),ii=1,nn)
  write(iotest,iostat=E_IO)(myvar03(ii),ii=1,nn)
  write(iotest,iostat=E_IO)(myvar04(ii),ii=1,nn)
  write(iotest,iostat=E_IO)(myvar05(ii),ii=1,nn)
  write(iotest,iostat=E_IO)(myvar06(ii),ii=1,nn)
  write(iotest,iostat=E_IO)(myvar07(ii),ii=1,nn)
  write(iotest,iostat=E_IO)(myvar08(ii),ii=1,nn)
  write(iotest,iostat=E_IO)(myvar09(ii),ii=1,nn)
  write(iotest,iostat=E_IO)(myvar10(ii),ii=1,nn)
  write(iotest,iostat=E_IO)(myvar11(ii),ii=1,nn)
  write(iotest,iostat=E_IO)(myvar12(ii),ii=1,nn)
  write(iotest,iostat=E_IO)(myvar13(ii),ii=1,nn)
  write(iotest,iostat=E_IO)(myvar14(ii),ii=1,nn)
  write(iotest,iostat=E_IO)(myvar15(ii),ii=1,nn)
  write(iotest,iostat=E_IO)(myvar16(ii),ii=1,nn)
  write(iotest,iostat=E_IO)(myvar17(ii),ii=1,nn)
  write(iotest,iostat=E_IO)(myvar18(ii),ii=1,nn)
  
  return
  
 end subroutine print_binary_pop_array
 
 subroutine print_binary_1d_array_int1(iotest,nn,myvar1d,E_IO)
 
!***********************************************************************
!     
!     LBsoft subroutine for writing a vector with double
!     precision in binary format in serial IO
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: iotest,nn
  integer(kind=1), intent(in), dimension(nn) :: myvar1d
  integer, intent(out) :: E_IO
  
  integer :: ii
  
  write(iotest,iostat=E_IO)(myvar1d(ii),ii=1,nn)
  
  return
  
 end subroutine print_binary_1d_array_int1
 
 subroutine close_file_binary(iotest,E_IO)
 
!***********************************************************************
!     
!     LBsoft subroutine for closing a binary stream file
!     in serial IO
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: iotest
  integer, intent(out) :: e_io
  
  close(unit=iotest,iostat=e_io)
  
  return
  
 endsubroutine close_file_binary
 
 subroutine minimage_onevec(imcons,dimms,dists)
 
!***********************************************************************
!     
!     LBsoft subroutine to compute athe minimum image convenction
!     on a vector
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     authors: F. Bonaccorso & M. Lauricella
!     last modification April 2020
!     
!*********************************************************************** 
 
  implicit none
  
  integer, intent(in) :: imcons
  real(kind=PRC), intent(in) :: dimms(3)
  real(kind=PRC) :: dists(3)

  integer :: i
  real(kind=PRC) aaa,bbb,ccc
  
  select case(imcons)
  case(0)
    return
  case(1)
    aaa=ONE/dimms(1)
      dists(1)=dists(1)-dimms(1)*nint(aaa*dists(1))
  case(2)
    bbb=ONE/dimms(2)
      dists(2)=dists(2)-dimms(2)*nint(bbb*dists(2))
  case(3)
    aaa=ONE/dimms(1)
    bbb=ONE/dimms(2)
      dists(1)=dists(1)-dimms(1)*nint(aaa*dists(1))
      dists(2)=dists(2)-dimms(2)*nint(bbb*dists(2))
  case(4)
    ccc=ONE/dimms(3)
      dists(3)=dists(3)-dimms(3)*nint(ccc*dists(3))
  case(5)
    aaa=ONE/dimms(1)
    ccc=ONE/dimms(3)
      dists(1)=dists(1)-dimms(1)*nint(aaa*dists(1))
      dists(3)=dists(3)-dimms(3)*nint(ccc*dists(3))
  case(6)
    bbb=ONE/dimms(2)
    ccc=ONE/dimms(3)
      dists(2)=dists(2)-dimms(2)*nint(bbb*dists(2))
      dists(3)=dists(3)-dimms(3)*nint(ccc*dists(3))
  case(7)
    aaa=ONE/dimms(1)
    bbb=ONE/dimms(2)
    ccc=ONE/dimms(3)
      dists(1)=dists(1)-dimms(1)*nint(aaa*dists(1))
      dists(2)=dists(2)-dimms(2)*nint(bbb*dists(2))
      dists(3)=dists(3)-dimms(3)*nint(ccc*dists(3))
  end select
  
 end subroutine minimage_onevec

 end module fluids_mod
