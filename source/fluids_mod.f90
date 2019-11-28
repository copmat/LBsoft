
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
 use aop_mod
 use utility_mod, only : Pi,modulvec,cross,dot,gauss,ibuffservice, &
                   allocate_array_ibuffservice,buffservice,conv_rad, &
                   allocate_array_buffservice,lbuffservice, &
                   allocate_array_lbuffservice,xcross,ycross,zcross, &
                   nbuffservice3d,buffservice3d,int_cube_sphere,fcut, &
                   rand_noseeded,linit_seed,gauss_noseeded,write_fmtnumb, &
                   openLogFile

 
 implicit none
 
 private
 
 integer, save, protected, public :: nx=0
 integer, save, protected, public :: ny=0
 integer, save, protected, public :: nz=0
 
#if LATTICE==319
 integer, parameter, public :: links=18
#endif

 integer, parameter, public :: nbcdir=6
 
 !max
 integer, save, public :: minx, maxx, miny, maxy, minz, maxz
 integer, save, public :: wminx, wmaxx, wminy, wmaxy, wminz, wmaxz
 integer, save, public :: ixpbc, iypbc, izpbc
 integer, save, public :: ix,iy,iz,mynx,myny,mynz
 type(REALPTR), dimension(0:links), public, protected :: aoptpR,aoptpB
 !max
 integer, save, protected, public :: LBintegrator=0
 
 integer, save, protected, public :: idistselect=0
 
 integer, save, protected, public :: ibctype=0
 
 integer, save, protected, public :: bc_type_east=0
 integer, save, protected, public :: bc_type_west=0
 integer, save, protected, public :: bc_type_rear=0
 integer, save, protected, public :: bc_type_front=0
 integer, save, protected, public :: bc_type_north=0
 integer, save, protected, public :: bc_type_south=0
 
 !timestep interval to rescale the total fluid mass
 integer, save, protected, public :: imass_rescale=100
 
 !selected wetting mode 0=average 1=fixed value (Benzi's style!)
 integer, save, protected, public :: iselwetting=0
 
 !interpolation mode of omega (0=phase field interp;1=Grunau interp)
 integer, parameter :: iselomegaint=1
 
 integer, save, protected, public, allocatable, dimension(:) :: &
  typeobjectliq
 
 integer, save :: nvar_managebc=0
 type(REALPTR_S), allocatable, dimension(:,:) :: rhoR_managebc,psiR_managebc
 type(REALPTR_S), allocatable, dimension(:,:) :: rhoB_managebc,psiB_managebc
 type(REALPTR_S), allocatable, dimension(:,:) :: u_managebc
 type(REALPTR_S), allocatable, dimension(:,:) :: v_managebc
 type(REALPTR_S), allocatable, dimension(:,:) :: w_managebc
 
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
 
 !is there a isfluid.dat file to read? 
 logical, save, protected, public :: lread_isfluid=.false.
 
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
 
 real(kind=PRC), save, protected, public :: bc_rhoR_east = ZERO
 real(kind=PRC), save, protected, public :: bc_rhoR_west = ZERO
 real(kind=PRC), save, protected, public :: bc_rhoR_rear= ZERO
 real(kind=PRC), save, protected, public :: bc_rhoR_front = ZERO
 real(kind=PRC), save, protected, public :: bc_rhoR_north= ZERO
 real(kind=PRC), save, protected, public :: bc_rhoR_south= ZERO
 
 real(kind=PRC), save, protected, public :: bc_rhoB_east = ZERO
 real(kind=PRC), save, protected, public :: bc_rhoB_west = ZERO
 real(kind=PRC), save, protected, public :: bc_rhoB_rear= ZERO
 real(kind=PRC), save, protected, public :: bc_rhoB_front = ZERO
 real(kind=PRC), save, protected, public :: bc_rhoB_north= ZERO
 real(kind=PRC), save, protected, public :: bc_rhoB_south= ZERO
 
 real(kind=PRC), save, protected, public :: bc_flow_east = ZERO
 real(kind=PRC), save, protected, public :: bc_flow_west = ZERO
 real(kind=PRC), save, protected, public :: bc_flow_rear= ZERO
 real(kind=PRC), save, protected, public :: bc_flow_front = ZERO
 real(kind=PRC), save, protected, public :: bc_flow_north= ZERO
 real(kind=PRC), save, protected, public :: bc_flow_south= ZERO
 
 real(kind=PRC), save, protected, public :: bc_u_east = ZERO
 real(kind=PRC), save, protected, public :: bc_u_west = ZERO
 real(kind=PRC), save, protected, public :: bc_u_rear= ZERO
 real(kind=PRC), save, protected, public :: bc_u_front = ZERO
 real(kind=PRC), save, protected, public :: bc_u_north= ZERO
 real(kind=PRC), save, protected, public :: bc_u_south= ZERO
 
 real(kind=PRC), save, protected, public :: bc_v_east = ZERO
 real(kind=PRC), save, protected, public :: bc_v_west = ZERO
 real(kind=PRC), save, protected, public :: bc_v_rear= ZERO
 real(kind=PRC), save, protected, public :: bc_v_front = ZERO
 real(kind=PRC), save, protected, public :: bc_v_north= ZERO
 real(kind=PRC), save, protected, public :: bc_v_south= ZERO
 
 real(kind=PRC), save, protected, public :: bc_w_east = ZERO
 real(kind=PRC), save, protected, public :: bc_w_west = ZERO
 real(kind=PRC), save, protected, public :: bc_w_rear= ZERO
 real(kind=PRC), save, protected, public :: bc_w_front = ZERO
 real(kind=PRC), save, protected, public :: bc_w_north= ZERO
 real(kind=PRC), save, protected, public :: bc_w_south= ZERO
 
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
 
 !fixed density to set the wall wetting as in Benzi's style
 real(kind=PRC), save, protected, public :: densR_wetting = ZERO
 real(kind=PRC), save, protected, public :: densB_wetting = ZERO
 
 integer(kind=1), save, protected, public, allocatable, dimension(:,:,:) :: isfluid
 integer(kind=1), save, protected, public, allocatable, dimension(:,:,:) :: new_isfluid
 integer(kind=1), save, protected, public, allocatable, dimension(:,:,:) :: bcfluid
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
 
 real(kind=PRC), save, protected, public, allocatable, dimension(:) :: bc_rhoR
 real(kind=PRC), save, protected, public, allocatable, dimension(:) :: bc_rhoB
 real(kind=PRC), save, protected, public, allocatable, dimension(:) :: bc_u
 real(kind=PRC), save, protected, public, allocatable, dimension(:) :: bc_v
 real(kind=PRC), save, protected, public, allocatable, dimension(:) :: bc_w
 real(kind=PRC), save, protected, public, allocatable, dimension(:) :: bc_flow
 real(kind=PRC), save, protected, public, allocatable, dimension(:) :: bc_omega
 
 real(kind=PRC), save, protected, public, allocatable, dimension(:,:) :: &
  objectdata
 
!real denoting the LOCAL fluid mass at actual time for red and blue
!(it changes if particles are presents)
 real(kind=PRC), save :: nfluidmassR = ZERO
 real(kind=PRC), save :: nfluidmassB = ZERO
!real denoting the TOTAL fluid mass at initial time for red and blue
 real(kind=PRC), save :: totnfluidmassR_init = ZERO
 real(kind=PRC), save :: totnfluidmassB_init = ZERO
 
!array containing the node list of the second belt 
 integer, save :: ndouble
 integer, save, allocatable, dimension(:) :: exdouble
 integer, save, allocatable, dimension(:) :: eydouble
 integer, save, allocatable, dimension(:) :: ezdouble
 
 integer, allocatable, save :: ibounce(:,:)
 integer, save :: nbounce0,nbounce6,nbounce7,nbounce8,nbounce9
 integer, dimension(0:nbcdir), save :: nbounce6dir,nbounce7dir, &
  nbounce8dir,nbounce9dir
 
 integer, allocatable, save :: isguards(:,:)
 integer, save :: nguards=0
 
 integer, allocatable, save :: isinglehalo(:,:)
 integer, save :: nsinglehalo=0

 integer(kind=1), allocatable, save :: ipoplistbc(:)
 integer, allocatable, save :: poplistbc(:,:)
 integer, save :: npoplistbc=0

 
#if LATTICE==319
 
 character(len=6), parameter, public :: latt_name="d3q19 "
 
 integer, parameter, public :: latt_dim=3
 
 integer, parameter, public :: nbuff=2
 
 real(kind=PRC), parameter, public :: cssq = ( ONE / THREE )
 
 real(kind=PRC), parameter :: p0 = ( ONE / THREE )
 real(kind=PRC), parameter :: p1 = ( ONE / EIGHTEEN )
 real(kind=PRC), parameter :: p2 = ( ONE / THIRTYSIX )
 real(kind=PRC), dimension(0:links), parameter, public :: &
  p = (/p0,p1,p1,p1,p1,p1,p1,p2,p2,p2,p2,p2,p2,p2,p2,p2,p2,p2,p2/)
 
 real(kind=PRC), parameter :: b_xi = real(2.d0,kind=PRC)
 real(kind=PRC), parameter :: b0 = -(TWO+TWO*b_xi)/(THREE*b_xi+TWELVE)
 real(kind=PRC), parameter :: b1 = (b_xi)/(SIX*b_xi+TWENTYFOUR)
 real(kind=PRC), parameter :: b2 = ONE/(SIX*b_xi+TWENTYFOUR)
 real(kind=PRC), dimension(0:links), parameter, public :: &
  b_l = (/b0,b1,b1,b1,b1,b1,b1,b2,b2,b2,b2,b2,b2,b2,b2,b2,b2,b2,b2/) 
  
 real(kind=PRC), parameter :: psi0 = -(FIVE/TWO)
 real(kind=PRC), parameter :: psi1 = -(ONE/SIX)
 real(kind=PRC), parameter :: psi2 = (ONE/TWENTYFOUR)
 real(kind=PRC), dimension(0:links), parameter, public :: &
  psi_CG = (/psi0,psi1,psi1,psi1,psi1,psi1,psi1,psi2,psi2,psi2,psi2, &
  psi2,psi2,psi2,psi2,psi2,psi2,psi2,psi2/) 
  
 real(kind=PRC), parameter :: xi0 = ZERO
 real(kind=PRC), parameter :: xi1 = (ONE/FOUR)
 real(kind=PRC), parameter :: xi2 = (ONE/EIGHT)
 real(kind=PRC), dimension(0:links), parameter, public :: &
  xi_CG = (/xi0,xi1,xi1,xi1,xi1,xi1,xi1,xi2,xi2,xi2,xi2,xi2,xi2,xi2, &
  xi2,xi2,xi2,xi2,xi2/) 
 
 real(kind=PRC), parameter :: phi0 = ZERO
 real(kind=PRC), parameter :: phi1 = (ONE/TWELVE)
 real(kind=PRC), parameter :: phi2 = (ONE/TWENTYFOUR)
 real(kind=PRC), dimension(0:links), parameter, public :: &
  phi_CG = (/phi0,phi1,phi1,phi1,phi1,phi1,phi1,phi2,phi2,phi2,phi2, &
  phi2,phi2,phi2,phi2,phi2,phi2,phi2,phi2/) 
 
 real(kind=PRC), parameter :: varphi0 = ONE
 real(kind=PRC), parameter :: varphi1 = -(ONE/TWELVE)
 real(kind=PRC), parameter :: varphi2 = -(ONE/TWENTYFOUR)
 real(kind=PRC), dimension(0:links), parameter, public :: &
  varphi_CG = (/varphi0,varphi1,varphi1,varphi1,varphi1,varphi1, &
  varphi1,varphi2,varphi2,varphi2,varphi2,varphi2,varphi2,varphi2, &
  varphi2,varphi2,varphi2,varphi2,varphi2/) 
  
 real(kind=PRC), dimension(0:links), parameter, public :: &
  a = (/ZERO,p1/cssq,p1/cssq,p1/cssq,p1/cssq,p1/cssq,p1/cssq,p2/cssq, &
  p2/cssq,p2/cssq,p2/cssq,p2/cssq,p2/cssq,p2/cssq,p2/cssq,p2/cssq, &
  p2/cssq,p2/cssq,p2/cssq/)
 
 !lattice vectors
 integer, dimension(0:links), parameter, public :: &
   !      0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18
  ex = (/ 0, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1, 0, 0, 0, 0/)
 integer, dimension(0:links), parameter, public :: &
  ey = (/ 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1/)
 integer, dimension(0:links), parameter, public :: &
  ez = (/ 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1/)
 integer, dimension(0:links), parameter, public :: &
  opp =(/ 0, 2, 1, 4, 3, 6, 5, 8, 7,10, 9,12,11,14,13,16,15,18,17/)
  
 real(kind=PRC), dimension(0:links), parameter, public :: &
  dex = real(ex,kind=PRC)
 real(kind=PRC), dimension(0:links), parameter, public :: &
  dey = real(ey,kind=PRC)
 real(kind=PRC), dimension(0:links), parameter, public :: &
  dez = real(ez,kind=PRC)
 
#endif 

! coeff for Colour Gradient
 real(kind=PRC), dimension(0:links), protected, public :: phiR_CG
 real(kind=PRC), dimension(0:links), protected, public :: phiB_CG

 integer, parameter, public :: linksd3q27=26
 
 !lattice vectors
 integer, dimension(0:linksd3q27), parameter, public :: &
   !           0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26
  exd3q27 = (/ 0, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1, 1,-1/)
 integer, dimension(0:linksd3q27), parameter, public :: &
  eyd3q27 = (/ 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1, 1,-1,-1, 1, 1,-1/)
 integer, dimension(0:linksd3q27), parameter, public :: &
  ezd3q27 = (/ 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1,-1, 1/)
 
 real(kind=PRC), dimension(0:linksd3q27), parameter, public :: &
  dexd3q27 = real(exd3q27,kind=PRC)
 real(kind=PRC), dimension(0:linksd3q27), parameter, public :: &
  deyd3q27 = real(eyd3q27,kind=PRC)
 real(kind=PRC), dimension(0:linksd3q27), parameter, public :: &
  dezd3q27 = real(ezd3q27,kind=PRC)
 
 real(kind=PRC), parameter :: p0d3q27 = ( TWO / THREE )**THREE
 real(kind=PRC), parameter :: p1d3q27 = (( TWO / THREE )**TWO ) * ( ONE/ SIX )
 real(kind=PRC), parameter :: p2d3q27 = (( ONE / SIX )**TWO) * ( TWO/ THREE)
 real(kind=PRC), parameter :: p3d3q27 = ( ONE / SIX )**THREE
 real(kind=PRC), dimension(0:linksd3q27), parameter, public :: &
  pd3q27 = (/p0d3q27,p1d3q27,p1d3q27,p1d3q27,p1d3q27,p1d3q27,p1d3q27, &
  p2d3q27,p2d3q27,p2d3q27,p2d3q27,p2d3q27,p2d3q27,p2d3q27,p2d3q27, &
  p2d3q27,p2d3q27,p2d3q27,p2d3q27,p3d3q27,p3d3q27,p3d3q27,p3d3q27, &
  p3d3q27,p3d3q27,p3d3q27,p3d3q27/)
  
 real(kind=PRC), dimension(0:linksd3q27), parameter, public :: &
  ad3q27 = (/ZERO,p1d3q27/cssq,p1d3q27/cssq,p1d3q27/cssq, &
  p1d3q27/cssq,p1d3q27/cssq,p1d3q27/cssq, &
  p2d3q27/cssq,p2d3q27/cssq,p2d3q27/cssq,p2d3q27/cssq,p2d3q27/cssq, &
  p2d3q27/cssq,p2d3q27/cssq,p2d3q27/cssq, &
  p2d3q27/cssq,p2d3q27/cssq,p2d3q27/cssq,p2d3q27/cssq,p3d3q27/cssq, &
  p3d3q27/cssq,p3d3q27/cssq,p3d3q27/cssq, &
  p3d3q27/cssq,p3d3q27/cssq,p3d3q27/cssq,p3d3q27/cssq/)
 
 real(kind=PRC), save, public, allocatable, target, dimension(:,:,:) :: & 
  f00R,f01R,f02R,f03R,f04R,f05R,f06R,f07R,f08R
 real(kind=PRC), save, public, allocatable, target, dimension(:,:,:) :: &
  f00B,f01B,f02B,f03B,f04B,f05B,f06B,f07B,f08B
 
 real(kind=PRC), save, public, allocatable, target, dimension(:,:,:) :: & 
  f09R,f10R,f11R,f12R,f13R,f14R,f15R,f16R,f17R,f18R
 real(kind=PRC), save, public, allocatable, target, dimension(:,:,:) :: &
  f09B,f10B,f11B,f12B,f13B,f14B,f15B,f16B,f17B,f18B


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
 public :: set_boundary_conditions_type
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
 public :: compute_densities_wall
 public :: compute_psi_sc
 public :: set_lsingle_fluid
 public :: set_lread_isfluid
 public :: driver_apply_bounceback_pop
 public :: probe_red_moments_in_node
 public :: probe_blue_moments_in_node
 public :: set_value_bc_east
 public :: set_value_bc_west
 public :: set_value_bc_rear
 public :: set_value_bc_front
 public :: set_value_bc_north
 public :: set_value_bc_south
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
 public :: particle_bounce_back, particle_bounce_back_phase2
 !public :: particle_moving_fluids
 public :: initialize_new_isfluid
 public :: update_isfluid
 public :: driver_bc_isfluid, driver_bc_new_isfluid
 public :: mapping_new_isfluid
 public :: particle_delete_fluids
 public :: particle_create_fluids
 public :: erase_fluids_in_particles
 public :: print_all_pops
 public :: print_all_pops_center
 public :: print_all_pops_area_shpere
 public :: pimage
 public :: omega_to_viscosity
 public :: compute_sc_particle_interact_phase1
 public :: compute_sc_particle_interact_phase2
 public :: print_all_pops2
 public :: driver_bc_pops
 public :: dumpHvar, restorePops, restoreHvar
 public :: driver_bc_psi
 public :: fixPops
 public :: set_cap_force
 public :: set_mass_rescale
 public :: rescale_fluid_mass
 public :: set_fluid_wall_mode
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
     !if(mod(nx,2)==0)mynx=mynx+1
     myny=maxy+nbuff
     mynz=maxz+nbuff

     istat=0
     
     allocate(isfluid(ix:mynx,iy:myny,iz:mynz),stat=istat(1))
     
     if(lparticles)then
       allocate(new_isfluid(ix:mynx,iy:myny,iz:mynz),stat=istat(98))
     endif
     
     allocate(bcfluid(ix:mynx,iy:myny,iz:mynz),stat=istat(99))
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

     allocate(f00R(ix:mynx,iy:myny,iz:mynz),stat=istat(30))
     allocate(f01R(ix:mynx,iy:myny,iz:mynz),stat=istat(31))
     allocate(f02R(ix:mynx,iy:myny,iz:mynz),stat=istat(32))
     allocate(f03R(ix:mynx,iy:myny,iz:mynz),stat=istat(33))
     allocate(f04R(ix:mynx,iy:myny,iz:mynz),stat=istat(34))
     allocate(f05R(ix:mynx,iy:myny,iz:mynz),stat=istat(35))
     allocate(f06R(ix:mynx,iy:myny,iz:mynz),stat=istat(36))
     allocate(f07R(ix:mynx,iy:myny,iz:mynz),stat=istat(37))
     allocate(f08R(ix:mynx,iy:myny,iz:mynz),stat=istat(38))
     allocate(f09R(ix:mynx,iy:myny,iz:mynz),stat=istat(39))
     allocate(f10R(ix:mynx,iy:myny,iz:mynz),stat=istat(40))
     allocate(f11R(ix:mynx,iy:myny,iz:mynz),stat=istat(41))
     allocate(f12R(ix:mynx,iy:myny,iz:mynz),stat=istat(42))
     allocate(f13R(ix:mynx,iy:myny,iz:mynz),stat=istat(43))
     allocate(f14R(ix:mynx,iy:myny,iz:mynz),stat=istat(44))
     allocate(f15R(ix:mynx,iy:myny,iz:mynz),stat=istat(45))
     allocate(f16R(ix:mynx,iy:myny,iz:mynz),stat=istat(46))
     allocate(f17R(ix:mynx,iy:myny,iz:mynz),stat=istat(47))
     allocate(f18R(ix:mynx,iy:myny,iz:mynz),stat=istat(48))

     aoptpR(0)%p => f00R
     aoptpR(1)%p => f01R
     aoptpR(2)%p => f02R
     aoptpR(3)%p => f03R
     aoptpR(4)%p => f04R
     aoptpR(5)%p => f05R
     aoptpR(6)%p => f06R
     aoptpR(7)%p => f07R
     aoptpR(8)%p => f08R
     aoptpR(9)%p => f09R
     aoptpR(10)%p => f10R
     aoptpR(11)%p => f11R
     aoptpR(12)%p => f12R
     aoptpR(13)%p => f13R
     aoptpR(14)%p => f14R
     aoptpR(15)%p => f15R
     aoptpR(16)%p => f16R
     aoptpR(17)%p => f17R
     aoptpR(18)%p => f18R
     
     

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

        allocate(f00B(ix:mynx,iy:myny,iz:mynz),stat=istat(60))
        allocate(f01B(ix:mynx,iy:myny,iz:mynz),stat=istat(61))
        allocate(f02B(ix:mynx,iy:myny,iz:mynz),stat=istat(62))
        allocate(f03B(ix:mynx,iy:myny,iz:mynz),stat=istat(63))
        allocate(f04B(ix:mynx,iy:myny,iz:mynz),stat=istat(64))
        allocate(f05B(ix:mynx,iy:myny,iz:mynz),stat=istat(65))
        allocate(f06B(ix:mynx,iy:myny,iz:mynz),stat=istat(66))
        allocate(f07B(ix:mynx,iy:myny,iz:mynz),stat=istat(67))
        allocate(f08B(ix:mynx,iy:myny,iz:mynz),stat=istat(68))
        allocate(f09B(ix:mynx,iy:myny,iz:mynz),stat=istat(69))
        allocate(f10B(ix:mynx,iy:myny,iz:mynz),stat=istat(60))
        allocate(f11B(ix:mynx,iy:myny,iz:mynz),stat=istat(61))
        allocate(f12B(ix:mynx,iy:myny,iz:mynz),stat=istat(62))
        allocate(f13B(ix:mynx,iy:myny,iz:mynz),stat=istat(63))
        allocate(f14B(ix:mynx,iy:myny,iz:mynz),stat=istat(64))
        allocate(f15B(ix:mynx,iy:myny,iz:mynz),stat=istat(65))
        allocate(f16B(ix:mynx,iy:myny,iz:mynz),stat=istat(66))
        allocate(f17B(ix:mynx,iy:myny,iz:mynz),stat=istat(67))
        allocate(f18B(ix:mynx,iy:myny,iz:mynz),stat=istat(68))

        aoptpB(0)%p => f00B
        aoptpB(1)%p => f01B
        aoptpB(2)%p => f02B
        aoptpB(3)%p => f03B
        aoptpB(4)%p => f04B
        aoptpB(5)%p => f05B
        aoptpB(6)%p => f06B
        aoptpB(7)%p => f07B
        aoptpB(8)%p => f08B
        aoptpB(9)%p => f09B
        aoptpB(10)%p => f10B
        aoptpB(11)%p => f11B
        aoptpB(12)%p => f12B
        aoptpB(13)%p => f13B
        aoptpB(14)%p => f14B
        aoptpB(15)%p => f15B
        aoptpB(16)%p => f16B
        aoptpB(17)%p => f17B
        aoptpB(18)%p => f18B
        
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
   
   isfluid(ix:mynx,iy:myny,iz:mynz)=3
   bcfluid(ix:mynx,iy:myny,iz:mynz)=0
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

   f00R(ix:mynx,iy:myny,iz:mynz)=ZERO
   f01R(ix:mynx,iy:myny,iz:mynz)=ZERO
   f02R(ix:mynx,iy:myny,iz:mynz)=ZERO
   f03R(ix:mynx,iy:myny,iz:mynz)=ZERO
   f04R(ix:mynx,iy:myny,iz:mynz)=ZERO
   f05R(ix:mynx,iy:myny,iz:mynz)=ZERO
   f06R(ix:mynx,iy:myny,iz:mynz)=ZERO
   f07R(ix:mynx,iy:myny,iz:mynz)=ZERO
   f08R(ix:mynx,iy:myny,iz:mynz)=ZERO
   f09R(ix:mynx,iy:myny,iz:mynz)=ZERO
   f10R(ix:mynx,iy:myny,iz:mynz)=ZERO
   f11R(ix:mynx,iy:myny,iz:mynz)=ZERO
   f12R(ix:mynx,iy:myny,iz:mynz)=ZERO
   f13R(ix:mynx,iy:myny,iz:mynz)=ZERO
   f14R(ix:mynx,iy:myny,iz:mynz)=ZERO
   f15R(ix:mynx,iy:myny,iz:mynz)=ZERO
   f16R(ix:mynx,iy:myny,iz:mynz)=ZERO
   f17R(ix:mynx,iy:myny,iz:mynz)=ZERO
   f18R(ix:mynx,iy:myny,iz:mynz)=ZERO
   
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

      f00B(ix:mynx,iy:myny,iz:mynz)=ZERO
      f01B(ix:mynx,iy:myny,iz:mynz)=ZERO
      f02B(ix:mynx,iy:myny,iz:mynz)=ZERO
      f03B(ix:mynx,iy:myny,iz:mynz)=ZERO
      f04B(ix:mynx,iy:myny,iz:mynz)=ZERO
      f05B(ix:mynx,iy:myny,iz:mynz)=ZERO
      f06B(ix:mynx,iy:myny,iz:mynz)=ZERO
      f07B(ix:mynx,iy:myny,iz:mynz)=ZERO
      f08B(ix:mynx,iy:myny,iz:mynz)=ZERO
      f09B(ix:mynx,iy:myny,iz:mynz)=ZERO
      f10B(ix:mynx,iy:myny,iz:mynz)=ZERO
      f11B(ix:mynx,iy:myny,iz:mynz)=ZERO
      f12B(ix:mynx,iy:myny,iz:mynz)=ZERO
      f13B(ix:mynx,iy:myny,iz:mynz)=ZERO
      f14B(ix:mynx,iy:myny,iz:mynz)=ZERO
      f15B(ix:mynx,iy:myny,iz:mynz)=ZERO
      f16B(ix:mynx,iy:myny,iz:mynz)=ZERO
      f17B(ix:mynx,iy:myny,iz:mynz)=ZERO
      f18B(ix:mynx,iy:myny,iz:mynz)=ZERO
   endif
   
   l=0
    do k=-2,2
      do j=-2,2
        do i=-2,2
          if(i==-2 .or. i==2 .or. j==-2 .or. j==2 .or. &
            k==-2 .or. k==2)then
             l=l+1
          endif
        enddo
      enddo
    enddo
    ndouble=l
    allocate(exdouble(ndouble),eydouble(ndouble),ezdouble(ndouble))
    l=0
    do k=-2,2
      do j=-2,2
        do i=-2,2
          if(i==-2 .or. i==2 .or. j==-2 .or. j==2 .or. &
            k==-2 .or. k==2)then
             l=l+1
             exdouble(l)=i
             eydouble(l)=j
             ezdouble(l)=k
          endif
        enddo
      enddo
    enddo

 end subroutine allocate_fluids
 
 subroutine set_boundary_conditions_type(itemp1,itemp2,itemp3)
 
  implicit none
  
  integer, intent(in) :: itemp1,itemp2,itemp3
  integer ::i,j,k
  integer, dimension(0:1,0:1,0:1), parameter :: iselbct = &
   reshape((/0,1,2,3,4,5,6,7/),(/2,2,2/))
 ! 0=F    1=T  periodic
 !     itemp1      itemp2      itemp3     ibctype
 !          0           0           0           0
 !          1           0           0           1
 !          0           1           0           2
 !          1           1           0           3
 !          0           0           1           4
 !          1           0           1           5
 !          0           1           1           6
 !          1           1           1           7
  
  
  if( itemp1 .ne. 0 .and. itemp1 .ne. 1)then
    call error(9)
  endif
  
  if( itemp2 .ne. 0 .and. itemp2 .ne. 1)then
    call error(9)
  endif
  
  if( itemp3 .ne. 0 .and. itemp3 .ne. 1)then
    call error(9)
  endif
  
  ibctype=iselbct(itemp1,itemp2,itemp3)
  
  return
  
 end subroutine set_boundary_conditions_type
 
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
  
 
  call initialiaze_manage_bc_hvar_selfcomm
  
  call initialiaze_manage_bc_pop_selfcomm

  
  return
  
 end subroutine driver_initialiaze_manage_bc_selfcomm
 
 subroutine initialize_isfluid_bcfluid(lvtkfilesub, lparticles)
 
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
  
  integer :: i,j,k,l,idir,ishift,jshift,kshift,ll
  integer :: itemp,jtemp,ktemp
  
  logical :: ltest(1)
  integer, parameter :: nistatmax=10
  integer, dimension(nistatmax) :: istat
  
  integer, parameter, dimension(3) :: myc=(/0,0,8/)
  real(8) :: dsum1,isum
 
  isfluid(:,:,:)=3
  bcfluid(:,:,:)=0
! set isfluid as you like (in future to be given as input file)
! all fluid
  
  forall(i=minx:maxx,j=miny:maxy,k=minz:maxz)isfluid(i,j,k)=1
  
! read isfluid.dat if necessary  
  call read_isfluid_dat(0)
  
  do k=minz-nbuff,maxz+nbuff
    do j=miny-nbuff,maxy+nbuff
      do i=minx-nbuff,maxx+nbuff
        if(i==0 .or. j==0 .or. k==0 .or. i==(nx+1) .or. j==(ny+1) .or. k==(nz+1))then
          if(ibctype==0)then ! 0 F F F
            if(i==(nx+1) .and. j==0 .and. k==(nz+1))then !north east rear corner 
              isfluid(i,j,k)=0
            elseif(i==(nx+1) .and. j==(ny+1) .and. k==(nz+1))then !north east front corner
              isfluid(i,j,k)=0
            elseif(i==0 .and. j==(ny+1) .and. k==(nz+1))then !north west front corner
              isfluid(i,j,k)=0
            elseif(i==0 .and. j==0 .and. k==(nz+1))then !north west rear corner
              isfluid(i,j,k)=0
            elseif(i==(nx+1) .and. j==0 .and. k==0)then !south east rear corner
              isfluid(i,j,k)=0
            elseif(i==0 .and. j==0 .and. k==0)then !south west rear corner
              isfluid(i,j,k)=0
            elseif(i==0 .and. j==(ny+1) .and. k==0)then !south west front corner
              isfluid(i,j,k)=0
            elseif(i==(nx+1) .and. j==(ny+1) .and. k==0)then !south east front corner
              isfluid(i,j,k)=0
            elseif(i==(nx+1) .and. j==0)then !rear east edge
              isfluid(i,j,k)=0
            elseif(i==0 .and. j==0)then !rear west edge
              isfluid(i,j,k)=0
            elseif(i==(nx+1) .and. k==(nz+1))then !north east edge
              isfluid(i,j,k)=0
            elseif(j==0 .and. k==(nz+1))then !north rear edge
              isfluid(i,j,k)=0
            elseif(j==(ny+1) .and. k==(nz+1))then !north front edge
              isfluid(i,j,k)=0
            elseif(i==0 .and. k==(nz+1))then !north west edge
              isfluid(i,j,k)=0
            elseif(i==(nx+1) .and. j==(ny+1))then !front east edge
              isfluid(i,j,k)=0
            elseif(i==0 .and. j==(ny+1))then !front west edge
              isfluid(i,j,k)=0
            elseif(i==(nx+1) .and. k==0)then !south east edge
              isfluid(i,j,k)=0
            elseif(j==0 .and. k==0)then !south east edge
              isfluid(i,j,k)=0
            elseif(j==(ny+1) .and. k==0)then !south front edge
              isfluid(i,j,k)=0
            elseif(i==0 .and. k==0)then !south west edge
              isfluid(i,j,k)=0
            elseif(i==(nx+1))then !east side
              if(bc_type_east/=0)then
                if(isfluid(i+ex(2),j,k)==1)then
                  isfluid(i,j,k)=bc_type_east+5
                  bcfluid(i,j,k)=2
                else
                  isfluid(i,j,k)=0
                endif
              else
                isfluid(i,j,k)=0
              endif
            elseif(i==0)then !west side
              if(bc_type_west/=0)then
                if(isfluid(i+ex(1),j,k)==1)then
                  isfluid(i,j,k)=bc_type_west+5
                  bcfluid(i,j,k)=1
                else
                  isfluid(i,j,k)=0
                endif
              else
                isfluid(i,j,k)=0
              endif
            elseif(k==(nz+1))then !north side
              if(bc_type_north/=0)then
                if(isfluid(i,j,k+ez(6))==1)then
                  isfluid(i,j,k)=bc_type_north+5
                  bcfluid(i,j,k)=6
                else
                  isfluid(i,j,k)=0
                endif
              else
                isfluid(i,j,k)=0
              endif
            elseif(k==0)then !south side
              if(bc_type_south/=0)then
                if(isfluid(i,j,k+ez(5))==1)then
                  isfluid(i,j,k)= bc_type_south+5
                  bcfluid(i,j,k)=5
                else
                  isfluid(i,j,k)=0
                endif
              else
                isfluid(i,j,k)=0
              endif
            elseif(j==0)then !rear side
              if(bc_type_rear/=0)then
                if(isfluid(i,j+ez(3),k)==1)then
                  isfluid(i,j,k)= bc_type_rear+5
                  bcfluid(i,j,k)=3
                else
                  isfluid(i,j,k)=0
                endif
              else
                isfluid(i,j,k)=0
              endif
            elseif(j==(ny+1))then !front side
              if(bc_type_front/=0)then
                if(isfluid(i,j+ez(4),k)==1)then
                  isfluid(i,j,k)= bc_type_front+5
                  bcfluid(i,j,k)=4
                else
                  isfluid(i,j,k)=0
                endif
              else
                isfluid(i,j,k)=0
              endif
            else
              isfluid(i,j,k)=0
            endif
          elseif(ibctype==1)then ! 1 T F F
            if(j==0 .or. j==(ny+1) .or. k==0 .or. k==(nz+1))then
              if(k==(nz+1) .and. j==0)then !north rear edge
                isfluid(i,j,k)=0
              elseif(k==(nz+1) .and. j==(ny+1))then !north front edge
                isfluid(i,j,k)=0
              elseif(k==0 .and. j==0)then !south rear edge
                isfluid(i,j,k)=0
              elseif(k==0 .and. j==(ny+1))then !south front edge
                isfluid(i,j,k)=0
              elseif(k==(nz+1))then !north side
                if(bc_type_north/=0)then
                  if(isfluid(i,j,k+ez(6))==1)then
                    isfluid(i,j,k)=bc_type_north+5
                    bcfluid(i,j,k)=6
                  else
                    isfluid(i,j,k)=0
                  endif
                else
                  isfluid(i,j,k)=0
                endif
              elseif(k==0)then !south side
                if(bc_type_south/=0)then
                  if(isfluid(i,j,k+ez(5))==1)then
                    isfluid(i,j,k)= bc_type_south+5
                    bcfluid(i,j,k)=5
                  else
                    isfluid(i,j,k)=0
                  endif
                else
                  isfluid(i,j,k)=0
                endif
              elseif(j==0)then !rear side
                if(bc_type_rear/=0)then
                  if(isfluid(i,j+ez(3),k)==1)then
                    isfluid(i,j,k)= bc_type_rear+5
                    bcfluid(i,j,k)=3
                  else
                    isfluid(i,j,k)=0
                  endif
                else
                  isfluid(i,j,k)=0
                endif
              elseif(j==(ny+1))then !front side
                if(bc_type_front/=0)then
                  if(isfluid(i,j+ez(4),k)==1)then
                    isfluid(i,j,k)= bc_type_front+5
                    bcfluid(i,j,k)=4
                  else
                    isfluid(i,j,k)=0
                  endif
                else
                  isfluid(i,j,k)=0
                endif
              else
                isfluid(i,j,k)=0
              endif
            endif
          elseif(ibctype==2)then ! 2 F T F
            if(i==0 .or. i==(nx+1) .or. k==0 .or. k==(nz+1))then
              if(i==(nx+1) .and. k==(nz+1))then !east north edge
                isfluid(i,j,k)=0
              elseif(i==(nx+1) .and. k==0)then !east south edge
                isfluid(i,j,k)=0
              elseif(i==0 .and. k==(nz+1))then !west north edge
                isfluid(i,j,k)=0
              elseif(i==0 .and. k==0)then !west south edge
                isfluid(i,j,k)=0
              elseif(i==(nx+1))then !east side
                if(bc_type_east/=0)then
                  if(isfluid(i+ex(2),j,k)==1)then
                    isfluid(i,j,k)=bc_type_east+5
                    bcfluid(i,j,k)=2
                  else
                    isfluid(i,j,k)=0
                  endif
                else
                  isfluid(i,j,k)=0
                endif
              elseif(i==0)then !west side
                if(bc_type_west/=0)then
                  if(isfluid(i+ex(1),j,k)==1)then
                    isfluid(i,j,k)=bc_type_west+5
                    bcfluid(i,j,k)=1
                  else
                    isfluid(i,j,k)=0
                  endif
                else
                  isfluid(i,j,k)=0
                endif
              elseif(k==(nz+1))then !north side
                if(bc_type_north/=0)then
                  if(isfluid(i,j,k+ez(6))==1)then
                    isfluid(i,j,k)=bc_type_north+5
                    bcfluid(i,j,k)=6
                  else
                    isfluid(i,j,k)=0
                  endif
                else
                  isfluid(i,j,k)=0
                endif
              elseif(k==0)then !south side
                if(bc_type_south/=0)then
                  if(isfluid(i,j,k+ez(5))==1)then
                    isfluid(i,j,k)= bc_type_south+5
                    bcfluid(i,j,k)=5
                  else
                    isfluid(i,j,k)=0
                  endif
                else
                  isfluid(i,j,k)=0
                endif
              else
                isfluid(i,j,k)=0
              endif
            endif
          elseif(ibctype==3)then ! 3 T T F
            if(k==0 .or. k==(nz+1))then
              if(k==(nz+1))then !north side
                if(bc_type_north/=0)then
                  if(isfluid(i,j,k+ez(6))==1)then
                    isfluid(i,j,k)=bc_type_north+5
                    bcfluid(i,j,k)=6
                  else
                    isfluid(i,j,k)=0
                  endif
                else
                  isfluid(i,j,k)=0
                endif
              elseif(k==0)then !south side
                if(bc_type_south/=0)then
                  if(isfluid(i,j,k+ez(5))==1)then
                    isfluid(i,j,k)= bc_type_south+5
                    bcfluid(i,j,k)=5
                  else
                    isfluid(i,j,k)=0
                  endif
                else
                  isfluid(i,j,k)=0
                endif
              else
                isfluid(i,j,k)=0
              endif
            endif
          elseif(ibctype==4)then ! 4 F F T
            if(i==0 .or. i==(nx+1) .or. j==0 .or. j==(ny+1))then
              if(i==(nx+1) .and. j==0)then !east rear edge
                isfluid(i,j,k)=0
              elseif(i==(nx+1) .and. j==(ny+1))then !east front edge
                isfluid(i,j,k)=0
              elseif(i==0 .and. j==0)then !west rear edge
                isfluid(i,j,k)=0
              elseif(i==0 .and. j==(ny+1))then !west front edge
                isfluid(i,j,k)=0
              elseif(i==(nx+1))then !east side
                if(bc_type_east/=0)then
                  if(isfluid(i+ex(2),j,k)==1)then
                    isfluid(i,j,k)=bc_type_east+5
                    bcfluid(i,j,k)=2
                  else
                    isfluid(i,j,k)=0
                  endif
                else
                  isfluid(i,j,k)=0
                endif
              elseif(i==0)then !west side
                if(bc_type_west/=0)then
                  if(isfluid(i+ex(1),j,k)==1)then
                    isfluid(i,j,k)=bc_type_west+5
                    bcfluid(i,j,k)=1
                  else
                    isfluid(i,j,k)=0
                  endif
                else
                  isfluid(i,j,k)=0
                endif
              elseif(j==0)then !rear side
                if(bc_type_rear/=0)then
                  if(isfluid(i,j+ez(3),k)==1)then
                    isfluid(i,j,k)= bc_type_rear+5
                    bcfluid(i,j,k)=3
                  else
                    isfluid(i,j,k)=0
                  endif
                else
                  isfluid(i,j,k)=0
                endif
              elseif(j==(ny+1))then !front side
                if(bc_type_front/=0)then
                  if(isfluid(i,j+ez(4),k)==1)then
                    isfluid(i,j,k)= bc_type_front+5
                    bcfluid(i,j,k)=4
                  else
                    isfluid(i,j,k)=0
                  endif
                else
                  isfluid(i,j,k)=0
                endif
              else
                isfluid(i,j,k)=0
              endif
            endif
          elseif(ibctype==5)then ! 5 T F T
            if(j==0 .or. j==(ny+1))then
              if(j==0)then !rear side
                if(bc_type_rear/=0)then
                  if(isfluid(i,j+ez(3),k)==1)then
                    isfluid(i,j,k)= bc_type_rear+5
                    bcfluid(i,j,k)=3
                  else
                    isfluid(i,j,k)=0
                  endif
                else
                  isfluid(i,j,k)=0
                endif
              elseif(j==(ny+1))then !front side
                if(bc_type_front/=0)then
                  if(isfluid(i,j+ez(4),k)==1)then
                    isfluid(i,j,k)= bc_type_front+5
                    bcfluid(i,j,k)=4
                  else
                    isfluid(i,j,k)=0
                  endif
                else
                  isfluid(i,j,k)=0
                endif
              else
                isfluid(i,j,k)=0
              endif
            endif
          elseif(ibctype==6)then ! 6 F T T
            if(i==0 .or. i==(nx+1))then
              if(i==(nx+1))then !east side
                if(bc_type_east/=0)then
                  if(isfluid(i+ex(2),j,k)==1)then
                    isfluid(i,j,k)=bc_type_east+5
                    bcfluid(i,j,k)=2
                  else
                    isfluid(i,j,k)=0
                  endif
                else
                  isfluid(i,j,k)=0
                endif
              elseif(i==0)then !west side
                if(bc_type_west/=0)then
                  if(isfluid(i+ex(1),j,k)==1)then
                    isfluid(i,j,k)=bc_type_west+5
                    bcfluid(i,j,k)=1
                  else
                    isfluid(i,j,k)=0
                  endif
                else
                  isfluid(i,j,k)=0
                endif
              else
                isfluid(i,j,k)=0
              endif
            endif
          endif 
        endif
      enddo
    enddo
  enddo
  
  
  !communicate isfluid over the processes applying the bc if necessary
  call comm_init_isfluid(isfluid)
  
  !communicate bcfluid over the processes applying the bc if necessary
  call comm_init_isfluid(bcfluid)
  
  !initialize the bc if necessary within the same process
  call initialiaze_manage_bc_isfluid_selfcomm
  
  !initialize the bc if necessary within the same process for a single halo
  call initialiaze_manage_bc_singlehalo_selfcomm
  
  !apply the bc if necessary within the same process
  call manage_bc_isfluid_selfcomm(isfluid)
  
  !apply the bc if necessary within the same process
  call manage_bc_isfluid_selfcomm(bcfluid)
  
  ltest=.false.
  do k=minz-1,maxz+1
    do j=miny-1,maxy+1
      do i=minx-1,maxx+1
        if(isfluid(i,j,k)>5 .and. bcfluid(i,j,k)==0)then
          !open boundary node without a direction is evil
          ltest=.true.
        endif
      enddo
    enddo
  enddo
  call or_world_larr(ltest,1)
  if(ltest(1))call error(16)
  
  if(lread_isfluid)then
    do l=1,links
      ishift=ex(l)
      jshift=ey(l)
      kshift=ez(l)
      do k=minz,maxz
        do j=miny,maxy
          do i=minx,maxx
            itemp=i+ishift
            jtemp=j+jshift
            ktemp=k+kshift
            if(itemp>=1 .and. itemp<=nx .and. jtemp>=1 .and. jtemp<=ny & 
             .and. ktemp>=1 .and. ktemp<=nz)then
              if(isfluid(i,j,k)==3 .and. isfluid(itemp,jtemp,ktemp)==1)then
                isfluid(i,j,k)=0 
              endif
            endif
          enddo
        enddo
      enddo
    enddo
    !communicate isfluid over the processes applying the bc if necessary
    call comm_init_isfluid(isfluid)
    !apply the bc if necessary within the same process
    call manage_bc_isfluid_selfcomm(isfluid)
  endif
  

#if 0
  open(unit=87,file='mio.xyz',status='replace')
  write(87,*)(nx+2)*(ny+2)*(nz+2)
  write(87,*)
  do k=minz-1,maxz+1
        do j=miny-1,maxy+1
          do i=minx-1,maxx+1
          write(87,'(i1,3f16.6)')isfluid(i,j,k),dble(i),dble(j),dble(k)
       enddo
       enddo
       enddo
  close(87)
#endif
  
  l=0
  do k=minz-1,maxz+1
    do j=miny-1,maxy+1
      do i=minx-1,maxx+1
        if(isfluid(i,j,k)==0)then
          l=l+1
        endif
      enddo
    enddo
  enddo
  nbounce0=l
  nbounce6dir(0)=l
  do idir=1,nbcdir
    do k=minz-1,maxz+1
      do j=miny-1,maxy+1
        do i=minx-1,maxx+1
          if(isfluid(i,j,k)==6 .and. bcfluid(i,j,k)==idir)then
            l=l+1
          endif
        enddo
      enddo
    enddo
    nbounce6dir(idir)=l
  enddo
  nbounce6=l
  
  nbounce7dir(0)=l
  do idir=1,nbcdir
    do k=minz-1,maxz+1
      do j=miny-1,maxy+1
        do i=minx-1,maxx+1
          if(isfluid(i,j,k)==7 .and. bcfluid(i,j,k)==idir)then
            l=l+1
          endif
        enddo
      enddo
    enddo
    nbounce7dir(idir)=l
  enddo
  nbounce7=l
  
  nbounce8dir(0)=l
  do idir=1,nbcdir
    do k=minz-1,maxz+1
      do j=miny-1,maxy+1
        do i=minx-1,maxx+1
          if(isfluid(i,j,k)==8 .and. bcfluid(i,j,k)==idir)then
            l=l+1
          endif
        enddo
      enddo
    enddo
    nbounce8dir(idir)=l
  enddo
  nbounce8=l
  
  nbounce9dir(0)=l
  do idir=1,nbcdir
    do k=minz-1,maxz+1
      do j=miny-1,maxy+1
        do i=minx-1,maxx+1
          if(isfluid(i,j,k)==9 .and. bcfluid(i,j,k)==idir)then
            l=l+1
          endif
        enddo
      enddo
    enddo
    nbounce9dir(idir)=l
  enddo
  nbounce9=l
  
  !allocate ibounce for inderect addressing 
  !NOTE that it is computationally less expensive than do a loop over all
  allocate(ibounce(3,nbounce9))
  
  l=0
  do k=minz-1,maxz+1
    do j=miny-1,maxy+1
      do i=minx-1,maxx+1
        if(isfluid(i,j,k)==0)then
          l=l+1
          ibounce(1,l)=i
          ibounce(2,l)=j
          ibounce(3,l)=k
        endif
      enddo
    enddo
  enddo
  
  !I apply an order to avoid an if in run time
  do idir=1,nbcdir
    do k=minz-1,maxz+1
      do j=miny-1,maxy+1
        do i=minx-1,maxx+1
          if(isfluid(i,j,k)==6 .and. bcfluid(i,j,k)==idir)then
            l=l+1
            ibounce(1,l)=i
            ibounce(2,l)=j
            ibounce(3,l)=k
          endif
        enddo
      enddo
    enddo
  enddo

  do idir=1,nbcdir
    do k=minz-1,maxz+1
      do j=miny-1,maxy+1
        do i=minx-1,maxx+1
          if(isfluid(i,j,k)==7 .and. bcfluid(i,j,k)==idir)then
            l=l+1
            ibounce(1,l)=i
            ibounce(2,l)=j
            ibounce(3,l)=k
          endif
        enddo
      enddo
    enddo
  enddo

  do idir=1,nbcdir
    do k=minz-1,maxz+1
      do j=miny-1,maxy+1
        do i=minx-1,maxx+1
          if(isfluid(i,j,k)==8 .and. bcfluid(i,j,k)==idir)then
            l=l+1
            ibounce(1,l)=i
            ibounce(2,l)=j
            ibounce(3,l)=k
          endif
        enddo
      enddo
    enddo
  enddo
  
  do idir=1,nbcdir
    do k=minz-1,maxz+1
      do j=miny-1,maxy+1
        do i=minx-1,maxx+1
          if(isfluid(i,j,k)==9 .and. bcfluid(i,j,k)==idir)then
            l=l+1
            ibounce(1,l)=i
            ibounce(2,l)=j
            ibounce(3,l)=k
          endif
        enddo
      enddo
    enddo
  enddo
  
  
  istat(1:nistatmax)=0
  allocate(bc_rhoR(nbounce0+1:nbounce9),stat=istat(1))
  allocate(bc_u(nbounce0+1:nbounce9),stat=istat(2))
  allocate(bc_v(nbounce0+1:nbounce9),stat=istat(3))
  allocate(bc_w(nbounce0+1:nbounce9),stat=istat(4))
  if(.not. lsingle_fluid)then
    allocate(bc_rhoB(nbounce0+1:nbounce9),stat=istat(5))
  endif
  allocate(bc_flow(nbounce0+1:nbounce9),stat=istat(6))
  allocate(bc_omega(nbounce0+1:nbounce9),stat=istat(7))

  !set the bc if necessary with their fixed values given in input
  call set_bc_hvar
  
  ltest=.false.
  if(any(istat.ne.0))then
    do i=1,nistatmax
      if(istat(i).ne.0)exit
    enddo
    call warning(2,dble(i))
     ltest=.true.
   endif
   call or_world_larr(ltest,1)
   if(ltest(1))call error(15)
   
   !certain bc need few hydrodynamic variables in the halo
   !check if it is the case and gives the order to communicate the hvar 
   !within the halo
   if(bc_type_east==1 .or. bc_type_west==1)then
     lexch_u=.true.
   endif
   if(bc_type_rear==1 .or. bc_type_front==1)then
     lexch_v=.true.
   endif
   if(bc_type_north==1 .or. bc_type_south==1)then
     lexch_w=.true.
   endif
   if(bc_type_east==2 .or. bc_type_west==2 .or. &
    bc_type_rear==2 .or. bc_type_front==2 .or. &
    bc_type_north==2 .or. bc_type_south==2)then
     lexch_dens=.true.
   endif
   
   if((lvtkfilesub .or. lparticles))then
     lexch_dens=.true.
     lexch_u=.true.
     lexch_v=.true.
     lexch_w=.true.
   endif

 end subroutine initialize_isfluid_bcfluid
 

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
  
  call manage_bc_hvar_selfcomm(rhoR_managebc)
  if(.not. lsingle_fluid)then
    call manage_bc_hvar_selfcomm(rhoB_managebc)
  endif
  
  call commwait_dens(rhoR,rhoB)
  
  
#ifdef DIAGNINIT
  call print_all_hvar(100,'inithvar_dopo',0,rhoR,u,v,w)
#endif
  
  call compute_densities_wall(0)
  
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
  call print_all_pops(100,'initpop',0,aoptpR)
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
  
    do k=1,nz
      do j=1,ny
        do i=1,nx
          rhoR(i,j,k)=meanR+stdevR*gauss()
        enddo
      enddo
    enddo
    
    if(lsingle_fluid)return
    
    do k=1,nz
      do j=1,ny
        do i=1,nx
          rhoB(i,j,k)=meanB+stdevB*gauss()
        enddo
      enddo
    enddo
  
  else
  
    forall (i=minx:maxx,j=miny:maxy,k=minz:maxz)
      rhoR(i,j,k)=meanR+stdevR*gauss_noseeded(i,j,k,1)
    end forall
    
    if(lsingle_fluid)return
    
    forall (i=minx:maxx,j=miny:maxy,k=minz:maxz)
      rhoB(i,j,k)=meanB+stdevB*gauss_noseeded(i,j,k,100)
    end forall
  
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
  
  if(lcolourG)then
  
    if(lunique_omega)then
      call set_initial_pop_fluids_CG_unique(alphaR_CG,rhoR,u,v,w,aoptpR)
      call set_initial_pop_fluids_CG_unique(alphaB_CG,rhoB,u,v,w,aoptpB)
    else
      call set_initial_pop_fluids_CG(alphaR_CG,rhoR,u,v,w,aoptpR)
      call set_initial_pop_fluids_CG(alphaB_CG,rhoB,u,v,w,aoptpB)
    endif
    
   else
  
    call set_initial_pop_fluids(rhoR,u,v,w,f00R,f01R,&
     f02R,f03R,f04R,f05R,f06R,f07R,f08R,f09R,f10R, &
     f11R,f12R,f13R,f14R,f15R,f16R,f17R,f18R,aoptpR)
  
    if(lsingle_fluid)return
  
    call set_initial_pop_fluids(rhoB,u,v,w,f00B,f01B,&
     f02B,f03B,f04B,f05B,f06B,f07B,f08B,f09B,f10B, &
     f11B,f12B,f13B,f14B,f15B,f16B,f17B,f18B,aoptpB)
  
  endif  
  
  return
  
 end subroutine driver_set_initial_pop_fluids
 
 subroutine set_initial_pop_fluids(rhosub,usub,vsub,wsub,f00sub,f01sub,&
  f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub,f09sub,f10sub, &
  f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub,aoptp)
 
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
   wsub,f00sub,f01sub,f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub,&
   f09sub,f10sub,f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub
  type(REALPTR), dimension(0:links):: aoptp
  
  integer :: i,j,k,l
  
  
  
  
  forall(i=minx-1:maxx+1,j=miny-1:maxy+1,k=minz-1:maxz+1)
    f00sub(i,j,k)=equil_pop00(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k), &
     wsub(i,j,k))
  end forall
  
  forall(i=minx-1:maxx+1,j=miny-1:maxy+1,k=minz-1:maxz+1)
    f01sub(i,j,k)=equil_pop01(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k), &
     wsub(i,j,k))
  end forall
  
  forall(i=minx-1:maxx+1,j=miny-1:maxy+1,k=minz-1:maxz+1)
    f02sub(i,j,k)=equil_pop02(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k), &
     wsub(i,j,k))
  end forall
  
  forall(i=minx-1:maxx+1,j=miny-1:maxy+1,k=minz-1:maxz+1)
    f03sub(i,j,k)=equil_pop03(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k), &
     wsub(i,j,k))
  end forall
  
  forall(i=minx-1:maxx+1,j=miny-1:maxy+1,k=minz-1:maxz+1)
    f04sub(i,j,k)=equil_pop04(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k), &
     wsub(i,j,k))
  end forall
  
  forall(i=minx-1:maxx+1,j=miny-1:maxy+1,k=minz-1:maxz+1)
    f05sub(i,j,k)=equil_pop05(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k), &
     wsub(i,j,k))
  end forall
  
  forall(i=minx-1:maxx+1,j=miny-1:maxy+1,k=minz-1:maxz+1)
    f06sub(i,j,k)=equil_pop06(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k), &
     wsub(i,j,k))
  end forall
  
  forall(i=minx-1:maxx+1,j=miny-1:maxy+1,k=minz-1:maxz+1)
    f07sub(i,j,k)=equil_pop07(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k), &
     wsub(i,j,k))
  end forall
  
  forall(i=minx-1:maxx+1,j=miny-1:maxy+1,k=minz-1:maxz+1)
    f08sub(i,j,k)=equil_pop08(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k), &
     wsub(i,j,k))
  end forall
  
  forall(i=minx-1:maxx+1,j=miny-1:maxy+1,k=minz-1:maxz+1)
    f09sub(i,j,k)=equil_pop09(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k), &
     wsub(i,j,k))
  end forall
  
  forall(i=minx-1:maxx+1,j=miny-1:maxy+1,k=minz-1:maxz+1)
    f10sub(i,j,k)=equil_pop10(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k), &
     wsub(i,j,k))
  end forall
  
  forall(i=minx-1:maxx+1,j=miny-1:maxy+1,k=minz-1:maxz+1)
    f11sub(i,j,k)=equil_pop11(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k), &
     wsub(i,j,k))
  end forall
  
  forall(i=minx-1:maxx+1,j=miny-1:maxy+1,k=minz-1:maxz+1)
    f12sub(i,j,k)=equil_pop12(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k), &
     wsub(i,j,k))
  end forall
  
  forall(i=minx-1:maxx+1,j=miny-1:maxy+1,k=minz-1:maxz+1)
    f13sub(i,j,k)=equil_pop13(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k), &
     wsub(i,j,k))
  end forall
  
  forall(i=minx-1:maxx+1,j=miny-1:maxy+1,k=minz-1:maxz+1)
    f14sub(i,j,k)=equil_pop14(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k), &
     wsub(i,j,k))
  end forall
  
  forall(i=minx-1:maxx+1,j=miny-1:maxy+1,k=minz-1:maxz+1)
    f15sub(i,j,k)=equil_pop15(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k), &
     wsub(i,j,k))
  end forall
  
  forall(i=minx-1:maxx+1,j=miny-1:maxy+1,k=minz-1:maxz+1)
    f16sub(i,j,k)=equil_pop16(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k), &
     wsub(i,j,k))
  end forall
  
  forall(i=minx-1:maxx+1,j=miny-1:maxy+1,k=minz-1:maxz+1)
    f17sub(i,j,k)=equil_pop17(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k), &
     wsub(i,j,k))
  end forall
  
  forall(i=minx-1:maxx+1,j=miny-1:maxy+1,k=minz-1:maxz+1)
    f18sub(i,j,k)=equil_pop18(rhosub(i,j,k),usub(i,j,k),vsub(i,j,k), &
     wsub(i,j,k))
  end forall
  
  
  return
  
 end subroutine set_initial_pop_fluids
 
 subroutine set_initial_pop_fluids_CG_unique(alpha_CG,rhosub,usub,vsub,wsub,aoptp)
 
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
  
  real(kind=PRC), intent(in) :: alpha_CG
  real(kind=PRC), allocatable, dimension(:,:,:)  :: rhosub,usub,vsub, &
   wsub
  type(REALPTR), dimension(0:links):: aoptp
  
  real(kind=PRC) :: grad_rhox,grad_rhoy,grad_rhoz,temp_omega
  integer :: i,j,k,l
  
  
  temp_omega=unique_omega
  
   do k=minz,maxz
     do j=miny,maxy
       do i=minx,maxx
  
         grad_rhox = ZERO
            grad_rhoy = ZERO
            grad_rhoz = ZERO

#ifdef GRADIENTD3Q27
         do l=1,linksd3q27
           grad_rhox=grad_rhox + ad3q27(l)*dexd3q27(l)* &
            (rhoR(i+exd3q27(l),j+eyd3q27(l),k+ezd3q27(l)))
           grad_rhoy=grad_rhoy + ad3q27(l)*deyd3q27(l)* &
            (rhoR(i+exd3q27(l),j+eyd3q27(l),k+ezd3q27(l)))
           grad_rhoz=grad_rhoz + ad3q27(l)*dezd3q27(l)* &
            (rhoR(i+exd3q27(l),j+eyd3q27(l),k+ezd3q27(l)))
         enddo
#else             
         do l=1,links
           grad_rhox=grad_rhox + a(l)*dex(l)* &
            (rhoR(i+ex(l),j+ey(l),k+ez(l)))
           grad_rhoy=grad_rhoy + a(l)*dey(l)* &
            (rhoR(i+ex(l),j+ey(l),k+ez(l)))
           grad_rhoz=grad_rhoz + a(l)*dez(l)* &
            (rhoR(i+ex(l),j+ey(l),k+ez(l)))
         enddo
#endif    

          do l=0,links
            aoptp(l)%p(i,j,k)= &
             equil_popCG(l,temp_omega,alpha_CG,rhosub(i,j,k), &
             usub(i,j,k),vsub(i,j,k),wsub(i,j,k), &
             grad_rhox,grad_rhoy,grad_rhoz)
          enddo
        enddo
      enddo
    enddo
  
  return
  
 end subroutine set_initial_pop_fluids_CG_unique
 
 subroutine set_initial_pop_fluids_CG(alpha_CG,rhosub,usub,vsub,wsub,aoptp)
 
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
  
  real(kind=PRC), intent(in) :: alpha_CG
  real(kind=PRC) :: grad_rhox,grad_rhoy,grad_rhoz,temp_omega
  integer :: i,j,k,l
  real(kind=PRC), allocatable, dimension(:,:,:)  :: rhosub,usub,vsub, &
   wsub
  type(REALPTR), dimension(0:links):: aoptp
  
  
   do k=minz,maxz
     do j=miny,maxy
       do i=minx,maxx
         
         temp_omega=omega(i,j,k)
         
         grad_rhox = ZERO
            grad_rhoy = ZERO
            grad_rhoz = ZERO

#ifdef GRADIENTD3Q27
         do l=1,linksd3q27
           grad_rhox=grad_rhox + ad3q27(l)*dexd3q27(l)* &
            (rhoR(i+exd3q27(l),j+eyd3q27(l),k+ezd3q27(l)))
           grad_rhoy=grad_rhoy + ad3q27(l)*deyd3q27(l)* &
            (rhoR(i+exd3q27(l),j+eyd3q27(l),k+ezd3q27(l)))
           grad_rhoz=grad_rhoz + ad3q27(l)*dezd3q27(l)* &
            (rhoR(i+exd3q27(l),j+eyd3q27(l),k+ezd3q27(l)))
         enddo
#else             
         do l=1,links
           grad_rhox=grad_rhox + a(l)*dex(l)* &
            (rhoR(i+ex(l),j+ey(l),k+ez(l)))
           grad_rhoy=grad_rhoy + a(l)*dey(l)* &
            (rhoR(i+ex(l),j+ey(l),k+ez(l)))
           grad_rhoz=grad_rhoz + a(l)*dez(l)* &
            (rhoR(i+ex(l),j+ey(l),k+ez(l)))
         enddo
#endif    

          do l=0,links
            aoptp(l)%p(i,j,k)= &
             equil_popCG(l,temp_omega,alpha_CG,rhosub(i,j,k), &
             usub(i,j,k),vsub(i,j,k),wsub(i,j,k), &
             grad_rhox,grad_rhoy,grad_rhoz)
          enddo
        enddo
      enddo
    enddo
  
  return
  
 end subroutine set_initial_pop_fluids_CG
 
 subroutine set_value_bc_east(itemp1,dtemp1,dtemp2,dtemp3,dtemp4,dtemp5, &
  dtemp6)
 
!***********************************************************************
!     
!     LBsoft subroutine for set the conditions of fluids 
!     at the east open boundary
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: itemp1
  real(kind=PRC), intent(in) :: dtemp1,dtemp2,dtemp3,dtemp4,dtemp5,dtemp6
  
  bc_type_east = itemp1
  bc_rhoR_east = dtemp1
  bc_rhoB_east = dtemp2
  bc_u_east = dtemp3
  bc_v_east = dtemp4
  bc_w_east = dtemp5
  bc_flow_east = dtemp6
  
  return
  
 end subroutine set_value_bc_east
 
 subroutine set_value_bc_west(itemp1,dtemp1,dtemp2,dtemp3,dtemp4,dtemp5, &
  dtemp6)
 
!***********************************************************************
!     
!     LBsoft subroutine for set the conditions of fluids 
!     at the west open boundary
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: itemp1
  real(kind=PRC), intent(in) :: dtemp1,dtemp2,dtemp3,dtemp4,dtemp5,dtemp6
  
  bc_type_west = itemp1
  bc_rhoR_west = dtemp1
  bc_rhoB_west = dtemp2
  bc_u_west = dtemp3
  bc_v_west = dtemp4
  bc_w_west = dtemp5
  bc_flow_west = dtemp6
  
  return
  
 end subroutine set_value_bc_west
 
 subroutine set_value_bc_rear(itemp1,dtemp1,dtemp2,dtemp3,dtemp4,dtemp5, &
  dtemp6)
 
!***********************************************************************
!     
!     LBsoft subroutine for set the conditions of fluids 
!     at the rear open boundary
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: itemp1
  real(kind=PRC), intent(in) :: dtemp1,dtemp2,dtemp3,dtemp4,dtemp5,dtemp6
  
  bc_type_rear = itemp1
  bc_rhoR_rear = dtemp1
  bc_rhoB_rear = dtemp2
  bc_u_rear = dtemp3
  bc_v_rear = dtemp4
  bc_w_rear = dtemp5
  bc_flow_rear = dtemp6
  
  return
  
 end subroutine set_value_bc_rear
 
 subroutine set_value_bc_front(itemp1,dtemp1,dtemp2,dtemp3,dtemp4,dtemp5, &
  dtemp6)
 
!***********************************************************************
!     
!     LBsoft subroutine for set the conditions of fluids 
!     at the front open boundary
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: itemp1
  real(kind=PRC), intent(in) :: dtemp1,dtemp2,dtemp3,dtemp4,dtemp5,dtemp6
  
  bc_type_front = itemp1
  bc_rhoR_front = dtemp1
  bc_rhoB_front = dtemp2
  bc_u_front = dtemp3
  bc_v_front = dtemp4
  bc_w_front = dtemp5
  bc_flow_front = dtemp6
  
  return
  
 end subroutine set_value_bc_front
 
 subroutine set_value_bc_north(itemp1,dtemp1,dtemp2,dtemp3,dtemp4,dtemp5, &
  dtemp6)
 
!***********************************************************************
!     
!     LBsoft subroutine for set the conditions of fluids 
!     at the north open boundary
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: itemp1
  real(kind=PRC), intent(in) :: dtemp1,dtemp2,dtemp3,dtemp4,dtemp5,dtemp6
  
  bc_type_north = itemp1
  bc_rhoR_north = dtemp1
  bc_rhoB_north = dtemp2
  bc_u_north = dtemp3
  bc_v_north = dtemp4
  if(bc_type_north==2)then
    bc_w_north = -dtemp5
  else
    bc_w_north = dtemp5
  endif
  bc_flow_north = dtemp6
  
  return
  
 end subroutine set_value_bc_north
 
 subroutine set_value_bc_south(itemp1,dtemp1,dtemp2,dtemp3,dtemp4,dtemp5, &
  dtemp6)
 
!***********************************************************************
!     
!     LBsoft subroutine for set the conditions of fluids 
!     at the south open boundary
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: itemp1
  real(kind=PRC), intent(in) :: dtemp1,dtemp2,dtemp3,dtemp4,dtemp5,dtemp6
  
  bc_type_south = itemp1
  bc_rhoR_south = dtemp1
  bc_rhoB_south = dtemp2
  bc_u_south = dtemp3
  bc_v_south = dtemp4
  if(bc_type_south==2)then
    bc_w_south = -dtemp5
  else
    bc_w_south = dtemp5
  endif
  bc_flow_south = dtemp6
  
  return
  
 end subroutine set_value_bc_south
 
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
 
 subroutine set_lread_isfluid(ltemp1)
 
!***********************************************************************
!     
!     LBsoft subroutine for set the value of lread_isfluid for reading 
!     the isfluid.dat input file
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
  
  implicit none
  
  logical, intent(in) :: ltemp1
  
  lread_isfluid=ltemp1
  
  return
  
 end subroutine set_lread_isfluid
 
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
  
#if LATTICE==319

  phiR_CG(0)=dtemp1
  phiR_CG(1:6)=(ONE-dtemp1) / TWELVE
  phiR_CG(7:18)=(ONE-dtemp1) / TWENTYFOUR
  
  phiB_CG(0)=dtemp2 
  phiB_CG(1:6)=(ONE-dtemp2) / TWELVE
  phiB_CG(7:18)=(ONE-dtemp2) / TWENTYFOUR
  
#else
  
  return
  
#endif
  
  return
  
 end subroutine initialize_phi_coeff_cg
 
 subroutine set_fluid_wall_mode(itemp1,dtemp1,dtemp2)
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the wetting of the fluid density 
!     at the wall by the Benzi's approach
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2019
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: itemp1
  real(kind=PRC), intent(in) :: dtemp1,dtemp2
  
  iselwetting = itemp1
  densR_wetting = dtemp1
  densB_wetting = dtemp2
  
  return
  
 end subroutine set_fluid_wall_mode
 
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
  
  real(kind=PRC), parameter :: myfactor = TWO / cssq
  
  real(kind=PRC) :: viscosity_to_omega
  
  viscosity_to_omega = TWO / ( myfactor * dtemp1 + ONE )
  
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
      call rescale_pops(rescalef,f00R,f01R,&
       f02R,f03R,f04R,f05R,f06R,f07R,f08R,f09R,f10R, &
       f11R,f12R,f13R,f14R,f15R,f16R,f17R,f18R)
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
      call rescale_pops(rescalef,f00R,f01R,&
       f02R,f03R,f04R,f05R,f06R,f07R,f08R,f09R,f10R, &
       f11R,f12R,f13R,f14R,f15R,f16R,f17R,f18R)
    endif
    !check the deviation
    if((abs(totnfluidmassB-totnfluidmassB_init)/totnfluidmassB_init)> &
     dmass_rescale)then
      rescalef=(totnfluidmassB_init*totnfluidnodes)/ &
       (totnfluidmassB*real(totnfluidnodes_init,kind=PRC))
      call rescale_pops(rescalef,f00B,f01B,&
       f02B,f03B,f04B,f05B,f06B,f07B,f08B,f09B,f10B, &
       f11B,f12B,f13B,f14B,f15B,f16B,f17B,f18B)
    endif
    
  endif
  
  return
  
 end subroutine rescale_fluid_mass
 
 subroutine rescale_pops(rescalesub,f00sub,f01sub,&
  f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub,f09sub,f10sub, &
  f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub)
 
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
  real(kind=PRC), allocatable, dimension(:,:,:)  :: &
   f00sub,f01sub,f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub,&
   f09sub,f10sub,f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub
  
  integer :: i,j,k
  
  do k=minz,maxz
    do j=miny,maxy
      do i=minx,maxx
        if (isfluid(i,j,k)/=1) cycle
        f00sub(i,j,k)=f00sub(i,j,k)*rescalesub
        f01sub(i,j,k)=f01sub(i,j,k)*rescalesub
        f02sub(i,j,k)=f02sub(i,j,k)*rescalesub
        f03sub(i,j,k)=f03sub(i,j,k)*rescalesub
        f04sub(i,j,k)=f04sub(i,j,k)*rescalesub
        f05sub(i,j,k)=f05sub(i,j,k)*rescalesub
        f06sub(i,j,k)=f06sub(i,j,k)*rescalesub
        f07sub(i,j,k)=f07sub(i,j,k)*rescalesub
        f08sub(i,j,k)=f08sub(i,j,k)*rescalesub
        f09sub(i,j,k)=f09sub(i,j,k)*rescalesub
        f10sub(i,j,k)=f10sub(i,j,k)*rescalesub
        f11sub(i,j,k)=f11sub(i,j,k)*rescalesub
        f12sub(i,j,k)=f12sub(i,j,k)*rescalesub
        f13sub(i,j,k)=f13sub(i,j,k)*rescalesub
        f14sub(i,j,k)=f14sub(i,j,k)*rescalesub
        f15sub(i,j,k)=f15sub(i,j,k)*rescalesub
        f16sub(i,j,k)=f16sub(i,j,k)*rescalesub
        f17sub(i,j,k)=f17sub(i,j,k)*rescalesub
        f18sub(i,j,k)=f18sub(i,j,k)*rescalesub
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
      call collision_fluids_BGK(rhoR,u,v,w,omega,f00R,f01R,&
       f02R,f03R,f04R,f05R,f06R,f07R,f08R,f09R,f10R, &
       f11R,f12R,f13R,f14R,f15R,f16R,f17R,f18R,aoptpR,nstep)
      if(lsingle_fluid)return
      call collision_fluids_BGK(rhoB,u,v,w,omega,f00B,f01B,&
       f02B,f03B,f04B,f05B,f06B,f07B,f08B,f09B,f10B, &
       f11B,f12B,f13B,f14B,f15B,f16B,f17B,f18B,aoptpB,nstep)

    case (1)
      ! done in collision_EDM) call convert_fluid_force_to_velshifted
      call collision_fluids_EDM(rhoR,u,v,w,fuR,fvR,fwR,omega,f00R,f01R,&
       f02R,f03R,f04R,f05R,f06R,f07R,f08R,f09R,f10R, &
       f11R,f12R,f13R,f14R,f15R,f16R,f17R,f18R,aoptpR,nstep)
      if(lsingle_fluid)return
      call collision_fluids_EDM(rhoB,u,v,w,fuB,fvB,fwB,omega,f00B,f01B,&
       f02B,f03B,f04B,f05B,f06B,f07B,f08B,f09B,f10B, &
       f11B,f12B,f13B,f14B,f15B,f16B,f17B,f18B,aoptpB,nstep)

    case default
      call error(14)
    end select
  
  else
    
    call driver_collision_fluids_CG(nstep)
    
  endif
  
  
 end subroutine driver_collision_fluids

 
 subroutine collision_fluids_BGK(rhosub,usub,vsub,wsub,omegas,f00sub,f01sub,&
  f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub,f09sub,f10sub, &
  f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub,aoptp,nstep)
 
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
   wsub,omegas,f00sub,f01sub,f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub,&
   f09sub,f10sub,f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub
  type(REALPTR), dimension(0:links):: aoptp
  integer :: l
  real(kind=PRC) :: locrho,locu,locv,locw, oneminusomega, temp_omega
#ifdef REGULARIZED
  real(kind=PRC) :: pxx,pyy,pzz,pxy,pxz,pyz
  real(kind=PRC), dimension(0:links) :: feq
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
          fneq = aoptp(l)%p(i,j,k) - feq(l)
			
	      !non equilibrium part of the momentum flux tensor
          pxx=pxx + (dex(l)*dex(l) - cssq)*fneq
          pyy=pyy + (dey(l)*dey(l) - cssq)*fneq
          pzz=pzz + (dez(l)*dez(l) - cssq)*fneq
          pxy=pxy +  dex(l)*dey(l)*fneq
          pxz=pxz +  dex(l)*dez(l)*fneq
          pyz=pyz +  dey(l)*dez(l)*fneq
        enddo
        
        do l=0,links
          aoptp(l)%p(i,j,k)= feq(l) + &
           ((HALF*p(l))/(cssq**TWO))*((dex(l)*dex(l) - cssq)*pxx + &
           (dey(l)*dey(l) - cssq)*pyy + &
           (dez(l)*dez(l) - cssq)*pzz + &
           TWO*dex(l)*dey(l)*pxy + &
           TWO*dex(l)*dez(l)*pxz + &
           TWO*dey(l)*dez(l)*pyz)           
        enddo
        f00sub(i,j,k)=f00sub(i,j,k)*oneminusomega + feq(0) * temp_omega
        f01sub(i,j,k)=f01sub(i,j,k)*oneminusomega + feq(1) * temp_omega
        f02sub(i,j,k)=f02sub(i,j,k)*oneminusomega + feq(2) * temp_omega
        f03sub(i,j,k)=f03sub(i,j,k)*oneminusomega + feq(3) * temp_omega
        f04sub(i,j,k)=f04sub(i,j,k)*oneminusomega + feq(4) * temp_omega
        f05sub(i,j,k)=f05sub(i,j,k)*oneminusomega + feq(5) * temp_omega
        f06sub(i,j,k)=f06sub(i,j,k)*oneminusomega + feq(6) * temp_omega
        f07sub(i,j,k)=f07sub(i,j,k)*oneminusomega + feq(7) * temp_omega
        f08sub(i,j,k)=f08sub(i,j,k)*oneminusomega + feq(8) * temp_omega
        f09sub(i,j,k)=f09sub(i,j,k)*oneminusomega + feq(9) * temp_omega
        f10sub(i,j,k)=f10sub(i,j,k)*oneminusomega + feq(10) * temp_omega
        f11sub(i,j,k)=f11sub(i,j,k)*oneminusomega + feq(11) * temp_omega
        f12sub(i,j,k)=f12sub(i,j,k)*oneminusomega + feq(12) * temp_omega
        f13sub(i,j,k)=f13sub(i,j,k)*oneminusomega + feq(13) * temp_omega
        f14sub(i,j,k)=f14sub(i,j,k)*oneminusomega + feq(14) * temp_omega
        f15sub(i,j,k)=f15sub(i,j,k)*oneminusomega + feq(15) * temp_omega
        f16sub(i,j,k)=f16sub(i,j,k)*oneminusomega + feq(16) * temp_omega
        f17sub(i,j,k)=f17sub(i,j,k)*oneminusomega + feq(17) * temp_omega
        f18sub(i,j,k)=f18sub(i,j,k)*oneminusomega + feq(18) * temp_omega
#else
        f00sub(i,j,k)=f00sub(i,j,k)*oneminusomega + &
         equil_pop00(locrho,locu,locv,locw) * temp_omega
        f01sub(i,j,k)=f01sub(i,j,k)*oneminusomega + &
         equil_pop01(locrho,locu,locv,locw) * temp_omega
        f02sub(i,j,k)=f02sub(i,j,k)*oneminusomega + &
         equil_pop02(locrho,locu,locv,locw) * temp_omega
        f03sub(i,j,k)=f03sub(i,j,k)*oneminusomega + &
         equil_pop03(locrho,locu,locv,locw) * temp_omega
        f04sub(i,j,k)=f04sub(i,j,k)*oneminusomega + &
         equil_pop04(locrho,locu,locv,locw) * temp_omega
        f05sub(i,j,k)=f05sub(i,j,k)*oneminusomega + &
         equil_pop05(locrho,locu,locv,locw) * temp_omega
        f06sub(i,j,k)=f06sub(i,j,k)*oneminusomega + &
         equil_pop06(locrho,locu,locv,locw) * temp_omega
        f07sub(i,j,k)=f07sub(i,j,k)*oneminusomega + &
         equil_pop07(locrho,locu,locv,locw) * temp_omega
        f08sub(i,j,k)=f08sub(i,j,k)*oneminusomega + &
         equil_pop08(locrho,locu,locv,locw) * temp_omega
        f09sub(i,j,k)=f09sub(i,j,k)*oneminusomega + &
         equil_pop09(locrho,locu,locv,locw) * temp_omega
        f10sub(i,j,k)=f10sub(i,j,k)*oneminusomega + &
         equil_pop10(locrho,locu,locv,locw) * temp_omega
        f11sub(i,j,k)=f11sub(i,j,k)*oneminusomega + &
         equil_pop11(locrho,locu,locv,locw) * temp_omega
        f12sub(i,j,k)=f12sub(i,j,k)*oneminusomega + &
         equil_pop12(locrho,locu,locv,locw) * temp_omega
        f13sub(i,j,k)=f13sub(i,j,k)*oneminusomega + &
         equil_pop13(locrho,locu,locv,locw) * temp_omega
        f14sub(i,j,k)=f14sub(i,j,k)*oneminusomega + &
         equil_pop14(locrho,locu,locv,locw) * temp_omega
        f15sub(i,j,k)=f15sub(i,j,k)*oneminusomega + &
         equil_pop15(locrho,locu,locv,locw) * temp_omega
        f16sub(i,j,k)=f16sub(i,j,k)*oneminusomega + &
         equil_pop16(locrho,locu,locv,locw) * temp_omega
        f17sub(i,j,k)=f17sub(i,j,k)*oneminusomega + &
         equil_pop17(locrho,locu,locv,locw) * temp_omega
        f18sub(i,j,k)=f18sub(i,j,k)*oneminusomega + &
         equil_pop18(locrho,locu,locv,locw) * temp_omega
#endif
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
          fneq = aoptp(l)%p(i,j,k) - feq(l)
			
	      !non equilibrium part of the momentum flux tensor
          pxx=pxx + (dex(l)*dex(l) - cssq)*fneq
          pyy=pyy + (dey(l)*dey(l) - cssq)*fneq
          pzz=pzz + (dez(l)*dez(l) - cssq)*fneq
          pxy=pxy +  dex(l)*dey(l)*fneq
          pxz=pxz +  dex(l)*dez(l)*fneq
          pyz=pyz +  dey(l)*dez(l)*fneq
        enddo
        
        do l=0,links
          aoptp(l)%p(i,j,k)= feq(l) + &
           ((HALF*p(l))/(cssq**TWO))*((dex(l)*dex(l) - cssq)*pxx + &
           (dey(l)*dey(l) - cssq)*pyy + &
           (dez(l)*dez(l) - cssq)*pzz + &
           TWO*dex(l)*dey(l)*pxy + &
           TWO*dex(l)*dez(l)*pxz + &
           TWO*dey(l)*dez(l)*pyz)           
        enddo
        f00sub(i,j,k)=f00sub(i,j,k)*oneminusomega + feq(0) * temp_omega
        f01sub(i,j,k)=f01sub(i,j,k)*oneminusomega + feq(1) * temp_omega
        f02sub(i,j,k)=f02sub(i,j,k)*oneminusomega + feq(2) * temp_omega
        f03sub(i,j,k)=f03sub(i,j,k)*oneminusomega + feq(3) * temp_omega
        f04sub(i,j,k)=f04sub(i,j,k)*oneminusomega + feq(4) * temp_omega
        f05sub(i,j,k)=f05sub(i,j,k)*oneminusomega + feq(5) * temp_omega
        f06sub(i,j,k)=f06sub(i,j,k)*oneminusomega + feq(6) * temp_omega
        f07sub(i,j,k)=f07sub(i,j,k)*oneminusomega + feq(7) * temp_omega
        f08sub(i,j,k)=f08sub(i,j,k)*oneminusomega + feq(8) * temp_omega
        f09sub(i,j,k)=f09sub(i,j,k)*oneminusomega + feq(9) * temp_omega
        f10sub(i,j,k)=f10sub(i,j,k)*oneminusomega + feq(10) * temp_omega
        f11sub(i,j,k)=f11sub(i,j,k)*oneminusomega + feq(11) * temp_omega
        f12sub(i,j,k)=f12sub(i,j,k)*oneminusomega + feq(12) * temp_omega
        f13sub(i,j,k)=f13sub(i,j,k)*oneminusomega + feq(13) * temp_omega
        f14sub(i,j,k)=f14sub(i,j,k)*oneminusomega + feq(14) * temp_omega
        f15sub(i,j,k)=f15sub(i,j,k)*oneminusomega + feq(15) * temp_omega
        f16sub(i,j,k)=f16sub(i,j,k)*oneminusomega + feq(16) * temp_omega
        f17sub(i,j,k)=f17sub(i,j,k)*oneminusomega + feq(17) * temp_omega
        f18sub(i,j,k)=f18sub(i,j,k)*oneminusomega + feq(18) * temp_omega
#else
        f00sub(i,j,k)=f00sub(i,j,k)*oneminusomega + &
         equil_pop00(locrho,locu,locv,locw) * temp_omega
        f01sub(i,j,k)=f01sub(i,j,k)*oneminusomega + &
         equil_pop01(locrho,locu,locv,locw) * temp_omega
        f02sub(i,j,k)=f02sub(i,j,k)*oneminusomega + &
         equil_pop02(locrho,locu,locv,locw) * temp_omega
        f03sub(i,j,k)=f03sub(i,j,k)*oneminusomega + &
         equil_pop03(locrho,locu,locv,locw) * temp_omega
        f04sub(i,j,k)=f04sub(i,j,k)*oneminusomega + &
         equil_pop04(locrho,locu,locv,locw) * temp_omega
        f05sub(i,j,k)=f05sub(i,j,k)*oneminusomega + &
         equil_pop05(locrho,locu,locv,locw) * temp_omega
        f06sub(i,j,k)=f06sub(i,j,k)*oneminusomega + &
         equil_pop06(locrho,locu,locv,locw) * temp_omega
        f07sub(i,j,k)=f07sub(i,j,k)*oneminusomega + &
         equil_pop07(locrho,locu,locv,locw) * temp_omega
        f08sub(i,j,k)=f08sub(i,j,k)*oneminusomega + &
         equil_pop08(locrho,locu,locv,locw) * temp_omega
        f09sub(i,j,k)=f09sub(i,j,k)*oneminusomega + &
         equil_pop09(locrho,locu,locv,locw) * temp_omega
        f10sub(i,j,k)=f10sub(i,j,k)*oneminusomega + &
         equil_pop10(locrho,locu,locv,locw) * temp_omega
        f11sub(i,j,k)=f11sub(i,j,k)*oneminusomega + &
         equil_pop11(locrho,locu,locv,locw) * temp_omega
        f12sub(i,j,k)=f12sub(i,j,k)*oneminusomega + &
         equil_pop12(locrho,locu,locv,locw) * temp_omega
        f13sub(i,j,k)=f13sub(i,j,k)*oneminusomega + &
         equil_pop13(locrho,locu,locv,locw) * temp_omega
        f14sub(i,j,k)=f14sub(i,j,k)*oneminusomega + &
         equil_pop14(locrho,locu,locv,locw) * temp_omega
        f15sub(i,j,k)=f15sub(i,j,k)*oneminusomega + &
         equil_pop15(locrho,locu,locv,locw) * temp_omega
        f16sub(i,j,k)=f16sub(i,j,k)*oneminusomega + &
         equil_pop16(locrho,locu,locv,locw) * temp_omega
        f17sub(i,j,k)=f17sub(i,j,k)*oneminusomega + &
         equil_pop17(locrho,locu,locv,locw) * temp_omega
        f18sub(i,j,k)=f18sub(i,j,k)*oneminusomega + &
         equil_pop18(locrho,locu,locv,locw) * temp_omega  
#endif

       enddo
   enddo
  enddo

  ! end forall
 end subroutine collision_fluids_BGK
 

 subroutine collision_fluids_EDM(rhosub,usub,vsub,wsub,fusub,fvsub, &
  fwsub,omegas,f00sub,f01sub,f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,&
  f08sub,f09sub,f10sub,f11sub,f12sub,f13sub,f14sub,f15sub,f16sub, &
  f17sub,f18sub,aoptp,nstep)
 
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
   wsub,fusub,fvsub,fwsub,omegas,f00sub,f01sub,f02sub,f03sub,f04sub, &
   f05sub,f06sub,f07sub,f08sub,f09sub,f10sub,f11sub,f12sub,f13sub, &
   f14sub,f15sub,f16sub,f17sub,f18sub
  type(REALPTR), dimension(0:links):: aoptp
  real(kind=PRC) :: locrho,locu,locv,locw, oneminusomega, &
        locfu,locfv,locfw,temp_omega,omegaminusone
  
  integer :: i,j,k,l
  
#ifdef REGULARIZED
  real(kind=PRC) :: pxx,pyy,pzz,pxy,pxz,pyz
  real(kind=PRC), dimension(0:links) :: feq
  real(kind=PRC) :: fneq
#endif
    
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
        
#ifdef SCANDREA
          locfu = fusub(i,j,k)*t_LB + locu
          locfv = fvsub(i,j,k)*t_LB + locv
          locfw = fwsub(i,j,k)*t_LB + locw
#else
          locfu = fusub(i,j,k)*t_LB / locrho + locu
          locfv = fvsub(i,j,k)*t_LB / locrho + locv
          locfw = fwsub(i,j,k)*t_LB / locrho + locw
#endif
          fusub(i,j,k) = locfu
          fvsub(i,j,k) = locfv
          fwsub(i,j,k) = locfw
          
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
            fneq = aoptp(l)%p(i,j,k) - feq(l)
		    
	        !non equilibrium part of the momentum flux tensor
            pxx=pxx + (dex(l)*dex(l) - cssq)*fneq
            pyy=pyy + (dey(l)*dey(l) - cssq)*fneq
            pzz=pzz + (dez(l)*dez(l) - cssq)*fneq
            pxy=pxy +  dex(l)*dey(l)*fneq
            pxz=pxz +  dex(l)*dez(l)*fneq
            pyz=pyz +  dey(l)*dez(l)*fneq
          enddo
        
          do l=0,links
            aoptp(l)%p(i,j,k)= feq(l) + &
             ((HALF*p(l))/(cssq**TWO))*((dex(l)*dex(l) - cssq)*pxx + &
             (dey(l)*dey(l) - cssq)*pyy + &
             (dez(l)*dez(l) - cssq)*pzz + &
             TWO*dex(l)*dey(l)*pxy + &
             TWO*dex(l)*dez(l)*pxz + &
             TWO*dey(l)*dez(l)*pyz)           
          enddo
           
          f00sub(i,j,k)=oneminusomega*f00sub(i,j,k)+omegaminusone* &
           feq(0)+ equil_pop00(locrho,locfu,locfv,locfw)
       
          f01sub(i,j,k)=oneminusomega*f01sub(i,j,k)+omegaminusone* &
           feq(1)+ equil_pop01(locrho,locfu,locfv,locfw)
    
          f02sub(i,j,k)=oneminusomega*f02sub(i,j,k)+omegaminusone* &
           feq(2)+ equil_pop02(locrho,locfu,locfv,locfw)
    
          f03sub(i,j,k)=oneminusomega*f03sub(i,j,k)+omegaminusone* &
           feq(3)+ equil_pop03(locrho,locfu,locfv,locfw)
       
          f04sub(i,j,k)=oneminusomega*f04sub(i,j,k)+omegaminusone* &
           feq(4)+ equil_pop04(locrho,locfu,locfv,locfw)
        
          f05sub(i,j,k)=oneminusomega*f05sub(i,j,k)+omegaminusone* &
           feq(5)+ equil_pop05(locrho,locfu,locfv,locfw)
        
          f06sub(i,j,k)=oneminusomega*f06sub(i,j,k)+omegaminusone* &
           feq(6)+ equil_pop06(locrho,locfu,locfv,locfw)
         
          f07sub(i,j,k)=oneminusomega*f07sub(i,j,k)+omegaminusone* &
           feq(7)+ equil_pop07(locrho,locfu,locfv,locfw)
        
          f08sub(i,j,k)=oneminusomega*f08sub(i,j,k)+omegaminusone* &
           feq(8)+ equil_pop08(locrho,locfu,locfv,locfw)
       
          f09sub(i,j,k)=oneminusomega*f09sub(i,j,k)+omegaminusone* &
           feq(9)+ equil_pop09(locrho,locfu,locfv,locfw)
        
          f10sub(i,j,k)=oneminusomega*f10sub(i,j,k)+omegaminusone* &
           feq(10)+ equil_pop10(locrho,locfu,locfv,locfw)
        
          f11sub(i,j,k)=oneminusomega*f11sub(i,j,k)+omegaminusone* &
           feq(11)+ equil_pop11(locrho,locfu,locfv,locfw)
        
          f12sub(i,j,k)=oneminusomega*f12sub(i,j,k)+omegaminusone* &
           feq(12)+ equil_pop12(locrho,locfu,locfv,locfw)
        
          f13sub(i,j,k)=oneminusomega*f13sub(i,j,k)+omegaminusone* &
           feq(13)+ equil_pop13(locrho,locfu,locfv,locfw)
        
          f14sub(i,j,k)=oneminusomega*f14sub(i,j,k)+omegaminusone* &
           feq(14)+ equil_pop14(locrho,locfu,locfv,locfw)
        
          f15sub(i,j,k)=oneminusomega*f15sub(i,j,k)+omegaminusone* &
           feq(15)+ equil_pop15(locrho,locfu,locfv,locfw)
          
          f16sub(i,j,k)=oneminusomega*f16sub(i,j,k)+omegaminusone* &
           feq(16)+ equil_pop16(locrho,locfu,locfv,locfw)
          
          f17sub(i,j,k)=oneminusomega*f17sub(i,j,k)+omegaminusone* &
           feq(17)+ equil_pop17(locrho,locfu,locfv,locfw)
          
          f18sub(i,j,k)=oneminusomega*f18sub(i,j,k)+omegaminusone* &
           feq(18)+ equil_pop18(locrho,locfu,locfv,locfw)
#else

          f00sub(i,j,k)=oneminusomega*f00sub(i,j,k)+omegaminusone* &
           equil_pop00(locrho,locu,locv,locw)+ &
           equil_pop00(locrho,locfu,locfv,locfw)
       
          f01sub(i,j,k)=oneminusomega*f01sub(i,j,k)+omegaminusone* &
           equil_pop01(locrho,locu,locv,locw)+ &
           equil_pop01(locrho,locfu,locfv,locfw)
    
          f02sub(i,j,k)=oneminusomega*f02sub(i,j,k)+omegaminusone* &
           equil_pop02(locrho,locu,locv,locw)+ &
           equil_pop02(locrho,locfu,locfv,locfw)
    
          f03sub(i,j,k)=oneminusomega*f03sub(i,j,k)+omegaminusone* &
           equil_pop03(locrho,locu,locv,locw)+ &
           equil_pop03(locrho,locfu,locfv,locfw)
       
          f04sub(i,j,k)=oneminusomega*f04sub(i,j,k)+omegaminusone* &
           equil_pop04(locrho,locu,locv,locw)+ &
           equil_pop04(locrho,locfu,locfv,locfw)
        
          f05sub(i,j,k)=oneminusomega*f05sub(i,j,k)+omegaminusone* &
           equil_pop05(locrho,locu,locv,locw)+ &
           equil_pop05(locrho,locfu,locfv,locfw)
        
          f06sub(i,j,k)=oneminusomega*f06sub(i,j,k)+omegaminusone* &
           equil_pop06(locrho,locu,locv,locw)+ &
           equil_pop06(locrho,locfu,locfv,locfw)
         
          f07sub(i,j,k)=oneminusomega*f07sub(i,j,k)+omegaminusone* &
           equil_pop07(locrho,locu,locv,locw)+ &
           equil_pop07(locrho,locfu,locfv,locfw)
        
          f08sub(i,j,k)=oneminusomega*f08sub(i,j,k)+omegaminusone* &
           equil_pop08(locrho,locu,locv,locw)+ &
           equil_pop08(locrho,locfu,locfv,locfw)
       
          f09sub(i,j,k)=oneminusomega*f09sub(i,j,k)+omegaminusone* &
           equil_pop09(locrho,locu,locv,locw)+ &
           equil_pop09(locrho,locfu,locfv,locfw)
        
          f10sub(i,j,k)=oneminusomega*f10sub(i,j,k)+omegaminusone* &
           equil_pop10(locrho,locu,locv,locw)+ &
           equil_pop10(locrho,locfu,locfv,locfw)
        
          f11sub(i,j,k)=oneminusomega*f11sub(i,j,k)+omegaminusone* &
           equil_pop11(locrho,locu,locv,locw)+ &
           equil_pop11(locrho,locfu,locfv,locfw)
        
          f12sub(i,j,k)=oneminusomega*f12sub(i,j,k)+omegaminusone* &
           equil_pop12(locrho,locu,locv,locw)+ &
           equil_pop12(locrho,locfu,locfv,locfw)
        
          f13sub(i,j,k)=oneminusomega*f13sub(i,j,k)+omegaminusone* &
           equil_pop13(locrho,locu,locv,locw)+ &
           equil_pop13(locrho,locfu,locfv,locfw)
        
          f14sub(i,j,k)=oneminusomega*f14sub(i,j,k)+omegaminusone* &
           equil_pop14(locrho,locu,locv,locw)+ &
           equil_pop14(locrho,locfu,locfv,locfw)
        
          f15sub(i,j,k)=oneminusomega*f15sub(i,j,k)+omegaminusone* &
           equil_pop15(locrho,locu,locv,locw)+ &
           equil_pop15(locrho,locfu,locfv,locfw)
          
          f16sub(i,j,k)=oneminusomega*f16sub(i,j,k)+omegaminusone* &
           equil_pop16(locrho,locu,locv,locw)+ &
           equil_pop16(locrho,locfu,locfv,locfw)
          
          f17sub(i,j,k)=oneminusomega*f17sub(i,j,k)+omegaminusone* &
           equil_pop17(locrho,locu,locv,locw)+ &
           equil_pop17(locrho,locfu,locfv,locfw)
          
          f18sub(i,j,k)=oneminusomega*f18sub(i,j,k)+omegaminusone* &
           equil_pop18(locrho,locu,locv,locw)+ &
           equil_pop18(locrho,locfu,locfv,locfw)
           
#endif
           
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
        
#ifdef SCANDREA
        locfu = fusub(i,j,k)*t_LB + locu
        locfv = fvsub(i,j,k)*t_LB + locv
        locfw = fwsub(i,j,k)*t_LB + locw
#else
        locfu = fusub(i,j,k)*t_LB / locrho + locu
        locfv = fvsub(i,j,k)*t_LB / locrho + locv
        locfw = fwsub(i,j,k)*t_LB / locrho + locw
#endif
        fusub(i,j,k) = locfu
        fvsub(i,j,k) = locfv
        fwsub(i,j,k) = locfw
        
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
          fneq = aoptp(l)%p(i,j,k) - feq(l)
          
	      !non equilibrium part of the momentum flux tensor
          pxx=pxx + (dex(l)*dex(l) - cssq)*fneq
          pyy=pyy + (dey(l)*dey(l) - cssq)*fneq
          pzz=pzz + (dez(l)*dez(l) - cssq)*fneq
          pxy=pxy +  dex(l)*dey(l)*fneq
          pxz=pxz +  dex(l)*dez(l)*fneq
          pyz=pyz +  dey(l)*dez(l)*fneq
        enddo
        
        do l=0,links
          aoptp(l)%p(i,j,k)= feq(l) + &
           ((HALF*p(l))/(cssq**TWO))*((dex(l)*dex(l) - cssq)*pxx + &
           (dey(l)*dey(l) - cssq)*pyy + &
           (dez(l)*dez(l) - cssq)*pzz + &
           TWO*dex(l)*dey(l)*pxy + &
           TWO*dex(l)*dez(l)*pxz + &
           TWO*dey(l)*dez(l)*pyz)           
        enddo
           
        f00sub(i,j,k)=oneminusomega*f00sub(i,j,k)+omegaminusone* &
         feq(0)+ equil_pop00(locrho,locfu,locfv,locfw)
        
        f01sub(i,j,k)=oneminusomega*f01sub(i,j,k)+omegaminusone* &
         feq(1)+ equil_pop01(locrho,locfu,locfv,locfw)
  
        f02sub(i,j,k)=oneminusomega*f02sub(i,j,k)+omegaminusone* &
         feq(2)+ equil_pop02(locrho,locfu,locfv,locfw)
    
        f03sub(i,j,k)=oneminusomega*f03sub(i,j,k)+omegaminusone* &
         feq(3)+ equil_pop03(locrho,locfu,locfv,locfw)
       
        f04sub(i,j,k)=oneminusomega*f04sub(i,j,k)+omegaminusone* &
         feq(4)+ equil_pop04(locrho,locfu,locfv,locfw)
       
        f05sub(i,j,k)=oneminusomega*f05sub(i,j,k)+omegaminusone* &
         feq(5)+ equil_pop05(locrho,locfu,locfv,locfw)
        
        f06sub(i,j,k)=oneminusomega*f06sub(i,j,k)+omegaminusone* &
         feq(6)+ equil_pop06(locrho,locfu,locfv,locfw)
         
        f07sub(i,j,k)=oneminusomega*f07sub(i,j,k)+omegaminusone* &
         feq(7)+ equil_pop07(locrho,locfu,locfv,locfw)
      
        f08sub(i,j,k)=oneminusomega*f08sub(i,j,k)+omegaminusone* &
         feq(8)+ equil_pop08(locrho,locfu,locfv,locfw)
       
        f09sub(i,j,k)=oneminusomega*f09sub(i,j,k)+omegaminusone* &
         feq(9)+ equil_pop09(locrho,locfu,locfv,locfw)
        
        f10sub(i,j,k)=oneminusomega*f10sub(i,j,k)+omegaminusone* &
         feq(10)+ equil_pop10(locrho,locfu,locfv,locfw)
        
        f11sub(i,j,k)=oneminusomega*f11sub(i,j,k)+omegaminusone* &
         feq(11)+ equil_pop11(locrho,locfu,locfv,locfw)
        
        f12sub(i,j,k)=oneminusomega*f12sub(i,j,k)+omegaminusone* &
         feq(12)+ equil_pop12(locrho,locfu,locfv,locfw)
        
        f13sub(i,j,k)=oneminusomega*f13sub(i,j,k)+omegaminusone* &
         feq(13)+ equil_pop13(locrho,locfu,locfv,locfw)
        
        f14sub(i,j,k)=oneminusomega*f14sub(i,j,k)+omegaminusone* &
         feq(14)+ equil_pop14(locrho,locfu,locfv,locfw)
        
        f15sub(i,j,k)=oneminusomega*f15sub(i,j,k)+omegaminusone* &
         feq(15)+ equil_pop15(locrho,locfu,locfv,locfw)
          
        f16sub(i,j,k)=oneminusomega*f16sub(i,j,k)+omegaminusone* &
         feq(16)+ equil_pop16(locrho,locfu,locfv,locfw)
          
        f17sub(i,j,k)=oneminusomega*f17sub(i,j,k)+omegaminusone* &
         feq(17)+ equil_pop17(locrho,locfu,locfv,locfw)
          
        f18sub(i,j,k)=oneminusomega*f18sub(i,j,k)+omegaminusone* &
         feq(18)+ equil_pop18(locrho,locfu,locfv,locfw)
#else

        f00sub(i,j,k)=oneminusomega*f00sub(i,j,k)+omegaminusone* &
         equil_pop00(locrho,locu,locv,locw)+ &
         equil_pop00(locrho,locfu,locfv,locfw)
       
        f01sub(i,j,k)=oneminusomega*f01sub(i,j,k)+omegaminusone* &
         equil_pop01(locrho,locu,locv,locw)+ &
         equil_pop01(locrho,locfu,locfv,locfw)
    
        f02sub(i,j,k)=oneminusomega*f02sub(i,j,k)+omegaminusone* &
         equil_pop02(locrho,locu,locv,locw)+ &
         equil_pop02(locrho,locfu,locfv,locfw)
     
        f03sub(i,j,k)=oneminusomega*f03sub(i,j,k)+omegaminusone* &
         equil_pop03(locrho,locu,locv,locw)+ &
         equil_pop03(locrho,locfu,locfv,locfw)
        
        f04sub(i,j,k)=oneminusomega*f04sub(i,j,k)+omegaminusone* &
         equil_pop04(locrho,locu,locv,locw)+ &
         equil_pop04(locrho,locfu,locfv,locfw)
        
        f05sub(i,j,k)=oneminusomega*f05sub(i,j,k)+omegaminusone* &
         equil_pop05(locrho,locu,locv,locw)+ &
         equil_pop05(locrho,locfu,locfv,locfw)
        
        f06sub(i,j,k)=oneminusomega*f06sub(i,j,k)+omegaminusone* &
         equil_pop06(locrho,locu,locv,locw)+ &
         equil_pop06(locrho,locfu,locfv,locfw)
        
        f07sub(i,j,k)=oneminusomega*f07sub(i,j,k)+omegaminusone* &
         equil_pop07(locrho,locu,locv,locw)+ &
         equil_pop07(locrho,locfu,locfv,locfw)
        
        f08sub(i,j,k)=oneminusomega*f08sub(i,j,k)+omegaminusone* &
         equil_pop08(locrho,locu,locv,locw)+ &
         equil_pop08(locrho,locfu,locfv,locfw)
        
        f09sub(i,j,k)=oneminusomega*f09sub(i,j,k)+omegaminusone* &
         equil_pop09(locrho,locu,locv,locw)+ &
         equil_pop09(locrho,locfu,locfv,locfw)
        
        f10sub(i,j,k)=oneminusomega*f10sub(i,j,k)+omegaminusone* &
         equil_pop10(locrho,locu,locv,locw)+ &
         equil_pop10(locrho,locfu,locfv,locfw)
        
        f11sub(i,j,k)=oneminusomega*f11sub(i,j,k)+omegaminusone* &
         equil_pop11(locrho,locu,locv,locw)+ &
         equil_pop11(locrho,locfu,locfv,locfw)
        
        f12sub(i,j,k)=oneminusomega*f12sub(i,j,k)+omegaminusone* &
         equil_pop12(locrho,locu,locv,locw)+ &
         equil_pop12(locrho,locfu,locfv,locfw)
        
        f13sub(i,j,k)=oneminusomega*f13sub(i,j,k)+omegaminusone* &
         equil_pop13(locrho,locu,locv,locw)+ &
         equil_pop13(locrho,locfu,locfv,locfw)
        
        f14sub(i,j,k)=oneminusomega*f14sub(i,j,k)+omegaminusone* &
         equil_pop14(locrho,locu,locv,locw)+ &
         equil_pop14(locrho,locfu,locfv,locfw)
        
        f15sub(i,j,k)=oneminusomega*f15sub(i,j,k)+omegaminusone* &
         equil_pop15(locrho,locu,locv,locw)+ &
         equil_pop15(locrho,locfu,locfv,locfw)
        
        f16sub(i,j,k)=oneminusomega*f16sub(i,j,k)+omegaminusone* &
         equil_pop16(locrho,locu,locv,locw)+ &
         equil_pop16(locrho,locfu,locfv,locfw)
        
        f17sub(i,j,k)=oneminusomega*f17sub(i,j,k)+omegaminusone* &
         equil_pop17(locrho,locu,locv,locw)+ &
         equil_pop17(locrho,locfu,locfv,locfw)
        
        f18sub(i,j,k)=oneminusomega*f18sub(i,j,k)+omegaminusone* &
         equil_pop18(locrho,locu,locv,locw)+ &
         equil_pop18(locrho,locfu,locfv,locfw)
         
#endif

       enddo
     enddo
  enddo
  ! end forall
  
  return
  
 end subroutine collision_fluids_EDM
 
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
#if 1
           locrhoR=rhoR(i,j,k)
           locrhoB=rhoB(i,j,k)
           omega(i,j,k)=viscosity_to_omega(ONE/ &
            ((locrhoR/meanR)/((locrhoR/meanR)+(locrhoB/meanB))*(ONE/viscR) + &
            (locrhoB/meanB)/((locrhoR/meanR)+(locrhoB/meanB))*(ONE/viscB)))
#else
          phis=(rhoR(i,j,k)/meanR-rhoB(i,j,k)/meanB)/ &
           (rhoR(i,j,k)/meanR+rhoB(i,j,k)/meanB)
          omega(i,j,k)=ONE/(HALF*(ONE+phis)*tauR+HALF*(ONE-phis)*TauB)
#endif
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
  
  real(kind=PRC) :: phis
  real(kind=PRC), parameter :: gradlim = real(1.d-16,kind=PRC)
  real(kind=PRC), parameter :: phislim = ONE - gradlim
  real(kind=PRC), parameter :: maxvel = real(0.2d0,kind=PRC)
  
  real(kind=PRC) :: psinorm,psinorm_sq,e_dot_psi,temp,rhoapp
  real(kind=PRC) :: psix,psiy,psiz,cosphi,fsum,acoeff
  real(kind=PRC), dimension(0:links) :: feq
  real(kind=PRC) :: oneminusomega
  real(kind=PRC) :: locrhoR,locrhoB,locu,locv,locw
  integer :: i,j,k,l
  
  real(kind=PRC) :: temp_omega,grad_rhoRx,grad_rhoRy,grad_rhoRz, &
   grad_rhoBx,grad_rhoBy,grad_rhoBz
   
  !real(kind=PRC) :: myalpha
  
#ifdef GRADIENTD3Q27
  real(kind=PRC), dimension(0:linksd3q27) :: rhodiff,rhosum
#else
  real(kind=PRC), dimension(0:links) :: rhodiff,rhosum
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
  
          locrhoR = rhoR(i,j,k)
          locrhoB = rhoB(i,j,k)
          locu = u(i,j,k)
          locv = v(i,j,k)
          locw = w(i,j,k)
          
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
            feqR(l) = &
             equil_popCG(l,temp_omega,alphaR_CG,locrhoR, &
             locu,locv,locw,grad_rhoRx,grad_rhoRy,grad_rhoRz)
            feqB(l) = &
             equil_popCG(l,temp_omega,alphaB_CG,locrhoB, &
             locu,locv,locw,grad_rhoBx,grad_rhoBy,grad_rhoBz)
             
            fneqR = aoptpR(l)%p(i,j,k) - feqR(l)
            fneqB = aoptpB(l)%p(i,j,k) - feqB(l)
			
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
            aoptpR(l)%p(i,j,k)= feqR(l) + &
             ((HALF*p(l))/(cssq**TWO))*((dex(l)*dex(l) - cssq)*pxxR + &
             (dey(l)*dey(l) - cssq)*pyyR + &
             (dez(l)*dez(l) - cssq)*pzzR + &
             TWO*dex(l)*dey(l)*pxyR + &
             TWO*dex(l)*dez(l)*pxzR + &
             TWO*dey(l)*dez(l)*pyzR)
            
            aoptpB(l)%p(i,j,k)= feqB(l) + &
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
            aoptpR(l)%p(i,j,k)=aoptpR(l)%p(i,j,k)*oneminusomega + &
             equil_popCG(l,temp_omega,alphaR_CG,locrhoR, &
             locu,locv,locw,grad_rhoRx,grad_rhoRy,grad_rhoRz) &
             * temp_omega
            aoptpB(l)%p(i,j,k)=aoptpB(l)%p(i,j,k)*oneminusomega + &
             equil_popCG(l,temp_omega,alphaB_CG,locrhoB, &
             locu,locv,locw,grad_rhoBx,grad_rhoBy,grad_rhoBz) &
             * temp_omega
          enddo
          
          
          !perturbation step
          phis=(locrhoR/meanR-locrhoB/meanB)/(locrhoR/meanR+locrhoB/meanB)
          !if(abs(phis)<phislim)then
            rhodiff(1:links) = ZERO
            rhosum(1:links) = ZERO
            psix = ZERO
            psiy = ZERO
            psiz = ZERO
#ifdef GRADIENTD3Q27
            do l=1,linksd3q27
              rhodiff(l)=(rhoR(i+exd3q27(l),j+eyd3q27(l),k+ezd3q27(l))/meanR - &
               rhoB(i+exd3q27(l),j+eyd3q27(l),k+ezd3q27(l))/meanB)
              rhosum(l)= (rhoR(i+exd3q27(l),j+eyd3q27(l),k+ezd3q27(l))/meanR + &
               rhoB(i+exd3q27(l),j+eyd3q27(l),k+ezd3q27(l))/meanB)
            enddo
            do l=1,linksd3q27
              psix=psix + ad3q27(l)*dexd3q27(l)*(rhodiff(l)/rhosum(l))
              psiy=psiy + ad3q27(l)*deyd3q27(l)*(rhodiff(l)/rhosum(l))
              psiz=psiz + ad3q27(l)*dezd3q27(l)*(rhodiff(l)/rhosum(l))
            enddo
#else             
            do l=1,links
              rhodiff(l)=(rhoR(i+ex(l),j+ey(l),k+ez(l))/meanR- &
               rhoB(i+ex(l),j+ey(l),k+ez(l))/meanB)
              rhosum(l)= (rhoR(i+ex(l),j+ey(l),k+ez(l))/meanR+ &
               rhoB(i+ex(l),j+ey(l),k+ez(l))/meanB)
            enddo
            do l=1,links
              psix=psix + a(l)*dex(l)*(rhodiff(l)/rhosum(l))
              psiy=psiy + a(l)*dey(l)*(rhodiff(l)/rhosum(l))
              psiz=psiz + a(l)*dez(l)*(rhodiff(l)/rhosum(l))
            enddo
#endif
            psinorm_sq = psix**TWO + psiy**TWO + psiz**TWO
            if(psinorm_sq<gradlim**TWO)cycle
            psinorm=sqrt(psinorm_sq)
            acoeff=( NINE / FOUR )*temp_omega*sigma_CG
            do l=0,links
              e_dot_psi=dex(l)*psix + dey(l)*psiy + dez(l)*psiz
              temp=psinorm*(p(l)*(e_dot_psi**TWO)/psinorm_sq - b_l(l))
              if(isnan(temp)) temp=ZERO
              aoptpR(l)%p(i,j,k)=aoptpR(l)%p(i,j,k) +(HALF*acoeff)*temp 
              aoptpB(l)%p(i,j,k)=aoptpB(l)%p(i,j,k) +(HALF*acoeff)*temp 
            enddo
            
            !recolouring step
            do l=0,links
              !compute the sum of the two equilibriums at zero velocity 
              !feq(l)=locrhoR*phiR_CG(l) + locrhoB*phib_CG(l)
              feq(l)=equil_popCG(l,temp_omega,alphaR_CG,locrhoR, &
             ZERO,ZERO,ZERO,grad_rhoRx,grad_rhoRy,grad_rhoRz)+ &
             equil_popCG(l,temp_omega,alphaB_CG,locrhoB, &
             ZERO,ZERO,ZERO,grad_rhoBx,grad_rhoBy,grad_rhoBz)
            enddo
            rhoapp=locrhoR+locrhoB
            do l=0,links
              e_dot_psi=dex(l)*psix + dey(l)*psiy + dez(l)*psiz
              temp=sqrt(dex(l)**TWO + dey(l)**TWO + dez(l)**TWO)*psinorm
              if (temp>ZERO) then
                 cosphi=e_dot_psi/temp
              else
                 cosphi=ZERO
              endif
              temp=beta_CG*locrhoR*locrhoB*cosphi/(rhoapp**TWO)
              fsum=aoptpR(l)%p(i,j,k) + aoptpB(l)%p(i,j,k)
              aoptpR(l)%p(i,j,k)=fsum*locrhoR/rhoapp + temp*feq(l)
              aoptpB(l)%p(i,j,k)=fsum*locrhoB/rhoapp - temp*feq(l)
            enddo
            
          !endif
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
  
  real(kind=PRC) :: phis
  real(kind=PRC), parameter :: gradlim = real(1.d-16,kind=PRC)
  real(kind=PRC), parameter :: phislim = ONE - gradlim
  real(kind=PRC), parameter :: maxvel = real(0.2d0,kind=PRC)
  
  real(kind=PRC) :: psinorm,psinorm_sq,e_dot_psi,temp,rhoapp
  real(kind=PRC) :: psix,psiy,psiz,cosphi,fsum,acoeff
  real(kind=PRC), dimension(0:links) :: feq
  real(kind=PRC) :: oneminusomega
  real(kind=PRC) :: locrhoR,locrhoB,locu,locv,locw
  integer :: i,j,k,l
  
  real(kind=PRC) :: temp_omega,grad_rhoRx,grad_rhoRy,grad_rhoRz, &
   grad_rhoBx,grad_rhoBy,grad_rhoBz
   
  !real(kind=PRC) :: myalpha
  
#ifdef GRADIENTD3Q27
  real(kind=PRC), dimension(0:linksd3q27) :: rhodiff,rhosum
#else
  real(kind=PRC), dimension(0:links) :: rhodiff,rhosum
#endif
  
#ifdef REGULARIZED
  real(kind=PRC) :: pxxR,pyyR,pzzR,pxyR,pxzR,pyzR
  real(kind=PRC) :: pxxB,pyyB,pzzB,pxyB,pxzB,pyzB
  real(kind=PRC), dimension(0:links) :: feqR,feqB
  real(kind=PRC) :: fneqR,fneqB
#endif
  
  
  
   do k=minz,maxz
     do j=miny,maxy
       do i=minx,maxx
          
          temp_omega=omega(i,j,k)
          oneminusomega = ONE-temp_omega
          
          locrhoR = rhoR(i,j,k)
          locrhoB = rhoB(i,j,k)
          locu = u(i,j,k)
          locv = v(i,j,k)
          locw = w(i,j,k)
          
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
            feqR(l) = &
             equil_popCG(l,temp_omega,alphaR_CG,locrhoR, &
             locu,locv,locw,grad_rhoRx,grad_rhoRy,grad_rhoRz)
            feqB(l) = &
             equil_popCG(l,temp_omega,alphaB_CG,locrhoB, &
             locu,locv,locw,grad_rhoBx,grad_rhoBy,grad_rhoBz)
             
            fneqR = aoptpR(l)%p(i,j,k) - feqR(l)
            fneqB = aoptpB(l)%p(i,j,k) - feqB(l)
			
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
            aoptpR(l)%p(i,j,k)= feqR(l) + &
             ((HALF*p(l))/(cssq**TWO))*((dex(l)*dex(l) - cssq)*pxxR + &
             (dey(l)*dey(l) - cssq)*pyyR + &
             (dez(l)*dez(l) - cssq)*pzzR + &
             TWO*dex(l)*dey(l)*pxyR + &
             TWO*dex(l)*dez(l)*pxzR + &
             TWO*dey(l)*dez(l)*pyzR)
            
            aoptpB(l)%p(i,j,k)= feqB(l) + &
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
            aoptpR(l)%p(i,j,k)=aoptpR(l)%p(i,j,k)*oneminusomega + &
             equil_popCG(l,temp_omega,alphaR_CG,locrhoR, &
             locu,locv,locw,grad_rhoRx,grad_rhoRy,grad_rhoRz) &
             * temp_omega
            aoptpB(l)%p(i,j,k)=aoptpB(l)%p(i,j,k)*oneminusomega + &
             equil_popCG(l,temp_omega,alphaB_CG,locrhoB, &
             locu,locv,locw,grad_rhoBx,grad_rhoBy,grad_rhoBz) &
             * temp_omega
          enddo
          
          
          !perturbation step
          phis=(locrhoR/meanR-locrhoB/meanB)/(locrhoR/meanR+locrhoB/meanB)
          !if(abs(phis)<phislim)then
            rhodiff(1:links) = ZERO
            rhosum(1:links) = ZERO
            psix = ZERO
            psiy = ZERO
            psiz = ZERO
#ifdef GRADIENTD3Q27
            do l=1,linksd3q27
              rhodiff(l)=(rhoR(i+exd3q27(l),j+eyd3q27(l),k+ezd3q27(l))/meanR - &
               rhoB(i+exd3q27(l),j+eyd3q27(l),k+ezd3q27(l))/meanB)
              rhosum(l)= (rhoR(i+exd3q27(l),j+eyd3q27(l),k+ezd3q27(l))/meanR + &
               rhoB(i+exd3q27(l),j+eyd3q27(l),k+ezd3q27(l))/meanB)
            enddo
            do l=1,linksd3q27
              psix=psix + ad3q27(l)*dexd3q27(l)*(rhodiff(l)/rhosum(l))
              psiy=psiy + ad3q27(l)*deyd3q27(l)*(rhodiff(l)/rhosum(l))
              psiz=psiz + ad3q27(l)*dezd3q27(l)*(rhodiff(l)/rhosum(l))
            enddo
#else             
            do l=1,links
              rhodiff(l)=(rhoR(i+ex(l),j+ey(l),k+ez(l))/meanR- &
               rhoB(i+ex(l),j+ey(l),k+ez(l))/meanB)
              rhosum(l)= (rhoR(i+ex(l),j+ey(l),k+ez(l))/meanR+ &
               rhoB(i+ex(l),j+ey(l),k+ez(l))/meanB)
            enddo
            do l=1,links
              psix=psix + a(l)*dex(l)*(rhodiff(l)/rhosum(l))
              psiy=psiy + a(l)*dey(l)*(rhodiff(l)/rhosum(l))
              psiz=psiz + a(l)*dez(l)*(rhodiff(l)/rhosum(l))
            enddo
#endif
            psinorm_sq = psix**TWO + psiy**TWO + psiz**TWO
            if(psinorm_sq<gradlim**TWO)cycle
            psinorm=sqrt(psinorm_sq)
            acoeff=( NINE / FOUR )*temp_omega*sigma_CG
            do l=0,links
              e_dot_psi=dex(l)*psix + dey(l)*psiy + dez(l)*psiz
              temp=psinorm*(p(l)*(e_dot_psi**TWO)/psinorm_sq - b_l(l))
              if(isnan(temp)) temp=ZERO
              aoptpR(l)%p(i,j,k)=aoptpR(l)%p(i,j,k) +(HALF*acoeff)*temp 
              aoptpB(l)%p(i,j,k)=aoptpB(l)%p(i,j,k) +(HALF*acoeff)*temp 
            enddo
            
            !recolouring step
            do l=0,links
              !compute the sum of the two equilibriums at zero velocity 
              !feq(l)=locrhoR*phiR_CG(l) + locrhoB*phib_CG(l)
              feq(l)=equil_popCG(l,temp_omega,alphaR_CG,locrhoR, &
             ZERO,ZERO,ZERO,grad_rhoRx,grad_rhoRy,grad_rhoRz)+ &
             equil_popCG(l,temp_omega,alphaB_CG,locrhoB, &
             ZERO,ZERO,ZERO,grad_rhoBx,grad_rhoBy,grad_rhoBz)
            enddo
            rhoapp=locrhoR+locrhoB
            do l=0,links
              e_dot_psi=dex(l)*psix + dey(l)*psiy + dez(l)*psiz
              temp=sqrt(dex(l)**TWO + dey(l)**TWO + dez(l)**TWO)*psinorm
              if (temp>ZERO) then
                 cosphi=e_dot_psi/temp
              else
                 cosphi=ZERO
              endif
              temp=beta_CG*locrhoR*locrhoB*cosphi/(rhoapp**TWO)
              fsum=aoptpR(l)%p(i,j,k) + aoptpB(l)%p(i,j,k)
              aoptpR(l)%p(i,j,k)=fsum*locrhoR/rhoapp + temp*feq(l)
              aoptpB(l)%p(i,j,k)=fsum*locrhoB/rhoapp - temp*feq(l)
            enddo
            
          !endif
        enddo
      enddo
    enddo
  
  return
  
 end subroutine collision_fluids_CG
 
 pure function equil_popCG(myl,myomega,myalpha, &
  myrho,myu,myv,myw,gradx,grady,gradz)
 
!***********************************************************************
!     
!     LBsoft function for Boltzmann equilibrium distribution 
!     for the populations in the Colour Gradient model
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification November 2019
!     
!***********************************************************************
 
  implicit none
  
  integer, parameter :: mylinks=links
  real(kind=PRC), parameter :: mycssq = cssq
  real(kind=PRC), parameter, dimension(0:mylinks) :: myp=p
  real(kind=PRC), parameter, dimension(0:mylinks) :: mydex=dex
  real(kind=PRC), parameter, dimension(0:mylinks) :: mydey=dey
  real(kind=PRC), parameter, dimension(0:mylinks) :: mydez=dez
  real(kind=PRC), parameter, dimension(0:mylinks) :: mypsi=psi_CG
  real(kind=PRC), parameter, dimension(0:mylinks) :: myxi=xi_CG
  real(kind=PRC), parameter, dimension(0:mylinks) :: myphi=phi_CG
  real(kind=PRC), parameter, dimension(0:mylinks) :: myvarphi=varphi_CG
  
  integer, intent(in) :: myl
  
  real(kind=PRC), intent(in) :: myalpha,myrho,myu,myv,myw, &
   gradx,grady,gradz,myomega
  
  real(kind=PRC) :: equil_popCG
  
  real(kind=PRC) :: uv,mycontraction,dotugrad,myviscos
  real(kind=PRC), dimension(3) :: myvel,mygrad
  real(kind=PRC), dimension(3,3) :: myG
  
  
  
  
  integer :: i,j
  
  
#if 1
  myvel(1)=myu
  myvel(2)=myv
  myvel(3)=myw
  mygrad(1)=gradx
  mygrad(2)=grady
  mygrad(3)=gradz
  do j=1,3
   do i=1,3
     myG(i,j)=myvel(i)*mygrad(j)+myvel(j)*mygrad(i)                          
   enddo
  enddo
  mycontraction=mydex(myl)*mydex(myl)*myG(1,1) + &
           mydey(myl)*mydey(myl)*myG(2,2) + &
           mydez(myl)*mydez(myl)*myG(3,3) + &
           mydex(myl)*mydey(myl)*myG(1,2) + &
           mydex(myl)*mydey(myl)*myG(2,1) + &
           mydex(myl)*mydez(myl)*myG(1,3) + &
           mydex(myl)*mydez(myl)*myG(3,1) + &
           mydey(myl)*mydez(myl)*myG(2,3) + &
           mydey(myl)*mydez(myl)*myG(3,2)
  dotugrad=dot_product(myvel,mygrad)
  myviscos=mycssq*(ONE / myomega - HALF)
  uv=(ONE/mycssq)*(myu*mydex(myl) + myv*mydey(myl) + myw*mydez(myl))
  equil_popCG=myrho*(myphi(myl)+myvarphi(myl)*myalpha+myp(myl)*(uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))) + &
   myviscos*(mypsi(myl)*dotugrad + myxi(myl)*mycontraction)
#else
  uv=(ONE/mycssq)*(myu*mydex(myl) + myv*mydey(myl) + myw*mydez(myl))
  equil_popCG=myrho*(myphi(myl)+myvarphi(myl)*myalpha+myp(myl)*(uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO)))
#endif
  
  return
  
 end function equil_popCG
 
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
  
  call manage_bc_pop_selfcomm(aoptpR,lparticles)
  
  if(lsingle_fluid)return
  
  call manage_bc_pop_selfcomm(aoptpB,lparticles)
 
 end subroutine driver_bc_pop_selfcomm
 


  subroutine stream_nocopy(f01sub,&
  f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub,f09sub,f10sub, &
  f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub)
    
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
  real(kind=PRC), allocatable, dimension(:,:,:)  :: f01sub,f02sub,f03sub,f04sub, &
   f05sub,f06sub,f07sub,f08sub,f09sub,f10sub,f11sub,f12sub,f13sub, &
   f14sub,f15sub,f16sub,f17sub,f18sub
  integer :: i,j,k
  logical, save :: isFirst = .true.


  if (isFirst) then
        isFirst = .false.
        if (idrank == 0) write(6,*) "STREAM: Using stream_nocopy"
  endif


      do i=maxx+1, minx, -1
         f01sub(i,:,:) =  f01sub(i-1,:,:)
      enddo
      do i=minx-1, maxx
         f02sub(i,:,:) =  f02sub(i+1,:,:)
      enddo

      do j=maxy+1, miny, -1
         f03sub(:,j,:) =  f03sub(:,j-1,:)
      enddo
      do j=miny-1, maxy
         f04sub(:,j,:) =  f04sub(:,j+1,:)
      enddo

      do k=maxz+1, minz, -1
         f05sub(:,:,k) =  f05sub(:,:,k-1)
      enddo
      do k=minz-1, maxz
         f06sub(:,:,k) =  f06sub(:,:,k+1)
      enddo

      do i=maxx+1, minx, -1
       do j=maxy+1, miny, -1
         f07sub(i,j,:) =  f07sub(i-1,j-1,:)
       enddo
      enddo
      do i=minx-1, maxx
       do j=miny-1, maxy
         f08sub(i,j,:) =  f08sub(i+1,j+1,:)
       enddo
      enddo

      do i=minx-1, maxx
       do j=maxy+1, miny, -1
         f09sub(i,j,:) =  f09sub(i+1,j-1,:)
       enddo
      enddo
      do i=maxx+1, minx, -1
       do j=miny-1, maxy
         f10sub(i,j,:) =  f10sub(i-1,j+1,:)
       enddo
      enddo

      do i=maxx+1, minx, -1
       do k=maxz+1, minz, -1
         f11sub(i,:,k) =  f11sub(i-1,:,k-1)
       enddo
      enddo
      do i=minx-1, maxx
       do k=minz-1, maxz
         f12sub(i,:,k) =  f12sub(i+1,:,k+1)
       enddo
      enddo

      do i=minx-1, maxx
       do k=maxz+1, minz, -1
         f13sub(i,:,k) =  f13sub(i+1,:,k-1)
       enddo
      enddo
      do i=maxx+1, minx, -1
       do k=minz-1, maxz
         f14sub(i,:,k) =  f14sub(i-1,:,k+1)
       enddo
      enddo

      do j=maxy+1, miny, -1
       do k=maxz+1, minz, -1
         f15sub(:,j,k) =  f15sub(:,j-1,k-1)
       enddo
      enddo
      do j=miny-1, maxy
       do k=minz-1, maxz
         f16sub(:,j,k) =  f16sub(:,j+1,k+1)
       enddo
      enddo

      do j=miny-1, maxy
       do k=maxz+1, minz, -1
         f17sub(:,j,k) =  f17sub(:,j+1,k-1)
       enddo
      enddo
      do j=maxy+1, miny, -1
       do k=minz-1, maxz
         f18sub(:,j,k) =  f18sub(:,j-1,k+1)
       enddo
      enddo

     end subroutine stream_nocopy
     
  subroutine stream_nocopy_isfluid(f01sub,&
  f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub,f09sub,f10sub, &
  f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub)
    
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
  real(kind=PRC), allocatable, dimension(:,:,:)  :: f01sub,f02sub,f03sub,f04sub, &
   f05sub,f06sub,f07sub,f08sub,f09sub,f10sub,f11sub,f12sub,f13sub, &
   f14sub,f15sub,f16sub,f17sub,f18sub
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
          f01sub(i,j,k) =  f01sub(i-1,j,k)
        endif
      enddo
    enddo
  enddo
  
  do k=minz-1, maxz+1
    do j=miny-1, maxy+1
      do i=minx-1, maxx
        if ( isfluid(i+1,j,k)<3 .or. isfluid(i+1,j,k)>4) then
          f02sub(i,j,k) =  f02sub(i+1,j,k)
        endif
      enddo
    enddo
  enddo
  
  do k=minz-1, maxz+1
    do j=maxy+1, miny, -1
      do i=minx-1, maxx+1
        if ( isfluid(i,j-1,k)<3 .or. isfluid(i,j-1,k)>4) then
          f03sub(i,j,k) =  f03sub(i,j-1,k)
        endif
      enddo
    enddo
  enddo
  
  do k=minz-1, maxz+1
    do j=miny-1, maxy
      do i=minx-1, maxx+1
        if ( isfluid(i,j+1,k)<3 .or. isfluid(i,j+1,k)>4) then
          f04sub(i,j,k) =  f04sub(i,j+1,k)
        endif
      enddo
    enddo
  enddo
  
  do k=maxz+1, minz, -1
    do j=miny-1, maxy+1
      do i=minx-1, maxx+1
        if ( isfluid(i,j,k-1)<3 .or. isfluid(i,j,k-1)>4) then
          f05sub(i,j,k) =  f05sub(i,j,k-1)
        endif
      enddo
    enddo
  enddo
  
  do k=minz-1, maxz
    do j=miny-1, maxy+1
      do i=minx-1, maxx+1
        if ( isfluid(i,j,k+1)<3 .or. isfluid(i,j,k+1)>4) then
          f06sub(i,j,k) =  f06sub(i,j,k+1)
        endif
      enddo
    enddo
  enddo
  
  do k=minz-1, maxz+1
    do i=maxx+1, minx, -1
      do j=maxy+1, miny, -1
        if ( isfluid(i-1,j-1,k)<3 .or. isfluid(i-1,j-1,k)>4) then
          f07sub(i,j,k) =  f07sub(i-1,j-1,k)
        endif
      enddo
    enddo
  enddo
  
  do k=minz-1, maxz+1
    do i=minx-1, maxx
      do j=miny-1, maxy
        if ( isfluid(i+1,j+1,k)<3 .or. isfluid(i+1,j+1,k)>4) then
          f08sub(i,j,k) =  f08sub(i+1,j+1,k)
        endif
      enddo
    enddo
  enddo
  
  do k=minz-1, maxz+1
    do i=minx-1, maxx
      do j=maxy+1, miny, -1
        if ( isfluid(i+1,j-1,k)<3 .or. isfluid(i+1,j-1,k)>4) then
          f09sub(i,j,k) =  f09sub(i+1,j-1,k)
        endif
      enddo
    enddo
  enddo
   
  do k=minz-1, maxz+1
    do i=maxx+1, minx, -1
      do j=miny-1, maxy
        if ( isfluid(i-1,j+1,k)<3 .or. isfluid(i-1,j+1,k)>4) then
          f10sub(i,j,k) =  f10sub(i-1,j+1,k)
        endif
      enddo
    enddo
  enddo
   
  do k=maxz+1, minz, -1
    do j=miny-1, maxy+1
      do i=maxx+1, minx, -1
        if ( isfluid(i-1,j,k-1)<3 .or. isfluid(i-1,j,k-1)>4) then
          f11sub(i,j,k) =  f11sub(i-1,j,k-1)
        endif
      enddo
    enddo
  enddo
    
  do k=minz-1, maxz
    do j=miny-1, maxy+1
      do i=minx-1, maxx
        if ( isfluid(i+1,j,k+1)<3 .or. isfluid(i+1,j,k+1)>4) then
          f12sub(i,j,k) =  f12sub(i+1,j,k+1)
        endif
      enddo
    enddo
  enddo
    
  do k=maxz+1, minz, -1
    do j=miny-1, maxy+1
      do i=minx-1, maxx
        if ( isfluid(i+1,j,k-1)<3 .or. isfluid(i+1,j,k-1)>4) then
          f13sub(i,j,k) =  f13sub(i+1,j,k-1)
        endif
       enddo
     enddo
  enddo
  
  do k=minz-1, maxz
    do j=miny-1, maxy+1
      do i=maxx+1, minx, -1
        if ( isfluid(i-1,j,k+1)<3 .or. isfluid(i-1,j,k+1)>4) then
          f14sub(i,j,k) =  f14sub(i-1,j,k+1)
        endif
      enddo
    enddo
  enddo
  
  do k=maxz+1, minz, -1
    do j=maxy+1, miny, -1
      do i=minx-1, maxx+1
        if ( isfluid(i,j-1,k-1)<3 .or. isfluid(i,j-1,k-1)>4) then
          f15sub(i,j,k) =  f15sub(i,j-1,k-1)
        endif
      enddo
    enddo
  enddo
      
  do k=minz-1, maxz
    do j=miny-1, maxy
      do i=minx-1, maxx+1
        if ( isfluid(i,j+1,k+1)<3 .or. isfluid(i,j+1,k+1)>4) then
          f16sub(i,j,k) =  f16sub(i,j+1,k+1)
        endif
      enddo
    enddo
  enddo
   
  do k=maxz+1, minz, -1
    do j=miny-1, maxy
      do i=minx-1, maxx+1
        if ( isfluid(i,j+1,k-1)<3 .or. isfluid(i,j+1,k-1)>4) then
          f17sub(i,j,k) =  f17sub(i,j+1,k-1)
        endif
      enddo
    enddo
  enddo
      
  do k=minz-1, maxz
    do j=maxy+1, miny, -1
      do i=minx-1, maxx+1
        if ( isfluid(i,j-1,k+1)<3 .or. isfluid(i,j-1,k+1)>4) then
          f18sub(i,j,k) =  f18sub(i,j-1,k+1)
        endif
      enddo
    enddo
  enddo
  
  return
  
  end subroutine stream_nocopy_isfluid

  subroutine stream_copy(fsub, l)
  implicit none
  real(kind=PRC), allocatable, dimension(:,:,:)  :: fsub
  integer, intent(in) :: l
  integer :: i,j,k,ishift,jshift,kshift
  logical, save :: isFirst = .true.


  if (isFirst) then
        isFirst = .false.
        if (idrank == 0) write(6,*) "STREAM: Using stream_copy"
  endif

    ishift=ex(l)
    jshift=ey(l)
    kshift=ez(l)
    do k=minz-1,maxz+1
     do j=miny-1,maxy+1
      do i=minx-1,maxx+1
    !forall(i=minx-1:maxx+1,j=miny-1:maxy+1,k=minz-1:maxz+1,isfluid(i,j,k)<3 .or. isfluid(i,j,k)>4)
        if ( isfluid(i,j,k)<3 .or. isfluid(i,j,k)>4) then
         buffservice3d(i+ishift,j+jshift,k+kshift) = fsub(i,j,k)
        endif
      enddo
     enddo
    enddo
    !end forall

    do k=minz-1,maxz+1
     do j=miny-1,maxy+1
      do i=minx-1,maxx+1
    !forall(i=minx-1:maxx+1,j=miny-1:maxy+1,k=minz-1:maxz+1,isfluid(i,j,k)<3 .or. isfluid(i,j,k)>4)
        if ( isfluid(i,j,k)<3 .or. isfluid(i,j,k)>4) then
         fsub(i,j,k) = buffservice3d(i,j,k)
        endif
      enddo
     enddo
    enddo
   !end forall
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
       call print_all_pops(100,'mioprima',iter,aoptpR)
     endif
   endif
#endif
    

  call commspop(aoptpR)




#ifdef STREAM_NOCOPY

  if(lparticles)then
  
  call stream_nocopy_isfluid(f01R,&
  f02R,f03R,f04R,f05R,f06R,f07R,f08R,f09R,f10R, &
  f11R,f12R,f13R,f14R,f15R,f16R,f17R,f18R)
  
  else
  
    if(lread_isfluid)then
      call stream_nocopy_isfluid(f01R,&
       f02R,f03R,f04R,f05R,f06R,f07R,f08R,f09R,f10R, &
       f11R,f12R,f13R,f14R,f15R,f16R,f17R,f18R)
    else
      call stream_nocopy(f01R,&
       f02R,f03R,f04R,f05R,f06R,f07R,f08R,f09R,f10R, &
       f11R,f12R,f13R,f14R,f15R,f16R,f17R,f18R)
    endif
  endif
#else
  call stream_copy(f01R, 1)
  call stream_copy(f02R, 2)
  call stream_copy(f03R, 3)
  call stream_copy(f04R, 4)
  call stream_copy(f05R, 5)
  call stream_copy(f06R, 6)
  call stream_copy(f07R, 7)
  call stream_copy(f08R, 8)
  call stream_copy(f09R, 9)
  call stream_copy(f10R, 10)
  call stream_copy(f11R, 11)
  call stream_copy(f12R, 12)
  call stream_copy(f13R, 13)
  call stream_copy(f14R, 14)
  call stream_copy(f15R, 15)
  call stream_copy(f16R, 16)
  call stream_copy(f17R, 17)
  call stream_copy(f18R, 18)
#endif /* STREAM_NOCOPY */





  call manage_bc_pop_selfcomm(aoptpR,lparticles)

  call commrpop(aoptpR,lparticles,isfluid)




#ifdef DIAGNSTREAM
  if(iter.eq.NDIAGNSTREAM) then
    if(allocated(ownern))then
       call print_all_pops(300,'miodopo',iter,aoptpR)
    endif
  endif
#endif


  
  if(lsingle_fluid)return


  call commspop(aoptpB)





#ifdef STREAM_NOCOPY
  if(lparticles)then
  
  call stream_nocopy_isfluid(f01B,&
  f02B,f03B,f04B,f05B,f06B,f07B,f08B,f09B,f10B, &
  f11B,f12B,f13B,f14B,f15B,f16B,f17B,f18B)
  
  else
  if(lread_isfluid)then
      call stream_nocopy_isfluid(f01B,&
       f02B,f03B,f04B,f05B,f06B,f07B,f08B,f09B,f10B, &
       f11B,f12B,f13B,f14B,f15B,f16B,f17B,f18B)
    else
      call stream_nocopy(f01B,&
       f02B,f03B,f04B,f05B,f06B,f07B,f08B,f09B,f10B, &
       f11B,f12B,f13B,f14B,f15B,f16B,f17B,f18B)
    endif
  endif
#else
  call stream_copy(f01B, 1)
  call stream_copy(f02B, 2)
  call stream_copy(f03B, 3)
  call stream_copy(f04B, 4)
  call stream_copy(f05B, 5)
  call stream_copy(f06B, 6)
  call stream_copy(f07B, 7)
  call stream_copy(f08B, 8)
  call stream_copy(f09B, 9)
  call stream_copy(f10B, 10)
  call stream_copy(f11B, 11)
  call stream_copy(f12B, 12)
  call stream_copy(f13B, 13)
  call stream_copy(f14B, 14)
  call stream_copy(f15B, 15)
  call stream_copy(f16B, 16)
  call stream_copy(f17B, 17)
  call stream_copy(f18B, 18)
#endif /* STREAM_NOCOPY */

  

  call manage_bc_pop_selfcomm(aoptpB,lparticles)

  call commrpop(aoptpB,lparticles,isfluid)


  


  
  
  return
  
 end subroutine driver_streaming_fluids
 
 subroutine streaming_fluids(lparticles,f00sub,f01sub,f02sub,f03sub,f04sub, &
   f05sub,f06sub,f07sub,f08sub,f09sub,f10sub,f11sub,f12sub,f13sub, &
   f14sub,f15sub,f16sub,f17sub,f18sub,aoptp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the streaming step
!     on the Boltzmann populations
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  logical, intent(in) :: lparticles
  real(kind=PRC), allocatable, dimension(:,:,:)  :: f00sub,f01sub, &
   f02sub,f03sub,f04sub,f05sub,f06sub,f07sub,f08sub,f09sub,f10sub, &
   f11sub,f12sub,f13sub,f14sub,f15sub,f16sub,f17sub,f18sub
   
  type(REALPTR), dimension(0:links):: aoptp
  
  integer :: i,j,k,l,ishift,jshift,kshift,itemp,jtemp,ktemp
  
  integer, save :: iter=0
  
  
  
  iter=iter+1

#if 0
   if(iter.eq.1) then
   do l=1,links
!      write((50+idrank),*)'in pop ',l
       do i=wminx,wmaxx
         do j=miny,maxy
           do k=wminz,wmaxz
               write((100*idrank)+250+l,*)i,j,k,aoptp(l)%p(i,j,k)
           enddo
         enddo
        enddo
     enddo
   endif
#endif

   call commspop(aoptp)
   

   do l=1,links

      ishift=ex(l)
      jshift=ey(l)
      kshift=ez(l)
      forall(i=wminx:wmaxx,j=wminy:wmaxy,k=wminz:wmaxz) buffservice3d(i,j,k) = -1
!      forall(i=wminx+1:wmaxx-1,j=wminy+1:wmaxy-1,k=wminz+1:wmaxz-1)
      forall(i=minx-1:maxx+1,j=miny-1:maxy+1,k=minz-1:maxz+1) 
         buffservice3d(i+ishift,j+jshift,k+kshift) = aoptp(l)%p(i,j,k)
      end forall
#if 1
      do i=wminx,wmaxx,wmaxx-wminx
         itemp=i+ishift
         if(itemp.le.1.or.itemp.ge.(nx)) then
            if(ixpbc.eq.1) then
               if(itemp.eq.0) then
                  itemp=nx
               endif
               if(itemp.eq.(nx+1)) then
                  itemp=1
               endif
            endif
         endif
!!         if(itemp.GE.wminx.AND.itemp.LE.wmaxx) then
!max            if(l.eq.3.AND.idrank.eq.1) then
!max               write(0,*)'2, wminy+1 ',wminy+1,'wmaxy-1 ',wmaxy-1
!max            endif   
           forall(j=wminy+1:wmaxy-1,k=wminz+1:wmaxz-1) buffservice3d(itemp,j+jshift,k+kshift) = aoptp(l)%p(i,j,k)
!!         endif
         do j=wminy,wmaxy,wmaxy-wminy
            jtemp=j+jshift
            if(jtemp.le.1.or.jtemp.ge.(ny)) then
               if(iypbc.eq.1) then
                  if(jtemp.eq.0) then
                     jtemp=ny
                  endif
                  if(jtemp.eq.(ny+1)) then
                     jtemp=1
                  endif
               endif
            endif
!!            if(itemp.GE.wminx.AND.itemp.LE.wmaxx.AND.jtemp.GE.wminy.AND.jtemp.LE.wmaxy) then
!max               if(l.eq.3.AND.idrank.eq.1) then
!max                  write(0,*)'4, jtemp ',jtemp
!max               endif
               forall(k=wminz+1:wmaxz-1) buffservice3d(itemp,jtemp,k+kshift) = aoptp(l)%p(i,j,k)
!!            endif
         enddo
         do k=wminz,wmaxz,wmaxz-wminz
            ktemp=k+kshift;
            if(ktemp.le.1.or.ktemp.ge.(nz)) then
               if(izpbc.eq.1) then
                  if(ktemp.eq.0) then
                     ktemp=nz
                  endif
                  if(ktemp.eq.(nz+1)) then
                     ktemp=1
                  endif
               endif
            endif
!!            if(itemp.GE.wminx.AND.itemp.LE.wmaxx.AND.ktemp.GE.wminz.AND.ktemp.LE.wmaxz) then
!max               if(l.eq.3.AND.idrank.eq.1) then
!max                     write(0,*)'6, wminy+1 ',wminy+1,'wmaxy-1 ',wmaxy-1
!max               endif
               forall(j=wminy+1:wmaxy-1) buffservice3d(itemp,j+jshift,ktemp) = aoptp(l)%p(i,j,k)
!!            endif
         enddo
      enddo
      do j=wminy,wmaxy,wmaxy-wminy
         jtemp=j+jshift
         if(jtemp.le.1.or.jtemp.ge.(ny)) then
            if(iypbc.eq.1) then
               if(jtemp.eq.0) then
                  jtemp=ny
               endif
               if(jtemp.eq.(ny+1)) then
                  jtemp=1
               endif
            endif
         endif
!!         if(jtemp.GE.wminy.AND.jtemp.LE.wmaxy) then
!max            if(l.eq.3.AND.idrank.eq.1) then
!max               write(0,*)'8, jtemp ',jtemp
!max            endif
            forall(i=wminx+1:wmaxx-1,k=wminz+1:wmaxz-1) buffservice3d(i+ishift,jtemp,k+kshift) = aoptp(l)%p(i,j,k)
!!         endif
         do i=wminx,wmaxx,wmaxx-wminx
            itemp=i+ishift
            if(itemp.le.1.or.itemp.ge.(nx)) then
               if(ixpbc.eq.1) then
                  if(itemp.eq.0) then
                     itemp=nx
                  endif
                  if(itemp.eq.(nx+1)) then
                     itemp=1
                  endif
               endif
            endif
!!            if(itemp.GE.wminx.AND.itemp.LE.wmaxx.AND.jtemp.GE.wminy.AND.jtemp.LE.wmaxy) then
!max               if(l.eq.3.AND.idrank.eq.1) then
!max                  write(0,*)'10, jtemp ',jtemp
!max               endif
               forall(k=wminz+1:wmaxz-1) buffservice3d(itemp,jtemp,k+kshift) = aoptp(l)%p(i,j,k)
!!            endif
         enddo
         do k=wminz,wmaxz,wmaxz-wminz
            ktemp=k+kshift;
            if(ktemp.le.1.or.ktemp.ge.(nz)) then
               if(izpbc.eq.1) then
                  if(ktemp.eq.0) then
                     ktemp=nz
                  endif
                  if(ktemp.eq.(nz+1)) then
                     ktemp=1
                  endif
               endif
            endif
!!            if(jtemp.GE.wminy.AND.jtemp.LE.wmaxy.AND.ktemp.GE.wminz.AND.ktemp.LE.wmaxz) then
!max               if(l.eq.3.AND.idrank.eq.1) then
!max                  write(0,*)'12, jtemp ',jtemp
!max               endif
               forall(i=wminx+1:wmaxx-1) buffservice3d(i+ishift,jtemp,ktemp) = aoptp(l)%p(i,j,k)
!!            endif
         enddo
      enddo
      do k=wminz,wmaxz,wmaxz-wminz
         ktemp=k+kshift
         if(ktemp.le.1.or.ktemp.ge.(nz)) then
            if(izpbc.eq.1) then
               if(ktemp.eq.0) then
                  ktemp=nz
               endif
               if(ktemp.eq.(nz+1)) then
                  ktemp=1
               endif
            endif
         endif
!!         if(ktemp.GE.wminz.AND.ktemp.LE.wmaxz) then
!max            if(l.eq.3.AND.idrank.eq.1) then
!max                  write(0,*)'14, wminy+1 ',wminy+1,'wmaxy-1 ',wmaxy-1
!max            endif
            forall(i=wminx+1:wmaxx-1,j=wminy+1:wmaxy-1) buffservice3d(i+ishift,j+jshift,ktemp) = aoptp(l)%p(i,j,k)
!!         endif
         do i=wminx,wmaxx,wmaxx-wminx
            itemp=i+ishift
            if(itemp.le.1.or.itemp.ge.(nx)) then
               if(ixpbc.eq.1) then
                  if(itemp.eq.0) then
                     itemp=nx
                  endif
                  if(itemp.eq.(nx+1)) then
                     itemp=1
                  endif
               endif
            endif
!!            if(itemp.GE.wminx.AND.itemp.LE.wmaxx.AND.ktemp.GE.wminz.AND.ktemp.LE.wmaxz) then
!max               if(l.eq.3.AND.idrank.eq.1) then
!max                  write(0,*)'16, wminy+1 ',wminy+1,'wmaxy-1 ',wmaxy-1
!max               endif
               forall(j=wminy+1:wmaxy-1) buffservice3d(itemp,j+jshift,ktemp) = aoptp(l)%p(i,j,k)
!!            endif
         enddo
         do j=wminy,wmaxy,wmaxy-wminy
            jtemp=j+jshift;
            if(jtemp.le.1.or.jtemp.ge.(ny)) then
               if(iypbc.eq.1) then
                  if(jtemp.eq.0) then
                     jtemp=ny
                  endif
                  if(jtemp.eq.(ny+1)) then
                     jtemp=1
                  endif
               endif
            endif
!!            if(jtemp.GE.wminy.AND.jtemp.LE.wmaxy.AND.ktemp.GE.wminz.AND.ktemp.LE.wmaxz) then
!max               if(l.eq.3.AND.idrank.eq.1) then
!max                   write(0,*)'18, jtemp ',jtemp
!max                  endif
               forall(i=wminx+1:wmaxx-1) buffservice3d(i+ishift,jtemp,ktemp) = aoptp(l)%p(i,j,k)
!!            endif
         enddo
      enddo
      do i=wminx,wmaxx,wmaxx-wminx
         itemp=i+ishift
         if(ixpbc.eq.1) then
            if(itemp.eq.0) then
               itemp=nx
            endif
            if(itemp.eq.(nx+1)) then
               itemp=1
            endif
         endif
         do j=wminy,wmaxy,wmaxy-wminy
            jtemp=j+jshift
            if(iypbc.eq.1) then
               if(jtemp.eq.0) then
                  jtemp=ny
               endif
               if(jtemp.eq.(ny+1)) then
                  jtemp=1
               endif
            endif
            do k=wminz,wmaxz,wmaxz-wminz
               ktemp=k+kshift;
               if(izpbc.eq.1) then
                  if(ktemp.eq.0) then
                     ktemp=nz
                  endif
                  if(ktemp.eq.(nz+1)) then
                     ktemp=1
                  endif
               endif
               buffservice3d(itemp,jtemp,ktemp) = aoptp(l)%p(i,j,k) 
            enddo
         enddo
      enddo
#endif    
!max   forall(i=wminx:wmaxx,j=wminy:wmaxy,k=wminz:wmaxz) aoptp(l)%p(i,j,k) = buffservice3d(i,j,k)
       forall(i=minx-1:maxx+1,j=miny-1:maxy+1,k=minz-1:maxz+1) aoptp(l)%p(i,j,k) = buffservice3d(i,j,k)
   enddo


   call commrpop(aoptp,lparticles,isfluid)

#if 0
  if(iter.eq.1) then
  do l=1,links
!     write(50+idrank,*)'out pop ',l
        do i=wminx,wmaxx
         do j=miny,maxy
            do k=wminz,wmaxz
               write((100*idrank)+270+l,*)i,j,k,aoptp(l)%p(i,j,k)
           enddo
         enddo
        enddo
     enddo
  endif
#endif
   



  
#if 0
  if(iter.eq.1) then
  do l=1,links
!     write(50+idrank,*)'out pop ',l
        do i=wminx,wmaxx
         do j=wminy,wmaxy
            do k=wminz,wmaxz
               write((100*idrank)+270+l,*)i,j,k,aoptp(l)%p(i,j,k)
           enddo
         enddo
        enddo
     enddo
  endif
#endif
   
  return
  
 end subroutine streaming_fluids
 
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
  real(kind=PRC) :: locrhor, locrhob, weight_RB
  real(kind=PRC) :: ddx,ddy,ddz,ddxB,ddyB,ddzB
  real(kind=8), parameter :: RHOEPS = 9.0E-15
  


  !compute density and accumulate mass flux
  
  !red fluid
  
  if(lsingle_fluid)then
    nfluidnodes=0
    nfluidmassR=ZERO
    ! forall(i=minx:maxx,j=miny:maxy,k=minz:maxz)
   do k=minz,maxz
    do j=miny,maxy
     do i=minx,maxx
      locrho = &
 f00R(i,j,k) + f01R(i,j,k) + f02R(i,j,k) + f03R(i,j,k) + f04R(i,j,k) + &
 f05R(i,j,k) + f06R(i,j,k) + f07R(i,j,k) + f08R(i,j,k) + f09R(i,j,k) + &
 f10R(i,j,k) + f11R(i,j,k) + f12R(i,j,k) + f13R(i,j,k) + f14R(i,j,k) + &
 f15R(i,j,k) + f16R(i,j,k) + f17R(i,j,k) + f18R(i,j,k)

      invrho = ONE / locrho

      locu    = invrho * ( &
 f01R(i,j,k) - f02R(i,j,k) + f07R(i,j,k) - f08R(i,j,k) - f09R(i,j,k) + &
 f10R(i,j,k) + f11R(i,j,k) - f12R(i,j,k) - f13R(i,j,k) + f14R(i,j,k) )

      locv    = invrho * ( &
 f03R(i,j,k) - f04R(i,j,k) + f07R(i,j,k) - f08R(i,j,k) + f09R(i,j,k) - &
 f10R(i,j,k) + f15R(i,j,k) - f16R(i,j,k) - f17R(i,j,k) + f18R(i,j,k) )

      locw    = invrho * ( &
 f05R(i,j,k) - f06R(i,j,k) + f11R(i,j,k) - f12R(i,j,k) + f13R(i,j,k) - &
 f14R(i,j,k) + f15R(i,j,k) - f16R(i,j,k) + f17R(i,j,k) - f18R(i,j,k) )

    if (isfluid(i,j,k)==1) then
      if (locrho < RHOEPS) locrho = RHOEPS
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
    !end forall

    
    return
  endif

#define NEW_MOMENTS_2FL
#ifdef NEW_MOMENTS_2FL
  factR = ONE/tauR
  factB = ONE/tauB
  nfluidnodes=0
  nfluidmassR=ZERO
  nfluidmassB=ZERO
   ! forall(i=minx:maxx,j=miny:maxy,k=minz:maxz)
   do k=minz,maxz
    do j=miny,maxy
     do i=minx,maxx
      locrhor = &
 f00R(i,j,k) + f01R(i,j,k) + f02R(i,j,k) + f03R(i,j,k) + f04R(i,j,k) + &
 f05R(i,j,k) + f06R(i,j,k) + f07R(i,j,k) + f08R(i,j,k) + f09R(i,j,k) + &
 f10R(i,j,k) + f11R(i,j,k) + f12R(i,j,k) + f13R(i,j,k) + f14R(i,j,k) + &
 f15R(i,j,k) + f16R(i,j,k) + f17R(i,j,k) + f18R(i,j,k)

      locrhob = &
 f00B(i,j,k) + f01B(i,j,k) + f02B(i,j,k) + f03B(i,j,k) + f04B(i,j,k) + &
 f05B(i,j,k) + f06B(i,j,k) + f07B(i,j,k) + f08B(i,j,k) + f09B(i,j,k) + &
 f10B(i,j,k) + f11B(i,j,k) + f12B(i,j,k) + f13B(i,j,k) + f14B(i,j,k) + &
 f15B(i,j,k) + f16B(i,j,k) + f17B(i,j,k) + f18B(i,j,k)

      locu  = factR * ( &
 f01R(i,j,k) - f02R(i,j,k) + f07R(i,j,k) - f08R(i,j,k) - f09R(i,j,k) + &
 f10R(i,j,k) + f11R(i,j,k) - f12R(i,j,k) - f13R(i,j,k) + f14R(i,j,k) ) + &
      factB * ( &
 f01B(i,j,k) - f02B(i,j,k) + f07B(i,j,k) - f08B(i,j,k) - f09B(i,j,k) + &
 f10B(i,j,k) + f11B(i,j,k) - f12B(i,j,k) - f13B(i,j,k) + f14B(i,j,k) )

      locv    = factR * ( &
 f03R(i,j,k) - f04R(i,j,k) + f07R(i,j,k) - f08R(i,j,k) + f09R(i,j,k) - &
 f10R(i,j,k) + f15R(i,j,k) - f16R(i,j,k) - f17R(i,j,k) + f18R(i,j,k) ) + &
        factB * ( &
 f03B(i,j,k) - f04B(i,j,k) + f07B(i,j,k) - f08B(i,j,k) + f09B(i,j,k) - &
 f10B(i,j,k) + f15B(i,j,k) - f16B(i,j,k) - f17B(i,j,k) + f18B(i,j,k) )

      locw    = factR * ( &
 f05R(i,j,k) - f06R(i,j,k) + f11R(i,j,k) - f12R(i,j,k) + f13R(i,j,k) - &
 f14R(i,j,k) + f15R(i,j,k) - f16R(i,j,k) + f17R(i,j,k) - f18R(i,j,k) ) + &
        factB * ( &
 f05B(i,j,k) - f06B(i,j,k) + f11B(i,j,k) - f12B(i,j,k) + f13B(i,j,k) - &
 f14B(i,j,k) + f15B(i,j,k) - f16B(i,j,k) + f17B(i,j,k) - f18B(i,j,k) )

    if (isfluid(i,j,k)==1) then
      if (locrhor < RHOEPS) locrhor = RHOEPS
      if (locrhob < RHOEPS) locrhob = RHOEPS
      rhoR(i,j,k) = locrhor
      rhoB(i,j,k) = locrhob
      weight_RB = ONE / (locrhor * factR + locrhob * factB)
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
    !end forall
  
  !compute speed from mass flux
  ! forall(i=minx:maxx,j=miny:maxy,k=minz:maxz,isfluid(i,j,k)==1)
  !   u(i,j,k) = u(i,j,k)/(rhoR(i,j,k)/tauR + rhoB(i,j,k)/tauB)
  !   v(i,j,k) = v(i,j,k)/(rhoR(i,j,k)/tauR + rhoB(i,j,k)/tauB)
  !   w(i,j,k) = w(i,j,k)/(rhoR(i,j,k)/tauR + rhoB(i,j,k)/tauB)
  !end forall
#else
  forall(i=minx:maxx,j=miny:maxy,k=minz:maxz)
    rhoR(i,j,k) = ZERO
    rhoB(i,j,k) = ZERO
    u(i,j,k)    = ZERO
    v(i,j,k)    = ZERO
    w(i,j,k)    = ZERO
  end forall

  !red and blue fluid

  l=0
  ddx=dex(l)/tauR
  ddy=dey(l)/tauR
  ddz=dez(l)/tauR
  ddxB=dex(l)/tauB
  ddyB=dey(l)/tauB
  ddzB=dez(l)/tauB
  forall(i=minx:maxx,j=miny:maxy,k=minz:maxz,isfluid(i,j,k)==1)
    rhoR(i,j,k) = f00R(i,j,k)
    rhoB(i,j,k) = f00B(i,j,k)
  end forall

  l=1
  ddx=dex(l)/tauR
  ddy=dey(l)/tauR
  ddz=dez(l)/tauR
  ddxB=dex(l)/tauB
  ddyB=dey(l)/tauB
  ddzB=dez(l)/tauB
  forall(i=minx:maxx,j=miny:maxy,k=minz:maxz,isfluid(i,j,k)==1)
    rhoR(i,j,k) = f01R(i,j,k) + rhoR(i,j,k)
    rhoB(i,j,k) = f01B(i,j,k) + rhoB(i,j,k)
    u(i,j,k)    = f01R(i,j,k)*ddx + f01B(i,j,k)*ddxB + u(i,j,k)
  end forall

  l=2
  ddx=dex(l)/tauR
  ddy=dey(l)/tauR
  ddz=dez(l)/tauR
  ddxB=dex(l)/tauB
  ddyB=dey(l)/tauB
  ddzB=dez(l)/tauB
  forall(i=minx:maxx,j=miny:maxy,k=minz:maxz,isfluid(i,j,k)==1)
    rhoR(i,j,k) = f02R(i,j,k) + rhoR(i,j,k)
    rhoB(i,j,k) = f02B(i,j,k) + rhoB(i,j,k)
    u(i,j,k)    = f02R(i,j,k)*ddx + f02B(i,j,k)*ddxB + u(i,j,k)
  end forall

  l=3
  ddx=dex(l)/tauR
  ddy=dey(l)/tauR
  ddz=dez(l)/tauR
  ddxB=dex(l)/tauB
  ddyB=dey(l)/tauB
  ddzB=dez(l)/tauB
  forall(i=minx:maxx,j=miny:maxy,k=minz:maxz,isfluid(i,j,k)==1)
    rhoR(i,j,k) = f03R(i,j,k) + rhoR(i,j,k)
    rhoB(i,j,k) = f03B(i,j,k) + rhoB(i,j,k)
    v(i,j,k)    = f03R(i,j,k)*ddy + f03B(i,j,k)*ddyB + v(i,j,k)
  end forall

  l=4
  ddx=dex(l)/tauR
  ddy=dey(l)/tauR
  ddz=dez(l)/tauR
  ddxB=dex(l)/tauB
  ddyB=dey(l)/tauB
  ddzB=dez(l)/tauB
  forall(i=minx:maxx,j=miny:maxy,k=minz:maxz,isfluid(i,j,k)==1)
    rhoR(i,j,k) = f04R(i,j,k) + rhoR(i,j,k)
    rhoB(i,j,k) = f04B(i,j,k) + rhoB(i,j,k)
    v(i,j,k)    = f04R(i,j,k)*ddy + f04B(i,j,k)*ddyB + v(i,j,k)
  end forall

  l=5
  ddx=dex(l)/tauR
  ddy=dey(l)/tauR
  ddz=dez(l)/tauR
  ddxB=dex(l)/tauB
  ddyB=dey(l)/tauB
  ddzB=dez(l)/tauB
  forall(i=minx:maxx,j=miny:maxy,k=minz:maxz,isfluid(i,j,k)==1)
    rhoR(i,j,k) = f05R(i,j,k) + rhoR(i,j,k)
    rhoB(i,j,k) = f05B(i,j,k) + rhoB(i,j,k)
    w(i,j,k)    = f05R(i,j,k)*ddz + f05B(i,j,k)*ddzB + w(i,j,k)
  end forall

  l=6
  ddx=dex(l)/tauR
  ddy=dey(l)/tauR
  ddz=dez(l)/tauR
  ddxB=dex(l)/tauB
  ddyB=dey(l)/tauB
  ddzB=dez(l)/tauB
  forall(i=minx:maxx,j=miny:maxy,k=minz:maxz,isfluid(i,j,k)==1)
    rhoR(i,j,k) = f06R(i,j,k) + rhoR(i,j,k)
    rhoB(i,j,k) = f06B(i,j,k) + rhoB(i,j,k)
    w(i,j,k)    = f06R(i,j,k)*ddz + f06B(i,j,k)*ddzB + w(i,j,k)
  end forall

  l=7
  ddx=dex(l)/tauR
  ddy=dey(l)/tauR
  ddz=dez(l)/tauR
  ddxB=dex(l)/tauB
  ddyB=dey(l)/tauB
  ddzB=dez(l)/tauB
  forall(i=minx:maxx,j=miny:maxy,k=minz:maxz,isfluid(i,j,k)==1)
    rhoR(i,j,k) = f07R(i,j,k) + rhoR(i,j,k)
    rhoB(i,j,k) = f07B(i,j,k) + rhoB(i,j,k)
    u(i,j,k)    = f07R(i,j,k)*ddx + f07B(i,j,k)*ddxB + u(i,j,k)
    v(i,j,k)    = f07R(i,j,k)*ddy + f07B(i,j,k)*ddyB + v(i,j,k)
  end forall

  l=8
  ddx=dex(l)/tauR
  ddy=dey(l)/tauR
  ddz=dez(l)/tauR
  ddxB=dex(l)/tauB
  ddyB=dey(l)/tauB
  ddzB=dez(l)/tauB
  forall(i=minx:maxx,j=miny:maxy,k=minz:maxz,isfluid(i,j,k)==1)
    rhoR(i,j,k) = f08R(i,j,k) + rhoR(i,j,k)
    rhoB(i,j,k) = f08B(i,j,k) + rhoB(i,j,k)
    u(i,j,k)    = f08R(i,j,k)*ddx + f08B(i,j,k)*ddxB + u(i,j,k)
    v(i,j,k)    = f08R(i,j,k)*ddy + f08B(i,j,k)*ddyB + v(i,j,k)
  end forall

  l=9
  ddx=dex(l)/tauR
  ddy=dey(l)/tauR
  ddz=dez(l)/tauR
  ddxB=dex(l)/tauB
  ddyB=dey(l)/tauB
  ddzB=dez(l)/tauB
  forall(i=minx:maxx,j=miny:maxy,k=minz:maxz,isfluid(i,j,k)==1)
    rhoR(i,j,k) = f09R(i,j,k) + rhoR(i,j,k)
    rhoB(i,j,k) = f09B(i,j,k) + rhoB(i,j,k)
    u(i,j,k)    = f09R(i,j,k)*ddx + f09B(i,j,k)*ddxB + u(i,j,k)
    v(i,j,k)    = f09R(i,j,k)*ddy + f09B(i,j,k)*ddyB + v(i,j,k)
  end forall

  l=10
  ddx=dex(l)/tauR
  ddy=dey(l)/tauR
  ddz=dez(l)/tauR
  ddxB=dex(l)/tauB
  ddyB=dey(l)/tauB
  ddzB=dez(l)/tauB
  forall(i=minx:maxx,j=miny:maxy,k=minz:maxz,isfluid(i,j,k)==1)
    rhoR(i,j,k) = f10R(i,j,k) + rhoR(i,j,k)
    rhoB(i,j,k) = f10B(i,j,k) + rhoB(i,j,k)
    u(i,j,k)    = f10R(i,j,k)*ddx + f10B(i,j,k)*ddxB + u(i,j,k)
    v(i,j,k)    = f10R(i,j,k)*ddy + f10B(i,j,k)*ddyB + v(i,j,k)
  end forall

  l=11
  ddx=dex(l)/tauR
  ddy=dey(l)/tauR
  ddz=dez(l)/tauR
  ddxB=dex(l)/tauB
  ddyB=dey(l)/tauB
  ddzB=dez(l)/tauB
  forall(i=minx:maxx,j=miny:maxy,k=minz:maxz,isfluid(i,j,k)==1)
    rhoR(i,j,k) = f11R(i,j,k) + rhoR(i,j,k)
    rhoB(i,j,k) = f11B(i,j,k) + rhoB(i,j,k)
    u(i,j,k)    = f11R(i,j,k)*ddx + f11B(i,j,k)*ddxB + u(i,j,k)
    w(i,j,k)    = f11R(i,j,k)*ddz + f11B(i,j,k)*ddzB + w(i,j,k)
  end forall

  l=12
  ddx=dex(l)/tauR
  ddy=dey(l)/tauR
  ddz=dez(l)/tauR
  ddxB=dex(l)/tauB
  ddyB=dey(l)/tauB
  ddzB=dez(l)/tauB
  forall(i=minx:maxx,j=miny:maxy,k=minz:maxz,isfluid(i,j,k)==1)
    rhoR(i,j,k) = f12R(i,j,k) + rhoR(i,j,k)
    rhoB(i,j,k) = f12B(i,j,k) + rhoB(i,j,k)
    u(i,j,k)    = f12R(i,j,k)*ddx + f12B(i,j,k)*ddxB + u(i,j,k)
    w(i,j,k)    = f12R(i,j,k)*ddz + f12B(i,j,k)*ddzB + w(i,j,k)
  end forall

  l=13
  ddx=dex(l)/tauR
  ddy=dey(l)/tauR
  ddz=dez(l)/tauR
  ddxB=dex(l)/tauB
  ddyB=dey(l)/tauB
  ddzB=dez(l)/tauB
  forall(i=minx:maxx,j=miny:maxy,k=minz:maxz,isfluid(i,j,k)==1)
    rhoR(i,j,k) = f13R(i,j,k) + rhoR(i,j,k)
    rhoB(i,j,k) = f13B(i,j,k) + rhoB(i,j,k)
    u(i,j,k)    = f13R(i,j,k)*ddx + f13B(i,j,k)*ddxB + u(i,j,k)
    w(i,j,k)    = f13R(i,j,k)*ddz + f13B(i,j,k)*ddzB + w(i,j,k)
  end forall

  l=14
  ddx=dex(l)/tauR
  ddy=dey(l)/tauR
  ddz=dez(l)/tauR
  ddxB=dex(l)/tauB
  ddyB=dey(l)/tauB
  ddzB=dez(l)/tauB
  forall(i=minx:maxx,j=miny:maxy,k=minz:maxz,isfluid(i,j,k)==1)
    rhoR(i,j,k) = f14R(i,j,k) + rhoR(i,j,k)
    rhoB(i,j,k) = f14B(i,j,k) + rhoB(i,j,k)
    u(i,j,k)    = f14R(i,j,k)*ddx + f14B(i,j,k)*ddxB + u(i,j,k)
    w(i,j,k)    = f14R(i,j,k)*ddz + f14B(i,j,k)*ddzB + w(i,j,k)
  end forall

  l=15
  ddx=dex(l)/tauR
  ddy=dey(l)/tauR
  ddz=dez(l)/tauR
  ddxB=dex(l)/tauB
  ddyB=dey(l)/tauB
  ddzB=dez(l)/tauB
  forall(i=minx:maxx,j=miny:maxy,k=minz:maxz,isfluid(i,j,k)==1)
    rhoR(i,j,k) = f15R(i,j,k) + rhoR(i,j,k)
    rhoB(i,j,k) = f15B(i,j,k) + rhoB(i,j,k)
    v(i,j,k)    = f15R(i,j,k)*ddy + f15B(i,j,k)*ddyB + v(i,j,k)
    w(i,j,k)    = f15R(i,j,k)*ddz + f15B(i,j,k)*ddzB + w(i,j,k)
  end forall

  l=16
  ddx=dex(l)/tauR
  ddy=dey(l)/tauR
  ddz=dez(l)/tauR
  ddxB=dex(l)/tauB
  ddyB=dey(l)/tauB
  ddzB=dez(l)/tauB
  forall(i=minx:maxx,j=miny:maxy,k=minz:maxz,isfluid(i,j,k)==1)
    rhoR(i,j,k) = f16R(i,j,k) + rhoR(i,j,k)
    rhoB(i,j,k) = f16B(i,j,k) + rhoB(i,j,k)
    v(i,j,k)    = f16R(i,j,k)*ddy + f16B(i,j,k)*ddyB + v(i,j,k)
    w(i,j,k)    = f16R(i,j,k)*ddz + f16B(i,j,k)*ddzB + w(i,j,k)
  end forall

  l=17
  ddx=dex(l)/tauR
  ddy=dey(l)/tauR
  ddz=dez(l)/tauR
  ddxB=dex(l)/tauB
  ddyB=dey(l)/tauB
  ddzB=dez(l)/tauB
  forall(i=minx:maxx,j=miny:maxy,k=minz:maxz,isfluid(i,j,k)==1)
    rhoR(i,j,k) = f17R(i,j,k) + rhoR(i,j,k)
    rhoB(i,j,k) = f17B(i,j,k) + rhoB(i,j,k)
    v(i,j,k)    = f17R(i,j,k)*ddy + f17B(i,j,k)*ddyB + v(i,j,k)
    w(i,j,k)    = f17R(i,j,k)*ddz + f17B(i,j,k)*ddzB + w(i,j,k)
  end forall

  l=18
  ddx=dex(l)/tauR
  ddy=dey(l)/tauR
  ddz=dez(l)/tauR
  ddxB=dex(l)/tauB
  ddyB=dey(l)/tauB
  ddzB=dez(l)/tauB
  forall(i=minx:maxx,j=miny:maxy,k=minz:maxz,isfluid(i,j,k)==1)
    rhoR(i,j,k) = f18R(i,j,k) + rhoR(i,j,k)
    rhoB(i,j,k) = f18B(i,j,k) + rhoB(i,j,k)
    v(i,j,k)    = f18R(i,j,k)*ddy + f18B(i,j,k)*ddyB + v(i,j,k)
    w(i,j,k)    = f18R(i,j,k)*ddz + f18B(i,j,k)*ddzB + w(i,j,k)
  end forall

  !compute speed from mass flux
  forall(i=minx:maxx,j=miny:maxy,k=minz:maxz,isfluid(i,j,k)==1)
    u(i,j,k) = u(i,j,k)/(rhoR(i,j,k)/tauR + rhoB(i,j,k)/tauB)
    v(i,j,k) = v(i,j,k)/(rhoR(i,j,k)/tauR + rhoB(i,j,k)/tauB)
    w(i,j,k) = w(i,j,k)/(rhoR(i,j,k)/tauR + rhoB(i,j,k)/tauB)
  end forall
#endif
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
  
  l=0
  dtemp1 = f00R(i,j,k) + dtemp1 
  dtemp2    = f00R(i,j,k)*dex(l) + dtemp2
  dtemp3    = f00R(i,j,k)*dey(l) + dtemp3
  dtemp4    = f00R(i,j,k)*dez(l) + dtemp4
  
  l=1
  dtemp1 = f01R(i,j,k) + dtemp1
  dtemp2    = f01R(i,j,k)*dex(l) + dtemp2
  dtemp3    = f01R(i,j,k)*dey(l) + dtemp3
  dtemp4    = f01R(i,j,k)*dez(l) + dtemp4
  
  l=2
  dtemp1 = f02R(i,j,k) + dtemp1
  dtemp2    = f02R(i,j,k)*dex(l) + dtemp2
  dtemp3    = f02R(i,j,k)*dey(l) + dtemp3
  dtemp4    = f02R(i,j,k)*dez(l) + dtemp4
  
  l=3
  dtemp1 = f03R(i,j,k) + dtemp1
  dtemp2    = f03R(i,j,k)*dex(l) + dtemp2
  dtemp3    = f03R(i,j,k)*dey(l) + dtemp3
  dtemp4    = f03R(i,j,k)*dez(l) + dtemp4
  
  l=4
  dtemp1 = f04R(i,j,k) + dtemp1
  dtemp2    = f04R(i,j,k)*dex(l) + dtemp2
  dtemp3    = f04R(i,j,k)*dey(l) + dtemp3
  dtemp4    = f04R(i,j,k)*dez(l) + dtemp4
  
  l=5
  dtemp1 = f05R(i,j,k) + dtemp1
  dtemp2    = f05R(i,j,k)*dex(l) + dtemp2
  dtemp3    = f05R(i,j,k)*dey(l) + dtemp3
  dtemp4    = f05R(i,j,k)*dez(l) + dtemp4
  
  l=6
  dtemp1 = f06R(i,j,k) + dtemp1
  dtemp2    = f06R(i,j,k)*dex(l) + dtemp2
  dtemp3    = f06R(i,j,k)*dey(l) + dtemp3
  dtemp4    = f06R(i,j,k)*dez(l) + dtemp4
  
  l=7
  dtemp1 = f07R(i,j,k) + dtemp1
  dtemp2    = f07R(i,j,k)*dex(l) + dtemp2
  dtemp3    = f07R(i,j,k)*dey(l) + dtemp3
  dtemp4    = f07R(i,j,k)*dez(l) + dtemp4
  
  l=8
  dtemp1 = f08R(i,j,k) + dtemp1
  dtemp2    = f08R(i,j,k)*dex(l) + dtemp2
  dtemp3    = f08R(i,j,k)*dey(l) + dtemp3
  dtemp4    = f08R(i,j,k)*dez(l) + dtemp4
  
  l=9
  dtemp1 = f09R(i,j,k) + dtemp1
  dtemp2    = f09R(i,j,k)*dex(l) + dtemp2
  dtemp3    = f09R(i,j,k)*dey(l) + dtemp3
  dtemp4    = f09R(i,j,k)*dez(l) + dtemp4
  
  l=10
  dtemp1 = f10R(i,j,k) + dtemp1
  dtemp2    = f10R(i,j,k)*dex(l) + dtemp2
  dtemp3    = f10R(i,j,k)*dey(l) + dtemp3
  dtemp4    = f10R(i,j,k)*dez(l) + dtemp4
  
  l=11
  dtemp1 = f11R(i,j,k) + dtemp1
  dtemp2    = f11R(i,j,k)*dex(l) + dtemp2
  dtemp3    = f11R(i,j,k)*dey(l) + dtemp3
  dtemp4    = f11R(i,j,k)*dez(l) + dtemp4
  
  l=12
  dtemp1 = f12R(i,j,k) + dtemp1
  dtemp2    = f12R(i,j,k)*dex(l) + dtemp2
  dtemp3    = f12R(i,j,k)*dey(l) + dtemp3
  dtemp4    = f12R(i,j,k)*dez(l) + dtemp4
  
  l=13
  dtemp1 = f13R(i,j,k) + dtemp1
  dtemp2    = f13R(i,j,k)*dex(l) + dtemp2
  dtemp3    = f13R(i,j,k)*dey(l) + dtemp3
  dtemp4    = f13R(i,j,k)*dez(l) + dtemp4
  
  l=14
  dtemp1 = f14R(i,j,k) + dtemp1
  dtemp2    = f14R(i,j,k)*dex(l) + dtemp2
  dtemp3    = f14R(i,j,k)*dey(l) + dtemp3
  dtemp4    = f14R(i,j,k)*dez(l) + dtemp4
  
  l=15
  dtemp1 = f15R(i,j,k) + dtemp1
  dtemp2    = f15R(i,j,k)*dex(l) + dtemp2
  dtemp3    = f15R(i,j,k)*dey(l) + dtemp3
  dtemp4    = f15R(i,j,k)*dez(l) + dtemp4
  
  l=16
  dtemp1 = f16R(i,j,k) + dtemp1
  dtemp2    = f16R(i,j,k)*dex(l) + dtemp2
  dtemp3    = f16R(i,j,k)*dey(l) + dtemp3
  dtemp4    = f16R(i,j,k)*dez(l) + dtemp4
  
  l=17
  dtemp1 = f17R(i,j,k) + dtemp1
  dtemp2    = f17R(i,j,k)*dex(l) + dtemp2
  dtemp3    = f17R(i,j,k)*dey(l) + dtemp3
  dtemp4    = f17R(i,j,k)*dez(l) + dtemp4
  
  l=18
  dtemp1 = f18R(i,j,k) + dtemp1
  dtemp2    = f18R(i,j,k)*dex(l) + dtemp2
  dtemp3    = f18R(i,j,k)*dey(l) + dtemp3
  dtemp4    = f18R(i,j,k)*dez(l) + dtemp4
  
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
  
  l=0
  dtemp1 = f00B(i,j,k) + dtemp1 
  dtemp2    = f00B(i,j,k)*dex(l) + dtemp2
  dtemp3    = f00B(i,j,k)*dey(l) + dtemp3
  dtemp4    = f00B(i,j,k)*dez(l) + dtemp4
  
  l=1
  dtemp1 = f01B(i,j,k) + dtemp1
  dtemp2    = f01B(i,j,k)*dex(l) + dtemp2
  dtemp3    = f01B(i,j,k)*dey(l) + dtemp3
  dtemp4    = f01B(i,j,k)*dez(l) + dtemp4
  
  l=2
  dtemp1 = f02B(i,j,k) + dtemp1
  dtemp2    = f02B(i,j,k)*dex(l) + dtemp2
  dtemp3    = f02B(i,j,k)*dey(l) + dtemp3
  dtemp4    = f02B(i,j,k)*dez(l) + dtemp4
  
  l=3
  dtemp1 = f03B(i,j,k) + dtemp1
  dtemp2    = f03B(i,j,k)*dex(l) + dtemp2
  dtemp3    = f03B(i,j,k)*dey(l) + dtemp3
  dtemp4    = f03B(i,j,k)*dez(l) + dtemp4
  
  l=4
  dtemp1 = f04B(i,j,k) + dtemp1
  dtemp2    = f04B(i,j,k)*dex(l) + dtemp2
  dtemp3    = f04B(i,j,k)*dey(l) + dtemp3
  dtemp4    = f04B(i,j,k)*dez(l) + dtemp4
  
  l=5
  dtemp1 = f05B(i,j,k) + dtemp1
  dtemp2    = f05B(i,j,k)*dex(l) + dtemp2
  dtemp3    = f05B(i,j,k)*dey(l) + dtemp3
  dtemp4    = f05B(i,j,k)*dez(l) + dtemp4
  
  l=6
  dtemp1 = f06B(i,j,k) + dtemp1
  dtemp2    = f06B(i,j,k)*dex(l) + dtemp2
  dtemp3    = f06B(i,j,k)*dey(l) + dtemp3
  dtemp4    = f06B(i,j,k)*dez(l) + dtemp4
  
  l=7
  dtemp1 = f07B(i,j,k) + dtemp1
  dtemp2    = f07B(i,j,k)*dex(l) + dtemp2
  dtemp3    = f07B(i,j,k)*dey(l) + dtemp3
  dtemp4    = f07B(i,j,k)*dez(l) + dtemp4
  
  l=8
  dtemp1 = f08B(i,j,k) + dtemp1
  dtemp2    = f08B(i,j,k)*dex(l) + dtemp2
  dtemp3    = f08B(i,j,k)*dey(l) + dtemp3
  dtemp4    = f08B(i,j,k)*dez(l) + dtemp4
  
  l=9
  dtemp1 = f09B(i,j,k) + dtemp1
  dtemp2    = f09B(i,j,k)*dex(l) + dtemp2
  dtemp3    = f09B(i,j,k)*dey(l) + dtemp3
  dtemp4    = f09B(i,j,k)*dez(l) + dtemp4
  
  l=10
  dtemp1 = f10B(i,j,k) + dtemp1
  dtemp2    = f10B(i,j,k)*dex(l) + dtemp2
  dtemp3    = f10B(i,j,k)*dey(l) + dtemp3
  dtemp4    = f10B(i,j,k)*dez(l) + dtemp4
  
  l=11
  dtemp1 = f11B(i,j,k) + dtemp1
  dtemp2    = f11B(i,j,k)*dex(l) + dtemp2
  dtemp3    = f11B(i,j,k)*dey(l) + dtemp3
  dtemp4    = f11B(i,j,k)*dez(l) + dtemp4
  
  l=12
  dtemp1 = f12B(i,j,k) + dtemp1
  dtemp2    = f12B(i,j,k)*dex(l) + dtemp2
  dtemp3    = f12B(i,j,k)*dey(l) + dtemp3
  dtemp4    = f12B(i,j,k)*dez(l) + dtemp4
  
  l=13
  dtemp1 = f13B(i,j,k) + dtemp1
  dtemp2    = f13B(i,j,k)*dex(l) + dtemp2
  dtemp3    = f13B(i,j,k)*dey(l) + dtemp3
  dtemp4    = f13B(i,j,k)*dez(l) + dtemp4
  
  l=14
  dtemp1 = f14B(i,j,k) + dtemp1
  dtemp2    = f14B(i,j,k)*dex(l) + dtemp2
  dtemp3    = f14B(i,j,k)*dey(l) + dtemp3
  dtemp4    = f14B(i,j,k)*dez(l) + dtemp4
  
  l=15
  dtemp1 = f15B(i,j,k) + dtemp1
  dtemp2    = f15B(i,j,k)*dex(l) + dtemp2
  dtemp3    = f15B(i,j,k)*dey(l) + dtemp3
  dtemp4    = f15B(i,j,k)*dez(l) + dtemp4
  
  l=16
  dtemp1 = f16B(i,j,k) + dtemp1
  dtemp2    = f16B(i,j,k)*dex(l) + dtemp2
  dtemp3    = f16B(i,j,k)*dey(l) + dtemp3
  dtemp4    = f16B(i,j,k)*dez(l) + dtemp4
  
  l=17
  dtemp1 = f17B(i,j,k) + dtemp1
  dtemp2    = f17B(i,j,k)*dex(l) + dtemp2
  dtemp3    = f17B(i,j,k)*dey(l) + dtemp3
  dtemp4    = f17B(i,j,k)*dez(l) + dtemp4
  
  l=18
  dtemp1 = f18B(i,j,k) + dtemp1
  dtemp2    = f18B(i,j,k)*dex(l) + dtemp2
  dtemp3    = f18B(i,j,k)*dey(l) + dtemp3
  dtemp4    = f18B(i,j,k)*dez(l) + dtemp4
  
  return
  
 end subroutine probe_blue_moments_in_node
 
 subroutine probe_pops_in_node(itemp,jtemp,ktemp,nstepsub,aoptp,rhosub)
  
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
  type(REALPTR), dimension(0:links):: aoptp
  real(kind=PRC), allocatable :: rhosub(:,:,:)
  integer(kind=IPRC) :: i4
  
  integer :: l
  
  i4=i4back(itemp,jtemp,ktemp)
  if(ownern(i4)==idrank)then
    write(6,*)'------------------------------------------------','id =',idrank, &
     rhosub(itemp,jtemp,ktemp)
    do l=0,18
    write(6,*)'l =',l,nstepsub,aoptp(l)%p(itemp,jtemp,ktemp)
    call flush(6)
    enddo
  endif
  return
  
 end subroutine probe_pops_in_node
 
!*************START PART TO DEFINE THE EQUILIBRIUM FUNCTS***************
 
 pure function equil_pop00(myrho,myu,myv,myw)
 
!***********************************************************************
!     
!     LBsoft function for Boltzmann equilibrium distribution 
!     for the population zero
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: myrho,myu,myv,myw
  
  real(kind=PRC) :: equil_pop00
  
  real(kind=PRC) :: uv
  
  integer, parameter :: myl = 0
  
  real(kind=PRC), parameter :: mycssq = cssq
  real(kind=PRC), parameter :: myp = p(myl)
  
  real(kind=PRC), parameter :: mydex = dex(myl)
  real(kind=PRC), parameter :: mydey = dey(myl)
  real(kind=PRC), parameter :: mydez = dez(myl)
  
  equil_pop00=myrho*myp*(ONE-(HALF/mycssq)*(myu**TWO+myv**TWO+myw**TWO))
 
  return
  
 end function equil_pop00
 
 pure function equil_pop01(myrho,myu,myv,myw)
 
!***********************************************************************
!     
!     LBsoft function for Boltzmann equilibrium distribution 
!     for the population one
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: myrho,myu,myv,myw
  
  real(kind=PRC) :: equil_pop01
  
  real(kind=PRC) :: uv
  
  integer, parameter :: myl = 1
  
  real(kind=PRC), parameter :: mycssq = cssq
  real(kind=PRC), parameter :: myp = p(myl)
  
  real(kind=PRC), parameter :: mydex = dex(myl)
  real(kind=PRC), parameter :: mydey = dey(myl)
  real(kind=PRC), parameter :: mydez = dez(myl)
  
#if LATTICE==319
  uv=(ONE/mycssq)*(myu*mydex)
  equil_pop01=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#else
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop01=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#endif
  
  return
  
 end function equil_pop01
 
 pure function equil_pop02(myrho,myu,myv,myw)
 
!***********************************************************************
!     
!     LBsoft function for Boltzmann equilibrium distribution 
!     for the population two
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: myrho,myu,myv,myw
  
  real(kind=PRC) :: equil_pop02
  
  real(kind=PRC) :: uv
  
  integer, parameter :: myl = 2
  
  real(kind=PRC), parameter :: mycssq = cssq
  real(kind=PRC), parameter :: myp = p(myl)
  
  real(kind=PRC), parameter :: mydex = dex(myl)
  real(kind=PRC), parameter :: mydey = dey(myl)
  real(kind=PRC), parameter :: mydez = dez(myl)
  
#if LATTICE==319
  uv=(ONE/mycssq)*(myu*mydex)
  equil_pop02=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#else
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop02=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#endif
  
  return
  
 end function equil_pop02
 
 pure function equil_pop03(myrho,myu,myv,myw)
 
!***********************************************************************
!     
!     LBsoft function for Boltzmann equilibrium distribution 
!     for the population three
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: myrho,myu,myv,myw
  
  real(kind=PRC) :: equil_pop03
  
  real(kind=PRC) :: uv
  
  integer, parameter :: myl = 3
  
  real(kind=PRC), parameter :: mycssq = cssq
  real(kind=PRC), parameter :: myp = p(myl)
  
  real(kind=PRC), parameter :: mydex = dex(myl)
  real(kind=PRC), parameter :: mydey = dey(myl)
  real(kind=PRC), parameter :: mydez = dez(myl)
  
#if LATTICE==319
  uv=(ONE/mycssq)*(myv*mydey)
  equil_pop03=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#else
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop03=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#endif
  
  return
  
 end function equil_pop03
 
 pure function equil_pop04(myrho,myu,myv,myw)
 
!***********************************************************************
!     
!     LBsoft function for Boltzmann equilibrium distribution 
!     for the population four
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: myrho,myu,myv,myw
  
  real(kind=PRC) :: equil_pop04
  
  real(kind=PRC) :: uv
  
  integer, parameter :: myl = 4
  
  real(kind=PRC), parameter :: mycssq = cssq
  real(kind=PRC), parameter :: myp = p(myl)
  
  real(kind=PRC), parameter :: mydex = dex(myl)
  real(kind=PRC), parameter :: mydey = dey(myl)
  real(kind=PRC), parameter :: mydez = dez(myl)
  
#if LATTICE==319
  uv=(ONE/mycssq)*(myv*mydey)
  equil_pop04=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#else
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop04=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#endif
  
  return
  
 end function equil_pop04
 
 pure function equil_pop05(myrho,myu,myv,myw)
 
!***********************************************************************
!     
!     LBsoft function for Boltzmann equilibrium distribution 
!     for the population five
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: myrho,myu,myv,myw
  
  real(kind=PRC) :: equil_pop05
  
  real(kind=PRC) :: uv
  
  integer, parameter :: myl = 5
  
  real(kind=PRC), parameter :: mycssq = cssq
  real(kind=PRC), parameter :: myp = p(myl)
  
  real(kind=PRC), parameter :: mydex = dex(myl)
  real(kind=PRC), parameter :: mydey = dey(myl)
  real(kind=PRC), parameter :: mydez = dez(myl)
  
#if LATTICE==319
  uv=(ONE/mycssq)*(myw*mydez)
  equil_pop05=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#else
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop05=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#endif
  
  return
  
 end function equil_pop05
 
 pure function equil_pop06(myrho,myu,myv,myw)
 
!***********************************************************************
!     
!     LBsoft function for Boltzmann equilibrium distribution 
!     for the population six
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: myrho,myu,myv,myw
  
  real(kind=PRC) :: equil_pop06
  
  real(kind=PRC) :: uv
  
  integer, parameter :: myl = 6
  
  real(kind=PRC), parameter :: mycssq = cssq
  real(kind=PRC), parameter :: myp = p(myl)
  
  real(kind=PRC), parameter :: mydex = dex(myl)
  real(kind=PRC), parameter :: mydey = dey(myl)
  real(kind=PRC), parameter :: mydez = dez(myl)
  
#if LATTICE==319
  uv=(ONE/mycssq)*(myw*mydez)
  equil_pop06=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#else
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop06=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#endif
  
  return
  
 end function equil_pop06
 
 pure function equil_pop07(myrho,myu,myv,myw)
 
!***********************************************************************
!     
!     LBsoft function for Boltzmann equilibrium distribution 
!     for the population seven
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: myrho,myu,myv,myw
  
  real(kind=PRC) :: equil_pop07
  
  real(kind=PRC) :: uv
  
  integer, parameter :: myl = 7
  
  real(kind=PRC), parameter :: mycssq = cssq
  real(kind=PRC), parameter :: myp = p(myl)
  
  real(kind=PRC), parameter :: mydex = dex(myl)
  real(kind=PRC), parameter :: mydey = dey(myl)
  real(kind=PRC), parameter :: mydez = dez(myl)
  
#if LATTICE==319
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey)
  equil_pop07=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#else
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop07=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#endif
  
  return
  
 end function equil_pop07
 
 pure function equil_pop08(myrho,myu,myv,myw)
 
!***********************************************************************
!     
!     LBsoft function for Boltzmann equilibrium distribution 
!     for the population eight
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: myrho,myu,myv,myw
  
  real(kind=PRC) :: equil_pop08
  
  real(kind=PRC) :: uv
  
  integer, parameter :: myl = 8
  
  real(kind=PRC), parameter :: mycssq = cssq
  real(kind=PRC), parameter :: myp = p(myl)
  
  real(kind=PRC), parameter :: mydex = dex(myl)
  real(kind=PRC), parameter :: mydey = dey(myl)
  real(kind=PRC), parameter :: mydez = dez(myl)
  
#if LATTICE==319
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey)
  equil_pop08=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#else
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop08=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#endif
  
  return
  
 end function equil_pop08
 
 pure function equil_pop09(myrho,myu,myv,myw)
 
!***********************************************************************
!     
!     LBsoft function for Boltzmann equilibrium distribution 
!     for the population nine
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: myrho,myu,myv,myw
  
  real(kind=PRC) :: equil_pop09
  
  real(kind=PRC) :: uv
  
  integer, parameter :: myl = 9
  
  real(kind=PRC), parameter :: mycssq = cssq
  real(kind=PRC), parameter :: myp = p(myl)
  
  real(kind=PRC), parameter :: mydex = dex(myl)
  real(kind=PRC), parameter :: mydey = dey(myl)
  real(kind=PRC), parameter :: mydez = dez(myl)
  
#if LATTICE==319
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey)
  equil_pop09=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#else
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop09=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#endif
  
  return
  
 end function equil_pop09
 
 pure function equil_pop10(myrho,myu,myv,myw)
 
!***********************************************************************
!     
!     LBsoft function for Boltzmann equilibrium distribution 
!     for the population ten
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: myrho,myu,myv,myw
  
  real(kind=PRC) :: equil_pop10
  
  real(kind=PRC) :: uv
  
  integer, parameter :: myl = 10
  
  real(kind=PRC), parameter :: mycssq = cssq
  real(kind=PRC), parameter :: myp = p(myl)
  
  real(kind=PRC), parameter :: mydex = dex(myl)
  real(kind=PRC), parameter :: mydey = dey(myl)
  real(kind=PRC), parameter :: mydez = dez(myl)
  
#if LATTICE==319
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey)
  equil_pop10=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#else
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop10=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#endif
  
  return
  
 end function equil_pop10
 
 pure function equil_pop11(myrho,myu,myv,myw)
 
!***********************************************************************
!     
!     LBsoft function for Boltzmann equilibrium distribution 
!     for the population eleven
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: myrho,myu,myv,myw
  
  real(kind=PRC) :: equil_pop11
  
  real(kind=PRC) :: uv
  
  integer, parameter :: myl = 11
  
  real(kind=PRC), parameter :: mycssq = cssq
  real(kind=PRC), parameter :: myp = p(myl)
  
  real(kind=PRC), parameter :: mydex = dex(myl)
  real(kind=PRC), parameter :: mydey = dey(myl)
  real(kind=PRC), parameter :: mydez = dez(myl)
  
#if LATTICE==319
  uv=(ONE/mycssq)*(myu*mydex + myw*mydez)
  equil_pop11=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#else
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop11=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#endif
  
  return
  
 end function equil_pop11
 
 pure function equil_pop12(myrho,myu,myv,myw)
 
!***********************************************************************
!     
!     LBsoft function for Boltzmann equilibrium distribution 
!     for the population twelve
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: myrho,myu,myv,myw
  
  real(kind=PRC) :: equil_pop12
  
  real(kind=PRC) :: uv
  
  integer, parameter :: myl = 12
  
  real(kind=PRC), parameter :: mycssq = cssq
  real(kind=PRC), parameter :: myp = p(myl)
  
  real(kind=PRC), parameter :: mydex = dex(myl)
  real(kind=PRC), parameter :: mydey = dey(myl)
  real(kind=PRC), parameter :: mydez = dez(myl)
  
#if LATTICE==319
  uv=(ONE/mycssq)*(myu*mydex  + myw*mydez)
  equil_pop12=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#else
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop12=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#endif
  
  return
  
 end function equil_pop12
 
 pure function equil_pop13(myrho,myu,myv,myw)
 
!***********************************************************************
!     
!     LBsoft function for Boltzmann equilibrium distribution 
!     for the population thirteen
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: myrho,myu,myv,myw
  
  real(kind=PRC) :: equil_pop13
  
  real(kind=PRC) :: uv
  
  integer, parameter :: myl = 13
  
  real(kind=PRC), parameter :: mycssq = cssq
  real(kind=PRC), parameter :: myp = p(myl)
  
  real(kind=PRC), parameter :: mydex = dex(myl)
  real(kind=PRC), parameter :: mydey = dey(myl)
  real(kind=PRC), parameter :: mydez = dez(myl)
  
#if LATTICE==319
  uv=(ONE/mycssq)*(myu*mydex + myw*mydez)
  equil_pop13=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#else
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop13=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#endif
  
  return
  
 end function equil_pop13
 
 pure function equil_pop14(myrho,myu,myv,myw)
 
!***********************************************************************
!     
!     LBsoft function for Boltzmann equilibrium distribution 
!     for the population fourteen
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: myrho,myu,myv,myw
  
  real(kind=PRC) :: equil_pop14
  
  real(kind=PRC) :: uv
  
  integer, parameter :: myl = 14
  
  real(kind=PRC), parameter :: mycssq = cssq
  real(kind=PRC), parameter :: myp = p(myl)
  
  real(kind=PRC), parameter :: mydex = dex(myl)
  real(kind=PRC), parameter :: mydey = dey(myl)
  real(kind=PRC), parameter :: mydez = dez(myl)
  
#if LATTICE==319
  uv=(ONE/mycssq)*(myu*mydex + myw*mydez)
  equil_pop14=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#else
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop14=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#endif
  
  return
  
 end function equil_pop14
 
 pure function equil_pop15(myrho,myu,myv,myw)
 
!***********************************************************************
!     
!     LBsoft function for Boltzmann equilibrium distribution 
!     for the population fiveteen
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: myrho,myu,myv,myw
  
  real(kind=PRC) :: equil_pop15
  
  real(kind=PRC) :: uv
  
  integer, parameter :: myl = 15
  
  real(kind=PRC), parameter :: mycssq = cssq
  real(kind=PRC), parameter :: myp = p(myl)
  
  real(kind=PRC), parameter :: mydex = dex(myl)
  real(kind=PRC), parameter :: mydey = dey(myl)
  real(kind=PRC), parameter :: mydez = dez(myl)
  
#if LATTICE==319
  uv=(ONE/mycssq)*(myv*mydey + myw*mydez)
  equil_pop15=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#else
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop15=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#endif
  
  return
  
 end function equil_pop15
 
 pure function equil_pop16(myrho,myu,myv,myw)
 
!***********************************************************************
!     
!     LBsoft function for Boltzmann equilibrium distribution 
!     for the population sixteen
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: myrho,myu,myv,myw
  
  real(kind=PRC) :: equil_pop16
  
  real(kind=PRC) :: uv
  
  integer, parameter :: myl = 16
  
  real(kind=PRC), parameter :: mycssq = cssq
  real(kind=PRC), parameter :: myp = p(myl)
  
  real(kind=PRC), parameter :: mydex = dex(myl)
  real(kind=PRC), parameter :: mydey = dey(myl)
  real(kind=PRC), parameter :: mydez = dez(myl)
  
#if LATTICE==319
  uv=(ONE/mycssq)*(myv*mydey + myw*mydez)
  equil_pop16=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#else
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop16=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#endif
  
  return
  
 end function equil_pop16
 
 pure function equil_pop17(myrho,myu,myv,myw)
 
!***********************************************************************
!     
!     LBsoft function for Boltzmann equilibrium distribution 
!     for the population seventeen
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: myrho,myu,myv,myw
  
  real(kind=PRC) :: equil_pop17
  
  real(kind=PRC) :: uv
  
  integer, parameter :: myl = 17
  
  real(kind=PRC), parameter :: mycssq = cssq
  real(kind=PRC), parameter :: myp = p(myl)
  
  real(kind=PRC), parameter :: mydex = dex(myl)
  real(kind=PRC), parameter :: mydey = dey(myl)
  real(kind=PRC), parameter :: mydez = dez(myl)
  
#if LATTICE==319
  uv=(ONE/mycssq)*(myv*mydey + myw*mydez)
  equil_pop17=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#else
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop17=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#endif
  
  return
  
 end function equil_pop17
 
 pure function equil_pop18(myrho,myu,myv,myw)
 
!***********************************************************************
!     
!     LBsoft function for Boltzmann equilibrium distribution 
!     for the population eighteen
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: myrho,myu,myv,myw
  
  real(kind=PRC) :: equil_pop18
  
  real(kind=PRC) :: uv
  
  integer, parameter :: myl = 18
  
  real(kind=PRC), parameter :: mycssq = cssq
  real(kind=PRC), parameter :: myp = p(myl)
  
  real(kind=PRC), parameter :: mydex = dex(myl)
  real(kind=PRC), parameter :: mydey = dey(myl)
  real(kind=PRC), parameter :: mydez = dez(myl)
  
#if LATTICE==319
  uv=(ONE/mycssq)*(myv*mydey + myw*mydez)
  equil_pop18=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#else
  uv=(ONE/mycssq)*(myu*mydex + myv*mydey + myw*mydez)
  equil_pop18=myrho*myp*(ONE+uv+HALF*(uv*uv)-(HALF/mycssq)* &
   (myu**TWO+myv**TWO+myw**TWO))
#endif
  
  return
  
 end function equil_pop18

!**************END PART TO DEFINE THE EQUILIBRIUM FUNCTS****************

!*****************START PART TO MANAGE THE PERIODIC BC******************

 subroutine driver_bc_isfluid
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the boundary conditions
!     to isfluid array
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
 
  implicit none
  Logical, save   :: isFirst = .true.


  

    call commexch_isfluid(isfluid)
  
    call manage_bc_isfluid_selfcomm(isfluid)
    
   
    call commwait_isfluid(isfluid)
   

 end subroutine driver_bc_isfluid
 
 subroutine driver_bc_new_isfluid
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the boundary conditions
!     to new_isfluid array
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
 
  implicit none
  Logical, save   :: isFirst = .true.


    call commexch_isfluid(new_isfluid)
    
    call manage_bc_isfluid_selfcomm(new_isfluid)
    
   
    call commwait_isfluid(new_isfluid)
  

 end subroutine driver_bc_new_isfluid


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
  
    call manage_bc_hvar_selfcomm(rhoR_managebc)
    if(.not. lsingle_fluid)then
      call manage_bc_hvar_selfcomm(rhoB_managebc)
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
  

    call manage_bc_hvar_selfcomm(psiR_managebc)
    if(.not. lsingle_fluid)then
      call manage_bc_hvar_selfcomm(psiB_managebc)
    endif

  
    call commwait_dens(psiR,psiB)
  
  
 end subroutine driver_bc_psi

 
 subroutine initialiaze_manage_bc_isfluid_selfcomm()
 
!***********************************************************************
!     
!     LBsoft subroutine to manage the buffer of isfluid nodes
!     within the same process for applying the boundary conditions
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
 
  implicit none
  
  integer :: i,j,k,l,itemp,jtemp,ktemp,itemp2,jtemp2,ktemp2
  integer(kind=IPRC) :: i4orig,i4
  
  logical, parameter :: lverbose=.false.
  
  nguards=0
  if(ixpbc.eq.0 .and. iypbc.eq.0 .and. izpbc.eq.0)return
  
  do k=minz-nbuff,maxz+nbuff
    do j=miny-nbuff,maxy+nbuff
      do i=minx-nbuff,maxx+nbuff
        if(i<minx.or.j<miny.or.k<minz.or.i>maxx.or.j>maxy.or.k>maxz)then
          if(ixpbc.eq.0 .and. iypbc.eq.0 .and. izpbc.eq.0)cycle
          itemp=i
          jtemp=j
          ktemp=k
          itemp2=i
          jtemp2=j
          ktemp2=k
          !apply periodic conditions if necessary
          if(ixpbc.eq.1) then
            if(itemp<1) then
              itemp=itemp+nx
            endif
            if(itemp>nx) then
              itemp=itemp-nx
            endif
          endif
          if(iypbc.eq.1) then
            if(jtemp<1) then
              jtemp=jtemp+ny
            endif
           if(jtemp>ny) then
              jtemp=jtemp-ny
            endif
          endif
          if(izpbc.eq.1) then
            if(ktemp<1) then
              ktemp=ktemp+nz
            endif
            if(ktemp>nz) then
              ktemp=ktemp-nz
            endif
          endif
          i4=i4back(itemp,jtemp,ktemp) 
          i4orig=i4back(itemp2,jtemp2,ktemp2) 
          if(ownern(i4).EQ.idrank.AND.ownern(i4orig).EQ.idrank) THEN
            nguards=nguards+1
          endif
        endif
      enddo
    enddo
  enddo
  
  if(lverbose)write(6,*)'id=',idrank,'nguards=',nguards
  
  allocate(isguards(6,nguards))
  
  nguards=0
  do k=minz-nbuff,maxz+nbuff
    do j=miny-nbuff,maxy+nbuff
      do i=minx-nbuff,maxx+nbuff
        if(i<minx.or.j<miny.or.k<minz.or.i>maxx.or.j>maxy.or.k>maxz)then
          if(ixpbc.eq.0 .and. iypbc.eq.0 .and. izpbc.eq.0)cycle
          itemp=i
          jtemp=j
          ktemp=k
          itemp2=i
          jtemp2=j
          ktemp2=k
          !apply periodic conditions if necessary
          if(ixpbc.eq.1) then
            if(itemp<1) then
              itemp=itemp+nx
            endif
            if(itemp>nx) then
              itemp=itemp-nx
            endif
          endif
          if(iypbc.eq.1) then
            if(jtemp<1) then
              jtemp=jtemp+ny
            endif
           if(jtemp>ny) then
              jtemp=jtemp-ny
            endif
          endif
          if(izpbc.eq.1) then
            if(ktemp<1) then
              ktemp=ktemp+nz
            endif
            if(ktemp>nz) then
              ktemp=ktemp-nz
            endif
          endif
          i4=i4back(itemp,jtemp,ktemp) 
          i4orig=i4back(itemp2,jtemp2,ktemp2) 
          if(ownern(i4).EQ.idrank.AND.ownern(i4orig).EQ.idrank) THEN
            nguards=nguards+1
            isguards(1,nguards)=itemp
            isguards(2,nguards)=jtemp
            isguards(3,nguards)=ktemp
            isguards(4,nguards)=itemp2
            isguards(5,nguards)=jtemp2
            isguards(6,nguards)=ktemp2
          endif
        endif
      enddo
    enddo
  enddo
   
  return
   
 end subroutine initialiaze_manage_bc_isfluid_selfcomm
 
 subroutine manage_bc_isfluid_selfcomm(istemp)
 
!***********************************************************************
!     
!     LBsoft subroutine to manage the buffer of isfluid nodes
!     within the same process for applying the boundary conditions
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer(kind=1), dimension(:,:,:), allocatable :: istemp
  
  integer :: i,j,k,l,itemp,jtemp,ktemp,itemp2,jtemp2,ktemp2
  integer(kind=IPRC) :: i4orig,i4
  
  if(nguards>=1)then
    forall(i=1:nguards)
      istemp(isguards(4,i),isguards(5,i),isguards(6,i))= &
       istemp(isguards(1,i),isguards(2,i),isguards(3,i))
    end forall
  endif
   
  return
   
 end subroutine manage_bc_isfluid_selfcomm
 
 subroutine initialiaze_manage_bc_singlehalo_selfcomm()
 
!***********************************************************************
!     
!     LBsoft subroutine to manage the buffer of isfluid nodes
!     within the same process for applying the boundary conditions
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
 
  implicit none
  
  integer :: i,j,k,l,itemp,jtemp,ktemp,itemp2,jtemp2,ktemp2
  integer(kind=IPRC) :: i4orig,i4
  
  logical, parameter :: lverbose=.false.
  
  nsinglehalo=0
  if(ixpbc.eq.0 .and. iypbc.eq.0 .and. izpbc.eq.0)return
  
  do k=minz-nbuff,maxz+nbuff
    do j=miny-nbuff,maxy+nbuff
      do i=minx-nbuff,maxx+nbuff
        if(i<minx.or.j<miny.or.k<minz.or.i>maxx.or.j>maxy.or.k>maxz)then
          if(ixpbc.eq.0 .and. iypbc.eq.0 .and. izpbc.eq.0)cycle
          itemp=i
          jtemp=j
          ktemp=k
          itemp2=i
          jtemp2=j
          ktemp2=k
          !apply periodic conditions if necessary
          if(ixpbc.eq.1) then
            if(itemp<1) then
              itemp=itemp+nx
            endif
            if(itemp>nx) then
              itemp=itemp-nx
            endif
          endif
          if(iypbc.eq.1) then
            if(jtemp<1) then
              jtemp=jtemp+ny
            endif
           if(jtemp>ny) then
              jtemp=jtemp-ny
            endif
          endif
          if(izpbc.eq.1) then
            if(ktemp<1) then
              ktemp=ktemp+nz
            endif
            if(ktemp>nz) then
              ktemp=ktemp-nz
            endif
          endif
          i4=i4back(itemp,jtemp,ktemp) 
          i4orig=i4back(itemp2,jtemp2,ktemp2) 
          if(ownern(i4).EQ.idrank.AND.ownern(i4orig).EQ.idrank) THEN
            nsinglehalo=nsinglehalo+1
          endif
        endif
      enddo
    enddo
  enddo
  
  if(lverbose)write(6,*)'id=',idrank,'nsinglehalo=',nsinglehalo
  
  allocate(isinglehalo(6,nsinglehalo))
  
  nsinglehalo=0
  do k=minz-nbuff,maxz+nbuff
    do j=miny-nbuff,maxy+nbuff
      do i=minx-nbuff,maxx+nbuff
        if(i<minx.or.j<miny.or.k<minz.or.i>maxx.or.j>maxy.or.k>maxz)then
          if(ixpbc.eq.0 .and. iypbc.eq.0 .and. izpbc.eq.0)cycle
          itemp=i
          jtemp=j
          ktemp=k
          itemp2=i
          jtemp2=j
          ktemp2=k
          !apply periodic conditions if necessary
          if(ixpbc.eq.1) then
            if(itemp<1) then
              itemp=itemp+nx
            endif
            if(itemp>nx) then
              itemp=itemp-nx
            endif
          endif
          if(iypbc.eq.1) then
            if(jtemp<1) then
              jtemp=jtemp+ny
            endif
           if(jtemp>ny) then
              jtemp=jtemp-ny
            endif
          endif
          if(izpbc.eq.1) then
            if(ktemp<1) then
              ktemp=ktemp+nz
            endif
            if(ktemp>nz) then
              ktemp=ktemp-nz
            endif
          endif
          i4=i4back(itemp,jtemp,ktemp) 
          i4orig=i4back(itemp2,jtemp2,ktemp2) 
          if(ownern(i4).EQ.idrank.AND.ownern(i4orig).EQ.idrank) THEN
            nsinglehalo=nsinglehalo+1
            isinglehalo(1,nsinglehalo)=itemp
            isinglehalo(2,nsinglehalo)=jtemp
            isinglehalo(3,nsinglehalo)=ktemp
            isinglehalo(4,nsinglehalo)=itemp2
            isinglehalo(5,nsinglehalo)=jtemp2
            isinglehalo(6,nsinglehalo)=ktemp2
          endif
        endif
      enddo
    enddo
  enddo
   
  return
   
 end subroutine initialiaze_manage_bc_singlehalo_selfcomm
 
 subroutine manage_bc_singlehalo_selfcomm(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine to manage the buffer of isfluid nodes
!     within the same process for applying the boundary conditions
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), dimension(:,:,:), allocatable :: dtemp
  
  integer :: i,j,k,l,itemp,jtemp,ktemp,itemp2,jtemp2,ktemp2
  integer(kind=IPRC) :: i4orig,i4
  
  if(nsinglehalo>=1)then
    do i=1,nsinglehalo
      dtemp(isinglehalo(4,i),isinglehalo(5,i),isinglehalo(6,i))= &
       dtemp(isinglehalo(1,i),isinglehalo(2,i),isinglehalo(3,i))
    enddo
  endif
   
  return
   
 end subroutine manage_bc_singlehalo_selfcomm
 
 subroutine initialiaze_manage_bc_hvar_selfcomm
 
!***********************************************************************
!     
!     LBsoft subroutine to create a list of nodes which should be
!     managed within the same process for applying the boundary 
!     conditions
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer :: i,j,k,l,itemp,jtemp,ktemp,itemp2,jtemp2,ktemp2
  integer(kind=IPRC) :: i4orig,i4
  
  logical, parameter :: lverbose=.false.
  
  nvar_managebc=0
  if(ixpbc.eq.0 .and. iypbc.eq.0 .and. izpbc.eq.0)return
  
  do k=minz-nbuff,maxz+nbuff
    do j=miny-nbuff,maxy+nbuff
      do i=minx-nbuff,maxx+nbuff
        if(i<minx.or.j<miny.or.k<minz.or.i>maxx.or.j>maxy.or.k>maxz)then
          
          itemp=i
          jtemp=j
          ktemp=k
          itemp2=i
          jtemp2=j
          ktemp2=k
          !apply periodic conditions if necessary
          if(ixpbc.eq.1) then
            if(itemp<1) then
              itemp=itemp+nx
            endif
            if(itemp>nx) then
              itemp=itemp-nx
            endif
          endif
          if(iypbc.eq.1) then
            if(jtemp<1) then
              jtemp=jtemp+ny
            endif
           if(jtemp>ny) then
              jtemp=jtemp-ny
            endif
          endif
          if(izpbc.eq.1) then
            if(ktemp<1) then
              ktemp=ktemp+nz
            endif
            if(ktemp>nz) then
              ktemp=ktemp-nz
            endif
          endif
          i4=i4back(itemp,jtemp,ktemp) 
          i4orig=i4back(itemp2,jtemp2,ktemp2) 
          if(ownern(i4).EQ.idrank.AND.ownern(i4orig).EQ.idrank) THEN
            nvar_managebc=nvar_managebc+1
          endif
        endif
      enddo
    enddo
  enddo
  if(lverbose)write(6,*)'id=',idrank,'nvar_managebc=',nvar_managebc
  
  allocate(rhoR_managebc(2,nvar_managebc))
  allocate(psiR_managebc(2,nvar_managebc))
  if(.not. lsingle_fluid)allocate(rhoB_managebc(2,nvar_managebc))
  if(.not. lsingle_fluid)allocate(psiB_managebc(2,nvar_managebc))
  
  if(lexch_u)then
    allocate(u_managebc(2,nvar_managebc))
  endif
  if(lexch_v)then
    allocate(v_managebc(2,nvar_managebc))
  endif
  if(lexch_w)then
    allocate(w_managebc(2,nvar_managebc))
  endif
  
  nvar_managebc=0
  do k=minz-nbuff,maxz+nbuff
    do j=miny-nbuff,maxy+nbuff
      do i=minx-nbuff,maxx+nbuff
        if(i<minx.or.j<miny.or.k<minz.or.i>maxx.or.j>maxy.or.k>maxz)then
          if(ixpbc.eq.0 .and. iypbc.eq.0 .and. izpbc.eq.0)cycle
          itemp=i
          jtemp=j
          ktemp=k
          itemp2=i
          jtemp2=j
          ktemp2=k
          !apply periodic conditions if necessary
          if(ixpbc.eq.1) then
            if(itemp<1) then
              itemp=itemp+nx
            endif
            if(itemp>nx) then
              itemp=itemp-nx
            endif
          endif
          if(iypbc.eq.1) then
            if(jtemp<1) then
              jtemp=jtemp+ny
            endif
           if(jtemp>ny) then
              jtemp=jtemp-ny
            endif
          endif
          if(izpbc.eq.1) then
            if(ktemp<1) then
              ktemp=ktemp+nz
            endif
            if(ktemp>nz) then
              ktemp=ktemp-nz
            endif
          endif
          i4=i4back(itemp,jtemp,ktemp) 
          i4orig=i4back(itemp2,jtemp2,ktemp2) 
          if(ownern(i4).EQ.idrank.AND.ownern(i4orig).EQ.idrank) THEN
            nvar_managebc=nvar_managebc+1
            rhoR_managebc(1,nvar_managebc)%p=>rhoR(itemp2,jtemp2,ktemp2)   !who must recevice is one
            rhoR_managebc(2,nvar_managebc)%p=>rhoR(itemp,jtemp,ktemp)
            if(lShanChen)then
              psiR_managebc(1,nvar_managebc)%p=>psiR(itemp2,jtemp2,ktemp2)   !who must recevice is one
              psiR_managebc(2,nvar_managebc)%p=>psiR(itemp,jtemp,ktemp)
            endif
            if(.not. lsingle_fluid)then
              rhoB_managebc(1,nvar_managebc)%p=>rhoB(itemp2,jtemp2,ktemp2) !who must recevice is one
              rhoB_managebc(2,nvar_managebc)%p=>rhoB(itemp,jtemp,ktemp)
              if(lShanChen)then
                psiB_managebc(1,nvar_managebc)%p=>psiB(itemp2,jtemp2,ktemp2) !who must recevice is one
                psiB_managebc(2,nvar_managebc)%p=>psiB(itemp,jtemp,ktemp)
              endif
            endif
            if(lexch_u)then
              u_managebc(1,nvar_managebc)%p=>u(itemp2,jtemp2,ktemp2)   !who must recevice is one
              u_managebc(2,nvar_managebc)%p=>u(itemp,jtemp,ktemp)
            endif
            if(lexch_v)then
              v_managebc(1,nvar_managebc)%p=>v(itemp2,jtemp2,ktemp2)   !who must recevice is one
              v_managebc(2,nvar_managebc)%p=>v(itemp,jtemp,ktemp)
            endif
            if(lexch_w)then
              w_managebc(1,nvar_managebc)%p=>w(itemp2,jtemp2,ktemp2)   !who must recevice is one
              w_managebc(2,nvar_managebc)%p=>w(itemp,jtemp,ktemp)
            endif
          endif
        endif
      enddo
    enddo
  enddo
   
  return
   
 end subroutine initialiaze_manage_bc_hvar_selfcomm
 
 subroutine manage_bc_hvar_selfcomm(dtemp)
 
!***********************************************************************
!     
!     LBsoft subroutine to manage the buffer fluid nodes
!     within the same process for applying the boundary conditions
!     using the node list created in subroutine
!     initialiaze_manage_bc_dens_selfcomm
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  type(REALPTR_S), allocatable, dimension(:,:) :: dtemp
  
  integer :: i
  
  forall(i=1:nvar_managebc)
    dtemp(1,i)%p=real(dtemp(2,i)%p,kind=PRC)
  end forall
   
  return
   
 end subroutine manage_bc_hvar_selfcomm
  
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

    call manage_bc_hvar_selfcomm(u_managebc)

    call commwait_vel_component(u)

  
  endif
  
  if(lexch_v)then
  

    call commexch_vel_component(v)

    call manage_bc_hvar_selfcomm(v_managebc)

    call commwait_vel_component(v)
  
  endif
  
  if(lexch_w)then
  

    call commexch_vel_component(w)

    call manage_bc_hvar_selfcomm(w_managebc)

    call commwait_vel_component(w)

  
  endif
  
  isFirst = .false.

  return

  
 end subroutine driver_bc_velocities
 
 subroutine initialiaze_manage_bc_pop_selfcomm()
 
!***********************************************************************
!     
!     LBsoft subroutine to create a list of nodes which should be
!     managed within the same process for applying the boundary 
!     conditions
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer :: i,j,k,l,itemp,jtemp,ktemp,itemp2,jtemp2,ktemp2
  integer :: ishift,jshift,kshift
  integer(kind=IPRC) :: i4orig,i4
  
  logical, parameter :: lverbose=.false.
  
  
  npoplistbc=0
  
  if(ixpbc.eq.0 .and. iypbc.eq.0 .and. izpbc.eq.0)return
  
  do l=1,links
    ishift=ex(l)
    jshift=ey(l)
    kshift=ez(l)
    do k=minz,maxz
      do j=miny,maxy
        do i=minx,maxx
          if(i==1.or.j==1.or.k==1.or.i==nx.or.j==ny.or.k==nz)then
            itemp=i+ishift
            jtemp=j+jshift
            ktemp=k+kshift
            itemp2=i+ishift
            jtemp2=j+jshift
            ktemp2=k+kshift
            if(itemp<minx.or.jtemp<miny.or.ktemp<minz.or. &
             itemp>maxx.or.jtemp>maxy.or.ktemp>maxz)then
              !apply periodic conditions if necessary
              itemp=pimage(ixpbc,itemp,nx)
              jtemp=pimage(iypbc,jtemp,ny)
              ktemp=pimage(izpbc,ktemp,nz)
              i4=i4back(itemp,jtemp,ktemp) 
              i4orig=i4back(itemp2,jtemp2,ktemp2)
              !if(ownern(i4).EQ.idrank.AND.ownern(i4orig).EQ.idrank &
              ! .AND. isfluid(i,j,k)/=3) THEN
              if(ownern(i4).EQ.idrank &
               .AND. isfluid(i,j,k)/=3 .AND. isfluid(i,j,k)/=0) THEN
                npoplistbc=npoplistbc+1
              endif
            endif
          endif
        enddo
      enddo
    enddo
  enddo
  if(lverbose)write(6,*)'id=',idrank,'npoplistbc=',npoplistbc
  
  if(allocated(poplistbc))deallocate(poplistbc)
  allocate(poplistbc(6,npoplistbc))
  allocate(ipoplistbc(npoplistbc))
  npoplistbc=0
  
  do l=1,links
    ishift=ex(l)
    jshift=ey(l)
    kshift=ez(l)
    do k=minz,maxz
      do j=miny,maxy
        do i=minx,maxx
          if(i==1.or.j==1.or.k==1.or.i==nx.or.j==ny.or.k==nz)then
            itemp=i+ishift
            jtemp=j+jshift
            ktemp=k+kshift
            itemp2=i+ishift
            jtemp2=j+jshift
            ktemp2=k+kshift
            if(itemp<minx.or.jtemp<miny.or.ktemp<minz.or. &
             itemp>maxx.or.jtemp>maxy.or.ktemp>maxz)then
              !apply periodic conditions if necessary
              itemp=pimage(ixpbc,itemp,nx)
              jtemp=pimage(iypbc,jtemp,ny)
              ktemp=pimage(izpbc,ktemp,nz)
              i4=i4back(itemp,jtemp,ktemp) 
              i4orig=i4back(itemp2,jtemp2,ktemp2)
              !if(ownern(i4).EQ.idrank.AND.ownern(i4orig).EQ.idrank &
              if(ownern(i4).EQ.idrank &
               .AND. isfluid(i,j,k)/=3 .AND. isfluid(i,j,k)/=0) THEN
                npoplistbc=npoplistbc+1
                ipoplistbc(npoplistbc)=l
                poplistbc(1,npoplistbc)=itemp
                poplistbc(2,npoplistbc)=jtemp
                poplistbc(3,npoplistbc)=ktemp
                poplistbc(4,npoplistbc)=itemp2
                poplistbc(5,npoplistbc)=jtemp2
                poplistbc(6,npoplistbc)=ktemp2
              endif
            endif
          endif
        enddo
      enddo
    enddo
  enddo
   
  return
   
 end subroutine initialiaze_manage_bc_pop_selfcomm
 
 subroutine manage_bc_pop_selfcomm(aoptp,lparticles)
 
!***********************************************************************
!     
!     LBsoft subroutine to manage the population buffer fluid nodes
!     within the same process for applying the boundary conditions
!     using the node list created in subroutine
!     initialiaze_manage_bc_pop_selfcomm
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  logical, intent(in) :: lparticles
  
  type(REALPTR), dimension(0:links):: aoptp
  
  integer :: i,itemp,jtemp,ktemp,itemp2,jtemp2,ktemp2,idir,iopp
  
  
  if(lparticles)then
    do i=1,npoplistbc
      itemp = poplistbc(1,i)
      jtemp = poplistbc(2,i)
      ktemp = poplistbc(3,i)
      itemp2 = poplistbc(4,i)
      jtemp2 = poplistbc(5,i)
      ktemp2 = poplistbc(6,i)
      idir=ipoplistbc(i)
      iopp=opp(idir)
      if(isfluid(itemp+ex(iopp),jtemp+ey(iopp),ktemp+ez(iopp))<2 .or. &
       isfluid(itemp+ex(iopp),jtemp+ey(iopp),ktemp+ez(iopp))>5)then
        aoptp(idir)%p(itemp,jtemp,ktemp)=aoptp(idir)%p(itemp2,jtemp2,ktemp2)
      endif
!      if(isfluid(itemp2+ex(iopp),jtemp2+ey(iopp),ktemp2+ez(iopp))==1)then
!        aoptp(idir)%p(itemp+ex(iopp),jtemp+ey(iopp),ktemp+ez(iopp))= &
!         aoptp(idir)%p(itemp2+ex(iopp),jtemp2+ey(iopp),ktemp2+ez(iopp))
!      endif
    enddo
  else
    forall(i=1:npoplistbc)
      aoptp(ipoplistbc(i))%p(poplistbc(1,i),poplistbc(2,i),poplistbc(3,i))= &
      real(aoptp(ipoplistbc(i))%p(poplistbc(4,i),poplistbc(5,i),poplistbc(6,i)),kind=PRC)
    end forall
  endif
  
  return
   
 end subroutine manage_bc_pop_selfcomm
 
 pure function pimage(ipbcsub,i,nssub)
 
!***********************************************************************
!     
!     LBsoft sfunction to impose the pbc 
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: ipbcsub,i,nssub
  integer :: pimage
  
  pimage=i
  
  if(ipbcsub==1)then
    if(i<1) then
      pimage=i+nssub
    endif
    if(i>nssub) then
      pimage=i-nssub
    endif
  endif
  
  return
  
 end function pimage

 
 subroutine driver_bc_pops()
  implicit none
 
 
 call commexch_single_halo(f00R)
 call manage_bc_singlehalo_selfcomm(f00R)
 call commwait_single_halo(f00R)
 
 call commexch_single_halo(f01R)
 call manage_bc_singlehalo_selfcomm(f01R)
 call commwait_single_halo(f01R)
 
 call commexch_single_halo(f02R)
 call manage_bc_singlehalo_selfcomm(f02R)
 call commwait_single_halo(f02R)
 
 call commexch_single_halo(f03R)
 call manage_bc_singlehalo_selfcomm(f03R)
 call commwait_single_halo(f03R)
 
 call commexch_single_halo(f04R)
 call manage_bc_singlehalo_selfcomm(f04R)
 call commwait_single_halo(f04R)
 
 call commexch_single_halo(f05R)
 call manage_bc_singlehalo_selfcomm(f05R)
 call commwait_single_halo(f05R)
 
 call commexch_single_halo(f06R)
 call manage_bc_singlehalo_selfcomm(f06R)
 call commwait_single_halo(f06R)
 
 call commexch_single_halo(f07R)
 call manage_bc_singlehalo_selfcomm(f07R)
 call commwait_single_halo(f07R)
 
 call commexch_single_halo(f08R)
 call manage_bc_singlehalo_selfcomm(f08R)
 call commwait_single_halo(f08R)
 
 call commexch_single_halo(f09R)
 call manage_bc_singlehalo_selfcomm(f09R)
 call commwait_single_halo(f09R)
 
 call commexch_single_halo(f10R)
 call manage_bc_singlehalo_selfcomm(f10R)
 call commwait_single_halo(f10R)
 
 call commexch_single_halo(f11R)
 call manage_bc_singlehalo_selfcomm(f11R)
 call commwait_single_halo(f11R)
 
 call commexch_single_halo(f12R)
 call manage_bc_singlehalo_selfcomm(f12R)
 call commwait_single_halo(f12R)
 
 call commexch_single_halo(f13R)
 call manage_bc_singlehalo_selfcomm(f13R)
 call commwait_single_halo(f13R)
 
 call commexch_single_halo(f14R)
 call manage_bc_singlehalo_selfcomm(f14R)
 call commwait_single_halo(f14R)
 
 call commexch_single_halo(f15R)
 call manage_bc_singlehalo_selfcomm(f15R)
 call commwait_single_halo(f15R)
 
 call commexch_single_halo(f16R)
 call manage_bc_singlehalo_selfcomm(f16R)
 call commwait_single_halo(f16R)
 
 call commexch_single_halo(f17R)
 call manage_bc_singlehalo_selfcomm(f17R)
 call commwait_single_halo(f17R)
 
 call commexch_single_halo(f18R)
 call manage_bc_singlehalo_selfcomm(f18R)
 call commwait_single_halo(f18R)
 
 if(lsingle_fluid)return
 
 call commexch_single_halo(f00B)
 call manage_bc_singlehalo_selfcomm(f00B)
 call commwait_single_halo(f00B)
 
 call commexch_single_halo(f01B)
 call manage_bc_singlehalo_selfcomm(f01B)
 call commwait_single_halo(f01B)
 
 call commexch_single_halo(f02B)
 call manage_bc_singlehalo_selfcomm(f02B)
 call commwait_single_halo(f02B)
 
 call commexch_single_halo(f03B)
 call manage_bc_singlehalo_selfcomm(f03B)
 call commwait_single_halo(f03B)
 
 call commexch_single_halo(f04B)
 call manage_bc_singlehalo_selfcomm(f04B)
 call commwait_single_halo(f04B)
 
 call commexch_single_halo(f05B)
 call manage_bc_singlehalo_selfcomm(f05B)
 call commwait_single_halo(f05B)
 
 call commexch_single_halo(f06B)
 call manage_bc_singlehalo_selfcomm(f06B)
 call commwait_single_halo(f06B)
 
 call commexch_single_halo(f07B)
 call manage_bc_singlehalo_selfcomm(f07B)
 call commwait_single_halo(f07B)
 
 call commexch_single_halo(f08B)
 call manage_bc_singlehalo_selfcomm(f08B)
 call commwait_single_halo(f08B)
 
 call commexch_single_halo(f09B)
 call manage_bc_singlehalo_selfcomm(f09B)
 call commwait_single_halo(f09B)
 
 call commexch_single_halo(f10B)
 call manage_bc_singlehalo_selfcomm(f10B)
 call commwait_single_halo(f10B)
 
 call commexch_single_halo(f11B)
 call manage_bc_singlehalo_selfcomm(f11B)
 call commwait_single_halo(f11B)
 
 call commexch_single_halo(f12B)
 call manage_bc_singlehalo_selfcomm(f12B)
 call commwait_single_halo(f12B)
 
 call commexch_single_halo(f13B)
 call manage_bc_singlehalo_selfcomm(f13B)
 call commwait_single_halo(f13B)
 
 call commexch_single_halo(f14B)
 call manage_bc_singlehalo_selfcomm(f14B)
 call commwait_single_halo(f14B)
 
 call commexch_single_halo(f15B)
 call manage_bc_singlehalo_selfcomm(f15B)
 call commwait_single_halo(f15B)
 
 call commexch_single_halo(f16B)
 call manage_bc_singlehalo_selfcomm(f16B)
 call commwait_single_halo(f16B)
 
 call commexch_single_halo(f17B)
 call manage_bc_singlehalo_selfcomm(f17B)
 call commwait_single_halo(f17B)
 
 call commexch_single_halo(f18B)
 call manage_bc_singlehalo_selfcomm(f18B)
 call commwait_single_halo(f18B)

  
 end subroutine driver_bc_pops




!******************END PART TO MANAGE THE PERIODIC BC*******************
 
!*****************START PART TO MANAGE THE BOUNCEBACK*******************
 
 subroutine driver_apply_bounceback_pop
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the bounce back of the fluid
!     populations if requested from the boundary conditions
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  integer ::i,j,k,l
  integer, save :: iter=0
  
  call error(41)
  
  iter=iter+1
  
  call set_bc_variable_hvar
  

  call apply_bounceback_pop(bc_rhoR,bc_u,bc_v,bc_w,aoptpR)


  
#ifdef DIAGNSTREAM
  if(iter==NDIAGNSTREAM)call print_all_pops(100,'miodopobounce',iter,aoptpR)
#endif
  
  if(lsingle_fluid)return
  
  call apply_bounceback_pop(bc_rhoB,bc_u,bc_v,bc_w,aoptpB)

  return
  
 end subroutine driver_apply_bounceback_pop
 
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
  
  call set_bc_variable_hvar
  
  
  call apply_bounceback_pop_halfway(nstep,bc_rhoR,bc_u,bc_v,bc_w,aoptpR, &
   alphaR_CG)
  
#ifdef DIAGNSTREAM
  if(iter==NDIAGNSTREAM)call print_all_pops(100,'miodopobounce',iter,aoptpR)
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
  
  call apply_bounceback_pop_halfway(nstep,bc_rhoB,bc_u,bc_v,bc_w,aoptpB, &
   alphaB_CG)
  
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
 
 subroutine set_bc_hvar
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the bc value of the hydrodynamic 
!     variable if requested at halfway grid point
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
 
  implicit none
  
  integer :: i,j,k,idir,inits,ends,ishift,jshift,kshift
  integer :: ishift2,jshift2,kshift2
  integer :: ishift3,jshift3,kshift3
  
  !INITIALIZE isfluid=6:9 ; bcfluid from 1 to 6
  inits=nbounce0+1
  ends=nbounce8
  do i=inits,ends
    bc_rhoR(i)=ZERO
    bc_u(i)=ZERO
    bc_v(i)=ZERO
    bc_w(i)=ZERO
  enddo
  if(.not.lsingle_fluid)then
    forall(i=inits:ends)
      bc_rhoB(i)=ZERO
    end forall
  endif
  
  !fixed bc value
  
  !isfluid=6 ; bcfluid from 1 to 6
  ! dirichlet condition

  inits=nbounce6dir(0)+1
  ends=nbounce6dir(1)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_rhoR(i)=bc_rhoR_west
    end forall
  endif
  
  inits=nbounce6dir(1)+1
  ends=nbounce6dir(2)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_rhoR(i)=bc_rhoR_east
    end forall
  endif
  
  inits=nbounce6dir(2)+1
  ends=nbounce6dir(3)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_rhoR(i)=bc_rhoR_rear
    end forall
  endif
  
  inits=nbounce6dir(3)+1
  ends=nbounce6dir(4)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_rhoR(i)=bc_rhoR_front
    end forall
  endif
  
  inits=nbounce6dir(4)+1
  ends=nbounce6dir(5)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_rhoR(i)=bc_rhoR_south
    end forall
  endif
  
  inits=nbounce6dir(5)+1
  ends=nbounce6dir(6)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_rhoR(i)=bc_rhoR_north
    end forall
  endif
    
  if(.not. lsingle_fluid)then
    
    !isfluid=6 ; bcfluid from 1 to 6
    ! dirichlet condition
    inits=nbounce6dir(0)+1
    ends=nbounce6dir(1)
    if(inits<=ends)then
      forall(i=inits:ends)
        bc_rhoB(i)=bc_rhoB_west !bc_rhoB_east
      end forall
    endif
    
    inits=nbounce6dir(1)+1
    ends=nbounce6dir(2)
    if(inits<=ends)then
      forall(i=inits:ends)
        bc_rhoB(i)=bc_rhoB_east !bc_rhoB_west
      end forall
    endif
    
    inits=nbounce6dir(2)+1
    ends=nbounce6dir(3)
    if(inits<=ends)then
      forall(i=inits:ends)
        bc_rhoB(i)=bc_rhoB_rear !bc_rhoB_front
      end forall
    endif
    
    inits=nbounce6dir(3)+1
    ends=nbounce6dir(4)
    if(inits<=ends)then
      forall(i=inits:ends)
        bc_rhoB(i)=bc_rhoB_front !bc_rhoB_rear
      end forall
    endif
     
    inits=nbounce6dir(4)+1
    ends=nbounce6dir(5)
    if(inits<=ends)then
      forall(i=inits:ends)
        bc_rhoB(i)=bc_rhoB_south !bc_rhoB_north
      end forall
    endif
    
    inits=nbounce6dir(5)+1
    ends=nbounce6dir(6)
    if(inits<=ends)then
      forall(i=inits:ends)
        bc_rhoB(i)=bc_rhoB_north !bc_rhoB_south
      end forall
    endif
  
  endif
  
  ! neumann condition
  
  inits=nbounce7dir(0)+1
  ends=nbounce7dir(1)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_u(i)=bc_u_west
      bc_v(i)=bc_v_west
      bc_w(i)=bc_w_west
    end forall
  endif
  
  inits=nbounce7dir(1)+1
  ends=nbounce7dir(2)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_u(i)=bc_u_east
      bc_v(i)=bc_v_east
      bc_w(i)=bc_w_east
    end forall
  endif
  
  inits=nbounce7dir(2)+1
  ends=nbounce7dir(3)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_u(i)=bc_u_rear
      bc_v(i)=bc_v_rear
      bc_w(i)=bc_w_rear
    end forall
  endif
  
  inits=nbounce7dir(3)+1
  ends=nbounce7dir(4)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_u(i)=bc_u_front
      bc_v(i)=bc_v_front
      bc_w(i)=bc_w_front
    end forall
  endif
  
  inits=nbounce7dir(4)+1
  ends=nbounce7dir(5)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_u(i)=bc_u_south
      bc_v(i)=bc_v_south
      bc_w(i)=bc_w_south
    end forall
  endif
  
  inits=nbounce7dir(5)+1
  ends=nbounce7dir(6)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_u(i)=bc_u_north
      bc_v(i)=bc_v_north
      bc_w(i)=bc_w_north
    end forall
  endif
  
  inits=nbounce8dir(0)+1
  ends=nbounce8dir(1)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_rhoR(i)=bc_rhoR_west
      bc_u(i)=bc_u_west
      bc_v(i)=bc_v_west
      bc_w(i)=bc_w_west
    end forall
  endif
  
  inits=nbounce8dir(1)+1
  ends=nbounce8dir(2)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_rhoR(i)=bc_rhoR_east
      bc_u(i)=bc_u_east
      bc_v(i)=bc_v_east
      bc_w(i)=bc_w_east
    end forall
  endif
  
  inits=nbounce8dir(2)+1
  ends=nbounce8dir(3)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_rhoB(i)=bc_rhoR_rear
      bc_u(i)=bc_u_rear
      bc_v(i)=bc_v_rear
      bc_w(i)=bc_w_rear
    end forall
  endif
  
  inits=nbounce8dir(3)+1
  ends=nbounce8dir(4)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_rhoR(i)=bc_rhoR_front
      bc_u(i)=bc_u_front
      bc_v(i)=bc_v_front
      bc_w(i)=bc_w_front
    end forall
  endif
  
  inits=nbounce8dir(4)+1
  ends=nbounce8dir(5)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_rhoR(i)=bc_rhoR_south
      bc_u(i)=bc_u_south
      bc_v(i)=bc_v_south
      bc_w(i)=bc_w_south
    end forall
  endif
  
  inits=nbounce8dir(5)+1
  ends=nbounce8dir(6)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_rhoR(i)=bc_rhoR_north
      bc_u(i)=bc_u_north
      bc_v(i)=bc_v_north
      bc_w(i)=bc_w_north
    end forall
  endif
  
  if(.not. lsingle_fluid)then
    
    !isfluid=8 ; bcfluid from 1 to 6
    ! dirichlet condition
    inits=nbounce8dir(0)+1
    ends=nbounce8dir(1)
    if(inits<=ends)then
      forall(i=inits:ends)
        
        bc_rhoB(i)=bc_rhoB_west
      end forall
    endif
    
    inits=nbounce8dir(1)+1
    ends=nbounce8dir(2)
    if(inits<=ends)then
      forall(i=inits:ends)
        bc_rhoB(i)=bc_rhoB_east
      end forall
    endif
    
    inits=nbounce8dir(2)+1
    ends=nbounce8dir(3)
    if(inits<=ends)then
      forall(i=inits:ends)
        bc_rhoB(i)=bc_rhoB_rear
      end forall
    endif
    
    inits=nbounce8dir(3)+1
    ends=nbounce8dir(4)
    if(inits<=ends)then
      forall(i=inits:ends)
        bc_rhoB(i)=bc_rhoB_front
      end forall
    endif
     
    inits=nbounce8dir(4)+1
    ends=nbounce8dir(5)
    if(inits<=ends)then
      forall(i=inits:ends)
        bc_rhoB(i)=bc_rhoB_south
      end forall
    endif
    
    inits=nbounce8dir(5)+1
    ends=nbounce8dir(6)
    if(inits<=ends)then
      forall(i=inits:ends)
        bc_rhoB(i)=bc_rhoB_north
      end forall
    endif
    
  endif
  
  !isfluid=9 ; bcfluid from 1 to 6
  ! flow condition
  inits=nbounce9dir(0)+1
  ends=nbounce9dir(1)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_flow(i)=bc_flow_west
    end forall
  endif
    
  inits=nbounce9dir(1)+1
  ends=nbounce9dir(2)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_flow(i)=bc_flow_east
    end forall
  endif
    
  inits=nbounce9dir(2)+1
  ends=nbounce9dir(3)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_flow(i)=bc_flow_rear
    end forall
  endif
    
  inits=nbounce9dir(3)+1
  ends=nbounce9dir(4)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_flow(i)=bc_flow_front
    end forall
  endif
     
  inits=nbounce9dir(4)+1
  ends=nbounce9dir(5)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_flow(i)=bc_flow_south
    end forall
  endif
    
  inits=nbounce9dir(5)+1
  ends=nbounce9dir(6)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_flow(i)=bc_flow_north
    end forall
  endif
  
  !not fixed bc value
  
  do idir=1,nbcdir
    inits=nbounce6dir(idir-1)+1
    ends=nbounce6dir(idir)
    if(inits>ends)cycle
    ishift=ex(idir)
    jshift=ey(idir)
    kshift=ez(idir)
    ishift2=ex(idir)*2
    jshift2=ey(idir)*2
    kshift2=ez(idir)*2
    ishift3=ex(idir)*3
    jshift3=ey(idir)*3
    kshift3=ez(idir)*3
    if(idir==1 .or. idir==2)then
      do i=inits,ends
        bc_u(i)=interpolation_order_2_hf( &
         u(ibounce(1,i)+ishift,ibounce(2,i)+jshift,ibounce(3,i)+kshift), &
         u(ibounce(1,i)+ishift2,ibounce(2,i)+jshift2,ibounce(3,i)+kshift2), &
         u(ibounce(1,i)+ishift3,ibounce(2,i)+jshift3,ibounce(3,i)+kshift3))
      enddo
      where(bc_u(inits:ends)>cssq*HALF)
        bc_u(inits:ends)=cssq*HALF
      elsewhere(bc_u(inits:ends)<-cssq*HALF)
        bc_u(inits:ends)=-cssq*HALF
      end where
    endif
    if(idir==3 .or. idir==4)then
      do i=inits,ends
        bc_v(i)=interpolation_order_2_hf( &
         v(ibounce(1,i)+ishift,ibounce(2,i)+jshift,ibounce(3,i)+kshift), &
         v(ibounce(1,i)+ishift2,ibounce(2,i)+jshift2,ibounce(3,i)+kshift2), &
         v(ibounce(1,i)+ishift3,ibounce(2,i)+jshift3,ibounce(3,i)+kshift3))
      end do
      where(bc_v(inits:ends)>cssq*HALF)
        bc_v(inits:ends)=cssq*HALF
      elsewhere(bc_v(inits:ends)<-cssq*HALF)
        bc_v(inits:ends)=-cssq*HALF
      end where
    endif
    if(idir==5 .or.idir==6)then
      do i=inits,ends
        bc_w(i)=interpolation_order_2_hf( &
         w(ibounce(1,i)+ishift,ibounce(2,i)+jshift,ibounce(3,i)+kshift), &
         w(ibounce(1,i)+ishift2,ibounce(2,i)+jshift2,ibounce(3,i)+kshift2), &
         w(ibounce(1,i)+ishift3,ibounce(2,i)+jshift3,ibounce(3,i)+kshift3))
      enddo
      where(bc_w(inits:ends)>cssq*HALF)
        bc_w(inits:ends)=cssq*HALF
      elsewhere(bc_w(inits:ends)<-cssq*HALF)
        bc_w(inits:ends)=-cssq*HALF
      end where
    endif
  enddo
  
  
  do idir=1,nbcdir
    inits=nbounce7dir(idir-1)+1
    ends=nbounce7dir(idir)
    if(inits>ends)cycle
    ishift=ex(idir)
    jshift=ey(idir)
    kshift=ez(idir)
    ishift2=ex(idir)*2
    jshift2=ey(idir)*2
    kshift2=ez(idir)*2
    ishift3=ex(idir)*3
    jshift3=ey(idir)*3
    kshift3=ez(idir)*3
    if(lsingle_fluid)then
      do i=inits,ends
        bc_rhoR(i)=interpolation_order_2_hf( &
         rhoR(ibounce(1,i)+ishift,ibounce(2,i)+jshift,ibounce(3,i)+kshift), &
         rhoR(ibounce(1,i)+ishift2,ibounce(2,i)+jshift2,ibounce(3,i)+kshift2), &
         rhoR(ibounce(1,i)+ishift3,ibounce(2,i)+jshift3,ibounce(3,i)+kshift3))
      enddo
    else
      do i=inits,ends
        bc_rhoR(i)=interpolation_order_2_hf( &
         rhoR(ibounce(1,i)+ishift,ibounce(2,i)+jshift,ibounce(3,i)+kshift), &
         rhoR(ibounce(1,i)+ishift2,ibounce(2,i)+jshift2,ibounce(3,i)+kshift2), &
         rhoR(ibounce(1,i)+ishift3,ibounce(2,i)+jshift3,ibounce(3,i)+kshift3))
        bc_rhoB(i)=interpolation_order_2_hf( &
         rhoB(ibounce(1,i)+ishift,ibounce(2,i)+jshift,ibounce(3,i)+kshift), &
         rhoB(ibounce(1,i)+ishift2,ibounce(2,i)+jshift2,ibounce(3,i)+kshift2), &
         rhoB(ibounce(1,i)+ishift3,ibounce(2,i)+jshift3,ibounce(3,i)+kshift3))
      enddo
    endif
  enddo
  
  do idir=1,nbcdir
    inits=nbounce8dir(idir-1)+1
    ends=nbounce8dir(idir)
    if(inits>ends)cycle
    ishift=ex(idir)
    jshift=ey(idir)
    kshift=ez(idir)
!    ishift2=ex(idir)*2
!    jshift2=ey(idir)*2
!    kshift2=ez(idir)*2
!    ishift3=ex(idir)*3
!    jshift3=ey(idir)*3
!    kshift3=ez(idir)*3
    if(lunique_omega)then
      do i=inits,ends
        bc_omega(i)=unique_omega
      enddo
    else
      do i=inits,ends
        bc_omega(i)=omega(ibounce(1,i)+ishift,ibounce(2,i)+jshift,ibounce(3,i)+kshift)
      enddo
    endif
  enddo
  
  do idir=1,nbcdir
    inits=nbounce9dir(idir-1)+1
    ends=nbounce9dir(idir)
    if(inits>ends)cycle
    ishift=ex(idir)
    jshift=ey(idir)
    kshift=ez(idir)
    ishift2=ex(idir)*2
    jshift2=ey(idir)*2
    kshift2=ez(idir)*2
    ishift3=ex(idir)*3
    jshift3=ey(idir)*3
    kshift3=ez(idir)*3
    if(lsingle_fluid)then
      do i=inits,ends
        bc_rhoR(i)=interpolation_order_2_hf( &
         rhoR(ibounce(1,i)+ishift,ibounce(2,i)+jshift,ibounce(3,i)+kshift), &
         rhoR(ibounce(1,i)+ishift2,ibounce(2,i)+jshift2,ibounce(3,i)+kshift2), &
         rhoR(ibounce(1,i)+ishift3,ibounce(2,i)+jshift3,ibounce(3,i)+kshift3))
      enddo
      if(idir==1 .or. idir==2)then
        do i=inits,ends
          bc_u(i)=bc_flow(i)/bc_rhoR(i)
        enddo
        where(bc_u(inits:ends)>cssq*HALF)
          bc_u(inits:ends)=cssq*HALF
        elsewhere(bc_u(inits:ends)<-cssq*HALF)
          bc_u(inits:ends)=-cssq*HALF
        end where
      endif
      if(idir==3 .or. idir==4)then
        do i=inits,ends
          bc_v(i)=bc_flow(i)/bc_rhoR(i)
        end do
        where(bc_v(inits:ends)>cssq*HALF)
          bc_v(inits:ends)=cssq*HALF
        elsewhere(bc_v(inits:ends)<-cssq*HALF)
          bc_v(inits:ends)=-cssq*HALF
        end where
      endif
      if(idir==5 .or.idir==6)then
        do i=inits,ends
          bc_w(i)=-bc_flow(i)/bc_rhoR(i)
        enddo
        where(bc_w(inits:ends)>cssq*HALF)
          bc_w(inits:ends)=cssq*HALF
        elsewhere(bc_w(inits:ends)<-cssq*HALF)
          bc_w(inits:ends)=-cssq*HALF
        end where
      endif
    else
      do i=inits,ends
        bc_rhoR(i)=interpolation_order_2_hf( &
         rhoR(ibounce(1,i)+ishift,ibounce(2,i)+jshift,ibounce(3,i)+kshift), &
         rhoR(ibounce(1,i)+ishift2,ibounce(2,i)+jshift2,ibounce(3,i)+kshift2), &
         rhoR(ibounce(1,i)+ishift3,ibounce(2,i)+jshift3,ibounce(3,i)+kshift3))
        bc_rhoB(i)=interpolation_order_2_hf( &
         rhoB(ibounce(1,i)+ishift,ibounce(2,i)+jshift,ibounce(3,i)+kshift), &
         rhoB(ibounce(1,i)+ishift2,ibounce(2,i)+jshift2,ibounce(3,i)+kshift2), &
         rhoB(ibounce(1,i)+ishift3,ibounce(2,i)+jshift3,ibounce(3,i)+kshift3))
      enddo
      if(idir==1 .or. idir==2)then
        do i=inits,ends
          bc_u(i)=bc_flow(i)/(bc_rhoR(i)+bc_rhoB(i))
        enddo
        where(bc_u(inits:ends)>cssq*HALF)
          bc_u(inits:ends)=cssq*HALF
        elsewhere(bc_u(inits:ends)<-cssq*HALF)
          bc_u(inits:ends)=-cssq*HALF
        end where
      endif
      if(idir==3 .or. idir==4)then
        do i=inits,ends
          bc_v(i)=bc_flow(i)/(bc_rhoR(i)+bc_rhoB(i))
        end do
        where(bc_v(inits:ends)>cssq*HALF)
          bc_v(inits:ends)=cssq*HALF
        elsewhere(bc_v(inits:ends)<-cssq*HALF)
          bc_v(inits:ends)=-cssq*HALF
        end where
      endif
      if(idir==5 .or.idir==6)then
        do i=inits,ends
          bc_w(i)=-bc_flow(i)/(bc_rhoR(i)+bc_rhoB(i))
        enddo
        where(bc_w(inits:ends)>cssq*HALF)
          bc_w(inits:ends)=cssq*HALF
        elsewhere(bc_w(inits:ends)<-cssq*HALF)
          bc_w(inits:ends)=-cssq*HALF
        end where
      endif
    endif
  enddo
  
  return
  
 end subroutine set_bc_hvar
 
 subroutine set_bc_fixed_hvar
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the fixed bc value of the hydrodynamic 
!     variable if requested (fixed = which were set at the beginning)
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification September 2018
!     
!***********************************************************************
 
  implicit none
  
  integer :: i,j,k,idir,inits,ends
  
  !fixed bc value
  
  !isfluid=6 ; bcfluid from 1 to 6
  ! dirichlet condition
  inits=nbounce6dir(0)+1
  ends=nbounce6dir(1)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_rhoR(i)=bc_rhoR_west
    end forall
  endif
  
  inits=nbounce6dir(1)+1
  ends=nbounce6dir(2)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_rhoR(i)=bc_rhoR_east
    end forall
  endif
  
  inits=nbounce6dir(2)+1
  ends=nbounce6dir(3)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_rhoR(i)=bc_rhoR_rear
    end forall
  endif
  
  inits=nbounce6dir(3)+1
  ends=nbounce6dir(4)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_rhoR(i)=bc_rhoR_front
    end forall
  endif
  
  inits=nbounce6dir(4)+1
  ends=nbounce6dir(5)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_rhoR(i)=bc_rhoR_south
    end forall
  endif
  
  inits=nbounce6dir(5)+1
  ends=nbounce6dir(6)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_rhoR(i)=bc_rhoR_north
    end forall
  endif
    
  if(.not. lsingle_fluid)then
    
    !isfluid=6 ; bcfluid from 1 to 6
    ! dirichlet condition
    inits=nbounce6dir(0)+1
    ends=nbounce6dir(1)
    if(inits<=ends)then
      forall(i=inits:ends)
        bc_rhoB(i)=bc_rhoB_west !bc_rhoB_east
      end forall
    endif
    
    inits=nbounce6dir(1)+1
    ends=nbounce6dir(2)
    if(inits<=ends)then
      forall(i=inits:ends)
        bc_rhoB(i)=bc_rhoB_east !bc_rhoB_west
      end forall
    endif
    
    inits=nbounce6dir(2)+1
    ends=nbounce6dir(3)
    if(inits<=ends)then
      forall(i=inits:ends)
        bc_rhoB(i)=bc_rhoB_rear !bc_rhoB_front
      end forall
    endif
    
    inits=nbounce6dir(3)+1
    ends=nbounce6dir(4)
    if(inits<=ends)then
      forall(i=inits:ends)
        bc_rhoB(i)=bc_rhoB_front !bc_rhoB_rear
      end forall
    endif
     
    inits=nbounce6dir(4)+1
    ends=nbounce6dir(5)
    if(inits<=ends)then
      forall(i=inits:ends)
        bc_rhoB(i)=bc_rhoB_south !bc_rhoB_north
      end forall
    endif
    
    inits=nbounce6dir(5)+1
    ends=nbounce6dir(6)
    if(inits<=ends)then
      forall(i=inits:ends)
        bc_rhoB(i)=bc_rhoB_north !bc_rhoB_south
      end forall
    endif
  
  endif
  
  ! neumann condition
  
  inits=nbounce7dir(0)+1
  ends=nbounce7dir(1)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_u(i)=bc_u_west
      bc_v(i)=bc_v_west
      bc_w(i)=bc_w_west
    end forall
  endif
  
  inits=nbounce7dir(1)+1
  ends=nbounce7dir(2)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_u(i)=bc_u_east
      bc_v(i)=bc_v_east
      bc_w(i)=bc_w_east
    end forall
  endif
  
  inits=nbounce7dir(2)+1
  ends=nbounce7dir(3)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_u(i)=bc_u_rear
      bc_v(i)=bc_v_rear
      bc_w(i)=bc_w_rear
    end forall
  endif
  
  inits=nbounce7dir(3)+1
  ends=nbounce7dir(4)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_u(i)=bc_u_front
      bc_v(i)=bc_v_front
      bc_w(i)=bc_w_front
    end forall
  endif
  
  inits=nbounce7dir(4)+1
  ends=nbounce7dir(5)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_u(i)=bc_u_south
      bc_v(i)=bc_v_south
      bc_w(i)=bc_w_south
    end forall
  endif
  
  inits=nbounce7dir(5)+1
  ends=nbounce7dir(6)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_u(i)=bc_u_north
      bc_v(i)=bc_v_north
      bc_w(i)=bc_w_north
    end forall
  endif
  
  inits=nbounce8dir(0)+1
  ends=nbounce8dir(1)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_rhoR(i)=bc_rhoR_west
      bc_u(i)=bc_u_west
      bc_v(i)=bc_v_west
      bc_w(i)=bc_w_west
    end forall
  endif
  
  inits=nbounce8dir(1)+1
  ends=nbounce8dir(2)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_rhoR(i)=bc_rhoR_east
      bc_u(i)=bc_u_east
      bc_v(i)=bc_v_east
      bc_w(i)=bc_w_east
    end forall
  endif
  
  inits=nbounce8dir(2)+1
  ends=nbounce8dir(3)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_rhoB(i)=bc_rhoR_rear
      bc_u(i)=bc_u_rear
      bc_v(i)=bc_v_rear
      bc_w(i)=bc_w_rear
    end forall
  endif
  
  inits=nbounce8dir(3)+1
  ends=nbounce8dir(4)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_rhoR(i)=bc_rhoR_front
      bc_u(i)=bc_u_front
      bc_v(i)=bc_v_front
      bc_w(i)=bc_w_front
    end forall
  endif
  
  inits=nbounce8dir(4)+1
  ends=nbounce8dir(5)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_rhoR(i)=bc_rhoR_south
      bc_u(i)=bc_u_south
      bc_v(i)=bc_v_south
      bc_w(i)=bc_w_south
    end forall
  endif
  
  inits=nbounce8dir(5)+1
  ends=nbounce8dir(6)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_rhoR(i)=bc_rhoR_north
      bc_u(i)=bc_u_north
      bc_v(i)=bc_v_north
      bc_w(i)=bc_w_north
    end forall
  endif
  
  if(.not. lsingle_fluid)then
    
    !isfluid=8 ; bcfluid from 1 to 6
    ! dirichlet condition
    inits=nbounce8dir(0)+1
    ends=nbounce8dir(1)
    if(inits<=ends)then
      forall(i=inits:ends)
        
        bc_rhoB(i)=bc_rhoB_west
      end forall
    endif
    
    inits=nbounce8dir(1)+1
    ends=nbounce8dir(2)
    if(inits<=ends)then
      forall(i=inits:ends)
        bc_rhoB(i)=bc_rhoB_east
      end forall
    endif
    
    inits=nbounce8dir(2)+1
    ends=nbounce8dir(3)
    if(inits<=ends)then
      forall(i=inits:ends)
        bc_rhoB(i)=bc_rhoB_rear
      end forall
    endif
    
    inits=nbounce8dir(3)+1
    ends=nbounce8dir(4)
    if(inits<=ends)then
      forall(i=inits:ends)
        bc_rhoB(i)=bc_rhoB_front
      end forall
    endif
     
    inits=nbounce8dir(4)+1
    ends=nbounce8dir(5)
    if(inits<=ends)then
      forall(i=inits:ends)
        bc_rhoB(i)=bc_rhoB_south
      end forall
    endif
    
    inits=nbounce8dir(5)+1
    ends=nbounce8dir(6)
    if(inits<=ends)then
      forall(i=inits:ends)
        bc_rhoB(i)=bc_rhoB_north
      end forall
    endif
    
  endif
  
  !isfluid=9 ; bcfluid from 1 to 6
  ! flow condition
  inits=nbounce9dir(0)+1
  ends=nbounce9dir(1)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_flow(i)=bc_flow_west
    end forall
  endif
    
  inits=nbounce9dir(1)+1
  ends=nbounce9dir(2)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_flow(i)=bc_flow_east
    end forall
  endif
    
  inits=nbounce9dir(2)+1
  ends=nbounce9dir(3)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_flow(i)=bc_flow_rear
    end forall
  endif
    
  inits=nbounce9dir(3)+1
  ends=nbounce9dir(4)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_flow(i)=bc_flow_front
    end forall
  endif
     
  inits=nbounce9dir(4)+1
  ends=nbounce9dir(5)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_flow(i)=bc_flow_south
    end forall
  endif
    
  inits=nbounce9dir(5)+1
  ends=nbounce9dir(6)
  if(inits<=ends)then
    forall(i=inits:ends)
      bc_flow(i)=bc_flow_north
    end forall
  endif
  
  return
  
 end subroutine set_bc_fixed_hvar
 
 subroutine set_bc_variable_hvar
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the variable bc value of the        
!     hydrodynamic variable if requested (variable = can be change      
!     along the run) at halfway grid point
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
 
  implicit none
  
  integer :: i,j,k,idir,inits,ends,ishift,jshift,kshift
  integer :: ishift2,jshift2,kshift2
  integer :: ishift3,jshift3,kshift3
  
  !not fixed bc value
  
  do idir=1,nbcdir
    inits=nbounce6dir(idir-1)+1
    ends=nbounce6dir(idir)
    if(inits>ends)cycle
    ishift=ex(idir)
    jshift=ey(idir)
    kshift=ez(idir)
    ishift2=ex(idir)*2
    jshift2=ey(idir)*2
    kshift2=ez(idir)*2
    ishift3=ex(idir)*3
    jshift3=ey(idir)*3
    kshift3=ez(idir)*3
    if(idir==1 .or. idir==2)then
      do i=inits,ends
        bc_u(i)=interpolation_order_2_hf( &
         u(ibounce(1,i)+ishift,ibounce(2,i)+jshift,ibounce(3,i)+kshift), &
         u(ibounce(1,i)+ishift2,ibounce(2,i)+jshift2,ibounce(3,i)+kshift2), &
         u(ibounce(1,i)+ishift3,ibounce(2,i)+jshift3,ibounce(3,i)+kshift3))
      enddo
      where(bc_u(inits:ends)>cssq*HALF)
        bc_u(inits:ends)=cssq*HALF
      elsewhere(bc_u(inits:ends)<-cssq*HALF)
        bc_u(inits:ends)=-cssq*HALF
      end where
    endif
    if(idir==3 .or. idir==4)then
      do i=inits,ends
        bc_v(i)=interpolation_order_2_hf( &
         v(ibounce(1,i)+ishift,ibounce(2,i)+jshift,ibounce(3,i)+kshift), &
         v(ibounce(1,i)+ishift2,ibounce(2,i)+jshift2,ibounce(3,i)+kshift2), &
         v(ibounce(1,i)+ishift3,ibounce(2,i)+jshift3,ibounce(3,i)+kshift3))
      enddo
      where(bc_v(inits:ends)>cssq*HALF)
        bc_v(inits:ends)=cssq*HALF
      elsewhere(bc_v(inits:ends)<-cssq*HALF)
        bc_v(inits:ends)=-cssq*HALF
      end where
    endif
    if(idir==5 .or.idir==6)then
      do i=inits,ends
        bc_w(i)=interpolation_order_2_hf( &
         w(ibounce(1,i)+ishift,ibounce(2,i)+jshift,ibounce(3,i)+kshift), &
         w(ibounce(1,i)+ishift2,ibounce(2,i)+jshift2,ibounce(3,i)+kshift2), &
         w(ibounce(1,i)+ishift3,ibounce(2,i)+jshift3,ibounce(3,i)+kshift3))
      enddo
      where(bc_w(inits:ends)>cssq*HALF)
        bc_w(inits:ends)=cssq*HALF
      elsewhere(bc_w(inits:ends)<-cssq*HALF)
        bc_w(inits:ends)=-cssq*HALF
      end where
    endif
  enddo
  
  
  do idir=1,nbcdir
    inits=nbounce7dir(idir-1)+1
    ends=nbounce7dir(idir)
    if(inits>ends)cycle
    ishift=ex(idir)
    jshift=ey(idir)
    kshift=ez(idir)
    ishift2=ex(idir)*2
    jshift2=ey(idir)*2
    kshift2=ez(idir)*2
    ishift3=ex(idir)*3
    jshift3=ey(idir)*3
    kshift3=ez(idir)*3
    if(lsingle_fluid)then
      do i=inits,ends
        bc_rhoR(i)=interpolation_order_2_hf( &
         rhoR(ibounce(1,i)+ishift,ibounce(2,i)+jshift,ibounce(3,i)+kshift), &
         rhoR(ibounce(1,i)+ishift2,ibounce(2,i)+jshift2,ibounce(3,i)+kshift2), &
         rhoR(ibounce(1,i)+ishift3,ibounce(2,i)+jshift3,ibounce(3,i)+kshift3))
      enddo
    else
      do i=inits,ends
        bc_rhoR(i)=interpolation_order_2_hf( &
         rhoR(ibounce(1,i)+ishift,ibounce(2,i)+jshift,ibounce(3,i)+kshift), &
         rhoR(ibounce(1,i)+ishift2,ibounce(2,i)+jshift2,ibounce(3,i)+kshift2), &
         rhoR(ibounce(1,i)+ishift3,ibounce(2,i)+jshift3,ibounce(3,i)+kshift3))
        bc_rhoB(i)=interpolation_order_2_hf( &
         rhoB(ibounce(1,i)+ishift,ibounce(2,i)+jshift,ibounce(3,i)+kshift), &
         rhoB(ibounce(1,i)+ishift2,ibounce(2,i)+jshift2,ibounce(3,i)+kshift2), &
         rhoB(ibounce(1,i)+ishift3,ibounce(2,i)+jshift3,ibounce(3,i)+kshift3))
      enddo
    endif
  enddo
  
  do idir=1,nbcdir
    inits=nbounce8dir(idir-1)+1
    ends=nbounce8dir(idir)
    if(inits>ends)cycle
    ishift=ex(idir)
    jshift=ey(idir)
    kshift=ez(idir)
!    ishift2=ex(idir)*2
!    jshift2=ey(idir)*2
!    kshift2=ez(idir)*2
!    ishift3=ex(idir)*3
!    jshift3=ey(idir)*3
!    kshift3=ez(idir)*3
    if(lunique_omega)then
      do i=inits,ends
        bc_omega(i)=unique_omega
      enddo
    else
      do i=inits,ends
        bc_omega(i)=omega(ibounce(1,i)+ishift,ibounce(2,i)+jshift,ibounce(3,i)+kshift)
      enddo
    endif
  enddo
  
  do idir=1,nbcdir
    inits=nbounce9dir(idir-1)+1
    ends=nbounce9dir(idir)
    if(inits>ends)cycle
    ishift=ex(idir)
    jshift=ey(idir)
    kshift=ez(idir)
    ishift2=ex(idir)*2
    jshift2=ey(idir)*2
    kshift2=ez(idir)*2
    ishift3=ex(idir)*3
    jshift3=ey(idir)*3
    kshift3=ez(idir)*3
    if(lsingle_fluid)then
      do i=inits,ends
        bc_rhoR(i)=interpolation_order_2_hf( &
         rhoR(ibounce(1,i)+ishift,ibounce(2,i)+jshift,ibounce(3,i)+kshift), &
         rhoR(ibounce(1,i)+ishift2,ibounce(2,i)+jshift2,ibounce(3,i)+kshift2), &
         rhoR(ibounce(1,i)+ishift3,ibounce(2,i)+jshift3,ibounce(3,i)+kshift3))
      enddo
      if(idir==1 .or. idir==2)then
        do i=inits,ends
          bc_u(i)=bc_flow(i)/bc_rhoR(i)
        enddo
        where(bc_u(inits:ends)>cssq*HALF)
          bc_u(inits:ends)=cssq*HALF
        elsewhere(bc_u(inits:ends)<-cssq*HALF)
          bc_u(inits:ends)=-cssq*HALF
        end where
      endif
      if(idir==3 .or. idir==4)then
        do i=inits,ends
          bc_v(i)=bc_flow(i)/bc_rhoR(i)
        end do
        where(bc_v(inits:ends)>cssq*HALF)
          bc_v(inits:ends)=cssq*HALF
        elsewhere(bc_v(inits:ends)<-cssq*HALF)
          bc_v(inits:ends)=-cssq*HALF
        end where
      endif
      if(idir==5)then
        do i=inits,ends
          bc_w(i)=-bc_flow(i)/bc_rhoR(i)
        enddo
        where(bc_w(inits:ends)>cssq*HALF)
          bc_w(inits:ends)=cssq*HALF
        elsewhere(bc_w(inits:ends)<-cssq*HALF)
          bc_w(inits:ends)=-cssq*HALF
        end where
      endif
      if(idir==5 .or.idir==6)then
        do i=inits,ends
          bc_w(i)=bc_flow(i)/bc_rhoR(i)
        enddo
        where(bc_w(inits:ends)>cssq*HALF)
          bc_w(inits:ends)=cssq*HALF
        elsewhere(bc_w(inits:ends)<-cssq*HALF)
          bc_w(inits:ends)=-cssq*HALF
        end where
      endif
    else
      do i=inits,ends
        bc_rhoR(i)=interpolation_order_2_hf( &
         rhoR(ibounce(1,i)+ishift,ibounce(2,i)+jshift,ibounce(3,i)+kshift), &
         rhoR(ibounce(1,i)+ishift2,ibounce(2,i)+jshift2,ibounce(3,i)+kshift2), &
         rhoR(ibounce(1,i)+ishift3,ibounce(2,i)+jshift3,ibounce(3,i)+kshift3))
        bc_rhoB(i)=interpolation_order_2_hf( &
         rhoB(ibounce(1,i)+ishift,ibounce(2,i)+jshift,ibounce(3,i)+kshift), &
         rhoB(ibounce(1,i)+ishift2,ibounce(2,i)+jshift2,ibounce(3,i)+kshift2), &
         rhoB(ibounce(1,i)+ishift3,ibounce(2,i)+jshift3,ibounce(3,i)+kshift3))
      enddo
      if(idir==1 .or. idir==2)then
        do i=inits,ends
          bc_u(i)=bc_flow(i)/(bc_rhoR(i)+bc_rhoB(i))
        enddo
        where(bc_u(inits:ends)>cssq*HALF)
          bc_u(inits:ends)=cssq*HALF
        elsewhere(bc_u(inits:ends)<-cssq*HALF)
          bc_u(inits:ends)=-cssq*HALF
        end where
      endif
      if(idir==3 .or. idir==4)then
        do i=inits,ends
          bc_v(i)=bc_flow(i)/(bc_rhoR(i)+bc_rhoB(i))
        end do
        where(bc_v(inits:ends)>cssq*HALF)
          bc_v(inits:ends)=cssq*HALF
        elsewhere(bc_v(inits:ends)<-cssq*HALF)
          bc_v(inits:ends)=-cssq*HALF
        end where
      endif
      if(idir==5 .or.idir==6)then
        do i=inits,ends
          bc_w(i)=bc_flow(i)/(bc_rhoR(i)+bc_rhoB(i))
        enddo
        where(bc_w(inits:ends)>cssq*HALF)
          bc_w(inits:ends)=cssq*HALF
        elsewhere(bc_w(inits:ends)<-cssq*HALF)
          bc_w(inits:ends)=-cssq*HALF
        end where
      endif
    endif
  enddo
  
  return
  
 end subroutine set_bc_variable_hvar
 
 subroutine apply_bounceback_pop(rho_s,u_s,v_s,w_s,aoptp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the bounceback 
!     to fluid populations if necessary
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification September 2018
!     
!***********************************************************************
 
  implicit none 
  
  real(kind=PRC), allocatable, dimension(:)  :: rho_s,u_s,v_s,w_s
  
  type(REALPTR), dimension(0:links):: aoptp
  
  integer :: i,j,k,l,sx,ex,sz,ez
  
  real(kind=PRC), parameter :: pref_bouzidi=TWO/cssq
  real(kind=PRC), parameter :: cssq2 = ( HALF / cssq )
  real(kind=PRC), parameter :: cssq4 = ( HALF / (cssq*cssq) )
  integer, parameter :: kk = 1
  
  
  if(nbounce0>=1)then
  forall(i=1:nbounce0)
    buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     aoptp((2))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))
    aoptp((2))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     aoptp(1)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))
    aoptp(1)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i))
      
     buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
      aoptp((4))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))
     aoptp((4))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
      aoptp(3)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))
     aoptp(3)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
      buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i))
      
     buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
      aoptp((6))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))
     aoptp((6))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
      aoptp(5)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))
     aoptp(5)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
      buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i))
      
     buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
      aoptp((8))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))
     aoptp((8))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
      aoptp(7)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))
     aoptp(7)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
      buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i))
      
     buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
      aoptp((10))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))
     aoptp((10))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
      aoptp(9)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))
     aoptp(9)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
      buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i))
      
     buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
      aoptp((12))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))
     aoptp((12))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
      aoptp(11)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))
     aoptp(11)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
      buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i))
      
     buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
      aoptp((14))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))
     aoptp((14))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
      aoptp(13)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))
     aoptp(13)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
      buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i))
      
     buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
      aoptp((16))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))
     aoptp((16))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
      aoptp(15)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))
     aoptp(15)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
      buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i))
      
     buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
      aoptp((18))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))
     aoptp((18))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
      aoptp(17)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))
     aoptp(17)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
      buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i))
  end forall
  
  endif
  if(nbounce6>=nbounce0+1)then
  
  forall(i=nbounce0+1:nbounce6)
    !dirichlet condition
    !anti-bounce-back approach
    !from page 200 Kruger's book "the lattice boltzmann method"
    !NOTE de[x,y,z]=zero eliminated
    
    buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     aoptp((2))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))
        
    aoptp((2))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -aoptp(1)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))+ &
     p(2)*TWO*rho_s(i)* &
     (ONE+(dex(2)*u_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    aoptp(1)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i))+ &
     p(1)*TWO*rho_s(i)* &
     (ONE+(dex(1)*u_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     aoptp((4))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))
    
    aoptp((4))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -aoptp(3)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))+ &
     p(4)*TWO*rho_s(i)* &
     (ONE+(dey(4)*v_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
       
    aoptp(3)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i))+ &
     p(3)*TWO*rho_s(i)* &
     (ONE+(dey(3)*v_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     aoptp((6))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))
    
    aoptp((6))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -aoptp(5)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))+ &
     p(6)*TWO*rho_s(i)* &
     (ONE+(dez(6)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
       
    aoptp(5)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i))+ &
     p(5)*TWO*rho_s(i)* &
     (ONE+(dez(5)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     aoptp((8))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))
    
    aoptp((8))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -aoptp(7)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))+ &
     p(8)*TWO*rho_s(i)* &
     (ONE+(dex(8)*u_s(i)+ &
     dey(8)*v_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
       
    aoptp(7)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i))+ &
     p(7)*TWO*rho_s(i)* &
     (ONE+(dex(7)*u_s(i)+ &
     dey(7)*v_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     aoptp((10))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))
    
    aoptp((10))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -aoptp(9)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))+ &
     p(10)*TWO*rho_s(i)* &
     (ONE+(dex(10)*u_s(i)+ &
     dey(10)*v_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
       
    aoptp(9)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i))+ &
     p(9)*TWO*rho_s(i)* &
     (ONE+(dex(9)*u_s(i)+ &
     dey(9)*v_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     aoptp((12))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))
    
    aoptp((12))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -aoptp(11)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))+ &
     p(12)*TWO*rho_s(i)* &
     (ONE+(dex(12)*u_s(i)+ &
     dez(12)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
       
    aoptp(11)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i))+ &
     p(11)*TWO*rho_s(i)* &
     (ONE+(dex(11)*u_s(i)+ &
     dez(11)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     aoptp((14))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))
        
    aoptp((14))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -aoptp(13)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))+ &
     p(14)*TWO*rho_s(i)* &
     (ONE+(dex(14)*u_s(i)+ &
     dez(14)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
       
    aoptp(13)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i))+ &
     p(13)*TWO*rho_s(i)* &
     (ONE+(dex(13)*u_s(i)+ &
     dez(13)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     aoptp((16))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))
        
    aoptp((16))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -aoptp(15)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))+ &
     p(16)*TWO*rho_s(i)* &
     (ONE+(dey(16)*v_s(i)+ &
     dez(16)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
       
    aoptp(15)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i))+ &
     p(15)*TWO*rho_s(i)* &
     (ONE+(dey(15)*v_s(i)+ &
     dez(15)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     aoptp((18))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))
        
    aoptp((18))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -aoptp(17)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))+ &
     p(18)*TWO*rho_s(i)* &
     (ONE+(dey(18)*v_s(i)+ &
     dez(18)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
       
    aoptp(17)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i))+ &
     p(17)*TWO*rho_s(i)* &
     (ONE+(dey(17)*v_s(i)+ &
     dez(17)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
        
  end forall
  
  endif
  if(nbounce7>=nbounce6+1)then
  
  forall(i=nbounce6+1:nbounce7)
    !neumann condition
    !moving walls bounce-back approach
    !from page 180 Kruger's book "the lattice boltzmann method"
    !NOTE de[x,y,z]=zero eliminated
    buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     aoptp((2))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))
     
    aoptp((2))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     aoptp(1)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))- &
     p(2)*pref_bouzidi*rho_s(i)* &
     (dex(2)*u_s(i))
    
    aoptp(1)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i))- &
     p(1)*pref_bouzidi*rho_s(i)* &
     (dex(1)*u_s(i))
     
    buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     aoptp((4))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))
     
    aoptp((4))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     aoptp(3)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))- &
     p(4)*pref_bouzidi*rho_s(i)* &
     (dey(4)*v_s(i))
    
    aoptp(3)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i))- &
     p(3)*pref_bouzidi*rho_s(i)* &
     (dey(3)*v_s(i))
     
    buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     aoptp((6))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))
     
    aoptp((6))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     aoptp(5)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))- &
     p(6)*pref_bouzidi*rho_s(i)* &
     (dez(6)*w_s(i))
    
    aoptp(5)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i))- &
     p(5)*pref_bouzidi*rho_s(i)* &
     (dez(5)*w_s(i))
     
    buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     aoptp((8))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))
     
    aoptp((8))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     aoptp(7)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))- &
     p(8)*pref_bouzidi*rho_s(i)* &
     (dex(8)*u_s(i)+ &
     dey(8)*v_s(i))
    
    aoptp(7)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i))- &
     p(7)*pref_bouzidi*rho_s(i)* &
     (dex(7)*u_s(i)+ &
     dey(7)*v_s(i))
     
    buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     aoptp((10))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))
     
    aoptp((10))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     aoptp(9)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))- &
     p(10)*pref_bouzidi*rho_s(i)* &
     (dex(10)*u_s(i)+ &
     dey(10)*v_s(i))
    
    aoptp(9)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i))- &
     p(9)*pref_bouzidi*rho_s(i)* &
     (dex(9)*u_s(i)+ &
     dey(9)*v_s(i))
     
    buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     aoptp((12))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))
     
    aoptp((12))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     aoptp(11)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))- &
     p(12)*pref_bouzidi*rho_s(i)* &
     (dex(12)*u_s(i)+ &
     dez(12)*w_s(i))
    
    aoptp(11)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i))- &
     p(11)*pref_bouzidi*rho_s(i)* &
     (dex(11)*u_s(i)+ &
     dez(11)*w_s(i))
     
    buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     aoptp((14))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))
  
    aoptp((14))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     aoptp(13)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))- &
     p(14)*pref_bouzidi*rho_s(i)* &
     (dex(14)*u_s(i)+ &
     dez(14)*w_s(i))
    
    aoptp(13)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i))- &
     p(13)*pref_bouzidi*rho_s(i)* &
     (dex(13)*u_s(i)+ &
     dez(13)*w_s(i))
     
    buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     aoptp((16))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))
     
    aoptp((16))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     aoptp(15)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))- &
     p(16)*pref_bouzidi*rho_s(i)* &
     (dey(16)*v_s(i)+ &
     dez(16)*w_s(i))
    
    aoptp(15)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i))- &
     p(15)*pref_bouzidi*rho_s(i)* &
     (dey(15)*v_s(i)+ &
     dez(15)*w_s(i))
     
    buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     aoptp((18))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))
     
    aoptp((18))%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     aoptp(17)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))- &
     p(18)*pref_bouzidi*rho_s(i)* &
     (dey(18)*v_s(i)+ &
     dez(18)*w_s(i))
    
    aoptp(17)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     buffservice3d(ibounce(1,i),ibounce(2,i),ibounce(3,i))- &
     p(17)*pref_bouzidi*rho_s(i)* &
     (dey(17)*v_s(i)+ &
     dez(17)*w_s(i))
  
  end forall
  
  endif
  if(nbounce8>=nbounce7+1)then
  
  forall(i=nbounce7+1:nbounce8)
    !robin conditions
    !equilibrium scheme
    !from page 191 Kruger's book "the lattice boltzmann method"
    aoptp(0)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))= &
     equil_pop00(rho_s(i),u_s(i),v_s(i),w_s(i))

    aoptp(1)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))= &
     equil_pop01(rho_s(i),u_s(i),v_s(i),w_s(i))

    aoptp(2)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))= &
     equil_pop02(rho_s(i),u_s(i),v_s(i),w_s(i))

    aoptp(3)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))= &
     equil_pop03(rho_s(i),u_s(i),v_s(i),w_s(i))

    aoptp(4)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))= &
     equil_pop04(rho_s(i),u_s(i),v_s(i),w_s(i))

    aoptp(5)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))= &
     equil_pop05(rho_s(i),u_s(i),v_s(i),w_s(i))

    aoptp(6)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))= &
     equil_pop06(rho_s(i),u_s(i),v_s(i),w_s(i))

    aoptp(7)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))= &
     equil_pop07(rho_s(i),u_s(i),v_s(i),w_s(i))

    aoptp(8)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))= &
     equil_pop08(rho_s(i),u_s(i),v_s(i),w_s(i))

    aoptp(9)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))= &
     equil_pop09(rho_s(i),u_s(i),v_s(i),w_s(i))

    aoptp(10)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))= &
     equil_pop10(rho_s(i),u_s(i),v_s(i),w_s(i))

    aoptp(11)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))= &
     equil_pop11(rho_s(i),u_s(i),v_s(i),w_s(i))

    aoptp(12)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))= &
     equil_pop12(rho_s(i),u_s(i),v_s(i),w_s(i))

    aoptp(13)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))= &
     equil_pop13(rho_s(i),u_s(i),v_s(i),w_s(i))

    aoptp(14)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))= &
     equil_pop14(rho_s(i),u_s(i),v_s(i),w_s(i))

    aoptp(15)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))= &
     equil_pop15(rho_s(i),u_s(i),v_s(i),w_s(i))

    aoptp(16)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))= &
     equil_pop16(rho_s(i),u_s(i),v_s(i),w_s(i))

    aoptp(17)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))= &
     equil_pop17(rho_s(i),u_s(i),v_s(i),w_s(i))

    aoptp(18)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))= &
     equil_pop18(rho_s(i),u_s(i),v_s(i),w_s(i))
     
  end forall
  
  endif
  
  return
  
 end subroutine apply_bounceback_pop
 
 subroutine apply_bounceback_pop_halfway(nstep,rho_s,u_s,v_s,w_s,aoptp, &
  myalpha)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the bounceback 
!     to fluid populations if necessary in halfway mode,
!     page 82 of book: "the lattice boltzmann equation", S.Succi, 2001.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification September 2018
!     
!***********************************************************************
 
  implicit none 
  
  integer, intent(in) :: nstep
  real(kind=PRC), allocatable, dimension(:)  :: rho_s,u_s,v_s,w_s
  
  type(REALPTR), dimension(0:links):: aoptp
  real(kind=PRC), intent(in) :: myalpha
  
  integer :: i,j,k,l,sx,sz,iloop,indlow,indhig,ii,jj,kk,io,jo,ko,iii
  
  real(kind=PRC), dimension(0:links) :: myphi_CG
  real(kind=PRC), parameter :: pref_bouzidi=TWO/cssq
  real(kind=PRC), parameter :: cssq2 = ( HALF / cssq )
  real(kind=PRC), parameter :: cssq4 = ( HALF / (cssq*cssq) )
  real(kind=PRC) :: grad_rhox,grad_rhoy,grad_rhoz,locrhoR,locrhoB
  
  if(nbounce0>=1)then
#if 1
  do iii=1,nbounce0
  
    i=ibounce(1,iii)
    j=ibounce(2,iii)
    k=ibounce(3,iii)
  
    do iloop = 1, 9
  
      indlow = iloop*2 - 1
      indhig = iloop*2
  
      ii=i+ex(indlow)
      jj=j+ey(indlow)
      kk=k+ez(indlow)
  
      if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
       kk>=minz .and. kk<=maxz)then
        if(isfluid(ii,jj,kk)==1)then
	      aoptp(indlow)%p(i,j,k) = aoptp(indhig)%p(ii,jj,kk)
       endif
     endif
  
     ii=i+ex(indhig)
     jj=j+ey(indhig)
     kk=k+ez(indhig)
  
     if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
      kk>=minz .and. kk<=maxz)then
       if(isfluid(ii,jj,kk)==1)then
	     aoptp(indhig)%p(i,j,k) = aoptp(indlow)%p(ii,jj,kk)
	   endif
     endif
   enddo
 enddo
#else
    do i=1,nbounce0
    
    aoptp(1)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(2)%p(ibounce(1,i)+ex(1),ibounce(2,i),ibounce(3,i)), &
     kind=PRC)
    aoptp(2)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(1)%p(ibounce(1,i)+ex(2),ibounce(2,i),ibounce(3,i)), &
     kind=PRC)
    
    aoptp(3)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(4)%p(ibounce(1,i),ibounce(2,i)+ey(3),ibounce(3,i)), &
     kind=PRC)
    aoptp(4)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(3)%p(ibounce(1,i),ibounce(2,i)+ey(4),ibounce(3,i)), &
     kind=PRC)
    
    aoptp(5)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(6)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)+ez(5)), &
     kind=PRC)
    aoptp(6)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(5)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)+ez(6)), &
     kind=PRC)
    
    aoptp(7)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(8)%p(ibounce(1,i)+ex(7),ibounce(2,i)+ey(7),ibounce(3,i)), &
     kind=PRC)
    aoptp(8)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(7)%p(ibounce(1,i)+ex(8),ibounce(2,i)+ey(8),ibounce(3,i)), &
     kind=PRC)
    
    aoptp(9)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(10)%p(ibounce(1,i)+ex(9),ibounce(2,i)+ey(9),ibounce(3,i)), &
     kind=PRC)
    aoptp(10)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(9)%p(ibounce(1,i)+ex(10),ibounce(2,i)+ey(10),ibounce(3,i)), &
     kind=PRC)
     
    aoptp(11)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(12)%p(ibounce(1,i)+ex(11),ibounce(2,i),ibounce(3,i)+ez(11)), &
     kind=PRC)
    aoptp(12)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(11)%p(ibounce(1,i)+ex(12),ibounce(2,i),ibounce(3,i)+ez(12)), &
     kind=PRC)
    
    aoptp(13)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(14)%p(ibounce(1,i)+ex(13),ibounce(2,i),ibounce(3,i)+ez(13)), &
     kind=PRC)
    aoptp(14)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(13)%p(ibounce(1,i)+ex(14),ibounce(2,i),ibounce(3,i)+ez(14)), &
     kind=PRC)
    
    aoptp(15)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(16)%p(ibounce(1,i),ibounce(2,i)+ey(15),ibounce(3,i)+ez(15)), &
     kind=PRC)
    aoptp(16)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(15)%p(ibounce(1,i),ibounce(2,i)+ey(16),ibounce(3,i)+ez(16)), &
     kind=PRC)
     
    aoptp(17)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(18)%p(ibounce(1,i),ibounce(2,i)+ey(17),ibounce(3,i)+ez(17)), &
     kind=PRC)
    aoptp(18)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(17)%p(ibounce(1,i),ibounce(2,i)+ey(18),ibounce(3,i)+ez(18)), &
     kind=PRC)
    
    enddo
    
#endif
  
  endif
  
  if(nbounce6>=nbounce0+1)then
  if(lColourG)then
  do i=nbounce0+1,nbounce6
    !dirichlet condition
    !anti-bounce-back approach
    !from page 200 Kruger's book "the lattice boltzmann method"
    !NOTE de[x,y,z]=zero eliminated
    
    myphi_CG(0:links)=phi_CG(0:links)+varphi_CG(0:links)*myalpha
    
    aoptp(1)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -real(aoptp(2)%p(ibounce(1,i)+ex(1),ibounce(2,i),ibounce(3,i)), &
     kind=PRC)+myphi_CG(1)*TWO*rho_s(i) + &
     p(1)*TWO*rho_s(i)* &
     ((dex(1)*u_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    aoptp(2)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -real(aoptp(1)%p(ibounce(1,i)+ex(2),ibounce(2,i),ibounce(3,i)), &
     kind=PRC)+myphi_CG(2)*TWO*rho_s(i) +  &
     p(2)*TWO*rho_s(i)* &
     ((dex(2)*u_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
    
    
    aoptp(3)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -real(aoptp(4)%p(ibounce(1,i),ibounce(2,i)+ey(3),ibounce(3,i)), &
     kind=PRC)+myphi_CG(3)*TWO*rho_s(i) +  &
     p(3)*TWO*rho_s(i)* &
     ((dey(3)*v_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    aoptp(4)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -real(aoptp(3)%p(ibounce(1,i),ibounce(2,i)+ey(4),ibounce(3,i)), &
     kind=PRC)+myphi_CG(4)*TWO*rho_s(i) +  &
     p(4)*TWO*rho_s(i)* &
     ((dey(4)*v_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
    
    
    aoptp(5)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -real(aoptp(6)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)+ez(5)), &
     kind=PRC)+myphi_CG(5)*TWO*rho_s(i) +  &
     p(5)*TWO*rho_s(i)* &
     ((dez(5)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    aoptp(6)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -real(aoptp(5)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)+ez(6)), &
     kind=PRC)+myphi_CG(6)*TWO*rho_s(i) +  &
     p(6)*TWO*rho_s(i)* &
     ((dez(6)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
    
    
    aoptp(7)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -real(aoptp(8)%p(ibounce(1,i)+ex(7),ibounce(2,i)+ey(7),ibounce(3,i)), &
     kind=PRC)+myphi_CG(7)*TWO*rho_s(i) +  &
     p(7)*TWO*rho_s(i)* &
     ((dex(7)*u_s(i)+ &
     dey(7)*v_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    aoptp(8)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -real(aoptp(7)%p(ibounce(1,i)+ex(8),ibounce(2,i)+ey(8),ibounce(3,i)), &
     kind=PRC)+myphi_CG(8)*TWO*rho_s(i) +  &
     p(8)*TWO*rho_s(i)* &
     ((dex(8)*u_s(i)+ &
     dey(8)*v_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
    
    
    aoptp(9)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -real(aoptp(10)%p(ibounce(1,i)+ex(9),ibounce(2,i)+ey(9),ibounce(3,i)), &
     kind=PRC)+myphi_CG(9)*TWO*rho_s(i) +  &
     p(9)*TWO*rho_s(i)* &
     ((dex(9)*u_s(i)+ &
     dey(9)*v_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    aoptp(10)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -real(aoptp(9)%p(ibounce(1,i)+ex(10),ibounce(2,i)+ey(10),ibounce(3,i)), &
     kind=PRC)+myphi_CG(10)*TWO*rho_s(i) +  &
     p(10)*TWO*rho_s(i)* &
     ((dex(10)*u_s(i)+ &
     dey(10)*v_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
     
    aoptp(11)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -real(aoptp(12)%p(ibounce(1,i)+ex(11),ibounce(2,i),ibounce(3,i)+ez(11)), &
     kind=PRC)+myphi_CG(11)*TWO*rho_s(i) +  &
     p(11)*TWO*rho_s(i)* &
     ((dex(11)*u_s(i)+ &
     dez(11)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    aoptp(12)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -real(aoptp(11)%p(ibounce(1,i)+ex(12),ibounce(2,i),ibounce(3,i)+ez(12)), &
     kind=PRC)+myphi_CG(12)*TWO*rho_s(i) +  &
     p(12)*TWO*rho_s(i)* &
     ((dex(12)*u_s(i)+ &
     dez(12)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    
    aoptp(13)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -real(aoptp(14)%p(ibounce(1,i)+ex(13),ibounce(2,i),ibounce(3,i)+ez(13)), &
     kind=PRC)+myphi_CG(13)*TWO*rho_s(i) +  &
     p(13)*TWO*rho_s(i)* &
     ((dex(13)*u_s(i)+ &
     dez(13)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    aoptp(14)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -real(aoptp(13)%p(ibounce(1,i)+ex(14),ibounce(2,i),ibounce(3,i)+ez(14)), &
     kind=PRC)+myphi_CG(14)*TWO*rho_s(i) +  &
     p(14)*TWO*rho_s(i)* &
     ((dex(14)*u_s(i)+ &
     dez(14)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
    
    
    aoptp(15)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -real(aoptp(16)%p(ibounce(1,i),ibounce(2,i)+ey(15),ibounce(3,i)+ez(15)), &
     kind=PRC)+myphi_CG(15)*TWO*rho_s(i) +  &
     p(15)*TWO*rho_s(i)* &
     ((dey(15)*v_s(i)+ &
     dez(15)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    aoptp(16)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -real(aoptp(15)%p(ibounce(1,i),ibounce(2,i)+ey(16),ibounce(3,i)+ez(16)), &
     kind=PRC)+myphi_CG(16)*TWO*rho_s(i) +  &
     p(16)*TWO*rho_s(i)* &
     ((dey(16)*v_s(i)+ &
     dez(16)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
     
    aoptp(17)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -real(aoptp(18)%p(ibounce(1,i),ibounce(2,i)+ey(17),ibounce(3,i)+ez(17)), &
     kind=PRC)+myphi_CG(17)*TWO*rho_s(i) +  &
     p(17)*TWO*rho_s(i)* &
     ((dey(17)*v_s(i)+ &
     dez(17)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    aoptp(18)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -real(aoptp(17)%p(ibounce(1,i),ibounce(2,i)+ey(18),ibounce(3,i)+ez(18)), &
     kind=PRC)+myphi_CG(18)*TWO*rho_s(i) +  &
     p(18)*TWO*rho_s(i)* &
     ((dey(18)*v_s(i)+ &
     dez(18)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
    
  enddo
  else
  forall(i=nbounce0+1:nbounce6)
    !dirichlet condition
    !anti-bounce-back approach
    !from page 200 Kruger's book "the lattice boltzmann method"
    !NOTE de[x,y,z]=zero eliminated
    
    aoptp(1)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -real(aoptp(2)%p(ibounce(1,i)+ex(1),ibounce(2,i),ibounce(3,i)), &
     kind=PRC)+ &
     p(1)*TWO*rho_s(i)* &
     (ONE+(dex(1)*u_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    aoptp(2)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -real(aoptp(1)%p(ibounce(1,i)+ex(2),ibounce(2,i),ibounce(3,i)), &
     kind=PRC)+ &
     p(2)*TWO*rho_s(i)* &
     (ONE+(dex(2)*u_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
    
    
    aoptp(3)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -real(aoptp(4)%p(ibounce(1,i),ibounce(2,i)+ey(3),ibounce(3,i)), &
     kind=PRC)+ &
     p(3)*TWO*rho_s(i)* &
     (ONE+(dey(3)*v_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    aoptp(4)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -real(aoptp(3)%p(ibounce(1,i),ibounce(2,i)+ey(4),ibounce(3,i)), &
     kind=PRC)+ &
     p(4)*TWO*rho_s(i)* &
     (ONE+(dey(4)*v_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
    
    
    aoptp(5)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -real(aoptp(6)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)+ez(5)), &
     kind=PRC)+ &
     p(5)*TWO*rho_s(i)* &
     (ONE+(dez(5)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    aoptp(6)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -real(aoptp(5)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)+ez(6)), &
     kind=PRC)+ &
     p(6)*TWO*rho_s(i)* &
     (ONE+(dez(6)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
    
    
    aoptp(7)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -real(aoptp(8)%p(ibounce(1,i)+ex(7),ibounce(2,i)+ey(7),ibounce(3,i)), &
     kind=PRC)+ &
     p(7)*TWO*rho_s(i)* &
     (ONE+(dex(7)*u_s(i)+ &
     dey(7)*v_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    aoptp(8)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -real(aoptp(7)%p(ibounce(1,i)+ex(8),ibounce(2,i)+ey(8),ibounce(3,i)), &
     kind=PRC)+ &
     p(8)*TWO*rho_s(i)* &
     (ONE+(dex(8)*u_s(i)+ &
     dey(8)*v_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
    
    
    aoptp(9)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -real(aoptp(10)%p(ibounce(1,i)+ex(9),ibounce(2,i)+ey(9),ibounce(3,i)), &
     kind=PRC)+ &
     p(9)*TWO*rho_s(i)* &
     (ONE+(dex(9)*u_s(i)+ &
     dey(9)*v_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    aoptp(10)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -real(aoptp(9)%p(ibounce(1,i)+ex(10),ibounce(2,i)+ey(10),ibounce(3,i)), &
     kind=PRC)+ &
     p(10)*TWO*rho_s(i)* &
     (ONE+(dex(10)*u_s(i)+ &
     dey(10)*v_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
     
    aoptp(11)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -real(aoptp(12)%p(ibounce(1,i)+ex(11),ibounce(2,i),ibounce(3,i)+ez(11)), &
     kind=PRC)+ &
     p(11)*TWO*rho_s(i)* &
     (ONE+(dex(11)*u_s(i)+ &
     dez(11)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    aoptp(12)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -real(aoptp(11)%p(ibounce(1,i)+ex(12),ibounce(2,i),ibounce(3,i)+ez(12)), &
     kind=PRC)+ &
     p(12)*TWO*rho_s(i)* &
     (ONE+(dex(12)*u_s(i)+ &
     dez(12)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    
    aoptp(13)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -real(aoptp(14)%p(ibounce(1,i)+ex(13),ibounce(2,i),ibounce(3,i)+ez(13)), &
     kind=PRC)+ &
     p(13)*TWO*rho_s(i)* &
     (ONE+(dex(13)*u_s(i)+ &
     dez(13)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    aoptp(14)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -real(aoptp(13)%p(ibounce(1,i)+ex(14),ibounce(2,i),ibounce(3,i)+ez(14)), &
     kind=PRC)+ &
     p(14)*TWO*rho_s(i)* &
     (ONE+(dex(14)*u_s(i)+ &
     dez(14)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
    
    
    aoptp(15)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -real(aoptp(16)%p(ibounce(1,i),ibounce(2,i)+ey(15),ibounce(3,i)+ez(15)), &
     kind=PRC)+ &
     p(15)*TWO*rho_s(i)* &
     (ONE+(dey(15)*v_s(i)+ &
     dez(15)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    aoptp(16)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -real(aoptp(15)%p(ibounce(1,i),ibounce(2,i)+ey(16),ibounce(3,i)+ez(16)), &
     kind=PRC)+ &
     p(16)*TWO*rho_s(i)* &
     (ONE+(dey(16)*v_s(i)+ &
     dez(16)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
     
    aoptp(17)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -real(aoptp(18)%p(ibounce(1,i),ibounce(2,i)+ey(17),ibounce(3,i)+ez(17)), &
     kind=PRC)+ &
     p(17)*TWO*rho_s(i)* &
     (ONE+(dey(17)*v_s(i)+ &
     dez(17)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    aoptp(18)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     -real(aoptp(17)%p(ibounce(1,i),ibounce(2,i)+ey(18),ibounce(3,i)+ez(18)), &
     kind=PRC)+ &
     p(18)*TWO*rho_s(i)* &
     (ONE+(dey(18)*v_s(i)+ &
     dez(18)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
    
  end forall
  endif
  
  endif
  if(nbounce7>=nbounce6+1)then
  
  forall(i=nbounce6+1:nbounce7)
    !neumann condition
    !moving walls bounce-back approach
    !from page 180 Kruger's book "the lattice boltzmann method"
    !NOTE de[x,y,z]=zero eliminated
    
    aoptp(1)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(2)%p(ibounce(1,i)+ex(1),ibounce(2,i),ibounce(3,i)), &
     kind=PRC)- &
     p(2)*pref_bouzidi*rho_s(i)* &
     (dex(2)*u_s(i))
    
    aoptp(2)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(1)%p(ibounce(1,i)+ex(2),ibounce(2,i),ibounce(3,i)), &
     kind=PRC)- &
     p(1)*pref_bouzidi*rho_s(i)* &
     (dex(1)*u_s(i))
    
    
    aoptp(3)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(4)%p(ibounce(1,i),ibounce(2,i)+ey(3),ibounce(3,i)), &
     kind=PRC)- &
     p(4)*pref_bouzidi*rho_s(i)* &
     (dey(4)*v_s(i))
    
    aoptp(4)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(3)%p(ibounce(1,i),ibounce(2,i)+ey(4),ibounce(3,i)), &
     kind=PRC)- &
     p(3)*pref_bouzidi*rho_s(i)* &
     (dey(3)*v_s(i))
    
    
    aoptp(5)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(6)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)+ez(5)), &
     kind=PRC)- &
     p(6)*pref_bouzidi*rho_s(i)* &
     (dez(6)*w_s(i))
    
    aoptp(6)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(5)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)+ez(6)), &
     kind=PRC)- &
     p(5)*pref_bouzidi*rho_s(i)* &
     (dez(5)*w_s(i))
    
    
    aoptp(7)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(8)%p(ibounce(1,i)+ex(7),ibounce(2,i)+ey(7),ibounce(3,i)), &
     kind=PRC)- &
     p(8)*pref_bouzidi*rho_s(i)* &
     (dex(8)*u_s(i)+ &
     dey(8)*v_s(i))
    
    aoptp(8)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(7)%p(ibounce(1,i)+ex(8),ibounce(2,i)+ey(8),ibounce(3,i)), &
     kind=PRC)- &
     p(7)*pref_bouzidi*rho_s(i)* &
     (dex(7)*u_s(i)+ &
     dey(7)*v_s(i))
    
    
    aoptp(9)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(10)%p(ibounce(1,i)+ex(9),ibounce(2,i)+ey(9),ibounce(3,i)), &
     kind=PRC)- &
     p(10)*pref_bouzidi*rho_s(i)* &
     (dex(10)*u_s(i)+ &
     dey(10)*v_s(i))
    
    aoptp(10)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(9)%p(ibounce(1,i)+ex(10),ibounce(2,i)+ey(10),ibounce(3,i)), &
     kind=PRC)- &
     p(9)*pref_bouzidi*rho_s(i)* &
     (dex(9)*u_s(i)+ &
     dey(9)*v_s(i))
    
    
    aoptp(11)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(12)%p(ibounce(1,i)+ex(11),ibounce(2,i),ibounce(3,i)+ez(11)), &
     kind=PRC)- &
     p(12)*pref_bouzidi*rho_s(i)* &
     (dex(12)*u_s(i)+ &
     dez(12)*w_s(i))
    
    aoptp(12)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(11)%p(ibounce(1,i)+ex(12),ibounce(2,i),ibounce(3,i)+ez(12)), &
     kind=PRC)- &
     p(11)*pref_bouzidi*rho_s(i)* &
     (dex(11)*u_s(i)+ &
     dez(11)*w_s(i))
    
    
    aoptp(13)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(14)%p(ibounce(1,i)+ex(13),ibounce(2,i),ibounce(3,i)+ez(13)), &
     kind=PRC)- &
     p(14)*pref_bouzidi*rho_s(i)* &
     (dex(14)*u_s(i)+ &
     dez(14)*w_s(i))
    
    aoptp(14)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(13)%p(ibounce(1,i)+ex(14),ibounce(2,i),ibounce(3,i)+ez(14)), &
     kind=PRC)- &
     p(13)*pref_bouzidi*rho_s(i)* &
     (dex(13)*u_s(i)+ &
     dez(13)*w_s(i))
    
    
    aoptp(15)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(16)%p(ibounce(1,i),ibounce(2,i)+ey(15),ibounce(3,i)+ez(15)), &
     kind=PRC)- &
     p(16)*pref_bouzidi*rho_s(i)* &
     (dey(16)*v_s(i)+ &
     dez(16)*w_s(i))
    
    aoptp(16)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(15)%p(ibounce(1,i),ibounce(2,i)+ey(16),ibounce(3,i)+ez(16)), &
     kind=PRC)- &
     p(15)*pref_bouzidi*rho_s(i)* &
     (dey(15)*v_s(i)+ &
     dez(15)*w_s(i))
    
    
    aoptp(17)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(18)%p(ibounce(1,i),ibounce(2,i)+ey(17),ibounce(3,i)+ez(17)), &
     kind=PRC)- &
     p(18)*pref_bouzidi*rho_s(i)* &
     (dey(18)*v_s(i)+ &
     dez(18)*w_s(i))
    
    aoptp(18)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(17)%p(ibounce(1,i),ibounce(2,i)+ey(18),ibounce(3,i)+ez(18)), &
     kind=PRC)- &
     p(17)*pref_bouzidi*rho_s(i)* &
     (dey(17)*v_s(i)+ &
     dez(17)*w_s(i))
    
  end forall
  
  endif
  if(nbounce8>=nbounce7+1)then
  
  if(lColourG)then
    do ii=nbounce7+1,nbounce8
      grad_rhox = ZERO
      grad_rhoy = ZERO
      grad_rhoz = ZERO
      i=ibounce(1,ii)
      j=ibounce(2,ii)
      k=ibounce(3,ii)
      locrhoR = rhoR(i,j,k)
      locrhoB = rhoB(i,j,k)
      
      do l=0,links
        aoptp(l)%p(ibounce(1,ii),ibounce(2,ii),ibounce(3,ii))= &
         equil_popCG(l,bc_omega(ii),myalpha,rho_s(ii), &
         u_s(ii),v_s(ii),w_s(ii),grad_rhox,grad_rhoy,grad_rhoz)
      enddo
    enddo
  else
  
  forall(i=nbounce7+1:nbounce8)
    !robin conditions
    !equilibrium scheme
    !from page 191 Kruger's book "the lattice boltzmann method"
    aoptp(0)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))= &
     equil_pop00(rho_s(i),u_s(i),v_s(i),w_s(i))

    aoptp(1)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))= &
     equil_pop01(rho_s(i),u_s(i),v_s(i),w_s(i))

    aoptp(2)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))= &
     equil_pop02(rho_s(i),u_s(i),v_s(i),w_s(i))

    aoptp(3)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))= &
     equil_pop03(rho_s(i),u_s(i),v_s(i),w_s(i))

    aoptp(4)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))= &
     equil_pop04(rho_s(i),u_s(i),v_s(i),w_s(i))

    aoptp(5)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))= &
     equil_pop05(rho_s(i),u_s(i),v_s(i),w_s(i))

    aoptp(6)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))= &
     equil_pop06(rho_s(i),u_s(i),v_s(i),w_s(i))

    aoptp(7)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))= &
     equil_pop07(rho_s(i),u_s(i),v_s(i),w_s(i))

    aoptp(8)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))= &
     equil_pop08(rho_s(i),u_s(i),v_s(i),w_s(i))

    aoptp(9)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))= &
     equil_pop09(rho_s(i),u_s(i),v_s(i),w_s(i))

    aoptp(10)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))= &
     equil_pop10(rho_s(i),u_s(i),v_s(i),w_s(i))

    aoptp(11)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))= &
     equil_pop11(rho_s(i),u_s(i),v_s(i),w_s(i))

    aoptp(12)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))= &
     equil_pop12(rho_s(i),u_s(i),v_s(i),w_s(i))

    aoptp(13)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))= &
     equil_pop13(rho_s(i),u_s(i),v_s(i),w_s(i))

    aoptp(14)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))= &
     equil_pop14(rho_s(i),u_s(i),v_s(i),w_s(i))

    aoptp(15)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))= &
     equil_pop15(rho_s(i),u_s(i),v_s(i),w_s(i))

    aoptp(16)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))= &
     equil_pop16(rho_s(i),u_s(i),v_s(i),w_s(i))

    aoptp(17)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))= &
     equil_pop17(rho_s(i),u_s(i),v_s(i),w_s(i))

    aoptp(18)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i))= &
     equil_pop18(rho_s(i),u_s(i),v_s(i),w_s(i))
     
  end forall
  
  endif
  
  endif
  
  if(nbounce9>=nbounce8+1)then
  
  forall(i=nbounce8+1:nbounce9)
    !neumann condition
    !moving walls bounce-back approach
    !from page 180 Kruger's book "the lattice boltzmann method"
    !NOTE de[x,y,z]=zero eliminated
    
    aoptp(1)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(2)%p(ibounce(1,i)+ex(1),ibounce(2,i),ibounce(3,i)), &
     kind=PRC)- &
     p(2)*pref_bouzidi*rho_s(i)* &
     (dex(2)*u_s(i))
    
    aoptp(2)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(1)%p(ibounce(1,i)+ex(2),ibounce(2,i),ibounce(3,i)), &
     kind=PRC)- &
     p(1)*pref_bouzidi*rho_s(i)* &
     (dex(1)*u_s(i))
    
    
    aoptp(3)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(4)%p(ibounce(1,i),ibounce(2,i)+ey(3),ibounce(3,i)), &
     kind=PRC)- &
     p(4)*pref_bouzidi*rho_s(i)* &
     (dey(4)*v_s(i))
    
    aoptp(4)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(3)%p(ibounce(1,i),ibounce(2,i)+ey(4),ibounce(3,i)), &
     kind=PRC)- &
     p(3)*pref_bouzidi*rho_s(i)* &
     (dey(3)*v_s(i))
    
    
    aoptp(5)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(6)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)+ez(5)), &
     kind=PRC)- &
     p(6)*pref_bouzidi*rho_s(i)* &
     (dez(6)*w_s(i))
    
    aoptp(6)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(5)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)+ez(6)), &
     kind=PRC)- &
     p(5)*pref_bouzidi*rho_s(i)* &
     (dez(5)*w_s(i))
    
    
    aoptp(7)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(8)%p(ibounce(1,i)+ex(7),ibounce(2,i)+ey(7),ibounce(3,i)), &
     kind=PRC)- &
     p(8)*pref_bouzidi*rho_s(i)* &
     (dex(8)*u_s(i)+ &
     dey(8)*v_s(i))
    
    aoptp(8)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(7)%p(ibounce(1,i)+ex(8),ibounce(2,i)+ey(8),ibounce(3,i)), &
     kind=PRC)- &
     p(7)*pref_bouzidi*rho_s(i)* &
     (dex(7)*u_s(i)+ &
     dey(7)*v_s(i))
    
    
    aoptp(9)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(10)%p(ibounce(1,i)+ex(9),ibounce(2,i)+ey(9),ibounce(3,i)), &
     kind=PRC)- &
     p(10)*pref_bouzidi*rho_s(i)* &
     (dex(10)*u_s(i)+ &
     dey(10)*v_s(i))
    
    aoptp(10)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(9)%p(ibounce(1,i)+ex(10),ibounce(2,i)+ey(10),ibounce(3,i)), &
     kind=PRC)- &
     p(9)*pref_bouzidi*rho_s(i)* &
     (dex(9)*u_s(i)+ &
     dey(9)*v_s(i))
    
    
    aoptp(11)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(12)%p(ibounce(1,i)+ex(11),ibounce(2,i),ibounce(3,i)+ez(11)), &
     kind=PRC)- &
     p(12)*pref_bouzidi*rho_s(i)* &
     (dex(12)*u_s(i)+ &
     dez(12)*w_s(i))
    
    aoptp(12)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(11)%p(ibounce(1,i)+ex(12),ibounce(2,i),ibounce(3,i)+ez(12)), &
     kind=PRC)- &
     p(11)*pref_bouzidi*rho_s(i)* &
     (dex(11)*u_s(i)+ &
     dez(11)*w_s(i))
    
    
    aoptp(13)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(14)%p(ibounce(1,i)+ex(13),ibounce(2,i),ibounce(3,i)+ez(13)), &
     kind=PRC)- &
     p(14)*pref_bouzidi*rho_s(i)* &
     (dex(14)*u_s(i)+ &
     dez(14)*w_s(i))
    
    aoptp(14)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(13)%p(ibounce(1,i)+ex(14),ibounce(2,i),ibounce(3,i)+ez(14)), &
     kind=PRC)- &
     p(13)*pref_bouzidi*rho_s(i)* &
     (dex(13)*u_s(i)+ &
     dez(13)*w_s(i))
    
    
    aoptp(15)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(16)%p(ibounce(1,i),ibounce(2,i)+ey(15),ibounce(3,i)+ez(15)), &
     kind=PRC)- &
     p(16)*pref_bouzidi*rho_s(i)* &
     (dey(16)*v_s(i)+ &
     dez(16)*w_s(i))
    
    aoptp(16)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(15)%p(ibounce(1,i),ibounce(2,i)+ey(16),ibounce(3,i)+ez(16)), &
     kind=PRC)- &
     p(15)*pref_bouzidi*rho_s(i)* &
     (dey(15)*v_s(i)+ &
     dez(15)*w_s(i))
    
    
    aoptp(17)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(18)%p(ibounce(1,i),ibounce(2,i)+ey(17),ibounce(3,i)+ez(17)), &
     kind=PRC)- &
     p(18)*pref_bouzidi*rho_s(i)* &
     (dey(18)*v_s(i)+ &
     dez(18)*w_s(i))
    
    aoptp(18)%p(ibounce(1,i),ibounce(2,i),ibounce(3,i)) = &
     real(aoptp(17)%p(ibounce(1,i),ibounce(2,i)+ey(18),ibounce(3,i)+ez(18)), &
     kind=PRC)- &
     p(17)*pref_bouzidi*rho_s(i)* &
     (dey(17)*v_s(i)+ &
     dez(17)*w_s(i))
    
  end forall
  
  endif
  
  return
  
 end subroutine apply_bounceback_pop_halfway
 
 pure function interpolation_order_2(arg0,arg1,arg2)
 
!***********************************************************************
!     
!     LBsoft function for applying the Maclaurin 2° order using  
!     the forward finite difference coeff for the derivatives
!     at fullway grid point
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: arg0,arg1,arg2
  
  real(kind=PRC) :: interpolation_order_2
  
  !forward finite difference coeff 
  !(2° order first derivative)
  !(1° order second derivative)
  real(kind=PRC), parameter, dimension(0:2) :: coeff1=(/-THREE*HALF,TWO,-HALF/)
  real(kind=PRC), parameter, dimension(0:2) :: coeff2=(/ONE,-TWO,ONE/)*HALF
  
  !derivates
!  der1=coeff1(0)*arg0+coeff1(1)*arg1+coeff1(2)*arg2
!  der2=coeff2(0)*arg0+coeff2(1)*arg1+coeff2(2)*arg2
  
  !apply Maclaurin 2° order
  interpolation_order_2=arg0- &
   (coeff1(0)*arg0+coeff1(1)*arg1+coeff1(2)*arg2)+ &
   (coeff2(0)*arg0+coeff2(1)*arg1+coeff2(2)*arg2)
  
  return
  
 end function interpolation_order_2
 
 pure function interpolation_order_2_hf(arg0,arg1,arg2)
 
!***********************************************************************
!     
!     LBsoft function for applying the Maclaurin 2° order using  
!     the forward finite difference coeff for the derivatives
!     at halfway grid point
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: arg0,arg1,arg2
  
  real(kind=PRC) :: interpolation_order_2_hf
  
  !forward finite difference coeff 
  !(2° order first derivative)
  !(1° order second derivative)
  real(kind=PRC), parameter, dimension(0:2) :: coeff1=(/-THREE*HALF,TWO,-HALF/)*HALF
  real(kind=PRC), parameter, dimension(0:2) :: coeff2=(/ONE,-TWO,ONE/)*ONE/EIGHT
  
  !derivates
!  der1=coeff1(0)*arg0+coeff1(1)*arg1+coeff1(2)*arg2
!  der2=coeff2(0)*arg0+coeff2(1)*arg1+coeff2(2)*arg2
  
  !apply Maclaurin 2° order
  interpolation_order_2_hf=arg0- &
   (coeff1(0)*arg0+coeff1(1)*arg1+coeff1(2)*arg2)+ &
   (coeff2(0)*arg0+coeff2(1)*arg1+coeff2(2)*arg2)
  
  return
  
 end function interpolation_order_2_hf
 
 subroutine bounceback_pop(aoptp)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying the bounceback 
!     to fluid populations if necessary
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Bernaschi
!     last modification September 2018
!     
!***********************************************************************
 
  implicit none 
  
  type(REALPTR), dimension(0:links):: aoptp
 
  integer :: i,j,k,l, ishift, jshift, kshift
  integer :: itemp, jtemp, ktemp
  integer, save :: iter=0, lminx, lminy, lminz, lmaxx, lmaxy, lmaxz
  
  
  if(iter.eq.0) then
     lmaxx=maxx+1
     lmaxy=maxy+1
     lmaxz=maxz+1
     lminx=minx-1
     lminy=miny-1
     lminz=minz-1
     if((ixpbc.ne.1.AND.iypbc.ne.1).OR.(ixpbc.ne.1.AND.izpbc.ne.1)) then
        lmaxx=maxx
        lminx=minx
     endif
    if((iypbc.ne.1.AND.izpbc.ne.1)) then
        lmaxy=maxy
        lminy=miny
    endif
  endif
  iter=iter+1

  do l=1,links
    ishift=ex(l)
    jshift=ey(l)
    kshift=ez(l)
    if(ixpbc.ne.1) then
      do i=wminx,wmaxx,wmaxx-wminx
          if(mod(l,2).eq.1) then
            forall(j=miny-1:maxy+1,k=minz-1:maxz+1) buffservice3d(i,j,k) = aoptp((l+1))%p(i,j,k)
            forall(j=miny-1:maxy+1,k=minz-1:maxz+1) aoptp((l+1))%p(i,j,k) = aoptp(l)%p(i,j,k)
            forall(j=miny-1:maxy+1,k=minz-1:maxz+1) aoptp(l)%p(i,j,k) = buffservice3d(i,j,k)
          endif
      enddo
    endif
    if(iypbc.ne.1) then      
      do j=wminy,wmaxy,wmaxy-wminy
          if(mod(l,2).eq.1) then
            forall(i=lminx:lmaxx,k=minz-1:maxz+1) buffservice3d(i,j,k) = aoptp((l+1))%p(i,j,k)
            forall(i=lminx:lmaxx,k=minz-1:maxz+1) aoptp((l+1))%p(i,j,k) = aoptp(l)%p(i,j,k)
            forall(i=lminx:lmaxx,k=minz-1:maxz+1) aoptp(l)%p(i,j,k) = buffservice3d(i,j,k)
          endif
      enddo
    endif
    if(izpbc.ne.1) then
       do k=wminz,wmaxz,wmaxz-wminz
          if(mod(l,2).eq.1) then
            forall(i=lminx:lmaxx,j=lminy:lmaxy) buffservice3d(i,j,k) = aoptp((l+1))%p(i,j,k)
            forall(i=lminx:lmaxx,j=lminy:lmaxy) aoptp((l+1))%p(i,j,k) = aoptp(l)%p(i,j,k)
            forall(i=lminx:lmaxx,j=lminy:lmaxy) aoptp(l)%p(i,j,k) = buffservice3d(i,j,k)
          endif
       enddo
    endif
 enddo
 
 return
  
 end subroutine bounceback_pop
 
!******************END PART TO MANAGE THE BOUNCEBACK********************

!******************START PART TO MANAGE THE COPY WALL*******************

 subroutine compute_densities_wall(nstep)
 
!***********************************************************************
!     
!     LBsoft subroutine for driving the reflection of the density 
!     variables if requested from the boundary conditions
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  integer, intent(in) :: nstep
  integer :: i,j,k,l,ishift,jshift,kshift, errorCount
  integer :: ishift2,jshift2,kshift2
  integer :: ishift3,jshift3,kshift3
  integer :: idir,inits,ends
  REAL(kind=PRC) :: dsum1,dsum2
  REAL(kind=PRC) :: isum
  logical :: ltest(1)=.false.
  
  
   errorCount = 0
   
   select case(iselwetting)
   case (1)
     if(lsingle_fluid)then
       do k=minz-1,maxz+1
         do j=miny-1,maxy+1
           do i=minx-1,maxx+1
             if(isfluid(i,j,k)==0 .or. isfluid(i,j,k)==2)then
               rhoR(i,j,k)=densR_wetting
             endif
           enddo
         enddo
       enddo
     else
       do k=minz-1,maxz+1
         do j=miny-1,maxy+1
           do i=minx-1,maxx+1
             if(isfluid(i,j,k)==0 .or. isfluid(i,j,k)==2)then
               rhoR(i,j,k)=densR_wetting
               rhoB(i,j,k)=densB_wetting
             endif
           enddo
         enddo
       enddo
     endif
   case default
     if(lsingle_fluid)then
       do k=minz-1,maxz+1
         do j=miny-1,maxy+1
           do i=minx-1,maxx+1
             if(isfluid(i,j,k)==0 .or. isfluid(i,j,k)==2)then
               dsum1=ZERO
               isum=ZERO
               do l=1,linksd3q27
                 ishift=i+exd3q27(l)
                 jshift=j+eyd3q27(l)
                 kshift=k+ezd3q27(l)
                 if(isfluid(ishift,jshift,kshift)==1)then
                   dsum1=dsum1+pd3q27(l)*rhoR(ishift,jshift,kshift)
                   isum=isum+pd3q27(l)
                 endif
               enddo
               if(isum==ZERO .and. isfluid(i,j,k)==2)then
                 dsum1=ZERO
                 isum=ZERO
                 do l=1,ndouble
                   ishift=i+exdouble(l)
                   jshift=j+eydouble(l)
                   kshift=k+ezdouble(l)
                   if(ishift<minx-nbuff)cycle
                   if(ishift>maxx+nbuff)cycle
                   if(jshift<miny-nbuff)cycle
                   if(jshift>maxy+nbuff)cycle
                   if(kshift<minz-nbuff)cycle
                   if(kshift>maxz+nbuff)cycle
                   if(isfluid(ishift,jshift,kshift)==1)then
                     dsum1=dsum1+rhoR(ishift,jshift,kshift)
                     isum=isum+ONE
                   endif
                 enddo
                 if(isum==ZERO)then
                   ltest(1)=.true.
                   errorCount = errorCount + 1
                   rhoR(i,j,k)= MINDENS
                 else
                   rhoR(i,j,k)=dsum1/isum
                 endif
               else
                 rhoR(i,j,k)=dsum1/isum
               endif
             endif
           enddo
         enddo
       enddo
     else
       do k=minz-1,maxz+1
         do j=miny-1,maxy+1
           do i=minx-1,maxx+1
             if(isfluid(i,j,k)==0 .or. isfluid(i,j,k)==2)then
               dsum1=ZERO
               dsum2=ZERO
               isum=ZERO
               do l=1,linksd3q27
                 ishift=i+exd3q27(l)
                 jshift=j+eyd3q27(l)
                 kshift=k+ezd3q27(l)
                 if(isfluid(ishift,jshift,kshift)==1)then
                   dsum1=dsum1+pd3q27(l)*rhoR(ishift,jshift,kshift)
                  dsum2=dsum2+pd3q27(l)*rhoB(ishift,jshift,kshift)
                  isum=isum+pd3q27(l)
                 endif
               enddo
               if(isum==ZERO .and. isfluid(i,j,k)==2)then
                 dsum1=ZERO
                 dsum2=ZERO
                 isum=ZERO
                 do l=1,ndouble
                   ishift=i+exdouble(l)
                   jshift=j+eydouble(l)
                   kshift=k+ezdouble(l)
                   if(ishift<minx-nbuff)cycle
                   if(ishift>maxx+nbuff)cycle
                   if(jshift<miny-nbuff)cycle
                   if(jshift>maxy+nbuff)cycle
                   if(kshift<minz-nbuff)cycle
                   if(kshift>maxz+nbuff)cycle
                   if(isfluid(ishift,jshift,kshift)==1)then
                     dsum1=dsum1+rhoR(ishift,jshift,kshift)
                     dsum2=dsum2+rhoB(ishift,jshift,kshift)
                     isum=isum+ONE
                   endif
                 enddo
                 if(isum==ZERO)then
                   ltest(1)=.true.
                   errorCount = errorCount + 1
                   rhoR(i,j,k)= MINDENS
                   rhoB(i,j,k)= MINDENS
                 else
                   if(isum>MINDENS)then
                     rhoR(i,j,k)=dsum1/isum
                     rhoB(i,j,k)=dsum2/isum
                   else
                     rhoR(i,j,k)= MINDENS
                     rhoB(i,j,k)= MINDENS
                   endif
                 endif
               else
                 if(isum>MINDENS)then
                   rhoR(i,j,k)=dsum1/isum
                   rhoB(i,j,k)=dsum2/isum
                 else
                   rhoR(i,j,k)= MINDENS
                   rhoB(i,j,k)= MINDENS
                 endif
               endif
             endif
           enddo
         enddo
       enddo
     endif
   end select
   
   do idir=1,nbcdir
    inits=nbounce6dir(idir-1)+1
    ends=nbounce6dir(idir)
    if(inits>ends)cycle
    ishift=ex(idir)
    jshift=ey(idir)
    kshift=ez(idir)
    ishift2=ex(idir)*2
    jshift2=ey(idir)*2
    kshift2=ez(idir)*2
    ishift3=ex(idir)*3
    jshift3=ey(idir)*3
    kshift3=ez(idir)*3
    if(lsingle_fluid)then
      do i=inits,ends
        rhoR(ibounce(1,i),ibounce(2,i),ibounce(3,i))=interpolation_order_2_hf( &
         rhoR(ibounce(1,i)+ishift,ibounce(2,i)+jshift,ibounce(3,i)+kshift), &
         rhoR(ibounce(1,i)+ishift2,ibounce(2,i)+jshift2,ibounce(3,i)+kshift2), &
         rhoR(ibounce(1,i)+ishift3,ibounce(2,i)+jshift3,ibounce(3,i)+kshift3))
      enddo
    else
      do i=inits,ends
        rhoR(ibounce(1,i),ibounce(2,i),ibounce(3,i))=interpolation_order_2_hf( &
         rhoR(ibounce(1,i)+ishift,ibounce(2,i)+jshift,ibounce(3,i)+kshift), &
         rhoR(ibounce(1,i)+ishift2,ibounce(2,i)+jshift2,ibounce(3,i)+kshift2), &
         rhoR(ibounce(1,i)+ishift3,ibounce(2,i)+jshift3,ibounce(3,i)+kshift3))
        rhoB(ibounce(1,i),ibounce(2,i),ibounce(3,i))=interpolation_order_2_hf( &
         rhoB(ibounce(1,i)+ishift,ibounce(2,i)+jshift,ibounce(3,i)+kshift), &
         rhoB(ibounce(1,i)+ishift2,ibounce(2,i)+jshift2,ibounce(3,i)+kshift2), &
         rhoB(ibounce(1,i)+ishift3,ibounce(2,i)+jshift3,ibounce(3,i)+kshift3))
      enddo
    endif
  enddo
  
  
  do idir=1,nbcdir
    inits=nbounce7dir(idir-1)+1
    ends=nbounce7dir(idir)
    if(inits>ends)cycle
    ishift=ex(idir)
    jshift=ey(idir)
    kshift=ez(idir)
    ishift2=ex(idir)*2
    jshift2=ey(idir)*2
    kshift2=ez(idir)*2
    ishift3=ex(idir)*3
    jshift3=ey(idir)*3
    kshift3=ez(idir)*3
    if(lsingle_fluid)then
      do i=inits,ends
        rhoR(ibounce(1,i),ibounce(2,i),ibounce(3,i))=interpolation_order_2_hf( &
         rhoR(ibounce(1,i)+ishift,ibounce(2,i)+jshift,ibounce(3,i)+kshift), &
         rhoR(ibounce(1,i)+ishift2,ibounce(2,i)+jshift2,ibounce(3,i)+kshift2), &
         rhoR(ibounce(1,i)+ishift3,ibounce(2,i)+jshift3,ibounce(3,i)+kshift3))
      enddo
    else
      do i=inits,ends
        rhoR(ibounce(1,i),ibounce(2,i),ibounce(3,i))=interpolation_order_2_hf( &
         rhoR(ibounce(1,i)+ishift,ibounce(2,i)+jshift,ibounce(3,i)+kshift), &
         rhoR(ibounce(1,i)+ishift2,ibounce(2,i)+jshift2,ibounce(3,i)+kshift2), &
         rhoR(ibounce(1,i)+ishift3,ibounce(2,i)+jshift3,ibounce(3,i)+kshift3))
        rhoB(ibounce(1,i),ibounce(2,i),ibounce(3,i))=interpolation_order_2_hf( &
         rhoB(ibounce(1,i)+ishift,ibounce(2,i)+jshift,ibounce(3,i)+kshift), &
         rhoB(ibounce(1,i)+ishift2,ibounce(2,i)+jshift2,ibounce(3,i)+kshift2), &
         rhoB(ibounce(1,i)+ishift3,ibounce(2,i)+jshift3,ibounce(3,i)+kshift3))
      enddo
    endif
  enddo
  
  do idir=1,nbcdir
    inits=nbounce8dir(idir-1)+1
    ends=nbounce8dir(idir)
    if(inits>ends)cycle
    ishift=ex(idir)
    jshift=ey(idir)
    kshift=ez(idir)
    ishift2=ex(idir)*2
    jshift2=ey(idir)*2
    kshift2=ez(idir)*2
    ishift3=ex(idir)*3
    jshift3=ey(idir)*3
    kshift3=ez(idir)*3
    if(lsingle_fluid)then
      do i=inits,ends
        rhoR(ibounce(1,i),ibounce(2,i),ibounce(3,i))=interpolation_order_2_hf( &
         rhoR(ibounce(1,i)+ishift,ibounce(2,i)+jshift,ibounce(3,i)+kshift), &
         rhoR(ibounce(1,i)+ishift2,ibounce(2,i)+jshift2,ibounce(3,i)+kshift2), &
         rhoR(ibounce(1,i)+ishift3,ibounce(2,i)+jshift3,ibounce(3,i)+kshift3))
      enddo
    else
      do i=inits,ends
        rhoR(ibounce(1,i),ibounce(2,i),ibounce(3,i))=interpolation_order_2_hf( &
         rhoR(ibounce(1,i)+ishift,ibounce(2,i)+jshift,ibounce(3,i)+kshift), &
         rhoR(ibounce(1,i)+ishift2,ibounce(2,i)+jshift2,ibounce(3,i)+kshift2), &
         rhoR(ibounce(1,i)+ishift3,ibounce(2,i)+jshift3,ibounce(3,i)+kshift3))
        rhoB(ibounce(1,i),ibounce(2,i),ibounce(3,i))=interpolation_order_2_hf( &
         rhoB(ibounce(1,i)+ishift,ibounce(2,i)+jshift,ibounce(3,i)+kshift), &
         rhoB(ibounce(1,i)+ishift2,ibounce(2,i)+jshift2,ibounce(3,i)+kshift2), &
         rhoB(ibounce(1,i)+ishift3,ibounce(2,i)+jshift3,ibounce(3,i)+kshift3))
      enddo
    endif
  enddo
  
  do idir=1,nbcdir
    inits=nbounce9dir(idir-1)+1
    ends=nbounce9dir(idir)
    if(inits>ends)cycle
    ishift=ex(idir)
    jshift=ey(idir)
    kshift=ez(idir)
    ishift2=ex(idir)*2
    jshift2=ey(idir)*2
    kshift2=ez(idir)*2
    ishift3=ex(idir)*3
    jshift3=ey(idir)*3
    kshift3=ez(idir)*3
    if(lsingle_fluid)then
      do i=inits,ends
        bc_rhoR(i)=interpolation_order_2_hf( &
         rhoR(ibounce(1,i)+ishift,ibounce(2,i)+jshift,ibounce(3,i)+kshift), &
         rhoR(ibounce(1,i)+ishift2,ibounce(2,i)+jshift2,ibounce(3,i)+kshift2), &
         rhoR(ibounce(1,i)+ishift3,ibounce(2,i)+jshift3,ibounce(3,i)+kshift3))
      enddo
    else
      do i=inits,ends
        bc_rhoR(i)=interpolation_order_2_hf( &
         rhoR(ibounce(1,i)+ishift,ibounce(2,i)+jshift,ibounce(3,i)+kshift), &
         rhoR(ibounce(1,i)+ishift2,ibounce(2,i)+jshift2,ibounce(3,i)+kshift2), &
         rhoR(ibounce(1,i)+ishift3,ibounce(2,i)+jshift3,ibounce(3,i)+kshift3))
        bc_rhoB(i)=interpolation_order_2_hf( &
         rhoB(ibounce(1,i)+ishift,ibounce(2,i)+jshift,ibounce(3,i)+kshift), &
         rhoB(ibounce(1,i)+ishift2,ibounce(2,i)+jshift2,ibounce(3,i)+kshift2), &
         rhoB(ibounce(1,i)+ishift3,ibounce(2,i)+jshift3,ibounce(3,i)+kshift3))
      enddo
    endif
  enddo
   
   
   !call or_world_larr(ltest,1)
   if(ltest(1)) write(216,fmt=1003) nstep, idrank, errorCount
1003 format ("step:", I6, " id:", I4,I10," fix in compute_densities_wall", X)
  
 end subroutine compute_densities_wall
 
!*******************END PART TO MANAGE THE COPY WALL********************

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
  
#ifdef SCANDREA
  !red fluid
  !forall(i=minx:maxx,j=miny:maxy,k=minz:maxz,isfluid(i,j,k)==1)
  do k=minz,maxz
   do j=miny,maxy
    do i=minx,maxx

        if (isfluid(i,j,k)/=1) cycle
    fuR(i,j,k) = fuR(i,j,k) - pair_SC*gradpsixB(i,j,k)
    fvR(i,j,k) = fvR(i,j,k) - pair_SC*gradpsiyB(i,j,k)
    fwR(i,j,k) = fwR(i,j,k) - pair_SC*gradpsizB(i,j,k)
      enddo
   enddo
  enddo
  !end forall
  !blue fluid
  !forall(i=minx:maxx,j=miny:maxy,k=minz:maxz,isfluid(i,j,k)==1)
  do k=minz,maxz
   do j=miny,maxy
    do i=minx,maxx

        if (isfluid(i,j,k)/=1) cycle
    fuB(i,j,k) = fuB(i,j,k) - pair_SC*gradpsixR(i,j,k)
    fvB(i,j,k) = fvB(i,j,k) - pair_SC*gradpsiyR(i,j,k)
    fwB(i,j,k) = fwB(i,j,k) - pair_SC*gradpsizR(i,j,k)
    enddo
   enddo
  enddo
  !end forall
#else
  !red fluid
  !forall(i=minx:maxx,j=miny:maxy,k=minz:maxz,isfluid(i,j,k)==1)
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
  !end forall
  !blue fluid
  !forall(i=minx:maxx,j=miny:maxy,k=minz:maxz,isfluid(i,j,k)==1)
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
  !end forall
#endif

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


 subroutine init_particle_2_isfluid(isub,jsub,ksub,nspheres, &
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
  
  integer, intent(in) :: isub,jsub,ksub,nspheres,nspheredeads
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
           otemp,vx,vy,vz,fx,fy,fz,tx,ty,tz,rhoR,aoptpR, debug,iatm, A, l)
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
        call particle_to_node_bounce_back_bc2_phase1(lrotate,nstep,i,j,k,aoptpR, debug,iatm)
      endif

      !the fluid bounce back is local so I have to do it
      if(i/=ii.or.j/=jj.or.k/=kk)then
       if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. kk>=minz .and. kk<=maxz)then
          call particle_to_node_bounce_back_bc2_phase1(lrotate,nstep,ii,jj,kk,aoptpR, debug,iatm)
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
           otemp,vx,vy,vz,fx,fy,fz,tx,ty,tz,rhoR,aoptpR, debug,iatm, A, l)
        call node_to_particle_bounce_back_bc2(lrotate,nstep,i,j,k,rtemp, &
           otemp,vx,vy,vz,fx,fy,fz,tx,ty,tz,rhoB,aoptpB, debug,iatm, A, l)

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
         aoptpR, debug,iatm)
        call particle_to_node_bounce_back_bc2_phase1(lrotate,nstep,i,j,k, &
         aoptpB, debug,iatm)
      endif

      !the fluid bounce back is local so I have to do it
      if(i/=ii.or.j/=jj.or.k/=kk)then
       if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. kk>=minz .and. kk<=maxz)then
          call particle_to_node_bounce_back_bc2_phase1(lrotate,nstep,ii,jj,kk, &
           aoptpR, debug,iatm)
          call particle_to_node_bounce_back_bc2_phase1(lrotate,nstep,ii,jj,kk, &
           aoptpB, debug,iatm)
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
         otemp,vx,vy,vz,rhoR,aoptpR, debug,iatm)
      endif

      !the fluid bounce back is local so I have to do it
      if(i/=ii.or.j/=jj.or.k/=kk)then
       if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. kk>=minz .and. kk<=maxz)then
          call particle_to_node_bounce_back_bc2_phase2(lrotate,nstep,ii,jj,kk,rtemp, &
           otemp,vx,vy,vz,rhoR,aoptpR, debug,iatm)
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
         otemp,vx,vy,vz,rhoR,aoptpR, debug,iatm)
        call particle_to_node_bounce_back_bc2_phase2(lrotate,nstep,i,j,k,rtemp, &
         otemp,vx,vy,vz,rhoB,aoptpB, debug,iatm)
      endif

      !the fluid bounce back is local so I have to do it
      if(i/=ii.or.j/=jj.or.k/=kk)then
       if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. kk>=minz .and. kk<=maxz)then
          call particle_to_node_bounce_back_bc2_phase2(lrotate,nstep,ii,jj,kk,rtemp, &
           otemp,vx,vy,vz,rhoR,aoptpR, debug,iatm)
          call particle_to_node_bounce_back_bc2_phase2(lrotate,nstep,ii,jj,kk,rtemp, &
           otemp,vx,vy,vz,rhoB,aoptpB, debug,iatm)
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
   otemp,vxs,vys,vzs,fx,fy,fz,tx,ty,tz,rhosub,aoptp, debug,iatm, A, l)
 
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
  type(REALPTR), dimension(0:links):: aoptp
  
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

      f2p = TWO * aoptp(indhig)%p(ii,jj,kk)- &
       p(indhig)*pref_bouzidi*rhosub(ii,jj,kk)*(dex(indhig)*vx+dey(indhig)*vy+dez(indhig)*vz)

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
      !       "pop",aoptp(indhig)%p(ii,jj,kk), "rho",rhosub(ii,jj,kk), &
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

      f2p = TWO * aoptp(indlow)%p(ii,jj,kk) - &
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
      !       "pop",aoptp(indlow)%p(ii,jj,kk), "rho",rhosub(ii,jj,kk), &
      !       "v", vx,vy,vz, vxs,vys,vzs
      ! endif
    endif
  endif
  endif
 enddo
  
 end subroutine node_to_particle_bounce_back_bc2
 
 subroutine node_to_particle_bounce_back_bc(lrotate,nstep,i,j,k,rversor, &
   otemp,vxs,vys,vzs,fx,fy,fz,tx,ty,tz,rhosub,aoptp)
 
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
  type(REALPTR), dimension(0:links):: aoptp
  
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
      f2p=TWO*real(aoptp(2)%p(ii,jj,kk),kind=PRC)- &
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
      f2p=TWO*real(aoptp(1)%p(ii,jj,kk),kind=PRC)- &
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
      f2p=TWO*real(aoptp(4)%p(ii,jj,kk),kind=PRC)- &
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
      f2p=TWO*real(aoptp(3)%p(ii,jj,kk),kind=PRC)- &
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
      f2p=TWO*real(aoptp(6)%p(ii,jj,kk),kind=PRC)- &
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
      f2p=TWO*real(aoptp(5)%p(ii,jj,kk),kind=PRC)- &
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
      f2p=TWO*real(aoptp(8)%p(ii,jj,kk),kind=PRC)- &
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
      f2p=TWO*real(aoptp(7)%p(ii,jj,kk),kind=PRC)- &
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
      f2p=TWO*real(aoptp(10)%p(ii,jj,kk),kind=PRC)- &
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
      f2p=TWO*real(aoptp(9)%p(ii,jj,kk),kind=PRC)- &
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
      f2p=TWO*real(aoptp(12)%p(ii,jj,kk),kind=PRC)- &
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
      f2p=TWO*real(aoptp(11)%p(ii,jj,kk),kind=PRC)- &
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
      f2p=TWO*real(aoptp(14)%p(ii,jj,kk),kind=PRC)- &
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
      f2p=TWO*real(aoptp(13)%p(ii,jj,kk),kind=PRC)- &
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
      f2p=TWO*real(aoptp(16)%p(ii,jj,kk),kind=PRC)- &
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
      f2p=TWO*real(aoptp(15)%p(ii,jj,kk),kind=PRC)- &
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
      f2p=TWO*real(aoptp(18)%p(ii,jj,kk),kind=PRC)- &
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
      f2p=TWO*real(aoptp(17)%p(ii,jj,kk),kind=PRC)- &
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

 subroutine particle_to_node_bounce_back_bc2_phase1(lrotate,nstep,i,j,k, aoptp, debug,iatm)
 
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
  type(REALPTR), dimension(0:links):: aoptp
  
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
	      aoptp(indlow)%p(i,j,k) = aoptp(indhig)%p(ii,jj,kk)
           ! if (debug) then
           !   write (iounit(3),*) __FILE__,__LINE__, ii,jj,kk, "pop=", indlow,indhig, &
           !       "pop now",aoptp(indhig)%p(ii,jj,kk), "@",i,j,k
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
	      aoptp(indhig)%p(i,j,k) = aoptp(indlow)%p(ii,jj,kk)
	      ! if (debug) then
              !   write (iounit(3),*) __FILE__,__LINE__, ii,jj,kk, "pop=", indhig,indlow, &
              !    "pop now",aoptp(indlow)%p(ii,jj,kk), "@",i,j,k
              ! endif
       endif
    endif
  endif
 enddo
  
 end subroutine particle_to_node_bounce_back_bc2_phase1

 subroutine particle_to_node_bounce_back_bc2_phase2(lrotate,nstep,i,j,k,rversor, &
   otemp,vxs,vys,vzs,rhosub,aoptp, debug,iatm)

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
  type(REALPTR), dimension(0:links):: aoptp

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
          aoptp(indlow)%p(i,j,k)= aoptp(indlow)%p(i,j,k) - &
           p(indhig)*pref_bouzidi*rhosub(ii,jj,kk)*(dex(indhig)*vx+dey(indhig)*vy+dez(indhig)*vz)
           ! if (debug) then
           !   write (iounit(3),*) __FILE__,__LINE__, ii,jj,kk, "pop=", indlow,indhig, &
           !       "rho",rhosub(ii,jj,kk), &
           !       "->pop low", aoptp(indlow)%p(i,j,k),i,j,k
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
          aoptp(indhig)%p(i,j,k) = aoptp(indhig)%p(i,j,k) - &
           p(indlow)*pref_bouzidi*rhosub(ii,jj,kk)*(dex(indlow)*vx+dey(indlow)*vy+dez(indlow)*vz)
          ! if (debug) then
          !    write (iounit(3),*) __FILE__,__LINE__, ii,jj,kk, "pop=", indhig,indlow, &
          !        "rho",rhosub(ii,jj,kk), &
          !        "->pop high", aoptp(indhig)%p(i,j,k),i,j,k
          !  endif
       endif
    endif
  endif
 enddo

 end subroutine particle_to_node_bounce_back_bc2_phase2

 subroutine particle_to_node_bounce_back_bc(lrotate,nstep,i,j,k,rversor, &
   otemp,vxs,vys,vzs,rhosub,aoptp)
 
 
  implicit none
  
  logical, intent(in) :: lrotate
  integer, intent(in) :: nstep,i,j,k
  real(kind=PRC), intent(in) :: vxs,vys,vzs
  real(kind=PRC), intent(in), dimension(3) :: rversor,otemp
  real(kind=PRC), allocatable, dimension(:,:,:)  :: rhosub
  type(REALPTR), dimension(0:links):: aoptp
  
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
      aoptp(1)%p(i,j,k)=real(aoptp(2)%p(ii,jj,kk),kind=PRC)- &
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
      aoptp(2)%p(i,j,k)=real(aoptp(1)%p(ii,jj,kk),kind=PRC)- &
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
      aoptp(3)%p(i,j,k)=real(aoptp(4)%p(ii,jj,kk),kind=PRC)- &
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
      aoptp(4)%p(i,j,k)=real(aoptp(3)%p(ii,jj,kk),kind=PRC)- &
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
      aoptp(5)%p(i,j,k)=real(aoptp(6)%p(ii,jj,kk),kind=PRC)- &
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
      aoptp(6)%p(i,j,k)=real(aoptp(5)%p(ii,jj,kk),kind=PRC)- &
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
      aoptp(7)%p(i,j,k)=real(aoptp(8)%p(ii,jj,kk),kind=PRC)- &
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
      aoptp(8)%p(i,j,k)=real(aoptp(7)%p(ii,jj,kk),kind=PRC)- &
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
      aoptp(9)%p(i,j,k)=real(aoptp(10)%p(ii,jj,kk),kind=PRC)- &
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
      aoptp(10)%p(i,j,k)=real(aoptp(9)%p(ii,jj,kk),kind=PRC)- &
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
      aoptp(11)%p(i,j,k)=real(aoptp(12)%p(ii,jj,kk),kind=PRC)- &
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
      aoptp(12)%p(i,j,k)=real(aoptp(11)%p(ii,jj,kk),kind=PRC)- &
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
      aoptp(13)%p(i,j,k)=real(aoptp(14)%p(ii,jj,kk),kind=PRC)- &
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
      aoptp(14)%p(i,j,k)=real(aoptp(13)%p(ii,jj,kk),kind=PRC)- &
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
      aoptp(15)%p(i,j,k)=real(aoptp(16)%p(ii,jj,kk),kind=PRC)- &
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
      aoptp(16)%p(i,j,k)=real(aoptp(15)%p(ii,jj,kk),kind=PRC)- &
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
      aoptp(17)%p(i,j,k)=real(aoptp(18)%p(ii,jj,kk),kind=PRC)- &
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
      aoptp(18)%p(i,j,k)=real(aoptp(17)%p(ii,jj,kk),kind=PRC)- &
       p(17)*pref_bouzidi*rhosub(ii,jj,kk)*(dey(17)*vy+dez(17)*vz)
    endif
  endif
  
  
  return
  
 end subroutine particle_to_node_bounce_back_bc

 subroutine fixPops(debug, nstep,isub,jsub,ksub, &
                 nspheres,spherelists)
  implicit none
  logical, intent(in) :: debug
  integer, intent(in) :: nstep,isub,jsub,ksub,nspheres
  integer, allocatable, dimension(:,:), intent(in) :: spherelists

  real(kind=PRC), parameter :: MINPOP = 1.0E-15
  integer :: i,j,k,l,ii,jj,kk


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
        if ( f01R(i,j,k)<MINPOP ) f01R(i,j,k) = MINPOP
        if ( f02R(i,j,k)<MINPOP ) f02R(i,j,k) = MINPOP
        if ( f03R(i,j,k)<MINPOP ) f03R(i,j,k) = MINPOP
        if ( f04R(i,j,k)<MINPOP ) f04R(i,j,k) = MINPOP
        if ( f05R(i,j,k)<MINPOP ) f05R(i,j,k) = MINPOP
        if ( f06R(i,j,k)<MINPOP ) f06R(i,j,k) = MINPOP
        if ( f07R(i,j,k)<MINPOP ) f07R(i,j,k) = MINPOP
        if ( f08R(i,j,k)<MINPOP ) f08R(i,j,k) = MINPOP
        if ( f09R(i,j,k)<MINPOP ) f09R(i,j,k) = MINPOP
        if ( f10R(i,j,k)<MINPOP ) f10R(i,j,k) = MINPOP
        if ( f11R(i,j,k)<MINPOP ) f11R(i,j,k) = MINPOP
        if ( f12R(i,j,k)<MINPOP ) f12R(i,j,k) = MINPOP
        if ( f13R(i,j,k)<MINPOP ) f13R(i,j,k) = MINPOP
        if ( f14R(i,j,k)<MINPOP ) f14R(i,j,k) = MINPOP
        if ( f15R(i,j,k)<MINPOP ) f15R(i,j,k) = MINPOP
        if ( f16R(i,j,k)<MINPOP ) f16R(i,j,k) = MINPOP
        if ( f17R(i,j,k)<MINPOP ) f17R(i,j,k) = MINPOP
        if ( f18R(i,j,k)<MINPOP ) f18R(i,j,k) = MINPOP
      endif

      !the fluid bounce back is local so I have to do it
      if(i/=ii.or.j/=jj.or.k/=kk)then
       if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. kk>=minz .and. kk<=maxz)then
        if ( f01R(ii,jj,kk)<MINPOP ) f01R(ii,jj,kk) = MINPOP
        if ( f02R(ii,jj,kk)<MINPOP ) f02R(ii,jj,kk) = MINPOP
        if ( f03R(ii,jj,kk)<MINPOP ) f03R(ii,jj,kk) = MINPOP
        if ( f04R(ii,jj,kk)<MINPOP ) f04R(ii,jj,kk) = MINPOP
        if ( f05R(ii,jj,kk)<MINPOP ) f05R(ii,jj,kk) = MINPOP
        if ( f06R(ii,jj,kk)<MINPOP ) f06R(ii,jj,kk) = MINPOP
        if ( f07R(ii,jj,kk)<MINPOP ) f07R(ii,jj,kk) = MINPOP
        if ( f08R(ii,jj,kk)<MINPOP ) f08R(ii,jj,kk) = MINPOP
        if ( f09R(ii,jj,kk)<MINPOP ) f09R(ii,jj,kk) = MINPOP
        if ( f10R(ii,jj,kk)<MINPOP ) f10R(ii,jj,kk) = MINPOP
        if ( f11R(ii,jj,kk)<MINPOP ) f11R(ii,jj,kk) = MINPOP
        if ( f12R(ii,jj,kk)<MINPOP ) f12R(ii,jj,kk) = MINPOP
        if ( f13R(ii,jj,kk)<MINPOP ) f13R(ii,jj,kk) = MINPOP
        if ( f14R(ii,jj,kk)<MINPOP ) f14R(ii,jj,kk) = MINPOP
        if ( f15R(ii,jj,kk)<MINPOP ) f15R(ii,jj,kk) = MINPOP
        if ( f16R(ii,jj,kk)<MINPOP ) f16R(ii,jj,kk) = MINPOP
        if ( f17R(ii,jj,kk)<MINPOP ) f17R(ii,jj,kk) = MINPOP
        if ( f18R(ii,jj,kk)<MINPOP ) f18R(ii,jj,kk) = MINPOP
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
        if ( f01R(i,j,k)<MINPOP ) f01R(i,j,k) = MINPOP
        if ( f02R(i,j,k)<MINPOP ) f02R(i,j,k) = MINPOP
        if ( f03R(i,j,k)<MINPOP ) f03R(i,j,k) = MINPOP
        if ( f04R(i,j,k)<MINPOP ) f04R(i,j,k) = MINPOP
        if ( f05R(i,j,k)<MINPOP ) f05R(i,j,k) = MINPOP
        if ( f06R(i,j,k)<MINPOP ) f06R(i,j,k) = MINPOP
        if ( f07R(i,j,k)<MINPOP ) f07R(i,j,k) = MINPOP
        if ( f08R(i,j,k)<MINPOP ) f08R(i,j,k) = MINPOP
        if ( f09R(i,j,k)<MINPOP ) f09R(i,j,k) = MINPOP
        if ( f10R(i,j,k)<MINPOP ) f10R(i,j,k) = MINPOP
        if ( f11R(i,j,k)<MINPOP ) f11R(i,j,k) = MINPOP
        if ( f12R(i,j,k)<MINPOP ) f12R(i,j,k) = MINPOP
        if ( f13R(i,j,k)<MINPOP ) f13R(i,j,k) = MINPOP
        if ( f14R(i,j,k)<MINPOP ) f14R(i,j,k) = MINPOP
        if ( f15R(i,j,k)<MINPOP ) f15R(i,j,k) = MINPOP
        if ( f16R(i,j,k)<MINPOP ) f16R(i,j,k) = MINPOP
        if ( f17R(i,j,k)<MINPOP ) f17R(i,j,k) = MINPOP
        if ( f18R(i,j,k)<MINPOP ) f18R(i,j,k) = MINPOP
        if ( f01B(i,j,k)<MINPOP ) f01B(i,j,k) = MINPOP
        if ( f02B(i,j,k)<MINPOP ) f02B(i,j,k) = MINPOP
        if ( f03B(i,j,k)<MINPOP ) f03B(i,j,k) = MINPOP
        if ( f04B(i,j,k)<MINPOP ) f04B(i,j,k) = MINPOP
        if ( f05B(i,j,k)<MINPOP ) f05B(i,j,k) = MINPOP
        if ( f06B(i,j,k)<MINPOP ) f06B(i,j,k) = MINPOP
        if ( f07B(i,j,k)<MINPOP ) f07B(i,j,k) = MINPOP
        if ( f08B(i,j,k)<MINPOP ) f08B(i,j,k) = MINPOP
        if ( f09B(i,j,k)<MINPOP ) f09B(i,j,k) = MINPOP
        if ( f10B(i,j,k)<MINPOP ) f10B(i,j,k) = MINPOP
        if ( f11B(i,j,k)<MINPOP ) f11B(i,j,k) = MINPOP
        if ( f12B(i,j,k)<MINPOP ) f12B(i,j,k) = MINPOP
        if ( f13B(i,j,k)<MINPOP ) f13B(i,j,k) = MINPOP
        if ( f14B(i,j,k)<MINPOP ) f14B(i,j,k) = MINPOP
        if ( f15B(i,j,k)<MINPOP ) f15B(i,j,k) = MINPOP
        if ( f16B(i,j,k)<MINPOP ) f16B(i,j,k) = MINPOP
        if ( f17B(i,j,k)<MINPOP ) f17B(i,j,k) = MINPOP
        if ( f18B(i,j,k)<MINPOP ) f18B(i,j,k) = MINPOP
      endif

      !the fluid bounce back is local so I have to do it
      if(i/=ii.or.j/=jj.or.k/=kk)then
       if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. kk>=minz .and. kk<=maxz)then
        if ( f01R(ii,jj,kk)<MINPOP ) f01R(ii,jj,kk) = MINPOP
        if ( f02R(ii,jj,kk)<MINPOP ) f02R(ii,jj,kk) = MINPOP
        if ( f03R(ii,jj,kk)<MINPOP ) f03R(ii,jj,kk) = MINPOP
        if ( f04R(ii,jj,kk)<MINPOP ) f04R(ii,jj,kk) = MINPOP
        if ( f05R(ii,jj,kk)<MINPOP ) f05R(ii,jj,kk) = MINPOP
        if ( f06R(ii,jj,kk)<MINPOP ) f06R(ii,jj,kk) = MINPOP
        if ( f07R(ii,jj,kk)<MINPOP ) f07R(ii,jj,kk) = MINPOP
        if ( f08R(ii,jj,kk)<MINPOP ) f08R(ii,jj,kk) = MINPOP
        if ( f09R(ii,jj,kk)<MINPOP ) f09R(ii,jj,kk) = MINPOP
        if ( f10R(ii,jj,kk)<MINPOP ) f10R(ii,jj,kk) = MINPOP
        if ( f11R(ii,jj,kk)<MINPOP ) f11R(ii,jj,kk) = MINPOP
        if ( f12R(ii,jj,kk)<MINPOP ) f12R(ii,jj,kk) = MINPOP
        if ( f13R(ii,jj,kk)<MINPOP ) f13R(ii,jj,kk) = MINPOP
        if ( f14R(ii,jj,kk)<MINPOP ) f14R(ii,jj,kk) = MINPOP
        if ( f15R(ii,jj,kk)<MINPOP ) f15R(ii,jj,kk) = MINPOP
        if ( f16R(ii,jj,kk)<MINPOP ) f16R(ii,jj,kk) = MINPOP
        if ( f17R(ii,jj,kk)<MINPOP ) f17R(ii,jj,kk) = MINPOP
        if ( f18R(ii,jj,kk)<MINPOP ) f18R(ii,jj,kk) = MINPOP
        if ( f01B(ii,jj,kk)<MINPOP ) f01B(ii,jj,kk) = MINPOP
        if ( f02B(ii,jj,kk)<MINPOP ) f02B(ii,jj,kk) = MINPOP
        if ( f03B(ii,jj,kk)<MINPOP ) f03B(ii,jj,kk) = MINPOP
        if ( f04B(ii,jj,kk)<MINPOP ) f04B(ii,jj,kk) = MINPOP
        if ( f05B(ii,jj,kk)<MINPOP ) f05B(ii,jj,kk) = MINPOP
        if ( f06B(ii,jj,kk)<MINPOP ) f06B(ii,jj,kk) = MINPOP
        if ( f07B(ii,jj,kk)<MINPOP ) f07B(ii,jj,kk) = MINPOP
        if ( f08B(ii,jj,kk)<MINPOP ) f08B(ii,jj,kk) = MINPOP
        if ( f09B(ii,jj,kk)<MINPOP ) f09B(ii,jj,kk) = MINPOP
        if ( f10B(ii,jj,kk)<MINPOP ) f10B(ii,jj,kk) = MINPOP
        if ( f11B(ii,jj,kk)<MINPOP ) f11B(ii,jj,kk) = MINPOP
        if ( f12B(ii,jj,kk)<MINPOP ) f12B(ii,jj,kk) = MINPOP
        if ( f13B(ii,jj,kk)<MINPOP ) f13B(ii,jj,kk) = MINPOP
        if ( f14B(ii,jj,kk)<MINPOP ) f14B(ii,jj,kk) = MINPOP
        if ( f15B(ii,jj,kk)<MINPOP ) f15B(ii,jj,kk) = MINPOP
        if ( f16B(ii,jj,kk)<MINPOP ) f16B(ii,jj,kk) = MINPOP
        if ( f17B(ii,jj,kk)<MINPOP ) f17B(ii,jj,kk) = MINPOP
        if ( f18B(ii,jj,kk)<MINPOP ) f18B(ii,jj,kk) = MINPOP
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
!     author: M. Lauricella
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
end subroutine compute_onebelt_density

 subroutine compute_onebelt_density_twofluids(i,j,k, Rsum, Bsum, Dsum)
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
end subroutine compute_onebelt_density_twofluids

 subroutine compute_secbelt_density(i,j,k, Rsum, Dsum)
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
end subroutine compute_secbelt_density

 subroutine compute_secbelt_density_twofluids(i,j,k, Rsum, Bsum, Dsum)
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
end subroutine compute_secbelt_density_twofluids


 subroutine fix_onebelt_density(i,j,k, Rsum, Dsum)
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
!     author: M. Lauricella
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
            call initialize_newnode_fluid(i,j,k,Rsum,myu,myv,myw,aoptpR)
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
            call initialize_newnode_fluid(i,j,k,Rsum,myu,myv,myw,aoptpR)
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
            call initialize_newnode_fluid(i,j,k,Rsum,myu,myv,myw,aoptpR)
            call initialize_newnode_fluid(i,j,k,Bsum,myu,myv,myw,aoptpB)
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
            call initialize_newnode_fluid(i,j,k,Rsum,myu,myv,myw,aoptpR)
            call initialize_newnode_fluid(i,j,k,Bsum,myu,myv,myw,aoptpB)
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
!     author: M. Lauricella
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

 
 subroutine fill_new_isfluid_inlist(isub,jsub,ksub, nspheres, spherelists, fvalue)
 implicit none
 integer, intent(in) :: isub,jsub,ksub,nspheres, fvalue
 integer, allocatable, dimension(:,:), intent(in) :: spherelists
  integer, save :: imin,imax,jmin,jmax,kmin,kmax
  logical, save :: lfirst=.true.
  integer :: l, i,j,k


  if(lfirst)then
    lfirst=.false.
    imin=minx-1
    imax=maxx+1
    jmin=miny-1
    jmax=maxy+1
    kmin=minz-1
    kmax=maxz+1
  endif

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

      new_isfluid(i,j,k)= fvalue
  enddo
 end subroutine fill_new_isfluid_inlist


 subroutine mapping_new_isfluid(natmssub,natms_ext,atmbook,nspheres,spherelists, &
   spheredists,nspheredeads,spherelistdeads,lmoved,xx,yy,zz)
  
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
  integer, intent(in) :: natmssub,natms_ext, nspheres,nspheredeads
  integer, allocatable, intent(in) :: atmbook(:)
  integer, allocatable, dimension(:,:), intent(in) :: spherelists, &
   spherelistdeads
  real(kind=PRC), allocatable, dimension(:), intent(in) :: spheredists
  logical(kind=1),allocatable, dimension(:), intent(in) :: lmoved
  real(kind=PRC), allocatable, dimension(:), intent(in) :: xx,yy,zz
  integer :: myi, isub,jsub,ksub, iatm, i,j,k, iofile(2)
  
  
  !delete fluid 
  do myi=1,natms_ext
    iatm=atmbook(myi)
    isub=nint(xx(iatm))
    jsub=nint(yy(iatm))
    ksub=nint(zz(iatm))

    call fill_new_isfluid_inlist(isub,jsub,ksub, nspheres, spherelists, 2)

    call fill_new_isfluid_inlist(isub,jsub,ksub, nspheredeads, spherelistdeads, 4)

    CYCLE_OUT_INTERVAL(isub, minx-nbuff, maxx+nbuff)
    CYCLE_OUT_INTERVAL(jsub, miny-nbuff, maxy+nbuff)
    CYCLE_OUT_INTERVAL(ksub, minz-nbuff, maxz+nbuff)

    ! Fill internal atoms
    new_isfluid(isub,jsub,ksub)=5
  enddo

 end subroutine mapping_new_isfluid

 
 subroutine initialize_new_isfluid()
 
!***********************************************************************
!     
!     LBsoft subroutine for initialize the new_isfluid
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
 
  implicit none
  
  where(isfluid(minx-nbuff:maxx+nbuff,miny-nbuff:maxy+nbuff,minz-nbuff:maxz+nbuff)==2 &
   .or. isfluid(minx-nbuff:maxx+nbuff,miny-nbuff:maxy+nbuff,minz-nbuff:maxz+nbuff)==4 &
   .or. isfluid(minx-nbuff:maxx+nbuff,miny-nbuff:maxy+nbuff,minz-nbuff:maxz+nbuff)==5)
    new_isfluid(minx-nbuff:maxx+nbuff,miny-nbuff:maxy+nbuff,minz-nbuff:maxz+nbuff)=1
  elsewhere
    new_isfluid(minx-nbuff:maxx+nbuff,miny-nbuff:maxy+nbuff,minz-nbuff:maxz+nbuff)= &
     isfluid(minx-nbuff:maxx+nbuff,miny-nbuff:maxy+nbuff,minz-nbuff:maxz+nbuff)
  end where
  
 end subroutine initialize_new_isfluid
 
 subroutine update_isfluid()
 
!***********************************************************************
!     
!     LBsoft subroutine for update isfluid from new_isfluid
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
 
  implicit none
  integer :: i,j,k
  
  
  ! forall(i=minx-nbuff:maxx+nbuff,j=miny-nbuff:maxy+nbuff,k=minz-nbuff:maxz+nbuff)
  do k=minz-nbuff,maxz+nbuff
    do j=miny-nbuff,maxy+nbuff
      do i=minx-nbuff,maxx+nbuff
         isfluid(i,j,k)=new_isfluid(i,j,k)
      enddo
    enddo
  enddo
  !end forall
 end subroutine update_isfluid
 

 subroutine initialize_newnode_fluid(i,j,k,rhosub,usub,vsub,wsub,aoptp)
 
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
  type(REALPTR), dimension(0:links):: aoptp
  
  !formula taken from eq. 19 of PRE 83, 046707 (2011)
    aoptp(0)%p(i,j,k)=equil_pop00(rhosub,usub,vsub, &
     wsub)
  
  
  
    aoptp(1)%p(i,j,k)=equil_pop01(rhosub,usub,vsub, &
     wsub)
  
  
  
    aoptp(2)%p(i,j,k)=equil_pop02(rhosub,usub,vsub, &
     wsub)
  
  
  
    aoptp(3)%p(i,j,k)=equil_pop03(rhosub,usub,vsub, &
     wsub)
  
  
  
    aoptp(4)%p(i,j,k)=equil_pop04(rhosub,usub,vsub, &
     wsub)
  
  
  
    aoptp(5)%p(i,j,k)=equil_pop05(rhosub,usub,vsub, &
     wsub)
  
  
  
    aoptp(6)%p(i,j,k)=equil_pop06(rhosub,usub,vsub, &
     wsub)
  
  
  
    aoptp(7)%p(i,j,k)=equil_pop07(rhosub,usub,vsub, &
     wsub)
  
  
  
    aoptp(8)%p(i,j,k)=equil_pop08(rhosub,usub,vsub, &
     wsub)
  
  
  
    aoptp(9)%p(i,j,k)=equil_pop09(rhosub,usub,vsub, &
     wsub)
  
  
  
    aoptp(10)%p(i,j,k)=equil_pop10(rhosub,usub,vsub, &
     wsub)
  
  
  
    aoptp(11)%p(i,j,k)=equil_pop11(rhosub,usub,vsub, &
     wsub)
  
  
  
    aoptp(12)%p(i,j,k)=equil_pop12(rhosub,usub,vsub, &
     wsub)
  
  
  
    aoptp(13)%p(i,j,k)=equil_pop13(rhosub,usub,vsub, &
     wsub)
  
  
  
    aoptp(14)%p(i,j,k)=equil_pop14(rhosub,usub,vsub, &
     wsub)
  
  
  
    aoptp(15)%p(i,j,k)=equil_pop15(rhosub,usub,vsub, &
     wsub)
  
  
  
    aoptp(16)%p(i,j,k)=equil_pop16(rhosub,usub,vsub, &
     wsub)
  
  
  
    aoptp(17)%p(i,j,k)=equil_pop17(rhosub,usub,vsub, &
     wsub)
  
  
  
    aoptp(18)%p(i,j,k)=equil_pop18(rhosub,usub,vsub, &
     wsub)
  
  
  return
  
 end subroutine initialize_newnode_fluid
 
!***********END PART TO MANAGE THE INTERACTION WITH PARTICLES***********

 subroutine write_test_map_pop(aoptp)
  
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
  
  type(REALPTR), dimension(0:links):: aoptp
  
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
            write(iomap,'(f20.10)')aoptp(l)%p(i,j,k) 
          enddo
        enddo
      enddo
    enddo

  
  close(iomap)
  
 end subroutine write_test_map_pop
 
 subroutine test_fake_dens(dtemp,ltestout)
 
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
 
 subroutine test_fake_pops(aoptp,ltestout)
 
  implicit none
  
  logical :: lwritesub
  type(REALPTR), dimension(0:links):: aoptp
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
              write(iotest)i,j,k,l,real(aoptp(l)%p(i,j,k),kind=PRC), &
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
              write(iotest,'(4i8,g20.10,4i4)')i,j,k,l,aoptp(l)%p(i,j,k), &
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
              if(real(aoptp(l)%p(i,j,k),kind=PRC)/=dtemp2(i,j,k,l))then
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
              if(dabs(real(aoptp(l)%p(i,j,k),kind=PRC)-dtemp2(i,j,k,l))>1.d-9)then
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
      write(6,*)'node = ',i,j,k,l,aoptp(l)%p(i,j,k),dtemp2(i,j,k,l)
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
 
 subroutine print_all_pops(iosub,filenam,itersub,aoptp)
 
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
  type(REALPTR), dimension(0:links):: aoptp
  
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
            write(iosub*idrank+23,*)i,j,k,l,aoptp(l)%p(i,j,k)
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
  rdimx,rdimy,rdimz,aoptp,xc,yc,zc,fmiosss,nspheres, &
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
  type(REALPTR), dimension(0:links):: aoptp
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
               aoptp(l)%p(ii,jj,kk)
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
               aoptp(opp(l))%p(i2,j2,k2)
      endif
    enddo
  enddo
  
  close(iosub*idrank+23)
  
  return
  
 end subroutine print_all_pops_center
 
 subroutine print_all_pops_area_shpere(iosub,filenam,itersub, &
  rdimx,rdimy,rdimz,aoptp,xc,yc,zc,fmiosss,nspheres, &
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
  type(REALPTR), dimension(0:links):: aoptp
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
          miapop(m,l)=aoptp(m)%p(itar,jtar,ktar)
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
            miapop(m,l)=aoptp(m)%p(itar,jtar,ktar)
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
            miapop(m,l)=aoptp(m)%p(itar,jtar,ktar)
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
       aoptpR(l)%p(i,j,k) = i*1000000 + j*10000 + k*100 + l
     end forall

     where(isfluid(minx-1:maxx+1,miny-1:maxy+1,minz-1:maxz+1)==0)
       aoptpR(l)%p(minx-1:maxx+1,miny-1:maxy+1,minz-1:maxz+1) = &
        - aoptpR(l)%p(minx-1:maxx+1,miny-1:maxy+1,minz-1:maxz+1)
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

          destVal = aoptpR(l)%p(i,j,k)

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
            write(iosub,*)i,j,k,l,aoptpR(l)%p(i,j,k),aoptpB(l)%p(i,j,k)
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
    read(133) aoptpR(l)%p(minx:maxx,miny:maxy,minz:maxz)
    read(133) aoptpB(l)%p(minx:maxx,miny:maxy,minz:maxz)
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
    reshape(f00R(1:nx,1:ny,1:nz),(/nn/)), &
    reshape(f01R(1:nx,1:ny,1:nz),(/nn/)), &
    reshape(f02R(1:nx,1:ny,1:nz),(/nn/)), &
    reshape(f03R(1:nx,1:ny,1:nz),(/nn/)), &
    reshape(f04R(1:nx,1:ny,1:nz),(/nn/)), &
    reshape(f05R(1:nx,1:ny,1:nz),(/nn/)), &
    reshape(f06R(1:nx,1:ny,1:nz),(/nn/)), &
    reshape(f07R(1:nx,1:ny,1:nz),(/nn/)), &
    reshape(f08R(1:nx,1:ny,1:nz),(/nn/)), &
    reshape(f09R(1:nx,1:ny,1:nz),(/nn/)), &
    reshape(f10R(1:nx,1:ny,1:nz),(/nn/)), &
    reshape(f11R(1:nx,1:ny,1:nz),(/nn/)), &
    reshape(f12R(1:nx,1:ny,1:nz),(/nn/)), &
    reshape(f13R(1:nx,1:ny,1:nz),(/nn/)), &
    reshape(f14R(1:nx,1:ny,1:nz),(/nn/)), &
    reshape(f15R(1:nx,1:ny,1:nz),(/nn/)), &
    reshape(f16R(1:nx,1:ny,1:nz),(/nn/)), &
    reshape(f17R(1:nx,1:ny,1:nz),(/nn/)), &
    reshape(f18R(1:nx,1:ny,1:nz),(/nn/)), &
    ierror)
   call close_file_binary(myo,ierror)
     
   if(.not. lsingle_fluid) then
     myo=myunit+2
     mynamefile=repeat(' ',120)
     mynamefile='dumpPopsB.' // mysave //'.dat'
     call open_file_binary(myo,trim(mynamefile),ierror)
     call print_binary_pop_array(myo,nn, &
      reshape(f00B(1:nx,1:ny,1:nz),(/nn/)), &
      reshape(f01B(1:nx,1:ny,1:nz),(/nn/)), &
      reshape(f02B(1:nx,1:ny,1:nz),(/nn/)), &
      reshape(f03B(1:nx,1:ny,1:nz),(/nn/)), &
      reshape(f04B(1:nx,1:ny,1:nz),(/nn/)), &
      reshape(f05B(1:nx,1:ny,1:nz),(/nn/)), &
      reshape(f06B(1:nx,1:ny,1:nz),(/nn/)), &
      reshape(f07B(1:nx,1:ny,1:nz),(/nn/)), &
      reshape(f08B(1:nx,1:ny,1:nz),(/nn/)), &
      reshape(f09B(1:nx,1:ny,1:nz),(/nn/)), &
      reshape(f10B(1:nx,1:ny,1:nz),(/nn/)), &
      reshape(f11B(1:nx,1:ny,1:nz),(/nn/)), &
      reshape(f12B(1:nx,1:ny,1:nz),(/nn/)), &
      reshape(f13B(1:nx,1:ny,1:nz),(/nn/)), &
      reshape(f14B(1:nx,1:ny,1:nz),(/nn/)), &
      reshape(f15B(1:nx,1:ny,1:nz),(/nn/)), &
      reshape(f16B(1:nx,1:ny,1:nz),(/nn/)), &
      reshape(f17B(1:nx,1:ny,1:nz),(/nn/)), &
      reshape(f18B(1:nx,1:ny,1:nz),(/nn/)), &
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

   call collective_writeFile_pops(nstep, mynamefile, &
     f00R,f01R,f02R,f03R,f04R, &
     f05R,f06R,f07R,f08R,f09R,f10R,f11R,f12R,f13R, &
     f14R,f15R,f16R,f17R,f18R, nbuff)

   if(.not. lsingle_fluid) then
    mynamefile=repeat(' ',120)
    mynamefile='dumpPopsB.' // mysave //'.dat'
    call collective_writeFile_pops(nstep, mynamefile, &
     f00B,f01B,f02B,f03B,f04B, &
     f05B,f06B,f07B,f08B,f09B,f10B,f11B,f12B,f13B, &
     f14B,f15B,f16B,f17B,f18B, nbuff)
   endif

   ! Conservative: dump also hvars & isfluid
   mynamefile=repeat(' ',120)
   mynamefile='dumpRhoR.' // mysave //'.dat'
   call collective_writeFile(nstep, mynamefile, rhoR, nbuff)

   if(.not. lsingle_fluid) then
    mynamefile=repeat(' ',120)
    mynamefile='dumpRhoB.' // mysave //'.dat'
    call collective_writeFile(nstep, mynamefile, rhoB, nbuff)
   endif

   mynamefile=repeat(' ',120)
   mynamefile='dumpU.' // mysave //'.dat'
   call collective_writeFile(nstep, mynamefile, u, nbuff)

   mynamefile=repeat(' ',120)
   mynamefile='dumpV.' // mysave //'.dat'
   call collective_writeFile(nstep, mynamefile, v, nbuff)

   mynamefile=repeat(' ',120)
   mynamefile='dumpW.' // mysave //'.dat'
   call collective_writeFile(nstep, mynamefile, w, nbuff)

   mynamefile=repeat(' ',120)
   mynamefile='dumpisFluid.' // mysave //'.dat'
   call collective_writeFile_int(nstep, mynamefile, isfluid, nbuff)


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
     call collective_readFile_pops(nstep, mynamefile, &
       f00R,f01R,f02R,f03R,f04R, &
       f05R,f06R,f07R,f08R,f09R,f10R,f11R,f12R,f13R, &
       f14R,f15R,f16R,f17R,f18R, nbuff)

     if(.not. lsingle_fluid) then
      mynamefile=repeat(' ',120)
      mynamefile='dumpPopsB.restart.dat'
      inquire(file=trim(mynamefile),exist=lexist)
      if(.not. lexist)then
        call warning(54,real(1.0,kind=PRC),trim(mynamefile))
        call error(36)
      endif
      call collective_readFile_pops(nstep, mynamefile, &
       f00B,f01B,f02B,f03B,f04B, &
       f05B,f06B,f07B,f08B,f09B,f10B,f11B,f12B,f13B, &
       f14B,f15B,f16B,f17B,f18B, nbuff)
     endif

     ! Conservative: read also hvars & isfluid
     mynamefile=repeat(' ',120)
     mynamefile='dumpRhoR.restart.dat'
     inquire(file=trim(mynamefile),exist=lexist)
     if(.not. lexist)then
       call warning(54,real(1.0,kind=PRC),trim(mynamefile))
       call error(36)
     endif
     call collective_readFile(nstep, mynamefile, rhoR, nbuff)

     if(.not. lsingle_fluid) then
      mynamefile=repeat(' ',120)
      mynamefile='dumpRhoB.restart.dat'
      inquire(file=trim(mynamefile),exist=lexist)
      if(.not. lexist)then
        call warning(54,real(1.0,kind=PRC),trim(mynamefile))
        call error(36)
      endif
      call collective_readFile(nstep, mynamefile, rhoB, nbuff)
     endif

     mynamefile=repeat(' ',120)
     mynamefile='dumpU.restart.dat'
     inquire(file=trim(mynamefile),exist=lexist)
     if(.not. lexist)then
       call warning(54,real(1.0,kind=PRC),trim(mynamefile))
       call error(36)
     endif
     call collective_readFile(nstep, mynamefile, u, nbuff)

     mynamefile=repeat(' ',120)
     mynamefile='dumpV.restart.dat'
     inquire(file=trim(mynamefile),exist=lexist)
     if(.not. lexist)then
       call warning(54,real(1.0,kind=PRC),trim(mynamefile))
       call error(36)
     endif
     call collective_readFile(nstep, mynamefile, v, nbuff)

     mynamefile=repeat(' ',120)
     mynamefile='dumpW.restart.dat'
     inquire(file=trim(mynamefile),exist=lexist)
     if(.not. lexist)then
       call warning(54,real(1.0,kind=PRC),trim(mynamefile))
       call error(36)
     endif
     call collective_readFile(nstep, mynamefile, w, nbuff)

     mynamefile=repeat(' ',120)
     mynamefile='dumpisFluid.restart.dat'
     inquire(file=trim(mynamefile),exist=lexist)
     if(.not. lexist)then
       call warning(54,real(1.0,kind=PRC),trim(mynamefile))
       call error(36)
     endif
     call collective_readFile_int(nstep, mynamefile, isfluid, nbuff)

     call driver_bc_densities
     call driver_bc_velocities
     call driver_bc_isfluid
  end subroutine restore_oneFile
  
 subroutine read_isfluid_dat(nstep)
  
!***********************************************************************
!     
!     LBsoft subroutine for reading the isfluid.dat file
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
  
  if(.not. lread_isfluid)return
  
  mynamefile=repeat(' ',120)
  mynamefile='isfluid.dat'
  inquire(file=trim(mynamefile),exist=lexist)
  if(.not. lexist)then
    call warning(63)
    call error(7)
  endif
  call collective_readFile_int(nstep, mynamefile, isfluid, nbuff)
     
  return
  
 end subroutine read_isfluid_dat

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

     call collective_writeFile(nstep, mynamefile, rhoR, nbuff)

     if(.not. lsingle_fluid) then
      mynamefile=repeat(' ',120)
      mynamefile='dumpBin/RhoB.' // write_fmtnumb(nstep) //'.raw'
      call collective_writeFile(nstep, mynamefile, rhoB, nbuff)
     endif

     mynamefile=repeat(' ',120)
     mynamefile='dumpBin/U.' // write_fmtnumb(nstep) //'.raw'
     call collective_writeFile(nstep, mynamefile, u, nbuff)

     mynamefile=repeat(' ',120)
     mynamefile='dumpBin/V.' // write_fmtnumb(nstep) //'.raw'
     call collective_writeFile(nstep, mynamefile, v, nbuff)

     mynamefile=repeat(' ',120)
     mynamefile='dumpBin/W.' // write_fmtnumb(nstep) //'.raw'
     call collective_writeFile(nstep, mynamefile, w, nbuff)

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
