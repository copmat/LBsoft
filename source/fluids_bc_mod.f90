
#include <default_macro.h>

 module fluids_bc_mod
 
!***********************************************************************
!     
!     LBsoft module containing variable and subroutines of the
!     fluid boundary conditions
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification May 2020
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
                   collective_writeFile_int,gminx,gmaxx, &
                   gminy,gmaxy,gminz,gmaxz
 use error_mod
 use utility_mod, only : Pi,modulvec,cross,dot,gauss,ibuffservice, &
                   allocate_array_ibuffservice,buffservice,conv_rad, &
                   allocate_array_buffservice,lbuffservice, &
                   allocate_array_lbuffservice,xcross,ycross,zcross, &
                   nbuffservice3d,buffservice3d,int_cube_sphere,fcut, &
                   rand_noseeded,linit_seed,gauss_noseeded,write_fmtnumb, &
                   openLogFile

 
 implicit none
 
 private
 
 integer, save, protected, public :: ibctype=0
 integer, save, public :: ixpbc, iypbc, izpbc
 integer(kind=1), save, public, allocatable, dimension(:,:,:) :: isfluid
 integer(kind=1), save, public, allocatable, dimension(:,:,:) :: new_isfluid
 integer(kind=1), save, public, allocatable, dimension(:,:,:) :: bcfluid
 
 integer, save :: ix,iy,iz,mynx,myny,mynz,nx,ny,nz
 integer, save :: minx, maxx, miny, maxy, minz, maxz
 
 !is there a isfluid.dat file to read? 
 logical, save, protected, public :: lread_isfluid=.false.
 
 !selected wetting mode 0=average 1=fixed value (Benzi's style!)
 integer, save, protected, public :: iselwetting=0
 
 !fixed density to set the wall wetting as in Benzi's style
 real(kind=PRC), save, protected, public :: densR_wetting = ZERO
 real(kind=PRC), save, protected, public :: densB_wetting = ZERO
 
 integer, save, protected, public :: bc_type_east=0
 integer, save, protected, public :: bc_type_west=0
 integer, save, protected, public :: bc_type_rear=0
 integer, save, protected, public :: bc_type_front=0
 integer, save, protected, public :: bc_type_north=0
 integer, save, protected, public :: bc_type_south=0
 
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
 
 integer, parameter, public :: nbcdir=6
 
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
 
 real(kind=PRC), save, protected, public, allocatable, dimension(:) :: bc_rhoR
 real(kind=PRC), save, protected, public, allocatable, dimension(:) :: bc_rhoB
 real(kind=PRC), save, protected, public, allocatable, dimension(:) :: bc_u
 real(kind=PRC), save, protected, public, allocatable, dimension(:) :: bc_v
 real(kind=PRC), save, protected, public, allocatable, dimension(:) :: bc_w
 real(kind=PRC), save, protected, public, allocatable, dimension(:) :: bc_flow
 real(kind=PRC), save, protected, public, allocatable, dimension(:) :: bc_omega
 
!array containing the node list of the second belt 
 integer, save, public :: ndouble
 integer, save, public, allocatable, dimension(:) :: exdouble
 integer, save, public, allocatable, dimension(:) :: eydouble
 integer, save, public, allocatable, dimension(:) :: ezdouble
 
 public :: set_boundary_conditions_type
 public :: allocate_isfluid
 public :: set_isfluid_bcfluid
 public :: set_double_range
 public :: set_lread_isfluid
 public :: set_fluid_wall_mode
 public :: set_value_bc_east
 public :: set_value_bc_west
 public :: set_value_bc_rear
 public :: set_value_bc_front
 public :: set_value_bc_north
 public :: set_value_bc_south
 public :: pimage
 public :: compute_densities_wall
 public :: driver_bc_new_isfluid
 public :: driver_bc_isfluid
 public :: mapping_new_isfluid
 public :: initialize_new_isfluid
 public :: update_isfluid
 public :: set_bc_variable_hvar
 public :: manage_bc_hvar_selfcomm
 public :: apply_bounceback_pop_halfway
 public :: manage_bc_singlehalo_selfcomm
 public :: manage_bc_pop_selfcomm
 public :: initialiaze_manage_bc_pop_selfcomm
 
 contains
 
 subroutine set_boundary_conditions_type(itemp1,itemp2,itemp3)
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the ibctype variable taking into
!     account the periodic bc
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification May 2020
!     
!***********************************************************************
   
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
 
 subroutine allocate_isfluid(lparticles)
 
!***********************************************************************
!     
!     LBsoft subroutine for allocating arrays which 
!     describe the fluid boundary conditions
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification June 2020
!     
!***********************************************************************

  implicit none
     
  logical, intent(in) :: lparticles
  integer, parameter :: nistatmax=10
  integer, dimension(nistatmax) :: istat
  integer :: i,j,k,l
  logical, dimension(1) :: ltest=.false.
  
  minx=gminx(idrank)
  maxx=gmaxx(idrank)
  miny=gminy(idrank)
  maxy=gmaxy(idrank)
  minz=gminz(idrank)
  maxz=gmaxz(idrank)
  
  nx=gmaxx(-1)
  ny=gmaxy(-1)
  nz=gmaxz(-1)
  
  ix=minx-nbuff
  iy=miny-nbuff
  iz=minz-nbuff
  mynx=maxx+nbuff
  myny=maxy+nbuff
  mynz=maxz+nbuff
  
  istat=0
  
  allocate(isfluid(ix:mynx,iy:myny,iz:mynz),stat=istat(1))
     
  if(lparticles)then
    allocate(new_isfluid(ix:mynx,iy:myny,iz:mynz),stat=istat(2))
  endif
     
  allocate(bcfluid(ix:mynx,iy:myny,iz:mynz),stat=istat(3))
  ltest=.false.
  if(any(istat.ne.0))then
    do i=1,nistatmax
      if(istat(i).ne.0)exit
    enddo
    call warning(2,real(i,kind=PRC))
    ltest=.true.
  endif

  call or_world_larr(ltest,1)
  if(ltest(1))call error(4)
     
  isfluid(ix:mynx,iy:myny,iz:mynz)=3
  bcfluid(ix:mynx,iy:myny,iz:mynz)=0
  
  return
  
 end subroutine allocate_isfluid
 
 subroutine set_isfluid_bcfluid(lvtkfilesub,lsingle_fluid,lparticles, &
  lunique_omega,lexch_dens,lexch_u,lexch_v,lexch_w,unique_omega,omega, &
  rhoR,rhoB,u,v,w)
 
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
  
  logical, intent(in) :: lvtkfilesub,lsingle_fluid,lparticles, &
   lunique_omega
  logical, intent(inout) :: lexch_dens
  logical, intent(inout) :: lexch_u
  logical, intent(inout) :: lexch_v
  logical, intent(inout) :: lexch_w
  real(kind=PRC), intent(in) :: unique_omega
  real(kind=PRC), allocatable, dimension(:,:,:)  :: omega,rhoR,rhoB,u,v,w
  
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
  
  isfluid(minx:maxx,miny:maxy,minz:maxz)=1
  
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
  call initialiaze_manage_bc_selfcomm
  
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
  call set_bc_hvar(lsingle_fluid,lunique_omega,unique_omega,omega, &
   rhoR,rhoB,u,v,w)
  
  ltest=.false.
  if(any(istat.ne.0))then
    do i=1,nistatmax
      if(istat(i).ne.0)exit
    enddo
    call warning(2,real(i,kind=PRC))
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

 end subroutine set_isfluid_bcfluid
 
 subroutine set_double_range()
 
!***********************************************************************
!     
!     LBsoft subroutine for initialize the double directions of 
!     over the lattice
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification June 2020
!     
!***********************************************************************
 
  implicit none
  
  integer :: i,j,k,l
  
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
    
    return
    
 end subroutine set_double_range
 
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
  call collective_readFile_int(nstep,mynamefile,isfluid,nbuff, &
   minx,maxx,miny,maxy,minz,maxz)
     
  return
  
 end subroutine read_isfluid_dat
 
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
  
  logical, save   :: isFirst = .true.

  call commexch_isfluid(isfluid)
  
  call manage_bc_isfluid_selfcomm(isfluid)
    
   
  call commwait_isfluid(isfluid)
  
  return 

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

 subroutine initialiaze_manage_bc_selfcomm()
 
!***********************************************************************
!     
!     LBsoft subroutine to initialize the list of nodes
!     within the same process for applying the boundary conditions
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification May 2020
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
   
 end subroutine initialiaze_manage_bc_selfcomm
 
 subroutine manage_bc_isfluid_selfcomm(istemp)
 
!***********************************************************************
!     
!     LBsoft subroutine to manage the buffer of isfluid nodes
!     within the same process for applying the boundary conditions
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification May 2020
!     
!***********************************************************************
 
  implicit none
  
  integer(kind=1), dimension(:,:,:), allocatable :: istemp
  
  integer :: i,j,k,l,itemp,jtemp,ktemp,itemp2,jtemp2,ktemp2
  integer(kind=IPRC) :: i4orig,i4
  
  do i=1,nguards
    istemp(isguards(4,i),isguards(5,i),isguards(6,i))= &
     istemp(isguards(1,i),isguards(2,i),isguards(3,i))
  enddo
   
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
 
 subroutine manage_bc_singlehalo_selfcomm(dtemp,nbuff,ix,fx,iy,fy,iz,fz)
 
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
  
  integer, intent(in) :: nbuff,ix,fx,iy,fy,iz,fz
  real(kind=PRC), &
   dimension(ix-nbuff:fx+nbuff,iy-nbuff:fy+nbuff,iz-nbuff:fz+nbuff) :: &
   dtemp
  
  integer :: i,j,k,l,itemp,jtemp,ktemp,itemp2,jtemp2,ktemp2
  integer(kind=IPRC) :: i4orig,i4
  
  do i=1,nsinglehalo
    dtemp(isinglehalo(4,i),isinglehalo(5,i),isinglehalo(6,i))= &
     dtemp(isinglehalo(1,i),isinglehalo(2,i),isinglehalo(3,i))
  enddo
   
  return
   
 end subroutine manage_bc_singlehalo_selfcomm
 
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
  
  real(kind=PRC), allocatable, dimension(:,:,:) :: dtemp
  
  integer :: i
  
  do i=1,nguards
    dtemp(isguards(4,i),isguards(5,i),isguards(6,i))= &
     dtemp(isguards(1,i),isguards(2,i),isguards(3,i))
  enddo
   
  return
   
 end subroutine manage_bc_hvar_selfcomm
 
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
 
 subroutine manage_bc_pop_selfcomm(fluidsub,lparticles)
 
!***********************************************************************
!     
!     LBsoft subroutine to manage the population buffer fluid nodes
!     within the same process for applying the boundary conditions
!     using the node list created in subroutine
!     initialiaze_manage_bc_pop_selfcomm
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification May 2020
!     
!***********************************************************************
 
  implicit none
  
  logical, intent(in) :: lparticles
  
  real(kind=PRC), allocatable, dimension(:,:,:,:)  :: fluidsub
  
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
        fluidsub(itemp,jtemp,ktemp,idir)=fluidsub(itemp2,jtemp2,ktemp2,idir)
      endif
    enddo
  else
    do i=1,npoplistbc
      idir=ipoplistbc(i)
      fluidsub(poplistbc(1,i),poplistbc(2,i),poplistbc(3,i),idir)= &
      fluidsub(poplistbc(4,i),poplistbc(5,i),poplistbc(6,i),idir)
    enddo
  endif
  
  return
   
 end subroutine manage_bc_pop_selfcomm
 
 subroutine set_bc_hvar(lsingle_fluid,lunique_omega,unique_omega,omega, &
  rhoR,rhoB,u,v,w)
 
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
  
  logical, intent(in) :: lsingle_fluid,lunique_omega
  real(kind=PRC), intent(in) :: unique_omega
  real(kind=PRC), allocatable, dimension(:,:,:)  :: omega,rhoR,rhoB,u,v,w
  
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
 
 subroutine set_bc_fixed_hvar(lsingle_fluid)
 
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
  
  logical, intent(in) :: lsingle_fluid
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
 
 subroutine set_bc_variable_hvar(lsingle_fluid,lunique_omega, &
  unique_omega,omega,rhoR,rhoB,u,v,w)
 
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
  
  logical, intent(in) :: lsingle_fluid,lunique_omega
  real(kind=PRC), intent(in) :: unique_omega
  real(kind=PRC), allocatable, dimension(:,:,:) :: omega,rhoR,rhoB,u,v,w
  
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
 
 subroutine apply_bounceback_pop_halfway(nstep,lColourG,rho_s,u_s,v_s,w_s, &
  rhoR,rhoB,fluidsub,alpha_CG)
 
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
  logical, intent(in) :: lColourG
  real(kind=PRC), allocatable, dimension(:)  :: rho_s,u_s,v_s,w_s
  
  real(kind=PRC), allocatable, dimension(:,:,:) :: rhoR,rhoB
  
  real(kind=PRC), allocatable :: fluidsub(:,:,:,:)
  real(kind=PRC), intent(in) :: alpha_CG
  
  integer :: i,j,k,l,sx,sz,iloop,indlow,indhig,ii,jj,kk,io,jo,ko,iii
  
  real(kind=PRC), dimension(0:links) :: myphi_CG
  real(kind=PRC), parameter :: pref_bouzidi=TWO/cssq
  real(kind=PRC), parameter :: cssq2 = ( HALF / cssq )
  real(kind=PRC), parameter :: cssq4 = ( HALF / (cssq*cssq) )
  real(kind=PRC) :: grad_rhox,grad_rhoy,grad_rhoz,locrhoR,locrhoB
  
  

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
	        fluidsub(i,j,k,indlow) = fluidsub(ii,jj,kk,indhig)
          endif
        endif
        
        ii=i+ex(indhig)
        jj=j+ey(indhig)
        kk=k+ez(indhig)
        
        if(ii>=minx .and. ii<=maxx .and. jj>=miny .and. jj<=maxy .and. &
         kk>=minz .and. kk<=maxz)then
          if(isfluid(ii,jj,kk)==1)then
  	        fluidsub(i,j,k,indhig) = fluidsub(ii,jj,kk,indlow)
  	      endif
        endif
      enddo
    enddo
  
  
  
  
  if(lColourG)then
  do i=nbounce0+1,nbounce6
    !dirichlet condition
    !anti-bounce-back approach
    !from page 200 Kruger's book "the lattice boltzmann method"
    !NOTE de[x,y,z]=zero eliminated
    
    myphi_CG(0:links)=phi_CG(0:links)+varphi_CG(0:links)*alpha_CG
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),1) = &
     -fluidsub(ibounce(1,i)+ex(1),ibounce(2,i),ibounce(3,i),2) &
     +myphi_CG(1)*TWO*rho_s(i) + &
     p(1)*TWO*rho_s(i)* &
     ((dex(1)*u_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),2) = &
     -fluidsub(ibounce(1,i)+ex(2),ibounce(2,i),ibounce(3,i),1) &
     +myphi_CG(2)*TWO*rho_s(i) +  &
     p(2)*TWO*rho_s(i)* &
     ((dex(2)*u_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
    
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),3) = &
     -fluidsub(ibounce(1,i),ibounce(2,i)+ey(3),ibounce(3,i),4) &
     +myphi_CG(3)*TWO*rho_s(i) +  &
     p(3)*TWO*rho_s(i)* &
     ((dey(3)*v_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),4) = &
     -fluidsub(ibounce(1,i),ibounce(2,i)+ey(4),ibounce(3,i),3) &
     +myphi_CG(4)*TWO*rho_s(i) +  &
     p(4)*TWO*rho_s(i)* &
     ((dey(4)*v_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
    
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),5) = &
     -fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i)+ez(5),6) &
     +myphi_CG(5)*TWO*rho_s(i) +  &
     p(5)*TWO*rho_s(i)* &
     ((dez(5)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),6) = &
     -fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i)+ez(6),5) &
     +myphi_CG(6)*TWO*rho_s(i) +  &
     p(6)*TWO*rho_s(i)* &
     ((dez(6)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
    
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),7) = &
     -fluidsub(ibounce(1,i)+ex(7),ibounce(2,i)+ey(7),ibounce(3,i),8) &
     +myphi_CG(7)*TWO*rho_s(i) +  &
     p(7)*TWO*rho_s(i)* &
     ((dex(7)*u_s(i)+ &
     dey(7)*v_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),8) = &
     -fluidsub(ibounce(1,i)+ex(8),ibounce(2,i)+ey(8),ibounce(3,i),7) &
     +myphi_CG(8)*TWO*rho_s(i) +  &
     p(8)*TWO*rho_s(i)* &
     ((dex(8)*u_s(i)+ &
     dey(8)*v_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
    
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),9) = &
     -fluidsub(ibounce(1,i)+ex(9),ibounce(2,i)+ey(9),ibounce(3,i),10) &
     +myphi_CG(9)*TWO*rho_s(i) +  &
     p(9)*TWO*rho_s(i)* &
     ((dex(9)*u_s(i)+ &
     dey(9)*v_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),10) = &
     -fluidsub(ibounce(1,i)+ex(10),ibounce(2,i)+ey(10),ibounce(3,i),9) &
     +myphi_CG(10)*TWO*rho_s(i) +  &
     p(10)*TWO*rho_s(i)* &
     ((dex(10)*u_s(i)+ &
     dey(10)*v_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
     
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),11) = &
     -fluidsub(ibounce(1,i)+ex(11),ibounce(2,i),ibounce(3,i)+ez(11),12) &
     +myphi_CG(11)*TWO*rho_s(i) +  &
     p(11)*TWO*rho_s(i)* &
     ((dex(11)*u_s(i)+ &
     dez(11)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),12) = &
     -fluidsub(ibounce(1,i)+ex(12),ibounce(2,i),ibounce(3,i)+ez(12),11) &
     +myphi_CG(12)*TWO*rho_s(i) +  &
     p(12)*TWO*rho_s(i)* &
     ((dex(12)*u_s(i)+ &
     dez(12)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),13) = &
     -fluidsub(ibounce(1,i)+ex(13),ibounce(2,i),ibounce(3,i)+ez(13),14) &
     +myphi_CG(13)*TWO*rho_s(i) +  &
     p(13)*TWO*rho_s(i)* &
     ((dex(13)*u_s(i)+ &
     dez(13)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),14) = &
     -fluidsub(ibounce(1,i)+ex(14),ibounce(2,i),ibounce(3,i)+ez(14),13) &
     +myphi_CG(14)*TWO*rho_s(i) +  &
     p(14)*TWO*rho_s(i)* &
     ((dex(14)*u_s(i)+ &
     dez(14)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
    
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),15) = &
     fluidsub(ibounce(1,i),ibounce(2,i)+ey(15),ibounce(3,i)+ez(15),16) &
     +myphi_CG(15)*TWO*rho_s(i) +  &
     p(15)*TWO*rho_s(i)* &
     ((dey(15)*v_s(i)+ &
     dez(15)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),16) = &
     -fluidsub(ibounce(1,i),ibounce(2,i)+ey(16),ibounce(3,i)+ez(16),15) &
     +myphi_CG(16)*TWO*rho_s(i) +  &
     p(16)*TWO*rho_s(i)* &
     ((dey(16)*v_s(i)+ &
     dez(16)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
     
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),17) = &
     -fluidsub(ibounce(1,i),ibounce(2,i)+ey(17),ibounce(3,i)+ez(17),18) &
     +myphi_CG(17)*TWO*rho_s(i) +  &
     p(17)*TWO*rho_s(i)* &
     ((dey(17)*v_s(i)+ &
     dez(17)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),18) = &
     -fluidsub(ibounce(1,i),ibounce(2,i)+ey(18),ibounce(3,i)+ez(18),17) &
     +myphi_CG(18)*TWO*rho_s(i) +  &
     p(18)*TWO*rho_s(i)* &
     ((dey(18)*v_s(i)+ &
     dez(18)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
    
  enddo
  else
   do i=nbounce0+1,nbounce6
    !dirichlet condition
    !anti-bounce-back approach
    !from page 200 Kruger's book "the lattice boltzmann method"
    !NOTE de[x,y,z]=zero eliminated
    
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),1) = &
     -fluidsub(ibounce(1,i)+ex(1),ibounce(2,i),ibounce(3,i),2) &
     +p(1)*TWO*rho_s(i)* &
     (ONE+(dex(1)*u_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),2) = &
     -fluidsub(ibounce(1,i)+ex(2),ibounce(2,i),ibounce(3,i),1) &
     +p(2)*TWO*rho_s(i)* &
     (ONE+(dex(2)*u_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
    
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),3) = &
     -fluidsub(ibounce(1,i),ibounce(2,i)+ey(3),ibounce(3,i),4) &
     +p(3)*TWO*rho_s(i)* &
     (ONE+(dey(3)*v_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),4) = &
     -fluidsub(ibounce(1,i),ibounce(2,i)+ey(4),ibounce(3,i),3) &
     +p(4)*TWO*rho_s(i)* &
     (ONE+(dey(4)*v_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
    
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),5) = &
     -fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i)+ez(5),6) &
     +p(5)*TWO*rho_s(i)* &
     (ONE+(dez(5)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),6) = &
     -fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i)+ez(6),5) &
     +p(6)*TWO*rho_s(i)* &
     (ONE+(dez(6)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
    
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),7) = &
     -fluidsub(ibounce(1,i)+ex(7),ibounce(2,i)+ey(7),ibounce(3,i),8) &
     +p(7)*TWO*rho_s(i)* &
     (ONE+(dex(7)*u_s(i)+ &
     dey(7)*v_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),8) = &
     -fluidsub(ibounce(1,i)+ex(8),ibounce(2,i)+ey(8),ibounce(3,i),7) &
     +p(8)*TWO*rho_s(i)* &
     (ONE+(dex(8)*u_s(i)+ &
     dey(8)*v_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
    
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),9) = &
     -fluidsub(ibounce(1,i)+ex(9),ibounce(2,i)+ey(9),ibounce(3,i),10) &
     +p(9)*TWO*rho_s(i)* &
     (ONE+(dex(9)*u_s(i)+ &
     dey(9)*v_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),10) = &
     -fluidsub(ibounce(1,i)+ex(10),ibounce(2,i)+ey(10),ibounce(3,i),9) &
     +p(10)*TWO*rho_s(i)* &
     (ONE+(dex(10)*u_s(i)+ &
     dey(10)*v_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
     
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),11) = &
     -fluidsub(ibounce(1,i)+ex(11),ibounce(2,i),ibounce(3,i)+ez(11),12) &
     +p(11)*TWO*rho_s(i)* &
     (ONE+(dex(11)*u_s(i)+ &
     dez(11)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),12) = &
     -fluidsub(ibounce(1,i)+ex(12),ibounce(2,i),ibounce(3,i)+ez(12),11) &
     +p(12)*TWO*rho_s(i)* &
     (ONE+(dex(12)*u_s(i)+ &
     dez(12)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),13) = &
     -fluidsub(ibounce(1,i)+ex(13),ibounce(2,i),ibounce(3,i)+ez(13),14) &
     +p(13)*TWO*rho_s(i)* &
     (ONE+(dex(13)*u_s(i)+ &
     dez(13)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),14) = &
     -fluidsub(ibounce(1,i)+ex(14),ibounce(2,i),ibounce(3,i)+ez(14),13) &
     +p(14)*TWO*rho_s(i)* &
     (ONE+(dex(14)*u_s(i)+ &
     dez(14)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
    
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),15) = &
     fluidsub(ibounce(1,i),ibounce(2,i)+ey(15),ibounce(3,i)+ez(15),16) &
     +p(15)*TWO*rho_s(i)* &
     (ONE+(dey(15)*v_s(i)+ &
     dez(15)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),16) = &
     -fluidsub(ibounce(1,i),ibounce(2,i)+ey(16),ibounce(3,i)+ez(16),15) &
     +p(16)*TWO*rho_s(i)* &
     (ONE+(dey(16)*v_s(i)+ &
     dez(16)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
     
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),17) = &
     -fluidsub(ibounce(1,i),ibounce(2,i)+ey(17),ibounce(3,i)+ez(17),18) &
     +p(17)*TWO*rho_s(i)* &
     (ONE+(dey(17)*v_s(i)+ &
     dez(17)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
     
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),18) = &
     -fluidsub(ibounce(1,i),ibounce(2,i)+ey(18),ibounce(3,i)+ez(18),17) &
     +p(18)*TWO*rho_s(i)* &
     (ONE+(dey(18)*v_s(i)+ &
     dez(18)*w_s(i))**TWO/cssq4 - &
     (u_s(i)**TWO+ &
     v_s(i)**TWO+ &
     w_s(i)**TWO)/cssq2)
    
   enddo
  endif
  

  
   do i=nbounce6+1,nbounce7
    !neumann condition
    !moving walls bounce-back approach
    !from page 180 Kruger's book "the lattice boltzmann method"
    !NOTE de[x,y,z]=zero eliminated
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),1) = &
     fluidsub(ibounce(1,i)+ex(1),ibounce(2,i),ibounce(3,i),2) &
     -p(2)*pref_bouzidi*rho_s(i)* &
     (dex(2)*u_s(i))
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),2) = &
     fluidsub(ibounce(1,i)+ex(2),ibounce(2,i),ibounce(3,i),1) &
     -p(1)*pref_bouzidi*rho_s(i)* &
     (dex(1)*u_s(i))
    
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),3) = &
     fluidsub(ibounce(1,i),ibounce(2,i)+ey(3),ibounce(3,i),4) &
     -p(4)*pref_bouzidi*rho_s(i)* &
     (dey(4)*v_s(i))
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),4) = &
     fluidsub(ibounce(1,i),ibounce(2,i)+ey(4),ibounce(3,i),3) &
     -p(3)*pref_bouzidi*rho_s(i)* &
     (dey(3)*v_s(i))
    
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),5) = &
     fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i)+ez(5),6) &
     -p(6)*pref_bouzidi*rho_s(i)* &
     (dez(6)*w_s(i))
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),6) = &
     fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i)+ez(6),5) &
     -p(5)*pref_bouzidi*rho_s(i)* &
     (dez(5)*w_s(i))
    
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),7) = &
     fluidsub(ibounce(1,i)+ex(7),ibounce(2,i)+ey(7),ibounce(3,i),8) &
     -p(8)*pref_bouzidi*rho_s(i)* &
     (dex(8)*u_s(i)+ &
     dey(8)*v_s(i))
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),8) = &
     fluidsub(ibounce(1,i)+ex(8),ibounce(2,i)+ey(8),ibounce(3,i),7) &
     -p(7)*pref_bouzidi*rho_s(i)* &
     (dex(7)*u_s(i)+ &
     dey(7)*v_s(i))
    
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),9) = &
     fluidsub(ibounce(1,i)+ex(9),ibounce(2,i)+ey(9),ibounce(3,i),10) &
     -p(10)*pref_bouzidi*rho_s(i)* &
     (dex(10)*u_s(i)+ &
     dey(10)*v_s(i))
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),10) = &
     fluidsub(ibounce(1,i)+ex(10),ibounce(2,i)+ey(10),ibounce(3,i),9) &
     -p(9)*pref_bouzidi*rho_s(i)* &
     (dex(9)*u_s(i)+ &
     dey(9)*v_s(i))
    
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),11) = &
     fluidsub(ibounce(1,i)+ex(11),ibounce(2,i),ibounce(3,i)+ez(11),12) &
     -p(12)*pref_bouzidi*rho_s(i)* &
     (dex(12)*u_s(i)+ &
     dez(12)*w_s(i))
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),12) = &
     fluidsub(ibounce(1,i)+ex(12),ibounce(2,i),ibounce(3,i)+ez(12),11) &
     -p(11)*pref_bouzidi*rho_s(i)* &
     (dex(11)*u_s(i)+ &
     dez(11)*w_s(i))
    
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),13) = &
     fluidsub(ibounce(1,i)+ex(13),ibounce(2,i),ibounce(3,i)+ez(13),14) &
     -p(14)*pref_bouzidi*rho_s(i)* &
     (dex(14)*u_s(i)+ &
     dez(14)*w_s(i))
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),14) = &
     fluidsub(ibounce(1,i)+ex(14),ibounce(2,i),ibounce(3,i)+ez(14),13) &
     -p(13)*pref_bouzidi*rho_s(i)* &
     (dex(13)*u_s(i)+ &
     dez(13)*w_s(i))
    
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),15) = &
     fluidsub(ibounce(1,i),ibounce(2,i)+ey(15),ibounce(3,i)+ez(15),16) &
     -p(16)*pref_bouzidi*rho_s(i)* &
     (dey(16)*v_s(i)+ &
     dez(16)*w_s(i))
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),16) = &
     fluidsub(ibounce(1,i),ibounce(2,i)+ey(16),ibounce(3,i)+ez(16),15) &
     -p(15)*pref_bouzidi*rho_s(i)* &
     (dey(15)*v_s(i)+ &
     dez(15)*w_s(i))
    
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),17) = &
     fluidsub(ibounce(1,i),ibounce(2,i)+ey(17),ibounce(3,i)+ez(17),18) &
     -p(18)*pref_bouzidi*rho_s(i)* &
     (dey(18)*v_s(i)+ &
     dez(18)*w_s(i))
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),18) = &
     fluidsub(ibounce(1,i),ibounce(2,i)+ey(18),ibounce(3,i)+ez(18),17) &
     -p(17)*pref_bouzidi*rho_s(i)* &
     (dey(17)*v_s(i)+ &
     dez(17)*w_s(i))
    
   enddo
  
  
  
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
        fluidsub(ibounce(1,ii),ibounce(2,ii),ibounce(3,ii),l)= &
         equil_popCG(l,bc_omega(ii),alpha_CG,rho_s(ii), &
         u_s(ii),v_s(ii),w_s(ii))
      enddo
    enddo
  else
  
  do i=nbounce7+1,nbounce8
    !robin conditions
    !equilibrium scheme
    !from page 191 Kruger's book "the lattice boltzmann method"
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),0)= &
     equil_pop00(rho_s(i),u_s(i),v_s(i),w_s(i))

    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),1)= &
     equil_pop01(rho_s(i),u_s(i),v_s(i),w_s(i))

    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),2)= &
     equil_pop02(rho_s(i),u_s(i),v_s(i),w_s(i))

    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),3)= &
     equil_pop03(rho_s(i),u_s(i),v_s(i),w_s(i))

    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),4)= &
     equil_pop04(rho_s(i),u_s(i),v_s(i),w_s(i))

    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),5)= &
     equil_pop05(rho_s(i),u_s(i),v_s(i),w_s(i))

    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),6)= &
     equil_pop06(rho_s(i),u_s(i),v_s(i),w_s(i))

    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),7)= &
     equil_pop07(rho_s(i),u_s(i),v_s(i),w_s(i))

    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),8)= &
     equil_pop08(rho_s(i),u_s(i),v_s(i),w_s(i))

    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),9)= &
     equil_pop09(rho_s(i),u_s(i),v_s(i),w_s(i))

    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),10)= &
     equil_pop10(rho_s(i),u_s(i),v_s(i),w_s(i))

    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),11)= &
     equil_pop11(rho_s(i),u_s(i),v_s(i),w_s(i))

    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),12)= &
     equil_pop12(rho_s(i),u_s(i),v_s(i),w_s(i))

    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),13)= &
     equil_pop13(rho_s(i),u_s(i),v_s(i),w_s(i))

    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),14)= &
     equil_pop14(rho_s(i),u_s(i),v_s(i),w_s(i))

    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),15)= &
     equil_pop15(rho_s(i),u_s(i),v_s(i),w_s(i))

    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),16)= &
     equil_pop16(rho_s(i),u_s(i),v_s(i),w_s(i))

    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),17)= &
     equil_pop17(rho_s(i),u_s(i),v_s(i),w_s(i))

    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),18)= &
     equil_pop18(rho_s(i),u_s(i),v_s(i),w_s(i))
     
  enddo
  
  endif
  
  
  

  
  do i=nbounce8+1,nbounce9
    !neumann condition
    !moving walls bounce-back approach
    !from page 180 Kruger's book "the lattice boltzmann method"
    !NOTE de[x,y,z]=zero eliminated
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),1) = &
     fluidsub(ibounce(1,i)+ex(1),ibounce(2,i),ibounce(3,i),2)- &
     p(2)*pref_bouzidi*rho_s(i)* &
     (dex(2)*u_s(i))
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),2) = &
     fluidsub(ibounce(1,i)+ex(2),ibounce(2,i),ibounce(3,i),1)- &
     p(1)*pref_bouzidi*rho_s(i)* &
     (dex(1)*u_s(i))
    
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),3) = &
     fluidsub(ibounce(1,i),ibounce(2,i)+ey(3),ibounce(3,i),4)- &
     p(4)*pref_bouzidi*rho_s(i)* &
     (dey(4)*v_s(i))
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),4) = &
     fluidsub(ibounce(1,i),ibounce(2,i)+ey(4),ibounce(3,i),3)- &
     p(3)*pref_bouzidi*rho_s(i)* &
     (dey(3)*v_s(i))
    
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),5) = &
     fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i)+ez(5),6)- &
     p(6)*pref_bouzidi*rho_s(i)* &
     (dez(6)*w_s(i))
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),6) = &
     fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i)+ez(6),5)- &
     p(5)*pref_bouzidi*rho_s(i)* &
     (dez(5)*w_s(i))
    
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),7) = &
     fluidsub(ibounce(1,i)+ex(7),ibounce(2,i)+ey(7),ibounce(3,i),8)- &
     p(8)*pref_bouzidi*rho_s(i)* &
     (dex(8)*u_s(i)+ &
     dey(8)*v_s(i))
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),8) = &
     fluidsub(ibounce(1,i)+ex(8),ibounce(2,i)+ey(8),ibounce(3,i),7)- &
     p(7)*pref_bouzidi*rho_s(i)* &
     (dex(7)*u_s(i)+ &
     dey(7)*v_s(i))
    
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),9) = &
     fluidsub(ibounce(1,i)+ex(9),ibounce(2,i)+ey(9),ibounce(3,i),10)- &
     p(10)*pref_bouzidi*rho_s(i)* &
     (dex(10)*u_s(i)+ &
     dey(10)*v_s(i))
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),10) = &
     fluidsub(ibounce(1,i)+ex(10),ibounce(2,i)+ey(10),ibounce(3,i),9)- &
     p(9)*pref_bouzidi*rho_s(i)* &
     (dex(9)*u_s(i)+ &
     dey(9)*v_s(i))
    
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),11) = &
     fluidsub(ibounce(1,i)+ex(11),ibounce(2,i),ibounce(3,i)+ez(11),12)- &
     p(12)*pref_bouzidi*rho_s(i)* &
     (dex(12)*u_s(i)+ &
     dez(12)*w_s(i))
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),12) = &
     fluidsub(ibounce(1,i)+ex(12),ibounce(2,i),ibounce(3,i)+ez(12),11)- &
     p(11)*pref_bouzidi*rho_s(i)* &
     (dex(11)*u_s(i)+ &
     dez(11)*w_s(i))
    
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),13) = &
     fluidsub(ibounce(1,i)+ex(13),ibounce(2,i),ibounce(3,i)+ez(13),14)- &
     p(14)*pref_bouzidi*rho_s(i)* &
     (dex(14)*u_s(i)+ &
     dez(14)*w_s(i))
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),14) = &
     fluidsub(ibounce(1,i)+ex(14),ibounce(2,i),ibounce(3,i)+ez(14),13)- &
     p(13)*pref_bouzidi*rho_s(i)* &
     (dex(13)*u_s(i)+ &
     dez(13)*w_s(i))
    
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),15) = &
     fluidsub(ibounce(1,i),ibounce(2,i)+ey(15),ibounce(3,i)+ez(15),16)- &
     p(16)*pref_bouzidi*rho_s(i)* &
     (dey(16)*v_s(i)+ &
     dez(16)*w_s(i))
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),16) = &
     fluidsub(ibounce(1,i),ibounce(2,i)+ey(16),ibounce(3,i)+ez(16),15)- &
     p(15)*pref_bouzidi*rho_s(i)* &
     (dey(15)*v_s(i)+ &
     dez(15)*w_s(i))
    
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),17) = &
     fluidsub(ibounce(1,i),ibounce(2,i)+ey(17),ibounce(3,i)+ez(17),18)- &
     p(18)*pref_bouzidi*rho_s(i)* &
     (dey(18)*v_s(i)+ &
     dez(18)*w_s(i))
    
    fluidsub(ibounce(1,i),ibounce(2,i),ibounce(3,i),18) = &
     fluidsub(ibounce(1,i),ibounce(2,i)+ey(18),ibounce(3,i)+ez(18),17)- &
     p(17)*pref_bouzidi*rho_s(i)* &
     (dey(17)*v_s(i)+ &
     dez(17)*w_s(i))
    
  enddo
  
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
 
 subroutine fill_new_isfluid_inlist(isub,jsub,ksub, nspheres, &
  spherelists, fvalue)
 
!***********************************************************************
!     
!     LBsoft subroutine to initialize the isfluid array values inside
!     the particle template
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     authors: F. Bonaccorso & M. Lauricella
!     last modification April 2020
!     
!*********************************************************************** 
 
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

 pure function linear(i,j,k)
 integer, intent(in) :: i,j,k      
 integer :: linear

 linear = i*100**2 + j*100 + k
end function linear

!******************START PART TO MANAGE THE COPY WALL*******************

 subroutine compute_densities_wall(nstep,lsingle_fluid,rhoR,rhoB)
 
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
  logical, intent(in) :: lsingle_fluid
  real(kind=PRC), allocatable, dimension(:,:,:)  :: rhoR,rhoB
  
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
 
 end module fluids_bc_mod
