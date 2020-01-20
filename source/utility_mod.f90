 
#include <default_macro.h>
 module utility_mod
 
!***********************************************************************
!     
!     LBsoft module containing generic supporting subroutines called 
!     by different modules 
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
 use version_mod, only : idrank,or_world_larr
 use error_mod
 
 implicit none
 
 private
 
 integer, public, save :: nlbuffservice=0
 integer, public, save :: nibuffservice=0
 integer, public, save :: nbuffservice=0
 integer, public, save :: nbdf=0
 integer, public, save, protected :: nbuffservice3d=0
 integer, public, save :: nbuffservice_x=0
 integer, public, save :: nbuffservice_y=0
 integer, public, save :: nbuffservice_z=0
 logical, public, save :: linit_seed=.false.
 logical, public, save :: ltest_mode=.false.
 logical, public, allocatable, save :: lbuffservice(:)
 integer, public, allocatable, save :: ibuffservice(:)
 real(kind=PRC), public, allocatable, save :: buffservice(:)
 real(kind=PRC), public, allocatable, save :: xdf(:),ydf(:),zdf(:)
 real(kind=PRC), public, allocatable, save :: buffservice3d(:,:,:)
  
 real(kind=PRC), public, parameter :: & 
  Pi=real(3.141592653589793238462643383279502884d0,kind=PRC)

 real(kind=PRC), public, parameter :: conv_rad=Pi/real(180.d0,kind=PRC)
 
 public :: allocate_array_lbuffservice
 public :: allocate_array_ibuffservice
 public :: allocate_array_buffservice
 public :: allocate_array_bdf
 public :: allocate_array_buffservice3d
 public :: init_random_seed
 public :: gauss
 public :: modulvec
 public :: dot
 public :: cross
 public :: xcross
 public :: ycross
 public :: zcross
 public :: sig
 public :: write_fmtnumb
 public :: write_fmtnumb8
 public :: get_prntime
 public :: rand_noseeded
 public :: gauss_noseeded
 public :: dcell
 public :: invert
 public :: int_cube_sphere
 public :: openLogFile
 public :: fcut
 public :: space_fmtnumb
 public :: space_fmtnumb12
 public :: test_little_endian
 
 contains
 
 subroutine allocate_array_lbuffservice(imiomax)

!***********************************************************************
!     
!     LBsoft subroutine for reallocating the service array buff 
!     which is used within this module
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification july 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: imiomax
  
  if(nlbuffservice/=0)then
    if(imiomax>nlbuffservice)then
      deallocate(lbuffservice)
      nlbuffservice=imiomax+100
      allocate(lbuffservice(0:nlbuffservice))
    endif
  else
    nlbuffservice=imiomax+100
    allocate(lbuffservice(0:nlbuffservice))
  endif
  
  return
  
 end subroutine allocate_array_lbuffservice
 
 subroutine allocate_array_ibuffservice(imiomax)

!***********************************************************************
!     
!     LBsoft subroutine for reallocating the service array buff 
!     which is used within this module
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: imiomax
  
  if(nibuffservice/=0)then
    if(imiomax>nibuffservice)then
      deallocate(ibuffservice)
      nibuffservice=imiomax+100
      allocate(ibuffservice(0:nibuffservice))
    endif
  else
    nibuffservice=imiomax+100
    allocate(ibuffservice(0:nibuffservice))
  endif
  
  return
  
 end subroutine allocate_array_ibuffservice
 
 subroutine allocate_array_buffservice(imiomax)

!***********************************************************************
!     
!     LBsoft subroutine for reallocating the service array buff 
!     which is used within this module
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: imiomax
  
  if(nbuffservice/=0)then
    if(imiomax>nbuffservice)then
      deallocate(buffservice)
      nbuffservice=imiomax+100
      allocate(buffservice(0:nbuffservice))
    endif
  else
    nbuffservice=imiomax+100
    allocate(buffservice(0:nbuffservice))
  endif
  
  return
  
 end subroutine allocate_array_buffservice
 
 subroutine allocate_array_bdf(imiomax)

!***********************************************************************
!     
!     LBsoft subroutine for reallocating the service array buff 
!     which is used outside this module for particles
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: imiomax
  
  if(nbdf/=0)then
    if(imiomax>nbdf)then
      deallocate(xdf,ydf,zdf)
      nbdf=imiomax+100
      allocate(xdf(1:nbdf),ydf(1:nbdf),zdf(1:nbdf))
    endif
  else
    nbdf=imiomax+100
    allocate(xdf(1:nbdf),ydf(1:nbdf),zdf(1:nbdf))
  endif
  
  return
  
 end subroutine allocate_array_bdf
 
 subroutine allocate_array_buffservice3d(imiomin_x,imiomax_x,imiomin_y, &
  imiomax_y,imiomin_z,imiomax_z)

!***********************************************************************
!     
!     LBsoft subroutine for reallocating the service array buff 
!     which is used within this module
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: imiomin_x,imiomax_x,imiomin_y, &
   imiomax_y,imiomin_z,imiomax_z
   
  integer :: ierr
  logical, dimension(1) :: ltest=.false.
  
  if(allocated(buffservice3d))then
    if(imiomax_x>nbuffservice_x .or. imiomax_y>nbuffservice_y .or. &
     imiomax_z>nbuffservice_z)then
      deallocate(buffservice3d)
      nbuffservice_x=imiomax_x+100
      nbuffservice_y=imiomax_y+100
      nbuffservice_z=imiomax_z+100
      allocate(buffservice3d(imiomin_x:nbuffservice_x,imiomin_y:nbuffservice_y, &
       imiomin_z:nbuffservice_z),stat=ierr)
      ltest=.false.
      if(ierr.ne.0)then
        call warning(2,dble(1))
        ltest=.true.
      endif
    endif
  else
    nbuffservice_x=imiomax_x+100
    nbuffservice_y=imiomax_y+100
    nbuffservice_z=imiomax_z+100
    allocate(buffservice3d(imiomin_x:nbuffservice_x,imiomin_y:nbuffservice_y, &
     imiomin_z:nbuffservice_z),stat=ierr)
    ltest=.false.
    if(ierr.ne.0)then
      call warning(2,dble(1))
      ltest=.true.
    endif
  endif
  
  nbuffservice3d=(nbuffservice_x-imiomin_x+1)* &
   (nbuffservice_y-imiomin_y+1)*(nbuffservice_z-imiomin_z+1)
  
  call or_world_larr(ltest,1)
  if(ltest(1))call error(8)
  
  return
  
 end subroutine allocate_array_buffservice3d
 
 subroutine init_random_seed(myseed)
 
!***********************************************************************
!     
!     LBsoft subroutine for initialising the random generator
!     by the seed given in input file or by a random seed
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  integer,intent(in),optional :: myseed
  integer :: i, n, clock
  
  integer, allocatable :: seed(:)
          
  call random_seed(size = n)
  
  allocate(seed(n))
  
  if(present(myseed))then
!   If the seed is given in input
    seed = myseed*(idrank+1) + 37 * (/ (i - 1, i = 1, n) /)
    
  else
!   If the seed is not given in input it is generated by the clock
    call system_clock(count=clock)
         
    seed = clock*(idrank+1) + 37 * (/ (i - 1, i = 1, n) /)
    
  endif
  
  call random_seed(put = seed)
       
  deallocate(seed)
  
  return
 
 end subroutine init_random_seed
 
 pure function rand_noseeded(i,j,k,l)

!***********************************************************************
!     
!     LBsoft random number generator based on the universal
!     random number generator of marsaglia et Al.
!     it is delivered in a form which is depending on a seed squence 
!     provided at each call in order to ensure the MPI independence
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************

  implicit none
      
  integer, intent(in) :: i,j,k,l
  integer :: isub,jsub,ksub,lsub,msub
  integer ::ii,jj
  real(4) :: s,t,u33,u97,csub,uni
  real(kind=PRC) :: rand_noseeded
  
  real(4), parameter :: c =  362436.0/16777216.0
  real(4), parameter :: cd= 7654321.0/16777216.0
  real(4), parameter :: cm=16777213.0/16777216.0
  
! initial values of i,j,k must be in range 1 to 178 (not all 1)
! initial value of l must be in range 0 to 168.
      
  isub=mod(i,178)
  isub=isub+1
  
  jsub=mod(j,178)
  jsub=jsub+1
  
  ksub=mod(k,178)
  ksub=ksub+1
  
  lsub=mod(l,169)
  
! initialization on fly
  
  ii=97
  s=0.0
  t=0.5
  do jj=1,24
    msub=mod(mod(isub*jsub,179)*ksub,179)
    isub=jsub
    jsub=ksub
    ksub=msub
    lsub=mod(53*lsub+1,169)
    if(mod(lsub*msub,64).ge.32)s=s+t
    t=0.5*t
  enddo
  u97=s
  
  ii=33
  s=0.0
  t=0.5
  do jj=1,24
    msub=mod(mod(isub*jsub,179)*ksub,179)
    isub=jsub
    jsub=ksub
    ksub=msub
    lsub=mod(53*lsub+1,169)
    if(mod(lsub*msub,64).ge.32)s=s+t
    t=0.5*t
  enddo
  u33=s
  uni=u97-u33
  if(uni.lt.0.0)uni=uni+1.0
  csub=c-cd
  if(csub.lt.0.0)csub=csub+cm
  uni=uni-csub
  if(uni.lt.0.0)uni=uni+1.0
  rand_noseeded=real(uni,kind=PRC)

  return
  
 end function rand_noseeded
 
 function gauss()
 
!***********************************************************************
!     
!     LBsoft subroutine for generating random number normally
!     distributed by the Box-Muller transformation
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC) :: gauss
  real(kind=PRC) :: dtemp1,dtemp2
  logical :: lredo
  
  call random_number(dtemp1)
  call random_number(dtemp2)
  
  lredo=.true.
  
! the number is extract again if it is nan
  do while(lredo)
    lredo=.false.
!   Box-Muller transformation
    gauss=sqrt(- TWO *log(dtemp1))*cos(2*pi*dtemp2)
    if(isnan(cos(gauss)))lredo=.true.
  enddo
  
 end function gauss
 
 pure function gauss_noseeded(i,j,k,l)
 
!***********************************************************************
!     
!     LBsoft subroutine for generating random number normally
!     distributed by the Box-Muller transformation without a seed
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  integer, intent(in) :: i,j,k,l
  integer :: kk
  real(kind=PRC) :: gauss_noseeded
  real(kind=PRC) :: dtemp1,dtemp2
  logical :: lredo
  
  lredo=.true.
  kk=0
! the number is extract again if it is nan
  do while(lredo)
    lredo=.false.
    dtemp1=rand_noseeded(i,j,k,l+kk)
    dtemp2=rand_noseeded(i,j,k,l+1+kk)
!   Box-Muller transformation
    gauss_noseeded=sqrt(- TWO *log(dtemp1))*cos(2*pi*dtemp2)
    if(isnan(cos(gauss_noseeded)))lredo=.true.
    kk=kk+1
  enddo
  
 end function gauss_noseeded
 
 pure function modulvec(a)
 
!***********************************************************************
!     
!     LBsoft function for computing the module of a vector
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC) :: modulvec
  real(kind=PRC), dimension(3), intent(in) :: a

  modulvec = dsqrt(a(1)** TWO + a(2)** TWO + a(3)** TWO )
  
  return
  
 end function modulvec
 
 pure function cross(a,b)
 
!***********************************************************************
!     
!     LBsoft function for computing the cross product of two vectors
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC), dimension(3) :: cross
  real(kind=PRC), dimension(3), intent(in) :: a, b

  cross(1) = a(2) * b(3) - a(3) * b(2)
  cross(2) = a(3) * b(1) - a(1) * b(3)
  cross(3) = a(1) * b(2) - a(2) * b(1)
  
  return
  
 end function cross
 
 pure function xcross(a,b)
 
!***********************************************************************
!     
!     LBsoft function for computing the cross product of two vectors
!     along x
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC) :: xcross
  real(kind=PRC), dimension(3), intent(in) :: a, b

  xcross = a(2) * b(3) - a(3) * b(2)
  
  return
  
 end function xcross
 
 pure function ycross(a,b)
 
!***********************************************************************
!     
!     LBsoft function for computing the cross product of two vectors
!     along y
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC) :: ycross
  real(kind=PRC), dimension(3), intent(in) :: a, b
  
  ycross = a(3) * b(1) - a(1) * b(3)
  
  return
  
 end function ycross
 
 pure function zcross(a,b)
 
!***********************************************************************
!     
!     LBsoft function for computing the cross product of two vectors
!     along z
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC) :: zcross
  real(kind=PRC), dimension(3), intent(in) :: a, b
  
  zcross = a(1) * b(2) - a(2) * b(1)
  
  return
  
 end function zcross
 
 pure function int_cube_sphere(rdim,xo,yo,zo,side,xc,yc,zc,ires)
 
!***********************************************************************
!     
!     LBsoft function for computing the volume intersection
!     between a cube and a sphere
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC) :: int_cube_sphere
  real(kind=PRC), intent(in) :: rdim,xo,yo,zo
  real(kind=PRC), intent(in) :: side,xc,yc,zc
  integer, intent(in) :: ires
  
  real(kind=PRC) :: myres,sqrdim,xmin,ymin,zmin,xx,yy,zz,rdist,mysum
  integer :: i,j,k,l
  
  sqrdim=rdim**TWO
  
  myres=side/real(ires,kind=PRC)
  xmin=xc+(myres-side)/TWO
  ymin=yc+(myres-side)/TWO
  zmin=zc+(myres-side)/TWO
  mysum=ZERO
  do k=1,ires
    zz=zmin+real(k-1,kind=PRC)*myres
    do j=1,ires
      yy=ymin+real(j-1,kind=PRC)*myres
      do i=1,ires
        xx=xmin+real(i-1,kind=PRC)*myres
        rdist=(xx-xo)**TWO+(yy-yo)**TWO+(zz-zo)**TWO
        if(rdist>sqrdim)then
          mysum=mysum+ONE
        endif
      enddo
    enddo
  enddo
  
  int_cube_sphere=mysum/real(ires*ires*ires,kind=PRC)
  
  return
  
 end function int_cube_sphere
 
 pure function dot(a,b)
 
!***********************************************************************
!     
!     LBsoft function for computing the dot product of two vectors
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  real(kind=PRC) :: dot
  real(kind=PRC), dimension(3), intent(in) :: a, b
  
  dot=dot_product(a,b)
  
  return
  
 end function dot
 
 pure function sig(num)
 
!***********************************************************************
!     
!     LBsoft function for returning the sign of a floating number
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in):: num
  
  real(kind=PRC) :: sig
  
  if(num>0)then
    sig= ONE
  elseif(num==0)then
    sig= ZERO
  else
    sig= - ONE
  endif
  
  return
 
 end function sig
 
 function dimenumb(inum)
 
!***********************************************************************
!     
!     LBsoft function for returning the number of digits
!     of an integer number
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************

  implicit none

  integer,intent(in) :: inum
  integer :: dimenumb
  integer :: i
  real(kind=PRC) :: tmp

  i=1
  tmp=real(inum,kind=PRC)
  do 
    if(tmp< TEN )exit
    i=i+1
    tmp=tmp/ TEN
  enddo

  dimenumb=i

  return

 end function dimenumb

 function write_fmtnumb(inum)
 
!***********************************************************************
!     
!     LBsoft function for returning the string of six characters 
!     with integer digits and leading zeros to the left
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none

  integer,intent(in) :: inum
  character(len=6) :: write_fmtnumb
  integer :: numdigit,irest
  real(kind=PRC) :: tmp
  character(len=22) :: cnumberlabel

  numdigit=dimenumb(inum)
  irest=6-numdigit
  if(irest>0)then
    write(cnumberlabel,"(a,i8,a,i8,a)")"(a",irest,",i",numdigit,")"
    write(write_fmtnumb,fmt=cnumberlabel)repeat('0',irest),inum
  else
    write(cnumberlabel,"(a,i8,a)")"(i",numdigit,")"
    write(write_fmtnumb,fmt=cnumberlabel)inum
  endif
  
  return

 end function write_fmtnumb
 
 function write_fmtnumb8(inum)
 
!***********************************************************************
!     
!     LBsoft function for returning the string of eight characters 
!     with integer digits and leading zeros to the left
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none

  integer,intent(in) :: inum
  character(len=8) :: write_fmtnumb8
  integer :: numdigit,irest
  real(kind=PRC) :: tmp
  character(len=22) :: cnumberlabel

  numdigit=dimenumb(inum)
  irest=8-numdigit
  if(irest>0)then
    write(cnumberlabel,"(a,i8,a,i8,a)")"(a",irest,",i",numdigit,")"
    write(write_fmtnumb8,fmt=cnumberlabel)repeat('0',irest),inum
  else
    write(cnumberlabel,"(a,i8,a)")"(i",numdigit,")"
    write(write_fmtnumb8,fmt=cnumberlabel)inum
  endif
  
  return

 end function write_fmtnumb8
 
 subroutine get_prntime(hms,timelp,prntim)
  
!***********************************************************************
!     
!     LBsoft subroutine for casting cpu elapsed time into days, hours,
!     minutes and seconds for printing (input timelp in seconds)
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  character(len=1), intent(out) :: hms
  real(kind=PRC), intent(in) :: timelp
  real(kind=PRC), intent(out) :: prntim
  
  if(timelp.ge. real(8.64d4,kind=PRC))then
    hms='d'
    prntim=timelp/real(8.64d4,kind=PRC)
  elseif(timelp.ge.real(3.6d3,kind=PRC))then
    hms='h'
    prntim=timelp/real(3.6d3,kind=PRC)
  elseif(timelp.ge.real(6.0d1,kind=PRC))then
    hms='m'
    prntim=timelp/real(6.0d1,kind=PRC)
  else
    hms='s'
    prntim=timelp
  endif
  
  return
  
 end subroutine get_prntime
 
 subroutine dcell(aaa,bbb)
      
!***********************************************************************
!     
!     LBsoft subroutine to calculate the dimensional properties of
!     a simulation cell specified by the input matrix aaa.
!     the results are returned in the array bbb, with :
!     
!     bbb(1 to 3) - lengths of cell vectors
!     bbb(4 to 6) - cosines of cell angles
!     bbb(7 to 9) - perpendicular cell widths
!     bbb(10)     - cell volume
!     
!     modified from DL Protein
!     original copyright - daresbury laboratory 1992
!     original author    - w. smith       april 1992
!     
!***********************************************************************
      
  real(kind=PRC), intent(in) :: aaa(9)
  real(kind=PRC), intent(out) :: bbb(10)
  real(kind=PRC) :: axb1,axb2,axb3,bxc1,bxc2,bxc3,cxa1,cxa2,cxa3
      
! calculate lengths of cell vectors
      
  bbb(1)=sqrt(aaa(1)*aaa(1)+aaa(2)*aaa(2)+aaa(3)*aaa(3))
  bbb(2)=sqrt(aaa(4)*aaa(4)+aaa(5)*aaa(5)+aaa(6)*aaa(6))
  bbb(3)=sqrt(aaa(7)*aaa(7)+aaa(8)*aaa(8)+aaa(9)*aaa(9))
      
! calculate cosines of cell angles
      
  bbb(4)=(aaa(1)*aaa(4)+aaa(2)*aaa(5)+aaa(3)*aaa(6))/(bbb(1)*bbb(2))
  bbb(5)=(aaa(1)*aaa(7)+aaa(2)*aaa(8)+aaa(3)*aaa(9))/(bbb(1)*bbb(3))
  bbb(6)=(aaa(4)*aaa(7)+aaa(5)*aaa(8)+aaa(6)*aaa(9))/(bbb(2)*bbb(3))
      
! calculate vector products of cell vectors
      
  axb1=aaa(2)*aaa(6)-aaa(3)*aaa(5)
  axb2=aaa(3)*aaa(4)-aaa(1)*aaa(6)
  axb3=aaa(1)*aaa(5)-aaa(2)*aaa(4)
  bxc1=aaa(5)*aaa(9)-aaa(6)*aaa(8)
  bxc2=aaa(6)*aaa(7)-aaa(4)*aaa(9)
  bxc3=aaa(4)*aaa(8)-aaa(5)*aaa(7)
  cxa1=aaa(8)*aaa(3)-aaa(2)*aaa(9)
  cxa2=aaa(1)*aaa(9)-aaa(3)*aaa(7)
  cxa3=aaa(2)*aaa(7)-aaa(1)*aaa(8)
      
! calculate volume of cell
      
  bbb(10)=abs(aaa(1)*bxc1+aaa(2)*bxc2+aaa(3)*bxc3)
      
!  calculate cell perpendicular widths
  bbb(7)=bbb(10)/sqrt(bxc1*bxc1+bxc2*bxc2+bxc3*bxc3)
  bbb(8)=bbb(10)/sqrt(cxa1*cxa1+cxa2*cxa2+cxa3*cxa3)
  bbb(9)=bbb(10)/sqrt(axb1*axb1+axb2*axb2+axb3*axb3)
      
  return
   
 end subroutine dcell
   
 subroutine invert(a,b,d)

!***********************************************************************
!     
!     LBsoft subroutine to invert a 3 * 3 matrix using cofactors
!     
!     modified from DL Protein
!     original copyright - daresbury laboratory 1992
!     original author    - w. smith       april 1992
!     
!***********************************************************************

  implicit none
  
  real(kind=PRC), intent(in) :: a(9)

  real(kind=PRC), intent(out) :: b(9),d
  real(kind=PRC) :: r


! calculate adjoint matrix
  b(1)=a(5)*a(9)-a(6)*a(8)
  b(2)=a(3)*a(8)-a(2)*a(9)
  b(3)=a(2)*a(6)-a(3)*a(5)
  b(4)=a(6)*a(7)-a(4)*a(9)
  b(5)=a(1)*a(9)-a(3)*a(7)
  b(6)=a(3)*a(4)-a(1)*a(6)
  b(7)=a(4)*a(8)-a(5)*a(7)
  b(8)=a(2)*a(7)-a(1)*a(8)
  b(9)=a(1)*a(5)-a(2)*a(4)

! calculate determinant
  d=a(1)*b(1)+a(4)*b(2)+a(7)*b(3)
  r=ZERO
  if(abs(d).gt.ZERO)r=ONE/d

! complete inverse matrix
  b(1)=r*b(1)
  b(2)=r*b(2)
  b(3)=r*b(3)
  b(4)=r*b(4)
  b(5)=r*b(5)
  b(6)=r*b(6)
  b(7)=r*b(7)
  b(8)=r*b(8)
  b(9)=r*b(9)

  return
      
 end subroutine invert
 
 subroutine openLogFile(nstep, msg, iounit)
  implicit none
  integer, intent(in) :: nstep,iounit
  character(len=*), intent(in) :: msg
  character(len=120) :: mynamefile


  mynamefile=repeat(' ',120)
  mynamefile=trim(msg)// '.iter'//write_fmtnumb(nstep)//'.'//write_fmtnumb(idrank)//'.dat'
  open(unit=iounit, file=trim(mynamefile), status='replace')
 end subroutine openLogFile

 pure function fcut(r,inner_cut,outer_cut)

!***********************************************************************
!
!     LBsoft function for fading an observable (r) within a given
!     interval
!
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification January 2018
!
!***********************************************************************

  implicit none

  real(kind=PRC), intent(in) :: r,inner_cut,outer_cut
  real(kind=PRC) :: fcut

  if ( r <= inner_cut ) then
    fcut = ONE
  elseif ( r > outer_cut ) then
      fcut = ZERO
  else
      fcut = HALF*cos((r-inner_cut)*Pi/(outer_cut-inner_cut))+HALF
  endif

  return

 end function fcut
 
 function space_fmtnumb(inum)
 
!***********************************************************************
!     
!     LBsoft function for returning the string of six characters 
!     with integer digits and leading spaces to the left
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none

  integer,intent(in) :: inum
  character(len=6) :: space_fmtnumb
  integer :: numdigit,irest
  real(kind=8) :: tmp
  character(len=22) :: cnumberlabel

  numdigit=dimenumb(inum)
  irest=6-numdigit
  if(irest>0)then
    write(cnumberlabel,"(a,i8,a,i8,a)")"(a",irest,",i",numdigit,")"
    write(space_fmtnumb,fmt=cnumberlabel)repeat(' ',irest),inum
  else
    write(cnumberlabel,"(a,i8,a)")"(i",numdigit,")"
    write(space_fmtnumb,fmt=cnumberlabel)inum
  endif
  
  return

 end function space_fmtnumb
 
 function space_fmtnumb12(inum)
 
!***********************************************************************
!     
!     LBsoft function for returning the string of six characters 
!     with integer digits and leading TWELVE spaces to the left
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none

  integer,intent(in) :: inum
  character(len=12) :: space_fmtnumb12
  integer :: numdigit,irest
  real(kind=8) :: tmp
  character(len=22) :: cnumberlabel

  numdigit=dimenumb(inum)
  irest=12-numdigit
  if(irest>0)then
    write(cnumberlabel,"(a,i8,a,i8,a)")"(a",irest,",i",numdigit,")"
    write(space_fmtnumb12,fmt=cnumberlabel)repeat(' ',irest),inum
  else
    write(cnumberlabel,"(a,i8,a)")"(i",numdigit,")"
    write(space_fmtnumb12,fmt=cnumberlabel)inum
  endif
  
  return

 end function space_fmtnumb12
 
 subroutine test_little_endian(ltest)
 
!***********************************************************************
!     
!     LBsoft subroutine for checking if the computing architecture
!     is working in little-endian or big-endian
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none 
  integer, parameter :: ik1 = selected_int_kind(2) 
  integer, parameter :: ik4 = selected_int_kind(9) 
   
  logical, intent(out) :: ltest
   
  if(btest(transfer(int((/1,0,0,0/),ik1),1_ik4),0)) then 
    !it is little endian
    ltest=.true.
  else 
    !it is big endian
    ltest=.false.
  end if 
   
  return
   
 end subroutine test_little_endian 

 end module utility_mod
