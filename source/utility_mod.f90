 
#include <default_macro.h>
 module utility_mod
 
!***********************************************************************
!     
!     LBsoft module containing generic supporting subroutines called 
!     by different modules 
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
 use version_mod, only : idrank
 
 implicit none
 
 private
 
 integer, public, save :: nlbuffservice=0
 integer, public, save :: nibuffservice=0
 integer, public, save :: nbuffservice=0
 integer, public, save :: nbuffservice_x=0
 integer, public, save :: nbuffservice_y=0
 integer, public, save :: nbuffservice_z=0
 logical, public, allocatable, save :: lbuffservice(:)
 integer, public, allocatable, save :: ibuffservice(:)
 real(kind=PRC), public, allocatable, save :: buffservice(:)
 real(kind=PRC), public, allocatable, save :: buffservice3d(:,:,:)
  
 real(kind=PRC), public, parameter :: & 
  Pi=real(3.141592653589793238462643383279502884d0,kind=PRC)
 
 public :: allocate_array_lbuffservice
 public :: allocate_array_ibuffservice
 public :: allocate_array_buffservice
 public :: allocate_array_buffservice3d
 public :: init_random_seed
 public :: gauss
 public :: modulvec
 public :: dot
 public :: cross
 public :: sig
 public :: write_fmtnumb
 public :: get_prntime
 
 contains
 
 subroutine allocate_array_lbuffservice(imiomax)

!***********************************************************************
!     
!     LBsoft subroutine for reallocating the service array buff 
!     which is used within this module
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
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
!     licensed under Open Software License v. 3.0 (OSL-3.0)
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
!     licensed under Open Software License v. 3.0 (OSL-3.0)
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
 
 subroutine allocate_array_buffservice3d(imiomax_x,imiomax_y,imiomax_z)

!***********************************************************************
!     
!     LBsoft subroutine for reallocating the service array buff 
!     which is used within this module
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: imiomax_x,imiomax_y,imiomax_z
  
  if(allocated(buffservice3d))then
    if(imiomax_x>nbuffservice_x .or. imiomax_y>nbuffservice_y .or. &
     imiomax_z>nbuffservice_z)then
      deallocate(buffservice3d)
      nbuffservice_x=imiomax_x+100
      nbuffservice_y=imiomax_y+100
      nbuffservice_z=imiomax_z+100
      allocate(buffservice3d(0:nbuffservice_x,0:nbuffservice_y, &
       0:nbuffservice_z))
    endif
  else
    nbuffservice_x=imiomax_x+100
    nbuffservice_y=imiomax_y+100
    nbuffservice_z=imiomax_z+100
    allocate(buffservice3d(0:nbuffservice_x,0:nbuffservice_y, &
     0:nbuffservice_z))
  endif
  
  return
  
 end subroutine allocate_array_buffservice3d
 
 subroutine init_random_seed(myseed)
 
!***********************************************************************
!     
!     LBsoft subroutine for initialising the random generator
!     by the seed given in input file or by a random seed
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
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
 
 function gauss()
 
!***********************************************************************
!     
!     LBsoft subroutine for generating random number normally
!     distributed by the Box-Muller transformation
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
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
 
 pure function modulvec(a)
 
!***********************************************************************
!     
!     LBsoft function for computing the module of a vector
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
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
!     licensed under Open Software License v. 3.0 (OSL-3.0)
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
 
 pure function dot(a,b)
 
!***********************************************************************
!     
!     LBsoft function for computing the dot product of two vectors
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
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
!     licensed under Open Software License v. 3.0 (OSL-3.0)
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
!     licensed under Open Software License v. 3.0 (OSL-3.0)
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
!     licensed under Open Software License v. 3.0 (OSL-3.0)
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
 
 subroutine get_prntime(hms,timelp,prntim)
  
!***********************************************************************
!     
!     LBsoft subroutine for casting cpu elapsed time into days, hours,
!     minutes and seconds for printing (input timelp in seconds)
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
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
 
 end module utility_mod


