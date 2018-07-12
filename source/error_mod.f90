 
#include <default_macro.h>
 module error_mod
 
!***********************************************************************
!     
!     LBsoft module containing subroutines which print warning and
!     close the software if an error is occurred
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification September 2017
!     
!***********************************************************************
 
 use version_mod, only : idrank,abort_world,finalize_world
 
 implicit none
 
 private
 
 public :: error
 public :: warning
 
 contains
 
 subroutine error(kode)
 
!***********************************************************************
!     
!     LBsoft subroutine for printing error banners and close the
!     program
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification September 2017
!     
!***********************************************************************
  
  implicit none
  
  integer,intent(in) :: kode
  
  integer,parameter :: outp=6
  character(len=*),parameter :: outf='(a)'
  character(len=*),parameter :: outf2='(2a)'
  
  if(idrank==0)then
     select case (kode)
      case (1)
        write(outp,outf)'ERROR - input file named input.dat not found.'
      case (2)
        write(outp,outf)'ERROR - input file opened with error.'
      case (3)
        write(outp,outf) &
        'ERROR - dimension of simulation box not found in input file!'
      case (4)
        write(outp,outf)'ERROR - bad allocation in subroutine allocate_fluids.'
      case (5)
        write(outp,outf)'ERROR - error in reading input file'
      case (6)
        write(outp,outf)'ERROR - unknown directive in input file.'
      case (7)
        write(outp,outf)'ERROR - incomplete input file.'
      case (8)
        write(outp,outf2)'ERROR - bad allocation in subroutine ', &
         'allocate_array_buffservice3d.'
      case (9)
        write(outp,outf)'ERROR - wrong boundary conditions read in input file.'
      case (10)
        write(outp,outf2)'ERROR - boundary conditions not found in input file! ', &
         'please set them in input.'
      case default
        write(outp,'(a,i18)')'unknown ERROR! code = ',kode
    end select
  endif
  
  call finalize_world()
  stop
  
 end subroutine error
 
 subroutine warning(kode,ddata,wstring)
 
!***********************************************************************
!     
!     LBsoft subroutine for printing warning banners
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification January 2017
!     
!***********************************************************************
  
  implicit none
 
  integer, intent(in) :: kode
  real(kind=PRC), intent(in), optional :: ddata
  character(len=*), intent(in), optional :: wstring
  
  integer,parameter :: outp=6
  character(len=*),parameter :: outf='(a)'
  character(len=10) :: r_char
  
  if(idrank/=0)return
  
  select case (kode)
    case (1)
      write (r_char,'(i10)')nint(ddata)
      write(outp,'(/,3a)')"WARNING - there is the following unknown directive at line ", &
       trim(adjustl(r_char))," in the input file:"
      write(outp,'(3a,/)')"'",trim(wstring),"'"
    case (2)
      write (r_char,'(i10)')nint(ddata)
      write(outp,'(/,3a,/)')"WARNING - the allocation number ", &
       trim(adjustl(r_char))," exits with error"
    case (3)
      write(outp,'(/,2a,/)')"WARNING - only full preiodic boundary ", &
      "conditions are actualy implemented"
    case (4)
      write(outp,'(/,2a,/)')"WARNING - viscosity or relaxation ", &
      "time tau have to be defined in input file"
    case default
      write(outp,'(/,a,i8,/)')"unknown WARNING! code = ",kode
  end select

  return

 end subroutine warning

 end module error_mod


