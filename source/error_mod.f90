 
#include <default_macro.h>
 module error_mod
 
!***********************************************************************
!     
!     JETSPIN module containing subroutines which print warning and
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
!     JETSPIN subroutine for printing error banners and close the
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
        write(outp,outf)'ERROR - ftype is wrong for this integrator.'
      case (5)
        write(outp,outf)'ERROR - error in reading input file'
        write(outp,outf2)'ERROR - none FINISH directive found ', &
         'in the input file.'
      case (6)
        write(outp,outf)'ERROR - unknown directive in input file.'
      case (7)
        write(outp,outf)'ERROR - incomplete input file.'
      case (8)
        write(outp,outf2)'ERROR - the resolution of jet ', &
         'discretization is too small.'
      case (9)
        write(outp,outf2)'ERROR - wrong selection of print dat ', &
        'style in input file.'
      case (10)
        write(outp,outf)'ERROR - '
      case (11)
        write(outp,outf)'ERROR - '
      case (12)
        write(outp,outf2)'ERROR - requested a spline fitting beyond ', &
         'the spline boundaries.'
      case (13)
        write(outp,outf2)'ERROR - in allocating jetptc in ', &
         'allocate_array_jetptc_mod.'
      case (14)
        write(outp,outf2)'ERROR - numerical instability.', &
         ' Please check the time step.'
      case (15)
        write(outp,outf)'ERROR - restart.dat file not found!'
      case (16)
        write(outp,outf2)'ERROR - dynamic refinement threshold too small!', &
         ' Please increase the dynamic refinement threshold.'
      case (17)
        write(outp,outf2)'ERROR - multiple step error', &
         ' Please check the input file.'
      case (18)
        write(outp,outf)'ERROR - the restart file is corrupted'
      case (19)
        write(outp,outf2)'ERROR - array yve not found in subroutine', &
         ' compute_coulomelec_driver!'
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
      write(outp,'(/,a,/,a,g20.10,/)')"WARNING - the 'resolution' value is not a submultiple of the 'initial length'", &
      "WARNING - the remaider is equal to",ddata
    case default
      write(outp,'(/,a,i8,/)')"unknown WARNING! code = ",kode
  end select

  return

 end subroutine warning

 end module error_mod


