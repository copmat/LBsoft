 
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
      case (11)
        write(outp,outf)'ERROR - wrong vtk subroutine selected! only SERIAL'
      case (12)
        write(outp,outf2)'ERROR - the actual boundary conditions are not supported', &
         ' with the current lattice scheme'
      case (13)
        write(outp,outf2)'ERROR - the actual OPEN boundary conditions are not', &
         ' implemented. Check the bc_type value selected'
      case (14)
        write(outp,outf)'ERROR - wrong LB integrator selected'
      case (15)
        write(outp,outf2)'ERROR - bad allocation in subroutine ', &
         'initialize_isfluid_bcfluid.'
      case (16)
        write(outp,outf2)'ERROR - in subroutine an open boundary node ', &
         'was set without a direction.'
      case (17)
        write(outp,outf)'ERROR - input file named input.xyz not found.'
      case (18)
        write(outp,outf)'ERROR - input.xyz file opened with error.'
      case (19)
        write(outp,outf)'ERROR - error in reading the input file input.xyz.'
      case (20)
        write(outp,outf)'ERROR - bad allocation in subroutine allocate_particles.'
      case (21)
        write(outp,outf)'ERROR - densvar is too low. Increase it in input file.'
      case (22)
        write(outp,outf)'ERROR - the quaternion check exits with error.'
      case (23)
        write(outp,outf)'ERROR - densvar is too low.'
      case (24)
        write(outp,outf)'ERROR - vertest array too small.'
      case (25)
        write(outp,outf)'ERROR - fail in allocate the vertest array.'
      case (26)
        write(outp,outf)'ERROR - fail in allocate the prmvdw array.'
      case (27)
        write(outp,outf)'ERROR - wrong ltpvdw in compute_inter_force'
      case (28)
        write(outp,outf)'ERROR - bad allocation in subroutine nve_lf.'
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
      write(outp,'(/,2a,/)')"WARNING - the requested boundary ", &
      "conditions is NOT actualy implemented"
    case (4)
      write(outp,'(/,2a,/)')"WARNING - viscosity or relaxation ", &
      "time tau have to be defined in input file"
    case (5)
      write(outp,'(/,2a,/)')"WARNING - nsteps is not defined in ", &
      "input file"
    case (6)
      write(outp,'(/,a)') &
      "WARNING - possible [keys] to be used are reported in the following table"
      write(outp,'(a)') &
      "WARNING - each [key] should be separated by a space character (e.g. t x y z rc)"
      write(outp,'(a)') &
      "WARNING - ********************************************************************************"
      write(outp,'(a)') &
      "WARNING - * [keys]               * [meanings]                                            *"
      write(outp,'(a)') &
      "WARNING - ********************************************************************************"
      write(outp,'(a)') &
      "WARNING - * t                    * time in lb units                                      *"
      write(outp,'(a)') &
      "WARNING - * dens1                * mean fluid density of first component in lb units     *"
      write(outp,'(a)') &
      "WARNING - * dens2                * mean fluid density of second component in lb units    *"
      write(outp,'(a)') &
      "WARNING - * maxd1                * max fluid density of first component in lb units      *"
      write(outp,'(a)') &
      "WARNING - * maxd2                * max fluid density of second component in lb units     *"
      write(outp,'(a)') &
      "WARNING - * mind1                * min fluid density of first component in lb units      *"
      write(outp,'(a)') &
      "WARNING - * mind2                * min fluid density of second component in lb units     *"
      write(outp,'(a)') &
      "WARNING - * maxvx                * max fluid velocity of fluid mix along x in lb units   *"
      write(outp,'(a)') &
      "WARNING - * minvx                * min fluid velocity of fluid mix along x in lb units   *"
      write(outp,'(a)') &
      "WARNING - * maxvy                * max fluid velocity of fluid mix along y in lb units   *"
      write(outp,'(a)') &
      "WARNING - * minvy                * min fluid velocity of fluid mix along y in lb units   *"
      write(outp,'(a)') &
      "WARNING - * maxvz                * max fluid velocity of fluid mix along z in lb units   *"
      write(outp,'(a)') &
      "WARNING - * minvz                * min fluid velocity of fluid mix along z in lb units   *"
      write(outp,'(a)') &
      "WARNING - * cpu                  * time for every print interval in seconds              *"
      write(outp,'(a)') &
      "WARNING - * cpur                 * remaining time to the end in seconds                  *"
      write(outp,'(a)') &
      "WARNING - * cpue                 * elapsed time in seconds                               *"
      write(outp,'(a,/)') &
      "WARNING - ********************************************************************************"
    case (7)
      write(outp,'(/,a)')"WARNING - 'print list' is not correctly specified"
      write(outp,'(a,/)')"WARNING - 'print list' should be specified as 'print list [keys]'"
    case (8)
      write (r_char,'(i10)')nint(ddata)
      write(outp,'(/,a)')"WARNING - the number fluid components should be 1 or 2."
       write(outp,'(2a,/)')"WARNING - the actual number of fluids is : ", &
       trim(adjustl(r_char))
    case (9)
        write(outp,'(/,2a)')"WARNING - the actual boundary conditions are not supported", &
         " with the compiled lattice scheme"
        write(outp,'(2a,/)')"WARNING - the actual lattice scheme is : ",trim(wstring)
    case (10)
        write(outp,'(/,2a,/)')"WARNING - the shanchen pair force cannot be", &
         " applyied with a single fluid component!"
    case (11)
      write (r_char,'(i10)')nint(ddata)
      write(outp,'(/,3a,/)')"WARNING - there is an error in at line ", &
       trim(adjustl(r_char))," in the input.xyz file:"
    case (12)
      write(outp,'(/,a)')"WARNING - 'read list' is not correctly specified"
      write(outp,'(a,/)')"WARNING - 'read list' should be specified as 'read list [keys]'"
    case (13)
      write(outp,'(/,a)') &
      "WARNING - possible [keys] to be used are reported in the following table"
      write(outp,'(a)') &
      "WARNING - each [key] should be separated by a space character (e.g. mass vx vy vz rad)"
      write(outp,'(a)') &
      "WARNING - ********************************************************************************"
      write(outp,'(a)') &
      "WARNING - * [keys]               * [meanings]                                            *"
      write(outp,'(a)') &
      "WARNING - ********************************************************************************"
!      write(outp,'(a)') &
!      "WARNING - * char                 * charge of the particle in lb units                    *"
      write(outp,'(a)') &
      "WARNING - * mass                 * mass of the particle in lb units                      *"
      write(outp,'(a)') &
      "WARNING - * vx                   * particle velocity along x in lb units                 *"
      write(outp,'(a)') &
      "WARNING - * vy                   * particle velocity along y in lb units                 *"
      write(outp,'(a)') &
      "WARNING - * vz                   * particle velocity along z in lb units                 *"
      write(outp,'(a)') &
      "WARNING - * ox                   * particle angular velocity along x in lb units         *"
      write(outp,'(a)') &
      "WARNING - * oy                   * particle angular velocity along y in lb units         *"
      write(outp,'(a)') &
      "WARNING - * oz                   * particle angular velocity along z in lb units         *"
      write(outp,'(a)') &
      "WARNING - * rad                  * radius of the spherical particle in lb units          *"
      write(outp,'(a)') &
      "WARNING - * radx                 * radius of the particle along x in lb units            *"
      write(outp,'(a)') &
      "WARNING - * rady                 * radius of the particle along x in lb units            *"
      write(outp,'(a)') &
      "WARNING - * radz                 * radius of the particle along x in lb units            *"
      write(outp,'(a)') &
      "WARNING - * txx                  * particle rotational matrix x term along x             *"
      write(outp,'(a)') &
      "WARNING - * txy                  * particle rotational matrix x term along y             *"
      write(outp,'(a)') &
      "WARNING - * txz                  * particle rotational matrix x term along z             *"
      write(outp,'(a)') &
      "WARNING - * tyx                  * particle rotational matrix y term along x             *"
      write(outp,'(a)') &
      "WARNING - * tyy                  * particle rotational matrix y term along y             *"
      write(outp,'(a)') &
      "WARNING - * tyz                  * particle rotational matrix y term along z             *"
      write(outp,'(a)') &
      "WARNING - * tzx                  * particle rotational matrix z term along x             *"
      write(outp,'(a)') &
      "WARNING - * tzy                  * particle rotational matrix z term along y             *"
      write(outp,'(a)') &
      "WARNING - * tzz                  * particle rotational matrix z term along z             *"
      write(outp,'(a)') &
      "WARNING - * qa                   * real part of quaternions                              *"
      write(outp,'(a)') &
      "WARNING - * qb                   * imaginary i part of quaternions                       *"
      write(outp,'(a)') &
      "WARNING - * qc                   * imaginary j part of quaternions                       *"
      write(outp,'(a)') &
      "WARNING - * qd                   * imaginary k part of quaternions                       *"
      write(outp,'(a,/)') &
      "WARNING - ********************************************************************************"
    case (14)
      write (r_char,'(f10.5)')ddata
      write(outp,'(/,a)')"WARNING - densvar should equal or greater than 1."
       write(outp,'(2a,/)')"WARNING - the actual value of densvar is : ", &
       trim(adjustl(r_char))
    case (15)
      write (r_char,'(f10.5)')ddata
      write(outp,'(/,a)')"WARNING - densvar is too low."
       write(outp,'(2a,/)')"WARNING - the actual value of densvar is : ", &
       trim(adjustl(r_char))
    case (16)
      write(outp,'(/,a)')"WARNING - 'read list' is not correctly specified"
      write(outp,'(a,/)')"WARNING - 'read list' contains radx rady radz but particles are spherical"
    case (17)
      write(outp,'(/,2a,/)')"WARNING - rcut (potential cut-off) is not defined in ", &
      "input file"
    case (18)
      write(outp,'(/,2a,/)')"WARNING - delr (width of verlet border) is not defined in ", &
      "input file"
    case (19)
      write(outp,'(/,2a)')"WARNING - particle orientation NOT specified in ", &
      "input file"
      write(outp,'(2a,/)')"an uniform random orientation has been set ", &
      "for each particle"
    case (20)
      write(outp,'(/,a)')"WARNING - 'read list' is not correctly specified"
      write(outp,'(a)')"WARNING - 'read list' both quaternions and rotation matrices"
      write(outp,'(a,/)')"WARNING -  only one of two can be used in input"
    case (21)
      write(outp,'(/,a)')"WARNING - link cell algorithm cannot be used"
      write(outp,'(a,/)')"WARNING - the sub cell size is too small compared to rcut"
    case (22)
      write (r_char,'(i10)')nint(ddata)
      write(outp,'(/,a)')"WARNING - mslistcell is too small"
      write(outp,'(2a,/)')"WARNING - the actual value of mslistcell is : ", &
       trim(adjustl(r_char))
    case (23)
      write (r_char,'(i10)')nint(ddata)
      write(outp,'(/,a)')"WARNING - mxlist is too small"
      write(outp,'(2a,/)')"WARNING - the actual value of mxlist is : ", &
       trim(adjustl(r_char))
    case (24)
      write (r_char,'(i10)')nint(ddata)
      write(outp,'(/,3a,/)')"WARNING - only ",trim(adjustl(r_char)), &
       " field pair types can be specified in input"
    case (25)
      write(outp,'(/,a)')"WARNING - 'read list' is not correctly specified"
      write(outp,'(2a,/)')"WARNING - particle masses not specified neither ", &
       "in input nor xyz files"
    case (26)
      write(outp,'(/,a)')"WARNING - bounce back is automatically switched to halfway mode"
      write(outp,'(a,/)')"WARNING - this is mandatory for the particle treatment"
    case default
      write(outp,'(/,a,i8,/)')"unknown WARNING! code = ",kode
  end select

  return

 end subroutine warning

 end module error_mod


