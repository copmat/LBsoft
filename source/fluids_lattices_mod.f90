
#include <default_macro.h>

 module fluids_lattices_mod
 
!***********************************************************************
!     
!     LBsoft module containing the lattice definition of the LB model
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification May 2020
!     
!***********************************************************************

 
 implicit none
 
 private
 
#if LATTICE==319
 integer, parameter, public :: links=18
  
 character(len=6), parameter, public :: latt_name="d3q19 "
 
 integer, parameter, public :: latt_dim=3
 
 integer, parameter, public :: nbuff=2
 
 real(kind=PRC), parameter, public :: cssq = ( ONE / THREE )
 
 real(kind=PRC), parameter :: p0 = ( ONE / THREE )
 real(kind=PRC), parameter :: p1 = ( ONE / EIGHTEEN )
 real(kind=PRC), parameter :: p2 = ( ONE / THIRTYSIX )
 real(kind=PRC), dimension(0:links), parameter, public :: &
  p = (/p0,p1,p1,p1,p1,p1,p1,p2,p2,p2,p2,p2,p2,p2,p2,p2,p2,p2,p2/)
 
 !real(kind=PRC), parameter :: b_xi = real(2.d0,kind=PRC)  !Liu, Valocchi
 real(kind=PRC), parameter :: b_xi = real(0.5d0,kind=PRC) !Leclaire
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
 
 end module fluids_lattices_mod
