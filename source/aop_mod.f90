
#include <default_macro.h>
MODULE aop_mod

!***********************************************************************
!     
!     LBsoft module containing definitions of pointer types in order
!     to create structures of arrays
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Bernaschi
!     last modification October 2019
!     
!***********************************************************************
#ifdef INTEL
TYPE REALPTR
   REAL(kind=PRC), contiguous, pointer:: p(:,:,:)
END TYPE REALPTR
#else
TYPE REALPTR
   REAL(kind=PRC), pointer:: p(:,:,:)
END TYPE REALPTR
#endif
TYPE REALPTR_S
   REAL(kind=PRC), pointer:: p
END TYPE REALPTR_S

END MODULE aop_mod
