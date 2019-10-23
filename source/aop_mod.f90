
#include <default_macro.h>
MODULE aop_mod

!***********************************************************************
!     
!     LBsoft module containing definitions of pointer types in order
!     to create structures of arrays
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Bernaschi
!     last modification October 2019
!     
!***********************************************************************

TYPE REALPTR
   REAL(kind=PRC), pointer:: p(:,:,:)
END TYPE REALPTR

TYPE REALPTR_S
   REAL(kind=PRC), pointer:: p
END TYPE REALPTR_S

END MODULE aop_mod
