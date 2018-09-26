
#include <default_macro.h>
MODULE aop_mod

TYPE REALPTR
   REAL(kind=PRC), pointer:: p(:,:,:)
END TYPE REALPTR

TYPE REALPTR_S
   REAL(kind=PRC), pointer:: p
END TYPE REALPTR_S

END MODULE aop_mod
