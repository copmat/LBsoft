
#include <default_macro.h>
MODULE aop_mod

#if LATTICE==319
 integer, parameter, private :: links=18
#endif

TYPE REALPTR
   REAL(kind=PRC), pointer:: p(:,:,:)
END TYPE REALPTR

END MODULE aop_mod
