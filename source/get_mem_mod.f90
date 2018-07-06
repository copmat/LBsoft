module check_mem
       use iso_c_binding
!
       implicit none
!
contains
       real(c_double) function get_mem() bind(c,name='get_mem')
       end function get_mem
!
end module check_mem
