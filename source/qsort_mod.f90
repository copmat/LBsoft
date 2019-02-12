#include <default_macro.h>

module qsort_mod
 
implicit none
 
type group
    real(kind=PRC) :: f(3), t(3)
    integer :: i,j,k       ! values to be sorted by
end type group

contains
 
      pure function getVal (A,nA, pos)
        implicit none
        integer, intent(in) :: nA,pos
        type (group), dimension(nA), intent(in) :: A
        integer :: i,j,k
        integer :: getVal

        i = A(pos)%i
        j = A(pos)%j
        k = A(pos)%k
        getVal = i + 100*j + 100*100*k
      end function getVal
 

      recursive subroutine QSort(a,na)
        implicit none
        integer, intent(in) :: nA
        type (group), dimension(nA), intent(in out) :: A
        ! LOCAL VARIABLES
        integer :: left, right
        real :: random
        real :: pivot
        type (group) :: temp
        integer :: marker
         
            if (nA > 1) then
         
                call random_number(random)
                pivot = getVal(A,nA, 1 )
                left = 0
                right = nA + 1
         
                do while (left < right)
                    right = right - 1
                    do while (getVal(A,nA,right) > pivot)
                        right = right - 1
                    end do
                    left = left + 1
                    do while (getVal(A,nA,left) < pivot)
                        left = left + 1
                    end do
                    if (left < right) then
                        temp = A(left)
                        A(left) = A(right)
                        A(right) = temp
                    end if
                end do
         
                if (left == right) then
                    marker = left + 1
                else
                    marker = left
                end if
         
                call QSort(A(:marker-1),marker-1)
                call QSort(A(marker:),nA-marker+1)
         
            end if
         
      end subroutine QSort
 
end module qsort_mod
