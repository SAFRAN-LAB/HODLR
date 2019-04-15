! C functions declaration
interface

    function createKernel_C(a, b) bind(C, name="createKernel")
        use iso_c_binding
        implicit none
        type(c_ptr) :: createKernel_C
        integer(c_int), value :: a
        integer(c_int), value :: b
    end function

    subroutine deleteKernel_C(K) bind(C, name="deleteKernel")
        use iso_c_binding
        implicit none
        type(c_ptr), value :: K
    end subroutine

    subroutine getMatrix_C(K, i, j, n_rows, n_cols) bind(C, name="getMatrix")
        use iso_c_binding
        implicit none
        type(c_ptr), value :: K
        integer(c_int), value :: i
        integer(c_int), value :: j
        integer(c_int), value :: n_rows
        integer(c_int), value :: n_cols
    end subroutine

end interface
