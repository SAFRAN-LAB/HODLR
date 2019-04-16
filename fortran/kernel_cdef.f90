! C functions declaration
interface

    function create_kernel_c(N) bind(C, name="create_kernel")
        use iso_c_binding
        implicit none
        type(c_ptr) :: create_kernel_c
        integer(c_int), value :: N
    end function

    subroutine delete_kernel_c(K) bind(C, name="delete_kernel")
        use iso_c_binding
        implicit none
        type(c_ptr), value :: K
    end subroutine

    function get_matrix_entry_c(K, i, j) bind(C, name="get_matrix_entry")
        use iso_c_binding
        implicit none
        real(c_double) :: get_matrix_entry_c
        type(c_ptr), intent(in), value :: K
        integer(c_int), value :: i
        integer(c_int), value :: j
    end function

end interface
