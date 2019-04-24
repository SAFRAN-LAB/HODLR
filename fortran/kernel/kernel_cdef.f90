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

    subroutine get_vector_x_c(K, x) bind(C, name="get_vector_x")
        use iso_c_binding
        implicit none
        type(c_ptr), intent(in), value :: K
        real(c_double) :: x(*)
    end subroutine

    ! function get_matrix_c(K, row_start, col_start, row_end, col_end) bind(C, name = "get_matrix")

    !     use iso_c_binding
    !     implicit none

    !     real(c_double), dimension(row_end - row_start, col_end - col_start) :: get_matrix_c
    !     type(c_ptr), intent(in), value :: K

    !     integer(c_int), value :: row_start
    !     integer(c_int), value :: col_start
    !     integer(c_int), value :: row_end
    !     integer(c_int), value :: col_end
    
    ! end function

end interface
